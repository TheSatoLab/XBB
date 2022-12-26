#!/usr/bin/env R

library(tidyverse)
library(data.table)
library(ggplot2)
library(rbin)
library(cmdstanr)
library(patchwork)
library(RColorBrewer)
library(maps)


args = commandArgs(trailingOnly=T)

setwd("/Users/jumpeiito/Desktop/analysis/Sato_analysis/SARS-CoV-2_variants/XBB/final_script/estimate_epidemic_freq")


##########args##########
#input

stan_f.name <- 'script/multinomial_independent.stan'
metadata.name <- '../data/metadata.tsv' #args[4]


##########parameters##########
#general
core.num <- 4

#period to be analyzed
date.end <- as.Date("2022-11-15")
date.start <- as.Date("2022-08-01")

#Transmissibility
bin.size <- 1
generation_time <- 2.1
min.count <- 20


#model
multi_nomial_model <- cmdstan_model(stan_f.name)

##########data preprocessing & QC##########
metadata <- fread(metadata.name,header=T,sep="\t",quote="",check.names=T)

metadata.filtered <- metadata %>%
  distinct(Accession.ID,.keep_all=T) %>%
  filter(Host == "Human",
         str_length(Collection.date) == 10,
         Pango.lineage != "",
         Pango.lineage != "None",
         Pango.lineage != "Unassigned",
         !str_detect(Additional.location.information,"[Qq]uarantine")
  )

metadata.filtered <- metadata.filtered %>%
  mutate(Collection.date = as.Date(Collection.date),
         region = str_split(Location," / ",simplify = T)[,1],
         country = str_split(Location," / ",simplify = T)[,2],
         state = str_split(Location," / ",simplify = T)[,3])

metadata.filtered <- metadata.filtered %>% mutate(Pango.lineage.mod = ifelse(str_detect(Pango.lineage,"BQ.1"),"BQ.1",
                                                           ifelse(str_detect(Pango.lineage,"XBB"),"XBB","others"
                                                           )))

metadata.filtered <- metadata.filtered %>% filter(Collection.date >= date.start)

count.country.df <- metadata.filtered %>% group_by(country) %>% summarize(count.country = n())
count.pango.country.df <- metadata.filtered %>% filter(Pango.lineage.mod %in% c("BQ.1","XBB")) %>% group_by(country,Pango.lineage.mod) %>% summarize(count.pango.country = n())

count.pango.country.df.filtered <- count.pango.country.df %>% filter(count.pango.country>=50)
count.country.df.filtered <- count.country.df %>% filter(count.country>=1000)

country.interest.v <- c(count.pango.country.df.filtered$country,count.country.df.filtered$country) %>% unique()

metadata.filtered.interest <- metadata.filtered %>% filter(country %in% country.interest.v)

metadata.filtered.interest <- metadata.filtered.interest %>% mutate(Pango.lineage.mod = factor(Pango.lineage.mod,levels=c("others","BQ.1","XBB")))





#hierarchical

metadata.filtered.interest <- metadata.filtered.interest %>%
  mutate(date.num = as.numeric(Collection.date) - min(as.numeric(Collection.date))  + 1,
        date.bin = cut(date.num,seq(0,max(date.num),bin.size)),
        date.bin.num = as.numeric(date.bin))

metadata.filtered.interest <- metadata.filtered.interest %>% filter(!is.na(date.bin))

out.name <- "metadata.used.lineage_freq.txt"
write.table(metadata.filtered.interest,out.name,col.names=T,row.names=F,sep="\t",quote=F)

metadata.filtered.interest.bin <- metadata.filtered.interest %>% group_by(date.bin.num,Pango.lineage.mod,country) %>% summarize(count = n()) %>% ungroup()


draw.df.theta.latest.sum <- data.frame()

for(country.interest in country.interest.v){
  
  metadata.filtered.interest.bin.interest <- metadata.filtered.interest.bin %>% filter(country == country.interest)
  count.pango.df.interest <- metadata.filtered.interest.bin.interest %>% group_by(Pango.lineage.mod) %>% summarize(count = sum(count))
  
  lineage.interest.v <- count.pango.df.interest %>% filter(count >= min.count) %>% pull(Pango.lineage.mod)
  
  if(length(lineage.interest.v) == 1) next
  
  metadata.filtered.interest.bin.interest <- metadata.filtered.interest.bin.interest %>% filter(Pango.lineage.mod %in% lineage.interest.v)
  
  metadata.filtered.interest.bin.interest.spread <- metadata.filtered.interest.bin.interest %>% spread(key = Pango.lineage.mod, value = count)
  metadata.filtered.interest.bin.interest.spread[is.na(metadata.filtered.interest.bin.interest.spread)] <- 0
  
  
  X <- as.matrix(data.frame(X0 = 1, X1 = metadata.filtered.interest.bin.interest.spread$date.bin.num))
  Y <- metadata.filtered.interest.bin.interest.spread %>% select(- date.bin.num,-country)
  
  count.group <- apply(Y,2,sum)
  count.total <- sum(count.group)
  prop.group <- count.group / count.total
  
  Y <- Y %>% as.matrix()
  
  group.df <- data.frame(group_Id = 1:ncol(Y), group = colnames(Y))
  
  Y_sum.v <- apply(Y,1,sum)
  
  data.stan <- list(K = ncol(Y),
                    D = 2,
                    N = nrow(Y),
                    X = X,
                    Y = Y,
                    generation_time = generation_time,
                    bin_size = bin.size,
                    Y_sum = c(Y_sum.v))
  
  
  fit.stan <- multi_nomial_model$sample(
    data=data.stan,
    iter_sampling=1000,
    iter_warmup=500,
    seed=1234,
    parallel_chains = 4,
    adapt_delta = 0.99,
    max_treedepth = 15,
    #pars=c('b_raw'),
    chains=4)
  
  draw.df.theta <- fit.stan$summary("theta") %>% as.data.frame()
  
  draw.df.theta <- draw.df.theta %>% mutate(date = str_match(variable,'theta\\[([0-9]+),[0-9]+\\]')[,2] %>% as.numeric(),
                                            group_Id = str_match(variable,'theta\\[[0-9]+,([0-9]+)\\]')[,2] %>% as.numeric())
  
  draw.df.theta <- draw.df.theta %>% inner_join(group.df,by="group_Id")
  draw.df.theta <- draw.df.theta %>% mutate(Collection.date = date.start + date - 1)
  
  draw.df.theta.latest <- draw.df.theta %>% filter(Collection.date <= as.Date("2022-11-15")) %>% group_by(group) %>% slice_max(date,n=1) %>% ungroup()
  draw.df.theta.latest <- draw.df.theta.latest %>% mutate(country = country.interest) %>% select(group,Collection.date,country,mean)
  
  draw.df.theta.latest.sum <- rbind(draw.df.theta.latest.sum,draw.df.theta.latest)
  
}


#plot
combination.df <- expand.grid(country.interest.v,c("others","BQ.1","XBB"))
colnames(combination.df) <- c('country','group')

draw.df.theta.latest.sum <- combination.df %>% left_join(draw.df.theta.latest.sum,by=c("country","group"))
draw.df.theta.latest.sum$mean[is.na(draw.df.theta.latest.sum$mean)] <- 0


draw.df.theta.latest.sum <- draw.df.theta.latest.sum %>% 
  mutate(freq.class = ifelse(mean >= 0.5,">50%",
                      ifelse(mean >= 0.2,">20%",
                      ifelse(mean >= 0.1,">10%",
                      ifelse(mean >= 0.01,">1%","<1%")))))

draw.df.theta.latest.sum <- draw.df.theta.latest.sum %>%
    mutate(country = ifelse(country == "United Kingdom","UK",as.character(country)))



world <- map_data("world")

draw.df.theta.latest.BQ.1 <- draw.df.theta.latest.sum %>% filter(group == "BQ.1")
draw.df.theta.latest.XBB <- draw.df.theta.latest.sum %>% filter(group == "XBB")


world.merged.BQ.1 <- world %>% left_join(draw.df.theta.latest.BQ.1 %>% select(region = country, freq.class),by="region")
world.merged.BQ.1 <- world.merged.BQ.1 %>% mutate(freq.class = factor(freq.class,levels=c('>50%','>20%','>10%','>1%','<1%')))

world.merged.XBB <- world %>% left_join(draw.df.theta.latest.XBB %>% select(region = country, freq.class),by="region")
world.merged.XBB <- world.merged.XBB %>% mutate(freq.class = factor(freq.class,levels=c('>50%','>20%','>10%','>1%','<1%')))


col.df <- data.frame(name = c('>50%','>20%','>10%','>1%','<1%'), col = c(rev(brewer.pal(9, "BuPu"))[c(1,3,5,7)],"gray60"))

g1 <- ggplot(world.merged.BQ.1,aes(x = long, y = lat, group = group, fill = freq.class))
g1 <- g1 + geom_polygon(colour = "white", size = 0.1)
g1 <- g1 + scale_fill_manual(values=col.df$col,breaks=col.df$name,na.value="gray75")
g1 <- g1 + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g1 <- g1 + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 strip.text = element_text(size=8))
g1 <- g1 + coord_quickmap() #+ coord_map() #coord_fixed(ratio = 1)

g2 <- ggplot(world.merged.XBB,aes(x = long, y = lat, group = group, fill = freq.class))
g2 <- g2 + geom_polygon(colour = "white", size = 0.1)
g2 <- g2 + scale_fill_manual(values=col.df$col,breaks=col.df$name,na.value="gray75")
g2 <- g2 + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g2 <- g2 + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 strip.text = element_text(size=8))
g2 <- g2 + coord_quickmap() #+ coord_map() #coord_fixed(ratio = 1)


pdf.name <- 'BQ.1_XBB_lineage_freq.pdf'
pdf(pdf.name,width=15,height=10)
plot(g1 / g2)
dev.off()


out.name <- 'BQ.1_XBB_lineage_freq.txt'
write.table(draw.df.theta.latest.sum,out.name,col.names=T,row.names=F,sep="\t",quote=F)

