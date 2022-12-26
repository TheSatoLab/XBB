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

##########args##########
#input

download.date <- "2022-11-21"
download.date <- as.Date(download.date)
metadata.name <- 'metadata.tsv' #args[4]

#output
out.prefix <- '2022_11_21.' #args[5]
pdf.observed.name <- paste(out.prefix,".method1.observed.state.pdf",sep="")
pdf.theta.name <- paste(out.prefix,".method1.theta.state.pdf",sep="")

pdf.growth.rate.name <- paste(out.prefix,".method1.growth_rate.state.pdf",sep="")
txt.growth.rate.name <- paste(out.prefix,".method1.growth_rate.state.txt",sep="")


##########parameters##########
#general
core.num <- 4
variant.ref <- "BA.5"

#period to be analyzed
day.delay <- 0
day.analyzed <- 120

date.end <- as.Date("2022-11-10") #download.date - day.delay
date.start <- as.Date("2022-08-01")


#min numbers
limit.count.analyzed <- 20

#Transmissibility
bin.size <- 1
generation_time <- 2.1


country.interest.v <- c('India','USA','United Kingdom')

#model
multi_nomial_model <- cmdstan_model(stan_f.name)


##########data preprocessing & QC##########
country.info <- read.table(country.info.name,header=T,sep="\t",quote="")
country.info <- country.info %>% select(-region)

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


metadata.filtered.BQ.1_XBB <- metadata.filtered %>% filter(str_detect(Pango.lineage,"BQ") | str_detect(Pango.lineage,"XBB"))


metadata.filtered <- metadata.filtered %>% mutate(Pango.lineage.mod = #ifelse(str_detect(Pango.lineage,"BQ.1") & Accession.ID %in% seq.with_R346T,"BQ.1.1",
                                                                      ifelse(str_detect(Pango.lineage,"BQ.1"),"BQ.1",
                                                                      #ifelse(str_detect(Pango.lineage,"XBB") & Accession.ID %in% seq.with_G252V,"XBB.1",
                                                                      ifelse(str_detect(Pango.lineage,"XBB"),"XBB",as.character(Pango.lineage)
                                                                      ))) #))

metadata.filtered <- metadata.filtered %>% filter(Collection.date >= date.start)

count.country.df <- metadata.filtered %>% group_by(country) %>% summarize(count.country = n())
count.pango.country.df <- metadata.filtered %>% filter(Pango.lineage.mod %in% c("BQ.1","XBB")) %>% group_by(country,Pango.lineage.mod) %>% summarize(count.pango.country = n())

count.pango.country.df.filtered <- count.pango.country.df %>% filter(count.pango.country>=50)
count.country.df.filtered <- count.country.df %>% filter(count.country>=1000)

country.interest.v <- c(count.pango.country.df.filtered$country,count.country.df.filtered$country) %>% unique()

metadata.filtered.interest <- metadata.filtered %>% filter(country %in% country.interest.v)

count.pango.country.df2 <- metadata.filtered.interest %>% group_by(country,Pango.lineage.mod) %>% summarize(count.pango.country = n())
count.pango.country.df2.mean.rank <- count.pango.country.df2 %>% ungroup() %>% group_by(Pango.lineage.mod) %>% summarize(count.pango.country.mean = mean(count.pango.country)) %>% arrange(desc(count.pango.country.mean))

metadata.filtered.interest <- metadata.filtered.interest %>% mutate(Pango.lineage.mod = ifelse(Pango.lineage.mod %in% c("BQ.1","XBB"),as.character(Pango.lineage.mod),"others"))
metadata.filtered.interest <- metadata.filtered.interest %>% mutate(Pango.lineage.mod = factor(Pango.lineage.mod,levels=c("others","BQ.1","XBB")))


#hierarchical

metadata.filtered.interest <- metadata.filtered.interest %>% mutate(date.num = as.numeric(Collection.date) - min(as.numeric(Collection.date))  + 1, date.bin = cut(date.num,seq(0,max(date.num),bin.size)), date.bin.num = as.numeric(date.bin))
metadata.filtered.interest <- metadata.filtered.interest %>% filter(!is.na(date.bin))


metadata.filtered.interest.bin <- metadata.filtered.interest %>% group_by(date.bin.num,Pango.lineage.mod,country) %>% summarize(count = n()) %>% ungroup()

metadata.filtered.interest.bin.spread <- metadata.filtered.interest.bin %>% spread(key=Pango.lineage.mod,value = count)

metadata.filtered.interest.bin.spread[is.na(metadata.filtered.interest.bin.spread)] <- 0
metadata.filtered.interest.bin.spread <- metadata.filtered.interest.bin.spread %>% mutate(country = factor(country))
country <- as.numeric(metadata.filtered.interest.bin.spread$country)

X <- metadata.filtered.interest.bin.spread$date.bin.num

Y <- metadata.filtered.interest.bin.spread %>% select(- date.bin.num,-country)

count.group <- apply(Y,2,sum)
count.total <- sum(count.group)
prop.group <- count.group / count.total

Y <- Y %>% as.matrix()

group.df <- data.frame(group_Id = 1:ncol(Y), group = colnames(Y))
country.df <- data.frame(country_Id = 1:length(levels(metadata.filtered.interest.bin.spread$country)), country = levels(metadata.filtered.interest.bin.spread$country))


Y_sum.v <- apply(Y,1,sum)

stan_f.name <- 'script/multinomial_time_series_hierarchical.stan'
multi_nomial_model <- cmdstan_model(stan_f.name)


data.stan <- list(K = ncol(Y),
                  N = nrow(Y),
                  projection = nrow(Y),
                  D = max(country),
                  X = X,
                  Y = Y,
                  country = country,
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
    max_treedepth = 20,
    #pars=c('b_raw'),
    chains=4)



data_Id.df <- data.frame(data_Id = 1:length(X), country_Id = country, date_Id = X, Y_sum = Y_sum.v, date = as.Date(X,origin=date.start))
data_Id.df <- data_Id.df %>% inner_join(country.df,by="country_Id")

draw.df.theta <- fit.stan$summary("theta") %>% as.data.frame()

draw.df.theta <- draw.df.theta %>% mutate(data_Id = str_match(variable,'theta\\[([0-9]+),[0-9]+\\]')[,2] %>% as.numeric(),
                                                    group_Id = str_match(variable,'theta\\[[0-9]+,([0-9]+)\\]')[,2] %>% as.numeric())

draw.df.theta <- draw.df.theta %>% inner_join(data_Id.df,by="data_Id") %>% inner_join(group.df,by="group_Id")

g <- ggplot(draw.df.theta,aes(x=date,y=mean,color=group))
g <- g + geom_line()
g <- g + facet_wrap(~country)
g

out.name <- "estimated_theta.txt"
write.table(draw.df.theta,out.name,col.names=T,row.names=F,sep="\t",quote=F)



	




draw.df.theta.latest <- draw.df.theta %>% filter(date <= as.Date("2022-11-10")) %>% group_by(country) %>% slice_max(date,n=1) %>% ungroup()
draw.df.theta.latest <- draw.df.theta.latest %>% 
                                           mutate(freq.class = ifelse(mean >= 0.5,">50%",
                                           ifelse(mean >= 0.2,">20%",
                                           ifelse(mean >= 0.1,">10%",
                                           ifelse(mean >= 0.01,">1%","<1%"))))) #)


world <- map_data("world")


draw.df.theta.latest <- draw.df.theta.latest %>% mutate(country = ifelse(country == "United Kingdom","UK",as.character(country)))


draw.df.theta.latest.BQ.1 <- draw.df.theta.latest %>% filter(group == "BQ.1")
draw.df.theta.latest.XBB <- draw.df.theta.latest %>% filter(group == "XBB")


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


pdf.name <- 'BQ.1_XBB_lineage_freq.11-10.advi.pdf'
pdf(pdf.name,width=15,height=10)
plot(g1 / g2)
dev.off()




