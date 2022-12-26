#!/usr/bin/env R

library(tidyverse)
library(data.table)
library(ggplot2)
library(rbin)
library(cmdstanr)
library(patchwork)
library(RColorBrewer)


args = commandArgs(trailingOnly=T)

setwd("/Users/jumpeiito/Desktop/analysis/Sato_analysis/SARS-CoV-2_variants/XBB/final_script/hierarchical_Re")


##########inputs##########
stan_f.name <- 'script/multinomial_time_series_hierarchical.stan'
metadata.name <- '../data/metadata.tsv'
mut.info.name <- '../data/metadata.mut_long.tsv'

pango.info.name <- "../data/lineage_report.2022-11-24.txt"
pango.info <- read.table(pango.info.name,header=T)

#model
multi_nomial_model <- cmdstan_model(stan_f.name)


##########parameters##########
#general
core.num <- 4
variant.ref <- "BA.5"

#period to be analyzed
date.end <- as.Date("2022-11-15")
date.start <- as.Date("2022-08-01")

#Transmissibility
bin.size <- 1
generation_time <- 2.1


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

metadata.filtered <- metadata.filtered %>% filter(Collection.date >= date.start)

metadata.filtered <- metadata.filtered %>% left_join(pango.info,by="Pango.lineage")
metadata.filtered <- metadata.filtered %>% mutate(Pango.lineage.class = ifelse(str_detect(Pango.lineage.original,"B.1.1.529.1"),"BA.1",
                               ifelse(str_detect(Pango.lineage.original,"B.1.1.529.2.75"),"BA.2.75",
                               ifelse(str_detect(Pango.lineage.original,"B.1.1.529.2"),"BA.2",
                               ifelse(str_detect(Pango.lineage.original,"B.1.1.529.3"),"BA.3",
                               ifelse(str_detect(Pango.lineage.original,"B.1.1.529.4"),"BA.4",
                               ifelse(str_detect(Pango.lineage.original,"B.1.1.529.5"),"BA.5","others")))))))


#mutation data
mut.info <- fread(mut.info.name,header=T,sep="\t",check.names=T)
mut.info <- mut.info %>% filter(Id %in% as.character(metadata.filtered$Accession.ID))
mut.info <- mut.info %>% mutate(mut.mod = gsub("[A-Za-z]+$","",mut))


mut.interest.BA.5.v <- c("Spike_R346","Spike_K444","Spike_N460")

mut.info.convergent_mut <- mut.info %>% filter(mut.mod %in% mut.interest.BA.5.v)

mut.info.R346T <- mut.info %>% filter(mut == "Spike_R346T")
mut.info.G252V <- mut.info %>% filter(mut == "Spike_G252V")

metadata.filtered <- metadata.filtered %>% mutate(Pango.lineage.mod = ifelse(str_detect(Pango.lineage,"BQ.1") & Accession.ID %in% as.character(mut.info.R346T$Id),"BQ.1.1",
                                                                      ifelse(str_detect(Pango.lineage,"BQ.1"),"BQ.1",
                                                                      ifelse(str_detect(Pango.lineage,"XBB") & Accession.ID %in% as.character(mut.info.G252V$Id),"XBB.1",
                                                                      ifelse(str_detect(Pango.lineage,"XBB"),"XBB",
                                                                      ifelse(Pango.lineage.class == "BA.5" & (!Accession.ID %in% as.character(mut.info.convergent_mut$Id)),"BA.5_wo_mut",as.character(Pango.lineage)
                                                                      ))))))


count.country.df <- metadata.filtered %>% group_by(country) %>% summarize(count.country = n())
count.pango.country.df <- metadata.filtered %>% filter(Pango.lineage.mod %in% c("XBB","XBB.1")) %>% group_by(country,Pango.lineage.mod) %>% summarize(count.pango.country = n())

count.pango.country.df.filtered <- count.pango.country.df %>% filter(count.pango.country>=100)
count.country.df.filtered <- count.country.df %>% filter(count.country>=2000)

country.interest1.v <- unique(count.pango.country.df.filtered$country)
country.interest2.v <- unique(count.country.df.filtered$country)

country.interest.v <- intersect(country.interest1.v,country.interest2.v)

metadata.filtered.interest <- metadata.filtered %>% filter(country %in% country.interest.v)


pango.interest.v <- c("BA.5_wo_mut","BA.2.75","BQ.1","BQ.1.1","XBB","XBB.1")

metadata.filtered.interest <- metadata.filtered.interest %>% filter(Pango.lineage.mod %in% pango.interest.v)
metadata.filtered.interest <- metadata.filtered.interest %>% mutate(Pango.lineage.mod = factor(Pango.lineage.mod,levels=pango.interest.v))


metadata.filtered.interest <- metadata.filtered.interest %>% mutate(date.num = as.numeric(Collection.date) - min(as.numeric(Collection.date))  + 1, date.bin = cut(date.num,seq(0,max(date.num),bin.size)), date.bin.num = as.numeric(date.bin))
metadata.filtered.interest <- metadata.filtered.interest %>% filter(!is.na(date.bin))

out.name <- "metadata.used_for_analysis.hierarchical.txt"
write.table(metadata.filtered.interest,out.name,col.names=T,row.names=F,sep="\t",quote=F)

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
country.df <- data.frame(country_Id = 1:length(levels(metadata.filtered.interest.bin.spread$country)),
                         country = levels(metadata.filtered.interest.bin.spread$country))


Y_sum.v <- apply(Y,1,sum)


##########fitting##########
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
    iter_sampling=2000,
    iter_warmup=1000,
    seed=1234,
    parallel_chains = 4,
    adapt_delta = 0.99,
    max_treedepth = 15,
    #pars=c('b_raw'),
    chains=4)


##########outputs##########
#growth rate mean
stat.info.mean <- fit.stan$summary("growth_rate_mean") %>% as.data.frame()
stat.info.mean.q <- fit.stan$summary("growth_rate_mean", ~quantile(.x, probs = c(0.025,0.975))) %>% as.data.frame()
stat.info.mean <- merge(stat.info.mean,stat.info.mean.q,by="variable")
stat.info.mean <- stat.info.mean %>% mutate(group_Id = str_match(variable,'growth_rate_mean\\[([0-9]+)\\]')[,2] %>% as.numeric() + 1)

stat.info.mean.merged <- merge(stat.info.mean,group.df,by="group_Id") %>% select(-group_Id,-variable)
stat.info.mean.merged <- stat.info.mean.merged %>% arrange(desc(mean))


#growth_rate each country
stat.info.each <- fit.stan$summary("growth_rate") %>% as.data.frame()
stat.info.each.q <- fit.stan$summary("growth_rate", ~quantile(.x, probs = c(0.025,0.975))) %>% as.data.frame()
stat.info.each <- merge(stat.info.each,stat.info.each.q,by="variable")

stat.info.each <- stat.info.each %>% mutate(country_Id = str_match(variable,'growth_rate\\[([0-9]+),[0-9]+\\]')[,2] %>% as.numeric(), group_Id = str_match(variable,'growth_rate\\[[0-9]+,([0-9]+)\\]')[,2] %>% as.numeric() + 1)
stat.info.each.merged <- stat.info.each %>% inner_join(group.df,by="group_Id") %>% inner_join(country.df,by="country_Id") %>% select(-group_Id,-country_Id,-variable)


out.name <- paste('growth_rate.mean.hierarchical.txt',sep="")
write.table(stat.info.mean.merged,out.name,col.names=T,row.names=F,sep="\t",quote=F)

out.name <- paste('growth_rate.each_country.hierarchical.txt',sep="")
write.table(stat.info.each.merged,out.name,col.names=T,row.names=F,sep="\t",quote=F)


###growth rate_mean
#growth rate
draw.df.growth_rate <- fit.stan$draws("growth_rate_mean", format = "df") %>% as.data.frame() %>% select(! contains('.'))
draw.df.growth_rate.long <- draw.df.growth_rate %>% gather(key = class, value = value)

draw.df.growth_rate.long <- draw.df.growth_rate.long %>% mutate(group_Id = str_match(draw.df.growth_rate.long$class,'growth_rate_mean\\[([0-9]+)\\]')[,2] %>% as.numeric() + 1)
draw.df.growth_rate.long <- merge(draw.df.growth_rate.long,group.df,by="group_Id") %>% select(value,group)
draw.df.growth_rate.long <- draw.df.growth_rate.long %>% group_by(group) %>% filter(value>=quantile(value,0.005),value<=quantile(value,0.995))
draw.df.growth_rate.long <- rbind(data.frame(group="BA.5_wo_mut",value=1),draw.df.growth_rate.long)
draw.df.growth_rate.long <- draw.df.growth_rate.long %>% mutate(group = factor(group,levels=pango.interest.v))

col.v <- brewer.pal(ncol(Y)+1, "Set1")[c(1:5,7)]
draw.df.growth_rate.long <- draw.df.growth_rate.long
g <- ggplot(draw.df.growth_rate.long,aes(x=group,y=value,color=group,fill=group))
g <- g + geom_hline(yintercept=1, linetype="dashed", alpha=0.5)
g <- g + geom_violin(alpha=0.6,scale="width")
g <- g + stat_summary(geom="pointrange",fun = median, fun.min = function(x) quantile(x,0.025), fun.max = function(x) quantile(x,0.975), size=0.5,fatten =1.5)
g <- g + scale_color_manual(values=col.v)
g <- g + scale_fill_manual(values=col.v)
g <- g + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               strip.text = element_text(size=8))
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
g <- g + xlab('') + ylab('Relative Re (/BA.5)')
g <- g + theme(legend.position = 'none')
g <- g + scale_y_continuous(limits=c(1,1.42),breaks=c(1,1.2,1.4))
g

pdf.name <- 'growth_rate.global_mean.pdf'

pdf(pdf.name,width=1.8,height=3.5)
plot(g)
dev.off()



###growth rate for each country
draw.df.growth_rate.each <- fit.stan$draws("growth_rate", format = "df") %>% as.data.frame() %>% select(! contains('.'))
draw.df.growth_rate.each.long <- draw.df.growth_rate.each %>% gather(key = class, value = value)

draw.df.growth_rate.each.long <- draw.df.growth_rate.each.long %>% mutate(country_Id = str_match(class,'growth_rate\\[([0-9]+),[0-9]+\\]')[,2] %>% as.numeric(), group_Id = str_match(class,'growth_rate\\[[0-9]+,([0-9]+)\\]')[,2] %>% as.numeric() + 1)
draw.df.growth_rate.each.long <- draw.df.growth_rate.each.long %>% inner_join(group.df,by="group_Id") %>% inner_join(country.df,by="country_Id") %>% select(-group_Id,-country_Id,-class)
draw.df.growth_rate.each.long <- rbind(data.frame(value=1,group="BA.5_wo_mut",country=country.df$country),draw.df.growth_rate.each.long)

draw.df.growth_rate.each.long <- draw.df.growth_rate.each.long
draw.df.growth_rate.each.long <- draw.df.growth_rate.each.long %>% group_by(group,country) %>% filter(value>=quantile(value,0.005),value<=quantile(value,0.995))

draw.df.growth_rate.each.long <- draw.df.growth_rate.each.long %>% mutate(group = factor(group,levels=pango.interest.v))



g <- ggplot(draw.df.growth_rate.each.long,aes(x=group,y=value,color=group,fill=group))
g <- g + geom_hline(yintercept=1, linetype="dashed", alpha=0.5)
g <- g + geom_violin(alpha=0.4,scale="width")
g <- g + stat_summary(geom="pointrange",fun = median, fun.min = function(x) quantile(x,0.025), fun.max = function(x) quantile(x,0.975), size=0.5,fatten =1)
g <- g + scale_color_manual(values=col.v)
g <- g + scale_fill_manual(values=col.v)
g <- g + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               strip.text = element_text(size=8))
g <- g + xlab('') + ylab('Relative Re (/BA.5)')
g <- g + theme(legend.position = 'none')
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

g <- g + facet_wrap(~country,ncol=5)
g

pdf.name <- 'reproduction_number.each_country.pdf'

pdf(pdf.name,width=6,height=5)
plot(g)
dev.off()





#theta
data_Id.df <- data.frame(data_Id = 1:length(X), country_Id = country, date_Id = X, Y_sum = Y_sum.v, date = as.Date(X,origin=date.start))


data.freq <- metadata.filtered.interest.bin %>% rename(group = Pango.lineage.mod) %>% filter(group != "others") %>% group_by(country,date.bin.num) %>% mutate(freq = count / sum(count))
data.freq <- merge(data.freq,data_Id.df,by.x="date.bin.num",by.y="date_Id")


draw.df.theta <- fit.stan$draws("theta", format = "df") %>% as.data.frame() %>% select(! contains('.'))
draw.df.theta.long <- draw.df.theta %>% gather(key = class, value = value)
draw.df.theta.long <- draw.df.theta.long %>% mutate(data_Id = str_match(class,'theta\\[([0-9]+),[0-9]+\\]')[,2] %>% as.numeric(),
                                                            group_Id = str_match(class,'theta\\[[0-9]+,([0-9]+)\\]')[,2] %>% as.numeric())


draw.df.theta.long <- draw.df.theta.long %>% inner_join(data_Id.df,by="data_Id")
draw.df.theta.long.sum <- draw.df.theta.long %>% group_by(group_Id, country_Id, date) %>% summarize(mean = mean(value),ymin = quantile(value,0.025),ymax = quantile(value,0.975))

draw.df.theta.long.sum <- draw.df.theta.long.sum %>% inner_join(group.df,by="group_Id") %>% inner_join(country.df,by="country_Id")
draw.df.theta.long.sum <- draw.df.theta.long.sum %>% inner_join(data.freq %>% select(group,count,freq,date,country),by=c("date","group","country"))

draw.df.theta.long.sum.filtered <- draw.df.theta.long.sum %>% mutate(group = factor(group,levels=pango.interest.v))

g <- ggplot(draw.df.theta.long.sum.filtered,aes(x=date, y = mean, fill=group, color = group))
g <- g + geom_ribbon(aes(ymin=ymin,ymax=ymax), color=NA,alpha=0.4)
g <- g + geom_line(size=0.3)
g <- g + scale_x_date(date_labels = "%y-%m", date_breaks = "1 months", date_minor_breaks = "1 month")
g <- g + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g <- g + theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text = element_text(size=8)
    )
g <- g + facet_wrap(~country,nrow=2)
g <- g + scale_color_manual(values = col.v)
g <- g + scale_fill_manual(values = col.v)
g <- g + scale_size_continuous(range = c(0.2, 4))
g


pdf.name <- 'theta.each_country.pdf'
pdf(pdf.name,width=10,height=4)
plot(g)
dev.off()


