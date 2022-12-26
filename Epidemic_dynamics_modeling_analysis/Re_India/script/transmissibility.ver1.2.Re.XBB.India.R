#!/usr/bin/env R

library(tidyverse)
library(data.table)
library(ggplot2)
library(rbin)
library(cmdstanr)
library(patchwork)
library(RColorBrewer)
library(scales)



args = commandArgs(trailingOnly=T)

setwd("/Users/jumpeiito/Desktop/analysis/Sato_analysis/SARS-CoV-2_variants/XBB/final_script/Re_India")


##########inputs##########
stan_f.name <- 'script//multinomial_independent.stan'
metadata.name <- '../data/metadata.tsv' #args[4]
mut.info.name <- '../data/metadata.mut_long.tsv'

pango.info.name <- "../data/lineage_report.2022-11-24.txt"
pango.info <- read.table(pango.info.name,header=T)


##########parameters##########
#general
core.num <- 4

date.end <- as.Date("2022-11-15")
date.start <- as.Date("2022-06-01")


#Transmissibility
bin.size <- 1
generation_time <- 2.1


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

metadata.filtered <- metadata.filtered %>% filter(Collection.date >= date.start, country == "India")

count.pango.df <- metadata.filtered %>% group_by(Pango.lineage) %>% summarize(count = n()) %>% arrange(desc(count))



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
mut.interest.XBB.v <- c("Spike_V83A","Spike_F486S", "Spike_F490S")
mut.interest.BA.2.75.v <- c("Spike_R346","Spike_F486","Spike_K444","Spike_L452")


mut.info.convergent_mut <- mut.info %>% filter(mut.mod %in% mut.interest.BA.5.v)
mut.info.omit.BA.2.10.1 <- mut.info %>% filter(mut %in% mut.interest.XBB.v)
mut.info.omit.BA.2.75 <- mut.info %>% filter(mut.mod %in% mut.interest.BA.2.75.v)



mut.info.R346T <- mut.info %>% filter(mut == "Spike_R346T")
mut.info.G252V <- mut.info %>% filter(mut == "Spike_G252V")

metadata.filtered <- metadata.filtered %>% mutate(Pango.lineage.mod = ifelse(str_detect(Pango.lineage,"BQ.1") & Accession.ID %in% as.character(mut.info.R346T$Id),"BQ.1.1",
                                                                      ifelse(str_detect(Pango.lineage,"BQ.1"),"BQ.1",
                                                                      ifelse(str_detect(Pango.lineage,"XBB") & Accession.ID %in% as.character(mut.info.G252V$Id),"XBB.1",
                                                                      ifelse(str_detect(Pango.lineage,"XBB"),"XBB",
                                                                      ifelse(Pango.lineage.class == "BA.5" & (!Accession.ID %in% as.character(mut.info.convergent_mut$Id)),"BA.5_wo_mut",as.character(Pango.lineage)
                                                                      ))))))



metadata.filtered <- metadata.filtered %>% filter(Is.complete. == "TRUE", !N.Content >0.05 | is.na(N.Content))
metadata.filtered <- metadata.filtered %>% filter(!(Pango.lineage.mod == "BA.2.10.1" & Accession.ID %in% mut.info.omit.BA.2.10.1$Id))
metadata.filtered <- metadata.filtered %>% filter(!(Pango.lineage.mod == "BA.2.75" & Accession.ID %in% mut.info.omit.BA.2.75$Id))

lineage.interest.v <- c("BA.2","BA.2.10","BA.5_wo_mut","BA.2.75","BM.1","BM.1.1","BM.1.1.1","BJ.1","XBB","XBB.1")

count.pango.df <- metadata.filtered %>% group_by(Pango.lineage.mod) %>% summarize(count = n()) %>% arrange(desc(count))
count.pango.df %>% filter(Pango.lineage.mod %in% lineage.interest.v)


metadata.filtered.interest <- metadata.filtered %>% filter(Pango.lineage.mod %in% lineage.interest.v)


metadata.filtered.interest <- metadata.filtered.interest %>% mutate(date.num = as.numeric(Collection.date) - min(as.numeric(Collection.date))  + 1, date.bin = cut(date.num,seq(0,max(date.num),bin.size)), date.bin.num = as.numeric(date.bin))
metadata.filtered.interest <- metadata.filtered.interest %>% filter(!is.na(date.bin))

out.name <- 'metadata.epidemic_freq.India.txt'
write.table(metadata.filtered.interest,out.name,col.names=T,row.names=F,sep="\t",quote=F)



metadata.filtered.interest.bin <- metadata.filtered.interest %>% group_by(date.bin.num,Pango.lineage.mod) %>% summarize(count = n()) %>% ungroup()

metadata.filtered.interest.bin.spread <- metadata.filtered.interest.bin %>% spread(key=Pango.lineage.mod,value = count)

metadata.filtered.interest.bin.spread[is.na(metadata.filtered.interest.bin.spread)] <- 0
metadata.filtered.interest.bin.spread <- metadata.filtered.interest.bin.spread



X <- as.matrix(data.frame(X0 = 1, X1 = metadata.filtered.interest.bin.spread$date.bin.num))
Y <- metadata.filtered.interest.bin.spread %>% select(- date.bin.num)


count.group <- apply(Y,2,sum)
count.total <- sum(count.group)
prop.group <- count.group / count.total

Y <- Y %>% as.matrix()
apply(Y,2,sum)


group.df <- data.frame(group_Id = 1:ncol(Y), group = colnames(Y))

Y_sum.v <- apply(Y,1,sum)


#########fitting##########
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
    chains=4)



##########outputs##########

stat.info <- fit.stan$summary("growth_rate") %>% as.data.frame()
stat.info$Pango.lineage <- colnames(Y)[2:ncol(Y)]

stat.info.q <- fit.stan$summary("growth_rate", ~quantile(.x, probs = c(0.025,0.975))) %>% as.data.frame() %>% rename(q2.5 = `2.5%`, q97.5 = `97.5%`)
stat.info <- stat.info %>% inner_join(stat.info.q,by="variable")

out.name <- 'growth_rate.India.2022-12-01.txt'
write.table(stat.info,out.name,col.names=T,row.names=F,sep="\t",quote=F)

stat.info.plot <- stat.info %>% select(mean,Pango.lineage)
stat.info.plot <- rbind(data.frame(mean=1,Pango.lineage="BA.2"),stat.info.plot)

stat.info.plot <-stat.info.plot %>% mutate(x=c(1,2,2,2,5,3,4,5,6,7),
                                           y=c(2,2,3,1,3,3,3,2,2,2)
                                           )


g2 <- ggplot(stat.info.plot,aes(x=x,y=y,fill=mean))
g2 <- g2 + geom_point(size=10, color = "black", shape = 21)
     
g2 <- g2 + scale_fill_gradientn(colors=c("white","gold","red"),limits=c(1,1.55),oob=squish)
g2 <- g2 + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g2 <- g2 + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               strip.text = element_text(size=8))
g2 <- g2 + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g2 <- g2 + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               strip.text = element_text(size=8))
g2 <- g2 + scale_y_continuous(limits=c(0.5,3.5))

pdf.name <- 'growth_rate.India.2022-12-01.dot.wyr.pdf'
pdf(pdf.name,width=5,height=2.5)
plot(g2)
dev.off()

