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
mut.info.name <- 'metadata.mut_long.tsv'


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

metadata.filtered <- metadata.filtered %>% filter(Collection.date >= date.start)

#mutation data
mut.info <- fread(mut.info.name,header=T,sep="\t",check.names=T)
mut.info <- mut.info %>% filter(Id %in% as.character(metadata.filtered$Accession.ID))

mut.info.R346T <- mut.info %>% filter(mut == "Spike_R346T")
mut.info.G252V <- mut.info %>% filter(mut == "Spike_G252V")


metadata.filtered <- metadata.filtered %>% mutate(Pango.lineage.mod = ifelse(str_detect(Pango.lineage,"BQ.1") & Accession.ID %in% as.character(mut.info.R346T$Id),"BQ.1.1",
                                                                      ifelse(str_detect(Pango.lineage,"BQ.1"),"BQ.1",
                                                                      ifelse(str_detect(Pango.lineage,"XBB") & Accession.ID %in% as.character(mut.info.G252V$Id),"XBB.1",
                                                                      ifelse(str_detect(Pango.lineage,"XBB"),"XBB",as.character(Pango.lineage)
                                                                      )))))


count.country.df <- metadata.filtered %>% group_by(country) %>% summarize(count.country = n())
count.pango.country.df <- metadata.filtered %>% filter(Pango.lineage.mod %in% c("BQ.1","BQ.1.1","XBB","XBB.1")) %>% group_by(country,Pango.lineage.mod) %>% summarize(count.pango.country = n())

count.pango.country.df.filtered <- count.pango.country.df %>% filter(count.pango.country>=200)
count.country.df.filtered <- count.country.df %>% filter(count.country>=5000)

country.interest1.v <- unique(count.pango.country.df.filtered$country)
country.interest2.v <- unique(count.country.df.filtered$country)

country.interest.v <- intersect(country.interest1.v,country.interest2.v)

metadata.filtered.interest <- metadata.filtered %>% filter(country %in% country.interest.v)


count.pango.country.df2 <- metadata.filtered.interest %>% group_by(country,Pango.lineage.mod) %>% summarize(count.pango.country = n())
count.pango.country.df2.mean.rank <- count.pango.country.df2 %>% ungroup() %>% group_by(Pango.lineage.mod) %>% summarize(count.pango.country.mean = mean(count.pango.country)) %>% arrange(desc(count.pango.country.mean))

pango.interest.v <- count.pango.country.df2.mean.rank %>% slice_max(count.pango.country.mean,n=5) %>% pull(Pango.lineage.mod)

pango.interest.v <- c(pango.interest.v,c("BQ.1","BQ.1.1","XBB","XBB.1")) %>% unique()

metadata.filtered.interest <- metadata.filtered.interest %>% filter(Pango.lineage.mod %in% pango.interest.v)
metadata.filtered.interest <- metadata.filtered.interest %>% mutate(Pango.lineage.mod = factor(Pango.lineage.mod,levels=pango.interest.v))


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
    max_treedepth = 15,
    #pars=c('b_raw'),
    chains=4)



data_Id.df <- data.frame(data_Id = 1:length(X), country_Id = country, date_Id = X, Y_sum = Y_sum.v, date = as.Date(X,origin=date.start))
data_Id.df <- data_Id.df %>% inner_join(country.df,by="country_Id")

draw.df.theta <- fit.stan$summary("theta") %>% as.data.frame()

draw.df.theta <- draw.df.theta %>% mutate(data_Id = str_match(variable,'theta\\[([0-9]+),[0-9]+\\]')[,2] %>% as.numeric(),
                                                    group_Id = str_match(variable,'theta\\[[0-9]+,([0-9]+)\\]')[,2] %>% as.numeric())

draw.df.theta <- draw.df.theta %>% inner_join(data_Id.df,by="data_Id") %>% inner_join(group.df,by="group_Id")

