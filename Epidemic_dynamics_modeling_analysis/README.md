# Epidemic dynamics modeling analysis

## Summary:

We modeled the epidemic dynamics of viral lineages based on the viral genomic surveillance data deposited in the GISAID database (https://www.gisaid.org/). In the present study, we performed three types of analyses: i) The estimation of the relative Re for lineages related to XBB in India, ii) The estimation of the epidemic frequencies of XBB and BQ.1 lineages in each country as of November 15, 2022, and iii) The estimation of the global and country-specific Re value of XBB and BQ.1 lineages in the countries where these variants circulated. For the three analyses, the metadata of viral sequences downloaded from the GISAID database on December 1st, 2022 was used. We excluded the sequence records with the following features: a lack of collection date information; sampling in animals other than humans; sampling by quarantine; or without the PANGO lineage information.

## System requirements:
* **Ubuntu** 20.04.4 LTS
* **R** v4.1.2
* **tidyverse** v1.3.1
* **data.table** v1.14.2
* **ggplot2** v3.3.6
* **rbin** v0.2.0
* **CmdStanr** v0.5.3
* **patchwork** v1.1.2
* **RColorBrewer** v1.1.3

## Note:
Please download the GISAID Metadata from the section "Download packages" (https://gisaid.org/) and save it in the **data** directory.
Subsequently, please run the following command to extract the amino acid information from the metadata.

```
cd <THIS DIRECTORY>

#extract mutation info

python3 script/summarize_mut_info.ver2.py \
  data/metadata.txt \
  data/metadata.mut_long.txt

```