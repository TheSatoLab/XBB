# The estimation of the epidemic frequencies of XBB and BQ.1 lineages in each country as of November 15, 2022

## Summary:

To estimate the epidemic frequencies of XBB and BQ.1 lineages in each country as of November 15, 2022, we analyzed the records for viral samples collected from August 1, 2022 to November 15th, 2022. In data for each country, we counted the daily lineage frequency of BQ.1 (including its decedent sublineages), XBB (including its decedent sublineages), and the other SARS-CoV-2 lineages (referred to as “Other lineages”). We analyzed the data only for countries with a total of ≥1000 samples or ≥50 samples of either the BQ.1 or XBB lineages. In this criterion, 56 countries remained. Subsequently, we fitted the multinomial logistic model described in the paragraph above to the daily lineage frequency data of each country separately, and the epidemic frequency of each viral lineage as of November 15, 2022 in each country was estimated. If the data for November 15, 2022 in a particular country are not available, the lineage frequencies at the latest date in the country were used instead. The estimated lineage frequencies for BQ.1 and XBB in each country were shown on the global map using The R library maps v3.4.1 (https://cran.r-project.org/web/packages/maps/index.html).

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



## Usages:
```

```