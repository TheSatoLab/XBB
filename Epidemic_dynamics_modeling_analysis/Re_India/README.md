# The estimation of the relative Re for lineages related to XBB in India

## Summary:
To estimate the relative Re for lineages related to XBB in India, we analyzed the records for samples from India from June 1, 2022 to November 15, 2022. We removed records with >5% undetermined (N) nucleotide sequences from the dataset. We first simplified the viral lineage classification based on the PANGO lineage. We renamed the sublineages of BA.5 as BA.5, and subsequently, we removed the BA.5 sequences harboring any of the convergent S substitutions, S:R346X, S:K444X, and S:N460X from our dataset in order to exclude the sequences belonging to the recent BA.5 sublineages exhibiting particularly higher Re such as BQ.1.1. Also, we removed the sequences of BA.2.75 harboring any of the convergent S substitutions, S:R346X, S:K444X, S:N460X, and S:F486X. Furthermore, since a part of BA.2.10.1 sequences harbor XBB-characteristic substitutions (S:V83A, S:F486S, and S:F490S) probably due to the misclassification of XBB, we removed the sequences of BA.2.10.1 harboring these XBB-characteristic substitutions. According to the modified viral lineages, we extracted records for viral lineages of interest: BA.2, BA.5, BA.2.75, BM.1, BM1.1, BM.1.1.1, BA.2.10, BJ.1, XBB, and XBB.1. Subsequently, we counted the daily frequency of each viral lineage. Relative Re value for each viral lineage was estimated according to the Bayesian multinomial logistic model, described in our previous study. Briefly, we estimated the logistic slope parameter $\beta_l$ for each viral lineage $l$ using the model and then calculated relative Re for each lineage $r_l$ as $r_l=exp(\gamma\beta_l)$ where $\gamma$ is the average viral generation time (2.1 days) (http://sonorouschocolate.com/covid19/index.php?title=Estimating_Generation_Time_Of_Omicron). Parameter estimation was performed via the MCMC approach implemented in CmdStan v2.30.1 (https://mc-stan.org) with CmdStanr v0.5.3 (https://mc-stan.org/cmdstanr/). Four independent MCMC chains were run with 500 and 1,000 steps in the warmup and sampling iterations, respectively.

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