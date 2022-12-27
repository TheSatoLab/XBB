# The estimation of the global and country-specific Re value of XBB and BQ.1 lineages in the countries where these variants circulated.

## Summary:

To estimate the global average and country-specific Re values for BQ.1 and XBB lineages, we analyzed the sequence records for viral samples collected from August 1, 2022 to November 15, 2022. We defined the sequences of BQ.1 (including its sublineages) harboring S:R346T as BQ.1.1 and the other BQ.1 sequences as BQ.1. Similarly, the sequences of XBB (including its sublineages) harboring S:G252V as XBB and the other XBB sequences as XBB. Subsequently, we extracted the sequence records of BQ.1, BQ.1.1, XBB, and XBB.1 in addition to BA.5 (including its sublineages) and BA.2.75 (including its sublineages), which are predominant lineages before the BQ.1 and XBB emergencies. Next, we counted the daily frequency of the lineages above in each country. We analyzed counties with a total of ≥1000 samples and ≥200 samples of either the BQ.1, BQ.1.1, XBB, or XBB.1 lineages. In this criterion, 11 countries (Australia, Austria, Denmark, India, Indonesia, Israel, Malaysia, Peru, Singapore, the UK, and the USA) remained. To estimate the global average Re values of the lineages above, we employed a hierarchal Bayesian multinomial logistic model, which we established in our previous studies. Briefly, this hierarchal model can estimate the global average and country-specific Re values of lineages of interest simultaneously according to the daily lineage frequency data from multiple countries. The relative Re of each viral lineage $l$ in each county $s$ ($r_{ls}$) was calculated according to the country-specific slope parameter, $\beta_{ls}$, as $r_{ls}=exp(\gamma\beta_{ls})$ where $\gamma$ is the average viral generation time (2.1 days). Similarly, the global average relative Re of each viral lineage was calculated according to the global average slope parameter, $\beta_l$, as $r_l=exp(\gamma\beta_l)$. For parameter estimation, the global average intercept and slope parameters of the BA.5 variant were fixed at 0. Consequently, the relative Re of BA.5 was fixed at 1, and those of the other lineages were estimated relative to that of BA.5. Parameter estimation was performed via the MCMC approach implemented in CmdStan v2.30.1 (https://mc-stan.org) with CmdStanr v0.5.3 (https://mc-stan.org/cmdstanr/). Four independent MCMC chains were run with 500 and 2,000 steps in the warmup and sampling iterations, respectively.

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