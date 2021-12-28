# Bayesian Instrumental Variable Analysis
This is a Github repo for the [Bayesian instrumental variable analysis paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4314427/) and its non-parametric Bayesian version (unpublished dissertation) written by Dr. Xuyang Lu and Dr. Gang Li.

## About the author

- [Gang Li](https://gangli.faculty.biostat.ucla.edu/)

  Dr. Gang Li is a professor at the [Department of Biostatistics, UCLA](https://www.biostat.ucla.edu/). His research interests include survival analysis, longitudinal data analysis with non-ignorable missing data, high dimensional data analysis, design and analysis of clinical trials, and evaluation and development of biomarkers.  He is currently serving as an Associate Editor of Biometrics and Journal of Nonparametric Statistics. He is a Fellow of the Institute of Mathematical Statistics, a Fellow of the American Statistical Association, an Elected Member of the International Statistics Institute, and an Elected Fellow of the Royal Statistical Society.  

  He collaborates extensively with investigators in biomedical research and public health, especially in cancer research.  He is the Director of UCLA’s Jonsson Comprehensive Cancer Center Biostatistics Shared Resource (BASE Unit). He currently also serves as the Director of Biostatistics for UCLA’s Center for Excellence in Pancreatic Diseases, Co-Director of Biostatistics for an NIH-funded P01 program grant on stem-cell engineering in men with melanoma, and Director of Biostatistics for an NIH-funded P01 program grant on targeting diet-induced promotion of Kras-initiated pancreatic adenocarcinoma.

- [Xuyang Lu]()

  Dr. Xuyang Lu got his PhD degree from the Department of Biostatistics, UCLA under the supervision of Dr. Gang Li. He joined the [Genentech](https://www.gene.com/) after his graduation and stayed there for 7 years (2014-2021). In 2021, Dr. Lu went back to China and joined [Adlai Nortye Biopharma Co., Ltd.](https://en.adlainortye.com/) as a senior director.
  
  His PhD dissertation focuses on building Bayesian models and developing [Markov Chain Monte Carlo (MCMC)](https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo) algorithms for [instrumental variable analysis](https://en.wikipedia.org/wiki/Instrumental_variables_estimation) in causal inference with time-to-event data. It is a poineering work in the related field to deal with time-to-event data with all four types of censoring within the framework of two stage least square estimation.

## Directories and Files

```diff
- *parametric IV model*: this folder contains two files, `IV_MH.R` and `IV_example.R`. They implement the parametric Bayesian IV model proposed in [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4314427/) published in *Statistics in Medicine*.
- *semiparametric IV model _ ARIC example*: this folder corresponds to [section 4.4 of Dr. Lu's PhD dissertation](https://escholarship.org/uc/item/8223z6fp).
- *semiparametric IV model _ simulation*: this folder contains codes to reproduce [section 4.3 of Dr. Lu's PhD dissertation](https://escholarship.org/uc/item/8223z6fp).
```
