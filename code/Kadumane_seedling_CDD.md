Density-dependence seedling mortality in Kadumane
================
Robert Bagchi and Ashwin Viswanathan
16 July 2024

- [set up](#set-up)
- [Data](#data)
- [Summary information for the
  paper](#summary-information-for-the-paper)
- [Models](#models)
  - [Scaled conspecific density
    models](#scaled-conspecific-density-models)
    - [Diagnostics](#diagnostics)
    - [Model inference](#model-inference)
    - [Take-homes:](#take-homes)
    - [Models split by categorical fragment
      size](#models-split-by-categorical-fragment-size)
    - [Take-homes](#take-homes-1)
    - [Graphics](#graphics)
    - [Plotting effects of fragment
      area](#plotting-effects-of-fragment-area)
  - [Species specific inferences.](#species-specific-inferences)
- [Session Information](#session-information)

## set up

``` r
library(tidyverse)
library(ggthemes)
library(knitr)
library(glmmTMB)
library(DHARMa)
library(broom.mixed)
library(ggeffects)
library(ggdist)
library(sjPlot)
library(patchwork)
```

``` r
theme_set(theme_tufte())
```

# Data

Load data

``` r
sp_codes <- read_csv("data/sp_codes.csv")
```

    ## Rows: 26 Columns: 4
    ## -- Column specification --------------------------------------------------------
    ## Delimiter: ","
    ## chr (4): code, genus, species, family
    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
site_dat <- read_rds("data/kadumane_site_metadata.rds")
sdls <- read_rds("data/kadumane_seedlings.rds")
plot_dat <- read_rds("data/kadumane_plot_metadata.rds")
```

1.  Join in site and plot meta-data.

2.  Remove rows for species with no seedlings at start of the census in
    plot.

3.  Select species that

- Are recorded in at least 5 plots.
- Vary in density among plots.
- Are groups of many (unidentified) species with the same code.

``` r
sdls |> summarise(n_sdls = sum(census.start), 
                  n_survs = sum(census.final),
                  n_species = n_distinct(species)) |> knitr::kable()
```

| n_sdls | n_survs | n_species |
|-------:|--------:|----------:|
|   7134 |    3054 |       102 |

``` r
sdls <-  sdls |> left_join(site_dat) |> left_join(plot_dat) |> 
  filter(census.start > 0)
```

    ## Joining with `by = join_by(site)`
    ## Joining with `by = join_by(site, location, group, plot)`

``` r
## find species with enough individuals and variation in density to allow 
## analysis
sp_list <- group_by(sdls, species) |> 
  summarise(abund = sum(census.start) ,  ## total abund
            n = sum(census.start > 0), ## number of plots withe species
            sd_dens = sd(census.start[census.start > 0])) |> ##var in density
  filter(n > 1) |> arrange(n)
## only lose 8 species by restricting to 5 or more occurrences (instead of 1)
## Seems reasonable. Also removing species that were unreliably identified.

sp_list <- filter(sp_list, n > 4, !(species %in% c("Palm", "Artoh", "SC")),
                  sd_dens > 0) 
sdls <- filter(sdls, species %in% sp_list$species)
dim(sdls) ## 968 columns
```

    ## [1] 968  19

``` r
# Number of plots= 37 locs * 3 groups * 5 plots
n_plots <- 37*3*5
n_sites <- 21
```

Rename columns, add total seedling density, scale and centre data.

``` r
## shorten names
sdls <- rename(sdls, 
               "trt_F" = "treatment.fungicide", 
               "trt_I" = "treatment.insecticide", 
               "Pr_m" = "proportion.mortality",
               "gr" = "group",## group causes problems with some helper funcs
               "loc" = "location") 
## categorical variable for treatment
sdls <- mutate(sdls, trt = case_when(
  trt_F == "0" & trt_I == "0" ~ "C",
  trt_F == "F" & trt_I == "0" ~ "F",
  trt_F == "0" & trt_I == "I" ~ "I",
  trt_F == "F" & trt_I == "I" ~ "FI"),
  trt = factor(trt, levels = c("C", "I", "F", "FI")))

## add total density
tot_dens <- sdls |> group_by(site, loc, gr, plot) |> 
  summarise(tot_dens = sum(census.start))
```

    ## `summarise()` has grouped output by 'site', 'loc', 'gr'. You can override using
    ## the `.groups` argument.

``` r
## add species mean density
## divide by total number of plots
sp_mean_dens <- sdls |> group_by(species) |> 
  summarise(sp_mean_condens = sum(census.start)/n_plots, 
            sp_mean_surv = sum(census.final)/sum(census.start))

sdls <- left_join(sdls, tot_dens, by = c("site", "loc", "gr", "plot")) |> 
  left_join(select(sp_mean_dens, - sp_mean_surv))
```

    ## Joining with `by = join_by(species)`

``` r
#scale density by mean, fix couple of NAs
sdls <- mutate(sdls, 
               slope.degrees = replace_na(slope.degrees, 5), 
               Pr_s = 1 - Pr_m,
               con_dens = census.start,
               con_dens_s = con_dens/sp_mean_condens,
               slope.degrees_s = as.vector(scale(slope.degrees)),
                trt_F = factor(trt_F, labels = c("0", "F")),
                trt_I = factor(trt_I, labels = c("0", "I"))
               )

## calculate total scaled density. Note this is the sum of the scaled densities
## of all species in the plot. Remember that the total_density contrast only 
## works when all conspecific densities sum to total density (replicating 
## a sum-to-zero contrast for conspecific and heterospecific densities.
sdls <- left_join(sdls, 
                  group_by(sdls, site, loc, gr, plot) |> 
                    summarise(tot_dens_s = sum(con_dens_s)))
```

    ## `summarise()` has grouped output by 'site', 'loc', 'gr'. You can override using
    ## the `.groups` argument.
    ## Joining with `by = join_by(site, loc, gr, plot)`

``` r
summary(sdls)
```

    ##       site     loc       gr      plot        species     census.start    
    ##  S3     :119   L1:549   G1:283   1:185   SR      :172   Min.   :  1.000  
    ##  S1     :109   L2:225   G2:351   2:210   Symp    :122   1st Qu.:  1.000  
    ##  S12    : 77   L3:161   G3:334   5:183   Cinam   : 86   Median :  1.000  
    ##  S13    : 75   L4: 33            6:199   Climber1: 77   Mean   :  6.413  
    ##  S39    : 74                     7:191   FLC     : 64   3rd Qu.:  3.000  
    ##  S5     : 63                             Litsea  : 61   Max.   :300.000  
    ##  (Other):451                             (Other) :386                    
    ##    census.mid       census.final         Pr_m        proportion.connectivity
    ##  Min.   :  0.000   Min.   : 0.000   Min.   :0.0000   Min.   :0.0000         
    ##  1st Qu.:  0.000   1st Qu.: 0.000   1st Qu.:0.0000   1st Qu.:0.0100         
    ##  Median :  1.000   Median : 1.000   Median :0.1667   Median :0.0200         
    ##  Mean   :  3.089   Mean   : 2.769   Mean   :0.3810   Mean   :0.1255         
    ##  3rd Qu.:  2.000   3rd Qu.: 2.000   3rd Qu.:1.0000   3rd Qu.:0.1500         
    ##  Max.   :109.000   Max.   :90.000   Max.   :1.0000   Max.   :0.6700         
    ##                                                                             
    ##  fragment.size    size.category      location.code   group.code 
    ##  Min.   :  1.10   Length:968         S14L1  : 54   S14L1G2: 23  
    ##  1st Qu.:  9.00   Class :character   S1L2   : 49   S12L1G2: 18  
    ##  Median : 46.00   Mode  :character   S12L2  : 40   S14L1G3: 17  
    ##  Mean   : 51.15                      S13L3  : 40   S1L1G3 : 17  
    ##  3rd Qu.: 64.50                      S12L1  : 37   S1L2G1 : 17  
    ##  Max.   :149.00                      S3L3   : 37   S1L2G2 : 17  
    ##                                      (Other):711   (Other):859  
    ##  proportion.rock  slope.degrees    treatment.control trt_F   trt_I     trt     
    ##  Min.   :0.0000   Min.   : 0.000   Min.   :0.0000    0:594   0:578   C   :395  
    ##  1st Qu.:0.0500   1st Qu.: 5.000   1st Qu.:0.0000    F:374   I:390   I   :  0  
    ##  Median :0.1000   Median : 5.000   Median :0.0000                    F   :  0  
    ##  Mean   :0.1888   Mean   : 9.199   Mean   :0.4081                    FI  :  0  
    ##  3rd Qu.:0.3000   3rd Qu.:10.000   3rd Qu.:1.0000                    NA's:573  
    ##  Max.   :0.9000   Max.   :40.000   Max.   :1.0000                              
    ##  NA's   :3                                                                     
    ##     tot_dens      sp_mean_condens        Pr_s           con_dens      
    ##  Min.   :  1.00   Min.   :0.01261   Min.   :0.0000   Min.   :  1.000  
    ##  1st Qu.:  3.00   1st Qu.:0.12793   1st Qu.:0.0000   1st Qu.:  1.000  
    ##  Median :  6.00   Median :0.24685   Median :0.8333   Median :  1.000  
    ##  Mean   : 15.79   Mean   :1.49738   Mean   :0.6190   Mean   :  6.413  
    ##  3rd Qu.: 14.00   3rd Qu.:1.14595   3rd Qu.:1.0000   3rd Qu.:  3.000  
    ##  Max.   :301.00   Max.   :6.48288   Max.   :1.0000   Max.   :300.000  
    ##                                                                       
    ##    con_dens_s       slope.degrees_s      tot_dens_s      
    ##  Min.   :  0.1543   Min.   :-1.08744   Min.   :  0.1543  
    ##  1st Qu.:  1.7732   1st Qu.:-0.49640   1st Qu.:  9.8715  
    ##  Median :  6.8841   Median :-0.49640   Median : 23.5926  
    ##  Mean   : 14.9070   Mean   : 0.00000   Mean   : 39.0650  
    ##  3rd Qu.: 13.2143   3rd Qu.: 0.09464   3rd Qu.: 47.7227  
    ##  Max.   :355.2000   Max.   : 3.64087   Max.   :357.8179  
    ## 

``` r
dim(sdls) ## 968 species x plot combinations
```

    ## [1] 968  27

# Summary information for the paper

``` r
sdls |> summarise(n_sdls = sum(census.start), 
                  n_survs = sum(census.final),
                  n_species = n_distinct(species)) |> knitr::kable()
```

| n_sdls | n_survs | n_species |
|-------:|--------:|----------:|
|   6208 |    2680 |        26 |

``` r
sdls <- droplevels(sdls)
table(sdls$species)[order(table(sdls$species))]
```

    ## 
    ##     Here       SN     Humb   Crypto     Holi      Mac   Rattan    Toona 
    ##        5        5        6        7        7        8        8        8 
    ##       SP  Beilsch   Acacia Mallotus      Mel    Trich     Poly    UID20 
    ##        9       12       15       16       19       21       26       31 
    ##     Olea  Climber  Psyflav     Dimo   Litsea      FLC Climber1    Cinam 
    ##       36       45       46       56       61       64       77       86 
    ##     Symp       SR 
    ##      122      172

# Models

## Scaled conspecific density models

The most abundant species initially will often have lower survival
(fecundity/ survival trade-off). This could generate what looks like a
density-dependent relationship when looking across species, even without
a relationship within species (i.e., Simpson’s paradox).

``` r
ggplot(sp_mean_dens, aes(x = log(sp_mean_condens), y = sp_mean_surv)) +
  geom_label(aes(label= species)) + geom_smooth(method = "lm") + theme_tufte()
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](Kadumane_seedling_CDD_files/figure-gfm/abundance_survival-1.png)<!-- -->

``` r
## not the clearest pattern, but worth accounting for.
(ggplot(sdls, aes(x = scale((con_dens_s)), colour = species)) + geom_density())/
(ggplot(sdls, aes(x = scale((con_dens)), colour = species)) + 
   geom_density()) + plot_layout(guides = "collect")
```

![](Kadumane_seedling_CDD_files/figure-gfm/abundance_survival-2.png)<!-- -->

To account for this, we can scale conspecific density by dividing by
mean density and refitting the models.

``` r
## Random intercept model
m_cdd_s_ri <- glmmTMB(Pr_s ~ slope.degrees_s + 
                        trt_I:trt_F +
                         (scale(tot_dens_s) + scale(con_dens_s)) *
                         (trt_I + trt_F) *
                         scale(fragment.size)  +
                        (1|species) +
                        (1|site/loc/gr/plot), 
                      weights = census.start, data = sdls, 
                      family=binomial)

## Random intercept and slope model for species specific effects
m_cdd_s_ris <- glmmTMB(Pr_s ~ slope.degrees_s + 
                         trt_I:trt_F +
                         (scale(tot_dens_s) + scale(con_dens_s)) *
                         (trt_I + trt_F) *
                         scale(fragment.size)  +
                         (scale(con_dens_s) + scale(tot_dens_s)|species) +
                         (1|site/loc/gr/plot), 
                       weights = census.start, 
                       data = sdls, 
                       family=binomial)
## note the slightly odd ordering of terms doesn't change the model structure,
## but does change the default ordering of terms in outputs and plots to a
## more convenient one for describing in paper (first CDD, then biocide effects
## then fragmentation effects).

anova(m_cdd_s_ri, m_cdd_s_ris) ## random slope *much* better
```

    ## Data: sdls
    ## Models:
    ## m_cdd_s_ri: Pr_s ~ slope.degrees_s + trt_I:trt_F + (scale(tot_dens_s) + scale(con_dens_s)) * , zi=~0, disp=~1
    ## m_cdd_s_ri:     (trt_I + trt_F) * scale(fragment.size) + (1 | species) + , zi=~0, disp=~1
    ## m_cdd_s_ri:     (1 | site/loc/gr/plot), zi=~0, disp=~1
    ## m_cdd_s_ris: Pr_s ~ slope.degrees_s + trt_I:trt_F + (scale(tot_dens_s) + scale(con_dens_s)) * , zi=~0, disp=~1
    ## m_cdd_s_ris:     (trt_I + trt_F) * scale(fragment.size) + (scale(con_dens_s) + , zi=~0, disp=~1
    ## m_cdd_s_ris:     scale(tot_dens_s) | species) + (1 | site/loc/gr/plot), zi=~0, disp=~1
    ##             Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)    
    ## m_cdd_s_ri  25 2232.1 2353.9 -1091.0   2182.1                            
    ## m_cdd_s_ris 30 2215.1 2361.3 -1077.5   2155.1 27.02      5  5.654e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Models fit without issues, even with the correlation among density
effects.

The improvement with the random slopes model suggests we need to look at
individual species. Proceeding with the random intercept and slope model
from here on.

### Diagnostics

``` r
res_s <- simulateResiduals(m_cdd_s_ris)
plot(res_s) ## ok - some deviation from ideal residual distribution, but 
```

![](Kadumane_seedling_CDD_files/figure-gfm/diag_scaled-1.png)<!-- -->

``` r
## acceptable.

## look at relationship with covariates
diag_dat <- data.frame(m_cdd_s_ris$frame, res = res_s$fittedResiduals)

diag_dat <-   rename_with(diag_dat, ~ str_replace(.x, "scale\\.", "")) |> 
  rename_with(~str_replace(.x, "\\.$", "")) 

ggplot(diag_dat, aes(x = con_dens_s, y = res)) +  
  geom_point(aes(colour = species)) + 
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_smooth(method="gam") ## no trend. 
```

    ## `geom_smooth()` using formula = 'y ~ s(x, bs = "cs")'

![](Kadumane_seedling_CDD_files/figure-gfm/diag_scaled-2.png)<!-- -->

``` r
ggplot(diag_dat, aes(x = fragment.size, y = res)) +  
  geom_point(aes(colour = species), position = "jitter") + 
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_smooth(method="gam") ## no trend. 
```

    ## `geom_smooth()` using formula = 'y ~ s(x, bs = "cs")'

![](Kadumane_seedling_CDD_files/figure-gfm/diag_scaled-3.png)<!-- -->

``` r
ggplot(diag_dat, aes(x = con_dens_s, y = res)) +  
  facet_wrap(~trt_F + trt_I ) +
  geom_point(aes(colour = species)) +
  geom_smooth(method="gam") ## no trend with treatment
```

    ## `geom_smooth()` using formula = 'y ~ s(x, bs = "cs")'

![](Kadumane_seedling_CDD_files/figure-gfm/diag_scaled-4.png)<!-- -->

The diagnostics aren’t perfect, but not particularly unusual for a
binomial model. None of the big problems (e.g., overdispersion, trends
with covariates) seem to apply here.

### Model inference

``` r
summary(m_cdd_s_ris)
```

    ##  Family: binomial  ( logit )
    ## Formula:          
    ## Pr_s ~ slope.degrees_s + trt_I:trt_F + (scale(tot_dens_s) + scale(con_dens_s)) *  
    ##     (trt_I + trt_F) * scale(fragment.size) + (scale(con_dens_s) +  
    ##     scale(tot_dens_s) | species) + (1 | site/loc/gr/plot)
    ## Data: sdls
    ## Weights: census.start
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   2215.0   2361.3  -1077.5   2155.0      938 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups           Name              Variance  Std.Dev.  Corr        
    ##  species          (Intercept)       1.027e+00 1.0131846             
    ##                   scale(con_dens_s) 3.887e-01 0.6234621  0.56       
    ##                   scale(tot_dens_s) 3.738e-02 0.1933384 -0.47 -0.71 
    ##  plot:gr:loc:site (Intercept)       2.788e-01 0.5280003             
    ##  gr:loc:site      (Intercept)       2.057e-01 0.4535170             
    ##  loc:site         (Intercept)       4.782e-01 0.6915250             
    ##  site             (Intercept)       2.188e-08 0.0001479             
    ## Number of obs: 968, groups:  
    ## species, 26; plot:gr:loc:site, 474; gr:loc:site, 110; loc:site, 37; site, 21
    ## 
    ## Conditional model:
    ##                                               Estimate Std. Error z value
    ## (Intercept)                                    0.66723    0.27771   2.403
    ## slope.degrees_s                               -0.18268    0.07133  -2.561
    ## scale(tot_dens_s)                             -0.10921    0.15595  -0.700
    ## scale(con_dens_s)                             -0.44496    0.18820  -2.364
    ## trt_II                                         0.55988    0.15137   3.699
    ## trt_FF                                         0.27661    0.13886   1.992
    ## scale(fragment.size)                          -0.25279    0.15514  -1.629
    ## trt_II:trt_FF                                 -0.20355    0.21738  -0.936
    ## scale(tot_dens_s):trt_II                       0.04971    0.19245   0.258
    ## scale(tot_dens_s):trt_FF                      -0.01093    0.14816  -0.074
    ## scale(con_dens_s):trt_II                      -0.25189    0.16952  -1.486
    ## scale(con_dens_s):trt_FF                       0.38152    0.15970   2.389
    ## scale(tot_dens_s):scale(fragment.size)         0.09861    0.11221   0.879
    ## scale(con_dens_s):scale(fragment.size)        -0.24883    0.16012  -1.554
    ## trt_II:scale(fragment.size)                    0.08701    0.12608   0.690
    ## trt_FF:scale(fragment.size)                    0.11108    0.11225   0.990
    ## scale(tot_dens_s):trt_II:scale(fragment.size)  0.07160    0.22351   0.320
    ## scale(tot_dens_s):trt_FF:scale(fragment.size) -0.17654    0.20314  -0.869
    ## scale(con_dens_s):trt_II:scale(fragment.size) -0.14412    0.19987  -0.721
    ## scale(con_dens_s):trt_FF:scale(fragment.size)  0.55367    0.19600   2.825
    ##                                               Pr(>|z|)    
    ## (Intercept)                                   0.016277 *  
    ## slope.degrees_s                               0.010437 *  
    ## scale(tot_dens_s)                             0.483771    
    ## scale(con_dens_s)                             0.018061 *  
    ## trt_II                                        0.000217 ***
    ## trt_FF                                        0.046363 *  
    ## scale(fragment.size)                          0.103231    
    ## trt_II:trt_FF                                 0.349094    
    ## scale(tot_dens_s):trt_II                      0.796159    
    ## scale(tot_dens_s):trt_FF                      0.941199    
    ## scale(con_dens_s):trt_II                      0.137312    
    ## scale(con_dens_s):trt_FF                      0.016894 *  
    ## scale(tot_dens_s):scale(fragment.size)        0.379537    
    ## scale(con_dens_s):scale(fragment.size)        0.120172    
    ## trt_II:scale(fragment.size)                   0.490132    
    ## trt_FF:scale(fragment.size)                   0.322352    
    ## scale(tot_dens_s):trt_II:scale(fragment.size) 0.748722    
    ## scale(tot_dens_s):trt_FF:scale(fragment.size) 0.384811    
    ## scale(con_dens_s):trt_II:scale(fragment.size) 0.470860    
    ## scale(con_dens_s):trt_FF:scale(fragment.size) 0.004731 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot_model(m_cdd_s_ris, show.values=TRUE) + ylim(c(0.2, 3)) 
```

    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

![](Kadumane_seedling_CDD_files/figure-gfm/inference_scaled-1.png)<!-- -->

``` r
plot_model(m_cdd_s_ris, show.values=TRUE) + ylim(c(0.2, 3)) 
```

    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

![](Kadumane_seedling_CDD_files/figure-gfm/inference_scaled-2.png)<!-- -->

``` r
## total density is never very important, so separating out it's effects
rmvars <- names(fixef(m_cdd_s_ris)$cond)
rmvars <- c(rmvars[grep("tot_dens_s", rmvars)], "slope.degrees_s")

## First confirm it isn't influential
plot_model(m_cdd_s_ris, show.values=TRUE, terms=rmvars) + ylim(c(0.2, 3)) 
```

    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

![](Kadumane_seedling_CDD_files/figure-gfm/inference_scaled-3.png)<!-- -->

``` r
# basically what we thought

## now decluttered version
plot_model(m_cdd_s_ris, show.values=TRUE, rm.terms=rmvars) + ylim(c(0.2, 3)) 
```

    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

![](Kadumane_seedling_CDD_files/figure-gfm/inference_scaled-4.png)<!-- -->

``` r
## compare to random intercept model

plot_models(m_cdd_s_ri, m_cdd_s_ris, rm.terms=rmvars, m.labels=c("ri", "ris"),
            p.shape = TRUE)
```

![](Kadumane_seedling_CDD_files/figure-gfm/inference_scaled-5.png)<!-- -->

``` r
summary(m_cdd_s_ri)
```

    ##  Family: binomial  ( logit )
    ## Formula:          
    ## Pr_s ~ slope.degrees_s + trt_I:trt_F + (scale(tot_dens_s) + scale(con_dens_s)) *  
    ##     (trt_I + trt_F) * scale(fragment.size) + (1 | species) +  
    ##     (1 | site/loc/gr/plot)
    ## Data: sdls
    ## Weights: census.start
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   2232.1   2353.9  -1091.0   2182.1      943 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups           Name        Variance  Std.Dev. 
    ##  species          (Intercept) 1.496e+00 1.2230929
    ##  plot:gr:loc:site (Intercept) 3.554e-01 0.5961589
    ##  gr:loc:site      (Intercept) 2.921e-01 0.5404754
    ##  loc:site         (Intercept) 4.210e-01 0.6488833
    ##  site             (Intercept) 1.522e-08 0.0001234
    ## Number of obs: 968, groups:  
    ## species, 26; plot:gr:loc:site, 474; gr:loc:site, 110; loc:site, 37; site, 21
    ## 
    ## Conditional model:
    ##                                               Estimate Std. Error z value
    ## (Intercept)                                    0.81111    0.30834   2.631
    ## slope.degrees_s                               -0.18928    0.07399  -2.558
    ## scale(tot_dens_s)                              0.06247    0.11818   0.529
    ## scale(con_dens_s)                             -0.57768    0.10619  -5.440
    ## trt_II                                         0.59141    0.15620   3.786
    ## trt_FF                                         0.25563    0.14449   1.769
    ## scale(fragment.size)                          -0.32059    0.15085  -2.125
    ## trt_II:trt_FF                                 -0.20244    0.22300  -0.908
    ## scale(tot_dens_s):trt_II                       0.04069    0.19865   0.205
    ## scale(tot_dens_s):trt_FF                      -0.07806    0.15129  -0.516
    ## scale(con_dens_s):trt_II                      -0.09714    0.14638  -0.664
    ## scale(con_dens_s):trt_FF                       0.30573    0.13274   2.303
    ## scale(tot_dens_s):scale(fragment.size)         0.04495    0.11292   0.398
    ## scale(con_dens_s):scale(fragment.size)        -0.43936    0.11188  -3.927
    ## trt_II:scale(fragment.size)                    0.12382    0.12955   0.956
    ## trt_FF:scale(fragment.size)                    0.08116    0.11529   0.704
    ## scale(tot_dens_s):trt_II:scale(fragment.size)  0.01694    0.22765   0.074
    ## scale(tot_dens_s):trt_FF:scale(fragment.size) -0.07885    0.19887  -0.396
    ## scale(con_dens_s):trt_II:scale(fragment.size)  0.12721    0.17947   0.709
    ## scale(con_dens_s):trt_FF:scale(fragment.size)  0.36934    0.17565   2.103
    ##                                               Pr(>|z|)    
    ## (Intercept)                                   0.008524 ** 
    ## slope.degrees_s                               0.010518 *  
    ## scale(tot_dens_s)                             0.597063    
    ## scale(con_dens_s)                             5.33e-08 ***
    ## trt_II                                        0.000153 ***
    ## trt_FF                                        0.076857 .  
    ## scale(fragment.size)                          0.033566 *  
    ## trt_II:trt_FF                                 0.363984    
    ## scale(tot_dens_s):trt_II                      0.837712    
    ## scale(tot_dens_s):trt_FF                      0.605885    
    ## scale(con_dens_s):trt_II                      0.506951    
    ## scale(con_dens_s):trt_FF                      0.021267 *  
    ## scale(tot_dens_s):scale(fragment.size)        0.690614    
    ## scale(con_dens_s):scale(fragment.size)        8.60e-05 ***
    ## trt_II:scale(fragment.size)                   0.339178    
    ## trt_FF:scale(fragment.size)                   0.481443    
    ## scale(tot_dens_s):trt_II:scale(fragment.size) 0.940699    
    ## scale(tot_dens_s):trt_FF:scale(fragment.size) 0.691760    
    ## scale(con_dens_s):trt_II:scale(fragment.size) 0.478457    
    ## scale(con_dens_s):trt_FF:scale(fragment.size) 0.035493 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Take-homes:

- Survival is negatively conspecific density dependent. This result is
  consistent across all permutations of the model tried.

- Insecticide and fungicide (to a lesser extent) increase survival
  significantly.

- Fungicide removes the CDD, and that effect strengthens with fragment
  area. However, that effect needs more exploration (see below).

- The interaction between CDD and fragment size is complex, but perhaps
  explained by species turnover. The CDD effect strengthens with
  fragment area, but not significantly (or even close) in the (better)
  random slope model. In the random intercept model, the relationship is
  significant. This is probably because of species turnover - perhaps
  species that prefer larger fragments happen to be ones that show
  stronger density dependence?

- Effects of fragment area and any interactions with density need more
  probing as the effects are sensitive to model structure.

### Models split by categorical fragment size

The effect of fragment

``` r
# Using 3 categories
sdls <- mutate(sdls, 
               frag_sizeclass3 = cut(fragment.size, 
                                     round(quantile(site_dat$fragment.size, 
                                              c(0, 0.33, 0.66, 1))), 
                                        include.lowest=TRUE))

## quick look at how many fragments and seedlings within each 
## fragment size category

sdls |> group_by(frag_sizeclass3) |> summarise(n_frag = n_distinct(site),
                                               n_sdl = n())
```

    ## # A tibble: 3 x 3
    ##   frag_sizeclass3 n_frag n_sdl
    ##   <fct>            <int> <int>
    ## 1 [1,5]                6   124
    ## 2 (5,31]               9   342
    ## 3 (31,149]             6   502

``` r
m_cdd_s_ris_frag3 <- lapply(split(sdls, f= sdls$frag_sizeclass3), function(d){
  var_d <- d |> group_by(species) |> 
    summarise(var_dens = var(con_dens_s)) |> 
    filter(!is.na(var_dens) & var_dens > 0) 
  
  glmmTMB(Pr_s ~ slope.degrees_s + 
            trt_I:trt_F +
            (scale(tot_dens_s) + scale(con_dens_s)) *
            (trt_I + trt_F) +
            (scale(con_dens_s) + scale(tot_dens_s)||species) +
            (1|site/loc/gr/plot), 
          weights = census.start, 
          data = filter(d, species %in% var_d$species), 
          family=binomial)})
```

    ## Warning in finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old): Model convergence
    ## problem; singular convergence (7). See vignette('troubleshooting'),
    ## help('diagnose')

``` r
## can probably ignore the warning, which comes from the random 
## effect structure and tiny variances. We don't want to remove it because
## that would hamper comparison.

#diagnose(m_cdd_s_ris_frag3[[2]])

term_nms <- names(fixef(m_cdd_s_ris_frag3[[1]])$cond) 
names(m_cdd_s_ris_frag3) <- c("1 - 4", "5 - 32", "33 - 150") ## renaming to

## be clearer and correct (from looking at observed boundaries)
plot_discrete <- plot_models(m_cdd_s_ris_frag3, 
            rm.terms=c("slope.degrees_s", 
            term_nms[str_detect(term_nms, "tot")]), ## declutter
            m.labels=names(m_cdd_s_ris_frag3), p.shape=TRUE, 
            show.values = TRUE) +
  labs(colour = "Fragment Area")

pl_discrete <- plot_discrete$data |> 
  mutate(term = str_replace_all(term, fixed("scale(con_dens_s)"),"ConDens"),
         term = str_replace_all(term, "trt_[FI]", ""),
         term = factor(term, levels = rev(unique(term))), 
         FragArea = factor(group, levels = c("1 - 4", "5 - 32", "33 - 150")))

pl_discrete <- ggplot(pl_discrete, aes(x = estimate, xmin = conf.low, xmax = conf.high, 
                         y = term, colour = FragArea)) + 
  geom_pointrange(position = position_dodge2(width = 0.3)) +
  geom_vline(xintercept=1, linetype = "dotted") +
  scale_x_continuous(trans = "log10", n.breaks = 7) +
  scale_colour_viridis_d(option="D", end=0.8) +
  labs(x = "log(odds ratio)", y = NULL, colour = "Fragment\nArea (ha)")

pl_discrete
```

![](Kadumane_seedling_CDD_files/figure-gfm/mod_by_fragsize-1.png)<!-- -->

``` r
ggsave(pl_discrete, file = "figures/CatFragPlot.png", height = 5, width = 5)

map(m_cdd_s_ris_frag3, summary)
```

    ## $`1 - 4`
    ##  Family: binomial  ( logit )
    ## Formula:          
    ## Pr_s ~ slope.degrees_s + trt_I:trt_F + (scale(tot_dens_s) + scale(con_dens_s)) *  
    ##     (trt_I + trt_F) + (scale(con_dens_s) + scale(tot_dens_s) ||  
    ##     species) + (1 | site/loc/gr/plot)
    ## Data: filter(d, species %in% var_d$species)
    ## Weights: census.start
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    319.9    368.2   -142.0    283.9       90 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups           Name              Variance Std.Dev. Corr      
    ##  species          (Intercept)       0.20970  0.4579             
    ##                   scale(con_dens_s) 0.21980  0.4688   0.00      
    ##                   scale(tot_dens_s) 0.09292  0.3048   0.00 0.00 
    ##  plot:gr:loc:site (Intercept)       0.26650  0.5162             
    ##  gr:loc:site      (Intercept)       0.32066  0.5663             
    ##  loc:site         (Intercept)       0.22937  0.4789             
    ##  site             (Intercept)       0.22937  0.4789             
    ## Number of obs: 108, groups:  
    ## species, 8; plot:gr:loc:site, 68; gr:loc:site, 17; loc:site, 6; site, 6
    ## 
    ## Conditional model:
    ##                          Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)               0.96258    0.48338   1.991   0.0464 *
    ## slope.degrees_s          -0.09686    0.23156  -0.418   0.6757  
    ## scale(tot_dens_s)        -0.09204    0.34452  -0.267   0.7893  
    ## scale(con_dens_s)        -0.75158    0.53049  -1.417   0.1566  
    ## trt_II                   -0.11596    0.39854  -0.291   0.7711  
    ## trt_FF                   -0.38553    0.40475  -0.953   0.3408  
    ## trt_II:trt_FF             0.23016    0.62914   0.366   0.7145  
    ## scale(tot_dens_s):trt_II  0.44546    0.40368   1.103   0.2698  
    ## scale(tot_dens_s):trt_FF -0.23049    0.71206  -0.324   0.7462  
    ## scale(con_dens_s):trt_II -0.44574    0.65129  -0.684   0.4937  
    ## scale(con_dens_s):trt_FF  0.26750    0.60297   0.444   0.6573  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $`5 - 32`
    ##  Family: binomial  ( logit )
    ## Formula:          
    ## Pr_s ~ slope.degrees_s + trt_I:trt_F + (scale(tot_dens_s) + scale(con_dens_s)) *  
    ##     (trt_I + trt_F) + (scale(con_dens_s) + scale(tot_dens_s) ||  
    ##     species) + (1 | site/loc/gr/plot)
    ## Data: filter(d, species %in% var_d$species)
    ## Weights: census.start
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    717.0    784.3   -340.5    681.0      293 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups           Name              Variance  Std.Dev.  Corr      
    ##  species          (Intercept)       1.344e+00 1.159e+00           
    ##                   scale(con_dens_s) 1.527e-01 3.907e-01 0.00      
    ##                   scale(tot_dens_s) 8.330e-71 9.127e-36 0.00 0.00 
    ##  plot:gr:loc:site (Intercept)       3.356e-01 5.793e-01           
    ##  gr:loc:site      (Intercept)       3.435e-01 5.861e-01           
    ##  loc:site         (Intercept)       5.088e-01 7.133e-01           
    ##  site             (Intercept)       1.248e-06 1.117e-03           
    ## Number of obs: 311, groups:  
    ## species, 19; plot:gr:loc:site, 153; gr:loc:site, 36; loc:site, 12; site, 9
    ## 
    ## Conditional model:
    ##                          Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)               0.39506    0.42938   0.920  0.35754   
    ## slope.degrees_s          -0.35259    0.12759  -2.764  0.00572 **
    ## scale(tot_dens_s)        -0.04612    0.25566  -0.180  0.85684   
    ## scale(con_dens_s)        -0.25801    0.27010  -0.955  0.33947   
    ## trt_II                    0.85949    0.28140   3.054  0.00226 **
    ## trt_FF                    0.37630    0.28207   1.334  0.18219   
    ## trt_II:trt_FF            -0.62306    0.40586  -1.535  0.12474   
    ## scale(tot_dens_s):trt_II -0.25100    0.24672  -1.017  0.30899   
    ## scale(tot_dens_s):trt_FF  0.17258    0.26432   0.653  0.51379   
    ## scale(con_dens_s):trt_II -0.09498    0.20487  -0.464  0.64292   
    ## scale(con_dens_s):trt_FF -0.12890    0.25318  -0.509  0.61066   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $`33 - 150`
    ##  Family: binomial  ( logit )
    ## Formula:          
    ## Pr_s ~ slope.degrees_s + trt_I:trt_F + (scale(tot_dens_s) + scale(con_dens_s)) *  
    ##     (trt_I + trt_F) + (scale(con_dens_s) + scale(tot_dens_s) ||  
    ##     species) + (1 | site/loc/gr/plot)
    ## Data: filter(d, species %in% var_d$species)
    ## Weights: census.start
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   1123.1   1197.8   -543.6   1087.1      449 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups           Name              Variance  Std.Dev.  Corr      
    ##  species          (Intercept)       8.871e-01 9.418e-01           
    ##                   scale(con_dens_s) 2.247e-01 4.741e-01 0.00      
    ##                   scale(tot_dens_s) 7.522e-11 8.673e-06 0.00 0.00 
    ##  plot:gr:loc:site (Intercept)       2.416e-01 4.915e-01           
    ##  gr:loc:site      (Intercept)       2.724e-01 5.219e-01           
    ##  loc:site         (Intercept)       3.565e-01 5.970e-01           
    ##  site             (Intercept)       1.271e-08 1.127e-04           
    ## Number of obs: 467, groups:  
    ## species, 17; plot:gr:loc:site, 236; gr:loc:site, 56; loc:site, 19; site, 6
    ## 
    ## Conditional model:
    ##                          Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)               0.65966    0.33804   1.951 0.051006 .  
    ## slope.degrees_s          -0.09605    0.10450  -0.919 0.358013    
    ## scale(tot_dens_s)        -0.03133    0.15311  -0.205 0.837843    
    ## scale(con_dens_s)        -0.53039    0.22085  -2.402 0.016325 *  
    ## trt_II                    0.66851    0.19030   3.513 0.000443 ***
    ## trt_FF                    0.37203    0.18583   2.002 0.045281 *  
    ## trt_II:trt_FF            -0.21407    0.28842  -0.742 0.457956    
    ## scale(tot_dens_s):trt_II  0.04324    0.26389   0.164 0.869845    
    ## scale(tot_dens_s):trt_FF  0.01724    0.17809   0.097 0.922882    
    ## scale(con_dens_s):trt_II -0.14922    0.22863  -0.653 0.513979    
    ## scale(con_dens_s):trt_FF  0.50369    0.20459   2.462 0.013819 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Take-homes

The model gets a bit complex when split like this, but it serves as a
useful exploration of what’s going on.

- Conspecific density reduces survival most clearly in the larger
  fragments.

- The results of fungicide:cdens make a bit more sense when the data are
  broken up like this. Fungicide only interacts with cdens in the
  largest fragments, which makes sense, because these are the fragments
  where cdens has an effect.

- insects seem to only reduce plant survival in the medium and large
  fragments

### Graphics

``` r
## refitting model to control order of coefficients
m_cdd_s_ris <- glmmTMB(Pr_s ~ slope.degrees_s +
                         (scale(tot_dens_s) + scale(con_dens_s)) +
                         scale(fragment.size) +
                         trt_I*trt_F +
                         (scale(tot_dens_s) + scale(con_dens_s)) *
                         (trt_I + trt_F) *
                         scale(fragment.size)  +
                         (scale(con_dens_s) + scale(tot_dens_s)|species) +
                         ## setting cor to 0 to converge
                         (1|site/loc/gr/plot),
                       weights = census.start,
                       data = sdls,
                       family=binomial)
fixef(m_cdd_s_ris)
```

    ## 
    ## Conditional model:
    ##                                   (Intercept)  
    ##                                       0.66723  
    ##                               slope.degrees_s  
    ##                                      -0.18268  
    ##                             scale(tot_dens_s)  
    ##                                      -0.10921  
    ##                             scale(con_dens_s)  
    ##                                      -0.44496  
    ##                          scale(fragment.size)  
    ##                                      -0.25279  
    ##                                        trt_II  
    ##                                       0.55988  
    ##                                        trt_FF  
    ##                                       0.27661  
    ##                                 trt_II:trt_FF  
    ##                                      -0.20355  
    ##                      scale(tot_dens_s):trt_II  
    ##                                       0.04971  
    ##                      scale(tot_dens_s):trt_FF  
    ##                                      -0.01093  
    ##                      scale(con_dens_s):trt_II  
    ##                                      -0.25189  
    ##                      scale(con_dens_s):trt_FF  
    ##                                       0.38152  
    ##        scale(tot_dens_s):scale(fragment.size)  
    ##                                       0.09861  
    ##        scale(con_dens_s):scale(fragment.size)  
    ##                                      -0.24883  
    ##                   scale(fragment.size):trt_II  
    ##                                       0.08701  
    ##                   scale(fragment.size):trt_FF  
    ##                                       0.11108  
    ## scale(tot_dens_s):scale(fragment.size):trt_II  
    ##                                       0.07160  
    ## scale(tot_dens_s):scale(fragment.size):trt_FF  
    ##                                      -0.17654  
    ## scale(con_dens_s):scale(fragment.size):trt_II  
    ##                                      -0.14412  
    ## scale(con_dens_s):scale(fragment.size):trt_FF  
    ##                                       0.55367

``` r
labs <-  c("ConDens","FragArea", "I", "F",  "I:F", 
           "ConDens:I", "ConDens:F", "FragArea:\n ConDens", 
           "FragArea:I" , "FragArea:F", 
           "FragArea:\n (ConDens:I)", "FragArea:\n (ConDens:F)")

tp_cdd_s_ris <- 
  tidy(m_cdd_s_ris, conf.int = TRUE) |> 
  filter(effect == "fixed", 
         !str_detect(term, "tot_dens"), 
         !term %in% c("(Intercept)", "slope.degrees_s")) |> 
  mutate(labels = factor(labs, levels = rev(labs)), 
         Biocide = case_when(
           str_detect(term, "trt_II") ~  "Insecticide",
           str_detect(term, "trt_FF") ~ "Fungicide",
           .default =  "Control"),
         Biocide = factor(ifelse(str_detect(term, "trt_II:trt_FF"),
                                 "Both", Biocide), 
                          levels = c("Control", "Insecticide", 
                                     "Fungicide", "Both"))) |> 
  ggplot(aes(y = labels, x = estimate, xmin = conf.low, xmax = conf.high, 
             colour = Biocide)) + 
  geom_pointrange(aes(shape = p.value < 0.05 )) +
  geom_vline(xintercept=0, linetype = "dotted") +
  scale_colour_brewer(palette="Set2") +
  scale_shape_manual(values=c(1, 16), guide = "none" ) +
  labs(y = NULL, x = "log(odds ratio)")  + 
  guides(colour = guide_legend(title= "Biocide:", position = "top", 
                               direction = "horizontal", 
                               title.position = "left", title.hjust = 0.5)) +
  theme(
    legend.margin = margin(0, 0, 0, 0), 
    legend.justification.top = "left",
    legend.location = "plot",
    plot.title.position = "plot"
  )                                        
tp_cdd_s_ris
```

![](Kadumane_seedling_CDD_files/figure-gfm/plots_scaled-1.png)<!-- -->

``` r
p_cdd_s_ris <- predict_response(m_cdd_s_ris, 
                                terms=c("con_dens_s[1:80, by = 1]", 
                                        "trt_I", "trt_F",
                                        "fragment.size[5, 65, 125]"))

pr_cdd_s_ris <- residualize_over_grid(p_cdd_s_ris, m_cdd_s_ri) |> 
  bind_cols(species = sdls$species) |> 
  mutate(trt = factor(
    case_when(
      group == "0" & facet == "0" ~ "Control",
      group == "0" & facet == "F" ~ "Fungicide",
      group == "I" & facet == "0" ~ "Insecticide",
      group == "I" & facet == "F" ~ "Both"),
    levels = c("Control", "Insecticide", "Fungicide", "Both")),
    facet_lab = "Fragment Area (ha)") |> 
  rename("con_dens" = "x", "frag_area" = "panel") |> 
  group_by(con_dens, trt, frag_area, facet_lab, species) |> 
  summarise(predicted = mean(predicted), n = n())
```

    ## `summarise()` has grouped output by 'con_dens', 'trt', 'frag_area',
    ## 'facet_lab'. You can override using the `.groups` argument.

``` r
p_cdd_s_ris <- as.data.frame(p_cdd_s_ris) |> 
  mutate(trt = factor(
    case_when(
      group == "0" & facet == "0" ~ "Control",
      group == "0" & facet == "F" ~ "Fungicide",
      group == "I" & facet == "0" ~ "Insecticide",
      group == "I" & facet == "F" ~ "Both"),
    levels = c("Control", "Insecticide", "Fungicide", "Both")),
    facet_lab = "Fragment Area (ha)") |> 
  rename("con_dens" = "x", "frag_area" = "panel")

# summary(p_cdd_s_ris)
# summary(pr_cdd_s_ris)

pl_cdd_s <- ggplot(p_cdd_s_ris, 
       aes(x=con_dens, y = predicted, colour = trt)) +
  ggh4x::facet_nested(~facet_lab + frag_area) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = trt), 
              colour = NA, alpha = 0.3) + 
  geom_line(linewidth = 1.2) +
  geom_point(data = pr_cdd_s_ris) +
  scale_colour_brewer(palette="Set2") + 
  scale_fill_brewer(palette="Set2", guide = "none") +
  labs(x = expression(Conspecific~density~(m^-2)), y = "Pr(Survival)", 
       colour = "Treatment") + theme(legend.position="none")

mod_plot <- (pl_cdd_s  | tp_cdd_s_ris )  + 
  plot_layout(widths=c(0.6, 0.4))
mod_plot
```

![](Kadumane_seedling_CDD_files/figure-gfm/plots_scaled-2.png)<!-- -->

``` r
# the predictions with the "both" treatment get really busy. Given the 
# lack of an interaction, trying without both

pl_cdd_s_2 <- ggplot(filter(p_cdd_s_ris, trt != "Both"), 
       aes(x=con_dens, y = predicted, colour = trt)) +
  ggh4x::facet_nested(~facet_lab + frag_area) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = trt), 
              colour = NA, alpha = 0.3) + 
  geom_line(linewidth = 1.2) +
  geom_point(data = filter(pr_cdd_s_ris, trt != "Both")) +
  scale_colour_brewer(palette="Set2") + 
  scale_fill_brewer(palette="Set2", guide = "none") +
  labs(x = expression(Conspecific~density~(m^-2)), y = "Pr(Survival)", 
       colour = "Treatment") + theme(legend.position="none")

mod_plot2 <- (pl_cdd_s_2  | tp_cdd_s_ris )  + 
  plot_layout(widths=c(0.6, 0.4)) + plot_annotation(tag_level = "A") ## better
mod_plot2
```

![](Kadumane_seedling_CDD_files/figure-gfm/plots_scaled-3.png)<!-- -->

``` r
ggsave(mod_plot2, file = "figures/model_results.png", height = 5.5, width = 7.5)
```

### Plotting effects of fragment area

The effects are complex and hinge on a non-significant interaction
between fragment area and the fungicide effect. To interrogate them
more, here are some plots (basically Neyman-Johnson interaction plots)

``` r
## plot the effect of density over fragment size
## helps to scale and log data outside the model formula to set this up
int_data <- sdls |> 
  mutate(con_dens_ss = scale(con_dens_s),
         tot_dens_ss = scale(tot_dens_s),
         fragsize_s = scale(fragment.size))

int_mod <- glmmTMB(Pr_s ~ slope.degrees_s + 
                     trt_I:trt_F +
                     (tot_dens_ss + con_dens_ss) *  
                     (trt_I + trt_F) * fragsize_s +
                     (1 + tot_dens_ss + con_dens_ss | species) +  
                     (1 | site/loc/gr/plot), 
                   weights = census.start, family = binomial, 
                   data = int_data)

preddat <- expand.grid(
  trt_F = levels(int_data$trt_F),
  trt_I = levels(int_data$trt_I),
  con_dens_ss = 1,
  fragsize_s = seq(-1, 2.1, length = 20)) |> 
  filter(!(trt_F == "F" & trt_I == "I"))

xmat <- model.matrix(~ con_dens_ss + con_dens_ss:(trt_F + trt_I + fragsize_s) +
                       con_dens_ss:trt_F:fragsize_s + 
                       con_dens_ss:trt_I:fragsize_s,
                     preddat)[ , -1]

preddat$int_hat <- as.vector(xmat %*% fixef(int_mod)$cond[colnames(xmat)])
vmat <- (vcov(int_mod)$cond[colnames(xmat), colnames(xmat)])
preddat$int_se <-  sqrt(diag(xmat %*% vmat %*% t(xmat)))

preddat <- preddat |> mutate(.lower = int_hat - 1.96*int_se, 
                             .upper = int_hat + 1.96*int_se, 
                             frag_area = fragsize_s*
                               sd(sdls$fragment.size) +
                               mean(sdls$fragment.size),
                             trt = factor(case_when(
                               trt_F == "0" & trt_I == "0" ~ "Control",
                               trt_F == "F" & trt_I == "0" ~ "Fungicide",
                               trt_F == "0" & trt_I == "I" ~ "Insecticide")))

int_plot <- ggplot(preddat, aes(x = frag_area)) + facet_wrap( ~ trt) + 
  geom_ribbon(aes(y = int_hat, ymin = .lower, ymax = .upper),
              alpha = 0.2, colour = NA) +
  geom_line(aes(y = int_hat)) + geom_hline(yintercept=0, linetype = "dotted") +
  labs(x =  "Fragment area (ha)", 
       y = "Effect of conspecific density on survival")


sdls <- mutate(sdls, trt = case_when(
  trt_F == "F" ~ "F",
  trt_I == "I" ~ "I",
  .default = "C"), 
  trt = factor(trt, labels = c("Control", "Fungicide", "Insecticide")))

int_plot <- int_plot + 
  ggdist::geom_dots(data = sdls, aes(x = fragment.size), y =  -2.8,
                    smooth = ggdist::smooth_unbounded(), layout = "swarm", 
                    side = "top", binwidth = 1,# alpha = 0.7, 
                    overflow = "compress")
int_plot
```

![](Kadumane_seedling_CDD_files/figure-gfm/fragarea_plots-1.png)<!-- -->

``` r
ggsave(int_plot, file = "figures/interaction_plot.png", width=6.6, height = 4)
```

## Species specific inferences.

Initial analyses suggested that some species were heavily influencing
the results.

The random slope should partly account for that.

``` r
sjPlot::plot_model(m_cdd_s_ris, type = "re", terms = "species", ri.nr = 1)
```

![](Kadumane_seedling_CDD_files/figure-gfm/species_effects-1.png)<!-- -->

``` r
## observable variation among species in CDD slopes, but not in tot_dens_s
sp_eff <- as.data.frame(ranef(m_cdd_s_ris, condVar = T)$cond$species)

sp_eff <- bind_cols(sp_eff, 
                    t(apply(attr(sp_eff, "condVar"), 3, function(x)
                      sqrt(diag(x)))))
```

    ## New names:
    ## * `` -> `...4`
    ## * `` -> `...5`
    ## * `` -> `...6`

``` r
names(sp_eff) <- c("Intercept", "con_dens", "tot_dens", 
                   "Intercept_se", "con_dens_se", "tot_dens_se")

## extract total effects
sp_eff <- mutate(sp_eff, sp_names = row.names(sp_eff),
                 Intercept = Intercept + fixef(m_cdd_s_ris)$cond[1],
                 con_dens = con_dens +
                   fixef(m_cdd_s_ris)$cond["scale(con_dens_s)"],
                 tot_dens = tot_dens +
                   fixef(m_cdd_s_ris)$cond["scale(tot_dens)"])

sp_eff <- left_join(sp_eff, sp_codes, by = c("sp_names" = "code"))
## wrangle species names
sp_eff <- sp_eff |> mutate(spbin = case_when(
  genus %in% c("Unidentified", "Meliac")   ~ paste(genus, species, sep = "~"),
  genus != "Unidentified" & str_starts("sp", species) ~ 
                    paste(paste0("italic('", genus, "')"), species, sep = "~"),
  .default =    paste0("italic('", paste(genus, species, sep = " "), "')")))

blup_plot <- sp_eff |> arrange(con_dens) |> 
  ggplot(aes(y = reorder(spbin, con_dens), x = con_dens, 
             xmin = con_dens - 2*con_dens_se, 
             xmax = con_dens + 2*con_dens_se)) +
    geom_vline(xintercept=0, linetype = "dotted") +
  geom_pointrange() + 
  scale_y_discrete(labels = scales::label_parse()) +
  labs(x = "Conspecific density effect", y = "Species")

blup_plot <- blup_plot + theme(axis.text.y = element_text(size = 7))
blup_plot
```

![](Kadumane_seedling_CDD_files/figure-gfm/species_effects-2.png)<!-- -->

``` r
# ggsave(blup_plot, file = "figures/cdd_blups.png", width = 4, height =5)

anova(m_cdd_s_ris, m_cdd_s_ri)
```

    ## Data: sdls
    ## Models:
    ## m_cdd_s_ri: Pr_s ~ slope.degrees_s + trt_I:trt_F + (scale(tot_dens_s) + scale(con_dens_s)) * , zi=~0, disp=~1
    ## m_cdd_s_ri:     (trt_I + trt_F) * scale(fragment.size) + (1 | species) + , zi=~0, disp=~1
    ## m_cdd_s_ri:     (1 | site/loc/gr/plot), zi=~0, disp=~1
    ## m_cdd_s_ris: Pr_s ~ slope.degrees_s + (scale(tot_dens_s) + scale(con_dens_s)) + , zi=~0, disp=~1
    ## m_cdd_s_ris:     scale(fragment.size) + trt_I * trt_F + (scale(tot_dens_s) + , zi=~0, disp=~1
    ## m_cdd_s_ris:     scale(con_dens_s)) * (trt_I + trt_F) * scale(fragment.size) + , zi=~0, disp=~1
    ## m_cdd_s_ris:     (scale(con_dens_s) + scale(tot_dens_s) | species) + (1 | , zi=~0, disp=~1
    ## m_cdd_s_ris:     site/loc/gr/plot), zi=~0, disp=~1
    ##             Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)    
    ## m_cdd_s_ri  25 2232.1 2353.9 -1091.0   2182.1                            
    ## m_cdd_s_ris 30 2215.1 2361.3 -1077.5   2155.1 27.02      5  5.654e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Substantial variation in cdd among species -

Perhaps remove most abundant species to check if patterns are robust.

``` r
sdls |> group_by(species) |> summarise(n = sum(census.start)) |> 
   arrange(desc(n)) |> print(n = 26)
```

    ## # A tibble: 26 x 2
    ##    species      n
    ##    <fct>    <int>
    ##  1 SR        3598
    ##  2 Symp       636
    ##  3 Climber1   627
    ##  4 FLC        313
    ##  5 Cinam      137
    ##  6 Psyflav    116
    ##  7 Dimo        99
    ##  8 SP          85
    ##  9 Litsea      72
    ## 10 Olea        71
    ## 11 UID20       65
    ## 12 Climber     59
    ## 13 Poly        54
    ## 14 Mel         43
    ## 15 Beilsch     42
    ## 16 Holi        34
    ## 17 Toona       28
    ## 18 SN          25
    ## 19 Trich       22
    ## 20 Mallotus    19
    ## 21 Acacia      16
    ## 22 Mac         12
    ## 23 Rattan      11
    ## 24 Humb         9
    ## 25 Crypto       8
    ## 26 Here         7

``` r
sp_common <- sdls |> group_by(species) |> summarise(n = sum(census.start)) |> 
   arrange(desc(n)) |> pull(species)
sp_mods <- lapply(c("none", as.character(sp_common[1:5])), function(i) {
  update(m_cdd_s_ris, data = filter(sdls, !species == i))})
```

    ## Warning in (function (start, objective, gradient = NULL, hessian = NULL, :
    ## NA/NaN function evaluation
    ## Warning in (function (start, objective, gradient = NULL, hessian = NULL, :
    ## NA/NaN function evaluation
    ## Warning in (function (start, objective, gradient = NULL, hessian = NULL, :
    ## NA/NaN function evaluation

    ## Warning in finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old): Model convergence
    ## problem; non-positive-definite Hessian matrix. See vignette('troubleshooting')

    ## Warning in finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old): Model convergence
    ## problem; singular convergence (7). See vignette('troubleshooting'),
    ## help('diagnose')

``` r
sp_codes <- mutate(sp_codes, spbin = paste(genus, species))
term_nms <- names(fixef(sp_mods[[1]])$cond)
names(sp_mods) <- c("None", as.character(sp_common[1:5]))
plot_models(sp_mods,
            rm.terms = c("slope.degrees_s", 
                         term_nms[str_detect(term_nms, "tot")]),
            m.labels = c("None", 
                         sp_codes$spbin[match(sp_common[1:5], sp_codes$code)]),
            p.shape=TRUE) + labs(colour = "Species excluded") +
  scale_color_viridis_d(end=0.9) ## getting rid of the awful yellow.
```

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

![](Kadumane_seedling_CDD_files/figure-gfm/species_influence-1.png)<!-- -->

``` r
map(sp_mods, summary)
```

    ## $None
    ##  Family: binomial  ( logit )
    ## Formula:          
    ## Pr_s ~ slope.degrees_s + (scale(tot_dens_s) + scale(con_dens_s)) +  
    ##     scale(fragment.size) + trt_I * trt_F + (scale(tot_dens_s) +  
    ##     scale(con_dens_s)) * (trt_I + trt_F) * scale(fragment.size) +  
    ##     (scale(con_dens_s) + scale(tot_dens_s) | species) + (1 |  
    ##     site/loc/gr/plot)
    ## Data: filter(sdls, !species == i)
    ## Weights: census.start
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   2215.0   2361.3  -1077.5   2155.0      938 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups           Name              Variance  Std.Dev.  Corr        
    ##  species          (Intercept)       1.027e+00 1.0131846             
    ##                   scale(con_dens_s) 3.887e-01 0.6234621  0.56       
    ##                   scale(tot_dens_s) 3.738e-02 0.1933384 -0.47 -0.71 
    ##  plot:gr:loc:site (Intercept)       2.788e-01 0.5280003             
    ##  gr:loc:site      (Intercept)       2.057e-01 0.4535170             
    ##  loc:site         (Intercept)       4.782e-01 0.6915250             
    ##  site             (Intercept)       2.188e-08 0.0001479             
    ## Number of obs: 968, groups:  
    ## species, 26; plot:gr:loc:site, 474; gr:loc:site, 110; loc:site, 37; site, 21
    ## 
    ## Conditional model:
    ##                                               Estimate Std. Error z value
    ## (Intercept)                                    0.66723    0.27771   2.403
    ## slope.degrees_s                               -0.18268    0.07133  -2.561
    ## scale(tot_dens_s)                             -0.10921    0.15595  -0.700
    ## scale(con_dens_s)                             -0.44496    0.18820  -2.364
    ## scale(fragment.size)                          -0.25279    0.15514  -1.629
    ## trt_II                                         0.55988    0.15137   3.699
    ## trt_FF                                         0.27661    0.13886   1.992
    ## trt_II:trt_FF                                 -0.20355    0.21738  -0.936
    ## scale(tot_dens_s):trt_II                       0.04971    0.19245   0.258
    ## scale(tot_dens_s):trt_FF                      -0.01093    0.14816  -0.074
    ## scale(con_dens_s):trt_II                      -0.25189    0.16952  -1.486
    ## scale(con_dens_s):trt_FF                       0.38152    0.15970   2.389
    ## scale(tot_dens_s):scale(fragment.size)         0.09861    0.11221   0.879
    ## scale(con_dens_s):scale(fragment.size)        -0.24883    0.16012  -1.554
    ## scale(fragment.size):trt_II                    0.08701    0.12608   0.690
    ## scale(fragment.size):trt_FF                    0.11108    0.11225   0.990
    ## scale(tot_dens_s):scale(fragment.size):trt_II  0.07160    0.22351   0.320
    ## scale(tot_dens_s):scale(fragment.size):trt_FF -0.17654    0.20314  -0.869
    ## scale(con_dens_s):scale(fragment.size):trt_II -0.14412    0.19987  -0.721
    ## scale(con_dens_s):scale(fragment.size):trt_FF  0.55367    0.19600   2.825
    ##                                               Pr(>|z|)    
    ## (Intercept)                                   0.016277 *  
    ## slope.degrees_s                               0.010437 *  
    ## scale(tot_dens_s)                             0.483771    
    ## scale(con_dens_s)                             0.018061 *  
    ## scale(fragment.size)                          0.103231    
    ## trt_II                                        0.000217 ***
    ## trt_FF                                        0.046363 *  
    ## trt_II:trt_FF                                 0.349094    
    ## scale(tot_dens_s):trt_II                      0.796159    
    ## scale(tot_dens_s):trt_FF                      0.941199    
    ## scale(con_dens_s):trt_II                      0.137312    
    ## scale(con_dens_s):trt_FF                      0.016894 *  
    ## scale(tot_dens_s):scale(fragment.size)        0.379537    
    ## scale(con_dens_s):scale(fragment.size)        0.120172    
    ## scale(fragment.size):trt_II                   0.490132    
    ## scale(fragment.size):trt_FF                   0.322352    
    ## scale(tot_dens_s):scale(fragment.size):trt_II 0.748722    
    ## scale(tot_dens_s):scale(fragment.size):trt_FF 0.384811    
    ## scale(con_dens_s):scale(fragment.size):trt_II 0.470860    
    ## scale(con_dens_s):scale(fragment.size):trt_FF 0.004731 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $SR
    ##  Family: binomial  ( logit )
    ## Formula:          
    ## Pr_s ~ slope.degrees_s + (scale(tot_dens_s) + scale(con_dens_s)) +  
    ##     scale(fragment.size) + trt_I * trt_F + (scale(tot_dens_s) +  
    ##     scale(con_dens_s)) * (trt_I + trt_F) * scale(fragment.size) +  
    ##     (scale(con_dens_s) + scale(tot_dens_s) | species) + (1 |  
    ##     site/loc/gr/plot)
    ## Data: filter(sdls, !species == i)
    ## Weights: census.start
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   1538.8   1679.2   -739.4   1478.8      766 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups           Name              Variance  Std.Dev. Corr        
    ##  species          (Intercept)       9.948e-01 0.997383             
    ##                   scale(con_dens_s) 1.666e-01 0.408201  0.49       
    ##                   scale(tot_dens_s) 5.029e-02 0.224248 -0.68 -0.95 
    ##  plot:gr:loc:site (Intercept)       9.819e-02 0.313348             
    ##  gr:loc:site      (Intercept)       1.978e-01 0.444692             
    ##  loc:site         (Intercept)       4.018e-01 0.633879             
    ##  site             (Intercept)       1.277e-08 0.000113             
    ## Number of obs: 796, groups:  
    ## species, 25; plot:gr:loc:site, 436; gr:loc:site, 110; loc:site, 37; site, 21
    ## 
    ## Conditional model:
    ##                                               Estimate Std. Error z value
    ## (Intercept)                                    0.63850    0.27529   2.319
    ## slope.degrees_s                               -0.02056    0.08536  -0.241
    ## scale(tot_dens_s)                             -0.17944    0.16530  -1.085
    ## scale(con_dens_s)                             -0.30396    0.15773  -1.927
    ## scale(fragment.size)                          -0.18107    0.15650  -1.157
    ## trt_II                                         0.49663    0.16145   3.076
    ## trt_FF                                         0.41105    0.16944   2.426
    ## trt_II:trt_FF                                 -0.25253    0.24472  -1.032
    ## scale(tot_dens_s):trt_II                      -0.08794    0.22453  -0.392
    ## scale(tot_dens_s):trt_FF                      -0.03447    0.19921  -0.173
    ## scale(con_dens_s):trt_II                       0.02762    0.18700   0.148
    ## scale(con_dens_s):trt_FF                       0.14500    0.18050   0.803
    ## scale(tot_dens_s):scale(fragment.size)        -0.06969    0.14652  -0.476
    ## scale(con_dens_s):scale(fragment.size)        -0.14534    0.15883  -0.915
    ## scale(fragment.size):trt_II                    0.21404    0.14875   1.439
    ## scale(fragment.size):trt_FF                    0.18394    0.15205   1.210
    ## scale(tot_dens_s):scale(fragment.size):trt_II  0.11592    0.27765   0.418
    ## scale(tot_dens_s):scale(fragment.size):trt_FF -0.10503    0.25202  -0.417
    ## scale(con_dens_s):scale(fragment.size):trt_II -0.04202    0.23315  -0.180
    ## scale(con_dens_s):scale(fragment.size):trt_FF  0.08224    0.23157   0.355
    ##                                               Pr(>|z|)   
    ## (Intercept)                                     0.0204 * 
    ## slope.degrees_s                                 0.8096   
    ## scale(tot_dens_s)                               0.2777   
    ## scale(con_dens_s)                               0.0540 . 
    ## scale(fragment.size)                            0.2473   
    ## trt_II                                          0.0021 **
    ## trt_FF                                          0.0153 * 
    ## trt_II:trt_FF                                   0.3021   
    ## scale(tot_dens_s):trt_II                        0.6953   
    ## scale(tot_dens_s):trt_FF                        0.8626   
    ## scale(con_dens_s):trt_II                        0.8826   
    ## scale(con_dens_s):trt_FF                        0.4218   
    ## scale(tot_dens_s):scale(fragment.size)          0.6343   
    ## scale(con_dens_s):scale(fragment.size)          0.3601   
    ## scale(fragment.size):trt_II                     0.1502   
    ## scale(fragment.size):trt_FF                     0.2264   
    ## scale(tot_dens_s):scale(fragment.size):trt_II   0.6763   
    ## scale(tot_dens_s):scale(fragment.size):trt_FF   0.6769   
    ## scale(con_dens_s):scale(fragment.size):trt_II   0.8570   
    ## scale(con_dens_s):scale(fragment.size):trt_FF   0.7225   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $Symp
    ##  Family: binomial  ( logit )
    ## Formula:          
    ## Pr_s ~ slope.degrees_s + (scale(tot_dens_s) + scale(con_dens_s)) +  
    ##     scale(fragment.size) + trt_I * trt_F + (scale(tot_dens_s) +  
    ##     scale(con_dens_s)) * (trt_I + trt_F) * scale(fragment.size) +  
    ##     (scale(con_dens_s) + scale(tot_dens_s) | species) + (1 |  
    ##     site/loc/gr/plot)
    ## Data: filter(sdls, !species == i)
    ## Weights: census.start
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##       NA       NA       NA       NA      816 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups           Name              Variance   Std.Dev.   Corr        
    ##  species          (Intercept)        1.197e+00  1.094e+00             
    ##                   scale(con_dens_s)  5.227e-01  7.230e-01  0.47       
    ##                   scale(tot_dens_s)  1.641e-02  1.281e-01 -0.36 -0.99 
    ##  plot:gr:loc:site (Intercept)        2.593e-01  5.093e-01             
    ##  gr:loc:site      (Intercept)        2.719e-01  5.214e-01             
    ##  loc:site         (Intercept)        4.625e-01  6.801e-01             
    ##  site             (Intercept)       1.189e-295 3.449e-148             
    ## Number of obs: 846, groups:  
    ## species, 25; plot:gr:loc:site, 443; gr:loc:site, 110; loc:site, 37; site, 21
    ## 
    ## Conditional model:
    ##                                                Estimate Std. Error z value
    ## (Intercept)                                    0.623087   0.299655   2.079
    ## slope.degrees_s                               -0.238159   0.078978  -3.016
    ## scale(tot_dens_s)                             -0.172296   0.157359  -1.095
    ## scale(con_dens_s)                             -0.431181   0.217797  -1.980
    ## scale(fragment.size)                          -0.204384   0.161807  -1.263
    ## trt_II                                         0.642882   0.164703   3.903
    ## trt_FF                                         0.288345   0.145095   1.987
    ## trt_II:trt_FF                                 -0.156364   0.231996  -0.674
    ## scale(tot_dens_s):trt_II                      -0.045501   0.201656  -0.226
    ## scale(tot_dens_s):trt_FF                      -0.003621   0.155956  -0.023
    ## scale(con_dens_s):trt_II                      -0.254715   0.190675  -1.336
    ## scale(con_dens_s):trt_FF                       0.465605   0.178066   2.615
    ## scale(tot_dens_s):scale(fragment.size)         0.098487   0.117733   0.837
    ## scale(con_dens_s):scale(fragment.size)        -0.211136   0.172688  -1.223
    ## scale(fragment.size):trt_II                   -0.065994   0.135574  -0.487
    ## scale(fragment.size):trt_FF                    0.090493   0.119279   0.759
    ## scale(tot_dens_s):scale(fragment.size):trt_II  0.108954   0.228006   0.478
    ## scale(tot_dens_s):scale(fragment.size):trt_FF -0.251520   0.211014  -1.192
    ## scale(con_dens_s):scale(fragment.size):trt_II -0.171561   0.220288  -0.779
    ## scale(con_dens_s):scale(fragment.size):trt_FF  0.650795   0.214434   3.035
    ##                                               Pr(>|z|)    
    ## (Intercept)                                    0.03759 *  
    ## slope.degrees_s                                0.00257 ** 
    ## scale(tot_dens_s)                              0.27355    
    ## scale(con_dens_s)                              0.04773 *  
    ## scale(fragment.size)                           0.20654    
    ## trt_II                                        9.49e-05 ***
    ## trt_FF                                         0.04689 *  
    ## trt_II:trt_FF                                  0.50032    
    ## scale(tot_dens_s):trt_II                       0.82149    
    ## scale(tot_dens_s):trt_FF                       0.98148    
    ## scale(con_dens_s):trt_II                       0.18160    
    ## scale(con_dens_s):trt_FF                       0.00893 ** 
    ## scale(tot_dens_s):scale(fragment.size)         0.40286    
    ## scale(con_dens_s):scale(fragment.size)         0.22146    
    ## scale(fragment.size):trt_II                    0.62642    
    ## scale(fragment.size):trt_FF                    0.44805    
    ## scale(tot_dens_s):scale(fragment.size):trt_II  0.63275    
    ## scale(tot_dens_s):scale(fragment.size):trt_FF  0.23328    
    ## scale(con_dens_s):scale(fragment.size):trt_II  0.43610    
    ## scale(con_dens_s):scale(fragment.size):trt_FF  0.00241 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $Climber1
    ##  Family: binomial  ( logit )
    ## Formula:          
    ## Pr_s ~ slope.degrees_s + (scale(tot_dens_s) + scale(con_dens_s)) +  
    ##     scale(fragment.size) + trt_I * trt_F + (scale(tot_dens_s) +  
    ##     scale(con_dens_s)) * (trt_I + trt_F) * scale(fragment.size) +  
    ##     (scale(con_dens_s) + scale(tot_dens_s) | species) + (1 |  
    ##     site/loc/gr/plot)
    ## Data: filter(sdls, !species == i)
    ## Weights: census.start
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   2004.3   2148.1   -972.1   1944.3      861 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups           Name              Variance  Std.Dev.  Corr        
    ##  species          (Intercept)       1.110e+00 1.0536851             
    ##                   scale(con_dens_s) 4.279e-01 0.6541436  0.58       
    ##                   scale(tot_dens_s) 3.649e-02 0.1910162 -0.54 -0.62 
    ##  plot:gr:loc:site (Intercept)       3.125e-01 0.5589869             
    ##  gr:loc:site      (Intercept)       2.224e-01 0.4716092             
    ##  loc:site         (Intercept)       4.840e-01 0.6957357             
    ##  site             (Intercept)       1.601e-08 0.0001265             
    ## Number of obs: 891, groups:  
    ## species, 25; plot:gr:loc:site, 468; gr:loc:site, 110; loc:site, 37; site, 21
    ## 
    ## Conditional model:
    ##                                               Estimate Std. Error z value
    ## (Intercept)                                    0.72418    0.29263   2.475
    ## slope.degrees_s                               -0.18938    0.07488  -2.529
    ## scale(tot_dens_s)                             -0.05951    0.17064  -0.349
    ## scale(con_dens_s)                             -0.49321    0.21218  -2.325
    ## scale(fragment.size)                          -0.27528    0.16120  -1.708
    ## trt_II                                         0.52847    0.16187   3.265
    ## trt_FF                                         0.33985    0.14975   2.269
    ## trt_II:trt_FF                                 -0.23080    0.23385  -0.987
    ## scale(tot_dens_s):trt_II                       0.02794    0.19797   0.141
    ## scale(tot_dens_s):trt_FF                      -0.02253    0.15547  -0.145
    ## scale(con_dens_s):trt_II                      -0.24455    0.19155  -1.277
    ## scale(con_dens_s):trt_FF                       0.42565    0.18361   2.318
    ## scale(tot_dens_s):scale(fragment.size)         0.13497    0.11919   1.132
    ## scale(con_dens_s):scale(fragment.size)        -0.25574    0.17360  -1.473
    ## scale(fragment.size):trt_II                    0.05831    0.13213   0.441
    ## scale(fragment.size):trt_FF                    0.10946    0.11909   0.919
    ## scale(tot_dens_s):scale(fragment.size):trt_II  0.01717    0.23536   0.073
    ## scale(tot_dens_s):scale(fragment.size):trt_FF -0.16414    0.22013  -0.746
    ## scale(con_dens_s):scale(fragment.size):trt_II -0.10702    0.22543  -0.475
    ## scale(con_dens_s):scale(fragment.size):trt_FF  0.57845    0.22089   2.619
    ##                                               Pr(>|z|)   
    ## (Intercept)                                    0.01333 * 
    ## slope.degrees_s                                0.01144 * 
    ## scale(tot_dens_s)                              0.72728   
    ## scale(con_dens_s)                              0.02010 * 
    ## scale(fragment.size)                           0.08770 . 
    ## trt_II                                         0.00110 **
    ## trt_FF                                         0.02324 * 
    ## trt_II:trt_FF                                  0.32367   
    ## scale(tot_dens_s):trt_II                       0.88777   
    ## scale(tot_dens_s):trt_FF                       0.88479   
    ## scale(con_dens_s):trt_II                       0.20171   
    ## scale(con_dens_s):trt_FF                       0.02044 * 
    ## scale(tot_dens_s):scale(fragment.size)         0.25748   
    ## scale(con_dens_s):scale(fragment.size)         0.14071   
    ## scale(fragment.size):trt_II                    0.65901   
    ## scale(fragment.size):trt_FF                    0.35802   
    ## scale(tot_dens_s):scale(fragment.size):trt_II  0.94184   
    ## scale(tot_dens_s):scale(fragment.size):trt_FF  0.45590   
    ## scale(con_dens_s):scale(fragment.size):trt_II  0.63498   
    ## scale(con_dens_s):scale(fragment.size):trt_FF  0.00882 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $FLC
    ##  Family: binomial  ( logit )
    ## Formula:          
    ## Pr_s ~ slope.degrees_s + (scale(tot_dens_s) + scale(con_dens_s)) +  
    ##     scale(fragment.size) + trt_I * trt_F + (scale(tot_dens_s) +  
    ##     scale(con_dens_s)) * (trt_I + trt_F) * scale(fragment.size) +  
    ##     (scale(con_dens_s) + scale(tot_dens_s) | species) + (1 |  
    ##     site/loc/gr/plot)
    ## Data: filter(sdls, !species == i)
    ## Weights: census.start
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   2081.8   2226.0  -1010.9   2021.8      874 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups           Name              Variance  Std.Dev.  Corr        
    ##  species          (Intercept)       1.054e+00 1.0267284             
    ##                   scale(con_dens_s) 5.047e-01 0.7104215  0.55       
    ##                   scale(tot_dens_s) 8.428e-02 0.2903098 -0.44 -0.72 
    ##  plot:gr:loc:site (Intercept)       2.697e-01 0.5193591             
    ##  gr:loc:site      (Intercept)       1.919e-01 0.4380878             
    ##  loc:site         (Intercept)       4.780e-01 0.6913691             
    ##  site             (Intercept)       2.726e-08 0.0001651             
    ## Number of obs: 904, groups:  
    ## species, 25; plot:gr:loc:site, 466; gr:loc:site, 110; loc:site, 37; site, 21
    ## 
    ## Conditional model:
    ##                                               Estimate Std. Error z value
    ## (Intercept)                                    0.60130    0.28614   2.101
    ## slope.degrees_s                               -0.19046    0.07263  -2.622
    ## scale(tot_dens_s)                             -0.24177    0.17930  -1.348
    ## scale(con_dens_s)                             -0.35898    0.21142  -1.698
    ## scale(fragment.size)                          -0.26189    0.15705  -1.668
    ## trt_II                                         0.54049    0.15528   3.481
    ## trt_FF                                         0.27584    0.14241   1.937
    ## trt_II:trt_FF                                 -0.19031    0.22314  -0.853
    ## scale(tot_dens_s):trt_II                       0.10583    0.20414   0.518
    ## scale(tot_dens_s):trt_FF                       0.04725    0.15841   0.298
    ## scale(con_dens_s):trt_II                      -0.25219    0.18033  -1.398
    ## scale(con_dens_s):trt_FF                       0.35289    0.16866   2.092
    ## scale(tot_dens_s):scale(fragment.size)         0.17530    0.11231   1.561
    ## scale(con_dens_s):scale(fragment.size)        -0.28594    0.17139  -1.668
    ## scale(fragment.size):trt_II                    0.04834    0.12808   0.377
    ## scale(fragment.size):trt_FF                    0.16345    0.11585   1.411
    ## scale(tot_dens_s):scale(fragment.size):trt_II -0.05304    0.23021  -0.230
    ## scale(tot_dens_s):scale(fragment.size):trt_FF -0.13798    0.21737  -0.635
    ## scale(con_dens_s):scale(fragment.size):trt_II -0.07390    0.21099  -0.350
    ## scale(con_dens_s):scale(fragment.size):trt_FF  0.51571    0.20874   2.471
    ##                                               Pr(>|z|)    
    ## (Intercept)                                    0.03560 *  
    ## slope.degrees_s                                0.00873 ** 
    ## scale(tot_dens_s)                              0.17753    
    ## scale(con_dens_s)                              0.08951 .  
    ## scale(fragment.size)                           0.09540 .  
    ## trt_II                                         0.00050 ***
    ## trt_FF                                         0.05274 .  
    ## trt_II:trt_FF                                  0.39374    
    ## scale(tot_dens_s):trt_II                       0.60416    
    ## scale(tot_dens_s):trt_FF                       0.76549    
    ## scale(con_dens_s):trt_II                       0.16196    
    ## scale(con_dens_s):trt_FF                       0.03641 *  
    ## scale(tot_dens_s):scale(fragment.size)         0.11854    
    ## scale(con_dens_s):scale(fragment.size)         0.09525 .  
    ## scale(fragment.size):trt_II                    0.70588    
    ## scale(fragment.size):trt_FF                    0.15828    
    ## scale(tot_dens_s):scale(fragment.size):trt_II  0.81778    
    ## scale(tot_dens_s):scale(fragment.size):trt_FF  0.52559    
    ## scale(con_dens_s):scale(fragment.size):trt_II  0.72616    
    ## scale(con_dens_s):scale(fragment.size):trt_FF  0.01349 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $Cinam
    ##  Family: binomial  ( logit )
    ## Formula:          
    ## Pr_s ~ slope.degrees_s + (scale(tot_dens_s) + scale(con_dens_s)) +  
    ##     scale(fragment.size) + trt_I * trt_F + (scale(tot_dens_s) +  
    ##     scale(con_dens_s)) * (trt_I + trt_F) * scale(fragment.size) +  
    ##     (scale(con_dens_s) + scale(tot_dens_s) | species) + (1 |  
    ##     site/loc/gr/plot)
    ## Data: filter(sdls, !species == i)
    ## Weights: census.start
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   2059.4   2202.9   -999.7   1999.4      852 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups           Name              Variance  Std.Dev.  Corr        
    ##  species          (Intercept)       9.639e-01 0.9817926             
    ##                   scale(con_dens_s) 4.032e-01 0.6349760  0.45       
    ##                   scale(tot_dens_s) 4.948e-02 0.2224359 -0.38 -0.79 
    ##  plot:gr:loc:site (Intercept)       2.906e-01 0.5390979             
    ##  gr:loc:site      (Intercept)       2.685e-01 0.5181503             
    ##  loc:site         (Intercept)       4.422e-01 0.6649552             
    ##  site             (Intercept)       2.109e-08 0.0001452             
    ## Number of obs: 882, groups:  
    ## species, 25; plot:gr:loc:site, 461; gr:loc:site, 110; loc:site, 37; site, 21
    ## 
    ## Conditional model:
    ##                                               Estimate Std. Error z value
    ## (Intercept)                                    0.70050    0.27933   2.508
    ## slope.degrees_s                               -0.19170    0.07382  -2.597
    ## scale(tot_dens_s)                             -0.13076    0.16765  -0.780
    ## scale(con_dens_s)                             -0.48358    0.19780  -2.445
    ## scale(fragment.size)                          -0.27909    0.15808  -1.765
    ## trt_II                                         0.59369    0.15775   3.763
    ## trt_FF                                         0.31667    0.14362   2.205
    ## trt_II:trt_FF                                 -0.24942    0.22545  -1.106
    ## scale(tot_dens_s):trt_II                       0.05476    0.20594   0.266
    ## scale(tot_dens_s):trt_FF                      -0.03671    0.15834  -0.232
    ## scale(con_dens_s):trt_II                      -0.27294    0.18369  -1.486
    ## scale(con_dens_s):trt_FF                       0.47156    0.17190   2.743
    ## scale(tot_dens_s):scale(fragment.size)         0.11020    0.12150   0.907
    ## scale(con_dens_s):scale(fragment.size)        -0.35449    0.16559  -2.141
    ## scale(fragment.size):trt_II                    0.09773    0.13213   0.740
    ## scale(fragment.size):trt_FF                    0.10274    0.11743   0.875
    ## scale(tot_dens_s):scale(fragment.size):trt_II  0.09621    0.23656   0.407
    ## scale(tot_dens_s):scale(fragment.size):trt_FF -0.23337    0.21662  -1.077
    ## scale(con_dens_s):scale(fragment.size):trt_II -0.17272    0.21472  -0.804
    ## scale(con_dens_s):scale(fragment.size):trt_FF  0.64615    0.21256   3.040
    ##                                               Pr(>|z|)    
    ## (Intercept)                                   0.012149 *  
    ## slope.degrees_s                               0.009406 ** 
    ## scale(tot_dens_s)                             0.435430    
    ## scale(con_dens_s)                             0.014492 *  
    ## scale(fragment.size)                          0.077487 .  
    ## trt_II                                        0.000168 ***
    ## trt_FF                                        0.027460 *  
    ## trt_II:trt_FF                                 0.268580    
    ## scale(tot_dens_s):trt_II                      0.790319    
    ## scale(tot_dens_s):trt_FF                      0.816674    
    ## scale(con_dens_s):trt_II                      0.137320    
    ## scale(con_dens_s):trt_FF                      0.006083 ** 
    ## scale(tot_dens_s):scale(fragment.size)        0.364392    
    ## scale(con_dens_s):scale(fragment.size)        0.032297 *  
    ## scale(fragment.size):trt_II                   0.459497    
    ## scale(fragment.size):trt_FF                   0.381631    
    ## scale(tot_dens_s):scale(fragment.size):trt_II 0.684220    
    ## scale(tot_dens_s):scale(fragment.size):trt_FF 0.281345    
    ## scale(con_dens_s):scale(fragment.size):trt_II 0.421187    
    ## scale(con_dens_s):scale(fragment.size):trt_FF 0.002366 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Only species that makes a difference is S. rubicundum - removing it
dampens interaction between density and fungicide (it is still there,
but marginally non-significant)

what about single species models?

``` r
## base model.
## remove species random effect
names(sp_common) <- sp_common
single_sp_mods <- map(sp_common[1:5], function(i) {
  glmmTMB(Pr_s ~  slope.degrees_s +
            (scale(tot_dens_s) + scale(con_dens_s)) +
            scale(fragment.size) +
            trt_I*trt_F +
            (scale(tot_dens_s) + scale(con_dens_s)) *
            (trt_I + trt_F) *
            scale(fragment.size)  +
            ## setting cor to 0 to converge
            (1|site/loc/gr/plot), 
          weights = census.start, 
          data = filter(sdls, species == i), 
          family=binomial)})
 
names(single_sp_mods) <- 
  sp_codes$spbin[match(names(single_sp_mods), sp_codes$code)]
single_sp_mods$All <- m_cdd_s_ris

map(single_sp_mods, summary)
```

    ## $`Syzygium rubicundum`
    ##  Family: binomial  ( logit )
    ## Formula:          
    ## Pr_s ~ slope.degrees_s + (scale(tot_dens_s) + scale(con_dens_s)) +  
    ##     scale(fragment.size) + trt_I * trt_F + (scale(tot_dens_s) +  
    ##     scale(con_dens_s)) * (trt_I + trt_F) * scale(fragment.size) +  
    ##     (1 | site/loc/gr/plot)
    ## Data: filter(sdls, species == i)
    ## Weights: census.start
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    664.8    740.4   -308.4    616.8      148 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups           Name        Variance  Std.Dev. 
    ##  plot:gr:loc:site (Intercept) 2.103e-01 0.4586119
    ##  gr:loc:site      (Intercept) 1.693e-01 0.4114793
    ##  loc:site         (Intercept) 6.510e-01 0.8068493
    ##  site             (Intercept) 1.922e-07 0.0004384
    ## Number of obs: 172, groups:  
    ## plot:gr:loc:site, 172; gr:loc:site, 61; loc:site, 30; site, 16
    ## 
    ## Conditional model:
    ##                                               Estimate Std. Error z value
    ## (Intercept)                                   -0.14989    0.23065  -0.650
    ## slope.degrees_s                               -0.31491    0.10682  -2.948
    ## scale(tot_dens_s)                              0.04320    0.14091   0.307
    ## scale(con_dens_s)                             -0.65441    0.16022  -4.084
    ## scale(fragment.size)                          -0.19506    0.23153  -0.842
    ## trt_II                                         0.69260    0.25888   2.675
    ## trt_FF                                        -0.01042    0.20949  -0.050
    ## trt_II:trt_FF                                 -0.24727    0.33821  -0.731
    ## scale(tot_dens_s):trt_II                      -0.04539    0.34754  -0.131
    ## scale(tot_dens_s):trt_FF                      -0.05204    0.18634  -0.279
    ## scale(con_dens_s):trt_II                       0.06969    0.24165   0.288
    ## scale(con_dens_s):trt_FF                       0.07628    0.21351   0.357
    ## scale(tot_dens_s):scale(fragment.size)         0.32513    0.12538   2.593
    ## scale(con_dens_s):scale(fragment.size)        -0.08381    0.11571  -0.724
    ## scale(fragment.size):trt_II                   -0.07710    0.24413  -0.316
    ## scale(fragment.size):trt_FF                   -0.28855    0.20967  -1.376
    ## scale(tot_dens_s):scale(fragment.size):trt_II -0.42050    0.33640  -1.250
    ## scale(tot_dens_s):scale(fragment.size):trt_FF -0.26047    0.33499  -0.778
    ## scale(con_dens_s):scale(fragment.size):trt_II -0.01368    0.17338  -0.079
    ## scale(con_dens_s):scale(fragment.size):trt_FF  0.37473    0.15094   2.483
    ##                                               Pr(>|z|)    
    ## (Intercept)                                    0.51578    
    ## slope.degrees_s                                0.00320 ** 
    ## scale(tot_dens_s)                              0.75919    
    ## scale(con_dens_s)                             4.42e-05 ***
    ## scale(fragment.size)                           0.39952    
    ## trt_II                                         0.00746 ** 
    ## trt_FF                                         0.96033    
    ## trt_II:trt_FF                                  0.46470    
    ## scale(tot_dens_s):trt_II                       0.89608    
    ## scale(tot_dens_s):trt_FF                       0.78005    
    ## scale(con_dens_s):trt_II                       0.77306    
    ## scale(con_dens_s):trt_FF                       0.72090    
    ## scale(tot_dens_s):scale(fragment.size)         0.00951 ** 
    ## scale(con_dens_s):scale(fragment.size)         0.46888    
    ## scale(fragment.size):trt_II                    0.75214    
    ## scale(fragment.size):trt_FF                    0.16877    
    ## scale(tot_dens_s):scale(fragment.size):trt_II  0.21130    
    ## scale(tot_dens_s):scale(fragment.size):trt_FF  0.43684    
    ## scale(con_dens_s):scale(fragment.size):trt_II  0.93710    
    ## scale(con_dens_s):scale(fragment.size):trt_FF  0.01304 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $`Symplocos racemosa`
    ##  Family: binomial  ( logit )
    ## Formula:          
    ## Pr_s ~ slope.degrees_s + (scale(tot_dens_s) + scale(con_dens_s)) +  
    ##     scale(fragment.size) + trt_I * trt_F + (scale(tot_dens_s) +  
    ##     scale(con_dens_s)) * (trt_I + trt_F) * scale(fragment.size) +  
    ##     (1 | site/loc/gr/plot)
    ## Data: filter(sdls, species == i)
    ## Weights: census.start
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    328.1    395.4   -140.0    280.1       98 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups           Name        Variance  Std.Dev. 
    ##  plot:gr:loc:site (Intercept) 1.625e-01 0.4030853
    ##  gr:loc:site      (Intercept) 2.503e-02 0.1582100
    ##  loc:site         (Intercept) 5.766e-01 0.7593371
    ##  site             (Intercept) 1.487e-08 0.0001219
    ## Number of obs: 122, groups:  
    ## plot:gr:loc:site, 122; gr:loc:site, 53; loc:site, 24; site, 17
    ## 
    ## Conditional model:
    ##                                               Estimate Std. Error z value
    ## (Intercept)                                    0.71220    0.34523   2.063
    ## slope.degrees_s                                0.19652    0.15944   1.233
    ## scale(tot_dens_s)                              0.30384    0.33337   0.911
    ## scale(con_dens_s)                              0.25698    0.61914   0.415
    ## scale(fragment.size)                          -0.22861    0.43806  -0.522
    ## trt_II                                        -0.26647    0.48631  -0.548
    ## trt_FF                                         1.60922    0.79724   2.018
    ## trt_II:trt_FF                                 -1.12097    0.69201  -1.620
    ## scale(tot_dens_s):trt_II                      -0.04651    0.73245  -0.064
    ## scale(tot_dens_s):trt_FF                       1.53285    1.47131   1.042
    ## scale(con_dens_s):trt_II                      -1.73483    1.47049  -1.180
    ## scale(con_dens_s):trt_FF                       0.89298    2.64710   0.337
    ## scale(tot_dens_s):scale(fragment.size)         0.10147    0.39163   0.259
    ## scale(con_dens_s):scale(fragment.size)         0.69530    0.85389   0.814
    ## scale(fragment.size):trt_II                   -0.33443    0.63470  -0.527
    ## scale(fragment.size):trt_FF                    1.42770    0.91043   1.568
    ## scale(tot_dens_s):scale(fragment.size):trt_II -0.79311    0.95984  -0.826
    ## scale(tot_dens_s):scale(fragment.size):trt_FF  2.43647    2.06836   1.178
    ## scale(con_dens_s):scale(fragment.size):trt_II -2.56759    1.99589  -1.286
    ## scale(con_dens_s):scale(fragment.size):trt_FF  1.12732    3.66319   0.308
    ##                                               Pr(>|z|)  
    ## (Intercept)                                     0.0391 *
    ## slope.degrees_s                                 0.2177  
    ## scale(tot_dens_s)                               0.3621  
    ## scale(con_dens_s)                               0.6781  
    ## scale(fragment.size)                            0.6018  
    ## trt_II                                          0.5837  
    ## trt_FF                                          0.0435 *
    ## trt_II:trt_FF                                   0.1053  
    ## scale(tot_dens_s):trt_II                        0.9494  
    ## scale(tot_dens_s):trt_FF                        0.2975  
    ## scale(con_dens_s):trt_II                        0.2381  
    ## scale(con_dens_s):trt_FF                        0.7359  
    ## scale(tot_dens_s):scale(fragment.size)          0.7956  
    ## scale(con_dens_s):scale(fragment.size)          0.4155  
    ## scale(fragment.size):trt_II                     0.5983  
    ## scale(fragment.size):trt_FF                     0.1168  
    ## scale(tot_dens_s):scale(fragment.size):trt_II   0.4086  
    ## scale(tot_dens_s):scale(fragment.size):trt_FF   0.2388  
    ## scale(con_dens_s):scale(fragment.size):trt_II   0.1983  
    ## scale(con_dens_s):scale(fragment.size):trt_FF   0.7583  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $`Spatholobus purpureus`
    ##  Family: binomial  ( logit )
    ## Formula:          
    ## Pr_s ~ slope.degrees_s + (scale(tot_dens_s) + scale(con_dens_s)) +  
    ##     scale(fragment.size) + trt_I * trt_F + (scale(tot_dens_s) +  
    ##     scale(con_dens_s)) * (trt_I + trt_F) * scale(fragment.size) +  
    ##     (1 | site/loc/gr/plot)
    ## Data: filter(sdls, species == i)
    ## Weights: census.start
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    221.1    277.4    -86.6    173.1       53 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups           Name        Variance  Std.Dev. 
    ##  plot:gr:loc:site (Intercept) 3.501e-11 5.917e-06
    ##  gr:loc:site      (Intercept) 6.978e-12 2.642e-06
    ##  loc:site         (Intercept) 1.120e-10 1.058e-05
    ##  site             (Intercept) 3.560e-11 5.967e-06
    ## Number of obs: 77, groups:  
    ## plot:gr:loc:site, 77; gr:loc:site, 30; loc:site, 15; site, 10
    ## 
    ## Conditional model:
    ##                                               Estimate Std. Error z value
    ## (Intercept)                                   -0.88691    0.24982  -3.550
    ## slope.degrees_s                               -0.40153    0.12585  -3.191
    ## scale(tot_dens_s)                             -1.01314    0.52171  -1.942
    ## scale(con_dens_s)                             -0.03815    0.13710  -0.278
    ## scale(fragment.size)                          -0.36087    0.24130  -1.496
    ## trt_II                                         1.22403    0.48781   2.509
    ## trt_FF                                         0.41325    0.37907   1.090
    ## trt_II:trt_FF                                 -0.31901    0.42278  -0.755
    ## scale(tot_dens_s):trt_II                       0.17701    0.88399   0.200
    ## scale(tot_dens_s):trt_FF                       0.78859    0.65603   1.202
    ## scale(con_dens_s):trt_II                      -0.20655    0.26658  -0.775
    ## scale(con_dens_s):trt_FF                      -0.18489    0.22084  -0.837
    ## scale(tot_dens_s):scale(fragment.size)        -0.77759    0.50864  -1.529
    ## scale(con_dens_s):scale(fragment.size)        -1.55347    0.86887  -1.788
    ## scale(fragment.size):trt_II                   -1.70031    1.44116  -1.180
    ## scale(fragment.size):trt_FF                    1.48540    0.63381   2.344
    ## scale(tot_dens_s):scale(fragment.size):trt_II  0.42592    0.79637   0.535
    ## scale(tot_dens_s):scale(fragment.size):trt_FF  0.40897    0.63428   0.645
    ## scale(con_dens_s):scale(fragment.size):trt_II -5.94422    4.07890  -1.457
    ## scale(con_dens_s):scale(fragment.size):trt_FF  3.21763    1.72268   1.868
    ##                                               Pr(>|z|)    
    ## (Intercept)                                   0.000385 ***
    ## slope.degrees_s                               0.001420 ** 
    ## scale(tot_dens_s)                             0.052141 .  
    ## scale(con_dens_s)                             0.780831    
    ## scale(fragment.size)                          0.134775    
    ## trt_II                                        0.012099 *  
    ## trt_FF                                        0.275639    
    ## trt_II:trt_FF                                 0.450521    
    ## scale(tot_dens_s):trt_II                      0.841293    
    ## scale(tot_dens_s):trt_FF                      0.229340    
    ## scale(con_dens_s):trt_II                      0.438433    
    ## scale(con_dens_s):trt_FF                      0.402457    
    ## scale(tot_dens_s):scale(fragment.size)        0.126322    
    ## scale(con_dens_s):scale(fragment.size)        0.073789 .  
    ## scale(fragment.size):trt_II                   0.238070    
    ## scale(fragment.size):trt_FF                   0.019098 *  
    ## scale(tot_dens_s):scale(fragment.size):trt_II 0.592770    
    ## scale(tot_dens_s):scale(fragment.size):trt_FF 0.519066    
    ## scale(con_dens_s):scale(fragment.size):trt_II 0.145030    
    ## scale(con_dens_s):scale(fragment.size):trt_FF 0.061790 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $`Ventilago madraspatana`
    ##  Family: binomial  ( logit )
    ## Formula:          
    ## Pr_s ~ slope.degrees_s + (scale(tot_dens_s) + scale(con_dens_s)) +  
    ##     scale(fragment.size) + trt_I * trt_F + (scale(tot_dens_s) +  
    ##     scale(con_dens_s)) * (trt_I + trt_F) * scale(fragment.size) +  
    ##     (1 | site/loc/gr/plot)
    ## Data: filter(sdls, species == i)
    ## Weights: census.start
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    140.5    192.4    -46.3     92.5       40 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups           Name        Variance  Std.Dev. 
    ##  plot:gr:loc:site (Intercept) 1.595e-10 1.263e-05
    ##  gr:loc:site      (Intercept) 3.794e-17 6.159e-09
    ##  loc:site         (Intercept) 3.726e-13 6.104e-07
    ##  site             (Intercept) 2.325e-13 4.821e-07
    ## Number of obs: 64, groups:  
    ## plot:gr:loc:site, 64; gr:loc:site, 27; loc:site, 16; site, 10
    ## 
    ## Conditional model:
    ##                                                Estimate Std. Error z value
    ## (Intercept)                                     1.08195    0.43536   2.485
    ## slope.degrees_s                                -0.33295    0.25156  -1.324
    ## scale(tot_dens_s)                               0.77152    0.99666   0.774
    ## scale(con_dens_s)                              -1.10041    0.66212  -1.662
    ## scale(fragment.size)                            0.32157    0.41979   0.766
    ## trt_II                                          5.03424    2.00733   2.508
    ## trt_FF                                          0.01163    0.61168   0.019
    ## trt_II:trt_FF                                  -2.26606    1.09780  -2.064
    ## scale(tot_dens_s):trt_II                       -3.30626    3.59749  -0.919
    ## scale(tot_dens_s):trt_FF                       -0.16767    0.93832  -0.179
    ## scale(con_dens_s):trt_II                       -7.36897    3.06914  -2.401
    ## scale(con_dens_s):trt_FF                        1.17060    1.27183   0.920
    ## scale(tot_dens_s):scale(fragment.size)         -1.60292    0.90233  -1.776
    ## scale(con_dens_s):scale(fragment.size)          1.32854    1.18201   1.124
    ## scale(fragment.size):trt_II                     4.12978    2.49481   1.655
    ## scale(fragment.size):trt_FF                    -1.07212    0.74249  -1.444
    ## scale(tot_dens_s):scale(fragment.size):trt_II  -3.00156    5.18441  -0.579
    ## scale(tot_dens_s):scale(fragment.size):trt_FF   0.20736    1.29481   0.160
    ## scale(con_dens_s):scale(fragment.size):trt_II -10.58223    4.14115  -2.555
    ## scale(con_dens_s):scale(fragment.size):trt_FF   0.54595    1.81312   0.301
    ##                                               Pr(>|z|)  
    ## (Intercept)                                     0.0129 *
    ## slope.degrees_s                                 0.1856  
    ## scale(tot_dens_s)                               0.4389  
    ## scale(con_dens_s)                               0.0965 .
    ## scale(fragment.size)                            0.4437  
    ## trt_II                                          0.0121 *
    ## trt_FF                                          0.9848  
    ## trt_II:trt_FF                                   0.0390 *
    ## scale(tot_dens_s):trt_II                        0.3581  
    ## scale(tot_dens_s):trt_FF                        0.8582  
    ## scale(con_dens_s):trt_II                        0.0164 *
    ## scale(con_dens_s):trt_FF                        0.3574  
    ## scale(tot_dens_s):scale(fragment.size)          0.0757 .
    ## scale(con_dens_s):scale(fragment.size)          0.2610  
    ## scale(fragment.size):trt_II                     0.0979 .
    ## scale(fragment.size):trt_FF                     0.1488  
    ## scale(tot_dens_s):scale(fragment.size):trt_II   0.5626  
    ## scale(tot_dens_s):scale(fragment.size):trt_FF   0.8728  
    ## scale(con_dens_s):scale(fragment.size):trt_II   0.0106 *
    ## scale(con_dens_s):scale(fragment.size):trt_FF   0.7633  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $`Cinnamomum sp.`
    ##  Family: binomial  ( logit )
    ## Formula:          
    ## Pr_s ~ slope.degrees_s + (scale(tot_dens_s) + scale(con_dens_s)) +  
    ##     scale(fragment.size) + trt_I * trt_F + (scale(tot_dens_s) +  
    ##     scale(con_dens_s)) * (trt_I + trt_F) * scale(fragment.size) +  
    ##     (1 | site/loc/gr/plot)
    ## Data: filter(sdls, species == i)
    ## Weights: census.start
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    158.7    217.6    -55.3    110.7       62 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups           Name        Variance  Std.Dev. 
    ##  plot:gr:loc:site (Intercept) 2.994e-09 5.472e-05
    ##  gr:loc:site      (Intercept) 9.293e-01 9.640e-01
    ##  loc:site         (Intercept) 3.086e-09 5.555e-05
    ##  site             (Intercept) 7.390e-01 8.597e-01
    ## Number of obs: 86, groups:  
    ## plot:gr:loc:site, 86; gr:loc:site, 42; loc:site, 23; site, 17
    ## 
    ## Conditional model:
    ##                                               Estimate Std. Error z value
    ## (Intercept)                                   -0.71081    0.57847  -1.229
    ## slope.degrees_s                               -0.21513    0.49865  -0.431
    ## scale(tot_dens_s)                             -0.82516    0.59356  -1.390
    ## scale(con_dens_s)                              0.37420    0.50413   0.742
    ## scale(fragment.size)                          -0.40154    0.59974  -0.669
    ## trt_II                                         0.41253    0.82832   0.498
    ## trt_FF                                         1.92466    1.11891   1.720
    ## trt_II:trt_FF                                 -1.35961    1.49140  -0.912
    ## scale(tot_dens_s):trt_II                      -0.18009    1.13034  -0.159
    ## scale(tot_dens_s):trt_FF                       3.47128    1.77665   1.954
    ## scale(con_dens_s):trt_II                      -0.38512    1.64299  -0.234
    ## scale(con_dens_s):trt_FF                      -0.33156    1.60517  -0.207
    ## scale(tot_dens_s):scale(fragment.size)         0.59999    0.40618   1.477
    ## scale(con_dens_s):scale(fragment.size)         0.36357    0.58390   0.623
    ## scale(fragment.size):trt_II                   -1.45627    1.03267  -1.410
    ## scale(fragment.size):trt_FF                    1.92983    1.33244   1.448
    ## scale(tot_dens_s):scale(fragment.size):trt_II -1.81786    1.65778  -1.097
    ## scale(tot_dens_s):scale(fragment.size):trt_FF  3.45466    2.25046   1.535
    ## scale(con_dens_s):scale(fragment.size):trt_II  0.07232    2.04916   0.035
    ## scale(con_dens_s):scale(fragment.size):trt_FF  0.98437    2.05481   0.479
    ##                                               Pr(>|z|)  
    ## (Intercept)                                     0.2192  
    ## slope.degrees_s                                 0.6662  
    ## scale(tot_dens_s)                               0.1645  
    ## scale(con_dens_s)                               0.4579  
    ## scale(fragment.size)                            0.5032  
    ## trt_II                                          0.6185  
    ## trt_FF                                          0.0854 .
    ## trt_II:trt_FF                                   0.3620  
    ## scale(tot_dens_s):trt_II                        0.8734  
    ## scale(tot_dens_s):trt_FF                        0.0507 .
    ## scale(con_dens_s):trt_II                        0.8147  
    ## scale(con_dens_s):trt_FF                        0.8364  
    ## scale(tot_dens_s):scale(fragment.size)          0.1396  
    ## scale(con_dens_s):scale(fragment.size)          0.5335  
    ## scale(fragment.size):trt_II                     0.1585  
    ## scale(fragment.size):trt_FF                     0.1475  
    ## scale(tot_dens_s):scale(fragment.size):trt_II   0.2728  
    ## scale(tot_dens_s):scale(fragment.size):trt_FF   0.1248  
    ## scale(con_dens_s):scale(fragment.size):trt_II   0.9718  
    ## scale(con_dens_s):scale(fragment.size):trt_FF   0.6319  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $All
    ##  Family: binomial  ( logit )
    ## Formula:          
    ## Pr_s ~ slope.degrees_s + (scale(tot_dens_s) + scale(con_dens_s)) +  
    ##     scale(fragment.size) + trt_I * trt_F + (scale(tot_dens_s) +  
    ##     scale(con_dens_s)) * (trt_I + trt_F) * scale(fragment.size) +  
    ##     (scale(con_dens_s) + scale(tot_dens_s) | species) + (1 |  
    ##     site/loc/gr/plot)
    ## Data: sdls
    ## Weights: census.start
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   2215.0   2361.3  -1077.5   2155.0      938 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups           Name              Variance  Std.Dev.  Corr        
    ##  species          (Intercept)       1.027e+00 1.0131846             
    ##                   scale(con_dens_s) 3.887e-01 0.6234621  0.56       
    ##                   scale(tot_dens_s) 3.738e-02 0.1933384 -0.47 -0.71 
    ##  plot:gr:loc:site (Intercept)       2.788e-01 0.5280003             
    ##  gr:loc:site      (Intercept)       2.057e-01 0.4535170             
    ##  loc:site         (Intercept)       4.782e-01 0.6915250             
    ##  site             (Intercept)       2.188e-08 0.0001479             
    ## Number of obs: 968, groups:  
    ## species, 26; plot:gr:loc:site, 474; gr:loc:site, 110; loc:site, 37; site, 21
    ## 
    ## Conditional model:
    ##                                               Estimate Std. Error z value
    ## (Intercept)                                    0.66723    0.27771   2.403
    ## slope.degrees_s                               -0.18268    0.07133  -2.561
    ## scale(tot_dens_s)                             -0.10921    0.15595  -0.700
    ## scale(con_dens_s)                             -0.44496    0.18820  -2.364
    ## scale(fragment.size)                          -0.25279    0.15514  -1.629
    ## trt_II                                         0.55988    0.15137   3.699
    ## trt_FF                                         0.27661    0.13886   1.992
    ## trt_II:trt_FF                                 -0.20355    0.21738  -0.936
    ## scale(tot_dens_s):trt_II                       0.04971    0.19245   0.258
    ## scale(tot_dens_s):trt_FF                      -0.01093    0.14816  -0.074
    ## scale(con_dens_s):trt_II                      -0.25189    0.16952  -1.486
    ## scale(con_dens_s):trt_FF                       0.38152    0.15970   2.389
    ## scale(tot_dens_s):scale(fragment.size)         0.09861    0.11221   0.879
    ## scale(con_dens_s):scale(fragment.size)        -0.24883    0.16012  -1.554
    ## scale(fragment.size):trt_II                    0.08701    0.12608   0.690
    ## scale(fragment.size):trt_FF                    0.11108    0.11225   0.990
    ## scale(tot_dens_s):scale(fragment.size):trt_II  0.07160    0.22351   0.320
    ## scale(tot_dens_s):scale(fragment.size):trt_FF -0.17654    0.20314  -0.869
    ## scale(con_dens_s):scale(fragment.size):trt_II -0.14412    0.19987  -0.721
    ## scale(con_dens_s):scale(fragment.size):trt_FF  0.55367    0.19600   2.825
    ##                                               Pr(>|z|)    
    ## (Intercept)                                   0.016277 *  
    ## slope.degrees_s                               0.010437 *  
    ## scale(tot_dens_s)                             0.483771    
    ## scale(con_dens_s)                             0.018061 *  
    ## scale(fragment.size)                          0.103231    
    ## trt_II                                        0.000217 ***
    ## trt_FF                                        0.046363 *  
    ## trt_II:trt_FF                                 0.349094    
    ## scale(tot_dens_s):trt_II                      0.796159    
    ## scale(tot_dens_s):trt_FF                      0.941199    
    ## scale(con_dens_s):trt_II                      0.137312    
    ## scale(con_dens_s):trt_FF                      0.016894 *  
    ## scale(tot_dens_s):scale(fragment.size)        0.379537    
    ## scale(con_dens_s):scale(fragment.size)        0.120172    
    ## scale(fragment.size):trt_II                   0.490132    
    ## scale(fragment.size):trt_FF                   0.322352    
    ## scale(tot_dens_s):scale(fragment.size):trt_II 0.748722    
    ## scale(tot_dens_s):scale(fragment.size):trt_FF 0.384811    
    ## scale(con_dens_s):scale(fragment.size):trt_II 0.470860    
    ## scale(con_dens_s):scale(fragment.size):trt_FF 0.004731 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
term_nms <- names(fixef(single_sp_mods[[1]])$cond)

plot_models(single_sp_mods, m.labels=names(single_sp_mods), p.shape=TRUE, 
                        rm.terms = c("slope.degrees_s", 
                         term_nms[str_detect(term_nms, "tot")])) +
  scale_colour_viridis_d(end=0.9, name="Species", option="H")
```

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

![](Kadumane_seedling_CDD_files/figure-gfm/single_sp-1.png)<!-- -->

``` r
## Ventilago distorts scale, so dropping it
plot_models(single_sp_mods[-4], m.labels=names(single_sp_mods)[-4],
            p.shape=TRUE, 
            rm.terms = c("slope.degrees_s", 
                         term_nms[str_detect(term_nms, "tot")])) +
  scale_colour_viridis_d(end=0.9, name="Species", option="H")
```

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

![](Kadumane_seedling_CDD_files/figure-gfm/single_sp-2.png)<!-- -->

``` r
## As convenient as the plot_models function is, it would be nice to 
## plot different species on different scales and simplify the 
## p-values

sp_effects <- map(single_sp_mods, \(d) {
  tidy(d, conf.int = TRUE ) |>  
  filter(effect == "fixed", str_detect(term, "tot", negate=TRUE),
         !(term %in% c("(Intercept)", "slope.degrees_s") )) 
}) |> bind_rows(.id = "Species")

labs <- data.frame(par = c("ConDens","FragArea", "I", "F",  "I:F", 
           "ConDens:I", "ConDens:F", "FragArea:\n ConDens", 
           "FragArea:I" , "FragArea:F", 
           "FragArea:\n (ConDens:I)", "FragArea:\n (ConDens:F)"),
           term = unique(sp_effects$term)) |> 
  mutate(par = factor(par, levels = rev(par)))


sp_effects <- left_join(sp_effects, labs, by = "term")

sp_effects <- mutate(sp_effects, 
                     Biocide = case_when(
                       str_detect(term, "trt_II") ~  "Insecticide",
                       str_detect(term, "trt_FF") ~ "Fungicide",
                       .default =  "Control"),
                     Biocide = factor(
                       ifelse(str_detect(term,"trt_II:trt_FF"), 
                              "Both", Biocide), 
                       levels = c("Control", "Insecticide", 
                                  "Fungicide", "Both")))
sp_effects <- sp_effects |> mutate(sp_effects, 
                                   Species = ifelse(Species == "All", "All", 
                                      paste0("italic('", Species, "')")))

pl_species <-
  filter(sp_effects, Species != "All", par != "I:F") |> 
  ggplot(aes(x = estimate, xmin = conf.low, xmax = conf.high, 
             y = par, colour = Biocide)) +
           facet_wrap(~Species, scales = "free_x", labeller = label_parsed) +
  geom_pointrange(aes(shape = p.value < 0.05 )) +
  geom_vline(xintercept=0, linetype = "dotted") +
  scale_colour_brewer(palette="Set2") +
  scale_shape_manual(values=c(21, 16), guide = "none" ) +
  labs(y = NULL, x = "log(odds ratio)") 

# lemon::reposition_legend(pl_species, position = 'top right',
#                                        panel='panel-3-2')


# ggsave(lemon::reposition_legend(pl_species, position = 'top right',
#                                        panel='panel-3-2'),
#        file = "figures/species_plot.png", height = 7, width = 7)

pl_species <- ((blup_plot) |
                  (pl_species + theme(legend.position = c(0.95, 0.05), 
                                      legend.justification=c(1, 0)))) + 
  plot_layout(widths=c(0.35, 0.65)) + plot_annotation(tag_levels = "A")
```

    ## Warning: A numeric `legend.position` argument in `theme()` was deprecated in ggplot2
    ## 3.5.0.
    ## i Please use the `legend.position.inside` argument of `theme()` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
pl_species
```

![](Kadumane_seedling_CDD_files/figure-gfm/single_sp-3.png)<!-- -->

``` r
ggsave(pl_species,
       file = "figures/species_plot.png", height = 7, width = 9)
```

# Session Information

``` r
sessionInfo()
```

    ## R version 4.4.0 (2024-04-24 ucrt)
    ## Platform: x86_64-w64-mingw32/x64
    ## Running under: Windows 10 x64 (build 17763)
    ## 
    ## Matrix products: default
    ## 
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.1252 
    ## [2] LC_CTYPE=English_United States.1252   
    ## [3] LC_MONETARY=English_United States.1252
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.1252    
    ## 
    ## time zone: America/New_York
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] patchwork_1.2.0     sjPlot_2.8.16       ggdist_3.3.2       
    ##  [4] ggeffects_1.6.0.2   broom.mixed_0.2.9.5 DHARMa_0.4.6       
    ##  [7] glmmTMB_1.1.9       knitr_1.47          ggthemes_5.1.0     
    ## [10] lubridate_1.9.3     forcats_1.0.0       stringr_1.5.1      
    ## [13] dplyr_1.1.4         purrr_1.0.2         readr_2.1.5        
    ## [16] tidyr_1.3.1         tibble_3.2.1        ggplot2_3.5.1      
    ## [19] tidyverse_2.0.0    
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3   rstudioapi_0.16.0    datawizard_0.10.0   
    ##   [4] magrittr_2.0.3       estimability_1.5.1   farver_2.1.2        
    ##   [7] nloptr_2.0.3         rmarkdown_2.27       ragg_1.3.2          
    ##  [10] vctrs_0.6.5          minqa_1.2.7          effectsize_0.8.8    
    ##  [13] htmltools_0.5.8.1    distributional_0.4.0 haven_2.5.4         
    ##  [16] broom_1.0.6          sjmisc_2.8.10        parallelly_1.37.1   
    ##  [19] StanHeaders_2.32.9   plyr_1.8.9           emmeans_1.10.2      
    ##  [22] TMB_1.9.11           mime_0.12            lifecycle_1.0.4     
    ##  [25] iterators_1.0.14     pkgconfig_2.0.3      gap_1.5-3           
    ##  [28] sjlabelled_1.2.0     Matrix_1.7-0         R6_2.5.1            
    ##  [31] fastmap_1.2.0        rbibutils_2.2.16     future_1.33.2       
    ##  [34] shiny_1.8.1.1        digest_0.6.35        numDeriv_2016.8-1.1 
    ##  [37] colorspace_2.1-0     furrr_0.3.1          textshaping_0.4.0   
    ##  [40] qgam_1.3.4           labeling_0.4.3       fansi_1.0.6         
    ##  [43] timechange_0.3.0     mgcv_1.9-1           compiler_4.4.0      
    ##  [46] bit64_4.0.5          withr_3.0.0          doParallel_1.0.17   
    ##  [49] backports_1.5.0      inline_0.3.19        performance_0.11.0  
    ##  [52] QuickJSR_1.2.0       pkgbuild_1.4.4       highr_0.11          
    ##  [55] MASS_7.3-60.2        sjstats_0.19.0       loo_2.7.0           
    ##  [58] tools_4.4.0          beeswarm_0.4.0       httpuv_1.6.15       
    ##  [61] glue_1.7.0           nlme_3.1-164         promises_1.3.0      
    ##  [64] grid_4.4.0           generics_0.1.3       gtable_0.3.5        
    ##  [67] tzdb_0.4.0           hms_1.1.3            utf8_1.2.4          
    ##  [70] foreach_1.5.2        pillar_1.9.0         vroom_1.6.5         
    ##  [73] later_1.3.2          splines_4.4.0        lattice_0.22-6      
    ##  [76] bit_4.0.5            tidyselect_1.2.1     gridExtra_2.3       
    ##  [79] stats4_4.4.0         xfun_0.44            matrixStats_1.3.0   
    ##  [82] rstan_2.32.6         stringi_1.8.4        yaml_2.3.8          
    ##  [85] boot_1.3-30          evaluate_0.23        codetools_0.2-20    
    ##  [88] cli_3.6.2            RcppParallel_5.1.7   xtable_1.8-4        
    ##  [91] parameters_0.21.7    systemfonts_1.1.0    Rdpack_2.6          
    ##  [94] munsell_0.5.1        Rcpp_1.0.12          globals_0.16.3      
    ##  [97] coda_0.19-4.1        parallel_4.4.0       ggh4x_0.2.8         
    ## [100] bayestestR_0.13.2    gap.datasets_0.0.6   lme4_1.1-35.3       
    ## [103] listenv_0.9.1        viridisLite_0.4.2    mvtnorm_1.2-5       
    ## [106] scales_1.3.0         insight_0.19.11      crayon_1.5.2        
    ## [109] rlang_1.1.3
