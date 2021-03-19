<!-- README.md is generated from README.Rmd. Please edit that file -->



# Cascade <img src="man/figures/logo.png" align="right" width="200"/>

# Cascade, Selection, Reverse-Engineering and Prediction in Cascade Networks
## Frédéric Bertrand and Myriam Maumy-Bertrand

<!-- badges: start -->
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/fbertran/Cascade/workflows/R-CMD-check/badge.svg)](https://github.com/fbertran/Cascade/actions)
[![Codecov test coverage](https://codecov.io/gh/fbertran/Cascade/branch/master/graph/badge.svg)](https://codecov.io/gh/fbertran/Cascade?branch=master)
[![CRAN status](https://www.r-pkg.org/badges/version/Cascade)](https://cran.r-project.org/package=Cascade)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/Cascade)](https://cran.r-project.org/package=Cascade)
[![GitHub Repo stars](https://img.shields.io/github/stars/fbertran/Cascade?style=social)](https://github.com/fbertran/Cascade)
[![DOI](https://zenodo.org/badge/18441815.svg)](https://zenodo.org/badge/latestdoi/18441815)
<!-- badges: end -->

Cascade is a modeling tool allowing gene selection, reverse engineering, and prediction in cascade networks. Jung, N., Bertrand, F., Bahram, S., Vallat, L., and Maumy-Bertrand, M. (2014) <https://doi.org/10.1093/bioinformatics/btt705>.


The package was presented at the [User2014!](http://user2014.r-project.org/) conference. Jung, N., Bertrand, F., Bahram, S., Vallat, L., and Maumy-Bertrand, M. (2014). "Cascade: a R-package to study, predict and simulate the diffusion of a signal through a temporal genenetwork", *book of abstracts*, User2014!, Los Angeles, page 153, <http://user2014.r-project.org/abstracts/posters/181_Jung.pdf>.


![Reverse-engineered network.](man/figures/figure_article_def.png)

![Simulation of an intervention on a gene of the network.](man/figures/figure_network_pred2.png)

This website and these examples were created by F. Bertrand and M. Maumy-Bertrand.

## Installation

You can install the released version of Cascade from [CRAN](https://CRAN.R-project.org) with:


```r
install.packages("Cascade")
```

You can install the development version of Cascade from [github](https://github.com) with:


```r
devtools::install_github("fbertran/Cascade")
```

## Examples

### Data management
Import Cascade Data (repeated measurements on several subjects) from the CascadeData package and turn them into a micro array object. The second line makes sure the CascadeData package is installed.

```r
library(Cascade)
if(!require(CascadeData)){install.packages("CascadeData")}
data(micro_US)
micro_US<-as.micro_array(micro_US,time=c(60,90,210,390),subject=6)
```

Get a summay and plots of the data:

```r
summary(micro_US)
#>    N1_US_T60        N1_US_T90        N1_US_T210    
#>  Min.   :   1.0   Min.   :   1.0   Min.   :   1.0  
#>  1st Qu.:  19.7   1st Qu.:  18.8   1st Qu.:  15.2  
#>  Median :  38.0   Median :  37.2   Median :  34.9  
#>  Mean   : 107.5   Mean   : 106.9   Mean   : 109.6  
#>  3rd Qu.:  80.6   3rd Qu.:  82.1   3rd Qu.:  82.8  
#>  Max.   :8587.9   Max.   :8311.7   Max.   :7930.3  
#>    N1_US_T390       N2_US_T60        N2_US_T90     
#>  Min.   :   1.0   Min.   :   1.0   Min.   :   1.0  
#>  1st Qu.:  20.9   1st Qu.:  18.5   1st Qu.:  17.1  
#>  Median :  40.2   Median :  36.9   Median :  36.7  
#>  Mean   : 105.7   Mean   : 110.6   Mean   : 102.1  
#>  3rd Qu.:  84.8   3rd Qu.:  85.3   3rd Qu.:  78.2  
#>  Max.   :7841.8   Max.   :7750.3   Max.   :8014.3  
#>    N2_US_T210       N2_US_T390       N3_US_T60     
#>  Min.   :   1.0   Min.   :   1.0   Min.   :   1.0  
#>  1st Qu.:  15.8   1st Qu.:  17.7   1st Qu.:  17.3  
#>  Median :  36.0   Median :  37.4   Median :  34.4  
#>  Mean   : 106.8   Mean   : 111.3   Mean   : 101.6  
#>  3rd Qu.:  83.5   3rd Qu.:  86.4   3rd Qu.:  75.4  
#>  Max.   :8028.6   Max.   :7498.4   Max.   :8072.2  
#>    N3_US_T90        N3_US_T210       N3_US_T390    
#>  Min.   :   1.0   Min.   :   1.0   Min.   :   1.0  
#>  1st Qu.:  19.5   1st Qu.:  16.4   1st Qu.:  20.9  
#>  Median :  38.2   Median :  34.7   Median :  41.0  
#>  Mean   : 107.1   Mean   : 100.3   Mean   : 113.9  
#>  3rd Qu.:  82.3   3rd Qu.:  76.3   3rd Qu.:  89.2  
#>  Max.   :7889.2   Max.   :8278.2   Max.   :6856.2  
#>    N4_US_T60        N4_US_T90        N4_US_T210    
#>  Min.   :   1.0   Min.   :   1.0   Min.   :   1.0  
#>  1st Qu.:  20.4   1st Qu.:  19.5   1st Qu.:  20.5  
#>  Median :  38.9   Median :  38.5   Median :  39.9  
#>  Mean   : 113.6   Mean   : 114.8   Mean   : 110.1  
#>  3rd Qu.:  84.6   3rd Qu.:  86.1   3rd Qu.:  86.8  
#>  Max.   :9502.3   Max.   :9193.4   Max.   :9436.0  
#>    N4_US_T390       N5_US_T60        N5_US_T90     
#>  Min.   :   1.0   Min.   :   1.0   Min.   :   1.0  
#>  1st Qu.:  19.9   1st Qu.:  16.8   1st Qu.:  18.8  
#>  Median :  38.8   Median :  34.5   Median :  36.9  
#>  Mean   : 111.7   Mean   : 111.3   Mean   : 108.0  
#>  3rd Qu.:  85.4   3rd Qu.:  82.0   3rd Qu.:  81.4  
#>  Max.   :8771.0   Max.   :8569.3   Max.   :7970.1  
#>    N5_US_T210       N5_US_T390       N6_US_T60     
#>  Min.   :   1.0   Min.   :   1.0   Min.   :   1.0  
#>  1st Qu.:  19.5   1st Qu.:  19.9   1st Qu.:  21.1  
#>  Median :  38.2   Median :  39.0   Median :  40.9  
#>  Mean   : 107.4   Mean   : 109.8   Mean   : 110.1  
#>  3rd Qu.:  82.4   3rd Qu.:  84.9   3rd Qu.:  86.3  
#>  Max.   :8371.0   Max.   :7686.5   Max.   :8241.0  
#>    N6_US_T90        N6_US_T210       N6_US_T390    
#>  Min.   :   1.0   Min.   :   1.0   Min.   :   1.0  
#>  1st Qu.:  21.5   1st Qu.:  19.9   1st Qu.:  20.2  
#>  Median :  40.8   Median :  39.1   Median :  39.4  
#>  Mean   : 108.5   Mean   : 112.0   Mean   : 109.5  
#>  3rd Qu.:  85.6   3rd Qu.:  86.3   3rd Qu.:  86.6  
#>  Max.   :8355.0   Max.   :8207.1   Max.   :9520.0
```

<img src="man/figures/README-plotmicroarrayclass-1.png" title="plot of chunk plotmicroarrayclass" alt="plot of chunk plotmicroarrayclass" width="100%" /><img src="man/figures/README-plotmicroarrayclass-2.png" title="plot of chunk plotmicroarrayclass" alt="plot of chunk plotmicroarrayclass" width="100%" />

### Gene selection
There are several functions to carry out gene selection before the inference. They are detailed in the two vignettes of the package. 

### Data simulation
Let's simulate some cascade data and then do some reverse engineering.

We first design the F matrix

```r
T<-4
F<-array(0,c(T-1,T-1,T*(T-1)/2))

for(i in 1:(T*(T-1)/2)){diag(F[,,i])<-1}
F[,,2]<-F[,,2]*0.2
F[2,1,2]<-1
F[3,2,2]<-1
F[,,4]<-F[,,2]*0.3
F[3,1,4]<-1
F[,,5]<-F[,,2]
```

We set the seed to make the results reproducible and draw a scale free random network.

```r
set.seed(1)
Net<-Cascade::network_random(
  nb=100,
  time_label=rep(1:4,each=25),
  exp=1,
  init=1,
  regul=round(rexp(100,1))+1,
  min_expr=0.1,
  max_expr=2,
  casc.level=0.4
)
Net@F<-F
```

We simulate gene expression according to the network that was previously drawn

```r
M <- Cascade::gene_expr_simulation(
  network=Net,
  time_label=rep(1:4,each=25),
  subject=5,
  level_peak=200)
```

Get a summay and plots of the simulated data:

```r
summary(M)
#>  log(S/US) : P1T1   log(S/US) : P1T2     log(S/US) : P1T3   
#>  Min.   :-759.882   Min.   :-2024.5979   Min.   :-1007.748  
#>  1st Qu.: -36.758   1st Qu.:  -22.5653   1st Qu.:  -68.054  
#>  Median :   6.265   Median :    0.5759   Median :   -4.192  
#>  Mean   :  10.613   Mean   :   -6.6230   Mean   :    3.085  
#>  3rd Qu.:  74.682   3rd Qu.:   78.2516   3rd Qu.:   74.866  
#>  Max.   : 647.643   Max.   :  870.7513   Max.   : 1155.413  
#>  log(S/US) : P1T4    log(S/US) : P2T1   log(S/US) : P2T2   
#>  Min.   :-1075.637   Min.   :-790.431   Min.   :-1505.543  
#>  1st Qu.:  -31.538   1st Qu.: -65.394   1st Qu.:  -59.833  
#>  Median :   -2.293   Median :   2.087   Median :   -1.262  
#>  Mean   :    9.055   Mean   :   7.791   Mean   :  -18.568  
#>  3rd Qu.:   75.316   3rd Qu.:  70.108   3rd Qu.:   76.908  
#>  Max.   :  556.449   Max.   : 669.203   Max.   : 1058.385  
#>  log(S/US) : P2T3   log(S/US) : P2T4   log(S/US) : P3T1    
#>  Min.   :-980.965   Min.   :-547.117   Min.   :-1278.6158  
#>  1st Qu.: -55.077   1st Qu.: -58.721   1st Qu.:  -42.1909  
#>  Median :  -7.144   Median :  -3.519   Median :    0.4064  
#>  Mean   : -35.647   Mean   : -24.277   Mean   :   -3.8860  
#>  3rd Qu.:  41.364   3rd Qu.:  37.621   3rd Qu.:   48.4275  
#>  Max.   :1114.897   Max.   : 270.423   Max.   :  527.3972  
#>  log(S/US) : P3T2   log(S/US) : P3T3    log(S/US) : P3T4   
#>  Min.   :-624.834   Min.   :-1018.897   Min.   :-2403.703  
#>  1st Qu.: -48.460   1st Qu.:  -52.456   1st Qu.:  -57.493  
#>  Median :  -2.505   Median :   -2.026   Median :   -4.529  
#>  Mean   : -11.950   Mean   :    5.893   Mean   :  -33.088  
#>  3rd Qu.:  33.781   3rd Qu.:   43.214   3rd Qu.:   51.808  
#>  Max.   : 576.141   Max.   : 1159.517   Max.   :  495.014  
#>  log(S/US) : P4T1   log(S/US) : P4T2     log(S/US) : P4T3  
#>  Min.   :-683.000   Min.   :-1957.3692   Min.   :-591.460  
#>  1st Qu.: -81.510   1st Qu.:  -39.3372   1st Qu.: -39.499  
#>  Median :   5.102   Median :   -0.0215   Median :   2.251  
#>  Mean   :  -2.034   Mean   :    2.9550   Mean   :  27.377  
#>  3rd Qu.:  74.738   3rd Qu.:   77.1869   3rd Qu.:  62.097  
#>  Max.   : 454.955   Max.   :  955.0680   Max.   :1341.859  
#>  log(S/US) : P4T4   log(S/US) : P5T1   log(S/US) : P5T2   
#>  Min.   :-577.069   Min.   :-436.986   Min.   :-647.1962  
#>  1st Qu.: -31.924   1st Qu.: -69.809   1st Qu.: -48.5156  
#>  Median :   2.456   Median :   2.156   Median :  -0.2949  
#>  Mean   :  29.675   Mean   :  -2.929   Mean   :   6.2300  
#>  3rd Qu.:  35.322   3rd Qu.:  47.462   3rd Qu.:  47.9558  
#>  Max.   :1577.042   Max.   : 651.596   Max.   :1359.9584  
#>  log(S/US) : P5T3   log(S/US) : P5T4  
#>  Min.   :-409.347   Min.   :-188.652  
#>  1st Qu.: -44.205   1st Qu.: -33.097  
#>  Median :  -1.056   Median :   1.873  
#>  Mean   :   6.190   Mean   :  26.740  
#>  3rd Qu.:  47.694   3rd Qu.:  67.879  
#>  Max.   : 434.178   Max.   : 743.820
```

<img src="man/figures/README-summarysimuldata-1.png" title="plot of chunk summarysimuldata" alt="plot of chunk summarysimuldata" width="100%" /><img src="man/figures/README-summarysimuldata-2.png" title="plot of chunk summarysimuldata" alt="plot of chunk summarysimuldata" width="100%" /><img src="man/figures/README-summarysimuldata-3.png" title="plot of chunk summarysimuldata" alt="plot of chunk summarysimuldata" width="100%" />


```r
plot(M)
```

<img src="man/figures/README-plotsimuldata-1.png" title="plot of chunk plotsimuldata" alt="plot of chunk plotsimuldata" width="100%" /><img src="man/figures/README-plotsimuldata-2.png" title="plot of chunk plotsimuldata" alt="plot of chunk plotsimuldata" width="100%" /><img src="man/figures/README-plotsimuldata-3.png" title="plot of chunk plotsimuldata" alt="plot of chunk plotsimuldata" width="100%" /><img src="man/figures/README-plotsimuldata-4.png" title="plot of chunk plotsimuldata" alt="plot of chunk plotsimuldata" width="100%" /><img src="man/figures/README-plotsimuldata-5.png" title="plot of chunk plotsimuldata" alt="plot of chunk plotsimuldata" width="100%" /><img src="man/figures/README-plotsimuldata-6.png" title="plot of chunk plotsimuldata" alt="plot of chunk plotsimuldata" width="100%" /><img src="man/figures/README-plotsimuldata-7.png" title="plot of chunk plotsimuldata" alt="plot of chunk plotsimuldata" width="100%" />

### Network inference
We infer the new network using subjectwise leave one out cross-validation (all measurement from the same subject are removed from the dataset)

```r
Net_inf_C <- Cascade::inference(M, cv.subjects=TRUE)
#> We are at step :  1
#> The convergence of the network is (L1 norm) : 0.0068
#> We are at step :  2
#> The convergence of the network is (L1 norm) : 0.00121
#> We are at step :  3
#> The convergence of the network is (L1 norm) : 0.00096
```

<img src="man/figures/README-netinf-1.png" title="plot of chunk netinf" alt="plot of chunk netinf" width="100%" /><img src="man/figures/README-netinf-2.png" title="plot of chunk netinf" alt="plot of chunk netinf" width="100%" />

Heatmap of the coefficients of the Omega matrix of the network

```r
stats::heatmap(Net_inf_C@network, Rowv=NA, Colv=NA, scale="none", revC=TRUE)
```

<img src="man/figures/README-heatresults-1.png" title="plot of chunk heatresults" alt="plot of chunk heatresults" width="100%" />

###Post inferrence network analysis
We switch to data that were derived from the inferrence of a real biological network and try to detect the optimal cutoff value: the best cutoff value for a network to fit a scale free network.


```r
data("network")
set.seed(1)
cutoff(network)
#> [1] "This calculation may be long"
#> [1] "1/10"
#> [1] "2/10"
#> [1] "3/10"
#> [1] "4/10"
#> [1] "5/10"
#> [1] "6/10"
#> [1] "7/10"
#> [1] "8/10"
#> [1] "9/10"
#> [1] "10/10"
#>  [1] 0.000 0.000 0.569 0.295 0.242 0.147 0.012 0.387 0.351 0.096
```

<img src="man/figures/README-cutoff-1.png" title="plot of chunk cutoff" alt="plot of chunk cutoff" width="100%" />

```
#> $p.value
#>  [1] 0.000 0.000 0.569 0.295 0.242 0.147 0.012 0.387 0.351 0.096
#> 
#> $p.value.inter
#>  [1] -0.0655747  0.2027477  0.3488525  0.3976355  0.2215038  0.1046155
#>  [7]  0.1257731  0.2560059  0.2561465  0.1560364
#> 
#> $sequence
#>  [1] 0.00000000 0.04444444 0.08888889 0.13333333 0.17777778 0.22222222
#>  [7] 0.26666667 0.31111111 0.35555556 0.40000000
```

Analyze the network with a cutoff set to the previouly found 0.14 optimal value.


```r
analyze_network(network,nv=0.14)
#> Error in h(simpleError(msg, call)): erreur d'évaluation de l'argument 'Omega' lors de la sélection d'une méthode pour la fonction 'analyze_network' : objet 'network' introuvable
```


