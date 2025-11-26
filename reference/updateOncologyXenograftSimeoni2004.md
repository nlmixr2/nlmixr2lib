# Update an oncology xenograft model based on Simeoni 2004

Update an oncology xenograft model based on Simeoni 2004

## Usage

``` r
updateOncologyXenograftSimeoni2004(
  object,
  ncmt,
  damagedCmtName = "damagedCells",
  drugEffectName = "drugEffectCyclingCells",
  undamagedCmtName = "cyclingCells",
  tumorVolName = "tumorVol",
  transitRateName = "damageTransit"
)
```

## Arguments

- object:

  Fitted object or function specifying the model.

- ncmt:

  The desired number of damaged cell compartments

- damagedCmtName, undamagedCmtName, tumorVolName:

  character string names for the compartments for damaged cells,
  undamaged cells, and the calculated tumor volume (the sum of undamaged
  and damaged cells)

- drugEffectName, transitRateName:

  character string names of the drug effect and transit rate (as used in
  the model block)

## Value

An updated model with the new number of compartments

## Examples

``` r
library(rxode2)
readModelDb("oncology_xenograft_simeoni_2004") %>%
  updateOncologyXenograftSimeoni2004(ncmt = 5)
#> ℹ You can modify the number of damaged cell compartments in the model using the function updateOncologyXenograftSimeoni2004(model, ncmt)
#>  
#>  
#> ℹ add covariate `damagedCells1`
#> ℹ add covariate `damagedCells2`
#> ℹ add covariate `damagedCells3`
#> ℹ add covariate `damagedCells4`
#> ℹ add covariate `damagedCells5`
#>  ── rxode2-based free-form 6-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>  ldamageTransit      ldrugSlope ltumorExpGrowth ltumorLinGrowth  tumorVolpropSd 
#>     -0.03252319     -7.37137930     -1.29828348     -0.20579491      0.20000000 
#>   tumorVoladdSd 
#>     30.00000000 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1     cyclingCells
#> 2                  2    damagedCells1
#> 3                  3    damagedCells2
#> 4                  4    damagedCells3
#> 5                  5    damagedCells4
#> 6                  6    damagedCells5
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     depends <- "Cc"
#>     description <- "Oncology tumor growth model in xenograft models"
#>     reference <- "Monica Simeoni, Paolo Magni, Cristiano Cammia, Giuseppe De Nicolao, Valter Croci, Enrico Pesenti, Massimiliano Germani, Italo Poggesi, Maurizio Rocchetti; Predictive Pharmacokinetic-Pharmacodynamic Modeling of Tumor Growth Kinetics in Xenograft Models after Administration of Anticancer Agents. Cancer Res 1 February 2004; 64 (3): 1094–1101. https://doi.org/10.1158/0008-5472.CAN-03-2524"
#>     units <- list(time = "day")
#>     ini({
#>         ldamageTransit <- c(-2.30258509299405, -0.0325231917055601, 
#>             2.30258509299405)
#>         label("Transit rate through damage (1/day)")
#>         ldrugSlope <- c(-11.5129254649702, -7.37137930126383, 
#>             -2.30258509299405)
#>         label("Linear drug effect on cycling cells (1/(day*ng/mL))")
#>         ltumorExpGrowth <- c(-6.90775527898214, -1.29828348379718, 
#>             0.693147180559945)
#>         label("Tumor exponential growth rate (1/day)")
#>         ltumorLinGrowth <- c(-4.60517018598809, -0.205794912979597, 
#>             1.6094379124341)
#>         label("Tumor linear growth rate (tumor volume/day)")
#>         tumorVolpropSd <- c(0, 0.2)
#>         label("Proportional residual error (fraction)")
#>         tumorVoladdSd <- c(0, 30)
#>         label("Additive residual error (tumor volume)")
#>     })
#>     model({
#>         damageTransit <- exp(ldamageTransit)
#>         drugSlope <- exp(ldrugSlope)
#>         tumorExpGrowth <- exp(ltumorExpGrowth)
#>         tumorLinGrowth <- exp(ltumorLinGrowth)
#>         cyclingCells(0) <- tumorVol0
#>         psi <- 20
#>         tumorVol <- cyclingCells + damagedCells1 + damagedCells2 + 
#>             damagedCells3 + damagedCells4 + damagedCells5
#>         drugEffectCyclingCells <- drugSlope * Cc
#>         d/dt(cyclingCells) <- tumorExpGrowth * cyclingCells/(1 + 
#>             (tumorExpGrowth/tumorLinGrowth * tumorVol)^psi)^(1/psi) - 
#>             drugEffectCyclingCells * cyclingCells
#>         tumorVol ~ prop(tumorVolpropSd) + add(tumorVoladdSd)
#>         d/dt(damagedCells1) <- drugEffectCyclingCells * cyclingCells - 
#>             damageTransit * damagedCells1
#>         d/dt(damagedCells2) <- damageTransit * (damagedCells1 - 
#>             damagedCells2)
#>         d/dt(damagedCells3) <- damageTransit * (damagedCells2 - 
#>             damagedCells3)
#>         d/dt(damagedCells4) <- damageTransit * (damagedCells3 - 
#>             damagedCells4)
#>         d/dt(damagedCells5) <- damageTransit * (damagedCells4 - 
#>             damagedCells5)
#>     })
#> }
```
