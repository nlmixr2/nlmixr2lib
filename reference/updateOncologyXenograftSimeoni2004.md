# Update an oncology xenograft model based on Simeoni 2004

Update an oncology xenograft model based on Simeoni 2004

## Usage

``` r
updateOncologyXenograftSimeoni2004(
  object,
  ncmt,
  damagedCmtName = "damaged_cells",
  drugEffectName = "drugEffectCyclingCells",
  undamagedCmtName = "cycling_cells",
  tumorVolName = "tumor_vol",
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
readModelDb("oncology_xenograft_simeoni_2004") |>
  updateOncologyXenograftSimeoni2004(ncmt = 5)
#> ℹ You can modify the number of damaged cell compartments in the model using the function updateOncologyXenograftSimeoni2004(model, ncmt)
#>  
#>  
#> ℹ add covariate `damaged_cells1`
#> ℹ add covariate `damaged_cells2`
#> ℹ add covariate `damaged_cells3`
#> ℹ add covariate `damaged_cells4`
#> ℹ add covariate `damaged_cells5`
#>  ── rxode2-based free-form 6-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>   ldamageTransit       ldrugSlope  ltumorExpGrowth  ltumorLinGrowth 
#>      -0.03252319      -7.37137930      -1.29828348      -0.20579491 
#> propSd_tumor_vol  addSd_tumor_vol 
#>       0.20000000      30.00000000 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1    cycling_cells
#> 2                  2   damaged_cells1
#> 3                  3   damaged_cells2
#> 4                  4   damaged_cells3
#> 5                  5   damaged_cells4
#> 6                  6   damaged_cells5
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     depends <- c("Cc", "tumor_vol0")
#>     description <- "Oncology tumor growth model in xenograft models"
#>     reference <- "Monica Simeoni, Paolo Magni, Cristiano Cammia, Giuseppe De Nicolao, Valter Croci, Enrico Pesenti, Massimiliano Germani, Italo Poggesi, Maurizio Rocchetti; Predictive Pharmacokinetic-Pharmacodynamic Modeling of Tumor Growth Kinetics in Xenograft Models after Administration of Anticancer Agents. Cancer Res 1 February 2004; 64 (3): 1094–1101. https://doi.org/10.1158/0008-5472.CAN-03-2524"
#>     units <- list(time = "day", dosing = "dose_unit", concentration = "conc_unit/vol_unit")
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
#>         propSd_tumor_vol <- c(0, 0.2)
#>         label("Proportional residual error (fraction)")
#>         addSd_tumor_vol <- c(0, 30)
#>         label("Additive residual error (tumor volume)")
#>     })
#>     model({
#>         damageTransit <- exp(ldamageTransit)
#>         drugSlope <- exp(ldrugSlope)
#>         tumorExpGrowth <- exp(ltumorExpGrowth)
#>         tumorLinGrowth <- exp(ltumorLinGrowth)
#>         cycling_cells(0) <- tumor_vol0
#>         psi <- 20
#>         tumor_vol <- cycling_cells + damaged_cells1 + damaged_cells2 + 
#>             damaged_cells3 + damaged_cells4 + damaged_cells5
#>         drugEffectCyclingCells <- drugSlope * Cc
#>         d/dt(cycling_cells) <- tumorExpGrowth * cycling_cells/(1 + 
#>             (tumorExpGrowth/tumorLinGrowth * tumor_vol)^psi)^(1/psi) - 
#>             drugEffectCyclingCells * cycling_cells
#>         tumor_vol ~ prop(propSd_tumor_vol) + add(addSd_tumor_vol)
#>         d/dt(damaged_cells1) <- drugEffectCyclingCells * cycling_cells - 
#>             damageTransit * damaged_cells1
#>         d/dt(damaged_cells2) <- damageTransit * (damaged_cells1 - 
#>             damaged_cells2)
#>         d/dt(damaged_cells3) <- damageTransit * (damaged_cells2 - 
#>             damaged_cells3)
#>         d/dt(damaged_cells4) <- damageTransit * (damaged_cells3 - 
#>             damaged_cells4)
#>         d/dt(damaged_cells5) <- damageTransit * (damaged_cells4 - 
#>             damaged_cells5)
#>     })
#> }
```
