# Add a peripheral compartment to a model

Adds a peripheral distribution compartment (\`peripheral1\`, then
\`peripheral2\`) to a central compartment, following the
\`PK_2cmt_des\`/\`PK_3cmt_des\` seed conventions: intercompartmental
clearance \`q\` (\`q2\` for the second) and peripheral volume \`vp\`
(\`vp2\`), with \`k12 \<- q/vc\`, \`k21 \<- q/vp\` (and \`k13\`/\`k31\`
for the second). At most two peripheral compartments are supported,
matching the seeds.

## Usage

``` r
addPeriph(ui, n = NULL, central = "central", model)
```

## Arguments

- ui:

  The model as a function (or something convertible to an rxUi object)

- n:

  which peripheral to add: 1 (\`peripheral1\`, \`q\`, \`vp\`) or 2
  (\`peripheral2\`, \`q2\`, \`vp2\`); defaults to the first missing one

- central:

  central compartment name

- model:

  Deprecated alias for `ui`. Supplying `model` instead of `ui` still
  works but emits a deprecation warning.

## Value

a model with the peripheral compartment added

## See also

Other distribution:
[`removePeriph()`](https://nlmixr2.github.io/nlmixr2lib/reference/removePeriph.md)

## Author

Matthew L. Fidler

## Examples

``` r

readModelDb("PK_1cmt_des") |> addPeriph()
#>  
#>  
#>  ── rxode2-based free-form 3-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lka    lcl    lvc propSd     lq    lvp 
#>   0.45   1.00   3.45   0.50   0.10   5.00 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2          central
#> 3                  3      peripheral1
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     compartmentData <- list(depot = list(analyte = "drug", units = NA_character_, 
#>         specimen = "administration site", verified = FALSE), 
#>         central = list(analyte = "drug", units = NA_character_, 
#>             specimen = "plasma", verified = FALSE))
#>     dosing <- c("central", "depot")
#>     reference <- "nlmixr2lib template"
#>     units <- list(time = "time_unit", dosing = "dose_unit", concentration = "conc_unit/vol_unit")
#>     ini({
#>         lka <- 0.45
#>         label("Absorption rate (Ka)")
#>         lcl <- 1
#>         label("Clearance (CL)")
#>         lvc <- 3.45
#>         label("Central volume of distribution (V)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         lq <- 0.1
#>         label("Intercompartmental clearance (q)")
#>         lvp <- 5
#>         label("Peripheral volume of distribution (vp)")
#>     })
#>     model({
#>         ka <- exp(lka)
#>         cl <- exp(lcl)
#>         vc <- exp(lvc)
#>         kel <- cl/vc
#>         d/dt(depot) <- -ka * depot
#>         q <- exp(lq)
#>         vp <- exp(lvp)
#>         k12 <- q/vc
#>         k21 <- q/vp
#>         d/dt(central) <- ka * depot - kel * central - k12 * central + 
#>             k21 * peripheral1
#>         d/dt(peripheral1) <- k12 * central - k21 * peripheral1
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }

readModelDb("PK_2cmt_des") |> addPeriph()
#>  
#>  
#>  ── rxode2-based free-form 4-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lka    lcl    lvc    lvp     lq propSd    lq2   lvp2 
#>   0.45   1.00   3.00   5.00   0.10   0.50   0.10   5.00 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2          central
#> 3                  3      peripheral2
#> 4                  4      peripheral1
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     compartmentData <- list(depot = list(analyte = "drug", units = NA_character_, 
#>         specimen = "administration site", verified = FALSE), 
#>         central = list(analyte = "drug", units = NA_character_, 
#>             specimen = "plasma", verified = FALSE), peripheral1 = list(analyte = "drug", 
#>             units = NA_character_, specimen = "plasma", verified = FALSE))
#>     reference <- "nlmixr2lib template"
#>     units <- list(time = "time_unit", dosing = "dose_unit", concentration = "conc_unit/vol_unit")
#>     ini({
#>         lka <- 0.45
#>         label("Absorption rate (Ka)")
#>         lcl <- 1
#>         label("Clearance (CL)")
#>         lvc <- 3
#>         label("Central volume of distribution (V)")
#>         lvp <- 5
#>         label("Peripheral volume of distribution (Vp)")
#>         lq <- 0.1
#>         label("Intercompartmental clearance (Q)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         lq2 <- 0.1
#>         label("Intercompartmental clearance (q2)")
#>         lvp2 <- 5
#>         label("Peripheral volume of distribution (vp2)")
#>     })
#>     model({
#>         ka <- exp(lka)
#>         cl <- exp(lcl)
#>         vc <- exp(lvc)
#>         vp <- exp(lvp)
#>         q <- exp(lq)
#>         kel <- cl/vc
#>         k12 <- q/vc
#>         k21 <- q/vp
#>         d/dt(depot) <- -ka * depot
#>         q2 <- exp(lq2)
#>         vp2 <- exp(lvp2)
#>         k13 <- q2/vc
#>         k31 <- q2/vp2
#>         d/dt(central) <- ka * depot - kel * central - k12 * central + 
#>             k21 * peripheral1 - k13 * central + k31 * peripheral2
#>         d/dt(peripheral2) <- k13 * central - k31 * peripheral2
#>         d/dt(peripheral1) <- k12 * central - k21 * peripheral1
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }
```
