# Add a property to a compartment

Add a property to a compartment

## Usage

``` r
addCmtProp(ui, prop = c("f", "lag", "dur", "rate", "ini"), cmt)

addBioavailability(ui, cmt, scale = c("log", "logit"), cmt2 = NULL, f = 0.8)

addLag(ui, cmt)

addDur(ui, cmt)

addRate(ui, cmt)

addIni(ui, cmt)
```

## Arguments

- ui:

  rxode2 ui object

- prop:

  property to add to a compartment:

  \- `F`: bioavailability

  \- `lag`: absorption lag time

  \- `dur`: modeled duration of infusion

  \- `rate`: modeled infusion rate

  \- `ini`: initial value of the compartment

- cmt:

  compartment to apply the property to

- scale:

  parameterization of the bioavailability: `"log"` (default) estimates
  it unboundedly on the log scale, while `"logit"` constrains it to
  (0,1) with a logit/expit parameterization (Monolix-style `F`, oral or
  subcutaneous)

- cmt2:

  optional second compartment for the logit dose-split form; see
  \[addBioavailabilityLogit()\]

- f:

  initial bioavailability fraction, in (0,1), for the logit form, or
  NULL to leave the initial estimate unset

## Value

rxode2 ui object with property applied

## Functions

- `addBioavailability()`: Adds the bioavailability to a compartment in
  the model

- `addLag()`: Adds the lag-time to a compartment in the model

- `addDur()`: Adds the modeled duration to a compartment in the model

- `addRate()`: Adds the modeled rate to a compartment in the model

- `addIni()`: Adds the initial value to the compartment

## Author

Matthew L. Fidler

## Examples

``` r

readModelDb("PK_3cmt_des") |> addCmtProp("f", "depot")
#>  
#>  
#>  ── rxode2-based free-form 4-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>     lka     lcl     lvc     lvp    lvp2      lq     lq2  propSd lfDepot 
#>    0.45    1.00    3.00    5.00    8.00    0.10    0.50    0.50    0.10 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2          central
#> 3                  3      peripheral1
#> 4                  4      peripheral2
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     compartmentData <- list(depot = list(analyte = "drug", units = NA_character_, 
#>         specimen = "administration site", verified = FALSE), 
#>         central = list(analyte = "drug", units = NA_character_, 
#>             specimen = "plasma", verified = FALSE), peripheral1 = list(analyte = "drug", 
#>             units = NA_character_, specimen = "plasma", verified = FALSE), 
#>         peripheral2 = list(analyte = "drug", units = NA_character_, 
#>             specimen = "plasma", verified = FALSE))
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
#>         lvp2 <- 8
#>         label("Second peripheral volume of distribution (Vp2)")
#>         lq <- 0.1
#>         label("Intercompartmental clearance (Q)")
#>         lq2 <- 0.5
#>         label("Second intercompartmental clearance (Q2)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         lfDepot <- 0.1
#>     })
#>     model({
#>         fDepot <- exp(lfDepot)
#>         ka <- exp(lka)
#>         cl <- exp(lcl)
#>         vc <- exp(lvc)
#>         vp <- exp(lvp)
#>         vp2 <- exp(lvp2)
#>         q <- exp(lq)
#>         q2 <- exp(lq2)
#>         kel <- cl/vc
#>         k12 <- q/vc
#>         k21 <- q/vp
#>         k13 <- q2/vc
#>         k31 <- q2/vp2
#>         d/dt(depot) <- -ka * depot
#>         f(depot) <- fDepot
#>         d/dt(central) <- ka * depot - kel * central - k12 * central + 
#>             k21 * peripheral1 - k13 * central + k31 * peripheral2
#>         d/dt(peripheral1) <- k12 * central - k21 * peripheral1
#>         d/dt(peripheral2) <- k13 * central - k31 * peripheral2
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }

readModelDb("PK_3cmt_des") |> addBioavailability(depot)
#>  
#>  
#>  
#>  
#>  ── rxode2-based free-form 4-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>        lka        lcl        lvc        lvp       lvp2         lq        lq2 
#>  0.4500000  1.0000000  3.0000000  5.0000000  8.0000000  0.1000000  0.5000000 
#>     propSd    lfDepot 
#>  0.5000000 -0.2231436 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2          central
#> 3                  3      peripheral1
#> 4                  4      peripheral2
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     compartmentData <- list(depot = list(analyte = "drug", units = NA_character_, 
#>         specimen = "administration site", verified = FALSE), 
#>         central = list(analyte = "drug", units = NA_character_, 
#>             specimen = "plasma", verified = FALSE), peripheral1 = list(analyte = "drug", 
#>             units = NA_character_, specimen = "plasma", verified = FALSE), 
#>         peripheral2 = list(analyte = "drug", units = NA_character_, 
#>             specimen = "plasma", verified = FALSE))
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
#>         lvp2 <- 8
#>         label("Second peripheral volume of distribution (Vp2)")
#>         lq <- 0.1
#>         label("Intercompartmental clearance (Q)")
#>         lq2 <- 0.5
#>         label("Second intercompartmental clearance (Q2)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         lfDepot <- -0.22314355131421
#>     })
#>     model({
#>         fDepot <- exp(lfDepot)
#>         ka <- exp(lka)
#>         cl <- exp(lcl)
#>         vc <- exp(lvc)
#>         vp <- exp(lvp)
#>         vp2 <- exp(lvp2)
#>         q <- exp(lq)
#>         q2 <- exp(lq2)
#>         kel <- cl/vc
#>         k12 <- q/vc
#>         k21 <- q/vp
#>         k13 <- q2/vc
#>         k31 <- q2/vp2
#>         d/dt(depot) <- -ka * depot
#>         f(depot) <- fDepot
#>         d/dt(central) <- ka * depot - kel * central - k12 * central + 
#>             k21 * peripheral1 - k13 * central + k31 * peripheral2
#>         d/dt(peripheral1) <- k12 * central - k21 * peripheral1
#>         d/dt(peripheral2) <- k13 * central - k31 * peripheral2
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }

readModelDb("PK_3cmt_des") |> addLag(depot)
#>  
#>  
#>  ── rxode2-based free-form 4-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>       lka       lcl       lvc       lvp      lvp2        lq       lq2    propSd 
#>      0.45      1.00      3.00      5.00      8.00      0.10      0.50      0.50 
#> llagDepot 
#>      0.10 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2          central
#> 3                  3      peripheral1
#> 4                  4      peripheral2
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     compartmentData <- list(depot = list(analyte = "drug", units = NA_character_, 
#>         specimen = "administration site", verified = FALSE), 
#>         central = list(analyte = "drug", units = NA_character_, 
#>             specimen = "plasma", verified = FALSE), peripheral1 = list(analyte = "drug", 
#>             units = NA_character_, specimen = "plasma", verified = FALSE), 
#>         peripheral2 = list(analyte = "drug", units = NA_character_, 
#>             specimen = "plasma", verified = FALSE))
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
#>         lvp2 <- 8
#>         label("Second peripheral volume of distribution (Vp2)")
#>         lq <- 0.1
#>         label("Intercompartmental clearance (Q)")
#>         lq2 <- 0.5
#>         label("Second intercompartmental clearance (Q2)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         llagDepot <- 0.1
#>     })
#>     model({
#>         lagDepot <- exp(llagDepot)
#>         ka <- exp(lka)
#>         cl <- exp(lcl)
#>         vc <- exp(lvc)
#>         vp <- exp(lvp)
#>         vp2 <- exp(lvp2)
#>         q <- exp(lq)
#>         q2 <- exp(lq2)
#>         kel <- cl/vc
#>         k12 <- q/vc
#>         k21 <- q/vp
#>         k13 <- q2/vc
#>         k31 <- q2/vp2
#>         d/dt(depot) <- -ka * depot
#>         lag(depot) <- lagDepot
#>         d/dt(central) <- ka * depot - kel * central - k12 * central + 
#>             k21 * peripheral1 - k13 * central + k31 * peripheral2
#>         d/dt(peripheral1) <- k12 * central - k21 * peripheral1
#>         d/dt(peripheral2) <- k13 * central - k31 * peripheral2
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }

readModelDb("PK_3cmt_des") |> addDur(depot)
#>  
#>  
#>  ── rxode2-based free-form 4-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>       lka       lcl       lvc       lvp      lvp2        lq       lq2    propSd 
#>      0.45      1.00      3.00      5.00      8.00      0.10      0.50      0.50 
#> ldurDepot 
#>      0.10 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2          central
#> 3                  3      peripheral1
#> 4                  4      peripheral2
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     compartmentData <- list(depot = list(analyte = "drug", units = NA_character_, 
#>         specimen = "administration site", verified = FALSE), 
#>         central = list(analyte = "drug", units = NA_character_, 
#>             specimen = "plasma", verified = FALSE), peripheral1 = list(analyte = "drug", 
#>             units = NA_character_, specimen = "plasma", verified = FALSE), 
#>         peripheral2 = list(analyte = "drug", units = NA_character_, 
#>             specimen = "plasma", verified = FALSE))
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
#>         lvp2 <- 8
#>         label("Second peripheral volume of distribution (Vp2)")
#>         lq <- 0.1
#>         label("Intercompartmental clearance (Q)")
#>         lq2 <- 0.5
#>         label("Second intercompartmental clearance (Q2)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         ldurDepot <- 0.1
#>     })
#>     model({
#>         durDepot <- exp(ldurDepot)
#>         ka <- exp(lka)
#>         cl <- exp(lcl)
#>         vc <- exp(lvc)
#>         vp <- exp(lvp)
#>         vp2 <- exp(lvp2)
#>         q <- exp(lq)
#>         q2 <- exp(lq2)
#>         kel <- cl/vc
#>         k12 <- q/vc
#>         k21 <- q/vp
#>         k13 <- q2/vc
#>         k31 <- q2/vp2
#>         d/dt(depot) <- -ka * depot
#>         dur(depot) <- durDepot
#>         d/dt(central) <- ka * depot - kel * central - k12 * central + 
#>             k21 * peripheral1 - k13 * central + k31 * peripheral2
#>         d/dt(peripheral1) <- k12 * central - k21 * peripheral1
#>         d/dt(peripheral2) <- k13 * central - k31 * peripheral2
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }

readModelDb("PK_3cmt_des") |> addRate(depot)
#>  
#>  
#>  ── rxode2-based free-form 4-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>        lka        lcl        lvc        lvp       lvp2         lq        lq2 
#>       0.45       1.00       3.00       5.00       8.00       0.10       0.50 
#>     propSd lrateDepot 
#>       0.50       0.10 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2          central
#> 3                  3      peripheral1
#> 4                  4      peripheral2
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     compartmentData <- list(depot = list(analyte = "drug", units = NA_character_, 
#>         specimen = "administration site", verified = FALSE), 
#>         central = list(analyte = "drug", units = NA_character_, 
#>             specimen = "plasma", verified = FALSE), peripheral1 = list(analyte = "drug", 
#>             units = NA_character_, specimen = "plasma", verified = FALSE), 
#>         peripheral2 = list(analyte = "drug", units = NA_character_, 
#>             specimen = "plasma", verified = FALSE))
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
#>         lvp2 <- 8
#>         label("Second peripheral volume of distribution (Vp2)")
#>         lq <- 0.1
#>         label("Intercompartmental clearance (Q)")
#>         lq2 <- 0.5
#>         label("Second intercompartmental clearance (Q2)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         lrateDepot <- 0.1
#>     })
#>     model({
#>         rateDepot <- exp(lrateDepot)
#>         ka <- exp(lka)
#>         cl <- exp(lcl)
#>         vc <- exp(lvc)
#>         vp <- exp(lvp)
#>         vp2 <- exp(lvp2)
#>         q <- exp(lq)
#>         q2 <- exp(lq2)
#>         kel <- cl/vc
#>         k12 <- q/vc
#>         k21 <- q/vp
#>         k13 <- q2/vc
#>         k31 <- q2/vp2
#>         d/dt(depot) <- -ka * depot
#>         rate(depot) <- rateDepot
#>         d/dt(central) <- ka * depot - kel * central - k12 * central + 
#>             k21 * peripheral1 - k13 * central + k31 * peripheral2
#>         d/dt(peripheral1) <- k12 * central - k21 * peripheral1
#>         d/dt(peripheral2) <- k13 * central - k31 * peripheral2
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }

readModelDb("PK_3cmt_des") |> addIni(depot)
#>  
#>  
#>  ── rxode2-based free-form 4-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>       lka       lcl       lvc       lvp      lvp2        lq       lq2    propSd 
#>      0.45      1.00      3.00      5.00      8.00      0.10      0.50      0.50 
#> liniDepot 
#>      0.10 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2          central
#> 3                  3      peripheral1
#> 4                  4      peripheral2
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     compartmentData <- list(depot = list(analyte = "drug", units = NA_character_, 
#>         specimen = "administration site", verified = FALSE), 
#>         central = list(analyte = "drug", units = NA_character_, 
#>             specimen = "plasma", verified = FALSE), peripheral1 = list(analyte = "drug", 
#>             units = NA_character_, specimen = "plasma", verified = FALSE), 
#>         peripheral2 = list(analyte = "drug", units = NA_character_, 
#>             specimen = "plasma", verified = FALSE))
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
#>         lvp2 <- 8
#>         label("Second peripheral volume of distribution (Vp2)")
#>         lq <- 0.1
#>         label("Intercompartmental clearance (Q)")
#>         lq2 <- 0.5
#>         label("Second intercompartmental clearance (Q2)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         liniDepot <- 0.1
#>     })
#>     model({
#>         iniDepot <- exp(liniDepot)
#>         ka <- exp(lka)
#>         cl <- exp(lcl)
#>         vc <- exp(lvc)
#>         vp <- exp(lvp)
#>         vp2 <- exp(lvp2)
#>         q <- exp(lq)
#>         q2 <- exp(lq2)
#>         kel <- cl/vc
#>         k12 <- q/vc
#>         k21 <- q/vp
#>         k13 <- q2/vc
#>         k31 <- q2/vp2
#>         d/dt(depot) <- -ka * depot
#>         depot(0) <- iniDepot
#>         d/dt(central) <- ka * depot - kel * central - k12 * central + 
#>             k21 * peripheral1 - k13 * central + k31 * peripheral2
#>         d/dt(peripheral1) <- k12 * central - k21 * peripheral1
#>         d/dt(peripheral2) <- k13 * central - k31 * peripheral2
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }

# bioavailability bounded to (0,1) with a logit parameterization
readModelDb("PK_1cmt_des") |> addBioavailability(depot, scale = "logit")
#>  
#>  
#>  ── rxode2-based free-form 2-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>         lka         lcl         lvc      propSd logitfDepot 
#>    0.450000    1.000000    3.450000    0.500000    1.386294 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2          central
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
#>         logitfDepot <- 1.38629436111989
#>         label("Bioavailability fraction (fDepot)")
#>     })
#>     model({
#>         fDepot <- expit(logitfDepot, 0, 1)
#>         ka <- exp(lka)
#>         cl <- exp(lcl)
#>         vc <- exp(lvc)
#>         kel <- cl/vc
#>         d/dt(depot) <- -ka * depot
#>         f(depot) <- fDepot
#>         d/dt(central) <- ka * depot - kel * central
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }

# dose split between two absorption paths, fraction named for depot
readModelDb("PK_1cmt_des") |>
  addDepot(depot = "depot2", ka = "ka2") |>
  addBioavailability(depot, depot2, scale = "logit")
#>  
#>  
#>  ── rxode2-based free-form 3-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>         lka         lcl         lvc      propSd        lka2 logitfDepot 
#>    0.450000    1.000000    3.450000    0.500000    0.100000    1.386294 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2           depot2
#> 3                  3          central
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
#>         lka2 <- 0.1
#>         label("First order absorption rate (ka2)")
#>         logitfDepot <- 1.38629436111989
#>         label("Fraction of dose absorbed from depot (fDepot)")
#>     })
#>     model({
#>         fDepot <- expit(logitfDepot, 0, 1)
#>         ka <- exp(lka)
#>         cl <- exp(lcl)
#>         vc <- exp(lvc)
#>         kel <- cl/vc
#>         d/dt(depot) <- -ka * depot
#>         f(depot) <- fDepot
#>         ka2 <- exp(lka2)
#>         d/dt(depot2) <- -ka2 * depot2
#>         f(depot2) <- 1 - fDepot
#>         d/dt(central) <- ka * depot - kel * central + ka2 * depot2
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }
```
