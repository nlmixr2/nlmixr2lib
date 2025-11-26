# Add a property to a compartment

Add a property to a compartment

## Usage

``` r
addCmtProp(ui, prop = c("f", "lag", "dur", "rate", "ini"), cmt)

addBioavailability(ui, cmt)

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
```
