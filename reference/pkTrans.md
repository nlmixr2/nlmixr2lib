# Change the transformation type for PK models

Change the transformation type for PK models

## Usage

``` r
pkTrans(
  ui,
  type = c("k", "k21", "vss", "aob", "alpha"),
  k13 = "k13",
  k31 = "k31",
  k12 = "k12",
  k21 = "k21",
  kel = "kel",
  vc = "vc",
  cl = "cl",
  vp = "vp",
  q = "q",
  vp2 = "vp2",
  q2 = "q2",
  vss = "vss",
  aob = "aob",
  alpha = "alpha",
  beta = "beta",
  gam = "gam",
  A = "A",
  B = "B",
  C = "C",
  s = "s",
  p = "p",
  tmp = "tmp",
  beforeCmt = c("depot", "central")
)
```

## Arguments

- ui:

  A model in terms of Clearance

- type:

  the type of PK transformation to make:

  \- `"k"`: Change to rate constants (kel, k12, k21, k13, k31)

  \- `"vss"`: Change to volume of distribution at steady state (cl, vc,
  q, vss)

  \- `"aob"`: Change to A/B ratio (aob, alpha, beta, vc)

  \- `"k21"`: Change to k21 constant (k21, alpha, beta, vc) or (k21,
  k31, alpha, beta, gam, vc)

  \- `"alpha"`: Change to macro constants (alpha, beta, gam, A, B, C,
  vc)

- k13:

  name of rate constant from central to periph2

- k31:

  name of rate constant from periph2 to central

- k12:

  name of rate constant from central to periph1

- k21:

  name of rate constant from periph1 to central

- kel:

  name of elimination rate constant

- vc:

  name of central compartment volume

- cl:

  name of clearance

- vp:

  name of volume of periph1

- q:

  name of intercompartmental clearance between central and periph1

- vp2:

  name of volume of periph2

- q2:

  name of intercompartmental clearance between central and periph2

- vss:

  name of volume of distribution at steady state

- aob:

  A/B ratio

- alpha:

  macro constant name for first exponential decay term

- beta:

  macro constant name for second exponential decay term

- gam:

  macro constant name for third exponential decay term

- A:

  macro coefficient for the first exponential decay term (corresponds
  with alpha)

- B:

  macro coefficient for the second exponential decay term (corresponds
  with beta)

- C:

  macro coefficient for the third exponential decay term (corresponds
  with gam)

- s:

  sum constant name for the k12 three compartment

- p:

  product constant name for the k12 three compartment

- tmp:

  name of temporary variable for the three compartment with \`A\`,
  \`B\`, \`C\`, \`alpha\`, \`beta\` and \`gam\`.

- beforeCmt:

  if the model is compartmental you can specify the preferred names
  where the estimates and extra lines are added before

## Value

ui with no PK parameters estimated

## Author

Matthew L. Fidler

## Examples

``` r
# \donttest{
# Three compartment model translations

readModelDb("PK_3cmt_des") |>
  pkTrans("k")
#>  
#>  
#>  ── rxode2-based free-form 4-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lka propSd   lk12   lk21   lk13   lk31   lkel    lvc 
#>   0.45   0.50   0.10   0.10   0.10   0.10   0.10   0.10 
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
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         lk12 <- 0.1
#>         label("Central->Periph1 constant (k12)")
#>         lk21 <- 0.1
#>         label("Periph1->Central constant (k21)")
#>         lk13 <- 0.1
#>         label("Central->Periph2 constant (k13)")
#>         lk31 <- 0.1
#>         label("Periph2->Central constant (k31)")
#>         lkel <- 0.1
#>         label("Elimination from central (kel)")
#>         lvc <- 0.1
#>         label("Central compartment volume (vc)")
#>     })
#>     model({
#>         ka <- exp(lka)
#>         k12 <- exp(lk12)
#>         k21 <- exp(lk21)
#>         k13 <- exp(lk13)
#>         k31 <- exp(lk31)
#>         kel <- exp(lkel)
#>         vc <- exp(lvc)
#>         d/dt(depot) <- -ka * depot
#>         d/dt(central) <- ka * depot - kel * central - k12 * central + 
#>             k21 * peripheral1 - k13 * central + k31 * peripheral2
#>         d/dt(peripheral1) <- k12 * central - k21 * peripheral1
#>         d/dt(peripheral2) <- k13 * central - k31 * peripheral2
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }

readModelDb("PK_3cmt_des") |>
  pkTrans("k21")
#>  
#>  
#>  ── rxode2-based free-form 4-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lka propSd   lk21   lk31 lalpha  lbeta   lgam    lvc 
#>   0.45   0.50   0.10   0.10   0.10   0.10   0.10   0.10 
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
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         lk21 <- 0.1
#>         label("Periph1->Central constant (k21)")
#>         lk31 <- 0.1
#>         label("Periph2->Central constant (k31)")
#>         lalpha <- 0.1
#>         label("alpha macro constant (alpha)")
#>         lbeta <- 0.1
#>         label("beta macro constant (beta)")
#>         lgam <- 0.1
#>         label("gam macro constant (gam)")
#>         lvc <- 0.1
#>         label("Volume of central compartment (vc)")
#>     })
#>     model({
#>         ka <- exp(lka)
#>         k21 <- exp(lk21)
#>         k31 <- exp(lk31)
#>         alpha <- exp(lalpha)
#>         beta <- exp(lbeta)
#>         gam <- exp(lgam)
#>         vc <- exp(lvc)
#>         kel <- alpha * beta * gam/(k21 * k31)
#>         s <- alpha + beta + gam
#>         p <- alpha * beta + alpha * gam + beta * gam
#>         k13 <- (p + k31 * k31 - k31 * s - kel * k21)/(k21 - k31)
#>         k12 <- s - kel - k13 - k21 - k31
#>         d/dt(depot) <- -ka * depot
#>         d/dt(central) <- ka * depot - kel * central - k12 * central + 
#>             k21 * peripheral1 - k13 * central + k31 * peripheral2
#>         d/dt(peripheral1) <- k12 * central - k21 * peripheral1
#>         d/dt(peripheral2) <- k13 * central - k31 * peripheral2
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }

readModelDb("PK_3cmt_des") |>
  pkTrans("alpha")
#>  
#>  
#>  ── rxode2-based free-form 4-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lka propSd lalpha  lbeta   lgam     lA     lB     lC 
#>   0.45   0.50   0.10   0.10   0.10   0.10   0.10   0.10 
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
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         lalpha <- 0.1
#>         label("alpha macro constant (alpha)")
#>         lbeta <- 0.1
#>         label("beta macro constant (beta)")
#>         lgam <- 0.1
#>         label("gam macro constant (gam)")
#>         lA <- 0.1
#>         label("A coefficient (A)")
#>         lB <- 0.1
#>         label("B coefficient (B)")
#>         lC <- 0.1
#>         label("C coefficent (C)")
#>     })
#>     model({
#>         ka <- exp(lka)
#>         alpha <- exp(lalpha)
#>         beta <- exp(lbeta)
#>         gam <- exp(lgam)
#>         A <- exp(lA)
#>         B <- exp(lB)
#>         C <- exp(lC)
#>         vc <- 1/(A + B + C)
#>         s <- -(alpha * C + alpha * B + gam * A + gam * B + beta * 
#>             A + beta * C) * vc
#>         p <- (alpha * beta * C + alpha * gam * B + beta * gam * 
#>             A) * vc
#>         tmp <- sqrt(p * p - 4 * s)
#>         k21 <- 0.5 * (-p + tmp)
#>         k31 <- 0.5 * (-p - tmp)
#>         kel <- alpha * beta * gam/(k21 * k31)
#>         k12 <- ((beta * gam + alpha * beta + alpha * gam) - k21 * 
#>             (alpha + beta + gam) - kel * k31 + k21 * k21)/(k31 - 
#>             k21)
#>         k13 <- alpha + beta + gam - (kel + k12 + k21 + k31)
#>         d/dt(depot) <- -ka * depot
#>         d/dt(central) <- ka * depot - kel * central - k12 * central + 
#>             k21 * peripheral1 - k13 * central + k31 * peripheral2
#>         d/dt(peripheral1) <- k12 * central - k21 * peripheral1
#>         d/dt(peripheral2) <- k13 * central - k31 * peripheral2
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }

# The most types of transformations are
# available for 2 compartment models

readModelDb("PK_2cmt_des") |>
  pkTrans("k")
#>  
#>  
#>  ── rxode2-based free-form 3-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lka propSd   lk12   lk21   lkel    lvc 
#>   0.45   0.50   0.10   0.10   0.10   0.10 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2          central
#> 3                  3      peripheral1
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     ini({
#>         lka <- 0.45
#>         label("Absorption rate (Ka)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         lk12 <- 0.1
#>         label("Central->Periph1 constant (k12)")
#>         lk21 <- 0.1
#>         label("Periph1->Central constant (k21)")
#>         lkel <- 0.1
#>         label("Elimination from central (kel)")
#>         lvc <- 0.1
#>         label("Central compartment volume (vc)")
#>     })
#>     model({
#>         ka <- exp(lka)
#>         k12 <- exp(lk12)
#>         k21 <- exp(lk21)
#>         kel <- exp(lkel)
#>         vc <- exp(lvc)
#>         d/dt(depot) <- -ka * depot
#>         d/dt(central) <- ka * depot - kel * central - k12 * central + 
#>             k21 * peripheral1
#>         d/dt(peripheral1) <- k12 * central - k21 * peripheral1
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }

readModelDb("PK_2cmt_des") |>
  pkTrans("vss")
#>  
#>  
#>  ── rxode2-based free-form 3-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lka propSd    lcl    lvc     lq   lvss 
#>   0.45   0.50   0.10   0.10   0.10   0.10 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2          central
#> 3                  3      peripheral1
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     ini({
#>         lka <- 0.45
#>         label("Absorption rate (Ka)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         lcl <- 0.1
#>         label("Clearance (cl)")
#>         lvc <- 0.1
#>         label("Central compartment volume (vc)")
#>         lq <- 0.1
#>         label("Periph1<->Central inter-compartmental clearance (q)")
#>         lvss <- 0.1
#>         label("Volume of distribution at steady state (vss)")
#>     })
#>     model({
#>         ka <- exp(lka)
#>         cl <- exp(lcl)
#>         vc <- exp(lvc)
#>         q <- exp(lq)
#>         vss <- exp(lvss)
#>         kel <- cl/vc
#>         k12 <- q/vc
#>         k21 <- q/(vss - vc)
#>         d/dt(depot) <- -ka * depot
#>         d/dt(central) <- ka * depot - kel * central - k12 * central + 
#>             k21 * peripheral1
#>         d/dt(peripheral1) <- k12 * central - k21 * peripheral1
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }

readModelDb("PK_2cmt_des") |>
  pkTrans("aob")
#>  
#>  
#>  ── rxode2-based free-form 3-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lka propSd   laob lalpha  lbeta    lvc 
#>   0.45   0.50   0.10   0.10   0.10   0.10 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2          central
#> 3                  3      peripheral1
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     ini({
#>         lka <- 0.45
#>         label("Absorption rate (Ka)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         laob <- 0.1
#>         label("A/B (aob)")
#>         lalpha <- 0.1
#>         label("alpha macro constant (alpha)")
#>         lbeta <- 0.1
#>         label("beta macro constant (beta)")
#>         lvc <- 0.1
#>         label("Volume of central compartment (vc)")
#>     })
#>     model({
#>         ka <- exp(lka)
#>         aob <- exp(laob)
#>         alpha <- exp(lalpha)
#>         beta <- exp(lbeta)
#>         vc <- exp(lvc)
#>         k21 <- (aob * beta + alpha)/(aob + 1)
#>         kel <- (alpha * beta)/k21
#>         k12 <- alpha + beta - k21 - kel
#>         d/dt(depot) <- -ka * depot
#>         d/dt(central) <- ka * depot - kel * central - k12 * central + 
#>             k21 * peripheral1
#>         d/dt(peripheral1) <- k12 * central - k21 * peripheral1
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }

readModelDb("PK_2cmt_des") |>
  pkTrans("k21")
#>  
#>  
#>  ── rxode2-based free-form 3-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lka propSd   lk21 lalpha  lbeta    lvc 
#>   0.45   0.50   0.10   0.10   0.10   0.10 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2          central
#> 3                  3      peripheral1
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     ini({
#>         lka <- 0.45
#>         label("Absorption rate (Ka)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         lk21 <- 0.1
#>         label("Periph1->Central constant (k21)")
#>         lalpha <- 0.1
#>         label("alpha macro constant (alpha)")
#>         lbeta <- 0.1
#>         label("beta macro constant (beta)")
#>         lvc <- 0.1
#>         label("Volume of central compartment (vc)")
#>     })
#>     model({
#>         ka <- exp(lka)
#>         k21 <- exp(lk21)
#>         alpha <- exp(lalpha)
#>         beta <- exp(lbeta)
#>         vc <- exp(lvc)
#>         kel <- alpha * beta/k21
#>         k12 <- alpha + beta - k21 - kel
#>         d/dt(depot) <- -ka * depot
#>         d/dt(central) <- ka * depot - kel * central - k12 * central + 
#>             k21 * peripheral1
#>         d/dt(peripheral1) <- k12 * central - k21 * peripheral1
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }

readModelDb("PK_2cmt_des") |>
  pkTrans("alpha")
#>  
#>  
#>  ── rxode2-based free-form 3-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lka propSd lalpha  lbeta     lA     lB 
#>   0.45   0.50   0.10   0.10   0.10   0.10 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2          central
#> 3                  3      peripheral1
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     ini({
#>         lka <- 0.45
#>         label("Absorption rate (Ka)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         lalpha <- 0.1
#>         label("alpha macro constant (alpha)")
#>         lbeta <- 0.1
#>         label("beta macro constant (beta)")
#>         lA <- 0.1
#>         label("A coefficient (A)")
#>         lB <- 0.1
#>         label("B coefficient (B)")
#>     })
#>     model({
#>         ka <- exp(lka)
#>         alpha <- exp(lalpha)
#>         beta <- exp(lbeta)
#>         A <- exp(lA)
#>         B <- exp(lB)
#>         vc <- 1/(A + B)
#>         k21 <- (A * beta + B * alpha) * vc
#>         kel <- alpha * beta/k21
#>         k12 <- alpha + beta - k21 - kel
#>         d/dt(depot) <- -ka * depot
#>         d/dt(central) <- ka * depot - kel * central - k12 * central + 
#>             k21 * peripheral1
#>         d/dt(peripheral1) <- k12 * central - k21 * peripheral1
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }

# One compartment transformations are also available:

readModelDb("PK_1cmt_des") |>
  pkTrans("k")
#>  
#>  
#>  ── rxode2-based free-form 2-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lka propSd   lkel    lvc 
#>   0.45   0.50   0.10   0.10 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2          central
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     dosing <- c("central", "depot")
#>     ini({
#>         lka <- 0.45
#>         label("Absorption rate (Ka)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         lkel <- 0.1
#>         label("Elimination from central (kel)")
#>         lvc <- 0.1
#>         label("Central compartment volume (vc)")
#>     })
#>     model({
#>         ka <- exp(lka)
#>         kel <- exp(lkel)
#>         vc <- exp(lvc)
#>         d/dt(depot) <- -ka * depot
#>         d/dt(central) <- ka * depot - kel * central
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }

readModelDb("PK_1cmt_des") |>
  pkTrans("alpha")
#>  
#>  
#>  ── rxode2-based free-form 2-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lka propSd lalpha     lA 
#>   0.45   0.50   0.10   0.10 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2          central
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     dosing <- c("central", "depot")
#>     ini({
#>         lka <- 0.45
#>         label("Absorption rate (Ka)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         lalpha <- 0.1
#>         label("alpha macro constant (alpha)")
#>         lA <- 0.1
#>         label("A coefficient (A)")
#>     })
#>     model({
#>         ka <- exp(lka)
#>         alpha <- exp(lalpha)
#>         A <- exp(lA)
#>         kel <- alpha
#>         vc <- 1/A
#>         d/dt(depot) <- -ka * depot
#>         d/dt(central) <- ka * depot - kel * central
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }

# also works without depot:

readModelDb("PK_3cmt_des") |>
  removeDepot() |>
  pkTrans("k")
#>  
#>  
#>  ── rxode2-based free-form 3-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#> propSd   lk12   lk21   lk13   lk31   lkel    lvc 
#>    0.5    0.1    0.1    0.1    0.1    0.1    0.1 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1          central
#> 2                  2      peripheral1
#> 3                  3      peripheral2
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     ini({
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         lk12 <- 0.1
#>         label("Central->Periph1 constant (k12)")
#>         lk21 <- 0.1
#>         label("Periph1->Central constant (k21)")
#>         lk13 <- 0.1
#>         label("Central->Periph2 constant (k13)")
#>         lk31 <- 0.1
#>         label("Periph2->Central constant (k31)")
#>         lkel <- 0.1
#>         label("Elimination from central (kel)")
#>         lvc <- 0.1
#>         label("Central compartment volume (vc)")
#>     })
#>     model({
#>         k12 <- exp(lk12)
#>         k21 <- exp(lk21)
#>         k13 <- exp(lk13)
#>         k31 <- exp(lk31)
#>         kel <- exp(lkel)
#>         vc <- exp(lvc)
#>         d/dt(central) <- -kel * central - k12 * central + k21 * 
#>             peripheral1 - k13 * central + k31 * peripheral2
#>         d/dt(peripheral1) <- k12 * central - k21 * peripheral1
#>         d/dt(peripheral2) <- k13 * central - k31 * peripheral2
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }

# }
```
