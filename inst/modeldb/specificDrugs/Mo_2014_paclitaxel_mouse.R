Mo_2014_paclitaxel_mouse <- function() {
  description <- "Preclinical (mouse, HCT116 human colon carcinoma xenograft). Lifespan-based tumor growth inhibition (LS TGI) model of intravenous paclitaxel with a cycle-specific drug effect (Mo 2014 Case Study 1, simultaneous fit of the control and 30 mg/kg Q4D x3 arms, Table 2). Tumor cells have a common lifespan td_div at the end of which they divide; the division efficiency p(w) falls from p0 towards 1 as the total tumor weight w passes the threshold wth (switch sharpness psi = 20). The tumor production rate is the recursive lifespan product kin(t) = kin0 * prod over past lifespans of (1 - E(C(t - i td_div))) * p(w(t - i td_div)), so the model is a delay differential equation with one delayed tumor weight and one delayed plasma concentration per elapsed lifespan; these are encoded with rxode2 delay() for 40 lifespans plus their copies shifted by the apoptosis duration td_apo, which makes the model valid for t < 40.5 * td_div (59 days at the published td_div). Paclitaxel acts only on cells reaching the end of their lifespan (sigmoid Emax on plasma concentration), routing them into a non-proliferating apoptotic pool that dies after td_apo. Following the authors' published MATLAB implementation, the number of past lifespans in each product is round(t / td_div), not the floor INT(t / td_div) printed in Equations 7, 30 and 31 (see the vignette). Two-compartment paclitaxel PK fixed to the Simeoni 2004 values quoted in the paper. Time is in days, doses in mg/kg into central, plasma concentration in ng/mL and tumor size in g. Solve with method = 'dop853' and maxsteps = 1e7."
  reference <- paste(
    "Mo G, Gibbons F, Schroeder P, Krzyzanski W. Lifespan based",
    "pharmacokinetic-pharmacodynamic model of tumor growth inhibition by",
    "anticancer therapeutics. PLoS One. 2014;9(10):e109747.",
    "doi:10.1371/journal.pone.0109747"
  )
  vignette <- "Mo_2014_lifespan_tgi"
  units <- list(time = "day", dosing = "mg/kg", concentration = "ng/mL")

  compartmentData <- list(
    central = list(analyte = "paclitaxel", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "paclitaxel", units = "mg/kg", specimen = "tissue", verified = TRUE),
    cycling_cells = list(analyte = "proliferating tumor cells M", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells1 = list(
      analyte = "apoptotic (non-proliferating) tumor cells A",
      units = "g",
      specimen = "tumor",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (HCT116 human colon carcinoma xenograft)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    disease_state = "subcutaneous HCT116 human colon carcinoma xenograft",
    dose_range = "paclitaxel 30 mg/kg intravenous every 4 days for 3 doses starting on day 8 after tumor inoculation, plus untreated controls",
    notes = "Mean tumor-weight data digitized by Mo 2014 from Simeoni et al. (Cancer Res 2004;64:1094-1101); the number of animals per arm is not restated. The model was fitted to the mean profiles with MATLAB nlinfit (Methods, Data Analysis); td_div and td_apo were found by grid search and then held fixed. Initial tumor weight 0.033 g at inoculation (Figure 7 caption)."
  )

  ini({
    # Lifespan TGI parameters: Mo 2014 Table 2 (simultaneous control +
    # paclitaxel fit, cycle-specific drug effect).
    ltd_div <- fixed(log(1.46)); label("Tumor cell lifespan (doubling time) T (day)") # Table 2: T = 1.46 day, fixed from the Table 1 grid search
    ltd_apo <- fixed(log(0.536)); label("Apoptosis duration TA (day)") # Table 2: TA = 0.536 day, fixed after grid search
    lp0 <- log(1.44); label("Cell division efficiency below the size threshold p0 (unitless)") # Table 2: p0 = 1.44
    lkin0 <- log(4.04e-2); label("Tumor production rate before the experiment kin0 (g/day)") # Table 2: kin0 = 4.04e-2 g/day
    lwth <- log(2.48); label("Tumor size threshold for falling division efficiency wth (g)") # Table 2: wth = 2.48 g
    psi <- fixed(20); label("Switch sharpness of the division-efficiency function psi (unitless)") # Model development, Eq. 13-14 and Figure 2B: psi = 20; supplement S5/S8 code psi = 20
    lemax <- fixed(log(1)); label("Maximum cycle-specific drug effect Emax (fraction)") # Table 2: Emax = 1, fixed for lack of escalating doses
    lec50 <- log(9.45); label("Paclitaxel plasma concentration at half-maximal effect EC50 (ng/mL)") # Table 2: EC50 = 9.45 ng/mL
    lrbase_tumor <- fixed(log(0.033)); label("Tumor weight at inoculation w0 (g)") # Figure 7 caption: initial tumor fixed at 0.033 g

    # Paclitaxel PK, fixed to the Simeoni 2004 two-compartment model
    # (Methods, Literature Data Acquisition). The paper reports the rate
    # constants in 1/h; they are converted to 1/day (x 24) for this
    # day-based model.
    lvc <- fixed(log(0.81)); label("Paclitaxel central volume V (L/kg)") # Methods: V = 0.81 L/kg
    lkel <- fixed(log(0.868 * 24)); label("Paclitaxel elimination rate constant kel (1/day)") # Methods: kel = 0.868 1/h
    lk12 <- fixed(log(0.006 * 24)); label("Paclitaxel central-to-peripheral rate constant k12 (1/day)") # Methods: k12 = 0.006 1/h
    lk21 <- fixed(log(0.0838 * 24)); label("Paclitaxel peripheral-to-central rate constant k21 (1/day)") # Methods: k21 = 0.0838 1/h

    # Mo 2014 fitted the mean profiles by unweighted nonlinear least
    # squares (MATLAB nlinfit) and reports no residual-error estimate.
    addSd_tumor_vol <- fixed(0); label("Additive residual error on tumor size (g; not reported)") # not reported (nlinfit least squares)
  })

  model({
    td_div <- exp(ltd_div)
    td_apo <- exp(ltd_apo)
    p0 <- exp(lp0)
    kin0 <- exp(lkin0)
    wth <- exp(lwth)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    rbase <- exp(lrbase_tumor)
    vc <- exp(lvc)
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)

    # Paclitaxel PK (dose in mg/kg into central; mg/L x 1000 = ng/mL)
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    Cc <- central / vc * 1000

    # Tumor weight w = M + A (Eq. 17); division efficiency p(w) (Eq. 13)
    tumor_vol <- cycling_cells + damaged_cells1
    pw <- 1 + (p0 - 1) / (1 + (tumor_vol / wth)^psi)^(1 / psi)
    # Cycle-specific Emax drug effect (Eq. 27)
    eff <- emax * Cc / (ec50 + Cc)

    # Delayed terms for i = 1..40 past lifespans: w_di = w(t - i T),
    # c_di = C(t - i T), p_di = p(w_di), s_di = 1 - E(c_di). Before the
    # experiment (t - i T < 0) the tumor is at w0, A = 0 and C = 0
    # (Eqs. 9, 25 and the supplement history functions S6 / S9). As in the
    # supplement S8 code, p_di enters once round(t / T) >= i, that is from
    # t = (i - 0.5) T, and is then p(w0) until t = i T.
    w_d1 <- rbase; c_d1 <- 0
    if (t >= td_div) {w_d1 <- delay(cycling_cells, td_div) + delay(damaged_cells1, td_div); c_d1 <- delay(central, td_div) / vc * 1000}
    p_d1 <- 1; if (t >= 0.5 * td_div) {p_d1 <- 1 + (p0 - 1) / (1 + (w_d1 / wth)^psi)^(1 / psi)}
    s_d1 <- 1 - emax * c_d1 / (ec50 + c_d1)
    w_d2 <- rbase; c_d2 <- 0
    if (t >= 2 * td_div) {w_d2 <- delay(cycling_cells, 2 * td_div) + delay(damaged_cells1, 2 * td_div); c_d2 <- delay(central, 2 * td_div) / vc * 1000}
    p_d2 <- 1; if (t >= 1.5 * td_div) {p_d2 <- 1 + (p0 - 1) / (1 + (w_d2 / wth)^psi)^(1 / psi)}
    s_d2 <- 1 - emax * c_d2 / (ec50 + c_d2)
    w_d3 <- rbase; c_d3 <- 0
    if (t >= 3 * td_div) {w_d3 <- delay(cycling_cells, 3 * td_div) + delay(damaged_cells1, 3 * td_div); c_d3 <- delay(central, 3 * td_div) / vc * 1000}
    p_d3 <- 1; if (t >= 2.5 * td_div) {p_d3 <- 1 + (p0 - 1) / (1 + (w_d3 / wth)^psi)^(1 / psi)}
    s_d3 <- 1 - emax * c_d3 / (ec50 + c_d3)
    w_d4 <- rbase; c_d4 <- 0
    if (t >= 4 * td_div) {w_d4 <- delay(cycling_cells, 4 * td_div) + delay(damaged_cells1, 4 * td_div); c_d4 <- delay(central, 4 * td_div) / vc * 1000}
    p_d4 <- 1; if (t >= 3.5 * td_div) {p_d4 <- 1 + (p0 - 1) / (1 + (w_d4 / wth)^psi)^(1 / psi)}
    s_d4 <- 1 - emax * c_d4 / (ec50 + c_d4)
    w_d5 <- rbase; c_d5 <- 0
    if (t >= 5 * td_div) {w_d5 <- delay(cycling_cells, 5 * td_div) + delay(damaged_cells1, 5 * td_div); c_d5 <- delay(central, 5 * td_div) / vc * 1000}
    p_d5 <- 1; if (t >= 4.5 * td_div) {p_d5 <- 1 + (p0 - 1) / (1 + (w_d5 / wth)^psi)^(1 / psi)}
    s_d5 <- 1 - emax * c_d5 / (ec50 + c_d5)
    w_d6 <- rbase; c_d6 <- 0
    if (t >= 6 * td_div) {w_d6 <- delay(cycling_cells, 6 * td_div) + delay(damaged_cells1, 6 * td_div); c_d6 <- delay(central, 6 * td_div) / vc * 1000}
    p_d6 <- 1; if (t >= 5.5 * td_div) {p_d6 <- 1 + (p0 - 1) / (1 + (w_d6 / wth)^psi)^(1 / psi)}
    s_d6 <- 1 - emax * c_d6 / (ec50 + c_d6)
    w_d7 <- rbase; c_d7 <- 0
    if (t >= 7 * td_div) {w_d7 <- delay(cycling_cells, 7 * td_div) + delay(damaged_cells1, 7 * td_div); c_d7 <- delay(central, 7 * td_div) / vc * 1000}
    p_d7 <- 1; if (t >= 6.5 * td_div) {p_d7 <- 1 + (p0 - 1) / (1 + (w_d7 / wth)^psi)^(1 / psi)}
    s_d7 <- 1 - emax * c_d7 / (ec50 + c_d7)
    w_d8 <- rbase; c_d8 <- 0
    if (t >= 8 * td_div) {w_d8 <- delay(cycling_cells, 8 * td_div) + delay(damaged_cells1, 8 * td_div); c_d8 <- delay(central, 8 * td_div) / vc * 1000}
    p_d8 <- 1; if (t >= 7.5 * td_div) {p_d8 <- 1 + (p0 - 1) / (1 + (w_d8 / wth)^psi)^(1 / psi)}
    s_d8 <- 1 - emax * c_d8 / (ec50 + c_d8)
    w_d9 <- rbase; c_d9 <- 0
    if (t >= 9 * td_div) {w_d9 <- delay(cycling_cells, 9 * td_div) + delay(damaged_cells1, 9 * td_div); c_d9 <- delay(central, 9 * td_div) / vc * 1000}
    p_d9 <- 1; if (t >= 8.5 * td_div) {p_d9 <- 1 + (p0 - 1) / (1 + (w_d9 / wth)^psi)^(1 / psi)}
    s_d9 <- 1 - emax * c_d9 / (ec50 + c_d9)
    w_d10 <- rbase; c_d10 <- 0
    if (t >= 10 * td_div) {w_d10 <- delay(cycling_cells, 10 * td_div) + delay(damaged_cells1, 10 * td_div); c_d10 <- delay(central, 10 * td_div) / vc * 1000}
    p_d10 <- 1; if (t >= 9.5 * td_div) {p_d10 <- 1 + (p0 - 1) / (1 + (w_d10 / wth)^psi)^(1 / psi)}
    s_d10 <- 1 - emax * c_d10 / (ec50 + c_d10)
    w_d11 <- rbase; c_d11 <- 0
    if (t >= 11 * td_div) {w_d11 <- delay(cycling_cells, 11 * td_div) + delay(damaged_cells1, 11 * td_div); c_d11 <- delay(central, 11 * td_div) / vc * 1000}
    p_d11 <- 1; if (t >= 10.5 * td_div) {p_d11 <- 1 + (p0 - 1) / (1 + (w_d11 / wth)^psi)^(1 / psi)}
    s_d11 <- 1 - emax * c_d11 / (ec50 + c_d11)
    w_d12 <- rbase; c_d12 <- 0
    if (t >= 12 * td_div) {w_d12 <- delay(cycling_cells, 12 * td_div) + delay(damaged_cells1, 12 * td_div); c_d12 <- delay(central, 12 * td_div) / vc * 1000}
    p_d12 <- 1; if (t >= 11.5 * td_div) {p_d12 <- 1 + (p0 - 1) / (1 + (w_d12 / wth)^psi)^(1 / psi)}
    s_d12 <- 1 - emax * c_d12 / (ec50 + c_d12)
    w_d13 <- rbase; c_d13 <- 0
    if (t >= 13 * td_div) {w_d13 <- delay(cycling_cells, 13 * td_div) + delay(damaged_cells1, 13 * td_div); c_d13 <- delay(central, 13 * td_div) / vc * 1000}
    p_d13 <- 1; if (t >= 12.5 * td_div) {p_d13 <- 1 + (p0 - 1) / (1 + (w_d13 / wth)^psi)^(1 / psi)}
    s_d13 <- 1 - emax * c_d13 / (ec50 + c_d13)
    w_d14 <- rbase; c_d14 <- 0
    if (t >= 14 * td_div) {w_d14 <- delay(cycling_cells, 14 * td_div) + delay(damaged_cells1, 14 * td_div); c_d14 <- delay(central, 14 * td_div) / vc * 1000}
    p_d14 <- 1; if (t >= 13.5 * td_div) {p_d14 <- 1 + (p0 - 1) / (1 + (w_d14 / wth)^psi)^(1 / psi)}
    s_d14 <- 1 - emax * c_d14 / (ec50 + c_d14)
    w_d15 <- rbase; c_d15 <- 0
    if (t >= 15 * td_div) {w_d15 <- delay(cycling_cells, 15 * td_div) + delay(damaged_cells1, 15 * td_div); c_d15 <- delay(central, 15 * td_div) / vc * 1000}
    p_d15 <- 1; if (t >= 14.5 * td_div) {p_d15 <- 1 + (p0 - 1) / (1 + (w_d15 / wth)^psi)^(1 / psi)}
    s_d15 <- 1 - emax * c_d15 / (ec50 + c_d15)
    w_d16 <- rbase; c_d16 <- 0
    if (t >= 16 * td_div) {w_d16 <- delay(cycling_cells, 16 * td_div) + delay(damaged_cells1, 16 * td_div); c_d16 <- delay(central, 16 * td_div) / vc * 1000}
    p_d16 <- 1; if (t >= 15.5 * td_div) {p_d16 <- 1 + (p0 - 1) / (1 + (w_d16 / wth)^psi)^(1 / psi)}
    s_d16 <- 1 - emax * c_d16 / (ec50 + c_d16)
    w_d17 <- rbase; c_d17 <- 0
    if (t >= 17 * td_div) {w_d17 <- delay(cycling_cells, 17 * td_div) + delay(damaged_cells1, 17 * td_div); c_d17 <- delay(central, 17 * td_div) / vc * 1000}
    p_d17 <- 1; if (t >= 16.5 * td_div) {p_d17 <- 1 + (p0 - 1) / (1 + (w_d17 / wth)^psi)^(1 / psi)}
    s_d17 <- 1 - emax * c_d17 / (ec50 + c_d17)
    w_d18 <- rbase; c_d18 <- 0
    if (t >= 18 * td_div) {w_d18 <- delay(cycling_cells, 18 * td_div) + delay(damaged_cells1, 18 * td_div); c_d18 <- delay(central, 18 * td_div) / vc * 1000}
    p_d18 <- 1; if (t >= 17.5 * td_div) {p_d18 <- 1 + (p0 - 1) / (1 + (w_d18 / wth)^psi)^(1 / psi)}
    s_d18 <- 1 - emax * c_d18 / (ec50 + c_d18)
    w_d19 <- rbase; c_d19 <- 0
    if (t >= 19 * td_div) {w_d19 <- delay(cycling_cells, 19 * td_div) + delay(damaged_cells1, 19 * td_div); c_d19 <- delay(central, 19 * td_div) / vc * 1000}
    p_d19 <- 1; if (t >= 18.5 * td_div) {p_d19 <- 1 + (p0 - 1) / (1 + (w_d19 / wth)^psi)^(1 / psi)}
    s_d19 <- 1 - emax * c_d19 / (ec50 + c_d19)
    w_d20 <- rbase; c_d20 <- 0
    if (t >= 20 * td_div) {w_d20 <- delay(cycling_cells, 20 * td_div) + delay(damaged_cells1, 20 * td_div); c_d20 <- delay(central, 20 * td_div) / vc * 1000}
    p_d20 <- 1; if (t >= 19.5 * td_div) {p_d20 <- 1 + (p0 - 1) / (1 + (w_d20 / wth)^psi)^(1 / psi)}
    s_d20 <- 1 - emax * c_d20 / (ec50 + c_d20)
    w_d21 <- rbase; c_d21 <- 0
    if (t >= 21 * td_div) {w_d21 <- delay(cycling_cells, 21 * td_div) + delay(damaged_cells1, 21 * td_div); c_d21 <- delay(central, 21 * td_div) / vc * 1000}
    p_d21 <- 1; if (t >= 20.5 * td_div) {p_d21 <- 1 + (p0 - 1) / (1 + (w_d21 / wth)^psi)^(1 / psi)}
    s_d21 <- 1 - emax * c_d21 / (ec50 + c_d21)
    w_d22 <- rbase; c_d22 <- 0
    if (t >= 22 * td_div) {w_d22 <- delay(cycling_cells, 22 * td_div) + delay(damaged_cells1, 22 * td_div); c_d22 <- delay(central, 22 * td_div) / vc * 1000}
    p_d22 <- 1; if (t >= 21.5 * td_div) {p_d22 <- 1 + (p0 - 1) / (1 + (w_d22 / wth)^psi)^(1 / psi)}
    s_d22 <- 1 - emax * c_d22 / (ec50 + c_d22)
    w_d23 <- rbase; c_d23 <- 0
    if (t >= 23 * td_div) {w_d23 <- delay(cycling_cells, 23 * td_div) + delay(damaged_cells1, 23 * td_div); c_d23 <- delay(central, 23 * td_div) / vc * 1000}
    p_d23 <- 1; if (t >= 22.5 * td_div) {p_d23 <- 1 + (p0 - 1) / (1 + (w_d23 / wth)^psi)^(1 / psi)}
    s_d23 <- 1 - emax * c_d23 / (ec50 + c_d23)
    w_d24 <- rbase; c_d24 <- 0
    if (t >= 24 * td_div) {w_d24 <- delay(cycling_cells, 24 * td_div) + delay(damaged_cells1, 24 * td_div); c_d24 <- delay(central, 24 * td_div) / vc * 1000}
    p_d24 <- 1; if (t >= 23.5 * td_div) {p_d24 <- 1 + (p0 - 1) / (1 + (w_d24 / wth)^psi)^(1 / psi)}
    s_d24 <- 1 - emax * c_d24 / (ec50 + c_d24)
    w_d25 <- rbase; c_d25 <- 0
    if (t >= 25 * td_div) {w_d25 <- delay(cycling_cells, 25 * td_div) + delay(damaged_cells1, 25 * td_div); c_d25 <- delay(central, 25 * td_div) / vc * 1000}
    p_d25 <- 1; if (t >= 24.5 * td_div) {p_d25 <- 1 + (p0 - 1) / (1 + (w_d25 / wth)^psi)^(1 / psi)}
    s_d25 <- 1 - emax * c_d25 / (ec50 + c_d25)
    w_d26 <- rbase; c_d26 <- 0
    if (t >= 26 * td_div) {w_d26 <- delay(cycling_cells, 26 * td_div) + delay(damaged_cells1, 26 * td_div); c_d26 <- delay(central, 26 * td_div) / vc * 1000}
    p_d26 <- 1; if (t >= 25.5 * td_div) {p_d26 <- 1 + (p0 - 1) / (1 + (w_d26 / wth)^psi)^(1 / psi)}
    s_d26 <- 1 - emax * c_d26 / (ec50 + c_d26)
    w_d27 <- rbase; c_d27 <- 0
    if (t >= 27 * td_div) {w_d27 <- delay(cycling_cells, 27 * td_div) + delay(damaged_cells1, 27 * td_div); c_d27 <- delay(central, 27 * td_div) / vc * 1000}
    p_d27 <- 1; if (t >= 26.5 * td_div) {p_d27 <- 1 + (p0 - 1) / (1 + (w_d27 / wth)^psi)^(1 / psi)}
    s_d27 <- 1 - emax * c_d27 / (ec50 + c_d27)
    w_d28 <- rbase; c_d28 <- 0
    if (t >= 28 * td_div) {w_d28 <- delay(cycling_cells, 28 * td_div) + delay(damaged_cells1, 28 * td_div); c_d28 <- delay(central, 28 * td_div) / vc * 1000}
    p_d28 <- 1; if (t >= 27.5 * td_div) {p_d28 <- 1 + (p0 - 1) / (1 + (w_d28 / wth)^psi)^(1 / psi)}
    s_d28 <- 1 - emax * c_d28 / (ec50 + c_d28)
    w_d29 <- rbase; c_d29 <- 0
    if (t >= 29 * td_div) {w_d29 <- delay(cycling_cells, 29 * td_div) + delay(damaged_cells1, 29 * td_div); c_d29 <- delay(central, 29 * td_div) / vc * 1000}
    p_d29 <- 1; if (t >= 28.5 * td_div) {p_d29 <- 1 + (p0 - 1) / (1 + (w_d29 / wth)^psi)^(1 / psi)}
    s_d29 <- 1 - emax * c_d29 / (ec50 + c_d29)
    w_d30 <- rbase; c_d30 <- 0
    if (t >= 30 * td_div) {w_d30 <- delay(cycling_cells, 30 * td_div) + delay(damaged_cells1, 30 * td_div); c_d30 <- delay(central, 30 * td_div) / vc * 1000}
    p_d30 <- 1; if (t >= 29.5 * td_div) {p_d30 <- 1 + (p0 - 1) / (1 + (w_d30 / wth)^psi)^(1 / psi)}
    s_d30 <- 1 - emax * c_d30 / (ec50 + c_d30)
    w_d31 <- rbase; c_d31 <- 0
    if (t >= 31 * td_div) {w_d31 <- delay(cycling_cells, 31 * td_div) + delay(damaged_cells1, 31 * td_div); c_d31 <- delay(central, 31 * td_div) / vc * 1000}
    p_d31 <- 1; if (t >= 30.5 * td_div) {p_d31 <- 1 + (p0 - 1) / (1 + (w_d31 / wth)^psi)^(1 / psi)}
    s_d31 <- 1 - emax * c_d31 / (ec50 + c_d31)
    w_d32 <- rbase; c_d32 <- 0
    if (t >= 32 * td_div) {w_d32 <- delay(cycling_cells, 32 * td_div) + delay(damaged_cells1, 32 * td_div); c_d32 <- delay(central, 32 * td_div) / vc * 1000}
    p_d32 <- 1; if (t >= 31.5 * td_div) {p_d32 <- 1 + (p0 - 1) / (1 + (w_d32 / wth)^psi)^(1 / psi)}
    s_d32 <- 1 - emax * c_d32 / (ec50 + c_d32)
    w_d33 <- rbase; c_d33 <- 0
    if (t >= 33 * td_div) {w_d33 <- delay(cycling_cells, 33 * td_div) + delay(damaged_cells1, 33 * td_div); c_d33 <- delay(central, 33 * td_div) / vc * 1000}
    p_d33 <- 1; if (t >= 32.5 * td_div) {p_d33 <- 1 + (p0 - 1) / (1 + (w_d33 / wth)^psi)^(1 / psi)}
    s_d33 <- 1 - emax * c_d33 / (ec50 + c_d33)
    w_d34 <- rbase; c_d34 <- 0
    if (t >= 34 * td_div) {w_d34 <- delay(cycling_cells, 34 * td_div) + delay(damaged_cells1, 34 * td_div); c_d34 <- delay(central, 34 * td_div) / vc * 1000}
    p_d34 <- 1; if (t >= 33.5 * td_div) {p_d34 <- 1 + (p0 - 1) / (1 + (w_d34 / wth)^psi)^(1 / psi)}
    s_d34 <- 1 - emax * c_d34 / (ec50 + c_d34)
    w_d35 <- rbase; c_d35 <- 0
    if (t >= 35 * td_div) {w_d35 <- delay(cycling_cells, 35 * td_div) + delay(damaged_cells1, 35 * td_div); c_d35 <- delay(central, 35 * td_div) / vc * 1000}
    p_d35 <- 1; if (t >= 34.5 * td_div) {p_d35 <- 1 + (p0 - 1) / (1 + (w_d35 / wth)^psi)^(1 / psi)}
    s_d35 <- 1 - emax * c_d35 / (ec50 + c_d35)
    w_d36 <- rbase; c_d36 <- 0
    if (t >= 36 * td_div) {w_d36 <- delay(cycling_cells, 36 * td_div) + delay(damaged_cells1, 36 * td_div); c_d36 <- delay(central, 36 * td_div) / vc * 1000}
    p_d36 <- 1; if (t >= 35.5 * td_div) {p_d36 <- 1 + (p0 - 1) / (1 + (w_d36 / wth)^psi)^(1 / psi)}
    s_d36 <- 1 - emax * c_d36 / (ec50 + c_d36)
    w_d37 <- rbase; c_d37 <- 0
    if (t >= 37 * td_div) {w_d37 <- delay(cycling_cells, 37 * td_div) + delay(damaged_cells1, 37 * td_div); c_d37 <- delay(central, 37 * td_div) / vc * 1000}
    p_d37 <- 1; if (t >= 36.5 * td_div) {p_d37 <- 1 + (p0 - 1) / (1 + (w_d37 / wth)^psi)^(1 / psi)}
    s_d37 <- 1 - emax * c_d37 / (ec50 + c_d37)
    w_d38 <- rbase; c_d38 <- 0
    if (t >= 38 * td_div) {w_d38 <- delay(cycling_cells, 38 * td_div) + delay(damaged_cells1, 38 * td_div); c_d38 <- delay(central, 38 * td_div) / vc * 1000}
    p_d38 <- 1; if (t >= 37.5 * td_div) {p_d38 <- 1 + (p0 - 1) / (1 + (w_d38 / wth)^psi)^(1 / psi)}
    s_d38 <- 1 - emax * c_d38 / (ec50 + c_d38)
    w_d39 <- rbase; c_d39 <- 0
    if (t >= 39 * td_div) {w_d39 <- delay(cycling_cells, 39 * td_div) + delay(damaged_cells1, 39 * td_div); c_d39 <- delay(central, 39 * td_div) / vc * 1000}
    p_d39 <- 1; if (t >= 38.5 * td_div) {p_d39 <- 1 + (p0 - 1) / (1 + (w_d39 / wth)^psi)^(1 / psi)}
    s_d39 <- 1 - emax * c_d39 / (ec50 + c_d39)
    w_d40 <- rbase; c_d40 <- 0
    if (t >= 40 * td_div) {w_d40 <- delay(cycling_cells, 40 * td_div) + delay(damaged_cells1, 40 * td_div); c_d40 <- delay(central, 40 * td_div) / vc * 1000}
    p_d40 <- 1; if (t >= 39.5 * td_div) {p_d40 <- 1 + (p0 - 1) / (1 + (w_d40 / wth)^psi)^(1 / psi)}
    s_d40 <- 1 - emax * c_d40 / (ec50 + c_d40)
    # kin(t) (Eq. 30 with round() for INT()) and kin(t - T) (Eq. 29):
    # kin(t - T) * (1 - E(C(t - T))) = kin(t) / p(w(t)) once t >= T, and
    # kin0 before (Eq. 5).
    kin <- kin0 * pw *
      p_d1 * p_d2 * p_d3 * p_d4 * p_d5 * p_d6 * p_d7 * p_d8 * p_d9 * p_d10 *
      p_d11 * p_d12 * p_d13 * p_d14 * p_d15 * p_d16 * p_d17 * p_d18 * p_d19 * p_d20 *
      p_d21 * p_d22 * p_d23 * p_d24 * p_d25 * p_d26 * p_d27 * p_d28 * p_d29 * p_d30 *
      p_d31 * p_d32 * p_d33 * p_d34 * p_d35 * p_d36 * p_d37 * p_d38 * p_d39 * p_d40 *
      s_d1 * s_d2 * s_d3 * s_d4 * s_d5 * s_d6 * s_d7 * s_d8 * s_d9 * s_d10 *
      s_d11 * s_d12 * s_d13 * s_d14 * s_d15 * s_d16 * s_d17 * s_d18 * s_d19 * s_d20 *
      s_d21 * s_d22 * s_d23 * s_d24 * s_d25 * s_d26 * s_d27 * s_d28 * s_d29 * s_d30 *
      s_d31 * s_d32 * s_d33 * s_d34 * s_d35 * s_d36 * s_d37 * s_d38 * s_d39 * s_d40
    kin_lag <- kin0
    if (t >= td_div) {kin_lag <- kin / pw}

    # kin(t - TA) for the apoptotic outflow: the same product evaluated
    # TA earlier (w_aj = w(t - TA - j T), c_aj = C(t - TA - j T)).
    w_a0 <- rbase; c_a0 <- 0
    if (t >= td_apo) {w_a0 <- delay(cycling_cells, td_apo) + delay(damaged_cells1, td_apo); c_a0 <- delay(central, td_apo) / vc * 1000}
    p_a0 <- 1 + (p0 - 1) / (1 + (w_a0 / wth)^psi)^(1 / psi)
    w_a1 <- rbase; c_a1 <- 0
    if (t >= td_apo + td_div) {w_a1 <- delay(cycling_cells, td_apo + td_div) + delay(damaged_cells1, td_apo + td_div); c_a1 <- delay(central, td_apo + td_div) / vc * 1000}
    p_a1 <- 1; if (t - td_apo >= 0.5 * td_div) {p_a1 <- 1 + (p0 - 1) / (1 + (w_a1 / wth)^psi)^(1 / psi)}
    s_a1 <- 1 - emax * c_a1 / (ec50 + c_a1)
    w_a2 <- rbase; c_a2 <- 0
    if (t >= td_apo + 2 * td_div) {w_a2 <- delay(cycling_cells, td_apo + 2 * td_div) + delay(damaged_cells1, td_apo + 2 * td_div); c_a2 <- delay(central, td_apo + 2 * td_div) / vc * 1000}
    p_a2 <- 1; if (t - td_apo >= 1.5 * td_div) {p_a2 <- 1 + (p0 - 1) / (1 + (w_a2 / wth)^psi)^(1 / psi)}
    s_a2 <- 1 - emax * c_a2 / (ec50 + c_a2)
    w_a3 <- rbase; c_a3 <- 0
    if (t >= td_apo + 3 * td_div) {w_a3 <- delay(cycling_cells, td_apo + 3 * td_div) + delay(damaged_cells1, td_apo + 3 * td_div); c_a3 <- delay(central, td_apo + 3 * td_div) / vc * 1000}
    p_a3 <- 1; if (t - td_apo >= 2.5 * td_div) {p_a3 <- 1 + (p0 - 1) / (1 + (w_a3 / wth)^psi)^(1 / psi)}
    s_a3 <- 1 - emax * c_a3 / (ec50 + c_a3)
    w_a4 <- rbase; c_a4 <- 0
    if (t >= td_apo + 4 * td_div) {w_a4 <- delay(cycling_cells, td_apo + 4 * td_div) + delay(damaged_cells1, td_apo + 4 * td_div); c_a4 <- delay(central, td_apo + 4 * td_div) / vc * 1000}
    p_a4 <- 1; if (t - td_apo >= 3.5 * td_div) {p_a4 <- 1 + (p0 - 1) / (1 + (w_a4 / wth)^psi)^(1 / psi)}
    s_a4 <- 1 - emax * c_a4 / (ec50 + c_a4)
    w_a5 <- rbase; c_a5 <- 0
    if (t >= td_apo + 5 * td_div) {w_a5 <- delay(cycling_cells, td_apo + 5 * td_div) + delay(damaged_cells1, td_apo + 5 * td_div); c_a5 <- delay(central, td_apo + 5 * td_div) / vc * 1000}
    p_a5 <- 1; if (t - td_apo >= 4.5 * td_div) {p_a5 <- 1 + (p0 - 1) / (1 + (w_a5 / wth)^psi)^(1 / psi)}
    s_a5 <- 1 - emax * c_a5 / (ec50 + c_a5)
    w_a6 <- rbase; c_a6 <- 0
    if (t >= td_apo + 6 * td_div) {w_a6 <- delay(cycling_cells, td_apo + 6 * td_div) + delay(damaged_cells1, td_apo + 6 * td_div); c_a6 <- delay(central, td_apo + 6 * td_div) / vc * 1000}
    p_a6 <- 1; if (t - td_apo >= 5.5 * td_div) {p_a6 <- 1 + (p0 - 1) / (1 + (w_a6 / wth)^psi)^(1 / psi)}
    s_a6 <- 1 - emax * c_a6 / (ec50 + c_a6)
    w_a7 <- rbase; c_a7 <- 0
    if (t >= td_apo + 7 * td_div) {w_a7 <- delay(cycling_cells, td_apo + 7 * td_div) + delay(damaged_cells1, td_apo + 7 * td_div); c_a7 <- delay(central, td_apo + 7 * td_div) / vc * 1000}
    p_a7 <- 1; if (t - td_apo >= 6.5 * td_div) {p_a7 <- 1 + (p0 - 1) / (1 + (w_a7 / wth)^psi)^(1 / psi)}
    s_a7 <- 1 - emax * c_a7 / (ec50 + c_a7)
    w_a8 <- rbase; c_a8 <- 0
    if (t >= td_apo + 8 * td_div) {w_a8 <- delay(cycling_cells, td_apo + 8 * td_div) + delay(damaged_cells1, td_apo + 8 * td_div); c_a8 <- delay(central, td_apo + 8 * td_div) / vc * 1000}
    p_a8 <- 1; if (t - td_apo >= 7.5 * td_div) {p_a8 <- 1 + (p0 - 1) / (1 + (w_a8 / wth)^psi)^(1 / psi)}
    s_a8 <- 1 - emax * c_a8 / (ec50 + c_a8)
    w_a9 <- rbase; c_a9 <- 0
    if (t >= td_apo + 9 * td_div) {w_a9 <- delay(cycling_cells, td_apo + 9 * td_div) + delay(damaged_cells1, td_apo + 9 * td_div); c_a9 <- delay(central, td_apo + 9 * td_div) / vc * 1000}
    p_a9 <- 1; if (t - td_apo >= 8.5 * td_div) {p_a9 <- 1 + (p0 - 1) / (1 + (w_a9 / wth)^psi)^(1 / psi)}
    s_a9 <- 1 - emax * c_a9 / (ec50 + c_a9)
    w_a10 <- rbase; c_a10 <- 0
    if (t >= td_apo + 10 * td_div) {w_a10 <- delay(cycling_cells, td_apo + 10 * td_div) + delay(damaged_cells1, td_apo + 10 * td_div); c_a10 <- delay(central, td_apo + 10 * td_div) / vc * 1000}
    p_a10 <- 1; if (t - td_apo >= 9.5 * td_div) {p_a10 <- 1 + (p0 - 1) / (1 + (w_a10 / wth)^psi)^(1 / psi)}
    s_a10 <- 1 - emax * c_a10 / (ec50 + c_a10)
    w_a11 <- rbase; c_a11 <- 0
    if (t >= td_apo + 11 * td_div) {w_a11 <- delay(cycling_cells, td_apo + 11 * td_div) + delay(damaged_cells1, td_apo + 11 * td_div); c_a11 <- delay(central, td_apo + 11 * td_div) / vc * 1000}
    p_a11 <- 1; if (t - td_apo >= 10.5 * td_div) {p_a11 <- 1 + (p0 - 1) / (1 + (w_a11 / wth)^psi)^(1 / psi)}
    s_a11 <- 1 - emax * c_a11 / (ec50 + c_a11)
    w_a12 <- rbase; c_a12 <- 0
    if (t >= td_apo + 12 * td_div) {w_a12 <- delay(cycling_cells, td_apo + 12 * td_div) + delay(damaged_cells1, td_apo + 12 * td_div); c_a12 <- delay(central, td_apo + 12 * td_div) / vc * 1000}
    p_a12 <- 1; if (t - td_apo >= 11.5 * td_div) {p_a12 <- 1 + (p0 - 1) / (1 + (w_a12 / wth)^psi)^(1 / psi)}
    s_a12 <- 1 - emax * c_a12 / (ec50 + c_a12)
    w_a13 <- rbase; c_a13 <- 0
    if (t >= td_apo + 13 * td_div) {w_a13 <- delay(cycling_cells, td_apo + 13 * td_div) + delay(damaged_cells1, td_apo + 13 * td_div); c_a13 <- delay(central, td_apo + 13 * td_div) / vc * 1000}
    p_a13 <- 1; if (t - td_apo >= 12.5 * td_div) {p_a13 <- 1 + (p0 - 1) / (1 + (w_a13 / wth)^psi)^(1 / psi)}
    s_a13 <- 1 - emax * c_a13 / (ec50 + c_a13)
    w_a14 <- rbase; c_a14 <- 0
    if (t >= td_apo + 14 * td_div) {w_a14 <- delay(cycling_cells, td_apo + 14 * td_div) + delay(damaged_cells1, td_apo + 14 * td_div); c_a14 <- delay(central, td_apo + 14 * td_div) / vc * 1000}
    p_a14 <- 1; if (t - td_apo >= 13.5 * td_div) {p_a14 <- 1 + (p0 - 1) / (1 + (w_a14 / wth)^psi)^(1 / psi)}
    s_a14 <- 1 - emax * c_a14 / (ec50 + c_a14)
    w_a15 <- rbase; c_a15 <- 0
    if (t >= td_apo + 15 * td_div) {w_a15 <- delay(cycling_cells, td_apo + 15 * td_div) + delay(damaged_cells1, td_apo + 15 * td_div); c_a15 <- delay(central, td_apo + 15 * td_div) / vc * 1000}
    p_a15 <- 1; if (t - td_apo >= 14.5 * td_div) {p_a15 <- 1 + (p0 - 1) / (1 + (w_a15 / wth)^psi)^(1 / psi)}
    s_a15 <- 1 - emax * c_a15 / (ec50 + c_a15)
    w_a16 <- rbase; c_a16 <- 0
    if (t >= td_apo + 16 * td_div) {w_a16 <- delay(cycling_cells, td_apo + 16 * td_div) + delay(damaged_cells1, td_apo + 16 * td_div); c_a16 <- delay(central, td_apo + 16 * td_div) / vc * 1000}
    p_a16 <- 1; if (t - td_apo >= 15.5 * td_div) {p_a16 <- 1 + (p0 - 1) / (1 + (w_a16 / wth)^psi)^(1 / psi)}
    s_a16 <- 1 - emax * c_a16 / (ec50 + c_a16)
    w_a17 <- rbase; c_a17 <- 0
    if (t >= td_apo + 17 * td_div) {w_a17 <- delay(cycling_cells, td_apo + 17 * td_div) + delay(damaged_cells1, td_apo + 17 * td_div); c_a17 <- delay(central, td_apo + 17 * td_div) / vc * 1000}
    p_a17 <- 1; if (t - td_apo >= 16.5 * td_div) {p_a17 <- 1 + (p0 - 1) / (1 + (w_a17 / wth)^psi)^(1 / psi)}
    s_a17 <- 1 - emax * c_a17 / (ec50 + c_a17)
    w_a18 <- rbase; c_a18 <- 0
    if (t >= td_apo + 18 * td_div) {w_a18 <- delay(cycling_cells, td_apo + 18 * td_div) + delay(damaged_cells1, td_apo + 18 * td_div); c_a18 <- delay(central, td_apo + 18 * td_div) / vc * 1000}
    p_a18 <- 1; if (t - td_apo >= 17.5 * td_div) {p_a18 <- 1 + (p0 - 1) / (1 + (w_a18 / wth)^psi)^(1 / psi)}
    s_a18 <- 1 - emax * c_a18 / (ec50 + c_a18)
    w_a19 <- rbase; c_a19 <- 0
    if (t >= td_apo + 19 * td_div) {w_a19 <- delay(cycling_cells, td_apo + 19 * td_div) + delay(damaged_cells1, td_apo + 19 * td_div); c_a19 <- delay(central, td_apo + 19 * td_div) / vc * 1000}
    p_a19 <- 1; if (t - td_apo >= 18.5 * td_div) {p_a19 <- 1 + (p0 - 1) / (1 + (w_a19 / wth)^psi)^(1 / psi)}
    s_a19 <- 1 - emax * c_a19 / (ec50 + c_a19)
    w_a20 <- rbase; c_a20 <- 0
    if (t >= td_apo + 20 * td_div) {w_a20 <- delay(cycling_cells, td_apo + 20 * td_div) + delay(damaged_cells1, td_apo + 20 * td_div); c_a20 <- delay(central, td_apo + 20 * td_div) / vc * 1000}
    p_a20 <- 1; if (t - td_apo >= 19.5 * td_div) {p_a20 <- 1 + (p0 - 1) / (1 + (w_a20 / wth)^psi)^(1 / psi)}
    s_a20 <- 1 - emax * c_a20 / (ec50 + c_a20)
    w_a21 <- rbase; c_a21 <- 0
    if (t >= td_apo + 21 * td_div) {w_a21 <- delay(cycling_cells, td_apo + 21 * td_div) + delay(damaged_cells1, td_apo + 21 * td_div); c_a21 <- delay(central, td_apo + 21 * td_div) / vc * 1000}
    p_a21 <- 1; if (t - td_apo >= 20.5 * td_div) {p_a21 <- 1 + (p0 - 1) / (1 + (w_a21 / wth)^psi)^(1 / psi)}
    s_a21 <- 1 - emax * c_a21 / (ec50 + c_a21)
    w_a22 <- rbase; c_a22 <- 0
    if (t >= td_apo + 22 * td_div) {w_a22 <- delay(cycling_cells, td_apo + 22 * td_div) + delay(damaged_cells1, td_apo + 22 * td_div); c_a22 <- delay(central, td_apo + 22 * td_div) / vc * 1000}
    p_a22 <- 1; if (t - td_apo >= 21.5 * td_div) {p_a22 <- 1 + (p0 - 1) / (1 + (w_a22 / wth)^psi)^(1 / psi)}
    s_a22 <- 1 - emax * c_a22 / (ec50 + c_a22)
    w_a23 <- rbase; c_a23 <- 0
    if (t >= td_apo + 23 * td_div) {w_a23 <- delay(cycling_cells, td_apo + 23 * td_div) + delay(damaged_cells1, td_apo + 23 * td_div); c_a23 <- delay(central, td_apo + 23 * td_div) / vc * 1000}
    p_a23 <- 1; if (t - td_apo >= 22.5 * td_div) {p_a23 <- 1 + (p0 - 1) / (1 + (w_a23 / wth)^psi)^(1 / psi)}
    s_a23 <- 1 - emax * c_a23 / (ec50 + c_a23)
    w_a24 <- rbase; c_a24 <- 0
    if (t >= td_apo + 24 * td_div) {w_a24 <- delay(cycling_cells, td_apo + 24 * td_div) + delay(damaged_cells1, td_apo + 24 * td_div); c_a24 <- delay(central, td_apo + 24 * td_div) / vc * 1000}
    p_a24 <- 1; if (t - td_apo >= 23.5 * td_div) {p_a24 <- 1 + (p0 - 1) / (1 + (w_a24 / wth)^psi)^(1 / psi)}
    s_a24 <- 1 - emax * c_a24 / (ec50 + c_a24)
    w_a25 <- rbase; c_a25 <- 0
    if (t >= td_apo + 25 * td_div) {w_a25 <- delay(cycling_cells, td_apo + 25 * td_div) + delay(damaged_cells1, td_apo + 25 * td_div); c_a25 <- delay(central, td_apo + 25 * td_div) / vc * 1000}
    p_a25 <- 1; if (t - td_apo >= 24.5 * td_div) {p_a25 <- 1 + (p0 - 1) / (1 + (w_a25 / wth)^psi)^(1 / psi)}
    s_a25 <- 1 - emax * c_a25 / (ec50 + c_a25)
    w_a26 <- rbase; c_a26 <- 0
    if (t >= td_apo + 26 * td_div) {w_a26 <- delay(cycling_cells, td_apo + 26 * td_div) + delay(damaged_cells1, td_apo + 26 * td_div); c_a26 <- delay(central, td_apo + 26 * td_div) / vc * 1000}
    p_a26 <- 1; if (t - td_apo >= 25.5 * td_div) {p_a26 <- 1 + (p0 - 1) / (1 + (w_a26 / wth)^psi)^(1 / psi)}
    s_a26 <- 1 - emax * c_a26 / (ec50 + c_a26)
    w_a27 <- rbase; c_a27 <- 0
    if (t >= td_apo + 27 * td_div) {w_a27 <- delay(cycling_cells, td_apo + 27 * td_div) + delay(damaged_cells1, td_apo + 27 * td_div); c_a27 <- delay(central, td_apo + 27 * td_div) / vc * 1000}
    p_a27 <- 1; if (t - td_apo >= 26.5 * td_div) {p_a27 <- 1 + (p0 - 1) / (1 + (w_a27 / wth)^psi)^(1 / psi)}
    s_a27 <- 1 - emax * c_a27 / (ec50 + c_a27)
    w_a28 <- rbase; c_a28 <- 0
    if (t >= td_apo + 28 * td_div) {w_a28 <- delay(cycling_cells, td_apo + 28 * td_div) + delay(damaged_cells1, td_apo + 28 * td_div); c_a28 <- delay(central, td_apo + 28 * td_div) / vc * 1000}
    p_a28 <- 1; if (t - td_apo >= 27.5 * td_div) {p_a28 <- 1 + (p0 - 1) / (1 + (w_a28 / wth)^psi)^(1 / psi)}
    s_a28 <- 1 - emax * c_a28 / (ec50 + c_a28)
    w_a29 <- rbase; c_a29 <- 0
    if (t >= td_apo + 29 * td_div) {w_a29 <- delay(cycling_cells, td_apo + 29 * td_div) + delay(damaged_cells1, td_apo + 29 * td_div); c_a29 <- delay(central, td_apo + 29 * td_div) / vc * 1000}
    p_a29 <- 1; if (t - td_apo >= 28.5 * td_div) {p_a29 <- 1 + (p0 - 1) / (1 + (w_a29 / wth)^psi)^(1 / psi)}
    s_a29 <- 1 - emax * c_a29 / (ec50 + c_a29)
    w_a30 <- rbase; c_a30 <- 0
    if (t >= td_apo + 30 * td_div) {w_a30 <- delay(cycling_cells, td_apo + 30 * td_div) + delay(damaged_cells1, td_apo + 30 * td_div); c_a30 <- delay(central, td_apo + 30 * td_div) / vc * 1000}
    p_a30 <- 1; if (t - td_apo >= 29.5 * td_div) {p_a30 <- 1 + (p0 - 1) / (1 + (w_a30 / wth)^psi)^(1 / psi)}
    s_a30 <- 1 - emax * c_a30 / (ec50 + c_a30)
    w_a31 <- rbase; c_a31 <- 0
    if (t >= td_apo + 31 * td_div) {w_a31 <- delay(cycling_cells, td_apo + 31 * td_div) + delay(damaged_cells1, td_apo + 31 * td_div); c_a31 <- delay(central, td_apo + 31 * td_div) / vc * 1000}
    p_a31 <- 1; if (t - td_apo >= 30.5 * td_div) {p_a31 <- 1 + (p0 - 1) / (1 + (w_a31 / wth)^psi)^(1 / psi)}
    s_a31 <- 1 - emax * c_a31 / (ec50 + c_a31)
    w_a32 <- rbase; c_a32 <- 0
    if (t >= td_apo + 32 * td_div) {w_a32 <- delay(cycling_cells, td_apo + 32 * td_div) + delay(damaged_cells1, td_apo + 32 * td_div); c_a32 <- delay(central, td_apo + 32 * td_div) / vc * 1000}
    p_a32 <- 1; if (t - td_apo >= 31.5 * td_div) {p_a32 <- 1 + (p0 - 1) / (1 + (w_a32 / wth)^psi)^(1 / psi)}
    s_a32 <- 1 - emax * c_a32 / (ec50 + c_a32)
    w_a33 <- rbase; c_a33 <- 0
    if (t >= td_apo + 33 * td_div) {w_a33 <- delay(cycling_cells, td_apo + 33 * td_div) + delay(damaged_cells1, td_apo + 33 * td_div); c_a33 <- delay(central, td_apo + 33 * td_div) / vc * 1000}
    p_a33 <- 1; if (t - td_apo >= 32.5 * td_div) {p_a33 <- 1 + (p0 - 1) / (1 + (w_a33 / wth)^psi)^(1 / psi)}
    s_a33 <- 1 - emax * c_a33 / (ec50 + c_a33)
    w_a34 <- rbase; c_a34 <- 0
    if (t >= td_apo + 34 * td_div) {w_a34 <- delay(cycling_cells, td_apo + 34 * td_div) + delay(damaged_cells1, td_apo + 34 * td_div); c_a34 <- delay(central, td_apo + 34 * td_div) / vc * 1000}
    p_a34 <- 1; if (t - td_apo >= 33.5 * td_div) {p_a34 <- 1 + (p0 - 1) / (1 + (w_a34 / wth)^psi)^(1 / psi)}
    s_a34 <- 1 - emax * c_a34 / (ec50 + c_a34)
    w_a35 <- rbase; c_a35 <- 0
    if (t >= td_apo + 35 * td_div) {w_a35 <- delay(cycling_cells, td_apo + 35 * td_div) + delay(damaged_cells1, td_apo + 35 * td_div); c_a35 <- delay(central, td_apo + 35 * td_div) / vc * 1000}
    p_a35 <- 1; if (t - td_apo >= 34.5 * td_div) {p_a35 <- 1 + (p0 - 1) / (1 + (w_a35 / wth)^psi)^(1 / psi)}
    s_a35 <- 1 - emax * c_a35 / (ec50 + c_a35)
    w_a36 <- rbase; c_a36 <- 0
    if (t >= td_apo + 36 * td_div) {w_a36 <- delay(cycling_cells, td_apo + 36 * td_div) + delay(damaged_cells1, td_apo + 36 * td_div); c_a36 <- delay(central, td_apo + 36 * td_div) / vc * 1000}
    p_a36 <- 1; if (t - td_apo >= 35.5 * td_div) {p_a36 <- 1 + (p0 - 1) / (1 + (w_a36 / wth)^psi)^(1 / psi)}
    s_a36 <- 1 - emax * c_a36 / (ec50 + c_a36)
    w_a37 <- rbase; c_a37 <- 0
    if (t >= td_apo + 37 * td_div) {w_a37 <- delay(cycling_cells, td_apo + 37 * td_div) + delay(damaged_cells1, td_apo + 37 * td_div); c_a37 <- delay(central, td_apo + 37 * td_div) / vc * 1000}
    p_a37 <- 1; if (t - td_apo >= 36.5 * td_div) {p_a37 <- 1 + (p0 - 1) / (1 + (w_a37 / wth)^psi)^(1 / psi)}
    s_a37 <- 1 - emax * c_a37 / (ec50 + c_a37)
    w_a38 <- rbase; c_a38 <- 0
    if (t >= td_apo + 38 * td_div) {w_a38 <- delay(cycling_cells, td_apo + 38 * td_div) + delay(damaged_cells1, td_apo + 38 * td_div); c_a38 <- delay(central, td_apo + 38 * td_div) / vc * 1000}
    p_a38 <- 1; if (t - td_apo >= 37.5 * td_div) {p_a38 <- 1 + (p0 - 1) / (1 + (w_a38 / wth)^psi)^(1 / psi)}
    s_a38 <- 1 - emax * c_a38 / (ec50 + c_a38)
    w_a39 <- rbase; c_a39 <- 0
    if (t >= td_apo + 39 * td_div) {w_a39 <- delay(cycling_cells, td_apo + 39 * td_div) + delay(damaged_cells1, td_apo + 39 * td_div); c_a39 <- delay(central, td_apo + 39 * td_div) / vc * 1000}
    p_a39 <- 1; if (t - td_apo >= 38.5 * td_div) {p_a39 <- 1 + (p0 - 1) / (1 + (w_a39 / wth)^psi)^(1 / psi)}
    s_a39 <- 1 - emax * c_a39 / (ec50 + c_a39)
    w_a40 <- rbase; c_a40 <- 0
    if (t >= td_apo + 40 * td_div) {w_a40 <- delay(cycling_cells, td_apo + 40 * td_div) + delay(damaged_cells1, td_apo + 40 * td_div); c_a40 <- delay(central, td_apo + 40 * td_div) / vc * 1000}
    p_a40 <- 1; if (t - td_apo >= 39.5 * td_div) {p_a40 <- 1 + (p0 - 1) / (1 + (w_a40 / wth)^psi)^(1 / psi)}
    s_a40 <- 1 - emax * c_a40 / (ec50 + c_a40)
    kin_apo <- kin0 *
      p_a0 * p_a1 * p_a2 * p_a3 * p_a4 * p_a5 * p_a6 * p_a7 * p_a8 * p_a9 *
      p_a10 * p_a11 * p_a12 * p_a13 * p_a14 * p_a15 * p_a16 * p_a17 * p_a18 * p_a19 *
      p_a20 * p_a21 * p_a22 * p_a23 * p_a24 * p_a25 * p_a26 * p_a27 * p_a28 * p_a29 *
      p_a30 * p_a31 * p_a32 * p_a33 * p_a34 * p_a35 * p_a36 * p_a37 * p_a38 * p_a39 *
      p_a40 * s_a1 * s_a2 * s_a3 * s_a4 * s_a5 * s_a6 * s_a7 * s_a8 * s_a9 *
      s_a10 * s_a11 * s_a12 * s_a13 * s_a14 * s_a15 * s_a16 * s_a17 * s_a18 * s_a19 *
      s_a20 * s_a21 * s_a22 * s_a23 * s_a24 * s_a25 * s_a26 * s_a27 * s_a28 * s_a29 *
      s_a30 * s_a31 * s_a32 * s_a33 * s_a34 * s_a35 * s_a36 * s_a37 * s_a38 * s_a39 *
      s_a40
    eff_apo <- emax * c_a0 / (ec50 + c_a0)
    out_apo <- 0
    if (t >= td_apo) {out_apo <- kin_apo * eff_apo}

    # Proliferating pool (Eq. 31 written as kin (1 - E) - kin(t - T)(1 - E(C(t - T))),
    # equivalently Eqs. 28-29) and apoptotic pool (Eq. 32)
    cycling_cells(0) <- rbase
    d/dt(cycling_cells) <- kin * (1 - eff) - kin_lag
    d/dt(damaged_cells1) <- kin * eff - out_apo

    tumor_vol ~ add(addSd_tumor_vol)
  })
}
