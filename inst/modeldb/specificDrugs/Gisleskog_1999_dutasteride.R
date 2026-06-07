Gisleskog_1999_dutasteride <- function() {
  description <- "Two-compartment population PK model for dutasteride (GI198745, a dual type-1/type-2 5-alpha-reductase inhibitor) in healthy male volunteers after single oral doses, with first-order absorption, an absorption lag-time, and parallel linear (CL_l) plus Michaelis-Menten (Vmax / Km) elimination from the central compartment (Gisleskog 1999). All volumes and clearances are apparent (oral, no IV reference); bioavailability is assumed dose-independent."
  reference <- "Olsson Gisleskog P, Hermann D, Hammarlund-Udenaes M, Karlsson MO. The pharmacokinetic modelling of GI198745 (dutasteride), a compound with parallel linear and nonlinear elimination. Br J Clin Pharmacol. 1999;47(1):53-58. doi:10.1046/j.1365-2125.1999.00843.x"
  vignette <- "Gisleskog_1999_dutasteride"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 32L,
    n_studies      = 1L,
    age_range      = "20-57 years (parent 48-subject cohort; same range applies to the 32 GI198745 recipients)",
    age_median     = "37 years (parent cohort median)",
    weight_range   = "56.3-102 kg (parent cohort)",
    weight_median  = "76.1 kg (parent cohort)",
    sex_female_pct = 0,
    race_ethnicity = "Not reported",
    disease_state  = "Healthy male volunteers in a single-dose escalation study.",
    dose_range     = "Single oral solution dose of 0.01, 0.1, 1, 2.5, 5, 10, 20, or 40 mg of GI198745 (4 subjects per dose level). Vehicle: PEG400 + TWEEN80 0.01% in 7.5 mL water (15 mL for the 40 mg dose), administered with 240 mL of water.",
    regions        = "United Kingdom (study approved by the institutional review board of Besselaar Ltd.).",
    notes          = "Source: Gisleskog 1999 Methods, p54. The parent study enrolled 48 healthy males aged 20-57 (median 37 y, weight 56.3-102 kg, median 76.1 kg) to receive GI198745, finasteride, or placebo; the 32 GI198745 recipients comprise the modelling analysis. Sampling: predose and at 0.5, 1, 2, 3, 4, 6, 8, 12, 16, 24 hours and 2, 3, 7, 14, 21, 28 days post-dose (and 56 days in three subjects). The 0.01 mg dose was below the LC/MS assay LLOQ (0.1 ng/mL) so produced no measurable concentrations; all other dose levels contributed observations. No covariates were tested -- the cohort was too small per the authors (Methods, Variability model, p55 and Discussion, p57)."
  )

  ini({
    # Structural PK parameters - Gisleskog 1999 Table 1 (p57), final model
    # (Model 3: 2-compartment with first-order absorption, absorption lag-time,
    # and parallel linear + saturable Michaelis-Menten elimination). All
    # volumes and clearances are apparent (the publication carries no IV
    # reference; bioavailability is assumed dose-independent, Methods p54).
    lka    <- log(2.41);   label("First-order absorption rate ka (1/h)")                          # Table 1: ka 2.41 /h (SE 22%)
    ltlag  <- log(0.316);  label("Absorption lag-time tlag (h)")                                  # Table 1: tlag 0.316 h (SE 22%)
    lcl    <- log(0.583);  label("Linear apparent clearance CL_l (L/h)")                          # Table 1: CL_l 0.583 L/h (SE 30%)
    lvc    <- log(173);    label("Apparent central volume of distribution Vc (L)")                # Table 1: Vc 173 L (SE 12%)
    lq     <- log(33.5);   label("Inter-compartmental clearance Q (L/h)")                         # Table 1: Q 33.5 L/h (SE 16%)
    lvp    <- log(338);    label("Apparent peripheral volume of distribution Vp (L)")             # Table 1: Vp 338 L (SE 13%). Vss = Vc + Vp = 511 L (calculated, footnote a).
    lvmax  <- log(5.91e-3);label("Saturable elimination capacity Vmax (mg/h)")                    # Table 1: Vmax 5.91 ug/h (SE 39%); the table's printed unit is mu-g h^-1 (micrograms/h), which pdftotext / docling render as 'mg' because pdfminer's ToUnicode CMap collapses the Symbol-font mu glyph to ASCII 'm'. The paper's prose Vmax/Km = 6.2 L/h confirms ug/h not mg/h: 5.91 ug/h / 0.957 ng/mL = 6.17 L/h. Encoded here as 0.00591 mg/h to keep the (mg) state derivative dimensionally consistent.
    lkm    <- log(0.957);  label("Michaelis-Menten constant Km (ng/mL)")                          # Table 1: Km 0.957 ng/mL (SE 47%)

    # IIV - Gisleskog 1999 Table 1, "Intersubject variability (CV%)" column.
    # No IIV is reported on Km; per Methods (Variability model, p55) variability
    # was initially applied to every parameter and dropped when it approached
    # zero, which is the reason Km has no eta. Off-diagonal omega elements
    # are not reported (diagonal OMEGA in NONMEM). Convert each CV% to a
    # log-normal variance via omega^2 = log(1 + CV^2):
    #   ka     CV 70% -> omega^2 = log(1 + 0.70^2) = 0.3988
    #   tlag   CV 38% -> omega^2 = log(1 + 0.38^2) = 0.1350
    #   CL_l   CV 69% -> omega^2 = log(1 + 0.69^2) = 0.3884
    #   Vmax   CV 43% -> omega^2 = log(1 + 0.43^2) = 0.1693
    #   Q      CV 62% -> omega^2 = log(1 + 0.62^2) = 0.3253
    #   Vc     CV 32% -> omega^2 = log(1 + 0.32^2) = 0.0975
    #   Vp     CV 41% -> omega^2 = log(1 + 0.41^2) = 0.1554
    etalka    ~ 0.3988                                                                            # Table 1: ka IIV 70% CV (SE of variability 77%)
    etaltlag  ~ 0.1350                                                                            # Table 1: tlag IIV 38% CV (SE of variability 53%)
    etalcl    ~ 0.3884                                                                            # Table 1: CL_l IIV 69% CV (SE of variability 30%)
    etalvmax  ~ 0.1693                                                                            # Table 1: Vmax IIV 43% CV (SE of variability 26%)
    etalq     ~ 0.3253                                                                            # Table 1: Q IIV 62% CV (SE of variability 16%)
    etalvc    ~ 0.0975                                                                            # Table 1: Vc IIV 32% CV (SE of variability 19%)
    etalvp    ~ 0.1554                                                                            # Table 1: Vp IIV 41% CV (SE of variability 19%)

    # Residual error - Gisleskog 1999 Eq. 8 (p55): C = C_hat * (1 + eps),
    # eps ~ N(0, sigma^2). Table 1 reports sigma = 0.13 and Results (p57)
    # describes this as "residual variability... at 13%", so 0.13 is the SD
    # on the proportional scale and maps directly to propSd = 0.13 in
    # nlmixr2's prop(...) error model.
    propSd <- 0.13;   label("Proportional residual error on Cc (fraction)")                       # Table 1: sigma 0.13 = 13% CV proportional
  })

  model({
    # Individual PK parameters (typical-value scale with log-normal eta).
    ka    <- exp(lka    + etalka)
    tlag  <- exp(ltlag  + etaltlag)
    cl    <- exp(lcl    + etalcl)
    vc    <- exp(lvc    + etalvc)
    q     <- exp(lq     + etalq)
    vp    <- exp(lvp    + etalvp)
    vmax  <- exp(lvmax  + etalvmax)
    km    <- exp(lkm)

    # Two-compartment with first-order oral absorption, an absorption lag-time,
    # and parallel linear plus saturable Michaelis-Menten elimination from the
    # central compartment (Gisleskog 1999 Eqs. 1-3 and reparameterisations
    # Eqs. 4-6, p54-55):
    #   dA1/dt = -ka * A1
    #   dA2/dt =  ka * A1 - k23 * A2 + k32 * A3 - k20 * A2 - Vmax * (A2/Vc) / (Km + A2/Vc)
    #   dA3/dt =  k23 * A2 - k32 * A3
    # with k20 = CL_l / Vc, k23 = Q / Vc, k32 = Q / Vp.
    #
    # Units. Dose is in mg so central / vc has units mg / L = ug/mL = 1000 ng/mL.
    # Gisleskog 1999 reports Cc and Km in ng/mL, so Cc = 1000 * central / vc.
    # With Cc and km both in ng/mL the Michaelis-Menten ratio is dimensionless
    # and vmax * Cc / (km + Cc) is in mg/h, matching the (mg) state derivative.
    Cc <- 1000 * central / vc

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          (cl / vc) * central -
                          vmax * Cc / (km + Cc) -
                          (q / vc) * central +
                          (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central -
                          (q / vp) * peripheral1

    alag(depot) <- tlag

    Cc ~ prop(propSd)
  })
}
