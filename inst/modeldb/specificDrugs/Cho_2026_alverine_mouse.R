Cho_2026_alverine_mouse <- function() {
  description <- "Preclinical (mouse). Integrated parent-metabolite PK model for alverine and its three major metabolites, 4-hydroxy alverine (M1), 4-hydroxy alverine glucuronide (M2) and N-desethyl alverine (M3), in male ICR mice after a single IV bolus (5 mg/kg alverine citrate) or oral gavage (75 mg/kg) dose (Cho 2026). Two-compartment parent; one-compartment M1, M2 and M3. Parent central loss (CL/Vp1) is partitioned into a residual elimination rate constant, formation of M1, and formation of M3 (derived as the remainder); M1 forms M2 sequentially. An oral dose is split by fixed and estimated fractions into three first-order input pathways: to parent central (fraction 0.04), and apparent presystemic inputs directly to M1 (0.92) and M3 (0.0399). Amounts in umol of alverine PER ANIMAL (dose in umol/kg x body weight; the source mice weighed 35-42 g); concentrations in umol/L. For an oral dose give the same amount to depot, depot_m1 and depot_m3 -- the bioavailability fractions split it. Typical-value model: the source reports no IIV or residual-error magnitudes."
  reference <- "Cho A, Jeong HC, Kim M, Cho I, Na HJ, Ko KC, Chae YJ, Lee KR. Development of Integrated Parent-Metabolite Pharmacokinetic Model With Apparent Presystemic Metabolism to Characterize the Disposition of Alverine. CPT Pharmacometrics Syst Pharmacol. 2026;15(9):e70342. doi:10.1002/psp4.70342"
  vignette <- "Cho_2026_alverine_mouse"
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")
  dosing <- c("central", "depot", "depot_m1", "depot_m3")
  covariateData <- list()
  compartmentData <- list(
    depot = list(analyte = "alverine", units = "umol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "alverine", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "alverine", units = "umol", specimen = "plasma", verified = TRUE),
    depot_m1 = list(
      analyte = "alverine (oral dose fraction appearing as M1 through apparent presystemic metabolism)",
      units = "umol",
      specimen = "administration site",
      verified = TRUE
    ),
    central_m1 = list(analyte = "4-hydroxy alverine", units = "umol", specimen = "plasma", verified = TRUE),
    central_m2 = list(analyte = "4-hydroxy alverine glucuronide", units = "umol", specimen = "plasma", verified = TRUE),
    depot_m3 = list(
      analyte = "alverine (oral dose fraction appearing as M3 through apparent presystemic metabolism)",
      units = "umol",
      specimen = "administration site",
      verified = TRUE
    ),
    central_m3 = list(analyte = "N-desethyl alverine", units = "umol", specimen = "plasma", verified = TRUE)
  )
  population <- list(
    species = "mouse (ICR)",
    n_subjects = 10L,
    n_studies = 1L,
    age_range = "8-9 weeks",
    weight_range = "35-42 g",
    sex_female_pct = 0,
    disease_state = "healthy, fasted 8 h before dosing",
    dose_range = "single IV bolus 5 mg/kg alverine citrate (10.56 umol/kg; n = 5) or single oral gavage 75 mg/kg alverine citrate (158.37 umol/kg; n = 5)",
    regions = "Republic of Korea (KRIBB)",
    notes = "Male ICR mice, serial sampling 0.083-24 h (IV) and 0.25-24 h (PO); plasma alverine, M1, M2 and M3 by LC-MS/MS, concentrations below the LLOQ excluded (Methods, In Vivo Pharmacokinetic Study Design and LC-MS/MS Analysis). Fitted in Monolix 2024R1 (SAEM)."
  )
  ini({
    # Parent (alverine) disposition -- Table 3 section 1
    lcl <- log(0.512); label("Alverine clearance from the central compartment CL (L/h)") # Table 3 'CL' = 0.512 L/h (RSE 1.993%)
    lvc <- log(0.443); label("Alverine central volume Vp1 (L)") # Table 3 'Vp1' = 0.443 L (RSE 5.556%)
    lvp <- log(0.29); label("Alverine peripheral volume Vp2 (L)") # Table 3 'Vp2' = 0.29 L (RSE 2.445%)
    lq <- log(0.372); label("Alverine intercompartmental clearance Q (L/h)") # Table 3 'Q' = 0.372 L/h (RSE not reported)

    # Parent elimination and metabolism -- Table 3 section 2
    lkel_other <- log(0.0249); label("Alverine residual (non-M1, non-M3) elimination rate constant kp (1/h)") # Table 3 'kp' = 0.0249 1/h (RSE not reported)
    lkmet_m1 <- log(0.582); label("Formation rate constant of M1 from parent central km1 (1/h)") # Table 3 'km1' = 0.582 1/h (RSE 12.313%)

    # Oral input pathways -- Table 3 section 3
    lka <- fixed(log(6.9)); label("Parent first-order absorption rate constant ka (1/h)") # Table 3 'ka' = 6.9 1/h, dagger = fixed; Methods, Parameter Estimation
    lka_m1 <- log(2.183); label("First-order presystemic input rate constant to M1 ka,m1 (1/h)") # Table 3 'kam1' = 2.183 1/h (RSE not reported)
    lka_m3 <- log(823.25); label("First-order presystemic input rate constant to M3 ka,m3 (1/h)") # Table 3 'kam3' = 823.25 1/h (RSE 8.82E3%)
    lfdepot <- fixed(log(0.04)); label("Fraction of the oral dose entering the parent pathway fp (unitless)") # Table 3 'fp' = 0.04, dagger = fixed to the NCA bioavailability; Methods, Parameter Estimation
    lfdepot_m3 <- log(0.0399); label("Fraction of the oral dose entering the M3 presystemic pathway fm3 (unitless)") # Table 3 'fm3' = 0.0399

    # M1 kinetics -- Table 3 section 4
    lvc_m1 <- log(47.905); label("Apparent volume of M1 Vm1 (L)") # Table 3 'Vm1' = 47.905 L (RSE 56.178%)
    lkel_m1 <- log(1.79e-8); label("Elimination rate constant of M1 other than conversion to M2 km1e (1/h)") # Table 3 'km1e' = 1.79E-8 1/h (RSE 1.09E10%)
    lkmet_m2 <- log(1.064); label("Formation rate constant of M2 from M1 km2 (1/h)") # Table 3 'km2' = 1.064 1/h (derived) = 'km1,total' 1.064 (RSE 48.2%) - 'km1e'; File S1 k_m2 = k_m1_total - k_m1e

    # M2 kinetics -- Table 3 section 5
    lvc_m2 <- log(8.1e-5); label("Apparent volume of M2 Vm2 (L)") # Table 3 'Vm2' = 8.1E-5 L (RSE 20.551%)
    lkel_m2 <- log(1874.08); label("Elimination rate constant of M2 km2,total (1/h)") # Table 3 'km2,total' = 1874.08 1/h (RSE 15.576%)

    # M3 kinetics -- Table 3 section 6
    lvc_m3 <- log(2.287); label("Apparent volume of M3 Vm3 (L)") # Table 3 'Vm3' = 2.287 L (RSE 55.264%)
    lkel_m3 <- log(1.683); label("Elimination rate constant of M3 km3,total (1/h)") # Table 3 'km3,total' = 1.683 1/h (RSE 7.444%)

    # Residual error: forms from Methods (combined additive + proportional for
    # the parent, proportional for each metabolite); magnitudes not reported.
    propSd <- fixed(0); label("Proportional residual error on alverine (fraction; not reported in source)") # Methods, Parameter Estimation: combined error for parent, no value reported
    addSd <- fixed(0); label("Additive residual error on alverine (umol/L; not reported in source)") # Methods, Parameter Estimation: combined error for parent, no value reported
    propSd_m1 <- fixed(0); label("Proportional residual error on M1 (fraction; not reported in source)") # Methods, Parameter Estimation: proportional error, no value reported
    propSd_m2 <- fixed(0); label("Proportional residual error on M2 (fraction; not reported in source)") # Methods, Parameter Estimation: proportional error, no value reported
    propSd_m3 <- fixed(0); label("Proportional residual error on M3 (fraction; not reported in source)") # Methods, Parameter Estimation: proportional error, no value reported
  })
  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)
    kel_other <- exp(lkel_other)
    kmet_m1 <- exp(lkmet_m1)
    ka <- exp(lka)
    ka_m1 <- exp(lka_m1)
    ka_m3 <- exp(lka_m3)
    fdepot <- exp(lfdepot)
    fdepot_m3 <- exp(lfdepot_m3)
    vc_m1 <- exp(lvc_m1)
    kel_m1 <- exp(lkel_m1)
    kmet_m2 <- exp(lkmet_m2)
    vc_m2 <- exp(lvc_m2)
    kel_m2 <- exp(lkel_m2)
    vc_m3 <- exp(lvc_m3)
    kel_m3 <- exp(lkel_m3)

    # Equations 2-4 and File S1: total parent central loss kp,total = CL/Vp1
    # is partitioned into kp + km1 + km3, with km3 the remainder.
    kel <- cl / vc
    kmet_m3 <- kel - kel_other - kmet_m1
    k12 <- q / vc
    k21 <- q / vp

    # Equation 1: fp + fm1 + fm3 = 1 (File S1: f_m1 = (1 - f_p) - f_m3)
    fdepot_m1 <- 1 - fdepot - fdepot_m3

    # Equations 6-13
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(depot_m1) <- -ka_m1 * depot_m1
    d/dt(central_m1) <- kmet_m1 * central + ka_m1 * depot_m1 - (kel_m1 + kmet_m2) * central_m1
    d/dt(central_m2) <- kmet_m2 * central_m1 - kel_m2 * central_m2
    d/dt(depot_m3) <- -ka_m3 * depot_m3
    d/dt(central_m3) <- kmet_m3 * central + ka_m3 * depot_m3 - kel_m3 * central_m3

    f(depot) <- fdepot
    f(depot_m1) <- fdepot_m1
    f(depot_m3) <- fdepot_m3

    # Equation 14
    Cc <- central / vc
    Cc_m1 <- central_m1 / vc_m1
    Cc_m2 <- central_m2 / vc_m2
    Cc_m3 <- central_m3 / vc_m3

    Cc ~ add(addSd) + prop(propSd)
    Cc_m1 ~ prop(propSd_m1)
    Cc_m2 ~ prop(propSd_m2)
    Cc_m3 ~ prop(propSd_m3)
  })
}
