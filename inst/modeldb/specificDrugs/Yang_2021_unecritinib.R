Yang_2021_unecritinib <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for the ALK/ROS1 kinase",
    "inhibitor unecritinib (TQ-B3101) and its active metabolite crizotinib",
    "(TQ-B3101M) in adults with advanced solid tumours and adolescents with",
    "relapsed or refractory ALK-positive anaplastic large cell lymphoma (Yang",
    "2021). Unecritinib is described by a one-compartment model with",
    "first-order absorption and first-order elimination; all of its",
    "clearance forms crizotinib by amide hydrolysis (fraction metabolised",
    "fixed to 1 on a molar basis for identifiability), which is described by",
    "a two-compartment model whose apparent clearance decreases exponentially",
    "with time since the first dose, CLm(t) = CLm0 * (1 + famp * (1 -",
    "exp(-k * t))) with famp = -0.41, to a steady-state value 41 percent",
    "below the first-dose value. All clearances and volumes are apparent (X/F",
    "for the parent, X/Fm for the metabolite, i.e. proportional to the true",
    "metabolite values by the unknown fraction metabolised). No covariate was",
    "retained. The paper names the analytes only by development code and",
    "withholds their molecular weights; the maintainers identified them from",
    "PubChem and apply the crizotinib / unecritinib molecular-weight ratio",
    "to the formation flux so that Cc_crizotinib is in crizotinib mass",
    "units. Below-quantitation-limit samples were discarded (M1 method) in",
    "the source analysis."
  )
  reference <- paste(
    "Yang F, Wu H, Bo Y, Lu Y, Pan H, Li S, Lu Q, Xie S, Liao H, Wang B.",
    "Population Pharmacokinetic Modeling and Simulation of TQ-B3101 to",
    "Inform Dosing in Pediatric Patients With Solid Tumors.",
    "Front Pharmacol. 2021;12:782518 (published 18 January 2022).",
    "doi:10.3389/fphar.2021.782518.",
    sep = " "
  )
  vignette <- "Yang_2021_unecritinib"

  units <- list(
    time = "h",
    dosing = "mg",
    concentration = paste(
      "ng/mL for Cc (unecritinib) and for Cc_crizotinib (crizotinib); the",
      "source fitted molar concentrations (LLOQ 1 nmol/L for both analytes,",
      "Yang 2021 Methods, Analytical Methods), converted here with PubChem",
      "molecular weights"
    )
  )

  # Issue #482: what each ODE state holds. Amounts are mg of unecritinib for
  # the parent states and mg of crizotinib for the metabolite states (the
  # formation flux carries the molecular-weight ratio, see model()).
  compartmentData <- list(
    depot = list(
      analyte = "unecritinib",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "unecritinib",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_crizotinib = list(
      analyte = "crizotinib",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_crizotinib = list(
      analyte = "crizotinib",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  # No covariate was retained in the final model (Yang 2021 Results,
  # Population PK Modeling: "the demographic covariates including body weight,
  # body mass index, BSA, age, gender, obesity, albumin, and markers of
  # hepatic and kidney functions, had no meaningful impact on the PK of
  # TQ-B3101 and TQ-B3101M"). The screened set is documented here; none of
  # these columns is referenced in model().
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Screened. Fixed allometric exponents (0.75 on clearances, 1 on",
        "volumes) were also tested and INCREASED the objective function by",
        "4.9 (parent) and 12.6 (metabolite) points, so weight was not",
        "retained (Yang 2021 Results). The paper's paediatric dosing",
        "simulations nevertheless applied those fixed exponents, scaled to",
        "the median adult weight; that scaling is a simulation device, not",
        "part of the fitted model (see the vignette)."
      )
    ),
    HT = list(
      description = "Body height.",
      units = "cm",
      type = "continuous",
      notes = "Screened; not retained (Yang 2021 Methods, covariate list)."
    ),
    BMI = list(
      description = "Body mass index.",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened; not retained (Yang 2021 Methods, covariate list)."
    ),
    BSA = list(
      description = "Body surface area.",
      units = "m^2",
      type = "continuous",
      notes = paste(
        "Screened; not retained (Yang 2021 Methods, covariate list). BSA",
        "defines the paper's paediatric dosing tiers (Table 5) but does not",
        "enter the fitted model."
      )
    ),
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      notes = "Screened; not retained (Yang 2021 Methods, covariate list)."
    ),
    ALB = list(
      description = "Serum albumin.",
      units = "g/L",
      type = "continuous",
      notes = "Screened; not retained (Yang 2021 Methods, covariate list)."
    ),
    CREAT = list(
      description = "Serum creatinine.",
      units = "umol/L",
      type = "continuous",
      notes = paste(
        "Screened. Significant in the univariate analysis on CL/F of",
        "TQ-B3101 but not retained because the association was positive and",
        "judged not physiologically meaningful (Yang 2021 Results)."
      )
    ),
    ALT = list(
      description = "Alanine aminotransferase.",
      units = "U/L",
      type = "continuous",
      notes = "Screened; not retained (Yang 2021 Methods, covariate list)."
    ),
    AST = list(
      description = "Aspartate aminotransferase.",
      units = "U/L",
      type = "continuous",
      notes = "Screened; not retained (Yang 2021 Methods, covariate list)."
    ),
    TBILI = list(
      description = "Total bilirubin.",
      units = "umol/L",
      type = "continuous",
      notes = paste(
        "Screened. Significant in the univariate analysis on V/F of",
        "TQ-B3101 and on Vcm/Fm of TQ-B3101M but not significant enough to be",
        "retained at the backward-elimination step (Yang 2021 Results)."
      )
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate.",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = "Screened; not retained (Yang 2021 Methods, covariate list)."
    ),
    DIS_OBESE = list(
      description = "Obesity indicator (1 = obese).",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not retained (Yang 2021 Methods, covariate list)."
    ),
    SEXF = list(
      description = "Sex indicator (1 = female).",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not retained (Yang 2021 Methods, covariate list)."
    ),
    ADOLESCENT = list(
      description = "Study population indicator (1 = adolescent, 0 = adult).",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened. Significant in the univariate analysis on CL/F of",
        "TQ-B3101 but not significant enough to be retained (Yang 2021",
        "Results)."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 40L,
    n_studies = 2L,
    n_observations = paste(
      "375 quantifiable TQ-B3101 and 658 quantifiable TQ-B3101M",
      "concentrations (340 and 42 BQL samples discarded)"
    ),
    age_range = "11-73 years (adults 28-73, adolescents 11-14)",
    age_median = "49.5 years",
    weight_range = "32.9-87.7 kg",
    weight_median = "58.3 kg (adults 59.0 kg, adolescents 41.8 kg)",
    sex_female_pct = 47.5,
    race_ethnicity = "Chinese patients (studies run at Chinese centres; race not tabulated)",
    disease_state = paste(
      "34 adults with advanced malignant solid tumours (Phase 1) or",
      "relapsed/refractory ALK-positive anaplastic large cell lymphoma",
      "(Phase 2), and 6 adolescents with relapsed/refractory ALK-positive",
      "ALCL (Phase 2)."
    ),
    dose_range = paste(
      "Oral TQ-B3101 under fasting conditions: single doses of 100 or 200 mg;",
      "100, 200 or 300 mg once daily for 28 days (Phase 1); 200, 250, 300 or",
      "350 mg twice daily for 28 days (Phase 2)."
    ),
    regions = "China",
    notes = paste(
      "Pooled Phase 1 dose-escalation study NCT03019276 and Phase 2",
      "single-arm study NCT04306887 (Yang 2021 Table 1). Baseline",
      "demographics: Yang 2021 Table 2 (21/40 = 52.5% male). LLOQ 1 nmol/L",
      "for both analytes; BQL samples were discarded (M1). Sequential fit in",
      "NONMEM 7.4 FOCE-I: parent parameters were estimated first and fixed in",
      "the combined parent-metabolite fit."
    )
  )

  ini({
    # -----------------------------------------------------------------
    # Unecritinib (TQ-B3101, parent). Yang 2021 Table 3, final oral PK model. The
    # half-life check ln(2) * V/F / (CL/F) = 0.693 * 4200 / 2850 = 1.02 h
    # reproduces the 1.0 h quoted in the Abstract and Results.
    # -----------------------------------------------------------------
    lcl <- log(2850); label("Unecritinib (TQ-B3101) apparent clearance CL/F (L/h)")                    # Yang 2021 Table 3: CL/F = 2,850 L/h (RSE 7%; bootstrap median 2854.3, 95% CI 2452.6-3290.2)
    lvc <- log(4200); label("Unecritinib (TQ-B3101) apparent volume of distribution V/F (L)")          # Yang 2021 Table 3: V/F = 4,200 L (RSE 9%; bootstrap median 4205.5, 95% CI 3506.4-5061.4)
    lka <- log(51.9); label("Unecritinib (TQ-B3101) first-order absorption rate constant Ka (1/h)")    # Yang 2021 Table 3: Ka = 51.9 1/h (RSE 67%; bootstrap median 46.0, 95% CI 15.2-165.7)

    # Fraction of TQ-B3101 clearance forming TQ-B3101M, fixed to 1 so that the
    # metabolite sub-model is identifiable; every metabolite parameter below is
    # therefore an apparent value X/Fm.
    fm <- fixed(1); label("Molar fraction of unecritinib clearance forming crizotinib (unitless)")  # Yang 2021 Methods: 'the fraction of TQ-B3101 to TQ-B3101M was fixed to one to obtain an identifiable model'; Figure 1 caption

    # -----------------------------------------------------------------
    # Crizotinib (TQ-B3101M, active metabolite). Yang 2021 Table 3.
    # -----------------------------------------------------------------
    lcl_crizotinib <- log(126);  label("Crizotinib (TQ-B3101M) apparent clearance at time 0, CLm0/Fm (L/h)")          # Yang 2021 Table 3: CLm0/Fm = 126 L/h (RSE 11%; bootstrap median 128.9, 95% CI 101.3-207.6)
    lvc_crizotinib <- log(2300); label("Crizotinib (TQ-B3101M) apparent central volume Vcm/Fm (L)")                   # Yang 2021 Table 3: Vcm/Fm = 2,300 L (RSE 9%; bootstrap median 2277.9, 95% CI 1925.4-2673.0)
    lq_crizotinib  <- log(113);  label("Crizotinib (TQ-B3101M) apparent inter-compartmental clearance Qm/Fm (L/h)")   # Yang 2021 Table 3: Qm/Fm = 113 L/h (RSE 18%; bootstrap median 106.6, 95% CI 70.9-201.0)
    lvp_crizotinib <- log(1480); label("Crizotinib (TQ-B3101M) apparent peripheral volume Vpm/Fm (L)")                # Yang 2021 Table 3: Vpm/Fm = 1,480 L (RSE 25%; bootstrap median 1498.0, 95% CI 564.2-2447.4)

    # Time-dependent clearance of TQ-B3101M (Yang 2021 Table 3 footnote a):
    # CLm/Fm = CLm0/Fm * [1 - TDPK * (1 - exp(-KTDPK * T))]. Encoded with the
    # signed fractional-amplitude token, cl_exp_famp = -TDPK (negative = decline).
    cl_exp_famp <- -0.41;      label("Fractional amplitude of the crizotinib (TQ-B3101M) clearance change, -TDPK (unitless)")  # Yang 2021 Table 3: TDPK on CLm/Fm = 0.41 (RSE 16%; bootstrap median 0.43, 95% CI 0.28-0.65); cl_exp_famp = -TDPK
    lcl_exp_kdes <- log(0.0363); label("First-order rate constant of the crizotinib (TQ-B3101M) clearance decrease, KTDPK (1/h)")  # Yang 2021 Table 3: KTDPK = 0.0363 1/h (RSE 26%; bootstrap median 0.037, 95% CI 0.024-0.12)

    # -----------------------------------------------------------------
    # Inter-individual variability. Yang 2021 Methods Eq. 1 (log-normal,
    # theta_i = theta_TV * exp(eta_i)); Table 3 reports each as %CV, converted
    # with omega^2 = log(1 + CV^2). No off-diagonal elements are reported, so
    # the etas are independent. Ka and Qm/Fm carry no IIV.
    # -----------------------------------------------------------------
    etalcl          ~ 0.0760   # Yang 2021 Table 3: eta CL/F 28.1 %CV; log(1 + 0.281^2) = 0.07600
    etalvc          ~ 0.1010   # Yang 2021 Table 3: eta V/F 32.6 %CV; log(1 + 0.326^2) = 0.10100
    etalcl_crizotinib ~ 0.1100   # Yang 2021 Table 3: eta CLm/Fm 34.1 %CV; log(1 + 0.341^2) = 0.11000
    etalvc_crizotinib ~ 0.2484   # Yang 2021 Table 3: eta Vcm/Fm 53.1 %CV; log(1 + 0.531^2) = 0.24839
    etalvp_crizotinib ~ 0.5300   # Yang 2021 Table 3: eta Vpm/Fm 83.6 %CV; log(1 + 0.836^2) = 0.52998

    # -----------------------------------------------------------------
    # Residual error. Table 3 reports one %CV per analyte (sigma_1p, sigma_1m),
    # i.e. a proportional error model.
    # -----------------------------------------------------------------
    propSd           <- 0.711; label("Unecritinib (TQ-B3101) proportional residual SD (fraction)")    # Yang 2021 Table 3: sigma_1p = 71.1 %CV for TQ-B3101 (RSE 3%; bootstrap median 70.7, 95% CI 66.1-75.8)
    propSd_crizotinib  <- 0.319; label("Crizotinib (TQ-B3101M) proportional residual SD (fraction)")   # Yang 2021 Table 3: sigma_1m = 31.9 %CV for TQ-B3101M (RSE 6%; bootstrap median 31.6, 95% CI 28.4-35.9)
  })

  model({
    # 1. Individual parameters (Yang 2021 Eq. 1).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    cl_crizotinib0 <- exp(lcl_crizotinib + etalcl_crizotinib)
    vc_crizotinib  <- exp(lvc_crizotinib + etalvc_crizotinib)
    q_crizotinib   <- exp(lq_crizotinib)
    vp_crizotinib  <- exp(lvp_crizotinib + etalvp_crizotinib)

    # 2. Time-dependent crizotinib clearance (Yang 2021 Table 3 footnote a).
    #    T is time since the first dose; here t, so simulations must start at
    #    the first dose. CLm falls from CLm0 at t = 0 to CLm0 * (1 +
    #    cl_exp_famp) = 74.3 L/h (typical) as t -> infinity, with half-time
    #    ln(2) / 0.0363 = 19.1 h.
    cl_exp_kdes <- exp(lcl_exp_kdes)
    cl_crizotinib <- cl_crizotinib0 * (1 + cl_exp_famp * (1 - exp(-cl_exp_kdes * t)))

    # 3. Micro-constants. All unecritinib clearance forms crizotinib (fm = 1
    #    on a molar basis; Figure 1 draws the (1 - FM) * CL/F elimination
    #    arm, which is zero).
    kel      <- cl / vc
    kform    <- fm * kel
    k12_m    <- q_crizotinib / vc_crizotinib
    k21_m    <- q_crizotinib / vp_crizotinib
    kel_m    <- cl_crizotinib / vc_crizotinib

    # Molecular weights (g/mol), not reported in Yang 2021 ("structure,
    # molecular formula and molecular weight cannot be provided"). PubChem:
    # unecritinib CID 71506874, C23H24Cl2FN5O2, 492.4 (synonym TQ-B3101);
    # crizotinib CID 11626560, C21H22Cl2FN5O, 450.3. Unecritinib is the
    # N-acetyl amide of crizotinib, and TQ-B3101 is converted to TQ-B3101M
    # "by amide hydrolysis" (Yang 2021 Discussion). The source fitted molar
    # data (LLOQ 1 nmol/L), so fm = 1 is one mole of crizotinib per mole of
    # unecritinib cleared; the ratio turns that into a mass flux. It matches
    # the paper's own 300 mg unecritinib ~ 275 mg crizotinib equivalence
    # (275 / 300 = 0.917 vs 450.3 / 492.4 = 0.914).
    mw_unecritinib <- 492.4
    mw_crizotinib  <- 450.3

    # 4. ODE system (Yang 2021 Figure 1); amounts in mg of each analyte.
    d/dt(depot)                  <- -ka * depot
    d/dt(central)                <-  ka * depot - kel * central
    d/dt(central_crizotinib)     <-  kform * central * mw_crizotinib / mw_unecritinib -
                                     (kel_m + k12_m) * central_crizotinib +
                                     k21_m * peripheral1_crizotinib
    d/dt(peripheral1_crizotinib) <-  k12_m * central_crizotinib - k21_m * peripheral1_crizotinib

    # 5. Observations. Dose in mg and volumes in L give mg/L; x 1000 -> ng/mL.
    Cc            <- 1000 * central / vc
    Cc_crizotinib <- 1000 * central_crizotinib / vc_crizotinib

    Cc          ~ prop(propSd)
    Cc_crizotinib ~ prop(propSd_crizotinib)
  })
}
