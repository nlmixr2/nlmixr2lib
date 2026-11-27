Nayak_2026_marstacimab <- function() {
  description <- "Two-compartment target-mediated drug disposition (quasi-steady-state approximation) population PK model for marstacimab, an anti-TFPI IgG1 monoclonal antibody, with first-order SC absorption after a lag, linear clearance, an additional Michaelis-Menten clearance of free drug, and turnover of total soluble TFPI, in healthy volunteers and adolescents and adults with hemophilia A or B (Nayak 2026). A sequential sigmoid Emax model links model-predicted free TFPI to peak thrombin in participants with hemophilia."
  reference <- "Nayak S, Suzuki A, Ravva P, Raje S. Clinical Pharmacokinetics and Pharmacodynamics of Marstacimab, an Anti-tissue Factor Pathway Inhibitor Monoclonal Antibody, in Adolescent and Adult Participants with Hemophilia. Clin Pharmacokinet. 2026;65(10):1591-1612. doi:10.1007/s40262-026-01695-5"
  vignette <- "Nayak_2026_marstacimab"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power model normalised to the median weight of the analysis set, 71.1 kg (Table 3, all studies combined; Supplementary eq 8.1 'WT/WTmedian'). Estimated exponents 1.14 on CL and 1.86 on Vc (Table 5). Weight also enters Q and Vp (Methods; Fig. 1 caption), but no exponent for them is reported; Q takes the CL exponent and Vp the Vc exponent.",
      source_name = "WT"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-volunteer indicator, 1 = healthy volunteer, 0 = participant with hemophilia",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (participant with hemophilia)",
      notes = "Linear effect on CL: CL * (1 + 0.199 * DIS_HEALTHY) (Table 5; Supplementary eq 7.1). Healthy volunteers come from studies NCT02531815 and NCT04832139 (n = 63).",
      source_name = "participant status (healthy vs hemophilia)"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator, 1 = Asian, 0 = non-Asian",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian)",
      notes = "Linear effect on CL: CL * (1 + 0.346 * RACE_ASIAN) (Table 5). Table 4 lists Native Hawaiian participants separately from Asian participants; they are non-Asian here.",
      source_name = "race (Asian vs non-Asian)"
    ),
    ADA_POS = list(
      description = "Anti-drug antibody status, 1 = ADA positive, 0 = ADA negative",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ADA negative)",
      notes = "Linear effect on CL: CL * (1 + 0.229 * ADA_POS) (Table 5). Coded per participant (Table 4, 'ADA status').",
      source_name = "ADA status"
    ),
    RENALIMP_MILD = list(
      description = "Mild renal impairment indicator, 1 = eGFR 60 to < 90 mL/min/1.73 m^2, 0 = normal renal function",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no renal impairment)",
      notes = "Linear effect on CL: CL * (1 - 0.168 * RENALIMP_MILD) (Table 5). Results define mild renal impairment as eGFR < 90 mL/min/1.73 m^2; the analysis set had no moderate or severe renal impairment (Table 4).",
      source_name = "baseline renal impairment (mild)"
    ),
    DIS_HEMOPHILIA_B = list(
      description = "Hemophilia B indicator, 1 = hemophilia B, 0 = hemophilia A or no hemophilia",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (hemophilia A; healthy volunteers also carry 0 and take their own DIS_HEALTHY effect)",
      notes = "Linear effect on CL: CL * (1 + 0.118 * DIS_HEMOPHILIA_B) (Table 5). The hemophilia A subgroup is the reference population (Table 5 footnote).",
      source_name = "hemophilia type (HB)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "marstacimab", units = "nmol", specimen = "administration site", verified = TRUE),
    central = list(
      analyte = "total (free + TFPI-bound) marstacimab",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(analyte = "free marstacimab", units = "nmol", specimen = "tissue", verified = TRUE),
    total_target = list(
      analyte = "total (free + marstacimab-bound) soluble TFPI",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 213L,
    n_studies = 6L,
    n_healthy = 63L,
    n_hemophilia = 150L,
    age_range = "13-66 years",
    age_median = "33 years",
    weight_range = "35-120 kg",
    weight_median = "71.1 kg",
    sex_female_pct = 0,
    race_ethnicity = c(White = 33.3, Black = 7.5, Asian = 32.4, `Native Hawaiian` = 26.3, Other = 0.5),
    disease_state = "Healthy male volunteers (n = 63) and male adolescents (n = 19) and adults with severe hemophilia A (n = 121) or moderately severe to severe hemophilia B (n = 29) without inhibitors (the phase Ib/II study also enrolled participants with inhibitors).",
    dose_range = "SC 30-450 mg single or once weekly (phase III: 300 mg SC loading dose then 150 mg SC QW, escalation to 300 mg QW allowed); IV 150 and 440 mg single dose in healthy volunteers.",
    regions = "Multinational",
    renal_function = "78.4% normal, 21.6% mild impairment (eGFR < 90 mL/min/1.73 m^2); none moderate or severe.",
    ada_positive_pct = 23.6,
    tfpi_baseline = "Total TFPI mean (SD) 132.5 (28.5) ng/mL (Table 3); typical model baseline 3.75 nM (Table 5).",
    notes = "Pooled from 6 studies: NCT02531815, NCT04832139, NCT04878731 (phase I), NCT02974855 (phase Ib/II), NCT03363321 (phase II) and NCT03938792 (phase III, non-inhibitor cohort, Day 360 data excluded). Baseline demographics in Tables 3 and 4. The peak-thrombin Emax model was fitted only to participants with hemophilia; its typical values do not describe healthy volunteers."
  )

  ini({
    # Disposition (Table 5; reference weight 71.1 kg, hemophilia A, non-Asian, ADA negative, normal renal function)
    lcl <- log(0.0188); label("Linear clearance CL (L/h)") # Table 5: CL = 0.0188 L/h (RSE 11.3%)
    lvc <- log(3.61); label("Central volume of distribution Vc (L)") # Table 5: Vc = 3.61 L (RSE 9.15%)
    lq <- log(0.00489); label("Inter-compartmental clearance Q (L/h)") # Table 5: Q = 0.00489 L/h (RSE 10.2%)
    lvp <- log(4.99); label("Peripheral volume of distribution Vp (L)") # Table 5: Vp = 4.99 L (RSE 10.5%)
    lvmax <- log(0.53); label("Maximum rate of the additional nonlinear clearance of free drug Vmax (nM/h)") # Table 5: Vmax = 0.53 nM/h (RSE 7.96%)
    lkm <- log(4.31); label("Free-drug concentration at half-maximal nonlinear clearance Km (nM)") # Table 5: Km = 4.31 nM (RSE 9.94%)

    # SC absorption (Table 5)
    lfdepot <- log(0.705); label("SC bioavailability F (fraction)") # Table 5: F = 0.705 (RSE 4.44%)
    lka <- log(0.00742); label("First-order SC absorption rate constant ka (1/h)") # Table 5: Ka = 0.00742 1/h (RSE 7.29%)
    ltlag <- log(2); label("SC absorption lag time (h)") # Table 5: Lag time = 2 h (RSE 2.25%)

    # Target (soluble TFPI) turnover and QSS binding (Table 5)
    lrbase_target <- log(3.75); label("Baseline total TFPI concentration (nM)") # Table 5: Baseline total TFPI = 3.75 nM (RSE 1.52%)
    lkdeg <- log(0.0161); label("Degradation rate constant of free soluble TFPI kdeg (1/h)") # Table 5: kdeg = 0.0161 1/h (RSE 5.72%)
    lkint <- log(0.000754); label("Elimination rate constant of the marstacimab-TFPI complex kcom (1/h)") # Table 5: kcom = 0.000754 1/h (RSE 12.3%)
    lkss <- log(72.2); label("Quasi-steady-state constant KSS (nM)") # Table 5: KSS = 72.2 (RSE 5.04%); nM, the concentration unit of the TMDD system

    # Body-weight exponents (Table 5); Q shares the CL exponent and Vp the Vc exponent (see vignette)
    e_wt_cl_q <- 1.14; label("Weight exponent on CL and Q (unitless)") # Table 5: Weight exponent on CL = 1.14 (RSE 32%)
    e_wt_vc_vp <- 1.86; label("Weight exponent on Vc and Vp (unitless)") # Table 5: Weight exponent on Vc = 1.86 (RSE 18.5%)

    # Categorical effects on CL, linear form CL * (1 + theta * indicator) (Supplementary eq 7.1)
    e_dis_healthy_cl <- 0.199; label("Fractional change in CL for healthy volunteers vs hemophilia (unitless)") # Table 5: 'Additional effect on CL for healthy population' = 0.199
    e_race_asian_cl <- 0.346; label("Fractional change in CL for Asian vs non-Asian participants (unitless)") # Table 5: 'Additional effect on CL for Asian participants' = 0.346
    e_ada_pos_cl <- 0.229; label("Fractional change in CL for ADA-positive participants (unitless)") # Table 5: 'Additional effect on CL for ADA +ve participants' = 0.229
    e_renalimp_mild_cl <- -0.168; label("Fractional change in CL for mild renal impairment (unitless)") # Table 5: 'Additional effect on CL for mild renal impairment' = -0.168
    e_dis_hemophilia_b_cl <- 0.118; label("Fractional change in CL for hemophilia B vs hemophilia A (unitless)") # Table 5: 'Additional effect on CL for HB participants' = 0.118

    # Peak thrombin sigmoid Emax model driven by free TFPI (Table 8; Supplementary eq 5)
    le0 <- log(37.8); label("Peak thrombin at high free TFPI, PBASE (nM)") # Table 8: PBASE = 37.8 nM (RSE 11.8%)
    lrmax_peak_thrombin <- log(91.4); label("Peak thrombin as free TFPI approaches zero, PMAX (nM)") # Table 8: PMAX = 91.4 nM (RSE 13.3%)
    lec50 <- log(2.92); label("Free TFPI concentration at half-maximal effect EC50 (nM)") # Table 8: EC50 = 2.92 nM (RSE 12.3%)
    hill <- fixed(5); label("Hill coefficient (unitless)") # Table 8: Hill coefficient (fixed) = 5

    # IIV (Table 5 and Table 8 report variances; the printed %CV equals sqrt(variance) * 100)
    etalcl + etalvc ~ c(0.262, -0.1299287, 0.232) # Table 5: IIV CL 0.262, IIV Vc 0.232, correlation -0.527; covariance = -0.527 * sqrt(0.262 * 0.232)
    etalq ~ fixed(0.0225) # Table 5: IIV on Q 0.0225 (15% CV, not estimated)
    etalvp ~ fixed(0.0225) # Table 5: IIV on Vp 0.0225 (15% CV, not estimated)
    etalkm ~ fixed(0.0225) # Table 5: IIV on Km 0.0225 (15% CV, not estimated)
    etalvmax ~ 0.101 # Table 5: IIV on Vmax 0.101 (31.8% CV)
    etalfdepot ~ 0.0332 # Table 5: IIV on F 0.0332 (18.2% CV)
    etalka ~ 0.278 # Table 5: IIV on Ka 0.278 (52.7% CV)
    etaltlag ~ fixed(0.0225) # Table 5: IIV on lag time 0.0225 (15% CV, not estimated)
    etalrbase_target ~ 0.0249 # Table 5: IIV on baseline total TFPI 0.0249 (15.8% CV)
    etalkdeg ~ 0.023 # Table 5: IIV on Kdeg 0.023 (15.2% CV)
    etalkint ~ fixed(0.0225) # Table 5: IIV on Kcom 0.0225 (15% CV, not estimated)
    etalkss ~ 0.105 # Table 5: IIV on KSS constant 0.105 (32.3% CV)
    etale0 ~ 0.104 # Table 8: IIV on PBASE 0.104 (32.2% CV)
    etalec50 ~ 0.0283 # Table 8: IIV on EC50 0.0283 (16.8% CV)
    etalrmax_peak_thrombin ~ 0.111 # Table 8: IIV on PMAX 0.111 (33.3% CV)

    # Residual error (Table 5: additive on the log scale; Table 8: additive on the natural scale)
    expSd <- 0.392; label("Log-scale residual SD of total marstacimab (unitless)") # Table 5: additive residual error on marstacimab (log scale) = 0.392
    expSd_total_target <- 0.158; label("Log-scale residual SD of total TFPI (unitless)") # Table 5: additive residual error on total TFPI (log scale) = 0.158
    addSd_peak_thrombin <- 26.3; label("Additive residual SD of peak thrombin (nM)") # Table 8: additive error peak TGA = 26.3 nM
  })

  model({
    # Marstacimab molecular weight converts mg doses to nmol and nM to ug/mL.
    # The paper does not state it; 146 kDa is the Hympavzi prescribing information.
    mw_marstacimab <- 146000 # g/mol
    nmol_per_mg <- 1e6 / mw_marstacimab

    wt_ratio <- WT / 71.1

    cov_cl <- (1 + e_dis_healthy_cl * DIS_HEALTHY) *
      (1 + e_race_asian_cl * RACE_ASIAN) *
      (1 + e_ada_pos_cl * ADA_POS) *
      (1 + e_renalimp_mild_cl * RENALIMP_MILD) *
      (1 + e_dis_hemophilia_b_cl * DIS_HEMOPHILIA_B)

    cl <- exp(lcl + etalcl) * wt_ratio^e_wt_cl_q * cov_cl
    vc <- exp(lvc + etalvc) * wt_ratio^e_wt_vc_vp
    q <- exp(lq + etalq) * wt_ratio^e_wt_cl_q
    vp <- exp(lvp + etalvp) * wt_ratio^e_wt_vc_vp
    vmax <- exp(lvmax + etalvmax)
    km <- exp(lkm + etalkm)
    fdepot <- exp(lfdepot + etalfdepot)
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag + etaltlag)

    rbase_target <- exp(lrbase_target + etalrbase_target)
    kdeg <- exp(lkdeg + etalkdeg)
    kint <- exp(lkint + etalkint)
    kss <- exp(lkss + etalkss)
    ksyn <- kdeg * rbase_target

    e0 <- exp(le0 + etale0)
    rmax_peak_thrombin <- exp(lrmax_peak_thrombin + etalrmax_peak_thrombin)
    ec50 <- exp(lec50 + etalec50)

    # Quasi-steady-state free drug (nM) from total drug and total target
    ctot <- central / vc
    qss_term <- ctot - total_target - kss
    cfree <- 0.5 * (qss_term + sqrt(qss_term^2 + 4 * kss * ctot))
    complex <- ctot - cfree
    free_target <- total_target - complex

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot -
      cl * cfree -
      q * cfree + q * peripheral1 / vp -
      kint * complex * vc -
      vmax * cfree / (km + cfree) * vc
    d/dt(peripheral1) <- q * cfree - q * peripheral1 / vp
    d/dt(total_target) <- ksyn - kdeg * free_target - kint * complex

    total_target(0) <- rbase_target

    f(depot) <- fdepot * nmol_per_mg
    alag(depot) <- tlag
    f(central) <- nmol_per_mg

    # Total marstacimab in ug/mL (nM * g/mol / 1e6)
    Cc <- ctot * mw_marstacimab / 1e6

    peak_thrombin <- e0 + (rmax_peak_thrombin - e0) / (1 + (free_target / ec50)^hill)

    Cc ~ lnorm(expSd)
    total_target ~ lnorm(expSd_total_target)
    peak_thrombin ~ add(addSd_peak_thrombin)
  })
}
