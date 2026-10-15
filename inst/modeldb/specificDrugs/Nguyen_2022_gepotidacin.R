Nguyen_2022_gepotidacin <- function() {
  description <- "Three-compartment intravenous population PK model for gepotidacin in healthy adults and adults with acute bacterial skin and skin structure infections, with allometric body-weight scaling of every clearance and volume and an estimated (modelled) infusion duration; extended to paediatric dose selection for plague by a fixed paracetamol-derived postmenstrual-age maturation function on clearance applied below 2 years of age. The paper's Simcyp whole-body PBPK model is not included."
  reference <- "Nguyen D, Shaik JS, Tai G, Tiffany C, Perry C, Dumont E, Gardiner D, Barth A, Singh R, Hossain M. Comparison between physiologically based pharmacokinetic and population pharmacokinetic modelling to select paediatric doses of gepotidacin in plague. Br J Clin Pharmacol. 2022;88(2):416-428. doi:10.1111/bcp.14996"
  vignette <- "Nguyen_2022_gepotidacin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline body weight. Power scaling of every clearance and volume with reference weight 70 kg: exponent 0.75 on CL and 1 on V1, Q2, V2, Q3 and V3 (Supporting Information Table S2 footnote equations).",
      source_name = "WT"
    ),
    PNA = list(
      description = "Postnatal age",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = "Used only for the paediatric maturation function on clearance, which applies when PNA < 24 months (2 years); at PNA >= 24 months the maturation factor is 1. Postmenstrual age is derived inside model() as PNA converted to weeks plus 40 weeks, assuming full-term birth (paper Section 2.2.2). Adult users must still supply the column (any value >= 24 months, e.g. AGE * 12).",
      source_name = "PNA"
    ),
    TINF = list(
      description = "Nominal duration of the intravenous infusion of the dose record",
      units = "h",
      type = "continuous",
      reference_category = NULL,
      notes = "Selects the modelled infusion duration D1: the estimated D1 for a nominal 1-h infusion (1.06 h) when TINF = 1, the estimated D1 for a nominal 2-h infusion (2.14 h) when TINF = 2, and the nominal TINF itself for any other duration (the paper estimated D1 only for the 1-h and 2-h infusions it fitted). Dose records must carry rate = -2 so that rxode2 uses the modelled dur(central).",
      source_name = "D1 (infusion duration category)"
    )
  )

  compartmentData <- list(
    central = list(analyte = "gepotidacin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "gepotidacin", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "gepotidacin", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 249,
    n_studies = 5,
    age_range = "adults; study mean ages 26.6-44.7 years (Table 1); paediatric simulation cohort 0.01-19.6 years",
    weight_range = "not reported for the fitted dataset; adult simulation cohort 49.5-117.9 kg; paediatric simulation cohort 2.5-79 kg",
    sex_female_pct = NA_real_,
    disease_state = "Healthy adult volunteers (four phase I studies, including Japanese subjects) and adults with acute bacterial skin and skin structure infections (phase II study BTZ116704); no PK difference between the healthy and patient populations.",
    dose_range = "Single IV 1-h infusions of 200-1800 mg and 2-h infusions of 1800 mg; repeat 2-h infusions of 400-1000 mg twice daily and 1000 mg three times daily; patients 750 or 1000 mg twice daily or 1000 mg three times daily for up to 10 days",
    regions = "not reported",
    notes = "Final base model pooled four phase I healthy-volunteer studies (N = 140: BTZ115198, BTZ115775, BTZ116666 and BTZ115774) and the phase II ABSSSI study BTZ116704 (N = 109) (Section 2.2, Table 1). Fitted in NONMEM 7.2. The paediatric maturation function was not estimated from gepotidacin data: its constants were taken from a published paracetamol model and applied for simulation only."
  )

  ini({
    # Supporting Information Table S2 'Gepotidacin PK Parameter Estimates of
    # the Final Covariate Population PK Model in Adults' (typical values for a
    # 70 kg adult).
    lcl <- log(38.0); label("Clearance, 70 kg adult (L/h)") # Table S2: CL 38.0 L/h (RSE 1.42%)
    lvc <- log(20.1); label("Central volume of distribution V1, 70 kg adult (L)") # Table S2: V1 20.1 L (RSE 10.5%)
    lq <- log(3.85); label("Intercompartmental clearance central-peripheral1 Q2, 70 kg adult (L/h)") # Table S2: Q2 3.85 L/h (RSE 6.86%)
    lvp <- log(67.1); label("Peripheral1 volume of distribution V2, 70 kg adult (L)") # Table S2: V2 67.1 L (RSE 5.25%)
    lq2 <- log(48.5); label("Intercompartmental clearance central-peripheral2 Q3, 70 kg adult (L/h)") # Table S2: Q3 48.5 L/h (RSE 4.21%)
    lvp2 <- log(70.3); label("Peripheral2 volume of distribution V3, 70 kg adult (L)") # Table S2: V3 70.3 L (RSE 1.88%)
    ld1_inf1h <- log(1.06); label("Modelled infusion duration D1 for a nominal 1-h infusion (h)") # Table S2: D1 1-h infusion 1.06 h (RSE 7.17%)
    ld1_inf2h <- log(2.14); label("Modelled infusion duration D1 for a nominal 2-h infusion (h)") # Table S2: D1 2-h infusion 2.14 h (RSE 0.79%)

    # Allometric exponents (Table S2 footnote: CL = 38.0*(WT/70)^0.75;
    # V1, Q2, V2, Q3 and V3 each scale linearly with WT/70). No uncertainty is
    # reported, so they are fixed.
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL (unitless)") # Table S2 footnote: (WT/70)^0.75
    e_wt_vc_q_vp_q2_vp2 <- fixed(1); label("Allometric exponent of body weight on V1, Q2, V2, Q3 and V3 (unitless)") # Table S2 footnote: (WT/70)

    # Paediatric maturation of clearance (Section 2.2.2 equation), constants
    # fixed from the paracetamol population PK model of the paper's ref 22.
    pma_tm50 <- fixed(52.2); label("Postmenstrual age at 50% of mature clearance TM50 (weeks)") # Section 2.2.2: TM50 52.2 weeks, fixed
    pma_hill <- fixed(3.43); label("Hill coefficient of the postmenstrual-age maturation of clearance (unitless)") # Section 2.2.2: Hill coefficient 3.43, fixed

    # IIV reported as CV%; converted with omega^2 = log(1 + CV^2).
    etalcl ~ 0.0302 # Table S2: IIV CL 17.5 CV% (shrinkage 6.7%)
    etalvc ~ 0.3727 # Table S2: IIV V1 67.2 CV% (shrinkage 18%)
    etalq ~ 0.1792 # Table S2: IIV Q2 44.3 CV% (shrinkage 28.6%)
    etalvp ~ 0.1554 # Table S2: IIV V2 41.0 CV% (shrinkage 28.6%)
    etalq2 ~ 0.0578 # Table S2: IIV Q3 24.4 CV% (shrinkage 35.4%)
    etalvp2 ~ 0.0098 # Table S2: IIV V3 9.92 CV% (shrinkage 47%)

    propSd <- 0.222; label("Proportional residual error (fraction)") # Table S2: residual error 22.2 CV% of proportional error
  })

  model({
    wt_ratio <- WT / 70

    # Postmenstrual age (weeks) from postnatal age (months), full-term birth
    # assumed: 30.4375 days per month / 7 days per week, plus 40 weeks.
    pma_wk <- PNA * 30.4375 / 7 + 40
    fmat <- pma_wk^pma_hill / (pma_tm50^pma_hill + pma_wk^pma_hill)
    # The maturation function was applied only below 2 years of age.
    if (PNA >= 24) {
      fmat <- 1
    }

    cl <- exp(lcl + etalcl) * wt_ratio^e_wt_cl * fmat
    vc <- exp(lvc + etalvc) * wt_ratio^e_wt_vc_q_vp_q2_vp2
    q <- exp(lq + etalq) * wt_ratio^e_wt_vc_q_vp_q2_vp2
    vp <- exp(lvp + etalvp) * wt_ratio^e_wt_vc_q_vp_q2_vp2
    q2 <- exp(lq2 + etalq2) * wt_ratio^e_wt_vc_q_vp_q2_vp2
    vp2 <- exp(lvp2 + etalvp2) * wt_ratio^e_wt_vc_q_vp_q2_vp2

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # Modelled infusion duration (NONMEM D1 with RATE = -2); dose records
    # must carry rate = -2. Durations other than the fitted 1 h and 2 h use
    # the nominal TINF.
    d1 <- TINF
    if (TINF == 1) {
      d1 <- exp(ld1_inf1h)
    }
    if (TINF == 2) {
      d1 <- exp(ld1_inf2h)
    }
    dur(central) <- d1

    # mg / L = ug/mL
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
