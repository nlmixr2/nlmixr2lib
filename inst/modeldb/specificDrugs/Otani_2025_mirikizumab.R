Otani_2025_mirikizumab <- function() {
  description <- "Two-compartment population PK model with first-order SC absorption and fixed allometric weight scaling for mirikizumab (anti-IL-23p19 mAb) in children and adolescents (2 to <18 years) with moderately to severely active ulcerative colitis (Otani 2025, phase II SHINE-1)"
  reference <- "Otani Y, Chua L, Komocsar WJ, Larkin A, Johns J, Zhang XC. Mirikizumab pharmacokinetics and exposure-response in pediatric patients with moderate-to-severe ulcerative colitis. CPT Pharmacometrics Syst Pharmacol. 2025;14(3):474-485. doi:10.1002/psp4.13286"
  vignette <- "Otani_2025_mirikizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(analyte = "mirikizumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "mirikizumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "mirikizumab", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Fixed allometric exponents 0.8 on CL and Q and 1.0 on Vc and Vp, normalised to 70.7 kg (the LUCENT 1+2 adult median; Table 2 footnotes d and e). Table 2 labels the effect 'baseline weight'; the control stream (Appendix S1) reads it from a column named BWTV.",
      source_name = "BWTV"
    ),
    ALB = list(
      description = "Serum albumin (time-varying)",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear effect on CL, CL multiplier = 1 + (-0.0218) * (ALB - 45) (Table 2 footnote d; control stream CLALBU). Results, Covariate effects: only time-varying albumin was retained.",
      source_name = "ALBU"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 26L,
    n_studies = 1L,
    n_observations = 226L,
    age_range = "2 to <18 years (3 aged 2-7, 4 aged 8-11, 19 aged 12-17)",
    age_mean = "11.8 years (SD 3.4)",
    weight_range = "about 16-80 kg (2 participants <= 20 kg, 13 > 20 to <= 40 kg, 11 > 40 kg)",
    weight_mean = "40.5 kg (SD 16.0)",
    bmi_mean = "18.3 kg/m^2 (SD 4.3)",
    sex_female_pct = 57.7,
    race_ethnicity = c(White = 68.0, Asian = 28.0, Other = 4.0),
    disease_state = "Moderately to severely active ulcerative colitis (modified Mayo score 4-9, endoscopic subscore >= 2) with intolerance, inadequate response or loss of response to corticosteroids, immunomodulators, biologics or JAK inhibitors.",
    dose_range = "Induction IV Q4W x 3 (infused over >= 40 min): 300 mg (> 40 kg), 5 or 10 mg/kg (> 20 to <= 40 kg), 5 mg/kg (<= 20 kg); extended IV induction for week-12 non-responders. Maintenance SC Q4W to week 48: 200 mg (> 40 kg), 100 mg (> 20 to <= 40 kg), 50 mg (<= 20 kg).",
    regions = "37 sites in 5 countries.",
    notes = "Phase II open-label SHINE-1 (NCT04004611); Otani 2025 Table 1. 232 serum samples; 6 BLQ troughs from patients who had stopped dosing were excluded (226 observations in the final fit, Appendix S1). Ka and Q fixed to the adult LUCENT 1+2 estimates (Chua 2023). Fit in NONMEM 7.5.2, FOCE-I."
  )

  ini({
    # Typical values at the reference subject (70.7 kg, albumin 45 g/L).
    # Rates are published per hour; x 24 converts to per day.
    lka <- fixed(log(0.00706 * 24)); label("SC absorption rate constant (1/day)") # Appendix S1 THETA2 KA 0.00706 1/h FIX; Table 2 footnote c: fixed to the LUCENT 1+2 estimate
    lcl <- log(0.0209 * 24); label("Clearance at 70.7 kg and albumin 45 g/L (L/day)") # Appendix S1 THETA3 CL 0.0209 L/h; Table 2: 0.021 L/h
    lvc <- log(3.48); label("Central volume of distribution at 70.7 kg (L)") # Appendix S1 THETA4 V2 3.48 L; Table 2
    lq <- fixed(log(0.0087 * 24)); label("Intercompartmental clearance at 70.7 kg (L/day)") # Appendix S1 THETA5 Q 0.0087 L/h FIX; Table 2 footnote c
    lvp <- log(1.40); label("Peripheral volume of distribution at 70.7 kg (L)") # Appendix S1 THETA6 V3 1.4 L; Table 2: 1.40 L
    logitfdepot <- log(0.498 / (1 - 0.498)); label("Logit of SC bioavailability (unitless)") # Appendix S1 THETA1 F1 0.498 (TVPHI = LOG(TVBIO/(1-TVBIO))); Table 2: 49.8%

    e_wt_cl_q <- fixed(0.8); label("Allometric exponent of weight on CL and Q (unitless)") # Appendix S1 THETA9 CLBWTV 0.8 FIX; Table 2
    e_wt_vc_vp <- fixed(1.0); label("Allometric exponent of weight on Vc and Vp (unitless)") # Appendix S1 THETA10 VBWTV 1 FIX; Table 2
    e_alb_cl <- -0.0218; label("Linear effect of albumin on CL (per g/L)") # Appendix S1 THETA11 CLALBU -0.0218; Table 2 footnote d

    # IIV variances from Appendix S1 OMEGA; Table 2 footnote b reports
    # 100 * sqrt(exp(OMEGA) - 1) for both.
    etalcl ~ 0.128 # Appendix S1 OMEGA3 CL 0.128; Table 2: 37.1%
    etalogitfdepot ~ 0.471 # Appendix S1 OMEGA1 F1 0.471, additive on the logit scale (PHI = TVPHI + ETA(1)); Table 2: 77.6%

    propSd <- 0.288; label("Proportional residual error (fraction)") # Appendix S1 THETA7 PROP 0.288 with SIGMA fixed to 1; additive THETA8 fixed to 0
  })
  model({
    wt_cl_q <- (WT / 70.7)^e_wt_cl_q
    wt_vc_vp <- (WT / 70.7)^e_wt_vc_vp
    alb_cl <- 1 + e_alb_cl * (ALB - 45)

    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * alb_cl * wt_cl_q
    vc <- exp(lvc) * wt_vc_vp
    q <- exp(lq) * wt_cl_q
    vp <- exp(lvp) * wt_vc_vp

    fdepot <- expit(logitfdepot + etalogitfdepot)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    # mg / L = ug/mL
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
