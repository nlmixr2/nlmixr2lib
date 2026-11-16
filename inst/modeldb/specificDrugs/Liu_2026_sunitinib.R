Liu_2026_sunitinib <- function() {
  description <- paste(
    "Joint population PK model for oral sunitinib and its active metabolite SU12662 in adults with",
    "metastatic renal cell carcinoma on 50 mg once daily 4-weeks-on/2-weeks-off (Liu 2026; 294",
    "patients pooled from four Pfizer trials). Each compound is two-compartment with its own",
    "first-order input depot: the sunitinib dose enters depot with bioavailability F (fixed to 1,",
    "+14.8% in women, IIV) and, in parallel, depot_su12662 with F * fm (fm fixed to 0.21, +24.3% in",
    "women, IIV), so the SU12662 clearance and volumes are apparent relative to F * fm. Apparent",
    "clearance and central and peripheral volumes decline hyperbolically with time since the first",
    "dose, CL0 / (1 + beta * week), with additive IIV on each beta shared between parent and",
    "metabolite through one random effect and a metabolite scaling factor. Baseline CL/F falls",
    "with age (exponential, centred 60 years), with higher albumin (power, 41 g/L) and in women;",
    "albumin and body weight (82.4 kg) scale Vc/F. Correlated parent-metabolite IIV on CL/F, Vc/F",
    "and Vp/F; proportional residual error on each analyte. Dose both depots with the full",
    "sunitinib dose; f() applies F and F * fm. Concentrations in ng/mL."
  )
  reference <- paste(
    "Liu H, Friberg LE. Time-Dependent Bias and Prognostic Confounding in Exposure-Response Analysis of",
    "Targeted Therapies: Lessons from Sunitinib in Metastatic Renal Cancer.",
    "Clin Pharmacol Ther. 2026. doi:10.1002/cpt.70447. PMID: 42587360; PMCID: PMC13468949.",
    "Final PK model: Supplementary Material S2 (time-varying CL/F and V/F equations and Table S2)."
  )
  vignette <- "Liu_2026_sunitinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")
  dosing <- c("depot", "depot_su12662")

  covariateData <- list(
    SEXF = list(
      description = "Sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Table S2 footnote 2: SEX = 0 for males and 1 for females, applied as a multiplicative",
        "factor (1 + theta_SEX * SEX) on F_sunitinib (+0.148), fm (+0.243), sunitinib CL/F",
        "(-0.0985) and SU12662 CL/F (-0.0864). The source column SEX already has the canonical",
        "SEXF orientation; no transformation."
      ),
      source_name = "SEX"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table S2 footnote 3: exponential effect exp((AGE - 60) * theta_AGE) on baseline CL/F,",
        "centred on the median age of 60 years; theta_AGE = -0.00590 (sunitinib) and -0.00793",
        "(SU12662) per year."
      ),
      source_name = "AGE"
    ),
    ALB = list(
      description = "Serum albumin at baseline",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table S2 footnote 4: power function (ALB / 41)^theta_ALB centred on the median albumin of",
        "41 g/L. One exponent shared by both analytes on CL/F (-0.507) and one on Vc/F (-0.741):",
        "Table S2 prints the same estimate and RSE in the parent and metabolite columns."
      ),
      source_name = "ALB"
    ),
    WT = list(
      description = "Body weight at baseline",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table S2 footnote 5: power function (WT / 82.4)^theta_WT centred on the median weight of",
        "82.4 kg, on sunitinib Vc/F (0.220), SU12662 CL/F (0.217) and SU12662 Vc/F (0.447). No",
        "weight effect on sunitinib CL/F."
      ),
      source_name = "WT"
    )
  )

  covariatesDataExcluded <- list(
    HGB = list(
      description = "Haemoglobin at baseline",
      units = "g/L",
      type = "continuous",
      notes = "Among the baseline laboratory covariates screened by stepwise covariate modelling (Methods; Supplementary Figure S1); not retained in the final PK model (Table S2)."
    ),
    ECOG_GE1 = list(
      description = "ECOG performance status >= 1",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Methods 'Baseline covariates'); not retained in the final PK model (Table S2)."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "sunitinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sunitinib (total)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "sunitinib (total)", units = "mg", specimen = "plasma", verified = TRUE),
    depot_su12662 = list(
      analyte = "SU12662 input (fraction F * fm of the sunitinib dose)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central_su12662 = list(analyte = "SU12662 (total)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_su12662 = list(analyte = "SU12662 (total)", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 294L,
    n_studies = 4L,
    age_range = "24-87 years (exposure-response subset, n = 165); PK-model median 60 years",
    age_mean = "57.1 years (SD 10.3; exposure-response subset)",
    weight_range = "45-150 kg (exposure-response subset); PK-model median 82.4 kg",
    sex_female_pct = 35,
    race_ethnicity = c(White = 92, NonWhite = 8),
    disease_state = "Metastatic renal cell carcinoma: cytokine-refractory (studies 1006 and 014) or treatment-naive (studies 1034 and 1065)",
    dose_range = "Sunitinib 50 mg orally once daily, 4 weeks on / 2 weeks off; dose reductions to 37.5 or 25 mg",
    regions = "Multinational (Pfizer trials NCT00077974, NCT00054886, NCT00083889, NCT00267748)",
    notes = paste(
      "PK analysis: 1688 sunitinib and 1715 SU12662 total plasma concentrations from 294 patients",
      "(105 from study 1006, 60 from 014, 42 from 1034, 88 from 1065; Table 1, Supplementary S2),",
      "trough samples in the on-treatment weeks and at the end of the 2-week break, up to 5 cycles.",
      "Demographics (age, weight, sex, race) are tabulated only for the 165-patient",
      "exposure-response subset (Supplementary Table S1); the PK-model centring values (age 60",
      "years, albumin 41 g/L, weight 82.4 kg) are the medians of the 294-patient PK dataset."
    )
  )

  ini({
    # Absorption: fixed to literature values (Methods; Table S2).
    lka <- fixed(log(0.224))
    label("Sunitinib first-order absorption rate constant ka (1/h)") # Table S2 Ka parent 0.224 (fixed)
    lka_su12662 <- fixed(log(0.132))
    label("SU12662 first-order input rate constant (1/h)") # Table S2 Ka metabolite 0.132 (fixed)

    # Sunitinib disposition at t = 0 (Table S2 'Parent' column).
    lcl <- log(35.3)
    label("Sunitinib baseline apparent clearance CL0/F (L/h)") # Table S2 CL0/F parent 35.3 (RSE 3.2%)
    lvc <- log(2250)
    label("Sunitinib baseline apparent central volume Vc0/F (L)") # Table S2 Vc0/F parent 2250 (RSE 2.9%)
    lq <- log(1.72)
    label("Sunitinib apparent intercompartmental clearance Q/F (L/h)") # Table S2 Q/F parent 1.72 (RSE 10.9%)
    lvp <- log(272)
    label("Sunitinib baseline apparent peripheral volume Vp0/F (L)") # Table S2 Vp0/F parent 272 (RSE 6.1%)

    # SU12662 disposition at t = 0 (Table S2 'Metabolite' column).
    lcl_su12662 <- log(18.1)
    label("SU12662 baseline apparent clearance CL0/(F*fm) (L/h)") # Table S2 CL0/F metabolite 18.1 (RSE 2.7%)
    lvc_su12662 <- log(2210)
    label("SU12662 baseline apparent central volume Vc0/(F*fm) (L)") # Table S2 Vc0/F metabolite 2210 (RSE 4.2%)
    lq_su12662 <- log(1.02)
    label("SU12662 apparent intercompartmental clearance Q/(F*fm) (L/h)") # Table S2 Q/F metabolite 1.02 (RSE 14.7%)
    lvp_su12662 <- log(384)
    label("SU12662 baseline apparent peripheral volume Vp0/(F*fm) (L)") # Table S2 Vp0/F metabolite 384 (RSE 15.8%)

    # Bioavailability and fraction metabolised (Methods: F_sunitinib fixed to 1,
    # fm fixed to 0.21, F_SU12662 = F_sunitinib * fm).
    lfdepot <- fixed(log(1))
    label("Sunitinib relative bioavailability F (unitless)") # Table S2 F_sunitinib 1 (fixed)
    lfm_su12662 <- fixed(log(0.21))
    label("Fraction of the sunitinib dose converted to SU12662 fm (unitless)") # Table S2 fm 0.21 (fixed)

    # Sex effects, multiplicative (1 + theta * SEXF) (Table S2 footnote 2).
    e_sexf_fdepot <- 0.148
    label("Fractional change in F for women (unitless)") # Table S2 SEX on F_sunitinib 0.148 (RSE 10.9%)
    e_sexf_fm_su12662 <- 0.243
    label("Fractional change in fm for women (unitless)") # Table S2 SEX on fm 0.243 (RSE 18.6%)
    e_sexf_cl <- -0.0985
    label("Fractional change in sunitinib CL/F for women (unitless)") # Table S2 SEX on CL/F parent -0.0985 (RSE 12.7%)
    e_sexf_cl_su12662 <- -0.0864
    label("Fractional change in SU12662 CL/F for women (unitless)") # Table S2 SEX on CL/F metabolite -0.0864 (RSE 17.4%)

    # Age, exp((AGE - 60) * theta) on CL/F (Table S2 footnote 3).
    e_age_cl <- -0.00590
    label("Exponential age coefficient on sunitinib CL/F (1/year)") # Table S2 AGE on CL/F parent -0.00590 (RSE 15.1%)
    e_age_cl_su12662 <- -0.00793
    label("Exponential age coefficient on SU12662 CL/F (1/year)") # Table S2 AGE on CL/F metabolite -0.00793 (RSE 14.2%)

    # Albumin, (ALB / 41)^theta, one exponent shared by both analytes
    # (Table S2 footnote 4).
    e_alb_cl <- -0.507
    label("Power exponent of ALB/41 on sunitinib and SU12662 CL/F (unitless)") # Table S2 ALB on CL/F -0.507 (RSE 19.0%), both columns
    e_alb_vc <- -0.741
    label("Power exponent of ALB/41 on sunitinib and SU12662 Vc/F (unitless)") # Table S2 ALB on Vc/F -0.741 (RSE 12.9%), both columns

    # Body weight, (WT / 82.4)^theta (Table S2 footnote 5).
    e_wt_cl_su12662 <- 0.217
    label("Power exponent of WT/82.4 on SU12662 CL/F (unitless)") # Table S2 WT on CL/F metabolite 0.217 (RSE 14.9%)
    e_wt_vc <- 0.220
    label("Power exponent of WT/82.4 on sunitinib Vc/F (unitless)") # Table S2 WT on Vc/F parent 0.220 (RSE 19.9%)
    e_wt_vc_su12662 <- 0.447
    label("Power exponent of WT/82.4 on SU12662 Vc/F (unitless)") # Table S2 WT on Vc/F metabolite 0.447 (RSE 14.9%)

    # Time-varying CL/F and V/F: parameter(t) = baseline / (1 + beta * week)
    # (Supplementary S2 equations; Table S2 footnote 6). beta_V acts on both
    # the central and the peripheral volume.
    cl_time_hyp <- 0.00365
    label("Hyperbolic time-decline rate of sunitinib CL/F, beta_CL (1/week)") # Table S2 beta_TIME,CL parent 0.00365 (RSE 13.3%)
    cl_time_hyp_su12662 <- 0.00801
    label("Hyperbolic time-decline rate of SU12662 CL/F, beta_CL (1/week)") # Table S2 beta_TIME,CL metabolite 0.00801 (RSE 11.6%)
    v_time_hyp <- 0.00168
    label("Hyperbolic time-decline rate of sunitinib Vc/F and Vp/F, beta_V (1/week)") # Table S2 beta_TIME,Vc parent 0.00168 (RSE 40.5%)
    v_time_hyp_su12662 <- 0.0125
    label("Hyperbolic time-decline rate of SU12662 Vc/F and Vp/F, beta_V (1/week)") # Table S2 beta_TIME,Vc metabolite 0.0125 (RSE 20.7%)

    # Divisors of the shared beta random effects for SU12662:
    # beta_SU12662,i = beta_SU12662 + eta_i / scale (Table S2 footnote 8).
    scale_etacl_time_hyp_su12662 <- 0.90
    label("Divisor of the shared beta_CL random effect for SU12662 (unitless)") # Table S2 IIV beta_TIME,CL/F scale 0.90 (RSE 24.0%)
    scale_etav_time_hyp_su12662 <- 1.38
    label("Divisor of the shared beta_V random effect for SU12662 (unitless)") # Table S2 IIV beta_TIME,Vc/F scale 1.38 (RSE 25.6%)

    # Exponential IIV. Table S2 footnote 7 reports each as sqrt(variance) x 100%,
    # so omega^2 = (value / 100)^2. Correlations: CL/F parent-metabolite 0.70;
    # V/F parent-metabolite 0.68, applied to the Vc/F pair and the Vp/F pair
    # (Supplementary S2 text; see the vignette assumptions).
    etalfdepot ~ 0.070225 # Table S2 IIV F_sunitinib 26.5 -> 0.265^2
    etalfm_su12662 ~ 0.051076 # Table S2 IIV fm 22.6 -> 0.226^2
    etalcl + etalcl_su12662 ~ c(0.00982081, 0.0104749, 0.022801) # Table S2 IIV CL/F 9.91 and 15.1, correlation 0.70
    etalvc + etalvc_su12662 ~ c(0.006889, 0.00660348, 0.013689) # Table S2 IIV Vc/F 8.30 and 11.7, correlation 0.68
    etalvp + etalvp_su12662 ~ c(0.068121, 0.105425, 0.352836) # Table S2 IIV Vp/F 26.1 and 59.4, correlation 0.68
    # Additive IIV on beta (Table S2 'Additive IIV and scales', 0.8 and 0.9
    # reported as sqrt(variance) x 100, in 1/week).
    etacl_time_hyp ~ 0.000064 # Table S2 IIV beta_TIME,CL/F 0.8 -> 0.008^2
    etav_time_hyp ~ 0.000081 # Table S2 IIV beta_TIME,Vc/F 0.9 -> 0.009^2

    # Proportional residual error (Table S2 RUV, reported as %CV).
    propSd <- 0.217
    label("Sunitinib proportional residual error (fraction)") # Table S2 RUV parent 21.7% (RSE 1.7%)
    propSd_su12662 <- 0.219
    label("SU12662 proportional residual error (fraction)") # Table S2 RUV metabolite 21.9% (RSE 2.1%)
  })

  model({
    # Time since the first dose in weeks (Supplementary S2: 'T is time in
    # weeks'); the first dose must be at t = 0.
    tweek <- t / 168

    # Individual time-decline rates; one random effect per parameter class is
    # shared by both analytes (Supplementary S2 equations).
    beta_cl <- cl_time_hyp + etacl_time_hyp
    beta_cl_su12662 <- cl_time_hyp_su12662 + etacl_time_hyp / scale_etacl_time_hyp_su12662
    beta_v <- v_time_hyp + etav_time_hyp
    beta_v_su12662 <- v_time_hyp_su12662 + etav_time_hyp / scale_etav_time_hyp_su12662

    # Bioavailability and fraction metabolised (Table S2 footnotes 1-2).
    fdepot <- exp(lfdepot + etalfdepot) * (1 + e_sexf_fdepot * SEXF)
    fm_su12662 <- exp(lfm_su12662 + etalfm_su12662) * (1 + e_sexf_fm_su12662 * SEXF)

    # Sunitinib.
    ka <- exp(lka)
    cl0 <- exp(lcl + etalcl) * (1 + e_sexf_cl * SEXF) * exp(e_age_cl * (AGE - 60)) *
      (ALB / 41)^e_alb_cl
    vc0 <- exp(lvc + etalvc) * (ALB / 41)^e_alb_vc * (WT / 82.4)^e_wt_vc
    vp0 <- exp(lvp + etalvp)
    q <- exp(lq)
    cl <- cl0 / (1 + beta_cl * tweek)
    vc <- vc0 / (1 + beta_v * tweek)
    vp <- vp0 / (1 + beta_v * tweek)

    # SU12662.
    ka_su12662 <- exp(lka_su12662)
    cl0_su12662 <- exp(lcl_su12662 + etalcl_su12662) * (1 + e_sexf_cl_su12662 * SEXF) *
      exp(e_age_cl_su12662 * (AGE - 60)) * (ALB / 41)^e_alb_cl * (WT / 82.4)^e_wt_cl_su12662
    vc0_su12662 <- exp(lvc_su12662 + etalvc_su12662) * (ALB / 41)^e_alb_vc *
      (WT / 82.4)^e_wt_vc_su12662
    vp0_su12662 <- exp(lvp_su12662 + etalvp_su12662)
    q_su12662 <- exp(lq_su12662)
    cl_su12662 <- cl0_su12662 / (1 + beta_cl_su12662 * tweek)
    vc_su12662 <- vc0_su12662 / (1 + beta_v_su12662 * tweek)
    vp_su12662 <- vp0_su12662 / (1 + beta_v_su12662 * tweek)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    kel_su12662 <- cl_su12662 / vc_su12662
    k12_su12662 <- q_su12662 / vc_su12662
    k21_su12662 <- q_su12662 / vp_su12662

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(depot_su12662) <- -ka_su12662 * depot_su12662
    d/dt(central_su12662) <- ka_su12662 * depot_su12662 - kel_su12662 * central_su12662 -
      k12_su12662 * central_su12662 + k21_su12662 * peripheral1_su12662
    d/dt(peripheral1_su12662) <- k12_su12662 * central_su12662 - k21_su12662 * peripheral1_su12662

    # Both depots receive the full sunitinib dose (Methods:
    # F_SU12662 = F_sunitinib * fm).
    f(depot) <- fdepot
    f(depot_su12662) <- fdepot * fm_su12662

    # mg/L -> ng/mL
    Cc <- 1000 * central / vc
    Cc_su12662 <- 1000 * central_su12662 / vc_su12662

    Cc ~ prop(propSd)
    Cc_su12662 ~ prop(propSd_su12662)
  })
}
