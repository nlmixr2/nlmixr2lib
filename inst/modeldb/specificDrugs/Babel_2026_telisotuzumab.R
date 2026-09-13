Babel_2026_telisotuzumab <- function() {
  description <- paste0(
    "Population pharmacokinetic model of the telisotuzumab vedotin ",
    "(Teliso-V) CONJUGATE (the c-Met-directed antibody-drug conjugate ",
    "measured by the total conjugate assay) in adults with c-Met ",
    "protein overexpressing advanced solid tumours, most of them ",
    "non-squamous non-small cell lung cancer (Babel 2026, n = 304 ",
    "pooled from a phase 1 dose-ranging study and the LUMINOSITY ",
    "phase 2 study). Linear two-compartment disposition after ",
    "intravenous infusion, with interindividual variability on ",
    "clearance and central volume and a combined proportional plus ",
    "additive residual error. Body weight enters as a shared power ",
    "exponent on CL and Q and a second shared exponent on Vc and Vp; ",
    "treatment-emergent anti-drug antibody status, baseline albumin ",
    "and race act on CL; age, baseline albumin and sex act on Vc. ",
    "Non-linear clearance was tested and rejected. The companion ",
    "unconjugated MMAE payload model, which Babel 2026 developed ",
    "INDEPENDENTLY of this one, is Babel_2026_telisotuzumab_mmae; the ",
    "four companion exposure-response logistic models are the ",
    "Babel_2026_telisotuzumab_orr / _neuropathy / _corneal / _teae ",
    "family."
  )
  reference <- paste(
    "Babel H, Brunsdon P, Engelhardt B, Schmitt V, Ratajczak C, Mensing S,",
    "Menon RM, Parikh A. Population pharmacokinetics and exposure-response",
    "analyses for telisotuzumab vedotin in patients with c-Met protein",
    "overexpressing tumors.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(1):e70219.",
    "doi:10.1002/psp4.70219. PMCID PMC12945708.",
    "All parameter values are from Data S1 (Supporting Information) Table S7,",
    "'Final Model Parameter Estimates and Variability of Teliso-V Conjugate",
    "and Unconjugated MMAE Payload Pharmacokinetics', Teliso-V Conjugate block.",
    sep = " "
  )
  vignette <- "Babel_2026_telisotuzumab"

  units <- list(
    time          = "day",
    dosing        = "mg",
    concentration = "ug/mL"
  )

  # CL and Q are in L/day and Vc and Vp in L (Table S7), so an amount in
  # mg divided by a volume in L gives mg/L == ug/mL, which is the scale
  # of the CavgADC axis in Babel 2026 Figures 2 and 4 (0 to 10 ug/mL).
  compartmentData <- list(
    central     = list(analyte = "telisotuzumab vedotin conjugate", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "telisotuzumab vedotin conjugate", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Normalised to the overall-population median of 68.9 kg inside model(). Babel 2026 Methods: 'Continuous covariates were normalized to the median of the overall population and incorporated into the model using a power function'. The median is read from Table S4 (All Participants, N = 304, Body Weight median 68.9 kg) and is corroborated by Table S3, whose forest-plot reference group is 'Body Weight <= 68.9 kg'. Range 36.0-144 kg. One exponent is shared by CL and Q and a second by Vc and Vp, exactly as Table S7 tabulates them ('Body Weight on CL and Q', 'Body Weight on Vc and Vp').",
      source_name        = "Body Weight"
    ),
    AGE = list(
      description        = "Age at baseline.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Normalised to the overall-population median of 65 years inside model() (Babel 2026 Table S4, All Participants median 65, range 30-87). Acts on Vc only in the conjugate model. Babel 2026 Table S3 describes the forest-plot reference stratum as 'Age 40-65 years', which is a binning of the same continuous covariate, not a different centring.",
      source_name        = "Age"
    ),
    ALB = list(
      description        = "Baseline serum albumin concentration.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Normalised to the overall-population median of 41.6 g/L inside model() (Babel 2026 Table S4, All Participants median 41.6 g/L, range 29.0-52.0; corroborated by the Table S3 reference group 'Baseline Albumin <= 41.6 g/L'). Enters as a power function on both CL and Vc, with negative exponents in each case, so a higher albumin lowers clearance and volume and therefore raises exposure - the direction shown in Babel 2026 Figure 1A, where albumin > 41.6 g/L versus <= 41.6 g/L gives an AUCtau ratio of 1.18.",
      source_name        = "Baseline Albumin"
    ),
    SEXF = list(
      description        = "Sex indicator; 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male; Babel 2026 Table S3 reference group 'Sex: Male')",
      notes              = "Multiplicative factor on Vc only. Babel 2026 Table S4: 115 of 304 (38%) female. The tabulated factor 0.925 is the female-versus-male value: a smaller central volume in women raises Cmax, matching the Figure 1A 'Female : Male' Cmax ratio of 1.07 while leaving AUCtau essentially unchanged at 1.01 (Vc does not enter AUCtau).",
      source_name        = "Sex"
    ),
    RACE_BLACK = list(
      description        = "Black or African American race indicator; 1 = Black or African American, 0 otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White; Babel 2026 Table S3 reference group 'Race: White')",
      notes              = "Multiplicative factor on CL. Babel 2026 Table S4: 8 of 304 (3%) Black or African American versus 198 of 304 (65%) White. The authors caution that the small n limits conclusions: the model-predicted AUCtau increase of about 56% relative to White patients (Figure 1A) sits inside the observed spread of White patients, and no dose adjustment by race is recommended.",
      source_name        = "Race: Black or African American"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator; 1 = Asian, 0 otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White; Babel 2026 Table S3 reference group 'Race: White')",
      notes              = "Multiplicative factor on CL. Babel 2026 Table S4: 98 of 304 (32%) Asian. The factor 0.887 implies an AUCtau ratio of 1/0.887 = 1.13 relative to White patients, matching the 1.10 shown in Figure 1A. RACE_BLACK and RACE_ASIAN are mutually exclusive; a White patient carries 0 for both.",
      source_name        = "Race: Asian"
    ),
    ADA_POS = list(
      description        = "Treatment-emergent anti-drug antibody status; 1 = positive, 0 = negative.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA negative; Babel 2026 Table S3 reference group 'Treatment-emergent ADA: Negative')",
      notes              = "Multiplicative factor on CL. Babel 2026 Table S2 footnote c defines treatment-emergent ADA as either (1) baseline ADA negative or missing with at least one post-baseline ADA positive, or (2) baseline ADA positive with at least one post-baseline positive result at least 3.3 x baseline in LUMINOSITY or 2 x baseline in NCT02099058, per each study's assay. Nominally time-varying (a patient becomes ADA positive at seroconversion); Babel 2026 does not state whether the covariate was implemented as a time-varying flag or as an ever-positive subject-level flag, so it is carried here as a subject-level indicator. The factor 1.17 implies an AUCtau ratio of 1/1.17 = 0.855 for ADA-positive patients, matching the 0.846 in Figure 1A. Figure 1A reports n = 61 positive versus n = 243 negative.",
      source_name        = "Treatment-emergent ADA status"
    )
  )

  covariatesDataExcluded <- list(
    RENALIMP = list(
      description = "Baseline renal-function category derived from Cockcroft-Gault creatinine clearance (normal, mild, moderate/severe).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened against conjugate CL in the Babel 2026 Table S2 covariate sweep but NOT retained in the final conjugate model (Table S7 lists no renal term for the conjugate). Renal function WAS retained on unconjugated MMAE payload clearance; see Babel_2026_telisotuzumab_mmae, which carries RENALIMP_MILD, RENALIMP_MOD and RENALIMP_SEV."
    ),
    HEPIMP_MILD = list(
      description = "Mild hepatic impairment indicator by National Cancer Institute criteria.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on conjugate CL and Vc (Babel 2026 Table S2) but not retained; the Discussion states mild hepatic impairment 'was not identified as a significant covariate on Teliso-V conjugate or free MMAE payload clearance'. No patient with moderate or severe hepatic impairment was enrolled, so those strata are unidentified."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 304L,
    n_studies      = 2L,
    age_range      = "median 65 years, range 30-87 (Babel 2026 Table S4, All Participants)",
    weight_range   = "median 68.9 kg, range 36.0-144 (Babel 2026 Table S4, All Participants)",
    sex_female_pct = 37.8,
    race_ethnicity = c(White = 65, Asian = 32, `Black or African American` = 3),
    disease_state  = "Advanced solid tumours likely to express c-Met (phase 1, NCT02099058, n = 35, all tumour types, monotherapy cohorts receiving Process II drug material only) and locally advanced or metastatic c-Met protein overexpressing non-small cell lung cancer (LUMINOSITY phase 2, NCT03539536, n = 269).",
    dose_range     = "Phase 1: 0.15-3.3 mg/kg every 3 weeks and 1.6-2.2 mg/kg every 2 weeks. LUMINOSITY: 1.6 or 1.9 mg/kg every 2 weeks. Approved regimen 1.9 mg/kg Q2W as an intravenous infusion, capped at 190 mg for patients weighing at least 100 kg.",
    regions        = "Europe 26%, North America 26%, Asia 28%, rest of world 20% (Babel 2026 Table S4)",
    notes          = paste0(
      "Baseline albumin median 41.6 g/L (range 29.0-52.0); baseline ",
      "renal function normal 36%, mild 43%, moderate 19%, severe 1%, ",
      "missing less than 1%; baseline hepatic function normal 88% with ",
      "the remainder mild and no moderate or severe (Babel 2026 Table ",
      "S4). Japanese 9%; Hispanic or Latino 3%. 6.03% of conjugate ",
      "records were below the limit of quantitation and a further ",
      "2.84% were excluded as outliers by a time-binned mixture-model ",
      "rule described in Data S1 Methods; observations with absolute ",
      "conditional weighted residual above 5 for the base model were ",
      "also excluded. Structural conjugate parameters were estimated ",
      "with relative standard errors of at most 23.1%."
    )
  )

  ini({
    # ==================================================================
    # Structural disposition. All four values are the population
    # estimates in Babel 2026 Data S1 Table S7, Teliso-V Conjugate
    # block. NONMEM 7.5.1; the table reports the estimates on the
    # natural (untransformed) scale together with %RSE and a 95% CI, so
    # each is wrapped in log() here to put it on the model's log scale.
    # ==================================================================
    lcl <- log(1.34); label("Clearance (L/day)")                        # Table S7 'CL (L/day)' 1.34, %RSE 2.98, 95% CI (1.26, 1.42)
    lvc <- log(3.44); label("Central volume of distribution (L)")       # Table S7 'Vc (L)' 3.44, %RSE 1.38, 95% CI (3.35, 3.54)
    lq  <- log(1.21); label("Intercompartmental clearance (L/day)")     # Table S7 'Q (L/day)' 1.21, %RSE 4.59, 95% CI (1.11, 1.33)
    lvp <- log(2.41); label("Peripheral volume of distribution (L)")    # Table S7 'Vp (L)' 2.41, %RSE 2.68, 95% CI (2.28, 2.54)

    # ==================================================================
    # Covariate effects. Continuous covariates are power functions of
    # the covariate divided by its overall-population median, and
    # categorical covariates are multiplicative factors relative to the
    # reference group (Babel 2026 Methods, Population Pharmacokinetic
    # Analyses). Table S7 prints one row per effect; the row label is
    # quoted in each trailing comment.
    # ==================================================================
    e_wt_cl_q   <- 0.405;  label("Power exponent on (WT/68.9 kg) shared by CL and Q (unitless)")     # Table S7 'Body Weight on CL and Q' 0.405, %RSE 18.9, 95% CI (0.255, 0.556)
    e_wt_vc_vp  <- 0.590;  label("Power exponent on (WT/68.9 kg) shared by Vc and Vp (unitless)")    # Table S7 'Body Weight on Vc and Vp' 0.590, %RSE 6.96, 95% CI (0.510, 0.671)
    e_ada_cl    <- 1.17;   label("Multiplicative factor on CL for treatment-emergent ADA positive versus negative (unitless)") # Table S7 'ADA on CL' 1.17, %RSE 3.64, 95% CI (1.09, 1.26)
    e_alb_cl    <- -0.970; label("Power exponent on (ALB/41.6 g/L) for CL (unitless)")               # Table S7 'Albumin on CL' -0.970, %RSE 18.7, 95% CI (-1.33, -0.614)
    e_black_cl  <- 0.763;  label("Multiplicative factor on CL for Black or African American versus White (unitless)")          # Table S7 'Black or African American vs. White on CL' 0.763, %RSE 7.17, 95% CI (0.663, 0.878)
    e_asian_cl  <- 0.887;  label("Multiplicative factor on CL for Asian versus White (unitless)")    # Table S7 'Asian vs. White on CL' 0.887, %RSE 3.48, 95% CI (0.828, 0.949)
    e_age_vc    <- 0.223;  label("Power exponent on (AGE/65 years) for Vc (unitless)")               # Table S7 'Age on Vc' 0.223, %RSE 21.7, 95% CI (0.128, 0.318)
    e_alb_vc    <- -0.464; label("Power exponent on (ALB/41.6 g/L) for Vc (unitless)")               # Table S7 'Albumin on Vc' -0.464, %RSE 23.1, 95% CI (-0.674, -0.254)
    e_sexf_vc   <- 0.925;  label("Multiplicative factor on Vc for female versus male (unitless)")    # Table S7 'Sex on Vc' 0.925, %RSE 2.32, 95% CI (0.884, 0.968)

    # ==================================================================
    # Interindividual variability. Table S7's 'Population Estimate'
    # column for the IIV rows is the VARIANCE of the log-normal random
    # effect, which the table self-pins: its footnote gives
    # '%CV was calculated as SQRT(exp(omega2)-1)*100', and
    # sqrt(exp(0.0946) - 1) * 100 = 31.5% and
    # sqrt(exp(0.0268) - 1) * 100 = 16.5%, reproducing the printed %CV
    # column exactly. So the numbers below are omega-squared, which is
    # what an nlmixr2 eta line takes.
    # ==================================================================
    etalcl ~ 0.0946  # Table S7 'IIV on CL' variance 0.0946, 31.5 %CV, 4.03% shrinkage
    etalvc ~ 0.0268  # Table S7 'IIV on Vc' variance 0.0268, 16.5 %CV, 13.0% shrinkage

    # ==================================================================
    # Residual unexplained variability. Babel 2026 Results: 'Residual
    # variability was best described by combined additive and
    # proportional error terms'. Table S7 reports both as VARIANCES, so
    # each is entered as its square root because nlmixr2 error terms
    # take standard deviations. The additive row is printed as
    # 'Additional Error (Variance)' in the conjugate block and as
    # 'Additive Error (Variance)' in the payload block; the two are the
    # same quantity and the conjugate wording is a typographical slip.
    # ==================================================================
    propSd <- sqrt(0.0528); label("Proportional residual error on Cc (fraction)")   # Table S7 'Proportional Error (Variance)' 0.0528, %RSE 1.13, 95% CI (0.0516, 0.0540)
    addSd  <- sqrt(0.0447); label("Additive residual error on Cc (ug/mL)")          # Table S7 'Additional Error (Variance)' 0.0447, %RSE 3.40, 95% CI (0.0417, 0.0476)
  })

  model({
    # ----- Individual parameters -----
    # Continuous covariates enter as power functions of covariate /
    # population median; categorical covariates as multiplicative
    # factors relative to the reference group. RACE_BLACK and
    # RACE_ASIAN are mutually exclusive indicators, so at most one of
    # the two race factors is active for any subject and a White
    # patient carries neither.
    cl <- exp(lcl + etalcl) *
      (WT / 68.9)^e_wt_cl_q *
      (ALB / 41.6)^e_alb_cl *
      e_ada_cl^ADA_POS *
      e_black_cl^RACE_BLACK *
      e_asian_cl^RACE_ASIAN
    vc <- exp(lvc + etalvc) *
      (WT / 68.9)^e_wt_vc_vp *
      (AGE / 65)^e_age_vc *
      (ALB / 41.6)^e_alb_vc *
      e_sexf_vc^SEXF
    q  <- exp(lq) * (WT / 68.9)^e_wt_cl_q
    vp <- exp(lvp) * (WT / 68.9)^e_wt_vc_vp

    # ----- Micro-constants -----
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ----- Two-compartment linear disposition (Babel 2026 Figure S1A) -----
    # Dosing is an intravenous infusion into central; there is no
    # absorption compartment and no non-linear elimination pathway
    # (non-linear CL was tested and did not improve the objective
    # function value).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ----- Observation -----
    # Cc is the Teliso-V conjugate serum concentration in ug/mL, the
    # analyte of the total conjugate assay.
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
