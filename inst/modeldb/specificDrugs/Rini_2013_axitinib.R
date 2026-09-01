Rini_2013_axitinib <- function() {
  description <- "Two-compartment population PK model for axitinib pooled across healthy volunteers and patients with metastatic renal cell carcinoma or other solid tumours (Rini 2013). First-order absorption with an estimated lag time; linear-proportional effects of age > 60 years, Japanese ethnicity and active smoking on systemic clearance; power-form effect of body weight on the central volume of distribution (reference 74.1 kg); a linear-proportional fasting effect on the absorption rate constant ka; a linear-proportional fasting effect on bioavailability F that applies only to crystal polymorph Form IV; and a linear-proportional reduction in F for the marketed crystal polymorph Form XLI relative to Form IV. Pooled data from 590 subjects (383 healthy volunteers, 181 metastatic RCC patients and 26 patients with other solid tumours) across 17 trials."
  reference   <- "Rini BI, Garrett M, Poland B, Dutcher JP, Rixe O, Wilding G, Stadler WM, Pithavala YK, Kim S, Tarazi J, Motzer RJ. Axitinib in metastatic renal cell carcinoma: results of a pharmacokinetic and pharmacodynamic analysis. J Clin Pharmacol. 2013;53(5):491-504. doi:10.1002/jcph.73"
  vignette    <- "Rini_2013_axitinib"
  units       <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Rini 2013 Methods 'Model development' names the states
  # explicitly (central and peripheral compartments of a linear 2-compartment
  # model with first-order absorption from an oral depot).
  compartmentData <- list(
    depot       = list(analyte = "axitinib", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "axitinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "axitinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight at screening.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on Vc centred on 74.1 kg; the exponent 0.778 was freely estimated rather than fixed at 1. Rini 2013 Results 'PK model' prints the typical-value equation Vc = 47.3 L * (weight / 74.1 kg)^0.778; the reported screening median body weight was 74 kg (range 37-136). Body weight was tested on CL and showed no relationship.",
      source_name        = "weight"
    ),
    AGE_GT60 = list(
      description        = "Indicator for subject age strictly greater than 60 years at screening.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (age <= 60 years).",
      notes              = "Rini 2013 Results 'Effects of covariates': age was first tested as a continuous covariate on CL (Supplementary Figure S2-A) but entered the final model as a binary threshold to avoid overestimating CL at younger ages given the non-uniform age distribution across the 17 pooled studies; several thresholds were tested and 60 years was the most significant (Supplementary Figure S2-B). Time-fixed per subject; derive as as.integer(AGE > 60) from the continuous age column.",
      source_name        = "Age>60"
    ),
    RACE_JAPANESE = list(
      description        = "Japanese-heritage indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Japanese: Caucasian, Black, Hispanic, other Asian, and 'Other/Not listed').",
      notes              = "Rini 2013 Results 'Effects of covariates' tested Asian and Japanese separately and in combination; only the Japanese indicator was significant on CL, so the non-Japanese reference deliberately pools other Asian subjects with Caucasian, Black, Hispanic and 'Other/Not listed'. The Discussion cautions that the Japanese effect may be confounded by the lower body weight and higher age of the Japanese subjects in the pooled dataset.",
      source_name        = "RaceJapanese"
    ),
    SMOKE = list(
      description        = "Active-smoker indicator at screening.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ex-smoker or non-smoker; Rini 2013 pools these two into a single reference group).",
      notes              = "Rini 2013 Results 'Effects of covariates': active smokers had 102% greater CL than ex-smokers or non-smokers, but only 19 of 590 subjects (3%) were active smokers and the estimate is imprecise (44% RSE, 95% CI 0.144 to 1.90), so the clinical significance is unclear and the Discussion calls for confirmation in a larger population. Set SMOKE = 0 to simulate the reference (ex-smoker / non-smoker) population. Two-level (current vs not-current) encoding, not the paired SMOKE_CURRENT / SMOKE_NEVER 3-level encoding, because the paper pools ex-smokers and never-smokers into one reference group.",
      source_name        = "Smokeractive"
    ),
    FED = list(
      description        = "Fed-vs-fasted indicator at the dose record.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (fed; the typical-value reference for ka and F in this model).",
      notes              = "The fasting effect is applied internally as (1 - FED), so FED = 1 (fed) leaves ka and F at their typical-value fed Form IV estimates and FED = 0 (fasted) activates the linear-proportional fasting increases (Rini 2013 Table 2 rows 'ka (hour-1): fed', 'Fasting effect on ka' and 'Fasting on F, Form IV'). The fasting effect on ka applies regardless of formulation; the fasting effect on F applies only to Form IV (Table 2 footnote g: 'there was no observed food effect with Form XLI'). Per-dose-record covariate; the crossover food-effect studies mean a single subject may carry both values.",
      source_name        = "FED"
    ),
    FORM_AXI_XLI = list(
      description        = "Crystal polymorph form indicator for axitinib (Form XLI marketed vs Form IV earlier).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Form IV; the typical-value F reference in Rini 2013 Table 2, F = 0.457 in the fed state).",
      notes              = "1 = Form XLI (marketed commercial crystal polymorph); 0 = Form IV (earlier Phase I crystal polymorph and the typical-value F reference). Linear-proportional effect on F only; no effect on ka, CL, Vc, Q or Vp was retained. Rini 2013 Table 2 footnote g gives F(Form XLI) = 0.457 * (1 - 0.121) = 0.402 and states that no food effect was observed with Form XLI, so the fasting effect on F is gated to Form IV records. Per-dose-record covariate.",
      source_name        = "Form"
    )
  )

  covariatesDataExcluded <- list(
    BSA = list(
      description = "Body surface area.",
      units       = "m^2",
      type        = "continuous",
      notes       = "Tested on CL in the stepwise covariate search but not retained in the final model (Rini 2013 Methods 'Model development' and Results 'Effects of covariates')."
    ),
    ALT = list(
      description = "Alanine aminotransferase.",
      units       = "U/L",
      type        = "continuous",
      notes       = "Tested on CL but not retained; median (range) at screening 21 U/L (5-188). The Discussion attributes the null result to the requirement for normal hepatic function at study entry (Rini 2013 Results 'Effects of covariates')."
    ),
    AST = list(
      description = "Aspartate aminotransferase.",
      units       = "U/L",
      type        = "continuous",
      notes       = "Tested on CL but not retained; median (range) at screening 22 U/L (9-154) (Rini 2013 Results 'Effects of covariates')."
    ),
    BILI = list(
      description = "Total bilirubin.",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Tested on CL but not retained; median (range) at screening 0.7 mg/dL (0.1-3) (Rini 2013 Results 'Effects of covariates')."
    ),
    CRCL = list(
      description = "Creatinine clearance estimated with the Cockcroft-Gault equation.",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Tested on CL but not retained; median (range) at screening 103 mL/min (8-214) (Rini 2013 Results 'Effects of covariates')."
    ),
    SEXF = list(
      description = "Female-sex indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = "The only categorical covariate tested on BOTH CL and Vc; not retained on either (Rini 2013 Methods 'Model development'). The pooled population was 85% male."
    ),
    ECOG_GE1 = list(
      description = "Baseline ECOG performance status of 1 or worse.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested on CL but not retained in the population PK model (Rini 2013 Methods 'Model development'). Baseline ECOG PS is separately a significant prognostic factor for progression-free and overall survival in the PK/PD Cox analysis (Table 5), which is not encoded here."
    ),
    DIS_HEALTHY = list(
      description = "Healthy-volunteer vs cancer-patient study-population indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested on CL but not retained: Rini 2013 Results 'Effects of covariates' reports no difference in axitinib CL between healthy volunteers and cancer patients, with clearance etas centred at zero in both groups (Supplementary Figure S3), which is what justifies pooling the two populations into a single model."
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 590L,
    n_studies       = 17L,
    age_range       = "18-85 years",
    age_median      = "42 years overall; 60 years (32-85) in patients and 32 years (18-69) in healthy volunteers",
    weight_range    = "37-136 kg",
    weight_median   = "74 kg",
    sex_female_pct  = 15,
    race_ethnicity  = c(Caucasian = 61, Other = 39),
    disease_state   = "Pooled healthy volunteers (n = 383) and patients with advanced solid tumours (n = 207, of whom 181 had metastatic renal cell carcinoma, cytokine- or sorafenib-refractory).",
    dose_range      = "Healthy volunteers received a single 5 mg oral dose (one study additionally gave 1 mg intravenously to estimate absolute bioavailability, and one Chinese study gave single 5, 7 and 10 mg oral doses). Patients received 5 mg orally twice daily with titration permitted up to 10 mg twice daily.",
    regions         = "United States, France, Germany, Japan, Belgium, Singapore and China.",
    renal_function  = "Creatinine clearance (Cockcroft-Gault) median 103 mL/min (range 8-214); normal renal function was a study-entry requirement.",
    hepatic_function = "AST median 22 U/L (9-154), ALT median 21 U/L (5-188), bilirubin median 0.7 mg/dL (0.1-3); normal hepatic function was a study-entry requirement. Sixteen subjects with Child-Pugh hepatic impairment (study 13) were excluded from the analysis dataset.",
    smoking_status  = "19 of 590 subjects (3%) were active smokers.",
    notes           = "Demographics from Rini 2013 Results 'Subject characteristics'; study designs from Table 1. Sixteen patients from dose-finding study 1 who started above the 5 mg bid maximum tolerated dose were also excluded. The PK/PD (exposure-efficacy) analyses used a 168-patient metastatic RCC subset described in Table 3; those Cox and logistic regression models are not encoded in this file (see the vignette Errata)."
  )

  ini({
    # Structural PK parameters -- Rini 2013 Table 2 final-model estimates.
    # Reference subject: a 74.1 kg, <= 60-year-old, non-Japanese, non-smoking
    # subject receiving 5 mg axitinib orally as crystal polymorph Form IV in
    # the fed state. Absolute (not apparent) CL and volumes: F was anchored by
    # the 1 mg intravenous arm of study 8 (Rini 2013 Methods 'Model development').
    lka     <- log(0.482) ; label("Absorption rate constant ka, fed reference (1/h)")                # Rini 2013 Table 2: ka (hour-1): fed = 0.482
    lcl     <- log(14.6)  ; label("Systemic clearance CL (L/h)")                                     # Rini 2013 Table 2: CL = 14.6 L/h
    lvc     <- log(47.3)  ; label("Central volume of distribution Vc at WT = 74.1 kg (L)")           # Rini 2013 Table 2: Vc = 47.3 L
    lq      <- log(4.00)  ; label("Inter-compartmental clearance Q (L/h)")                           # Rini 2013 Table 2: Q = 4.00 L/h
    lvp     <- log(393)   ; label("Peripheral volume of distribution Vp (L)")                        # Rini 2013 Table 2: Vp = 393 L
    lfdepot <- log(0.457) ; label("Absolute oral bioavailability F, Form IV fed reference (unitless)")  # Rini 2013 Table 2: F: fed/Form IV = 0.457
    ltlag   <- log(0.454) ; label("Absorption lag time (h)")                                         # Rini 2013 Table 2: tlag = 0.454 h

    # Covariate effects -- Rini 2013 Table 2 final-model estimates. The three
    # CL effects and the two F effects are linear-proportional shifts in the
    # paper's own notation theta1 * (1 + theta2 * indicator) (Methods
    # 'Model development'); weight on Vc is a power function centred on the
    # population median.
    e_age_cl   <- -0.213 ; label("Linear-proportional effect of age > 60 years on CL")                          # Rini 2013 Table 2 'Age >60-yr effect on CL' = -0.213 and the Results 'PK model' equation CL = 14.6 * (1 - 0.213 * Age>60) * ...
    e_jpn_cl   <- -0.249 ; label("Linear-proportional effect of Japanese ethnicity on CL")                      # Rini 2013 Table 2 'Japanese ethnicity effect on CL' = -0.249 and the Results 'PK model' equation ... * (1 - 0.249 * RaceJapanese) * ...
    e_smoke_cl <-  1.02  ; label("Linear-proportional effect of active smoking on CL")                          # Rini 2013 Table 2 'Smoking status on CL' = 1.02 and the Results 'PK model' equation ... * (1 + 1.02 * Smokeractive)
    e_wt_vc    <-  0.778 ; label("Power exponent of (WT / 74.1 kg) on Vc (unitless)")                           # Rini 2013 Table 2 'Weight effect on Vc' = 0.778 and the Results 'PK model' equation Vc = 47.3 * (weight / 74.1 kg)^0.778
    e_fast_ka  <-  1.97  ; label("Linear-proportional fasting effect on ka (relative to the fed reference)")    # Rini 2013 Table 2 'Fasting effect on ka' = 1.97; footnote e gives ka_fasted = 0.482 * (1 + 1.97) = 1.43 /h
    e_fast_f   <-  0.330 ; label("Linear-proportional fasting effect on F, Form IV only")                       # Rini 2013 Table 2 'Fasting on F, Form IV' = 0.330; footnote f gives F_fasted,IV = 0.457 * (1 + 0.33) = 0.608
    e_xli_f    <- -0.121 ; label("Linear-proportional Form XLI effect on F (relative to the Form IV reference)") # Rini 2013 Table 2 'Form XLI on F' = -0.121; footnote g gives F_XLI = 0.457 * (1 - 0.121) = 0.402

    # IIV -- exponential on each structural PK parameter (Rini 2013 Methods
    # 'Model development'). Table 2 footnote a reports the interindividual
    # variability as %CV without stating the back-transformation. The
    # convention is pinned by the companion healthy-volunteer analysis from
    # the same modelling group (Garrett 2014, BJCP 77:480-492), whose Methods
    # print the definition verbatim -- "%IIV = sqrt(omega^2) x 100" -- and
    # whose Table 3 satisfies it numerically: omega^2(Vc) = 0.0949 gives
    # sqrt(0.0949) = 30.8%, exactly the 30.8% Vc IIV quoted in that paper's
    # Discussion. The exact-log-normal reading, sqrt(exp(omega^2) - 1), would
    # give 31.6% and does not match. So %CV is sqrt(omega^2) * 100 here and
    # omega^2 = (%CV / 100)^2. No correlations between CL, Vc and ka were
    # reported. Table 2 footnote d states that Q and Vp were modelled as 100%
    # correlated with the same interindividual variability, which is encoded
    # exactly as the single shared random effect etalq_lvp (a 2x2 block with
    # unit correlation would be singular).
    etalcl    ~ 0.358801 ; label("IIV variance on log CL (59.9% CV)")            # Rini 2013 Table 2: CL 14.6 (59.9) -> 0.599^2
    etalvc    ~ 0.157609 ; label("IIV variance on log Vc (39.7% CV)")            # Rini 2013 Table 2: Vc 47.3 (39.7) -> 0.397^2
    etalq_lvp ~ 0.753424 ; label("IIV variance shared by log Q and log Vp (86.8% CV)")  # Rini 2013 Table 2: Q 4.00 (86.8) with footnote d (Q and Vp 100% correlated, same IIV) -> 0.868^2
    etalka    ~ 0.592900 ; label("IIV variance on log ka (77% CV)")              # Rini 2013 Table 2: ka 0.482 (77) -> 0.77^2

    # Residual error -- Rini 2013 Methods 'Model development' modelled
    # log-transformed concentrations with a proportional residual, i.e. an
    # additive residual on the log scale, which is a log-normal residual on
    # the linear concentration scale. Table 2 reports the residual as a
    # percentage separately for the oral (58.2%) and intravenous (33.5%)
    # routes; only the oral residual is encoded here because a single rxode2
    # endpoint carries a single error model and the library use of this model
    # is oral dosing. See the vignette Errata for the intravenous value.
    expSd <- 0.582 ; label("Log-normal residual error SD (oral; log-scale fraction)")  # Rini 2013 Table 2: Residual error (%), oral = 58.2
  })

  model({
    # Reference covariate value for the Vc power function (Rini 2013 Results
    # 'PK model' equation; the screening median body weight was 74 kg).
    ref_wt <- 74.1

    # Individual PK parameters. The three CL covariates and the Vc weight
    # term reproduce the two typical-value equations printed in Rini 2013
    # Results 'PK model':
    #   CL = 14.6 L/h * (1 - 0.213 * Age>60) * (1 - 0.249 * RaceJapanese) * (1 + 1.02 * Smokeractive)
    #   Vc = 47.3 L   * (weight / 74.1 kg)^0.778
    ka <- exp(lka + etalka) * (1 + e_fast_ka * (1 - FED))
    cl <- exp(lcl + etalcl) *
      (1 + e_age_cl * AGE_GT60) *
      (1 + e_jpn_cl * RACE_JAPANESE) *
      (1 + e_smoke_cl * SMOKE)
    vc <- exp(lvc + etalvc) * (WT / ref_wt)^e_wt_vc
    q  <- exp(lq  + etalq_lvp)
    vp <- exp(lvp + etalq_lvp)

    # Absolute oral bioavailability. The fasting effect is gated to Form IV
    # records because Rini 2013 Table 2 footnote g reports no observed food
    # effect with Form XLI, and the Form XLI effect is estimated in the fed
    # state (Results 'Effects of covariates': "When comparing the 2
    # formulations in the fed state, F for Form XLI was 12.1% lower than
    # Form IV"). F = 1 is anchored by the 1 mg intravenous arm of study 8.
    fdepot <- exp(lfdepot) *
      (1 + e_fast_f * (1 - FED) * (1 - FORM_AXI_XLI)) *
      (1 + e_xli_f * FORM_AXI_XLI)

    tlag <- exp(ltlag)

    # Micro-constants for the explicit two-compartment ODE system
    # (NONMEM ADVAN4/TRANS4 equivalent).
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    alag(depot) <- tlag
    f(depot)    <- fdepot

    # Concentrations are expressed in ng/mL: the axitinib dose is entered in
    # mg so central / vc is mg/L, which is scaled by 1000 to ng/mL.
    Cc <- central / vc * 1000
    Cc ~ lnorm(expSd)
  })
}
