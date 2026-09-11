Lee_2024_gentamicin_teigen <- function() {
  description <- paste(
    "Teigen one-compartment population PK model for intravenous gentamicin in adults with end-stage",
    "renal disease receiving thrice-weekly intermittent hemodialysis, as transcribed and used by",
    "Lee 2024 as the a priori prior for MAP Bayesian (NONMEM POSTHOC) estimation of individual PK",
    "parameters in a single obese hemodialysis patient. Elimination is linear and switches between",
    "two clearance regimes: a non-hemodialysis clearance CL_NHD that scales linearly with",
    "Cockcroft-Gault creatinine clearance computed on ideal body weight, and a total on-dialysis",
    "clearance CL_HD that REPLACES CL_NHD while a session is running (gated by the time-varying",
    "RRT_HEMODIAL_ACTIVE regressor). Volume of distribution carries no covariate. All population",
    "parameters are FIXED priors; Lee 2024 estimated only the individual random effects.",
    sep = " "
  )
  reference <- paste(
    "Lee H, Yoon S, Chung JY. Individual pharmacokinetic parameter estimation of gentamicin in an",
    "obese hemodialysis patient using non-linear mixed effect model. Transl Clin Pharmacol.",
    "2024;32(3):150-158. doi:10.12793/tcp.2024.32.e14 (Methods, 'Individual PK parameter estimation",
    "using NONMEM and Monolix'; Results; Table 2; Figs. 1-2).",
    "The structural model and every population parameter value originate from Teigen MM, Duffull S,",
    "Dang L, Johnson DW. Dosing of gentamicin in patients with end-stage renal disease receiving",
    "hemodialysis. J Clin Pharmacol. 2006;46(11):1259-1267. doi:10.1177/0091270006292987, which is",
    "closed access and was NOT available on disk; the values encoded here are those Lee 2024 prints.",
    sep = " "
  )
  vignette <- "Lee_2024_gentamicin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "gentamicin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance, Cockcroft-Gault computed on IDEAL body weight, raw (NOT body-surface-area normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Lee 2024 Methods: 'The CrCL was calculated using the Cockcroft and Gault formula, which used",
        "the ideal body weight.' The source model equation is written with CrCL in L/h --",
        "'CLNHD ... 0.453 x CrCL/0.53 L/h' -- so the reference value 0.53 L/h is 8.833 mL/min and the",
        "conversion is applied inside model(). This is the raw un-normalized variant of the CRCL",
        "canonical (same convention as Delattre_2010_amikacin.R, Dohmann_2025_piperacillin.R and",
        "Takada_2025_vancomycin.R); supplying a BSA-normalized value here would silently rescale the",
        "renal term. Teigen development-cohort values quoted by Lee 2024: group MEAN 0.53 L/h",
        "(8.83 mL/min, the reference) and group MAXIMUM 1.24 L/h (20.7 mL/min). Lee 2024's patient",
        "had a mean CrCL of 2.44 L/h (40.7 mL/min; range 1.87-3.24 L/h = 31.2-54.0 mL/min), which the",
        "Discussion argues is an OVERestimate of true renal function because prolonged hospitalization",
        "and chronic illness had reduced her muscle mass and hence her serum creatinine. Enters as the",
        "plain linear ratio (CRCL / 8.833), i.e. an exponent of exactly 1 -- Lee 2024 prints no power",
        "term. Because the term is unbounded and multiplies clearance directly, extrapolating far above",
        "the Teigen group maximum is exactly the failure mode Lee 2024 was written to demonstrate.",
        sep = " "
      ),
      source_name        = "CrCL"
    ),
    RRT_HEMODIAL_ACTIVE = list(
      description        = "Hemodialysis-active indicator (1 during a dialysis session, 0 otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (interdialytic / no dialysis running)",
      notes              = paste(
        "Time-varying within subject. Lee 2024 Methods: 'The model contained HD as a covariate for CL",
        "during HD (CLHD)'. Unlike the additive dialysis-arm idiom used by the sibling hemodialysis",
        "models in this library (Veinstein_2013_gentamicin.R, Dohmann_2025_piperacillin.R,",
        "Liesenfeld_2013_dabigatran.R), CL_HD here REPLACES CL_NHD rather than adding to it -- see the",
        "cl_total line in model() and the vignette Errata for the quantitative evidence.",
        "Lee 2024's patient underwent eleven 3-4 h sessions over the 569 h record (Table 1) on an",
        "FX CorDiax 60 high-flux dialyzer at a blood flow of 15 L/h, with a mean session length of",
        "3.7 h and interdialytic intervals of 18.6-68.5 h. The Teigen development cohort was dialysed",
        "thrice weekly. Set to 0 throughout for a subject who is not dialysed, which reduces the model",
        "to a one-compartment model with the CrCL-scaled CL_NHD alone.",
        sep = " "
      ),
      source_name        = "HD"
    )
  )

  covariatesDataExcluded <- list(
    IBW = list(
      description = "Ideal body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Not a model input. Lee 2024 uses IBW only as an INPUT to the Cockcroft-Gault CrCL that is",
        "supplied as the CRCL column; the model itself never references weight of any kind. Recorded",
        "here so the provenance of the CRCL column is not lost. Lee 2024's patient was 158 cm and",
        "66.9 kg with an ideal body weight greater than 1.25 x total body weight (Methods, Patient",
        "and ethics), which is the criterion the paper uses to classify her as obese for",
        "aminoglycoside-dosing purposes.",
        sep = " "
      )
    ),
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened but NOT in the model. Lee 2024's Discussion notes that 'both the CL and Vd of",
        "gentamicin increase with total body weight in obese populations' and that the a priori",
        "Teigen model carries NO total-body-weight term -- CrCL enters on ideal body weight only.",
        "That absence is the paper's central finding, not an omission in this transcription.",
        sep = " "
      )
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 1L,
    n_studies        = 1L,
    age_range        = "53 years (single patient)",
    weight_range     = "66.9 kg (single patient)",
    sex_female_pct   = 100,
    race_ethnicity   = "Korean (single patient)",
    disease_state    = paste(
      "APPLICATION population, Lee 2024: one 53-year-old obese Korean woman with end-stage renal",
      "disease on thrice-weekly intermittent hemodialysis, treated with gentamicin for a",
      "carbapenem-resistant Pseudomonas aeruginosa surgical-site infection. Height 158 cm, weight",
      "66.9 kg, BMI 26.8 kg/m2 (obese by the WHO Asia-Pacific adult cut-off of 25), ideal body weight",
      "> 1.25 x total body weight. Extensive critical-illness history: ST-elevation myocardial",
      "infarction with percutaneous coronary intervention, cardiogenic shock, post-infarction",
      "ventricular septal defect, two heart transplants with extracorporeal membrane oxygenation",
      "between them, continuous renal replacement therapy progressing to end-stage renal disease,",
      "and repeated surgical-site infections. On dialysis approximately 4.4 months at the time of",
      "the gentamicin course.",
      sep = " "
    ),
    dose_range       = paste(
      "Four intravenous infusions over 569 h (Lee 2024 Table 1): 140 mg at 140 mg/h (t = 0 h),",
      "110 mg at 88 mg/h (t = 72.72 h), 100 mg at 100 mg/h (t = 483.53 h), and 100 mg at 100 mg/h",
      "(t = 531.32 h). Three serum samples were drawn, at t = 72.72, 530.78 and 568.80 h; per the",
      "Discussion only samples taken at least 2 h after the end of an infusion AND at least 2 h after",
      "the end of a dialysis session were used, to avoid the delayed distribution and post-dialysis",
      "redistribution known to affect aminoglycosides in renal failure.",
      sep = " "
    ),
    regions          = "Republic of Korea (Seoul National University Bundang Hospital)",
    renal_function   = paste(
      "End-stage renal disease on thrice-weekly intermittent hemodialysis (FX CorDiax 60 high-flux",
      "dialyzer, blood flow 15 L/h, session length 3.7 h, interdialytic interval 18.6-68.5 h).",
      "Cockcroft-Gault CrCL on ideal body weight: mean 2.44 L/h (range 1.87-3.24), i.e. roughly",
      "twice the Teigen development cohort's group maximum of 1.24 L/h.",
      sep = " "
    ),
    notes            = paste(
      "IMPORTANT -- two distinct populations. The `population` block above describes Lee 2024's",
      "APPLICATION cohort (n = 1), i.e. the patient this model was USED to fit, not the cohort the",
      "model was ESTIMATED from. The DEVELOPMENT population is Teigen 2006 (closed access, not on",
      "disk): adults with end-stage renal disease receiving intermittent hemodialysis who had been",
      "dialysed for at least one month, modelled in NONMEM version 5. Everything Lee 2024 reports",
      "about that cohort is its group MEAN creatinine clearance (0.53 L/h) and group MAXIMUM",
      "(1.24 L/h); no subject count, demographics, sampling design or dialyzer list is quoted.",
      "Lee 2024's own contribution is the individual (MAP Bayesian / POSTHOC) estimates of Table 2,",
      "obtained twice -- once with the patient's measured CrCL and once with CrCL forced to the",
      "Teigen group mean -- in both NONMEM 7.4.4 and Monolix 2024R1. Those four parameter sets are",
      "reproduced in the validation vignette, not in this file, because they are individual",
      "realizations of this population model rather than models in their own right.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # ALL population parameters below are FIXED priors. Lee 2024 Methods:
    # "we employed the POSTHOC option using NONMEM 7.4.4", i.e. a MAP Bayesian
    # run in which the population parameters are held at the a priori values
    # and only the individual random effects are estimated. None of these
    # values was estimated in Lee 2024; all originate from Teigen 2006 and are
    # quoted by Lee 2024 in Methods, "Individual PK parameter estimation using
    # NONMEM and Monolix". Encoded with fixed() accordingly, following the
    # convention of Tong_2026_vancomycin_hughes.R and its siblings.
    #
    # Lee 2024 prints the trio as:
    #   "The population parameter values (standard error %) of CLNHD, CLHD, and
    #    Vd were 0.453 x CrCL/0.53 (0.9) L/h and 4.69 (0.86) L/h and,
    #    23.5 (0.91) L, respectively"
    # The parenthesized numbers are standard errors in per cent, not values.
    # ------------------------------------------------------------------------
    lcl <- fixed(log(0.453))
    label("Non-hemodialysis clearance CL_NHD at the reference CrCL of 8.833 mL/min = 0.53 L/h (L/h)")
    # Lee 2024 Methods: CLNHD = 0.453 x CrCL/0.53 L/h (SE 0.9%). Originates from Teigen 2006.

    lcl_hemodialysis <- fixed(log(4.69))
    label("TOTAL clearance while a hemodialysis session is running CL_HD (L/h)")
    # Lee 2024 Methods: CLHD = 4.69 L/h (SE 0.86%). Originates from Teigen 2006.
    # NOTE this is the total on-dialysis clearance, NOT an increment added to
    # CL_NHD -- see the cl_total comment in model().

    lvc <- fixed(log(23.5))
    label("Volume of distribution Vd (L)")
    # Lee 2024 Methods: Vd = 23.5 L (SE 0.91%). Originates from Teigen 2006.
    # No covariate is reported on Vd.

    # ------------------------------------------------------------------------
    # Inter-individual variability. Lee 2024 Methods states the a priori model
    # carried "log-normal between-subject variability for non-HD CL (CLNHD)
    # and Vd", and the Discussion quotes ONE magnitude:
    #   "According to the a priori information, the interindividual variability
    #    (%CV) of CLNHD was 51%".
    # omega^2 = log(CV^2 + 1) = log(0.51^2 + 1) = 0.23129, the house convention
    # (cf. Veinstein_2013_gentamicin.R). Under the alternative NONMEM shorthand
    # omega^2 = CV^2 the variance would be 0.2601 and the realized CV 54.6%
    # instead of 51.0%; the two readings are discussed in the vignette Errata.
    # The paper's own sanity check -- "the 95% confidence interval (CI) spans
    # from 0% to 200% of the population's typical value" -- is 1 +/- 1.96 x 0.51
    # on the linear scale and does not discriminate between them.
    #
    # The IIV on Vd is stated to EXIST but its magnitude is never printed, in
    # Lee 2024 or in any on-disk source; Teigen 2006 is closed access. It is
    # encoded as fixed(0) rather than invented -- simulations from this model
    # therefore carry between-subject variability on clearance only.
    # ------------------------------------------------------------------------
    etalcl ~ fixed(0.23129)  # Lee 2024 Discussion: %CV of CLNHD = 51%
    etalvc ~ fixed(0)        # Lee 2024 Methods: log-normal BSV on Vd exists; magnitude NOT published

    # ------------------------------------------------------------------------
    # Residual variability. Lee 2024 Methods states the a priori model used
    # "a combined error model", and that for the Monolix cross-check "we
    # utilized the sigma values for additive and proportional residual errors
    # from NONMEM" -- so both components exist, but NEITHER magnitude is
    # printed anywhere in the paper. Encoded as zero rather than invented,
    # following Takada_2025_vancomycin.R and Setiawan_2023_sulbactam.R.
    # Simulations from this model are residual-error-free.
    # ------------------------------------------------------------------------
    propSd <- fixed(0); label("Proportional residual SD (fraction; 0 -- not reported in the source)")
    addSd  <- fixed(0); label("Additive residual SD (mg/L; 0 -- not reported in the source)")
  })

  model({
    # 1. Unit alignment. The canonical CRCL column is in mL/min; the source
    #    equation is written in L/h with a reference of 0.53 L/h. The ratio is
    #    dimensionless, so the reference is expressed in mL/min directly:
    #    0.53 L/h = 530 mL/h = 8.8333 mL/min.
    #
    # 2. Individual PK parameters, and the clearance-regime switch.
    #
    #    Lee 2024 Methods gives the covariate model as the plain linear ratio
    #    "0.453 x CrCL/0.53" -- no power term is printed, so the exponent is
    #    exactly 1.
    #
    #    Lee 2024 describes "HD as a covariate for CL during HD (CLHD)", and
    #    CL_HD is the TOTAL clearance while a session runs -- it REPLACES
    #    CL_NHD rather than adding to it. This differs from the additive
    #    dialysis-arm idiom used elsewhere in this library and was settled
    #    against the paper's own published predictions rather than assumed:
    #    digitising the 310 plotted points of Fig. 1B (the panel in which CrCL
    #    is forced to the 0.53 L/h group mean, so the covariate factor is
    #    exactly 1 and no quantity is free to absorb the difference) and
    #    re-simulating this model at the Table 2 estimates gives a median
    #    absolute error of 1.1% for the switch form against 6.5% for the
    #    additive form. See vignette Errata for the full comparison.
    #
    #    IMPORTANT -- `cl` must be the clearance ACTUALLY IN EFFECT, not the
    #    interdialytic arm. rxode2 5.1.7 recognises the joint presence of `cl`
    #    and `vc` and solves the one-compartment system with its analytic
    #    linear-compartment kernel driven by those two variables, bypassing the
    #    d/dt() right-hand side. Assigning the interdialytic arm to `cl` and
    #    carrying the switch in a separate variable (the shape used by the
    #    sibling hemodialysis models in this library) therefore produces a
    #    model whose dialysis arm is silently INERT -- the reported cl/kel
    #    output columns still switch, but the solved concentrations do not.
    #    Writing the switch directly into `cl` is what makes it load-bearing.
    cl_hemodialysis <- exp(lcl_hemodialysis)  # L/h, total CL while dialysing
    cl <- (1 - RRT_HEMODIAL_ACTIVE) * exp(lcl + etalcl) * (CRCL / 8.8333) +
      RRT_HEMODIAL_ACTIVE * cl_hemodialysis   # L/h
    vc <- exp(lvc + etalvc)                   # L

    kel <- cl / vc

    # 4. One-compartment intravenous disposition. Gentamicin was given as a
    #    zero-order intravenous infusion (Lee 2024 Table 1 supplies an infusion
    #    rate for every dose), so the rate is set by the event record and there
    #    is no depot.
    d/dt(central) <- -kel * central

    # 5. Observation. Dose in mg over volume in L gives mg/L, the units of
    #    Figs. 1-2 ("Gentamicin serum concentration (mg/L)"). The assay was a
    #    chemiluminescent microparticle immunoassay with an LLOQ of
    #    0.3 ug/mL == 0.3 mg/L and a validated range of 0.3-10.00 ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
