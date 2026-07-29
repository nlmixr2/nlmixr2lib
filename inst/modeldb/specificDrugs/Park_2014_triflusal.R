Park_2014_triflusal <- function() {
  description <- paste(
    "One-compartment population PK with first-order metabolite-formation",
    "kinetics for the active triflusal metabolite hydroxy-4-(trifluoromethyl)",
    "benzoic acid (HTB) in healthy Korean male volunteers, with a binary",
    "probability PD model for inhibition of platelet aggregation (IPA).",
    "Triflusal is an antiplatelet prodrug; only HTB is measured analytically.",
    "NONMEM ADVAN2 TRANS2 is used by the source paper -- the canonical depot",
    "compartment carries triflusal and the canonical first-order rate constant",
    "(here `lka`) plays the role of the paper's HTB formation rate constant kf",
    "(0.341 1/h). Apparent oral clearance CL/F (0.200 L/h at 71.65 kg) and",
    "apparent oral volume V/F (8.300 L at 71.65 kg) describe HTB disposition;",
    "F absorbs the unknown fraction of triflusal converted to HTB. Body weight",
    "is the only retained covariate and enters as a power on CL/F (exponent",
    "0.845) and direct proportionality on V/F (exponent fixed to 1). PD",
    "endpoint is binary IPA = 1 when platelet aggregation < 74% else 0; the",
    "instantaneous probability of IPA is a sigmoid Hill function of HTB",
    "concentration, prob_ipa = Cc^gamma / (EC50^gamma + Cc^gamma), with EC50",
    "= 84.9 ug/mL and gamma = 19.2 (BSV on gamma fixed to 0). The Hill",
    "exponent is very steep (quantal-like concentration-response). Parameter",
    "values from Park 2014 Table 2 Estimates column."
  )
  reference <- paste(
    "Park SM, Lee J, Seong SJ, Park JG, Gwon MR, Lim MS, Lee HW, Yoon YR,",
    "Yang DH, Kwon KI, Han S. Population pharmacokinetic and pharmacodynamic",
    "modeling of transformed binary effect data of triflusal in healthy",
    "Korean male volunteers: a randomized, open-label, multiple dose,",
    "crossover study. BMC Pharmacology and Toxicology 2014;15:75.",
    "doi:10.1186/2050-6511-15-75."
  )
  vignette <- "Park_2014_triflusal"
  units    <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed baseline. Reference 71.65 kg is the cohort median",
        "(Park 2014 Methods page 4 / Table 2 footnote -- 'subject whose",
        "weight = 71.65 kg'). Power exponent 0.845 on CL/F (Park 2014",
        "Table 2 theta_4); direct proportionality (exponent fixed to 1)",
        "on V/F (Park 2014 Table 2 V/F equation V_d/F = theta_2 *",
        "(weight/71.65)). Cohort weights ranged 53.3 - 89.7 kg with mean",
        "70.8 +/- 9.0 kg (Table 1)."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 34L,
    n_studies      = 1L,
    age_range      = "21 - 28 years",
    age_median     = "24.1 years (SD 1.7)",
    weight_range   = "53.3 - 89.7 kg",
    weight_median  = "70.8 kg (SD 9.0)",
    height_range   = "167.1 - 184.2 cm (mean 176.1, SD 4.9)",
    sex_female_pct = 0,
    race_ethnicity = c(Asian = 100),
    disease_state  = paste(
      "Healthy adult male volunteers; subjects were within 20% of their",
      "ideal body weight (IBW = (height_cm - 100) * 0.9) and passed",
      "physical examination plus routine laboratory tests (blood",
      "hematology, biochemistry, prothrombin time, bleeding time,",
      "urinalysis)."
    ),
    dose_range     = paste(
      "Oral triflusal 900 mg loading dose on Day 1 followed by a 600",
      "mg/day maintenance dose (given as two 300 mg capsules once daily)",
      "from Day 2 to Day 9 inclusive (Park 2014 Methods: Study drugs and",
      "dosage regimen). The Figure 1 caption lists only the PK sampling",
      "days (2, 3, 5, 7, 8, 9), not the dosing days; daily dosing across",
      "Day 2 - Day 9 is the authoritative regimen and matches the",
      "paper's reported accumulation factor 2.26 and Css,min 103.5",
      "ug/mL."
    ),
    regions        = paste(
      "South Korea -- Kyungpook National University Hospital Clinical",
      "Trial Center (KNUH CTC), Daegu."
    ),
    cris_registration = "KCT0001299 (Clinical Research Information Service, Korea)",
    notes          = paste(
      "38 volunteers enrolled August - September 2008; 34 completed and",
      "were included in the population PK/PD analysis (4 dropped out:",
      "severe dental pain, medication administration error, missing",
      "follow-up urinalysis, platelet aggregation 1% of baseline in",
      "second period). 476 HTB plasma concentration data points and 340",
      "platelet aggregation data points were used for model",
      "construction. HTB was quantified by HPLC-MS/MS with calibration",
      "range 1 - 300 ug/mL (LLOQ = 1 ug/mL; intra-day and inter-day CV",
      "< 15%; Park 2014 Methods: HTB plasma concentration measurements).",
      "Platelet aggregation was measured by light transmission",
      "aggregometry on platelet-rich plasma with 0.5 uL arachidonic acid",
      "as agonist (Chronolog 490-4D); PD observations were dichotomised",
      "at 74% aggregation per a published reference range for normal",
      "platelet-rich plasma. Age, height, weight, albumin, creatinine,",
      "AST, ALT, and Cockcroft-Gault CLCR were screened as covariates;",
      "only weight was retained for CL/F and V/F. Creatinine clearance",
      "was not selected even though HTB is > 60% renally excreted -- the",
      "authors attribute this to confounding with weight."
    )
  )

  ini({
    # -----------------------------------------------------------------
    # Structural PK parameters -- Park 2014 Table 2 (Estimates column).
    # NONMEM ADVAN2 TRANS2 (1-cmt first-order absorption). The paper's
    # k_f (first-order formation rate constant of HTB from triflusal)
    # is mathematically identical to ADVAN2 k_a and is encoded as the
    # canonical `lka`. Reference body weight 71.65 kg is the cohort
    # median used in the source paper's covariate equations.
    # -----------------------------------------------------------------
    lka <- log(0.341); label("First-order HTB formation rate from triflusal (kf == ka, 1/h)") # Park 2014 Table 2 theta_3 (k_f = 0.341 1/h, RSE 15.1%)
    lcl <- log(0.200); label("Apparent oral clearance of HTB at WT = 71.65 kg (CL/F, L/h)")    # Park 2014 Table 2 theta_1 (CL/F = 0.200 L/h, RSE 2.7%)
    lvc <- log(8.300); label("Apparent oral volume of distribution of HTB at WT = 71.65 kg (V/F, L)") # Park 2014 Table 2 theta_2 (V_d/F = 8.300 L, RSE 2.7%)

    # -----------------------------------------------------------------
    # Covariate effects on PK -- Park 2014 Table 2 footnotes.
    #   CL/F = theta_1 * (WT/71.65)^theta_4   (power form)
    #   V/F  = theta_2 * (WT/71.65)           (linear; exponent 1)
    # The V/F exponent is reported in the equation as no exponent
    # (i.e. linear proportionality, equivalent to exponent = 1) and
    # is not separately estimated -- encode it as fixed(1).
    # -----------------------------------------------------------------
    e_wt_cl <- 0.845;        label("Allometric power exponent of WT on CL/F (unitless)") # Park 2014 Table 2 theta_4 = 0.845 (RSE 17.4%)
    e_wt_vc <- fixed(1.000); label("Allometric power exponent of WT on V/F (unitless, fixed at 1 per V/F = theta_2 * (WT/71.65))") # Park 2014 Table 2 V/F equation (no separate exponent reported)

    # -----------------------------------------------------------------
    # Structural PD parameters -- Park 2014 Table 2 (Estimates column).
    # Binary probability PD: prob_ipa = Cc^gamma / (EC50^gamma + Cc^gamma).
    # gamma is the canonical Hill coefficient (`lhill`); EC50 in ug/mL.
    # The shape is quantal-like (gamma ~ 19) -- the concentration-
    # response transitions sharply around EC50.
    # -----------------------------------------------------------------
    lec50 <- log(84.9); label("HTB concentration at which P(IPA) = 0.5 (EC50, ug/mL)")        # Park 2014 Table 2 theta_5 (EC50 = 84.9 ug/mL, RSE 4.0%)
    lhill <- log(19.2); label("Hill / shape exponent of IPA probability vs HTB concentration (unitless)") # Park 2014 Table 2 theta_6 (gamma = 19.2, RSE 22.4%)

    # -----------------------------------------------------------------
    # BSV (between-subject variability) -- Park 2014 Table 2 % CV.
    # Log-normal etas, omega^2 = log(CV^2 + 1) with CV as a fraction:
    #   CL/F  CV 14.9% -> omega^2 = log(1 + 0.149^2) = 0.021960
    #   V/F   CV  9.5% -> omega^2 = log(1 + 0.095^2) = 0.008984
    #   k_f   CV 88.0% -> omega^2 = log(1 + 0.880^2) = 0.573537
    #   EC50  CV 21.8% -> omega^2 = log(1 + 0.218^2) = 0.046430
    #   gamma 0 (Fixed -- no IIV on gamma; not declared below)
    # -----------------------------------------------------------------
    etalcl   ~ 0.021960 # Park 2014 Table 2 omega_1^2 BSV CL/F 14.9% CV (RSE 21.0%; shrinkage 0.6%)
    etalvc   ~ 0.008984 # Park 2014 Table 2 omega_2^2 BSV V/F  9.5% CV (RSE 57.8%; shrinkage 40.2%)
    etalka   ~ 0.573537 # Park 2014 Table 2 omega_3^2 BSV k_f 88.0% CV (RSE 26.5%; shrinkage 10.5%)
    etalec50 ~ 0.046430 # Park 2014 Table 2 omega_4^2 BSV EC50 21.8% CV (RSE 28.4%; shrinkage 8.3%)

    # -----------------------------------------------------------------
    # Residual error -- Park 2014 Table 2 sigma^2 = 0.098 (proportional
    # only; additive error was not selected per Table 2 footnote b).
    # NONMEM proportional EPS reports the variance directly; nlmixr2's
    # `prop(propSd)` carries the SD on the linear concentration scale,
    # so propSd = sqrt(sigma^2) = sqrt(0.098) ~ 0.3130 (31.3% CV).
    # -----------------------------------------------------------------
    propSd <- 0.3130; label("Proportional residual error on HTB Cc (fraction)") # Park 2014 Table 2 sigma^2 = 0.098 (RSE 6.5%; shrinkage 11.5%); propSd = sqrt(0.098)
  })

  model({
    # -----------------------------------------------------------------
    # 1. Individual PK parameters.
    #    CL/F carries a power weight effect; V/F carries linear weight
    #    proportionality (Park 2014 Table 2 footnotes).
    # -----------------------------------------------------------------
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 71.65)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 71.65)^e_wt_vc

    kel <- cl / vc

    # -----------------------------------------------------------------
    # 2. ODE system: 1-cmt first-order metabolite formation. The
    #    depot compartment receives triflusal dose; it converts to HTB
    #    in central at rate `ka` (== paper k_f). HTB is eliminated from
    #    central with rate kel = CL/V.
    # -----------------------------------------------------------------
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Dose in mg / V in L gives mg/L = ug/mL (matches paper units).
    Cc <- central / vc

    # -----------------------------------------------------------------
    # 3. PD: binary probability of IPA. prob_ipa is a deterministic
    #    transform of Cc (the paper's binary probability model).
    #    Users simulating binary outcomes draw rbinom(1, 1, prob_ipa)
    #    in post-processing -- see the validation vignette.
    # -----------------------------------------------------------------
    ec50 <- exp(lec50 + etalec50)
    hill <- exp(lhill)

    prob_ipa <- Cc^hill / (ec50^hill + Cc^hill)

    Cc ~ prop(propSd)
  })
}
