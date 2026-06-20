Wittau_2015_meropenem <- function() {
  description <- "Two-compartment intravenous population PK model for meropenem in morbidly obese adults (Wittau 2015). Allometric scaling on fat-free mass with a reference FFM of 53 kg. Unbound meropenem concentrations in subcutaneous tissue and peritoneal fluid are described as the plasma concentration multiplied by the site-to-plasma AUC ratios FSC and FPF (the final model assumed very rapid equilibration with plasma, so SC and PF are not carried as separate ODE states)."
  reference <- paste(
    "Wittau M, Scheele J, Kurlbaum M, Brockschmidt C, Wolf AM,",
    "Hemper E, Henne-Bruns D, Bulitta JB.",
    "Population pharmacokinetics and target attainment of meropenem",
    "in plasma and tissue of morbidly obese patients after laparoscopic",
    "intraperitoneal surgery.",
    "Antimicrob Agents Chemother. 2015;59(10):6241-6247.",
    "doi:10.1128/AAC.00259-15. PMID 26248373.",
    sep = " "
  )
  vignette <- "Wittau_2015_meropenem"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass; derived from total body weight, height, and sex via the Janmahasatian et al. (2005) formula. Used for allometric scaling on CL, CLd, V1, and V2 with reference FFM = 53 kg.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Reference FFM 53 kg per Table 2 footnote b. Cohort range 52.3 - 94.0 kg; geometric mean 70.3 kg, 27.4% CV used for Monte Carlo simulations (Methods 'Monte Carlo simulations').",
      source_name        = "FFM"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Considered alongside FFM for allometric scaling but not retained in the final model: 'we considered FFM (15) to account for the significantly altered body composition of morbidly obese patients' (Methods 'Parameter variability model and covariate effects'). Discussion notes the cohort lacked the WT-vs-FFM contrast needed to rank them definitively; FFM was preferred on physiological grounds. Cohort range 116 - 203 kg."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Recorded in Table 1 (mean 40 +/- 7.87 years, range 31 - 49) but not screened as a covariate; the paper explicitly avoided empirical covariate model building given the small sample size (n = 5)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Cohort included 3 of 5 female (Table 1) but sex was not screened as a covariate; small-sample exclusion applies."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Inclusion criterion (BMI >= 40) and cohort descriptor (mean 54.2, range 47.6 - 62.3 kg/m^2; Table 1) but not used as a model covariate."
    ),
    HEIGHT = list(
      description = "Standing height",
      units       = "cm",
      type        = "continuous",
      notes       = "Used internally to derive FFM via Janmahasatian; not used as a model covariate. Cohort mean 170 +/- 14.5 cm."
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Inclusion criterion (CrCL > 30 mL/min) excluded severe renal insufficiency; not screened as a covariate. Serum creatinine reported in Table 1 (median 75 umol/L, range 64 - 80)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 5L,
    n_studies      = 1L,
    age_range      = "31 - 49 years (median 39)",
    age_median     = "39 years",
    weight_range   = "116 - 203 kg total body weight (median 163)",
    weight_median  = "163 kg",
    ffm_range      = "52.3 - 94.0 kg (median 63.8); reference 53 kg for allometric scaling",
    bmi_range      = "47.6 - 62.3 kg/m^2 (median 51.9)",
    sex_female_pct = 60,
    disease_state  = "Morbidly obese (BMI >= 40), noninfected adults undergoing elective laparoscopic intraperitoneal surgery (5 sleeve gastrectomies, 1 abdominal-wall hernia repair; one excluded from analysis). Severe renal insufficiency (CrCL <= 30 mL/min), severe hepatic disease, and concurrent valproic acid were exclusion criteria.",
    dose_range     = "1 g meropenem every 8 h as a 15-min IV infusion; PK sampled after the fourth dose (1 g q8h, IV; intraoperative + post-operative).",
    regions        = "Germany (University of Ulm, Department of Visceral Surgery; August 2012 - January 2013).",
    notes          = "ClinicalTrials.gov NCT01407965. Rich sampling at 0.5, 1, 2, 3, 5, and 8 h after the fourth dose in plasma, subcutaneous adipose interstitial space fluid (microdialysis), and peritoneal fluid (microdialysis). All concentrations analyzed simultaneously via population modeling in S-ADAPT (importance sampling, pmethod = 4). Demographic summary from Table 1; modelling details from Materials and Methods 'Population modeling and Monte Carlo simulations'."
  )

  ini({
    # =========================================================================
    # Structural disposition parameters (Table 2 column 'Population mean';
    # Table 2 footnote b: 'Estimates represent a patient with normal body size
    # (i.e., 53 kg fat-free mass) and are based on an allometric body size
    # model.').
    # =========================================================================
    lcl  <- log(18.7); label("Total clearance CL (L/h) at FFM = 53 kg")                              # Table 2 row 'Total clearance'
    lvc  <- log(21.5); label("Central volume of distribution V1 (L) at FFM = 53 kg")                 # Table 2 row 'Volume of distribution of central compartment'
    lvp  <- log(6.16); label("Peripheral volume of distribution V2 (L) at FFM = 53 kg")              # Table 2 row 'Volume of distribution of peripheral compartment'
    lq   <- log(29.4); label("Distribution clearance CLd (L/h) at FFM = 53 kg")                      # Table 2 row 'Distribution clearance between central and peripheral compartments'

    # =========================================================================
    # Site-to-plasma AUC ratios for unbound meropenem in subcutaneous tissue
    # and peritoneal fluid. The final model multiplies the central plasma
    # concentration by FSC and FPF (Materials and Methods 'Structural model'
    # paragraph 1; Table 2 footnote d: 'The half-lives of equilibration ...
    # were rapid (equilibration half-life, <0.5 h). The final model assumed a
    # very rapid equilibration between the respective peripheral site and
    # plasma.'). FSC and FPF therefore enter the model as multiplicative
    # scalars, not as additional ODE states. Naming follows the precedent set
    # by Landersdorfer_2009_moxifloxacin (lfcortical / lfcancellous and
    # Ccortical / Ccancellous outputs), which uses the same ratio-scaling
    # construction.
    # =========================================================================
    lfsc <- log(0.721); label("Subcutaneous-tissue:plasma AUC ratio FSC (unitless)")                 # Table 2 row 'Ratio of AUC values in subcutaneous tissue and plasma'
    lfpf <- log(0.943); label("Peritoneal-fluid:plasma AUC ratio FPF (unitless)")                    # Table 2 row 'Ratio of AUC values in peritoneal fluid and plasma'

    # =========================================================================
    # Allometric exponents on FFM. The paper states 'standard allometric
    # scaling' (Methods 'Parameter variability model and covariate effects')
    # without printing the exponent values; the Discussion contrasts a 15%
    # difference between allometric and linear over the 1.8-fold cohort FFM
    # range, consistent with the canonical 0.75 / 1.0 (1.8^0.75 / 1.8^1.0 = 0.86,
    # a 14% difference). Encoded as fixed() because the paper holds them
    # constant and reports no uncertainty.
    # =========================================================================
    e_ffm_cl_q  <- fixed(0.75); label("Allometric exponent on CL and CLd (unitless; fixed)")        # Methods 'Parameter variability model and covariate effects'; Discussion paragraph 4
    e_ffm_vc_vp <- fixed(1.00); label("Allometric exponent on V1 and V2 (unitless; fixed)")         # Methods 'Parameter variability model and covariate effects'; Discussion paragraph 4

    # =========================================================================
    # Between-subject variability. The model is exponential (Methods
    # 'Parameter variability model and covariate effects'), with a major-
    # diagonal omega matrix (Results paragraph 4: 'The model included a
    # major-diagonal variance-covariance matrix'). Table 2 footnote c clarifies
    # that the 'Between-subject variability' column reports the apparent CV of
    # the log-normal distribution; we convert via omega^2 = log(CV^2 + 1).
    # =========================================================================
    etalcl  ~ 0.001489    # CV  3.86% -> log(0.0386^2 + 1) = 0.001489 - Table 2 'Total clearance' BSV
    etalq   ~ 1.43607     # CV   179% -> log(1.79^2   + 1) = 1.43607  - Table 2 'Distribution clearance' BSV
    etalvc  ~ 0.010758    # CV  10.4% -> log(0.104^2  + 1) = 0.010758 - Table 2 'Volume of distribution of central compartment' BSV
    etalvp  ~ 0.001779    # CV  4.22% -> log(0.0422^2 + 1) = 0.001779 - Table 2 'Volume of distribution of peripheral compartment' BSV
    etalfsc ~ 0.81277     # CV   112% -> log(1.12^2   + 1) = 0.81277  - Table 2 'Ratio of AUC values in subcutaneous tissue and plasma' BSV
    etalfpf ~ 0.09229     # CV  31.1% -> log(0.311^2  + 1) = 0.09229  - Table 2 'Ratio of AUC values in peritoneal fluid and plasma' BSV

    # =========================================================================
    # Residual error - combined additive + proportional per observation site
    # (Table 2 footnote a: 'The standard deviations of the additive and
    # proportional residual errors were 0.0235 mg/liter and 21.6% in plasma,
    # 1.26 mg/liter and 9.58% in SC tissue, and 0.617 mg/liter and 10.3% in
    # peritoneal fluid.'). Concentrations carry units of mg/L; proportional
    # SDs are reported as fractions.
    # =========================================================================
    propSd     <- 0.216;  label("Proportional residual SD on plasma Cc (fraction)")                  # Table 2 footnote a
    addSd      <- 0.0235; label("Additive residual SD on plasma Cc (mg/L)")                          # Table 2 footnote a
    propSd_Csc <- 0.0958; label("Proportional residual SD on subcutaneous-tissue Csc (fraction)")    # Table 2 footnote a
    addSd_Csc  <- 1.26;   label("Additive residual SD on subcutaneous-tissue Csc (mg/L)")            # Table 2 footnote a
    propSd_Cpf <- 0.103;  label("Proportional residual SD on peritoneal-fluid Cpf (fraction)")       # Table 2 footnote a
    addSd_Cpf  <- 0.617;  label("Additive residual SD on peritoneal-fluid Cpf (mg/L)")               # Table 2 footnote a
  })

  model({
    # -----------------------------------------------------------------------
    # 1. Allometric scaling of disposition parameters on FFM (reference 53 kg;
    #    Table 2 footnote b). Standard exponents from the paper's 'standard
    #    allometric scaling' wording: 0.75 on CL and CLd, 1.0 on V1 and V2.
    # -----------------------------------------------------------------------
    ffm_ratio <- FFM / 53

    cl <- exp(lcl  + etalcl)  * ffm_ratio^e_ffm_cl_q
    q  <- exp(lq   + etalq)   * ffm_ratio^e_ffm_cl_q
    vc <- exp(lvc  + etalvc)  * ffm_ratio^e_ffm_vc_vp
    vp <- exp(lvp  + etalvp)  * ffm_ratio^e_ffm_vc_vp

    # -----------------------------------------------------------------------
    # 2. Site-to-plasma AUC ratios (no body-size dependence in the source).
    # -----------------------------------------------------------------------
    fsc <- exp(lfsc + etalfsc)
    fpf <- exp(lfpf + etalfpf)

    # -----------------------------------------------------------------------
    # 3. Micro-rate constants.
    # -----------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # -----------------------------------------------------------------------
    # 4. Two-compartment IV disposition (Fig 2). State variables carry drug
    #    amount; doses (mg) and volumes (L) give plasma Cc directly in mg/L.
    # -----------------------------------------------------------------------
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # -----------------------------------------------------------------------
    # 5. Observation variables. Cc is plasma; Csc and Cpf are the unbound
    #    site-to-plasma scaled concentrations (Materials and Methods
    #    'Structural model' paragraph 1; Table 2 footnote d).
    # -----------------------------------------------------------------------
    Cc  <- central / vc
    Csc <- fsc * Cc
    Cpf <- fpf * Cc

    Cc  ~ add(addSd)     + prop(propSd)
    Csc ~ add(addSd_Csc) + prop(propSd_Csc)
    Cpf ~ add(addSd_Cpf) + prop(propSd_Cpf)
  })
}
