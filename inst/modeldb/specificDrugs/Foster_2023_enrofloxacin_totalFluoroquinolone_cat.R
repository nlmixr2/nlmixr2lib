Foster_2023_enrofloxacin_totalFluoroquinolone_cat <- function() {
  description <- paste(
    "Veterinary (domestic cat). Two-compartment population PK model for TOTAL",
    "FLUOROQUINOLONE (the arithmetic sum of plasma enrofloxacin and its active",
    "metabolite ciprofloxacin, which Foster 2023 treats as one analyte because the",
    "two compounds' antibacterial effects are additive) after a single 5 mg/kg",
    "intravenous enrofloxacin dose infused over 30 minutes in 34 client-owned cats",
    "hospitalised for clinical illness and spanning normal to markedly reduced kidney",
    "function. Sparse sampling (3 samples per cat, 98 samples total) with nonlinear",
    "mixed-effects modelling in Phoenix NLME. Body weight, centred at the population",
    "value of 3.8 kg, is the only retained covariate and acts as a power effect on",
    "both the central volume and the elimination clearance; serum creatinine,",
    "symmetric dimethylarginine, blood urea nitrogen, age and sex were screened and",
    "none was significant, which is the paper's central finding (no dose adjustment",
    "is indicated for azotemic cats). Foster 2023 Table 1 labels the volumes 'mL/kg'",
    "and the clearances 'mL/kg/h', but those parameters are ABSOLUTE (mL and mL/h),",
    "not weight-normalised: only the absolute reading reproduces the paper's own",
    "Table 1 AUC0-24 of 46.0 ug*h/mL and the concentration range of Figure 1 -- see",
    "the lvc comment and the vignette Errata for the arithmetic. Between-subject",
    "variances and the residual-error magnitude are not reported anywhere in the",
    "paper or its supplement, so every eta and propSd is fixed(0) and the model",
    "simulates typical-value profiles unless the user supplies variances.",
    sep = " "
  )
  reference <- paste(
    "Foster JD, Abouraya M, Papich MG, Muma NA. (2023).",
    "Population pharmacokinetic analysis of enrofloxacin and its active metabolite",
    "ciprofloxacin after intravenous injection to cats with reduced kidney function.",
    "Journal of Veterinary Internal Medicine 37(6):2230-2240.",
    "doi:10.1111/jvim.16866.",
    sep = " "
  )
  vignette <- "Foster_2023_enrofloxacin_cat"
  units <- list(
    time = "h",
    dosing = "ug",
    concentration = "ug/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. The analyte is the summed enrofloxacin + ciprofloxacin
  # plasma concentration that Foster 2023 calls "total fluoroquinolone"
  # (Methods, "Therefore, the total fluroquinolone concentration (sum of
  # enrofloxacin and ciprofloxacin concentrations) was used for
  # pharmacokinetic analysis").
  compartmentData <- list(
    central = list(
      analyte = "enrofloxacin + ciprofloxacin (total fluoroquinolone)",
      units = "ug", specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "enrofloxacin + ciprofloxacin (total fluoroquinolone)",
      units = "ug", specimen = "plasma", verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed baseline body weight. Foster 2023 Results section 3.2: 'In this",
        "model body weight was centered around the mean population body weight, 3.8",
        "kg.' Enters as a power-of-ratio effect (WT / 3.8)^exponent on both the",
        "central volume V1 and the elimination clearance Cl1; the peripheral volume",
        "V2 and the intercompartmental clearance Cl2 carry no covariate ('Only body",
        "weight was found to significantly affect the volume and clearance of the",
        "first compartment'). The centred power form is Phoenix NLME's standard",
        "continuous-covariate model; the alternative centred-exponential reading",
        "exp(1.63 * (WT - 3.8)) would multiply V1 by 3463-fold at the cohort's",
        "largest cat (8.8 kg) and is physically impossible. Cohort median 3.8 kg,",
        "range 1.8-8.8 kg (Results section 3).",
        sep = " "
      ),
      source_name        = "weight (Foster 2023 Table 1 covariate names dVd1weight / dCl1weight)"
    )
  )

  # Screened in the stepwise forward-addition / backward-deletion covariate
  # search (Foster 2023 Methods paragraph 8, Figures 2 and 3) but NOT retained
  # in the final model: 'All other covariates were found to not significantly
  # affect volume or clearance of either compartment.' No point estimates are
  # published for any of them, so they are documented here rather than encoded.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age", units = "years", type = "continuous",
      notes = "Screened as a continuous covariate on Cl1, Cl2, V1 and V2 (Figure 3); not retained. Cohort mean 11.4 +/- 4.3 years."
    ),
    SEXF = list(
      description = "Female sex indicator", units = "(binary)", type = "binary",
      reference_category = "male (castrated)",
      notes = "Foster 2023 coded sex as female = 0 and male = 1 (Methods paragraph 8, the inverse of the canonical SEXF orientation); screened on Cl1, Cl2, V1 and V2 (Figure 2) and not retained. Cohort 18 spayed females / 16 castrated males."
    ),
    CREAT = list(
      description = "Serum creatinine concentration", units = "mg/dL", type = "continuous",
      notes = "Screened as a marker of kidney function; not retained. Cohort mean 7.5 +/- 6.0 mg/dL. Also used to stratify cats into normal (<= 2.4), moderate (2.5-10) and severe (> 10) kidney dysfunction for the noncompartmental comparison."
    ),
    BUN = list(
      description = "Blood urea nitrogen concentration", units = "mg/dL", type = "continuous",
      notes = "Screened as a marker of kidney function; not retained in THIS (total fluoroquinolone) model. Cohort mean 137.6 +/- 103.0 mg/dL. BUN IS retained on the metabolic formation clearance in the companion parent + metabolite model, Foster_2023_enrofloxacin_ciprofloxacin_cat."
    ),
    SDMA = list(
      description = "Symmetric dimethylarginine concentration", units = "ug/dL", type = "continuous",
      notes = "Screened as a marker of glomerular filtration rate; not retained. Cohort mean 39.8 +/- 27.0 ug/dL. SDMA has no canonical entry in inst/references/covariate-columns.md; none is proposed here because the covariate is documentation only and carries no published point estimate."
    )
  )

  population <- list(
    species        = "cat (domestic, client-owned)",
    n_subjects     = 34L,
    n_studies      = 1L,
    age_range      = "mean 11.4 +/- 4.3 years",
    weight_range   = "1.8-8.8 kg (median 3.8 kg)",
    weight_median  = "3.8 kg",
    sex_female_pct = 52.9,
    disease_state  = paste(
      "Hospitalised for clinical illness with suspected bacterial infection and",
      "variable kidney function: 9 cats with normal kidney function (median serum",
      "creatinine 1.6 mg/dL, range 0.9-2.2), 15 with moderate kidney dysfunction",
      "(median 5.6, range 2.7-8.6) and 10 with severe kidney dysfunction (median",
      "13.3, range 11.9-21.5). Cats with pre-existing retinal lesions or blindness,",
      "anemia (hematocrit < 20%), hypoxia, congestive heart failure, hypotension or",
      "hypovolemia were excluded, as were cats given any fluoroquinolone in the",
      "preceding 7 days.",
      sep = " "
    ),
    renal_function = "Serum creatinine mean 7.5 +/- 6.0 mg/dL; BUN mean 137.6 +/- 103.0 mg/dL; SDMA mean 39.8 +/- 27.0 ug/dL",
    dose_range     = "Single 5 mg/kg enrofloxacin (Baytril 2.27%, Elanco) diluted 1:1 with sterile saline and infused intravenously over 30 minutes through a peripheral catheter; mean total dose 20.5 +/- 7.4 mg",
    sampling       = "Sparse: each cat scheduled for 3 plasma samples in the 24 h after dosing, randomly assigned to one of 7 optimised schedules (Table S1); 98 samples analysed from 34 cats (4 cats contributed only 2 usable samples). Plasma enrofloxacin and ciprofloxacin by validated HPLC.",
    regions        = "United States (single referral hospital, enrolment 2019-01-01 to 2021-09-01)",
    notes          = paste(
      "Baseline demographics are in Foster 2023 Results section 3 (first paragraph).",
      "The 2-compartment structure was selected over 1-compartment on mean AIC",
      "(-23.61 vs -8.97) and the multiplicative residual model over proportional and",
      "additive alternatives (AIC 207.14 vs 209.61 vs 230.24). Internal validation",
      "was a 300-replicate bootstrap (Table S2); the bootstrap flagged V2 as poorly",
      "estimated. All cats continued clinically indicated antibiotic therapy after",
      "the sampling window, so the model describes first-dose kinetics only.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Structural disposition parameters -- Foster 2023 Table 1
    # =====================================================================
    # UNITS. Table 1 prints "theta Vd1 (mL/kg) 4794.7" and "theta Cl1
    # (mL/kg/h) 254.3", i.e. it labels the primary parameters as
    # weight-normalised. They are NOT. Reading them as absolute (mL and
    # mL/h) with the dose entered as the TOTAL amount in ug is the only
    # reading that reproduces the paper's own numbers:
    #
    #   absolute reading, mean total dose 20.5 mg = 20500 ug
    #     AUC(0-inf) = 20500 / 254.3                      = 80.6 ug*h/mL
    #     AUC(0-24)  = 80.6 * 0.569 (2-cmt fraction)      = 45.9 ug*h/mL
    #     -> Table 1 reports AUC(0-24) = 46.0 ug*h/mL     MATCH
    #     C(0)       = 20500 / 4794.7                     =  4.3 ug/mL
    #     -> Figure 1 spaghetti plot spans ~0.2-8 ug/mL   MATCH
    #     Dose / AUC(0-24) = 20500 / 45.9                 = 447 mL/h
    #     -> Results 3.1 "observed clearance" median 440.8 MATCH
    #
    #   weight-normalised reading, dose 5 mg/kg = 5000 ug/kg
    #     AUC(0-inf) = 5000 / 254.3                       = 19.7 ug*h/mL
    #     AUC(0-24)                                       = 11.2 ug*h/mL
    #     -> vs Table 1's 46.0 and the observed NCA median 41.44  4x LOW
    #     C(0)       = 5000 / 4794.7                       = 1.0 ug/mL
    #     -> vs Figure 1's early concentrations of 3-8 ug/mL      4x LOW
    #
    # The absolute reading is therefore used throughout, doses are entered
    # in ug, and the "/kg" in the Table 1 row labels is treated as a
    # printing error. It also makes the fitted weight exponents
    # interpretable: an exponent near 1 on an ABSOLUTE volume is ordinary
    # body-size scaling, whereas the same exponent stacked on top of an
    # already-per-kg volume would imply total scaling with WT^2.63. See
    # vignette Errata.
    lvc <- log(4794.7); label("Central volume of distribution V1 for a 3.8 kg cat (mL)")  # Table 1: theta Vd1 = 4794.7 (SE 359.6, CV% 7.22)
    lvp <- log(2.9e4); label("Peripheral volume of distribution V2 (mL)")                 # Table 1: Vd2 = 2.9 x 10^4 (SE 5.1 x 10^5); the paper flags this parameter as poorly estimated (Results 3.3, Discussion paragraph 3)
    lcl <- log(254.3); label("Elimination clearance Cl1 for a 3.8 kg cat (mL/h)")         # Table 1: theta Cl1 = 254.3 (SE 31.8, CV% 12.49)
    lq <- log(132.1); label("Intercompartmental clearance Cl2 (mL/h)")                    # Table 1: Cl2 = 132.1 (SE 27.4, CV% 20.71)

    # =====================================================================
    # Body-weight covariate -- Foster 2023 Abstract and Table 1
    # =====================================================================
    # Values come from the Abstract, which prints three significant figures
    # ("effect 1.63, SE 0.19" for volume and "effect 1.63, SE 0.27" for
    # clearance); Table 1 rounds both to 1.6. The Table 1 CV% for the
    # clearance exponent (16.38) reproduces exactly as 0.27 / 1.63, which
    # confirms 1.63 rather than 1.6 is the underlying estimate.
    e_wt_vc <- 1.63; label("Body-weight power exponent on V1 (unitless)")   # Abstract: "volume of distribution (effect 1.63, SE 0.19, P < .01)"; Table 1 dVd1weight = 1.6 (SE 0.2)
    e_wt_cl <- 1.63; label("Body-weight power exponent on Cl1 (unitless)")  # Abstract: "clearance (effect 1.63, SE 0.27, P < .01)"; Table 1 dCl1weight = 1.6 (SE 0.3)

    # =====================================================================
    # Between-subject variability
    # =====================================================================
    # Foster 2023 Methods paragraph 7 states that interindividual
    # variability was exponential, P_i = P_pop * exp(eta_iP), with eta ~
    # N(0, omega^2), and Figure 2 / Figure 3 plot the empirical Bayes etas
    # for all four structural parameters (nCl, nCl2, nV, nV2 -- so all four
    # carried an eta). NO omega, omega^2, SD or CV of any eta is reported
    # in the paper, in Table S2 (bootstrap, which lists only the six
    # structural / covariate parameters) or in any other supplement; the
    # Table 1 CV% column is the relative standard error of the typical
    # value, not between-subject variability. Per the standing policy for
    # unreported IIV with structural values present, each eta is fixed(0)
    # rather than invented. Users with variance estimates can unfix them.
    etalvc ~ fixed(0)
    etalvp ~ fixed(0)
    etalcl ~ fixed(0)
    etalq ~ fixed(0)

    # =====================================================================
    # Residual unexplained variability
    # =====================================================================
    # Methods paragraph 7: "A multiplicative model was chosen (among
    # additive, log-additive, power, and mixed error models) to describe
    # the residual random variability", Cobs_ij = Cpred_ij * (1 + eps_ij),
    # which is the proportional form in nlmixr2. Results 3.1 confirms the
    # choice on AIC. The magnitude of sigma is not reported anywhere, so
    # it is fixed(0) on the same policy as the etas above; simulations are
    # noise-free typical-value profiles.
    propSd <- fixed(0); label("Proportional residual error (fraction; not published)")  # Methods paragraph 7 states the multiplicative form; no sigma value is given
  })

  model({
    # 1. Individual parameters. Body weight enters V1 and Cl1 as a power of
    #    the ratio to the 3.8 kg centring value (Results 3.2). V2 and Cl2
    #    carry no covariate.
    vc <- exp(lvc + etalvc) * (WT / 3.8)^e_wt_vc
    vp <- exp(lvp + etalvp)
    cl <- exp(lcl + etalcl) * (WT / 3.8)^e_wt_cl
    q <- exp(lq + etalq)

    # 2. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. ODE system. The enrofloxacin dose is given intravenously (30-minute
    #    infusion), so the user doses `central` directly with a rate.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 4. Observation. Amounts are in ug and volumes in mL, so Cc comes out
    #    in ug/mL, the unit of Figure 1's y-axis.
    Cc <- central / vc

    Cc ~ prop(propSd)
  })
}
