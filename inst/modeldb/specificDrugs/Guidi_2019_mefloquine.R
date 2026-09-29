Guidi_2019_mefloquine <- function() {
  description <- paste(
    "Two-compartment population PK model of mefloquine (MQ) with first-order",
    "absorption in 48 Kenyan children aged 6-59 months with uncomplicated",
    "Plasmodium falciparum malaria, treated with the fixed-dose",
    "artesunate-mefloquine (ASMQ) dispersible tablet (Guidi 2019).",
    "Relative bioavailability F1 is fixed at 1 and carries between-subject",
    "variability. Body weight scales both clearances allometrically",
    "(exponent 0.75) and both volumes linearly (exponent 1). The absorption",
    "rate constant is estimated separately for the first treatment day and",
    "for the second and third days, the authors attributing the increase to",
    "resolution of malaria-related gastrointestinal disturbance after the",
    "first artesunate dose, and decreases with age. The model was externally",
    "validated against sparse data from 378 further children across six",
    "African centres. The companion artesunate / dihydroartemisinin model",
    "from the same trial is modellib('Guidi_2019_artesunate')."
  )
  reference <- paste(
    "Guidi M, Mercier T, Aouri M, Decosterd LA, Csajka C, Ogutu B, Carn G,",
    "Kiechel JR. Population pharmacokinetics and pharmacodynamics of the",
    "artesunate-mefloquine fixed dose combination for the treatment of",
    "uncomplicated falciparum malaria in African children.",
    "Malaria Journal 2019;18:139. doi:10.1186/s12936-019-2754-6"
  )
  vignette <- "Guidi_2019_artesunate_mefloquine"

  # Mefloquine concentrations are reported in ng/mL throughout the paper
  # (Table 4 Cmax 2874 ng/mL, LOQ 2.5 ng/mL) but exposures are reported in
  # mg/L*h (Table 4 AUC0-inf 650; Results AUC0-day63 725). With doses in mg
  # and volumes in L, `central / vc` is mg/L (= ug/mL); multiply by 1000 to
  # compare against the paper's ng/mL concentrations.
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L"
  )

  compartmentData <- list(
    depot = list(analyte = "mefloquine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "mefloquine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "mefloquine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling of all four disposition parameters on the median",
        "population body weight MBW = 12.2 kg, given verbatim in the Table 3",
        "footnote d: CL_ind = CL * (BW/MBW)^0.75, Q_ind = Q * (BW/MBW)^0.75,",
        "VC_ind = VC * BW/MBW and VP_ind = VP * BW/MBW. Both exponents were",
        "FIXED, not estimated (Methods 'Covariate analysis': 'PWR the",
        "function power fixed to 0.75 for clearances and 1 for volumes of",
        "distribution'). The weight effect only became detectable after a",
        "sensitivity analysis removed one outlying patient with extremely",
        "low concentrations after the second and third doses, which had",
        "'masked the real impact of BW on clearances and volumes of",
        "distribution' (Results 'Covariate analysis').",
        "Table 1 rounds the model-building cohort median weight to 12 kg;",
        "the allometric reference of 12.2 kg is the unrounded value used by",
        "the authors and is the one encoded here."
      ),
      source_name = "BW"
    ),
    AGE = list(
      description = "Chronological age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Linear effect on the absorption rate constant Ka, centred on the",
        "median population age MAGE = 2.6 years, in the fractional-deviation",
        "form given in the Table 3 footnote: '(1 + theta_AGE_Ka",
        "(AGE-MAGE)/MAGE) with MAGE = 2.6 years'. (The printed footnote",
        "renders this as '(AGE-AGE)/MAGE'; the numerator's second term is",
        "MAGE, as in the identically-structured Table 2 footnote for the age",
        "effect on artesunate F1.) With theta = -0.67 a child of twice the",
        "median age has covariate value 1 and so Ka is reduced by 67%. The",
        "authors interpret this as a feeding effect: younger children are",
        "breastfed and so receive food closer to dosing, and food is",
        "reported to increase mefloquine Ka (Discussion).",
        "Multivariate analysis showed age accounted for an apparent effect",
        "of sex on Ka, so sex was not retained."
      ),
      source_name = "AGE"
    ),
    DAY2 = list(
      description = "Day-after-first-dose landmark indicator for the 3-day ASMQ course",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the first treatment day)",
      notes = paste(
        "1 = the record falls on the second or third treatment day, 0 = the",
        "first treatment day. Time-varying within subject. Selects between",
        "the two separately estimated absorption rate constants reported in",
        "Table 3 as 'Ka (h-1) DAY = 1' (0.17 /h) and 'Ka (h-1) DAY > 1'",
        "(0.40 /h); note the paper's Table 3 indexes the treatment day from",
        "1 whereas its Methods text indexes from 0, but both describe the",
        "same contrast of the first dosing day against the two that follow.",
        "The authors attribute the 2.4-fold increase to the dramatic drop in",
        "parasite load after the first artesunate dose improving the",
        "patient's state of health and resolving gastrointestinal",
        "disturbance (Discussion)."
      ),
      source_name = "DAY"
    )
  )

  # Screened in the covariate analysis but NOT retained in the final model.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Female children had 55% lower Ka than male children in univariate",
        "analysis (dOFV <= -7.0, p < 0.05), but 'multivariate analysis",
        "showed that age accounted for the effect of sex on Ka', so sex was",
        "not retained (Results 'Covariate analysis')."
      )
    ),
    PARA = list(
      description = "Plasmodium falciparum parasitaemia (asexual parasites/uL)",
      units = "parasites/uL",
      type = "continuous",
      notes = paste(
        "Tested on mefloquine relative bioavailability as a log10-transformed",
        "per-day baseline count; the univariate analyses 'showed no",
        "association between the covariates tested and MQ bioavailability'."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 48L,
    n_subjects_external_validation = 378L,
    n_studies = 1L,
    n_observations = list(model_building = 216L, external_validation = 538L),
    age_range = "0.6-5.0 years (median 2.6)",
    weight_range = "7-17 kg (median 12; allometric reference 12.2)",
    weight_median = "12.2 kg (MBW, the allometric reference in Table 3 footnote d)",
    age_median = "2.6 years (MAGE, the centring value for the age effect on Ka)",
    sex_female_pct = 60,
    race_ethnicity = "Black African (Kenyan, Tanzanian and Burkinabe)",
    disease_state = "Uncomplicated Plasmodium falciparum malaria (parasite density 2000-200,000 asexual parasites/uL, fever >= 37.5 C)",
    dose_range = "One or two dispersible tablets of 25 mg artesunate / 55 mg mefloquine once daily for 3 days: one tablet for ages 6-11 months, two tablets for ages 12-59 months",
    regions = "Kenya (model-building); Kenya, Tanzania and Burkina Faso, six centres (external validation)",
    notes = paste(
      "Demographics from Guidi 2019 Table 1. The structural model was",
      "developed on the 48 intensively sampled Kenyan children of the",
      "'Model-building dataset (n = 48)' column (216 MQ concentrations, none",
      "below the 2.5 ng/mL limit of quantification, median 5 samples per",
      "child). It was then externally validated on 538 sparse concentrations",
      "from the 378 children of the 'MQ validation dataset (n = 378)'",
      "column, giving a negligible individual-level bias of 0%",
      "(CI95 -2 to 1%) with 16% precision; per-site biases were",
      "non-significant or small (absolute value <= 6%, Table 5). The",
      "pharmacokinetic-pharmacodynamic dataset (n = 451) is a third,",
      "larger cohort used only for the exploratory logistic regression of",
      "recrudescence against predicted AUC, which found no association and",
      "is not part of this model (see the vignette)."
    )
  )

  ini({
    # ---- Relative bioavailability -------------------------------------
    # F1 fixed to 1 with an estimated BSV (Table 3, row 'F1': '1 FIX').
    # Its inclusion 'significantly decreased the OFV whilst explaining all
    # the BSV associated to VC' (Results), which is why no eta sits on vc.
    lfdepot <- fixed(log(1))
    label("Relative bioavailability of mefloquine, F1 (fraction)") # Guidi 2019 Table 3: F1 = 1 FIX

    # ---- Disposition ---------------------------------------------------
    # Typical values at the allometric reference weight MBW = 12.2 kg
    # (Table 3 footnote d).
    lcl <- log(0.45)
    label("Apparent mefloquine clearance at 12.2 kg, CL/F (L/h)") # Guidi 2019 Table 3: CL = 0.45 L/h (RSE 7%)
    lvc <- log(95)
    label("Apparent mefloquine central volume of distribution at 12.2 kg, VC/F (L)") # Guidi 2019 Table 3: VC = 95 L (RSE 7%)
    lq <- log(0.35)
    label("Apparent mefloquine inter-compartmental clearance at 12.2 kg, Q/F (L/h)") # Guidi 2019 Table 3: Q = 0.35 L/h (RSE 28%)
    lvp <- log(60)
    label("Apparent mefloquine peripheral volume of distribution at 12.2 kg, VP/F (L)") # Guidi 2019 Table 3: VP = 60 L (RSE 9%)

    # ---- Absorption ----------------------------------------------------
    # The authors estimated Ka SEPARATELY in each of two treatment-day
    # strata within the single joint fit rather than estimating a
    # reference value plus a multiplicative offset, so both strata carry an
    # explicit suffix and neither keeps the bare canonical `lka` (see
    # references/parameter-names.md, 'Stratum-suffixed parameters').
    lka_day1 <- log(0.17)
    label("Mefloquine absorption rate constant on the first treatment day, Ka (1/h)") # Guidi 2019 Table 3: Ka DAY = 1 -> 0.17 /h (RSE 17%)
    lka_day2 <- log(0.40)
    label("Mefloquine absorption rate constant on treatment days 2-3, Ka (1/h)") # Guidi 2019 Table 3: Ka DAY > 1 -> 0.40 /h (RSE 22%)

    # Age effect on Ka. Table 3 footnote form:
    #   Ka = Ka_day * (1 + theta_AGE_Ka * (AGE - MAGE) / MAGE), MAGE = 2.6 y
    # Table 3 gives -0.67 (a 67% fall at double the median age). The
    # Results text quotes 74%; that figure comes from the univariate step
    # before multivariate refinement, and the final Table 3 estimate is
    # the one encoded here.
    e_age_ka <- -0.67
    label("Effect of age on mefloquine Ka (fractional deviation from 2.6 years)") # Guidi 2019 Table 3: theta_AGE Ka = -0.67 (RSE 18%)

    # Allometric exponents, both FIXED by the authors rather than estimated
    # (Methods 'Covariate analysis', citing Holford 1996 [19]). One shared
    # exponent covers CL and Q, another covers VC and VP, matching the four
    # equations in the Table 3 footnote d.
    e_wt_cl_q <- fixed(0.75)
    label("Allometric exponent shared by mefloquine CL and Q (unitless)") # Guidi 2019 Table 3 footnote d; Methods 'Covariate analysis'
    e_wt_vc_vp <- fixed(1)
    label("Allometric exponent shared by mefloquine VC and VP (unitless)") # Guidi 2019 Table 3 footnote d; Methods 'Covariate analysis'

    # ---- Between-subject variability ----------------------------------
    # Methods: 'Exponential errors were assumed to capture BSV in all the
    # pharmacokinetic parameters'; Table 3 reports BSV as a percentage.
    # Converted with omega^2 = log(CV^2 + 1):
    #   F1  CV 28% -> log(1 + 0.28^2) = 0.0754785
    #   CL  CV 39% -> log(1 + 0.39^2) = 0.1415864
    #   Ka  CV 91% -> log(1 + 0.91^2) = 0.6032772
    # This CV reading is confirmed quantitatively by Table 4: propagating
    # the F1 and CL variances through AUC = F * Dose / CL gives a 95%
    # prediction-interval upper bound of 650 * exp(1.96 * sqrt(0.0754785 +
    # 0.1415864)) = 1620 mg/L*h against the published 1619. Reading the
    # same percentages as omega directly (sqrt(omega^2) = CV) would give
    # 1666 and is therefore excluded. See the vignette.
    # Table 3 shows BSV on F1, CL and Ka only: 'No improvement of the model
    # fit was observed associating BSV on Q or VP' and the F1 BSV absorbed
    # the VC variability.
    etalfdepot ~ 0.0754785 # Guidi 2019 Table 3: BSV on F1 = 28% (RSE 15%), bootstrap 27% (CI95 18 to 34)
    etalcl ~ 0.1415864 # Guidi 2019 Table 3: BSV on CL = 39% (RSE 17%), bootstrap 38% (CI95 24 to 51)
    etalka ~ 0.6032772 # Guidi 2019 Table 3: BSV on Ka = 91% (RSE 12%), bootstrap 87% (CI95 64 to 110)

    # ---- Residual variability -----------------------------------------
    # 'Finally, a proportional model was retained to describe the
    # intra-patient variability' (Results 'Structural and statistical
    # model'); Table 3 reports it as a CV%.
    propSd <- 0.21
    label("Proportional residual SD for mefloquine (fraction)") # Guidi 2019 Table 3: sigma_prop = 21% CV (RSE 13%)
  })

  model({
    # 1. Absorption rate constant: pick the treatment-day stratum, then
    #    apply the age effect and the between-subject variability. The
    #    (1 - DAY2) / DAY2 weighting selects exactly one of the two
    #    published estimates on each record.
    lka_day <- (1 - DAY2) * lka_day1 + DAY2 * lka_day2
    ka <- exp(lka_day + etalka) * (1 + e_age_ka * (AGE - 2.6) / 2.6)

    # 2. Individual parameters, allometrically scaled on MBW = 12.2 kg.
    fdepot <- exp(lfdepot + etalfdepot)
    cl <- exp(lcl + etalcl) * (WT / 12.2)^e_wt_cl_q
    vc <- exp(lvc) * (WT / 12.2)^e_wt_vc_vp
    q <- exp(lq) * (WT / 12.2)^e_wt_cl_q
    vp <- exp(lvp) * (WT / 12.2)^e_wt_vc_vp

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Bioavailability
    f(depot) <- fdepot

    # 6. Observation. mg/L; multiply by 1000 for the paper's ng/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
