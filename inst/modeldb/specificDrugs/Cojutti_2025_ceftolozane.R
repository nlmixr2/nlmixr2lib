Cojutti_2025_ceftolozane <- function() {
  description <- "One-compartment IV population PK model for ceftolozane given as ceftolozane/tazobactam continuous infusion in adults with suspected or documented Pseudomonas aeruginosa infection (retrospective therapeutic drug monitoring cohort, Bologna, Italy). Clearance is a power function of the 2021 race-free CKD-EPI eGFR de-indexed to the individual Du Bois body surface area (mL/min); the central volume is fixed at 32.1 L from the registration population PK model. Supports continuous-infusion probability-of-target-attainment analysis against the aggressive fCss/MIC >= 4 target."
  reference <- "Cojutti PG, Pai MP, Gatti M, Rinaldi M, Tonetti T, Siniscalchi A, Viale P, Pea F. Targeting Pseudomonas aeruginosa resistance at the exposure frontier: a population PK/PD blueprint for ceftolozane/tazobactam continuous infusion. Antimicrob Agents Chemother. 2025;69(11):e01042-25. doi:10.1128/aac.01042-25"
  vignette <- "Cojutti_2025_ceftolozane"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Doses are ceftolozane amounts in mg (the ceftolozane component of the 2:1
  # ceftolozane/tazobactam product; e.g. a '6 g/3 g daily' regimen is 6000 mg
  # ceftolozane per 24 h). Volumes in L, so central / vc is mg/L.
  compartmentData <- list(
    central = list(analyte = "ceftolozane", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "2021 race-free CKD-EPI creatinine eGFR de-indexed to the individual Du Bois body surface area (eGFR x BSA_DuBois / 1.73); raw mL/min, NOT BSA-normalized",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source column ckdepi2021mlmin_bsa_dd. Supplement Table S3 model 2 (bold, the selected",
        "covariate model, AIC 2711.36) uses 'beta_ckdepi2021_bsa_dd'; Supplement Table S4 defines",
        "ckdepi2021mlmin_bsa_dd as the 2021 CKD-EPI eGFR 'normalized to individual body surface",
        "area - Dubois Dubois equation', unit mL/min, computed in the supplement STATA code as",
        "CKDEPI2021 * bsa_dd / 1.73 with bsa_dd = ht^0.725 * 0.007184 * wt^0.425 (ht in cm, wt in",
        "kg) and serum creatinine in mg/dL. The main text calls this 'the 2021 non-race-based",
        "CKD-EPI equation with body surface area (BSA) indexation'; the operative quantity is the",
        "de-indexed mL/min value, which is also the unit the Monte Carlo renal classes are stated",
        "in. Table 1 reports only the BSA-indexed CKD-EPI eGFR, median (IQR) 75.1 (32.5-105.7)",
        "mL/min/1.73 m^2, range 4.5-165.8. Reference value 60 mL/min (Supplement Table S3",
        "footnote a 'eGFR/60'); see the e_crcl_cl comment in ini() for why 60 rather than the 100",
        "printed under Table 2. Patients on renal replacement therapy were excluded."
      ),
      source_name = "ckdepi2021mlmin_bsa_dd"
    )
  )

  # Supplement Table S2 univariate screen and Table S3 candidate models: the
  # other covariates tested on CLTOL but not retained in the final model.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened (Supplement Table S2, coefficient -0.1, P = 0.514) but not retained. Cohort median (IQR) 65.5 (58.0-74.3) years, range 19-90 (Table 1)."
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened together with ideal, adjusted, lean and dosing body weight, BMI and BSA (Supplement Table S2, all P > 0.08) but not retained. Cohort median (IQR) 72.5 (64.8-81.5) kg, range 37-119 (Table 1)."
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Supplement Table S2, P = 0.077) but not retained. 40/88 female (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 88L,
    n_studies = 1L,
    n_observations = 171L,
    age_range = "Adults; median (IQR) 65.5 (58.0-74.3) years, range 19-90 (Table 1, Results)",
    weight_range = "Median (IQR) 72.5 (64.8-81.5) kg, range 37-119 (Table 1, Results)",
    bmi_range = "Median (IQR) 25.5 (23.2-28.6) kg/m^2 (Table 1)",
    sex_female_pct = 45.5,
    race_ethnicity = "Not reported (single Italian tertiary university hospital)",
    disease_state = paste(
      "Adults receiving continuous-infusion ceftolozane/tazobactam for suspected or documented",
      "Pseudomonas aeruginosa infection: hospital- or ventilator-acquired pneumonia 32 (36.4%),",
      "bloodstream infection 18 (20.5%), skin and soft tissue 8 (9.1%), urinary tract 6 (6.8%),",
      "bone and joint 4 (4.5%), malignant otitis externa 3 (3.4%), intra-abdominal 2 (2.3%),",
      "infective endocarditis 1 (1.1%); empirical treatment 14 (15.9%). Mechanical ventilation 39",
      "(44.3%), vasopressors 27 (30.7%); median (IQR) SOFA score 7 (3-9) (Table 1, Results)."
    ),
    renal_function = "BSA-indexed 2021 CKD-EPI eGFR median (IQR) 75.1 (32.5-105.7) mL/min/1.73 m^2, range 4.5-165.8; augmented renal clearance 1 (1.1%); renal replacement therapy excluded (Table 1, Methods)",
    dose_range = paste(
      "2 g/1 g ceftolozane/tazobactam loading dose over 1 h, then continuous infusion initially by",
      "eGFR (6 g/3 g daily if eGFR > 50, 3 g/1.5 g daily if 30-50, 1.5 g/0.75 g daily if 15-29, 0.9",
      "g/0.45 g daily if < 15 mL/min/1.73 m^2) and subsequently TDM-adjusted. Median (range) daily",
      "dose 4 g/2 g (2 g/1 g - 6 g/3 g); treatment length median (IQR) 10 (8-14.3) days (Methods,",
      "Table 1)."
    ),
    regions = "Italy (IRCCS Azienda Ospedaliero-Universitaria di Bologna)",
    notes = paste(
      "Retrospective TDM cohort, November 2023 to October 2024. 171 steady-state total ceftolozane",
      "concentrations, median (range) 42.9 (2.1-216.8) mg/L; tazobactam Css median (IQR) 8.5",
      "(4.9-15.8) mg/L was measured but not modelled. LC-MS/MS assay, LLOQ 0.2 mg/L ceftolozane.",
      "Fitted with SAEM in Monolix 2024R1."
    )
  )

  ini({
    # Volume of the single compartment. Methods: 'Volumes of distribution of
    # ceftolozane ... was fixed based on a prior population pharmacokinetic
    # model used to establish the dosing recommendations ... namely 32.1 L'
    # (the Merck registration model, Supplement Figure S6). Continuous-infusion
    # steady-state data carry almost no information on V (eta shrinkage 75.9%).
    lvc <- fixed(log(32.1)); label("Central volume V1 (L)") # Table 2: V1 = 32.1 L, RSE '-' (fixed per Methods)

    # Typical clearance at the reference eGFR of 60 mL/min.
    lcl <- log(3.08); label("Clearance CLTOL at CRCL = 60 mL/min (L/h)") # Table 2: CLTOL = 3.08 L/h (RSE 7.2%)

    # Power exponent on eGFR. REFERENCE VALUE: Table 2 footnote and the
    # Discussion both print 'CLTOL = 3.08 x (eGFR/100)^0.71', but every number
    # the paper derives from the model requires eGFR/60:
    #  - Supplement Table S3 footnote a: the selected model 2 is '(eGFR/60)'.
    #  - Discussion: central-tendency CL 2.71 L/h at eGFR 50 (the value Zhang
    #    2021 gives 2.98 L/h at). 3.08 x (50/60)^0.71 = 2.71; (50/100) gives 1.88.
    #  - Discussion: 5.90 L/h at 150 mL/min. 3.08 x (150/60)^0.71 = 5.90;
    #    (150/100) gives 4.11.
    #  - Results: median individual CL 3.69 L/h, vs 3.08 x (~80/60)^0.71 ~ 3.8
    #    at the cohort median de-indexed eGFR ((80/100) gives 2.6).
    #  - Figure 1 PTA at MIC 8 mg/L (fCss/MIC >= 4) is reproduced with /60 and
    #    over-predicted by 15-25 points with /100 (see the vignette).
    # The parameter-table side wins over the display equation; the vignette
    # Assumptions and deviations section shows the /100 sensitivity.
    e_crcl_cl <- 0.71; label("Power exponent of (CRCL/60) on CL (unitless)") # Table 2: beta_CLTOL = 0.71 (RSE 13.9%)

    # IIV. Monolix reports omega as the SD of the log-normal random effect with
    # the CV in parentheses. Variances below use omega^2 = log(1 + CV^2) from
    # the four-significant-figure CVs (0.58^2 = 0.3364 and 0.42^2 = 0.1764 from
    # the rounded SDs agree).
    etalcl ~ log(1 + 0.6265^2) # Table 2: omega CLTOL = 0.58 (CV 62.65%), RSE 10.3%
    etalvc ~ fixed(log(1 + 0.4392^2)) # Table 2: omega V1 = 0.42 (CV 43.92%), RSE '-' (held at its value together with V1; no SE in Table 2 or Table S1)

    # Residual error. Supplement Table S3: the final structure uses Monolix
    # 'Combined2' error, sd = sqrt(a^2 + (b * f)^2), which is nlmixr2's
    # add() + prop() with combined2().
    addSd <- 12.36; label("Additive residual error a (mg/L)") # Table 2: a = 12.36 (RSE 9.4%)
    propSd <- 0.21; label("Proportional residual error b (fraction)") # Table 2: b = 0.21 (RSE 15.5%)
  })

  model({
    # Individual parameters. Clearance is a power function of the
    # BSA-de-indexed CKD-EPI 2021 eGFR centred at 60 mL/min.
    cl <- exp(lcl + etalcl) * (CRCL / 60)^e_crcl_cl
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    # One-compartment disposition with zero-order input supplied through the
    # event table (rate or dur on the dosing records): a 1 h loading infusion
    # followed by a continuous infusion into the central compartment.
    d/dt(central) <- -kel * central

    # Total ceftolozane plasma concentration (mg/L). The paper's PK/PD target
    # uses the free concentration 0.81 x Css (19% protein binding, Methods,
    # Monte Carlo simulation); that conversion is applied in the vignette.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined2()
  })
}
