Sanches_2022_piperacillin <- function() {
  description <- "Two-compartment IV-infusion population PK model for piperacillin in critically ill Brazilian adults, with clearance proportional to Cockcroft-Gault creatinine clearance normalised to 60 mL/min/1.73 m^2 and distribution written as the central-to-peripheral and peripheral-to-central rate constants KCP and KPC. Estimated with the Pmetrics non-parametric adaptive grid (NPAG) algorithm; the NPAG marginal means are the typical values and the reported %CV values are carried as log-normal between-subject variability. The fitted model also estimated per-subject initial conditions for the day-5 sampling interval, which are not reported, so a simulation from treatment start reaches concentrations well above the paper's observed day-5 data (see the vignette Assumptions and deviations) (Sanches 2022)"
  reference <- "Sanches C, Alves GCS, Farkas A, da Silva SD, de Castro WV, Chequer FMD, Beraldi-Magalhaes F, Magalhaes IRdS, Baldoni AdO, Chatfield MD, Lipman J, Roberts JA, Parker SL. Population Pharmacokinetic Model of Piperacillin in Critically Ill Patients and Describing Interethnic Variation Using External Validation. Antibiotics (Basel). 2022;11(4):434. doi:10.3390/antibiotics11040434"
  vignette <- "Sanches_2022_piperacillin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Sanches 2022 section 4.2: total piperacillin was assayed in plasma by
  # HPLC-UV; doses and amounts are in mg and volumes in L, giving mg/L.
  compartmentData <- list(
    central = list(analyte = "piperacillin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "piperacillin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault equation, normalised to 1.73 m^2 body surface area",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Sanches 2022 section 2 (Results): 'The inclusion of creatinine clearance (CRCL) normalized to 60 mL/min/1.73 m 2 as covariate ... by the equation CL = TVCL *(CRCL/60)'. Linear, proportional, with no intercept and no estimated exponent. Section 3 (Discussion) and section 4.5 identify the estimate as Cockcroft-Gault. Cohort median 60 mL/min/1.73 m^2, IQR 47-83 (Table 1). Time-fixed per subject.",
      source_name = "CRCL"
    )
  )

  # Screened in the covariate analysis (section 4.3 step iii) but not retained
  # in the final model. SAPS 3, MODS and the CKD-EPI eGFR were also screened;
  # they have no register entry and are listed only in the vignette.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (section 4.3 step iii), not retained. Cohort median 72 years, IQR 57-78 (Table 1)."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (section 4.3 step iii), not retained. Not summarised in Table 1."
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (section 4.3 step iii), not retained. Cohort median 69 kg, IQR 57-77 (Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened as 'sex' (section 4.3 step iii), not retained. 9 of 24 patients male (Table 1)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (section 4.3 step iii), not retained. Cohort median 22 kg/m^2, IQR 21-31 (Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened as 'creatinine' (section 4.3 step iii), not retained. Patients with serum creatinine > 2 mg/dL were excluded (section 4.1)."
    ),
    DIS_SEPSIS = list(
      description = "Sepsis indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no sepsis)",
      notes = "Screened as 'presence of sepsis' (section 4.3 step iii), not retained. 12 of 24 patients (Table 1)."
    ),
    SOFA = list(
      description = "Sequential Organ Failure Assessment score at the time of sampling",
      units = "points",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (section 4.3 step iii), not retained. Cohort median 5, IQR 4-7 (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 24L,
    n_studies = 1L,
    age_range = "IQR 57-78 years",
    age_median = "72 years",
    weight_range = "IQR 57-77 kg",
    weight_median = "69 kg",
    sex_female_pct = 62.5,
    race_ethnicity = "Brazilian (admixed European, African and Amerindian ancestry; ethnicity not recorded per patient)",
    disease_state = "Critically ill ICU adults with confirmed or suspected infection treated with piperacillin/tazobactam. Sepsis 50%, microbiologically confirmed infection 58%, vasoactive drugs 29%, ICU mortality 33%. SAPS 3 median 53 (IQR 45-63), SOFA 5 (4-7), MODS 3 (2-4). Patients with serum creatinine > 2 mg/dL or a rise above twice baseline were excluded.",
    dose_range = "Piperacillin 4 g every 8 h as an intermittent infusion over about 30 min (empiric arm), or 2 g every 6 h or 3.3 g every 4 h (individually designed optimum dosing strategy arm), from a randomised trial.",
    regions = "Brazil (ICU of a medium-sized hospital in the Midwest region of Minas Gerais)",
    renal_function = "Cockcroft-Gault creatinine clearance median 60 mL/min/1.73 m^2 (IQR 47-83)",
    notes = "Sanches 2022 Table 1 and sections 2, 4.1-4.3. Plasma sampled on day 5 of treatment, one or two samples per patient within a dosing interval (up to three: predose and 1 and 3 h after the start of infusion). Total piperacillin by HPLC-UV, range 2.5-100 mg/L; no concentration exceeded 100 mg/L. Fitted with Pmetrics 1.5.0 NPAG. External validation (not used for fitting): 20 Australian ICU patients (Udy 2015) and 10 Indigenous Australian ICU patients (Tsai 2016)."
  )

  ini({
    # Structural parameters: Sanches 2022 Table 2, final covariate model.
    # Pmetrics NPAG reports the mean (SD), median and %CV of each marginal of
    # the non-parametric joint density. The MEANS are used as typical values:
    # the Abstract and Results quote them ('Clearance and volume of
    # distribution were (mean +/- SD) 3.33 +/- 1.24 L h-1 and 10.69 +/- 4.50
    # L'), and the %CV column is SD/mean (1.24/3.33 = 37%, 4.50/10.69 = 42%,
    # 0.15/1.15 = 13%). Neither column reproduces the paper's Figure 2 PTA
    # (see the vignette Assumptions and deviations).
    lcl <- log(3.33)
    label("Clearance at CRCL = 60 mL/min/1.73 m^2 (L/h)") # Table 2, 'CL (L/h)' mean 3.33 (SD 1.24), median 3.01, CV 37%
    lvc <- log(10.69)
    label("Central volume of distribution (L)") # Table 2, 'V (L)' mean 10.69 (SD 4.50), median 9.03, CV 42%
    lk12 <- log(1.15)
    label("Central-to-peripheral rate constant KCP (1/h)") # Table 2, 'KCP (h-1)' mean 1.15 (SD 0.15), median 1.21, CV 13%
    lk21 <- log(0.08)
    label("Peripheral-to-central rate constant KPC (1/h)") # Table 2, 'KPC (h-1)' mean 0.08 (SD 0.09), median 0.03, CV 120%

    # Between-subject variability. NPAG estimates a discrete joint density,
    # not an omega matrix. Each marginal is approximated by a log-normal with
    # the reported %CV, omega^2 = log(1 + CV^2); no correlations are reported,
    # so the etas are independent.
    etalcl ~ 0.128393 # Table 2, CL CV 37% -> log(1 + 0.37^2)
    etalvc ~ 0.162472 # Table 2, V CV 42% -> log(1 + 0.42^2)
    etalk12 ~ 0.016759 # Table 2, KCP CV 13% -> log(1 + 0.13^2)
    etalk21 ~ 0.891998 # Table 2, KPC CV 120% -> log(1 + 1.20^2); the rounded mean and SD give 0.09/0.08 = 113%

    # Residual error: the Pmetrics gamma model, error = SD * gamma, with the
    # assay SD polynomial SD = C0 + C1 * Y (section 4.3 step ii). Results
    # section 2: 'Residual error was modelled as gamma * (1 + 0.1*concentration),
    # value = 5', i.e. C0 = 1 mg/L and C1 = 0.1 held fixed, gamma = 5 fitted.
    addSd <- fixed(1)
    label("Assay SD polynomial intercept C0 (mg/L)") # Results section 2, 'gamma * (1 + 0.1*concentration)'
    propSd <- fixed(0.1)
    label("Assay SD polynomial linear coefficient C1 (fraction)") # Results section 2, 'gamma * (1 + 0.1*concentration)'
    gammaSd <- 5
    label("Pmetrics gamma multiplier on the assay SD (unitless)") # Results section 2, 'value = 5'
  })
  model({
    # Individual parameters. Results section 2: CL = TVCL * (CRCL / 60).
    cl <- exp(lcl + etalcl) * (CRCL / 60)
    vc <- exp(lvc + etalvc)
    kcp <- exp(lk12 + etalk12)
    kpc <- exp(lk21 + etalk21)

    # The paper parameterises distribution with KCP and KPC. The ODEs are
    # driven through the exactly equivalent q = KCP * V and vp = q / KPC: a
    # system written straight from stored micro-constants can be rewritten by
    # rxSolve()'s default linear-compartment conversion into a one-compartment
    # solve that silently drops peripheral1.
    q <- kcp * vc
    vp <- q / kpc

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Total plasma piperacillin (mg/L).
    Cc <- central / vc

    # Pmetrics evaluates the SD polynomial on the observed concentration;
    # nlmixr2 evaluates it on the prediction.
    sdCc <- gammaSd * (addSd + propSd * Cc)
    Cc ~ add(sdCc)
  })
}
