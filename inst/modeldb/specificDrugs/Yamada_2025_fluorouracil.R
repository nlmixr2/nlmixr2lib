Yamada_2025_fluorouracil <- function() {
  description <- "One-compartment population PK model of fluorouracil (5-FU) with zero-order IV input and first-order elimination in patients with locally advanced unresectable or metastatic HER2-negative, CLDN18.2-positive gastric or gastroesophageal junction (G/GEJ) adenocarcinoma receiving mFOLFOX6 with or without concomitant zolbetuximab (Yamada 2025, ILUSTRO Cohort 2). Zolbetuximab coadministration was tested as a covariate on CL and on V and was not statistically significant in either case, so the best base model is also the final model and the extracted model carries no covariate effects."
  reference <- paste(
    "Yamada A, Yang J, Bonate PL, Heo N, Poondru S. Population",
    "pharmacokinetic analysis of fluorouracil and oxaliplatin in the",
    "absence or presence of zolbetuximab in locally advanced unresectable",
    "or metastatic gastric or gastroesophageal junction adenocarcinoma.",
    "Cancer Chemother Pharmacol. 2025;95:89.",
    "doi:10.1007/s00280-025-04808-2 (PMCID PMC12449379).",
    "Parameter estimates are from Supplementary Table 1 ('Parameter",
    "estimates of best base model for 5-FU and models with zolbetuximab",
    "effect on 5-FU PK'), 'Best base model' column, in the online",
    "Supplementary Material (280_2025_4808_MOESM1_ESM.pdf).",
    "Companion oxaliplatin model from the same paper:",
    "modellib('Yamada_2025_oxaliplatin').",
    sep = " "
  )
  vignette <- "Yamada_2025_fluorouracil_oxaliplatin_zolbetuximab"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. The paper reports CL in L/h and V in L, so with dose
  # in mg the natural concentration unit is mg/L = ug/mL. The paper's own
  # tables are printed in ng/mL (a factor of 1000 larger); see the vignette
  # source-trace table.
  compartmentData <- list(
    central = list(analyte = "fluorouracil", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    CONMED_ZOLBETUXIMAB = list(
      description        = "Concomitant zolbetuximab coadministration (1 = 5-FU given in the presence of zolbetuximab, 0 = 5-FU alone)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (5-FU alone; ILUSTRO Cohort 2 Cycle 1, in which zolbetuximab was deliberately delayed to Day 3 so that 5-FU PK could be sampled without it)",
      notes              = "Screened as a covariate on CL and on V using the fractional model P = theta_p * (1 + theta_zolbetuximab * X_zolbetuximab) (Methods, 'Evaluation of zolbetuximab effect as a covariate'; Article Equation a). Neither effect reached the forward-selection criterion: dOFV was -0.787 for the effect on CL and -2.279 for the effect on V, both well short of the 10.83 required at alpha = 0.001 with 1 degree of freedom (Supplementary Table 1; Results, 'Evaluation of zolbetuximab effect on 5-FU PK'). The rejected point estimates were -0.0531 (RSE 110.5%) on CL and 0.248 (RSE 74.6%) on V. Retained here for provenance only; the final 5-FU model has no covariates. The effect IS retained on oxaliplatin V1-V3 -- see modellib('Yamada_2025_oxaliplatin').",
      source_name        = "X_zolbetuximab"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 21L,
    n_studies      = 1L,
    disease_state  = "Adults with previously untreated, locally advanced unresectable or metastatic HER2-negative, CLDN18.2-positive gastric or gastroesophageal junction (G/GEJ) adenocarcinoma (CLDN18.2 positivity defined as >= 75% of tumor cells with moderate to strong membranous CLDN18 staining intensity). Cohort 2 of the phase 2 ILUSTRO study (NCT03505320).",
    dose_range     = "mFOLFOX6 on Days 1, 15 and 29 of 42-day cycles for up to 4 cycles: 5-FU 400 mg/m^2 IV bolus over 5-15 min, followed by 5-FU 2400 mg/m^2 continuous IV infusion over 46-48 h (given with oxaliplatin 85 mg/m^2 and folinic acid 400 mg/m^2 infused over 2 h). Zolbetuximab 800 mg/m^2 loading dose then 600 mg/m^2, each a minimum 2-h IV infusion, on Days 1 and 22 of each cycle except Cycle 1, where it was given on Day 3 so that 5-FU PK could be sampled in its absence.",
    notes          = "Baseline clinical and demographic characteristics are not tabulated in this paper; it cites the previously published ILUSTRO Cohort 2 report (reference 4) for them, which was not available on disk at extraction time. Analysis dataset: of the original 206 non-BQL postdose 5-FU concentrations (108 without zolbetuximab, 98 with zolbetuximab), 12 non-BQL postdose concentrations, 2 non-BQL predose concentrations and 35 BQL postdose concentrations were excluded during base model development for implausible PK profiles or dosing records. PK sampling: predose and 0.5, 1, 2, 5, 24 and 48 h after the start of dosing in Cycles 1 (chemotherapy alone) and 2 (chemotherapy plus zolbetuximab). Assay: validated LC-MS/MS, LLOQ 5.5 ng/mL. Estimation in NONMEM 7.5.0; condition number of the final model 2.32."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Supplementary Table 1, 'Best base model
    # parameter estimate with (%RSE) or [shrinkage]' column. This column
    # IS the final model: Results, 'Evaluation of zolbetuximab effect on
    # 5-FU PK' concludes "Therefore, the best base model was the final
    # model for 5-FU."
    #
    # Consistency check (not a source value): CL/V = 179/65.6 = 2.729 /h
    # gives a terminal half-life of 0.254 h, matching the 0.25 h quoted
    # in the Discussion.
    # ------------------------------------------------------------------
    lcl <- log(179);  label("Systemic clearance CL (L/h)")                            # Supp Table 1, best base model: CL 179 L/h (RSE 8.7%)
    lvc <- log(65.6); label("Central volume of distribution V (L)")                   # Supp Table 1, best base model: V 65.6 L (RSE 12.7%)

    # ------------------------------------------------------------------
    # Inter-individual variability -- Supplementary Table 1 'IIV
    # [shrinkage]' rows, reported as CV%. Converted to a log-normal
    # variance with the library convention omega^2 = log(CV^2 + 1).
    # Reported shrinkage: 9.5% on CL, 32.9% on V.
    # ------------------------------------------------------------------
    etalcl ~ 0.111842                                                                 # Supp Table 1: omega_CL = 34.4 % CV -> log(1 + 0.344^2) = 0.111842
    etalvc ~ 0.161745                                                                 # Supp Table 1: omega_V  = 41.9 % CV -> log(1 + 0.419^2) = 0.161745

    # ------------------------------------------------------------------
    # Residual variability -- Supplementary Table 1 'Residual variability
    # (%RSE)' row, 'Proportional error [CV%]' = 0.583 (RSE 7.0%).
    #
    # The row header reads [CV%] but the value is printed WITHOUT a
    # percent sign, unlike every other percentage in the same table (the
    # IIV rows are printed as "34.4%" and "41.9%"). Read literally as a
    # percentage it would be 0.583% CV, which is implausible for these
    # data. It is therefore read as the CV expressed as a fraction, i.e.
    # 58.3% CV, which is consistent with the paper's own account of the
    # 5-FU residuals (skewed / heavy-tailed CWRES Q-Q plot,
    # underestimation of higher concentrations; Results, 'Predictive
    # performance of final model for 5-FU'). The competing reading -- the
    # value being a NONMEM variance, giving sqrt(0.583) = 76.4% CV -- is
    # documented in the vignette Errata.
    #
    # A mixed (additive + proportional) residual model had a lower OFV but
    # was rejected because the RSE on the IIV of V was 286%; a purely
    # proportional model was selected (Results, 'Base model for 5-FU').
    # ------------------------------------------------------------------
    propSd <- 0.583; label("Proportional residual error on 5-FU plasma concentration (fraction)")  # Supp Table 1: proportional error 0.583 (RSE 7.0%)
  })

  model({
    # 1. Individual parameters
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    # 2. Micro-constant
    kel <- cl / vc

    # 3. ODE system. Zero-order input: 5-FU is given as a short IV
    #    "bolus" (5-15 min) followed by a 46-48 h continuous infusion,
    #    both encoded as infusions on the dosing records (rate / dur),
    #    so no absorption compartment is required.
    d/dt(central) <- -kel * central

    # 4. Observation and error
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
