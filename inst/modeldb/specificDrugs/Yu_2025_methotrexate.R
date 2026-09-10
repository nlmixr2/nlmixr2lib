Yu_2025_methotrexate <- function() {
  description <- paste(
    "Two-compartment IV-infusion population PK model for high-dose",
    "methotrexate (3 or 5 g/m^2 over 24 h) in 214 Chinese children with",
    "acute lymphoblastic leukaemia (Yu 2025; 1,672 plasma concentrations).",
    "Clearance uses an age-cutoff structure at 1 year: a separate typical",
    "clearance is estimated in each age stratum (4.46 L/h for age > 1 year,",
    "1.69 L/h for age <= 1 year), with covariate effects shared across both",
    "strata -- a bedside-Schwartz eGFR power effect (exponent 0.537,",
    "reference 160 mL/min/1.73 m^2), a body-weight power effect (exponent",
    "0.45, reference 20 kg) and a blood-urea-nitrogen power effect (exponent",
    "-0.0823, reference 3 mmol/L). The paper's stated novelty is this pair:",
    "the 1-year clearance cutoff, and the use of BUN alongside eGFR as a",
    "second, non-collinear renal marker. Central volume scales with body",
    "surface area (15.9 L at 0.77 m^2, exponent 1.10); intercompartmental",
    "clearance and peripheral volume carry no covariates. Exponential",
    "between-subject variability on clearance, intercompartmental clearance",
    "and peripheral volume (none on central volume), with a 44.2%",
    "proportional residual error.",
    sep = " "
  )
  reference <- paste(
    "Yu B, Wan Y, Mei K, Zhan D, Tang Q, Hu X, Ji W, Cai H (2025).",
    "Population Pharmacokinetics and Covariate Analysis of Methotrexate in",
    "Pediatric Acute Lymphoblastic Leukemia.",
    "Drug Des Devel Ther 19:8473-8486. doi:10.2147/DDDT.S545368.",
    sep = " "
  )
  vignette <- "Yu_2025_methotrexate"
  units    <- list(time = "h", dosing = "umol", concentration = "umol/L")

  # Issue #482. Verified against Yu 2025 Methods "Population Pharmacokinetic
  # Model Development" (NONMEM ADVAN3 TRANS4, i.e. two-compartment IV) and the
  # final-model equations on p. 8480, which parameterise exactly CL, V1, Q and
  # V2. Amount units are umol because every methotrexate concentration in the
  # paper is reported in umol/L (LLOQ 0.17 umol/L; steady-state targets 26-60
  # umol/L for low-risk and 52-100 umol/L for intermediate/high-risk patients).
  compartmentData <- list(
    central     = list(analyte = "methotrexate", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "methotrexate", units = "umol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Yu 2025 Table 1 (Modeling Data column): median 5 years, range 0.65-14, mean 5.52. Used ONLY as the 1-year cutoff that selects between the two clearance strata of the two final-model equations on p. 8480 ('If age > 1 years old CL (L/h) = ...' / 'If age <= 1 years old CL (L/h) = ...'). It does not enter the model as a continuous term: the paper screened AGE continuously first and reports that the piecewise 1-year cutoff fit substantially better (dOFV -55.825 for the cutoff versus -12.535 for the continuous term; Discussion). The cutoff is the paper's stated novelty, motivated by immature hepatic and renal function below 1 year of age.",
      source_name        = "age"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate by the bedside Schwartz equation (BSA-normalized)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Yu 2025 Methods: 'The eGFR was calculated using the bedside Schwartz equation', i.e. 0.413 * height (cm) / Scr (mg/dL), reported normalized to 1.73 m^2. Table 1 (Modeling Data): median 160 mL/min/1.73 m^2, range 21.90-405.65, mean 162.40 -- a supranormal-skewed paediatric cohort, so the model carries little information about frank renal impairment even though the paper's own Monte Carlo simulations extrapolate down to 80 mL/min/1.73 m^2. The reference value 160 in the CL equations is the cohort median. Enters CL as (eGFR/160)^0.537, shared across both age strata; it was the single most influential covariate (dOFV -170.45). Stored under canonical CRCL, which explicitly admits the paediatric bedside-Schwartz variant (precedent Jung_2024_vancomycin.R); the estimating equation is documented here per the register's per-model requirement.",
      source_name        = "eGFR"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Yu 2025 Table 1 (Modeling Data): median 19.50 kg, range 7-62.50, mean 22.43. The reference weight 20 kg in the CL equations is the cohort median rounded by the authors. Enters CL as (WT/20)^0.45, shared across both age strata. Note the exponent is estimated (0.45, RSE 11%), not fixed at the allometric 0.75, and that body SIZE enters clearance through weight but volume through body surface area -- the paper models the two size descriptors on different parameters rather than using one throughout.",
      source_name        = "WT"
    ),
    BUN = list(
      description        = "Blood urea nitrogen",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Yu 2025 Table 1 (Modeling Data): median 2.95 mmol/L, range 0.40-16.70, mean 3.15. The reference value 3 in the CL equations is the cohort median rounded by the authors. Enters CL as (BUN/3)^(-0.0823), shared across both age strata; the NEGATIVE exponent means a higher urea burden lowers methotrexate clearance. Retaining BUN alongside eGFR is one of the paper's two stated novelties: the authors checked for collinearity between the two renal markers and found effectively none (Spearman rho = -0.1192, variance inflation factor 1.007457), and argue BUN adds tubular / volume-status information that glomerular eGFR misses. The effect is the weakest term in the model (RSE 24%, bootstrap 95% CI -0.1553 to -0.0093, the only interval that approaches zero), so it moves clearance by only about 6% across the interquartile range of BUN.",
      source_name        = "BUN"
    ),
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Yu 2025 Table 1 (Modeling Data): median 0.77 m^2, range 0.36-1.71. The reference value 0.77 in the V1 equation is the cohort median, and is also the BSA of the 5-year-old, 19 kg virtual child used for every Monte Carlo simulation in the paper. Enters the CENTRAL volume only, as (BSA/0.77)^1.10; the peripheral volume V2 and the intercompartmental clearance Q carry no covariates. BSA is separately the basis of the dose itself (3 g/m^2 for low-risk, 5 g/m^2 for intermediate/high-risk patients), so it acts on both sides of the exposure calculation.",
      source_name        = "BSA"
    )
  )

  # Covariates the paper screened but did not retain in the final model. No
  # point estimate was published for any of them, so none can be encoded.
  covariatesDataExcluded <- list(
    SEXF = list(
      description        = "Female sex indicator. Screened as a categorical covariate and rejected.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Yu 2025 Methods 'Covariate Model' lists sex among the categorical variables evaluated; it does not appear in the final model of Table 2. Table 1 (Modeling Data): 105 male / 66 female. No point estimate published."
    ),
    CREAT = list(
      description        = "Serum creatinine. Screened as a renal-function covariate and rejected in favour of the derived eGFR.",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Yu 2025 Methods lists Scr among the continuous variables evaluated; it does not appear in Table 2. Table 1 (Modeling Data): median 0.29 mg/dL, range 0.07-2.09, mean 0.32. Serum creatinine is an input to the bedside-Schwartz CRCL that IS retained, but does not enter the model separately. No point estimate published."
    ),
    ALT = list(
      description        = "Alanine aminotransferase. Screened as a hepatic-function covariate and rejected.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Yu 2025 Methods lists ALT among the continuous variables evaluated; it does not appear in Table 2. Table 1 (Modeling Data): median 29.50 U/L, range 3.30-1268.6. No point estimate published."
    ),
    AST = list(
      description        = "Aspartate aminotransferase. Screened as a hepatic-function covariate and rejected.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Yu 2025 Methods lists AST among the continuous variables evaluated; it does not appear in Table 2. Table 1 (Modeling Data): median 36.30 U/L, range 10.20-1258. No point estimate published."
    ),
    MTXNUM = list(
      description        = "Number of prior high-dose methotrexate chemotherapy cycles. Screened as a categorical covariate and rejected.",
      units              = "(count, treated as categorical by the source)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Yu 2025 Methods 'Covariate Model': 'The categorical variables were sex and the number of methotrexate chemotherapy cycles (MTXNUM)'. It does not appear in the final model of Table 2, and the paper never states how many levels MTXNUM took or how they were coded, so the column cannot be reconstructed even for a screening replication. No point estimate published."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 214L,
    n_studies        = 1L,
    age_range        = "0.65-14 years",
    age_median       = "5 years (mean 5.52) in the 171-patient model-building set; 6 years (mean 6.35, range 0.86-14) in the 43-patient external validation set",
    weight_range     = "7-62.50 kg",
    weight_median    = "19.50 kg (mean 22.43)",
    sex_female_pct   = 38.6,
    race_ethnicity   = "Not reported by the authors as a covariate. The Discussion describes the cohort as 'primarily Han Chinese children treated at our institution' and names this ethnic homogeneity as a limitation on external generalizability.",
    disease_state    = "Children (age <= 18 years) with pathologically confirmed acute lymphoblastic leukaemia receiving high-dose methotrexate consolidation, risk-stratified into low-risk (LR) and intermediate/high-risk (I/HR) groups under the CCLG-ALL-2018 protocol (before October 2021) or the CCCG-ALL-2020 protocol (from October 2021). Inclusion required at least one plasma methotrexate concentration.",
    dose_range       = "3 g/m^2 (low-risk) or 5 g/m^2 (intermediate/high-risk) as a 24-hour intravenous infusion, given as a loading-dose strategy: 10% of the total dose over 0.5 h followed by the remaining 90% over 23.5 h.",
    regions          = "China: Anhui Provincial Children's Hospital, Hefei, Anhui Province. Single centre, May 2021 to November 2024.",
    renal_function   = "Bedside-Schwartz eGFR median 160 mL/min/1.73 m^2 (range 21.90-405.65, mean 162.40) -- supranormal-skewed, as is usual in children. Serum creatinine median 0.29 mg/dL (range 0.07-2.09); blood urea nitrogen median 2.95 mmol/L (range 0.40-16.70).",
    n_concentrations = 1672L,
    notes            = "Baseline demographics from Yu 2025 Table 1. Of the 214 patients contributing 1,672 concentrations, 171 patients / 1,342 concentrations formed the model-building set and a randomly selected 20% (43 patients / 330 concentrations) formed the external validation set; the population fields above quote the model-building column except where stated. Retrospective, single-centre. Samples were assayed by enzyme-multiplied immunoassay (EMIT, Viva-ProE, Siemens) calibrated over 0.3-2600 umol/L, with an LLOQ of 0.17 umol/L; values below the LLOQ were excluded. Routine monitoring at 20-24 h, 44-48 h and 68-72 h after the start of infusion, continuing until methotrexate fell to <= 0.2 umol/L. Fit in NONMEM 7.4 with PsN 4.6.0 using FOCE-I. Internal evaluation by GOF plots, prediction-corrected VPC and a bootstrap (93.6% success rate, all estimates within 16% of the final model); external evaluation gave MPE -3.99%, MAPE 22.4%, F20 46.36% and F30 64.55%. The paper reports NO pharmacogenetic data (SLCO1B1, ABCC2, MTHFR), which the authors name as a limitation. The Monte Carlo dosing analysis (Tables 3-4, Figure 3) is a simulation from this same model, not a separate model."
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural parameters. Yu 2025 Table 2, "NONMEM Estimate" column; the
    # final-model equations are the five display equations on p. 8480:
    #
    #   if age >  1 year  CL (L/h) = 4.46 * (eGFR/160)^0.537 * (WT/20)^0.45 * (BUN/3)^-0.0823
    #   if age <= 1 year  CL (L/h) = 1.69 * (eGFR/160)^0.537 * (WT/20)^0.45 * (BUN/3)^-0.0823
    #                     V1 (L)   = 15.90 * (BSA/0.77)^1.10
    #                     Q  (L/h) = 0.149
    #                     V2 (L)   = 7.23
    #
    # NOTE ON SOURCING: those five equations are embedded in the PDF as
    # VECTOR GRAPHICS, not text. `pdftotext` and the preprocessed
    # `_trimmed.md` both emit the bare stubs "If age >1 years old CL (L/h) ="
    # with nothing after them, and the generic covariate forms (Eq. 1-2 on
    # p. 8477) are lost the same way. The equations above were read by
    # rendering p. 3 and p. 6 to PNG and reading them visually. Anyone
    # re-verifying this file from a text extraction will find the reference
    # values 160 / 20 / 3 / 0.77 MISSING; they are not in any table.
    #
    # This is ONE jointly-fit model: only the typical clearance is
    # stratum-specific, so it alone carries an explicit `_agele1` / `_agegt1`
    # suffix (symmetric stratum-suffix scheme; neither stratum is the
    # reference). Every covariate effect, every IIV term and the residual
    # error are shared, and so keep their bare canonical names.
    # ---------------------------------------------------------------------

    lcl_agegt1 <- log(4.46); label("Typical clearance for age > 1 year at eGFR = 160 mL/min/1.73 m^2, WT = 20 kg, BUN = 3 mmol/L (L/h)")   # Yu 2025 Table 2: theta(CL age>1) = 4.46 L/h (RSE 3%, bootstrap median 4.63, 95% CI 4.23-4.86)
    lcl_agele1 <- log(1.69); label("Typical clearance for age <= 1 year at eGFR = 160 mL/min/1.73 m^2, WT = 20 kg, BUN = 3 mmol/L (L/h)")  # Yu 2025 Table 2: theta(CL age<=1) = 1.69 L/h (RSE 10%, bootstrap median 1.76, 95% CI 0.65-2.73)
    lvc        <- log(15.90); label("Central volume of distribution at BSA = 0.77 m^2 (L)")                                                # Yu 2025 Table 2: theta(V1) = 15.90 L (RSE 4%, bootstrap median 16.01, 95% CI 13.84-17.89)
    lq         <- log(0.149); label("Intercompartmental clearance (L/h)")                                                                  # Yu 2025 Table 2: theta(Q) = 0.149 L/h (RSE 5%, bootstrap median 0.147, 95% CI 0.124-0.174)
    lvp        <- log(7.23); label("Peripheral volume of distribution (L)")                                                                # Yu 2025 Table 2: theta(V2) = 7.23 L (RSE 9%, bootstrap median 7.27, 95% CI 4.89-9.57)

    # ---------------------------------------------------------------------
    # Covariate effects. All are shared across the two age strata: the strata
    # differ only in the multiplicative typical clearance above. Each is a
    # power of a median-normalised ratio, per Eq. (1) on p. 8477
    # (Pi = theta_p * (COV/COV_m)^theta). The Methods call Eq. 1 "the
    # exponential model", but the printed equation is unambiguously a POWER
    # model, and so are all five final-model equations; the equations win.
    # ---------------------------------------------------------------------

    e_crcl_cl <- 0.537;   label("Power exponent on (CRCL/160) for CL, both age strata (unitless)")  # Yu 2025 Table 2: theta(GFR) = 0.537 (RSE 6%, bootstrap median 0.512, 95% CI 0.354-0.721)
    e_wt_cl   <- 0.45;    label("Power exponent on (WT/20) for CL, both age strata (unitless)")     # Yu 2025 Table 2: theta(WT) = 0.45 (RSE 11%, bootstrap median 0.45, 95% CI 0.35-0.55)
    e_bun_cl  <- -0.0823; label("Power exponent on (BUN/3) for CL, both age strata (unitless)")     # Yu 2025 Table 2: theta(BUN) = -0.0823 (RSE 24%, bootstrap median -0.0779, 95% CI -0.1553 to -0.0093)
    e_bsa_vc  <- 1.10;    label("Power exponent on (BSA/0.77) for V1 (unitless)")                   # Yu 2025 Table 2: theta(BSA) = 1.10 (RSE 10%, bootstrap median 1.09, 95% CI 0.76-1.45)

    # ---------------------------------------------------------------------
    # Inter-individual variability. Exponential model (Methods: "Exponential
    # models were used to describe inter-individual variability").
    #
    # Table 2 reports these three rows as "Interindividual variability
    # omega CL (%)" = 21.24, "omega Q (%)" = 27.66 and "omega V2 (%)" = 41.47.
    # The row symbol is omega (an SD), not omega^2, and the values are
    # quoted as PERCENTAGES, which a variance cannot be. So these are
    # SD-scale values: var(etalcl) = 0.2124^2, and likewise for Q and V2.
    #
    # The variance reading is not merely implausible, it is falsified by the
    # paper's own Table 3. Re-simulating the paper's virtual 5-year-old
    # (19 kg, BSA 0.77) under the SD reading reproduces the published
    # incidence of delayed excretion across all six eGFR x risk-group cells
    # (28.6 / 72.5 / 85.7 and 52.4 / 88.0 / 95.4 against the published
    # 25.4 / 67.0 / 84.3 and 50.1 / 86.5 / 99.5); reading the same numbers as
    # variances gives 39.8 / 60.6 / 70.2 and 52.7 / 70.4 / 79.5, which is
    # wrong in BOTH directions and misses the 99.5% cell by 20 points. See
    # the vignette for the full check.
    #
    # There is deliberately NO etalvc: Table 2 reports IIV on CL, Q and V2
    # only, and the paper never mentions estimating or removing IIV on V1.
    # ---------------------------------------------------------------------

    etalcl ~ 0.04511376  # 0.2124^2  # Yu 2025 Table 2: omega(CL) = 21.24% (RSE 14%, bootstrap median 20.73, 95% CI 17.61-24.33)
    etalq  ~ 0.07650756  # 0.2766^2  # Yu 2025 Table 2: omega(Q)  = 27.66% (RSE 31%, bootstrap median 27.60, 95% CI 9.38-38)
    etalvp ~ 0.17197609  # 0.4147^2  # Yu 2025 Table 2: omega(V2) = 41.47% (RSE 27%, bootstrap median 40.24, 95% CI 27.39-51.86)

    # ---------------------------------------------------------------------
    # Residual error: proportional (Methods: "proportional models were
    # applied to characterize residual variability").
    # ---------------------------------------------------------------------

    propSd <- 0.4416; label("Proportional residual error (fraction)")  # Yu 2025 Table 2: residual unexplained variability = 44.16% (RSE 5%, bootstrap median 43.82, 95% CI 40.74-47.32)
  })

  model({
    # ---- 1. Age-stratum indicator (Yu 2025 p. 8480, cutoff at 1 year) ----
    is_agele1 <- (AGE <= 1.0)

    # ---- 2. Individual PK parameters ----
    # Only the typical clearance is stratum-specific; the eGFR, weight and
    # BUN power effects and the exponential IIV are shared across strata.
    cl_age <- exp(lcl_agele1) * is_agele1 + exp(lcl_agegt1) * (1.0 - is_agele1)
    cl <- cl_age * (CRCL / 160.0)^e_crcl_cl * (WT / 20.0)^e_wt_cl *
      (BUN / 3.0)^e_bun_cl * exp(etalcl)

    vc <- exp(lvc) * (BSA / 0.77)^e_bsa_vc
    q  <- exp(lq) * exp(etalq)
    vp <- exp(lvp) * exp(etalvp)

    # ---- 3. Micro-constants ----
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- 4. ODE system (two compartments, first-order elimination from the
    #         central compartment; NONMEM ADVAN3 TRANS4). Methotrexate is
    #         given as a 24-hour IV infusion, so dosing goes straight into
    #         `central` and there is no absorption compartment -- the
    #         Methods phrase "with first-order absorption and elimination"
    #         is inconsistent with the ADVAN3 TRANS4 it names in the same
    #         sentence, with the all-intravenous regimen, and with Table 2,
    #         which reports no absorption parameter. See vignette Errata. ----
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ---- 5. Observation and proportional residual error ----
    # Dose in umol, vc in L -> central/vc has units umol/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
