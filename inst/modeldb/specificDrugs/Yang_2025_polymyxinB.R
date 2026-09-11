Yang_2025_polymyxinB <- function() {
  description <- "One-compartment intravenous population PK model for polymyxin B in critically ill adults with carbapenem-resistant organism infections, built from paired steady-state trough and peak plasma concentrations (Yang 2025). Creatinine clearance and platelet count are power covariates on clearance, both normalized to the modeling-set median (CrCL 75.99 mL/min, PLT 163.50 x 10^9/L) with exponents 0.26 and -0.14. Volume of distribution carries no covariate. Combined proportional plus additive residual error."
  reference <- paste(
    "Yang J, Yu M, Gan Y, Cheng L, Yang G, Xiong L, Liu F, Chen Y.",
    "Population pharmacokinetics of polymyxin B in critically ill patients",
    "with carbapenem-resistant organisms infections: insights from",
    "steady-state trough and peak plasma concentration.",
    "Front Pharmacol. 2025;16:1511088.",
    "doi:10.3389/fphar.2025.1511088. PMCID PMC11936910.",
    sep = " "
  )
  vignette <- "Yang_2025_polymyxinB"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Yang 2025: the assay measures total
  # polymyxin B (PMB1 + PMB2 summed on a molar basis) in plasma, and the
  # single disposition compartment is the plasma/central compartment of the
  # one-compartment model of Equations 6-7.
  compartmentData <- list(
    central = list(analyte = "polymyxinB", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, NOT BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Yang 2025 Methods 'Data collection' item (4): 'liver and renal function indices, with CrCL calculated using the Cockcroft-Gault equation'. Raw mL/min, not normalized to 1.73 m^2 BSA; stored under the canonical CRCL column per inst/references/covariate-columns.md, which accepts raw mL/min provided the assay form is documented per model (precedent: Delattre_2010_amikacin.R, Chen_2023_nemonoxacin.R, Valade_2015_emtricitabine.R). Applied as a power covariate on CL, CL = 2.03 * (CRCL / 75.99)^0.26 * ... (Equation 6). The normalizing constant is printed only as the symbol 'CrCLmedian' in Equation 6; 75.99 mL/min is the modeling-set median from Table 1 (75.99, IQR 38.46-130.58). See the vignette source-trace section for the Table 3 arithmetic that confirms this reading. Time-fixed per subject in this analysis (a single baseline CrCL per patient).",
      source_name        = "CrCL"
    ),
    PLT = list(
      description        = "Platelet (thrombocyte) count from the routine complete blood count",
      units              = "10^9/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Yang 2025 Table 1 'PLT (10^9/L)'. Applied as a power covariate on CL, CL = ... * (PLT / 163.50)^(-0.14) (Equation 6); a HIGHER platelet count gives a LOWER clearance. The normalizing constant is printed only as the symbol 'PLTmedian' in Equation 6; 163.50 x 10^9/L is the modeling-set median from Table 1 (163.50, IQR 84.5-266.25). Yang 2025 Discussion states this is the first report of platelet count as a covariate on polymyxin B PK, and the Limitations paragraph notes that patients with elevated platelet counts were excluded from modeling, so the relationship is supported only over roughly 85-266 x 10^9/L (the modeling-set interquartile range, which is also the range used for the Monte Carlo simulations of Table 3). Time-fixed per subject in this analysis.",
      source_name        = "PLT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened per Yang 2025 Methods ('Covariates included age, weight, APACHE II score, the presence of sepsis, and all laboratory parameters listed in Table 1') using Spearman correlation against empirical Bayes estimates followed by stepwise forward addition / backward elimination. Not retained: only CrCL and PLT met the inclusion criteria, and no covariate was retained on V."
    ),
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened (Yang 2025 Methods) and not retained. Notable because Table 4 shows that two of the five comparator polymyxin B popPK studies (Manchandani et al., Crass et al.) did retain total body weight on CL."
    ),
    APACHEII = list(
      description = "Acute Physiology and Chronic Health Evaluation II severity-of-illness score",
      units       = "(points)",
      type        = "continuous",
      notes       = "Screened (Yang 2025 Methods) and not retained. Modeling-set value 28 +/- 9.25 (Table 1)."
    ),
    SEPSIS = list(
      description = "Presence of sepsis",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened (Yang 2025 Methods, 'the presence of sepsis') and not retained."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened as one of the Table 1 laboratory parameters and not retained. Yang 2025 Limitations notes that only total (not unbound) polymyxin B was measured, so an albumin-driven protein-binding effect could not be resolved."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened as one of the Table 1 laboratory parameters and not retained; the derived Cockcroft-Gault CRCL was retained instead."
    ),
    OTHER_TABLE1_LABS = list(
      description = "Remaining Table 1 laboratory parameters screened as candidate covariates",
      units       = "(various)",
      type        = "continuous",
      notes       = "BUN, eGFR, ALT, AST, total protein, total bilirubin, hemoglobin, white blood cell count, IL-6, INR, APTT and fibrinogen were all screened per Yang 2025 Methods ('all laboratory parameters listed in Table 1') and none was retained in the final model. Grouped into a single entry because the paper reports the screen collectively rather than per analyte."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 80L,
    n_studies      = 1L,
    age_median     = "60 years (IQR 47-74)",
    weight_median  = "63 kg (IQR 55-74)",
    sex_female_pct = 25,
    race_ethnicity = "Not reported (single-center Chinese ICU cohort)",
    disease_state  = "Critically ill adult ICU patients with microbiologically confirmed carbapenem-resistant organism (CRO) infection: carbapenem-resistant Acinetobacter baumannii (70 cases), carbapenem-resistant Enterobacterales (56 cases) and carbapenem-resistant Pseudomonas aeruginosa (15 cases); lung was the dominant infection site (81.25%). APACHE II 28 +/- 9.25. Patients receiving any form of renal replacement therapy during polymyxin B treatment were excluded.",
    dose_range     = "Intravenous polymyxin B, loading dose in 97.5% of the modeling set followed by a maintenance dose of 1.31 +/- 0.25 mg/kg per administration; Monte Carlo simulations covered 50, 75, 100 and 125 mg q12h given as 1-hour infusions. Median treatment duration 13 days (IQR 9-16).",
    regions        = "China (intensive care unit, First Affiliated Hospital of Army Medical University, Chong Qing)",
    renal_function = "Cockcroft-Gault CrCL median 75.99 mL/min (IQR 38.46-130.58); serum creatinine median 79.70 umol/L (IQR 51.60-158.90). Renal replacement therapy was an exclusion criterion.",
    n_observations = "184 polymyxin B plasma concentrations from 80 patients in the modeling set (11 patients had repeated sampling). Two samples per occasion: one immediately before an infusion (steady-state trough) and one immediately after (steady-state peak), collected after at least 48 h of therapy. An additional 15 patients / 30 samples formed a chronologically separated external validation set.",
    notes          = "Single-center prospective study, August 2021 to July 2024; ethics approval No. (A) KY2021064. Baseline demographics and laboratory parameters in Yang 2025 Table 1 (modeling set column). Concentrations measured by validated UPLC-MS/MS with polymyxin E2 as internal standard; total polymyxin B was computed by summing the molar contributions of PMB1 and PMB2. External validation gave MPE% 2.69, MAPE% 28.45, F20 36.67% and F30 73.33%."
  )

  ini({
    # Structural parameters -- Yang 2025 Table 2 (final model), typical values
    # at the modeling-set median CrCL (75.99 mL/min) and PLT (163.50 x 10^9/L).
    lcl <- log(2.03); label("Clearance CL (L/h) at CRCL = 75.99 mL/min and PLT = 163.5 x 10^9/L") # Yang 2025 Table 2: CL = 2.03 L/h (RSE 5.30%; bootstrap median 2.02, 95% CI 1.81-2.24)
    lvc <- log(18);   label("Volume of distribution V (L)")                                       # Yang 2025 Table 2: V = 18 L (RSE 5.30%; bootstrap median 17.98, 95% CI 16.30-19.80)

    # Covariate effects on CL -- Yang 2025 Equation 6:
    #   CL = 2.03 * (CrCL/CrCLmedian)^0.26 * (PLT/PLTmedian)^(-0.14) * exp(eta_CL)
    # Table 2 names the two exponents dCLdCrCL and dCLdPLT. Both were
    # ESTIMATED (each carries an RSE and a bootstrap 95% CI), so neither is
    # wrapped in fixed(); the Table 2 footnote wording "fixed parameter
    # coefficient" is the NONMEM sense of a fixed (population) effect as
    # opposed to a random effect, not a held-constant value.
    e_crcl_cl <- 0.26;  label("Power exponent on (CRCL / 75.99 mL/min) for CL (unitless)")        # Yang 2025 Table 2: dCLdCrCL = 0.26 (RSE 20.40%; bootstrap 95% CI 0.15-0.36)
    e_plt_cl  <- -0.14; label("Power exponent on (PLT / 163.5 x 10^9/L) for CL (unitless)")       # Yang 2025 Table 2: dCLdPLT = -0.14 (RSE 24.20%; bootstrap 95% CI -0.22 to -0.066)

    # Inter-individual variability. Yang 2025 Table 2 reports a single IIV
    # term, "etaCL (%) 38.50" (RSE 10.70%; bootstrap median 38.03, 95% CI
    # 29.70-46.60), under an exponential random-effects model (Results,
    # "PopPK model analysis and validation"). The row is read as an
    # SD-scale/CV percentage, omega = 0.385, giving variance 0.385^2 =
    # 0.148225 -- NOT as a variance omega^2 = 0.385 (which would imply
    # omega = 0.62). The Table 2 abbreviation footnote glosses eta as
    # "variance of inter-individual variability", but three independent
    # checks refute that gloss; see the vignette source-trace section:
    #   (1) the "(%)" in the row label is meaningless for a variance;
    #   (2) the model-based AUCss,24h IQR of Table 1 (55.81-94.07, i.e. a
    #       1.686-fold spread) implies a total SD of log(CL) of about 0.387,
    #       which already CAPS omega below 0.62 before any covariate
    #       contribution is removed;
    #   (3) reconstructing the Table 3 AUCss,24h band percentages from
    #       Equation 6 requires omega ~ 0.385 (e.g. 75 mg q12h, CrCL 90,
    #       PLT 85 gives P(AUC < 50) = 25.4% against the published 25.8%,
    #       and P(AUC > 100) = 12.8% against the published 14.2%).
    etalcl ~ 0.148225 # Yang 2025 Table 2: etaCL = 38.50% -> omega = 0.385 -> variance 0.148225

    # Equation 7 is written V = 18 L * exp(eta_V), but Table 2 reports NO
    # variance for eta_V (the "Inter-individual variability" block contains
    # only the etaCL row) and no supplement on disk supplies one. Per the
    # standing unreported-IIV policy the term is encoded as an explicit
    # zero-variance random effect rather than invented; see the vignette
    # Assumptions and deviations section.
    etalvc ~ fixed(0) # Yang 2025 Equation 7 has exp(eta_V) but Table 2 reports no eta_V variance

    # Residual variability -- Yang 2025 Table 2, combined proportional plus
    # additive error (Results: "residual variability was described using both
    # proportional and additive error models"). Both rows are read on the SD
    # scale, not as NONMEM $SIGMA variances: at variance scale the residual CV
    # would be about 55% and the external-validation metrics of the Results
    # (MAPE% 28.45, F30 73.33%) could not be attained. See the vignette
    # source-trace section for the arithmetic.
    propSd <- 0.30; label("Proportional residual error (fraction)") # Yang 2025 Table 2: proportional error = 0.30 (RSE 8.20%; bootstrap 95% CI 0.25-0.35)
    addSd  <- 0.21; label("Additive residual error (mg/L)")         # Yang 2025 Table 2: additive error = 0.21 (RSE 28.40%; bootstrap 95% CI 0.070-0.35)
  })
  model({
    # Individual PK parameters -- Yang 2025 Equations 6 and 7. The two
    # normalizing constants are the modeling-set medians of Table 1
    # (CrCL 75.99 mL/min, PLT 163.50 x 10^9/L); Equation 6 prints them only
    # as the symbols CrCLmedian and PLTmedian.
    cl <- exp(lcl + etalcl) * (CRCL / 75.99)^e_crcl_cl * (PLT / 163.5)^e_plt_cl
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    # One-compartment disposition with first-order elimination; polymyxin B is
    # given as an intravenous infusion, so drug enters the central compartment
    # directly (Yang 2025 Results, "A one-compartment model with first-order
    # elimination best fit the population data").
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
