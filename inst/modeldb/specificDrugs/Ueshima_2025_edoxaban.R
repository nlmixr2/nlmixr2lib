Ueshima_2025_edoxaban <- function() {
  description <- "One-compartment population PK model for oral edoxaban in Japanese adults with non-valvular atrial fibrillation, fitted to sparse real-world therapeutic-drug-monitoring data (one steady-state sample per patient, 9-47 h after the last dose). Absorption is not modelled: with no data on the absorption phase the authors used NONMEM ADVAN1 TRANS2, so a dose enters the central compartment directly. Apparent oral clearance carries a power effect of raw Cockcroft-Gault creatinine clearance; the apparent volume of distribution is fixed to the ENGAGE AF-TIMI 48 Asian value with a body-weight exponent fixed at 1. CYP3A5*3 and the ABCB1 1236C>T / 2677G>T,A / 3435C>T polymorphisms (and the 2677/3435 haplotype) were screened and none affected edoxaban pharmacokinetics (Ueshima 2025)."
  reference   <- "Ueshima S, Hira D, Matsuda S, Michihata R, Tabuchi Y, Ozawa T, Itoh H, Iguchi M, Akao M, Aizawa T, Kashiwa A, Shizuta S, Makiyama T, Nakagawa Y, Horie M, Terada T, Katsura T. Population pharmacokinetics and pharmacogenomics of edoxaban in Japanese adults with atrial fibrillation. J Pharm Health Care Sci. 2025;11:46. doi:10.1186/s40780-025-00453-2"
  vignette    <- "Ueshima_2025_edoxaban"
  units       <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance estimated by the Cockcroft-Gault equation, raw mL/min, NOT body-surface-area normalized",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power scaling (CRCL/61.8)^0.692 on CL/F; 61.8 mL/min is the cohort median (Ueshima 2025 Table 1),",
        "and the same constant is printed inside the final-model equation (Eq. 10) and inside the Table 2",
        "parameter header 'CL/F (L/h) = theta1 x (CLcr/61.8)^theta3'. Computed with the Cockcroft-Gault",
        "equation per Ueshima 2025 Methods (reference [23], Cockcroft & Gault 1976); serum creatinine was",
        "read from the routine clinical record, so CRCL is time-fixed at the value contemporaneous with the",
        "single PK sample. Observed range 20.0-135.5 mL/min. The paper's own Monte Carlo simulations",
        "(Fig. 3) extrapolate the term to 30, 60, 90 and 120 mL/min, which stays inside the fitted range.",
        "Two structural notes carried from the source. First, the authors decomposed CL/F into apparent",
        "renal (CLR/F) and apparent non-renal (CLNR/F) arms (Eq. 4) and fitted the CLcr effect on the renal",
        "arm; the population mean of CLNR/F converged to 0.01 L/h, which they judged negligible, so the",
        "final model (Eq. 10) is the single power term encoded here and carries no separate non-renal",
        "intercept. Second, the source itself flags that this is likely an underestimate of the true",
        "non-renal arm -- an intravenous study puts renal clearance at 49.1% of total -- because no urine or",
        "metabolite samples were collected and each patient contributed only one concentration.",
        "Raw (non-BSA-normalized) Cockcroft-Gault CRCL follows the precedent of Nakai_2025_tranexamicAcid.R,",
        "Wada_2023_sparsentan.R and Chen_2023_nemonoxacin.R.",
        sep = " "
      ),
      source_name        = "CLcr"
    ),
    WT = list(
      description        = "Actual body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Linear scaling (WT/70)^1 on Vd/F per Ueshima 2025 Eq. 3, Vd/F = 336.4 x (BW/70). BOTH parts of",
        "this term are fixed, not estimated: the exponent was fixed at 1 'due to a lack of individual data'",
        "(Methods, PPK modeling) and the 336.4 L intercept was taken from the ENGAGE AF-TIMI 48 analysis",
        "(reference [22], Krekels 2016) as 193 L central volume x 1.743 Asian fold-change. Consequently",
        "body weight does NOT act on clearance in this model. Cohort weight median 61.4 kg, range",
        "37.3-97.3 kg (Table 1); the 70 kg reference is therefore above the cohort median and is inherited",
        "from the upstream ENGAGE AF-TIMI 48 parameterisation rather than being a cohort-centred value.",
        sep = " "
      ),
      source_name        = "BW"
    )
  )

  # Covariates the source screened but did not retain in the final model. None
  # of these is referenced in model(); they are recorded so the provenance of
  # the paper's covariate screen -- which is the whole point of the paper, since
  # its headline finding is a set of NEGATIVE pharmacogenomic results -- is not
  # lost. Ueshima 2025 Results, "PPK modeling": stepwise forward inclusion
  # (dOBJ > 3.84) and backward elimination (dOBJ > 6.63), plus a 95%-confidence-
  # interval rule on the effect size, retained CLcr alone.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a power (Eq. 7) or linear (Eq. 8) effect on CLR/F and CLNR/F; not retained (Ueshima 2025 Results, PPK modeling). Cohort median 72.2 years, range 35.2-92.6."
    ),
    SEXF = list(
      description = "Sex, 1 = female",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a dichotomous effect (Eq. 9); not retained (Ueshima 2025 Results, PPK modeling). Cohort 83 male / 48 female."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened as a power (Eq. 7) or linear (Eq. 8) effect on CLR/F and CLNR/F; not retained. Cohort median 19 IU/L, range 5-77. The source reads this as consistent with the clinical finding that mild-to-moderate hepatic impairment does not change edoxaban pharmacokinetics (Discussion, reference [29])."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened as a power (Eq. 7) or linear (Eq. 8) effect on CLR/F and CLNR/F; not retained. Cohort median 25 IU/L, range 11-92."
    ),
    CONMED_CYP3A4_PGP_INH = list(
      description = "Concomitant CYP3A4 and/or P-glycoprotein inhibitor, 1 = on an inhibitor",
      units       = "(binary)",
      type        = "binary",
      notes       = "Pooled indicator over amiodarone (n = 5), diltiazem (n = 6) and verapamil (n = 5); Ueshima 2025 Table 1. Screened as a dichotomous effect (Eq. 9); not retained. The source flags the small exposed n as a limitation and notes the result runs contrary to dedicated DDI studies and to the ENGAGE AF-TIMI 48 popPK model, which found P-gp-inhibitor effects on both CL/F and bioavailability -- effects this study could not detect because it has no absorption-phase data."
    ),
    SNP_CYP3A5_RS776746 = list(
      description = "CYP3A5 rs776746 (CYP3A5*3) genotype",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened as a dichotomous effect (Eq. 9); not retained (Ueshima 2025 Results, PPK modeling). Genotype counts *1/*1, *1/*3, *3/*3 = 2, 38, 91 (Table 1); consistent with Hardy-Weinberg equilibrium. Note the register also carries the derived-phenotype canonical CYP3A5_EXPR; the raw-genotype SNP name is used here because the source reports and stratifies on genotype, matching Gu_2025_rivaroxaban.R."
    ),
    SNP_ABCB1_RS1128503 = list(
      description = "ABCB1 rs1128503 (1236C>T) genotype",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened as a dichotomous effect (Eq. 9); not retained. Genotype counts C/C, C/T, T/T = 56, 56, 19 (Ueshima 2025 Table 1)."
    ),
    SNP_ABCB1_RS2032582 = list(
      description = "ABCB1 rs2032582 (2677G>T/A) genotype",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened as a dichotomous effect (Eq. 9); not retained. Genotype counts G/G, G/T, G/A, A/A, T/A, T/T = 18, 20, 48, 23, 19, 3 (Ueshima 2025 Table 1)."
    ),
    SNP_ABCB1_RS1045642 = list(
      description = "ABCB1 rs1045642 (3435C>T) genotype",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened as a dichotomous effect (Eq. 9); not retained. Genotype counts C/C, C/T, T/T = 43, 65, 23 (Ueshima 2025 Table 1). The 2677G>T,A / 3435C>T HAPLOTYPE was screened separately -- the two loci are in strong linkage disequilibrium -- and was likewise not retained."
    )
  )

  compartmentData <- list(
    central = list(analyte = "edoxaban", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 131L,
    n_studies      = 1L,
    n_observations = 131L,
    age_range      = "35.2-92.6 years",
    age_median     = "72.2 years",
    weight_range   = "37.3-97.3 kg",
    weight_median  = "61.4 kg",
    sex_female_pct = 36.6,
    race_ethnicity = "Japanese",
    disease_state  = "Adults with atrial fibrillation on chronic once-daily oral edoxaban (Lixiana) for prevention of cardioembolic stroke",
    dose_range     = "15-60 mg once daily (15 mg n = 3, 30 mg n = 70, 60 mg n = 58)",
    renal_function = "Cockcroft-Gault creatinine clearance median 61.8 mL/min, range 20.0-135.5; serum creatinine median 0.89 mg/dL, range 0.56-2.00",
    hepatic_function = "AST median 25 IU/L (range 11-92); ALT median 19 IU/L (range 5-77)",
    co_medication  = "CYP3A4 and/or P-glycoprotein inhibitors in 16 of 131 patients: amiodarone 5, diltiazem 6, verapamil 5",
    regions        = "Three centres in Japan: Shiga University of Medical Science Hospital, National Hospital Organization Kyoto Medical Center, and Kyoto University Hospital (January 2017 to August 2019)",
    notes          = paste(
      "Retrospective real-world therapeutic-drug-monitoring cohort. ONE blood sample per patient, drawn at",
      "steady state 9-47 h after the last dose, so the dataset has 131 observations from 131 patients and",
      "no absorption-phase information whatsoever -- this is what forces the no-absorption structure and the",
      "fixed Vd/F. Patients with poor drug compliance or with no compliance record were excluded. Edoxaban",
      "measured by LC-ESI-MS/MS, lower limit of quantification 1 ng/mL, inter-day variation < 4.4%; observed",
      "plasma concentrations median 19.5 ng/mL, range 1.7-152.0. Genotyping by real-time PCR. Estimation in",
      "NONMEM 7.3.0 with FOCE-I; evaluated by goodness-of-fit plots, a 1000-replicate prediction-corrected",
      "VPC, and a 1000-replicate non-parametric bootstrap (PsN 4.7.0).",
      sep = " "
    )
  )

  ini({
    # Structural parameters. Ueshima 2025 Table 2, "Original data / Mean"
    # column, reproduced verbatim as the final model in Results Eq. 10:
    #   CL/F (L/h) = 28.2 x (CLcr/61.8)^0.692
    #   Vd/F (L)   = 336.4 x (BW/70)
    lcl <- log(28.2); label("Apparent oral clearance CL/F at CRCL 61.8 mL/min (L/h)") # Ueshima 2025 Table 2, theta1 (95% CI 26.9-29.5; bootstrap median 28.1)

    # Vd/F was NOT estimated. Ueshima 2025 Methods, PPK modeling: with a single
    # sample per patient on the elimination phase only, "the Vd/F was fixed at
    # the reported value for Asian in the clinical trial [22]" -- 193 L central
    # volume x 1.743 Asian fold-change = 336.4 L. Table 2 prints it as
    # "336.4 fixed" with no confidence interval. Taken from Krekels 2016
    # (ENGAGE AF-TIMI 48), not fitted here.
    lvc <- fixed(log(336.4)); label("Apparent volume of distribution Vd/F at WT 70 kg, taken from Krekels 2016 (L)") # Ueshima 2025 Table 2, theta2 / Eq. 3

    # Covariate effects.
    e_crcl_cl <- 0.692; label("Power exponent on (CRCL/61.8) for CL/F (unitless)") # Ueshima 2025 Table 2, theta3 (95% CI 0.582-0.801; bootstrap median 0.693). Estimated: the CI excludes 0, which is the source's own retention rule.
    # The body-weight exponent on Vd/F was fixed at 1 (Ueshima 2025 Methods,
    # PPK modeling: "with the exponent fixed to 1 due to a lack of individual
    # data"), so Eq. 3 is linear in weight rather than allometric. It is kept
    # as an explicit parameter rather than folded into model() so the
    # fixed-vs-estimated status stays machine-readable.
    e_wt_vc <- fixed(1); label("Power exponent on (WT/70) for Vd/F (unitless)") # Ueshima 2025 Methods, PPK modeling / Eq. 3

    # Inter-individual variability. Ueshima 2025 Eq. 2 puts a single
    # exponential eta on CL/F: (CL/F)i = theta1 x exp(eta1i). Table 2 reports
    # "omega1 (CV%) 26.4 (27.4)", where footnote b states the value is the
    # coefficient of variation of the inter-individual variability for CL/F and
    # footnote c states the parenthesised 27.4 is the eta SHRINKAGE in percent,
    # not an RSE. The reported quantity is omega1 itself expressed as a
    # percentage -- the standard NONMEM %CV = 100 x sqrt(omega^2) convention --
    # so omega1 = 0.264 and the eta variance is 0.264^2 = 0.0697.
    etalcl ~ 0.0697 # Ueshima 2025 Table 2, omega1 = 26.4 CV% (95% CI 18.8-34.0; bootstrap median 25.9); variance = 0.264^2

    # Residual error. Ueshima 2025 Eq. 1 is a pure exponential (log-normal)
    # error model, Cobs = Cpred x exp(eps), which maps to nlmixr2's lnorm()
    # with the SD on the log scale. Table 2 reports "sigma (CV%) 58.7 fixed";
    # Methods explains the value was NOT estimated here but fixed at the
    # ENGAGE AF-TIMI 48 value (reference [22]) built as the healthy-subject
    # residual 14.6% plus the atrial-fibrillation-patient increment 44.1%.
    # Same %CV convention as omega1 above, so sigma = 0.587 on the log scale.
    expSd <- fixed(0.587); label("Exponential (log-normal) residual SD, taken from Krekels 2016 (log-scale)") # Ueshima 2025 Table 2, sigma = 58.7 CV% / Methods Eq. 1
  })

  model({
    # Individual parameters. Creatinine clearance acts on CL/F only and body
    # weight on Vd/F only (Ueshima 2025 Eqs. 3 and 10); there is no allometric
    # weight term on clearance in this model.
    cl <- exp(lcl + etalcl) * (CRCL / 61.8)^e_crcl_cl
    vc <- exp(lvc) * (WT / 70)^e_wt_vc

    # One compartment, no absorption. Ueshima 2025 Methods, PPK modeling: "As
    # there were no data on the absorption phase, a 1-compartment model without
    # first-order absorption was employed ... (ADVAN1 TRANS2)". NONMEM ADVAN1 is
    # the one-compartment linear model whose dosing compartment IS the central
    # compartment, so an oral dose is placed directly into `central` with no
    # depot, no ka and no separately identified bioavailability term (F is
    # folded into the apparent CL/F and Vd/F). Event tables must therefore dose
    # into `central`.
    d/dt(central) <- -(cl / vc) * central

    # central is in mg and vc in L, so central/vc is mg/L = ug/mL; multiply by
    # 1000 to express Cc in the ng/mL used throughout Ueshima 2025.
    Cc <- central / vc * 1000
    Cc ~ lnorm(expSd)
  })
}
