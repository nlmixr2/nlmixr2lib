Liu_2025_tranexamicAcid <- function() {
  description <- "Two-compartment population PK model for intravenous tranexamic acid (TXA) with first-order elimination, in Chinese adults undergoing cardiac surgery with cardiopulmonary bypass; allometric body weight on all four disposition parameters with exponents fixed at 0.75 (clearances) and 1 (volumes) (Liu 2025)."
  reference   <- "Liu Y, Zhou C, Lv H, Tian L, Jiang J, Shi J. Population Pharmacokinetics of Tranexamic Acid in Chinese Population Undergoing Cardiac Surgery with Cardiopulmonary Bypass. Drug Des Devel Ther. 2025;19:4343-4353. doi:10.2147/DDDT.S493485"
  vignette    <- "Liu_2025_tranexamicAcid"
  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Actual body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The sole covariate retained in the final model. Enters as an allometric power term",
        "centred on 70 kg -- (WT/70)^0.75 on both clearances and (WT/70)^1 on both volumes",
        "(Liu 2025 Table 3, final row). The display equation typeset between the Methods",
        "'Covariate Models' paragraphs gives the general form as V1 = theta1 x (WT/70)^theta2.",
        "Liu 2025 Methods states 70 kg was chosen because the cohort mean body weight was close",
        "to 70 kg (74.3 +/- 19.9 kg high-dose, 68.2 +/- 12.2 kg low-dose; Table 1) and because",
        "the same reference had been used by Grassin-Delyle et al. The exponents were NOT",
        "estimated: Methods states 'exponent theta2 was set as a fixed value, while the exponent",
        "for the clearance parameter was typically 0.75 and the exponent for the volume of",
        "distribution parameter was typically 1', and Table 4 reports no %RSE for either.",
        "Time-fixed at the pre-operative value.",
        sep = " "
      ),
      source_name        = "BW"
    )
  )

  # Screened during covariate model building (Liu 2025 Methods, "Covariate
  # Models"; 73 stepwise covariate models were constructed) but NOT retained in
  # the final model, so none of these is referenced in model(). The best
  # stepwise model (the 23rd) kept sex on CL2 and age on CL1, but Liu 2025
  # Results record that it "was more complex but did not explain additional
  # variability" and it scored worse than the body-weight model on every
  # criterion (-2LL 1437 vs 1391, AIC 1455 vs 1409, BIC 1485 vs 1439; Table 3),
  # so the authors discarded it. Its coefficients are recorded in the notes
  # below for provenance only -- they belong to a rejected candidate model and
  # must not be combined with the final-model estimates in ini().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened but not retained. The rejected 23rd stepwise model carried age on CL1 as",
        "CL1 = CL1_typical x (AGE/55)^0.75 (Liu 2025 Table 3, middle row). Cohort 51.4 +/- 11.3",
        "years (high-dose) and 59.3 +/- 9.9 years (low-dose); inclusion criteria required 18-70",
        "years (Methods).",
        sep = " "
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "categorical",
      notes       = paste(
        "Screened but not retained. The rejected 23rd stepwise model carried sex on CL2 as a",
        "proportional effect CL2 = CL2_typical x (1 + 0.75 x Sex) (Liu 2025 Table 3, middle row);",
        "the paper does not state which sex the indicator codes as 1, which is a further reason",
        "the term is not reproduced here. Cohort 4/7 female (high-dose) and 1/9 female",
        "(low-dose), i.e. 5/16 = 31.3% overall (Table 1).",
        sep = " "
      )
    ),
    CPB_ON = list(
      description = "Cardiopulmonary bypass phase indicator",
      units       = "(binary)",
      type        = "categorical",
      notes       = paste(
        "Screened but not retained. Liu 2025 followed the method of Dowd et al in testing the",
        "state of CPB as a dichotomous covariate, and the Discussion reports that no significant",
        "effect of CPB on any PK parameter was detected -- consistent with Grassin-Delyle et al.",
        "The authors attribute the null result partly to the pump prime dose, which offsets the",
        "haemodilution that would otherwise mark the onset of bypass, and partly to sparse",
        "sampling during the bypass window itself.",
        sep = " "
      )
    ),
    T_CPB = list(
      description = "Total cardiopulmonary bypass duration",
      units       = "min",
      type        = "continuous",
      notes       = paste(
        "Screened as a continuous covariate but not retained (Liu 2025 Discussion). Cohort",
        "129.71 +/- 36.91 min (high-dose) and 116.44 +/- 50.80 min (low-dose); Table 2.",
        sep = " "
      )
    ),
    BODYTEMP = list(
      description = "Minimum rectal temperature during cardiopulmonary bypass",
      units       = "degC",
      type        = "continuous",
      notes       = paste(
        "Screened as a continuous covariate but not retained (Liu 2025 Methods 'Covariate Models'",
        "and Discussion). This is the intra-operative nadir rectal temperature during the bypass",
        "run, not an admission body temperature; Liu 2025 does not tabulate its distribution.",
        sep = " "
      )
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened but not retained. Cohort 9.32 +/- 3.88 (high-dose) and 13.73 +/- 7.92 umol/L (low-dose); Liu 2025 Table 2."
    ),
    DBIL = list(
      description = "Direct (conjugated) bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened but not retained. Cohort 2.45 +/- 0.81 (high-dose) and 5.28 +/- 3.42 umol/L (low-dose); Liu 2025 Table 2."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained. Cohort 29.71 +/- 22.47 (high-dose) and 23.67 +/- 9.90 IU/L (low-dose); Liu 2025 Table 2."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained. Cohort 29.29 +/- 11.54 (high-dose) and 25.44 +/- 4.95 IU/L (low-dose); Liu 2025 Table 2."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = paste(
        "Screened but not retained. Cohort 78.15 +/- 13.89 (high-dose) and 91.53 +/- 30.01 umol/L",
        "(low-dose); Liu 2025 Table 2. Note that the sibling model Nakai_2025_tranexamicAcid.R,",
        "fitted to a Japanese cohort in the same clinical setting, DID retain renal function",
        "(Cockcroft-Gault creatinine clearance) on CL; Liu 2025 tested serum creatinine directly",
        "rather than a derived clearance, in a cohort with a narrower creatinine range.",
        sep = " "
      )
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened but not retained. Cohort 5.55 +/- 2.05 (high-dose) and 5.82 +/- 1.50 mmol/L (low-dose); Liu 2025 Table 2."
    )
  )

  compartmentData <- list(
    central     = list(analyte = "tranexamic acid", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tranexamic acid", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 16L,
    n_studies      = 1L,
    n_observations = 224L,
    age_range      = "18-70 years by protocol; 51.4 +/- 11.3 years (high-dose) and 59.3 +/- 9.9 years (low-dose)",
    weight_range   = "74.3 +/- 19.9 kg (high-dose) and 68.2 +/- 12.2 kg (low-dose)",
    bmi_range      = "25.7 +/- 5.4 (high-dose) and 23.9 +/- 3.3 kg/m^2 (low-dose)",
    sex_female_pct = 31.3,
    race_ethnicity = "Chinese",
    disease_state  = "Adults undergoing cardiac surgery with cardiopulmonary bypass; NYHA class I-III, no class IV. Hypertension in 7/16, diabetes 2/16, hyperlipidaemia 5/16 (Liu 2025 Table 1).",
    renal_function = "Serum creatinine 78.15 +/- 13.89 umol/L (high-dose) and 91.53 +/- 30.01 umol/L (low-dose); blood urea nitrogen 5.55 +/- 2.05 and 5.82 +/- 1.50 mmol/L. End-stage disease with expected survival under 3 months was an exclusion criterion; no dialysis patients are described.",
    dose_range     = "Randomised to a high-dose arm (n = 7; 30 mg/kg loading dose infused over 20 min after induction of anaesthesia, 16 mg/kg/h maintenance infusion until the end of the operation, and a 2 mg/kg pump prime dose added to the CPB priming solution) or a low-dose arm (n = 9; 10 mg/kg loading, 2 mg/kg/h maintenance, 1 mg/kg pump prime).",
    regions        = "Single centre, Fuwai Hospital, Chinese Academy of Medical Sciences, Beijing, China (Ethics No. 2022-1866)",
    notes          = paste(
      "Prospective randomised study, 16 participants, all with the full 14-timepoint sampling",
      "schedule and no observations below the 1 ug/mL lower limit of quantification and no",
      "missing values or outliers (Liu 2025 Results). Sampling: pre-dose; 10 min after the start",
      "of the loading dose; end of the loading dose; 30 min, 1 h and 2 h after the start of the",
      "maintenance dose; end of the maintenance dose; and 15 min, 30 min, 1 h, 2 h, 3 h, 4 h and",
      "6 h after the end of the maintenance dose. Plasma TXA by LC-MS/MS. Surgery duration",
      "285.29 +/- 66.05 min (high-dose) and 272.44 +/- 61.28 min (low-dose); CPB duration",
      "129.71 +/- 36.91 and 116.44 +/- 50.80 min (Table 2). Estimation in Phoenix NLME 8.3 with",
      "FOCE-ELS. The final model was confirmed by a 1000-replicate bootstrap (95.8% convergence)",
      "and a 1000-replicate visual predictive check. This is the first published population PK",
      "model of tranexamic acid in a Chinese population.",
      sep = " "
    )
  )

  ini({
    # Structural parameters. Liu 2025 Table 4 "Estimate (% RSE)" column, which
    # reports the typical value for a 70 kg reference subject; the same four
    # numbers are restated in the Abstract Results and in the Discussion
    # ("CL1 = 4.7 L/h, V1 = 4.9 L, CL2 = 17.0 L/h and V2 = 11.1 L"), and all
    # four fall inside the bootstrap 95% CI of Table 5.
    lcl <- log(4.7);  label("Clearance CL1 at WT 70 kg (L/h)")                          # Liu 2025 Table 4 (%RSE 6.89; bootstrap median 4.8, 95% CI 4.2-5.4)
    lvc <- log(4.9);  label("Central volume of distribution V1 at WT 70 kg (L)")        # Liu 2025 Table 4 (%RSE 9.86; bootstrap median 4.8, 95% CI 0.1-5.8)
    lq  <- log(17.0); label("Inter-compartmental clearance CL2 at WT 70 kg (L/h)")      # Liu 2025 Table 4 (%RSE 21.36; bootstrap median 17.0, 95% CI 12.18-48.8)
    lvp <- log(11.1); label("Peripheral volume of distribution V2 at WT 70 kg (L)")     # Liu 2025 Table 4 (%RSE 6.83; bootstrap median 11.1, 95% CI 9.8-15.7)

    # Allometric exponents. Liu 2025 Methods, "Covariate Models": "exponent
    # theta2 was set as a fixed value, while the exponent for the clearance
    # parameter was typically 0.75 and the exponent for the volume of
    # distribution parameter was typically 1". Table 3's final row and Table 4's
    # "Covariate effect" column print them as literal exponents with no %RSE, so
    # both are encoded as fixed().
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent on (WT/70) for CL1 and CL2 (unitless)")  # Liu 2025 Table 4 "Covariate effect" (BW/70)^0.75
    e_wt_vc_vp <- fixed(1);    label("Allometric exponent on (WT/70) for V1 and V2 (unitless)")    # Liu 2025 Table 4 "Covariate effect" (BW/70)^1

    # Between-subject variability. Liu 2025 Methods, "Random Effect Model":
    # "The between-subject variability (BSV) [...] in random effects [was]
    # represented by [an] exponential [model]", i.e. P_i = TV(P) x exp(eta_i).
    #
    # SCALE OF THE TABLE 4 BSV COLUMN. Table 4 heads the column "BSV (% RSE)
    # (shrinkage)" without stating whether the number is the variance omega^2 or
    # the standard deviation omega, and Phoenix NLME can print either. The
    # values are read here as STANDARD DEVIATIONS on the log scale (so
    # CV ~= omega), on the strength of the paper's own Figure 1: at 140 min the
    # seven high-dose profiles read 210.7, 199.9, 191.0, 182.0, 164.3, 117.6 and
    # 108.6 mg/L, a between-subject SD of log concentration of 0.263, and the
    # nine low-dose profiles span only 33-51 mg/L. That 0.263 is the TOTAL
    # observed dispersion, so it also contains the residual error; the model's
    # BSV alone therefore cannot exceed it. Simulating the same read-off with the
    # residual at zero gives a BSV-only log-SD of 0.191 under the SD reading
    # (compatible, leaving an 18% residual budget) against 0.335 under the
    # variance reading (omega = 0.529, 0.663, 0.548, 0.469), which would require a
    # negative residual variance and is therefore impossible. The vignette runs
    # this comparison as an executable gate. Weighed against it, the reported
    # %RSEs of 37-48% on 16 subjects sit
    # closer to the asymptotic RSE of a variance (sqrt(2/n) = 35%) than of an SD
    # (1/sqrt(2n) = 18%), so the reading is an inference from the paper's data
    # rather than a printed statement -- see the vignette Errata.
    etalcl ~ 0.28^2  # Liu 2025 Table 4 BSV for CL1: 0.28 (%RSE 37, shrinkage 0.03), read as a log-scale SD
    etalq  ~ 0.44^2  # Liu 2025 Table 4 BSV for CL2: 0.44 (%RSE 44, shrinkage 0.21), read as a log-scale SD
    etalvc ~ 0.30^2  # Liu 2025 Table 4 BSV for V1:  0.30 (%RSE 48, shrinkage 0.16), read as a log-scale SD
    etalvp ~ 0.22^2  # Liu 2025 Table 4 BSV for V2:  0.22 (%RSE 45, shrinkage 0.12), read as a log-scale SD

    # Residual error. Liu 2025 Methods, "Random Effect Model", declares that the
    # within-subject variability (WSV) "was represented by [a] proportional
    # [model]", and the Results confirm "the WSV was suitable for the
    # proportional model". No magnitude is reported anywhere in the paper:
    # Table 4 lists four rows only (CL1, CL2, V1, V2) with no residual row,
    # Table 5 reports bootstrap results for the same four structural parameters,
    # and the two supplementary items the paper cites (Figures S1 and S2) are
    # concentration-time profiles rather than a parameter table. The declared
    # structure is therefore recorded with a fixed(0) magnitude rather than an
    # invented one; set propSd to a plausible value before using this model for
    # any simulation that needs residual scatter.
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - declared by the source but not reported)")  # Liu 2025 Methods, Random Effect Model (WSV proportional; magnitude never reported)
  })

  model({
    # Allometric body-size scaling, centred on 70 kg. Liu 2025 Table 3 final row
    # applies it to all four disposition parameters:
    #   CL1 = CL1_typical x (BW/70)^0.75    CL2 = CL2_typical x (BW/70)^0.75
    #   V1  = V1_typical  x (BW/70)^1       V2  = V2_typical  x (BW/70)^1
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q
    q  <- exp(lq  + etalq)  * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp

    # Micro-constants for the explicit two-compartment system.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Intravenous administration: the loading infusion, the maintenance infusion
    # and the pump prime dose all enter the central compartment directly
    # (Liu 2025 Methods, "Study Subjects and Dosing Regimen"). No depot and no
    # bioavailability term.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Dose in mg and vc in L give mg/L, the unit used on the Liu 2025 Figure 1
    # and Figure 3 concentration axes.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
