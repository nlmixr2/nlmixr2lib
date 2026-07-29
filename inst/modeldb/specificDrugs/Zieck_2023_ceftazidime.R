Zieck_2023_ceftazidime <- function() {
  description <- "One-compartment IV population PK model for ceftazidime in adult general-ward patients spanning adequate to severely impaired renal function (Zieck 2023); clearance is a power function of CKD-EPI estimated glomerular filtration rate and is 1.56-fold higher with concomitant use of another antibiotic."
  reference <- paste(
    "Zieck SE, de Vroom SL, Mulder FP, van Twillert G, Mathot RAA,",
    "Geerlings SE, van Hest RM. Pharmacokinetic/pharmacodynamic target",
    "attainment of ceftazidime in adult patients on general wards with",
    "different degrees of renal function: a prospective observational",
    "bicenter cohort study.",
    "Antibiotics (Basel). 2023;12(3):469.",
    "doi:10.3390/antibiotics12030469.",
    sep = " "
  )
  vignette <- "Zieck_2023_ceftazidime"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "CKD-EPI estimated glomerular filtration rate, BSA-normalized",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column CKDEPI (Supplementary File S1 $INPUT). Renal function",
        "was estimated with the CKD-EPI (Chronic Kidney Disease Epidemiology",
        "Collaboration) creatinine equation on the day of inclusion",
        "(Zieck 2023 Table 2 footnote a). Power-form effect on CL centred on",
        "76.85 mL/min/1.73 m^2: (CRCL / 76.85)^0.75, per the Supplementary",
        "File S1 $PK line 'TVCL=THETA(3)*(CKDEPI/76.85)**THETA(4)*",
        "THETA(6)**FLAG1'. The Table 3 footnote prints the same normaliser as",
        "76.86; the control-stream value 76.85 is used here because it is the",
        "value the estimated model actually executed (the 0.013% difference is",
        "numerically irrelevant). The exponent 0.75 was ESTIMATED, not fixed",
        "(RSE 13.9%, bootstrap 95% CI 0.56-0.93; the $THETA comment reads",
        "'est exponent TVCL effect CKDEPI'), so it is not the canonical",
        "allometric 3/4 power. Cohort median eGFR 73.5 (IQR 34.3-111.4)",
        "mL/min/1.73 m^2, spanning 10.6 to 124.8 across the three renal",
        "function strata (Table 2). Also tested but not retained: MDRD and",
        "Cockcroft-Gault eGFR, and raw serum creatinine (Methods 4.6)."
      ),
      source_name        = "CKDEPI"
    ),
    CONMED_ABX = list(
      description        = "Concomitant use of another (non-ceftazidime) antibiotic",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant other antibiotic)",
      notes              = paste(
        "Source column COMED (Supplementary File S1 $INPUT), entering as",
        "FLAG1 via 'IF(COMED.EQ.1)FLAG1=1'. Multiplicative power-form effect",
        "on CL: 1.56^CONMED_ABX (Zieck 2023 Table 3 final model, RSE 12.4%,",
        "bootstrap 1.57, 95% CI 1.20-1.94), i.e. the multiplier is 1 when",
        "CONMED_ABX = 0 and 1.56 when CONMED_ABX = 1. 27 of 40 patients (68%)",
        "received a concomitant other antibiotic (Table 2). The source paper",
        "never names the co-administered agents, so this is an umbrella",
        "any-other-antibiotic indicator rather than a per-agent covariate.",
        "Inclusion dropped IIV on CL from 37.6% to 31.3% CV (16.8% of the IIV",
        "in CL; Results 2.2). The authors flag the association as",
        "physiologically unexplained and possibly coincidental (Discussion),",
        "so it should be read as an empirical marker of a sicker / more",
        "heavily co-treated subset rather than a mechanistic drug-drug",
        "interaction."
      ),
      source_name        = "COMED"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 40L,
    n_studies      = 1L,
    n_observations = 119L,
    age_range      = "Adults (>=18 years by inclusion criteria); oldest subject 86 years (Discussion)",
    age_median     = "62.0 years (IQR 47.0-72.0)",
    weight_range   = "Median 79.6 kg (IQR 69.7-92.3); stratum IQRs span 57.1-140.6 kg",
    weight_median  = "79.6 kg",
    height_median  = "175.5 cm (IQR 167.0-185.0)",
    bmi_median     = "25.0 kg/m^2 (IQR 22.0-29.0)",
    sex_female_pct = 42.5,
    race_ethnicity = c(Caucasian = 80, `African American` = 10, Asian = 7.5, Hispanic = 2.5),
    disease_state  = paste(
      "Adults admitted to a general (non-ICU) hospital ward receiving",
      "therapeutic ceftazidime as standard care for a suspected or proven",
      "Gram-negative infection. Patients in the ICU, on renal replacement",
      "therapy, with cystic fibrosis, or with severe burns were excluded",
      "because of known altered pharmacokinetics (Methods 4.2)."
    ),
    renal_function = paste(
      "CKD-EPI eGFR median 73.5 (IQR 34.3-111.4) mL/min/1.73 m^2. Three",
      "prespecified strata: adequate (eGFR >= 50, n = 25, median 102.8),",
      "moderate impairment (eGFR 30-50, n = 10, median 34.3), and severe",
      "impairment (eGFR < 30, n = 5, median 18.6). Serum creatinine median",
      "100.0 (IQR 67.3-135.3) umol/L (Table 2)."
    ),
    co_medication  = "27/40 (68%) on a concomitant other antibiotic (Table 2)",
    dose_range     = paste(
      "Guideline-recommended, renal-function-stratified IV regimens infused",
      "over 0.5 h (Table 1): 2000 mg q8h for adequate renal function,",
      "1000 mg q12h for moderate impairment, and 1000 mg q24h for severe",
      "impairment. One severely impaired patient received 2000 mg q24h",
      "instead of the recommended 1000 mg q24h and was retained for the",
      "primary outcome only (Results 2.1)."
    ),
    regions        = "The Netherlands (bicenter: Amsterdam UMC location AMC, and Noordwest Ziekenhuisgroep location Alkmaar)",
    notes          = paste(
      "Prospective observational cohort study run October 2019 to December",
      "2021. Baseline demographics per Zieck 2023 Table 2 (medians with",
      "interquartile ranges). Three samples per patient (one trough plus two",
      "random) were drawn within 72 h of treatment start; 121 samples were",
      "collected, 2 excluded, leaving 119 for NONMEM 7.5 estimation, of which",
      "1 was below the 0.1 mg/L LLOQ and none above the 40 mg/L ULOQ. Five",
      "patients contributed only two samples. The model was internally",
      "validated with a 1000-run bootstrap and a prediction-corrected VPC",
      "(Figure 1)."
    )
  )

  ini({
    # ===== Structural PK (Zieck 2023 Table 3 "Final Model" column; ==========
    # ===== values confirmed against Supplementary File S1 $THETA) ==========
    # One-compartment IV disposition, CL/V parameterization
    # ($SUBROUTINES ADVAN1 TRANS2, S1 = V).
    lcl <- log(3.74); label("Typical CL at CKD-EPI eGFR 76.85 mL/min/1.73 m^2 with no concomitant antibiotic (L/h)")  # Table 3 Final Model CL = 3.74 L/h (RSE 9.80%, bootstrap 3.74, 95% CI 3.03-4.41); File S1 $THETA(3)
    lvc <- log(21.8); label("Volume of distribution (L)")                                                             # Table 3 Final Model V  = 21.8 L   (RSE 7.90%, bootstrap 22.1, 95% CI 19.1-25.1); File S1 $THETA(5)

    # ===== Covariate effects on CL (Zieck 2023 Table 3 footnote) ===========
    # CL (L/h) = 3.74 * (CKDEPI / 76.85)^0.75 * 1.56^flag
    # Both coefficients were estimated (neither is fixed): the exponent
    # carries RSE 13.9% and bootstrap CI 0.56-0.93, so it is NOT the
    # canonical allometric 0.75 held constant -- it coincidentally
    # rounds to the same value.
    # File S1 $THETA(4) is commented as an estimated (not fixed) exponent;
    # $THETA(6) is commented as the concomitant-antibiotic (AB) effect.
    e_crcl_cl       <- 0.75; label("Power exponent on (CRCL / 76.85 mL/min/1.73 m^2) for CL (unitless)")                          # Table 3 Final Model eGFR (CKD-EPI) on CL = 0.75 (RSE 13.9%, bootstrap 0.74, 95% CI 0.56-0.93); File S1 $THETA(4)
    e_conmed_abx_cl <- 1.56; label("Multiplicative factor on CL for concomitant other-antibiotic use (CONMED_ABX = 1, unitless)")  # Table 3 Final Model concomitant antibiotic use on CL = 1.56 (RSE 12.4%, bootstrap 1.57, 95% CI 1.20-1.94); File S1 $THETA(6)

    # ===== IIV (exponential / log-normal on CL and V) =====================
    # File S1: CL = TVCL * EXP(ETA(1)); V = THETA(5) * EXP(ETA(2)).
    # Table 3 reports IIV as %CV; omega^2 = log(1 + CV^2).
    #
    # etalvc uses the Table 3 final-model 40.2% CV rather than the File S1
    # $OMEGA(2) value of 0.157 (which implies 41.2% CV and matches none of
    # the three V figures the paper reports: structural 40.5%, final 40.2%,
    # bootstrap 40.9%). Control-stream $OMEGA entries are commonly stale
    # initial estimates, and this block also carries a leftover editing
    # comment. By contrast $OMEGA(1) = 0.0936 for CL reproduces Table 3
    # exactly, so it is used verbatim (one more significant digit than the
    # rounded table value). Operator-confirmed: sidecar oare_PMC10044023
    # request-001 q2. See the vignette Errata.
    etalcl ~ 0.0936  # File S1 $OMEGA(1) = 0.0936 -> CV = sqrt(exp(0.0936) - 1) = 31.3%, matching Table 3 Final Model IIV CL 31.3% CV (RSE 29.6%, shrinkage 7%)
    etalvc ~ 0.1498  # Table 3 Final Model IIV V = 40.2% CV (RSE 61.8%, shrinkage 18%) -> omega^2 = log(1 + 0.402^2) = 0.1498

    # ===== Residual error (proportional only) =============================
    # File S1 $ERROR: W = IPRED*THETA(1) + THETA(2); Y = IPRED + W*ERR(1),
    # with $THETA(2) = (0 fix) (additive term switched off) and $SIGMA 1 FIX.
    # The proportional SD is therefore THETA(1) directly.
    propSd <- 0.186; label("Proportional residual error (fraction)")  # Table 3 Final Model proportional error = 18.6% (RSE 15.9%, bootstrap 18.7, 95% CI 13.9-23.6); File S1 $THETA(1) = 0.186
  })

  model({
    # ----- Covariate factors on CL -------------------------------------
    # Zieck 2023 Table 3 footnote / File S1 $PK:
    #   FLAG1 = 0; IF(COMED.EQ.1) FLAG1 = 1
    #   TVCL  = THETA(3) * (CKDEPI / 76.85)**THETA(4) * THETA(6)**FLAG1
    renal_factor <- (CRCL / 76.85)^e_crcl_cl
    abx_factor   <- e_conmed_abx_cl^CONMED_ABX

    # ----- Individual PK parameters ------------------------------------
    cl <- exp(lcl + etalcl) * renal_factor * abx_factor
    vc <- exp(lvc + etalvc)

    # ----- Micro-constants ---------------------------------------------
    kel <- cl / vc

    # ----- ODE system ---------------------------------------------------
    # Ceftazidime is given as an intermittent IV infusion (0.5 h) directly
    # into the central compartment; infusion duration/rate comes from the
    # event table (File S1 carries RATE in $INPUT).
    d/dt(central) <- -kel * central

    # ----- Output --------------------------------------------------------
    # Total ceftazidime plasma concentration: dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
