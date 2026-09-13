Jia_2026_rivaroxaban <- function() {
  description <- paste0(
    "One-compartment population PK model for rivaroxaban in 38 adults after ",
    "transjugular intrahepatic portosystemic shunt (TIPS) placement, with ",
    "sequential zero-order (D1 = 0.831 h) then first-order (Ka = 0.140/h) ",
    "absorption, an absorption lag time (1.23 h) and linear elimination. No ",
    "covariate reached significance, so the model is covariate-free despite a ",
    "29-variable screen. CL/F is 7.48 L/h and V/F only 4.75 L, the latter ",
    "markedly below the 21.7-101 L reported in non-TIPS populations; because ",
    "Ka (0.140/h) is far smaller than kel (CL/F divided by V/F = 1.57/h) the ",
    "disposition is flip-flop, so the terminal phase is absorption-rate-limited ",
    "and Cmax is set by Dose times Ka divided by CL/F almost independently of ",
    "V/F (Jia 2026)"
  )
  reference <- paste(
    "Jia M, Chai Y, Gao Y, Jing C, Zhu K, Zhu T, Wang L, Sun A, Yang J, Zhu Y,",
    "Feng Y, Cao Y, Li J. Population pharmacokinetics of rivaroxaban after",
    "transjugular intrahepatic portosystemic shunt. Eur J Clin Pharmacol. 2026.",
    "doi:10.1007/s00228-026-04034-6.",
    "De-identified concentration-time data deposited by the authors at",
    "doi:10.5281/zenodo.17035573.",
    sep = " "
  )
  vignette <- "Jia_2026_rivaroxaban"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Rivaroxaban was given orally, 5 or 10 mg once daily
  # (Methods "Study design, population and treatment"), and measured as total
  # drug in plasma by a validated LC-MS/MS assay with a 1-1000 ug/L calibration
  # range (Methods "Laboratory analysis").
  compartmentData <- list(
    depot   = list(analyte = "rivaroxaban", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "rivaroxaban", units = "mg", specimen = "plasma", verified = TRUE)
  )

  # The final model carries NO covariates, so there is no `covariateData`.
  #
  # This is a substantive negative result rather than an omission: the authors
  # screened 29 variables by forward inclusion (p < 0.01, dOFV > 6.63) then
  # backward elimination (p < 0.001, dOFV > 10.83) and retained none. Results
  # "Covariate analysis": "A nominal associated effect of ALT on CL/F (P < 0.01)
  # was detected in the univariate screen, and the variance component dropped to
  # 0.0091. However, further forward inclusion did not reach the predetermined
  # significance level (dOFV < 6.63; p < 0.01). No other demographic, renal, or
  # hepatic variables significantly decreased IIV on any structural parameter;
  # hence no covariate was included in the final model."
  #
  # Every entry below therefore has NO point estimate in the paper and cannot be
  # encoded. The distribution quoted for each is Table 1, as median (min-max)
  # for continuous variables and n (%) for categorical ones; those ranges are
  # what bound the screen, and the paper repeatedly attributes the null results
  # to their narrowness. Canonical register names are used where one exists;
  # entries marked "not a register canonical" carry the source dataset's own
  # column name (from the deposited Zenodo dataset) because no canonical covers
  # the concept. These are documentation only -- `checkModelConventions()` does
  # not resolve `covariatesDataExcluded` names against the register, and none is
  # referenced in `model()`.
  covariatesDataExcluded <- list(
    ALT = list(
        description = "Alanine aminotransferase",
        units       = "U/L",
        type        = "continuous",
        notes       = paste(
          "THE NEAR MISS, and the only covariate the paper singles out. Significant on",
          "CL/F in the univariate screen at p < 0.01, dropping the CL/F variance",
          "component to 0.0091, but it failed to clear the forward-inclusion threshold",
          "(dOFV < 6.63) and was not retained (Results 'Covariate analysis'). No",
          "coefficient is printed, so the effect cannot be encoded. Median 36 U/L",
          "(range 7-289), Table 1."
        )
      ),
      CRCL = list(
        description = "Cockcroft-Gault creatinine clearance",
        units       = "mL/min",
        type        = "continuous",
        notes       = paste(
          "THE NOTABLE NEGATIVE. Creatinine clearance is a retained covariate on CL/F in",
          "essentially every non-TIPS rivaroxaban popPK model, and its absence here is",
          "discussed at length: 'Whereas prior PopPK studies in non-shunted populations",
          "have consistently observed creatinine clearance (CrCl) as an important",
          "predictor of rivaroxaban clearance, no such effect was noted among our",
          "post-TIPS subjects'. The authors give three reasons: the CrCl range was narrow",
          "(no patient below 40 mL/min), creatinine-based estimates overstate renal",
          "function in cirrhosis because sarcopenia lowers creatinine production, and",
          "portosystemic shunting raises F so that a rise in intrinsic clearance is",
          "offset when clearance is expressed as CL/F. Computed with the",
          "Cockcroft-Gault equation from serum creatinine in umol/L (Methods 'Data",
          "collection'). Median 127.1 mL/min (range 50.5-341.0), Table 1 -- note this is",
          "RAW mL/min, not BSA-normalized.",
          "Compare the retained CRCL power effect in `modellib('Lai_2026_rivaroxaban')`."
        )
      ),
      ASCITES = list(
        description        = "Baseline ascites severity (none / mild / moderate-to-severe)",
        units              = "(categorical)",
        type               = "categorical",
        reference_category = "None",
        notes              = paste(
          "Explicitly explored on V/F and explicitly rejected, which matters because the",
          "paper's headline finding is a 78-95% reduction in V/F: 'Regarding ascites, we",
          "explored baseline ascites status (none/mild/moderate-severe; see Table 1) as a",
          "covariate on Vd/F. Ascites was not supported as a significant covariate",
          "(p > 0.05) and was therefore not retained in the final model, suggesting that,",
          "within the range of ascites severity represented in our cohort, ascites alone",
          "did not explain the lower Vd/F.' The authors reason that rivaroxaban is highly",
          "protein-bound and not primarily distributed into extracellular fluid, so",
          "ascites may proxy broader decompensation rather than add distribution volume,",
          "while cautioning that power was limited. None 11 (28.9%), mild 16 (42.1%),",
          "moderate-to-severe 11 (28.9%), Table 1. Not a register canonical; the",
          "deposited dataset carries no ascites column at all."
        )
      ),
      AGE = list(
        description = "Age",
        units       = "years",
        type        = "continuous",
        notes       = "Screened, not retained. Median 57 years (range 32-76), Table 1; eligibility was 18-70 years (Methods)."
      ),
      WT = list(
        description = "Actual body weight",
        units       = "kg",
        type        = "continuous",
        notes       = paste(
          "Screened, not retained; the covariate list names it 'weight (actual body",
          "weight)' (Methods 'Covariate model'). No allometric term appears in the final",
          "model. Median 62 kg (range 47.5-100), Table 1. Still needed indirectly as an",
          "input to the Cockcroft-Gault CrCl."
        )
      ),
      HT = list(
        description = "Body height",
        units       = "cm",
        type        = "continuous",
        notes       = "Screened, not retained. Median 168 cm (range 155-190), Table 1."
      ),
      SEXF = list(
        description        = "Female sex indicator",
        units              = "(binary)",
        type               = "binary",
        reference_category = "0 (male)",
        notes              = paste(
          "Screened as a categorical covariate by analysis of variance (Methods",
          "'Covariate model'), not retained. 14 of 38 female (36.8%), Table 1. The",
          "deposited dataset codes SEX with 1 = female for the two subjects shown, so a",
          "user mapping that column must confirm the direction before use."
        )
      ),
      EGFR = list(
        description = "Estimated glomerular filtration rate",
        units       = "mL/min/1.73 m^2",
        type        = "continuous",
        notes       = paste(
          "Screened alongside CrCl as the second renal descriptor, not retained. Median",
          "108.8 (range 62.9-148.5), Table 1. Table 1 mislabels the unit as 'mL/min/L';",
          "eGFR is conventionally BSA-normalized to mL/min/1.73 m^2, which is the unit",
          "recorded here. See the CRCL entry for why no renal effect was detectable."
        )
      ),
      AST = list(
        description = "Aspartate aminotransferase",
        units       = "U/L",
        type        = "continuous",
        notes       = "Screened, not retained. Median 35 U/L (range 9-230), Table 1."
      ),
      GGT = list(
        description = "Gamma-glutamyl transpeptidase",
        units       = "U/L",
        type        = "continuous",
        notes       = "Screened, not retained. Median 34 U/L (range 11-144), Table 1."
      ),
      ALP = list(
        description = "Alkaline phosphatase",
        units       = "U/L",
        type        = "continuous",
        notes       = "Screened, not retained. Median 88 U/L (range 31-291), Table 1."
      ),
      TBA = list(
        description = "Total bile acids",
        units       = "umol/L",
        type        = "continuous",
        notes       = "Screened, not retained. Median 61.2 umol/L (range 17.5-237.9), Table 1."
      ),
      TBILI = list(
        description = "Total bilirubin",
        units       = "umol/L",
        type        = "continuous",
        notes       = "Screened, not retained. Median 25.1 umol/L (range 10.6-65), Table 1."
      ),
      DBIL = list(
        description = "Direct (conjugated) bilirubin",
        units       = "umol/L",
        type        = "continuous",
        notes       = "Screened, not retained. Median 11.9 umol/L (range 5.7-45.9), Table 1."
      ),
      TPRO = list(
        description = "Total serum protein",
        units       = "g/L",
        type        = "continuous",
        notes       = "Screened, not retained; the paper's column is 'TP'. Median 56.9 g/L (range 37.6-76.8), Table 1."
      ),
      ALB = list(
        description = "Serum albumin",
        units       = "g/L",
        type        = "continuous",
        notes       = paste(
          "Screened, not retained. Median 29.7 g/L (range 25.2-38.7), Table 1 -- uniformly",
          "low, which the Discussion uses as evidence of advanced liver disease and hence",
          "of why creatinine-based renal markers mislead in this cohort."
        )
      ),
      GLB = list(
        description = "Serum globulin",
        units       = "g/L",
        type        = "continuous",
        notes       = "Screened, not retained. Median 26.1 g/L (range 10.7-45.7), Table 1. Not a register canonical; source column 'GLB'."
      ),
      PT_SEC = list(
        description = "Prothrombin time",
        units       = "s",
        type        = "continuous",
        notes       = "Screened, not retained; the paper's column is 'PT'. Median 11.2 s (range 8.9-15.4), Table 1."
      ),
      PTA = list(
        description = "Prothrombin activity",
        units       = "%",
        type        = "continuous",
        notes       = "Screened, not retained. Median 65.4% (range 41.6-89.6), Table 1. Not a register canonical; source column 'PTA'."
      ),
      PTR = list(
        description = "Prothrombin time ratio",
        units       = "(ratio)",
        type        = "continuous",
        notes       = "Screened, not retained. Median 1.23 (range 0.98-1.69), Table 1."
      ),
      INR = list(
        description = "International Normalized Ratio",
        units       = "(ratio)",
        type        = "continuous",
        notes       = "Screened, not retained. Median 1.2 (range 0.98-1.59), Table 1. Not a register canonical; source column 'INR'."
      ),
      APTT = list(
        description = "Activated partial thromboplastin time",
        units       = "s",
        type        = "continuous",
        notes       = "Screened, not retained. Median 43.1 s (range 30.7-71.1), Table 1. Not a register canonical; source column 'APTT'."
      ),
      APTTR = list(
        description = "Activated partial thromboplastin time ratio",
        units       = "(ratio)",
        type        = "continuous",
        notes       = "Screened, not retained. Median 1.35 (range 0.96-2.22), Table 1. Not a register canonical; source column 'APTTR'."
      ),
      FIB = list(
        description = "Plasma fibrinogen",
        units       = "g/L",
        type        = "continuous",
        notes       = "Screened, not retained. Median 2.1 g/L (range 1.02-3.99), Table 1."
      ),
      TT = list(
        description = "Thrombin time",
        units       = "s",
        type        = "continuous",
        notes       = "Screened, not retained. Median 20.3 s (range 17.1-23.7), Table 1. Not a register canonical; source column 'TT'."
      ),
      DDIMER = list(
        description = "Plasma D-dimer",
        units       = "ug/L",
        type        = "continuous",
        notes       = paste(
          "Screened, not retained; the paper's column is 'DD'. Median 3.1 (range 0.3-8.9),",
          "Table 1. Table 1 gives the unit as ug/L, but values of 0.3-8.9 are far below any",
          "plausible ug/L D-dimer and match the mg/L FEU scale used clinically in China, so",
          "the tabulated unit is very likely mis-stated; a user supplying this column",
          "should confirm the scale against their own assay."
        )
      ),
      FDP = list(
        description = "Fibrin degradation products",
        units       = "mg/L",
        type        = "continuous",
        notes       = "Screened, not retained. Median 13.4 mg/L (range 1.07-39.2), Table 1. Not a register canonical; source column 'FDP'."
      ),
      ATA = list(
        description = "Antithrombin activity",
        units       = "%",
        type        = "continuous",
        notes       = paste(
          "Screened, not retained. Median 52.2% (range 29.1-95.4), Table 1. Not a register",
          "canonical; source column 'ATA'. Distinct from the registered AT_BL_UDL, which is",
          "a per-subject BASELINE antithrombin activity on a U/dL scale."
        )
      ),
      CTP = list(
        description        = "Child-Turcotte-Pugh classification (A or B)",
        units              = "(categorical)",
        type               = "categorical",
        reference_category = "A",
        notes              = paste(
          "Screened as a categorical covariate, not retained. A 15 (39.5%), B 23 (60.5%),",
          "Table 1. NO patient was class C -- a stated limitation ('The cohort lacked",
          "representation of patients with advanced hepatic dysfunction (CTP class C),",
          "limiting our ability to characterize pharmacokinetics in this most severely",
          "impaired subgroup'), so the model carries no information about decompensated",
          "class-C disease. Not a register canonical; the deposited dataset codes CTP as",
          "1 = A, 2 = B."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 38,
    n_studies      = 1,
    n_observations = 131,
    age_range      = "32-76 years (Table 1 median 57); eligibility 18-70 years",
    age_median     = "57 years",
    weight_range   = "47.5-100 kg",
    weight_median  = "62 kg",
    sex_female_pct = 36.8,
    race_ethnicity = c(Asian = 100),
    disease_state  = paste(
      "Portal hypertension treated by transjugular intrahepatic portosystemic",
      "shunt (TIPS), receiving post-operative rivaroxaban for stent-thrombosis",
      "prophylaxis. Child-Turcotte-Pugh class A 39.5% / B 60.5%, no class C.",
      "Ascites: none 28.9%, mild 42.1%, moderate-to-severe 28.9%."
    ),
    dose_range     = "5 mg once daily (30 patients) or 10 mg once daily (8 patients), orally, from post-operative day 3",
    regions        = "China (single center: Beijing Youan Hospital, Capital Medical University, Beijing)",
    renal_function = "Cockcroft-Gault CrCl median 127.1 mL/min (range 50.5-341.0); eGFR median 108.8 (62.9-148.5). No patient had CrCl below 40 mL/min, which the Discussion cites as the reason no renal effect on CL/F was detectable.",
    hepatic_function = "Cirrhotic, Child-Turcotte-Pugh A/B only. Albumin median 29.7 g/L, total bilirubin median 25.1 umol/L, ALT median 36 U/L.",
    co_medication  = "None relevant: 'No co-administered medications acting on P450 enzymes or P-gp inducing/inhibiting drugs were identified among the enrolled patients' (Methods 'Covariate model').",
    notes          = paste(
      "Prospective single-center study run July 2023 to March 2025; Chinese",
      "Clinical Trial Registry ChiCTR2300073784. Table 1 baseline demographics.",
      "39 patients enrolled and 136 samples collected; 3 samples were invalid",
      "(protocol deviations or below the 1 ug/L limit of quantification) leaving",
      "133 from 38 patients, and a sensitivity analysis then excluded 2",
      "influential points, giving the 131 observations from 38 patients used to",
      "fit the final model (Results 'Patient characteristics' and 'Population PK",
      "model'). Sampling is sparse by design because post-TIPS patients are",
      "clinically fragile: 2, 4 and 24 h (+/- 0.5 h) after the first",
      "pharmacist-observed dose for all subjects, with an additional 8 h sample",
      "for the last 23 participants after a protocol revision. The authors warn",
      "that this 'inevitably reduced the precision of absorption-phase estimates',",
      "and that 'the sparse early sampling design limits precise identification",
      "of true Cmax'. Assay: validated LC-MS/MS, verapamil internal standard,",
      "1-1000 ug/L calibration range, intra- and inter-day accuracy and precision",
      "below 15%. Model built in NONMEM 7.5.0 with PsN 5.2.6 using FOCE-I.",
      "Bootstrap converged in only 897 of 1000 replicates (89.7%), and the",
      "authors caution that V/F is 'the least reliably estimated parameter'.",
      "No pharmacodynamic endpoint (bleeding, thrombosis) and no anti-Xa",
      "calibration were collected, and ABCB1 / ABCG2 / CYP3A4 genotypes were not",
      "analysed."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural parameters, all from Table 2 column "Final Estimate".
    # Time in h, dose in mg, CL/F in L/h, V/F in L. The Results text repeats
    # every one of these five values with its RSE, and the Abstract repeats
    # CL/F and V/F ("The apparent clearance and volume of distribution were
    # 7.48 L/h and 4.75 L, respectively").
    #
    # ABSORPTION IS SEQUENTIAL ZERO-ORDER THEN FIRST-ORDER. The dose is
    # delivered into `depot` at a constant rate over D1 hours starting at
    # `tlag`, then transfers to `central` with first-order rate ka. The
    # authors' deposited dataset (doi:10.5281/zenodo.17035573) confirms this
    # encoding independently of the prose: dose records carry CMT = 1 and
    # RATE = -2 (NONMEM's flag for a duration modelled by D1) while
    # observations carry CMT = 2.
    #
    # NOTE ON FLIP-FLOP KINETICS. kel = CL/F / (V/F) = 7.48 / 4.75 = 1.57/h
    # (t1/2 = 0.44 h) is an order of magnitude FASTER than ka = 0.140/h
    # (t1/2 = 4.95 h), so absorption is rate-limiting and the observed
    # terminal slope reflects ka, not elimination. A practical consequence is
    # that Cmax is approximately Dose * ka / (CL/F) and therefore nearly
    # independent of V/F, which is why the strikingly small V/F does not
    # translate into implausible peak concentrations.
    lka   <- log(0.140); label("Absorption rate constant (ka, 1/h)")                            # Table 2 row 'Ka (1/h)' 0.140 (RSE 5.97%)
    lcl   <- log(7.48);  label("Apparent clearance (CL/F, L/h)")                                # Table 2 row 'CL/F (L/h)' 7.48 (RSE 9.52%)
    lvc   <- log(4.75);  label("Apparent central volume of distribution (Vd/F, L)")             # Table 2 row 'Vd/F (L)' 4.75 (RSE 37.4%)
    ld1   <- log(0.831); label("Duration of the zero-order absorption input (D1, h)")           # Table 2 row 'D1 (h)' 0.831 (RSE 29.3%)
    ltlag <- log(1.23);  label("Absorption lag time (ALAG1, h)")                                # Table 2 row 'ALAG1 (h)' 1.23 (RSE 20.0%)

    # ---------------------------------------------------------------------
    # IIV. Table 2 prints four IIV rows as percentages and the Results text
    # calls them coefficients of variation ("Inter-individual variability
    # (IIV), expressed as the coefficient of variation, was 61.5% (RSE 24.6%)
    # for CL/F and 83.4% (RSE 51.4%) for Vd/F"). Two conventions could produce
    # such a percentage, and they differ materially here:
    #
    #   (a) omega^2 = (CV/100)^2                  -> etalcl = 0.378
    #   (b) omega^2 = log((CV/100)^2 + 1)         -> etalcl = 0.321
    #
    # Table 2's own algebra cannot settle it (the bootstrap-median column is
    # wildly different from the Final Estimate column, and several CIs have
    # negative lower bounds). It was settled instead by RE-FITTING this exact
    # structural model to the authors' own deposited dataset
    # (doi:10.5281/zenodo.17035573) with nlmixr2 FOCEI, which reproduced the
    # published fit closely -- ka 0.1398 vs 0.140, CL/F 7.44 vs 7.48, Vd/F 4.65
    # vs 4.75, D1 0.844 vs 0.831, ALAG1 1.2300 vs 1.23 -- and returned omega^2
    # of 0.3785 / 0.7118 / 0.4501 / 0.0409. Those are the SQUARES of the
    # printed percentages (0.615^2 = 0.378, 0.834^2 = 0.696, 0.658^2 = 0.433,
    # 0.204^2 = 0.0416), not the log-normal form (b), which would have given
    # 0.321 / 0.471 / 0.353 / 0.0411. Convention (a) is therefore used, and the
    # values below are the exact squares of Table 2.
    #
    # There is NO IIV on ka: Table 2 has no IIV_Ka row. No correlations between
    # etas are reported, so the matrix is diagonal.
    etalcl   ~ 0.378225  # Table 2 row 'IIV_CL/F (%)' 61.5 (RSE 24.6%); 0.615^2
    etalvc   ~ 0.695556  # Table 2 row 'IIV_Vd/F (%)' 83.4 (RSE 51.4%); 0.834^2
    etald1   ~ 0.432964  # Table 2 row 'IIV_D1 (%)' 65.8 (RSE 69.4%); 0.658^2
    etaltlag ~ 0.041616  # Table 2 row 'IIV_ALAG1 (%)' 20.4 (RSE 96.2%); 0.204^2

    # ---------------------------------------------------------------------
    # Residual error: combined proportional plus additive. Table 2 prints
    # 0.0825 (unitless) and 4.08 (labelled ug/L) without saying whether either
    # is a variance or an SD. The same re-fit settles it: nlmixr2 returned
    # propSd 0.2837 and addSd 2.010, i.e. the SQUARE ROOTS of the printed
    # numbers (sqrt(0.0825) = 0.2872, sqrt(4.08) = 2.0199). Both Table 2 rows
    # are therefore NONMEM $SIGMA VARIANCES, and Table 2's "ug/L" unit label on
    # the additive row is loose -- as a variance that quantity is in ug^2/L^2.
    # Taking 0.0825 at face value as an SD would understate proportional
    # residual error more than threefold (8.25% instead of 28.7% CV).
    propSd <- 0.2872; label("Proportional residual error (fraction)")  # sqrt of Table 2 row 'Proportional residual error (-)' 0.0825 (RSE 24.2%)
    addSd  <- 2.02;   label("Additive residual error (ug/L)")          # sqrt of Table 2 row 'Additional residual error (ug/L)' 4.08 (RSE 52.0%)
  })

  model({
    # Individual parameters. Log-normal IIV on CL/F, V/F, D1 and the lag time;
    # ka has no IIV (see ini()).
    ka   <- exp(lka)
    cl   <- exp(lcl + etalcl)
    vc   <- exp(lvc + etalvc)
    d1   <- exp(ld1 + etald1)
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Sequential zero-order then first-order absorption. `dur(depot) <- d1`
    # spreads each dose uniformly over d1 hours, and `alag(depot) <- tlag`
    # delays its start. IMPORTANT: dose records must carry `rate = -2` so
    # rxode2 uses the modelled duration; without that flag the dose enters
    # `depot` as an instantaneous bolus and the zero-order phase is skipped.
    dur(depot)  <- d1
    alag(depot) <- tlag

    # Dose in mg and vc in L give mg/L; x1000 converts to the ug/L used for
    # concentrations throughout the paper (the assay range is 1-1000 ug/L and
    # the safety threshold is Cmax,ss <= 140 ug/L). Note the paper reports
    # AUCss,24 in mg*h/L, so an AUC computed from this ug/L output must be
    # divided by 1000 to compare against the 1.77 mg*h/L threshold.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd) + add(addSd)
  })
}
