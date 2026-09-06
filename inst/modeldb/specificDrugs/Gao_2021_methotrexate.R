Gao_2021_methotrexate <- function() {
  description <- paste(
    "Three-compartment population PK model for high-dose intravenous",
    "methotrexate (1-5 g/m^2) in Chinese children with acute lymphoblastic",
    "leukaemia (Gao 2021; n = 311 patients, 4,517 concentrations). Linear",
    "first-order elimination from the central compartment. Body weight enters",
    "every clearance and every volume as a FIXED allometric function",
    "(exponents 0.75 and 1.0) normalised to a 19 kg typical child, and serum",
    "creatinine lowers clearance through a linear term (-0.97% per umol/L",
    "above 26 umol/L). Exponential between-subject variability on CL and on",
    "the first intercompartmental clearance, and an additive-on-log-scale",
    "residual error equivalent to a 35.4% proportional error. Identified for",
    "external evaluation by Wang 2023, which found this model the least",
    "biased of the six it tested (MPE 3.22% population / 1.33% individual)",
    "but still imprecise (RMSE 72.39% / 63.96%).",
    sep = " "
  )
  reference <- paste(
    "Gao X, Qian XW, Zhu XH, Yu Y, Miao H, Meng JH, Jiang JY, Wang HS,",
    "Zhai XW (2021). Population Pharmacokinetics of High-Dose Methotrexate in",
    "Chinese Pediatric Patients With Acute Lymphoblastic Leukemia.",
    "Front Pharmacol 12:701452. doi:10.3389/fphar.2021.701452.",
    "Extracted as part of the Wang 2023 external-evaluation set",
    "(Wang S, Yin Q, Yang M, Cheng Z, Xie F (2023). Pharmaceutics 15(2):569.",
    "doi:10.3390/pharmaceutics15020569); see",
    "vignette('Wang_2023_methotrexate'). Parameter values here come from the",
    "Gao 2021 primary, NOT from Wang 2023 Table 2, which mis-states two of",
    "them -- see the ini() comments on etalcl and etalq.",
    sep = " "
  )
  vignette <- "Wang_2023_methotrexate"
  units    <- list(time = "h", dosing = "umol", concentration = "umol/L")

  # Issue #482. Verified against Gao 2021 Results ("the current study based on
  # 311 pediatric patients enabled to fit the MTX concentration-time data by
  # using a three-compartment disposition model") and its Table 2, which names
  # the states V_C, V_P1 and V_P2. Amount units are umol because every
  # methotrexate concentration in both Gao 2021 and Wang 2023 is reported in
  # umol/L (Wang's therapeutic targets are quoted as 65 / 33 / 16 umol/L).
  compartmentData <- list(
    central     = list(analyte = "methotrexate", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "methotrexate", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "methotrexate", units = "umol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column Weight. Allometric power terms on ALL six disposition",
        "parameters, normalised to a 19 kg typical child. Gao 2021 Table 2",
        "footnote: 'The parameters in Table 2 are given for a \"typical\" child",
        "with body weight of 19 kg. Body weight was implemented as a fixed",
        "allometric function on all clearance and volume of distribution",
        "parameters using exponent of 0.75 and 1.0, respectively.' Both",
        "exponents are FIXED, not estimated, and are wrapped in fixed() in",
        "ini() accordingly. The 19 kg reference is the cohort median (Wang",
        "2023 Table 1 records median 19.0 kg, range 4.5-113.0 for this study)."
      ),
      source_name        = "Weight"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column Scr, in SI units (umol/L), centred at 26 umol/L. Enters",
        "CL as a LINEAR (not power) term. Gao 2021 Abstract: 'The serum",
        "creatinine significantly affected the MTX clearance, with a 0.97%",
        "decrease in clearance per 1 umol/L of serum creatinine', and Gao",
        "Table 2 prints the coefficient as 'SCr on CL (%) -0.97'. Gao's own",
        "Table 2 footnote writes the equation as",
        "'[CL = CL_typical x ((SCr-26) x 0.0097)]', which is garbled -- it",
        "drops both the '1 +' and the minus sign, and would make clearance",
        "zero at the reference. The coherent form, used here and printed by",
        "Wang 2023 Table 2, is CL = CL_typical x (1 + (Scr - 26) x (-0.0097)).",
        "The 26 umol/L centre is consistent with the cohort's median serum",
        "creatinine of 0.3 mg/dL (= 26.5 umol/L; Wang 2023 Table 1).",
        "CAUTION -- this linear term is NOT bounded below: it reaches zero at",
        "SCr = 129.1 umol/L and turns NEGATIVE above it, which lies inside",
        "Gao's own reported range (0.1-1.5 mg/dL = 8.8-132.6 umol/L). The",
        "published equation is reproduced as printed rather than clamped;",
        "simulations must keep CREAT below ~129 umol/L. See the vignette",
        "Errata."
      ),
      source_name        = "Scr"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Tested and not retained. Gao 2021 Results: 'Inclusion of age-related",
        "maturation effect on CL did not show a significant improvement in",
        "model fit further.' Wang 2023 Table 1 records median 5.0 years",
        "(range 0.75-15.2) for this cohort."
      )
    ),
    TBILI = list(
      description = "Total serum bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = paste(
        "Significant but deliberately dropped. Gao 2021 Results: 'Adding TBIL",
        "on clearance (delta OFV = -30.975) ... improved model fit",
        "significantly, with slope estimates of -0.0044 ... However, these two",
        "covariates did not substantially reduce either the inter-individual",
        "or residual variability further (i.e., <5%).'"
      )
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = paste(
        "Significant but deliberately dropped alongside TBIL. Gao 2021",
        "Results: 'albumin in the central volume of distribution",
        "(delta OFV = -36.722) using linear function improved model fit",
        "significantly, with slope estimates of ... -0.070', then dropped for",
        "the same <5% variability-reduction rule."
      )
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Tested and not retained. Gao 2021 Results: 'Other covariates (e.g.,",
        "sex, AST, and ALT) did not significantly affect MTX PK properties.'",
        "Wang 2023 Table 1 records 197 male / 114 female for this cohort."
      )
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Tested and not retained (same Gao 2021 Results sentence as SEXF)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Tested and not retained (same Gao 2021 Results sentence as SEXF)."
    ),
    CONMED_OMEPRAZOLE = list(
      description = "Concomitant omeprazole indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Statistically significant but not retained. Gao 2021 Results: 'Their",
        "inclusion on clearance improved model fit significantly",
        "(delta OFV = -64.331 and -42.874, respectively); however, the",
        "reduction in either inter-individual or residual variability was",
        "minimal (<1.1%).' Co-medicated in 7.7% of patients."
      )
    ),
    CONMED_NSAID = list(
      description = "Concomitant non-steroidal anti-inflammatory drug indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Significant but not retained, alongside omeprazole (same Gao 2021 Results sentence)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 311L,
    n_studies      = 1L,
    age_range      = "0.75 to 15.2 years",
    age_median     = "5.0 years",
    weight_range   = "4.5 to 113.0 kg",
    weight_median  = "19.0 kg",
    height_range   = "67 to 175 cm (median 112)",
    sex_female_pct = 100 * 114 / 311,
    disease_state  = "Childhood acute lymphoblastic leukaemia (ALL) receiving high-dose methotrexate consolidation.",
    renal_function = paste(
      "Serum creatinine median 0.3 mg/dL (range 0.1-1.5), i.e. about",
      "26 umol/L (range 8.8-132.6). Note the model's linear SCr term on",
      "clearance goes negative above 129 umol/L, inside this range."
    ),
    hepatic_function = "ALT median 16.0 U/L (range 2.0-390.0); AST median 26.0 U/L (range 8.0-135.0).",
    dose_range     = "1 to 5 g/m^2 intravenous high-dose methotrexate.",
    regions        = "China (Children's Hospital of Fudan University, Shanghai).",
    notes          = paste(
      "Demographics from Wang 2023 Table 1 (the external-evaluation paper that",
      "tabulates all six evaluated cohorts side by side); model structure and",
      "parameter values from the Gao 2021 primary, Table 2 and Results.",
      "4,517 concentration measurements. Wang 2023 notes this is the largest",
      "of the six evaluated datasets and that, despite that, its central",
      "volume (20.7 L) is higher than the other models', producing 'an obvious",
      "underprediction for the peak concentrations of MTX'. Uncertainty in",
      "Gao 2021 was derived by sampling importance resampling (2,000 samples,",
      "1,000 resamples)."
    )
  )

  ini({
    # Structural parameters -- Gao 2021 Table 2 ('NONMEM estimates' column,
    # %RSE in parentheses), for a typical 19 kg child at SCr = 26 umol/L:
    #
    #   CL (L/h)   6.9   (2.5)      V_C (L)   20.7  (4.9)
    #   Q1 (L/h)   0.255 (7.4)      V_P1 (L)  41.0  (11.4)
    #   Q2 (L/h)   0.217 (8.7)      V_P2 (L)   3.17 (9.8)
    #
    # Wang 2023 Table 2 reproduces all six of these values correctly.
    lcl  <- log(6.9);   label("Clearance for a 19 kg child at SCr 26 umol/L (L/h)")                       # Gao 2021 Table 2 CL = 6.9 (RSE 2.5; SIR median 6.9, 95% CI 6.62-7.19)
    lvc  <- log(20.7);  label("Central volume of distribution for a 19 kg child (L)")                     # Gao 2021 Table 2 V_C = 20.7 (RSE 4.9; SIR median 20.5, 95% CI 18.5-22.4)
    lq   <- log(0.255); label("Intercompartmental clearance to peripheral1 for a 19 kg child (L/h)")      # Gao 2021 Table 2 Q1 = 0.255 (RSE 7.4; SIR median 0.258, 95% CI 0.232-0.285)
    lvp  <- log(41.0);  label("First peripheral volume of distribution for a 19 kg child (L)")            # Gao 2021 Table 2 V_P1 = 41.0 (RSE 11.4; SIR median 42.1, 95% CI 34.4-51.3)
    lq2  <- log(0.217); label("Intercompartmental clearance to peripheral2 for a 19 kg child (L/h)")      # Gao 2021 Table 2 Q2 = 0.217 (RSE 8.7; SIR median 0.224, 95% CI 0.193-0.260)
    lvp2 <- log(3.17);  label("Second peripheral volume of distribution for a 19 kg child (L)")           # Gao 2021 Table 2 V_P2 = 3.17 (RSE 9.8; SIR median 3.28, 95% CI 2.75-3.88)

    # Allometric exponents. Gao 2021 Table 2 footnote: 'Body weight was
    # implemented as a FIXED allometric function on all clearance and volume of
    # distribution parameters using exponent of 0.75 and 1.0, respectively.'
    # Fixed, not estimated -- Table 2 gives them no row, no RSE and no CI.
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on all clearances (unitless)")  # Gao 2021 Table 2 footnote
    e_wt_vc <- fixed(1);    label("Allometric exponent of body weight on all volumes (unitless)")     # Gao 2021 Table 2 footnote

    # Serum creatinine on clearance, linear.
    e_creat_cl <- -0.0097; label("Linear coefficient of serum creatinine on CL (per umol/L)")  # Gao 2021 Table 2 'SCr on CL (%) -0.97' (RSE 4.7; SIR median -0.96, 95% CI -0.91 to -1.00); Abstract '0.97% decrease in clearance per 1 umol/L'

    # Between-subject variability, exponential (Gao 2021 Eq. 1:
    # theta_i = theta . exp(eta_i,theta)). Gao's Table 2 reports a 'CV for IIV'
    # column that is populated on EXACTLY TWO rows -- CL and Q1 -- and is a long
    # dash on all four others.
    #
    # ---------------------------------------------------------------------
    # DEVIATION FROM Wang 2023 Table 2 -- Wang mis-transcribes both of these.
    #
    #   * Wang prints the CL IIV as 17.9%. Gao's Table 2 prints 17.5%
    #     (RSE 5.8; SIR median 17.6, 95% CI 16.0-19.7).
    #
    #   * Wang puts the second eta on V2 in its equation column and calls it
    #     V1 in its IIV column -- disagreeing with itself. Gao puts it on Q1,
    #     and says so twice: Table 2 populates 'CV for IIV' on the Q1 row
    #     (26.2, RSE 17.2, SIR 95% CI 18.4-33.6) and leaves V_P1/V_C blank, and
    #     Gao's Results states 'Inclusion of SCr in a linear function on
    #     clearance significantly improved model fit ... and reduced the
    #     variability of Q1 by 29.2% (CV% of IIV from 37.7 to 26.7%)'.
    #
    # The Gao primary is followed here. See the vignette Errata.
    # ---------------------------------------------------------------------
    #
    # Converted from CV% with the log-normal identity omega^2 = log(1 + CV^2).
    # (Gao's Table 2 footnote states the reverse conversion as
    # '100 x (e^variance)^1/2', which cannot be right -- it is < 1 for every
    # positive variance -- so the standard identity is used.)
    etalcl ~ 0.030165  # Gao 2021 Table 2 'CV for IIV' on the CL row = 17.5% -> log(1 + 0.175^2)
    etalq  ~ 0.066391  # Gao 2021 Table 2 'CV for IIV' on the Q1 row = 26.2% -> log(1 + 0.262^2)

    # Residual error. Gao 2021 Methods: 'The residual unexplained variability
    # ... was modeled with an additive error on the log-transformed
    # concentrations', which is proportional error in nlmixr2's linear space.
    # Wang 2023 Table 2 marks this row 'Prop * = 35.4%' with the footnote
    # '* The exponential residual error model was used, in which the estimate
    # is approximately equal to the proportional error model.'
    propSd <- 0.354; label("Proportional residual error (fraction)")  # Gao 2021 Table 2 'sigma 0.354' (RSE 4.3; SIR median 0.350, 95% CI 0.336-0.366)
  })

  model({
    # Individual PK parameters. Every clearance carries the 0.75 exponent and
    # every volume the 1.0 exponent on (WT/19); CL additionally carries the
    # linear serum-creatinine term.
    cl  <- exp(lcl + etalcl) * (WT / 19)^e_wt_cl * (1 + e_creat_cl * (CREAT - 26))
    vc  <- exp(lvc)          * (WT / 19)^e_wt_vc
    q   <- exp(lq  + etalq)  * (WT / 19)^e_wt_cl
    vp  <- exp(lvp)          * (WT / 19)^e_wt_vc
    q2  <- exp(lq2)          * (WT / 19)^e_wt_cl
    vp2 <- exp(lvp2)         * (WT / 19)^e_wt_vc

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Three-compartment disposition. Methotrexate is infused intravenously
    # straight into the central compartment, so there is no absorption phase
    # and no bioavailability term.
    d/dt(central)     <- -kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # Plasma methotrexate concentration in umol/L (amounts in umol, volumes L).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
