Gomes_2017_gentamicin_standard <- function() {
  description <- "One-compartment intravenous population PK model for gentamicin: the MWPharm 'Standard model' (Evers and Vinks 1995; van Lent-Evers 1999) as tabulated and externally validated in adults with infective endocarditis by Gomes 2017. Clearance is the sum of a metabolic arm scaled linearly by total body weight and a renal arm proportional to Cockcroft-Gault creatinine clearance; volume is a per-kg coefficient on the MWPharm fat-corrected lean body mass. Log-normal between-subject variability on all three parameters; residual error is the EMIT assay-error polynomial Gomes 2017 used."
  reference <- paste(
    "Gomes A, van der Wijk L, Proost JH, Sinha B, Touw DJ.",
    "Pharmacokinetic modeling of gentamicin in treatment of infective",
    "endocarditis: Model development and validation of existing models.",
    "PLoS One. 2017;12(5):e0177324. doi:10.1371/journal.pone.0177324.",
    "PMCID: PMC5419648. Parameters originate from Evers N, Vinks A.",
    "Development of gentamicin population models for ICU and non-ICU",
    "patients using a non parametric algorithm (NPEM). Ther Drug Monit.",
    "1995;17(4):424 (abstract), cited by Gomes 2017 as references 19 and 24.",
    sep = " "
  )
  vignette <- "Gomes_2017_gentamicin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters twice. (1) Linearly on the metabolic clearance arm as",
        "cl_nonren = CLm * (WT / 70) (Methods, 'CL = CLm*(BW/70) + fr*CLcr');",
        "70 kg is the unit denominator of the MWPharm CLm (L/h/70kgBW), not",
        "a cohort statistic. (2) Through the fat-corrected lean body mass",
        "LBMc = LBM + 0.4 * (WT - LBM) that scales the volume (see HT).",
        "Gomes 2017 endocarditis modeling cohort mean 76.2 kg, range 46-121 kg",
        "(Table 1); the source population of this model is not described."
      ),
      source_name = "BW"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Input to the lean-body-mass formula. The paper cites Chennavasin",
        "and Brater (reference 27) for LBM and Schwartz et al. (reference 28)",
        "for the fat correction but prints neither formula; both were",
        "back-solved from the per-patient LBM and LBMc columns of the S1 Data",
        "workbook (sheet 'The Hague', all 65 modeling-cohort patients) and",
        "reproduce every tabulated value to its 0.1 kg rounding:",
        "LBM = min(WT, 50 + 2.3 * (HT / 2.54 - 60)) for men and",
        "LBM = min(WT, 45.5 + 2.3 * (HT / 2.54 - 60)) for women (the Devine",
        "ideal body weight with height in inches, capped at total body",
        "weight), then LBMc = LBM + 0.4 * (WT - LBM). Gomes 2017 modeling cohort mean",
        "173.9 cm, range 149-193 cm (Table 1)."
      ),
      source_name = "Height"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Selects the intercept of the lean-body-mass formula only (50 kg for",
        "men, 45.5 kg for women). Sex has no separate effect on clearance or",
        "on the volume coefficient (its effect on creatinine clearance is",
        "already inside the supplied CRCL). S1 Data codes sex as 'M' / 'V'",
        "(Dutch 'vrouw' = female). Gomes 2017 modeling cohort 21 F / 44 M (Table 1)."
      ),
      source_name = "Gender"
    ),
    CRCL = list(
      description = paste(
        "Creatinine clearance by Cockcroft-Gault in ABSOLUTE mL/min, NOT",
        "BSA-normalized. Stored under the canonical CRCL column, which",
        "accepts raw mL/min when the source applies no BSA normalization."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The renal arm is fr * CLcr with CLcr in L/h; model() converts the",
        "mL/min column by 0.06 (60 min/h / 1000 mL/L). The absolute (not",
        "/1.73 m^2) scale is fixed by the paper's own Figure 1: its ICU-model",
        "typical curve for the 'standard patient' (male, 70 y, 76.2 kg,",
        "174 cm, CLcr 58 mL/min/1.73 m^2) decays at 0.156 1/h, which the",
        "ICU-model parameters reproduce only with the de-normalized",
        "58 * BSA / 1.73 = 64.0 mL/min (Du Bois BSA 1.91 m^2); the",
        "normalized 58 mL/min would give 0.142 1/h. The S1 Data workbook",
        "carries both CLcr (mL/min) and CLcr_n (/1.73 m^2) per patient.",
        "Gomes 2017 modeling cohort mean 64.3 mL/min, range 8.7-157.5 mL/min (Table 1).",
        "Dialysis patients were excluded, and fewer than 10% of the cohort had",
        "CLcr below 15 mL/min (Discussion)."
      ),
      source_name = "CLcr"
    )
  )

  compartmentData <- list(
    central = list(analyte = "gentamicin", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    disease_state = "Adult non-ICU (general ward) patients receiving gentamicin (source population of the MWPharm Standard model)",
    regions = "Netherlands",
    notes = paste(
      "Gomes 2017 does not describe the population this model was built",
      "from; it tabulates the parameters (Table 3) as one of the two",
      "MWPharm gentamicin models in Dutch clinical use and cites Evers",
      "and Vinks 1995 (an NPEM conference abstract) and van Lent-Evers 1999.",
      "Gomes 2017 validated it on 14 University Medical Center Groningen",
      "endocarditis patients / 30 serum samples (Table 4: MDPE -0.9%, MDAPE 4.8%),",
      "and used the ICU model for its Figure 1 and 2 simulations."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Gomes 2017 Table 3, row 'Standard model' (rounded) and the S2 Data
    # workbook, sheet 'Table 2 and 3', rows 'Standard Pop / SD' (unrounded).
    # Population mean and interindividual SD on the NATURAL scale; the
    # MWPharm parameters are log-normal (Methods, Phase 1), so the mean
    # is carried as the typical value and the SD converted by
    # omega^2 = log(1 + (SD / mean)^2). No covariances are reported, so the
    # block is diagonal.
    # ---------------------------------------------------------------------

    # Metabolic (non-renal) clearance per 70 kg total body weight.
    lcl_nonren <- log(0.2254)
    label("Metabolic clearance CLm at 70 kg total body weight (L/h/70 kg)") # Table 3 'Standard model' CLm mean 0.225, SD 0.113; S2 Data CLm 0.2254 SD 0.1127

    # fr, renal gentamicin clearance as a fraction of creatinine clearance
    # (both in L/h), on the log scale to carry its variability.
    le_crcl_cl_renal <- log(0.7272)
    label("Renal gentamicin clearance as a fraction fr of creatinine clearance (unitless)") # Table 3 'Standard model' fr mean 0.727, SD 0.471; S2 Data fr 0.7272 SD 0.4709

    # Volume coefficient per kg of fat-corrected lean body mass.
    lvd <- log(0.2729)
    label("Volume of distribution per kg fat-corrected lean body mass (L/kg LBMc)") # Table 3 'Standard model' Vd mean 0.273, SD 0.082 L/kgLBMc; S2 Data V 0.2729 SD 0.0822

    # MWPharm fat-distribution fraction in LBMc = LBM + fd * (WT - LBM),
    # back-solved from the Gomes 2017 S1 Data LBM / LBMc columns (median
    # 0.4000 across all patients with WT > LBM). Shared by all MWPharm
    # gentamicin models in Gomes 2017, which state every volume in
    # L/kgLBMc (Table 3); the paper's Figure 1 ICU-model peak is reproduced
    # with this value (see the vignette).
    e_wt_vd <- fixed(0.4)
    label("Fat-distribution fraction fd weighting (WT - LBM) inside LBMc (unitless)") # back-solved from S1 Data sheet 'The Hague' LBM / LBMc / Weight columns; Methods cite reference 28 for the fat correction

    # IIV: log(1 + CV^2) of the S2 Data natural-scale SDs.
    etalcl_nonren ~ 0.223144 # CLm SD 0.1127 / mean 0.2254 -> CV 50.0 %
    etale_crcl_cl_renal ~ 0.350181 # fr SD 0.4709 / mean 0.7272 -> CV 64.8 %
    etalvd ~ 0.086844 # Vd SD 0.0822 / mean 0.2729 -> CV 30.1 %

    # Residual error: the EMIT assay-error polynomial used for modeling
    # (Methods, Gentamicin assay): SD = 0.0766 + 0.0006 * C + 0.0064 * C^2
    # mg/L. Assay characteristics, not fitted, hence fixed (the S2 Data
    # bootstrap sheet holds S0, S1, S2 constant at these values).
    addSd <- fixed(0.0766)
    label("Assay SD polynomial intercept S0 (mg/L)") # Methods 'SD = 0.0766 + 0.0006 C + 0.0064 C2'
    propSd <- fixed(0.0006)
    label("Assay SD polynomial linear coefficient S1 (fraction)") # Methods 'SD = 0.0766 + 0.0006 C + 0.0064 C2'
    quadSd <- fixed(0.0064)
    label("Assay SD polynomial quadratic coefficient S2 (L/mg)") # Methods 'SD = 0.0766 + 0.0006 C + 0.0064 C2'
  })

  model({
    # Lean body mass (kg): Devine ideal body weight with height in inches,
    # capped at total body weight; then the MWPharm fat correction. Both
    # back-solved from the S1 Data per-patient LBM / LBMc columns.
    ibw <- (50 - 4.5 * SEXF) + 2.3 * (HT / 2.54 - 60)
    lbm <- ibw
    if (WT < ibw) lbm <- WT
    lbmc <- lbm + e_wt_vd * (WT - lbm)

    # Cockcroft-Gault creatinine clearance from mL/min to L/h.
    clcr <- 0.06 * CRCL

    # Methods: CL = CLm * (BW / 70) + fr * CLcr.
    cl_nonren <- exp(lcl_nonren + etalcl_nonren) * (WT / 70)
    fr <- exp(le_crcl_cl_renal + etale_crcl_cl_renal)
    cl_renal <- fr * clcr
    cl <- cl_nonren + cl_renal

    # Volume = Vd (L/kg LBMc) * LBMc.
    vd <- exp(lvd + etalvd)
    vc <- vd * lbmc

    kel <- cl / vc

    # One compartment, first-order elimination (Methods, Phase 1).
    d/dt(central) <- -kel * central

    Cc <- central / vc
    sdCc <- addSd + propSd * Cc + quadSd * Cc^2
    Cc ~ add(sdCc)
  })
}
