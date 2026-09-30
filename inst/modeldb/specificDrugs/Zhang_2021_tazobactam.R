Zhang_2021_tazobactam <- function() {
  description <- paste(
    "Two-compartment population PK model for tazobactam given as a 1-hour",
    "intravenous infusion (alone or with ceftolozane), fitted to 5,679 plasma",
    "concentrations from 835 adults pooled across 16 studies, including 305",
    "ventilated patients with hospital-acquired or ventilator-associated",
    "pneumonia from the phase 3 ASPECT-NP trial (Zhang 2021). Elimination is",
    "first order. Baseline Cockcroft-Gault creatinine clearance acts on",
    "clearance through a power function, with an additional end-stage renal",
    "disease factor on clearance and central volume; body weight acts on both",
    "volumes through estimated power exponents; infection type (cIAI,",
    "pneumonia, other infection on Vc; cUTI, cIAI, pneumonia on Vp) shifts",
    "the volumes against a healthy-participant reference. Correlated IIV on",
    "CL and Vc, independent IIV on Vp. A hypothetical epithelial-lining-fluid",
    "(ELF) link compartment without mass transfer (influx K1E from the",
    "central amount, elimination KE0, ELF volume equal to Vc) predicts ELF",
    "concentrations; pneumonia lowers K1E and KE0 by one shared factor and",
    "IIV on K1E differs between healthy participants and pneumonia patients.",
    "Companion model to Zhang_2021_ceftolozane, fitted separately in the",
    "same paper."
  )
  reference <- paste(
    "Zhang Z, Patel YT, Fiedler-Kelly J, Feng HP, Bruno CJ, Gao W.",
    "Population Pharmacokinetic Analysis for Plasma and Epithelial Lining",
    "Fluid Ceftolozane/Tazobactam Concentrations in Patients With Ventilated",
    "Nosocomial Pneumonia. J Clin Pharmacol. 2021;61(2):254-268.",
    "doi:10.1002/jcph.1733.",
    "Plasma fixed effects, IIV and residual error: Table 2 and the",
    "parameter-covariate equations beneath it. ELF link parameters and ELF",
    "residual error: Table 3 and its equations. ELF model structure:",
    "Supplemental Figure S1, model 3.",
    sep = " "
  )
  vignette <- "Zhang_2021_ceftolozane_tazobactam"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Doses in mg, concentrations in ug/mL (= mg/L). The ELF state is a
  # hypothetical link compartment scaled by Vc (Supplemental Figure S1: 'V3
  # ... assumed to equal V1') that does not deplete the central compartment.
  compartmentData <- list(
    central = list(analyte = "tazobactam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tazobactam", units = "mg", specimen = "plasma", verified = TRUE),
    elf = list(analyte = "tazobactam", units = "mg", specimen = "epithelial lining fluid", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Baseline creatinine clearance estimated by the Cockcroft-Gault",
        "formula; raw mL/min, NOT body-surface-area normalized"
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table 2 row 'Exponent of (CrCL/100) for CL' = 0.623 (RSE 6.42%);",
        "equation 'CL = 16.6 x (CrCL/100)^0.623 x 0.626^ESRD'. Reference",
        "100 mL/min. Cockcroft-Gault, baseline, time-fixed. Table 1: mean",
        "111.9 mL/min (SD 59.5), range 6.3-531.3."
      ),
      source_name = "CrCL"
    ),
    WT = list(
      description = "Actual total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table 2 rows 'Exponent of (WTKG/70) for Vc' = 0.629 (RSE 10.7%) and",
        "'Exponent of (WTKG/70) for Vp' = 0.530 (RSE 14.4%). Reference 70 kg.",
        "Table 1: mean 74.5 kg, range 33.5-150.1."
      ),
      source_name = "WTKG"
    ),
    RENALIMP_ESRD = list(
      description = "End-stage renal disease indicator (1 = ESRD)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no ESRD)",
      notes = paste(
        "Table 2 rows 'Fold-change in CL for ESRD' = 0.626 (RSE 14.6%) and",
        "'Fold-change in Vc for ESRD' = 0.749 (RSE 11.0%); Results: 'ESRD did",
        "not significantly impact Vp'. Multiplies the continuous CrCL power",
        "term ('reduced by an additional 37.4%, in addition to the decrease",
        "associated with reduced CrCL'). The stratum is 6 hemodialysis",
        "patients sampled after dialysis; supply from the diagnosis, not by",
        "thresholding CRCL."
      ),
      source_name = "ESRD"
    ),
    DIS_CUTI = list(
      description = "Complicated urinary tract infection cohort indicator (1 = cUTI)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy participant, including renal-impairment subjects without infection; all infection indicators 0)",
      notes = paste(
        "Table 2 row 'Fold-change in Vp for cUTI' = 1.25 (RSE 4.27%); the",
        "tazobactam CL and Vc cUTI rows are '-' (not in the model). 103 of 835",
        "subjects (12.3%)."
      ),
      source_name = "cUTI"
    ),
    DIS_CIAI = list(
      description = "Complicated intra-abdominal infection cohort indicator (1 = cIAI)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy participant; all infection indicators 0)",
      notes = paste(
        "Table 2 rows 'Fold-change in Vc for cIAI' = 1.49 (RSE 3.70%) and",
        "'Fold-change in Vp for cIAI' = 1.34 (RSE 4.57%). 174 of 835 (20.8%)."
      ),
      source_name = "cIAI"
    ),
    DIS_PNEUMONIA = list(
      description = "Pneumonia infection-type indicator (1 = pneumonia)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy participant; all infection indicators 0)",
      notes = paste(
        "Table 2 rows 'Fold-change in Vc for pneumonia' = 2.17 (RSE 4.46%)",
        "and 'Fold-change in Vp for pneumonia' = 2.06 (RSE 8.89%). Table 3:",
        "pneumonia also scales BOTH ELF rate constants by 0.479 ('TVK1E =",
        "0.262 x 0.479^Pneu', 'TVKE0 = 0.691 x 0.479^Pneu') and switches the",
        "K1E IIV to the pneumonia-patient variance. Pools ASPECT-NP HABP/VABP",
        "and CXA-ICU-14-01 critically ill pneumonia patients (331 of 835,",
        "39.6%) without separating HABP from VABP."
      ),
      source_name = "Pneumonia"
    ),
    DIS_OTHER_INFECT = list(
      description = "Other-infection cohort indicator (1 = the paper's 'other infection' type)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy participant; all infection indicators 0)",
      notes = paste(
        "Table 2 row 'Fold-change in Vc for other infections' = 2.49 (RSE",
        "16.5%); equation term '2.49^OINF'. The 10 critically ill patients of",
        "the augmented-renal-clearance cohort of study MK-7625A-007 (1.2% of",
        "835). A residual infection-type category, not an ARC flag."
      ),
      source_name = "OINF"
    )
  )

  covariatesDataExcluded <- list(
    RACE_JAPANESE = list(
      description = "Japanese race indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Race was tested as white / Japanese / other and rejected. Results:",
        "'No significant differences in the PK of tazobactam were observed",
        "based on race.' Japanese 25.1%, other 12.6% (Table 1)."
      )
    ),
    RACE_OTHER = list(
      description = "Other (non-white, non-Japanese) race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened with RACE_JAPANESE and rejected; see that entry."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 835L,
    n_studies = 16L,
    n_concentrations = "5679 plasma; 42 ELF (one bronchoalveolar lavage sample per participant)",
    age_range = "18-98 years",
    age_median = "mean 53.9 years (SD 19.6); ASPECT-NP mean 59.8",
    weight_range = "33.5-150.1 kg",
    weight_median = "mean 74.5 kg (SD 18.0); ASPECT-NP mean 81.4",
    sex_female_pct = 38.9,
    race_ethnicity = c(White = 62.3, Japanese = 25.1, Other = 12.6),
    disease_state = paste(
      "Pooled healthy participants (including renal impairment without",
      "infection, 26.0%), cUTI (12.3%), cIAI (20.8%), pneumonia (39.6%) and",
      "other infection (1.2%) (Table 1)."
    ),
    dose_range = paste(
      "Tazobactam 0.25-1.5 g as 1-hour intravenous infusions, alone or with",
      "ceftolozane. ASPECT-NP: 1 g tazobactam (3 g C/T) q8h for CrCL > 50",
      "mL/min, 0.5 g for CrCL 30-50, 0.25 g for CrCL 15-29."
    ),
    renal_function = paste(
      "Cockcroft-Gault CrCL mean 111.9 mL/min, range 6.3-531.3; > 50 mL/min",
      "89.9%, 30-50 7.2%, 15-29 2.0%, < 15 0.8% (Table 1)."
    ),
    elf_subpopulation = paste(
      "ELF from CXA-ELF-10-03 (healthy) and CXA-ICU-14-01 (critically ill",
      "pneumonia); 42 quantifiable ELF observations from 42 participants."
    ),
    regions = "Multinational (16 studies; Supplemental Tables S1-S2)",
    notes = paste(
      "NONMEM 7.3.0 (Cognigen). 23.3% of postdose tazobactam plasma samples",
      "were below the LLOQ and omitted. The ELF link parameters come from a",
      "joint plasma-ELF fit whose plasma component is not reported; the",
      "paper states the final ELF models 'were the plasma models described",
      "above with a hypothetical ELF compartment linked to the plasma",
      "compartment', which is how the two are combined here."
    )
  )

  ini({
    # =====================================================================
    # PLASMA STRUCTURAL PARAMETERS -- Table 2, tazobactam 'Typical Value'.
    # Reference subject: CrCL 100 mL/min, 70 kg, healthy, no ESRD.
    # =====================================================================
    lcl <- log(16.6); label("Tazobactam clearance at CrCL 100 mL/min, no ESRD (L/h)") # Table 2: Systemic CL = 16.6 L/h (RSE 1.80%)
    lvc <- log(13.1); label("Tazobactam central volume at 70 kg, healthy (L)") # Table 2: Vc = 13.1 L (RSE 1.96%)
    lq <- log(4.05); label("Tazobactam intercompartmental clearance (L/h)") # Table 2: Q = 4.05 L/h (RSE 3.81%)
    lvp <- log(4.89); label("Tazobactam peripheral volume at 70 kg, healthy (L)") # Table 2: Vp = 4.89 L (RSE 2.59%)

    # =====================================================================
    # COVARIATE EFFECTS -- Table 2 equations:
    #   CL = 16.6 x (CrCL/100)^0.623 x 0.626^ESRD
    #   Vc = 13.1 x (WTKG/70)^0.629 x 1.49^cIAI x 2.17^Pneumonia x 2.49^OINF
    #        x 0.749^ESRD
    #   Vp = 4.89 x (WTKG/70)^0.530 x 1.25^cUTI x 1.34^cIAI x 2.06^Pneumonia
    # (the typeset Vp equation drops the 'x' before 2.06^Pneumonia; Table 2
    # lists the 2.06 row, and the Conclusions give 'Vp for tazobactam ...
    # 106% higher'.)
    # =====================================================================
    e_crcl_cl <- 0.623; label("Power exponent of CrCL/100 on tazobactam CL (unitless)") # Table 2: 0.623 (RSE 6.42%)
    e_wt_vc <- 0.629; label("Power exponent of WT/70 on tazobactam Vc (unitless)") # Table 2: 0.629 (RSE 10.7%)
    e_wt_vp <- 0.530; label("Power exponent of WT/70 on tazobactam Vp (unitless)") # Table 2: 0.530 (RSE 14.4%)
    e_esrd_cl <- log(0.626); label("ESRD effect on tazobactam CL (log fold-change)") # Table 2: 0.626 (RSE 14.6%)
    e_ciai_vc <- log(1.49); label("cIAI effect on tazobactam Vc (log fold-change)") # Table 2: 1.49 (RSE 3.70%)
    e_pneumonia_vc <- log(2.17); label("Pneumonia effect on tazobactam Vc (log fold-change)") # Table 2: 2.17 (RSE 4.46%)
    e_other_infect_vc <- log(2.49); label("Other-infection effect on tazobactam Vc (log fold-change)") # Table 2: 2.49 (RSE 16.5%)
    e_esrd_vc <- log(0.749); label("ESRD effect on tazobactam Vc (log fold-change)") # Table 2: 0.749 (RSE 11.0%)
    e_cuti_vp <- log(1.25); label("cUTI effect on tazobactam Vp (log fold-change)") # Table 2: 1.25 (RSE 4.27%)
    e_ciai_vp <- log(1.34); label("cIAI effect on tazobactam Vp (log fold-change)") # Table 2: 1.34 (RSE 4.57%)
    e_pneumonia_vp <- log(2.06); label("Pneumonia effect on tazobactam Vp (log fold-change)") # Table 2: 2.06 (RSE 8.89%)

    # =====================================================================
    # ELF LINK PARAMETERS -- Table 3, tazobactam. Reference: healthy.
    #   TVK1E = 0.262 x 0.479^Pneu ; TVKE0 = 0.691 x 0.479^Pneu
    # =====================================================================
    lk_central_elf <- log(0.262); label("Tazobactam plasma-to-ELF rate constant K1E, healthy (1/h)") # Table 3: K1E = 0.262 /h (RSE 30.6%)
    lke0 <- log(0.691); label("Tazobactam ELF elimination rate constant KE0, healthy (1/h)") # Table 3: KE0 = 0.691 /h (RSE 25.1%)
    e_pneumonia_elf <- log(0.479); label("Pneumonia effect on tazobactam K1E and KE0, shared (log fold-change)") # Table 3: 0.479 (RSE 44.0%)

    # =====================================================================
    # BETWEEN-SUBJECT VARIABILITY -- omega^2 = CV^2 (the %CV column is
    # sqrt(omega^2) x 100; see Zhang_2021_ceftolozane, where the printed
    # CL-Vc correlation r = 0.474 fixes the scale for the paper's tables).
    # =====================================================================
    etalcl + etalvc ~ c(
      0.281961,
      0.071, 0.149769
    ) # Table 2: CL 53.1% CV -> 0.531^2; Vc 38.7% CV -> 0.387^2; cov 0.071 (RSE 33.3%)
    etalvp ~ 0.037636 # Table 2: Vp 19.4% CV -> 0.194^2 (RSE 54.7%)
    etalk_central_elf ~ 0.429025 # Table 3: K1E IIV, healthy participants 65.5% CV -> 0.655^2
    etalk_central_elf_pneumonia ~ 0.712336 # Table 3: K1E IIV, pneumonia patients 84.4% CV -> 0.844^2

    # =====================================================================
    # RESIDUAL VARIABILITY -- proportional only; printed values are sigma^2
    # (Table 2: 0.081 -> '28.5% CV' = sqrt(0.081); Table 3: 0.055 ->
    # '23.4% CV' = sqrt(0.055)).
    # =====================================================================
    propSd <- 0.284605; label("Plasma proportional residual SD (fraction)") # Table 2: sigma^2 prop = 0.081 -> 28.5% CV
    propSd_Celf <- 0.234521; label("ELF proportional residual SD (fraction)") # Table 3: sigma^2 prop = 0.055 -> 23.4% CV
  })

  model({
    # 1. Individual plasma parameters (Table 2 equations).
    cl <- exp(lcl + etalcl + e_esrd_cl * RENALIMP_ESRD) * (CRCL / 100)^e_crcl_cl
    vc <- exp(
      lvc + etalvc + e_ciai_vc * DIS_CIAI + e_pneumonia_vc * DIS_PNEUMONIA +
        e_other_infect_vc * DIS_OTHER_INFECT + e_esrd_vc * RENALIMP_ESRD
    ) * (WT / 70)^e_wt_vc
    q <- exp(lq)
    vp <- exp(
      lvp + etalvp + e_cuti_vp * DIS_CUTI + e_ciai_vp * DIS_CIAI +
        e_pneumonia_vp * DIS_PNEUMONIA
    ) * (WT / 70)^e_wt_vp

    # 2. ELF link parameters (Table 3 equations), separate K1E IIV for
    #    healthy participants and pneumonia patients; no IIV on KE0.
    k_central_elf <- exp(
      lk_central_elf + e_pneumonia_elf * DIS_PNEUMONIA +
        etalk_central_elf * (1 - DIS_PNEUMONIA) +
        etalk_central_elf_pneumonia * DIS_PNEUMONIA
    )
    ke0 <- exp(lke0 + e_pneumonia_elf * DIS_PNEUMONIA)

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODEs (Supplemental Figure S1 model 3; link compartment does not
    #    deplete plasma).
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(elf) <- k_central_elf * central - ke0 * elf

    # 5. Observations; ELF volume equals Vc.
    Cc <- central / vc
    Celf <- elf / vc
    Cc ~ prop(propSd)
    Celf ~ prop(propSd_Celf)
  })
}
