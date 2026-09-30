Zhang_2021_ceftolozane <- function() {
  description <- paste(
    "Two-compartment population PK model for ceftolozane given as a 1-hour",
    "intravenous infusion (alone or with tazobactam), fitted to 8,330 plasma",
    "concentrations from 968 adults pooled across 16 studies, including 305",
    "ventilated patients with hospital-acquired or ventilator-associated",
    "pneumonia from the phase 3 ASPECT-NP trial (Zhang 2021). Elimination is",
    "first order. Baseline Cockcroft-Gault creatinine clearance acts on",
    "clearance through a power function, with an additional end-stage renal",
    "disease factor; body weight acts on both volumes through estimated power",
    "exponents; infection type (cUTI, cIAI, pneumonia, other infection) shifts",
    "clearance and central volume against a healthy-participant reference.",
    "Correlated IIV on CL and Vc, independent IIV on Vp. A hypothetical",
    "epithelial-lining-fluid (ELF) link compartment without mass transfer",
    "(influx K1E from the central amount, elimination KE0, ELF volume equal to",
    "Vc) predicts ELF concentrations; pneumonia lowers K1E and KE0 by one",
    "shared factor and IIV on K1E differs between healthy participants and",
    "pneumonia patients. Companion model to Zhang_2021_tazobactam, fitted",
    "separately in the same paper."
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

  # Doses are in mg and the paper reports concentrations in ug/mL, which is
  # numerically identical to mg/L, so amount(mg) / volume(L) lands directly in
  # the reported unit. The ELF state is a hypothetical link compartment: its
  # amount is scaled by Vc (Supplemental Figure S1: 'V3 ... assumed to equal
  # V1') and is not removed from the central compartment.
  compartmentData <- list(
    central = list(analyte = "ceftolozane", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ceftolozane", units = "mg", specimen = "plasma", verified = TRUE),
    elf = list(analyte = "ceftolozane", units = "mg", specimen = "epithelial lining fluid", verified = TRUE)
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
        "Table 2 row 'Exponent of (CrCL/100) for CL' = 0.701 (RSE 4.32%);",
        "equation 'CL = 4.84 x (CrCL/100)^0.701 x ...'. Reference 100 mL/min.",
        "Results: renal function 'estimated as CrCL (calculated using the",
        "Cockcroft-Gault formula)'. Baseline value, time-fixed (Discussion:",
        "'covariates were assumed to remain at their baseline values').",
        "Table 1: mean 109.9 mL/min (SD 56.7), range 6.3-531.3 overall;",
        "ASPECT-NP mean 124.1 (range 14.9-531.3)."
      ),
      source_name = "CrCL"
    ),
    WT = list(
      description = "Actual total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table 2 rows 'Exponent of (WTKG/70) for Vc' = 0.684 (RSE 9.69%) and",
        "'Exponent of (WTKG/70) for Vp' = 0.484 (RSE 16.3%). Reference 70 kg",
        "(Results: 'relative to a 70-kg adult'). Not a covariate on CL in",
        "the final model. Table 1: mean 74.7 kg, range 33.5-173.0."
      ),
      source_name = "WTKG"
    ),
    RENALIMP_ESRD = list(
      description = "End-stage renal disease indicator (1 = ESRD)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no ESRD)",
      notes = paste(
        "Table 2 rows 'Fold-change in CL for ESRD' = 0.320 (RSE 10.9%) and",
        "'Fold-change in Vc for ESRD' = 1.30 (RSE 11.7%). MULTIPLIES the",
        "continuous CrCL power term rather than replacing it -- Results:",
        "'ceftolozane CL was reduced by an additional 68%, in addition to the",
        "decrease associated with reduced CrCL'. The ESRD stratum is the 6",
        "patients with ESRD on hemodialysis whose PK was collected after (not",
        "during) the dialysis session (Methods); the effect of dialysis",
        "itself is not in this model. Supply the column from the diagnosis,",
        "not by thresholding CRCL."
      ),
      source_name = "ESRD"
    ),
    DIS_CUTI = list(
      description = "Complicated urinary tract infection cohort indicator (1 = cUTI)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy participant, including renal-impairment subjects without infection; all infection indicators 0)",
      notes = paste(
        "Table 2 rows 'Fold-change in CL for cUTI' = 1.18 (RSE 2.99%) and",
        "'Fold-change in Vc for cUTI' = 1.25 (RSE 3.52%), entered as",
        "1.18^cUTI and 1.25^cUTI. One of the paper's mutually exclusive",
        "infection types (cUTI, cIAI, pneumonia, other infection); the",
        "all-zero reference is the 'Healthy' row of Table 1, whose footnote",
        "says it 'Included patients with renal impairment without",
        "infection'. 176 of 968 subjects (18.2%)."
      ),
      source_name = "cUTI"
    ),
    DIS_CIAI = list(
      description = "Complicated intra-abdominal infection cohort indicator (1 = cIAI)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy participant; all infection indicators 0)",
      notes = paste(
        "Table 2 rows 'Fold-change in CL for cIAI' = 1.43 (RSE 2.99%) and",
        "'Fold-change in Vc for cIAI' = 1.59 (RSE 4.97%). Mutually exclusive",
        "with the other infection indicators. 174 of 968 subjects (18.0%)."
      ),
      source_name = "cIAI"
    ),
    DIS_PNEUMONIA = list(
      description = "Pneumonia infection-type indicator (1 = pneumonia)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy participant; all infection indicators 0)",
      notes = paste(
        "Table 2 row 'Fold-change in Vc for pneumonia' = 2.00 (RSE 3.29%);",
        "no pneumonia term on CL or Vp. Table 3: pneumonia also scales BOTH",
        "ELF rate constants by 0.0339 (equations 'TVK1E = 0.808 x",
        "0.0339^Pneu', 'TVKE0 = 1.56 x 0.0339^Pneu') and switches the K1E",
        "IIV from the healthy-participant to the pneumonia-patient variance.",
        "The paper's pneumonia stratum pools the 305 ventilated HABP/VABP",
        "patients of ASPECT-NP with the critically ill patients with",
        "confirmed or suspected pneumonia of study CXA-ICU-14-01 (331 of 968",
        "subjects, 34.2%) and never separates hospital-acquired from",
        "ventilator-associated pneumonia, so the unqualified pneumonia flag",
        "is used rather than DIS_HABP / DIS_VABP."
      ),
      source_name = "Pneumonia"
    ),
    DIS_OTHER_INFECT = list(
      description = "Other-infection cohort indicator (1 = the paper's 'other infection' type)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy participant; all infection indicators 0)",
      notes = paste(
        "Table 2 row 'Fold-change in Vc for other infections' = 2.14 (RSE",
        "14.8%); equation term '2.14^OINF'. Methods: 'The patients",
        "categorized as other infection were from the augmented renal",
        "clearance cohort of study MK-7625A-007 ... in which CrCL was",
        "calculated using the Cockcroft-Gault equation'; Table 2 note:",
        "'Other infections included critically ill adult patients from study",
        "MK-7625A-007'. 10 of 968 subjects (1.0%). This is a residual",
        "infection-type category of the paper's mutually exclusive scheme,",
        "NOT an augmented-renal-clearance flag: ASPECT-NP patients with CrCL",
        "> 150 mL/min are pneumonia patients here, not other-infection."
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
        "Race was tested as white / Japanese / other (Methods) and rejected.",
        "Results: 'No significant differences in the PK of ceftolozane were",
        "observed based on race.' Japanese 21.7%, other 12.2% (Table 1)."
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
    n_subjects = 968L,
    n_studies = 16L,
    n_concentrations = "8330 plasma; 47 ELF (one bronchoalveolar lavage sample per participant)",
    age_range = "18-98 years",
    age_median = "mean 53.2 years (SD 19.7); ASPECT-NP mean 59.8",
    weight_range = "33.5-173.0 kg",
    weight_median = "mean 74.7 kg (SD 17.8); ASPECT-NP mean 81.5",
    sex_female_pct = 39.8,
    race_ethnicity = c(White = 66.1, Japanese = 21.7, Other = 12.2),
    disease_state = paste(
      "Pooled healthy participants (including renal impairment without",
      "infection, 28.6%), cUTI (18.2%), cIAI (18.0%), pneumonia (34.2%;",
      "305 ventilated HABP/VABP patients from ASPECT-NP plus critically ill",
      "patients with confirmed or suspected pneumonia) and other infection",
      "(1.0%, critically ill augmented-renal-clearance cohort) (Table 1)."
    ),
    dose_range = paste(
      "Ceftolozane 0.25-3.0 g as 1-hour intravenous infusions, alone or with",
      "tazobactam, single dose or every 6, 8, 12 or 24 h. ASPECT-NP: 2 g",
      "ceftolozane (3 g C/T) q8h for CrCL > 50 mL/min, 1 g for CrCL 30-50,",
      "0.5 g for CrCL 15-29 (Introduction; Supplemental Table S1)."
    ),
    renal_function = paste(
      "Cockcroft-Gault CrCL mean 109.9 mL/min, range 6.3-531.3; > 50 mL/min",
      "91.3%, 30-50 6.2%, 15-29 1.8%, < 15 0.7% (Table 1). Six ESRD patients",
      "on hemodialysis contribute post-dialysis data."
    ),
    elf_subpopulation = paste(
      "ELF data from 2 phase 1 studies: CXA-ELF-10-03 (healthy volunteers)",
      "and CXA-ICU-14-01 (critically ill, confirmed or suspected pneumonia);",
      "51 participants, 53.2% healthy, sampled once at 1, 2, 4, 6 or 8 h",
      "after the start of infusion (Methods, Results)."
    ),
    regions = "Multinational (16 studies; Supplemental Tables S1-S2)",
    notes = paste(
      "NONMEM 7.3.0 (Cognigen). Postdose samples below the LLOQ omitted",
      "(5.2% of ceftolozane plasma samples). Outliers (|CWRES| >= 6) excluded",
      "from the final model. The ELF link parameters were estimated in a",
      "joint plasma-ELF fit whose plasma component used the pooled data",
      "without the phase 3 studies; the paper reports only the ELF",
      "parameters of that fit and states that 'The final ceftolozane and",
      "tazobactam ELF disposition models were the plasma models described",
      "above with a hypothetical ELF compartment linked to the plasma",
      "compartment', which is how the two are combined here."
    )
  )

  ini({
    # =====================================================================
    # PLASMA STRUCTURAL PARAMETERS -- Table 2, ceftolozane 'Typical Value'.
    # Reference subject: CrCL 100 mL/min, 70 kg, healthy, no ESRD.
    # =====================================================================
    lcl <- log(4.84); label("Ceftolozane clearance at CrCL 100 mL/min, healthy (L/h)") # Table 2: Systemic CL = 4.84 L/h (RSE 1.66%)
    lvc <- log(9.23); label("Ceftolozane central volume at 70 kg, healthy (L)") # Table 2: Vc = 9.23 L (RSE 1.83%)
    lq <- log(3.13); label("Ceftolozane intercompartmental clearance (L/h)") # Table 2: Q = 3.13 L/h (RSE 6.93%)
    lvp <- log(4.78); label("Ceftolozane peripheral volume at 70 kg (L)") # Table 2: Vp = 4.78 L (RSE 3.24%)

    # =====================================================================
    # COVARIATE EFFECTS -- Table 2 and the equations beneath it:
    #   CL = 4.84 x (CrCL/100)^0.701 x 1.18^cUTI x 1.43^cIAI x 0.32^ESRD
    #   Vc = 9.23 x (WTKG/70)^0.684 x 1.25^cUTI x 1.59^cIAI x 2.00^Pneumonia
    #        x 2.14^OINF x 1.30^ESRD
    #   Vp = 4.78 x (WTKG/70)^0.484
    # Categorical effects are printed as fold-changes raised to the 0/1
    # indicator; stored here as their logs and exponentiated in model().
    # =====================================================================
    e_crcl_cl <- 0.701; label("Power exponent of CrCL/100 on ceftolozane CL (unitless)") # Table 2: 0.701 (RSE 4.32%)
    e_wt_vc <- 0.684; label("Power exponent of WT/70 on ceftolozane Vc (unitless)") # Table 2: 0.684 (RSE 9.69%)
    e_wt_vp <- 0.484; label("Power exponent of WT/70 on ceftolozane Vp (unitless)") # Table 2: 0.484 (RSE 16.3%)
    e_esrd_cl <- log(0.320); label("ESRD effect on ceftolozane CL (log fold-change)") # Table 2: 0.320 (RSE 10.9%)
    e_cuti_cl <- log(1.18); label("cUTI effect on ceftolozane CL (log fold-change)") # Table 2: 1.18 (RSE 2.99%)
    e_ciai_cl <- log(1.43); label("cIAI effect on ceftolozane CL (log fold-change)") # Table 2: 1.43 (RSE 2.99%)
    e_cuti_vc <- log(1.25); label("cUTI effect on ceftolozane Vc (log fold-change)") # Table 2: 1.25 (RSE 3.52%)
    e_ciai_vc <- log(1.59); label("cIAI effect on ceftolozane Vc (log fold-change)") # Table 2: 1.59 (RSE 4.97%)
    e_pneumonia_vc <- log(2.00); label("Pneumonia effect on ceftolozane Vc (log fold-change)") # Table 2: 2.00 (RSE 3.29%)
    e_other_infect_vc <- log(2.14); label("Other-infection effect on ceftolozane Vc (log fold-change)") # Table 2: 2.14 (RSE 14.8%)
    e_esrd_vc <- log(1.30); label("ESRD effect on ceftolozane Vc (log fold-change)") # Table 2: 1.30 (RSE 11.7%)

    # =====================================================================
    # ELF LINK PARAMETERS -- Table 3, ceftolozane. Reference: healthy.
    #   TVK1E = 0.808 x 0.0339^Pneu ; TVKE0 = 1.56 x 0.0339^Pneu
    # Table 3 prints the shared pneumonia factor as 0.034 (RSE 44.0%); the
    # equations beneath it carry the extra digit, 0.0339, used here.
    # =====================================================================
    lk_central_elf <- log(0.808); label("Ceftolozane plasma-to-ELF rate constant K1E, healthy (1/h)") # Table 3: K1E = 0.808 /h (RSE 11.6%)
    lke0 <- log(1.56); label("Ceftolozane ELF elimination rate constant KE0, healthy (1/h)") # Table 3: KE0 = 1.56 /h (RSE 8.97%)
    e_pneumonia_elf <- log(0.0339); label("Pneumonia effect on ceftolozane K1E and KE0, shared (log fold-change)") # Table 3 equations: 0.0339^Pneu on both

    # =====================================================================
    # BETWEEN-SUBJECT VARIABILITY -- Table 2 / Table 3 'Magnitude'.
    # SCALE: the %CV column is sqrt(omega^2) x 100, not
    # sqrt(exp(omega^2) - 1). Table 2 footnote c prints the CL-Vc
    # correlation as r = 0.474 alongside cov = 0.073; with
    # omega^2 = CV^2 that gives 0.073 / sqrt(0.361^2 x 0.429^2) = 0.471,
    # whereas omega^2 = log(1 + CV^2) gives 0.507, outside the rounding of
    # every printed digit. omega^2 = CV^2 throughout.
    # =====================================================================
    etalcl + etalvc ~ c(
      0.130321,
      0.073, 0.184041
    ) # Table 2: CL 36.1% CV -> 0.361^2; Vc 42.9% CV -> 0.429^2; cov 0.073 (RSE 17.1%)
    etalvp ~ 0.022801 # Table 2: Vp 15.1% CV -> 0.151^2 (RSE 57.8%, shrinkage 55.8%)
    etalk_central_elf ~ 0.156816 # Table 3: K1E IIV, healthy participants 39.6% CV -> 0.396^2
    etalk_central_elf_pneumonia ~ 0.659344 # Table 3: K1E IIV, pneumonia patients 81.2% CV -> 0.812^2

    # =====================================================================
    # RESIDUAL VARIABILITY -- combined proportional + additive on the
    # variance scale. Table 2 footnote d: 'RV (%CV) = sqrt(F^2 x 0.0248 +
    # 0.00984) / F x 100', i.e. the printed 0.025 / 0.010 are sigma^2. That
    # reproduces the printed '100%-15.7% CV' at F = 0.1 and 200 ug/mL and the
    # '0.099 SD' additive term. Table 3 footnote d gives the ELF-model pair
    # 0.0249 / 0.00834 ('92.7%-15.8% CV' at F = 0.1-100 ug/mL).
    # =====================================================================
    propSd <- 0.157480; label("Plasma proportional residual SD (fraction)") # Table 2: sigma^2 prop = 0.0248 -> sqrt = 0.1575 (15.7% CV)
    addSd <- 0.099197; label("Plasma additive residual SD (ug/mL)") # Table 2: sigma^2 add = 0.00984 -> sqrt = 0.0992 ('0.099 SD')
    propSd_Celf <- 0.157797; label("ELF proportional residual SD (fraction)") # Table 3: sigma^2 prop = 0.0249 -> sqrt = 0.1578 (15.8% CV)
    addSd_Celf <- 0.091324; label("ELF additive residual SD (ug/mL)") # Table 3: sigma^2 add = 0.00834 -> sqrt = 0.0913
  })

  model({
    # 1. Individual plasma parameters (Table 2 equations).
    cl <- exp(
      lcl + etalcl + e_cuti_cl * DIS_CUTI + e_ciai_cl * DIS_CIAI +
        e_esrd_cl * RENALIMP_ESRD
    ) * (CRCL / 100)^e_crcl_cl
    vc <- exp(
      lvc + etalvc + e_cuti_vc * DIS_CUTI + e_ciai_vc * DIS_CIAI +
        e_pneumonia_vc * DIS_PNEUMONIA + e_other_infect_vc * DIS_OTHER_INFECT +
        e_esrd_vc * RENALIMP_ESRD
    ) * (WT / 70)^e_wt_vc
    q <- exp(lq)
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vp

    # 2. ELF link parameters (Table 3 equations). K1E carries a separate
    #    IIV for healthy participants and pneumonia patients (Methods:
    #    'estimating IIV with K1E separately for healthy participants and
    #    patients with pneumonia'); KE0 has no IIV ('NE').
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

    # 4. ODEs. Zero-order 1-hour infusion into central (event table).
    #    Supplemental Figure S1 model 3: dA3/dt = K1E x A1 - KE0 x A3, with
    #    no A3 term in dA1/dt -- the link compartment does not deplete
    #    plasma ('without mass transfer', Methods).
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(elf) <- k_central_elf * central - ke0 * elf

    # 5. Observations. ELF volume V3 equals V1 (Supplemental Figure S1).
    Cc <- central / vc
    Celf <- elf / vc
    Cc ~ add(addSd) + prop(propSd)
    Celf ~ add(addSd_Celf) + prop(propSd_Celf)
  })
}
