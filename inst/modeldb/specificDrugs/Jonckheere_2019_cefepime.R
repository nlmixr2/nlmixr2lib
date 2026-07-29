Jonckheere_2019_cefepime <- function() {
  description <- "Two-compartment population PK model for IV cefepime in critically ill ICU patients (Jonckheere 2019), updated by simultaneously fitting plasma + urine PK from the original Jonckheere 2017 pilot (STDY1) and the Jonckheere 2019 target-controlled-infusion cohort (STDY2). Total clearance is the sum of an estimated-creatinine-clearance-driven renal arm (CL_renal = 2.29 * (eCrCL/60)^0.943 L/h per 70 kg) and a covariate-free non-renal arm (CL_nonren = 0.795 L/h per 70 kg); all PK parameters are scaled allometrically with body weight (reference 70 kg, exponent 3/4 for clearances, 1 for volumes). The structural form encodes the non-dialysis patient (paper Equations 1-4); a separate CL_dialysis = 4.48 L/h applied during intermittent hemodialysis sessions in the source dataset is documented in the vignette but not enabled in this model file."
  reference <- paste(
    "Jonckheere S, De Neve N, Verbeke J, De Decker K, Brandt I, Boel A,",
    "Van Bocxlaer J, Struys MMRF, Colin PJ. (2020).",
    "Target-Controlled Infusion of Cefepime in Critically Ill Patients.",
    "Antimicrob Agents Chemother 64(1):e01552-19.",
    "doi:10.1128/AAC.01552-19"
  )
  vignette <- "Jonckheere_2019_cefepime"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on CL, Q (exponent 3/4) and Vc, Vp (exponent 1) with reference weight 70 kg per Jonckheere 2019 Equations 1-4 and Discussion paragraph 4.",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Estimated creatinine clearance by the Cockcroft-Gault formula (raw, NOT BSA-normalized).",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying covariate; carried as the raw Cockcroft-Gault value in mL/min, NOT the canonical BSA-normalized mL/min/1.73 m^2 form (the same deviation documented in Frey_2010_tocilizumab.R for the same paper-reported quantity). Reference value 60 mL/min per Jonckheere 2019 Table 2 / Equation 1, approximating the cohort median (50.4 mL/min, IQR 29.1-100). The paper interpolates eCrCL between observation times by constant backward-fill (Methods, Update of the previously reported population pharmacokinetic model paragraph). Methods state 'eCrCL was calculated according to the Cockcroft-Gault equation'; the CKD-EPI value (Table 1) is reported for context only and is NOT the covariate driving CL_renal.",
      source_name        = "CRCL"
    )
  )

  population <- list(
    n_subjects     = 21L,
    n_studies      = 2L,
    age_range      = "IQR 72-78 years (median 76)",
    weight_range   = "IQR 67-86 kg (median 76)",
    sex_female_pct = 24.0,
    bmi_range      = "IQR 23.5-27.8 kg/m^2 (median 26.3)",
    bsa_range      = "IQR 1.77-2.03 m^2 (median 1.88)",
    sofa_range     = "IQR 3-8 (median 7) at ICU inclusion",
    crcl_cg_range  = "IQR 29.1-100 mL/min (median 50.4) Cockcroft-Gault",
    crcl_mdrd_range = "IQR 30.1-106 mL/min/1.73 m^2 (median 42.3) MDRD",
    crcl_ckdepi_range = "IQR 26.8-83.3 mL/min/1.73 m^2 (median 38.8) CKD-EPI",
    serum_creat_range = "IQR 0.66-2.14 mg/dL (median 1.49)",
    crp_range      = "IQR 95.7-287 mg/L at study inclusion (median 197); declined over 96 h",
    icu_los_range  = "IQR 5-10 days (median 7)",
    hospital_mortality_pct = 24.0,
    n_observations = 201L,
    obs_per_subject = "median 10 (IQR 9-11) plasma samples per patient",
    disease_state  = "Critically ill ICU patients (n=21) on continuous-infusion cefepime via target-controlled infusion targeting 16 mg/L; indications were respiratory infection (86%), abdominal infection (5%), combined respiratory + abdominal (5%), and unknown origin (5%). Pathogens isolated in 76% of patients (Klebsiella spp., E. coli, Citrobacter spp., Proteus, Pseudomonas, Morganella, Enterobacter, S. aureus, H. influenzae).",
    dose_range     = "Median daily cefepime dose 1.3-1.8 g/day delivered as continuous IV infusion via TCI; treatment median 4.0 days (IQR 2.0-5.0). Maximum infusion rate capped at 4 g/h.",
    regions        = "Single-centre Belgian cohort: OLV Hospital, Aalst (intensive care unit). Patients enrolled May 2016 - August 2017. ClinicalTrials.gov NCT02688582.",
    notes          = "Baseline demographics from Jonckheere 2019 Table 1. The model file encodes the non-dialysis population PK; intermittent hemodialysis (IHD) was a special case in the source data with CL_dialysis = 4.48 L/h applied during dialysis sessions and CL_renal = 0 between sessions (renal clearance was assumed absent in IHD patients). To simulate a dialysis patient, set CRCL = 0 and add the CL_dialysis term manually -- see vignette Errata. The model was developed by simultaneously fitting STDY1 (the Jonckheere 2017 pilot, doi:10.1128/AAC.00756-16) plus STDY2 (this study); plasma residual error from STDY2 (12.8% CV) is used here, see vignette Assumptions and deviations."
  )

  ini({
    # Structural parameters -- typical values for a non-dialysis patient at
    # 70 kg body weight and eCrCL = 60 mL/min (Cockcroft-Gault). All seven
    # estimates come from Jonckheere 2019 Table 2.
    lcl_renal <- log(2.29);   label("Renal clearance for a 70 kg patient at eCrCL = 60 mL/min (CL_renal, L/h)")  # Jonckheere 2019 Table 2 (theta_1)
    e_crcl_cl_renal <- 0.943; label("Power exponent of (eCrCL / 60) on CL_renal (unitless)")                     # Jonckheere 2019 Table 2 (theta_2)
    lcl_nonren <- log(0.795); label("Non-renal clearance for a 70 kg patient (CL_other in the paper, L/h)")       # Jonckheere 2019 Table 2 (CL_other)
    lvc <- log(10.7);         label("Central volume of distribution for a 70 kg patient (V_1, L)")                # Jonckheere 2019 Table 2 (V_1)
    lvp <- log(12.2);         label("Peripheral volume of distribution for a 70 kg patient (V_2, L)")             # Jonckheere 2019 Table 2 (V_2)
    lq  <- log(11.0);         label("Inter-compartmental clearance for a 70 kg patient (Q_2, L/h)")               # Jonckheere 2019 Table 2 (Q_2)

    # Allometric exponents (fixed by Jonckheere 2019 Discussion paragraph 4:
    # "scaling of all PK parameters with body weight according to allometric
    # theory"; 3/4 for clearance terms, 1 for volume terms).
    e_wt_cl_q  <- fixed(0.75); label("Allometric (WT) exponent on total CL and Q (unitless)")  # Jonckheere 2019 Discussion (allometric theory)
    e_wt_vc_vp <- fixed(1.00); label("Allometric (WT) exponent on V_1 and V_2 (unitless)")     # Jonckheere 2019 Discussion (allometric theory)

    # Inter-individual variability. Jonckheere 2019 Table 2 (footnote b)
    # reports BSV as CV%; the table footnote defines CV% via
    # sqrt(omega^2) * 100 where omega^2 is the NONMEM variance estimate, so
    # the back-transform omega^2 = (CV)^2 (NOT log(1 + CV^2)). The paper
    # treats omega as a standard deviation on the natural-log-transformed
    # individual parameters, so a reported "CV%" is omega itself expressed
    # as a percentage, not a CV-on-the-arithmetic-scale.
    #   CL_renal  : 24.6% CV -> omega^2 = 0.246^2  = 0.060516
    #   V_1       : 45.7% CV -> omega^2 = 0.457^2  = 0.208849
    #   CL_nonren : 69.4% CV -> omega^2 = 0.694^2  = 0.481636
    # No correlations between etas are reported, so each is independent.
    etalcl_renal  ~ 0.060516   # Jonckheere 2019 Table 2 (BSV CL_renal, 24.6% CV)
    etalvc        ~ 0.208849   # Jonckheere 2019 Table 2 (BSV V_1, 45.7% CV)
    etalcl_nonren ~ 0.481636   # Jonckheere 2019 Table 2 (BSV CL_other, 69.4% CV)

    # Residual error -- proportional model on plasma. Jonckheere 2019
    # Table 2 reports two plasma residual-error magnitudes from the joint
    # fit: STDY1 (the 2017 pilot, 31.8% CV) and STDY2 (the current TCI
    # cohort, 12.8% CV). The TCI use case the paper develops corresponds
    # to STDY2, so STDY2 is taken as the default plasma residual error
    # here; the STDY1 alternative is documented in vignette deviations.
    propSd <- 0.128; label("Proportional residual error on plasma (fraction)")  # Jonckheere 2019 Table 2 (Plasma_STDY2, 12.8% CV)
  })
  model({
    # Individual clearance components. Renal clearance scales with eCrCL
    # via a power covariate (Jonckheere 2019 Equation 1); the non-renal
    # arm is covariate-free at the typical level.
    cl_renal  <- exp(lcl_renal  + etalcl_renal)  * (CRCL / 60)^e_crcl_cl_renal
    cl_nonren <- exp(lcl_nonren + etalcl_nonren)
    cl <- (cl_renal + cl_nonren) * (WT / 70)^e_wt_cl_q

    # Volumes and inter-compartmental clearance: allometric on body weight.
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    vp <- exp(lvp)          * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 70)^e_wt_cl_q

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Plasma concentration; dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
