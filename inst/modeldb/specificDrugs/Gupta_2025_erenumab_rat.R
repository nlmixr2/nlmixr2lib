Gupta_2025_erenumab_rat <- function() {
  description <- paste(
    "Preclinical (rat). QSP. Integrated PK / semi-mechanistic immune-cell",
    "model of anti-drug antibody (ADA) formation and immune tolerance",
    "induction against the human IgG2 monoclonal antibody erenumab in",
    "Sprague-Dawley rats (Gupta 2025): one-compartment first-order",
    "absorption erenumab PK with aggregation-state effects on CL/F and ka and",
    "an ADA-magnitude power effect on CL/F, driving dendritic-cell",
    "maturation and CD4+ T-helper / T-regulatory cell activation and",
    "proliferation; tacrolimus (inhibits T-helper proliferation) and",
    "methotrexate (inhibits T-helper activation) carried as coupled",
    "one-compartment PK states."
  )
  reference <- paste(
    "Gupta P, Ryman JT, Jawa V, Meibohm B.",
    "A Semi-Mechanistic Mathematical Model of Immune Tolerance Induction",
    "to Support Preclinical Studies of Human Monoclonal Antibodies in Rats.",
    "Pharmaceutics. 2025;17(7):845. doi:10.3390/pharmaceutics17070845.",
    "PK parameters from Table 1; immune-cell parameters from Table 2 and",
    "Supplementary Table S1; equations from Equations 1-12 and",
    "Supplementary Equations S1-S15.",
    sep = " "
  )
  vignette <- "Gupta_2025_erenumab_rat"
  units <- list(
    time = "day",
    dosing = "ug (erenumab, absolute amount); mg/kg (tacrolimus and methotrexate doses into depot_tacrolimus / depot_methotrexate)",
    concentration = "ug/mL"
  )
  dosing <- c("depot", "depot_methotrexate", "depot_tacrolimus")

  # Immune-cell states use the paper's own symbols (Figure 2, Equations
  # S4-S9); the tacrolimus and methotrexate PK states are the immunosuppressant
  # concentration drivers of Supplementary Equations S12-S15.
  paper_specific_compartments <- c(
    "ND",
    "MD",
    "NThlp",
    "AThlp",
    "NTreg",
    "ATreg",
    "depot_tacrolimus",
    "central_tacrolimus",
    "depot_methotrexate",
    "central_methotrexate"
  )
  # Study-specific variance terms: Gupta 2025 estimated IIV on Vd/F and the
  # combined residual error separately for Study 1 and Study 2 (Section 2.2.1,
  # Table 1).
  paper_specific_etas <- c("etalvc_s2")
  paper_specific_residual_sds <- c("addSd_s2", "propSd_s2")

  covariateData <- list(
    FORM_ERENUMAB_AGGREGATE = list(
      description = "Erenumab aggregation state of the administered dose: 1 = aggregated erenumab, 0 = monomeric erenumab",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Multiplies CL/F by theta_aggregation,CL = 2.40 and ka by",
        "theta_aggregation,Ka = 0.12 (Equation 3, Table 1), and switches the",
        "dendritic-cell activation EC50 from the monomer value 9.85 ug/mL to",
        "the aggregate value 2.02 ug/mL (Table 2). In the source data all",
        "Study 2 animals received aggregates and all Study 1 animals",
        "received monomer, so this column and STUDY_ERENUMAB_S2 are",
        "identical in the fitted data set."
      ),
      source_name = "aggregation state (monomer = 0, aggregate = 1)"
    ),
    STUDY_ERENUMAB_S2 = list(
      description = "Gupta 2025 / Ryman immune-tolerance Study 2 indicator: 1 = Study 2 (1 mg/kg weekly erenumab aggregates), 0 = Study 1 (10 mg/kg weekly erenumab monomer)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Selects the study-specific IIV on Vd/F (86.9% Study 1, 21.2% Study",
        "2), the study-specific combined residual error (Table 1), and the",
        "ADA-positivity S/N cut-off used in the ADA effect on CL/F (S/N > 2",
        "Study 1, S/N > 1.5 Study 2; Section 2.1). The paper attributes the",
        "separate variance terms to different assays, instruments, erenumab",
        "batches and animal populations in studies run 15 months apart."
      ),
      source_name = "Study 1 / Study 2"
    ),
    OCC = list(
      description = "Dosing occasion (1-8) for the inter-occasion variability on ka",
      units = "(integer)",
      type = "count",
      reference_category = NULL,
      notes = paste(
        "Each of the eight weekly erenumab injections (weeks 1-4 induction,",
        "weeks 9-12 rechallenge) is its own occasion (Section 2.2.1), so OCC",
        "= k from the k-th erenumab dose until the next. OCC outside 1-8",
        "applies no IOV term."
      ),
      source_name = "occasion"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "erenumab", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "erenumab", units = "ug", specimen = "serum", verified = TRUE),
    depot_tacrolimus = list(
      analyte = "tacrolimus",
      units = "mg/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    central_tacrolimus = list(analyte = "tacrolimus", units = "mg/kg", specimen = "plasma", verified = TRUE),
    depot_methotrexate = list(
      analyte = "methotrexate",
      units = "mg/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    central_methotrexate = list(analyte = "methotrexate", units = "mg/kg", specimen = "plasma", verified = TRUE),
    ND = list(analyte = "naive dendritic cells", units = "cells/uL", specimen = "blood cell", verified = TRUE),
    MD = list(analyte = "mature dendritic cells", units = "cells/uL", specimen = "blood cell", verified = TRUE),
    NThlp = list(analyte = "naive CD4+ T-helper cells", units = "cells/uL", specimen = "blood cell", verified = TRUE),
    AThlp = list(
      analyte = "activated CD4+ T-helper cells",
      units = "cells/uL",
      specimen = "blood cell",
      verified = TRUE
    ),
    NTreg = list(analyte = "naive T-regulatory cells", units = "cells/uL", specimen = "blood cell", verified = TRUE),
    ATreg = list(analyte = "activated T-regulatory cells", units = "cells/uL", specimen = "blood cell", verified = TRUE)
  )

  population <- list(
    species = "rat (Sprague-Dawley)",
    n_subjects = 52,
    n_studies = 2,
    weight_range = "approximately 270-290 g at study start, 380-410 g at week 12 (Discussion)",
    disease_state = "Healthy animals in two 12-week erenumab immune-tolerance studies (4-week induction, 4-week washout, 4-week rechallenge).",
    dose_range = paste(
      "Erenumab SC once weekly in weeks 1-4 and 9-12: 10 mg/kg monomer",
      "(Study 1) or 1 mg/kg aggregates (Study 2). Immunosuppression in the",
      "induction phase: methotrexate 5 mg/kg IP weekly weeks 1-7 (Group",
      "1.2), methotrexate 3 mg/kg IP on days 1-3 (Group 1.3), or sirolimus",
      "6 mg/kg PO once on day 1 plus tacrolimus 2 mg/kg SC weekly weeks 1-7",
      "(Groups 1.4 and 2.2)."
    ),
    regions = "United States (University of Tennessee Health Science Center)",
    notes = paste(
      "Study 1: four groups of n = 8 (control, two methotrexate regimens,",
      "tacrolimus/sirolimus); Study 2: two groups of n = 10 (control,",
      "tacrolimus/sirolimus). 1924 erenumab serum concentrations and 624 ADA",
      "S/N measurements; 14 animals ADA-positive. Sex not reported in the",
      "source. Section 2.1 and Section 3."
    )
  )

  ini({
    # --- Erenumab population PK (Gupta 2025 Table 1, Monolix SAEM) ---
    lka <- log(0.91); label("Absorption rate constant for monomeric erenumab (1/day)") # Table 1 'Ka, day-1 0.91 (7.94)'
    lcl <- log(1.83); label("Apparent clearance CL/F for monomeric erenumab, ADA-negative (mL/day)") # Table 1 'Base CL/F, mL/day 1.83 (3.82)'
    lvc <- log(32.6); label("Apparent volume of distribution Vd/F (mL)") # Table 1 'Base Vd/F, mL 32.6 (3.74)'
    e_aggregate_cl <- 2.40; label("Multiplicative effect of erenumab aggregates on CL/F (unitless)") # Table 1 'theta aggregation, CL 2.40 (5.60)'
    e_aggregate_ka <- 0.12; label("Multiplicative effect of erenumab aggregates on ka (unitless)") # Table 1 'theta aggregation, Ka 0.12 (13.6)'
    e_ada_cl <- 0.37; label("Exponent of the ADA S/N effect on CL/F (unitless)") # Table 1 'theta ADA Exponent in ADA effect on CL/F 0.37 (3.27)'
    sn_cutoff <- fixed(2); label("ADA-positivity S/N cut-off, Study 1 (unitless)") # Section 2.1 'S/N ratio of greater than 2 ... for Study 1'; Equation 4 SNcutoff
    sn_cutoff_s2 <- fixed(1.5); label("ADA-positivity S/N cut-off, Study 2 (unitless)") # Section 2.1 '... and 1.5 ... for Study 2'; Equation 4 SNcutoff

    # --- Immune-cell model: literature values (Table 2 part C, Table S1), fixed ---
    betaND <- fixed(0.0924); label("Death rate of naive dendritic cells (1/day)") # Table 2 C 'betaND 0.0924'
    ND0 <- fixed(3700); label("Initial number of naive dendritic cells (cells/uL)") # Table 2 C 'ND0 3700'
    deltaND <- fixed(1.5); label("Maximum activation rate of naive dendritic cells (1/day)") # Table 2 C 'deltaND 1.5'
    ec50 <- fixed(9.85); label("Erenumab monomer concentration for half-maximal naive DC activation (ug/mL)") # Table 2 C 'EC50monomer 9.85'
    betaMD <- fixed(0.2310); label("Death rate of mature dendritic cells (1/day)") # Table 2 C 'betaMD 0.2310'
    NThlp0 <- fixed(723); label("Initial number of naive CD4+ T-helper cells (cells/uL)") # Table 2 C 'NThlp,0 723'
    betaNThlp <- fixed(0.0056); label("Death rate of naive CD4+ T-helper cells (1/day)") # Table 2 C 'betaNThlp 0.0056'
    deltaNThlp <- fixed(1.5); label("Maximum activation rate of naive CD4+ T-helper cells (1/day)") # Table 2 C 'deltaNThlp 1.5'
    NTreg0 <- fixed(62); label("Initial number of naive T-regulatory cells (cells/uL)") # Table 2 C 'NTreg,0 62'
    betaATreg <- fixed(0.18); label("Death rate of activated T-regulatory cells (1/day)") # Table 2 C 'betaATreg 0.18'
    betaNTreg <- fixed(0.0056); label("Death rate of naive T-regulatory cells (1/day)") # Table 2 C 'betaNTreg 0.0056'
    deltaNTreg <- fixed(1.5); label("Maximum activation rate of naive T-regulatory cells (1/day)") # Table 2 C 'deltaNTreg 1.5'
    alphaADA <- fixed(77); label("ADA secretion rate (nmol/day)") # Table 2 C 'alpha 77'; Table S1 unit nmol/day
    kADA <- fixed(0.138); label("ADA elimination rate constant (1/day)") # Table 2 C 'KADA 0.138'

    # --- Immune-cell model: calibrated values (Table 2 part A, Table S1), fixed ---
    ec50_aggregate <- fixed(2.02); label("Erenumab aggregate concentration for half-maximal naive DC activation (ug/mL)") # Table 2 A 'EC50aggregate 2.02'
    rhoAThlp <- fixed(2.45); label("Maximum proliferation rate of activated CD4+ T-helper cells (1/day)") # Table 2 A 'rhoAThlp 2.45'
    betaAThlp <- fixed(0.05); label("Death rate of activated CD4+ T-helper cells (1/day)") # Table 2 A 'betaAThlp 0.05'
    kADAprop <- fixed(500); label("ADA proportionality constant K (cells/day)") # Table S1 'K Constant Cells/day 500 Model calibration'
    nADA <- fixed(1000); label("ADA-to-S/N scaling constant N (nmol/uL)") # Table S1 'N Constant nmol/uL 1000 Model calibration'

    # --- Immune-cell model: estimated values (Table 2 part B) ---
    lrhoATreg <- log(3.02); label("Maximum proliferation rate of activated T-regulatory cells (1/day)") # Table 2 B 'rhoATreg 3.02 (3.01)'
    lATreg50 <- log(1000); label("Activated T-regulatory cells for half-maximal suppression of T-helper proliferation (cells/uL)") # Table 2 B 'ATreg50 1000' (no RSE printed)
    ltlag_ada <- log(28.3); label("Lag time for ADA formation (day)") # Table 2 B 'Tlag 28.3 (2.19)'

    # --- Tacrolimus PK (Table S1, Equation S12), fixed literature values ---
    f_tacrolimus <- fixed(0.2); label("Tacrolimus SC bioavailability (fraction)") # Table S1 'FTAC 0.2'
    ka_tacrolimus <- fixed(7.2); label("Tacrolimus absorption rate constant (1/day)") # Table S1 'KaTAC 7.2'
    kel_tacrolimus <- fixed(3.6); label("Tacrolimus elimination rate constant (1/day)") # Table S1 'KeTAC 3.6'
    vc_tacrolimus <- fixed(10); label("Tacrolimus volume of distribution (L/kg)") # Table S1 'VdTAC 10'
    ic50_tacrolimus <- fixed(0.2); label("Tacrolimus IC50 for inhibition of T-helper proliferation (nM)") # Table S1 'IC50TAC 0.2'; Section 2.2.2 'IC50 of 0.2 nM'
    mw_tacrolimus <- fixed(804.02); label("Tacrolimus molecular weight (g/mol)") # C44H69NO12; not printed in the source, needed to compare CTAC (mg/L) with IC50TAC (nM)

    # --- Methotrexate PK (Table S1, Equation S14), fixed literature values ---
    f_methotrexate <- fixed(1); label("Methotrexate IP bioavailability (fraction)") # Table S1 'FMTX 1'
    ka_methotrexate <- fixed(88.56); label("Methotrexate absorption rate constant (1/day)") # Table S1 'KaMTX 88.56'
    kel_methotrexate <- fixed(56.41); label("Methotrexate elimination rate constant (1/day)") # Table S1 'KeMTX 56.41'
    vc_methotrexate <- fixed(0.174); label("Methotrexate volume of distribution (L/kg)") # Table S1 'VdMTX 0.174'
    ic50_methotrexate <- fixed(1); label("Methotrexate IC50 for inhibition of T-helper activation (nM)") # Table S1 'IC50MTX 1'; Section 2.2.2 'IC50MTX of 1 nM'
    mw_methotrexate <- fixed(454.44); label("Methotrexate molecular weight (g/mol)") # C20H22N8O5; not printed in the source, needed to compare CMTX (mg/L) with IC50MTX (nM)

    # --- PK random effects (Table 1). Monolix reports omega (SD of eta); variance = omega^2 ---
    etalcl ~ 0.031684 # Table 1 IIV CL/F 17.8% -> 0.178^2
    etalvc ~ 0.755161 # Table 1 IIV Vd/F Study 1 86.9% -> 0.869^2
    etalvc_s2 ~ 0.044944 # Table 1 IIV Vd/F Study 2 21.2% -> 0.212^2
    etaiov_lka_1 ~ 3.3489 # Table 1 IOV Ka 183.0% -> 1.83^2; one occasion per erenumab dose
    etaiov_lka_2 ~ fixed(3.3489) # shared IOV variance, occasion 2
    etaiov_lka_3 ~ fixed(3.3489) # shared IOV variance, occasion 3
    etaiov_lka_4 ~ fixed(3.3489) # shared IOV variance, occasion 4
    etaiov_lka_5 ~ fixed(3.3489) # shared IOV variance, occasion 5
    etaiov_lka_6 ~ fixed(3.3489) # shared IOV variance, occasion 6
    etaiov_lka_7 ~ fixed(3.3489) # shared IOV variance, occasion 7
    etaiov_lka_8 ~ fixed(3.3489) # shared IOV variance, occasion 8

    # --- Immune-model random effects (Table 2 part B); variance = omega^2 ---
    etalrhoATreg ~ 0.023409 # Table 2 B IIV rhoATreg 15.3% -> 0.153^2
    etalATreg50 ~ 0.001156 # Table 2 B IIV ATreg50 3.4% -> 0.034^2
    etaltlag_ada ~ 0.2809 # Table 2 B IIV Tlag 53% -> 0.53^2

    # --- Residual error on erenumab concentrations (Table 1), combined additive + proportional ---
    addSd <- 0.09; label("Additive residual error, Study 1 (ug/mL)") # Table 1 'Additive residual error (Study 1) 0.09 (22.2)'
    propSd <- 0.29; label("Proportional residual error, Study 1 (fraction)") # Table 1 'Proportional residual error (Study 1) 0.29 (2.40)'
    addSd_s2 <- 0.12; label("Additive residual error, Study 2 (ug/mL)") # Table 1 'Additive residual error (Study 2) 0.12 (12.4)'
    propSd_s2 <- 0.13; label("Proportional residual error, Study 2 (fraction)") # Table 1 'Proportional residual error (Study 2) 0.13 (3.45)'
  })

  model({
    # 1. Inter-occasion variability on ka, one occasion per erenumab dose.
    iov_ka <- (OCC == 1) * etaiov_lka_1 + (OCC == 2) * etaiov_lka_2 +
      (OCC == 3) * etaiov_lka_3 + (OCC == 4) * etaiov_lka_4 +
      (OCC == 5) * etaiov_lka_5 + (OCC == 6) * etaiov_lka_6 +
      (OCC == 7) * etaiov_lka_7 + (OCC == 8) * etaiov_lka_8

    # 2. Individual erenumab PK parameters (Equations 1 and 3). Vd/F IIV is
    #    study-specific.
    ka <- exp(lka + iov_ka) * e_aggregate_ka^FORM_ERENUMAB_AGGREGATE
    cl_base <- exp(lcl + etalcl) * e_aggregate_cl^FORM_ERENUMAB_AGGREGATE
    vc <- exp(lvc + etalvc * (1 - STUDY_ERENUMAB_S2) + etalvc_s2 * STUDY_ERENUMAB_S2)

    # 3. Individual immune-model parameters.
    rhoATreg <- exp(lrhoATreg + etalrhoATreg)
    ATreg50 <- exp(lATreg50 + etalATreg50)
    tlag_ada <- exp(ltlag_ada + etaltlag_ada)
    ec50_ind <- ec50 * (1 - FORM_ERENUMAB_AGGREGATE) + ec50_aggregate * FORM_ERENUMAB_AGGREGATE

    # 4. ADA magnitude (Equations 10 and S10-S11): zero before Tlag, then
    #    proportional to the activated T-helper count. t is time since the
    #    first erenumab dose.
    ada <- alphaADA * AThlp / kADAprop * (1 - exp(-kADA * (time - tlag_ada))) * (time >= tlag_ada)
    ada_sn <- ada / nADA

    # 5. ADA effect on CL/F (Equation 4): power function of S/N relative to
    #    the study's ADA-positivity cut-off, applied only above the cut-off.
    #    The model-predicted S/N closes the PK / ADA loop of Figure 2.
    sn_cut <- sn_cutoff * (1 - STUDY_ERENUMAB_S2) + sn_cutoff_s2 * STUDY_ERENUMAB_S2
    ada_fac <- 1
    if (ada_sn > sn_cut) ada_fac <- (ada_sn / sn_cut)^e_ada_cl
    cl <- cl_base * ada_fac
    kel <- cl / vc

    # 6. Immunosuppressant concentrations (Equations S12 and S14 as
    #    single-dose superposition), converted from mg/L to nM.
    Cc_tacrolimus <- central_tacrolimus / vc_tacrolimus
    Cc_methotrexate <- central_methotrexate / vc_methotrexate
    c_tacrolimus_nM <- Cc_tacrolimus / mw_tacrolimus * 1e6
    c_methotrexate_nM <- Cc_methotrexate / mw_methotrexate * 1e6
    inh_tacrolimus <- ic50_tacrolimus / (c_tacrolimus_nM + ic50_tacrolimus)
    inh_methotrexate <- ic50_methotrexate / (c_methotrexate_nM + ic50_methotrexate)

    # 7. Erenumab concentration = maturation signal (Equations 5 and S3).
    Cc <- central / vc

    # Fraction of mature DCs available to co-operate with T cells (Equation 9).
    D <- MD / (MD + NThlp + AThlp + NTreg + ATreg)
    act_nd <- deltaND * Cc / (ec50_ind + Cc)

    # 8. ODE system.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    d/dt(depot_tacrolimus) <- -ka_tacrolimus * depot_tacrolimus
    d/dt(central_tacrolimus) <- ka_tacrolimus * depot_tacrolimus - kel_tacrolimus * central_tacrolimus
    d/dt(depot_methotrexate) <- -ka_methotrexate * depot_methotrexate
    d/dt(central_methotrexate) <- ka_methotrexate * depot_methotrexate - kel_methotrexate * central_methotrexate

    # Dendritic cells (Equations S4-S5).
    d/dt(ND) <- betaND * ND0 - betaND * ND - act_nd * ND
    d/dt(MD) <- act_nd * ND - betaMD * MD
    # Naive and activated CD4+ T-helper cells (Equations S6-S7). Methotrexate
    # inhibits the activation term (Equation S15) and tacrolimus the
    # proliferation term (Equation S13); each factor is 1 when its drug is
    # absent, so one equation covers every treatment group.
    d/dt(NThlp) <- betaNThlp * NThlp0 - betaNThlp * NThlp - deltaNThlp * D * NThlp
    d/dt(AThlp) <- deltaNThlp * D * NThlp * inh_methotrexate +
      rhoAThlp * AThlp * D * ATreg50 / (ATreg50 + ATreg) * inh_tacrolimus -
      betaAThlp * AThlp
    # Naive and activated T-regulatory cells (Equations S8-S9).
    d/dt(NTreg) <- betaNTreg * NTreg0 - betaNTreg * NTreg - deltaNTreg * D * NTreg
    d/dt(ATreg) <- deltaNTreg * D * NTreg + rhoATreg * ATreg * D - betaATreg * ATreg

    # 9. Bioavailability of the immunosuppressants (Equations S12, S14).
    f(depot_tacrolimus) <- f_tacrolimus
    f(depot_methotrexate) <- f_methotrexate

    # 10. Initial conditions: homeostatic naive pools; mature DCs and
    #     activated T cells start at zero.
    ND(0) <- ND0
    NThlp(0) <- NThlp0
    NTreg(0) <- NTreg0

    # 11. Observation and residual error (combined, study-specific).
    addSd_study <- addSd * (1 - STUDY_ERENUMAB_S2) + addSd_s2 * STUDY_ERENUMAB_S2
    propSd_study <- propSd * (1 - STUDY_ERENUMAB_S2) + propSd_s2 * STUDY_ERENUMAB_S2
    Cc ~ add(addSd_study) + prop(propSd_study)
  })
}
