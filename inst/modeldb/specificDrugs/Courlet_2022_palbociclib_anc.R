Courlet_2022_palbociclib_anc <- function() {
  description <- "Final population PK/PD model of palbociclib-induced neutropenia in women with advanced breast cancer followed in routine care (Courlet 2022, Table 2; PK and PD parameters estimated jointly). PK: two-compartment model with first-order absorption, an absorption lag time and first-order elimination; apparent clearance takes a separately estimated value when palbociclib is taken under fasting conditions (or with a light meal) together with a proton-pump inhibitor. PD: Friberg semi-mechanistic myelosuppression model (proliferating pool, three maturation transit compartments, circulating neutrophils, (Base / ANC)^gamma feedback) with an Emax palbociclib effect inhibiting proliferation (EC50 fixed to the literature value). IIV on ka, CL/F, Vc/F, Base, MTT, Emax and EC50; proportional residual error on palbociclib and additive error on log-transformed ANC. The PK-only model of the same paper is Courlet_2022_palbociclib."
  reference <- paste(
    "Courlet P, Cardoso E, Bandiera C, Stravodimou A, Zurcher JP, Chtioui H, Locatelli I, Decosterd LA,",
    "Darnaud L, Blanchet B, Alexandre J, Wagner AD, Zaman K, Schneider MP, Guidi M, Csajka C. (2022).",
    "Population Pharmacokinetics of Palbociclib and Its Correlation with Clinical Efficacy and Safety in",
    "Patients with Advanced Breast Cancer. Pharmaceutics 14(7):1317. doi:10.3390/pharmaceutics14071317.",
    sep = " "
  )
  vignette <- "Courlet_2022_palbociclib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL", anc = "10^9/L")
  # Unit notes.
  # 1. PK. Table 2 reports CL/F and Q in L/h, Vc/F and Vp/F in L, ka in 1/h
  #    and ALAG in h; doses are in mg, so central / vc is in mg/L and the 1000
  #    factor gives Cc in ng/mL (the assay unit, Section 2.2).
  # 2. PD. EC50 is in ng/mL (Table 2), matching Cc. MTT is in hours (Table 2)
  #    and ANC (Base) in G/L = 10^9 cells/L.

  compartmentData <- list(
    depot = list(analyte = "palbociclib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "palbociclib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "palbociclib", units = "mg", specimen = "tissue", verified = TRUE),
    prol = list(analyte = "neutrophil progenitor cells", units = "10^9/L", specimen = "tissue", verified = TRUE),
    transit1 = list(analyte = "maturing neutrophils", units = "10^9/L", specimen = "tissue", verified = TRUE),
    transit2 = list(analyte = "maturing neutrophils", units = "10^9/L", specimen = "tissue", verified = TRUE),
    transit3 = list(analyte = "maturing neutrophils", units = "10^9/L", specimen = "tissue", verified = TRUE),
    circ = list(analyte = "neutrophils", units = "10^9/L", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    CONMED_PPI = list(
      description = "Concomitant proton-pump inhibitor co-administered with the palbociclib dose",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no PPI co-administration)",
      notes = "Per-concentration (time-varying) indicator; 78 of 255 concentrations (31%) were measured under PPI co-administration (Table 1). PPI use alone had no effect on any PK parameter (Section 3.2.1); it acts only in combination with fasting, through the interaction CONMED_PPI * (1 - FED) that selects the separately estimated clearance lcl_ppifasted.",
      source_name = "PPI"
    ),
    FED = list(
      description = "Palbociclib taken with a meal (1) versus under fasting conditions or with a light meal (0)",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (taken with a meal, as labelled)",
      notes = "Per-concentration indicator. The source flags 'administration under fasting conditions' (54 of 255 concentrations, 21%, Table 1), and Section 3.2.1 groups 'fasting conditions (or with a light meal)' together, so FED = 1 - fasting flag. Fasting alone had no retained effect; it acts only together with a PPI (CONMED_PPI * (1 - FED)).",
      source_name = "fasting conditions"
    )
  )

  covariatesDataExcluded <- list(
    GGT = list(
      description = "Gamma-glutamyltransferase at baseline",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Significant on Base in the univariate step (Section 3.2.2, dOFV < -6.3) but removed at backward deletion. Table 1 median 29 [18-51] U/L.",
      source_name = "GGT"
    ),
    CRCL = list(
      description = "BSA-normalized renal function (Cockcroft-Gault estimate per the source)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "The paper's eGFR (Cockcroft-Gault) was significant on Emax in the univariate step (Section 3.2.2) but removed at backward deletion; also tested on CL/F. Table 1 reports it as 71 [54-98] in mL/min/1.73 m^2.",
      source_name = "eGFR"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on the PK and PD parameters (Section 2.3.3); not retained. Table 1 median 65 [55-75] years.",
      source_name = "Age"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL/F, Vc/F and the PD parameters (Section 2.3.3); not retained. Table 1 median 67 [61-80] kg.",
      source_name = "Body weight"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL/F and the PD parameters (Section 2.3.3); not retained. Table 1 median 43 [41-45] g/L.",
      source_name = "Albumin"
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL/F and the PD parameters (Section 2.3.3); not retained. Table 1 median 61 [49-81] U/L.",
      source_name = "ALK"
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL/F and the PD parameters (Section 2.3.3); not retained. Table 1 median 5 [4-7] umol/L.",
      source_name = "BILT"
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL/F and the PD parameters (Section 2.3.3); not retained.",
      source_name = "AST"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL/F and the PD parameters (Section 2.3.3); not retained.",
      source_name = "ALT"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 44L,
    n_studies = 1L,
    age_range = "median 65 years (IQR 55-75)",
    weight_range = "median 67 kg (IQR 61-80)",
    sex_female_pct = 100,
    disease_state = "Women treated with palbociclib for advanced (metastatic) breast cancer in routine care (OpTAT study, NCT04484064); 70% of ANC records under concomitant fulvestrant, 22% after chemotherapy in the previous line.",
    dose_range = "Oral palbociclib 75-125 mg once daily (median 100 mg), 21 days on / 7 days off (28-day cycle).",
    regions = "Switzerland (Lausanne University Hospital)",
    observations = "255 palbociclib plasma concentrations (2.0-159.0 ng/mL) and 1174 absolute neutrophil counts (40 before palbociclib initiation); ANC within 7 days of G-CSF administration excluded.",
    notes = "Baseline characteristics from Table 1. Estimated in NONMEM 7.4.3 with FOCEI (ADVAN6); PK and PD parameters were estimated together in the final model (Section 2.3.2)."
  )

  ini({
    # ---- PK: Table 2, 'Pharmacokinetics' rows of the final PK/PD model ----
    lka <- log(0.8); label("ka: first-order absorption rate constant (1/h)") # Table 2: ka = 0.8 1/h (RSE 30%)
    ltlag <- log(2.0); label("ALAG: absorption lag time (h)") # Table 2: ALAG = 2.0 h (RSE 8%)
    lcl <- log(67); label("CL/F: apparent clearance, fed or without PPI (L/h)") # Table 2: CL/F = 67 L/h (RSE 5%)
    lcl_ppifasted <- log(131); label("CL/F under fasting conditions with a co-administered PPI (L/h)") # Table 2: CL/F PPI,no food = 131 L/h (RSE 4%)
    lvc <- log(2800); label("Vc/F: apparent central volume of distribution (L)") # Table 2: Vc/F = 2800 L (RSE 7%)
    lq <- log(7); label("Q/F: apparent intercompartmental clearance (L/h)") # Table 2: Q = 7 L/h (RSE 31%)
    lvp <- log(704); label("Vp/F: apparent peripheral volume of distribution (L)") # Table 2: Vp/F = 704 L (RSE 9%)

    # ---- PD: Table 2, 'Pharmacodynamics' rows ----
    lcirc0 <- log(4.1); label("Base: baseline circulating neutrophil count (10^9/L)") # Table 2: Base = 4.1 G/L (RSE 7%)
    lmtt <- log(122); label("MTT: mean transit time through the maturation chain (h)") # Table 2: MTT = 122 h (RSE 5%)
    lemax <- log(0.22); label("Emax: maximum fractional inhibition of proliferation by palbociclib (unitless)") # Table 2: Emax = 0.22 (RSE 7%)
    # Section 3.2.2: the estimate (8.8 ng/mL) was unreliable, so EC50 was fixed
    # to the literature value of reference 19.
    lec50 <- fixed(log(40.1)); label("EC50: palbociclib concentration giving half of Emax (ng/mL)") # Table 2: EC50 = 40.1 ng/mL, FIX
    lgamma <- log(0.13); label("gamma: feedback exponent on (Base / ANC) (unitless)") # Table 2: gamma = 0.13 (RSE 9%)

    # ---- IIV: exponential (Section 2.3.1); CV% converted with omega^2 = log(CV^2 + 1) ----
    etalka ~ 0.94098 # Table 2: omega ka = 125% CV -> log(1.25^2 + 1) = 0.94098
    etalcl ~ 0.080750 # Table 2: omega CL = 29% CV -> log(0.29^2 + 1) = 0.080750
    etalvc ~ 0.097490 # Table 2: omega Vc = 32% CV -> log(0.32^2 + 1) = 0.097490
    etalcirc0 ~ 0.115567 # Table 2: omega base = 35% CV -> log(0.35^2 + 1) = 0.115567
    etalmtt ~ 0.014297 # Table 2: omega MTT = 12% CV -> log(0.12^2 + 1) = 0.014297
    etalemax ~ 0.022250 # Table 2: omega Emax = 15% CV -> log(0.15^2 + 1) = 0.022250
    etalec50 ~ 0.623207 # Table 2: omega EC50 = 93% CV -> log(0.93^2 + 1) = 0.623207

    # ---- Residual error ----
    propSd <- 0.18; label("Proportional residual error on palbociclib Cc (fraction)") # Table 2: proportional residual error = 18% (RSE 9%)
    # Section 2.3.2: ANC log-transformed, additive error on the log scale.
    expSd_ANC <- 0.31; label("Log-scale (exponential) residual error on ANC (SD, log units)") # Table 2: additive residual error (log scale) = 0.31 (RSE 9%)
  })

  model({
    # ---- 1. PK individual parameters ----
    # Section 3.2.1: CL/F differs only when palbociclib is taken fasting (or
    # with a light meal) AND with a PPI; PPI use alone has no effect.
    ppi_fasted <- CONMED_PPI * (1 - FED)
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag)
    cl <- exp((1 - ppi_fasted) * lcl + ppi_fasted * lcl_ppifasted + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- 2. PK ODEs: two compartments, first-order absorption with lag ----
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    alag(depot) <- tlag
    # central / vc is mg/L; 1 mg/L = 1000 ng/mL.
    Cc <- 1000 * central / vc

    # ---- 3. PD individual parameters ----
    circ0 <- exp(lcirc0 + etalcirc0)
    mtt <- exp(lmtt + etalmtt)
    emax <- exp(lemax + etalemax)
    ec50 <- exp(lec50 + etalec50)
    gamma <- exp(lgamma)
    # Section 2.3.2 and Figure 1: ktr = (n + 1) / MTT with n = 3 transit
    # compartments, and kprol = ktr = kcirc.
    ktr <- 4 / mtt

    # ---- 4. Friberg myelosuppression chain (Figure 1) ----
    # Equation 2: E_palbo = Emax * C / (EC50 + C), acting on the proliferation
    # compartment; feedback (Base / Circ_t)^gamma.
    edrug <- emax * Cc / (ec50 + Cc)
    feed <- (circ0 / circ)^gamma
    d/dt(prol) <- ktr * prol * (1 - edrug) * feed - ktr * prol
    d/dt(transit1) <- ktr * prol - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(circ) <- ktr * transit3 - ktr * circ

    # Friberg model: before any palbociclib, all five compartments equal Base.
    prol(0) <- circ0
    transit1(0) <- circ0
    transit2(0) <- circ0
    transit3(0) <- circ0
    circ(0) <- circ0

    # ---- 5. Observations and error ----
    ANC <- circ
    Cc ~ prop(propSd)
    ANC ~ lnorm(expSd_ANC)
  })
}
