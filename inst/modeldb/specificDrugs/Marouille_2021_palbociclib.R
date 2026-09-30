Marouille_2021_palbociclib <- function() {
  description <- "Population PK/PD model of palbociclib-induced neutropenia in women with breast cancer followed in routine care (Le Marouille 2021). PK: one-compartment model with first-order absorption, a fixed absorption lag time, a fixed apparent volume and first-order elimination; apparent oral clearance increases with Cockcroft-Gault creatinine clearance and decreases with serum alkaline phosphatase (power models centred on the cohort medians), IIV on CL/F only and additive residual error. PD: Friberg semi-mechanistic myelosuppression model (proliferating pool, three maturation transit compartments, circulating neutrophils, (Base / ANC)^gamma feedback) with a linear palbociclib effect Slope * C on proliferation; baseline ANC increases with age (power model). IIV on Base, Slope and MTT; exponential (log-normal) residual error on ANC."
  reference <- paste(
    "Le Marouille A, Petit E, Kaderbhai C, Desmoulins I, Hennequin A, Mayeur D, Fumet JD, Ladoire S,",
    "Tharin Z, Ayati S, Ilie S, Royer B, Schmitt A. (2021).",
    "Pharmacokinetic/Pharmacodynamic Model of Neutropenia in Real-Life Palbociclib-Treated Patients.",
    "Pharmaceutics 13(10):1708. doi:10.3390/pharmaceutics13101708.",
    "V/F and Tlag were fixed to the values of Royer B et al. (2021) Population Pharmacokinetics of",
    "Palbociclib in a Real-World Situation. Pharmaceuticals 14(3):181. doi:10.3390/ph14030181.",
    sep = " "
  )
  vignette <- "Marouille_2021_palbociclib"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L", anc = "10^9/L")
  # Unit notes.
  # 1. PK. Table 3 reports CL/F in L/h, V/F in L, ka in 1/h and Tlag in h;
  #    doses are in mg (Section 3.1: 75, 100 and 125 mg per day), so
  #    central / vc is in mg/L. Concentrations are reported in ug/L throughout
  #    (LLOQ 6 ng/mL, Section 2.2; additive error in ug/L, Table 3), hence the
  #    1000 ug/L per mg/L factor on Cc.
  # 2. PD. The drug-effect Slope is in L/ug (Table 4), i.e. it multiplies the
  #    concentration in ug/L. MTT is reported in days (Table 4) and is
  #    converted inline to hours so the whole model runs on one clock.
  #    ANC is in G/L (10^9 cells/L).

  compartmentData <- list(
    depot = list(analyte = "palbociclib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "palbociclib", units = "mg", specimen = "plasma", verified = TRUE),
    prol = list(analyte = "neutrophil progenitor cells", units = "10^9/L", specimen = "tissue", verified = TRUE),
    transit1 = list(analyte = "maturing neutrophils", units = "10^9/L", specimen = "tissue", verified = TRUE),
    transit2 = list(analyte = "maturing neutrophils", units = "10^9/L", specimen = "tissue", verified = TRUE),
    transit3 = list(analyte = "maturing neutrophils", units = "10^9/L", specimen = "tissue", verified = TRUE),
    circ = list(analyte = "neutrophils", units = "10^9/L", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance estimated with the Cockcroft-Gault formula (raw, NOT BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL/F centred on the popPK cohort median 71.6 mL/min (Table 3 row 'Cl cr (med = 71.6 mL/mn) on Cl/F'; Table 1 median 71.6, range 22.1-282.3). No BSA normalization is described, so values are raw mL/min. Missing values were imputed with the cohort median (Section 2.6).",
      source_name = "Clcr"
    ),
    ALP = list(
      description = "Serum alkaline phosphatase activity",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL/F centred on 88.6 U/L as printed in Table 3 ('ALP (med = 88.6 UI/L) on Cl/F'); Table 1 prints the popPK-cohort median as 89 (range 11-819) U/L. Higher ALP lowers CL/F. Missing values were imputed with the cohort median (Section 2.6).",
      source_name = "ALP"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the baseline ANC (Base) centred on 63.7 years as printed in Table 4 ('Age on Base (med = 63.7 years)'); Table 1 prints the PK/PD-cohort median as 63 (range 40-92) years. Older patients have a higher baseline ANC.",
      source_name = "Age"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Section 3.2 / Discussion: tested on CL/F (graphical approach / Wald test) and significant alone, but not retained after backward elimination once CRCL was in the model. Table 1 median 67 (37-140) kg.",
      source_name = "Weight"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Section 3.2 / Discussion: tested on CL/F and significant alone, but not retained after backward elimination. Table 1 median 68.0 (31.0-301.8) umol/L.",
      source_name = "Serum creatinine"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Section 3.3: significant on Base in the forward step but not in the backward step (its effect is thought to be absorbed by age, with which it is negatively correlated). Table 1 PK/PD-cohort median 39.0 (20.0-48.0) g/L.",
      source_name = "Albumin"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 143L,
    n_studies = 1L,
    age_range = "40-92 years (median 69 in the popPK cohort, 63 in the PK/PD cohort)",
    weight_range = "37-140 kg (median 67)",
    sex_female_pct = 100,
    disease_state = "Women treated with palbociclib for HR+/HER2- (or HER2+ non-amplified) advanced or metastatic breast cancer in routine care, combined with fulvestrant (33.6%), letrozole (49.6%) or another hormonotherapy (9.8%); 75.5% with metastases.",
    dose_range = "Oral palbociclib 125 mg (n = 127 samples), 100 mg (n = 38) or 75 mg (n = 16) once daily, 21 days on / 7 days off (28-day cycle).",
    regions = "France (Centre Georges-Francois Leclerc, Dijon)",
    observations = "PK: 181 plasma concentrations from 143 patients (routine therapeutic drug monitoring, 28 October 2018 to 15 March 2021), 0.9-197.25 h after the last dose, 6-229 ug/L (mean 77); multiple samples from one patient were treated as independent. PK/PD: 1508 absolute neutrophil counts from 128 patients (3-41 per patient, mean 11.8) over up to one year (13 cycles), range 0.29-8 G/L.",
    notes = "Baseline characteristics from Tables 1 and 2. Estimated in Monolix 2020R1 (SAEM). The PK/PD model was fitted sequentially: each patient's individual popPK parameters entered the PK/PD fit as regressors (Section 2.4)."
  )

  ini({
    # ---- PK: Table 3, 'Final Population PK Model with Covariates' ----
    lka <- log(0.187); label("ka: first-order absorption rate constant (1/h)") # Table 3: ka = 0.187 /h (RSE 22.2%)
    lcl <- log(57.13); label("CL/F: apparent oral clearance at CRCL = 71.6 mL/min and ALP = 88.6 U/L (L/h)") # Table 3: Cl/F = 57.13 L/h (RSE 2.8%)
    # V/F and Tlag were not estimable (4% of samples within 2 h of dosing) and
    # were fixed to the Royer 2021 one-compartment palbociclib model (Section 3.2).
    lvc <- fixed(log(1580)); label("V/F: apparent volume of distribution (L)") # Table 3: V/F = 1580 L (no RSE; fixed per Section 3.2)
    ltlag <- fixed(log(0.658)); label("Tlag: absorption lag time (h)") # Table 3: 'T lag (h) (fix)' = 0.658 h

    # Covariate effects, Equation 1: Param_i = Param_pop * (COV_i / COV_med)^beta
    e_crcl_cl <- 0.44; label("Exponent of the (CRCL / 71.6 mL/min) power effect on CL/F (unitless)") # Table 3: 'Cl cr (med = 71.6 mL/mn) on Cl/F' = 0.44 (RSE 15.1%)
    e_alp_cl <- -0.14; label("Exponent of the (ALP / 88.6 U/L) power effect on CL/F (unitless)") # Table 3: 'ALP (med = 88.6 UI/L) on Cl/F' = -0.14 (RSE 34.1%)

    # ---- PD: Table 4, 'Final PK/PD Model with Covariates' ----
    lcirc0 <- log(2.92); label("Base: baseline circulating neutrophil count at age 63.7 years (10^9/L)") # Table 4: Base = 2.92 G/L (RSE 3.55%)
    lslope <- log(0.0011); label("Slope: linear palbociclib effect on proliferation (L/ug)") # Table 4: Slope = 0.0011 L/ug (RSE 6.46%)
    lmtt <- log(5.29 * 24); label("MTT: mean transit time through the maturation chain (h)") # Table 4: MTT = 5.29 days (RSE 4.93%) = 126.96 h
    lgamma <- log(0.103); label("gamma: feedback exponent on (Base / ANC) (unitless)") # Table 4: Gamma = 0.103 (RSE 7.01%); the Abstract and Discussion quote 0.102
    e_age_circ0 <- 0.465; label("Exponent of the (AGE / 63.7 years) power effect on Base (unitless)") # Table 4: 'Age on Base (med = 63.7 years)' = 0.465 (RSE 34.00%)

    # ---- Inter-individual variability ----
    # Log-normal IIV (Section 2.5); Tables 3 and 4 give CV% = sqrt(exp(omega^2) - 1),
    # so omega^2 = log(CV^2 + 1).
    etalcl ~ 0.10101 # Table 3: Cl/F IIV = 32.6% -> log(0.326^2 + 1) = 0.10101
    etalcirc0 ~ 0.08393 # Table 4: Base IIV = 29.6% -> log(0.296^2 + 1) = 0.08393
    etalslope ~ 0.07965 # Table 4: Slope IIV = 28.8% -> log(0.288^2 + 1) = 0.07965
    etalmtt ~ 0.03153 # Table 4: MTT IIV = 17.9% -> log(0.179^2 + 1) = 0.03153

    # ---- Residual error ----
    addSd <- 13.84; label("Additive residual error on palbociclib Cc (ug/L)") # Table 3: Additive error = 13.84 ug/L (RSE 20.5%)
    # Section 2.4: ANC log-transformed with an additive error on the log scale.
    expSd_ANC <- 0.34; label("Log-scale (exponential) residual error on ANC (SD, log units)") # Table 4: Exponential error = 0.34 (RSE 2.02%)
  })

  model({
    # ---- 1. PK individual parameters ----
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (CRCL / 71.6)^e_crcl_cl * (ALP / 88.6)^e_alp_cl
    vc <- exp(lvc)
    tlag <- exp(ltlag)
    kel <- cl / vc

    # ---- 2. PK ODEs: one compartment, first-order absorption with lag ----
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag
    # central / vc is mg/L; 1 mg/L = 1000 ug/L.
    Cc <- 1000 * central / vc

    # ---- 3. PD individual parameters ----
    circ0 <- exp(lcirc0 + etalcirc0) * (AGE / 63.7)^e_age_circ0
    slope <- exp(lslope + etalslope)
    mtt <- exp(lmtt + etalmtt)
    gamma <- exp(lgamma)
    # Section 2.4: ktr = (3 + 1) / MTT and kprol = ktr = kcirc.
    ktr <- 4 / mtt

    # ---- 4. Friberg myelosuppression chain (Figure 1) ----
    # Linear drug effect E_D = C * Slope on proliferation; feedback (CIRC0 / CIRCt)^gamma.
    edrug <- slope * Cc
    feed <- (circ0 / circ)^gamma
    d/dt(prol) <- ktr * prol * (1 - edrug) * feed - ktr * prol
    d/dt(transit1) <- ktr * prol - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(circ) <- ktr * transit3 - ktr * circ

    # Section 2.4: at baseline all five compartments equal Base.
    prol(0) <- circ0
    transit1(0) <- circ0
    transit2(0) <- circ0
    transit3(0) <- circ0
    circ(0) <- circ0

    # ---- 5. Observations and error ----
    ANC <- circ
    Cc ~ add(addSd)
    ANC ~ lnorm(expSd_ANC)
  })
}
