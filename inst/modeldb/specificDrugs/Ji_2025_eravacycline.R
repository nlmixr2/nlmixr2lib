Ji_2025_eravacycline <- function() {
  description <- paste(
    "Three-compartment intravenous population PK model for eravacycline in healthy adult volunteers",
    "(Ji 2025), fitted jointly to total plasma concentrations and to urea-corrected epithelial lining",
    "fluid (ELF) concentrations derived from bronchoalveolar lavage. ELF is not a kinetic compartment:",
    "under the authors' quasi-equilibrium simplification the unbound ELF concentration is the unbound",
    "plasma concentration times an estimated ELF distribution ratio of 8.26. Unbound plasma",
    "concentration follows a fixed nonlinear power function of total concentration",
    "(Cu = 0.141 * Ctot^0.873, mg/L) carried from published binding data. Body weight enters",
    "allometrically (fixed exponents 0.75 on clearances and 1 on volumes, referenced to the 81.4 kg",
    "cohort median); no other covariate was retained.",
    sep = " "
  )
  reference <- paste(
    "Ji X-w, Mak WY, Xue F, Yang W-y, Kuan IH-S, Xiang X-q, Li Y, Zhu X.",
    "Population pharmacokinetics and pulmonary modeling of eravacycline and the determination of",
    "microbiological breakpoint and cutoff of PK/PD. Antimicrob Agents Chemother. 2025;69(3):e01065-24.",
    "doi:10.1128/aac.01065-24",
    sep = " "
  )
  vignette <- "Ji_2025_eravacycline"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    # The three ODE states. ELF is NOT a compartment in the final model: the
    # quasi-equilibrium simplification (supplement eq. 10-12) sets dA_ELF/dt = 0,
    # so the two ELF flux terms of the central ODE (supplement eq. 5) cancel and
    # Celf becomes an algebraic rescaling of the unbound plasma concentration.
    central = list(analyte = "eravacycline", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "eravacycline", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "eravacycline", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling on all six disposition parameters, referenced to the cohort median",
        "81.4 kg (Ji 2025 Table 5; range 49.4-105 kg). Exponents are fixed at 0.75 on CL, Q1, Q2",
        "and 1 on Vc, Vp1, Vp2 (Table 1 footnote a). Included a priori as part of the structural",
        "model (supplement 'Covariate model')."
      ),
      source_name = "WT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      notes = "Screened in the stepwise covariate search and not retained (supplement Table S2 note)."
    ),
    ALB = list(
      description = "Serum albumin.",
      units = "not reported (no albumin row in Table 5)",
      type = "continuous",
      notes = "Screened in the stepwise covariate search and not retained (supplement Table S2 note)."
    ),
    ALT = list(
      description = "Alanine aminotransferase.",
      units = "U/L",
      type = "continuous",
      notes = paste(
        "ALT on CL was included in forward selection round 2 (dOFV -10.77) and removed in backward",
        "elimination (supplement Table S2)."
      )
    ),
    CRCL = list(
      description = "Creatinine clearance, not BSA-normalized (Table 5 median 125.7 mL/min).",
      units = "mL/min",
      type = "continuous",
      notes = "Screened in the stepwise covariate search and not retained (supplement Table S2 note)."
    ),
    TPRO = list(
      description = "Serum total protein.",
      units = "g/dL",
      type = "continuous",
      notes = "Screened in the stepwise covariate search and not retained (supplement Table S2 note)."
    ),
    RBC = list(
      description = "Red blood cell count.",
      units = "10^12/L",
      type = "continuous",
      notes = paste(
        "RBC on Q2 and on Vp met the forward-inclusion criterion in round 2 but was not chosen",
        "(supplement Table S2)."
      )
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Sex on Q was included in forward selection round 1 (dOFV -9.18) and removed in backward",
        "elimination round 2 (supplement Table S2)."
      )
    ),
    WBC = list(
      description = "White blood cell count.",
      units = "10^9/L",
      type = "continuous",
      notes = "Screened in the stepwise covariate search and not retained (supplement Table S2 note)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 20,
    n_studies = 1,
    age_range = "19-47 years",
    age_median = "27.5 years",
    weight_range = "49.4-105 kg",
    weight_median = "81.4 kg",
    sex_female_pct = 35,
    race_ethnicity = c(Caucasian = 75, Other = 25),
    disease_state = "healthy volunteers",
    dose_range = "1 mg/kg IV over 1 h every 12 h for 7 doses",
    regions = "United States (Phase I study TP-434-006, NCT01989949)",
    sampling = paste(
      "Plasma before and 1, 2, 4, 6 and 12 h after the last dose; one bronchoalveolar lavage per",
      "subject at 1, 2, 4, 6 or 12 h after the last dose (20 BAL samples, 5 below the quantification limit)."
    ),
    notes = "Ji 2025 Table 5 (demographics), Methods, and supplement section 'Human PK data collection'."
  )

  ini({
    # Structural parameters - typical values at WT = 81.4 kg (Table 1)
    lcl <- log(16.3); label("Clearance (L/h)")                                    # Table 1: CL = 16.3 L/h (RSE 5%; bootstrap 16.2, 15.0-17.4)
    lvc <- log(3.88); label("Central volume of distribution (L)")                 # Table 1: Vc = 3.88 L (RSE 25%; bootstrap 3.5, 1.8-7.2)
    lq <- log(33.9); label("Intercompartmental clearance to peripheral1 (L/h)")   # Table 1: Q1 = 33.9 L/h (RSE 11%; bootstrap 32.7, 24.7-47.7)
    lvp <- log(39); label("First peripheral volume of distribution (L)")          # Table 1: Vp1 = 39 L (RSE 20%; bootstrap 37.6, 29.2-54.5)
    lq2 <- log(23.5); label("Intercompartmental clearance to peripheral2 (L/h)")  # Table 1: Q2 = 23.5 L/h (RSE 16%; bootstrap 23.3, 19.2-29.2)
    lvp2 <- log(122); label("Second peripheral volume of distribution (L)")       # Table 1: Vp2 = 122 L (RSE 15%; bootstrap 121, 109.9-134.9)

    # ELF:plasma unbound concentration ratio, Ratio = Cfree,ELF / Cfree,central
    lfelf <- log(8.26); label("ELF:plasma unbound concentration ratio (unitless)")  # Table 1: Ratio = 8.26 (RSE 12%; bootstrap 8.3, 6.8-9.8); supplement eq. 12

    # Nonlinear plasma protein binding, Cu = 0.141 * Ctot^0.873 with Ctot and Cu
    # in mg/L, carried from published binding data (supplement refs 1-2) and
    # not estimated in this analysis.
    fucoef <- fixed(0.141); label("Unbound plasma concentration at total 1 mg/L (mg/L)")  # Supplement eq. 1 and 9: coefficient 0.141
    fupow <- fixed(0.873); label("Power exponent of total concentration in the unbound-concentration relationship (unitless)")  # Supplement eq. 1 and 9: exponent 0.873

    # Allometric exponents (Table 1 footnote a; fixed, no uncertainty reported)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of WT on CL (unitless)")    # Table 1 footnote a: CL = 16.3 x (WT/81.4)^0.75
    e_wt_vc <- fixed(1); label("Allometric exponent of WT on Vc (unitless)")       # Table 1 footnote a: Vc = 3.88 x (WT/81.4)
    e_wt_q <- fixed(0.75); label("Allometric exponent of WT on Q1 (unitless)")     # Table 1 footnote a: Q1 = 33.9 x (WT/81.4)^0.75
    e_wt_vp <- fixed(1); label("Allometric exponent of WT on Vp1 (unitless)")      # Table 1 footnote a: Vp1 = 39 x (WT/81.4)
    e_wt_q2 <- fixed(0.75); label("Allometric exponent of WT on Q2 (unitless)")    # Table 1 footnote a: Q2 = 23.5 x (WT/81.4)^0.75
    e_wt_vp2 <- fixed(1); label("Allometric exponent of WT on Vp2 (unitless)")     # Table 1 footnote a: Vp2 = 122 x (WT/81.4)

    # IIV - exponential (supplement eq. 13). Table 1 prints bare percentages,
    # read as lognormal CVs: omega^2 = log(1 + CV^2).
    etalcl ~ 0.01832703  # log(1 + 0.136^2); Table 1 IIV CL = 13.6% (RSE 69%, shrinkage 0%)
    etalvc ~ 0.8525412   # log(1 + 1.16^2); Table 1 IIV Vc = 116% (RSE 49%, shrinkage 19%)
    etalvp ~ 0.2642855   # log(1 + 0.55^2); Table 1 IIV Vp1 = 55% (RSE 56%, shrinkage 8%)

    # Residual error - plasma (combined) and ELF (proportional)
    propSd <- 0.0704; label("Proportional residual SD for plasma Cc (fraction)")  # Table 1 serum: proportional error = 7.04% (RSE 21%)
    addSd <- 8.73 / 1000; label("Additive residual SD for plasma Cc (mg/L)")       # Table 1 serum: additive error = 8.73 (no unit printed; read as ng/mL, /1000 to mg/L) (RSE 48%)
    propSd_Celf <- 0.423; label("Proportional residual SD for ELF Celf (fraction)") # Table 1 BAL: proportional error = 42.3% (RSE 28%)
  })

  model({
    # Reference weight: cohort median (Table 5), the normaliser of Table 1 footnote a
    wt_ref <- 81.4

    cl <- exp(lcl + etalcl) * (WT / wt_ref)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / wt_ref)^e_wt_vc
    q <- exp(lq) * (WT / wt_ref)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / wt_ref)^e_wt_vp
    q2 <- exp(lq2) * (WT / wt_ref)^e_wt_q2
    vp2 <- exp(lvp2) * (WT / wt_ref)^e_wt_vp2
    felf <- exp(lfelf)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Supplement eq. 5-7 with the ELF flux terms removed: at quasi-equilibrium
    # (eq. 10) k_cl * A1,free = k_lc * A_ELF, so the two terms cancel exactly.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # Total plasma concentration (mg/L)
    Cc <- central / vc

    # Unbound plasma concentration, supplement eq. 1 / 9 (mg/L)
    Cu <- fucoef * Cc^fupow

    # Unbound ELF concentration, supplement eq. 12. Observed ELF concentrations
    # are urea-corrected BAL concentrations, taken as unbound (supplement
    # 'ELF concentration determined from BAL').
    Celf <- felf * Cu

    Cc ~ add(addSd) + prop(propSd)
    Celf ~ prop(propSd_Celf)
  })
}
