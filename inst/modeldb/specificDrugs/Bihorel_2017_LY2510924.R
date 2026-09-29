Bihorel_2017_LY2510924 <- function() {
  description <- "Coupled population PK/PD model for LY2510924, a cyclic-peptide CXCR4 antagonist given as once-daily subcutaneous injections to 227 patients with advanced cancer across three studies. PK is a two-compartment model with first-order absorption and an apparent clearance that DECREASES with the administered daily dose along a decreasing sigmoid (a high-dose asymptote CLmin/F plus a span CLdelta/F removed with half-effect dose50, truncated below 1 mg), with allometric body-weight scaling on CL/F and V2/F. The PD layer is a precursor-dependent indirect-response model for blood CD34+ cell counts: a precursor pool exchanges cells reversibly with the circulating compartment, and LY2510924 stimulates the pool-to-blood mobilisation rate through a saturable Emax term. A second, LY2510924-independent stimulation driven by an empirical first-order signal build-up applies only to patients in Study CXAC, standing in for the etoposide/carboplatin standard of care and concomitant G-CSF/EPO."
  reference <- paste(
    "Bihorel S, Raddad E, Fiedler-Kelly J, Stille JR, Hing J, Ludwig E.",
    "Population Pharmacokinetic and Pharmacodynamic Modeling of LY2510924",
    "in Patients With Advanced Cancer.",
    "CPT Pharmacometrics Syst Pharmacol. 2017;6(9):614-624.",
    "doi:10.1002/psp4.12221.",
    sep = " "
  )
  vignette <- "Bihorel_2017_LY2510924"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL",
    cd34 = "cells/uL"
  )

  # `signal` is the paper's empirical latent state S (Eq. 5): a dimensionless
  # first-order build-up saturating at 1, standing in for the LY2510924-
  # independent CD34+ rise seen in Study CXAC. It is not a mass compartment and
  # has no canonical counterpart.
  paper_specific_compartments <- c("signal")

  compartmentData <- list(
    depot = list(analyte = "LY2510924", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "LY2510924", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "LY2510924", units = "mg", specimen = "plasma", verified = TRUE),
    precursor1 = list(
      analyte = "CD34+ progenitor cells",
      units = "cells/uL",
      specimen = "not applicable",
      verified = TRUE
    ),
    circ = list(analyte = "CD34+ progenitor cells", units = "cells/uL", specimen = "whole blood", verified = TRUE),
    signal = list(
      analyte = "empirical stimulatory signal",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling of CL/F (Eq. 1) and V2/F (Eq. 2), each normalised to 80.1 kg, the median body weight of the analysis population (Table 1 overall mean 84.40 kg, range 39.6-167.8 kg).",
      source_name = "WTKG"
    ),
    DOSE_LY2510924_MGD = list(
      description = "Administered once-daily subcutaneous LY2510924 dose",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives the decreasing-sigmoid dose dependence of CL/F in Eq. 1. This is the DOSE LEVEL as a covariate, not the amount of the current dosing record: CL/F is a function of the daily dose the patient is assigned to, so the column is constant within a dosing regimen and changes only when the assigned daily dose changes. Eq. 1 is truncated at 1 mg (the lowest dose studied, giving the highest clearance), so values below 1 mg/day give the same CL/F as 1 mg/day. Studied levels: 1, 2.5, 5, 10, 20 and 30 mg/day (Study CXAA dose escalation); 2.5 or 20 mg/day (CXAA dose confirmation); 20 mg/day (Studies CXAB and CXAC).",
      source_name = "dose"
    ),
    STUDY_CXAA = list(
      description = "Enrolment in phase 1 Study I2V-MC-CXAA",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (phase 2 Studies I2V-MC-CXAB and I2V-MC-CXAC)",
      notes = "Selects which of the two LY2510924 residual-error models applies. The paper estimated separate additive-plus-constant-CV residual-error models for phase 1 and phase 2 studies (Results, PK model section; Table 2 rows 'Phase 1 RV' and 'Phase 2 RV'). Study CXAA is the only phase 1 study, so this indicator is exactly the phase 1 / phase 2 split. Not a structural covariate: no PK or PD parameter depends on it.",
      source_name = "study"
    ),
    STUDY_CXAC = list(
      description = "Enrolment in phase 2 Study I2V-MC-CXAC (extensive-stage small cell lung carcinoma, etoposide/carboplatin standard of care)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Studies I2V-MC-CXAA and I2V-MC-CXAB)",
      notes = "The paper's FCXAC indicator (Eq. 5), which gates the empirical signal build-up that drives the LY2510924-independent CD34+ rise. FCXAC is 1 after the first dose for Study CXAC patients and 0 otherwise; because `signal` starts at 0 and the model's time origin is the first dose, the plain indicator reproduces it. Both arms of Study CXAC carry the signal, as the effect could not be separated from the standard of care, cancer type, or study (Methods).",
      source_name = "FCXAC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 227L,
    n_studies = 3L,
    age_range = "29-85 years",
    age_median = "64.3 years (mean; SD 9.8)",
    weight_range = "39.6-167.8 kg",
    weight_median = "80.1 kg (median of the analysis population, used as the allometric reference; Table 1 overall mean 84.40 kg, SD 22.00)",
    sex_female_pct = 44.5,
    race_ethnicity = c(
      Caucasian = 90.7,
      `Black/African American` = 7.5,
      `American Indian/Alaskan Native` = 0.4,
      Unknown = 1.3
    ),
    disease_state = "Patients with advanced and/or metastatic cancer. Study CXAA (phase 1, n = 39) enrolled a broad mix of tumour types, mostly gastrointestinal (46.2%), lung (15.4%) and genitourinary (12.8%). Study CXAB (phase 2, n = 100) enrolled treatment-naive metastatic renal cell carcinoma, randomised 2:1 to LY2510924 plus sunitinib or sunitinib alone. Study CXAC (phase 2, n = 88) enrolled treatment-naive extensive-stage small cell lung carcinoma, randomised 1:1 to LY2510924 plus etoposide/carboplatin or etoposide/carboplatin alone.",
    dose_range = "Once-daily subcutaneous LY2510924. Study CXAA: 1, 2.5, 5, 10, 20 or 30 mg/day for 28-day cycles (dose escalation), then 2.5 or 20 mg/day (dose confirmation). Studies CXAB and CXAC: 20 mg/day.",
    co_medication = "Study CXAB: sunitinib 50 mg once daily orally for the first 4 weeks of each 6-week cycle. Study CXAC: carboplatin targeting AUC 5 mg/mL/min intravenously on day 1 of each 21-day cycle plus etoposide 100 mg/m2 intravenously on days 1-3. G-CSF was used by 12.3% and erythropoietin by 2.6% of patients overall, all of them in Studies CXAB and CXAC (Table 1).",
    notes = "Demographics from Table 1. 767 LY2510924 plasma samples from 147 patients entered the PK analysis and 1,042 CD34+ cell counts from 227 patients entered the PK/PD analysis; patients with only a baseline CD34+ measurement were included in the PK/PD but not the PK analysis. No patient developed anti-drug antibodies. Baseline CD34+ cell count overall mean 1.8 cells/uL (SD 1.6, range 0-13, n = 211)."
  )

  ini({
    # ---------------------------------------------------------------- PK ----
    # Eq. 1 -- apparent elimination clearance is a DECREASING sigmoid function
    # of the administered daily dose: CL/F falls from a maximum of
    # (CLmin/F + CLdelta/F) at doses <= 1 mg towards the asymptote CLmin/F at
    # high doses, with half the span removed at dose50. Truncated at 1 mg, the
    # lowest dose studied (Discussion).
    lcl_dosemin <- log(4.75)
    label("Apparent elimination clearance at high dose, the asymptote TVCLmin/F (L/h)") # Table 2, TVCLmin/F = 4.75 L/h (%SEM 7.09)
    lcl_dosespan <- log(12.9)
    label("Span of apparent elimination clearance removed by increasing dose, TVCLdelta/F (L/h)") # Table 2, TVCLdelta/F = 12.9 L/h (%SEM 11.1)
    lcl_dose50 <- fixed(log(3.60))
    label("Daily dose removing half the clearance span, dose50 (mg)") # Table 2, dose50 = 3.60 mg, FIXED (fixed to the Study CXAA base-model estimate; Results)
    e_wt_cl <- 0.870
    label("Allometric exponent of body weight on CL/F (unitless)") # Table 2, bCL = 0.870 (%SEM 14.8)

    lvc <- log(35.0)
    label("Apparent central volume of distribution TVV2/F (L)") # Table 2, TVV2/F = 35.0 L (%SEM 4.89)
    e_wt_vc <- 0.948
    label("Allometric exponent of body weight on V2/F (unitless)") # Table 2, bV = 0.948 (%SEM 14.8)
    lq <- log(3.74)
    label("Apparent distribution clearance TVQ/F (L/h)") # Table 2, TVQ/F = 3.74 L/h (%SEM 18.8)
    lvp <- log(21.9)
    label("Apparent peripheral volume of distribution TVV3/F (L)") # Table 2, TVV3/F = 21.9 L (%SEM 6.67)
    lka <- fixed(log(10.0))
    label("First-order subcutaneous absorption rate TVKA (1/h)") # Table 2, TVKA = 10.0 1/h, FIXED (estimates were > 40 1/h and poorly estimated; Discussion)

    # ---------------------------------------------------------------- PD ----
    # Eq. 5-6 -- precursor-dependent indirect response with REVERSIBLE transfer
    # between the precursor pool and the circulating compartment.
    lrbase <- log(1.45)
    label("Baseline blood CD34+ cell count CD34_0 (cells/uL)") # Table 2, CD34_0 = 1.45 cells/uL (%SEM 5.66)
    lkpc_cell <- log(42.4e-6)
    label("First-order transfer rate of CD34+ cells from the precursor pool to blood, Kpc (1/h)") # Table 2, Kpc = 42.4 in units of 1000000/h, i.e. 42.4e-6 1/h (%SEM 25.1); see the units note in model()
    lkcp_cell <- log(0.185)
    label("First-order transfer rate of CD34+ cells from blood back to the precursor pool, Kcp (1/h)") # Table 2, Kcp = 0.185 1/h (%SEM 20.4)
    lkout <- log(0.0104)
    label("First-order elimination rate of circulating CD34+ cells, Kout (1/h)") # Table 2, Kout = 0.0104 1/h (%SEM 51.1)
    lemax <- log(13.1)
    label("Maximum fractional stimulation of CD34+ cell mobilisation by LY2510924, Smax (unitless)") # Table 2, Smax = 13.1 (%SEM 14.4)
    lec50 <- log(6.87)
    label("LY2510924 concentration giving half-maximal stimulation of mobilisation, SC50 (ng/mL)") # Table 2, SC50 = 6.87 ng/mL (%SEM 45.1)
    lksig <- log(0.00804)
    label("First-order build-up rate of the empirical Study CXAC signal, Kt (1/h)") # Table 2, Kt = 0.00804 1/h (%SEM 15.5)
    lsstim_signal_mob <- log(5.63)
    label("Stimulatory factor of the empirical signal on CD34+ mobilisation, alpha (unitless)") # Table 2, alpha = 5.63 (%SEM 18.2)

    # ------------------------------------------------- Interindividual ------
    # Exponential IIV (Eqs. 3-4 for PK; 'exponential variability models ... in a
    # simple diagonal matrix form', Methods). Table 2 reports magnitudes as %CV;
    # the variances below use the exact log-normal relation
    # omega^2 = log(1 + CV^2).
    etalcl ~ 0.0851 # Table 2, CL/F IIV 29.8%CV (shrinkage 17.7%) -> log(1 + 0.298^2)
    etalvc ~ 0.0684 # Table 2, V2/F IIV 26.6%CV (shrinkage 36.0%) -> log(1 + 0.266^2)
    etalrbase ~ 0.4521 # Table 2, CD34_0 IIV 75.6%CV (shrinkage 7.06%) -> log(1 + 0.756^2)
    etalkout ~ 1.7647 # Table 2, Kout IIV 220%CV (shrinkage 63.3%) -> log(1 + 2.20^2)
    etalemax ~ 0.1919 # Table 2, Smax IIV 46.0%CV (shrinkage 48.0%) -> log(1 + 0.460^2)
    etalec50 ~ 1.0280 # Table 2, SC50 IIV 134%CV (shrinkage 62.1%) -> log(1 + 1.34^2)

    # ------------------------------------------------- Residual error -------
    # Table 2 reports each additive-plus-constant-CV residual error as the two
    # endpoints of the %CV over the concentration range given in the table
    # footnote. Both endpoints are printed values; the (addSd, propSd) pair
    # below is their exact algebraic inversion of
    #   %CV(C) = 100 * sqrt(addSd^2 + (propSd * C)^2) / C
    # solved at the two stated concentrations. See vignette Errata.
    propSd <- 0.2221
    label("LY2510924 proportional residual error, phase 1 Study CXAA (fraction)") # Table 2, Phase 1 RV '249 - 22.3%CV' over 0.2-25 ng/mL (footnote a)
    addSd <- 0.4960
    label("LY2510924 additive residual error, phase 1 Study CXAA (ng/mL)") # Table 2, Phase 1 RV '249 - 22.3%CV' over 0.2-25 ng/mL (footnote a)
    propSd_phase2 <- 0.4058
    label("LY2510924 proportional residual error, phase 2 Studies CXAB and CXAC (fraction)") # Table 2, Phase 2 RV '730 - 41.0%CV' over 0.2-25 ng/mL (footnote a)
    addSd_phase2 <- 1.4577
    label("LY2510924 additive residual error, phase 2 Studies CXAB and CXAC (ng/mL)") # Table 2, Phase 2 RV '730 - 41.0%CV' over 0.2-25 ng/mL (footnote a)

    propSd_CD34 <- 0.4669
    label("CD34+ cell count proportional residual error (fraction)") # Table 2, PD RV '396 - 46.7%CV' over 0.6-200 cells/uL (footnote b)
    addSd_CD34 <- 2.3594
    label("CD34+ cell count additive residual error (cells/uL)") # Table 2, PD RV '396 - 46.7%CV' over 0.6-200 cells/uL (footnote b)
  })

  model({
    # ---- 1. Dose-dependent apparent clearance, Eq. 1 -----------------------
    # Eq. 1 is printed as a two-branch function of the daily dose:
    #   dose <= 1 mg : CLmin/F + CLdelta/F
    #   dose >  1 mg : CLmin/F + CLdelta/F - CLdelta/F*(dose-1)/((dose50-1)+(dose-1))
    # `doseAbove1` is the truncated excess over 1 mg, so the single expression
    # below reproduces BOTH branches exactly (it collapses to the first when
    # doseAbove1 is 0).
    cl_dosemin <- exp(lcl_dosemin)
    cl_dosespan <- exp(lcl_dosespan)
    cl_dose50 <- exp(lcl_dose50)
    doseAbove1 <- (DOSE_LY2510924_MGD - 1) * (DOSE_LY2510924_MGD > 1)
    clTypical <- (WT / 80.1)^e_wt_cl *
      (cl_dosemin + cl_dosespan -
        cl_dosespan * doseAbove1 / ((cl_dose50 - 1) + doseAbove1))

    # ---- 2. Individual PK parameters, Eqs. 2-4 -----------------------------
    cl <- clTypical * exp(etalcl) # Eq. 3
    vc <- exp(lvc + etalvc) * (WT / 80.1)^e_wt_vc # Eqs. 2 and 4
    q <- exp(lq)
    vp <- exp(lvp)
    ka <- exp(lka)

    # ---- 3. Micro-constants ------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- 4. PK ODEs (Figure 1) ---------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Doses are in mg and volumes in L, so central/vc is mg/L == ug/mL;
    # multiply by 1000 to give ng/mL, the unit of the PK observations and of
    # SC50 (units$concentration).
    Cc <- (central / vc) * 1000

    # ---- 5. Individual PD parameters ---------------------------------------
    rbase <- exp(lrbase + etalrbase)
    kout <- exp(lkout + etalkout)
    emax <- exp(lemax + etalemax)
    ec50 <- exp(lec50 + etalec50)
    ksig <- exp(lksig)
    sstim_signal_mob <- exp(lsstim_signal_mob)

    # Table 2 prints the unit of Kpc as '1000000/h', i.e. the tabulated 42.4 is
    # the rate multiplied by 1e6, so Kpc = 42.4e-6 1/h. The scaling is forced by
    # the model's own behaviour and not a free reading: at 42.4e-6 1/h the
    # printed initial condition puts about 6,700 cells/uL in the precursor pool,
    # whose turnover time Kin/P(0) is far longer than a treatment cycle, which
    # is what produces the slowly developing, small-magnitude tolerance the
    # paper reports. At 42.4 1/h the pool would hold 0.0067 cells/uL and empty
    # within the hour, giving instantaneous and complete tolerance and no
    # sustained response at all -- the opposite of every result in the paper.
    kpc_cell <- exp(lkpc_cell)
    kcp_cell <- exp(lkcp_cell)

    kin <- rbase * kout # Eq. 6

    # ---- 6. PD ODEs (Eq. 5) ------------------------------------------------
    # `mobilisation` is the stimulated pool-to-blood transfer rate: the baseline
    # Kpc multiplied by 1 plus the saturable LY2510924 effect plus the empirical
    # Study CXAC signal term.
    mobilisation <- kpc_cell * (1 + emax * Cc / (ec50 + Cc) + sstim_signal_mob * signal)

    d/dt(precursor1) <- kin - mobilisation * precursor1 + kcp_cell * circ
    # Eq. 5 as typeset prints '+ (Kcp + Kout) x R' on this line. That sign is a
    # typographical error in the published equation: cells leave the circulating
    # compartment both by recycling to the pool (Kcp) and by elimination (Kout),
    # so the term is a loss. The minus is forced independently by the paper's own
    # printed initial conditions together with Eq. 6 -- P(0) and R(0) below are
    # exactly the drug-free steady state of the minus form, and with a plus the
    # system has no steady state at all and R diverges. Minus signs render
    # correctly elsewhere in the same display equation ('Kin - Kpc', '1 - S'),
    # so this is the paper's typo and not a glyph artefact. See vignette Errata.
    d/dt(circ) <- mobilisation * precursor1 - (kcp_cell + kout) * circ
    d/dt(signal) <- ksig * STUDY_CXAC * (1 - signal)

    # Drug-free steady state assumed before the first dose (Eq. 5 initial
    # conditions). Note P(0) uses the UNSTIMULATED Kpc.
    precursor1(0) <- rbase * (kout + kcp_cell) / kpc_cell
    circ(0) <- rbase
    signal(0) <- 0

    # ---- 7. Observations and residual error --------------------------------
    CD34 <- circ

    # Separate additive-plus-constant-CV residual-error models were estimated
    # for the phase 1 study (CXAA) and the phase 2 studies (CXAB, CXAC)
    # (Results, PK model section). STUDY_CXAA selects between them.
    propSdPk <- propSd * STUDY_CXAA + propSd_phase2 * (1 - STUDY_CXAA)
    addSdPk <- addSd * STUDY_CXAA + addSd_phase2 * (1 - STUDY_CXAA)

    Cc ~ add(addSdPk) + prop(propSdPk)
    CD34 ~ add(addSd_CD34) + prop(propSd_CD34)
  })
}
