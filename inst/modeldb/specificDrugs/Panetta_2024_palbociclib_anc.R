Panetta_2024_palbociclib_anc <- function() {
  description <- "Semi-mechanistic PK/PD model of palbociclib-induced neutropenia in children and young adults (4.9-21.6 years) with recurrent, progressive or refractory brain tumors treated on the Pediatric Brain Tumor Consortium phase 1 study PBTC-042 (Panetta 2024, Table 3A). The palbociclib PK layer is the one-compartment first-order-absorption model with an absorption lag time and an AST power effect on CL/F reported in Table 2 (also available on its own as modellib('Panetta_2024_palbociclib')), converted from hours to days. Plasma palbociclib inhibits the proliferation rate of a bone-marrow pool of neutrophil precursors through a fractional Imax term IC50^n / (IC50^n + C^n); the pool feeds three maturation transit compartments and then the circulating neutrophil compartment, which exerts G-CSF-like negative feedback on proliferation through a KM / (KM + ANC) term. KM and the drug-free initial conditions are not estimated: they are fixed by the paper's steady-state requirement (Equations 6 and 7). The transit rate constant is parameterised as the estimated fraction fk_bp of the proliferation rate kin, on a logit scale, which guarantees kin > k_bp for every simulated subject. The circulating-neutrophil elimination rate kout and the Hill coefficient n are fixed. Volumes and clearances are per m^2 of body-surface area and doses are supplied in mg/m^2, so BSA cancels and no BSA covariate is needed. The platelet counterpart is modellib('Panetta_2024_palbociclib_plt')."
  reference <- paste(
    "Panetta JC, Selvo NS, Van Mater D, Stewart CF. (2024).",
    "Population Pharmacokinetic and Pharmacodynamic Study of Palbociclib in Children and Young Adults",
    "with Recurrent, Progressive, or Refractory Brain Tumors.",
    "Pharmaceutics 16(12):1528. doi:10.3390/pharmaceutics16121528.",
    "The myelosuppression structure was previously described in",
    "Panetta JC, Schaiquevich P, Santana VM, Stewart CF. (2008) Using pharmacokinetic and",
    "pharmacodynamic modeling and simulation to evaluate importance of schedule in topotecan therapy",
    "for pediatric neuroblastoma. Clin Cancer Res 14(1):318-325. doi:10.1158/1078-0432.CCR-07-1243,",
    "itself adapted from a temozolomide myelosuppression model.",
    "The underlying phase 1 trial is PBTC-042 (NCT02255461).",
    sep = " "
  )
  vignette <- "Panetta_2024_palbociclib"
  units <- list(
    time          = "day",
    dosing        = "mg/m^2",
    concentration = "ng/mL",
    anc           = "10^3/uL"
  )
  # Unit notes.
  # 1. TIME. Panetta 2024 reports the PK parameters per hour (Table 2) and the
  #    PD parameters per day (Table 3A). This model runs in DAYS, so the PK
  #    parameters are converted inline (ka and CL/F multiplied by 24, Tlag
  #    divided by 24) with the arithmetic written out so it is auditable.
  #    Log-scale variances are invariant to that rescaling and are unchanged.
  # 2. DOSE AND VOLUME. V/F is L/m^2 and CL/F is L/h/m^2 (Table 2); section 2.3
  #    states the BSA-normalized dosage was the model input, so doses are in
  #    mg/m^2, BSA cancels, and central / vc is in mg/L. Cc multiplies by the
  #    1000 ng/mL per mg/L factor to reach the ng/mL that Panetta 2024 reports.
  # 3. PK-PD LINK. Equation 1 is driven by Cc in ng/mL. Panetta 2024 labels
  #    IC50 and C_palbociclib "nM" (Table 3A), but that reading contradicts the
  #    paper's own Supplementary Figure S2, Figure 4 and Table 4; see the IC50
  #    entry in ini() for the three checks that settle it.

  compartmentData <- list(
    depot      = list(analyte = "palbociclib", units = "mg/m^2", specimen = "administration site", verified = TRUE),
    central    = list(analyte = "palbociclib", units = "mg/m^2", specimen = "plasma", verified = TRUE),
    precursor1 = list(analyte = "neutrophils", units = "10^3/uL", specimen = "tissue", verified = TRUE),
    precursor2 = list(analyte = "neutrophils", units = "10^3/uL", specimen = "tissue", verified = TRUE),
    precursor3 = list(analyte = "neutrophils", units = "10^3/uL", specimen = "tissue", verified = TRUE),
    precursor4 = list(analyte = "neutrophils", units = "10^3/uL", specimen = "tissue", verified = TRUE),
    circ       = list(analyte = "neutrophils", units = "10^3/uL", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    AST = list(
      description        = "Baseline serum aspartate aminotransferase activity.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL/F centred on the cohort median AST of 25 U/L (Table 1), the only covariate retained by the pharmacokinetic covariate analysis (Table 2, AST column; Table S1). See modellib('Panetta_2024_palbociclib') for the full covariate-screen description.",
      source_name        = "AST"
    ),
    OCC = list(
      description        = "Treatment-course index used by the inter-occasion random effects.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Panetta 2024 estimated the pharmacodynamic parameters on courses 1 and 2 (section 3.4), so OCC = 1 is course 1 and OCC = 2 is course 2. The paper reports inter-occasion variability on kin, k_bp and IC50 (Table 3A) but never prints an occasion column name; OCC is the nlmixr2lib canonical. The CL/F inter-occasion term of the PK layer is keyed to the same OCC here: in the PK analysis its occasions were the course 1 day 1 and course 1 day 21 sampling visits, roughly three weeks apart and therefore comparable in span to a treatment course. For a single-occasion simulation set OCC = 1 throughout.",
      source_name        = NA_character_
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 31L,
    n_studies      = 1L,
    age_range      = "4.9-21.6 years (median 12.8)",
    weight_range   = "23.8-110.4 kg (median 52.5)",
    bsa_range      = "0.8-2.4 m^2 (median 1.5)",
    sex_female_pct = 36,
    race_ethnicity = c(Caucasian = 58, Black = 12, Other = 30),
    disease_state  = "Recurrent, progressive or refractory pediatric brain tumors. Stratum I (not heavily pretreated) versus stratum II (heavily pretreated); simulations at 75 mg/m^2 predicted grade 3 or greater neutropenia in 75% of stratum II versus 35% of stratum I patients (Figure S4).",
    dose_range     = "Oral palbociclib 50, 75 or 95 mg/m^2 once daily for 21 days of each 28-day course.",
    regions        = "United States (Pediatric Brain Tumor Consortium)",
    observations   = "190 absolute neutrophil counts from 31 patients during course 1 and 99 during course 2; 54 further counts from six patients during courses 3 and 4 were held out to test prediction. Of the 33 patients in the PK analysis, 2 were excluded from the PD analysis because they had only a single pre-treatment ANC observation. A complete blood count with differential was drawn pre-therapy and then weekly, averaging 11 ANC observations per patient.",
    notes              = "Estimated in Monolix 2023R1 by SAEM with the palbociclib pharmacokinetics fixed to each individual's post hoc estimates (section 2.4.2). Observed outcomes from the individual model estimates over courses 1 and 2: median (5th-95th percentile) time below 1.5, 1.0 and 0.5 x 10^3/uL was 15 (0, 35), 4 (0, 21) and 0 (0, 11) days per course, and the median nadir was 1.0 (0.4, 3.3) x 10^3/uL. Courses 3 and 4 were predicted without bias (median residual -6.1%, p = 0.22), i.e. the neutropenia was reversible and not cumulative."
  )

  ini({
    # =========================================================================
    # PALBOCICLIB PK LAYER. Panetta 2024 Table 2, AST column (the final model).
    # Section 2.4.2: "The palbociclib pharmacokinetics were fixed to each
    # individual's post hoc estimated parameters", so these are the same
    # population values as modellib('Panetta_2024_palbociclib'), converted from
    # hours to days.
    # =========================================================================
    ltlag <- log(0.8 / 24); label("Tlag: absorption lag time (day)")                   # Table 2, AST column: Tlag = 0.8 h  = 0.033333 day
    lka   <- log(0.46 * 24); label("ka: first-order absorption rate constant (1/day)") # Table 2, AST column: ka = 0.46 /h  = 11.04 /day
    lvc   <- log(653); label("V/F: apparent oral volume per BSA (L/m^2)")              # Table 2, AST column: V/F = 653 L/m^2
    lcl   <- log(36.5 * 24); label("CL/F: apparent oral clearance per BSA (L/day/m^2)") # Table 2, AST column: CL/F = 36.5 L/h/m^2 = 876 L/day/m^2

    e_ast_cl <- -0.3; label("Exponent of the (AST / 25 U/L) power effect on CL/F (unitless)")  # Table 2, AST column: AST on CL/F = -0.3 (RSE 39.3%); p = 0.0066

    # =========================================================================
    # NEUTROPHIL PD LAYER. Panetta 2024 Table 3A.
    # =========================================================================
    lrbase <- log(3.95); label("Ncirc0: drug-free steady-state circulating neutrophil count (10^3/uL)")  # Table 3A: Ncirc0 = 3.95 x 10^3/uL (RSE 8.7%)
    lkin   <- log(0.90); label("kin: proliferation rate of the bone-marrow neutrophil pool (1/day)")     # Table 3A: kin = 0.90 /day (RSE 8.6%)
    # IC50 UNITS. Panetta 2024 labels IC50 and C_palbociclib "nM" (Table 3A and
    # the Equation 1-5 glossary), but the value is only self-consistent with
    # the paper's own PK model and simulations when read as ng/mL -- the unit
    # the PK model actually produces. Three independent checks, all failed by
    # the nM reading (which needs the molecular weight 447.54 g/mol, i.e.
    # 1 ng/mL = 2.2344 nM):
    # Note the direction: because the model is driven by Cc in ng/mL, reading
    # the tabulated number as nM makes the drug MORE potent on that scale
    # (383.5 nM = 383.5 / 2.2344 = 171.6 ng/mL), so it predicts DEEPER
    # depletion.
    #   1. Supplementary Figure S2 plots the model-estimated concentration on a
    #      log axis against a left axis labelled nM, with the curves peaking
    #      around 10^2. This model's steady-state Cmax at 75 mg/m^2 is
    #      121 ng/mL but 271 nM, so the plotted numbers are the ng/mL values.
    #   2. This model's own typical ANC nadirs at 50, 75 and 95 mg/m^2 are
    #      1.14, 0.59 and 0.36 x 10^3/uL under the ng/mL reading but 0.24,
    #      0.07 and 0.03 under the nM reading, i.e. the circulating pool is
    #      essentially wiped out. Table 4 reports only 19% of individuals
    #      below 1.0 x 10^3/uL at 50 mg/m^2 and only 3.1, 13 and 32% below
    #      0.5 x 10^3/uL across the three dosages, which the nM reading
    #      cannot produce. The platelet model rejects it independently
    #      (typical nadirs 191, 163 and 143 x 10^9/L under the ng/mL reading
    #      against roughly 205, 175 and 160 read off the Figure 4D-F median
    #      curves, versus 128, 89 and 67 under the nM reading).
    #   3. Table 4 reports only 6, 8 and 12% of individuals dropping below
    #      100 x 10^9/L at 50, 75 and 95 mg/m^2, which requires a typical
    #      platelet nadir well above 100; the nM reading puts the typical
    #      subject below it at every dosage.
    # The value below is therefore the tabulated 383.5 with the units read as
    # ng/mL. This is a reporting error in the source, recorded in the
    # vignette's assumptions section.
    lic50  <- log(383.5); label("IC50: palbociclib concentration halving marrow proliferation (ng/mL)")  # Table 3A: IC50 = 383.5 (RSE 17.8%), printed as nM; see the units note above

    # k_bp is not estimated directly. Section 2.4.1: "For this relationship to
    # be positive, kin > k_bp. To ensure this requirement, we defined
    # k_bp = fk_bp * kin where fk_bp in [0, 1] and estimated fk_bp", and
    # section 2.4.2 states fk_bp was given a logit distribution rather than the
    # log-normal used for every other parameter. Table 3A tabulates the derived
    # typical k_bp = 0.78 /day, so the typical fraction is
    #   fk_bp = 0.78 / 0.90 = 0.866667,  logit(0.866667) = log(6.5) = 1.8718.
    # Cross-check: the Discussion quotes a mean transit time of 5.12 days and
    # section 2.4.1 defines it as 4 / k_bp; 4 / 0.78 = 5.128 days.
    logitfktr <- qlogis(0.78 / 0.90)
    label("logit(fk_bp): logit of the transit rate constant expressed as a fraction of kin (unitless)")  # Table 3A: k_bp = 0.78 /day (RSE 3.5%) with kin = 0.90 /day

    lkout <- fixed(log(4.80)); label("kout: elimination rate of circulating neutrophils (1/day)")  # Table 3A: kout = 4.80 /day, reported as Fixed
    lhill <- fixed(log(1));    label("n: Hill coefficient of the palbociclib marrow effect (unitless)")  # Table 3A: n = 1, reported as Fixed

    # =========================================================================
    # Inter-individual and inter-occasion variability. Table 3A reports these
    # as CV%, and section 2.4.2 states the random effects are log-normal (all
    # parameters except fk_bp), so omega^2 = log(CV^2 + 1). The inversion is
    # written out inline so it is auditable.
    # =========================================================================
    etalrbase ~ log(0.470^2 + 1)  # Table 3A: IIV Ncirc0 = 47.0% (RSE 13.9%) -> omega^2 = 0.19958
    etalkin   ~ log(0.294^2 + 1)  # Table 3A: IIV kin    = 29.4% (RSE 50.2%) -> omega^2 = 0.08290
    etalic50  ~ log(0.256^2 + 1)  # Table 3A: IIV IC50   = 25.6% (RSE 53.0%) -> omega^2 = 0.06346

    # fk_bp is logit-distributed, so its variability is a logit-scale variance
    # rather than log(CV^2 + 1). Table 3A's CV is converted with the delta
    # method for the logit transform, d logit(p) / dp = 1 / (p (1 - p)), which
    # for a CV expressed relative to p gives sd_logit = CV / (1 - p) with
    # p = 0.866667:
    #   IIV: 0.015 / 0.133333 = 0.1125,  IOV: 0.058 / 0.133333 = 0.4350.
    # Both are checked by round-trip: expit(1.8718 +/- 0.1125) gives k_bp
    # between 0.767 and 0.792 /day, a 1.5% CV, and expit(1.8718 +/- 0.4350)
    # gives 0.727 to 0.819 /day, a 5.8% CV.
    etalogitfktr ~ 0.1125^2  # Table 3A: IIV k_bp = 1.5% (RSE 203.4%)

    etaiov_kin_1  ~ log(0.143^2 + 1)
    label("Inter-occasion variability in kin, occasion 1 (log-scale variance)")   # Table 3A: IOV kin  = 14.3% (RSE 54.3%) -> 0.02024
    etaiov_kin_2  ~ fix(log(0.143^2 + 1))                                         # same IOV magnitude on occasion 2
    etaiov_ic50_1 ~ log(0.242^2 + 1)
    label("Inter-occasion variability in IC50, occasion 1 (log-scale variance)")  # Table 3A: IOV IC50 = 24.2% (RSE 149.4%) -> 0.05692
    etaiov_ic50_2 ~ fix(log(0.242^2 + 1))                                         # same IOV magnitude on occasion 2
    etaiov_fktr_1 ~ 0.4350^2
    label("Inter-occasion variability in logit(fk_bp), occasion 1 (logit-scale variance)")  # Table 3A: IOV k_bp = 5.8% (RSE 123.7%)
    etaiov_fktr_2 ~ fix(0.4350^2)                                                  # same IOV magnitude on occasion 2

    # CL/F inter-occasion variability carried over from the PK layer.
    etaltlag ~ 0.67^2  # Table 2, AST column: IIV Tlag = 0.67 (reported as a standard deviation)
    etalka   ~ 0.87^2  # Table 2, AST column: IIV ka   = 0.87
    etalvc   ~ 0.04^2  # Table 2, AST column: IIV V/F  = 0.04
    etalcl   ~ 0.26^2  # Table 2, AST column: IIV CL/F = 0.26
    etaiov_cl_1 ~ 0.12^2
    label("Inter-occasion variability in CL/F, occasion 1 (log-scale variance)")  # Table 2, AST column: IOV CL/F = 0.12
    etaiov_cl_2 ~ fix(0.12^2)                                                     # same IOV magnitude on occasion 2

    # Residual error. Section 2.4.2: "A proportional residual error model was
    # used". Table 3A reports it as b = 0.29. The palbociclib concentrations
    # are not an endpoint of this fit -- the PK was fixed to individual post hoc
    # estimates -- so only the ANC carries an error model here; the plasma
    # residual error lives in modellib('Panetta_2024_palbociclib').
    propSd_ANC <- 0.29; label("Proportional residual error on ANC (fraction)")  # Table 3A: residual error b = 0.29 (RSE 6.0%)
  })

  model({
    # ---- 1. Occasion indicators ----
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_cl   <- oc1 * etaiov_cl_1   + oc2 * etaiov_cl_2
    iov_kin  <- oc1 * etaiov_kin_1  + oc2 * etaiov_kin_2
    iov_ic50 <- oc1 * etaiov_ic50_1 + oc2 * etaiov_ic50_2
    iov_fktr <- oc1 * etaiov_fktr_1 + oc2 * etaiov_fktr_2

    # ---- 2. Individual parameters ----
    tlag <- exp(ltlag + etaltlag)
    ka   <- exp(lka + etalka)
    vc   <- exp(lvc + etalvc)
    cl   <- exp(lcl + etalcl + iov_cl) * (AST / 25)^e_ast_cl

    rbase <- exp(lrbase + etalrbase)
    kin   <- exp(lkin + etalkin + iov_kin)
    ic50  <- exp(lic50 + etalic50 + iov_ic50)
    kout  <- exp(lkout)
    hill  <- exp(lhill)
    # k_bp = fk_bp * kin with fk_bp in (0, 1), so kin > k_bp always holds and
    # the feedback constant km_fb below is always positive.
    fktr <- expit(logitfktr + etalogitfktr + iov_fktr)
    ktr  <- fktr * kin

    # ---- 3. Palbociclib PK: one compartment, first-order absorption ----
    kel <- cl / vc
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central
    alag(depot) <- tlag
    # central / vc is mg/L because doses are mg/m^2 and vc is L/m^2;
    # 1 mg/L = 1000 ng/mL.
    Cc <- 1000 * central / vc

    # ---- 4. PK-PD link ----
    # Equation 1 is driven directly by Cc, the plasma palbociclib concentration
    # in ng/mL. See the IC50 entry in ini() for why no nM conversion is applied
    # despite Panetta 2024 labelling IC50 and C_palbociclib "nM".

    # ---- 5. Steady-state-derived constants (Equation 6) ----
    # KM = k_bp / (kin - k_bp) * Ncirc0. Substituting k_bp = fk_bp * kin makes
    # kin cancel: KM = fk_bp / (1 - fk_bp) * Ncirc0. At the typical values that
    # is 0.866667 / 0.133333 * 3.95 = 25.675 x 10^3/uL. Written in the
    # cancelled form because it is numerically safer and algebraically
    # identical.
    km_fb <- fktr / (1 - fktr) * rbase

    # ---- 6. Neutrophil dynamics: Equations 1-5 ----
    # Equation 1 as printed omits the trailing N_PBM factor on the bracketed
    # proliferation term; that it belongs there is forced by the paper's own
    # Equations 6 and 7. Setting dN_PBM/dt = 0 in the drug-free steady state
    # with the factor present gives kin KM / (KM + Ncirc0) = k_bp, which
    # rearranges to exactly Equation 6, and the resulting initial conditions
    # are exactly Equation 7. Without the factor neither identity holds.
    feedback <- km_fb / (km_fb + circ)
    inhib    <- ic50^hill / (ic50^hill + Cc^hill)
    d/dt(precursor1) <- kin * feedback * inhib * precursor1 - ktr * precursor1  # Equation 1, N_PBM
    d/dt(precursor2) <- ktr * (precursor1 - precursor2)                        # Equation 2, N_T1
    d/dt(precursor3) <- ktr * (precursor2 - precursor3)                        # Equation 3, N_T2
    d/dt(precursor4) <- ktr * (precursor3 - precursor4)                        # Equation 4, N_T3
    d/dt(circ)       <- ktr * precursor4 - kout * circ                         # Equation 5, N_circ

    # ---- 7. Drug-free steady-state initial conditions (Equation 7) ----
    precursor1(0) <- kout / ktr * rbase
    precursor2(0) <- kout / ktr * rbase
    precursor3(0) <- kout / ktr * rbase
    precursor4(0) <- kout / ktr * rbase
    circ(0)       <- rbase

    # ---- 8. Observation and error ----
    ANC <- circ
    ANC ~ prop(propSd_ANC)
  })
}
