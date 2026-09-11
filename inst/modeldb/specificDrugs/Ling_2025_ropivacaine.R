Ling_2025_ropivacaine <- function() {
  description <- paste(
    "Two-compartment population PK model for ropivacaine given as a single ultrasound-guided",
    "superficial serratus anterior plane block (SAPB) at 3 mg/kg in adults undergoing",
    "video-assisted thoracoscopic lobectomy (Ling 2025). Absorption from the fascial plane is a",
    "parallel mixed-order process: a fraction frel = 72.6 % of the dose enters a depot and is",
    "absorbed first-order at rate ka, while the complementary 27.4 % enters the central",
    "compartment directly as a zero-order input of duration D2 = 0.015 h beginning after a lag",
    "ALAG2 = 0.49 h. Disposition is parameterised on the rate constant k rather than on a",
    "clearance, with an apparent central volume Vc/F = 125 L, apparent inter-compartmental",
    "clearance Q/F = 14.7 L/h and apparent peripheral volume Vp/F = 197 L; the implied apparent",
    "clearance k * Vc/F = 7.48 L/h is a derived quantity the paper quotes but does not estimate.",
    "Two covariates were retained. The concentration of the injected ropivacaine solution acts on",
    "ka, which was estimated separately in each concentration stratum (32.0, 19.4 and 14.4 1/h for",
    "the 0.25 %, 0.5 % and 0.75 % w/v solutions), and platelet count acts on Vc/F. IMPORTANT: the",
    "paper prints the platelet coefficient (-0.438) but never prints the covariate equation or its",
    "centring value; a median-normalised power form referenced to 200 x 10^9/L is used here. See",
    "the vignette Errata for that and for the three other reading decisions this extraction had to",
    "make."
  )
  reference <- paste(
    "Ling J, Xu C, Tang L, Qiu L, Hu N. Comparison of the pharmacokinetic variations of different",
    "concentrations of ropivacaine used for serratus anterior plane block in patients undergoing",
    "thoracoscopic lobectomy: a population pharmacokinetics analysis.",
    "Front Pharmacol. 2025;16:1540606. doi:10.3389/fphar.2025.1540606 (Table 1, Table 2, Table 3,",
    "Table 4 and the Population pharmacokinetic modeling / Simulation sections)."
  )
  vignette <- "Ling_2025_ropivacaine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    PLT = list(
      description        = "Preoperative platelet count from the routine complete blood count",
      units              = "10^9 cells/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Ling 2025 Table 1 reports the platelet count as median (range) per concentration group:",
        "182 (43-341) in the 0.25 % group, 213 (133-344) in the 0.5 % group and 187 (132-240) in",
        "the 0.75 % group (p = 0.208 across groups). The n-weighted mean of those three medians is",
        "194 x 10^9/L. Time-fixed: a single preoperative value per patient. Retained in the final",
        "model on the apparent central volume Vc/F with theta = -0.438 (RSE 29.7 %, Table 3), so a",
        "higher platelet count gives a SMALLER apparent central volume. The paper never prints the",
        "covariate equation, the transformation or the centring constant, so this file uses the",
        "median-normalised power form (PLT / 200)^-0.438 that is standard for a continuous",
        "covariate on a volume and that the covariate register records for PLT",
        "(Stitt_2026_tranexamicAcid.R uses (PLT/196)^0.468 on clearance). 200 x 10^9/L is the",
        "rounded clinical standard and sits within 3 % of the cohort's own 194 x 10^9/L centre;",
        "under a power form the reference choice is nearly inconsequential, since moving it from",
        "195 to 200 rescales Vc/F by only 1.1 %. See the vignette Errata for the alternative",
        "readings that were considered and rejected."
      ),
      source_name        = "Platelet count (x10^9 L^-1)"
    ),
    FORM_ROPI_SOLN05 = list(
      description        = paste(
        "Ropivacaine injectate-concentration indicator: 1 = the block was performed with the 0.5 %",
        "w/v (5 mg/mL) ropivacaine solution, 0 = a different solution strength was used"
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the 0.25 % w/v solution when FORM_ROPI_SOLN075 is also 0)",
      notes              = paste(
        "Time-fixed per patient: each patient received one single-shot block from one solution",
        "strength. Ling 2025 randomised patients to 0.25 %, 0.5 % and 0.75 % w/v ropivacaine at a",
        "constant 3 mg/kg dose, so the injected mass is the same in every arm and only the injected",
        "volume differs. The concentration of the solution was the only covariate retained on the",
        "first-order absorption rate constant, which the authors estimated separately in each",
        "stratum (Table 3: ka = 32.0, 19.4 and 14.4 1/h for 0.25 %, 0.5 % and 0.75 %), so this",
        "indicator selects the 0.5 % stratum estimate rather than scaling a reference value. Pair",
        "with FORM_ROPI_SOLN075; both zero selects the 0.25 % reference stratum. The two",
        "indicators are mutually exclusive."
      ),
      source_name        = "The concentration of ropivacaine (0.5 %)"
    ),
    FORM_ROPI_SOLN075 = list(
      description        = paste(
        "Ropivacaine injectate-concentration indicator: 1 = the block was performed with the",
        "0.75 % w/v (7.5 mg/mL) ropivacaine solution, 0 = a different solution strength was used"
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the 0.25 % w/v solution when FORM_ROPI_SOLN05 is also 0)",
      notes              = paste(
        "Time-fixed per patient; see FORM_ROPI_SOLN05. Selects the 0.75 % stratum estimate",
        "ka = 14.4 1/h (Ling 2025 Table 3). Ling 2025 also enrolled two patients who received the",
        "0.375 % solution; those two were held out for external validation and are NOT part of the",
        "fitted model, so no indicator is defined for that strength. The paper interpolates their",
        "ka post hoc from a log-linear regression through the three fitted values (Results, Model",
        "validation) -- that regression is a validation device, not part of the final model, and is",
        "not encoded here."
      ),
      source_name        = "The concentration of ropivacaine (0.75 %)"
    )
  )

  # Screened in the stepwise covariate search (Ling 2025, Population
  # pharmacokinetic modeling) but NOT retained in the final model, so no
  # coefficient is published for any of them. Recorded here for provenance
  # only; none is referenced in model().
  covariatesDataExcluded <- list(
    AGE   = list(description = "Age", units = "years", type = "continuous",
                 notes = "Table 1 medians 60.5 / 58 / 59 y (range 31-75); screened, not retained."),
    SEXF  = list(description = "Female sex indicator", units = "(binary)", type = "binary",
                 notes = "Screened, not retained. Ling 2025 Table 1 does not report the sex distribution of the cohort at all."),
    WT    = list(description = "Body weight", units = "kg", type = "continuous",
                 notes = "Table 1 medians 57.3 / 60.5 / 60 kg (range 50-81). Screened, not retained; the Discussion attributes this to the narrow body-size range. Body weight still determines the administered dose (3 mg/kg)."),
    WBC   = list(description = "White blood cell count", units = "10^9 cells/L", type = "continuous",
                 notes = "Table 1 medians 5.26 / 6.97 / 6.31; screened, not retained."),
    RBC   = list(description = "Red blood cell count", units = "10^12 cells/L", type = "continuous",
                 notes = "Table 1 medians 4.09 / 4.33 / 4.17, printed with the units x10^9 L^-1 which is a typo for x10^12 L^-1; screened, not retained."),
    HGB   = list(description = "Hemoglobin", units = "g/L", type = "continuous",
                 notes = "Named in the screening list; no summary statistics printed. Not retained."),
    HCT   = list(description = "Hematocrit", units = "(fraction)", type = "continuous",
                 notes = "Named in the screening list; no summary statistics printed. Not retained."),
    ALB   = list(description = "Serum albumin", units = "g/L", type = "continuous",
                 notes = "Named in the screening list; no summary statistics printed. Not retained. The Discussion notes that protein binding is reported elsewhere to affect ropivacaine PK."),
    TBIL  = list(description = "Total bilirubin", units = "umol/L", type = "continuous",
                 notes = "Named in the screening list; no summary statistics printed. Not retained."),
    ALP   = list(description = "Alkaline phosphatase", units = "U/L", type = "continuous",
                 notes = "Named in the screening list; no summary statistics printed. Not retained."),
    ALT   = list(description = "Alanine aminotransferase", units = "U/L", type = "continuous",
                 notes = "Table 1 medians 15.3 / 16.6 / 20.1 U/L; screened, not retained."),
    AST   = list(description = "Aspartate aminotransferase", units = "U/L", type = "continuous",
                 notes = "Table 1 medians 20.6 / 20.9 / 24.4 U/L; screened, not retained."),
    CREAT = list(description = "Serum creatinine", units = "umol/L", type = "continuous",
                 notes = "Table 1 medians 65 / 57.5 / 65, printed with the units mmol/L which is a typo for umol/L (65 mmol/L is not a survivable creatinine). Screened, not retained."),
    TBA   = list(description = "Serum total bile acid", units = "umol/L", type = "continuous",
                 notes = "Named in the screening list; no summary statistics printed. Not retained."),
    EGFR  = list(description = "Glomerular filtration rate", units = "mL/min/1.73m^2", type = "continuous",
                 notes = "Named in the screening list; the estimating equation is not stated and no summary statistics are printed. Not retained."),
    URATE = list(description = "Serum uric acid", units = "umol/L", type = "continuous",
                 notes = "Named in the screening list; no summary statistics printed. Not retained."),
    CONMED_PROPOFOL  = list(description = "Concomitant propofol", units = "(binary)", type = "binary",
                            notes = "Table 1: used in 25.0 / 71.4 / 66.7 % of the three groups (p = 0.035). Screened, not retained."),
    CONMED_LIDOCAINE = list(description = "Concomitant lidocaine", units = "(binary)", type = "binary",
                            notes = "Table 1: used in 58.3 / 57.1 / 53.3 % of the three groups (p = 0.962). Screened, not retained."),
    CONMED_DYCLONINE = list(description = "Concomitant dyclonine mucilage", units = "(binary)", type = "binary",
                            notes = "Table 1: used in 58.3 / 42.8 / 93.3 % of the three groups (p = 0.013). Screened, not retained. Only comedications used by more than 5 % of patients were tested.")
  )

  compartmentData <- list(
    # Ling 2025 sampled ARTERIAL plasma throughout (Methods, Patients); the
    # register's specimen vocabulary does not distinguish arterial from venous
    # plasma, so "plasma" is used and the arterial sampling site is recorded
    # here and in population$notes. The depot is the serratus anterior fascial
    # plane into which the block is injected.
    depot       = list(analyte = "ropivacaine", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "ropivacaine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ropivacaine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 41L,
    n_studies      = 1L,
    age_range      = "31-75 years (group medians 60.5, 58 and 59 years)",
    weight_range   = "50-81 kg (group medians 57.3, 60.5 and 60 kg)",
    sex_female_pct = NA_real_,
    disease_state  = "Adults undergoing primary elective video-assisted thoracoscopic lung resection (lobectomy), ASA physical status I-III, without chronic pain, cognitive dysfunction or overt organ dysfunction",
    dose_range     = "Single 3 mg/kg ropivacaine superficial serratus anterior plane block (150-243 mg observed) as a 0.25 %, 0.5 % or 0.75 % w/v solution",
    regions        = "China (The First People's Hospital of Changzhou / The Third Affiliated Hospital of Soochow University)",
    notes          = paste(
      "43 patients were enrolled between April and December 2023 and randomised by random-number",
      "table to 0.25 % (n = 12), 0.375 % (n = 2), 0.5 % (n = 14) and 0.75 % (n = 15) ropivacaine.",
      "388 arterial plasma concentrations from the 41 patients in the 0.25 / 0.5 / 0.75 % arms",
      "built the model; the 18 concentrations from the two 0.375 % patients were held out for",
      "external validation (mean prediction error -0.18 %, mean absolute prediction error 9.42 %,",
      "Table 4). Arterial blood was drawn at 1, 15, 30 and 45 min and 1, 2, 4, 8, 12 and 24 h after",
      "the block. Ropivacaine was assayed by LC-MS/MS with a lower limit of quantification of",
      "4 ng/mL. Fitted in NONMEM by FOCE with interaction; evaluated by goodness-of-fit plots, a",
      "1000-replicate VPC and NPDE. Ling 2025 does not report the sex distribution of the cohort",
      "even though sex was screened as a covariate."
    )
  )

  ini({
    # ====================================================================
    # Structural parameters -- Ling 2025 Table 3, "Estimate (RSE%)".
    #
    # Disposition is parameterised on the elimination RATE CONSTANT k with
    # an explicit apparent volume Vc/F, exactly as the paper reports it,
    # and is deliberately NOT reparameterised to CL/F + Vc/F: separate
    # random effects sit on k and on Vc/F (omega_k = 68.3 %, omega_Vc/F =
    # 24.0 %), so the clearance form would require a correlated eta block
    # the authors never fitted. The paper's quoted CL/F = 7.475 L/h is
    # simply k * Vc/F = 0.0598 * 125 and is not an estimated parameter.
    # ====================================================================

    # The first-order absorption rate constant is estimated SEPARATELY in
    # each injectate-concentration stratum (Table 3 splits the "ka" row
    # into three sub-rows), so all three carry an explicit stratum suffix
    # and none keeps the bare canonical name. A single shared IIV
    # (etalka) applies across the strata.
    lka_soln025 <- log(32.0);  label("First-order absorption rate constant ka with the 0.25 % w/v solution (1/h)")  # Ling 2025 Table 3: ka (0.25% ropivacaine) = 32.0 1/h, RSE 22.1%
    lka_soln05  <- log(19.4);  label("First-order absorption rate constant ka with the 0.5 % w/v solution (1/h)")   # Ling 2025 Table 3: ka (0.50% ropivacaine) = 19.4 1/h, RSE 19.9%
    lka_soln075 <- log(14.4);  label("First-order absorption rate constant ka with the 0.75 % w/v solution (1/h)")  # Ling 2025 Table 3: ka (0.75% ropivacaine) = 14.4 1/h, RSE 18.3%

    lkel <- log(0.0598); label("First-order elimination rate constant k from the central compartment (1/h)")        # Ling 2025 Table 3: k = 0.0598 1/h, RSE 13.4%
    lvc  <- log(125);    label("Apparent central volume of distribution Vc/F at PLT = 200 x 10^9/L (L)")            # Ling 2025 Table 3: Vc/F = 125 L, RSE 4.8%
    lq   <- log(14.7);   label("Apparent inter-compartmental clearance Q/F (L/h)")                                  # Ling 2025 Table 3: Q/F = 14.7 L/h, RSE 31.5%
    lvp  <- log(197);    label("Apparent peripheral volume of distribution Vp/F (L)")                               # Ling 2025 Table 3: Vp/F = 197 L, RSE 10.5%

    # ====================================================================
    # Parallel mixed-order absorption. F1 = 72.6 % of the dose takes the
    # first-order (depot) route; the complementary 27.4 % is delivered
    # into central as a zero-order input of duration D2 = 0.015 h starting
    # after the lag ALAG2 = 0.49 h. Both routes are apparent (every
    # disposition parameter is divided by the true bioavailability F), so
    # frel and 1 - frel describe only how the absorbed dose is split
    # between the two routes, not how much of it is absorbed -- which is
    # why the canonical is logitfrel (a release-process split) rather than
    # logitfdepot (a bioavailability).
    #
    # frel is held on the LOGIT scale because it carries IIV: an
    # exponential (log-scale) IIV of 21.1 % on a typical value of 0.726
    # puts 6.5 % of simulated subjects above 1, which would make
    # f(central) = 1 - frel negative. The point estimate is preserved
    # exactly, since expit(log(0.726 / 0.274)) = 0.726.
    # ====================================================================
    logitfrel <- log(0.726 / (1 - 0.726)); label("Logit of the fraction F1 of the dose absorbed by the first-order (depot) route (logit units)")  # Ling 2025 Table 3: F1 = 72.6%, RSE 8.0%; the Abstract gives the complement, "the proportion of zero-order absorption was 27.4%"; logit(0.726) = 0.9738
    ld2       <- log(0.015);                label("Duration D2 of the zero-order absorption input into the central compartment (h)")               # Ling 2025 Table 3: D2 = 0.015 h, RSE 12.5%
    ltlag     <- log(0.49);                 label("Lag time ALAG2 before the zero-order absorption input begins (h)")                              # Ling 2025 Table 3: ALAG2 = 0.49 h, RSE 0.4%

    # ====================================================================
    # Covariate effect. Ling 2025 Table 3 prints theta(PLT-Vc/F) = -0.438
    # (RSE 29.7%) but the paper NEVER prints the covariate equation, the
    # transformation, or the centring constant. Encoded here as the
    # median-normalised power form
    #     Vc/F = 125 * (PLT / 200)^(-0.438)
    # which is the standard shape for a continuous covariate on a volume,
    # is the shape the covariate register records for PLT, is strictly
    # positive at every platelet count, and is identical to first order
    # around the reference to a median-normalised linear form. Reading
    # -0.438 as a bare per-unit linear coefficient is arithmetically
    # impossible (it drives Vc/F negative 2.3 x 10^9/L above the
    # reference). See the vignette Errata.
    # ====================================================================
    e_plt_vc <- -0.438; label("Power exponent on (PLT / 200) for the apparent central volume Vc/F (unitless)")  # Ling 2025 Table 3: theta(PLT-Vc/F) = -0.438, RSE 29.7%

    # ====================================================================
    # Inter-individual variability. Ling 2025 Methods define the IIV as
    # exponential, P_j = P_hat * exp(eta_j), with "eta_j ... a random
    # variable distributed with a mean of zero and variance of omega^2".
    # Table 3 lists the symbol OMEGA itself (not omega^2, and not a %CV)
    # for each parameter, in percent -- so each tabulated value divided by
    # 100 IS the log-scale standard deviation, and the variance below is
    # its square. There is no CV-to-variance conversion to do.
    #
    # The IIV on the absorption fraction is carried on the logit scale
    # (see logitfrel above). Converting the reported log-scale SD to the
    # logit scale by the delta method preserves SD(F) exactly to first
    # order: SD(logit F) = SD(ln F) / (1 - F) = 0.211 / 0.274 = 0.770,
    # and 0.726 * 0.274 * 0.770 = 0.1532 = 0.726 * 0.211.
    # ====================================================================
    etalogitfrel ~ 0.770^2   # Ling 2025 Table 3: omega_F1   = 21.1% (RSE 33.4%), rescaled from the log scale to the logit scale by the delta method as described above
    etalka       ~ 0.582^2   # Ling 2025 Table 3: omega_ka   = 58.2% (RSE 36.0%); one shared eta across the three concentration strata
    etalkel      ~ 0.683^2   # Ling 2025 Table 3: omega_k    = 68.3% (RSE 23.3%)
    etalvc       ~ 0.240^2   # Ling 2025 Table 3: omega_Vc/F = 24.0% (RSE 27.2%)
    etalvp       ~ 1.054^2   # Ling 2025 Table 3: omega_Vp/F = 105.4% (RSE 29.2%)

    # ====================================================================
    # Residual unexplained variability. Ling 2025 Methods give the
    # combined form C_ij = Chat_ij * (1 + eps_1) + eps_2 on the observed
    # ng/mL scale, which is nlmixr2's prop() + add().
    # ====================================================================
    propSd <- 0.145; label("Proportional residual error on the plasma ropivacaine concentration (fraction)")  # Ling 2025 Table 3: delta_prop = 14.5%, RSE 21.7%
    addSd  <- 80.1;  label("Additive residual error on the plasma ropivacaine concentration (ng/mL)")         # Ling 2025 Table 3: delta_add = 80.1 ng/mL, RSE 30.8%
  })

  model({
    # ====================================================================
    # Individual parameters.
    #
    # The typical first-order absorption rate constant is selected by the
    # two injectate-concentration indicators; both zero gives the 0.25 %
    # reference stratum. The fixed effects and the shared eta are
    # collected on one line so the eta stays mu-referenced.
    # ====================================================================
    lka_ind <- lka_soln025 +
      (lka_soln05  - lka_soln025) * FORM_ROPI_SOLN05 +
      (lka_soln075 - lka_soln025) * FORM_ROPI_SOLN075 +
      etalka
    ka <- exp(lka_ind)

    kel  <- exp(lkel + etalkel)
    vc   <- exp(lvc  + etalvc) * (PLT / 200)^e_plt_vc
    q    <- exp(lq)
    vp   <- exp(lvp  + etalvp)
    d2   <- exp(ld2)
    tlag <- exp(ltlag)

    logitfrel_ind <- logitfrel + etalogitfrel
    frel <- expit(logitfrel_ind)

    # ====================================================================
    # Two-compartment disposition with parallel mixed-order absorption
    # from the serratus anterior fascial plane:
    #
    #   depot   : fraction frel = 0.726 of the dose, absorbed first-order
    #             at rate ka, whose value depends on the concentration of
    #             the injected solution.
    #   central : fraction 1 - frel = 0.274, delivered as a zero-order
    #             input of duration D2 = 0.015 h starting after the lag
    #             ALAG2 = 0.49 h. Ling 2025 describes the fascial plane as
    #             acting like a drug depot pump from which the local
    #             anaesthetic diffuses into the surrounding tissue.
    #
    # A user therefore encodes each block as TWO dose records at the same
    # time carrying the same full dose amount: one with cmt = "depot" and
    # one with cmt = "central" and rate = -2 (modelled duration). The f()
    # multipliers split the dose between the two routes; dur(central)
    # imposes the zero-order duration and alag(central) the lag.
    # ====================================================================
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central -
                          (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1

    f(depot)      <- frel
    f(central)    <- 1 - frel
    dur(central)  <- d2
    alag(central) <- tlag

    # Arterial plasma ropivacaine concentration. Doses are in mg and vc is
    # in L, so central / vc is mg/L = ug/mL; the factor 1000 converts to
    # the ng/mL units in which Ling 2025 reports every concentration,
    # Cmax, AUC and toxicity threshold.
    Cc <- (central / vc) * 1000

    Cc ~ add(addSd) + prop(propSd)
  })
}
