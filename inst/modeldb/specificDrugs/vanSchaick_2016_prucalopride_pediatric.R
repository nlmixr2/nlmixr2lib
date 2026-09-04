vanSchaick_2016_prucalopride_pediatric <- function() {
  description <- "Two-compartment oral population PK model for prucalopride in children with functional constipation (van Schaick 2016), jointly fit to a richly sampled single-dose phase 1 study (PRU-USA-12) and a sparsely sampled multiple-dose phase 3 study (SPD555-303). Absorption is a dual sequential first-order process: a slow rate applies before a fixed cut-off time after each dose and a fast rate applies after it. CL and Q are allometrically scaled to body weight (fixed exponent 0.75) with a fixed Rhodin 2009 postmenstrual-age renal-maturation function on CL; Vc and Vp are allometrically scaled (fixed exponent 1.0). Typical CL, CL interindividual variability and residual error are study-specific."
  reference   <- "van Schaick E, Benninga MA, Levine A, Magnusson M, Troy S. Development of a population pharmacokinetic model of prucalopride in children with functional constipation. Pharmacol Res Perspect. 2016;4(4):e00236. doi:10.1002/prp2.236"
  vignette    <- "vanSchaick_2016_prucalopride_pediatric"
  units       <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Non-canonical IIV and residual-error names declared so
  # checkModelConventions() can distinguish deliberate study-specific
  # parameterisation from naming drift. Both CL (typical value and IIV)
  # and the residual error are estimated separately per source study
  # (van Schaick 2016 Table 3), and the two absorption rate constants
  # carry separate fixed IIVs inherited from the adult model.
  paper_specific_etas         <- c("etalcl_pru", "etalcl_spd", "etalka_early", "etalka_late")
  paper_specific_residual_sds <- c("expSdPRU", "expSdSPD")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot       = list(analyte = "prucalopride", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "prucalopride", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "prucalopride", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used for standard allometric scaling on CL and Q (fixed exponent 3/4, van Schaick 2016 equation 2) and on Vc and Vp (fixed exponent 1, equation 3), reference 70 kg. Cohort medians 27.9 kg (PRU-USA-12) and 24.0 kg (SPD555-303); pooled range 11.0-110.0 kg (Table 2).",
      source_name        = "WT"
    ),
    PAGE = list(
      description        = "Postmenstrual age",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "van Schaick 2016 defines postmenstrual age as calculated age in weeks at start of treatment plus 40 weeks gestational age (Methods 'Structural model components', Table 2 footnote). Drives the fixed Rhodin 2009 renal-maturation Hill function on CL (equation 4), which is defined in WEEKS of PMA; the canonical PAGE unit is months, so model() converts as PMA_weeks = PAGE * 4.35 (same handling as LlanosPaez_2020_gentamicin.R). Cohort medians 9.3 years (PRU-USA-12) and 8.6 years (SPD555-303); pooled range approximately 2.4-18.8 years PMA (Table 2).",
      source_name        = "PMA"
    ),
    STUDY_SPD555303 = list(
      description        = "Source-study cohort indicator (1 = SPD555-303 phase 3 multiple-dose study, 0 = PRU-USA-12 phase 1 single-dose study)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (PRU-USA-12 phase 1 single-dose study)",
      notes              = "Switches between the study-specific typical clearances (CL_PRU-USA-12 = 22.9 vs CL_SPD555-303 = 20.1 L/h/70 kg), their study-specific IIVs (omega^2 0.0151 vs 0.1191) and the study-specific log-scale residual SDs (0.142 vs 0.35); see van Schaick 2016 Table 3. Set STUDY_SPD555303 = 1 to reproduce the paper's own dosing simulations, which the paper states were 'performed with clearance based on SPD555-303' (Results, 'Simulated plasma concentration-time profiles').",
      source_name        = "derived"
    )
  )

  # Screened during model development but NOT retained in the final model:
  # van Schaick 2016 equation 5 contains only body weight and the PMA
  # maturation term. Documented here to preserve the provenance of the
  # covariate screen without carrying convention warnings.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Calculated (chronological) age",
      units       = "years",
      type        = "continuous",
      notes       = "Listed in Methods 'Structural model components' among the covariates 'considered in the analysis because they were expected to impact the pharmacokinetics of prucalopride', and plotted against the random effects in Figure 4A / Figure S2. Not retained: the paper concludes 'the allometric relationships between CL, V2, and V3 and body weight adequately accounted for body size and age in this pediatric population'. Age enters the final model only indirectly, through PAGE in the maturation term."
    ),
    CRCL = list(
      description = "Creatinine clearance (Schwartz 1976 as adjusted for body weight by Rowland and Tozer 2010)",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened via van Schaick 2016 equation 1, CrCL[mL/min] = 42.5 * Height[cm] / SerumCreatinine[umol/L] * (WT[kg]/70)^0.7 -- raw, NOT BSA-normalized. Cohort means 82.8 (PRU-USA-12) and 71.9 mL/min (SPD555-303), pooled range 27.9-180.0 (Table 2). Plotted against the random effects in Figure 4C / Figure S2 but not retained in the final model; no point estimate for a CrCL effect is reported anywhere in the paper."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 175L,
    n_studies      = 2L,
    age_range      = "4.0-12.0 years (PRU-USA-12); 1.7-18.0 years (SPD555-303)",
    age_median     = "8.5 years (PRU-USA-12); 7.9 years (SPD555-303)",
    weight_range   = "15.0-61.0 kg (PRU-USA-12); 11.0-110.0 kg (SPD555-303)",
    weight_median  = "27.9 kg (PRU-USA-12); 24.0 kg (SPD555-303)",
    sex_female_pct = 52.6,
    race_ethnicity = "Not reported",
    disease_state  = "Functional constipation",
    dose_range     = "PRU-USA-12: single oral dose of prucalopride 0.03 mg/kg (0.02 mg/kg in one patient) as a 0.2 mg/mL oral solution. SPD555-303: prucalopride 0.04 mg/kg once daily for body weight <= 50 kg (maximum 2 mg), titratable to 0.06 mg/kg for insufficient response or down to 0.02 mg/kg for tolerability, as a 0.4 mg/mL oral solution or 2 mg tablet.",
    regions        = "Not reported (PRU-USA-12 conducted in the USA; SPD555-303 was a multinational phase 3 trial, NCT01330381)",
    postmenstrual_age_range = "Approximately 2.4-18.8 years PMA across both studies (calculated age plus 40 weeks gestational age; Table 2)",
    samples_plasma = "PRU-USA-12: 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 18, 24, 48 and 72 h post-dose, mean 12.6 samples per patient. SPD555-303: one sample 1-3 h after the first dose plus two steady-state samples 14-26 h post-dose at weeks 8 and 24.",
    notes          = "Demographics from van Schaick 2016 Table 2. n_subjects is the pooled PK dataset (38 in PRU-USA-12 + 137 in SPD555-303) per Table 1 'Number of patients in PK dataset', Table 2 and Table 4; note that the Results section 'Available data' instead describes the final analysis dataset as '481 records from 38 of 38 patients in PRU-USA-12, and 244 records from 106 of 107 randomized patients in SPD555-303', an internal inconsistency in the paper that is not resolved by any other reported figure. sex_female_pct is the pooled value (13 of 38 = 34.2% in PRU-USA-12 and 79 of 137 = 57.7% in SPD555-303; 92 of 175 = 52.6% pooled). Assay was radioimmunoassay with LLOQ 0.1 ng/mL in PRU-USA-12 and LC-MS/MS with LLOQ 0.2 ng/mL in SPD555-303; below-LLOQ observations were excluded rather than modelled. Estimation by FOCE with interaction. The absorption parameters (Ka1, Ka2, MTIME, their IIVs) and the relative bioavailability F1 = 0.858 were fixed to previously estimated adult values ('data on file'), not estimated here."
  )

  ini({
    # -------- Structural disposition (van Schaick 2016 Table 3, upper
    # block; all population estimates are scaled to a body weight of
    # 70 kg per the Table 3 footnote). --------
    lcl_pru <- log(22.9); label("Apparent clearance CL/F standardised to 70 kg, PRU-USA-12 (L/h/70 kg)")     # Table 3: CL_PRU-USA-12 = 22.9 L/h (RSE 2.4%, 95% CI 21.9-24)
    lcl_spd <- log(20.1); label("Apparent clearance CL/F standardised to 70 kg, SPD555-303 (L/h/70 kg)")     # Table 3: CL_SPD555-303 = 20.1 L/h (RSE 3.6%, 95% CI 18.6-21.5)
    lvc     <- log(446);  label("Apparent central volume Vc/F standardised to 70 kg (L/70 kg)")              # Table 3: V2 = 446 L (RSE 3.3%, 95% CI 417-475)
    lq      <- log(16.9); label("Apparent intercompartmental clearance Q/F standardised to 70 kg (L/h/70 kg)") # Table 3: Q = 16.9 L/h (RSE 15%, 95% CI 11.7-22)
    lvp     <- log(248);  label("Apparent peripheral volume Vp/F standardised to 70 kg (L/70 kg)")           # Table 3: V3 = 248 L (RSE 7.9%, 95% CI 210-286)

    # -------- Dual sequential first-order absorption. van Schaick 2016
    # Methods 'Structural model components': "one absorption rate (Ka1)
    # was applied before an estimated cut-off time (MTIME) and a second
    # absorption rate (Ka2) applied after this cut-off time", and the
    # Table 3 legend defines "Ka1, absorption rate at time less than
    # MTIME; Ka2: absorption rate at time greater than MTIME". All three
    # values, and the relative bioavailability, were "fixed to the value
    # previously estimated in adults" (Methods; F1 = 85.8%, data on
    # file), so all four carry fixed(). --------
    lka_early <- fixed(log(0.792)); label("First-order absorption rate constant before the cut-off time (1/h)") # Table 3: Ka1 = 0.792 1/h, fixed
    lka_late  <- fixed(log(3.87));  label("First-order absorption rate constant after the cut-off time (1/h)")  # Table 3: Ka2 = 3.87 1/h, fixed
    tkacut    <- fixed(0.734);      label("Time after dose at which absorption switches from the early to the late rate (h)") # Table 3: MTIME = 0.734 h, fixed
    lfdepot   <- fixed(log(0.858)); label("Relative bioavailability of the oral dose (fraction)")               # Table 3: F1 = 0.858, fixed (Methods: 85.8%, data on file)

    # -------- Allometric exponents, fixed at the standard 3/4 and 1
    # power model. van Schaick 2016 Methods 'Allometric scaling and
    # age-related maturation': "Q, CL, V2, and V3 parameters were scaled
    # across the ages from 6 months to 18 years by body size, using
    # standard allometric scaling equations 2 and 3" (Anderson and
    # Holford 2008), reference 70 kg. --------
    e_wt_cl_q  <- fixed(0.75); label("Shared allometric exponent on CL and Q (unitless)")   # Equation 2: (WT/70)^(3/4) on CL and Q
    e_wt_vc_vp <- fixed(1.00); label("Shared allometric exponent on Vc and Vp (unitless)")  # Equation 3: (WT/70)^1 on V2 and V3

    # -------- Renal-maturation Hill function on CL, fixed at the Rhodin
    # 2009 GFR-maturation parameters. van Schaick 2016 equation 4:
    # MaturationGFR = PMA^3.4 / (47.7^3.4 + PMA^3.4), PMA in weeks.
    # Same fixed pair as Padari_2018_penicillin_G.R. Self-check reported
    # in the paper: "the fractional GFR is expected to be over 97.5% for
    # a 24-month-old child" -- at PMA = 24 months postnatal + 40 weeks
    # gestational = 144.3 weeks, equation 4 gives 97.7%. --------
    tmat50   <- fixed(47.7); label("Postmenstrual age at 50% renal maturation (weeks)")   # Equation 4 (Rhodin et al. 2009)
    hill_mat <- fixed(3.4);  label("Hill coefficient for renal maturation (unitless)")    # Equation 4 (Rhodin et al. 2009)

    # -------- Interindividual variability (van Schaick 2016 Table 3,
    # lower block, "IIV estimate (% CV)" column). Methods: IIV "was
    # described using an exponential error model", so these are
    # variances on the log scale. The four ESTIMATED rows are variances:
    # sqrt(0.0151) = 12.3%, sqrt(0.1191) = 34.5%, sqrt(0.0202) = 14.2%
    # and sqrt(0.158) = 39.8% reproduce the printed 12%, 35%, 14% and
    # 40% parentheticals. (The two residual-error rows in the same
    # column are on a DIFFERENT scale -- see the residual-error block
    # below and the vignette Errata.) The Ka1 and Ka2 IIVs were fixed to
    # the adult values along with the rate constants themselves.
    # eta-shrinkage was 6% (CL PRU-USA-12) and 24% (CL SPD555-303) but
    # large at 58% (V2) and 70% (V3), which the paper attributes to the
    # sparse SPD555-303 sampling. --------
    etalcl_pru   ~ 0.0151        # Table 3: IIV CL_PRU-USA-12 = 0.0151 (12%), RSE 28%, eta-shrinkage 6%
    etalcl_spd   ~ 0.1191        # Table 3: IIV CL_SPD555-303 = 0.1191 (35%), RSE 33%, eta-shrinkage 24%
    etalvc       ~ 0.0202        # Table 3: IIV V2 = 0.0202 (14%), RSE 47%, eta-shrinkage 58%
    etalvp       ~ 0.158         # Table 3: IIV V3 = 0.158 (40%), RSE 68%, eta-shrinkage 70%
    etalka_early ~ fixed(0.794)  # Table 3: IIV Ka1 = 0.794, from the adult model
    etalka_late  ~ fixed(0.507)  # Table 3: IIV Ka2 = 0.507, from the adult model

    # -------- Residual error. van Schaick 2016 Methods: "The residual
    # variability was explained with an additive error on
    # log-transformed data", which is the log-normal (exponential)
    # residual family in nlmixr2, so the parameter is a log-scale SD.
    # Unlike the IIV rows above, these two Table 3 rows are SDs and not
    # variances: 0.142 and 0.35 reproduce the printed 14% and 35%
    # parentheticals directly as value x 100, and the Discussion states
    # independently that "the residual error (e) in SPD555-303 was much
    # larger than that in PRU-USA-12 (35% in SPD555-303 vs. 14% in
    # PRU-USA-12)". A variance reading would give sqrt(0.142) = 37.7%
    # and sqrt(0.35) = 59.2%, contradicting both the table's own
    # parentheticals and the prose. --------
    expSdPRU <- 0.142; label("Log-scale residual SD, PRU-USA-12 (unitless)")  # Table 3: Residual error PRU-USA-12 = 0.142 (14%), RSE 9.8%
    expSdSPD <- 0.35;  label("Log-scale residual SD, SPD555-303 (unitless)")  # Table 3: Residual error SPD555-303 = 0.35 (35%), RSE 16%
  })

  model({
    # 1. Derived covariate terms.
    #
    # Renal-function maturation (van Schaick 2016 equation 4, Rhodin
    # 2009). Equation 4 is defined on postmenstrual age in WEEKS; the
    # canonical PAGE column carries months, so convert here.
    pma_wks <- PAGE * 4.35
    fmat    <- pma_wks^hill_mat / (tmat50^hill_mat + pma_wks^hill_mat)

    # Study-specific typical clearance and clearance IIV (Table 3). Both
    # etas are declared for every subject but only the one selected by
    # STUDY_SPD555303 reaches cl, matching the paper's joint fit in which
    # each subject belongs to exactly one study.
    cl_pru <- exp(lcl_pru + etalcl_pru)
    cl_spd <- exp(lcl_spd + etalcl_spd)
    cl_tv  <- cl_spd * STUDY_SPD555303 + cl_pru * (1 - STUDY_SPD555303)

    # 2. Individual PK parameters. Allometric size scaling to the 70 kg
    # reference (equations 2 and 3) with the PMA maturation factor on CL
    # only (equation 5: CL_i = CL_TV * (WT_i/70)^(3/4) * MaturationGFR).
    cl <- cl_tv             * (WT / 70)^e_wt_cl_q  * fmat
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp

    # Dual sequential first-order absorption. The early rate applies
    # before tkacut hours after the dose and the late rate after it;
    # tad(depot) resets at every dose, so the sequence repeats for each
    # administration (the reading confirmed by the operator, recorded in
    # the vignette Errata -- the paper does not state whether MTIME is
    # measured from the dose or from the start of the record).
    ka <- exp(lka_early + etalka_early)
    if (tad(depot) >= tkacut) ka <- exp(lka_late + etalka_late)

    # 3. Micro-constants for the two-compartment disposition.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 4. ODE system: first-order oral absorption into a two-compartment
    # linear disposition model.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central -
                          k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. Relative bioavailability of the oral dose (Table 3: F1 = 0.858,
    # fixed from the adult model).
    f(depot) <- exp(lfdepot)

    # 6. Observation. Dose in mg and volumes in L give mg/L, i.e.
    # ug/mL; multiply by 1000 to report ng/mL as in the paper. Unit
    # check against Table 4 (SPD555-303 post hoc means): AUC =
    # F1 * Dose / CL = 0.858 * (0.04 mg/kg * 32.4 kg) / 11.3 L/h =
    # 0.0984 mg*h/L = 98.4 ng*h/mL versus the published 100.3, and
    # Css = AUC / 24 h = 4.18 versus the published 4.18 ng/mL.
    Cc <- (central / vc) * 1000

    # Study-specific log-scale residual SD (Table 3), switched on the
    # same study indicator as clearance.
    expSd <- expSdSPD * STUDY_SPD555303 + expSdPRU * (1 - STUDY_SPD555303)
    Cc ~ lnorm(expSd)
  })
}
