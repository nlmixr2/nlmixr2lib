Bellanti_2015_deferoxamine <- function() {
  description <- "Indirect-response disease model for serum ferritin in chronic transfusional iron overload (beta-thalassaemia major) with proportional deferoxamine effect on ferritin degradation rate. Two-compartment 8-h SC-infusion deferoxamine PK (literature-derived) enters as a time-varying steady-state concentration covariate; the ferritin compartment captures the baseline turnover (Kin, Kout), the disease-status-modulated transfusion-driven production (CRT), and the chelator effect (1 + DFO) on Kout."
  reference <- "Bellanti F, Del Vecchio GC, Putti MC, Cosmi C, Fotzi I, Bakshi SD, Danhof M, Della Pasqua O. Model-Based Optimisation of Deferoxamine Chelation Therapy. Pharm Res. 2016 Feb;33(2):498-509. doi:10.1007/s11095-015-1805-0 (published online 10 November 2015)."
  vignette <- "Bellanti_2015_deferoxamine"

  units <- list(
    time = "hour",
    dosing = "mg",
    concentration = "ug/L"
  )

  covariateData <- list(
    CSS_DFO = list(
      description = "Average steady-state plasma deferoxamine concentration driving the chelation effect on ferritin degradation",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying input. Bellanti 2015 generates this externally from a literature-derived two-compartment, zero-order-absorption (8-h SC infusion), first-order-elimination DFO PK model (CL/F = 19.3 L/h, Q/F = 17.6 L/h, V/F = 77.4 L, Vp/F = 238 L at an adult 70-kg reference; allometric exponents 0.75 on clearances and 1.00 on volumes for paediatric / adolescent scaling) and feeds it into the ferritin ODE as the drug-exposure driver in the linear concentration-effect term DFO = SLP * CSS_DFO on Kout (Bellanti 2015 Eq. 3 and Eq. 4). For new simulations: supply a constant typical-Css value, or a time-varying column that switches to 0 during drug holidays. The paper's compliance-corrected effective concentration TCss_AV = SCss_AV * (1 - CMPL) (Eq. 7 and Eq. 8) collapses into a single time-varying CSS_DFO by precomputing the (1 - CMPL) reduction in the input data. Typical-Css reference values: ~3.5 ug/mL for 30 mg/kg/day, ~5.5 ug/mL for 45 mg/kg/day, ~7.5 ug/mL for 60 mg/kg/day at 45 kg body weight on the 5-days-per-week 8-h SC infusion schedule (Bellanti 2015 Fig 3).",
      source_name = "SCssAV / TCssAV"
    ),
    FERRITIN_BL = list(
      description = "Per-subject baseline serum ferritin concentration at t = 0 (initial condition for the ferritin compartment)",
      units = "ug/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Subject-level initial condition. Bellanti 2015 fits to retrospective clinical data where each patient enters with an observed pre-treatment ferritin value; the cohort median was 2260 ug/L and the range 393-8500 ug/L across n=27 transfusion-dependent beta-thalassaemia major patients (Bellanti 2015 Table I). For simulation use: supply the desired starting ferritin per virtual subject -- the disease-progression and chelation trajectories are highly sensitive to the starting state because the SCL / SHP feedback terms (Bellanti 2015 Eq. 5 / Eq. 6) are non-linear functions of the current ferritin value.",
      source_name = "FERRITIN at baseline"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 27,
    n_studies = 1,
    age_range = "6.8-19.9 years",
    age_median = "14.6 years",
    weight_range = "17.5-71 kg",
    weight_median = "46 kg",
    sex_female_pct = NA_real_,
    disease_state = "Transfusion-dependent beta-thalassaemia major receiving deferoxamine monotherapy for iron chelation; baseline serum ferritin median 2260 ug/L (range 393-8500); chronic-transfusion-dependent haemoglobinopathy",
    dose_range = "Most prevalent regimen 40 mg/kg/day deferoxamine as 8-h SC infusion 5 days/week; recommended range 20-60 mg/kg/day. Doses were case-specifically adjusted over the up-to-10-year follow-up.",
    regions = "Italy (three centres: A.O. Universitaria Consorziale Policlinico di Bari, A.O. Universitaria Policlinico di Sassari, A.O. di Padova)",
    notes = "Baseline characteristics from Bellanti 2015 Table I (median, range): height 154 cm (111-173); TSH 2.34 mIU/L (0.58-83.2); FT4 1.05 ng/dL (0.73-1.43); AST 33 U/L (7-159); ALT 56 U/L (9-372); glucose 91 mg/dL (52-444); creatinine 0.6 mg/dL (0.2-1.12); ejection fraction 64 percent (35-77); ferritin 2260 ug/L (393-8500). Retrospective clinical data spanning up to 10 years; patients contributed a mean of 40.2 ferritin observations each (sd 17) with a minimum of 4 samples per year, sampled approximately every 2-3 months. NONMEM 7.2.0 FOCE estimation; bootstrap 1000 samples in PsN v3.5.3. The paediatric / adolescent cohort body weights span 17.5-71 kg; the paper applies fixed allometric scaling (0.75 on clearances, 1.00 on volumes) to extrapolate the adult literature-derived DFO PK to this body-weight range when generating CssAV externally."
  )

  ini({
    # ----- Baseline disease-model parameters (Bellanti 2015 Table III) -----
    # All four were estimated in an unpublished prior analysis and held FIX in this paper;
    # values reproduced inline from Table III. The unpublished provenance is documented in
    # the validation vignette's Assumptions and deviations section.

    # Units in Table III for Kin and Kout: ug/h and 1/h -- the ug/h label is a publication
    # typo; for dimensional consistency in the ODE d(F)/dt = Kin + CRT - Kout * F with F in ug/L
    # and Kout in 1/h, Kin must have units of ug/L/h. This interpretation reproduces a plausible
    # drug-free baseline ferritin Kin/Kout = 0.0002 / 4.5e-6 = 44.4 ug/L (~ lower end of healthy
    # range), matching the model interpretation of CRT as the disease-state perturbation.
    lkin    <- fixed(log(0.0002))
    label("Baseline ferritin zero-order production rate Kin (ug/L/h)")  # Bellanti 2015 Table III (FIX, prior unpublished disease model)
    lkout   <- fixed(log(4.5e-6))
    label("Baseline ferritin first-order degradation rate Kout (1/h)")   # Bellanti 2015 Table III (FIX, prior unpublished disease model)

    # SHP units in Table III are reported as 1/h; for dimensional consistency with the
    # exp(-SHP * FERRITIN) term in Eq. 2 (FERRITIN in ug/L), SHP must have units of L/ug.
    # SHP_ref * FERRITIN_median = 0.00026 * 2260 = 0.59, a sensible dimensionless exponent.
    lsclref <- fixed(log(0.383))
    label("Reference scaling factor SCL_ref on transfusion-driven production (ug/L/h)")  # Bellanti 2015 Table III (FIX, prior unpublished disease model)
    lshpref <- fixed(log(0.00026))
    label("Reference shape factor SHP_ref in exp(-SHP * FERRITIN) (L/ug)")  # Bellanti 2015 Table III (FIX, prior unpublished disease model)

    # ----- Disease-status feedback exponents (estimated; Bellanti 2015 Table III Eq. 5/Eq. 6) -----
    e_dis_scl <- 0.845
    label("Exponent on (FERRITIN / FERRITIN_MED) feedback on SCL (Eq. 5)")  # Bellanti 2015 Table III row DIS-exp-on-SCL
    e_dis_shp <- 1.29
    label("Exponent on (FERRITIN / FERRITIN_MED) feedback on SHP (Eq. 6)")  # Bellanti 2015 Table III row DIS-exp-on-SHP

    # ----- Drug effect (estimated; Bellanti 2015 Table III Eq. 4) -----
    # SLP units: 1/(ug/mL) so DFO = SLP * CSS_DFO is dimensionless inside Kout * (1 + DFO).
    lslp <- log(4.81)
    label("Slope of deferoxamine concentration-effect on Kout (per ug/mL)")  # Bellanti 2015 Table III row Slope

    # ----- Between-subject variability -----
    etalslp ~ 0.082   # Bellanti 2015 Table III row IIV-on-Slope (variance, log-normal IIV)

    # Bellanti 2015 reports inter-occasion variability (IOV) of 0.252 (variance) on the
    # conversion rate CRT. nlmixr2lib has no canonical pattern for occasion-keyed IOV without
    # an OCC column, so this is encoded as ordinary IIV on a fixed-at-1 multiplicative CRT
    # factor lcrt_mult (typical value exp(log(1)) = 1; the IIV variance is the paper's IOV
    # variance). The model loses the within-subject between-occasion drift the paper estimates
    # but preserves the population-level CRT spread. See vignette Assumptions and deviations.
    lcrt_mult <- fixed(log(1))
    label("Multiplicative CRT factor anchor (log scale, fixed at log(1) = 0; carries IIV originally reported as IOV)")  # structural anchor for the IOV-encoded-as-IIV pattern
    etalcrt_mult ~ 0.252   # Bellanti 2015 Table III row IOV-on-CRT (variance; encoded as IIV on the log-1 anchor)

    # ----- Residual error (proportional on serum ferritin) -----
    # Bellanti 2015 Table III lists the proportional error magnitude as -0.173. The negative
    # sign on a NONMEM proportional-error THETA is a well-known sign artefact -- the residual
    # epsilon distribution is symmetric around zero, so only |W| is operative. Encoded here
    # as the absolute magnitude.
    propSd <- 0.173
    label("Proportional residual error on serum ferritin (fraction)")  # Bellanti 2015 Table III row Error-proportional, abs value
  })

  model({
    # ----- Mechanism constant -----
    # Population median ferritin used as the normalising denominator in the disease-status
    # feedback equations (Bellanti 2015 Eq. 5 / Eq. 6). Value from Bellanti 2015 Table I.
    FERRITIN_MED <- 2260   # ug/L

    # ----- Typical-value disease-model parameters -----
    kin    <- exp(lkin)
    kout   <- exp(lkout)
    sclref <- exp(lsclref)
    shpref <- exp(lshpref)

    # ----- Individual drug-effect slope (log-normal IIV) -----
    slp <- exp(lslp + etalslp)

    # ----- Multiplicative CRT factor (typical value 1; carries IIV-encoded-as-IOV) -----
    crt_mult <- exp(lcrt_mult + etalcrt_mult)

    # ----- Disease-status feedback (Bellanti 2015 Eq. 5 / Eq. 6) -----
    # SCL_i and SHP_i depend on the current ferritin value; the disease state continuously
    # modulates its own production dynamics through these non-linear feedback terms.
    fr <- central / FERRITIN_MED
    scl <- sclref * fr^e_dis_scl
    shp <- shpref * fr^e_dis_shp

    # ----- Transfusion-driven production CRT (Bellanti 2015 Eq. 2) -----
    crt <- scl * exp(-shp * central) * crt_mult

    # ----- Deferoxamine effect on ferritin degradation (Bellanti 2015 Eq. 3 / Eq. 4) -----
    # Linear concentration-effect on Kout via DFO = SLP * CSS_DFO. The paper drives this with
    # an externally-computed average steady-state concentration CssAV; this implementation
    # accepts CSS_DFO as a time-varying covariate so a user can encode any dose schedule
    # (constant typical-Css, intermittent dosing, drug holidays) directly.
    dfo <- slp * CSS_DFO

    # ----- Ferritin ODE in canonical 'central' compartment (Bellanti 2015 Eq. 3) -----
    # 'central' here is the plasma ferritin pool (ug/L); compartment name follows the
    # nlmixr2lib endogenous-biomarker convention (see comparable model Bisaso_2014_albumin
    # which uses 'central' for the plasma albumin pool).
    # d(F)/dt = Kin + CRT - Kout * F * (1 + DFO)
    # Units: [ug/L/h] = [ug/L/h] + [ug/L/h] - [1/h] * [ug/L] * [-]
    d/dt(central) <- kin + crt - kout * central * (1 + dfo)

    # ----- Initial condition: per-subject baseline ferritin -----
    central(0) <- FERRITIN_BL

    # ----- Observation -----
    # The single observation in Bellanti 2015 is serum ferritin (the disease biomarker); the
    # deferoxamine concentration enters as a covariate input (CSS_DFO) rather than as a
    # PK output. Cc here is the plasma-pool ferritin concentration in ug/L (single-output
    # naming convention per nlmixr2lib).
    Cc <- central
    Cc ~ prop(propSd)
  })
}
