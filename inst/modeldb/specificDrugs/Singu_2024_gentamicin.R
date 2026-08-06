Singu_2024_gentamicin <- function() {
  description <- "One-compartment IV population PK model for gentamicin in 50 PK-evaluable Namibian neonates (82.7% preterm) with suspected or confirmed sepsis receiving 5 mg/kg every 24 h, developed from a truncated (7.5 h) two-sample informative-block randomised sampling design (Singu 2024). Clearance carries an allometric birth-weight effect (exponent 1.30, reference 1.57 kg), a power effect of log-transformed white blood cell count (exponent -0.560, reference log-WBC 2.39), and a logistic postnatal-age maturation function FMAT = PNA^GAMMA / (PNA^GAMMA + PNA50^GAMMA) with GAMMA and PNA50 fixed at 0.551 and 0.0332 years; volume carries an allometric birth-weight effect (exponent 1.76, reference 1.57 kg). Residual variability is a log-scale combination (additive + proportional) model with two concentration-dependent magnitudes switched at a predicted concentration of 2264.25 ng/mL (log(concentration) = 7.725), an approach the authors take from Hussein 2005; see modellib('Hussein_2005_levocetirizine')."
  reference <- paste(
    "Singu BS, Verbeeck RK, Pieper CH, Ette EI.",
    "Confirming the suitability of a gentamicin dosing strategy in neonates",
    "using the population pharmacokinetic approach with truncated sampling",
    "duration.",
    "Children (Basel). 2024;11(8):898.",
    "doi:10.3390/children11080898.",
    sep = " "
  )
  vignette <- "Singu_2024_gentamicin"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Residual-error parameters that do not follow the bare canonical
  # propSd / addSd / expSd naming because the source stratifies residual
  # variability by predicted concentration (see ini() and model()).
  paper_specific_residual_sds <- c("propSdLow", "addSdLow", "propSdHigh", "addSdHigh")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central = list(analyte = "gentamicin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT_BIRTH = list(
      description        = "Birth weight (time-fixed per subject).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric (power) effect on both CL and Vc, centred on the cohort",
        "median birth weight 1.57 kg (Table 1, range 0.90-3.92 kg). Table 3",
        "footnote: 'CLTV and VTV are the clearance and volume of distribution",
        "for the typical (i.e., average or reference) subject in the dataset",
        "weighing 1.57 kg'. Estimated exponents are CL_WT = 1.30 and",
        "V_WT = 1.76 (Table 3). The paper notes the 90% bootstrap CI for the",
        "birth-weight effect on CL included 1.0 but the term was retained",
        "'for a biological reason (i.e., dosing based on allometry)'",
        "(Results Section 3.4). The covariate is BIRTH weight, not a",
        "time-varying current weight."
      ),
      source_name        = "WT"
    ),
    PNA = list(
      description        = "Postnatal age (chronological time since birth).",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Singu 2024 reports postnatal age in DAYS in Table 1 (median 4.0,",
        "range 1.0-17) and in YEARS inside the maturation function (the",
        "typical neonate has PNA 4 days = 0.011 years; PNA50 = 0.0332",
        "years = 12.1 days). The canonical PNA column is in MONTHS, so",
        "model() converts with PNA / 12 before evaluating the maturation",
        "function, keeping lpna50 on the paper's year scale. Drives the",
        "logistic maturation function FMAT = PNA^GAMMA /",
        "(PNA^GAMMA + PNA50^GAMMA) on CL (Results Section 3.4 and Table 2",
        "footnote). The first day of life is defined as day 1, so PNA is",
        "strictly positive and FMAT is well defined."
      ),
      source_name        = "PNA"
    ),
    WBC = list(
      description        = "White blood cell count (marker of sepsis).",
      units              = "10^9 cells/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters CL as a power effect on the ratio of NATURAL-LOG-transformed",
        "values: (log(WBC) / 2.39)^-0.560. Equation 3 states 'WBC was",
        "log-transformed' and Results Section 3.4 defines the typical neonate",
        "by a 'median log-transformed WBC of 2.39'; the cohort median WBC is",
        "11.0 x 10^9/L (Table 1, range 1.67-37.4) and log(11.0) = 2.3979,",
        "confirming natural-log transformation and that 2.39 is the reference",
        "on the LOG scale (not on the raw count scale). The negative exponent",
        "makes higher WBC lower CL, which the Discussion attributes to",
        "sepsis-associated inflammation causing acute kidney injury and a",
        "decline in eGFR. This is the first report of WBC as a covariate on",
        "gentamicin CL in neonates."
      ),
      source_name        = "WBC"
    )
  )

  covariatesDataExcluded <- list(
    GA = list(
      description = "Gestational age at birth (weeks). Screened during covariate model building.",
      units       = "weeks",
      type        = "continuous",
      notes       = paste(
        "Screened; not retained. Cohort median 32 weeks (range 24-40;",
        "Table 1); 43/52 (82.7%) were preterm (< 37 weeks). Methods Section",
        "2.4.4 lists gestational age among the covariates tested",
        "('Birthweight, GA, PNA, white blood cell count (WBC), and the",
        "reciprocal of serum creatinine (1/SCr) were tested for inclusion')."
      )
    ),
    CREAT = list(
      description = "Serum creatinine (mg/dL); tested as its reciprocal 1/SCr. Screened during covariate model building.",
      units       = "mg/dL",
      type        = "continuous",
      notes       = paste(
        "Screened as the reciprocal 1/SCr to approximate kidney function and",
        "'to reduce collinearity in the covariate vector' (Results Section",
        "3.4); not retained. Cohort median 0.72 mg/dL (range 0.20-1.66;",
        "Table 1). Results Section 3.4: 'the 1/SCr showed no apparent",
        "relationship with these parameters'; Table 2 run 6 (Cl~RSCR) gave",
        "LLD -2.533 on 1 df, below the 3.84 critical value, so not",
        "significant."
      )
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 50L,
    n_enrolled      = 52L,
    n_studies       = 1L,
    n_observations  = 100L,
    age_range       = "GA 24-40 weeks (median 32) at birth; PNA 1.0-17 days (median 4.0) at sampling",
    age_median      = "GA 32 weeks; PNA 4.0 days (0.011 years)",
    weight_range    = "birth weight 0.90-3.92 kg",
    weight_median   = "birth weight 1.57 kg",
    sex_female_pct  = 44.2,
    race_ethnicity  = "Not reported (single-centre Namibian neonatal unit).",
    ga_range        = "24-40 weeks (median 32); 43/52 (82.7%) preterm at < 37 weeks GA",
    disease_state   = paste(
      "Neonates admitted to the Neonatal Unit of the Maternity Ward,",
      "Windhoek Central Hospital, with suspected or confirmed sepsis.",
      "Sepsis was suspected on clinical symptoms (temperature spike,",
      "tachycardia, lethargy, increased respiratory rate) and confirmed",
      "with elevated inflammatory markers (C-reactive protein and white",
      "blood cell count; combined sensitivity 90.3%). Neonates with anemia",
      "or congenital anomalies were excluded. Renal impairment was not an",
      "issue in this population (Discussion Section 4.2). Median WBC",
      "11.0 x 10^9/L (range 1.67-37.4); median serum creatinine 0.72 mg/dL",
      "(range 0.20-1.66); median height 41 cm (range 30-53)."
    ),
    dose_range      = paste(
      "Gentamicin 5 mg/kg as an intravenous bolus over 3-5 s every 24 h",
      "(median dose 7.9 mg, range 4.0-17 mg), co-administered with",
      "benzylpenicillin 100,000 IU/kg every 12 h or ampicillin 50 mg/kg",
      "every 8 h."
    ),
    co_medication   = "Benzylpenicillin 100,000 IU/kg q12h or ampicillin 50 mg/kg q8h.",
    regions         = "Single centre (Windhoek Central Hospital, Windhoek, Namibia).",
    notes           = paste(
      "Prospective, non-randomised observational study. 52 neonates",
      "enrolled, 50 PK-evaluable, each contributing exactly two serum",
      "samples for 100 observations total. Sampling used the informative",
      "PK profile (block) randomised design of Jones 1996 with a TRUNCATED",
      "7.5 h duration: 8% of samples from the 0.08-0.14 h block, 68% from",
      "the 0.14-4.2 h block, and 24% from the 4.2 h and later block; 38%",
      "of samples were drawn less than 1.0 h after the dose and the median",
      "(range) interval between a subject's two samples was 3.17 h",
      "(0.17-6.5 h). Gentamicin assayed on an Indiko Plus autoanalyser",
      "(LLOQ 0.3 ug/mL); serum creatinine by the kinetic alkaline picrate",
      "Jaffe method on a Cobas 6000. Fit in NONMEM 7.4.1 with FOCE-I on",
      "natural-log-transformed concentrations. The paper's two",
      "log-concentration scales differ by log(1000): Figure 1 plots the",
      "data on a log mg/L axis (spanning about -2.2 to 4.1) while the",
      "residual-error threshold is quoted as the natural log of the ng/mL",
      "value (log(2264.25 ng/mL) = 7.725). The physical threshold of",
      "2264.25 ng/mL = 2.264 mg/L is unambiguous and is what this model",
      "uses. Final model OFV",
      "-22.151 (base model 15.076). Shrinkage 3.19% (eta CL), 3.18%",
      "(eta V), 1.26% (eta PNA50). Baseline demographics from Table 1;",
      "parameter estimates from Table 3; covariate model-building run log",
      "in Table 2. The reported 90% bootstrap confidence intervals in",
      "Table 3 are internally inconsistent: the omega^2 PNA50 row gives a",
      "point estimate of 6.23 with a CI of 9.10-13.4, which excludes its",
      "own point estimate, and the omega^2 CL and omega^2 CL:V rows report",
      "negative lower bounds (-0.289 and -0.221) for quantities that",
      "cannot be negative. That column is therefore not used here to judge",
      "parameter identifiability. Table 4 (simulated exposure by postnatal",
      "age and birth weight) is likewise inconsistent with Table 3: its",
      "median Cmax is flat across birth-weight strata (6.83-7.22 mg/L),",
      "which requires V proportional to weight with an exponent of 1.0,",
      "whereas Table 3 reports V_WT = 1.76; and its median peak:trough",
      "ratio of 5.48 over 24 h implies a half-life of 9.8 h against the",
      "4.2 h the paper states. See the vignette for the full analysis."
    )
  )

  ini({
    # ----- Structural PK parameters (Singu 2024 Table 3, final model) -----
    # CL_TV,REF is the clearance BEFORE the maturation multiplier FMAT is
    # applied, at the reference birth weight 1.57 kg and reference
    # log-transformed WBC 2.39. Cross-check of the whole CL chain against
    # the paper's own typical-value statement (Results Section 3.4: 'For
    # the typical (average) neonate in the PK dataset (median weight 1.57
    # kg, median postnatal age of 4 days (0.011 years), median
    # log-transformed WBC of 2.39), the predicted CL and V are 0.069 L/h
    # and 0.417 L'):
    #   FMAT = 0.011^0.551 / (0.011^0.551 + 0.0332^0.551) = 0.35236
    #   CL   = 0.196 * 1 * 1 * 0.35236 = 0.06906 L/h  -> matches 0.069 L/h
    # This also falsifies the alternative reading in which the denominator
    # carries a bare PNA50 rather than PNA50^GAMMA, which would give
    # FMAT = 0.71510 and CL = 0.140 L/h.
    lcl <- log(0.196); label("Clearance before maturation, at WT_BIRTH = 1.57 kg and log(WBC) = 2.39 (L/h)") # Singu 2024 Table 3 'CL (L/h)' = 0.196
    lvc <- log(0.417); label("Volume of distribution at WT_BIRTH = 1.57 kg (L)")                             # Singu 2024 Table 3 'V (L)' = 0.417

    # ----- Allometric birth-weight exponents (Singu 2024 Table 3) -----
    # Both estimated (not fixed). The CL exponent's 90% bootstrap CI
    # (0.558, 1.68) includes 1.0; Results Section 3.4 states the term was
    # nonetheless retained 'for a biological reason (i.e., dosing based on
    # allometry, Table 3)'.
    e_wt_birth_cl <- 1.30; label("Power exponent of (WT_BIRTH / 1.57 kg) on CL (unitless)") # Singu 2024 Table 3 'CL_WT' = 1.30
    e_wt_birth_vc <- 1.76; label("Power exponent of (WT_BIRTH / 1.57 kg) on Vc (unitless)") # Singu 2024 Table 3 'V_WT'  = 1.76

    # ----- White blood cell effect on CL (Singu 2024 Table 3) -----
    # Applied to the ratio of NATURAL-LOG-transformed WBC against the
    # reference log-WBC of 2.39 (Equation 3 with 'WBC was log-transformed';
    # Results Section 3.4 'median log-transformed WBC of 2.39').
    e_wbc_cl <- -0.560; label("Power exponent of (log(WBC) / 2.39) on CL (unitless)") # Singu 2024 Table 3 'CL_WBC' = -0.560

    # ----- Postnatal-age maturation of CL (Singu 2024 Table 3) -----
    # FMAT = PNA^GAMMA / (PNA^GAMMA + PNA50^GAMMA) with PNA in YEARS
    # (Results Section 3.4 and the Table 2 footnote). Table 3 reports both
    # GAMMA and PNA50 with '90% Bootstrap Confidence Interval = Fixed' and
    # the Table 3 footnote states 'GAMMA and PNA50 were initially estimated
    # and later fixed in the final model', so both are wrapped in fixed().
    # PNA50 nonetheless retains an ESTIMATED inter-individual variance (see
    # etalpna50 below).
    e_pna_cl_hill <- fixed(0.551);      label("Steepness (GAMMA) of the postnatal-age maturation function on CL (unitless)") # Singu 2024 Table 3 'GAMMA' = 0.551, Fixed
    lpna50        <- fixed(log(0.0332)); label("Postnatal age at 50% of mature CL (years)")                                   # Singu 2024 Table 3 'PNA 50 (yr)' = 0.0332, Fixed

    # ----- Inter-individual variability (Singu 2024 Table 3) -----
    # Reported directly as VARIANCES on the log scale (Table 3 footnote a:
    # 'Intersubject in CL, V, and PNA50 ... estimated as their variances,
    # and the numbers in parentheses are the respective intersubject
    # expressed as percentages'). The parenthesised percentages are simply
    # 100 * sqrt(variance): sqrt(1.00e-5) = 0.316% (printed 0.31%),
    # sqrt(0.501) = 70.78% (printed 70.8%), sqrt(6.23) = 249.6% (printed
    # 250%). No CV%-to-variance back-transformation is therefore needed.
    #
    # omega^2 CL is essentially zero. Discussion Section 4.1: 'The
    # infinitely small variability in CL obtained was because the spread of
    # the PK samples in the region of the PK profile containing information
    # for CL estimation was inadequate for its estimation. ... The estimate
    # was not fixed because it was important to estimate the covariance
    # between CL and V at the individual subject level, which was essential
    # for the efficient simulation of gentamicin concentrations.' It is
    # therefore ESTIMATED (not fixed) and the CL-V covariance is retained.
    # The 2x2 block is positive definite as printed:
    #   det = 1.00e-5 * 0.501 - 0.00024^2 = 4.95e-6 > 0
    #   corr = 0.00024 / sqrt(1.00e-5 * 0.501) = 0.107
    etalcl + etalvc ~ c(1.00e-5,
                        0.00024, 0.501) # Singu 2024 Table 3: omega^2 CL = 1.00e-5, omega^2 CL:V = 0.00024, omega^2 V = 0.501

    # IIV on PNA50. Discussion Section 4.1: 'The large variability
    # associated with PNA50 results from not having an adequate spread of
    # samples in the region of the profile that contained information on
    # the fixed effect parameter', but 'the parameter estimate was
    # estimated with good precision and negligible bias, as indicated by a
    # minimal shrinkage of 1.26%'. No covariance with etalcl / etalvc is
    # reported, so it is a separate diagonal element.
    etalpna50 ~ 6.23 # Singu 2024 Table 3 'omega^2 PNA 50' = 6.23 (250%)

    # ----- Residual variability (Singu 2024 Table 3, footnote b) -----
    # Footnote b: 'Residual errors are estimated as their standard
    # deviations.' The concentration data were natural-log-transformed and
    # 'the error model was appropriately modified to the logarithm scale'
    # (Methods Section 2.4.3), i.e. a log-transform-both-sides combination
    # model whose log-scale SD is
    #   W = sqrt(propSd^2 + (addSd / Cc)^2).
    #
    # Results Section 3.4 splits that model at a predicted concentration of
    # log(concentration) = 7.725, 'i.e., 2264.25 ng/mL' (exp(7.725) =
    # 2264.25 exactly, confirming natural log and ng/mL for the modelling
    # dataset). The threshold 'was selected based on sensitivity analysis'
    # and 'had no clinical or therapeutic implications'; the approach is
    # taken from Hussein 2005 (reference 35), which is itself in this
    # library as Hussein_2005_levocetirizine and encodes the same
    # concentration-stratified switch.
    #
    # Tier assignment: subscript 1 is the LOW stratum, on three grounds.
    #  (a) Table 3 lists ADDITIVE1 / PROPORTIONAL1 before ADDITIVE2 /
    #      PROPORTIONAL2, matching the order in which Results Section 3.4
    #      names the two strata ('log(concentration) <= 7.725 and >7.725'),
    #      so subscript 1 is the first-named (low) stratum.
    #  (b) The additive terms point the same way. ADDITIVE1 = 4.96 ng/mL is
    #      a magnitude at which an additive term can matter; ADDITIVE2 =
    #      1.00e-5 ng/mL is effectively zero, which is only sensible in the
    #      high-concentration stratum where the additive term is irrelevant.
    #  (c) Figure 1 plots the observed data as log(concentration) against
    #      time with a y-axis spanning about -2.2 to 4.1 -- i.e. natural log
    #      of mg/L (exp(4.1) = 60 mg/L, exp(-2.2) = 0.11 mg/L; a ng/mL axis
    #      would have put every point between 5.7 and 11). The 2264.25 ng/mL
    #      threshold is 2.264 mg/L, which on that axis is log(2.264) = 0.817.
    #      Essentially every plotted observation lies above it: only one
    #      profile (running from about -0.85 to -2.15) and the tail of one
    #      other fall below. The high stratum therefore carries almost all
    #      100 observations, so the tight, well-determined 0.0205 belongs to
    #      it, and the large 0.99 magnitude belongs to the sparse low tail.
    # The one piece of contrary evidence is that the '2' terms carry
    # uninformative bootstrap intervals (PROPORTIONAL2's CI is -1.95, 1.96),
    # which would normally mark them as the data-poor pair. That column is
    # not used to adjudicate the assignment because it is internally
    # inconsistent throughout: the omega^2 PNA50 CI excludes its own point
    # estimate and the omega^2 CL / omega^2 CL:V rows report negative lower
    # bounds for quantities that cannot be negative (see population$notes).
    # The assignment changes only the residual-error magnitude -- the
    # structural model, the covariate model and the IIV are identical either
    # way, so a user simulating Cc (IPRED) is unaffected.
    #
    # Sign. PROPORTIONAL1 is printed as -0.99. A standard deviation is
    # identified only up to sign (the likelihood depends on the variance),
    # so the magnitude 0.99 is used.
    #
    # Units. Table 3's additive terms are on the modelling dataset's ng/mL
    # scale; this model reports Cc in mg/L, so both are divided by 1000.
    # Both are numerically negligible against their own stratum
    # (4.96 ng/mL is 0.22% of the 2264.25 ng/mL threshold), but they are
    # retained because Results Section 3.4 states 'The two components of
    # the residual error were estimated without fixing the additive
    # component to zero. The model could not minimize successfully without
    # it.'
    propSdLow  <- 0.99;    label("Proportional (log-scale) residual SD when predicted Cc <= 2.26425 mg/L (fraction)") # Singu 2024 Table 3 'Res PROPORTIONAL1' = -0.99 (magnitude)
    addSdLow   <- 0.00496; label("Additive residual SD when predicted Cc <= 2.26425 mg/L (mg/L)")                     # Singu 2024 Table 3 'Res ADDITIVE1' = 4.96 ng/mL = 0.00496 mg/L
    propSdHigh <- 0.0205;  label("Proportional (log-scale) residual SD when predicted Cc > 2.26425 mg/L (fraction)")  # Singu 2024 Table 3 'Res PROPORTIONAL2' = 0.0205
    addSdHigh  <- 1e-8;    label("Additive residual SD when predicted Cc > 2.26425 mg/L (mg/L)")                      # Singu 2024 Table 3 'Res ADDITIVE2' = 1.00e-5 ng/mL = 1e-8 mg/L
  })

  model({
    # ----- 1. Derived covariate terms -----
    # Postnatal-age maturation of CL. The canonical PNA column is in months
    # (see covariateData); the paper's maturation function is in years, so
    # convert first. FMAT = PNA^GAMMA / (PNA^GAMMA + PNA50^GAMMA) returns
    # 0.5 at PNA = PNA50 and asymptotes to 1 with increasing PNA
    # (Singu 2024 Equation 2 with theta_eff2 = PNA50^GAMMA, Results Section
    # 3.4, and the Table 2 footnote).
    pna_yr        <- PNA / 12
    pna50_yr      <- exp(lpna50 + etalpna50)
    maturation_cl <- pna_yr ^ e_pna_cl_hill /
      (pna_yr ^ e_pna_cl_hill + pna50_yr ^ e_pna_cl_hill)

    # White blood cell effect on CL, applied to the ratio of the
    # natural-log-transformed count against the reference log-WBC of 2.39.
    f_wbc_cl <- (log(WBC) / 2.39) ^ e_wbc_cl

    # ----- 2. Individual PK parameters (Singu 2024 Equations 3 and 4) -----
    cl <- exp(lcl + etalcl) * f_wbc_cl *
      (WT_BIRTH / 1.57) ^ e_wt_birth_cl * maturation_cl
    vc <- exp(lvc + etalvc) * (WT_BIRTH / 1.57) ^ e_wt_birth_vc

    # ----- 3. Micro-constants -----
    kel <- cl / vc

    # ----- 4. ODE system: one-compartment IV with first-order elimination -----
    # Gentamicin was given as an intravenous bolus over 3-5 s, so no
    # infusion duration is hard-coded; users supply rate / dur per dose in
    # their event table if they wish to model a finite infusion.
    d/dt(central) <- -kel * central

    # ----- 5. Observation and error -----
    Cc <- central / vc

    # Concentration-stratified log-scale combination residual error.
    # Predicted Cc <= 2.26425 mg/L (= 2264.25 ng/mL = exp(7.725) ng/mL)
    # uses the Low pair; Cc > 2.26425 mg/L uses the High pair. The
    # arithmetic switch is equivalent to ifelse() but stays inside the
    # rxode2 expression grammar, following the encoding already used by
    # Hussein_2005_levocetirizine. The combined log-scale SD
    # sqrt(prop^2 + (add/Cc)^2) is the standard log-transform-both-sides
    # form; it is consumed by lnorm(), for which the supplied SD is on the
    # natural-log scale. Cc must be strictly positive for a log-scale
    # residual to be defined, which holds at every post-dose time for this
    # bolus model.
    low_bin    <- (Cc <= 2.26425)
    propSd_eff <- propSdLow * low_bin + propSdHigh * (1 - low_bin)
    addSd_eff  <- addSdLow  * low_bin + addSdHigh  * (1 - low_bin)
    expSd      <- sqrt(propSd_eff ^ 2 + (addSd_eff / Cc) ^ 2)
    Cc ~ lnorm(expSd)
  })
}
