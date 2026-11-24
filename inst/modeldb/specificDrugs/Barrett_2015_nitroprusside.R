Barrett_2015_nitroprusside <- function() {
  description <- "Population K-PD model of the mean arterial pressure (MAP) response to intravenous sodium nitroprusside (SNP) infusion in children aged 0-17 years under anaesthesia needing deliberate hypotension or controlled normotension (Barrett 2015). SNP plasma concentrations were not measurable, so dose enters a single virtual effect compartment (depot_kpd) whose volume is fixed at 1 L/70 kg; its clearance and volume scale allometrically with weight (exponents 0.75 and 1). MAP is an inhibitory sigmoid Emax function of the effect-compartment concentration subtracted from an age-dependent baseline (power of age, reference 11.33 years) plus a linear drift in time. A two-class mixture on EC50 (high EC50 460 ug/L, probability 0.699; low EC50 138 ug/L) is carried by the binary covariate MIX_LOW_EC50. Emax is held as a logit-transformed fraction of each subject's baseline (typical 22.65 of 76 mmHg); the drift slope has a proportional random effect that allows either sign. Proportional residual error on MAP."
  reference <- paste(
    "Barrett JS, Hirankarn S, Holford N, Hammer GB, Drover DR, Cohane CA,",
    "Anderson B, Dombrowski E, Reece T, Zajicek A, Schulman SR.",
    "A hemodynamic model to guide blood pressure control during deliberate",
    "hypotension with sodium nitroprusside in children.",
    "Front Pharmacol. 2015;6:151.",
    "doi:10.3389/fphar.2015.00151.",
    sep = " "
  )
  vignette <- "Barrett_2015_nitroprusside"

  # Dose in ug and time in hours, so an infusion of R ug/kg/min enters as a
  # rate of R * WT * 60 ug/h. The effect-compartment concentration is ug/L and
  # MAP is mmHg.
  units <- list(time = "h", dosing = "ug", concentration = "ug/L", map = "mmHg")

  compartmentData <- list(
    depot_kpd = list(
      analyte = "sodium nitroprusside",
      units = "ug",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Theory-based allometry on the effect-compartment clearance (exponent 0.75) and volume (exponent 1), reference 70 kg (Equations 3-4). Table 1 weight 2.8-112.2 kg, median 31.6 kg.",
      source_name = "WT"
    ),
    AGE = list(
      description = "Post-natal age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power function on the baseline MAP S0, reference 11.33 years (Equation 5; the cohort median, Table 1 reports 11.3). Table 1 age 0.01-17.01 years. All neonates were full term, so post-natal age is the whole of the subject's age.",
      source_name = "Age"
    ),
    MIX_LOW_EC50 = list(
      description = "Latent mixture-model class indicator for SNP sensitivity: 1 = low-EC50 (more sensitive) subpopulation, 0 = high-EC50 (less sensitive) subpopulation",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (high-EC50 subpopulation; the majority class, probability 0.699)",
      notes = "Not a measured covariate: the class is the NONMEM mixture assignment on EC50 (Methods 'PD model'; Equation 1). Table 2 'Probability of high EC50 subpopulation' = 0.699, so P(MIX_LOW_EC50 = 1) = 0.301. The paper found no subject characteristic that predicts the class (Discussion). For population simulation draw MIX_LOW_EC50 ~ Bernoulli(0.301) per subject; for a typical subject of the majority class set MIX_LOW_EC50 = 0.",
      source_name = "NONMEM mixture subpopulation (low vs high EC50)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 202,
    n_studies = 1,
    age_range = "0.01-17.01 years (median 11.3); 4 neonates, 51 infants, 47 children, 100 adolescents",
    weight_range = "2.8-112.2 kg (median 31.6)",
    sex_female_pct = 64,
    race_ethnicity = c(
      White = 81.2,
      Black = 7.92,
      `American Indian` = 6.4,
      Asian = 3.5,
      `Pacific Islander` = 0.5,
      Other = 0.5
    ),
    disease_state = "Children under general anaesthesia (surgery) who required deliberate hypotension or controlled normotension for an anticipated duration of at least 2 h.",
    dose_range = "Blinded phase: randomised to 0.3, 1.0, 2.0 or 3.0 ug/kg/min continuous IV infusion for up to 30 min with titration allowed; observed dose rates 0.09-3.04 ug/kg/min (median 0.66) over 1-30 min infusion segments (median 5).",
    regions = "Multi-centre trial NICHD-2003-09-DR under the US Best Pharmaceuticals for Children Act.",
    notes = "3038 MAP measurements from 202 of 203 enrolled children in the blinded dose-ranging phase (including baseline values up to 6 min before the infusion) were used for model building; 10303 MAP values from 197 children in the open-label titration phase were used for external evaluation. Hispanic ethnicity 23.3%. Table 1."
  )

  ini({
    # Effect-compartment kinetics (Table 2 'Model' column; Equations 3-4).
    lcl <- log(3.12); label("Effect-compartment clearance for a 70 kg subject (L/h)") # Table 2 'CL (L/h/70 kg)' Model = 3.12 (bootstrap 3.13, RSE 9.6%)
    lvc <- fixed(log(1)); label("Effect-compartment volume for a 70 kg subject (L)") # Table 2 'V (L/70 kg)' = 1 (FIXED); Methods 'nominal volume of distribution (V) of 1 L/70 kg'
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL (unitless)") # Equation 3 exponent 0.75, theory-based allometry
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on V (unitless)") # Equation 4 linear (WT/70), theory-based allometry

    # Baseline MAP (Equations 2 and 5).
    lrbase <- log(76); label("Baseline MAP S0 at the reference age of 11.33 years (mmHg)") # Table 2 'S0 (mm Hg)' Model = 76 (RSE 1.0%)
    e_age_rbase <- 0.0333; label("Power exponent of (AGE/11.33) on baseline MAP (unitless)") # Table 2 'Power of age on S0' Model = 0.0333 (Equation 5 prints 0.0338)

    # Inhibitory sigmoid Emax (Equation 1).
    # Emax is estimated as a logit-transformed fraction of the baseline (Methods
    # 'PD model'; Table 2 footnote e). Table 2 prints it in mmHg at the typical
    # S0, so the typical fraction is 22.65 / 76 and its logit is
    # log(22.65 / (76 - 22.65)).
    logitemax <- log(22.65 / (76 - 22.65)); label("Logit of the typical maximal MAP decrease as a fraction of baseline, Emax/S0 (unitless)") # Table 2 'Emax (mm Hg)' Model = 22.65 (RSE 7.6%) with 'S0 (mm Hg)' = 76
    lec50_high <- log(460); label("EC50 in the high-EC50 subpopulation (ug/L)") # Table 2 'EC50 (ug/L) - High' Model = 460 (RSE 10.3%)
    lec50_low <- log(138); label("EC50 in the low-EC50 subpopulation (ug/L)") # Table 2 'EC50 (ug/L) - Low' Model = 138 (bootstrap average 104.3, RSE 21.6%)
    lhill <- log(7.23); label("Hill coefficient of the sigmoid Emax function (unitless)") # Table 2 'Hill coefficient' Model = 7.23 (RSE 11.2%)

    # Linear drift of MAP in time ('disease progression', Equation 1 alpha).
    slope <- 12.3; label("Linear slope of MAP in time (mmHg/h)") # Table 2 'Disease-progression slope' Model = 12.3; per-hour unit from Figures 4-5

    # Between-subject variability. Table 2 footnote d: CV% = sqrt(omega^2) x 100,
    # so omega^2 = (CV/100)^2. Only the bootstrap averages are printed for BSV.
    etalcl ~ 0.393129 # Table 2 CL BSV 62.7%; 0.627^2
    etalvc ~ 0.6724 # Table 2 V BSV 82.0%; 0.820^2
    etalrbase ~ 0.017424 # Table 2 S0 BSV 13.2%; 0.132^2
    etalogitemax ~ 0.585225 # Table 2 Emax BSV 76.5% on the logit scale (footnote e); 0.765^2
    etalec50_high ~ 0.011664 # Table 2 EC50-High BSV 10.8%; 0.108^2
    etalec50_low ~ 0.048841 # Table 2 EC50-Low BSV 22.1%; 0.221^2
    etaslope ~ 4.4521 # Table 2 slope BSV 211%, proportional random effect; 2.11^2

    propSd <- 0.329; label("Proportional residual error on MAP (fraction)") # Table 2 'Residual variability (%CV)' Model = 32.9
  })

  model({
    # Effect-compartment clearance and volume with theory-based allometry
    # (Equations 3-4), reference 70 kg.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    kel <- cl / vc

    d/dt(depot_kpd) <- -kel * depot_kpd

    # Effect-compartment SNP concentration (ug/L).
    Ce <- depot_kpd / vc

    # Baseline MAP: power of post-natal age, reference 11.33 years (Equation 5).
    rbase <- exp(lrbase + etalrbase) * (AGE / 11.33)^e_age_rbase

    # Emax as a fraction of the subject's baseline, logit-transformed so the
    # random effect keeps it between 0 and the baseline (Methods 'PD model';
    # Table 2 footnote e).
    logitfemax <- logitemax + etalogitemax
    emax <- rbase / (1 + exp(-logitfemax))

    # Two-class mixture on EC50 (Equation 1).
    ec50_high <- exp(lec50_high + etalec50_high)
    ec50_low <- exp(lec50_low + etalec50_low)
    ec50 <- ec50_high * (1 - MIX_LOW_EC50) + ec50_low * MIX_LOW_EC50
    hill <- exp(lhill)

    # Proportional random effect on the drift slope, so individual slopes may
    # take either sign (Methods 'PD model').
    slope_i <- slope * (1 + etaslope)

    # Time zero is the start of the SNP infusion.
    map <- rbase - emax * Ce^hill / (ec50^hill + Ce^hill) + slope_i * t

    map ~ prop(propSd)
  })
}
