Farrell_2013_conestatAlfa <- function() {
  description <- "One-compartment population PK model with Michaelis-Menten elimination for intravenous recombinant human C1 inhibitor (rhC1INH; conestat alfa; Ruconest) in healthy volunteers and adolescent / adult patients with hereditary angioedema (Farrell 2013). Total functional plasma C1INH is modelled as the sum of an estimated endogenous baseline (separate baselines for healthy volunteers and HAE patients) plus exogenously administered rhC1INH, with the endogenous production rate derived from the Michaelis-Menten elimination at baseline so the no-dose steady state is preserved. Allometric power scaling of central volume on body weight (exponent 0.612)."
  reference   <- "Farrell C, Hayes S, Relan A, van Amersfoort ES, Pijpstra R, Hack CE. Population pharmacokinetics of recombinant human C1 inhibitor in patients with hereditary angioedema. Br J Clin Pharmacol. 2013;76(6):897-907. doi:10.1111/bcp.12132"
  vignette    <- "Farrell_2013_conestatAlfa"
  units       <- list(time = "h", dosing = "U", concentration = "U/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight (baseline)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Power scaling (WT/70)^0.612 on the central volume (Farrell 2013 Table 2 row 'Effect of bodyweight on volume', 95% CI 0.351-0.873). Pooled-cohort range 45-128 kg, median 72 kg.",
      source_name        = "WT"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy participant indicator (1 = healthy volunteer, 0 = patient with hereditary angioedema)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (hereditary angioedema patient)",
      notes              = "Time-fixed per subject. Pooled across asymptomatic and symptomatic HAE for the HAE-patient reference (Farrell 2013 Methods 'Population (healthy volunteers, symptomatic patients and asymptomatic patients) ... were evaluated as covariates'). DIS_HEALTHY gates which endogenous baseline parameter and which baseline-IIV random effect apply for each subject: lrbase_hv (0.901 U/mL) with etalrbase_hv (12.7% CV) for healthy volunteers; lrbase_hae (0.176 U/mL) with etalrbase_hae (54.4% CV) for HAE patients. The cohort-specific baseline IIV magnitudes are reproduced verbatim from Farrell 2013 Table 2.",
      source_name        = "Population (healthy volunteer vs HAE patient)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age at the time of dosing",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in the Farrell 2013 covariate analysis but not retained in the final model (Methods 'Population ..., weight, age, gender, dose and study were evaluated as covariates ... No other relationships were evident between the individual random effects and covariates'). Cohort range 13-66 years, median 33."
    ),
    SEXF = list(
      description = "Female-sex indicator (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened in the Farrell 2013 covariate analysis but not retained in the final model (Methods 'gender ... evaluated as covariate ... No other relationships were evident'). Cohort breakdown: 86 / 133 female (64.7%)."
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 133L,
    n_studies       = 6L,
    n_healthy       = 14L,
    n_patients_asym = 12L,
    n_patients_sym  = 107L,
    age_range       = "13 to 66 years",
    age_median      = "33 years",
    weight_range    = "45 to 128 kg",
    weight_median   = "72 kg",
    sex_female_pct  = round(100 * 86 / 133, 1),
    race_ethnicity  = "Not reported",
    disease_state   = "Pooled adolescent and adult patients with hereditary angioedema (HAE) together with healthy volunteers. The HAE patients comprise asymptomatic individuals dosed between attacks (Study C1 1101-01, n = 12) and patients treated for acute attacks of angioedema (Studies C1 1202-01, C1 1203-01, C1 1205-01, C1 1304-01; n = 107).",
    dose_range      = "Intravenous rhC1INH 6.25-121 U/kg infused over 2 to 30 min; clinical regimen 50 U/kg up to a body weight of 84 kg or a fixed dose of 4200 U (= 2 x 2100 U vials) above 84 kg, with the option of a second 50 U/kg or 2100-4200 U dose within 4 h of the first if the attack had not resolved.",
    regions         = "Multinational clinical-trial programme (six studies pooled).",
    notes           = "Pooled data from six rhC1INH studies (Farrell 2013 Table 1): C1 1101-01 (12 asymptomatic HAE, 6.25-100 U/kg), C1 1106-02 (14 healthy volunteers, 100 U/kg), C1 1202-01 (4 symptomatic HAE, 100 U/kg), C1 1203-01 (10 symptomatic HAE, 100 U/kg), C1 1205-01 (52 symptomatic HAE, 50 and 100 U/kg), C1 1304-01 (41 symptomatic HAE, up to 3 x 2100 U). 294 administrations / 1362 quantifiable concentrations plus 656 below-quantification samples handled with the M4 likelihood method (BLQ included). Chromogenic functional-C1INH assay LLOQ 0.07 U/mL (specialized laboratory for Study 1101-01) or 0.28 U/mL (central contract laboratory for all other studies)."
  )

  ini({
    # ---- Structural PK parameters (Farrell 2013 Table 2 final-model column) ----
    # Reference subject: 70 kg HAE patient (DIS_HEALTHY = 0). The central
    # volume is reported in litres in Table 2; this file stores vc in mL so
    # the dose unit U divided by vc gives concentration directly in U/mL,
    # matching the units used throughout the paper and consistent with the
    # Cc ~ prop(propSd) + add(addSd) error model.
    lvc        <- log(2.86 * 1000); label("Central volume of distribution at 70 kg (V, mL)")               # Farrell 2013 Table 2: V = 2.86 L, 95% CI 2.68-3.03; converted to mL
    lvmax      <- log(1.63);        label("Maximal Michaelis-Menten elimination rate (Vmax, U/mL/h)")       # Farrell 2013 Table 2: Vmax = 1.63 U/mL/h in healthy volunteers and symptomatic HAE patients
    lkm        <- log(1.60);        label("Michaelis-Menten concentration for half-Vmax elimination (Km, U/mL)")  # Farrell 2013 Table 2: Km = 1.60 U/mL, 95% CI 1.14-2.24
    lrbase_hv  <- log(0.901);       label("Endogenous baseline functional C1INH in healthy volunteers (U/mL)")     # Farrell 2013 Table 2: baseline HV = 0.901 U/mL, 95% CI 0.839-0.968
    lrbase_hae <- log(0.176);       label("Endogenous baseline functional C1INH in HAE patients (U/mL)")           # Farrell 2013 Table 2: baseline HAE = 0.176 U/mL, 95% CI 0.154-0.200

    # ---- Covariate effects (Farrell 2013 Table 2) ----
    # The paper applies the power form V = TVV * (WT/70)^e_wt_vc; reference
    # WT = 70 kg is stated in the Methods discussion ('For subjects with the
    # lowest and highest bodyweights in the analysis (45 and 128 kg), volume
    # would be 24% lower and 45% higher compared with a 70 kg subject').
    e_wt_vc <- 0.612; label("Allometric exponent of WT on V (unitless)")                                   # Farrell 2013 Table 2: 'Effect of bodyweight on volume' = 0.612, 95% CI 0.351-0.873

    # ---- IIV (variance on the log-scale ~ omega^2) ----
    # Farrell 2013 Table 2 footnote: 'The %CV for both intersubject and
    # residual variability is an approximation taken as the square root of
    # the variance for the proportional error term * 100.' This is the
    # standard small-variance approximation CV ~= sqrt(omega^2); to preserve
    # the paper's verbatim values, omega^2 is encoded as (CV/100)^2 rather
    # than the exact log(CV^2 + 1) conversion. The lognormal/linear gap is
    # negligible for the four small CVs (12.7, 16.2, 29.2%) and becomes
    # noticeable only for the largest CV (54.4%, HAE-patient baseline), where
    # the linear form (CV/100)^2 = 0.2959 is what the paper meant by '54.4%
    # CV' under the footnote convention.
    etalvc        ~ 0.02624    # 16.2% CV; Farrell 2013 Table 2 'Interindividual variability in volume'
    etalvmax      ~ 0.08526    # 29.2% CV; Farrell 2013 Table 2 'Interindividual variability in Vmax'
    etalrbase_hv  ~ 0.01613    # 12.7% CV; Farrell 2013 Table 2 'Interindividual variability in healthy volunteer baseline levels'
    etalrbase_hae ~ 0.29594    # 54.4% CV; Farrell 2013 Table 2 'Interindividual variability in HAE patient baseline levels'

    # ---- Residual error ----
    # Farrell 2013 Table 2 reports three proportional residual errors: 10.5%
    # for Studies C1 1101-01 / 1106-02 / 1202-01 / 1203-01 (treated as a
    # single common value in the final model), 23.6% for Study C1 1205-01,
    # and 53.6% for Study C1 1304-01. This file registers the 10.5% value as
    # the default propSd because it applies to four of the six studies and
    # to all of the assay-laboratory configurations except the two pivotal
    # phase III studies. Study-specific residual error is not encoded as a
    # covariate effect; the deviation is documented in the validation
    # vignette's 'Assumptions and deviations' section.
    propSd <- 0.105;  label("Proportional residual error on Cc (fraction; 'all other studies' value)")  # Farrell 2013 Table 2: propRV = 10.5% for Studies 1101-01, 1106-02, 1202-01, 1203-01
    addSd  <- 0.0556; label("Additive residual standard deviation on Cc (U/mL)")                        # Farrell 2013 Table 2: addRV = 0.0556 U/mL
  })

  model({
    # ---- Individual PK parameters ----
    # Each cohort-specific baseline is built as a clean exp(lrbase_X +
    # etalrbase_X) expression so the mu-referencing rules apply to both
    # etas; the binary DIS_HEALTHY indicator then routes each subject to
    # exactly one of the two exponentials (the other contributes zero).
    # Pattern is the mu-referenceable variant of the paired-cohort form in
    # Klunder_2017_upadacitinib.R.
    vc    <- exp(lvc   + etalvc)   * (WT / 70)^e_wt_vc
    vmax  <- exp(lvmax + etalvmax)
    km    <- exp(lkm)
    rbase_hv  <- exp(lrbase_hv  + etalrbase_hv)
    rbase_hae <- exp(lrbase_hae + etalrbase_hae)
    rbase     <- rbase_hv * DIS_HEALTHY + rbase_hae * (1 - DIS_HEALTHY)

    # ---- Endogenous-production rate derived from the steady-state assumption ----
    # Farrell 2013 Methods, p. 4: 'For nonlinear models, baseline
    # concentration was estimated, and the rate of production was a function
    # of the nonlinear clearance and baseline concentration.' At the no-dose
    # steady state the Michaelis-Menten elimination rate equals the
    # endogenous production rate, giving kin = Vmax * rbase / (Km + rbase)
    # (concentration-per-time units). Multiplying by vc converts to an
    # amount-per-time production rate so the ODE on central (in U) balances.
    kin_amt <- vmax * rbase / (km + rbase) * vc

    # ---- ODE system (one-compartment, MM elimination, endogenous input) ----
    # `central` carries the total functional C1INH amount in plasma
    # (endogenous + exogenous, in U). The initial condition is the
    # cohort-specific endogenous steady state (rbase * vc). Exogenous IV
    # doses add to `central` directly via dosing events (cmt = central).
    Cc <- central / vc
    d/dt(central) <- kin_amt - vmax * Cc / (km + Cc) * vc
    central(0)    <- rbase * vc

    # ---- Observation and residual error ----
    Cc ~ add(addSd) + prop(propSd)
  })
}
