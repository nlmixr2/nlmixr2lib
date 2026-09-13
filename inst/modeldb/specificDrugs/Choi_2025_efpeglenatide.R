Choi_2025_efpeglenatide <- function() {
  description <- "Two-compartment population PK model for subcutaneous efpeglenatide (HM11260C, a long-acting Fc-fusion GLP-1 receptor agonist) with dual parallel absorption into a single subcutaneous depot: a first-order bolus fraction plus a delayed fraction routed through a Savic 2007 transit-compartment chain, reproducing the double absorption peak. Pooled across one phase 1 and five phase 2 studies in adults with type 2 diabetes or non-diabetic obesity (Choi 2025). Body weight acts on ka and CL/F; disease status (T2DM vs obesity) acts on CL/F."
  reference <- paste(
    "Choi S, Seo J, Park S, Kim NY, Kim H, Lim H-S.",
    "Population pharmacokinetics of efpeglenatide in individuals with obesity and with type 2 diabetes.",
    "Front Pharmacol. 2025;16:1715585.",
    "doi:10.3389/fphar.2025.1715585.",
    sep = " "
  )
  vignette <- "Choi_2025_efpeglenatide"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot       = list(analyte = "efpeglenatide", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "efpeglenatide", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "efpeglenatide", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effect normalized to 92 kg on ka (exponent -0.927) and on CL/F",
        "(exponent 0.964) per Choi 2025 Table 3 and the Discussion. The 92 kg",
        "reference is the 'approximate median of the population' called for by",
        "Eq. 11; the pooled cohort median in Table 2 is 93.6 kg. Baseline (not",
        "time-varying) weight was used -- Choi 2025 lists the absence of",
        "longitudinal weight as a study limitation."
      ),
      source_name        = "WT"
    ),
    DIS_DIAB = list(
      description        = "Type 2 diabetes disease status",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-diabetic obesity, the HM-EXC-205 cohort)",
      notes              = paste(
        "1 = participant with type 2 diabetes mellitus (n = 293; studies",
        "HM-EXC-102, -201, -202, -203, -204); 0 = non-diabetic participant with",
        "obesity (n = 205; study HM-EXC-205) per Choi 2025 Table 2. Choi 2025",
        "screened 'obesity status'; because the pooled cohort is exactly",
        "partitioned into T2DM and non-diabetic obesity, the covariate is",
        "recorded on the canonical DIS_DIAB column with obesity as the",
        "reference level, matching the Overgaard 2019 semaglutide precedent.",
        "The canonical column does not separate Type 1 from Type 2; this",
        "cohort is entirely Type 2. Multiplicative effect 1.375 on CL/F",
        "(0.044 L/h in T2DM vs 0.032 L/h in obesity). Choi 2025 notes that the",
        "obesity data came from a single study, so disease status is fully",
        "confounded with study in this dataset and the effect cannot be",
        "separated from a study effect."
      ),
      source_name        = "Obesity status / disease status (T2DM vs obese)"
    )
  )

  # Covariates that Choi 2025 screened in the stepwise covariate-modeling
  # search (Methods 2.6) but did NOT retain in the final model. Documented so
  # the provenance of the covariate screen survives; not referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Age at baseline. Screened on all PK parameters in the SCM forward-inclusion step (Methods 2.6) and not retained in the final model (Table 3).",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Simulated age subgroups (<=34, 34-52, >=52 years) showed exposure GMRs within 0.80-1.25 of the reference (Figure 6), i.e. no clinically meaningful effect."
    ),
    SEXF = list(
      description        = "Sex. Screened on all PK parameters (Methods 2.6) and not retained in the final model (Table 3).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Cohort was 275/498 female (55.22%, Table 2). Simulated male/female exposure GMRs fell within 0.80-1.25 (Figure 6)."
    ),
    BMI = list(
      description        = "Body mass index at baseline. Screened on all PK parameters (Methods 2.6) and not retained; body weight was the retained size descriptor.",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Pooled median 33.4 kg/m^2 (range 19.2-57.7, Table 2)."
    ),
    LBM = list(
      description        = "Lean body mass. Screened on all PK parameters (Methods 2.6) and not retained; total body weight was the retained size descriptor.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Choi 2025 Table 2 reports this as 'LBW' (lean body weight), a registered alias of the canonical LBM column; pooled median 55.1 kg (range 34.7-89.4). The body-composition formula is not stated in the paper."
    ),
    RACE_WHITE = list(
      description        = "White / Caucasian race indicator. Race was screened as a five-level categorical covariate on all PK parameters (Methods 2.6) and not retained in the final model (Table 3).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not Caucasian)",
      notes              = "Choi 2025 Table 2: 398/498 Caucasian (79.92%). Simulated race-subgroup exposure GMRs fell within 0.80-1.25 (Figure 6)."
    ),
    RACE_BLACK = list(
      description        = "Black / African American race indicator. Part of the same five-level race screen; not retained.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not Black)",
      notes              = "Choi 2025 Table 2: 65/498 Black (13.05%)."
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator. Part of the same five-level race screen; not retained.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not Asian)",
      notes              = "Choi 2025 Table 2: 23/498 Asian (4.62%)."
    ),
    RACE_OTHER = list(
      description        = "Residual race-category indicator. Part of the same five-level race screen; not retained.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not in the residual category)",
      notes              = paste(
        "Choi 2025 Table 2 reports two residual levels: 'Others' (13/498) and",
        "'Native Hawaiian or Pacific Islander' (2/498). Both are folded onto",
        "this single canonical column because no NHPI canonical exists and the",
        "whole race covariate was rejected by the SCM -- registering a new",
        "canonical for a screened-and-discarded two-subject stratum would add",
        "an entry no model references. See the vignette Errata."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 498L,                                   # Choi 2025 Table 2 Total column
    n_studies      = 6L,                                     # Choi 2025 Table 1: HM-EXC-102 (phase 1) and -201/-202/-203/-204/-205 (phase 2)
    n_observations = 3316L,                                  # Choi 2025 Methods 2.2: post-dose concentrations above LLOQ used in estimation (3,596 collected)
    age_range      = "18-64 years",                          # Choi 2025 Table 2 Total column; see notes -- per-study ranges reach 75 years
    age_median     = "52 years",                             # Choi 2025 Table 2 Total column
    weight_range   = "49.7-191 kg",                          # Choi 2025 Table 2 Total column
    weight_median  = "93.6 kg",                              # Choi 2025 Table 2 Total column
    bmi_range      = "19.2-57.7 kg/m^2",                     # Choi 2025 Table 2 Total column
    bmi_median     = "33.4 kg/m^2",                          # Choi 2025 Table 2 Total column
    height_range   = "142-201 cm",                           # Choi 2025 Table 2 Total column
    sex_female_pct = 55.22,                                  # Choi 2025 Table 2: 275/498 female
    race_ethnicity = c(Caucasian = 79.92, Black = 13.05, Asian = 4.62, `Native Hawaiian or Pacific Islander` = 0.4, Other = 1.92),  # Choi 2025 Table 2 Total column
    disease_state  = "Type 2 diabetes mellitus (293 participants, 58.84%) or non-diabetic obesity (205 participants, 41.16%)",  # Choi 2025 Table 2
    dose_range     = "Subcutaneous efpeglenatide. Single ascending doses 2-100 ug/kg; multiple doses 0.3-6 mg once weekly (QW), 6-8 mg every two weeks (Q2W), and 8-16 mg once monthly (QM)",  # Choi 2025 Table 1
    regions        = "Multicenter, international",           # Choi 2025 Table 1 (HM-EXC-203 described as multicenter, international)
    notes          = paste(
      "Pooled analysis of one phase 1 (HM-EXC-102) and five phase 2",
      "(HM-EXC-201, -202, -203, -204, -205) studies; liraglutide comparator",
      "arms in HM-EXC-102 and HM-EXC-203 were excluded (Methods 2.1).",
      "Demographics per Table 2. Non-diabetic obesity data came from",
      "HM-EXC-205 alone, so disease status is fully confounded with study.",
      "Two Table 2 internal inconsistencies are reproduced as printed and",
      "flagged in the vignette Errata: the Total age range (18-64 years) is",
      "narrower than the HM-EXC-201 per-study range (47-75 years), and the",
      "race counts sum to 501 rather than 498. BLQ handling: 117 pre-dose",
      "(3.25%) and 153 post-dose (4.0%) samples were BLQ and treated as",
      "missing; all pre-first-dose records were excluded."
    )
  )

  # Implementation notes (see the vignette 'Assumptions and deviations /
  # Errata' section for the full justification):
  #
  # * Dual absorption into ONE depot. Choi 2025 Eq. 4 splits each
  #   administered dose between (i) a bolus fraction BIOA1 deposited
  #   directly into the subcutaneous depot Asc, and (ii) a fraction
  #   BIOA2 = 1 - BIOA1 routed through a Savic 2007 transit chain that
  #   discharges into the SAME depot. Both fractions then absorb into
  #   central with the SAME first-order rate constant ka (Discussion:
  #   'both pathways were described using the same first-order rate
  #   constant (ka = 0.006 h-1)'). This differs from the sibling
  #   Lee_2015_sumatriptan dual-pathway model, which uses two depots with
  #   two distinct ka values; here a single `depot` state is correct.
  #
  # * Dose-record convention. The user supplies ONE dose record targeting
  #   `depot` carrying the full administered amount. `f(depot) <- ffo`
  #   deposits only the bolus fraction, while `transit(nn, mtt, 1 - ffo)`
  #   supplies the delayed fraction. rxode2's transit() reads podo() as the
  #   UNSCALED record amount (verified: with f(depot) = 0.248 and
  #   transit(..., 0.752), total mass delivered equals the record amt
  #   exactly), so the two arms sum to the whole dose with no double
  #   counting. transit() is compartment-aware and fires only because the
  #   dose event targets the same compartment whose d/dt() carries the call.
  #
  # * Transit input rate carries the leading Ktr. Choi 2025 Eq. 4 writes the
  #   transit term as Dose*BIOA2*(Ktr*t)^n*exp(-Ktr*t)/n!, which is an
  #   AMOUNT, not a rate, and integrates to Dose*BIOA2/Ktr rather than
  #   Dose*BIOA2. The paper's own Eqs. 1-2 supply the correct form: Eq. 1
  #   gives the outflow from the last transit compartment as Ktr*a_n and
  #   Eq. 2 gives a_n(t), so the input rate is Ktr*a_n(t) =
  #   Dose*BIOA2*Ktr*(Ktr*t)^n*exp(-Ktr*t)/n! -- exactly Savic 2007 and
  #   exactly what rxode2's transit() computes. The dropped Ktr in Eq. 4 is
  #   a transcription slip; encoding it would break mass balance by a
  #   factor of Ktr = 2.43. Ktr = (n+1)/MTT = 6.52/2.68 = 2.433 /h is
  #   computed internally by transit() from nn and mtt.
  #
  # * Exact gamma instead of Stirling. Choi 2025 Eq. 3 evaluates n! for
  #   non-integer n with Stirling's approximation
  #   n! ~ sqrt(2*pi)*n^(n+0.5)*exp(-n). rxode2's transit() uses
  #   lgamma(n+1) directly. At n = 5.52 Stirling gives log(n!) = 5.6834 vs
  #   the exact 5.6985, so the paper's kernel amplitude is 1.5% high and
  #   its transit arm delivers 0.763 rather than 0.752 of the dose. The
  #   exact form is used here because it conserves mass exactly (the
  #   vignette gates on this); 1.5% is an order of magnitude below the
  #   residual error. Same choice as Lee_2015_sumatriptan.
  #
  # * Table 3 definition-cell typo. The row 'Covariate effect (theta) of
  #   body weight on CL/F' carries the definition text
  #   '(V_C/F) x (WT/92)^theta', copied from the ka row's template. The row
  #   label itself, the Results ('body weight influenced both ka and CL/F
  #   ... exponents of -0.927 and 0.964, respectively') and the Discussion
  #   ('a significant covariate on both ka and CL/F ... normalized to
  #   92 kg; the exponents of -0.927 and 0.964, respectively') all place
  #   the 0.964 exponent on CL/F. Three statements against one cell: the
  #   exponent is on CL/F.
  #
  # * No IOV, no IIV correlations. 'IOV was evaluated but not found to be
  #   significant and was therefore not included in the final model' and
  #   'no statistically significant covariance identified between IIV
  #   terms' (Results 3.1), so omega is diagonal.

  ini({
    # ---- Structural parameters. All from Choi 2025 Table 3 ("Final
    # population pharmacokinetic parameter estimates of efpeglenatide").
    # Each typical value is log-transformed per the nlmixr2lib convention;
    # Choi 2025 Eq. 8 parameterises IIV as PTV*exp(eta), which is the same
    # exponential form.
    lka  <- log(0.006); label("Absorption rate constant from depot to central, at 92 kg (1/h)")  # Table 3 ka = 0.006 1/h (RSE 4.33%)
    lcl  <- log(0.032); label("Apparent clearance CL/F in non-diabetic obesity, at 92 kg (L/h)")  # Table 3 'CL/F in obesity' = 0.032 L/h (RSE 1.79%); reference level of the disease-status covariate
    lvc  <- log(2.80);  label("Apparent central volume of distribution Vc/F (L)")                 # Table 3 Vc/F = 2.80 L (RSE 7.39%)
    lvp  <- log(3.96);  label("Apparent peripheral volume of distribution Vp/F (L)")              # Table 3 Vp/F = 3.96 L (RSE 6.57%); Vss = Vc + Vp = 6.76 L, cited as 'Vss ~ 6.8 L' in the Discussion
    lq   <- log(0.073); label("Apparent inter-compartmental clearance Q/F (L/h)")                 # Table 3 Q/F = 0.073 L/h (RSE 11.14%)

    # ---- Dual absorption pathway. Choi 2025 Eq. 7 holds the bolus
    # fraction on the logit scale so it stays in (0, 1) and can carry IIV:
    # BIOA1 = exp(BIOF)/(1 + exp(BIOF)), BIOA2 = 1 - BIOA1.
    logitffo <- -1.11;   label("Logit of BIOA1, the dose fraction absorbed via the first-order bolus pathway (logit units)")  # Table 3 BIOF = -1.11 (RSE 27.59%); expit(-1.11) = 0.2479, matching the Table 3 derived row 'BIO A1 = 0.248' and its footnote a. The complement BIOA2 = 0.752 (Table 3 footnote b) is the transit-chain fraction.
    lmtt     <- log(2.680); label("Mean transit time of the delayed-absorption transit chain (h)")                            # Table 3 MTT = 2.680 h (RSE 13.47%)
    lnn      <- log(5.520); label("Number of transit compartments N (continuous, unitless)")                                  # Table 3 N = 5.520 (RSE 22.28%); estimated as a continuous quantity via the analytical Savic 2007 input form, so a non-integer value is expected

    # ---- Covariate effects. Continuous covariates enter as a power
    # function normalized to the approximate population median (Eq. 11);
    # the categorical disease-status effect is the ratio of the two typical
    # CL/F values Choi 2025 reports (Eq. 12 switches between them, which is
    # numerically identical to a multiplier for a binary covariate).
    e_wt_ka   <- -0.927; label("Power exponent of (WT / 92 kg) on ka (unitless)")     # Table 3 'Covariate effect of body weight on ka' = -0.927 (RSE 14.24%); form (ka)*(WT/92)^theta
    e_wt_cl   <- 0.964;  label("Power exponent of (WT / 92 kg) on CL/F (unitless)")   # Table 3 'Covariate effect of body weight on CL/F' = 0.964 (RSE 6.75%); confirmed on CL/F by Results 3.1 and the Discussion -- the row's definition cell mistakenly prints V_C/F, see the implementation notes
    e_diab_cl <- 1.375;  label("Multiplier on CL/F for type 2 diabetes (1.375^DIS_DIAB)")  # Table 3: CL/F = 0.044 L/h in T2DM vs 0.032 L/h in obesity -> 0.044/0.032 = 1.375, the '38% higher clearance in T2DM' stated in the Discussion

    # ---- Inter-individual variability. Choi 2025 Table 3 reports IIV as
    # %CV (table footnote) computed with Eq. 9,
    # %CV = sqrt(exp(omega^2) - 1)*100, so the internal variance is
    # inverted as omega^2 = log(1 + CV^2). Table 3 gives no IIV on Q/F or
    # N, and Results 3.1 reports no significant covariance between IIV
    # terms, so omega is diagonal with those two etas omitted.
    etalka       ~ 0.045312  # log(1 + 0.2153^2) -- Table 3 ka IIV 21.53% CV (RSE 28.48%)
    etalcl       ~ 0.051592  # log(1 + 0.2301^2) -- Table 3 'CL/F in obesity' IIV 23.01% CV (RSE 5.60%); the single CL eta applies to both disease strata
    etalvc       ~ 0.382964  # log(1 + 0.6831^2) -- Table 3 Vc/F IIV 68.31% CV (RSE 18.28%)
    etalvp       ~ 0.060484  # log(1 + 0.2497^2) -- Table 3 Vp/F IIV 24.97% CV (RSE 27.11%)
    etalmtt      ~ 0.193996  # log(1 + 0.4627^2) -- Table 3 MTT IIV 46.27% CV (RSE 37.14% per Table 3; Results 3.1 prints 35.88% for the same term)
    etalogitffo  ~ 0.370157  # log(1 + 0.6693^2) -- Table 3 'BIO A1' IIV 66.93% CV (RSE 14.6%); additive on the logit scale, which is the stated purpose of the Eq. 7 transformation ('to allow implementation of interindividual variability')

    # ---- Residual error. Choi 2025 Eq. 10 is purely proportional,
    # Y = IPRED*(1 + eps_prop), and Table 3 reports eps_prop = 0.15.
    # Read as the proportional SD (15%), not as a variance: a variance of
    # 0.15 implies an SD of 38.7%, whose residual-only 2.5th percentile is
    # 1 - 1.96*0.387 = 0.24 of the prediction, already far below the lower
    # 95% prediction-interval edge of the paper's own pcVPC (Figure 3,
    # panel C: lower edge / median ~ 0.53) BEFORE any IIV is added. The
    # vignette reproduces that comparison quantitatively.
    propSd <- 0.15; label("Proportional residual error SD (fraction)")  # Table 3 'eps prop' = 0.15 (RSE 0.48%)
  })

  model({
    # ---- Individual parameters. Covariate forms per Choi 2025 Eq. 11
    # (continuous, power, median-normalized) and Eq. 12 (categorical).
    # Reference weight 92 kg per Table 3 and the Discussion.
    ka  <- exp(lka + etalka) * (WT / 92)^e_wt_ka
    cl  <- exp(lcl + etalcl) * (WT / 92)^e_wt_cl * e_diab_cl^DIS_DIAB
    vc  <- exp(lvc + etalvc)
    vp  <- exp(lvp + etalvp)
    q   <- exp(lq)
    mtt <- exp(lmtt + etalmtt)
    nn  <- exp(lnn)

    # Bolus (first-order) dose fraction BIOA1; the transit chain carries
    # the complement BIOA2 = 1 - ffo (Choi 2025 Eq. 7).
    ffo <- expit(logitffo + etalogitffo)

    # ---- Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- ODE system (Choi 2025 Eqs. 4-6, Figure 1).
    # depot: the subcutaneous depot Asc. Receives the bolus fraction
    #   directly from the dose record (scaled by f(depot) below) plus the
    #   Savic 2007 analytical transit input for the delayed fraction, and
    #   drains into central at the shared rate ka.
    # central / peripheral1: two-compartment linear disposition with
    #   first-order elimination.
    d/dt(depot)       <- transit(nn, mtt, 1 - ffo) - ka * depot
    d/dt(central)     <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ---- Dose routing. One dose record targets `depot` with the full
    # administered amount; f(depot) admits only the bolus fraction as a
    # bolus, and transit() supplies (1 - ffo) of the same record amount
    # through the chain. The two arms sum to the whole dose.
    f(depot) <- ffo

    # ---- Observation. central is in mg and vc in L, giving mg/L; x1000
    # converts to the ng/mL of Choi 2025 Figures 2-3 and Methods 2.2.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
