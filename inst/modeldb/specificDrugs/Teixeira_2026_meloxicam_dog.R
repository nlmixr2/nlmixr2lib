Teixeira_2026_meloxicam_dog <- function() {
  description <- paste(
    "Preclinical (dog). Two-compartment population PK model for meloxicam in 18",
    "female mixed-breed dogs undergoing elective ovariohysterectomy, after a single",
    "0.2 mg/kg oral dose of either a free meloxicam solution or meloxicam-loaded",
    "poly(epsilon-caprolactone) polymeric nanocapsules (Teixeira 2026). Absorption is a",
    "double extravascular process: a fraction F1 of the dose enters depot and is",
    "absorbed with a slow first-order rate ka1, while the remaining 1 - F1 enters depot2",
    "and is absorbed with a much faster first-order rate ka2 after a lag time Tlag2 -",
    "the structure the authors used to describe the secondary plasma peak seen in",
    "several animals. The nanocapsule formulation is a covariate on Tlag2 (2.1-fold",
    "longer) and on the apparent peripheral volume V2 (3.0-fold larger); clearance,",
    "central volume and both absorption rate constants are unaffected.",
    sep = " "
  )
  reference <- paste(
    "Teixeira FEG, Lock GA, Giacomeli R, Pacheco CO, Maciel TR, Funghetto-Ribeiro AP,",
    "Lugoch G, Beckmann DV, de Oliveira MT, Haas SE.",
    "Evaluating orally administered meloxicam-loaded polymeric nanocapsules in female",
    "dogs: a population pharmacokinetic modeling study.",
    "Pharmaceuticals (Basel). 2026;19(3):412. doi:10.3390/ph19030412.",
    sep = " "
  )
  vignette <- "Teixeira_2026_meloxicam_dog"

  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Confirmed against Teixeira 2026 Figure 1 (final model
  # structure diagram) and Section 4.4 (oral dosing, jugular plasma sampling).
  compartmentData <- list(
    depot       = list(analyte = "meloxicam", units = "mg", specimen = "administration site", verified = TRUE),
    depot2      = list(analyte = "meloxicam", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "meloxicam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "meloxicam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per dog. Weight is NOT a covariate effect in Teixeira 2026 - it was",
        "screened (supplementary Table S2 footnote b) and not retained. It appears here",
        "purely as a unit conversion: Table 2 reports CL_pop in mL/min/kg and V1_pop /",
        "V2_pop in L/kg, i.e. already normalised to body weight, while Q_pop is reported",
        "in absolute L/h. Multiplying the per-kg clearance and volumes by WT puts every",
        "disposition parameter on the same absolute (L, L/h) scale as Q so the",
        "micro-constants Q/V1 and Q/V2 are dimensionally correct. Dose the model in mg",
        "(i.e. 0.2 * WT for the published regimen). Study means: 14.2 kg (free MLX arm)",
        "and 12.9 kg (NC-MLX arm); overall range 10.5-16.6 kg."
      ),
      source_name        = "Body Weight"
    ),
    FORM_MLX_NANOCAPSULE = list(
      description        = "Meloxicam-loaded polymeric nanocapsule formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = free meloxicam solution (MLX in 60% v/v PEG-400 diluted with water)",
      notes              = paste(
        "1 = NC-MLX, the poly(epsilon-caprolactone) nanocapsule formulation prepared by",
        "interfacial polymer deposition (mean diameter 326 +/- 13 nm, zeta potential",
        "-26.2 +/- 6.4 mV, SPAN 1.10, drug content 99.47%, encapsulation efficiency",
        "99.71%; Section 2.1). Time-fixed per dog: parallel-arm design, 9 dogs per arm.",
        "Both arms received the same 0.2 mg/kg meloxicam dose, so the indicator marks the",
        "delivery system and not the amount administered. Retained on Tlag2 and on V2",
        "only (Section 2.3 and Table 2); it was not retained on CL, V1, Q, ka1, ka2 or F1."
      ),
      source_name        = "NC-MLX"
    )
  )

  # Screened by forward inclusion / backward elimination (Section 4.6.2 and
  # supplementary Table S2 footnote a) but not retained in the final model, so
  # no coefficient is reported and AGE is never referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at study entry",
      units       = "months",
      type        = "continuous",
      notes       = paste(
        "Tested as a covariate on the one-compartment and double-extravascular-absorption",
        "candidate models (Table S2 footnote a) but absent from the final model in Table 2.",
        "Cohort: 9-48 months overall; 26.1 +/- 18.3 months (free MLX) and 19.5 +/- 14.4",
        "months (NC-MLX) per supplementary Table S1. Recorded here to preserve the",
        "covariate screen; no point estimate is published."
      )
    )
  )

  population <- list(
    species        = "dog (female, mixed breed)",
    n_subjects     = 18L,
    n_studies      = 1L,
    age_range      = "9-48 months",
    age_median     = "26.1 +/- 18.3 months (free MLX arm); 19.5 +/- 14.4 months (NC-MLX arm), mean +/- SD (Table S1)",
    weight_range   = "10.5-16.6 kg",
    weight_median  = "14.2 +/- 1.6 kg (free MLX arm); 12.9 +/- 1.8 kg (NC-MLX arm), mean +/- SD (Table S1)",
    sex_female_pct = 100,
    race_ethnicity = "Not applicable (mixed-breed dogs; the authors note that breed heterogeneity is itself a source of metabolic variability, Discussion)",
    disease_state  = "Healthy female dogs undergoing elective ovariohysterectomy; drug given 4 h before surgery",
    dose_range     = paste(
      "Single oral dose of 0.2 mg/kg meloxicam, given 4 h before ovariohysterectomy.",
      "NC-MLX dogs received 2.8 +/- 0.3 mL of the nanocapsule suspension (1 mg/mL);",
      "free-MLX dogs received 5.2 +/- 0.4 mL of the PEG-400 solution (0.5 mg/mL)."
    ),
    regions        = "Brazil (UNIPAMPA Veterinary Hospital, Uruguaiana, Rio Grande do Sul; single centre)",
    n_observations = 196L,
    notes          = paste(
      "196 plasma meloxicam observations from 18 dogs (9 per arm), all above the",
      "quantification limit. Sampling at 0.5, 1, 2, 4, 6, 8, 12, 24, 36, 48 and 60 h",
      "post-dose; plasma meloxicam by HPLC-PDA with piroxicam internal standard.",
      "Fit in Monolix 2021R2 by SAEM; model selection assisted by Sycomore 2021R2;",
      "nonparametric bootstrap (1000 replicates) via Rsmlx. Internal evaluation by",
      "goodness-of-fit, NPDE and a 1000-replicate VPC; external evaluation against two",
      "digitised literature datasets (28 and 66 mean-profile observations from five",
      "published canine oral-meloxicam studies) gave MPE 1.01% / 0.63% and RMSE 7.04 /",
      "6.05. Anaesthesia lasted 100 +/- 16 min (free MLX) and 109 +/- 31 min (NC-MLX);",
      "surgery 55 +/- 21 and 60 +/- 22 min respectively (Table S1)."
    )
  )

  ini({
    # =================================================================
    # Structural parameters - Teixeira 2026 Table 2 ("Estimation of the
    # PopPK model parameters"), Estimate (RSE%) column.
    #
    # UNITS NOTE. The units column of Table 2 is internally inconsistent:
    # the rows Ka2_pop, F1_pop and Tlag2_pop carry a cyclic one-row shift
    # of their unit labels (Ka2 is labelled "(h)", F1 is labelled
    # "(h-1)", Tlag2 is left blank). The table footnote defines Ka2 as a
    # "first-order absorption rate constant" and F1 as a "fraction
    # absorbed", and Section 2.3 states "Tlag2_pop_Free MLX = 1.22 h", so
    # the correct assignment is ka1 [1/h], ka2 [1/h], F1 [unitless],
    # Tlag2 [h]. Applied here; see the vignette Errata.
    # =================================================================

    # ---- Absorption (double extravascular, Figure 1) ----------------
    lka   <- log(0.086); label("First-order absorption rate constant ka1 of the first (unlagged) absorption process (1/h)")  # Table 2: Ka1_pop = 0.086 (RSE 20.1%); bootstrap median 0.099 [0.049, 0.397]
    lka2  <- log(1.82);  label("First-order absorption rate constant ka2 of the second (lagged) absorption process (1/h)")   # Table 2: Ka2_pop = 1.82 (RSE 51.3%); bootstrap median 2.722 [0.626, 16.602]

    # F1 is a fraction of the dose and is therefore carried on the logit
    # scale so that F1 and 1 - F1 both stay in (0, 1) for every
    # individual. Section 4.6.2 describes IIV generically as exponential
    # (Equation 1), but an exponential IIV of omega = 0.67 on a typical
    # value of 0.85 would put roughly 40% of dogs above F1 = 1 and give
    # them a negative dose in depot2, which is not simulatable. The
    # logit-normal form is Monolix's standard distribution for a bounded
    # fraction and is the only physically valid reading; documented in
    # the vignette Errata. Precedent: Bjornsson_2023_buprenorphine.R.
    logitfdepot <- log(0.85 / (1 - 0.85)); label("Logit of the fraction F1 of the dose entering the first absorption depot (unitless)")  # Table 2: F1_pop = 0.85 (RSE 5.87%); logit(0.85) = 1.7346; bootstrap median 0.81 [0.652, 0.944]
    ltlag2      <- log(1.22);              label("Lag time Tlag2 of the second absorption process, free-MLX reference (h)") # Table 2: Tlag2_pop = 1.22 (RSE 16.3%); Section 2.3 "Tlag2_pop_Free MLX = 1.22 h"

    # ---- Disposition ------------------------------------------------
    # CL, V1 and V2 are reported normalised to body weight; Q is
    # reported in absolute L/h. See covariateData$WT$notes.
    lcl <- log(0.1 * 60 / 1000); label("Apparent clearance CL/F per kg body weight (L/h/kg)")                        # Table 2: CL_pop = 0.1 mL/min (RSE 8.56%) = 0.1 mL/min/kg per Table 1 CL/F; 0.1 * 60 / 1000 = 0.006 L/h/kg
    lvc <- log(0.049);           label("Apparent central volume V1/F per kg body weight (L/kg)")                     # Table 2: V1_pop = 0.049 L/kg (RSE 56.0%); bootstrap median 0.151 [0.023, 0.493]
    lq  <- log(0.24);            label("Intercompartmental clearance Q (L/h, absolute - not weight-normalised)")     # Table 2: Q_pop = 0.24 L/h (RSE 24.1%); bootstrap median 0.271 [0.09, 1.046]
    lvp <- log(0.134);           label("Apparent peripheral volume V2/F per kg body weight, free-MLX reference (L/kg)") # Table 2: V2_pop = 0.134 L/kg (RSE 27.8%); bootstrap median 1.394 [0.296, 2.207]

    # ---- Covariate effects (formulation) ----------------------------
    # Monolix categorical-covariate coefficients on log-normally
    # distributed parameters enter additively on the log scale. Both are
    # confirmed by back-substitution against the values the paper prints
    # in prose: 1.22 * exp(0.74) = 2.56 h (paper: 2.55 h) and
    # 0.134 * exp(1.11) = 0.4066 L/kg (paper Section 2.3 / Discussion:
    # 0.406 L/kg).
    e_form_nanocapsule_tlag2 <- 0.74; label("Additive shift on log(Tlag2) for the NC-MLX nanocapsule formulation vs free MLX (unitless)") # Table 2: beta_Tlag2_NC-MLX = 0.74 (RSE 28.7%); bootstrap median 0.926 [-0.782, 7.591]
    e_form_nanocapsule_vp    <- 1.11; label("Additive shift on log(V2) for the NC-MLX nanocapsule formulation vs free MLX (unitless)")    # Table 2: beta_V2_NC-MLX = 1.11 (RSE 31.4%); bootstrap median 1.081 [0.229, 2.836]

    # =================================================================
    # Inter-individual variability - Teixeira 2026 Table 2, "IIV" block.
    # The tabulated omega values are standard deviations of the random
    # effect on its own transformed scale (Monolix reports lower-case
    # omega_<param> as an SD; the "(%)" column header and the "Omega:
    # variance" footnote are boilerplate that contradict each other).
    # Two independent checks support the SD reading: (i) five of the
    # seven RSEs are below sqrt(2/N) = 33% for N = 18, which an
    # asymptotic variance estimate cannot achieve; (ii) reading omega_F1
    # = 0.67 as a logit-scale SD reproduces the F1 bootstrap interval
    # almost exactly - expit(logit(0.85) +/- 1.645 * 0.67) = [0.653,
    # 0.945] vs the tabulated [0.652, 0.944]. nlmixr2 takes variances,
    # so each entry below is (paper omega)^2. Q carries no IIV
    # (Section 2.3: "Interindividual variability was maintained for all
    # PK parameters except Q"); Table 2 leaves its shrinkage cell blank.
    # No IIV correlation structure is reported.
    # =================================================================
    etalka         ~ 0.1444  # Table 2: omega_Ka1   = 0.38 (RSE 21.1%) -> 0.38^2
    etalka2        ~ 1.6384  # Table 2: omega_Ka2   = 1.28 (RSE 35.3%) -> 1.28^2
    etalogitfdepot ~ 0.4489  # Table 2: omega_F1    = 0.67 (RSE 53.8%) -> 0.67^2 (logit scale)
    etaltlag2      ~ 0.1444  # Table 2: omega_Tlag2 = 0.38 (RSE 20.7%) -> 0.38^2
    etalcl         ~ 0.1024  # Table 2: omega_CL    = 0.32 (RSE 19.4%) -> 0.32^2
    etalvc         ~ 2.0736  # Table 2: omega_V1    = 1.44 (RSE 29.8%) -> 1.44^2
    etalvp         ~ 0.4225  # Table 2: omega_V2    = 0.65 (RSE 20.3%) -> 0.65^2

    # =================================================================
    # Residual error - Teixeira 2026 Table 2 "Residual variability".
    # Section 2.3: "The residual error was described as proportional."
    # Monolix's proportional model is y = f * (1 + b * eps) with
    # b = 0.19, which is nlmixr2's prop(propSd) with propSd = 0.19.
    # =================================================================
    propSd <- 0.19; label("Proportional residual SD on plasma meloxicam concentration (fraction)")  # Table 2: proportional model = 0.19 (RSE 6.83%); bootstrap median 0.192 [0.154, 0.24]
  })

  model({
    # ---- Individual parameters --------------------------------------
    # CL, V1 and V2 are per-kg quantities in Table 2 and are multiplied
    # by WT to reach the absolute (L/h, L) scale on which Q is reported.
    ka    <- exp(lka  + etalka)
    ka2   <- exp(lka2 + etalka2)
    fdepot <- expit(logitfdepot + etalogitfdepot)
    tlag2 <- exp(ltlag2 + e_form_nanocapsule_tlag2 * FORM_MLX_NANOCAPSULE + etaltlag2)
    cl    <- exp(lcl + etalcl) * WT
    vc    <- exp(lvc + etalvc) * WT
    vp    <- exp(lvp + e_form_nanocapsule_vp * FORM_MLX_NANOCAPSULE + etalvp) * WT
    q     <- exp(lq)

    # ---- Concentration ----------------------------------------------
    # central in mg, vc in L -> mg/L, which equals ug/mL, the unit of
    # Teixeira 2026 Table 1 (Cmax in ug/mL, AUC in ug h/mL).
    Cc <- central / vc

    # ---- ODE system (Figure 1) --------------------------------------
    # The dose is split between two parallel first-order absorption
    # routes feeding the same central compartment: fraction F1 into
    # depot (rate ka1, no lag) and 1 - F1 into depot2 (rate ka2, lag
    # Tlag2). Central exchanges with a single peripheral compartment via
    # Q and is eliminated linearly by CL.
    d/dt(depot)       <- -ka  * depot
    d/dt(depot2)      <- -ka2 * depot2
    d/dt(central)     <-  ka  * depot + ka2 * depot2 -
                          cl * Cc -
                          q  * (Cc - peripheral1 / vp)
    d/dt(peripheral1) <-  q  * (Cc - peripheral1 / vp)

    # ---- Dose split and absorption lag ------------------------------
    # Dose BOTH depots with the same total amount; f() partitions it so
    # the absorbed mass sums to the administered dose.
    f(depot)     <- fdepot
    f(depot2)    <- 1 - fdepot
    alag(depot2) <- tlag2

    # ---- Residual error ---------------------------------------------
    Cc ~ prop(propSd)
  })
}
