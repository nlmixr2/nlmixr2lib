Ahmed_2015_topiramate <- function() {
  description <- "Population PK/PD model of topiramate (TPM) and its acute effect on phonemic generative fluency (Controlled Oral Word Association, COWA) in healthy adult volunteers given single oral or intravenous doses of 50-100 mg (Ahmed 2015). Two-compartment popPK with first-order absorption and elimination, oral bioavailability ~108%, allometric body-weight scaling on CL/Q (fixed 3/4) and Vc/Vp (fixed 1); separate proportional residual errors for oral and IV cohorts. PD: COWA = baseline * practice_factor * exp(-KE * Cc), where the practice factor inflates baseline by 12% beginning with the fourth (and subsequent) COWA test administration and KE = 0.157 L/mg gives a 14.5% drop in COWA per 1 mg/L rise in plasma TPM."
  reference <- "Ahmed GF, Marino SE, Brundage RC, Pakhomov SVS, Leppik IE, Cloyd JC, Clark A, Birnbaum AK. Pharmacokinetic-pharmacodynamic modelling of intravenous and oral topiramate and its effect on phonemic fluency in adult healthy volunteers. Br J Clin Pharmacol. 2015;79(5):820-830. doi:10.1111/bcp.12556"
  paper_specific_residual_sds <- c("propSdOral", "propSdIv")
  vignette <- "Ahmed_2015_topiramate"
  units <- list(
    time          = "hour",
    dosing        = "mg",
    concentration = "mg/L (TPM plasma); words per three 60-second trials (COWA)"
  )

  covariateData <- list(
    WT = list(
      description        = "Actual body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on CL and Q (exponent 0.75, fixed) and on Vc and Vp (exponent 1, fixed) centred at 70 kg (Ahmed 2015 Table 2 footer equations). Studied range 54.73-112.30 kg.",
      source_name        = "WT"
    ),
    OCC = list(
      description        = "Cumulative count of COWA test administrations for the subject by the current observation row (1 for the first administration, 2 for the second, ...). Used by the PD model as a derived binary practice-effect threshold (OCC >= 4 triggers the baseline-inflation factor); PK observations may carry any OCC value as the practice-effect term only enters the COWA observation equation.",
      units              = "(count)",
      type               = "count",
      reference_category = NULL,
      notes              = "Time-varying within subject. Ahmed 2015 reports the 12% baseline inflation as a covariate effect on COWA baseline keyed to NCOWA >= 4 (page 5 / Table 3, parameter theta_NCOWA>=4); the canonical OCC integer covariate is reused here for the COWA test-administration count.",
      source_name        = "NCOWA"
    ),
    ROUTE_IV = list(
      description        = "Indicator for intravenous (IV) administration of the stable-labelled topiramate formulation",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (oral tablet)",
      notes              = "Per-observation dosing-route indicator (1 = IV infusion, 0 = oral tablet). Ahmed 2015 retained no route covariate effect on the structural PK parameters (Results page 5: 'age, sex, race and TPM formulation had insignificant effects on CL') but fit a separate proportional residual error magnitude for each formulation (Table 2: oral %CV = 18.4, IV %CV = 7.2). The model body selects between propSdOral and propSdIv via ROUTE_IV. Distinct from the rxode2 cmt event column (cmt = central for IV doses, cmt = depot for oral doses).",
      source_name        = "FLAG"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 32L,
    n_studies      = 3L,
    age_range      = "19-55 years",
    age_median     = "26.5 years",
    weight_range   = "54.73-112.30 kg",
    weight_median  = "77.27 kg",
    sex_female_pct = 37.5,
    race_ethnicity = c(Caucasian = 75.0, AfricanAmerican = 15.625, Other = 6.25, Unknown = 3.125),
    disease_state  = "Healthy adult volunteers; normal renal function; no medications known to interact with TPM or alter cognitive function",
    dose_range     = "50-100 mg single oral dose; 50-100 mg single IV infusion over 15 min (stable-labelled formulation)",
    regions        = "USA (University of Minnesota; University of Florida)",
    notes          = "Pooled across three randomized crossover studies (Ahmed 2015 Table 1). Study I: n=12, two PK / PD visits with rich PK sampling (15 timepoints over 120 h, three COWA timepoints 0.25/2.5/6 h post-dose). Study II: n=11 (placebo-controlled), one sparse PK sample and one COWA test 2-3 h post-dose. Study III: n=9, single sparse PK / single COWA timepoint (third-period 2 mg lorazepam arm excluded from the modelling dataset). External validation cohort (n=9, single 200 mg oral dose) referenced but not used for model fitting."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural PK parameters at the 70 kg reference weight, full
    # maturity (healthy adults only). All values from Ahmed 2015 Table 2.
    # ------------------------------------------------------------------
    lcl     <- log(1.21);  label("Topiramate plasma clearance CL at 70 kg reference weight (L/h)")          # Ahmed 2015 Table 2 (CL = 1.21 * (WT/70)^0.75)
    lvc     <- log(59.3);  label("Central volume of distribution Vc at 70 kg reference weight (L)")         # Ahmed 2015 Table 2 (Vc = 59.3 * (WT/70)^1.0)
    lq      <- log(1.02);  label("Intercompartmental clearance Q at 70 kg reference weight (L/h)")          # Ahmed 2015 Table 2 (Q  = 1.02 * (WT/70)^0.75)
    lvp     <- log(12.1);  label("Peripheral volume of distribution Vp at 70 kg reference weight (L)")      # Ahmed 2015 Table 2 (Vp = 12.1 * (WT/70)^1.0)
    lka     <- log(2.38);  label("First-order absorption rate constant Ka (1/h)")                            # Ahmed 2015 Table 2 (Ka = 2.38)
    lfdepot <- log(1.08);  label("Oral bioavailability F (fraction)")                                        # Ahmed 2015 Table 2 (F = 1.08; ~108% per Discussion page 7)

    # Allometric exponents (fixed in the final model; CL exponent fixed
    # at the theoretical 0.75 after the 95% CI [-0.16, 0.99] of the
    # data-estimated value included 0.75; volume exponents fixed at 1).
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent on body weight for CL and Q (unitless)")          # Ahmed 2015 Table 2 footer / Results page 5
    e_wt_vc_vp <- fixed(1.0);  label("Allometric exponent on body weight for Vc and Vp (unitless)")         # Ahmed 2015 Table 2 footer / Results page 5

    # ------------------------------------------------------------------
    # PD structural parameters (continuous-data analysis column of
    # Ahmed 2015 Table 3). The exponential decline link is
    #   COWA_ij = BL_i * practice_factor * exp(-KE * Cc) + eps_ij,
    # where practice_factor = 1 if NCOWA < 4 (typical pretreatment /
    # early in-trial baseline) and = theta_NCOWA>=4 = 1.12 once the
    # subject has been administered the COWA test four or more times.
    # ------------------------------------------------------------------
    lrbase     <- log(42.5);   label("Typical baseline COWA score BL (words per three 60-s trials)")        # Ahmed 2015 Table 3 (BL = 42.5)
    lke_cowa   <- log(0.157);  label("Exponential decline constant KE on COWA per plasma TPM (L/mg)")       # Ahmed 2015 Table 3 (KE = 0.157)
    e_practice_cowa <- 0.12;   label("Fractional increase in baseline COWA after >= 4 test administrations") # Ahmed 2015 Table 3 (theta_NCOWA>=4 = 1.12, encoded here as the fractional shift +0.12)

    # ------------------------------------------------------------------
    # Inter-individual variability (exponential model, NONMEM
    # OMEGA = log(CV^2 + 1)). Reported as %CV in Tables 2 / 3.
    # Eta on Q, Vp, F, and KE were fixed to zero in the published
    # model because either the estimate of the variance was
    # unrealistically small (F) or the model failed to converge with
    # the variance estimated (Q, Vp, KE); they are omitted entirely
    # here, matching the published OMEGA structure (diagonal, four
    # estimated etas: lcl, lvc, lka, lrbase).
    # ------------------------------------------------------------------
    etalcl    ~ 0.036572  # CV = 19.3% -> log(1 + 0.193^2); Ahmed 2015 Table 2
    etalvc    ~ 0.058291  # CV = 24.5% -> log(1 + 0.245^2); Ahmed 2015 Table 2
    etalka    ~ 0.250197  # CV = 53.3% -> log(1 + 0.533^2); Ahmed 2015 Table 2
    etalrbase ~ 0.028493  # CV = 17.0% -> log(1 + 0.170^2); Ahmed 2015 Table 3 IIV of BL

    # ------------------------------------------------------------------
    # Residual error. Per-formulation proportional residual errors on
    # plasma TPM (Ahmed 2015 Table 2): oral %CV = 18.4, IV %CV = 7.2;
    # the model body selects between them via the ROUTE_IV covariate.
    # PD: additive on the COWA score, reported as variance 7.1 words^2
    # per three trials -> SD ~ 2.665 words.
    # ------------------------------------------------------------------
    propSdOral <- 0.184;     label("Proportional residual error on plasma TPM, oral formulation (fraction)") # Ahmed 2015 Table 2 (oral TPM %CV = 18.4)
    propSdIv   <- 0.072;     label("Proportional residual error on plasma TPM, IV formulation (fraction)")   # Ahmed 2015 Table 2 (IV TPM %CV = 7.2)
    addSd_COWA <- sqrt(7.1); label("Additive residual error on COWA score (words per three trials)")         # Ahmed 2015 Table 3 (sigma_eps^2 = 7.1)
  })

  model({
    # 1. Individual parameters with allometric scaling on body weight
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp)          * (WT / 70)^e_wt_vc_vp
    ka <- exp(lka + etalka)
    fdepot <- exp(lfdepot)

    # 2. Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 3. ODE system. Dose enters depot (mg, oral) or central (mg, IV)
    # via the event table; central / peripheral1 carry mg.
    d/dt(depot)        <- -ka * depot
    d/dt(central)      <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1)  <-  k12 * central - k21 * peripheral1

    # 4. Bioavailability for the oral depot (IV doses bypass depot via
    # the event table by dosing directly into central, so this scaling
    # only affects oral-route input).
    f(depot) <- fdepot

    # 5. Plasma TPM concentration in mg/L (dose in mg, volume in L)
    Cc <- central / vc

    # 6. PD: exponential decline of COWA from a subject-specific baseline,
    # with a 12% baseline inflation once the subject has been
    # administered the COWA test four or more times (Ahmed 2015 Results
    # page 6 and Table 3). OCC is the canonical integer covariate for
    # the test-administration count; the threshold is OCC >= 4.
    rbase    <- exp(lrbase + etalrbase)
    ke_cowa  <- exp(lke_cowa)
    practice <- 1 + e_practice_cowa * (OCC >= 4)
    COWA     <- rbase * practice * exp(-ke_cowa * Cc)

    # 7. Route-conditional proportional residual error on plasma TPM
    # (Ahmed 2015 Table 2): propSd collapses to propSdOral when
    # ROUTE_IV = 0 and to propSdIv when ROUTE_IV = 1.
    propSd <- propSdOral + (propSdIv - propSdOral) * ROUTE_IV

    # 8. Observation models
    Cc   ~ prop(propSd)
    COWA ~ add(addSd_COWA)
  })
}
