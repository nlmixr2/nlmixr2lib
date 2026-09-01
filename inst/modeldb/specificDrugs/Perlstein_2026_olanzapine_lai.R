Perlstein_2026_olanzapine_lai <- function() {
  description <- "Two-compartment population PK model for olanzapine following subcutaneous TV-44749, an extended-release copolymer-based long-acting injectable, in adults with schizophrenia or schizoaffective disorder and in healthy participants (Perlstein 2026). Absorption is a convolution-based prescribed input function: a double Weibull release profile splitting the dose between a rapid first release process and a delayed sustained second release process, implemented as two parallel depot compartments each emptying with its own Weibull hazard. Allometric body weight on CL/F and V/F with estimated exponents, and a linear dose effect on the second-process release time. Carries the published Emax dopamine D2 receptor occupancy layer (Mamo 2008) as an algebraic observable."
  reference <- paste(
    "Perlstein I, Cherniakov I, Elgart A, Gomeni R, Gutman D, Merenlender Wagner A, Singh R.",
    "Population Pharmacokinetic Model-Based Dose Selection of Extended-Release Injectable Olanzapine",
    "(TV-44749) for Subcutaneous Use in Phase 3 Clinical Trial in Adults with Schizophrenia.",
    "The Journal of Clinical Pharmacology. 2026;66(1):e70144. doi:10.1002/jcph.70144.",
    "D2 receptor occupancy layer from Mamo D, Kapur S, Keshavan M, et al.",
    "D2 receptor occupancy of olanzapine pamoate depot using positron emission tomography:",
    "an open-label study in patients with schizophrenia.",
    "Neuropsychopharmacology. 2008;33(2):298-304.",
    sep = " "
  )
  vignette <- "Perlstein_2026_olanzapine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")
  # Declared explicitly because buildModelDb()'s fallback heuristic recognises
  # only the literal names "depot" and "central" and would report
  # `depot,central` -- omitting depot2, which every injection MUST also target
  # (the release fraction is applied as a bioavailability split across the two
  # parallel depots), and naming central, which this subcutaneous-only model
  # never doses.
  dosing <- c("depot", "depot2")

  compartmentData <- list(
    depot       = list(analyte = "olanzapine", units = "mg", specimen = "administration site", verified = TRUE),
    depot2      = list(analyte = "olanzapine", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "olanzapine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "olanzapine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling on CL/F and V/F with EXPONENTS ESTIMATED by the paper",
        "(Table 1: CL_P = 0.31, RSE 86.8%; V_P = 0.364, RSE 110.7%), not fixed at the",
        "standard 0.75 / 1. The Methods text on p. 6 describes a candidate model with the",
        "exponents fixed at 0.75 and 1, but that model is contradicted by Table 1, which is",
        "headed 'Final Model Parameter Estimates' and reports an RSE for each exponent -- a",
        "fixed exponent has no RSE. See the vignette Errata for the numeric test that settles",
        "it against Table 2. The paper prints NO reference (centring) weight anywhere, so the",
        "standard 70 kg is used; the cohort mean weight was 85.8 kg (Results, Data).",
        "Baseline value, not time-varying."
      ),
      source_name        = "WT"
    ),
    DOSE_TV44749_MG = list(
      description        = "Administered TV-44749 subcutaneous injection dose, as olanzapine milligrams",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-dose-record covariate carrying the milligram dose of the TV-44749 injection.",
        "Enters as a linear additive effect on the SECOND-process Weibull release TIME:",
        "TD1 = TD1_0 + TD1_1 * DOSE (Perlstein 2026 Table 1 footnote: 'TD_0 and TD1_0, the",
        "parameter values at the dose = 0; TD_1 and TD1_1, dose dependent changes').",
        "Because this file parameterises the release processes by their rate scalers ra / ra2",
        "(= 1 / release time), the effect is applied on the reciprocal:",
        "ra2 = 1 / (1 / ra2_0 + e_dose_ra2 * DOSE_TV44749_MG), which is algebraically identical",
        "to the paper's printed additive form and preserves the printed 0.0939 h/mg exactly.",
        "Doses studied: 70 and 105 mg (healthy, subtherapeutic), 283, 318, 425, 531 and 566 mg",
        "(patients). Must be supplied on every record, dose and observation alike, because it",
        "is read inside model()."
      ),
      source_name        = "DOSE"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in the stepwise forward/backward covariate search (Methods, Covariate Analysis) but NOT retained in the final model; no coefficient is reported."
    ),
    BMI = list(
      description = "Body mass index at baseline",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened as a covariate because it is well characterised for olanzapine in the oral prescribing information (Discussion), but NOT retained in the final model; no coefficient is reported."
    ),
    HT = list(
      description = "Height at baseline",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate search but NOT retained in the final model; no coefficient is reported."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened in the stepwise covariate search but NOT retained in the final model; no coefficient is reported. Cohort was 73% male (Results, Data)."
    ),
    RACE_BLACK = list(
      description = "Black or African American race indicator (1 = yes, 0 = no)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened in the stepwise covariate search but NOT retained in the final model; no coefficient is reported. Cohort was 79% Black or African American (Results, Data)."
    ),
    SMOKER = list(
      description = "Current smoker indicator (1 = smoker, 0 = non-smoker)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a covariate because smoking induces CYP1A2-mediated olanzapine clearance and is discussed in the oral prescribing information (Discussion), but NOT retained in the final TV-44749 model; no coefficient is reported. Cohort was 69% smokers (Results, Data)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 89,
    n_studies      = 1,
    age_range      = "18-60 years (healthy participants) and 18-65 years (patients); cohort mean 45.6 years (SD 10.3)",
    weight_mean    = "85.8 kg (SD 15.7)",
    bmi_mean       = "28.3 kg/m^2 (SD 4.7)",
    sex_female_pct = 27,
    race_ethnicity = c(Black = 79, White = 29),
    disease_state  = "Healthy participants, and adults with a DSM-5-confirmed diagnosis of schizophrenia or schizoaffective disorder on a clinically stable oral olanzapine regimen (no dose change in the preceding 4 weeks, no other concurrent antipsychotic).",
    dose_range     = "Single subcutaneous TV-44749 injections of 70 or 105 mg (healthy, subtherapeutic) or 318, 425 or 531 mg (patients); three consecutive once-monthly injections of 283 or 566 mg (patients, multiple-dose cohort).",
    regions        = "United States",
    notes          = paste(
      "Phase 1 study TV-44749-SAD-10154, open-label, single- and multiple-dose.",
      "3530 olanzapine plasma PK samples from 89 TV-44749 recipients entered the TV-44749",
      "model; one further participant was excluded after taking an oral olanzapine dose",
      "following the TV-44749 injection (Results, Data). Race percentages are quoted as",
      "printed in the paper and sum to more than 100. Full per-cohort demographics are in",
      "Supplemental Table 1. Estimation in NONMEM 7.4."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Double Weibull release input function (Perlstein 2026 Methods p. 4):
    #   r(t) = FF*exp(-(t/TD)^SS) + (1 - FF)*exp(-(t/TD1)^SS1)
    # r(t) is the fraction of the dose still UNRELEASED, so the input into
    # the central compartment is Dose*(-dr/dt). Written here with the
    # registered Weibull stems ra / gam1 (and their second-process partners
    # ra2 / gam2), for which the equivalent form is
    #   r(t) = frel*exp(-(ra*t)^gam1) + (1 - frel)*exp(-(ra2*t)^gam2)
    # so ra = 1/TD and ra2 = 1/TD1 exactly. The reciprocal is spelled out in
    # each comment below because it means the ini() value is not the number
    # printed in Table 1 (operator sidecar oare_PMC12775547 request-001,
    # answer q1 = B).
    # ---------------------------------------------------------------------
    lra   <- log(1 / 117);  label("Weibull rate-scaling parameter of the first (rapid) release process (RA = 1/TD, 1/h)")     # Table 1: TD = 117 h (RSE 3.8%); ra = 1/117 = 0.008547 1/h
    lgam1 <- log(1.4);      label("Weibull shape / sigmoidicity of the first (rapid) release process (GAM1, unitless)")       # Table 1: SS = 1.4 (RSE 1.9%)
    lra2  <- log(1 / 323);  label("Weibull rate-scaling parameter of the second (sustained) release process at zero dose (RA2 = 1/TD1_0, 1/h)")  # Table 1: TD1_0 = 323 h (RSE 5.3%); ra2 = 1/323 = 0.003096 1/h
    lgam2 <- log(3.25);     label("Weibull shape / sigmoidicity of the second (sustained) release process (GAM2, unitless)")  # Table 1: SS1 = 3.25 (RSE 4.4%)

    # Fraction of the dose entering the first release process. Held on the
    # logit scale so it stays inside (0, 1) for every eta draw; a log-scale
    # encoding of a bounded fraction can leak above 1 (operator sidecar
    # oare_PMC12775547 answer q1 = B). logit(0.236) = -1.17457.
    logitfrel <- log(0.236 / (1 - 0.236)); label("Logit of the fraction of the dose entering the first (rapid) release process (FF, unitless)")  # Table 1: FF = 0.236 (RSE 6.2%)

    # Disposition -- Table 1. The paper parameterises distribution by the
    # micro-constants k12 and k21 rather than by Q and Vp, so they are carried
    # as primary parameters here.
    lcl  <- log(16.6);     label("Apparent clearance CL/F at the 70 kg reference weight (L/h)")                       # Table 1: CL/F = 16.6 L/h (RSE 8.2%)
    lvc  <- log(15.4);     label("Apparent central volume of distribution V/F at the 70 kg reference weight (L)")     # Table 1: V/F = 15.4 L (RSE 14.5%)
    lk12 <- log(1.61);     label("First-order transfer rate constant central to peripheral1 (k12, 1/h)")              # Table 1: k12 = 1.61 1/h (RSE 12.9%)
    lk21 <- log(0.00313);  label("First-order transfer rate constant peripheral1 to central (k21, 1/h)")              # Table 1: k21 = 0.00313 1/h (RSE 11.9%)

    # Allometric body-weight exponents. ESTIMATED in the final model, not
    # fixed at 0.75 / 1 -- Table 1 reports an RSE for each, and a fixed
    # exponent has no RSE. See covariateData$WT and the vignette Errata.
    e_wt_cl <- 0.31;   label("Body-weight allometric exponent on CL/F (unitless)")   # Table 1: CL_P = 0.31 (RSE 86.8%)
    e_wt_vc <- 0.364;  label("Body-weight allometric exponent on V/F (unitless)")    # Table 1: V_P = 0.364 (RSE 110.7%)

    # Linear dose effect on the second-process release time (see
    # covariateData$DOSE_TV44749_MG for the reciprocal bookkeeping).
    e_dose_ra2 <- 0.0939; label("Linear effect of TV-44749 dose on the second-process release time TD1 = TD1_0 + e_dose_ra2*DOSE (h/mg)")  # Table 1: TD1_1 = 0.0939 (RSE 41.4%)

    # Published Emax dopamine D2 receptor occupancy layer, inherited from
    # Mamo 2008 (the paper's reference 13) rather than estimated here, so both
    # parameters are fixed. Perlstein 2026 Methods, PK/D2RO Model:
    #   D2RO = ROmax * Cp / (EC50 + Cp)
    emax  <- fixed(100);      label("Maximal attainable dopamine D2 receptor occupancy ROmax (%)")            # Methods PK/D2RO Model: ROmax fixed to 100% (Mamo 2008)
    lec50 <- fixed(log(11));  label("Plasma olanzapine concentration giving 50% D2 receptor occupancy (EC50, ng/mL)")  # Methods PK/D2RO Model: EC50 = 11 ng/mL, estimated by Mamo 2008

    # ---------------------------------------------------------------------
    # Inter-individual variability -- Table 1, Random effect block. The rows
    # are VARIANCES, not SDs: for a variance the RSE cannot fall much below
    # sqrt(2/N), which is 15.0% at N = 89 subjects, and every row here is at
    # or above that bound. The TD row also makes no sense as an SD (0.0038
    # would be a 0.38% CV carrying 74.5% shrinkage); as a variance it is a
    # 6.2% CV. Log-normal (exponential) IIV is assumed on every parameter --
    # the paper does not state the transform. See the vignette Errata.
    #
    # ra = 1/TD, so var(log ra) = var(-log TD) = var(log TD): the published
    # omega carries onto the reciprocal parameterisation unchanged.
    # ---------------------------------------------------------------------
    etalra        ~ 0.0038  # Table 1 Random effect TD = 0.0038 (RSE 95.0%, shrinkage 74.5%); 6.2% CV
    etalgam1      ~ 0.0121  # Table 1 Random effect SS = 0.0121 (RSE 21.8%, shrinkage 7.23%); 11.0% CV
    etalra2       ~ 0.0128  # Table 1 Random effect TD1_0 = 0.0128 (RSE 22.5%, shrinkage 17.0%); 11.3% CV
    etalgam2      ~ 0.0581  # Table 1 Random effect SS1 = 0.0581 (RSE 30.8%, shrinkage 22.9%); 24.4% CV
    etalogitfrel  ~ 0.0669  # Table 1 Random effect FF = 0.0669 (RSE 24.4%, shrinkage 23.8%); applied on the logit scale, see Errata
    etalcl        ~ 0.165   # Table 1 Random effect CL/F = 0.165 (RSE 24.7%, shrinkage 11.1%); 42.4% CV
    etalvc        ~ 0.131   # Table 1 Random effect V/F = 0.131 (RSE 58.9%, shrinkage 36.7%); 37.5% CV
    etalk12       ~ 0.214   # Table 1 Random effect k12 = 0.214 (RSE 48.6%, shrinkage 24.9%); 48.8% CV
    etalk21       ~ 0.445   # Table 1 Random effect k21 = 0.445 (RSE 27.4%, shrinkage 15.0%); 74.0% CV

    # Residual error -- Table 1. These rows are SDs, not variances: for a
    # variance the RSE cannot fall much below sqrt(2/N) = 2.38% at N = 3530
    # observations, and the rows report 1.2% and 0.5%.
    addSd  <- 0.632;  label("Additive residual error (ng/mL)")           # Table 1: Additive RSE model = 0.632 (RSE 1.2%)
    propSd <- 0.233;  label("Proportional residual error (fraction)")    # Table 1: Proportional RSE model = 0.233 (RSE 0.5%)
  })

  model({
    # 1. Release input function. ra and ra2 are the reciprocals of the
    #    paper's release times TD and TD1; the dose effect is additive on the
    #    TIME scale, so it is applied to 1/ra2 (Table 1 footnote).
    ra   <- exp(lra   + etalra)
    gam1 <- exp(lgam1 + etalgam1)
    gam2 <- exp(lgam2 + etalgam2)
    ra2  <- 1 / (1 / exp(lra2 + etalra2) + e_dose_ra2 * DOSE_TV44749_MG)

    # Fraction into the first release process, back-transformed from the
    # logit scale so that IIV is additive on the logit and frel is bounded.
    # The fixed effect and the eta are collected on their own line so the
    # term stays in a mu-referenced position (rxode2 otherwise warns that
    # etalogitfrel defaulted to non-mu-referenced); same shape as the
    # logitfm_ind idiom registered in parameter-names.md.
    logitfrel_ind <- logitfrel + etalogitfrel
    frel <- 1 / (1 + exp(-logitfrel_ind))

    # 2. Individual disposition parameters, allometrically scaled to the
    #    70 kg reference weight.
    cl  <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc  <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)

    # 3. Micro-constants
    kel <- cl / vc

    # 4. Weibull release hazards. Each depot empties with the hazard of its
    #    own Weibull survival function S(t) = exp(-(rate*t)^shape), whose
    #    hazard is h(t) = shape*rate*(rate*t)^(shape - 1). Two depots emptying
    #    this way is an exact ODE representation of the paper's Dose*(-dr/dt)
    #    input, and unlike an explicit -dr/dt term it superposes correctly
    #    across repeated injections. Time after dose is floored at a small
    #    positive number because (rate*t)^(shape - 1) is not finite at t = 0
    #    for shape < 1 and is evaluated on every record.
    tr1 <- max(tad(depot),  1e-6)
    tr2 <- max(tad(depot2), 1e-6)
    h1  <- gam1 * ra  * (ra  * tr1)^(gam1 - 1)
    h2  <- gam2 * ra2 * (ra2 * tr2)^(gam2 - 1)

    # 5. ODE system. The injected dose is recorded on both depots and split
    #    between them by the release fraction.
    d/dt(depot)       <- -h1 * depot
    d/dt(depot2)      <- -h2 * depot2
    d/dt(central)     <-  h1 * depot + h2 * depot2 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 6. Dose split across the two parallel release processes
    f(depot)  <- frel
    f(depot2) <- 1 - frel

    # 7. Observation. Dose in mg over volume in L gives mg/L; the factor of
    #    1000 converts to the ng/mL the paper reports.
    Cc <- 1000 * central / vc

    # Published Emax D2 receptor occupancy (percent), Mamo 2008 via
    # Perlstein 2026 Methods. Reported without residual error -- it is a
    # deterministic transform of the plasma concentration, not a measured
    # endpoint of this study.
    ec50 <- exp(lec50)
    D2RO <- emax * Cc / (ec50 + Cc)

    Cc ~ add(addSd) + prop(propSd)
  })
}
