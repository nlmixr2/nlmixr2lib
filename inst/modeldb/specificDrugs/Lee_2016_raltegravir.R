Lee_2016_raltegravir <- function() {
  description <- paste(
    "Population PK model for oral raltegravir (a UGT1A1 phenotyping probe)",
    "and its glucuronide metabolite in 24 East Asian patients with advanced",
    "solid tumours receiving FOLFIRI chemotherapy (Lee 2016). Raltegravir",
    "absorption is described with a depot, a single transit compartment",
    "(the paper estimates a non-integer NN = 1.07 in the Savic 2007",
    "transit-chain framework; the packaged model approximates this with",
    "one explicit transit compartment), and a one-compartment central",
    "compartment with first-order elimination (CL/F, V/F). Raltegravir",
    "glucuronide is described by a one-compartment metabolite compartment",
    "(central_gluc) with V_GLU fixed at 1 L (a structural identifiability",
    "anchor) and a first-order metabolite clearance CL_GLU. The formation",
    "rate constant kmet maps to the source paper's FMET, which the authors",
    "define as the formation rate of glucuronide divided by V_GLU; with",
    "V_GLU fixed at 1 L, kmet has units 1/h and drives dA_gluc / dt =",
    "kmet * V_GLU * C_RAL_central - CL_GLU * C_gluc. Bioavailability F",
    "is fixed at 1 (single oral dose; absolute F not identifiable). IIV",
    "is reported on CL/F, MTT, F, V/F, kmet (FMET), and CL_GLU with a",
    "single off-diagonal covariance between CL/F and V/F (correlation",
    "0.567). The residual error was reported as additive on log-transformed",
    "observations for both raltegravir and glucuronide, which maps to a",
    "proportional residual on the linear-concentration scale. No baseline",
    "covariates (age, sex, weight, body surface area, serum albumin /",
    "creatinine / bilirubin / liver enzymes, ethnicity, or UGT1A1 * 6 /",
    "* 28 / * 60 and CYP3A5 * 3 genotypes) were retained in the final",
    "model.")
  reference <- paste(
    "Lee LS, Seng KY, Wang LZ, Yong WP, Hee KH, Soh TI, Wong A, Cheong PF,",
    "Soong R, Sapari NS, Soo R, Fan L, Lee SC, Goh BC.",
    "Phenotyping of UGT1A1 Activity Using Raltegravir Predicts Pharmacokinetics",
    "and Toxicity of Irinotecan in FOLFIRI.",
    "PLoS One. 2016 Jan 25;11(1):e0147681.",
    "doi:10.1371/journal.pone.0147681.",
    sep = " "
  )
  vignette <- "Lee_2016_raltegravir"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 24L,
    n_studies      = 1L,
    age_range      = "39-79 years",
    age_median     = "59 years",
    weight_range   = "42.4-81.1 kg",
    weight_median  = "55 kg",
    sex_female_pct = 21,
    race_ethnicity = "Asian: 75% Chinese, 21% Malay, 4% Indian (Singapore-recruited cohort; Table 1)",
    disease_state  = paste(
      "Asian patients with advanced-stage gastro-intestinal cancers requiring",
      "FOLFIRI chemotherapy (folinic acid, 5-fluorouracil, irinotecan).",
      "Inclusion required Karnofsky performance status > 70% and adequate",
      "haematologic / renal / hepatic function. Recruited in Singapore",
      "between 2009 and 2011.",
      sep = " "
    ),
    dose_range     = paste(
      "Single oral dose of 400 mg raltegravir as a UGT1A1 phenotyping probe,",
      "administered fasted one day before FOLFIRI cycle 1. Serial blood",
      "samples were collected at baseline and at 0.5, 1, 2, 4, 6, 8, and",
      "24 h post-dose.",
      sep = " "
    ),
    regions        = "Singapore (single centre, National University Health System)",
    co_medication  = paste(
      "Intravenous midazolam 1 mg co-administered as a CYP3A4 probe on the",
      "same day; FOLFIRI started the day after raltegravir. The midazolam",
      "PK model was not refit in this paper but reused from the authors'",
      "prior publication (Hee 2015, reference 16 in the source).",
      sep = " "
    ),
    trial_id       = "ClinicalTrials.gov NCT00808184",
    notes          = paste(
      "Baseline demographics from Table 1 of Lee 2016. The eta-shrinkage",
      "values reported for the final raltegravir parameters are CL_RAL/F",
      "24.1%, V_RAL/F 23.7%, MTT 15.5%, F 1.4% (parent model) and CL_GLU",
      "30.0%, FMET 28.3% (glucuronide model). Genotype frequencies in the",
      "cohort: CYP3A5 * 3 wt/het/var = 6/7/11; UGT1A1 * 6 wt/het/var = 17/7/0;",
      "UGT1A1 * 28 wt/het/var = 17/7/0; UGT1A1 * 60 wt/het/var = 12/9/3.",
      "No CYP3A4 variants were detected. The paper tested age, sex, body",
      "weight, body surface area, serum albumin, serum creatinine, total",
      "bilirubin, ALP, ALT, AST, race / ethnicity, and the UGT1A1 / CYP3A5",
      "genotypes as candidate covariates on the raltegravir-glucuronide",
      "model; none was retained in the final model reported in Table 2.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Raltegravir parent absorption parameters (Lee 2016 Table 2)
    # ------------------------------------------------------------------
    # The source uses a Savic 2007 transit-compartment chain with the
    # estimated number of (hypothetical) transit compartments NN = 1.07
    # and mean transit time MTT = 1.04 h, plus a separate first-order
    # absorption rate ka = 4.23 1/h (Fig 2 of Lee 2016 and the analytic
    # transit-chain formulation in reference [19] of the source). nlmixr2lib
    # does not implement the Savic 2007 non-integer transit chain; the
    # packaged model approximates the published NN = 1.07 with one explicit
    # transit compartment between depot and central. The transit rate
    # constant is computed inside model() as ktr = (NN + 1) / MTT with
    # NN approximated as 1, so the IIV applied to log(MTT) directly
    # propagates to the transit rate. See vignette Assumptions and
    # deviations for the discretisation rationale.

    lka     <- log(4.23);     label("First-order absorption rate constant ka (1/h)")             # Table 2: ka = 4.23 1/h (RSE 19.3%)
    lmtt    <- log(1.04);     label("Mean transit time MTT through the transit chain (h)")       # Table 2: MTT = 1.04 h (RSE 25.8%)
    lfdepot <- fixed(log(1)); label("Oral bioavailability F (unitless, fixed at 1)")             # Table 2: F = 1 FIXED (no IV reference data; absolute F not identifiable)

    # ------------------------------------------------------------------
    # Raltegravir parent disposition parameters (Lee 2016 Table 2)
    # ------------------------------------------------------------------
    lcl <- log(41.7);   label("Apparent oral clearance CL_RAL/F (L/h)")                          # Table 2: CL_RAL/F = 41.7 L/h (RSE 24.9%)
    lvc <- log(157);    label("Apparent oral central volume of distribution V_RAL/F (L)")        # Table 2: V_RAL/F = 157 L (RSE 32.4%)

    # ------------------------------------------------------------------
    # Glucuronide metabolite parameters (Lee 2016 Table 2)
    # ------------------------------------------------------------------
    # The metabolite compartment is parametrised with V_GLU fixed at 1 L
    # (a structural identifiability anchor; absolute glucuronide volume is
    # not identifiable from plasma data without a directly administered
    # metabolite reference). The source paper names the formation parameter
    # FMET and defines it as "the ratio of the formation rate of glucuronide
    # to V_GLU" -- i.e. a formation rate constant in 1/h when applied to
    # the parent central concentration. This maps to the nlmixr2lib paper-
    # named parameter `kmet` (formation rate constant from parent central
    # to metabolite); see R/conventions.R::paperNamedParams.
    lkmet_gluc <- log(0.0324); label("Glucuronide formation rate constant kmet (= FMET; 1/h)")    # Table 2: FMET = 0.0324 (RSE 10.4%)
    lcl_gluc   <- log(0.715);  label("Glucuronide first-order clearance CL_GLU (L/h)")            # Table 2: CL_GLU = 0.715 L/h (RSE 10.2%)
    lvc_gluc   <- fixed(log(1)); label("Glucuronide volume of distribution V_GLU (L, fixed)")     # Table 2: V_GLU = 1 L FIXED

    # ------------------------------------------------------------------
    # Inter-individual variability
    # ------------------------------------------------------------------
    # The paper reports IIV as %CV on the lognormal scale; the internal
    # variance on the log-scale is recovered as omega^2 = log(CV^2 + 1).
    # CL_RAL/F and V_RAL/F have a reported correlation coefficient of
    # 0.567 (Table 2: "Correlation CL_RAL/F, V_RAL/F = 0.57"); the
    # corresponding covariance is cov = 0.567 * SD_lcl * SD_lvc where
    # SD_lcl = sqrt(omega^2_lcl) and SD_lvc = sqrt(omega^2_lvc).
    #
    #   CV CL_RAL/F = 30.5%  -> omega^2 = log(0.305^2 + 1) = 0.0890
    #   CV V_RAL/F  = 81.4%  -> omega^2 = log(0.814^2 + 1) = 0.5085
    #   cov(lcl, lvc)        = 0.567 * sqrt(0.0890 * 0.5085) = 0.1206
    #   CV MTT      = 107.2% -> omega^2 = log(1.072^2 + 1) = 0.7651
    #   CV F        = 123.7% -> omega^2 = log(1.237^2 + 1) = 0.9283
    #   CV FMET     = 37.7%  -> omega^2 = log(0.377^2 + 1) = 0.1329
    #   CV CL_GLU   = 13.6%  -> omega^2 = log(0.136^2 + 1) = 0.0183
    etalcl + etalvc ~ c(0.0890,
                        0.1206, 0.5085)                                                          # Table 2: omega(CL_RAL/F)=30.5%CV, omega(V_RAL/F)=81.4%CV, correlation = 0.567
    etalmtt        ~ 0.7651                                                                       # Table 2: omega(MTT)    = 107.2 %CV
    etalfdepot     ~ 0.9283                                                                       # Table 2: omega(F)      = 123.7 %CV
    etalkmet_gluc  ~ 0.1329                                                                       # Table 2: omega(FMET)   = 37.7 %CV
    etalcl_gluc    ~ 0.0183                                                                       # Table 2: omega(CL_GLU) = 13.6 %CV

    # ------------------------------------------------------------------
    # Residual error -- additive on log-transformed observations
    # ------------------------------------------------------------------
    # The paper Methods (p. 5) state the residual was best described by
    # an additive model on log-transformed data for both raltegravir and
    # glucuronide. A NONMEM additive-on-log-scale residual maps to an
    # nlmixr2 proportional residual on the linear-concentration scale
    # (see references/parameter-names.md "Residual error"). Table 2 gives
    # SD values directly: sigma_RAL = 0.15, sigma_GLU = 0.18.
    propSd      <- 0.15; label("Proportional residual SD for raltegravir parent (fraction; SD on log scale)") # Table 2: sigma_RAL = 0.15 (RSE 3.2%); paper Methods: additive on log-transformed
    propSd_gluc <- 0.18; label("Proportional residual SD for raltegravir glucuronide (fraction; SD on log scale)") # Table 2: sigma_GLU = 0.18 (RSE 2.6%); paper Methods: additive on log-transformed
  })

  model({
    # ---- Individual raltegravir parameters ----
    ka       <- exp(lka)
    mtt      <- exp(lmtt + etalmtt)
    fdepot   <- exp(lfdepot + etalfdepot)
    cl       <- exp(lcl + etalcl)
    vc       <- exp(lvc + etalvc)

    # ---- Individual glucuronide parameters ----
    kmet_gluc <- exp(lkmet_gluc + etalkmet_gluc)
    cl_gluc   <- exp(lcl_gluc + etalcl_gluc)
    vc_gluc   <- exp(lvc_gluc)

    # ---- Transit-chain rate constant ----
    # With NN approximated as one explicit transit compartment between
    # depot and central, the chain has NN + 1 = 2 transitions (depot ->
    # transit1 and transit1 -> central). The mean transit time of the
    # entire chain is then MTT = (NN + 1) / ktr, so ktr = 2 / MTT. ka is
    # retained as a separate first-order absorption rate constant from
    # depot (Lee 2016 Fig 2: ka labels the depot -> plasma absorption
    # rate and ktr labels the inter-transit rate; here ka is applied to
    # depot -> transit1 and ktr to transit1 -> central).
    ktr <- 2 / mtt

    # ---- Micro-constants ----
    kel      <- cl / vc
    kel_gluc <- cl_gluc / vc_gluc

    # ---- ODE system: depot -> transit1 -> central (raltegravir),
    #                  central -> central_gluc (glucuronide) ----
    d/dt(depot)        <- -ka * depot
    d/dt(transit1)     <-  ka * depot - ktr * transit1
    d/dt(central)      <-  ktr * transit1 - kel * central
    # Glucuronide formation is driven by raltegravir parent concentration
    # scaled to the metabolite distribution volume: formation rate per
    # V_GLU is kmet, so the total formation rate (mass / time) into the
    # metabolite compartment is kmet * V_GLU * C_RAL_parent. With V_GLU
    # fixed at 1 L, the simulated glucuronide "concentration" lives in a
    # scaled frame and is not directly comparable to physical plasma
    # glucuronide concentrations; see vignette Assumptions and deviations.
    d/dt(central_gluc) <-  kmet_gluc * vc_gluc * (central / vc) - kel_gluc * central_gluc

    # ---- Bioavailability ----
    f(depot) <- fdepot

    # ---- Observations ----
    Cc      <- central / vc
    Cc_gluc <- central_gluc / vc_gluc

    Cc      ~ prop(propSd)
    Cc_gluc ~ prop(propSd_gluc)
  })
}
