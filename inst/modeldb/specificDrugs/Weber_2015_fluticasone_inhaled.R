Weber_2015_fluticasone_inhaled <- function() {
  description <- paste(
    "Semi-mechanistic. Population PK model for inhaled fluticasone propionate (FP)",
    "in healthy adult volunteers (Weber 2015), used for Monte-Carlo simulation of",
    "PK-based bioequivalence trials. Separate central (LC1 -> LC2) and peripheral",
    "(LP1 -> LP2) lung deposition compartments hold undissolved drug particles",
    "(LC1, LP1) and dissolved drug (LC2, LP2); mucociliary clearance kmuc removes",
    "undissolved particles from central lung regions only; dissolved drug is",
    "absorbed into a two-compartment systemic disposition with",
    "central-to-peripheral rate constants k12 and k21. Each administration",
    "splits across LC1 (bioavailability flung * fc) and LP1",
    "(bioavailability flung * (1 - fc)); the remaining (1 - flung) fraction is",
    "assumed to have negligible oral bioavailability. F_Lung and F_C are",
    "logit-transformed; all other parameters are log-transformed.",
    "Structural parameters and BSV were taken from the validated FP inhalation",
    "model of Weber and Hochhaus 2013 (reference 13 of Weber 2015); BOV on",
    "F_Lung, F_C, and kmuc described as a paper-specific extension for",
    "crossover-trial simulation is NOT encoded in this model file (see vignette",
    "Assumptions and deviations)."
  )
  reference <- paste(
    "Weber B, Hochhaus G. (2015).",
    "A Systematic Analysis of the Sensitivity of Plasma Pharmacokinetics to",
    "Detect Differences in the Pulmonary Performance of Inhaled Fluticasone",
    "Propionate Products Using a Model-Based Simulation Approach.",
    "AAPS J 17(4):999-1010. doi:10.1208/s12248-015-9768-y."
  )
  vignette <- "Weber_2015_fluticasone_inhaled"
  paper_specific_compartments <- c("LC1", "LC2", "LP1", "LP2")
  dosing <- c("LC1", "LP1")

  units <- list(time = "hour", dosing = "ug", concentration = "ng/mL")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = NA_integer_,
    age_range      = "Healthy adult volunteers (specifics not reported in Weber 2015 Methods)",
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = NA,
    disease_state  = "Healthy adult volunteers without respiratory disease.",
    dose_range     = paste(
      "Single inhaled dose; absolute dose level not specified.",
      "The PK model is scale-invariant in the emitted dose (no saturable",
      "elimination); reproducing the paper's reported AUC and Cmax T/R",
      "ratios is independent of the chosen dose magnitude.",
      "Bioequivalence trial simulations in Weber 2015 used a single-dose",
      "two-period crossover design with sample sizes of 10-200 subjects",
      "(Series 1: 10, 20, 30, 40, 50, 60, 70; Series 2 and 3: 50, 100, 200)."
    ),
    administration = paste(
      "Orally inhaled dry-powder fluticasone propionate (FP).",
      "The reference (R) product represents the FP Diskus inhaler;",
      "test (T) products differ from the reference in F_Lung, F_C,",
      "and/or k_diss."
    ),
    notes          = paste(
      "Structural model, typical parameter values, and BSV (Weber 2015",
      "Table I) are taken from the previously validated FP inhalation",
      "model of Weber and Hochhaus 2013 Mol Pharm 10(8):2873-85",
      "(reference 13 of Weber 2015). That upstream paper was the source",
      "of the BSV (30% CV across all parameters) and the proportional",
      "residual variability (20%); the current paper (Weber 2015) adds",
      "between-occasion variability (BOV; 30% CV) on F_Lung, F_C, and",
      "kmuc as the contribution that enables crossover-trial simulation",
      "for bioequivalence assessment. The current paper does not fit",
      "observed plasma concentrations."
    )
  )

  ini({
    # Lung-deposition fractions on the logit scale. Weber 2015 Methods
    # ("Pharmacokinetic Model" section, paragraph following Fig. 1)
    # states: "a logit transformation was applied to model parameters
    # that are naturally bound between 0 and 1 (F_Lung and F_C)".
    # Typical values from Table I (FP Diskus reference product).
    logitflung <- log(0.16 / (1 - 0.16)); label("Logit of fraction of emitted dose deposited in lung F_Lung (unitless)") # Weber 2015 Table I: F_Lung = 0.16
    logitfc    <- log(0.50 / (1 - 0.50)); label("Logit of fraction of lung dose deposited in central lung regions F_C (unitless)") # Weber 2015 Table I: F_C = 0.50

    # Lung-deposition rate constants on the log scale (Weber 2015 Table I,
    # FP Diskus reference product). All k values are first-order
    # (Methods, "Pharmacokinetic Model" section).
    lkdiss <- log(0.189); label("Log pulmonary dissolution rate kdiss (1/h)")                     # Weber 2015 Table I: k_diss  = 0.189 1/h
    lkmuc  <- log(0.938); label("Log mucociliary clearance rate kmuc, central lung only (1/h)")    # Weber 2015 Table I: k_muc   = 0.938 1/h
    lkpulc <- log(10);    label("Log pulmonary absorption rate from central lung kpulc (1/h)")     # Weber 2015 Table I: k_pul,C = 10 1/h
    lkpulp <- log(20);    label("Log pulmonary absorption rate from peripheral lung kpulp (1/h)")  # Weber 2015 Table I: k_pul,P = 20 1/h

    # Systemic two-compartment body PK (Weber 2015 Table I). The paper
    # parameterises the body 2-compartment with rate constants k12 and
    # k21 directly (rather than the canonical q / vp pair), so the
    # primary ini parameters use the rate-constant form lk12 / lk21
    # to preserve the paper's independent-eta structure on each rate.
    lcl  <- log(73);   label("Log systemic clearance CL (L/h)")                              # Weber 2015 Table I: CL  = 73 L/h
    lvc  <- log(31);   label("Log central body volume Vc (L)")                               # Weber 2015 Table I: V_C = 31 L
    lk12 <- log(1.78); label("Log central-to-peripheral body distribution rate k12 (1/h)")    # Weber 2015 Table I: k12 = 1.78 1/h
    lk21 <- log(0.09); label("Log peripheral-to-central body distribution rate k21 (1/h)")    # Weber 2015 Table I: k21 = 0.09 1/h

    # Between-subject variability (BSV; Weber 2015 Table I). Reported as
    # 30% CV for every parameter. For log-normal parameters the log-scale
    # variance is omega^2 = log(1 + CV^2) = log(1 + 0.30^2) = log(1.09)
    # ~= 0.0862. For the two logit-transformed parameters (logitflung,
    # logitfc) the same omega^2 is used on the logit scale; the resulting
    # back-transformed CV is approximate and depends on the typical value
    # (see vignette Assumptions and deviations for the delta-method
    # check). BSV terms are diagonal (no inter-parameter correlation
    # reported in Weber 2015 Table I).
    etalogitflung ~ 0.0862                # Weber 2015 Table I: BSV F_Lung  30% CV -> log(1 + 0.30^2) = 0.0862
    etalogitfc    ~ 0.0862                # Weber 2015 Table I: BSV F_C     30% CV -> log(1 + 0.30^2) = 0.0862
    etalkdiss     ~ 0.0862                # Weber 2015 Table I: BSV k_diss  30% CV -> log(1 + 0.30^2) = 0.0862
    etalkmuc      ~ 0.0862                # Weber 2015 Table I: BSV k_muc   30% CV -> log(1 + 0.30^2) = 0.0862
    etalkpulc     ~ 0.0862                # Weber 2015 Table I: BSV k_pul,C 30% CV -> log(1 + 0.30^2) = 0.0862
    etalkpulp     ~ 0.0862                # Weber 2015 Table I: BSV k_pul,P 30% CV -> log(1 + 0.30^2) = 0.0862
    etalcl        ~ 0.0862                # Weber 2015 Table I: BSV CL      30% CV -> log(1 + 0.30^2) = 0.0862
    etalvc        ~ 0.0862                # Weber 2015 Table I: BSV V_C     30% CV -> log(1 + 0.30^2) = 0.0862
    etalk12       ~ 0.0862                # Weber 2015 Table I: BSV k12     30% CV -> log(1 + 0.30^2) = 0.0862
    etalk21       ~ 0.0862                # Weber 2015 Table I: BSV k21     30% CV -> log(1 + 0.30^2) = 0.0862

    # Residual error. Weber 2015 Table I footer states "Residual
    # variability was 20%". Encoded as proportional residual SD; the
    # paper does not report a combined additive + proportional form.
    propSd <- 0.20; label("Proportional residual SD (fraction)") # Weber 2015 Table I footer: RV = 20%
  })
  model({
    # Individual parameters. logitflung / logitfc back-transform to [0, 1]
    # via the logistic; all other parameters back-transform via exp().
    # The mu-referenced logit sum sits on its own line so nlmixr2's parser
    # recognises the eta as mu-referenced (Aguiar 2021 ustekinumab pattern).
    logit_flung_i <- logitflung + etalogitflung
    logit_fc_i    <- logitfc    + etalogitfc
    flung <- 1 / (1 + exp(-logit_flung_i))
    fc    <- 1 / (1 + exp(-logit_fc_i))
    kdiss <- exp(lkdiss + etalkdiss)
    kmuc  <- exp(lkmuc  + etalkmuc)
    kpulc <- exp(lkpulc + etalkpulc)
    kpulp <- exp(lkpulp + etalkpulp)
    cl    <- exp(lcl    + etalcl)
    vc    <- exp(lvc    + etalvc)
    k12   <- exp(lk12   + etalk12)
    k21   <- exp(lk21   + etalk21)

    kel <- cl / vc

    # Lung kinetics (Weber 2015 Fig. 1 and caption). LC1 and LP1 hold
    # undissolved FP particles in central and peripheral lung regions,
    # respectively; LC2 and LP2 hold dissolved FP ready for pulmonary
    # absorption. Mucociliary clearance kmuc removes undissolved
    # particles from central lung regions only (no kmuc on LP1).
    d/dt(LC1) <- -kdiss * LC1 - kmuc * LC1
    d/dt(LC2) <-  kdiss * LC1 - kpulc * LC2
    d/dt(LP1) <- -kdiss * LP1
    d/dt(LP2) <-  kdiss * LP1 - kpulp * LP2

    # Systemic two-compartment disposition fed by pulmonary absorption
    # from both lung regions (Weber 2015 Fig. 1 caption).
    d/dt(central)     <-  kpulc * LC2 + kpulp * LP2 - kel * central -
                          k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Dose distribution. Each administration enters the event table as
    # two simultaneous dose rows (cmt = LC1 and cmt = LP1, each with
    # amt = total emitted dose); the bioavailabilities below split the
    # emitted dose into the central-lung region (flung * fc) and the
    # peripheral-lung region (flung * (1 - fc)). The remaining
    # (1 - flung) fraction (oropharyngeal deposition and swallowed
    # drug) is assumed to have negligible systemic bioavailability per
    # Weber 2015 Methods and Discussion.
    f(LC1) <- flung * fc
    f(LP1) <- flung * (1 - fc)

    # Plasma concentration. Dose enters in ug, vc is in L, so
    # central / vc yields ug/L = ng/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
