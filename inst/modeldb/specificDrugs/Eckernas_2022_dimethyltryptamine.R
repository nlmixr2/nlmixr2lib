Eckernas_2022_dimethyltryptamine <- function() {
  description <- "Population PK/PD model of intravenous N,N-dimethyltryptamine (DMT) in healthy adults (Eckernas 2022). Plasma DMT is described by a two-compartment model with first-order elimination from the central compartment; all DMT elimination forms the major metabolite indole-3-acetic acid (IAA, metabolic fraction fixed to 1), whose plasma increment over endogenous baseline follows a one-compartment model with first-order elimination. The subjective intensity of the psychedelic experience (0-10 rating) is driven by the effect-site concentration through a sigmoidal Emax function with zero baseline and Emax fixed to 10, the top of the rating scale. Between-subject variability is on DMT clearance, IAA volume, the effect-site EC50 and the Hill coefficient. Residual error is proportional for DMT and IAA and additive on the logit (0-10) scale for the intensity rating."
  reference <- paste(
    "Eckernas E, Timmermann C, Carhart-Harris R, Roshammar D, Ashton M.",
    "Population pharmacokinetic/pharmacodynamic modeling of the psychedelic",
    "experience induced by N,N-dimethyltryptamine - Implications for dose",
    "considerations. Clin Transl Sci. 2022;15(12):2928-2937.",
    "doi:10.1111/cts.13410.",
    sep = " "
  )
  vignette <- "Eckernas_2022_dimethyltryptamine"
  units <- list(time = "min", dosing = "nmol", concentration = "nmol/L")

  compartmentData <- list(
    central = list(analyte = "DMT", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "DMT", units = "nmol", specimen = "not applicable", verified = TRUE),
    central_iaa = list(
      analyte = "IAA formed from DMT (increment over endogenous baseline)",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    effect = list(analyte = "DMT", units = "nmol/L", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 13L, # Methods (Clinical trial): 13 healthy subjects
    n_studies = 1L, # single placebo-controlled, fixed-order study at the NIHR Imperial Clinical Research Facility
    age_range = "22-48 years", # Methods (Clinical trial)
    age_median = "33 years", # Methods (Clinical trial)
    sex_female_pct = 46.2, # Methods: 7 men of 13 subjects; (13 - 7) / 13 = 46.2%
    disease_state = "healthy volunteers",
    dose_range = "Placebo at visit 1, then one week later a single intravenous bolus of 7 mg (n = 3), 14 mg (n = 4), 18 mg (n = 1) or 20 mg (n = 5) DMT fumarate; doses were converted to nmol of DMT base for modelling.", # Methods (Clinical trial; Modeling approach)
    regions = "United Kingdom (London).",
    n_observations = "93 DMT and 87 IAA plasma concentrations (nine samples per subject up to 60 min) and 273 subjective intensity ratings (one per minute, 0-20 min).", # Results, first paragraph
    notes = "Intensity was rated on a 0-10 scale (0 = no effect, 10 = the most intense experience imaginable). IAA was measured as the change from the high endogenous baseline. Placebo-visit data contained no measurable DMT and no subjective effect and were not modelled. No covariate analysis was performed because of the small sample size."
  )

  ini({
    # ---- DMT plasma PK (Table 1) --------------------------------------------
    lcl <- log(26.0);  label("DMT clearance CL (L/min)")                                # Table 1: CL = 26.0 (20.6-33.6) L/min, %RSE 15.1
    lvc <- log(221);   label("DMT central volume of distribution Vc (L)")               # Table 1: Vc = 221 (181-273) L, %RSE 12.1
    lq  <- log(2.99);  label("DMT intercompartmental clearance Q (L/min)")              # Table 1: Q = 2.99 (1.87-5.44) L/min, %RSE 36.6
    lvp <- log(59.0);  label("DMT peripheral volume of distribution Vp (L)")            # Table 1: Vp = 59.0 (48.0-82.7) L, %RSE 18.5

    # ---- IAA metabolite PK (Table 1) ----------------------------------------
    fm         <- fixed(1);     label("Fraction of DMT elimination forming IAA (unitless)")   # Methods (PK model development): 'fixed metabolic fraction of 1'
    lcl_iaa    <- log(0.093);   label("Apparent IAA clearance CL(m) (L/min)")                 # Table 1: CL(m) = 0.093 (0.076-0.11) L/min, %RSE 10.6
    lvc_iaa    <- log(9.55);    label("Apparent IAA volume of distribution V(m) (L)")         # Table 1: V(m) = 9.55 (8.07-11.1) L, %RSE 9.8

    # ---- Intensity-rating PD (Table 1) --------------------------------------
    emax  <- fixed(10);   label("Maximum intensity rating Emax (rating units, 0-10 scale)")   # Table 1: Emax = 10 FIX; Results: highest allowed rating
    lec50 <- log(94.7);   label("Effect-site concentration for half-maximal rating EC50,e (nmol/L)")  # Table 1: EC50,e = 94.7 (75.4-123) nM, %RSE 14.9
    lke0  <- log(1.38);   label("Effect-compartment equilibration rate constant ke0 (1/min)")  # Table 1: ke0 = 1.38 (1.16-1.94) 1/min, %RSE 17.5
    lhill <- log(2.87);   label("Hill coefficient gamma (unitless)")                          # Table 1: gamma = 2.87 (2.63-3.02), %RSE 4.2

    # ---- Between-subject variability ----------------------------------------
    # Table 1 reports BSV as %CV = sqrt(omega^2) * 100: the same authors'
    # 2023 EEG analysis carries this BSV CL forward as `$OMEGA 0.224 FIX`
    # = 0.473^2 (doi:10.1002/psp4.12933, Appendix S1). The variances below are
    # therefore (CV/100)^2.
    etalcl     ~ 0.2237   # Table 1: BSV CL = 47.3 %CV; 0.473^2
    etalvc_iaa ~ 0.0841   # Table 1: BSV V(m) = 29.0 %CV; 0.290^2
    etalec50   ~ 0.1490   # Table 1: BSV EC50,e = 38.6 %CV; 0.386^2
    etalhill   ~ 0.5975   # Table 1: BSV gamma = 77.3 %CV; 0.773^2

    # ---- Residual error -----------------------------------------------------
    propSd                 <- 0.500; label("Proportional residual error on plasma DMT (fraction)")   # Table 1: residual error DMT = 50.0 (43.9-58.9) %CV
    propSd_iaa             <- 0.190; label("Proportional residual error on plasma IAA (fraction)")   # Table 1: residual error IAA = 19.0 (16.8-21.8) %CV
    addSd_psychedelicintensity <- 0.82;  label("Additive residual SD of the intensity rating on the logit scale")  # Table 1: residual error intensity ratings = 0.82 (0.55-0.90) SD; Eq 1
  })

  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    cl_iaa <- exp(lcl_iaa)
    vc_iaa <- exp(lvc_iaa + etalvc_iaa)

    ec50 <- exp(lec50 + etalec50)
    ke0  <- exp(lke0)
    hill <- exp(lhill + etalhill)

    kel     <- cl / vc
    k12     <- q / vc
    k21     <- q / vp
    kel_iaa <- cl_iaa / vc_iaa

    # Doses are nmol of DMT base, so central / vc is in nmol/L (= nM). IAA is
    # formed at the DMT elimination rate (fm = 1), mole for mole; its state is
    # the increment over the endogenous IAA baseline, which is how IAA was
    # measured.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(central_iaa) <-  fm * kel * central - kel_iaa * central_iaa

    Cc     <- central / vc
    Cc_iaa <- central_iaa / vc_iaa

    # Effect compartment (biophase): dCe/dt = ke0 * (Cp - Ce).
    d/dt(effect) <- ke0 * (Cc - effect)

    # Sigmoidal Emax response with zero baseline (all subjects rated 0 before
    # dosing). Eq 1 places the residual error on the logit scale bounded by
    # 0 and 10, which is rxode2's logitNorm(sd, 0, 10).
    psychedelicintensity <- emax * effect^hill / (ec50^hill + effect^hill)

    Cc              ~ prop(propSd)
    Cc_iaa          ~ prop(propSd_iaa)
    psychedelicintensity ~ logitNorm(addSd_psychedelicintensity, 0, 10)
  })
}
