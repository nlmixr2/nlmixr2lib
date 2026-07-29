Muller_2010_clindamycin <- function() {
  description <- paste(
    "Three-compartment intravenous popPK model for clindamycin in pregnant",
    "women during the peripartum period (Muller 2010). Fit to 175 maternal",
    "venous serum concentrations from 7 women receiving either 600 mg over",
    "20 min every 6 h (endocarditis prophylaxis) or 900 mg over 30 min",
    "every 8 h (group B streptococcal disease prophylaxis). No covariates",
    "were retained in the final model; demographic and laboratory screens",
    "(maternal age, gestational age, BMI, weight, edema, temperature,",
    "creatinine, ALP, AST, ALT, mode of delivery) are documented in",
    "covariatesDataExcluded. Proportional residual error with a per-subject",
    "log-normal scaling eta on the residual error magnitude (NONMEM",
    "omega-sigma interaction with an extra ETA on epsilon)."
  )
  reference <- paste(
    "Muller AE, Mouton JW, Oostvogel PM, Dorr PJ, Voskuyl RA, DeJongh J,",
    "Steegers EA, Danhof M.",
    "Pharmacokinetics of clindamycin in pregnant women in the peripartum period.",
    "Antimicrob Agents Chemother. 2010;54(5):2175-2181.",
    "doi:10.1128/AAC.01017-09.",
    sep = " "
  )
  vignette <- "Muller_2010_clindamycin"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  # Screened in the exploratory covariate analysis (paper Methods,
  # "Pharmacokinetic analysis"; Results: "None of the covariates could
  # improve the model fit"). Recorded here for provenance so that
  # downstream users can see what was tested. None of these are referenced
  # in model(); checkModelConventions() treats covariatesDataExcluded as
  # documentation only.
  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Maternal age at the time of clindamycin administration.",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Cohort summary (Table 1): mean 36.1 years, SD 4.24, range 31.3-41.8 (n=7). Screened in the exploratory covariate analysis but not retained in the final model.",
      source_name        = "Maternal age"
    ),
    GA = list(
      description        = "Gestational age at the time of clindamycin administration (close to gestational age at delivery for the 6 women in labor; identical to GA at birth for the 6 singleton-pregnancy labor cohort, ambiguous for the 1 preterm-PROM twin pregnancy where antibiotic preceded delivery).",
      units              = "week",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Cohort summary (Table 1): mean 38.3 weeks, SD 3.01, range 34-42.3 (n=7). The canonical GA register entry is 'gestational age at birth'; in this study the reported value is gestational age at the antibiotic infusion, which equals GA at birth for the 6 women treated in labor.",
      source_name        = "Gestational age"
    ),
    BMI = list(
      description        = "Maternal body mass index at the time of clindamycin administration.",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Cohort summary (Table 1): mean 32.1, SD 5.36, range 22.1-39.1 (n=7).",
      source_name        = "Body mass index"
    ),
    WT = list(
      description        = "Maternal body weight at the time of clindamycin administration.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Cohort summary (Table 1): mean 86.1 kg, SD 14.2, range 59.5-104.8 (n=7).",
      source_name        = "Wt"
    ),
    BODYTEMP = list(
      description        = "Maternal oral body temperature at the time of clindamycin administration.",
      units              = "degC",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Cohort summary (Table 1): mean 36.9 degC, SD 0.23, range 36.7-37.4 (n=7).",
      source_name        = "Temp"
    ),
    CREAT = list(
      description        = "Maternal serum creatinine concentration at the time of clindamycin administration.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Cohort summary (Table 1): mean 55.9 umol/L, SD 20.3, range 38-100 (n=7).",
      source_name        = "Creatinine concn"
    ),
    ALP = list(
      description        = "Maternal serum alkaline phosphatase activity at the time of clindamycin administration.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Cohort summary (Table 1): mean 495 U/L, SD 601, range 168-1794 (n=7). The wide range reflects the physiological rise of placental ALP during late pregnancy.",
      source_name        = "ALP"
    ),
    AST = list(
      description        = "Maternal serum aspartate aminotransferase activity at the time of clindamycin administration.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Cohort summary (Table 1): mean 22.6 U/L, SD 6.78, range 15-34 (n=7).",
      source_name        = "AST"
    ),
    ALT = list(
      description        = "Maternal serum alanine aminotransferase activity at the time of clindamycin administration.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Cohort summary (Table 1): mean 10.9 U/L, SD 5.27, range 5-18 (n=7).",
      source_name        = "ALT"
    )
  )
  # Edema (none / ankle / up-to-knee) and mode of delivery
  # (vaginal / secondary cesarean) were also tabulated in Table 1 and
  # screened. There is no canonical register entry for either, so they
  # are documented here in prose only rather than as covariatesDataExcluded
  # list entries to avoid registering uncanonicalised covariate names.

  population <- list(
    species         = "human",
    n_subjects      = 7L,
    n_studies       = 1L,
    age_range       = "31.3-41.8 years (mean 36.1, SD 4.24)",
    weight_range    = "59.5-104.8 kg (mean 86.1, SD 14.2)",
    sex_female_pct  = 100,
    disease_state   = paste(
      "Pregnant women in the peripartum period requiring intravenous",
      "clindamycin for either group B streptococcal disease prophylaxis",
      "(n=4), endocarditis prophylaxis (n=2), or both (n=1). Six women",
      "had singleton pregnancies in labor; one had a preterm premature",
      "rupture of membranes with a twin pregnancy. None of the patients",
      "was clinically ill. Underlying cardiac indications: minor mitral",
      "stenosis (n=1), prosthetic aortic + prosthetic mitral valves with",
      "tricuspid valve surgical history (n=1), aortic dysfunction with",
      "autoimmune thrombocytopenic purpura + prior splenectomy (n=1)."
    ),
    dose_range      = paste(
      "600 mg over 20 min (12 mg/mL in 0.9% NaCl) every 6 h, or 900 mg",
      "over 30 min (9 mg/mL in 0.9% NaCl) every 8 h, intravenous infusion.",
      "Three women received additional postpartum doses. Average dose 814 mg",
      "and infusion rate 1743 mg/h were used by the authors for the VPC",
      "simulations (Figure 3 caption)."
    ),
    ga_range        = "34-42.3 weeks (mean 38.3, SD 3.01) gestational age at antibiotic administration",
    bmi_range       = "22.1-39.1 kg/m^2 (mean 32.1, SD 5.36)",
    mode_of_delivery = "4 vaginal / 3 secondary cesarean",
    regions         = "Medical Centre Haaglanden, The Hague, Netherlands (single-centre, Feb 2005 - Mar 2007)",
    n_observations  = 175L,
    notes           = paste(
      "Maternal venous serum concentrations only; arterial and venous",
      "umbilical cord blood concentrations were collected (Table 2) but",
      "could not be incorporated into the PK model due to the limited",
      "number of cord samples (the paper reports a venous-cord/maternal",
      "ratio range 0.22-0.89 with one extreme of 1.59). 177 samples",
      "collected; 2 excluded for being below the 0.1 mg/L LLOQ. NONMEM",
      "version VI, release 2 with ADVAN11 subroutine, FOCE-I estimation.",
      "Bootstrap validation: 598 of 1000 runs converged successfully."
    )
  )

  # The paper's residual-error model is a proportional error with an
  # additional per-subject log-normal scaling factor on the residual
  # magnitude (Table 3 "IIV in error" row). Encoded here via a fixed
  # anchor lrv and an estimated etalrv that scale propSdBase in model();
  # the same idiom is used in Dogterom_2018_asenapine.R for an IIV-on-rV
  # term.
  paper_specific_residual_sds <- c("propSdBase")

  ini({
    # ============================================================
    # Structural PK parameters (Muller 2010 Table 3, "Structural
    # model" block). NONMEM ADVAN11 (3-compartment IV) reported on
    # the linear scale; log-transformed here per nlmixr2lib convention.
    # ============================================================
    lcl  <- log(10.0);  label("CL: clearance from central compartment (L/h)")                       # Table 3: CL = 10.0 L/h (RSE 37.5%, 95% CI 2.65-17.35)
    lvc  <- log(12.4);  label("V1: central volume of distribution (L)")                              # Table 3: V1 = 12.4 L (RSE 10.6%, 95% CI 9.83-14.9)
    lvp  <- log(52.2);  label("V2: first peripheral volume of distribution (L)")                     # Table 3: V2 = 52.2 L (RSE 6.25%, 95% CI 45.8-58.6)
    lvp2 <- log(6250);  label("V3: second peripheral volume of distribution (L)")                    # Table 3: V3 = 6250 L (RSE 22.9%, 95% CI 3447-9052)
    lq   <- log(137);   label("Q1: intercompartmental clearance V1 <-> V2 (L/h)")                    # Table 3: Q1 = 137 L/h (RSE 6.32%, 95% CI 120-154)
    lq2  <- log(21.1);  label("Q2: intercompartmental clearance V1 <-> V3 (L/h)")                    # Table 3: Q2 = 21.1 L/h (RSE 9.95%, 95% CI 17.0-25.2)

    # ============================================================
    # Inter-individual variability (Muller 2010 Table 3, "Variance
    # model" block). The paper assumed individual PK parameters
    # follow a log-normal distribution ("an exponential distribution
    # model was used to account for interindividual variability"),
    # so the NONMEM omega values are variances on the log scale:
    # P_i = TVP * exp(eta_i), eta_i ~ N(0, omega^2). IIV on CL and V3
    # were not statistically significant in the original analysis
    # (the 95% CIs for both include zero; the paper attributes this
    # to the small cohort, n=7) but are retained because they are
    # part of the published final model. No IIV was reported on V1,
    # V2, Q1, or Q2.
    # ============================================================
    etalcl  ~ 0.292    # Table 3 IIV in CL = 0.292 (RSE 56.2%, 95% CI -0.0294 to 0.613)
    etalvp2 ~ 0.00185  # Table 3 IIV in V3 = 0.00185 (RSE 156%, 95% CI -0.00379 to 0.00749) -- effectively negligible

    # ============================================================
    # IIV-on-residual-error scaling (Muller 2010 Table 3
    # "IIV in error" = 0.306). The paper uses NONMEM's omega-sigma
    # interaction option under FOCE with an additional ETA that
    # scales the residual-error magnitude per subject. Encoded as in
    # Dogterom_2018_asenapine.R: a fixed log-anchor lrv pairs the eta
    # with a typical-value fixed effect, and etalrv carries the
    # reported variance.
    # ============================================================
    lrv     <- fixed(log(1));  label("Residual-variability scaling anchor (fixed log(1))")           # structural anchor; pairs etalrv with a typical-value fixed effect
    etalrv  ~ 0.306             # Table 3 IIV in error = 0.306 (RSE 50.0%, 95% CI 0.00612-0.606)

    # ============================================================
    # Residual variability (Muller 2010 Table 3 "Residual variability"
    # = 0.0424). The paper reports the proportional EPS variance on
    # the linear scale (sigma^2 = 0.0424), so the proportional SD is
    # sqrt(0.0424) = 0.2059. The per-subject SD picks up the etalrv
    # scaling in model(): propSd = propSdBase * exp(lrv + etalrv).
    # ============================================================
    propSdBase <- 0.2059;  label("Proportional residual error baseline magnitude (fraction; sqrt(0.0424))")  # Table 3: residual variability sigma^2 = 0.0424 (RSE 42.5%, 95% CI 0.00712-0.0777); propSdBase = sqrt(0.0424)
  })

  model({
    # ---- Individual PK parameters (exponential IIV on log-scale TVs) ----
    cl  <- exp(lcl + etalcl)
    vc  <- exp(lvc)
    vp  <- exp(lvp)
    vp2 <- exp(lvp2 + etalvp2)
    q   <- exp(lq)
    q2  <- exp(lq2)

    # ---- Micro-constants (h^-1) ----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # ---- 3-compartment IV PK (Muller 2010 NONMEM ADVAN11; Table 3) ----
    # Doses go directly to the central compartment as IV infusions
    # (600 mg over 20 min q6h or 900 mg over 30 min q8h); infusion
    # rate/duration is specified per dose in the user's event table
    # rather than hard-coded here.
    d/dt(central)     <- -(kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # ---- Observation: plasma clindamycin concentration ----
    # Dose in mg, vc in L -> mg/L (the paper's reporting unit; LLOQ
    # 0.1 mg/L; assay linear to 50 mg/L).
    Cc <- central / vc

    # ---- Per-subject proportional residual SD ----
    # Muller 2010 'IIV in error' (Table 3) scales propSdBase via a
    # log-normal eta: propSd = propSdBase * exp(lrv + etalrv). lrv is
    # a fixed log-anchor pairing etalrv with a typical-value fixed
    # effect (matches the Dogterom_2018_asenapine.R precedent).
    propSd <- propSdBase * exp(lrv + etalrv)
    Cc ~ prop(propSd)
  })
}
