Germovsek_2018_meropenem <- function() {
  description <- "One-compartment plasma + CSF (two-state) IV population PK model for meropenem in neonates and young infants (<=90 days) with late-onset sepsis and/or meningitis (Germovsek 2018; NeoMero-1 and NeoMero-2 studies). Plasma CL and Vc are allometrically scaled to body weight (fixed exponent 0.632 on CL, 1.0 on Vc) with a fixed Rhodin-style postmenstrual-age maturation Hill function on CL and a power covariate of (CREAT_REF / CREAT) on CL; an additional CSF compartment with fixed Vcsf = 0.15 L/70 kg and estimated inter-compartmental clearance CL_CSF carries a logit-scale CSF penetration fraction (typical 8.4 %) modulated by CSF total protein concentration."
  reference <- "Germovsek E, Lutsar I, Kipper K, Karlsson MO, Planche T, Chazallon C, Meyer L, Trafojer UMT, Metsvaht T, Fournier I, Sharland M, Heath P, Standing JF; NeoMero Consortium. Plasma and CSF pharmacokinetics of meropenem in neonates and young infants: results from the NeoMero studies. J Antimicrob Chemother. 2018;73(7):1908-1916. doi:10.1093/jac/dky128"
  vignette <- "Germovsek_2018_meropenem"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight at enrolment",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Allometric scaling on CL (fixed exponent 0.632), Vc (fixed exponent 1), CL_CSF (allometric exponent 0.75) and Vcsf (linear) with reference 70 kg per Germovsek 2018 Methods 'PK modelling'.",
      source_name        = "WT"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age + postnatal age)",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Drives the fixed Rhodin-style GFR-maturation Hill function on CL (Tmat50 = 47.7 weeks PMA, Hill = 3.4) per Germovsek 2018 Methods reference 27 (Rhodin et al. 2009).",
      source_name        = "PMA"
    ),
    CREAT = list(
      description        = "Patient measured serum creatinine concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Used together with CREAT_REF in the renal-function factor (CREAT_REF / CREAT)^0.40 on CL. Germovsek 2018 Table 2 reports the exponent with sign convention theta_creatinine = -0.40 applied to (CREAT / CREAT_REF); the encoding here uses the algebraically identical (CREAT_REF / CREAT)^0.40 to match the nlmixr2lib convention used in Hennig 2013 / Llanos-Paez 2020.",
      source_name        = "SCR"
    ),
    CREAT_REF = list(
      description        = "Postmenstrual-age-expected normal-mean serum creatinine",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Externally-computed reference SCR for the individual (denoted standardised SCR in Germovsek 2018; Methods reference 28). For the typical infant in the study (PMA 37.4 weeks, raw SCR 32 umol/L) the standardised SCR was reported as 60 umol/L (Discussion). Users must compute CREAT_REF from the chosen PMA-stratified reference table before passing to the model; for normal renal function set CREAT_REF = CREAT so the renal-function factor evaluates to 1.",
      source_name        = "SCR_standardised"
    ),
    CSF_TPRO = list(
      description        = "Cerebrospinal-fluid total protein concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying CSF total protein measured contemporaneously with the CSF PK sample. Reference value 1.2 g/L (typical infant; Germovsek 2018 Discussion). Drives an additive deviation on the logit-scale CSF uptake parameter (theta_CSFproteins = -0.17; Table 2). Missing values were imputed to the median during model fitting (Germovsek 2018 Methods 'PK modelling').",
      source_name        = "CSF_protein"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 167L,
    n_studies      = 2L,
    age_range      = "PNA 1-90 days; PMA 22.6-51.3 weeks (gestational + postnatal)",
    age_median     = "PNA 13 days; PMA 37.4 weeks",
    weight_range   = "0.48-6.32 kg",
    weight_median  = "2.12 kg",
    sex_female_pct = 46.7,
    race_ethnicity = "Not reported (multicentre European neonatal ICU cohort)",
    disease_state  = "Neonates and young infants (<=90 days) with suspected or confirmed late-onset sepsis (NeoMero-1, n = 123 including 5 transferred to NeoMero-2) or bacterial meningitis (NeoMero-2, n = 49 including the 5 transferred). Exclusion criteria: renal failure, severe congenital malformations, causative pathogen resistant to meropenem, known intolerance or contraindications.",
    dose_range     = "Meropenem 20 mg/kg q12h for LOS in <32 weeks GA and <2 weeks PNA, 20 mg/kg q8h for LOS in all others, 40 mg/kg q8h or q12h for meningitis. 30-minute IV infusion in all subjects.",
    gestational_age_range = "22.6-41.9 weeks GA at birth (median 33.3 weeks)",
    postmenstrual_age_range = "22.6-51.3 weeks PMA (median 37.4 weeks)",
    samples_plasma = "401 plasma samples (median 2.4 per patient; 109 patients contributed 3 optimally timed samples, 44 patients contributed a single trough)",
    samples_csf    = "78 CSF samples from 56 patients (median 0.47 per patient; collected opportunistically at lumbar puncture, median 5.27 h post-dose)",
    regions        = "Europe (multicentre: UK, Italy, Estonia, France, plus sites listed in NeoMero Consortium members)",
    notes          = "Demographics from Germovsek 2018 Table 1. Five infants switched from NeoMero-1 to NeoMero-2 on later meningitis diagnosis. Eleven peak plasma samples below 10 mg/L were excluded as suspected data-entry errors and one biologically implausible CSF protein concentration (102 g/L) was excluded. Median raw SCR 32 umol/L (range 3.54-197.4); median CSF protein 1.2 g/L; median CSF lactate 1.8 mmol/L."
  )

  ini({
    # Structural plasma parameters (Germovsek 2018 Table 2 final-model column).
    lcl  <- log(16.7); label("Plasma clearance standardised to 70 kg (L/h)")            # Germovsek 2018 Table 2: CL = 16.7 L/h/70 kg
    lvc  <- log(38.6); label("Central volume of distribution standardised to 70 kg (L)") # Germovsek 2018 Table 2: V  = 38.6 L/70 kg

    # CSF-compartment parameters.
    lcl_csf <- log(0.017); label("Inter-compartmental clearance plasma <-> CSF standardised to 70 kg (L/h)") # Germovsek 2018 Table 2: CL_CSF = 0.017 L/h/70 kg

    # CSF uptake on logit scale. logituptake = 2.39 corresponds to a
    # barrier fraction of expit(2.39) = 0.916 and therefore a CSF
    # penetration fraction p = 1 - 0.916 = 0.084 (= 8.4 % typical;
    # Germovsek 2018 Abstract and Discussion).
    logituptake <- 2.39; label("Logit of the CSF barrier fraction (1 - penetration); typical penetration p = 1 - expit(logituptake)")  # Germovsek 2018 Table 2: CSF uptake = 2.39 (logit scale)

    # Allometric exponents (Germovsek 2018 Methods 'PK modelling': "weight
    # was included with allometric scaling (with exponents fixed to 1 for
    # central volume and 0.632 for CL)").
    e_wt_cl  <- fixed(0.632); label("Allometric exponent on CL (unitless)")                    # Germovsek 2018 Methods: fixed at 0.632
    e_wt_vc  <- fixed(1.000); label("Allometric exponent on Vc (unitless)")                    # Germovsek 2018 Methods: fixed at 1
    e_wt_clcsf <- fixed(0.75); label("Allometric exponent on CL_CSF (unitless)")               # standard small-molecule renal-CL allometry assumed
    e_wt_vcsf  <- fixed(1.000); label("Allometric exponent on Vcsf (unitless)")                # standard volume allometry

    # Rhodin renal-maturation Hill function (Germovsek 2018 Methods
    # reference 27): fmat = PMA^Hill / (Tmat50^Hill + PMA^Hill). Both
    # parameters fixed to the cited reference values (Rhodin et al. 2009
    # Pediatr Nephrol 24:67-76; widely used for renal-clearance maturation
    # in neonates).
    tmat50 <- fixed(47.7); label("PMA at 50 % renal maturation (weeks)")  # Germovsek 2018 Methods reference 27 (Rhodin 2009)
    hill_mat <- fixed(3.4); label("Hill coefficient for renal maturation (unitless)")  # Germovsek 2018 Methods reference 27 (Rhodin 2009)

    # Renal-function covariate on CL (Germovsek 2018 Table 2). Paper
    # writes the exponent as -0.40 with convention (CREAT / CREAT_REF)^theta;
    # encoded here as +0.40 on the equivalent (CREAT_REF / CREAT) ratio to
    # match the Hennig 2013 / Llanos-Paez 2020 nlmixr2lib pattern.
    e_creat_cl <- 0.40; label("Exponent on (CREAT_REF / CREAT) ratio for CL (unitless)")  # Germovsek 2018 Table 2: theta_creatinine = -0.40 on (CREAT/CREAT_REF)

    # CSF protein covariate on the logit barrier (Germovsek 2018 Table 2).
    # Enters additively on the logit scale; sign convention is that
    # higher CSF protein lowers the barrier (i.e. increases penetration)
    # so the encoded value carries the published -0.17 magnitude.
    e_csftpro_uptake <- -0.17; label("CSF protein effect on logit CSF barrier (per g/L deviation from 1.2)")  # Germovsek 2018 Table 2: theta_CSFproteins = -0.17 (logit scale)

    # Inter-individual variability. Germovsek 2018 Table 2 reports
    # variances (omega^2) directly (the "IIV on CL = 0.255" / "IIV on
    # V = 0.153" / "Cov IIV CL-V = 0.167" entries are variance / covariance
    # on the log-normal eta scale).
    etalcl + etalvc ~ c(0.255,
                        0.167, 0.153)  # Germovsek 2018 Table 2: var(CL) = 0.255, var(V) = 0.153, cov(CL,V) = 0.167

    # Residual error. Germovsek 2018 used the Dosne/Karlsson Box-Cox
    # dynamic-transform-both-sides (dTBS) residual error model with
    # estimated shape (lambda) and scedasticity (delta) parameters:
    #   Plasma: lambda = 0.280, delta = -0.174, RUV SD = 0.679
    #   CSF   : lambda = 0.285, delta = 0 (fixed), RUV SD = 1.19
    # nlmixr2 has no native dTBS residual-error specifier; this
    # implementation encodes the residual as a standard proportional
    # error using the paper's RUV SDs. The Box-Cox dTBS shape (lambda)
    # close to zero approximates the log-transform-both-sides limit, so
    # the proportional form is a reasonable but not identical surrogate.
    # See the vignette's "Assumptions and deviations" section for the
    # implication for stochastic VPCs.
    propSd      <- 0.679; label("Plasma proportional residual SD (Box-Cox dTBS SD; see Assumptions)")  # Germovsek 2018 Table 2: RUV_plasma = 0.679
    propSd_Ccsf <- 1.19;  label("CSF proportional residual SD (Box-Cox dTBS SD; see Assumptions)")     # Germovsek 2018 Table 2: RUV_CSF    = 1.19
  })

  model({
    # Declare named compartments for both ODE states and the algebraic
    # PK observables (Cc, Ccsf) so event tables can reference them by
    # name; without this, rxode2's auto-injected cmt(Cc)/cmt(Ccsf) push
    # central/csf into observation slots and the solver complains
    # "csf is required for solving".
    cmt(central)
    cmt(csf)
    cmt(Cc)
    cmt(Ccsf)

    # Renal-function maturation (Rhodin et al. 2009 Hill function on PMA).
    fmat <- PAGE^hill_mat / (tmat50^hill_mat + PAGE^hill_mat)

    # SCR ratio: paper standardises raw SCR by PMA-expected reference and
    # then enters the standardised value into a power form on CL.
    # Equivalent to (CREAT_REF / CREAT)^e_creat_cl with e_creat_cl > 0.
    scr_factor <- (CREAT_REF / CREAT)^e_creat_cl

    # Plasma individual PK parameters (allometric size + maturation + SCR).
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * fmat * scr_factor
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    # CSF compartment parameters (Germovsek 2018 Methods: CSF volume fixed
    # to 0.15 L/70 kg per cited reference 32; CL_CSF estimated).
    cl_csf <- exp(lcl_csf) * (WT / 70)^e_wt_clcsf
    vcsf   <- 0.15 * (WT / 70)^e_wt_vcsf

    # CSF barrier / penetration fraction. The paper estimates the logit
    # of the BARRIER fraction (1 - p); additive CSF-protein deviation
    # from the 1.2 g/L reference modulates this logit. p = 1 - barrier
    # is the steady-state CSF / plasma concentration ratio.
    logit_barr <- logituptake + e_csftpro_uptake * (CSF_TPRO - 1.2)
    barrier    <- 1 / (1 + exp(-logit_barr))
    penetration <- 1 - barrier

    # Plasma + CSF ODE system. CL_CSF * Cc enters the CSF compartment;
    # the CSF efflux rate constant is scaled by 1 / penetration so that
    # the steady-state Ccsf / Cc ratio equals the penetration fraction
    # (Germovsek 2018 Methods 'PK modelling').
    kel <- cl / vc
    d/dt(central) <- -kel * central -
                       (cl_csf / vc) * central +
                       (cl_csf / (penetration * vcsf)) * csf
    d/dt(csf)     <-   (cl_csf / vc) * central -
                       (cl_csf / (penetration * vcsf)) * csf

    # Observation variables. Dose in mg, V in L -> concentrations in mg/L.
    Cc   <- central / vc
    Ccsf <- csf     / vcsf

    Cc   ~ prop(propSd)
    Ccsf ~ prop(propSd_Ccsf)
  })
}
