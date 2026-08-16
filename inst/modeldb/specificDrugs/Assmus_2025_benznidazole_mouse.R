Assmus_2025_benznidazole_mouse <- function() {
  description <- paste(
    "Preclinical (mouse). One-compartment population PK model with",
    "first-order absorption for benznidazole in female BALB/c mice, with",
    "a priori allometric scaling on CL/F (exponent 0.75) and Vc/F",
    "(exponent 1.0) centred on 19.4 g, and a power-function effect of the",
    "administered mg/kg dose on the absorption rate constant (KA falls",
    "from 5.11 to 0.86 1/h between 10 and 100 mg/kg). The model carries",
    "two PK/PD index accumulator states -- cumulative plasma AUC and",
    "cumulative time above the protein-binding-corrected in vitro IC90",
    "against Trypanosoma cruzi amastigotes (6.427 ug/mL) -- which drive",
    "the paper's binary univariate logistic exposure-response models for",
    "sterile parasitological cure in chronically infected mice, measured",
    "by in vivo and ex vivo bioluminescence imaging.",
    sep = " "
  )
  reference <- paste(
    "Assmus F, Adehin A, Hoglund RM, Fortes Francisco A, Lewis MD,",
    "Kelly JM, Charman SA, White KL, Shackleford DM, Escudie F,",
    "Chatelain E, Scandale I, Tarning J.",
    "Pharmacokinetic-pharmacodynamic modeling of benznidazole and its",
    "antitrypanosomal activity in a murine model of chronic Chagas",
    "disease. PLoS Negl Trop Dis. 2025;19(5):e0012968.",
    "doi:10.1371/journal.pntd.0012968.",
    sep = " "
  )
  vignette <- "Assmus_2025_benznidazole_mouse"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ug/mL"
  )

  # Volumes are carried in L and clearances in L/h (the paper reports
  # 19.7 mL and 10.6 mL/h; see the ini() comments for the 1000-fold
  # conversion) so that a dose in mg gives a concentration in mg/L, which
  # is numerically identical to the ug/mL used throughout the paper.

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline weight at the start of dosing. Reference weight 0.0194 kg",
        "(19.4 g), the median weight of the 52 satellite PK mice (Results,",
        "'Pharmacokinetics of benznidazole'); the S1 Code $PK block writes",
        "this as (WT/0.0194) with WT in kg. Allometric exponents were fixed",
        "a priori at 0.75 on CL/F and 1.0 on Vc/F rather than estimated,",
        "because the satellite-PK weight range (18.1-21.6 g) is too narrow",
        "to support estimation (delta-OFV of 1.1); the authors retained them",
        "to allow simulation at the heavier 25 g median weight of the",
        "efficacy cohort.",
        sep = " "
      ),
      source_name        = "WT"
    ),
    DOSE_BZN_MGKG = list(
      description        = "Administered benznidazole dose per administration, in mg per kg body weight",
      units              = "mg/kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Carried on every record (not only dose rows) so the power",
        "covariate function on KA can be evaluated at observation times.",
        "Enters as (DOSE_BZN_MGKG / 30)^e_dose_ka with the reference set to",
        "the 30 mg/kg median dose of the satellite PK study (Methods,",
        "'Population pharmacokinetic analysis'). Dose levels studied:",
        "10, 30 and 100 mg/kg in the satellite PK study; 10, 20, 30, 50 and",
        "100 mg/kg per administration in the efficacy studies.",
        sep = " "
      ),
      source_name        = "DOSE_PER_KG"
    )
  )

  compartmentData <- list(
    depot        = list(analyte = "benznidazole", units = "mg", specimen = "administration site", verified = TRUE),
    central      = list(analyte = "benznidazole", units = "mg", specimen = "plasma", verified = TRUE),
    auc_total    = list(analyte = "benznidazole", units = NA_character_, specimen = "not applicable", verified = TRUE),
    t_above_mic  = list(analyte = "benznidazole", units = NA_character_, specimen = "not applicable", verified = TRUE)
  )

  # Bookkeeping accumulators for the two PK/PD indices the exposure-
  # response analysis is built on; same construction as the
  # Lallemand_2023_benzylpenicillin_horse.R auc_free / t_above_mic pair.
  paper_specific_compartments <- c("auc_total", "t_above_mic")

  population <- list(
    species        = "mouse (BALB/c, female)",
    n_subjects     = 52L,
    n_studies      = 1L,
    age_range      = "7-8 weeks old at infection (efficacy cohort)",
    weight_range   = "18.1-21.6 g (satellite PK cohort)",
    weight_median  = "19.4 g (satellite PK cohort); 25 g (efficacy cohort)",
    sex_female_pct = 100,
    disease_state  = paste(
      "Population PK was estimated in 52 UNINFECTED female BALB/c satellite",
      "mice. The exposure-response layer was fitted in a separate cohort of",
      "118 female BALB/c mice with chronic-stage Trypanosoma cruzi CL Brener",
      "infection (1000 bioluminescent blood trypomastigotes, chronic stage",
      "reached 50-70 days post-infection), weighing 20.5-28.1 g (median 25 g)",
      "at the start of dosing.",
      sep = " "
    ),
    dose_range     = paste(
      "Satellite PK: single oral gavage doses of 10 mg/kg (n = 16),",
      "30 mg/kg (n = 18) and 100 mg/kg (n = 18). Efficacy studies:",
      "10, 20, 30, 50 and 100 mg/kg per administration, once daily",
      "(one 50 mg/kg twice-daily regimen), for 5, 10 or 20 days",
      "(10 regimens, Table 1).",
      sep = " "
    ),
    regions        = "Not applicable (preclinical; London School of Hygiene and Tropical Medicine efficacy studies, Monash University bioanalysis)",
    notes          = paste(
      "PK sampling was sparse (1-3 samples per mouse) at 0.25, 0.5, 1, 2, 4,",
      "4.5, 6, 8, 10, 12, 14, 16 and 24 h post-dose; 110 samples collected,",
      "2 outliers excluded, data censored beyond 8 h (10 mg/kg) and 12 h",
      "(30 mg/kg) because all later levels were below the 5 ng/mL LLOQ,",
      "leaving 93 samples of which 90 were above LLOQ and entered the final",
      "fit. Estimated in NONMEM 7.4 with FOCE-I on natural-log-transformed",
      "concentrations; parameter precision from sampling importance",
      "resampling. Efficacy cohort: 83 of 118 mice cured (Table 1);",
      "parasitological cure defined as absence of bioluminescence on both",
      "in vivo and ex vivo imaging after cyclophosphamide immunosuppression.",
      sep = " "
    )
  )

  ini({
    # -------------------------------------------------------------------
    # Structural PK. Population estimates are Table 2 and are reproduced
    # verbatim in the S1 Code $THETA block, which confirms the
    # parameterisation:
    #   TVCL = THETA(1) * (WT/0.0194)**0.75
    #   TVV2 = THETA(2) * (WT/0.0194)**1.00
    #   TVKA = THETA(3) * (DOSE_PER_KG/30)**THETA(5)
    #   F1   = THETA(4)                                   [1 FIX]
    # Volumes/clearances are converted from the paper's mL to L so that a
    # dose in mg yields a concentration in mg/L = ug/mL.
    # -------------------------------------------------------------------
    lka <- log(2.18);    label("Absorption rate constant at the 30 mg/kg reference dose (1/h)")  # Table 2: K_A = 2.18 1/h (RSE 15.9 pct, 95 pct CI 1.65-3.03); S1 Code $THETA 3
    lcl <- log(0.0106);  label("Apparent clearance CL/F at 0.0194 kg body weight (L/h)")         # Table 2: CL/F = 10.6 mL/h (RSE 5.6 pct, 95 pct CI 9.62-11.88) = 0.0106 L/h; S1 Code $THETA 1
    lvc <- log(0.0197);  label("Apparent central volume Vc/F at 0.0194 kg body weight (L)")      # Table 2: V_C/F = 19.7 mL (RSE 8.2 pct, 95 pct CI 17.05-23.18) = 0.0197 L; S1 Code $THETA 2

    # Relative bioavailability was fixed to unity for the population so
    # that CL and Vc are apparent (CL/F, Vc/F) values.
    lfdepot <- fixed(log(1)); label("Relative oral bioavailability (fraction)")                  # S1 Code $THETA 4: "1 FIX"; Methods: "Relative bioavailability (F, fixed to unity for the population)"

    # -------------------------------------------------------------------
    # Allometry, applied a priori (not estimated). Methods: "Body weight
    # was implemented a priori as an allometric function on clearance
    # (exponent 0.75) and volume of distribution (exponent 1), centered on
    # the median weight of mice in the satellite PK study (19.4 g)."
    # -------------------------------------------------------------------
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL/F (unitless)")   # Methods, 'Population pharmacokinetic analysis'; S1 Code TVCL exponent 0.75
    e_wt_vc <- fixed(1.0);  label("Allometric exponent on Vc/F (unitless)")   # Methods, 'Population pharmacokinetic analysis'; S1 Code TVV2 exponent 1.00

    # -------------------------------------------------------------------
    # Dose effect on absorption (the only covariate retained in the final
    # model). Power function on the mg/kg dose centred on the 30 mg/kg
    # median. The negative exponent reproduces the paper's finding that a
    # higher dose slows absorption: KA = 5.11, 2.18 and 0.86 1/h at 10, 30
    # and 100 mg/kg respectively (Results, 'Pharmacokinetics of
    # benznidazole'). A dose effect on CL/F was screened but NOT retained.
    # -------------------------------------------------------------------
    e_dose_ka <- -0.775; label("Power exponent of (dose / 30 mg/kg) on the absorption rate constant (unitless)")  # Table 2: theta_Dose = -0.775 (RSE 13.3 pct, 95 pct CI -0.977 to -0.583); S1 Code $THETA 5

    # -------------------------------------------------------------------
    # Inter-individual variability. S1 Code $OMEGA gives the variances
    # directly: 0.0298 (CL) and 0.0271 (Vc). These back-transform to the
    # %CV reported in Table 2 via CV = sqrt(exp(omega^2) - 1):
    #   sqrt(exp(0.0298) - 1) = 0.1739 -> 17.4 %CV  (Table 2: 17.4)
    #   sqrt(exp(0.0271) - 1) = 0.1657 -> 16.6 %CV  (Table 2: 16.6)
    # IIV on KA and on F was estimated at <10 %CV and fixed to zero in the
    # final model ("estimated IIV below 10% was fixed to zero in the final
    # model"; S1 Code $OMEGA entries 3 and 4 are "0 FIX"). Those two etas
    # are omitted here rather than written as `~ fixed(0)`, because a
    # zero-variance diagonal makes OMEGA singular and rxode2 then fails to
    # Cholesky-decompose it at simulation time; a fixed-zero variance is
    # mathematically identical to having no eta at all.
    # -------------------------------------------------------------------
    etalcl ~ 0.0298  # S1 Code $OMEGA 1 (IIV_CL); Table 2: 17.4 %CV (RSE 17.0 pct, 95 pct CI 10.5-22.5)
    etalvc ~ 0.0271  # S1 Code $OMEGA 2 (IIV_V2); Table 2: 16.6 %CV (RSE 29.7 pct, 95 pct CI 6.9-25.6)

    # -------------------------------------------------------------------
    # Residual error. Methods: "Residual unexplained variability was
    # implemented as an additive error on the log-transformed observed
    # concentrations (equivalent to an exponential residual error on an
    # arithmetic scale)", i.e. S1 Code $ERROR Y = LOG(IPRED) + EPS(1).
    # NONMEM's $SIGMA holds a VARIANCE, so the 0.122 in Table 2 (labelled
    # "sigma") and in S1 Code $SIGMA is sigma^2; the log-scale SD that
    # nlmixr2's lnorm() error model expects is sqrt(0.122) = 0.3493.
    # -------------------------------------------------------------------
    expSd <- sqrt(0.122); label("Residual SD on the natural-log concentration scale (log ug/mL)")  # Table 2: sigma = 0.122 (RSE 10.9 pct, 95 pct CI 0.087-0.188); S1 Code $SIGMA 0.122 (a variance)

    # -------------------------------------------------------------------
    # Antitrypanosomal target concentration. In vitro IC90 of benznidazole
    # against T. cruzi amastigotes (Tulahuen strain) in 3T3 host cells,
    # corrected for protein binding to the assay medium and to mouse
    # plasma (S1 Text). Held fixed: it is a measured in vitro property,
    # not an estimated model parameter. The paper's sensitivity analysis
    # sweeps it over a 10-fold range (S4/S5 Tables), which a user can
    # reproduce with params = c(mic = ...).
    # -------------------------------------------------------------------
    mic <- fixed(6.427); label("Protein-binding-corrected plasma IC90 against T. cruzi amastigotes (ug/mL)")  # Table 1 footnote: "T>IC90, Time above IC90 in plasma (6.427 ug/mL...)"; Methods gives 24.7 uM = 6.43 ug/mL

    # -------------------------------------------------------------------
    # Binary univariate logistic exposure-response for sterile
    # parasitological cure (Equation 1: logit(p) = beta0 + beta1 * x).
    # Fitted in R by binary logistic regression of per-mouse cure status
    # (83/118 cured) on the regimen-level MEDIAN simulated exposure. Table 4
    # reports six such regressions; two of them are ODE-native and are
    # encoded here. The other four (cumulative dose, CMAX, AUC12, AUC24) are
    # design or summary NCA quantities with no natural ODE state, and are
    # reproduced in the validation vignette from the same Table 4
    # coefficients.
    # -------------------------------------------------------------------
    logite0_cure_auc <- -2.170;  label("Logit of the cure probability at zero cumulative AUC (unitless)")            # Table 4, AUCinf column: Intercept = -2.170 (SE 0.566)
    e_auc_cure       <- 0.0033;  label("Log-odds of parasitological cure per ug*h/mL of cumulative plasma AUC")      # Table 4, AUCinf column: LogOdds = 0.0033 (SE 0.0006); odds increase 0.33 pct per unit
    logite0_cure_tmic <- -1.859; label("Logit of the cure probability at zero time above IC90 (unitless)")           # Table 4, T>IC90 column: Intercept = -1.859 (SE 0.521)
    e_tmic_cure      <- 1.497;   label("Log-odds of parasitological cure per day above the plasma IC90")             # Table 4, T>IC90 column: LogOdds = 1.497 (SE 0.289); odds increase 347 pct per day
  })

  model({
    # 1. Individual parameters. Allometry is centred on the 19.4 g median
    #    satellite-PK mouse; the dose effect on KA is centred on the
    #    30 mg/kg median satellite-PK dose.
    ka <- exp(lka) * (DOSE_BZN_MGKG / 30)^e_dose_ka
    cl <- exp(lcl + etalcl) * (WT / 0.0194)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 0.0194)^e_wt_vc

    # 2. Micro-constant
    kel <- cl / vc

    # 3. One-compartment disposition with first-order absorption
    #    (S1 Code: ADVAN5 TRANS1, COMP(1) absorption -> COMP(2) central,
    #    K12 = KA, K20 = CL/V2).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 4. Bioavailability (fixed to unity for the population)
    f(depot) <- exp(lfdepot)

    # 5. Observation
    Cc <- central / vc

    # 6. PK/PD index accumulators. `auc_total` integrates the total plasma
    #    concentration, so at the end of a treatment course it equals the
    #    paper's cumulative AUCinf; `t_above_mic` accumulates the time the
    #    concentration is at or above the plasma IC90 and is divided by 24
    #    so it is reported in DAYS, matching Table 1 and Table 4.
    d/dt(auc_total)   <- Cc
    d/dt(t_above_mic) <- (Cc >= mic) / 24

    # 7. Exposure-response (Equation 1). Two parallel univariate logistic
    #    models for the probability of sterile parasitological cure.
    pcure_auc  <- 1 / (1 + exp(-(logite0_cure_auc  + e_auc_cure  * auc_total)))
    pcure_tmic <- 1 / (1 + exp(-(logite0_cure_tmic + e_tmic_cure * t_above_mic)))

    # 8. Residual error: additive on the natural-log concentration scale.
    Cc ~ lnorm(expSd)
  })
}
