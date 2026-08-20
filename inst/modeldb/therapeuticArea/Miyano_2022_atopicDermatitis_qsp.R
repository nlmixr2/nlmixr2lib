Miyano_2022_atopicDermatitis_qsp <- function() {
  description <- paste0(
    "QSP + MBMA. Systems-level mathematical model of atopic dermatitis (AD) ",
    "pathogenesis with drug-effect modules for nine biologics (dupilumab, ",
    "lebrikizumab, tralokinumab, secukinumab, fezakinumab, nemolizumab, ",
    "tezepelumab, GBR 830, and recombinant IFN-gamma). Fourteen ODEs describe ",
    "skin barrier integrity, infiltrated pathogens, four helper T-cell subsets ",
    "(Th1, Th2, Th17, Th22), seven cytokines (IL-4, IL-13, IL-17A, IL-22, ",
    "IL-31, IFN-gamma, TSLP), and OX40L; the EASI efficacy endpoint is an ",
    "algebraic function of skin barrier and pathogens. 51 log-normally ",
    "distributed patient-specific parameters define virtual-patient ",
    "heterogeneity; model-based meta-analysis across 9 drug trials tuned the ",
    "51 mu and 51 sigma distribution parameters. Drug effects enter as ",
    "multiplicative inhibition of the effective cytokine concentration ",
    "(r_inhibit = 0.99 for antibodies, or 0.99 x e_a2 = 0.436 for tralokinumab) ",
    "or as an additive rIFNg concentration (c_rifng = 210). Time in weeks."
  )
  reference <- paste(
    "Miyano T, Irvine AD, Tanaka RJ.",
    "A mathematical model to identify optimal combinations of drug targets",
    "for dupilumab poor responders in atopic dermatitis.",
    "Allergy. 2022;77(2):582-594.",
    "doi:10.1111/all.14870.",
    "PMID: 33894014.",
    "Model equations and parameter values (mu.csv, sigma.csv) reproduced",
    "from the authors' MIT-licensed reference implementation at",
    "https://github.com/Tanaka-Group/AD_QSP_model (cited in paper Section 4.6)."
  )
  vignette <- "Miyano_2022_atopicDermatitis_qsp"

  # 14 ODE states, all paper-mechanistic (no canonical PK compartment role).
  # AD_QSP_model.py diff_eq() docstring lines 68-83.
  paper_specific_compartments <- c(
    "barrier", "pathogens",
    "th1", "th2", "th17", "th22",
    "il4", "il13", "il17", "il22", "il31", "ifng", "tslp", "ox40l"
  )

  units <- list(
    time          = "week",
    dosing        = "n/a (drug effects enter algebraically via inhibition fractions and additive rIFNg concentration; no explicit PK compartment)",
    concentration = "fold-change vs healthy skin (cytokines, OX40L); count fold-change (T cells); unitless [0,1] fraction (skin barrier); unitless fold-change (pathogens); 0-72 (EASI score)"
  )

  # No subject-level covariates. Drug-scenario configuration is exposed as
  # ini() parameters (r_placebo, r_il4, ... c_rifng, e_a2) that a user
  # overrides at rxSolve() time via params = c(...); virtual-patient
  # heterogeneity is captured entirely by the 51 etas on the mechanistic
  # rate constants.
  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    barrier   = list(analyte = "skin barrier integrity", units = NA_character_, specimen = "not applicable", verified = FALSE),
    pathogens = list(analyte = "infiltrated pathogens", units = NA_character_, specimen = "not applicable", verified = FALSE),
    th1       = list(analyte = "Th1 helper T-cells", units = NA_character_, specimen = "blood cell", verified = FALSE),
    th2       = list(analyte = "Th2 helper T-cells", units = NA_character_, specimen = "blood cell", verified = FALSE),
    th17      = list(analyte = "Th17 helper T-cells", units = NA_character_, specimen = "blood cell", verified = FALSE),
    th22      = list(analyte = "Th22 helper T-cells", units = NA_character_, specimen = "blood cell", verified = FALSE),
    il4       = list(analyte = "IL-4 cytokine", units = NA_character_, specimen = "serum", verified = FALSE),
    il13      = list(analyte = "IL-13 cytokine", units = NA_character_, specimen = "serum", verified = FALSE),
    il17      = list(analyte = "IL-17A cytokine", units = NA_character_, specimen = "serum", verified = FALSE),
    il22      = list(analyte = "IL-22 cytokine", units = NA_character_, specimen = "serum", verified = FALSE),
    il31      = list(analyte = "IL-31 cytokine", units = NA_character_, specimen = "serum", verified = FALSE),
    ifng      = list(analyte = "IFN-gamma cytokine", units = NA_character_, specimen = "serum", verified = FALSE),
    tslp      = list(analyte = "TSLP cytokine", units = NA_character_, specimen = "serum", verified = FALSE),
    ox40l     = list(analyte = "OX40L cytokine", units = NA_character_, specimen = "serum", verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = "1000 virtual patients per simulation; MBMA pool 663 placebo + 848 drug-arm patients across 9 published trials (Miyano 2022 Table 1)",
    n_studies      = 9L,
    age_range      = "adults with moderate-to-severe atopic dermatitis (MBMA source trials, e.g. Simpson 2016 dupilumab Ph3 adults >=18 years)",
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = "Moderate-to-severe atopic dermatitis (MBMA of biologics-treatment trials; virtual patients defined only by pathophysiologic parameter heterogeneity, not demographics)",
    dose_range     = paste0(
      "Highest-dose arm per drug (Miyano 2022 Table 1): dupilumab 300 mg qw SC + TCS; ",
      "lebrikizumab 250 mg q2w SC + TCS; tralokinumab 300 mg q2w SC + TCS; ",
      "secukinumab 300 mg qw x4 then 300 mg q4w; fezakinumab 600 mg d0 then 300 mg q2w IV; ",
      "nemolizumab 60 mg q4w SC; tezepelumab 280 mg q2w SC + TCS; ",
      "GBR 830 10 mg/kg q4w IV; rIFNg 50 ug/m^2 qd SC + TCS."
    ),
    regions        = "multi-national (source-trial dependent)",
    notes          = paste0(
      "Model-based meta-analysis pools 9 randomised placebo-controlled trials ",
      "of investigational biologics for AD (Miyano 2022 Table 1). The QSP layer ",
      "describes system-level AD pathogenesis independently of demographics; ",
      "the MBMA layer tunes 51 log-normal parameter distributions to reproduce ",
      "the mean EASI score and EASI-75 responder rate at 24 weeks for each drug. ",
      "The model reproduces reported baseline levels of biological factors ",
      "(Miyano 2022 Figure 3A) as steady-state solutions and reproduces the ",
      "mean EASI trajectories and EASI-75 rates for the 9 drugs ",
      "(Miyano 2022 Figure 3B/3C, RMSE 2.1 for mean EASI and 7.4% for EASI-75)."
    )
  )

  # ------------------------------------------------------------------------
  # PARAMETER PROVENANCE
  # ------------------------------------------------------------------------
  # All 51 mu values (log-scale means) come from Tanaka-Group/AD_QSP_model
  # mu.csv (indices 0..50 in row order). All 51 sigma values (log-scale
  # standard deviations) come from sigma.csv (same indexing). Parameter
  # name -> index -> mu -> sigma mapping is preserved in the trailing
  # in-file comments so a downstream reader can source-trace every value
  # in one line.
  #
  # The 11 elimination-rate parameters (d1 skin barrier, d8 pathogens,
  # d9 Th cells, d10 IL-4, d11 IL-13, d12 IL-17A, d13 IL-22, d14 IL-31,
  # d15 IFN-gamma, d16 TSLP, d17 OX40L) were determined from published in
  # vivo (serum) human half-lives (Miyano 2022 Methods section 2.4, Table
  # S5 -- S5 is in the paper supplement which is not on disk; the numerical
  # values here are reproduced from the authors' companion mu.csv).
  #
  # The 40 remaining parameters were tuned to reproduce mean and CV of
  # baseline biological factor levels (Table S2) plus mean EASI and EASI-75
  # of each drug (Miyano 2022 Figure 1). See paper Methods section 2.4.
  #
  # ini() encodes the 51 parameters as l<name> log-scale THETAs at mu, with
  # eta<lname> ~ sigma^2 log-normal IIV. Drug-scenario overrides are exposed
  # as separate ini() parameters at neutral defaults (see block after the
  # structural parameters).
  ini({
    # === Skin barrier integrity (dc0/dt) -- Eq. from AD_QSP_model.py L155 ===
    lk1 <- 0.537   ; label("Log recovery rate of skin barrier via skin turnover (log 1/week)")               # mu.csv row 1
    lk2 <- -1.5    ; label("Log strength of IL-22 promotion of skin barrier recovery (log unitless)")        # mu.csv row 2
    lk3 <- 2.95    ; label("Log placebo/trial-participation skin-barrier recovery rate (log 1/week)")        # mu.csv row 3
    lb1 <- -8.69   ; label("Log strength of IL-4 inhibition of skin barrier recovery (log 1/cytokine)")     # mu.csv row 4
    lb2 <- -3.54   ; label("Log strength of IL-13 inhibition of skin barrier recovery (log 1/cytokine)")    # mu.csv row 5
    lb3 <- -3.92   ; label("Log strength of IL-17A inhibition of skin barrier recovery (log 1/cytokine)")   # mu.csv row 6
    lb4 <- -0.595  ; label("Log strength of IL-22 inhibition of skin barrier recovery (log 1/cytokine)")    # mu.csv row 7
    lb5 <- -2.09   ; label("Log strength of IL-31 inhibition of skin barrier recovery (log 1/cytokine)")    # mu.csv row 8
    ld1 <- -1.599  ; label("Log skin barrier self-degradation rate (log 1/week; elimination-rate parameter)")   # mu.csv row 9 (Table S5)
    ld2 <- -2.64   ; label("Log strength of IL-31 promotion of skin barrier degradation (log 1/cytokine/week)")  # mu.csv row 10
    ld3 <- 1.25    ; label("Log strength of infiltrated pathogens on skin barrier degradation (log 1/pathogen)") # mu.csv row 11

    # === Infiltrated pathogens (dc1/dt) -- Eq. from AD_QSP_model.py L158 ===
    lb6 <- 0.395   ; label("Log strength of skin barrier repression of pathogen infiltration (log 1/barrier)")   # mu.csv row 12
    ld4 <- -1      ; label("Log strength of pathogen self-repression of pathogen loss (log 1/pathogen)")         # mu.csv row 13
    ld5 <- -5.1    ; label("Log strength of IL-17A promotion of pathogen loss (log 1/cytokine)")                 # mu.csv row 14
    ld6 <- -5.01   ; label("Log strength of IL-22 promotion of pathogen loss (log 1/cytokine)")                  # mu.csv row 15
    ld7 <- -8.7    ; label("Log strength of IFN-gamma promotion of pathogen loss (log 1/cytokine)")              # mu.csv row 16
    lb7 <- -7.99   ; label("Log strength of IL-4 inhibition of pathogen loss (log 1/cytokine)")                  # mu.csv row 17
    lb8 <- -3.5    ; label("Log strength of IL-13 inhibition of pathogen loss (log 1/cytokine)")                 # mu.csv row 18
    ld8 <- -1.599  ; label("Log pathogen self-loss rate = pathogen production rate k4 (log 1/week; elim-rate)")  # mu.csv row 19 (Table S5)

    # === Th cell subsets (dc2..dc5/dt) -- Eq. from AD_QSP_model.py L161-164 ===
    lk5  <- 2.95   ; label("Log rate constant for pathogen-driven Th1 differentiation (log 1/week/pathogen)")   # mu.csv row 20
    lk9  <- -3.32  ; label("Log strength of IFN-gamma promotion of Th1 differentiation (log 1/cytokine)")       # mu.csv row 21
    ld9  <- -0.5482; label("Log Th-cell elimination rate (shared across Th1/Th2/Th17/Th22; log 1/week; elim)")  # mu.csv row 22 (Table S5)
    lb9  <- -2.66  ; label("Log strength of OX40L inhibition of Th cell elimination (log 1/OX40L)")             # mu.csv row 23
    lk6  <- 3.89   ; label("Log rate constant for pathogen-driven Th2 differentiation (log 1/week/pathogen)")   # mu.csv row 24
    lk10 <- -5.82  ; label("Log strength of IL-4 promotion of Th2 differentiation (log 1/cytokine)")            # mu.csv row 25
    lk7  <- 2.6    ; label("Log rate constant for pathogen-driven Th17 differentiation (log 1/week/pathogen)")  # mu.csv row 26
    lk8  <- 4.83   ; label("Log rate constant for pathogen-driven Th22 differentiation (log 1/week/pathogen)")  # mu.csv row 27

    # === Cytokines (dc6..dc13/dt) -- Eq. from AD_QSP_model.py L167-174 ===
    # IL-4 (from Th2)
    lk11 <- 5.5    ; label("Log Th2-driven IL-4 production rate (log 1/week/Th2)")                    # mu.csv row 28
    lk12 <- 9.25   ; label("Log constitutive IL-4 production rate (log fold-change/week)")            # mu.csv row 29
    ld10 <- 5.908  ; label("Log IL-4 elimination rate (log 1/week; elim-rate parameter)")             # mu.csv row 30 (Table S5)

    # IL-13 (from Th2)
    lk13 <- 6.65   ; label("Log Th2-driven IL-13 production rate (log 1/week/Th2)")                   # mu.csv row 31
    lk14 <- 8.84   ; label("Log constitutive IL-13 production rate (log fold-change/week)")           # mu.csv row 32
    ld11 <- 5.908  ; label("Log IL-13 elimination rate (log 1/week; elim-rate parameter)")            # mu.csv row 33 (Table S5)

    # IL-17A (from Th17)
    lk15 <- 4.01   ; label("Log Th17-driven IL-17A production rate (log 1/week/Th17)")                # mu.csv row 34
    lk16 <- 2.64   ; label("Log constitutive IL-17A production rate (log fold-change/week)")          # mu.csv row 35
    ld12 <- 3.332  ; label("Log IL-17A elimination rate (log 1/week; elim-rate parameter)")           # mu.csv row 36 (Table S5)

    # IL-22 (from Th22)
    lk17 <- 1.19   ; label("Log Th22-driven IL-22 production rate (log 1/week/Th22)")                 # mu.csv row 37
    lk18 <- 1.46   ; label("Log constitutive IL-22 production rate (log fold-change/week)")           # mu.csv row 38
    ld13 <- 3.332  ; label("Log IL-22 elimination rate (log 1/week; elim-rate parameter)")            # mu.csv row 39 (Table S5)

    # IL-31 (from Th2)
    lk19 <- 1.54   ; label("Log Th2-driven IL-31 production rate (log 1/week/Th2)")                   # mu.csv row 40
    lk20 <- 1.99   ; label("Log constitutive IL-31 production rate (log fold-change/week)")           # mu.csv row 41
    ld14 <- 3.332  ; label("Log IL-31 elimination rate (log 1/week; elim-rate parameter)")            # mu.csv row 42 (Table S5)

    # IFN-gamma (from Th1)
    lk21 <- 0.459  ; label("Log Th1-driven IFN-gamma production rate (log 1/week/Th1)")               # mu.csv row 43
    lk22 <- 2.62   ; label("Log constitutive IFN-gamma production rate (log fold-change/week)")       # mu.csv row 44
    ld15 <- 2.681  ; label("Log IFN-gamma elimination rate (log 1/week; elim-rate parameter)")        # mu.csv row 45 (Table S5)

    # TSLP (from keratinocytes; driven by infiltrated pathogens)
    lk23 <- 4      ; label("Log pathogen-driven TSLP production rate (log 1/week/pathogen)")          # mu.csv row 46
    lk24 <- 4.43   ; label("Log constitutive TSLP production rate (log fold-change/week)")            # mu.csv row 47
    ld16 <- 3.332  ; label("Log TSLP elimination rate (log 1/week; elim-rate parameter)")             # mu.csv row 48 (Table S5)

    # OX40L (driven by effective TSLP)
    lk25 <- 0.421  ; label("Log TSLP-driven OX40L production rate (log 1/week/TSLP)")                 # mu.csv row 49
    lk26 <- 1.68   ; label("Log constitutive OX40L production rate (log fold-change/week)")           # mu.csv row 50
    ld17 <- 0.4824 ; label("Log OX40L elimination rate (log 1/week; elim-rate parameter)")            # mu.csv row 51 (Table S5)

    # ============================================================
    # Inter-patient variability (log-normal). One eta per parameter;
    # all diagonal (no reported off-diagonal covariance in the source).
    # variance = sigma^2, where sigma is the log-scale SD from sigma.csv.
    # ============================================================
    # Skin barrier terms
    etalk1 ~ 0.628849  ; # sigma.csv row 1  (0.793)
    etalk2 ~ 0.152100  ; # sigma.csv row 2  (0.39)
    etalk3 ~ 2.016400  ; # sigma.csv row 3  (1.42)
    etalb1 ~ 0.330625  ; # sigma.csv row 4  (0.575)
    etalb2 ~ 2.496400  ; # sigma.csv row 5  (1.58)
    etalb3 ~ 0.121801  ; # sigma.csv row 6  (0.349)
    etalb4 ~ 0.636804  ; # sigma.csv row 7  (0.798)
    etalb5 ~ 0.815409  ; # sigma.csv row 8  (0.903)
    etald1 ~ 2.340900  ; # sigma.csv row 9  (1.53)
    etald2 ~ 0.321489  ; # sigma.csv row 10 (0.567)
    etald3 ~ 5.198400  ; # sigma.csv row 11 (2.28)
    # Pathogen terms
    etalb6 ~ 0.197136  ; # sigma.csv row 12 (0.444)
    etald4 ~ 0.224676  ; # sigma.csv row 13 (0.474)
    etald5 ~ 0.198025  ; # sigma.csv row 14 (0.445)
    etald6 ~ 0.151321  ; # sigma.csv row 15 (0.389)
    etald7 ~ 0.384400  ; # sigma.csv row 16 (0.62)
    etalb7 ~ 0.008464  ; # sigma.csv row 17 (0.092)
    etalb8 ~ 0.002304  ; # sigma.csv row 18 (0.048)
    etald8 ~ 0.001274  ; # sigma.csv row 19 (0.0357)
    # Th cells
    etalk5 ~ 0.001436  ; # sigma.csv row 20 (0.0379)
    etalk9 ~ 0.988036  ; # sigma.csv row 21 (0.994)
    etald9 ~ 0.007797  ; # sigma.csv row 22 (0.0883)
    etalb9 ~ 0.090000  ; # sigma.csv row 23 (0.3)
    etalk6 ~ 0.000009  ; # sigma.csv row 24 (0.00306)
    etalk10 ~ 0.004556 ; # sigma.csv row 25 (0.0675)
    etalk7 ~ 0.004858  ; # sigma.csv row 26 (0.0697)
    etalk8 ~ 0.225625  ; # sigma.csv row 27 (0.475)
    # Cytokines: IL-4
    etalk11 ~ 0.026569 ; # sigma.csv row 28 (0.163)
    etalk12 ~ 0.256036 ; # sigma.csv row 29 (0.506)
    etald10 ~ 0.047961 ; # sigma.csv row 30 (0.219)
    # IL-13
    etalk13 ~ 0.136900 ; # sigma.csv row 31 (0.37)
    etalk14 ~ 0.206116 ; # sigma.csv row 32 (0.454)
    etald11 ~ 0.012321 ; # sigma.csv row 33 (0.111)
    # IL-17
    etalk15 ~ 0.071824 ; # sigma.csv row 34 (0.268)
    etalk16 ~ 0.722500 ; # sigma.csv row 35 (0.85)
    etald12 ~ 0.326041 ; # sigma.csv row 36 (0.571)
    # IL-22
    etalk17 ~ 0.042436 ; # sigma.csv row 37 (0.206)
    etalk18 ~ 0.000767 ; # sigma.csv row 38 (0.0277)
    etald13 ~ 0.314721 ; # sigma.csv row 39 (0.561)
    # IL-31
    etalk19 ~ 0.034969 ; # sigma.csv row 40 (0.187)
    etalk20 ~ 0.331776 ; # sigma.csv row 41 (0.576)
    etald14 ~ 0.010000 ; # sigma.csv row 42 (0.1)
    # IFN-gamma
    etalk21 ~ 1.000000 ; # sigma.csv row 43 (1)
    etalk22 ~ 0.124609 ; # sigma.csv row 44 (0.353)
    etald15 ~ 0.021316 ; # sigma.csv row 45 (0.146)
    # TSLP
    etalk23 ~ 0.272484 ; # sigma.csv row 46 (0.522)
    etalk24 ~ 0.421201 ; # sigma.csv row 47 (0.649)
    etald16 ~ 0.041616 ; # sigma.csv row 48 (0.204)
    # OX40L
    etalk25 ~ 0.231361 ; # sigma.csv row 49 (0.481)
    etalk26 ~ 0.608400 ; # sigma.csv row 50 (0.78)
    etald17 ~ 0.055696 ; # sigma.csv row 51 (0.236)

    # ============================================================
    # Drug-scenario overrides (defaults reproduce steady-state / no drug).
    # Set at rxSolve() time via params = c(r_il13 = 0.99, ...) to simulate
    # a specific drug arm. All are dimensionless.
    # Source: AD_QSP_model.py simulate_one() drug_effect vector L201-237.
    # ============================================================
    r_placebo <- fixed(0)   ; label("Placebo / trial-participation flag (0=steady-state simulation; 1=any trial arm activates k3)")   # AD_QSP_model.py L91 (k3 = min(x[2], de[0]))
    r_il4     <- fixed(0)   ; label("IL-4 inhibition fraction by drug (0=no inhibition; 0.99=full antibody inhibition)")             # paper Sec 2.3, r_inhibit
    r_il13    <- fixed(0)   ; label("IL-13 inhibition fraction by drug (0=no inhibition; 0.99=full antibody inhibition)")            # paper Sec 2.3, r_inhibit
    r_il17    <- fixed(0)   ; label("IL-17A inhibition fraction by drug (0=no inhibition; 0.99=full antibody inhibition)")           # paper Sec 2.3, r_inhibit
    r_il22    <- fixed(0)   ; label("IL-22 inhibition fraction by drug (0=no inhibition; 0.99=full antibody inhibition)")            # paper Sec 2.3, r_inhibit
    r_il31    <- fixed(0)   ; label("IL-31 inhibition fraction by drug (0=no inhibition; 0.99=full antibody inhibition)")            # paper Sec 2.3, r_inhibit
    r_tslp    <- fixed(0)   ; label("TSLP inhibition fraction by drug (0=no inhibition; 0.99=full antibody inhibition)")             # paper Sec 2.3, r_inhibit
    r_ox40l   <- fixed(0)   ; label("OX40L inhibition fraction by drug (0=no inhibition; 0.99=full antibody inhibition)")            # paper Sec 2.3, r_inhibit
    c_rifng   <- fixed(0)   ; label("Additive rIFNg concentration in skin (0=no rIFNg; 210=rIFNg 50 ug/m^2 qd SC arm)")               # paper Sec 2.3 Eq, c_rIFNg
    e_a2      <- fixed(1)   ; label("Tralokinumab effective-inhibition factor multiplier (1.0 for lebrikizumab/dupilumab; 0.4396 for tralokinumab)")   # paper Sec 4.2 (e_a2 = 0.44); AD_QSP_model.py L140
  })

  model({
    # ==== Individual-patient parameter draws (exp(l<name> + eta<lname>)) ====
    k1 <- exp(lk1 + etalk1)
    k2 <- exp(lk2 + etalk2)
    k3_raw <- exp(lk3 + etalk3)
    b1 <- exp(lb1 + etalb1)
    b2 <- exp(lb2 + etalb2)
    b3 <- exp(lb3 + etalb3)
    b4 <- exp(lb4 + etalb4)
    b5 <- exp(lb5 + etalb5)
    d1 <- exp(ld1 + etald1)
    d2 <- exp(ld2 + etald2)
    d3 <- exp(ld3 + etald3)
    b6 <- exp(lb6 + etalb6)
    d4 <- exp(ld4 + etald4)
    d5 <- exp(ld5 + etald5)
    d6 <- exp(ld6 + etald6)
    d7 <- exp(ld7 + etald7)
    b7 <- exp(lb7 + etalb7)
    b8 <- exp(lb8 + etalb8)
    d8 <- exp(ld8 + etald8)
    k5 <- exp(lk5 + etalk5)
    k9 <- exp(lk9 + etalk9)
    d9 <- exp(ld9 + etald9)
    b9 <- exp(lb9 + etalb9)
    k6 <- exp(lk6 + etalk6)
    k10 <- exp(lk10 + etalk10)
    k7 <- exp(lk7 + etalk7)
    k8 <- exp(lk8 + etalk8)
    k11 <- exp(lk11 + etalk11)
    k12 <- exp(lk12 + etalk12)
    d10 <- exp(ld10 + etald10)
    k13 <- exp(lk13 + etalk13)
    k14 <- exp(lk14 + etalk14)
    d11 <- exp(ld11 + etald11)
    k15 <- exp(lk15 + etalk15)
    k16 <- exp(lk16 + etalk16)
    d12 <- exp(ld12 + etald12)
    k17 <- exp(lk17 + etalk17)
    k18 <- exp(lk18 + etalk18)
    d13 <- exp(ld13 + etald13)
    k19 <- exp(lk19 + etalk19)
    k20 <- exp(lk20 + etalk20)
    d14 <- exp(ld14 + etald14)
    k21 <- exp(lk21 + etalk21)
    k22 <- exp(lk22 + etalk22)
    d15 <- exp(ld15 + etald15)
    k23 <- exp(lk23 + etalk23)
    k24 <- exp(lk24 + etalk24)
    d16 <- exp(ld16 + etald16)
    k25 <- exp(lk25 + etalk25)
    k26 <- exp(lk26 + etalk26)
    d17 <- exp(ld17 + etald17)

    # k4 = d8 by design (pathogen production rate is set equal to the
    # pathogen self-loss constant). AD_QSP_model.py L141: `k4 = d8`.
    k4 <- d8

    # Placebo / trial-participation activation of the skin-barrier k3 term.
    # AD_QSP_model.py L91: k3 = min(x[2], de[0]); with de[0] in {0, 1e20}
    # this collapses to k3 = 0 (steady-state simulation) or k3 = x[2]
    # (any trial arm -- placebo, dupilumab, ...). Encoded here as a
    # multiplicative 0/1 switch driven by r_placebo.
    k3 <- k3_raw * r_placebo

    # === Drug-modified effective cytokine concentrations ===
    # Antibody drugs: c_effective = c_actual * (1 - r_inhibit).
    # rIFNg drug:     c_effective = c_actual + c_rifng (additive skin conc.).
    # Tralokinumab:   r_il13 is applied as r_il13 * e_a2 (e_a2 = 0.4396,
    # other IL-13 antibodies use e_a2 = 1). AD_QSP_model.py L140,144-151.
    IL4_eff  <- (1 - r_il4)             * il4
    IL13_eff <- (1 - r_il13 * e_a2)     * il13
    IL17_eff <- (1 - r_il17)            * il17
    IL22_eff <- (1 - r_il22)            * il22
    IL31_eff <- (1 - r_il31)            * il31
    TSLP_eff <- (1 - r_tslp)            * tslp
    OX40_eff <- (1 - r_ox40l)           * ox40l
    IFNg_eff <- ifng + c_rifng

    # =========================================================
    # 14 ODEs -- AD_QSP_model.py L155-174
    # =========================================================

    # Skin barrier integrity (0 <= barrier <= 1 at steady state):
    # dbarrier/dt = (1 - barrier) * (k1 + k2*IL22_eff + k3)
    #             / [(1 + b1*IL4)(1 + b2*IL13)(1 + b3*IL17)(1 + b4*IL22)(1 + b5*IL31)]
    #             - barrier * (d1*(1 + d3*pathogens) + d2*IL31)
    d/dt(barrier) <-
      (1 - barrier) * (k1 + k2 * IL22_eff + k3) /
        ((1 + b1 * IL4_eff) * (1 + b2 * IL13_eff) * (1 + b3 * IL17_eff) *
         (1 + b4 * IL22_eff) * (1 + b5 * IL31_eff)) -
      barrier * (d1 * (1 + d3 * pathogens) + d2 * IL31_eff)

    # Infiltrated pathogens:
    # dpathogens/dt = k4 / (1 + b6*barrier)
    #               - pathogens * [((1+d4*pathogens)(1+d5*IL17)(1+d6*IL22)(1+d7*IFNg))
    #                              / ((1+b7*IL4)(1+b8*IL13))
    #                              + d8]
    d/dt(pathogens) <-
      k4 / (1 + b6 * barrier) -
      pathogens *
        (((1 + d4 * pathogens) * (1 + d5 * IL17_eff) * (1 + d6 * IL22_eff) *
          (1 + d7 * IFNg_eff)) /
         ((1 + b7 * IL4_eff) * (1 + b8 * IL13_eff)) + d8)

    # Th cells: shared denominator (4 + k9*IFNg + k10*IL4) captures competitive
    # T-cell differentiation from a common naive-cell pool; d9 is shared
    # elimination modulated by OX40L (b9).
    d/dt(th1)  <- k5 * pathogens * (1 + k9 * IFNg_eff) /
                  (4 + k9 * IFNg_eff + k10 * IL4_eff) -
                  d9 * th1 / (1 + b9 * OX40_eff)
    d/dt(th2)  <- k6 * pathogens * (1 + k10 * IL4_eff) /
                  (4 + k9 * IFNg_eff + k10 * IL4_eff) -
                  d9 * th2 / (1 + b9 * OX40_eff)
    d/dt(th17) <- k7 * pathogens /
                  (4 + k9 * IFNg_eff + k10 * IL4_eff) -
                  d9 * th17 / (1 + b9 * OX40_eff)
    d/dt(th22) <- k8 * pathogens /
                  (4 + k9 * IFNg_eff + k10 * IL4_eff) -
                  d9 * th22 / (1 + b9 * OX40_eff)

    # Cytokines -- first-order production from Th sources + constitutive
    # production - first-order elimination. Note that dc/dt uses the raw
    # ODE-state cytokine level (not the drug-modified effective level);
    # drug inhibition acts algebraically at the "cytokine sees receptor"
    # step above, not on cytokine synthesis or clearance.
    d/dt(il4)   <- k11 * th2 + k12 - d10 * il4
    d/dt(il13)  <- k13 * th2 + k14 - d11 * il13
    d/dt(il17)  <- k15 * th17 + k16 - d12 * il17
    d/dt(il22)  <- k17 * th22 + k18 - d13 * il22
    d/dt(il31)  <- k19 * th2 + k20 - d14 * il31
    d/dt(ifng)  <- k21 * th1 + k22 - d15 * ifng
    d/dt(tslp)  <- k23 * pathogens + k24 - d16 * tslp
    d/dt(ox40l) <- k25 * TSLP_eff + k26 - d17 * ox40l

    # === Initial conditions (paper seed values; AD_QSP_model.py L183-198). ===
    # These are the Table-S2-informed baseline biological-factor levels used
    # to seed the 1000-week steady-state run. They are NOT the model's fitted
    # steady state -- users simulating a drug arm should first equilibrate for
    # >= ~1000 weeks and then apply drug. See vignette Steady-state check.
    barrier(0)   <- 0.5931
    pathogens(0) <- 0.4069
    th1(0)   <- 3.1
    th2(0)   <- 8.7
    th17(0)  <- 2.0
    th22(0)  <- 21.0
    il4(0)   <- 38.0
    il13(0)  <- 40.5
    il17(0)  <- 5.4
    il22(0)  <- 3.0
    il31(0)  <- 2.0
    ifng(0)  <- 1.5
    tslp(0)  <- 4.4
    ox40l(0) <- 9.7

    # === EASI algebraic observation ===
    # EASI = 72 * (2*pathogens + 2*(1 - barrier)) / 4
    #      = 36 * (pathogens + 1 - barrier)
    # AD_QSP_model.py EASI() L61-65. Ranges 0-72 (72 = max EASI).
    EASI <- 36 * (pathogens + 1 - barrier)
  })
}
