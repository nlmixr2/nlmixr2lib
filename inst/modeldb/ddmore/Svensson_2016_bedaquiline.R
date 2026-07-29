Svensson_2016_bedaquiline <- function() {
  description <- "Three-compartment population PK model for the antimycobacterial bedaquiline (BDQ) and its one-compartment N-desmethyl metabolite M2 in adult patients with multidrug-resistant tuberculosis (MDR-TB), with two-transit-compartment first-order oral absorption, time-varying body-weight allometric scaling, time-varying serum-albumin power effects on disposition and metabolite formation/elimination, and additional Black-race and linear age covariate effects on bedaquiline and M2 clearance."
  reference <- paste(
    "Svensson E. M., Dosne A.-G., Karlsson M. O. (2016).",
    "Population Pharmacokinetics of Bedaquiline and Metabolite M2 in",
    "Patients With Drug-Resistant Tuberculosis: The Effect of",
    "Time-Varying Weight and Albumin.",
    "CPT Pharmacometrics Syst Pharmacol 5(12):682-691.",
    "doi:10.1002/psp4.12147.",
    "DDMORE Foundation Model Repository: DDMODEL00000219.",
    sep = " "
  )
  vignette <- "Svensson_2016_bedaquiline"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")
  ddmore_id    <- "DDMODEL00000219"
  replicate_of <- NULL

  covariateData <- list(
    WT = list(
      description        = "Body weight (time-varying; supplied per observation)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying body weight applied as allometric power on clearances and volumes around a reference of 70 kg: `cl = cl_typ * (WT/70)^e_wt_cl` (estimated exponent 0.181) and `vc = vc_typ * (WT/70)^e_wt_vc` (exponent fixed to 1.0). The same `(WT/70)^e_wt_cl` factor is shared across CL, Q1, Q2, CLM2 and the same `(WT/70)^e_wt_vc` factor across Vc, Vp1, Vp2, Vc_m2 (Svensson 2016 Methods, equations after Table 3). The publication itself fits a semi-physiological linear weight trajectory (WT0 -> WT120 over 120 weeks) and uses *individually predicted* time-varying weight inside the PK model; the nlmixr2lib implementation expects the user to supply WT either as a constant baseline column or as the time-varying trajectory of choice (see vignette Errata).",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Serum albumin concentration (time-varying; supplied per observation)",
      units              = "g/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying serum albumin in g/dL applied as a power effect around the typical population steady-state value Ass = 4.04 g/dL: `cl = cl_typ * (ALB/4.04)^e_alb_cl` with estimated exponent 1.64 (Svensson 2016 Table 3, 'Individual time varying effect of albumin CL/CLM2/fbm'). The same exponent enters CLM2 twice — once directly and once via the inverse of the albumin power on the bedaquiline-to-M2 fraction-metabolised fbm — so the net albumin power on CLM2 is `(ALB/4.04)^(2*1.64)`. Volumes carry an additional fixed-coefficient unbound-fraction-driven adjustment `(4.04/ALB)^e_alb_vc` with `e_alb_vc` fixed to 1.0 (THETA(21) FIX in the .mod and 'Time varying effect of protein binding on disposition - 1 Fix' in Table 3). The publication itself fits a semi-physiological self-limiting logistic albumin model (A0 = 3.65 g/dL recovering to Ass = 4.04 g/dL with T1/2_return = 20.4 weeks per Table 2) and uses individually predicted time-varying albumin inside the PK model; the nlmixr2lib implementation expects the user to supply ALB either as a constant baseline column or as the time-varying trajectory of choice (see vignette Errata).",
      source_name        = "ALB"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear deviation around the median of 32 years applied to clearances: `cl = cl_typ * (1 + e_age_cl * (32 - AGE))` with estimated coefficient 0.00881 1/year (Svensson 2016 Table 3 footnote 'Age effect on CL/CLM2'). The reference age is the cohort median (32 years) reported in Table 1. Younger patients (AGE < 32) have higher CL; older patients have lower CL.",
      source_name        = "AGE"
    ),
    RACE_BLACK = list(
      description        = "Black-race indicator (1 = Black, 0 = non-Black)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Black; reference category includes White, Asian, Hispanic, and Other in the Svensson 2016 cohort).",
      notes              = "Multiplicative effect on clearances: `cl = cl_typ * (1 + e_race_black_cl * RACE_BLACK)` with estimated coefficient 0.84 (Svensson 2016 Table 3, 'Effect of black race on CL/CM2'). Black patients have ~84% higher BDQ and M2 CL relative to the non-Black reference. The .mod source column is `RACE` with the value 2 mapped to Black via `IF (RACE.EQ.2) BLACK = 1`; supply the canonical 0/1-coded `RACE_BLACK` column directly.",
      source_name        = "RACE"
    )
  )

  population <- list(
    n_subjects     = 335L,
    n_studies      = 2L,
    age_range      = "18-68 years (Table 1)",
    age_median     = "32 years",
    weight_range   = "Stage 1: 37-81 kg (median 55 kg) at start of treatment; stage 2: 30-113 kg (median 57 kg)",
    weight_median  = "55-57 kg (study-stage dependent; Table 1)",
    sex_female_pct = NA_real_,
    race_ethnicity = c(
      Black = 38.2,
      Asian = 38.2,
      Other = 23.6
    ),
    disease_state  = "Adult patients with pulmonary multidrug-resistant tuberculosis (MDR-TB) or pre-extensively / extensively drug-resistant tuberculosis (Pre-XDR / XDR-TB), enrolled in the C208 (placebo-controlled, two-stage) and C209 (open-label) phase II trials. Subjects with HIV co-infection were eligible if CD4+ counts were >250 cells/microL.",
    dose_range     = "Oral bedaquiline 400 mg once daily for the first 2 weeks (loading), followed by 200 mg three times weekly for 22 weeks, on top of an individualized second-line MDR-TB background regimen. Total treatment: 24 weeks of bedaquiline. Followed by 96 weeks of off-bedaquiline follow-up.",
    regions        = "Multicenter international (C208 and C209). Specific regional breakdown is not reported in the main text.",
    notes          = "Pooled analysis of two phase II trials (C208, two-stage placebo-controlled; C209, open-label) reported in Svensson 2016 Methods. Race breakdown is taken from Table 1 study-stage columns; the precise overall percentages depend on weighting C208 stage 1 (n=21 active), C208 stage 2 (n=79 active), and C209 (n=233 active) — Table 1 reports per-stage percentages and absolute counts. The race mix above represents the pooled C208 stage 2 + C209 mix and should be treated as approximate; consult Table 1 for cohort-specific values."
  )

  ini({
    # Final estimates are taken from the .mod $THETA / $OMEGA / $SIGMA initial-value
    # block of DDMODEL00000219's `Executable_BDQ_M2_PK_plus_WT_ALB_in_MDR-TB_patients.mod`
    # rather than from the bundle's `Output_real_*.lst`. The .lst is an
    # `ESTIMATION STEP OMITTED: YES` evaluation-only run with `R MATRIX
    # ALGORITHMICALLY SINGULAR` / `COVARIANCE STEP ABORTED`, and its INITIAL
    # ESTIMATE block disagrees with the bundled .mod $THETA for THETA(3),
    # THETA(12), THETA(14), THETA(16), THETA(18), and THETA(24) by factors of
    # 100-1000 (likely a different .mod with internal rescaling generated the
    # .lst). The .mod $THETA values used here match Svensson 2016 Table 3 (the
    # final fixed-effects estimates) and Table 2 (the albumin/weight
    # longitudinal model) exactly; see the validation vignette's Assumptions
    # and deviations section.

    # Structural absorption: two-transit-compartment first-order chain
    # parameterized via the mean absorption time MAT (hours) and the
    # fraction-of-MAT residing in the transit delay FR (paper symbols).
    # KA and KTR are derived in model() from MAT and FR.
    lmat     <- log(0.6620 * 6)              ; label("Mean absorption time (hours, on log scale)") # THETA(9) typical value 0.662 fraction of 6 h ≡ 3.97 h; Table 3 'MAT, fraction of 6 hours = 0.66'
    logitfr  <- log(0.4664 / (1 - 0.4664))   ; label("Fraction of MAT in transit delay (logit)")   # THETA(10) FR = 0.466; Table 3 'FR = 0.47'

    # Structural disposition: 3-compartment bedaquiline (parent) +
    # 1-compartment M2 metabolite. Apparent volumes and clearances (i.e.,
    # CL/F, V/F for the parent and CL/(F*fm), V/(F*fm) for the metabolite).
    lcl     <- log(2.616)   ; label("Apparent linear clearance of bedaquiline CL/F (L/h)")               # THETA(11); Table 3 'CL/F = 2.62 L/h'
    lvc     <- log(198.34)  ; label("Apparent central volume of bedaquiline Vc/F (L)")                   # THETA(12); Table 3 'V/F = 198 L/70 kg'
    lq      <- log(3.658)   ; label("Apparent first-peripheral inter-compartmental clearance Q1/F (L/h)")  # THETA(13); Table 3 'Q1/F = 3.66 L/h'
    lvp     <- log(8549.06) ; label("Apparent first-peripheral volume Vp1/F (L)")                          # THETA(14); Table 3 'VP1/F = 8550 L/70 kg'
    lq2     <- log(7.335)   ; label("Apparent second-peripheral inter-compartmental clearance Q2/F (L/h)") # THETA(15); Table 3 'Q2/F = 7.34 L/h'
    lvp2    <- log(2690.91) ; label("Apparent second-peripheral volume Vp2/F (L)")                         # THETA(16); Table 3 'VP2/F = 2690 L/70 kg'
    lcl_m2  <- log(10.0496) ; label("Apparent M2 metabolite clearance CLM2/(F*fm) (L/h)")                  # THETA(17); Table 3 'CLM2/(F*fm) = 10.0 L/h'
    lvc_m2  <- log(2203.71) ; label("Apparent M2 metabolite central volume Vc_m2/(F*fm) (L)")              # THETA(18); Table 3 'VM2/(F*fm) = 2200 L/70 kg'

    # Covariate effects.
    e_wt_cl        <- 0.1809                       ; label("Allometric body-weight exponent on clearances (CL, Q1, Q2, CLM2; reference 70 kg)")   # THETA(19); Table 3 'Allometric scaling clearances = 0.18'
    e_wt_vc        <- fixed(1.0)                   ; label("Allometric body-weight exponent on volumes (Vc, Vp1, Vp2, Vc_m2; reference 70 kg)")    # THETA(20) FIX; Table 3 footnote 'Allometric scaling of volumes - coefficients fixed to 1'
    e_alb_vc       <- fixed(1.0)                   ; label("Albumin-driven unbound-fraction power exponent applied to apparent volumes (reference Ass = 4.04 g/dL)")  # THETA(21) FIX; Table 3 'Time varying effect of protein binding on disposition - 1 Fix'
    e_alb_cl       <- 1.640                        ; label("Albumin power exponent on bedaquiline CL and M2 fbm (reference Ass = 4.04 g/dL)")   # THETA(22); Table 3 'Individual time varying effect of albumin CL/CLM2/fbm = 1.64'
    e_race_black_cl <- 0.8387                      ; label("Multiplicative Black-race effect on CL and CLM2 (vs non-Black reference)")           # THETA(23); Table 3 'Effect of black race on CL/CM2 = 0.84'
    e_age_cl       <- 0.008808                     ; label("Linear age effect on CL and CLM2 (1/year, around 32 years)")                          # THETA(24); Table 3 'Age effect on CL/CLM2 = 0.0088'

    # Bioavailability is fixed to F = 1 because CL and V are reported as
    # apparent F-relative values (CL/F, V/F). The .mod's `F1 = 1.8002 *
    # EXP(BOVF + BSVF)` carries an explicit unit-conversion factor (mg dose
    # to nmol/mL plasma concentration) — that factor is dropped here because
    # nlmixr2lib reports concentration in mg/L; only the BSV part is retained
    # below. See vignette Errata.
    lfdepot <- fixed(log(1))                       ; label("Bioavailability F (fixed at 1 because CL and V are apparent F-relative values)")

    # Inter-individual variability. Final estimates (variances) from the .mod
    # $OMEGA blocks. The full .mod also models BSV on baseline albumin / Ass /
    # T1/2_return / WT0 / WT120 (Box-Cox-transformed), BOV on F and MAT
    # between-occasion, and BSV on residual error (BSVRUV1 / BSVRUV2). The
    # nlmixr2lib implementation drops the albumin-and-weight semi-physiological
    # model (and therefore those BSVs / Box-Cox shapes) and drops the BOV and
    # residual-error etas; see the validation vignette's Assumptions and
    # deviations section.
    etalcl + etalcl_m2 ~ c(0.1528, 0.1349, 0.2121)    # .mod $OMEGA BLOCK(2) blocks 11-12; Table 3 BSV CL = 40.7 % CV, BSV CLM2 = 48.6 % CV, correlation 75.0 %
    etalvc      ~ 0.1719                              # .mod $OMEGA block 13; Table 3 BSV V = 43.3 % CV
    etalq       ~ 0.1812                              # .mod $OMEGA block 14; Table 3 BSV Q1 = 44.5 % CV
    etalvc_m2   ~ 0.1502                              # .mod $OMEGA block 15; Table 3 BSV VM2 = 40.2 % CV
    etalfdepot  ~ 0.0803                              # .mod $OMEGA block 10 'BSV F'; Table 3 BSV F = 28.9 % CV (BOV F = 19.7 % CV is dropped)
    etalmat     ~ 1.1620                              # .mod $OMEGA block 8 'BOV MAT' value (also approximates the absent dedicated BSV-MAT block); Table 3 BOV MAT = 148 % CV

    # Residual error.  The .mod codes residuals as additive on log-transformed
    # concentrations (`Y = log(C) + EPS` with `var(EPS) = sigma^2`), which is
    # equivalent to a proportional residual in linear space with
    # SD(C/C_pred) ~= sqrt(sigma^2). The .mod's correlated $SIGMA BLOCK(2)
    # between BDQ and M2 residual errors (correlation 43.3 % per Table 3) is
    # rendered as independent proportional residuals here because nlmixr2lib
    # has no idiomatic encoding for cross-output residual correlation; see the
    # vignette's Assumptions and deviations section.
    propSd     <- sqrt(0.05182)   ; label("Bedaquiline residual error (proportional CV)")   # .mod $SIGMA block 3 (var = 0.0518); Table 3 'Proportional residual error BDQ = 23.1 % CV'
    propSd_m2  <- sqrt(0.03668)   ; label("M2 metabolite residual error (proportional CV)") # .mod $SIGMA block 4 diagonal (var = 0.0367); Table 3 'Proportional residual error M2 = 19.3 % CV'
  })

  model({
    # Covariate effects (multiplicative, around their respective reference
    # categories). `RACE_BLACK`, `AGE`, `WT`, and `ALB` are user-supplied
    # columns from the dataset (canonical names per
    # `inst/references/covariate-columns.md`).
    blackcl <- 1 + e_race_black_cl * RACE_BLACK
    agecl   <- 1 + e_age_cl * (32 - AGE)
    allcl   <- (WT / 70)^e_wt_cl
    allv    <- (WT / 70)^e_wt_vc
    albcl   <- (ALB / 4.04)^e_alb_cl
    fmcorr  <- (ALB / 4.04)^e_alb_cl   # = 1 / (ALB/4.04)^(-e_alb_cl); the M2-side albumin correction via 1/fbm (paper's `1/FM` term)
    albvfu  <- (4.04 / ALB)^e_alb_vc   # albumin-driven free-fraction power on apparent volumes (THETA(21) = 1 FIX)

    # Absorption derived from MAT (hours) and FR (unitless fraction in (0, 1)).
    # MAT inherits its IIV via etalmat; FR has no IIV in the source model.
    mat <- exp(lmat + etalmat)
    fr  <- 1 / (1 + exp(-logitfr))
    mtt <- mat * fr
    ka  <- log(2) / (mat * (1 - fr) / 3.3)   # absorption rate (1/h); the 3.3-fold delay-to-absorption-half-life ratio is the paper's parameterization (Svensson 2016 Methods, transit description)
    ktr <- 2 / mtt                            # transit-chain rate (1/h) for two transit compartments

    # Bedaquiline (parent) PK parameters.
    cl  <- exp(lcl  + etalcl)    * blackcl * agecl * albcl * allcl
    vc  <- exp(lvc  + etalvc)    * allv
    q   <- exp(lq   + etalq)     * allcl
    vp  <- exp(lvp)              * allv
    q2  <- exp(lq2)              * allcl
    vp2 <- exp(lvp2)             * allv

    # M2 metabolite PK parameters. The metabolite clearance carries both the
    # albumin-on-CL effect (`albcl`) and the inverse albumin-on-fbm correction
    # (`fmcorr`), so the net albumin power on cl_m2 is `(ALB/4.04)^(2*e_alb_cl)`
    # — matching the .mod's `CLM2 = CLM2B/FM*COVALBI*ALLCL` after substituting
    # `FM = (ALB/Ass)^(-e_alb_cl)` and `COVALBI = (ALB/Ass)^e_alb_cl`. The
    # metabolite volume carries the single-power albumin-via-fbm correction.
    cl_m2 <- exp(lcl_m2 + etalcl_m2) * blackcl * agecl * albcl * fmcorr * allcl
    vc_m2 <- exp(lvc_m2 + etalvc_m2) * fmcorr * allv

    # Effective volumes for the observation equations only — i.e., the
    # albumin-driven free-fraction correction (`albvfu`) the .mod applies via
    # `VE = VB*ALLVE*COVFUIP` in $ERROR while leaving `V` unchanged in $PK.
    # The ODE clearance terms below use the un-corrected `vc`, matching the
    # source model's mass-balance bookkeeping; the observation equations
    # divide by these `*_obs` volumes.
    vc_obs    <- vc    * albvfu
    vc_m2_obs <- vc_m2 * albvfu

    # Two-transit-compartment first-order absorption chain feeding bedaquiline
    # central. depot and transit1 share the rate ktr; transit2 -> central uses
    # ka. Mirrors the .mod $DES (compartments DEPOT/TRANSI1/TRANSI2 -> BDQC).
    d/dt(depot)        <- -ktr * depot
    d/dt(transit1)     <-  ktr * depot - ktr * transit1
    d/dt(transit2)     <-  ktr * transit1 - ka * transit2
    d/dt(central)      <-  ka * transit2 - cl  * central / vc -
                           q  * central / vc + q  * peripheral1 / vp -
                           q2 * central / vc + q2 * peripheral2 / vp2
    d/dt(peripheral1)  <-  q  * central / vc - q  * peripheral1 / vp
    d/dt(peripheral2)  <-  q2 * central / vc - q2 * peripheral2 / vp2
    d/dt(central_m2)   <-  cl * central / vc - cl_m2 * central_m2 / vc_m2

    # Bioavailability (etalfdepot carries BSV on F; the unit-conversion factor
    # 1.8002 in the .mod's `F1 = 1.8002 * EXP(BOVF + BSVF)` is dropped because
    # nlmixr2lib reports concentration in mg/L rather than the .mod's
    # nmol/mL = umol/L; see vignette Errata).
    f(depot) <- exp(lfdepot + etalfdepot)

    # Observed plasma concentrations. The albumin-driven unbound-fraction
    # correction enters here only (matches the .mod's $ERROR block).
    Cc    <- central    / vc_obs
    Cc_m2 <- central_m2 / vc_m2_obs

    Cc    ~ prop(propSd)
    Cc_m2 ~ prop(propSd_m2)
  })
}
