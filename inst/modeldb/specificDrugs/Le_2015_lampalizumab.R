Le_2015_lampalizumab <- function() {
  description <- "Combined ocular-serum target-mediated drug-disposition (TMDD) model with quasi-steady-state binding approximation for intravitreally administered lampalizumab (anti-complement factor D Fab) and total complement factor D (CFD) in adults with geographic atrophy secondary to age-related macular degeneration. Vitreous humor is the dosing compartment (depot) and the site of drug-target binding; aqueous humor lampalizumab and aqueous humor total CFD observations are derived from vitreous via constant partition coefficients; serum lampalizumab is the central elimination compartment with linear first-order clearance. Age and female sex modify ocular and systemic elimination rates respectively (Le 2015 Table 1, Eq. 1-7)."
  reference <- "Le KN, Gibiansky L, van Lookeren Campagne M, Good J, Davancaze T, Loyet KM, Morimoto A, Strauss EC, Jin JY. Population Pharmacokinetics and Pharmacodynamics of Lampalizumab Administered Intravitreally to Patients With Geographic Atrophy. CPT Pharmacometrics Syst Pharmacol. 2015;4(10):595-604. doi:10.1002/psp4.12031. PMID: 26535160."
  vignette <- "Le_2015_lampalizumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L (equivalent to ug/mL) for lampalizumab and total CFD")

  covariateData <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Power-form effect on the ocular elimination rate (kout) and the systemic elimination rate (k), normalised by reference age 80 years per Le 2015 Table 1 footnote 'For a typical 80-year-old male patient'. Exponents -0.770 (kout) and -1.63 (k); older patients have slower elimination.",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Female sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Time-fixed per subject. Multiplicative effect on the systemic elimination rate k; females have k multiplied by 0.739 relative to males (Le 2015 Table 1 row 'Multiplier for sex effect on k'). Reference category is male.",
      source_name        = "SEX (derived to SEXF: SEXF = as.integer(SEX == 'F'))"
    )
  )

  population <- list(
    species           = "human",
    n_subjects        = 117L,
    n_studies         = 2L,
    age_range         = "Reference 80 years (Le 2015 Table 1 footnote; per-subject age statistics live in Supplemental Table 1 which was not bundled with the main PDF)",
    age_median        = "Typical 80-year-old (Le 2015 Table 1 reference patient)",
    weight_range      = NULL,
    weight_median     = NULL,
    sex_female_pct    = NULL,
    race_ethnicity    = NULL,
    disease_state     = "Adults with geographic atrophy (GA) secondary to age-related macular degeneration (AMD).",
    dose_range        = "Single ITV doses 0.1, 0.5, 1, 2, 5, 10 mg per eye (phase Ia, CFD4711g; n=18). Multiple ITV doses 10 mg per eye every 4 weeks or every 8 weeks for up to 18 months (phase Ib/II MAHALO, CFD4870g; n=99).",
    regions           = NULL,
    n_observations    = "697 serum lampalizumab + 24 aqueous humor lampalizumab + 62 aqueous humor total CFD concentrations; ocular PK/PD subset n=21.",
    trial_identifiers = "CFD4711g (NCT00973011); CFD4870g MAHALO (NCT01229215).",
    notes             = "Detailed baseline demographics (weight, race, exact age distribution) reside in Supplemental Table 1 of the source paper, which was not on disk with the main PDF. The systemic PK dataset spans all 117 subjects; the ocular subset has 21 subjects with aqueous humor samples."
  )

  ini({
    # ------------------------------------------------------------------
    # Ocular structural parameters (typical values for 80-yr-old male).
    # Le 2015 Table 1.
    # ------------------------------------------------------------------
    lkout  <- log(0.117);          label("Ocular elimination rate constant of unbound drug kout (1/day; 80-yr-old male reference)")    # Le 2015 Table 1, row 1
    lkoutC <- log(0.135);          label("Ocular elimination rate constant of drug-target complex koutC (1/day)")                       # Le 2015 Table 1, row 2
    lkinT  <- log(0.364);          label("Ocular target influx/synthesis rate kinT (mg/L/day, equivalent to ug/mL/day)")                # Le 2015 Table 1, row 3 (paper Discussion p.600 also cites 0.32 ug/mL/day; table point estimate 0.364 used)
    lkA    <- log(2.23);           label("Correction factor kA for drug:target molar weight and assay differences (unitless)")          # Le 2015 Table 1, row 4
    lVVITR <- log(0.00309);        label("Vitreous humor volume of distribution VVITR (L; equivalent to 3.09 mL)")                       # Le 2015 Table 1, row 5 (3.09 mL)
    lkoutT <- log(0.27);           label("Ocular degradation/elimination rate constant of target koutT (1/day)")                         # Le 2015 Table 1, row 6
    lkAQ   <- log(13.0);           label("Vitreous-to-aqueous partition coefficient of drug kAQ (unitless)")                             # Le 2015 Table 1, row 7

    # ------------------------------------------------------------------
    # Systemic structural parameters.
    # ------------------------------------------------------------------
    lVc    <- log(2.41);           label("Serum (central) volume of distribution Vc (L; equivalent to 2410 mL)")                          # Le 2015 Table 1, row 8 (2410 mL)
    lk     <- log(1.89);           label("Systemic first-order elimination rate constant k (1/day; 80-yr-old male reference)")            # Le 2015 Table 1, row 9

    # ------------------------------------------------------------------
    # Quasi-steady-state binding constant. Fixed in the source paper at
    # 20 pM (~ in-vitro KD of 19.7 pM), equivalent to ~0.96e-3 mg/L for
    # a 48 kDa Fab. Le 2015 Table 1 row 'Kss' footnote d: "Fixed to 20 pM,
    # approximately equal to in vitro KD value (19.7 pM), as the estimated
    # elimination rate of the drug-target complex (koutC) is much smaller
    # than dissociation rate of the drug-target complex." log() goes
    # inside fixed() per naming-conventions.md.
    # ------------------------------------------------------------------
    lKss   <- fixed(log(0.96e-3)); label("Quasi-steady-state binding constant Kss (mg/L = ug/mL); fixed at ~20 pM ~ in-vitro KD")        # Le 2015 Table 1, row 10 (FIXED)

    # ------------------------------------------------------------------
    # Covariate effects on elimination rates (Le 2015 Table 1, covariates
    # block). Power-form on AGE/80 for kout and k; multiplicative on k
    # for female sex (reference = male).
    # ------------------------------------------------------------------
    e_age_kout <- -0.770;          label("Power exponent for (AGE/80) on log-kout (unitless)")                                            # Le 2015 Table 1, row 'Power for age effect on kout'
    e_age_k    <- -1.63;           label("Power exponent for (AGE/80) on log-k (unitless)")                                              # Le 2015 Table 1, row 'Power for age effect on k'
    e_sexf_k   <-  0.739;          label("Multiplicative effect on k for female sex (unitless; reference = male)")                       # Le 2015 Table 1, row 'Multiplier for sex effect on k'

    # ------------------------------------------------------------------
    # Inter-subject variability (log-normal). Paper reports CV% for the
    # three retained random effects. Mapping CV% -> variance of the log
    # eta: omega^2 = CV^2 (treating CV as the SD of log eta on the small-
    # variance approximation that the source paper uses; the equivalent
    # log(1 + CV^2) form is within ~10% for these CV magnitudes). The
    # other parameters had IIV fixed to 1% CV during backward elimination
    # (Le 2015 Methods, p.596) and are not carried as etas here.
    # ------------------------------------------------------------------
    etalkout ~ 0.0745   # Le 2015 Table 1: CV 27.3% on kout -> 0.273^2 = 0.0745
    etalkA   ~ 0.349    # Le 2015 Table 1: CV 59.1% on kA -> 0.591^2 = 0.349
    etalk    ~ 0.0756   # Le 2015 Table 1: CV 27.5% on k -> 0.275^2 = 0.0756

    # ------------------------------------------------------------------
    # Residual error. Paper reports a SINGLE proportional residual for
    # aqueous-humor measurements (25.8% CV) shared between aqueous
    # lampalizumab and aqueous CFD, plus a separate proportional residual
    # for serum lampalizumab (32.9% CV). nlmixr2's residual-error syntax
    # requires one endpoint parameter per output, so the shared aqueous
    # residual is encoded here as two parameters (CAQpropSd, RAQpropSd)
    # initialised to the same Table 1 value. Re-estimation would treat
    # them as independent; documented in the vignette's Assumptions and
    # deviations.
    # ------------------------------------------------------------------
    propSd     <- 0.329;           label("Proportional residual SD for serum lampalizumab Cc (fraction)")                                  # Le 2015 Table 1, row 'Residual error for serum measurement' (32.9%)
    propSd_CAQ <- 0.258;           label("Proportional residual SD for aqueous humor lampalizumab (fraction; shared with RAQ in source)")  # Le 2015 Table 1, row 'Residual error for aqueous measurements' (25.8%; one residual for both aqueous outputs)
    propSd_RAQ <- 0.258;           label("Proportional residual SD for aqueous humor total CFD (fraction; shared with CAQ in source)")     # Le 2015 Table 1, row 'Residual error for aqueous measurements' (25.8%; one residual for both aqueous outputs)
  })

  model({
    # ---- Individual parameters (Le 2015 Methods + Table 1 covariate forms) ----
    kout  <- exp(lkout + etalkout) * (AGE / 80)^e_age_kout
    koutC <- exp(lkoutC)
    kinT  <- exp(lkinT)
    kA    <- exp(lkA + etalkA)
    VVITR <- exp(lVVITR)
    koutT <- exp(lkoutT)
    kAQ   <- exp(lkAQ)
    Vc    <- exp(lVc)
    k     <- exp(lk + etalk) * (AGE / 80)^e_age_k * (e_sexf_k^SEXF)
    Kss   <- exp(lKss)

    # ---- Derived: target vitreous-to-aqueous partition coefficient (Eq. 7) ----
    kTAQ <- kAQ * kA

    # ---- Quasi-steady-state algebra (Le 2015 Eq. 4): Cunbound from total drug and total target ----
    # depot holds the total (free + bound to CFD) lampalizumab amount in vitreous (mg).
    # total_target holds the total (free + bound to drug) CFD CONCENTRATION in vitreous (mg/L).
    CVITR    <- depot / VVITR
    qss_disc <- (CVITR - total_target - Kss)^2 + 4 * Kss * CVITR
    Cunbound <- 0.5 * ((CVITR - total_target - Kss) + sqrt(qss_disc))

    # ---- Pre-dose steady-state baseline of target (set kinT - koutT * Rss = 0) ----
    total_target(0) <- kinT / koutT

    # ---- ODE system (Le 2015 Eq. 1-3) ----
    # Eq. 1: drug in vitreous (depot) loses unbound drug at rate kout * Cunbound * VVITR
    #        and complex at rate koutC * (R * Cunbound / (Kss + Cunbound)) * VVITR.
    d/dt(depot)        <- -kout  * Cunbound * VVITR -
                           koutC * (total_target * Cunbound * VVITR) / (Kss + Cunbound)

    # Eq. 2: total target in vitreous: zero-order influx kinT minus linear degradation of free
    #        target (koutT) minus net additional loss of bound target via complex elimination.
    d/dt(total_target) <-  kinT -
                           koutT * total_target -
                           (koutC - koutT) * (total_target * Cunbound) / (Kss + Cunbound)

    # Eq. 3: serum drug gains the egressing unbound drug and complex from vitreous and is
    #        eliminated linearly at rate k.
    d/dt(central)      <-  kout  * Cunbound * VVITR +
                           koutC * (total_target * Cunbound * VVITR) / (Kss + Cunbound) -
                           k * central

    # ---- Observations (Le 2015 Eq. 5 and 6 plus serum concentration) ----
    Cc  <- central / Vc                # Serum lampalizumab (mg/L = ug/mL)
    CAQ <- CVITR / kAQ                 # Aqueous humor total lampalizumab (Le 2015 Eq. 5)
    RAQ <- total_target / kTAQ         # Aqueous humor total CFD (Le 2015 Eq. 6)

    Cc  ~ prop(propSd)
    CAQ ~ prop(propSd_CAQ)
    RAQ ~ prop(propSd_RAQ)
  })
}
