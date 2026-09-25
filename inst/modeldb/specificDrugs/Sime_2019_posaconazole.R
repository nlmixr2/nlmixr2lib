Sime_2019_posaconazole <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for posaconazole in critically ill adults,",
    "fitted SIMULTANEOUSLY to total and unbound plasma concentrations with explicit",
    "capacity-limited (Michaelis-Menten) albumin binding. The central compartment carries two",
    "states: a free posaconazole pool and an albumin-bound pool exchanging by second-order",
    "association (kon) and first-order dissociation (koff), giving the equilibrium isotherm",
    "Cbound = Bmax * Cfree / (KD + Cfree) with KD = koff / kon. Elimination and intercompartmental",
    "distribution both act on the UNBOUND concentration. Serum albumin sets the binding capacity",
    "Bmax = ALB * N * (MW_posaconazole / MW_albumin) * 1000 with N fixed to 1 binding site per",
    "albumin molecule, and body mass index scales the central volume linearly (V = Vtheta * BMI/24).",
    "Fitted NON-PARAMETRICALLY with the NPAG algorithm in Pmetrics; the Table 2 means are encoded",
    "as lognormal medians and the tabulated CV percentages as independent lognormal marginal",
    "variances. Residual unexplained variability is carried as fixed(0) because the selected",
    "Pmetrics error model was never published.",
    sep = " "
  )
  reference <- paste(
    "Sime FB, Byrne CJ, Parker S, Stuart J, Butler J, Starr T, Pandey S, Wallis SC, Lipman J,",
    "Roberts JA. Population pharmacokinetics of total and unbound concentrations of intravenous",
    "posaconazole in adult critically ill patients. Crit Care. 2019;23(1):205.",
    "doi:10.1186/s13054-019-2483-9. PMCID PMC6554926.",
    "Structural model from Fig. 1 and Methods Eqs. 2-5; parameter estimates from Table 2;",
    "covariate model from Results 'Pharmacokinetic model building'.",
    sep = " "
  )
  vignette <- "Sime_2019_posaconazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central = list(
      analyte = "posaconazole (unbound)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE,
      notes = paste(
        "Holds the FREE (protein-unbound) posaconazole amount in the central compartment,",
        "Cf(t) in Sime 2019 Fig. 1. Dividing by vc gives the unbound plasma concentration",
        "directly; the measured total concentration is the sum of this state and the bound",
        "state. The intravenous dose is administered into this state.",
        sep = " "
      )
    ),
    complex = list(
      analyte = "posaconazole bound to serum albumin",
      units = "mg",
      specimen = "plasma",
      verified = TRUE,
      notes = paste(
        "Holds the ALBUMIN-BOUND posaconazole amount in the central compartment, Cb(t) in",
        "Sime 2019 Fig. 1 (the dashed inner box). Reuses the registered TMDD `complex` state:",
        "albumin is the binding partner, so this is a reversibly bound drug-protein species in",
        "mass-action exchange with the free pool, exactly the role the register assigns to",
        "`complex`. Free albumin is not carried as a `target` state because the model tracks",
        "only the aggregate binding capacity Bmax, whose unoccupied part is the algebraic",
        "difference Bmax - Cbound.",
        sep = " "
      )
    ),
    peripheral1 = list(
      analyte = "posaconazole",
      units = "mg",
      specimen = "plasma",
      verified = TRUE,
      notes = paste(
        "Cp(t) in Sime 2019 Fig. 1. Only unbound posaconazole distributes in and out (Fig. 1",
        "legend: kcp and kpc are 'rate constant[s] for distribution of unbound posaconazole'),",
        "so no binding is modelled in this compartment.",
        sep = " "
      )
    )
  )

  covariateData <- list(
    ALB = list(
      description = "Serum albumin concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Sets the maximum albumin binding capacity in Sime 2019 Eq. 3,",
        "Bmax = Alb * N * (Mposa / MAlb) * 1000, with N fixed to 1 and the factor 1000",
        "converting g/L to mg/L. Albumin therefore enters as a binding-capacity driver, NOT",
        "as a covariate effect with a reference value and an e_<cov>_<param> coefficient --",
        "the same kind of use the register records for AAG in Said 2025 imatinib.",
        "Cohort median 20 g/L (IQR 18-24), Sime 2019 Table 1; that median reproduces the",
        "paper's reported mean unbound fraction of 0.65 percent to three significant figures.",
        "Results 'Pharmacokinetic model building': 'The only covariates that improved the",
        "goodness of fit and significantly reduced the objective function were BMI for volume",
        "of distribution (V) and albumin for Bmax.'",
        sep = " "
      ),
      source_name = "Albumin"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Scales the central volume linearly and through the origin,",
        "V = Vtheta * BMI / 24 (Sime 2019 Results 'Pharmacokinetic model building':",
        "'BMI was best related to V linearly and normalised to 24 (i.e. V = V x BMI/24, where",
        "V is typical value of V and 24 is the median BMI of study patients)').",
        "NOTE a discrepancy in the source: Table 1 reports the cohort median BMI as",
        "22.6 kg/m^2 (IQR 20.2-29.7), not 24. The normalising constant 24 is used here",
        "because it is the value written into the covariate equation; see vignette Errata.",
        "The dosing simulations of Tables 3-8 span BMI 17-38 kg/m^2.",
        sep = " "
      ),
      source_name = "BMI"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on V and CL but not retained (Sime 2019 Methods 'Development of the covariate model'). Cohort median 46 years (IQR 40-51), Table 1."
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Screened as 'gender' on V and CL but not retained (Sime 2019 Methods). Cohort 7 male / 1 female, Table 1."
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      notes = "Screened on V and CL but not retained (Sime 2019 Methods). Not tabulated separately; enters the retained BMI covariate."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened on V and CL but not retained (Sime 2019 Methods); BMI was retained on V instead. Cohort median 68 kg (IQR 65-82), Table 1."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Screened on V and CL but not retained (Sime 2019 Methods). Cohort median 106 umol/L (IQR 78-197), Table 1."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened on V and CL but not retained (Sime 2019 Methods). Cohort median 53 (IQR 28-60), Table 1, where the unit is printed as IU/mL; IU/L is intended."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened on V and CL but not retained (Sime 2019 Methods). Cohort median 47 (IQR 38-130), Table 1, where the unit is printed as IU/mL; IU/L is intended."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units = "U/L",
      type = "continuous",
      notes = "Screened on V and CL but not retained (Sime 2019 Methods). Cohort median 75 (IQR 63-108), Table 1, where the unit is printed as IU/mL; IU/L is intended."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened on V and CL but not retained (Sime 2019 Methods). Cohort median 11 umol/L (IQR 10-20), Table 1."
    ),
    GGT = list(
      description = "Gamma-glutamyltransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened on V and CL but not retained (Sime 2019 Methods). Not tabulated in Table 1."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 8L,
    n_studies = 1L,
    n_observations = 93L,
    age_median = "46 years (IQR 40-51)",
    weight_median = "68 kg (IQR 65-82)",
    bmi_median = "22.6 kg/m^2 (IQR 20.2-29.7)",
    sex_female_pct = 12.5,
    race_ethnicity = "not reported in the source paper",
    disease_state = paste(
      "Critically ill adults admitted to a quaternary referral intensive care unit with presumed",
      "or confirmed invasive fungal infection requiring systemic antifungal therapy. Median",
      "APACHE II score at ICU admission 17 (IQR 17-24); median SOFA score 5 (IQR 3-6) on day 1",
      "and 3 (IQR 2-4) on day 2. Four of eight patients (50 percent) had a positive fungal",
      "culture (Candida albicans, C. dubliniensis, C. parapsilosis, C. glabrata complex,",
      "Candida spp.). All patients were receiving another antifungal as usual care",
      "(fluconazole 75 percent, voriconazole 25 percent, caspofungin 37 percent,",
      "amphotericin 12 percent). Marked hypoalbuminaemia: median albumin 20 g/L (IQR 18-24).",
      sep = " "
    ),
    renal_function = "Serum creatinine median 106 umol/L (IQR 78-197); measured urinary creatinine clearance median 74 mL/min (IQR 53-109)",
    hepatic_function = "ALT median 53, AST median 47, ALP median 75 (units printed as IU/mL in Table 1); total bilirubin median 11 umol/L (IQR 10-20)",
    dose_range = "Single 300 mg intravenous posaconazole given as a 90-minute infusion through a central venous catheter",
    regions = "Australia (Royal Brisbane and Women's Hospital intensive care unit, Brisbane)",
    notes = paste(
      "Prospective observational population PK study; baseline demographics from Sime 2019",
      "Table 1. Each patient received one 300 mg intravenous posaconazole dose as an add-on to",
      "their prescribed antifungal therapy, with 14 blood samples over 48 h (pre-dose; 15, 45,",
      "75 and 90 min during the infusion; then 3, 5, 8, 12, 18, 24, 30, 36 and 48 h after the",
      "start of infusion). Total and unbound concentrations were measured by UHPLC-MS/MS",
      "(total calibration range 0.02-5 mg/L; unbound 0.0005-0.1 mg/L after ultracentrifugation),",
      "yielding 93 paired total / unbound observations. The observed unbound fraction was",
      "median 0.55 percent (IQR 0.36-1.9), mean 0.65 percent (SD 0.39), CV 58.5 percent.",
      "The full screened covariate set in Methods 'Development of the covariate model' was age,",
      "gender, height, weight, body mass index, albumin, serum creatinine, urinary creatinine",
      "clearance, presence of renal replacement therapy, alanine aminotransferase, aspartate",
      "aminotransferase, alkaline phosphatase, bilirubin, gamma-glutamyltransferase and SOFA",
      "score; only BMI (on V) and albumin (on Bmax) were retained. Urinary creatinine clearance,",
      "renal replacement therapy and the SOFA score are recorded here rather than in",
      "covariatesDataExcluded because no canonical covariate column matches them as the paper",
      "reports them (the creatinine clearance is a measured urinary value in plain mL/min, not",
      "BSA-normalised; the renal replacement modality is not stated).",
      "The paper notes the small sample size and narrow covariate spread as its main limitation.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # Structural parameters. All six values are the MEAN of the NPAG
    # non-parametric parameter distribution reported in Sime 2019 Table 2
    # ('Pharmacokinetic parameter estimates for the final covariate model').
    # NPAG estimates a discrete joint density over support points rather than a
    # point estimate plus a parametric variance; following the repository
    # convention for Pmetrics fits (see Hughes 2024 vancomycin, Setiawan 2023
    # sulbactam), the tabulated mean is encoded as the MEDIAN of a lognormal
    # marginal. See vignette Errata for the mean-versus-median consequence.
    #
    # ke, kcp and kpc are the micro-constants of Fig. 1. ke and kcp multiply the
    # UNBOUND central concentration, which is why their magnitudes are so large:
    # scaled by the unbound fraction of about 0.65 percent they give a clearance
    # of roughly 20 L/h, a steady-state volume of roughly 500 L and a terminal
    # half-life of roughly 19 h, all consistent with published posaconazole
    # disposition. See the vignette Source trace for the arithmetic.
    # ------------------------------------------------------------------------
    lkel <- log(42.07); label("Elimination rate constant acting on the unbound central concentration (1/h)")
    # Table 2, row 'Ke (h-1)' mean = 42.07 (SD 23.68, CV 56)
    lvc <- log(72.19); label("Central volume of distribution at the reference BMI of 24 kg/m2 (L)")
    # Table 2, row 'V theta (l)' mean = 72.19 (SD 43.14, CV 60); Abstract rounds to 72 (43) L
    lk12 <- log(334.27); label("Distribution rate constant central to peripheral1, acting on the unbound central concentration (1/h)")
    # Table 2, row 'Kcp (h-1)' mean = 334.27 (SD 236.42, CV 71)
    lk21 <- log(0.37); label("Distribution rate constant peripheral1 to central (1/h)")
    # Table 2, row 'Kpc (h-1)' mean = 0.37 (SD 0.11, CV 31)

    # ------------------------------------------------------------------------
    # Albumin-binding rate constants (Sime 2019 Eq. 4, KD = 1/KA = koff/kon).
    # Their ratio KD = 3897.40 / 2820.35 = 1.3819 mg/L is the only combination
    # the observable kinetics are sensitive to, because both constants are four
    # to five orders of magnitude faster than every PK rate constant in the
    # model: the free-bound exchange is at quasi-equilibrium on any observable
    # timescale. They are carried separately anyway because that is how the
    # authors parameterised and reported the model.
    # ------------------------------------------------------------------------
    lkon <- log(2820.35); label("Second-order association rate constant for posaconazole binding to albumin (L/mg/h)")
    # Table 2, row 'Kon (L/mg/h)' mean = 2820.35 (SD 671.99, CV 24)
    lkoff <- log(3897.40); label("First-order dissociation rate constant for posaconazole from albumin (1/h)")
    # Table 2, row 'Koff (h-1)' mean = 3897.40 (SD 596.46, CV 15)

    # ------------------------------------------------------------------------
    # Binding stoichiometry. Methods 'Structural base model and binding model',
    # final sentence: 'N was assumed to be 1.' Not estimated.
    # ------------------------------------------------------------------------
    nalb <- fixed(1); label("Number of posaconazole binding sites per albumin molecule (unitless)")
    # Methods, Eq. 3 definition list: 'N is the number of posaconazole binding sites per molecule of albumin. N was assumed to be 1.'

    # ------------------------------------------------------------------------
    # Inter-individual variability. NPAG estimates a discrete non-parametric
    # distribution rather than a parametric omega; Sime 2019 Table 2 reports a
    # mean, an SD and a CV percentage for each parameter, and the CV column is
    # SD/mean over the support points (23.68/42.07 = 56 percent and so on for
    # every row). The CV column is carried here as a LOG-NORMAL approximation
    # using omega^2 = log(CV^2 + 1), which reproduces each reported CV exactly.
    # Off-diagonal covariances are not reported, so the marginals are encoded as
    # independent -- an approximation imposed on a non-parametric joint density.
    # See vignette Assumptions and deviations.
    # ------------------------------------------------------------------------
    # Table 2 Ke CV = 56 -> log(0.56^2 + 1)
    etalkel ~ 0.27277146
    # Table 2 V theta CV = 60 -> log(0.60^2 + 1)
    etalvc ~ 0.30748470
    # Table 2 Kcp CV = 71 -> log(0.71^2 + 1)
    etalk12 ~ 0.40819471
    # Table 2 Kpc CV = 31 -> log(0.31^2 + 1)
    etalk21 ~ 0.09175843
    # Table 2 Kon CV = 24 -> log(0.24^2 + 1)
    etalkon ~ 0.05600219
    # Table 2 Koff CV = 15 -> log(0.15^2 + 1)
    etalkoff ~ 0.02225061

    # ------------------------------------------------------------------------
    # Residual unexplained variability is NOT reported. Methods 'Error model'
    # states only that a multiplicative (Error = SD*gamma) and an additive
    # (Error = [SD^2 + lambda^2]^0.5) error model were tested, and that assay
    # error was modelled as a linear polynomial Error = C0 + C1*[obs] 'starting
    # with a generic set of coefficients, followed by iterative optimization'.
    # Neither the selected model nor any of gamma, lambda, C0 or C1 appears in
    # the paper, and there is no supplement on disk carrying them. Carried as
    # fixed(0) on both outputs rather than invented -- see vignette Errata.
    # ------------------------------------------------------------------------
    propSd <- fixed(0); label("Proportional residual SD on the total concentration (fraction; 0 because not reported in the source)")
    propSd_Cunbound <- fixed(0); label("Proportional residual SD on the unbound concentration (fraction; 0 because not reported in the source)")
  })

  model({
    # 1. Physical constants used by the binding-capacity equation (Sime 2019
    #    Eq. 3). The paper defines Mposa and MAlb symbolically without printing
    #    values. Posaconazole (C37H42F2N8O4) has a molecular weight of
    #    700.8 g/mol and human serum albumin of 66500 g/mol -- the same albumin
    #    molecular weight already used by Fauchet_2015_lopinavir_unbound.R and
    #    recorded in the ALB covariate register entry. That pair is confirmed
    #    arithmetically by the paper's own numbers: at the cohort median albumin
    #    of 20 g/L they give Bmax = 210.77 mg/L, and with
    #    KD = koff/kon = 1.3819 mg/L the low-concentration unbound fraction
    #    fu = KD / (Bmax + KD) = 0.651 percent reproduces the reported mean
    #    unbound fraction of 0.65 percent (Results, 'Plasma protein binding').
    MW_POSA <- 700.8
    MW_ALB <- 66500

    # 2. Individual parameters. BMI scales the central volume linearly through
    #    the origin, normalised to 24 kg/m2 (Results: V = V x BMI/24).
    vc <- exp(lvc + etalvc) * (BMI / 24)
    kel <- exp(lkel + etalkel)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    kon <- exp(lkon + etalkon)
    koff <- exp(lkoff + etalkoff)

    # 3. Maximum albumin binding capacity, Sime 2019 Eq. 3:
    #      Bmax = Alb * N * (Mposa / MAlb) * 1000
    #    The 1000 converts the albumin concentration from g/L to mg/L so that
    #    Bmax is a posaconazole concentration in mg/L, matching KD.
    bmax <- ALB * nalb * (MW_POSA / MW_ALB) * 1000

    # 4. Plasma concentrations of the two central species.
    Cunbound <- central / vc
    Cbound <- complex / vc

    # 5. Mass-action albumin binding (Sime 2019 Fig. 1, kon / koff arrows). The
    #    unoccupied binding capacity is bmax - Cbound, so the net association
    #    flux in mg/h is
    #      kon * Cunbound * (bmax - Cbound) * vc - koff * complex
    #    whose stationary point is exactly the Eq. 2 isotherm
    #      Cbound = bmax * Cunbound / (KD + Cunbound),  KD = koff / kon.
    bindflux <- kon * Cunbound * (bmax - Cbound) * vc - koff * complex

    # 6. ODE system. Elimination (ke) and distribution to the peripheral
    #    compartment (kcp) act on the UNBOUND central concentration, i.e. on the
    #    free central amount, per the Fig. 1 legend; return from the peripheral
    #    compartment is first-order in the peripheral amount. The intravenous
    #    dose enters the free pool.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 - bindflux
    d/dt(complex) <- bindflux
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 7. Observations. Sime 2019 Eq. 5 gives Cfree = Ctotal - Cbound, i.e. the
    #    measured total concentration is the sum of the free and bound species.
    Cc <- Cunbound + Cbound

    Cc ~ prop(propSd)
    Cunbound ~ prop(propSd_Cunbound)
  })
}
