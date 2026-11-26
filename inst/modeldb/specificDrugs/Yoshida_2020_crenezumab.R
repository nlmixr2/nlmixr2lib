Yoshida_2020_crenezumab <- function() {
  description <- "Population PK/PD model for crenezumab (humanized anti-amyloid-beta IgG4) in adults with mild-to-moderate Alzheimer's disease: two-compartment crenezumab PK with first-order SC absorption, joined by a Michaelis-Menten-approximated target-mediated drug disposition (TMDD) model to two endogenous targets, plasma monomeric amyloid-beta(1-40) and amyloid-beta(1-42), each tracked as a total (free + crenezumab-bound) pool with its own peripheral distribution compartment."
  reference <- "Yoshida K, Moein A, Bittner T, Ostrowitzki S, Lin H, Honigberg L, Jin JY, Quartino A. Pharmacokinetics and pharmacodynamic effect of crenezumab on plasma and cerebrospinal fluid beta-amyloid in patients with mild-to-moderate Alzheimer's disease. Alzheimers Res Ther. 2020;12:16. doi:10.1186/s13195-020-0580-2"
  vignette <- "Yoshida_2020_crenezumab"
  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "ug/mL",
    concentration_Ctotal_ab40 = "pM",
    concentration_Ctotal_ab42 = "pM"
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (WT / 72.3)^e on crenezumab CLel, CLint and Vcent. 72.3 kg is the weight of the typical patient in Figure 7 (population median used by the PsN stepwise covariate search); 90% of the observed range was 52.4-95.9 kg.",
      source_name = "BWT"
    ),
    AGE = list(
      description = "Age at enrollment",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear effect (1 + e * (AGE - 72)) on the amyloid-beta(1-40) baseline only. 72 years is the typical patient of Figure 7; 90% of the observed range was 54-82 years.",
      source_name = "Age"
    ),
    CRCL = list(
      description = "Glomerular filtration rate, BSA-normalized",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear effect (1 + e * (CRCL - 72.5)) on both amyloid-beta baselines. 72.5 mL/min/1.73 m^2 is the typical patient of Figure 7; 90% of the observed range was 46.5-94.0. The paper does not name the estimating equation behind its GFR.",
      source_name = "GFR"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (female)",
      notes = "The paper estimates a male deviation on Kdeg ('Sex (male) on Kdeg' = -0.16), shared by amyloid-beta(1-40) and (1-42). Encoded as (1 + e_sex_kdeg * (1 - SEXF)) so that the printed coefficient is used unchanged and females are the reference.",
      source_name = "Sex"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "crenezumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "crenezumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "crenezumab", units = "mg", specimen = "tissue", verified = TRUE),
    total_ab40 = list(
      analyte = "amyloid-beta(1-40), total (free + crenezumab-bound)",
      units = "pmol",
      specimen = "plasma",
      verified = TRUE
    ),
    total_ab40_peripheral = list(
      analyte = "amyloid-beta(1-40), total",
      units = "pmol",
      specimen = "tissue",
      verified = TRUE
    ),
    total_ab42 = list(
      analyte = "amyloid-beta(1-42), total (free + crenezumab-bound)",
      units = "pmol",
      specimen = "plasma",
      verified = TRUE
    ),
    total_ab42_peripheral = list(
      analyte = "amyloid-beta(1-42), total",
      units = "pmol",
      specimen = "tissue",
      verified = TRUE
    )
  )

  # Total (free + crenezumab-bound) amyloid-beta pools of the Michaelis-Menten
  # TMDD approximation and their peripheral distribution compartments. The
  # ab40 / ab42 analyte tokens follow Ramakrishnan_2023_amyloid_qsp.R.
  paper_specific_compartments <- c(
    "total_ab40",
    "total_ab40_peripheral",
    "total_ab42",
    "total_ab42_peripheral"
  )

  population <- list(
    species = "human",
    n_subjects = 131L,
    n_studies = 3L,
    age_range = "50-90 years (study inclusion criteria); 90% of observed range 54-82 years",
    age_median = "72 years (typical patient of Figure 7)",
    weight_range = "90% of observed range 52.4-95.9 kg",
    weight_median = "72.3 kg (typical patient of Figure 7)",
    renal_function = "GFR 90% of observed range 46.5-94.0 mL/min/1.73 m^2; typical 72.5",
    disease_state = "Mild-to-moderate Alzheimer's disease",
    dose_range = "Crenezumab 300 mg SC every 2 weeks, or 15, 30, 45, 60 or 120 mg/kg IV every 4 weeks",
    regions = "Multinational (phase II ABBY and BLAZE, phase Ib GN29632)",
    notes = "Pooled crenezumab-treated patients from ABBY (NCT01343966; amyloid-beta data from the safety run-in cohort only), BLAZE (NCT01397578) and GN29632 (NCT02353598); 1332 serum PK and 2203 plasma PD samples (Yoshida 2020 Methods, PK/PD analysis). The paper prints no baseline-demographics table; ranges are those quoted for the covariate analysis (Results; Figure 7). Sex distribution not reported."
  )

  ini({
    # Crenezumab PK -- typical values for a 72.3 kg patient (Yoshida 2020 Table 2)
    lcl <- log(0.159); label("Elimination clearance of free crenezumab, CLel (L/day)") # Table 2 'CLel (L/day)' 0.159
    lclint <- log(1.01); label("Intrinsic clearance of the crenezumab-amyloid-beta complex, CLint (L/day)") # Table 2 'CLint (L/day)' 1.01
    lvc <- log(2.89); label("Central volume of distribution, Vcent (L)") # Table 2 'Vcent (L)' 2.89
    lvp <- log(1.84); label("Peripheral volume of distribution of crenezumab, Vperiph (L)") # Table 2 crenezumab 'Vperiph (L)' 1.84
    lq <- log(0.142); label("Intercompartmental clearance of crenezumab, Q (L/day)") # Table 2 crenezumab 'Q (L/day)' 0.142
    lfdepot <- log(0.662); label("Subcutaneous bioavailability, Fsc (fraction)") # Table 2 'Fsc (%)' 66.2
    lka <- log(0.161); label("First-order SC absorption rate constant, Ka (1/day)") # Table 2 'Ka (/day)' 0.161

    # Binding (Michaelis-Menten / quasi-equilibrium constant)
    lkd_ab40 <- log(12.0); label("Crenezumab-amyloid-beta(1-40) equilibrium dissociation constant, Kd (nM)") # Table 2 'Kd Abeta(1-40) (nM)' 12.0
    lkd_ab42 <- log(9.37); label("Crenezumab-amyloid-beta(1-42) equilibrium dissociation constant, Kd (nM)") # Table 2 'Kd Abeta(1-42) (nM)' 9.37

    # Endogenous amyloid-beta turnover and distribution
    lrbase_ab40 <- log(142); label("Baseline plasma amyloid-beta(1-40) for a 72-year-old with GFR 72.5 (pM)") # Table 2 Abeta(1-40) 'Baseline (pM)' 142
    lkdeg_ab40 <- log(114); label("Degradation rate constant of free amyloid-beta(1-40), female (1/day)") # Table 2 Abeta(1-40) 'Kdeg (/day)' 114
    lq_ab40 <- log(8.72); label("Intercompartmental clearance of amyloid-beta(1-40) (L/day)") # Table 2 Abeta(1-40) 'Q (L/day)' 8.72
    lvp_ab40 <- log(3.77); label("Peripheral distribution volume of amyloid-beta(1-40) (L)") # Table 2 Abeta(1-40) 'Vperiph (L)' 3.77
    lrbase_ab42 <- log(5.98); label("Baseline plasma amyloid-beta(1-42) for GFR 72.5 (pM)") # Table 2 Abeta(1-42) 'Baseline (pM)' 5.98
    lkdeg_ab42 <- log(287); label("Degradation rate constant of free amyloid-beta(1-42), female (1/day)") # Table 2 Abeta(1-42) 'Kdeg (/day)' 287
    lq_ab42 <- log(229); label("Intercompartmental clearance of amyloid-beta(1-42) (L/day)") # Table 2 Abeta(1-42) 'Q (L/day)' 229
    lvp_ab42 <- log(8.19); label("Peripheral distribution volume of amyloid-beta(1-42) (L)") # Table 2 Abeta(1-42) 'Vperiph (L)' 8.19

    # Covariate effects: power form on crenezumab parameters, linear form on amyloid-beta parameters (Table 2 footnote)
    e_wt_cl <- 0.835; label("Power exponent of body weight on CLel (unitless)") # Table 2 'BWT (kg) on CLel' 0.835
    e_wt_clint <- 0.474; label("Power exponent of body weight on CLint (unitless)") # Table 2 'BWT (kg) on CLint' 0.474
    e_wt_vc <- 0.725; label("Power exponent of body weight on Vcent (unitless)") # Table 2 'BWT (kg) on Vcent' 0.725
    e_age_rbase_ab40 <- 0.0037; label("Linear age effect on amyloid-beta(1-40) baseline (1/year)") # Table 2 Abeta(1-40) 'Age (year) on baseline' 0.0037
    e_crcl_rbase_ab40 <- -0.0035; label("Linear GFR effect on amyloid-beta(1-40) baseline (per mL/min/1.73 m^2)") # Table 2 Abeta(1-40) 'GFR on baseline' -0.0035
    e_crcl_rbase_ab42 <- -0.0048; label("Linear GFR effect on amyloid-beta(1-42) baseline (per mL/min/1.73 m^2)") # Table 2 Abeta(1-42) 'GFR on baseline' -0.0048
    e_sex_kdeg <- -0.16; label("Male fractional change in Kdeg, shared by amyloid-beta(1-40) and (1-42) (unitless)") # Table 2 'Sex (male) on Kdeg' -0.16, footnote b

    # Unit conversion for the binding terms (Kd in nM, crenezumab in mg/L)
    mw <- fixed(150000); label("Crenezumab molecular weight used for mass-to-molar conversion (g/mol)") # not printed in the paper; nominal IgG molecular weight

    # IIV (Table 2 'BSV (CV%)', converted with omega^2 = log(1 + CV^2))
    etalcl ~ 0.03583 # Table 2 CLel BSV 19.1%
    etalclint ~ 0.02314 # Table 2 CLint BSV 15.3%
    etalvc ~ 0.03294 # Table 2 Vcent BSV 18.3%
    etalvp ~ 0.4482 # Table 2 crenezumab Vperiph BSV 75.2%
    etalkd ~ 0.01833 # Table 2 Kd BSV 13.6%, one eta shared by both Kd values (footnote a)
    etalrbase_ab40 ~ 0.006866 # Table 2 Abeta(1-40) baseline BSV 8.3%
    etalrbase_ab42 ~ 0.01886 # Table 2 Abeta(1-42) baseline BSV 13.8%
    etalkdeg ~ 0.03153 # Table 2 Kdeg BSV 17.9%, one eta shared by both Kdeg values (footnote a)

    # Residual error: not reported in the paper
    propSd <- fixed(0); label("Proportional residual error on serum crenezumab (fraction; not reported in source)") # not reported in the source
    propSd_Ctotal_ab40 <- fixed(0); label("Proportional residual error on total plasma amyloid-beta(1-40) (fraction; not reported in source)") # not reported in the source
    propSd_Ctotal_ab42 <- fixed(0); label("Proportional residual error on total plasma amyloid-beta(1-42) (fraction; not reported in source)") # not reported in the source
  })

  model({
    # Crenezumab individual parameters (power covariate form, reference 72.3 kg)
    cl <- exp(lcl + etalcl) * (WT / 72.3)^e_wt_cl
    clint <- exp(lclint + etalclint) * (WT / 72.3)^e_wt_clint
    vc <- exp(lvc + etalvc) * (WT / 72.3)^e_wt_vc
    vp <- exp(lvp + etalvp)
    q <- exp(lq)
    ka <- exp(lka)
    fdepot <- exp(lfdepot)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    # Complex elimination rate constant: CLint acts on the bound amyloid-beta
    # carried in the central (serum/plasma) volume shared with crenezumab.
    kint <- clint / vc

    # Amyloid-beta individual parameters (linear covariate form)
    kd_ab40 <- exp(lkd_ab40 + etalkd)
    kd_ab42 <- exp(lkd_ab42 + etalkd)
    rbase_ab40 <- exp(lrbase_ab40 + etalrbase_ab40) *
      (1 + e_age_rbase_ab40 * (AGE - 72)) *
      (1 + e_crcl_rbase_ab40 * (CRCL - 72.5))
    rbase_ab42 <- exp(lrbase_ab42 + etalrbase_ab42) *
      (1 + e_crcl_rbase_ab42 * (CRCL - 72.5))
    sex_kdeg <- 1 + e_sex_kdeg * (1 - SEXF)
    kdeg_ab40 <- exp(lkdeg_ab40 + etalkdeg) * sex_kdeg
    kdeg_ab42 <- exp(lkdeg_ab42 + etalkdeg) * sex_kdeg
    q_ab40 <- exp(lq_ab40)
    vp_ab40 <- exp(lvp_ab40)
    q_ab42 <- exp(lq_ab42)
    vp_ab42 <- exp(lvp_ab42)

    # Zero-order synthesis that holds each total pool at its baseline before
    # dosing (no loss from the peripheral compartments)
    ksyn_ab40 <- kdeg_ab40 * rbase_ab40 * vc
    ksyn_ab42 <- kdeg_ab42 * rbase_ab42 * vc

    # Crenezumab concentration (mg/L = ug/mL) and its molar equivalent (nM)
    Cc <- central / vc
    Cc_nM <- Cc * 1e6 / mw

    # Michaelis-Menten approximation (Figure 1): the free fraction of each
    # amyloid-beta species is Kd / (Kd + Ccren) and the bound fraction
    # Ccren / (Kd + Ccren)
    ffree_ab40 <- kd_ab40 / (kd_ab40 + Cc_nM)
    ffree_ab42 <- kd_ab42 / (kd_ab42 + Cc_nM)

    # Target-mediated elimination of the complex (pmol/day); the matching loss
    # of crenezumab mass is (pmol/day) * mw * 1e-9 mg/pmol
    tmde_ab40 <- kint * total_ab40 * (1 - ffree_ab40)
    tmde_ab42 <- kint * total_ab42 * (1 - ffree_ab42)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1 -
      (tmde_ab40 + tmde_ab42) * mw * 1e-9
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(total_ab40) <- ksyn_ab40 - kdeg_ab40 * total_ab40 * ffree_ab40 - tmde_ab40 -
      q_ab40 / vc * total_ab40 + q_ab40 / vp_ab40 * total_ab40_peripheral
    d/dt(total_ab40_peripheral) <- q_ab40 / vc * total_ab40 - q_ab40 / vp_ab40 * total_ab40_peripheral
    d/dt(total_ab42) <- ksyn_ab42 - kdeg_ab42 * total_ab42 * ffree_ab42 - tmde_ab42 -
      q_ab42 / vc * total_ab42 + q_ab42 / vp_ab42 * total_ab42_peripheral
    d/dt(total_ab42_peripheral) <- q_ab42 / vc * total_ab42 - q_ab42 / vp_ab42 * total_ab42_peripheral

    total_ab40(0) <- rbase_ab40 * vc
    total_ab40_peripheral(0) <- rbase_ab40 * vp_ab40
    total_ab42(0) <- rbase_ab42 * vc
    total_ab42_peripheral(0) <- rbase_ab42 * vp_ab42

    f(depot) <- fdepot

    # Total (measured) and free (model-inferred) plasma amyloid-beta, pM
    Ctotal_ab40 <- total_ab40 / vc
    Ctotal_ab42 <- total_ab42 / vc
    Cfree_ab40 <- Ctotal_ab40 * ffree_ab40
    Cfree_ab42 <- Ctotal_ab42 * ffree_ab42

    Cc ~ prop(propSd)
    Ctotal_ab40 ~ prop(propSd_Ctotal_ab40)
    Ctotal_ab42 ~ prop(propSd_Ctotal_ab42)
  })
}
