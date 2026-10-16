Ferl_2018_mosunetuzumab_cyno <- function() {
  description <- "Preclinical (cynomolgus monkey). Two-compartment population PK model of the anti-CD20/CD3 T-cell-dependent bispecific antibody mosunetuzumab (BTCT4465A) with first-order subcutaneous absorption and three elimination terms: a constant linear clearance, a rapidly decaying linear clearance (cl_exp_component * exp(-cl_exp_kdes * t), t1/2 1.6 h) and a slowly decaying Michaelis-Menten elimination (vmax * exp(-vmax_exp_kdes * t), t1/2 4.8 days); base Model VI without covariates (Ferl 2018)."
  reference <- "Ferl GZ, Reyes A, Sun LL, Cheu M, Oldendorp A, Ramanujan S, Stefanich EG. A Preclinical Population Pharmacokinetic Model for Anti-CD20/CD3 T-Cell-Dependent Bispecific Antibodies. Clin Transl Sci. 2018;11(3):296-304. doi:10.1111/cts.12535"
  vignette <- "Ferl_2018_mosunetuzumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "mosunetuzumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "mosunetuzumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "mosunetuzumab", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "cynomolgus monkey (Chinese origin)",
    n_subjects = NA,
    n_studies = 7,
    weight_median = "3.12 kg",
    disease_state = "Healthy, purpose-bred, treatment-naive cynomolgus monkeys (Charles River Laboratories, Reno, NV).",
    dose_range = "Single IV slow bolus 0.001-1 mg/kg (3-4 animals per group); repeat dosing once weekly x 4 at 0.01-1 mg/kg as IV slow bolus or 1-h IV infusion, or 1 mg/kg SC (3-8 animals per group).",
    notes = "Seven NHP safety / PK studies pooled (Ferl 2018 Methods 'PK/PD studies in nonhuman primates'); two of the seven were previously reported in Sun 2015 (Sci Transl Med 7:287ra70). Two cohorts received CHO-derived material, the rest E. coli-derived material; the single CHO study had n = 44. The total number of animals is not reported. ADA-positive samples were excluded. The molecule is called only BTCT4465A in the paper, which predates the INN; BTCT4465A is the Genentech development code for mosunetuzumab (RG7828). A naive-pooled murine fit (human CD20 / human CD3 double-transgenic mice, Figure 2a) is not included because its V1, V2 and CLd are not reported."
  )

  ini({
    # Ferl 2018 Table 2 (Model VI, no covariates) reports CL1, CL2_0, CLd, V1 and V2 per kg,
    # normalised by the 3.12 kg median body weight 'in a post hoc step' (Table 2 footnote), so the
    # NONMEM fit itself was in absolute units and carried no body-weight covariate. The values below
    # undo that normalisation for the typical 3.12 kg animal and convert mL to L (dose mg, Cc mg/L =
    # ug/mL). Vmax_0 is not in the post-hoc list and is already absolute.
    lka <- log(1.33); label("SC absorption rate constant Ka (1/day)") # Table 2 Ka = 1.33 1/day
    logitfdepot <- 1.66; label("Logit of SC bioavailability, F = 1/(1 + exp(-Z)) (unitless)") # Table 2 Z = 1.66 (F = 0.84)
    lcl_exp_inf <- log(7.05 * 3.12 / 1000); label("Constant nonsaturable clearance CL1 (L/day)") # Table 2 CL1 = 7.05 mL/day/kg x 3.12 kg (Figure S6: 22 mL/day)
    lcl_exp_component <- log(63.8 * 3.12 / 1000); label("Initial value of the time-varying linear clearance CL2_0 (L/day)") # Table 2 CL2_0 = 63.8 mL/day/kg x 3.12 kg
    lcl_exp_kdes <- log(10.1); label("Decay constant lambda2 of the time-varying linear clearance (1/day)") # Table 2 lambda2 = 10.1 1/day
    lq <- log(31.3 * 3.12 / 1000); label("Distribution clearance CLd (L/day)") # Table 2 CLd = 31.3 mL/day/kg x 3.12 kg
    lvc <- log(44.4 * 3.12 / 1000); label("Central volume V1 (L)") # Table 2 V1 = 44.4 mL/kg x 3.12 kg
    lvp <- log(47.3 * 3.12 / 1000); label("Peripheral volume V2 (L)") # Table 2 V2 = 47.3 mL/kg x 3.12 kg
    lvmax <- log(1280 / 1000); label("Initial Michaelis-Menten maximum elimination rate Vmax_0 (mg/day)") # Table 2 Vmax_0 = 1280 (units printed ug/mL; ug/day per Results text and Table S3)
    lvmax_exp_kdes <- log(0.144); label("Decay constant lambda1 of the time-varying Vmax (1/day)") # Table 2 lambda1 = 0.144 1/day
    lkm <- log(19.6); label("Michaelis-Menten constant KM (ug/mL)") # Table 2 KM = 19.6 ug/mL

    # IIV: omega^2 = log(1 + CV^2) from the Table 2 apparent %CVs (Eq. 5). Table S2, which is
    # labelled as the Model VI variance-covariance matrix, reproduces the Table S3 covariate-model
    # %CVs instead (see Ferl_2018_mosunetuzumab_cyno_covariates), so the base model carries the
    # Table 2 variances with no covariances.
    etalcl_exp_component ~ 0.5143 # Table 2 CL2_0 %CV = 82
    etalvmax ~ 0.1554 # Table 2 Vmax_0 %CV = 41
    etalcl_exp_inf ~ 0.0606 # Table 2 CL1 %CV = 25
    etalq ~ 0.3894 # Table 2 CLd %CV = 69
    etalvc ~ 0.1554 # Table 2 V1 %CV = 41
    etalvp ~ 0.0918 # Table 2 V2 %CV = 31

    # Supplement section B.9: NONMEM 7.3.0 with a proportional error model; its magnitude is not
    # reported anywhere in the paper or supplement.
    propSd <- fixed(0); label("Proportional residual error (fraction; not reported)")
  })
  model({
    ka <- exp(lka)
    fdepot <- expit(logitfdepot)
    cl_exp_inf <- exp(lcl_exp_inf + etalcl_exp_inf)
    cl_exp_component <- exp(lcl_exp_component + etalcl_exp_component)
    cl_exp_kdes <- exp(lcl_exp_kdes)
    q <- exp(lq + etalq)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    vmax <- exp(lvmax + etalvmax)
    vmax_exp_kdes <- exp(lvmax_exp_kdes)
    km <- exp(lkm)

    # Eq. (1) and Figure 1: CL(t) = CL1 + CL2_0 * exp(-lambda2 * t) and
    # Vmax(t) = Vmax_0 * exp(-lambda1 * t), with t the time since the first dose.
    cl <- cl_exp_inf + cl_exp_component * exp(-cl_exp_kdes * time)
    vmax_t <- vmax * exp(-vmax_exp_kdes * time)

    Cc <- central / vc

    # Eqs. (1)-(3)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (cl + q) / vc * central + q / vp * peripheral1 -
      vmax_t * Cc / (Cc + km)
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1
    f(depot) <- fdepot

    Cc ~ prop(propSd)
  })
}
