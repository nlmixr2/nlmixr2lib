Ferl_2018_mosunetuzumab_human <- function() {
  description <- "Allometric human (70 kg) projection of the cynomolgus-monkey two-compartment PK model of the anti-CD20/CD3 T-cell-dependent bispecific antibody mosunetuzumab (BTCT4465A): constant linear clearance, a decaying linear clearance fixed to the rituximab NHL decay rate, and a decaying Michaelis-Menten elimination scaled from the NHP covariate model; typical-value IV model reproducing Figure 4 (Ferl 2018)."
  reference <- "Ferl GZ, Reyes A, Sun LL, Cheu M, Oldendorp A, Ramanujan S, Stefanich EG. A Preclinical Population Pharmacokinetic Model for Anti-CD20/CD3 T-Cell-Dependent Bispecific Antibodies. Clin Transl Sci. 2018;11(3):296-304. doi:10.1111/cts.12535"
  vignette <- "Ferl_2018_mosunetuzumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "mosunetuzumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "mosunetuzumab", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = NA,
    weight_median = "70 kg (nominal)",
    disease_state = "Prospective projection for patients with B-cell malignancies; no human data were fitted.",
    dose_range = "Figure 4 simulates 1, 10 and 100 mg flat IV doses every 28 days for five cycles.",
    notes = "Projected from the cynomolgus-monkey final model with covariates (Ferl_2018_mosunetuzumab_cyno_covariates; Ferl 2018 Table S3, E. coli-derived reference material at the median baseline CD4+ count, 3.12 kg median animal) by Eq. (8), p_h = p_c * (BW_h / BW_c)^b, with b = 0.85 on CL1 and CLd, b = 1 on V1, V2 and Vmax_0, KM unscaled, lambda1 fixed to the NHP value, and CL2_0 and lambda2 fixed to the rituximab NHL estimates of Ferl 2018 Table 3 (Li 2012) as an approximation of CD20+ tumour-mediated clearance (Methods 'Allometric scaling'; Results 'Model-predicted human plasma PK'). Typical-value projection only: no between-subject variability was projected, and SC absorption was not projected. The authors caution that the nonlinear and time-varying parameters should not be extrapolated to humans without clinical validation (Discussion). BTCT4465A is the Genentech development code for mosunetuzumab (RG7828)."
  )

  ini({
    # Eq. (8) from the Table S3 cyno values (E. coli reference, median CD4+) at 3.12 kg to a
    # 70 kg human; mL converted to L. The Results text and the Figure 4 caption both name the
    # final model with covariates, and with these values the model reproduces the Figure 4 inset
    # (Cmin, AUC and steady-state Cmax at all three doses); the Table 2 base values do not.
    lcl_exp_inf <- log(6.95 * 3.12 * (70 / 3.12)^0.85 / 1000); label("Constant nonsaturable clearance CL1 (L/day)") # Table S3 CL1 = 6.95 mL/day/kg; Methods b = 0.85 (Table 3: 4.4 mL/day/kg)
    # Table 3 prints the rituximab NHL CL2_0 as 18.0 mL/day/kg. Figure 4 is reproduced only when
    # 18.0 is applied as mL/day for the whole body (0.018 L/day); the per-kg reading
    # (18.0 x 70 = 1.26 L/day) gives cycle-1 AUCs about 2.4-fold below the inset.
    lcl_exp_component <- fixed(log(18.0 / 1000)); label("Initial value of the time-varying linear clearance CL2_0 (L/day)") # Table 3 rituximab NHL CL2_0 = 18.0, applied as mL/day (see comment above)
    lcl_exp_kdes <- fixed(log(0.023)); label("Decay constant lambda2 of the time-varying linear clearance (1/day)") # Table 3 rituximab NHL lambda2 = 0.023 1/day
    lq <- log(34.3 * 3.12 * (70 / 3.12)^0.85 / 1000); label("Distribution clearance CLd (L/day)") # Table S3 CLd,E.coli = 34.3 mL/day/kg; Methods b = 0.85
    lvc <- log(43.6 * 70 / 1000); label("Central volume V1 (L)") # Table S3 V1,E.coli = 43.6 mL/kg; Methods b = 1
    lvp <- log(40.4 * 70 / 1000); label("Peripheral volume V2 (L)") # Table S3 V2 = 40.4 mL/kg; Methods b = 1
    lvmax <- log(1760 * (70 / 3.12) / 1000); label("Initial Michaelis-Menten maximum elimination rate Vmax_0 (mg/day)") # Table S3 Vmax_0 = 1760 ug/day; Methods b = 1
    lvmax_exp_kdes <- fixed(log(0.167)); label("Decay constant lambda1 of the time-varying Vmax (1/day)") # Table S3 lambda1 = 0.167 1/day, fixed to the NHP value per Methods
    lkm <- log(31.1); label("Michaelis-Menten constant KM (ug/mL)") # Table S3 KM = 31.1 ug/mL, not scaled

    propSd <- fixed(0); label("Proportional residual error (fraction; not projected)")
  })
  model({
    cl_exp_inf <- exp(lcl_exp_inf)
    cl_exp_component <- exp(lcl_exp_component)
    cl_exp_kdes <- exp(lcl_exp_kdes)
    q <- exp(lq)
    vc <- exp(lvc)
    vp <- exp(lvp)
    vmax <- exp(lvmax)
    vmax_exp_kdes <- exp(lvmax_exp_kdes)
    km <- exp(lkm)

    # Eq. (1) and Figure 1, t = time since the first dose
    cl <- cl_exp_inf + cl_exp_component * exp(-cl_exp_kdes * time)
    vmax_t <- vmax * exp(-vmax_exp_kdes * time)

    Cc <- central / vc

    # Eqs. (1)-(2), IV dosing only
    d/dt(central) <- -(cl + q) / vc * central + q / vp * peripheral1 -
      vmax_t * Cc / (Cc + km)
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1

    Cc ~ prop(propSd)
  })
}
