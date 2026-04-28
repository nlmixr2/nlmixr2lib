LeTilly_2021_trastuzumab <- function() {
  description <- "Two-compartment serum/CSF population PK model for trastuzumab after intrathecal and intravenous administration in adults with HER2+ breast cancer leptomeningeal metastases (Le Tilly 2021); zero-order serum-to-CSF transfer plus first-order CSF-to-serum return, with a Friberg-style chain of latent target (HER2) transit compartments and irreversible binding-driven elimination of trastuzumab in the CSF compartment."
  reference <- "Le Tilly O, Azzopardi N, Bonneau C, Ohresser M, Ternant D, Thomas K, Olivier F, Trouillas I, Etcheverry M, Demarquay C, Garcia M, Paintaud G, Goupille O. Antigen Mass May Influence Trastuzumab Concentrations in Cerebrospinal Fluid After Intrathecal Administration. Clin Pharmacol Ther. 2021;110(1):210-219. doi:10.1002/cpt.2188"
  vignette <- "LeTilly_2021_trastuzumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    # Le Tilly 2021 tested body weight, glycorrhachia, and presence of an
    # intrathecal drug delivery device as candidate covariates (Methods,
    # "Influence of covariates"). None were retained in the final model
    # (Table 2 reports no covariate effects), so the implementation has no
    # covariate inputs.
  )

  population <- list(
    n_subjects     = 21L,
    n_studies      = 1L,
    n_observations = 304L,
    age_range      = "24-66 years",
    age_median     = "52 years",
    weight_range   = "38-90 kg",
    weight_median  = "65 kg",
    sex_female_pct = NA_real_,
    disease_state  = "HER2-positive breast cancer with leptomeningeal carcinomatosis (LMC).",
    dose_range     = "Weekly intrathecal trastuzumab 30, 60, 100, or 150 mg for up to 8 doses (n=21). 13/21 (62%) also received concurrent intravenous trastuzumab 6 mg/kg every 3 weeks; the IV infusion (30 minutes) was always given after the IT dose on shared days.",
    regions        = "Multicentric phase I/II clinical trial in France (NCT01373710).",
    administration_routes = "Intrathecal route via lumbar puncture (n=7), Ommaya reservoir (n=9), or indwelling intrathecal drug delivery device (n=11); some patients used multiple routes across the 8-week treatment period (Le Tilly 2021 Table 1).",
    samples        = "150 CSF samples and 154 serum samples; 39 CSF (26%) excluded after a Grubbs test on the CSF:serum ratio identified washout-induced overestimates from indwelling devices. Concentrations below LLOQ (0.144 mg/L CSF; 0.241 mg/L serum) were retained as censored values (Le Tilly 2021 Table 1).",
    notes          = "21 adult patients (presumed predominantly or all female given HER2+ breast cancer LMC; sex not tabulated in Le Tilly 2021 Table 1). Baseline demographics from Le Tilly 2021 Table 1. All patients received systematic corticosteroid prophylaxis (>=20 mg/day prednisolone or equivalent for at least 3 days before each IT injection plus 25 mg IT hydrocortisone hemisuccinate immediately before each IT trastuzumab dose). Software: MONOLIX 2019R1 SAEM."
  )

  ini({
    # Structural PK parameters (Le Tilly 2021 Table 2, final model).
    # Apparent serum and CSF volumes are close to physiological volumes.
    lvc   <- log(3.25);   label("Apparent serum (V1) volume of distribution (L)")                     # Le Tilly 2021 Table 2: V1 = 3.25 L
    lcl   <- log(0.139);  label("Linear serum clearance CL (L/day)")                                  # Le Tilly 2021 Table 2: CL = 0.139 L/day
    lvp   <- log(0.644);  label("Apparent CSF (V2) volume of distribution (L)")                       # Le Tilly 2021 Table 2: V2 = 0.644 L
    lk_f2s <- log(0.311); label("First-order CSF-to-serum transfer rate (1/day)")                     # Le Tilly 2021 Table 2: k21 = 0.311 day^-1
    lk_s2f <- log(0.264); label("Zero-order serum-to-CSF transfer rate (mg/day); saturated flow approximation") # Le Tilly 2021 Table 2: k12 = 0.264 mg/day

    # Latent target (HER2) production / transit / binding parameters.
    # The latent target is a phenomenological variable whose dynamics absorb
    # tumor burden, HER2 release on cell lysis, and distribution-volume
    # variations (Le Tilly 2021 Discussion); kdeg uses nmol^-1 as written
    # in Le Tilly 2021 Table 2, applied identically in dCSF/dt and dL/dt
    # equations (3) and (7) of the paper - mass-balance stoichiometry is
    # not strictly preserved because L is a phenomenological latent species,
    # not a measured antigen mass.
    lkin   <- log(11.88); label("Zero-order latent HER2 production rate (nmol/day)")                  # Le Tilly 2021 Table 2: kin = 11.88 nmol/day
    lktr   <- log(0.325); label("Latent HER2 transit / output rate (1/day); ktr = kout (constrained equal to avoid overparameterization)") # Le Tilly 2021 Table 2: ktr = 0.325 day^-1
    lkdeg  <- log(0.0116); label("Second-order trastuzumab-HER2 binding/degradation rate kdeg (nmol^-1 day^-1)") # Le Tilly 2021 Table 2: kdeg = 0.0116 nmol^-1 day^-1

    # Inter-individual variability. MONOLIX exponential IIV: theta_i =
    # theta_pop * exp(eta_i), with eta_i ~ N(0, omega^2). Le Tilly 2021
    # Table 2 reports SDs of the random effects directly (omega), so the
    # variances below are the squared values.
    etalvc   ~ 0.469225  # omega_V1  = 0.685, var = 0.685^2 -- Le Tilly 2021 Table 2
    etalcl   ~ 0.114921  # omega_CL  = 0.339, var = 0.339^2 -- Le Tilly 2021 Table 2
    etalk_f2s ~ 0.665856 # omega_k21 = 0.816, var = 0.816^2 -- Le Tilly 2021 Table 2
    etalkin  ~ 0.620944  # omega_kin = 0.788, var = 0.788^2 -- Le Tilly 2021 Table 2
    etalktr  ~ 1.168561  # omega_ktr = 1.081, var = 1.081^2 -- Le Tilly 2021 Table 2

    # Residual error. Le Tilly 2021 reports proportional error for both
    # serum and CSF concentrations, "Best error model was of proportional
    # type for both serum and CSF trastuzumab concentrations" (Results),
    # with values from Table 2 in % CV. Stored as fractions in nlmixr2's
    # linear-space proportional parameterization.
    CcpropSd   <- 0.2094; label("Serum proportional residual error (fraction)")  # Le Tilly 2021 Table 2: sigma_prop,serum = 20.94%
    CcsfpropSd <- 0.5409; label("CSF proportional residual error (fraction)")    # Le Tilly 2021 Table 2: sigma_prop,CSF   = 54.09%
  })

  model({
    # ---- Individual PK parameters --------------------------------------
    vc      <- exp(lvc + etalvc)
    cl      <- exp(lcl + etalcl)
    vp      <- exp(lvp)                          # no IIV retained on V2 per Le Tilly 2021 Table 2
    k_f2s   <- exp(lk_f2s + etalk_f2s)           # paper's k21 (first-order, 1/day)
    k_s2f   <- exp(lk_s2f)                       # paper's k12 (zero-order, mg/day)
    kin     <- exp(lkin + etalkin)
    ktr     <- exp(lktr + etalktr)
    kout    <- ktr                               # ktr and kout fixed equal per Le Tilly 2021 Methods
    kdeg    <- exp(lkdeg)
    kel     <- cl / vc                           # serum elimination rate constant (1/day)

    # Latent target initial condition (steady state of dL0/dt = kin - ktr*L0
    # at baseline gives L0 = kin/ktr; the paper sets every transit and the
    # final compartment to the same baseline kin/kout = kin/ktr).
    lat_init <- kin / kout

    lat0(0) <- lat_init
    lat1(0) <- lat_init
    lat2(0) <- lat_init
    lat(0)  <- lat_init

    # Negative feedback on latent target production (Le Tilly 2021
    # equation (8)): Kin(t) = kin * (L(t=0) / L(t))^gamma, with gamma = 1
    # (fixed after sensitivity analysis - "the power parameter gamma was
    # not identifiable and was therefore fixed to 1").
    kin_t <- kin * lat_init / lat

    # ---- ODE system ----------------------------------------------------
    # central tracks serum amount (mg); csf tracks CSF amount (mg). Le Tilly
    # 2021 equations (1) and (3) (serum / CSF), (5)-(7) (latent target
    # chain). The zero-order serum-to-CSF flux k_s2f is a saturated-flow
    # approximation (Discussion: "may suggest the presence of a saturable
    # phenomenon (saturated at the concentrations studied)"); it remains
    # constant in mg/day and is intended for use only at simulated trough
    # concentrations achieved with the studied IT/IV dosing.
    d/dt(central) <- -k_s2f + k_f2s * csf - kel * central
    d/dt(csf)     <-  k_s2f - k_f2s * csf - kdeg * csf * lat

    d/dt(lat0)    <-  kin_t - ktr  * lat0
    d/dt(lat1)    <-  ktr   * lat0 - kout * lat1
    d/dt(lat2)    <-  ktr   * lat1 - kout * lat2
    d/dt(lat)     <-  ktr   * lat2 - kout * lat - kdeg * lat * csf

    # ---- Observations and residual error -------------------------------
    Cc   <- central / vc                         # serum concentration (mg/L = ug/mL)
    Ccsf <- csf     / vp                         # CSF concentration   (mg/L = ug/mL)

    Cc   ~ prop(CcpropSd)
    Ccsf ~ prop(CcsfpropSd)
  })
}
