Ferron_2013_cabazitaxel <- function() {
  description <- "Three-compartment population PK model for intravenous cabazitaxel in patients with advanced solid tumors (Ferron 2013)"
  reference <- "Ferron GM, Dai Y, Semiond D. Population pharmacokinetics of cabazitaxel in patients with advanced solid tumors. Cancer Chemother Pharmacol. 2013;71(3):681-692. doi:10.1007/s00280-012-2058-9"
  vignette <- "Ferron_2013_cabazitaxel"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear scaling on CL with reference 1.84 m^2 (median of pooled cohort, Ferron 2013 Table 3 footnote and Equation 9). Ferron 2013 does not report which BSA formula was used; assume DuBois (commonly applied in oncology).",
      source_name        = "BSA"
    ),
    TUMTP_BC = list(
      description        = "Tumor-type indicator for breast cancer (TT1 in the source)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (all other tumor types: prostate, gastrointestinal, other)",
      notes              = "Multiplicative reduction on CL for breast-cancer patients (Ferron 2013 Equation 9: CL = 48.5 * BSA/1.84 * (1 - 0.543 * TT1)). Source column TT1 is renamed to canonical TUMTP_BC per inst/references/covariate-columns.md. The breast-cancer effect is confounded with study (34 of 37 breast-cancer patients came from a single Phase II study, ARD6191); see vignette Assumptions and deviations.",
      source_name        = "TT1"
    )
  )

  population <- list(
    n_subjects     = 170L,
    n_studies      = 5L,
    age_range      = "25-83 years",
    age_median     = "63 years",
    weight_range   = "35-133 kg",
    weight_median  = "74.8 kg",
    bsa_range      = "1.30-2.53 m^2",
    bsa_median     = "1.84 m^2",
    bmi_range      = "13.3-40.5 kg/m^2",
    bmi_median     = "26.2 kg/m^2",
    sex_female_pct = 35.3,
    race_ethnicity = c(Caucasian = 84.7, Black = 2.4, Oriental = 5.3, Hispanic = 4.1, Other = 3.5),
    disease_state  = "Advanced solid tumors (prostate 45.3%, breast 21.8%, gastrointestinal 13.5%, other 19.4%); includes the Phase III TROPIC trial in metastatic castration-resistant prostate cancer refractory to docetaxel.",
    dose_range     = "10-30 mg/m^2 cabazitaxel as a 1-h IV infusion every 3 weeks or once weekly for the first 4 weeks of a 5-week treatment cycle",
    regions        = "Multinational (TROPIC trial NCT00417079 plus four supporting studies)",
    renal_function = "CrCl median 89.2 mL/min; 59 mild renal impairment (50 <= CrCl <= 80), 14 moderate (30 <= CrCl < 50), 1 severe (CrCl < 30).",
    co_medication  = "Concomitant CYP inducers (mainly prednisone/prednisolone): 0% in cycle 1; up to 94% in subsequent cycles of the EFC6193 (TROPIC) study.",
    sampling       = "2,322 measurable cabazitaxel plasma concentrations from 170 patients with 4-50 sampling points across 1-3 cycles per patient (Ferron 2013 Tables 1, 2, 3 and Methods). LLOQ 1.0 ng/mL; concentrations below LLOQ excluded from the analysis.",
    notes          = "Pooled five studies (TED6188 Phase I, TED6189 Phase I, TED6190 Phase I, ARD6191 Phase II in taxane-resistant breast cancer, EFC6193/TROPIC Phase III in mCRPC). Baseline demographics from Ferron 2013 Table 3."
  )

  ini({
    # Structural PK -- final population PK model (Ferron 2013 Table 4).
    # Parameterization: ADVAN11/TRANS1-style hybrid where CL and V1 are
    # estimated alongside the four intercompartmental first-order rate
    # constants K12, K21, K13, K31 (h^-1). Derived elimination rate
    # K = K10 = CL / V1. Time unit hours; concentration ng/mL.
    lcl  <- log(48.5);   label("Clearance for non-breast-cancer at BSA 1.84 m^2 (CL, L/h)") # Ferron 2013 Table 4, theta1; Eq. 9
    lvc  <- log(26.0);   label("Central volume of distribution (V1, L)")                    # Ferron 2013 Table 4, theta2
    lk12 <- log(2.48);   label("Rate constant central -> peripheral1 (K12, 1/h)")           # Ferron 2013 Table 4, theta3
    lk21 <- log(0.604);  label("Rate constant peripheral1 -> central (K21, 1/h)")           # Ferron 2013 Table 4, theta4
    lk13 <- log(4.84);   label("Rate constant central -> peripheral2 (K13, 1/h)")           # Ferron 2013 Table 4, theta5
    lk31 <- log(0.0266); label("Rate constant peripheral2 -> central (K31, 1/h)")           # Ferron 2013 Table 4, theta6

    # Covariate effect on CL: CL = TVCL * BSA/1.84 * (1 - 0.543 * TUMTP_BC)
    # The 95% CI on the breast-cancer fractional reduction was 0.217-0.869
    # (Ferron 2013 Table 4, theta7); the authors note this finding is
    # confounded with study because 34 of 37 breast-cancer patients came
    # from the single Phase II study ARD6191 (Discussion).
    e_tumtp_bc_cl <- 0.543; label("Fractional reduction in CL for breast-cancer patients (unitless)") # Ferron 2013 Table 4, theta7; Eq. 9

    # IIV (log-normal) on CL, V1, K12, K13, K31 (no IIV on K21).
    # Convert CV% from Table 4 to log-scale variance via omega^2 = log(1 + CV^2).
    etalcl  ~ 0.13987 # Ferron 2013 Table 4: IIV CL 38.8% CV -> log(1 + 0.388^2)
    etalvc  ~ 0.62427 # Ferron 2013 Table 4: IIV V1 93.4% CV -> log(1 + 0.934^2)
    etalk12 ~ 0.53063 # Ferron 2013 Table 4: IIV K12 84.0% CV -> log(1 + 0.840^2)
    etalk13 ~ 0.34203 # Ferron 2013 Table 4: IIV K13 64.2% CV -> log(1 + 0.642^2)
    etalk31 ~ 0.07647 # Ferron 2013 Table 4: IIV K31 28.2% CV -> log(1 + 0.282^2)

    # Proportional residual error 27.8% CV (Ferron 2013 Table 4 Eq. 2: Cobs = Cpred * (1 + eps))
    propSd <- 0.278; label("Proportional residual error (fraction)") # Ferron 2013 Table 4
  })
  model({
    # Individual PK parameters (Ferron 2013 Eq. 9 for CL covariate effect)
    cl  <- exp(lcl + etalcl) * (BSA / 1.84) * (1 - e_tumtp_bc_cl * TUMTP_BC)
    vc  <- exp(lvc + etalvc)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21)
    k13 <- exp(lk13 + etalk13)
    k31 <- exp(lk31 + etalk31)

    kel <- cl / vc

    # 3-compartment IV model (no depot; cabazitaxel is given as a 1-h IV
    # infusion administered into the central compartment via the rate column
    # on dose rows).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-                    k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-                                                      k13 * central - k31 * peripheral2

    # Concentration: dose in mg, V1 in L -> mg/L = ug/mL = 1000 ng/mL.
    # Multiply by 1000 to express the central concentration in ng/mL, the
    # unit reported by the paper (LLOQ 1.0 ng/mL; Table 4 residual error
    # interpreted in ng/mL).
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
