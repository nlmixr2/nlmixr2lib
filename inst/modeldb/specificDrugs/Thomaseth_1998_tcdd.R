Thomaseth_1998_tcdd <- function() {
  description <- paste(
    "Minimal physiologic toxicokinetic (MPTK) model of",
    "2,3,7,8-tetrachlorodibenzo-p-dioxin (TCDD) in adult humans",
    "(Thomaseth and Salvan 1998, Environ Health Perspect). A single",
    "whole-body TCDD burden compartment is eliminated by hepatic",
    "degradation at a rate proportional to the fraction of total body lipid",
    "held in the liver, so the apparent elimination rate (and half-life)",
    "tracks body composition. Adipose, liver and other-tissue lipid volumes",
    "are anthropometric functions of body-mass index (BMI) and height, and",
    "BMI drifts with age, giving a time-varying distribution volume. A daily",
    "TCDD intake proportional to body weight (background plus optional",
    "occupational exposure) drives the burden. The observable is",
    "lipid-adjusted serum TCDD (ppt = pg TCDD per g lipid). The model was",
    "fit to repeated serum-dioxin measures in Operation Ranch Hand Vietnam",
    "veterans and used to estimate occupational exposure in a NIOSH cohort",
    "of chemical-plant workers. Deterministic structure with additive",
    "between-subject variability on the hepatic elimination constant and",
    "proportional residual error on log-transformed concentrations."
  )
  reference <- paste(
    "Thomaseth K, Salvan A. Estimation of occupational exposure to",
    "2,3,7,8-tetrachlorodibenzo-p-dioxin using a minimal physiologic",
    "toxicokinetic model. Environ Health Perspect. 1998;106(Suppl",
    "2):743-753. doi:10.1289/ehp.98106743. Structural parent model:",
    "Dankovic DA, Andersen ME, Salvan A, Stayner L. A simplified model",
    "describing the kinetics of TCDD in humans [abstract]. Toxicologist.",
    "1995;15:272."
  )
  vignette <- "Thomaseth_1998_tcdd"

  units <- list(
    time = "day",
    dosing = "pg",
    concentration = "ppt",
    amount = "pg",
    weight = "kg"
  )

  covariateData <- list(
    HT = list(
      description = "Body height (constant over the adult follow-up)",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters the adipose-tissue formula (Knapik et al. 1983) via BMI and",
        "converts BMI to body weight (bw = BMI * (HT/100)^2). The paper's",
        "worked examples (Tables 3 and 7, Figure 6) use HT = 170 cm."
      ),
      source_name = "H"
    ),
    BMI = list(
      description = "Baseline body-mass index at the start of the simulation (initial condition for the time-varying BMI)",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "BMI at the first modelled time. It then drifts with age per the",
        "population BMI sub-model (Equations 20-21). NIOSH workers at hire:",
        "mean 24.6 kg/m^2 (Table 2)."
      ),
      source_name = "BMI"
    ),
    AGE_BASE = list(
      description = "Subject age at the start of the simulation",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline age paired with the time-varying age = AGE_BASE + t/365.25.",
        "Only used by the BMI age-drift sub-model. NIOSH workers at hire:",
        "mean 28.5 years (Table 2)."
      ),
      source_name = "age"
    )
  )

  # Single lumped whole-body TCDD burden state; not a canonical compartment name.
  paper_specific_compartments <- c("tcdd")

  compartmentData <- list(
    tcdd = list(
      analyte = "2,3,7,8-TCDD",
      units = "pg",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 277L,
    n_studies = 2L,
    age_range = "25-75 years",
    disease_state = "healthy adult males (occupational / environmental TCDD exposure)",
    sex_female_pct = 0,
    dose_range = "environmental background to occupational TCDD intake (pg/kg/day)",
    regions = "USA (Operation Ranch Hand veterans; NIOSH chemical-plant workers)",
    notes = paste(
      "Kinetic parameters (kf, background input) were estimated from repeated",
      "serum-TCDD measures in Operation Ranch Hand Vietnam veterans: 277",
      "exposed veterans with three samples (1982/1987/1992) plus 65 with a",
      "1987 level, and 33 unexposed reference-group observations (Table 1).",
      "The occupational-exposure parameter was then estimated by linear",
      "regression on a NIOSH subcohort of 253 male chemical-plant workers",
      "with a single serum-TCDD measure (Table 2); the full NIOSH cohort is",
      "5172 workers at 12 US plants. All subjects male."
    )
  )

  ini({
    # -------------------------------------------------------------------
    # Estimated kinetic parameters (nonlinear mixed-effects, log-transformed
    # data: Table 4, 'NCME and log transformation' row).
    # -------------------------------------------------------------------
    kf <- 0.02199
    label("Hepatic (liver-lipid) TCDD elimination rate constant kf (1/day)") # Table 4 NCME row = 0.02199 (95% CI 0.020-0.024)

    input <- 0.1251
    label("Background daily TCDD intake per kg body weight (pg/kg/day)") # Table 4 NCME row = 0.1251 (95% CI 0.071-0.179); NIOSH scenario re-fit to 0.293 (Table 6, to hold 7 ppt background)

    exposure <- 0
    label("Occupational daily TCDD intake per kg body weight (pg/kg/day); 0 = background only") # Table 6 NIOSH occupational estimate = 232.7 (95% CI 192-273); default 0 for a background individual

    # -------------------------------------------------------------------
    # Population BMI age-drift sub-model (Equation 20, dBMI/dt = alphaBMI*age
    # + betaBMI); estimates from the NIOSH subcohort (Figure 3 text / Fig 4).
    # -------------------------------------------------------------------
    alphaBMI <- -3.755e-3
    label("BMI age-drift quadratic coefficient alpha (kg/m^2/year^2)") # Figure 3 text: alphaBMI = -3.755e-3 +/- 0.9e-3 (SE)

    betaBMI <- 0.26907
    label("BMI age-drift linear coefficient beta (kg/m^2/year)") # Figure 3 text: betaBMI = 0.26907 +/- 0.04 (SE)

    # -------------------------------------------------------------------
    # Fixed anatomical / lipid constants.
    # -------------------------------------------------------------------
    kadipose1 <- fixed(0.01264)
    label("Adipose-mass BMI slope, g adipose per kg bw per unit BMI (unitless)") # page 746: v_adipose = 1000*(0.01264*BMI - 0.13305); Knapik et al. 1983

    kadipose0 <- fixed(0.13305)
    label("Adipose-mass intercept (unitless)") # page 746: v_adipose = 1000*(0.01264*BMI - 0.13305)

    fliverlean <- fixed(0.0311)
    label("Liver mass as a fraction of lean body mass (unitless)") # page 745: V_liver = 3.11% of lean body mass (ICRP)

    flipadipose <- fixed(0.8)
    label("Lipid fraction of adipose tissue (unitless)") # page 745: LV_adipose = 0.8 * V_adipose (80% of adipose mass)

    flipliver <- fixed(0.069)
    label("Lipid fraction of liver tissue (unitless)") # page 745: LV_liver = 0.069 * V_liver (6.9%)

    flipother <- fixed(0.022)
    label("Lipid fraction of other tissue (unitless)") # page 745: LV_other = 0.022 * V_other (2.2%)

    ladj0 <- fixed(7)
    label("Assumed lipid-adjusted serum TCDD at the start of the simulation (ppt)") # page 745/749: initial condition fixed to 7 ppt (mean unexposed level)

    # -------------------------------------------------------------------
    # Between-subject variability on kf: ADDITIVE normal (Equation 18,
    # kf(i) = kf + e_kf, e_kf ~ N(0, sigma2_kf)), NOT log-normal.
    # -------------------------------------------------------------------
    etakf ~ 2.9e-5 # Table 4 footnote: 'sigma2_kf = 2.9e-5' (additive variance of kf; SD 0.00539, ~24% of kf)

    # -------------------------------------------------------------------
    # Residual error: additive on the log-transformed concentration
    # (Equation 19), equivalent to proportional in linear space.
    # -------------------------------------------------------------------
    propSd <- 0.259
    label("Proportional residual error (fraction)") # Table 4 footnote: 'sigma2_eps = 0.0671' on log scale; sqrt = 0.259 (25.9% CV, stated page 749)
  })

  model({
    # Current age (years) and the time-varying BMI (closed form of Eq 20-21,
    # integrated from AGE_BASE). t is the ODE time in days.
    age <- AGE_BASE + t / 365.25
    bmi <- BMI + alphaBMI / 2 * (age^2 - AGE_BASE^2) + betaBMI * (age - AGE_BASE)
    bw <- bmi * (HT / 100)^2 # body weight (kg)

    # Lipid volumes per kg body weight (g lipid / kg bw). The elimination
    # term uses the liver-to-total lipid ratio, in which body weight cancels.
    vadipose <- 1000 * (kadipose1 * bmi - kadipose0) # g adipose / kg bw
    vliver <- fliverlean * (1000 - vadipose) # g liver / kg bw
    vother <- 1000 - vadipose - vliver # g other / kg bw
    lvliver <- flipliver * vliver # g liver lipid / kg bw
    tlv <- flipadipose * vadipose + flipliver * vliver + flipother * vother # g total lipid / kg bw

    kfi <- kf + etakf # additive BSV on kf (Eq 18)

    # Initial whole-body burden X(t0) = ladj(t0) * TLV(t0) (Eq 2), using BMI
    # at t0 (= the BMI covariate).
    vadip0 <- 1000 * (kadipose1 * BMI - kadipose0)
    vliv0 <- fliverlean * (1000 - vadip0)
    voth0 <- 1000 - vadip0 - vliv0
    tlv0 <- flipadipose * vadip0 + flipliver * vliv0 + flipother * voth0
    bw0 <- BMI * (HT / 100)^2
    tcdd(0) <- ladj0 * tlv0 * bw0

    # Whole-body TCDD burden (pg), Equation 1. Elimination = kf * (liver
    # lipid / total lipid) * X; input = (background + occupational) * bw.
    d/dt(tcdd) <- -kfi * (lvliver / tlv) * tcdd + (input + exposure) * bw

    # Observable: lipid-adjusted serum TCDD (ppt = pg / g lipid) = X / TLV,
    # TLV = total lipid volume (g) = tlv * bw. The paper calls it ladj; Cc is
    # the canonical concentration output carrying the residual error.
    TLV <- tlv * bw
    ladj <- tcdd / TLV
    Cc <- ladj
    Cc ~ prop(propSd)
  })
}
