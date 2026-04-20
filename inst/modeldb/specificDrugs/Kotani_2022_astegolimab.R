Kotani_2022_astegolimab <- function() {
  description <- "Two-compartment population PK model for astegolimab (anti-ST2 IgG2) in adults with severe asthma (Kotani 2022)"
  reference <- "Kotani N, Dolton M, Svensson RJ, et al. Population Pharmacokinetics and Exposure-Response Relationships of Astegolimab in Patients With Severe Asthma. J Clin Pharmacol. 2022;62(7):905-917. doi:10.1002/jcph.2021"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline-only (time-fixed per subject in Kotani 2022). Shared allometric exponent on CL and Q; separate shared exponent on Vc and Vp. Reference 79 kg (median in the PK analysis data set, Table 1).",
      source_name        = "BWT"
    ),
    eGFR = list(
      description        = "Baseline estimated glomerular filtration rate",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline-only (time-fixed per subject). Power effect on CL. Reference 87.9 mL/min/1.73 m^2 (median in the PK analysis data set, Table 1).",
      source_name        = "BEGFR"
    ),
    BEOS = list(
      description        = "Baseline blood eosinophil count",
      units              = "cells/uL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline-only (time-fixed per subject). Power effect on CL. Reference 180 cells/uL (median in the PK analysis data set, Table 1).",
      source_name        = "BEOS"
    ),
    DOSE_70MG = list(
      description        = "Indicator for the 70 mg SC Q4W dose group",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (210 mg or 490 mg Q4W regimen)",
      notes              = "Subject-level indicator in Kotani 2022: 1 = 70 mg dose arm, 0 = 210 mg or 490 mg dose arm (reference). Multiplicative -15.3% relative change on relative bioavailability.",
      source_name        = "Dose70mg"
    )
  )

  population <- list(
    n_subjects     = 368,
    n_studies      = 1,
    study_name     = "Zenyatta (NCT02918019)",
    age_range      = "18-75 years",
    age_median     = "53 years",
    weight_range   = "43-130 kg",
    weight_median  = "79 kg",
    sex_female_pct = 66.8,
    race_ethnicity = c(
      White             = 84.0,
      Black             = 5.4,
      Asian             = 4.6,
      `Native American` = 4.6,
      Multiple          = 1.4
    ),
    disease_state  = "Severe asthma on medium- or high-dose inhaled corticosteroid plus at least one additional controller, with at least one asthma exacerbation in the 12 months before screening.",
    dose_range     = "70, 210, or 490 mg SC every 4 weeks for 52 weeks",
    regions        = "Multinational Phase 2b trial (Zenyatta).",
    ada_pos_pct    = 7.3,
    notes          = "Baseline demographics and covariates from Kotani 2022 Table 1 (n = 368 patients with at least one quantifiable postdose PK observation out of 502 randomized in Zenyatta)."
  )

  ini({
    # Structural parameters - reference subject: 79 kg, eGFR 87.9 mL/min/1.73 m^2,
    # blood eosinophil 180 cells/uL, 210 or 490 mg Q4W dose (reference dose group).
    lka <- log(0.0437); label("Absorption rate constant (ka, 1/day)")                                 # Table 2 (ka = 0.0437 day^-1)
    lcl <- log(0.244);  label("Apparent clearance for reference subject (CL/F, L/day)")              # Table 2 (CL/F = 0.244 L/day)
    lvc <- log(0.614);  label("Apparent central volume of distribution for reference subject (Vc/F, L)")   # Table 2 (Vc/F = 0.614 L)
    lvp <- log(2.74);   label("Apparent peripheral volume of distribution for reference subject (Vp/F, L)") # Table 2 (Vp/F = 2.74 L)
    lq  <- log(0.171);  label("Apparent intercompartmental clearance for reference subject (Q/F, L/day)")   # Table 2 (Q/F = 0.171 L/day)

    # Shared body-weight allometric exponents (Kotani 2022 constrained one coefficient
    # on both CL and Q, and another on both Vc and Vp; Table 2 footnote b).
    allo_clq <- 0.986; label("Allometric exponent on CL and Q (unitless)")  # Table 2 (BWT on CL and Q = 0.986)
    allo_v   <- 1.02;  label("Allometric exponent on Vc and Vp (unitless)") # Table 2 (BWT on Vtot = 1.02)

    # Power effects on CL
    e_egfr_cl <- 0.431;  label("Power exponent for baseline eGFR on CL (unitless)")                    # Table 2 (BEGFR on CL = 0.431)
    e_eos_cl  <- 0.0905; label("Power exponent for baseline blood eosinophil count on CL (unitless)")  # Table 2 (BEOS on CL = 0.0905)

    # Relative bioavailability effect for the 70 mg dose arm (fractional change).
    e_dose70_frel <- -0.153; label("Relative change in Frel for the 70 mg dose group (fraction)")     # Table 2 (Dose70mg on F = -0.153)

    # Box-Cox transformation shape parameter applied to the IIV on Frel (Eq. 2).
    boxcox_frel <- -2.81; label("Box-Cox shape parameter for IIV on Frel (unitless)")                  # Table 2 (BoxCoxIIV,Frel = -2.81)

    # Inter-individual variability (reported on the CV scale for log-normal ka and CL;
    # omega^2 = log(CV^2 + 1)). IIV on Frel is reported as the SD of the pre-Box-Cox eta;
    # variance is eta SD^2 (the Box-Cox transformation is applied in model()).
    etalka ~ 0.20507  # CV 0.477 -> log(1 + 0.477^2)                                                  # Table 2 (IIV_Ka CV = 0.477)
    etalcl ~ 0.04896  # CV 0.224 -> log(1 + 0.224^2)                                                  # Table 2 (IIV_CL CV = 0.224)
    etafrel ~ 0.05905 # SD 0.243 (pre-Box-Cox eta); omega^2 = 0.243^2                                 # Table 2 (IIV_Frel SD = 0.243)

    # Residual error (combined additive + proportional; CV% for the proportional term,
    # concentration units for the additive term).
    propSd <- 0.198; label("Proportional residual error (fraction)")                                   # Table 2 (RUV proportional CV = 0.198)
    addSd  <- 0.603; label("Additive residual error (ug/mL)")                                          # Table 2 (RUV additive = 0.603 ug/mL)
  })
  model({
    # Covariate-adjusted multipliers on apparent clearance and volume parameters.
    # Reference values (medians, Table 1): WT = 79 kg, eGFR = 87.9 mL/min/1.73 m^2,
    # BEOS = 180 cells/uL.
    cl_cov <- (WT / 79)^allo_clq * (eGFR / 87.9)^e_egfr_cl * (BEOS / 180)^e_eos_cl
    q_cov  <- (WT / 79)^allo_clq
    v_cov  <- (WT / 79)^allo_v

    # Box-Cox-distributed IIV on Frel (Eq. 2 in Kotani 2022): the normal eta is
    # transformed before being exponentiated. Reference Frel is 1; the 70 mg
    # dose group incurs a fractional -0.153 change per Eq. 4.
    boxcox_eta_frel <- ((exp(etafrel))^boxcox_frel - 1) / boxcox_frel
    frel <- (1 + e_dose70_frel * DOSE_70MG) * exp(boxcox_eta_frel)

    # Individual PK parameters
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * cl_cov
    vc <- exp(lvc)          * v_cov
    vp <- exp(lvp)          * v_cov
    q  <- exp(lq)           * q_cov

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- frel

    # Concentration: dose in mg, volumes in L -> mg/L = ug/mL
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
