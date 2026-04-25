Yang_2021_cemiplimab <- function() {
  description <- "Two-compartment population PK model for cemiplimab (anti-PD-1 IgG4) with time-varying clearance (sigmoid Emax) in adults with advanced solid tumors including cutaneous squamous cell carcinoma (Yang 2021)"
  reference <- "Yang F, Paccaly AJ, Rippley RK, Davis JD, DiCioccio AT. Population pharmacokinetic characteristics of cemiplimab in patients with advanced malignancies. J Pharmacokinet Pharmacodyn. 2021 Aug;48(4):479-494. doi:10.1007/s10928-021-09739-y"
  vignette <- "Yang_2021_cemiplimab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on shared CL/Q with reference 76.2 kg and on shared V2/V3 with reference 76.2 kg (Yang 2021 Table 2 median weight; Eqs. for CL_i, Q_i, V2_i, V3_i in the Final PopPK model section).",
      source_name        = "WGTBL"
    ),
    ALB = list(
      description        = "Baseline serum albumin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on shared CL/Q with reference 38 g/L (Yang 2021 Table 2 median albumin; Eqs. for CL_i and Q_i).",
      source_name        = "ALBBL"
    ),
    IGG = list(
      description        = "Baseline endogenous serum immunoglobulin G concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on shared CL/Q with reference 9.65 g/L (Yang 2021 Table 2 median IgG; Eqs. for CL_i and Q_i).",
      source_name        = "IGGBL"
    ),
    ALT = list(
      description        = "Baseline serum alanine aminotransferase activity",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on shared CL/Q with reference 20 IU/L (Yang 2021 Table 2 median ALT; Eqs. for CL_i and Q_i).",
      source_name        = "ALTBL"
    ),
    BMI = list(
      description        = "Baseline body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on shared V2/V3 with reference 26.5 kg/m^2 (Yang 2021 Table 2 median BMI; Eqs. for V2_i and V3_i).",
      source_name        = "BMIBL"
    )
  )

  population <- list(
    n_subjects     = 548L,
    n_studies      = 2L,
    age_range      = "27-96 years (median 65)",
    age_median     = "65 years",
    weight_range   = "30.9-172 kg (median 76.2)",
    weight_median  = "76.2 kg",
    sex_female_pct = 39.6,
    race_ethnicity = c(White = 90.9, Black = 3.6, Asian = 1.6, Other = 3.8),
    disease_state  = "Advanced solid tumors (any), including advanced cutaneous squamous cell carcinoma (CSCC; metastatic or locally advanced)",
    dose_range     = "1, 3, or 10 mg/kg Q2W; 3 mg/kg Q3W; 200 mg Q2W; 350 mg Q3W (all 30-min IV infusion)",
    regions        = "Multinational (not enumerated in the paper)",
    notes          = "Pooled from Study 1423 (NCT02383212; first-in-human in advanced malignancies) and Study 1540 (NCT02760498; phase 2 in advanced CSCC). Final analysis set 11,178 PK observations from 548 patients (Yang 2021 Tables 1-2). Median labs: albumin 38 g/L, IgG 9.65 g/L, ALT 20 IU/L, BMI 26.5 kg/m^2 (Yang 2021 Table 2)."
  )

  ini({
    # Structural parameters - reference values for a typical patient with median
    # baseline covariates: WT 76.2 kg, ALB 38 g/L, IGG 9.65 g/L, ALT 20 IU/L,
    # BMI 26.5 kg/m^2 (Yang 2021 Final PopPK model section + Table 2 medians).
    lcl  <- log(0.290); label("Baseline clearance CL_BASE,REF for the paper's reference covariates (L/day)") # Yang 2021 Table 3: TVCL = 0.290 L/day
    lvc  <- log(3.32);  label("Central volume of distribution V2_REF (L)")                                   # Yang 2021 Table 3: TVV2 = 3.32 L
    lq   <- log(0.638); label("Inter-compartmental clearance Q_REF (L/day)")                                 # Yang 2021 Table 3: TVQ  = 0.638 L/day
    lvp  <- log(1.65);  label("Peripheral volume of distribution V3_REF (L)")                                # Yang 2021 Table 3: TVV3 = 1.65 L

    # Time-varying clearance: sigmoid Emax of time since first dose
    # (Yang 2021 Final PopPK model: CL_i = CL_BASE,REF * exp(Emax_i * T^HILL /
    # (T50_i^HILL + T^HILL)) * covariate-terms * exp(eta_clq)). Emax is the
    # log-fold maximal change; with Emax = -0.410 the steady-state CL approaches
    # exp(-0.410) = 0.664 of baseline (a 33.6% reduction), consistent with the
    # paper's reported ~35.9% mean reduction within 16 weeks of treatment
    # (Yang 2021 Base model section). Emax can be negative, so it is not log-
    # transformed in the model file; T50 is positive and is log-transformed.
    emax  <- -0.410;     label("Maximal log-fold change in CL EMAX_REF (unitless)")                         # Yang 2021 Table 3: Emax = -0.410
    lt50  <- log(28.9);  label("Time at which the change in CL is 50%% of EMAX, T50_REF (days)")            # Yang 2021 Table 3: T50  = 28.9 days
    hill  <- 2.79;       label("Hill / sigmoidicity exponent of time on CL (unitless)")                     # Yang 2021 Table 3: HILL = 2.79

    # Covariate effects on shared CL/Q (Yang 2021 Final PopPK model
    # Eqs. for CL_i and Q_i; all power-form). Same exponent on CL and Q.
    e_wt_clq  <-  0.477;   label("Power exponent of WT on shared CL/Q (unitless)")                          # Yang 2021 Table 3: WGT_ON_CLQ = 0.477
    e_alb_clq <- -0.926;   label("Power exponent of ALB on shared CL/Q (unitless)")                         # Yang 2021 Table 3: ALB_ON_CLQ = -0.926
    e_igg_clq <-  0.184;   label("Power exponent of IGG on shared CL/Q (unitless)")                         # Yang 2021 Table 3: IGG_ON_CLQ = 0.184
    e_alt_clq <- -0.0795;  label("Power exponent of ALT on shared CL/Q (unitless)")                         # Yang 2021 Table 3: ALT_ON_CLQ = -0.0795

    # Covariate effects on shared V2/V3 (Yang 2021 Final PopPK model
    # Eqs. for V2_i and V3_i; all power-form). Same exponents on V2 and V3.
    e_wt_vss  <-  0.970;   label("Power exponent of WT on shared V2/V3 (unitless)")                         # Yang 2021 Table 3: WGT_ON_VSS = 0.970
    e_bmi_vss <- -0.560;   label("Power exponent of BMI on shared V2/V3 (unitless)")                        # Yang 2021 Table 3: BMI_ON_VSS = -0.560

    # Inter-individual variability (Yang 2021 Table 3). The paper estimated a
    # SHARED IIV on CL and Q (etalcl, used in the eta of both CL_i and Q_i)
    # and a SHARED IIV on V2 and V3 (etalvc, used in the eta of both V2_i and
    # V3_i), with off-diagonal covariance between the CLQ and V_ss etas. IIV
    # on Emax and T50 are independent.
    etalcl + etalvc ~ c(0.0870,
                        0.0422, 0.0432)                                                                    # Yang 2021 Table 3: omega^2_CLQ = 0.0870, cov_CLQ:VSS = 0.0422, omega^2_VSS = 0.0432
    etaemax ~ 0.228                                                                                        # Yang 2021 Table 3: omega^2_Emax = 0.228 (Emax_i = EMAX_REF * exp(eta))
    etalt50 ~ 0.610                                                                                        # Yang 2021 Table 3: omega^2_T50  = 0.610 (T50_i  = T50_REF  * exp(eta))

    # Residual error (Yang 2021 Final PopPK model: log-transformed
    # additive-plus-proportional, Y = F + F*ERR(1) + ERR(2)). RUVCV maps
    # to nlmixr2 propSd; RUVSD maps to addSd in the same units as Cc.
    propSd <- 0.188;  label("Proportional residual error (fraction)")                                       # Yang 2021 Table 3: RUVCV = 0.188
    addSd  <- 1.48;   label("Additive residual error (mg/L)")                                               # Yang 2021 Table 3: RUVSD = 1.48 mg/L
  })
  model({
    # Individual baseline parameters with covariate adjustments (Yang 2021
    # Final PopPK model Eqs.). The same etalcl is used on CL and Q (shared
    # eta on CL/Q) and the same etalvc is used on V2 and V3 (shared eta on
    # V2/V3). Reference values (denominators below) come from Yang 2021
    # Table 2 medians: WT 76.2 kg, ALB 38 g/L, IGG 9.65 g/L, ALT 20 IU/L,
    # BMI 26.5 kg/m^2.
    cl_base <- exp(lcl + etalcl) *
      (WT  / 76.2)^e_wt_clq *
      (ALB / 38)^e_alb_clq *
      (IGG / 9.65)^e_igg_clq *
      (ALT / 20)^e_alt_clq

    q <- exp(lq + etalcl) *
      (WT  / 76.2)^e_wt_clq *
      (ALB / 38)^e_alb_clq *
      (IGG / 9.65)^e_igg_clq *
      (ALT / 20)^e_alt_clq

    vc <- exp(lvc + etalvc) *
      (WT  / 76.2)^e_wt_vss *
      (BMI / 26.5)^e_bmi_vss

    vp <- exp(lvp + etalvc) *
      (WT  / 76.2)^e_wt_vss *
      (BMI / 26.5)^e_bmi_vss

    # Time-varying clearance (Yang 2021 Final PopPK model Eq. for CL_i).
    # Multiplicative IIV on Emax (Emax_i = EMAX_REF * exp(eta)) and on T50
    # (T50_i = T50_REF * exp(eta)).
    emax_i <- emax * exp(etaemax)
    t50_i  <- exp(lt50 + etalt50)
    cl <- cl_base * exp(emax_i * t^hill / (t50_i^hill + t^hill))

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg and volumes in L -> central/vc has units mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
