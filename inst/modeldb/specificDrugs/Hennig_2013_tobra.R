Hennig_2013_tobra <- function() {
  description <- "Two-compartment intravenous population PK model for tobramycin in adults and children with and without cystic fibrosis (Hennig 2013); fat-free mass allometric scaling on CL/Q (estimated exponent) and on V1/V2 (linear), sex-specific reference CL and V1, piecewise-linear age effect on CL with breakpoint at 18 years, and a power effect of the SCR_mean/SCR ratio on CL."
  reference <- "Hennig S, Standing JF, Staatz CE, Thomson AH. Population pharmacokinetics of tobramycin in patients with and without cystic fibrosis. Clin Pharmacokinet. 2013;52(4):289-301. doi:10.1007/s40262-013-0036-y"
  vignette <- "Hennig_2013_tobra"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass (Janmahasatian formula)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Estimated via the Janmahasatian formula (paper reference 20). Reference 70 kg. Allometric exponent on CL and Q is estimated (theta_FFM = 0.952); on V1 and V2 the exponent is fixed at 1 (linear scaling).",
      source_name        = "FFM"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at study entry. Piecewise-linear effect on CL with breakpoint at 18 years: slope theta_AGE = -0.021 per year for AGE <= 18 and -0.010 per year for AGE > 18. Reference age is 18 years (f_age = 1 at AGE = 18).",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Hennig 2013 reports separate typical population values for CL and V1 in females (8.1 L/h/70 kg, 20.1 L/70 kg) and males (9.4 L/h/70 kg, 25.1 L/70 kg). The model selects the sex-specific reference value based on SEXF; the canonical SEXF (1 = female) maps directly to the paper's 'female sex' indicator (Table 2 footnote 'F female').",
      source_name        = "SEXF"
    ),
    CREAT = list(
      description        = "Measured serum creatinine concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Patient's measured SCR. Per Hennig 2013 Methods 'Covariate Models', SCR was floored at 60 umol/L during covariate building (CREAT < 60 was set to 60 to avoid overestimating CL in patients with sub-reference SCR). Missing SCR was set to the dataset average 62.1 umol/L (4.2 % of patients). The effect on CL is f_scr = (CREAT_REF / CREAT)^theta_SCR.",
      source_name        = "SCR"
    ),
    CREAT_REF = list(
      description        = "Sex-, age- and size-adjusted normal-mean serum creatinine for the individual",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Externally-computed reference SCR for the patient (denoted SCR_mean in Hennig 2013 Eq. 5 / Eq. 6). Per the paper this was 'derived according to the relationships previously described' citing Ceriotti 2008 (Clin Chem) and Junge 2004 (Clin Chim Acta) for paediatrics/adults plus Johansson 2011 (Ther Drug Monit) for the algorithmic aggregation; the paper does not state the exact formula. Users must compute CREAT_REF from age, sex and (where the chosen formula uses it) body size before passing to the model. When no covariate value can be derived, set CREAT_REF = CREAT so the renal-function factor f_scr = 1 (matches the paper's 'covariate set to 1 for missing data' rule).",
      source_name        = "SCR_mean"
    )
  )

  population <- list(
    n_subjects     = 732L,
    n_studies      = "8 centres pooled (5 prior published tobramycin studies plus retrospective TDM data from Royal Children's Hospital Brisbane, Cincinnati Children's Hospital, and Gartnavel General Hospital Glasgow)",
    n_observations = 5605L,
    age_range      = "0.01-85 years (paediatric 0.5 weeks to 17.9 years; adults 18-85 years)",
    age_median     = "paediatric 7.68 years (CF 11.1, non-CF 5.0); adult 31.7 years (CF 24.3, non-CF 52.0)",
    weight_range   = "3.3-120.0 kg",
    weight_median  = "paediatric 25.5 kg (CF 31.8, non-CF 18.6); adult 58.0 kg (CF 53.9, non-CF 67.0)",
    ffm_range      = "3.0-65.1 kg",
    ffm_median     = "paediatric 19.8 kg; adult 43.5 kg",
    sex_female_pct = "paediatric 53 % (207/391 known); adult 48 % (99/208)",
    disease_state  = "Mixed: 465 patients with cystic fibrosis (351 children, 114 adults) and 267 without cystic fibrosis (173 children including febrile-neutropenia oncology patients, 94 adults from a heterogeneous adult cohort)",
    dose_range     = "Paediatric median 10.0 mg/kg/day (range 1.7-28.8); adult median 5.2 mg/kg/day (range 0.9-12.0). Mixed once-daily, twice-daily and three-times-daily regimens (Table 1).",
    administration = "Intravenous bolus injection (97 patients) or short intravenous infusion (635 patients); paediatric infusions via burette include a fixed lag time representing transit through the administration device (paper Methods 2.3.2).",
    regions        = "Australia, USA, UK",
    notes          = "Demographics from Hennig 2013 Table 1. SCR_CR median (mL/min) varies by sub-cohort (paediatric 84.1, adult 71.5); SCR median (umol/L) paediatric 44.7, adult 68.3 (Table 1)."
  )

  ini({
    # Reference (female) typical values for a 70 kg-FFM adult, 18 years,
    # CREAT = CREAT_REF (Hennig 2013 Table 2, "Final model" column). The
    # paper reports separate female / male point estimates for CL and V1;
    # here female is the reference, with a log-scale male multiplicative
    # effect e_male_cl / e_male_vc that recovers the published male values
    # (CL_male = exp(lcl + e_male_cl) = 9.4; V1_male = exp(lvc + e_male_vc) = 25.1).
    lcl <- log(8.1);  label("Typical CL for a 70 kg-FFM female adult at 18 y (L/h)")  # Hennig 2013 Table 2: theta_CL,female = 8.1 L/h/70 kg
    lvc <- log(20.1); label("Typical V1 for a 70 kg-FFM female adult (L)")            # Hennig 2013 Table 2: theta_V1,female = 20.1 L/70 kg
    lq  <- log(1.5);  label("Typical inter-compartmental clearance Q for a 70 kg-FFM adult (L/h)") # Hennig 2013 Table 2: theta_Q2 = 1.5 L/h/70 kg
    lvp <- log(10.0); label("Typical peripheral volume V2 for a 70 kg-FFM adult (L)") # Hennig 2013 Table 2: theta_V2 = 10.0 L/70 kg

    # Sex effects (male relative to the female reference; applied via (1 - SEXF))
    e_male_cl <- log(9.4 / 8.1);   label("Log-scale male sex effect on CL (unitless; applied as (1 - SEXF))")  # derived from Hennig 2013 Table 2: log(theta_CL,male / theta_CL,female) = log(9.4 / 8.1)
    e_male_vc <- log(25.1 / 20.1); label("Log-scale male sex effect on V1 (unitless; applied as (1 - SEXF))")  # derived from Hennig 2013 Table 2: log(theta_V1,male / theta_V1,female) = log(25.1 / 20.1)

    # Allometric / covariate parameters
    e_ffm_cl_q <- 0.952;  label("Allometric exponent of FFM/70 on CL and Q (unitless)")          # Hennig 2013 Table 2: theta_FFM = 0.952
    e_age_le18 <- -0.021; label("Slope of (AGE - 18) on f_age for AGE <= 18 (per year)")         # Hennig 2013 Table 2: theta_AGE (<18 years) = -0.021
    e_age_gt18 <- -0.010; label("Slope of (AGE - 18) on f_age for AGE >  18 (per year)")         # Hennig 2013 Table 2: theta_AGE (>18 years) = -0.010
    e_scr_cl   <- 0.222;  label("Power exponent of (CREAT_REF / CREAT) on CL (unitless)")        # Hennig 2013 Table 2: theta_SCR = 0.222

    # Inter-individual variability. Hennig 2013 Table 2 reports BSV CV% on
    # CL = 25.9 %, V1 = 15.2 %, Q2 = 41.8 %, V2 = 58.5 % with correlations
    # CL-V1 = 65.8 %, CL-Q2 = 71.1 %, V1-Q2 = 47.5 %. V2 has no listed
    # correlation with the other parameters and is treated as independent.
    # Variances on the log scale: omega^2 = log(CV^2 + 1).
    #   CL : log(0.259^2 + 1) = 0.06492
    #   V1 : log(0.152^2 + 1) = 0.02283
    #   Q  : log(0.418^2 + 1) = 0.16104
    #   V2 : log(0.585^2 + 1) = 0.29443
    # Covariances: cov(X,Y) = r(X,Y) * sqrt(var(X) * var(Y))
    #   cov(CL,V1) = 0.658 * sqrt(0.06492 * 0.02283) = 0.02534
    #   cov(CL,Q ) = 0.711 * sqrt(0.06492 * 0.16104) = 0.07270
    #   cov(V1,Q ) = 0.475 * sqrt(0.02283 * 0.16104) = 0.02881
    etalcl + etalvc + etalq ~ c(
      0.06492,
      0.02534, 0.02283,
      0.07270, 0.02881, 0.16104
    )                                                                                            # Hennig 2013 Table 2: BSV CV% and correlations on CL, V1, Q2
    etalvp ~ 0.29443                                                                             # Hennig 2013 Table 2: BSV CV% on V2 (independent)

    # Residual error: total proportional 20.4 % (Hennig 2013 Table 2 "Final
    # model"). The 8.4 % "within-sample" component is an L2 (duplicate-sample)
    # term specific to the one centre that recorded duplicate measurements
    # and is not portable to a generic simulation, so only the total
    # proportional RUV is carried forward here.
    propSd <- 0.204; label("Proportional residual error SD (fraction)")                          # Hennig 2013 Table 2: Prop RUV = 20.4 %
  })
  model({
    # Piecewise-linear age effect on CL with breakpoint at 18 years.
    # f_age = 1 + theta_AGE_<=18 * min(0, AGE - 18)
    #           + theta_AGE_>18  * max(0, AGE - 18)
    # AGE = 2  -> f_age = 1 + (-0.021)*(-16)         = 1.336
    # AGE = 18 -> f_age = 1
    # AGE = 70 -> f_age = 1 + (-0.010)*( 52)         = 0.480
    f_age <- 1 + e_age_le18 * min(0, AGE - 18) + e_age_gt18 * max(0, AGE - 18)

    # Renal-function factor on CL (Hennig 2013 Eq. 5).
    f_scr <- (CREAT_REF / CREAT)^e_scr_cl

    # Individual PK parameters (Hennig 2013 Eqs. 6-9). Female is the
    # reference for CL and V1; male values are recovered via the
    # exp(e_male_*) multiplier on (1 - SEXF).
    cl <- exp(lcl + (1 - SEXF) * e_male_cl + etalcl) * (FFM / 70)^e_ffm_cl_q * f_age * f_scr
    vc <- exp(lvc + (1 - SEXF) * e_male_vc + etalvc) * (FFM / 70)
    q  <- exp(lq  + etalq)                           * (FFM / 70)^e_ffm_cl_q
    vp <- exp(lvp + etalvp)                          * (FFM / 70)

    # Two-compartment IV micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg, V in L -> central / vc has units mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
