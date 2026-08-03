Steichert_2025_enalapril_enalaprilat_pediatric <- function() {
  description <- "Simultaneous parent + active-metabolite population PK model for oral enalapril (ODMT) and enalaprilat in ACEi-naive children with heart failure (Steichert 2025, LENA studies). Combined one-compartment model for enalapril (first-order absorption with a lag) coupled with a one-compartment model for enalaprilat via a fixed fraction metabolised fm = 0.7. Allometric scaling (fixed exponents 0.75 on CL, 1 on V) referenced to 5 kg body weight. Covariate effects retained in the final model: age and serum creatinine on the apparent clearance of enalaprilat, and modified Ross score on the apparent volume of distribution of enalaprilat."
  reference   <- "Steichert M, Cawello W, Laeer S; LENA Consortium. Population Pharmacokinetic Analysis of Enalapril and Enalaprilat in Newly Treated Children with Heart Failure: Implications for Safe Dosing of Enalapril (LENA Studies). Clin Pharmacokinet. 2025;64(7):1103-1118. doi:10.1007/s40262-025-01520-5"
  vignette    <- "Steichert_2025_enalapril_enalaprilat_pediatric"
  units       <- list(time = "h", dosing = "ug", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight, used as the allometric size descriptor on the apparent clearances and apparent volumes of distribution of both enalapril and enalaprilat. Reference weight is the weighted median (5 kg) computed by Perl-speaks-NONMEM across the LENA ACEi-naive cohort.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying; the Steichert 2025 dataset carries body weight per PK-sample record and applies last-observation-carried-forward when a scheduled weight is missing (Section 2.2.1). Range 2.52-11.3 kg in the analysed cohort (Table 1).",
      source_name        = "Weight"
    ),
    AGE = list(
      description        = "Subject age, entered as a power-form covariate on the apparent clearance of enalaprilat: `(AGE / 0.34)^0.311`. Reference age is the weighted median (0.34 years) computed by Perl-speaks-NONMEM across the LENA ACEi-naive cohort.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying (Section 2.2.2). Range 0.07-2.09 years (25 days to 2.1 years) in the analysed cohort (Table 1). The paper tested a power form initially and confirmed no exponential-form improvement in the backward step.",
      source_name        = "Age"
    ),
    CREAT = list(
      description        = "Serum creatinine concentration, entered as an exponential-form covariate on the apparent clearance of enalaprilat: `exp(-0.0141 * (CREAT - 23.37))`. Reference is the weighted median 23.37 umol/L computed by Perl-speaks-NONMEM across the LENA ACEi-naive cohort.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Range 12-68 umol/L in the analysed cohort (Table 1). The paper tested a power form first and found the exponential form superior (Section 2.2.2).",
      source_name        = "Serum creatinine"
    ),
    SCORE_ROSS = list(
      description        = "Modified Ross score, a paediatric heart-failure severity instrument, entered as an exponential-form covariate on the apparent volume of distribution of enalaprilat: `exp(-0.15 * (SCORE_ROSS - 4))`. Reference is the weighted median score 4 across the LENA ACEi-naive cohort.",
      units              = "(score; integer 0-12)",
      type               = "count",
      reference_category = NULL,
      notes              = "Time-varying, determined by the investigator (Section 2.2.2). Range 0-9 observed in the analysed cohort (Table 1). Exponential form required because the score can legitimately be zero, so a power form is undefined at that boundary. Modified Reithmann / Ross / Connolly form as used across the LENA paediatric heart-failure programme.",
      source_name        = "Ross score"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex, coded 0 = male, 1 = female in Steichert 2025 (their equation TV = theta1 * (1 + theta2 * sex) with female = 0, male = 1, which the register maps by inverting the sign of the effect on the SEXF = 1 female indicator).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a covariate on both CL and V of enalapril and enalaprilat in the stepwise search (Section 2.2.2). Not retained in the final model."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 34,
    n_studies      = 2,
    age_range      = "0.07-2.09 years (25 days to 2.1 years)",
    age_median     = "0.3 years (unweighted) / 0.34 years (weighted median used as covariate reference)",
    weight_range   = "2.52-11.3 kg",
    weight_median  = "4.47 kg (unweighted) / 5 kg (weighted median used as allometric reference)",
    sex_female_pct = 52.9,
    disease_state  = "ACEi-naive children with heart failure due to congenital heart disease (91.2%) or dilated cardiomyopathy (8.8%). Modified Ross score 0-9 (median 5). Serum creatinine 12-68 umol/L (median 27).",
    dose_range     = "0.03-0.29 mg/kg/day enalapril maleate administered as orodispersible mini-tablets (Aqumeldi 0.25 mg or 1 mg). Titration doses were age- and weight-dependent per the LENA dosing regimen (Supplementary Table 1).",
    regions        = "Austria, Germany, Hungary, the Netherlands, Serbia (LENA multicentre phase II/III PK-bridging studies EudraCT 2015-002335-17 for DCM and EudraCT 2015-002396-18 for CHD).",
    n_observations = "173 quantifiable enalapril + 268 quantifiable enalaprilat serum concentrations. Samples below LOQ were excluded (M1 method); LOQ 0.195 ug/L for enalapril and 0.180 ug/L for enalaprilat.",
    notes          = "Demographics and dosing summary from Steichert 2025 Table 1 and Sections 2.1.1 / 3.1. Sampling regimen: full PK profile (predose + 1, 2, 4, 6, 12 h) at first dose or steady state, plus single trough samples at titration / dose-confirmation / study-control / end-of-study visits."
  )

  ini({
    # Structural PK parameters from Steichert 2025 Table 2 ("Final model" column).
    # Reference-subject values are given for a 5 kg body weight (the population
    # weighted median). ka and fm are fixed (see below); lcl / lvc / lcl_enaat /
    # lvc_enaat / ltlag are estimated point estimates.
    lka       <- fixed(log(0.6))
    label("Absorption rate constant of enalapril (1/h), from a prior LENA analysis")            # Table 2: ka = 0.6 1/h, Fixed (paper text Section 2.2.1)
    lcl       <- log(4.61)
    label("Apparent clearance of enalapril CL_ENA/F at 5 kg (L/h)")                                    # Table 2: CL_ENA/F = 4.61 L/h at BW = 5 kg (RSE 12.6%)
    lvc       <- log(4.98)
    label("Apparent volume of distribution of enalapril Vd_ENA/F at 5 kg (L)")                         # Table 2: Vd_ENA/F = 4.98 L at BW = 5 kg (RSE 18.1%)
    lcl_enaat <- log(1.55)
    label("Apparent clearance of enalaprilat CL_ENAAT/F at 5 kg (L/h)")                                # Table 2: CL_ENAAT/F = 1.55 L/h at BW = 5 kg (RSE 7.1%)
    lvc_enaat <- log(34.1)
    label("Apparent volume of distribution of enalaprilat Vd_ENAAT/F at 5 kg (L)")                     # Table 2: Vd_ENAAT/F = 34.1 L at BW = 5 kg (RSE 15.7%)
    ltlag     <- log(0.515)
    label("Absorption lag time on enalapril depot (h)")                                                # Table 2: tlag = 0.515 h (RSE 2.8%)

    # Fixed structural constants (Steichert 2025 Section 2.2.1).
    lfm       <- fixed(log(0.7))
    label("Fraction of enalapril metabolised to enalaprilat")                                   # Section 2.2.1: fm = 0.7, fixed from literature [ref 9]

    # Fixed allometric exponents (Steichert 2025 Section 2.2.1; Table 2).
    e_wt_cl       <- fixed(0.75)
    label("Allometric exponent on CL_ENA/F")                                                    # Table 2: Weight on CL_ENA/F = 0.75, Fixed
    e_wt_cl_enaat <- fixed(0.75)
    label("Allometric exponent on CL_ENAAT/F")                                                  # Table 2: Weight on CL_ENAAT/F = 0.75, Fixed
    e_wt_vc       <- fixed(1)
    label("Allometric exponent on Vd_ENA/F")                                                    # Table 2: Weight on Vd_ENA/F = 1, Fixed
    e_wt_vc_enaat <- fixed(1)
    label("Allometric exponent on Vd_ENAAT/F")                                                  # Table 2: Weight on Vd_ENAAT/F = 1, Fixed

    # Estimated covariate effects retained in the final model.
    # Age enters as a power exponent (CL_ENAAT/F multiplied by (AGE/0.34)^0.311).
    # Serum creatinine enters as an exponential coefficient (CL_ENAAT/F multiplied
    # by exp(e_creat_cl_enaat * (CREAT - 23.37))). SCORE_ROSS enters as an
    # exponential coefficient on Vd_ENAAT/F (Vd_ENAAT/F multiplied by
    # exp(e_score_ross_vc_enaat * (SCORE_ROSS - 4))).
    e_age_cl_enaat        <- 0.311
    label("Power exponent for AGE on CL_ENAAT/F (referenced to 0.34 y)")                                # Table 2: Age on CL_ENAAT/F = 0.311 (RSE 29.3%); Section 3.2.2 equation
    e_creat_cl_enaat      <- -0.0141
    label("Exponential coefficient for CREAT on CL_ENAAT/F (1/(umol/L); referenced to 23.37 umol/L)")   # Table 2: Serum creatinine on CL_ENAAT/F = -0.0141 (RSE 33.3%); Section 3.2.2 equation
    e_score_ross_vc_enaat <- -0.15
    label("Exponential coefficient for SCORE_ROSS on Vd_ENAAT/F (1/point; referenced to 4)")            # Table 2: Ross score on Vd_ENAAT/F = -0.15 (RSE 26.3%); Section 3.2.2 equation

    # IIV variances (log-normal). Steichert 2025 Table 2 reports IIV as CV%
    # calculated as `sqrt(omega^2) * 100%` (Table 2 footnote a), so omega^2 =
    # (CV/100)^2 rather than the log-normal exact form `log(1 + CV^2)`.
    #   65.3% CV -> omega^2 = 0.653^2 = 0.4264
    #   80.2% CV -> omega^2 = 0.802^2 = 0.6432
    #   37.7% CV -> omega^2 = 0.377^2 = 0.1421
    #   90.4% CV -> omega^2 = 0.904^2 = 0.8172
    etalcl       ~ 0.4264
    etalvc       ~ 0.6432
    etalcl_enaat ~ 0.1421
    etalvc_enaat ~ 0.8172

    # Residual variability. Steichert 2025 Table 2 footnotes b/c give:
    #   proportional CV% = sqrt(sigma^2) * 100, so propSd = CV/100 as the
    #   proportional-SD parameter (linear-space, applied via `prop(propSd)`);
    #   additive SD = sqrt(sigma^2) in the reported concentration units.
    propSd          <- 0.535
    label("Proportional residual error on enalapril (fraction)")                                        # Table 2: Residual variability enalapril, Proportional error 53.5% CV
    addSd           <- 1.34
    label("Additive residual error on enalapril (ug/L)")                                                # Table 2: Residual variability enalapril, Additive error 1.34 ug/L
    propSd_enaat <- 0.395
    label("Proportional residual error on enalaprilat (fraction)")                                      # Table 2: Residual variability enalaprilat, Proportional error 39.5% CV
  })

  model({
    # Reference values (population weighted medians from Perl-speaks-NONMEM;
    # Steichert 2025 Section 3.2.2 and Table 2 footnote).
    wt_ref    <- 5      # kg, weighted median of the LENA ACEi-naive cohort
    age_ref   <- 0.34   # years, weighted median
    creat_ref <- 23.37  # umol/L, weighted median
    ross_ref  <- 4      # score, weighted median

    # Fixed structural constants.
    fm    <- exp(lfm)     # 0.7, fraction metabolised (Section 2.2.1)

    # Absorption.
    ka    <- exp(lka)     # 0.6 1/h, fixed
    tlag  <- exp(ltlag)

    # Individual PK parameters. Steichert 2025 Section 3.2.2 equations:
    #   CL_ENA/F     = 4.61 * (WT/5)^0.75 * exp(eta1)
    #   Vd_ENA/F     = 4.98 * (WT/5)^1    * exp(eta3)
    #   CL_ENAAT/F   = 1.55 * (WT/5)^0.75 * (AGE/0.34)^0.311
    #                       * exp(-0.0141 * (CREAT - 23.37)) * exp(eta2)
    #   Vd_ENAAT/F   = 34.1 * (WT/5)^1
    #                       * exp(-0.15 * (SCORE_ROSS - 4)) * exp(eta4)
    cl       <- exp(lcl       + etalcl)       * (WT / wt_ref)^e_wt_cl
    vc       <- exp(lvc       + etalvc)       * (WT / wt_ref)^e_wt_vc
    cl_enaat <- exp(lcl_enaat + etalcl_enaat) * (WT / wt_ref)^e_wt_cl_enaat *
                (AGE / age_ref)^e_age_cl_enaat *
                exp(e_creat_cl_enaat * (CREAT - creat_ref))
    vc_enaat <- exp(lvc_enaat + etalvc_enaat) * (WT / wt_ref)^e_wt_vc_enaat *
                exp(e_score_ross_vc_enaat * (SCORE_ROSS - ross_ref))

    # Elimination rate constants (Steichert 2025 Fig 2 caption transfer-rate
    # equations: k20 = (1-fm) * CL_ENA / Vd_ENA to urine; k23 = fm * CL_ENA /
    # Vd_ENA to enaat; k30 = CL_ENAAT / Vd_ENAAT to urine). The two parent-side
    # rates sum to CL_ENA/Vd_ENA, so the parent-central ODE loss is written as
    # a single first-order elimination term and the metabolite formation flux
    # is a scaled fraction (fm) of the same parent-loss flux.
    kel_ena   <- cl       / vc
    kel_enaat <- cl_enaat / vc_enaat

    # ODEs (Steichert 2025 Fig 2). depot -> central [enalapril] -> central_enaat
    # [enalaprilat]; the (1-fm) fraction of parent flux is excreted directly and
    # is not carried in an explicit urine compartment (Fig 2 shows k20 to urine
    # but no urine sampling was fit).
    d/dt(depot)         <- -ka * depot
    d/dt(central)       <-  ka * depot - kel_ena * central
    d/dt(central_enaat) <-  fm * kel_ena * central - kel_enaat * central_enaat

    # Absorption lag on the depot (Steichert 2025 Section 2.2.1 and Table 2).
    alag(depot) <- tlag

    # Observations and residual error. Cc = enalapril, Cc_enaat = enalaprilat.
    Cc       <- central       / vc
    Cc_enaat <- central_enaat / vc_enaat

    Cc       ~ add(addSd) + prop(propSd)
    Cc_enaat ~ prop(propSd_enaat)
  })
}
