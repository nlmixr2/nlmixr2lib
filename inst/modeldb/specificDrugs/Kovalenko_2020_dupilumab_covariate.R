Kovalenko_2020_dupilumab_covariate <- function() {
  description <- "Dupilumab primary covariate population PK model from Kovalenko 2020 (Model 4): 2-compartment with parallel linear + Michaelis-Menten elimination and a 3-transit-compartment SC absorption chain; fit to Phase 3 atopic-dermatitis data with body weight + albumin on Vc and BMI + EASI + race (White) on the linear elimination rate."
  reference <- "Kovalenko P, Davis JD, Li M, et al. Base and Covariate Population Pharmacokinetic Analyses of Dupilumab Using Phase 3 Data. Clinical Pharmacology in Drug Development. 2020;9(6):756-767. doi:10.1002/cpdd.780"
  vignette <- "Kovalenko_2020_dupilumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")
  # Model 4 (primary COVARIATE model) from Kovalenko 2020 Table 1, Table 3,
  # and Supplementary Tables S3 and S4 (Model 4 / "Parameterized Using ke"
  # columns).  Phase 3 atopic-dermatitis cohort only (R668-AD-1334 / SOLO 1,
  # R668-AD-1416 / SOLO 2, R668-AD-1224 / CHRONOS).  The structural model is
  # identical to the base model in `Kovalenko_2020_dupilumab_base.R`
  # (Model 3); the difference is that albumin enters Vc multiplicatively and
  # BMI, EASI, and race (White) enter ke multiplicatively.  All fixed
  # parameters (kcp, kpc, ka, MTT, Vm, Km, F) are carried forward from the
  # rich-data fits (Models 1-2) as in Model 3.
  #
  # The paper notes that while ADA, albumin, race, BMI, and EASI score were
  # statistically significant, ONLY body weight had a notable effect on Vc
  # explaining interindividual variability.  This file implements every
  # significant covariate as published; users who want a simpler model can
  # fall back to `Kovalenko_2020_dupilumab_base.R`.
  #
  # Covariate equation forms:
  #
  #   * Continuous covariates use the power form
  #     Y(cov) = Y * (cov / cov_ref)^theta, as stated in the paper's Methods.
  #     Reference values (cov_ref) are the cohort medians; the paper states
  #     they are taken at "median or another selected level of covariate"
  #     but does not tabulate the chosen reference values.  The values
  #     adopted here are the closest published medians from sibling
  #     analyses of the same Phase 3 dupilumab AD trials:
  #
  #       WT_ref  = 75   kg     (Kovalenko 2016 doi:10.1002/psp4.12136 Eq. 1)
  #       ALB_ref = 44   g/L    (Zhang 2021 doi:10.1002/psp4.12667 pooled
  #                              dupilumab cohort median; closest published
  #                              dupilumab popPK median in g/L)
  #       BMI_ref = 26   kg/m^2 (LIBERTY AD SOLO 1/2 mean baseline BMI
  #                              ~26.4 kg/m^2 per Simpson 2016 NEJM
  #                              doi:10.1056/NEJMoa1610020 Table 1)
  #       EASI_ref = 32         (LIBERTY AD SOLO 1/2 mean baseline EASI
  #                              ~33 per Simpson 2016 Table 1; rounded to
  #                              32 to match CHRONOS median; Blauvelt 2017
  #                              doi:10.1016/S0140-6736(17)31191-1 Table 1)
  #
  #     All reference values are documented in the vignette's Errata so a
  #     reader can re-scale if newer cohort medians become available.
  #
  #   * Dichotomous race (White) uses the multiplicative form
  #     Y(RACE_WHITE) = Y * (1 + theta * RACE_WHITE), consistent with the
  #     Zhang 2021 dupilumab popPK ADA implementation (same drug, same
  #     authorship group), giving Y(non-White) = Y and
  #     Y(White) = Y * (1 + theta).  For ke at theta = -0.123 the result is
  #     ke_White = 0.877 * ke_non-White (i.e., White subjects clear
  #     dupilumab via the linear route ~12% more slowly than non-White
  #     subjects).
  #
  # IIV: Supplementary Table S3 reports SD(ln Vc) = 0.206,
  # SD(ln ke) = 0.293, Corr(ln(ke), ln(Vc)) = -0.450.  Variances and
  # covariance for the correlated (etalvc, etalkel) block in ini() are
  # computed inline so the reader can trace the SDs to the source table.
  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on central volume Vc with reference 75 kg.  Reference weight is not restated in Kovalenko 2020; the 75 kg value is inherited from Kovalenko 2016 (doi:10.1002/psp4.12136, Eq. 1) and matches the sibling Kovalenko_2020_dupilumab.R (Model 1) and Kovalenko_2020_dupilumab_base.R (Model 3) files.",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Serum albumin (baseline)",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on central volume Vc with reference 44 g/L.  Reference value is not restated in Kovalenko 2020; the 44 g/L value is the published median from Zhang 2021 (doi:10.1002/psp4.12667) which analysed an overlapping dupilumab popPK cohort.  ALB in g/L is the SI canonical; multiply by 0.1 to convert US-convention g/dL values to g/L on data ingestion.",
      source_name        = "ALB"
    ),
    BMI = list(
      description        = "Body mass index (baseline)",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on the linear elimination rate ke with reference 26 kg/m^2 (approx. median baseline BMI of the LIBERTY AD SOLO 1/2 pooled cohorts per Simpson 2016 doi:10.1056/NEJMoa1610020 Table 1).  Reference value is not restated in Kovalenko 2020.",
      source_name        = "BMI"
    ),
    SCORE_EASI = list(
      description        = "Eczema Area and Severity Index (baseline)",
      units              = "(score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on the linear elimination rate ke with reference 32 (approx. median baseline EASI of the Phase 3 dupilumab AD cohort; Simpson 2016 reports mean baseline EASI ~33 in SOLO 1/2 Table 1; CHRONOS Blauvelt 2017 Table 1 reports a similar median).  Reference value is not restated in Kovalenko 2020.  Canonical column name is SCORE_EASI; the prior EASI alias is preserved for compatibility.",
      source_name        = "EASI"
    ),
    RACE_WHITE = list(
      description        = "White race indicator (1 = White, 0 = non-White)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-White; pools Black/African American, Asian, American Indian/Alaska Native, Native Hawaiian/Pacific Islander, Other / Not reported in the Phase 3 dupilumab AD trials per Simpson 2016)",
      notes              = "Multiplicative effect on the linear elimination rate ke as ke * (1 + e_white_kel * RACE_WHITE); the paper's covariate label is 'ke ~ race (White)' with the theta = -0.123 indicating White subjects clear dupilumab via the linear route ~12% more slowly than the non-White reference.  Implementation form follows the Zhang 2021 dupilumab ADA precedent (same drug, same authorship group).",
      source_name        = "RACE"
    )
  )

  population <- list(
    n_subjects     = "Phase 3 cohort only.  The article reports the pooled total across all 16 studies as N = 2115 on study and 2041 on active treatment, with 18,243 of 20,809 samples included.  Phase 3 SOLO 1 (R668-AD-1334) N = 447, SOLO 2 (R668-AD-1416) N = 472, CHRONOS (R668-AD-1224) N = 424 per Supplementary Table S1 (PK analysis set).",
    n_studies      = 3L,
    age_range      = "Adults with moderate-to-severe atopic dermatitis (detailed age breakdown not reported in the main text).",
    age_median     = "Not reported in the main text.",
    weight_range   = "Not reported in the main text.",
    weight_median  = "Not reported in the main text; reference weight 75 kg is the Kovalenko 2016 inheritance.",
    sex_female_pct = "Not reported in the main text.",
    race_ethnicity = "Race was a tested covariate and retained in the final model as a multiplicative effect of White vs non-White on the linear elimination rate.  Phase 3 AD trials are predominantly White, with Black, Asian, and Other categories represented.",
    disease_state  = "Adults with moderate-to-severe atopic dermatitis (SOLO 1, SOLO 2 monotherapy; CHRONOS concomitant topical corticosteroids).",
    dose_range     = "600 mg SC loading dose on day 1 followed by 300 mg SC qw or q2w for 15 weeks (SOLO 1, SOLO 2) or 51 weeks (CHRONOS).",
    regions        = "Multi-regional Phase 3 programme; see Supplementary Table S1 for per-study geographic coverage.",
    notes          = "Primary covariate model (Model 4) for regulatory submissions.  The structural parameters kcp, kpc, ka, MTT, Vm, Km, and F are FIXED to values obtained from Models 1 and 2 fits on rich Phase 1/2 data (per the stepwise modelling strategy described in the paper's Methods).  Only weight had a notable effect on Vc explaining interindividual variability; albumin / BMI / EASI / race retained as statistically significant (P < 1e-8 for all) but with smaller effect sizes."
  )

  ini({
    # Estimated structural parameters on Phase 3 data (Supplementary Table S3, Model 4 column)
    lvc  <- log(2.74);   label("central volume at the reference covariate values (L)")  # Supp. Table S3 Model 4: Vc = 2.74 L (SE 0.021)
    lkel <- log(0.0477); label("linear elimination rate at the reference covariate values (1/d)") # Supp. Table S3 Model 4: ke = 0.0477 1/d (SE 0.00078)

    # Fixed structural parameters (same as Model 3; carried forward from Models 1 and 2)
    lkcp    <- fixed(log(0.211));  label("central-to-peripheral rate kcp (1/d)")        # Supp. Table S3 Model 4 / Table 1: kcp = 0.211 (fixed)
    lkpc    <- fixed(log(0.310));  label("peripheral-to-central rate kpc (1/d)")        # Supp. Table S3 Model 4 / Table 1: kpc = 0.310 (fixed)
    lka     <- fixed(log(0.306));  label("absorption rate ka (1/d)")                     # Supp. Table S3 Model 4 / Table 1: ka  = 0.306 (fixed)
    lmtt    <- fixed(log(0.105));  label("mean transit time MTT (d)")                    # Supp. Table S3 Model 4 / Table 1: MTT = 0.105 (fixed)
    lvmax   <- fixed(log(1.07));   label("maximum target-mediated rate of elimination Vmax (mg/L/d)") # Supp. Table S3 Model 4 / Table 1: Vm = 1.07 (fixed)
    Km      <- fixed(0.01);        label("Michaelis-Menten constant Km (mg/L)")          # Supp. Table S3 Model 4 / Table 1: Km = 0.01 (fixed; carried over from Kovalenko 2016)
    lfdepot <- fixed(log(0.642));  label("subcutaneous bioavailability F (fraction)")    # Supp. Table S3 Model 4 / Table 1: F  = 0.642 (fixed)

    # Covariate effects on Vc (power form, multiplicative)
    e_wt_vc  <-  0.817; label("Power exponent of WT/75 on Vc (unitless)")               # Supp. Table S3 Model 4 / Table 3: Vc ~ weight  =  0.817 (SE 0.031)
    e_alb_vc <- -0.653; label("Power exponent of ALB/44 on Vc (unitless)")              # Supp. Table S3 Model 4 / Table 3: Vc ~ albumin = -0.653 (SE 0.072)

    # Covariate effects on ke (power form for continuous, additive multiplier for binary)
    e_bmi_kel        <-  0.368; label("Power exponent of BMI/26 on ke (unitless)")     # Supp. Table S3 Model 4 / Table 3: ke ~ BMI         =  0.368 (SE 0.053)
    e_score_easi_kel <-  0.143; label("Power exponent of SCORE_EASI/32 on ke (unitless)") # Supp. Table S3 Model 4 / Table 3: ke ~ EASI        =  0.143 (SE 0.021)
    e_race_white_kel <- -0.123; label("Fractional effect of RACE_WHITE on ke as ke * (1 + e_race_white_kel * RACE_WHITE) (unitless)") # Supp. Table S3 Model 4 / Table 3: ke ~ race (White) = -0.123 (SE 0.018)

    # Inter-individual variability: correlated (ln Vc, ln ke) block.
    # Supp. Table S3 Model 4 reports omegas as SDs on the log scale:
    #   SD(ln Vc) = 0.206, SD(ln ke) = 0.293
    #   Corr(ln(ke), ln(Vc)) = -0.450
    # Variances: 0.206^2 = 0.042436, 0.293^2 = 0.085849
    # Covariance: -0.450 * 0.206 * 0.293 = -0.0271611
    etalvc + etalkel ~ c(0.042436,
                         -0.0271611, 0.085849)

    # Residual error (combined proportional + additive on Cc in mg/L)
    propSd <- 0.125; label("Proportional residual error (fraction)")                    # Supp. Table S3 Model 4: sigma_prop (CV%) = 12.5 (SE 0.18)
    addSd  <- 6.06;  label("Additive residual error (mg/L)")                            # Supp. Table S3 Model 4: sigma_add (mg/L) = 6.06 (SE 0.23)
  })
  model({
    # Individual structural parameters with the Kovalenko 2020 Model 4 covariate equations.
    # Typical reference patient: 75 kg body weight, 44 g/L serum albumin, 26 kg/m^2 BMI,
    # SCORE_EASI = 32, non-White race (RACE_WHITE = 0).
    vc   <- exp(lvc  + etalvc)  *
            (WT  / 75)^e_wt_vc *
            (ALB / 44)^e_alb_vc
    kel  <- exp(lkel + etalkel) *
            (BMI        / 26)^e_bmi_kel *
            (SCORE_EASI / 32)^e_score_easi_kel *
            (1 + e_race_white_kel * RACE_WHITE)
    kcp  <- exp(lkcp)
    kpc  <- exp(lkpc)
    ka   <- exp(lka)
    MTT  <- exp(lmtt)
    vmax <- exp(lvmax)

    # Transit-chain rate for 3 transit compartments
    ktr <- (3 + 1) / MTT

    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * (depot    - transit1)
    d/dt(transit2)    <-  ktr * (transit1 - transit2)
    d/dt(transit3)    <-  ktr * transit2 - ka * transit3
    d/dt(central)     <-  ka * transit3 -
                          kel * central -
                          kcp * central + kpc * peripheral1 -
                          central * (vmax / (Km + central / vc))
    d/dt(peripheral1) <-  kcp * central - kpc * peripheral1

    # Subcutaneous bioavailability on the depot; IV doses bypass the depot via the event record
    f(depot) <- exp(lfdepot)

    # mg dosed / L central volume gives mg/L concentration; no unit conversion required
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
