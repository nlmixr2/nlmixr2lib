Yin_2021_trastuzumabDeruxtecan <- function() {
  description <- "Two-compartment population PK model for intact trastuzumab deruxtecan (T-DXd, DS-8201, anti-HER2 antibody-drug conjugate) with linear elimination and covariate effects of body weight, albumin, baseline tumor size, sex, and Japan-country indicator in patients with HER2-positive breast cancer or other HER2-expressing solid tumors (Yin 2021)"
  reference <- "Yin O, Iwata H, Lin C-C, Tamura K, Watanabe J, Wada R, Kastrissios H, Garimella T, Lee C, Zhang L, Shahidi J, Fujisaki Y, LaCreta F. Population Pharmacokinetics of Trastuzumab Deruxtecan in Patients With HER2-Positive Breast Cancer and Other Solid Tumors. Clin Pharmacol Ther. 2021;109(5):1314-1325. doi:10.1002/cpt.2096"
  vignette <- "Yin_2021_trastuzumabDeruxtecan"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effect on CL_intact (exponent 0.370) and V1_intact (exponent 0.489); reference 57.8 kg per Yin 2021 final model equations (Results: 'Development of intact trastuzumab deruxtecan model').",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Baseline serum albumin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value (SI units; not g/dL). Power effect on CL_intact (exponent -0.533); reference 40 g/L per Yin 2021 final model equations.",
      source_name        = "ALB"
    ),
    TUMSZ = list(
      description        = "Baseline tumor size (sum of diameters of target lesions per RECIST)",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effect on CL_intact (exponent 0.0710); reference 57 mm per Yin 2021 final model equations.",
      source_name        = "Tumor size"
    ),
    SEXF = list(
      description        = "Biological sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male) in the canonical column. The paper's own reference category is female (see notes).",
      notes              = "Yin 2021 final model encodes sex as a male-indicator (1 = male, 0 = female) with female as the reference category (Results: 'CL_intact = ... x (1.174, if male)' and 'V1,intact = ... x (1.197, if male)'). To store under the canonical SEXF (1 = female, 0 = male), the effect is applied in model() as (1 + e_male_cl * (1 - SEXF)) and (1 + e_male_v1 * (1 - SEXF)) so SEXF = 1 yields factor 1 (paper female reference) and SEXF = 0 yields the paper's male multiplier (1.174 on CL, 1.197 on V1).",
      source_name        = "Sex"
    ),
    REGION_JAPAN = list(
      description        = "Japan enrollment-country indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Japan country)",
      notes              = "Yin 2021 retained Country (Japan vs non-Japan) over Race because the two were highly confounded (correlation -0.81) and Country was more significant on all PK parameters. Multiplicative fractional effect of 0.903 on CL_intact and 0.738 on V2_intact for REGION_JAPAN = 1, applied as (1 + e_japan_cl * REGION_JAPAN) and (1 + e_japan_v2 * REGION_JAPAN).",
      source_name        = "Country"
    )
  )

  population <- list(
    n_subjects     = 639L,
    n_studies      = 5L,
    phase_mix      = "4 phase I studies (J101, J102, A103, A104) and 1 phase II study (DESTINY-Breast01)",
    n_observations_intact = 11434L,
    age_range      = "median 57 years (per Figure 5 caption: 'A typical patient is defined as a 57-year-old female...')",
    weight_range   = "5th-95th percentile not reported in main text; median 57.8 kg (per final model reference values and Figure 4 caption)",
    weight_median  = "57.8 kg",
    sex_female_pct = NA_real_,
    race_ethnicity = "Multi-regional cohort with sufficient Japanese enrollment for a Japan-vs-non-Japan effect to be supported. Race and Country were highly confounded (correlation -0.81); Race was dropped in favour of Country in the final model.",
    disease_state  = "HER2-positive (or HER2-expressing) advanced or metastatic solid tumors. 512 / 639 (80.1%) patients had unresectable / metastatic breast cancer; 445 / 639 (69.6%) had HER2-positive unresectable / metastatic breast cancer. Other tumor types include HER2-expressing gastric, lung, colorectal, salivary-gland and other solid tumors.",
    dose_range     = "0.8 - 8.0 mg/kg IV every 3 weeks (q3w), 21-day treatment cycles",
    regions        = "Multi-regional (US, Japan, EU). Japanese phase I studies J101 and J102 contributed substantially to the Japan-country covariate signal.",
    reference_subject = "Female, non-Japan country, 57.8 kg body weight, 40 g/L albumin, 57 mm baseline tumor size (Figure 4 caption: 'A typical patient is defined as a female from a non-Japan country with body weight 57.8 kg, albumin 40 g/L, and baseline tumor size 57 mm.').",
    notes          = "NONMEM 7.3, FOCE-INTER. Released-drug PK (1-compartment with time-varying release-rate constant and DAR/molar-mass-adjusted input from intact T-DXd) is described in the paper but is NOT implemented in this model file: it requires (a) a dosing-cycle index derived from time and (b) a DAR x molar-mass conversion factor that the paper does not give numerically. See the validation vignette 'Assumptions and deviations' section. Final analysis data set: 11,434 intact T-DXd serum concentrations from 639 patients."
  )

  ini({
    # Structural parameters - typical values for the Yin 2021 reference patient
    # (female, non-Japan, WT 57.8 kg, ALB 40 g/L, TUMSZ 57 mm) per Table 1.
    # NONMEM Table 1 reports CL_intact and Q_intact directly in L/day; V1 and V2
    # are reported in L. Time unit kept as day to match the paper's tabulation.
    lcl <- log(0.421); label("Clearance of intact T-DXd CL_intact at reference covariates (L/day)") # Yin 2021 Table 1: CL_intact = 0.421 L/day
    lvc <- log(2.77);  label("Central volume of distribution V1_intact at reference (L)")           # Yin 2021 Table 1: V1_intact = 2.77 L
    lq  <- log(0.199); label("Distributional clearance Q_intact (L/day)")                           # Yin 2021 Table 1: Q_intact = 0.199 L/day
    lvp <- log(5.16);  label("Peripheral volume of distribution V2_intact at reference (L)")        # Yin 2021 Table 1: V2_intact = 5.16 L

    # Covariate effects (Yin 2021 Results 'Development of intact trastuzumab
    # deruxtecan model' final-model equations and Table 1 'Final intact T-DXd
    # model' rows). Continuous covariates enter as (cov / ref)^exponent. The
    # categorical covariates (Japan country, male sex) enter as fractional
    # multipliers (1 + theta * indicator) so that theta_male_cl = 0.174 yields
    # the paper-reported 1.174 multiplier on CL when male, and theta_japan_cl
    # = -0.0970 yields the paper-reported 0.903 multiplier when Japanese.
    e_wt_cl     <-  0.370;   label("Power exponent of WT on CL_intact (unitless)")                                       # Yin 2021 Table 1: Body weight on CL_intact = 0.370
    e_alb_cl    <- -0.533;   label("Power exponent of ALB on CL_intact (unitless)")                                      # Yin 2021 Table 1: Albumin on CL_intact = -0.533
    e_tumsz_cl  <-  0.0710;  label("Power exponent of TUMSZ on CL_intact (unitless)")                                    # Yin 2021 Table 1: Tumor size on CL_intact = 0.0710
    e_japan_cl  <- -0.0970;  label("Fractional multiplicative effect of REGION_JAPAN on CL_intact (unitless)")            # Yin 2021 Table 1: Country (Japan) on CL_intact = -0.0970 -> multiplier 0.903
    e_male_cl   <-  0.174;   label("Fractional multiplicative effect of male sex on CL_intact (unitless; on (1 - SEXF))") # Yin 2021 Table 1: Sex (male) on CL_intact = 0.174 -> multiplier 1.174
    e_wt_v1     <-  0.489;   label("Power exponent of WT on V1_intact (unitless)")                                       # Yin 2021 Table 1: Body weight on V1_intact = 0.489
    e_male_v1   <-  0.197;   label("Fractional multiplicative effect of male sex on V1_intact (unitless; on (1 - SEXF))") # Yin 2021 Table 1: Sex (male) on V1_intact = 0.197 -> multiplier 1.197
    e_japan_v2  <- -0.262;   label("Fractional multiplicative effect of REGION_JAPAN on V2_intact (unitless)")            # Yin 2021 Table 1: Country (Japan) on V2_intact = -0.262 -> multiplier 0.738

    # Inter-individual variability. Yin 2021 Table 1 'Between-patient
    # variability' block reports variances directly: Var(CL_intact) = 0.0630,
    # Var(V1_intact) = 0.0250, Var(Q_intact) = 0.0900, Var(V2_intact) = 0.430,
    # Cov(CL_intact, V1_intact) = 0.0210. The reported 'Magnitude (%CV)'
    # column matches sqrt(omega^2) (the small-CV approximation), not the exact
    # log-normal %CV = sqrt(exp(omega^2) - 1).
    etalcl + etalvc ~ c(0.0630,
                        0.0210, 0.0250)  # Yin 2021 Table 1: Var(CL_intact)=0.0630, Cov(CL,V1)=0.0210, Var(V1_intact)=0.0250
    etalq  ~ 0.0900                       # Yin 2021 Table 1: Var(Q_intact) = 0.0900
    etalvp ~ 0.430                        # Yin 2021 Table 1: Var(V2_intact) = 0.430

    # Residual error (Yin 2021 Table 1 'Residual variability' block).
    # Combined proportional + additive on intact T-DXd serum concentration.
    # Additive SD reported in ng/mL; converted to ug/mL (the model's
    # concentration unit) via 1 ng/mL = 0.001 ug/mL -> 1181 ng/mL = 1.181 ug/mL.
    propSd <- 0.163; label("Proportional residual error SD (fraction)")                # Yin 2021 Table 1: proportional residual error SD = 0.163
    addSd  <- 1.181; label("Additive residual error SD on Cc (ug/mL)")                  # Yin 2021 Table 1: additive residual error SD = 1,181 ng/mL = 1.181 ug/mL
  })

  model({
    # Derived sex term. Yin 2021 encodes sex as a male-indicator with female as
    # the reference category, so (1 - SEXF) reproduces the paper's male = 1
    # column while keeping SEXF (1 = female) as the canonical storage form.
    sex_male <- 1 - SEXF

    # Individual PK parameters with covariate adjustments
    # (Yin 2021 final-model equations, Results section). Reference values:
    # WT 57.8 kg, ALB 40 g/L, TUMSZ 57 mm, female (SEXF = 1), non-Japan
    # (REGION_JAPAN = 0). Continuous covariates use power-of-ratio scaling;
    # categorical covariates use fractional multipliers (1 + theta * indicator).
    cl <- exp(lcl + etalcl) *
      (WT    / 57.8)^e_wt_cl *
      (ALB   / 40  )^e_alb_cl *
      (TUMSZ / 57  )^e_tumsz_cl *
      (1 + e_japan_cl * REGION_JAPAN) *
      (1 + e_male_cl  * sex_male)

    vc <- exp(lvc + etalvc) *
      (WT / 57.8)^e_wt_v1 *
      (1 + e_male_v1 * sex_male)

    q  <- exp(lq + etalq)
    vp <- exp(lvp + etalvp) *
      (1 + e_japan_v2 * REGION_JAPAN)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L -> central / vc has units mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
