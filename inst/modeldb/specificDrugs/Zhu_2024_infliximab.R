Zhu_2024_infliximab <- function() {
  description <- "One-compartment IV population PK model of infliximab in Chinese patients with Crohn's disease, with fat-free mass on clearance and volume and serum albumin and antibody-to-infliximab titre on clearance (Zhu 2024)"
  reference <- "Zhu K, Ding X, Xue L, Liu L, Wang Y, Li Y, Xi Q, Pang X, Chen W, Miao L. Optimising infliximab induction dosing to achieve clinical remission in Chinese patients with Crohn's disease. Front Pharmacol. 2024;15:1430120. doi:10.3389/fphar.2024.1430120"
  vignette <- "Zhu_2024_infliximab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "infliximab", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    FFM = list(
      description = "Fat-free mass computed from body weight, height and sex by the Janmahasatian et al. (2005) formula",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Zhu 2024 Methods ('popPK modelling'): 'Fat free mass (FFM) incorporated sex, height and bodyweight (Janmahasatian et al., 2005)'. Janmahasatian: FFM = 9270 * WT / (6680 + 216 * BMI) for men and 9270 * WT / (8780 + 244 * BMI) for women, BMI = WT / HT^2 (HT in metres). Normalised to 56.1 kg in Eqs. 3-4: linear on V, exponent 0.75 on CL. Weight was screened as both a baseline and a time-varying covariate; the paper does not say whether the final FFM column is baseline or time-varying.",
      source_name = "FFM"
    ),
    ALB = list(
      description = "Serum albumin concentration (time-varying)",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL normalised to 40 g/L, exponent -0.48 (Zhu 2024 Eq. 4; Supplementary Table S5). Screened as a time-varying covariate (Methods, 'popPK modelling').",
      source_name = "ALB"
    ),
    ADA_TITER = list(
      description = "Antibody-to-infliximab (ATI) titre as a reciprocal dilution at the matched sample time; 0 = ATI absent.",
      units = "(reciprocal dilution; 0 = ATI absent)",
      type = "continuous",
      reference_category = "0 (ATI absent); titre 1:20 (ADA_TITER = 20) shares the reference clearance.",
      notes = "Zhu 2024 reports ATI in four titre classes (Table 1): absent, 1:20, 1:60, >=1:180, a three-fold dilution series. Eq. 4 applies a step function on CL: theta = 1 for ATI <= 1:20 (absent or 1:20), 1.17 for 1:60 and 1.38 for >=1:180. model() derives the class indicators from the continuous column (60 <= ADA_TITER < 180 and ADA_TITER >= 180), so users supply the reciprocal dilution, with 0 for ATI-absent samples. Screened as a time-varying covariate.",
      source_name = "ATI"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on CL and V; not retained (Zhu 2024 Results 'Modelling'; Supplementary Table S3)."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL and V as a standalone covariate and not retained (Supplementary Table S3); sex still enters the model through the FFM formula."
    )
  )

  population <- list(
    n_subjects = 140L,
    n_studies = 1L,
    n_pk_samples = 526L,
    age_range = "15-58 years (median 34, IQR 24-44)",
    age_median = "34 years",
    weight_range = "39-80 kg (mean 59.9, SD 10.8)",
    height = "150-188 cm (mean 158.5, SD 4.95, as printed in Table 1)",
    sex_female_pct = 24.3,
    race_ethnicity = "Chinese (single-centre cohort, Suzhou, China).",
    disease_state = "Crohn's disease treated with intravenous infliximab; 43 of the 140 had moderate-to-severe disease (CDAI >= 220) and formed the week-14 clinical-remission cohort.",
    dose_range = "200, 300 or 400 mg IV (19 / 105 / 16 patients) at weeks 0, 2 and 6, then every 8 weeks.",
    albumin = "27.5-53.1 g/L (median 38.3, IQR 34.4-43.1) before the first infusion",
    ati = "ATI present in 234 of 526 samples (1:20 in 96, 1:60 in 86, >=1:180 in 52).",
    regions = "China (First Affiliated Hospital of Soochow University, October 2020 - December 2022).",
    notes = "Demographics from Zhu 2024 Table 1. 66 of 526 samples (12.55%) were below the 0.50 ug/mL lower limit of quantification and were handled with the M5 method; 96.87% of them were ATI-positive. Treatment duration 2-258 weeks (median 38.6). Height and time-varying erythrocyte sedimentation rate were also screened and not retained. Concentrations are reported in ug/mL, numerically equal to the mg/L used here."
  )

  ini({
    # Typical values for the reference patient: FFM 56.1 kg, albumin 40 g/L,
    # ATI titre <= 1:20 (Zhu 2024 Eqs. 3-4; Supplementary Table S5).
    lcl <- log(0.45); label("Clearance for the reference patient (L/day)") # Supplementary Table S5 'CL (L/d)' = 0.45; Eq. 4
    lvc <- log(9.65); label("Central volume of distribution for the reference patient (L)") # Supplementary Table S5 'V (L)' = 9.65; Eq. 3

    # FFM exponents are not estimated: Eq. 3 prints V linear in FFM / 56.1
    # and Eq. 4 prints the 0.75 exponent on CL; neither appears in Table S5.
    e_ffm_cl <- fixed(0.75); label("Power exponent of FFM on CL (unitless)") # Eq. 4 '(FFM/56.1)^0.75'
    e_ffm_vc <- fixed(1); label("Power exponent of FFM on V (unitless)") # Eq. 3 'V = 9.65 x (FFM/56.1)'

    e_alb_cl <- -0.48; label("Power exponent of serum albumin on CL (unitless)") # Supplementary Table S5 'Albumin on CL' = -0.48; Eq. 4 '(ALB/40)^-0.48'

    # Step-function ATI effect on CL (Eq. 2 categorical form; Eq. 4 theta).
    e_ada60_cl <- 1.17; label("Multiplicative factor on CL for ATI titre 1:60 (unitless)") # Supplementary Table S5 'ATI (1:60) on CL' = 1.17; Eq. 4
    e_ada180_cl <- 1.38; label("Multiplicative factor on CL for ATI titre >= 1:180 (unitless)") # Supplementary Table S5 'ATI (>=1:180) on CL' = 1.38; Eq. 4

    # Exponential IIV on CL only; Table S5 'IIV_CL(%)' = 25.1 read as a CV:
    # omega^2 = log(0.251^2 + 1) = 0.06110.
    etalcl ~ 0.06110 # Supplementary Table S5 'IIV_CL(%)' = 25.1

    # Combined proportional + additive residual error (Methods).
    propSd <- 0.266; label("Proportional residual error (fraction)") # Supplementary Table S5 'Proportional error (%)' = 26.6
    addSd <- 0.15; label("Additive residual error (mg/L)") # Supplementary Table S5 'Additive error, ug/mL' = 0.15
  })
  model({
    # ATI titre class indicators (three-fold dilution series 20 / 60 / 180).
    ada60 <- (ADA_TITER >= 60) * (ADA_TITER < 180)
    ada180 <- (ADA_TITER >= 180)
    ada_cl <- e_ada60_cl^ada60 * e_ada180_cl^ada180

    cl <- exp(lcl + etalcl) * (FFM / 56.1)^e_ffm_cl * (ALB / 40)^e_alb_cl * ada_cl
    vc <- exp(lvc) * (FFM / 56.1)^e_ffm_vc

    kel <- cl / vc
    d/dt(central) <- -kel * central

    # Dose in mg, volume in L -> mg/L (= ug/mL, the paper's unit).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
