Bilal_2024_rifampicin_sexwt <- function() {
  description <- paste(
    "One-compartment population PK model for a single 600 mg oral dose of",
    "rifampicin in healthy Caucasian adult volunteers (Bilal 2024), the",
    "alternative sex + body weight covariate model that the authors ranked",
    "second-best. Absorption is a zero-order input of duration Tk0",
    "(0.470 h) into the central compartment after a lag time (0.340 h);",
    "elimination is Michaelis-Menten (Vmax/F 157 mg/h in females, Km",
    "20.2 mg/L) from a central volume V/F of 33.2 L in females. Male sex",
    "multiplies V/F by exp(0.190) and Vmax/F by exp(0.420); body weight",
    "scales V/F linearly (exponent fixed to 1), normalised to the cohort",
    "mean weight of about 72.6 kg. Between-subject variability is on",
    "Vmax/F only (19.7%); inter-occasion variability over the two crossover",
    "periods is on lag time, Tk0, V/F and Vmax/F, with correlated lag-time",
    "and Tk0 occasion effects. Residual error is the Monolix combined1 form",
    "(additive 0.060 mg/L plus 11.0% proportional, summed as standard",
    "deviations). Dose records into central must carry rate = -2 so the",
    "modelled zero-order duration applies."
  )
  reference <- paste(
    "Bilal M, Ullah S, Jaehde U, Trueck C, Zaremba D, Wachall B, Wargenau M,",
    "Scheidel B, Wiesen MHJ, Gazzaz M, Chen C, Buesker S, Fuhr U, Taubert M,",
    "Dokos C. Assessment of body mass-related covariates for rifampicin",
    "pharmacokinetics in healthy Caucasian volunteers. Eur J Clin Pharmacol.",
    "2024;80(9):1271-1283. doi:10.1007/s00228-024-03697-3.",
    "Parameter estimates from Table 3 (sex + BW covariate model column); the",
    "structural model and random-effect structure from the Monolix code of",
    "the base model in the supplementary material.",
    sep = " "
  )
  vignette <- "Bilal_2024_rifampicin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(
      analyte = "rifampicin",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power model on V/F with the exponent fixed to 1 (Table 3), normalised",
        "to the observation-weighted mean body weight (Methods, Covariate",
        "analysis; Eq. 1). The weighted-mean reference is not printed. Every",
        "subject contributed 38 observations (912 / 24; 19 post-dose samples",
        "in each of two periods), so it equals the arithmetic mean body",
        "weight, which is also not printed. 72.6 kg is the",
        "sex-composition-weighted average of the per-sex Table 1 medians,",
        "(13 * 63.8 + 11 * 82.9) / 24; the overall cohort median is 68 kg",
        "(Results). Median (range) 63.8 (48.5-73.1) kg in females and 82.9",
        "(57.8-91.0) kg in males (Table 1)."
      ),
      source_name = "BW"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (female)",
      notes = paste(
        "Females are the reference category (Methods, Covariate analysis), so",
        "the male effects enter on the indicator (1 - SEXF). The paper calls",
        "the sex model 'a linear model' without printing it. Monolix, which",
        "fitted the model, enters a categorical covariate on a log-normal",
        "parameter linearly on the log scale, log(V) = log(Vpop) + beta *",
        "male, so male V/F is exp(0.190) = 1.21 and male Vmax/F exp(0.420) =",
        "1.52 times the female values. The (1 + beta) reading would give 1.19",
        "and 1.42; see the vignette. 13 of 24 subjects female (Table 1)."
      ),
      source_name = "SEX"
    ),
    OCC = list(
      description = "Crossover period (1 or 2), the occasion for the inter-occasion variability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "The two periods of the bioequivalence crossover were treated as two",
        "occasions (Methods, Population PK analysis); each period is one",
        "single 600 mg dose and its 24-h sampling. Set OCC = 1 for the first",
        "period and OCC = 2 for the second. The Monolix code declares the",
        "variability levels id and id*occ."
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    FFM = list(
      description = "Fat-free mass by the Janmahasatian (2005) formula.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Not a covariate of this model. FFM on V/F and Vmax/F was the best covariate model (dOFV -56.3 to -56.4; Table 4, model 7) and is the sibling Bilal_2024_rifampicin_ffm."
    ),
    BSA = list(
      description = "Body surface area by the Mosteller formula.",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Replacing BW by BSA on V/F (with sex on Vmax/F) lowered the OFV by a further 3.24 points but failed the backward-elimination criterion (Results; Table 4, model 4). Median (range) 1.66 (1.47-1.85) m^2 in females and 2.05 (1.62-2.18) m^2 in males (Table 1)."
    ),
    BMI = list(
      description = "Body mass index.",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Significant on Tlag (dOFV -4.02) and V/F (dOFV -5.11) in forward inclusion but failed backward elimination (Results). Median (range) 23.09 (18.7-26.7) in females and 24.7 (20.1-29.7) in males (Table 1)."
    ),
    HT = list(
      description = "Body height.",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "BH on V/F and Vmax/F lowered the OFV by 42.9 points (Table 4, model 3), less than sex + BW. Reported in metres: median (range) 1.65 (1.55-1.76) m in females and 1.81 (1.63-1.89) m in males (Table 1)."
    ),
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Methods, Covariate analysis) and not retained. Median (range) 37.0 (21.0-58.0) years in females and 43.0 (22.0-64.0) years in males (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 24L,
    n_studies = 1L,
    age_range = "21-64 years (median 39.5)",
    weight_range = "48.5-91.0 kg (median 68)",
    sex_female_pct = 54.2,
    race_ethnicity = "Caucasian (100%)",
    disease_state = "Healthy volunteers, BMI 18.5-30 kg/m^2, non-smokers",
    dose_range = "Single oral 600 mg rifampicin tablet in each of two periods (test tablet and the EREMFAT reference tablet), separated by a washout of at least 9 days",
    regions = "Germany (Cologne)",
    notes = paste(
      "Phase I/IV randomised, two-period crossover, open-label",
      "bioequivalence study (EudraCT 2017-004418-24). 25 enrolled, one",
      "dropped out before the first dose; 24 completed (11 men, 13 women).",
      "Plasma sampled predose and at 0.16, 0.33, 0.5, 0.75, 1, 1.25, 1.5,",
      "1.75, 2, 2.25, 2.5, 3, 3.5, 4, 6, 9, 12, 16 and 24 h. 912",
      "concentrations, 99 (10.8%) below the 0.1 mg/L LLOQ and treated as",
      "interval-censored. The two tablets were bioequivalent by NCA and the",
      "preparation had no significant effect on any parameter, so the data",
      "were pooled. Fit in Monolix 2023R1 (SAEM). Demographics are Table 1."
    )
  )

  ini({
    # Structural parameters: Table 3, 'Sex + BW covariate model' estimates
    # column. V/F and Vmax/F are the female typical values at the reference
    # body weight.
    ltlag <- log(0.340)
    label("Absorption lag time Tlag (h)") # Table 3 sex + BW model Tlag = 0.340 h (RSE 5.49%)
    ld1 <- log(0.470)
    label("Duration of the zero-order absorption Tk0 (h)") # Table 3 sex + BW model Tk0 = 0.470 h (RSE 9.02%)
    lvc <- log(33.2)
    label("Apparent volume of distribution V/F in females at the reference weight of 72.6 kg (L)") # Table 3 sex + BW model V/F = 33.2 L (RSE 3.11%)
    lvmax <- log(157)
    label("Apparent maximum elimination rate Vmax/F in females (mg/h)") # Table 3 sex + BW model Vmax/F = 157 mg/h (RSE 6.70%)
    lkm <- log(20.2)
    label("Michaelis-Menten constant Km (mg/L)") # Table 3 sex + BW model Km = 20.2 mg/L (RSE 4.70%)

    # Covariate effects. The sex coefficients are on the log scale (Monolix
    # categorical covariate on a log-normal parameter).
    e_sexm_vc <- 0.190
    label("Log-scale effect of male sex on V/F (unitless)") # Table 3 'beta sex on V/F' = 0.190 (RSE 24.6%)
    e_wt_vc <- fixed(1)
    label("Power exponent of body weight on V/F (unitless)") # Table 3 'beta BW on V/F' = 1.00 (Fixed)
    e_sexm_vmax <- 0.420
    label("Log-scale effect of male sex on Vmax/F (unitless)") # Table 3 'beta sex on Vmax/F' = 0.420 (RSE 20.3%)

    # Monolix reports random effects as omega, the SD of the log-normal eta;
    # Table 3 prints them x 100 under a 'CV %' heading. Variances are omega^2.
    etalvmax ~ 0.038809 # Table 3 sex + BW model IIV Vmax/F = 19.7; 0.197^2

    # Inter-occasion variability, two occasions (crossover periods). Monolix
    # gives one gamma per parameter, shared by both occasions; occasion 2
    # repeats the occasion-1 value as fixed. Lag-time and Tk0 occasion
    # effects are correlated: cov = 0.420 * 0.388 * 0.657 = 0.10706472.
    etaiov_tlag_1 + etaiov_d1_1 ~ c(
      0.150544,
      0.10706472,
      0.431649
    ) # Table 3 sex + BW model IOV Tlag = 38.8, IOV Tk0 = 65.7, Corr. IOV Tlag and IOV Tk0 = 0.420
    etaiov_tlag_2 + etaiov_d1_2 ~ fixed(c(
      0.150544,
      0.10706472,
      0.431649
    ))
    etaiov_vc_1 ~ 0.023716 # Table 3 sex + BW model IOV V = 15.4; 0.154^2
    etaiov_vc_2 ~ fixed(0.023716)
    etaiov_vmax_1 ~ 0.00868624 # Table 3 sex + BW model IOV Vmax = 9.32; 0.0932^2
    etaiov_vmax_2 ~ fixed(0.00868624)

    # Residual error: Monolix combined1, SD = a + b * f.
    addSd <- 0.060
    label("Additive residual error (mg/L)") # Table 3 sex + BW model additive residual error = 0.060 mg/L (RSE 7.85%)
    propSd <- 0.110
    label("Proportional residual error (fraction)") # Table 3 sex + BW model proportional residual error = 11.0% (RSE 3.33%)
  })
  model({
    # Occasion indicators: crossover period 1 or 2.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_tlag <- oc1 * etaiov_tlag_1 + oc2 * etaiov_tlag_2
    iov_d1 <- oc1 * etaiov_d1_1 + oc2 * etaiov_d1_2
    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2
    iov_vmax <- oc1 * etaiov_vmax_1 + oc2 * etaiov_vmax_2

    # Female reference; Eq. 1 power model for body weight centred on the
    # cohort mean weight (see covariateData).
    tlag <- exp(ltlag + iov_tlag)
    d1 <- exp(ld1 + iov_d1)
    vc <- exp(lvc + e_sexm_vc * (1 - SEXF) + iov_vc) * (WT / 72.6)^e_wt_vc
    vmax <- exp(lvmax + e_sexm_vmax * (1 - SEXF) + etalvmax + iov_vmax)
    km <- exp(lkm)

    # Figure 1: one compartment with Michaelis-Menten elimination
    # Vmax * Cp / (Km + Cp); the Monolix pkmodel(Tlag, Tk0, V, Vm, Km).
    Cc <- central / vc
    d/dt(central) <- -vmax * Cc / (km + Cc)

    # Zero-order input of duration Tk0 starting Tlag after the dose; dose
    # records need rate = -2.
    alag(central) <- tlag
    dur(central) <- d1

    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
