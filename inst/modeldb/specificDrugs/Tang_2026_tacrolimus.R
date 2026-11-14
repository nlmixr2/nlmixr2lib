Tang_2026_tacrolimus <- function() {
  description <- "Sequential population PK/PD model of once-daily modified-release oral tacrolimus (Graceptor) in adult liver transplant recipients during the first post-transplant year (Tang 2026). PK: one-compartment model with first-order absorption (ka fixed at 0.19 1/h) fitted to whole-blood troughs; apparent clearance CL/F rises with post-operative day up to day 14, (min(POD, 14)/14)^1.68, and falls with hematocrit, (HCT/30)^-0.59; correlated exponential IIV on CL/F and V/F; additive residual error. PD: direct linear model with no intercept relating the cumulative tacrolimus AUC over a 3-month interval (cAUC3mon, mg*h/L) to the change in the estimated glomerular filtration rate over the same interval; the slope (-0.59 mL/min/1.73 m^2 per mg*h/L at a pre-transplant eGFR >= 90) is modified by the pre-transplant eGFR category, with additive IIV on the slope and additive residual error. The cumulative AUC follows the paper's definition, the running sum of daily dose / individual CL/F (state cauc, for once-daily dosing); the integral of the concentration is carried separately in auc_central."
  reference <- paste(
    "Tang N, Maeda S, Okuda M.",
    "Population Pharmacokinetic/Pharmacodynamic Model of the Impact of",
    "Cumulative Tacrolimus Exposure on Kidney Function After Liver Transplant.",
    "Ther Drug Monit. 2026;48(4):473-480.",
    "doi:10.1097/FTD.0000000000001418.",
    "Parameter values are the final-model estimates in Table 2; the covariate",
    "equations are the Table 2 footnote.",
    sep = " "
  )
  vignette <- "Tang_2026_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE),
    auc_central = list(analyte = "tacrolimus", units = "mg*h/L", specimen = "not applicable", verified = TRUE),
    cauc = list(analyte = "tacrolimus", units = "mg*h/L", specimen = "not applicable", verified = TRUE)
  )

  # `cauc` is the paper's cumulative AUC, the running sum of the daily
  # AUC0-24h computed as daily dose / individual CL/F (Tang 2026 Methods
  # 'Population PK/PD Model'). It is not the integral of the concentration
  # (that is `auc_central`), so it is not the canonical accumulator.
  paper_specific_compartments <- c("cauc")

  covariateData <- list(
    POD = list(
      description = "Post-operative day: days elapsed since liver transplantation",
      units = "days",
      type = "continuous",
      reference_category = NULL,
      notes = "TIME-VARYING within subject. Enters CL/F through POD14X = min(POD, 14) as the power term (POD14X / 14)^1.68 (Tang 2026 Table 2 footnote: 'if POD <= 14, POD14X = POD, else POD14X = 14'), so CL/F rises steeply over the first two post-operative weeks and is constant afterwards. Cutoffs of 7, 14 and 30 days were tested and 14 gave the best fit (Methods 'Population PK Model'). The term is zero at POD = 0, which makes CL/F zero; supply POD >= 1.",
      source_name = "POD"
    ),
    HCT = list(
      description = "Hematocrit, packed red blood cell volume fraction times 100",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying in the source data (routine laboratory values over the first post-transplant year). Enters CL/F as (HCT / 30)^-0.59 (Tang 2026 Table 2 footnote). Cohort median 28.7 % (range 13.8-50, Table 1). Higher hematocrit lowers whole-blood apparent clearance because tacrolimus partitions into erythrocytes.",
      source_name = "HCT"
    ),
    RENALIMP_MILD = list(
      description = "Pre-transplant baseline eGFR 60 to < 90 mL/min/1.73 m^2 (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (pre-transplant baseline eGFR >= 90 mL/min/1.73 m^2 when RENALIMP_MOD and RENALIMP_SEV are also 0)",
      notes = "Time-fixed. One of three mutually exclusive indicators decomposing the paper's four-level baseline-eGFR category (>= 90 reference, 60 to < 90, 30 to < 60, < 30 mL/min/1.73 m^2; Results 'Population PK/PD Model'). The baseline eGFR is the mean of the last three eGFR values before transplantation, computed with the Japanese equation 194 * age^-0.287 * SCr^-1.094 (* 0.739 if female) (Methods and Table 1 footnote). Multiplies the eGFR slope by (1 - 0.68).",
      source_name = "eGFR baseline category"
    ),
    RENALIMP_MOD = list(
      description = "Pre-transplant baseline eGFR 30 to < 60 mL/min/1.73 m^2 (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (pre-transplant baseline eGFR >= 90 mL/min/1.73 m^2 when RENALIMP_MILD and RENALIMP_SEV are also 0)",
      notes = "Time-fixed; see RENALIMP_MILD for the category definition. Multiplies the eGFR slope by (1 - 1.14), which reverses its sign: in this stratum the predicted eGFR rises slightly with tacrolimus exposure (Results; the authors read it as a data artefact, not a protective effect).",
      source_name = "eGFR baseline category"
    ),
    RENALIMP_SEV = list(
      description = "Pre-transplant baseline eGFR < 30 mL/min/1.73 m^2 (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (pre-transplant baseline eGFR >= 90 mL/min/1.73 m^2 when RENALIMP_MILD and RENALIMP_MOD are also 0)",
      notes = "Time-fixed; see RENALIMP_MILD for the category definition. This stratum includes any end-stage values below 15 mL/min/1.73 m^2; the paper does not split them out. Multiplies the eGFR slope by (1 - 0.92).",
      source_name = "eGFR baseline category"
    )
  )

  # Screened in the stepwise covariate search but not retained (Tang 2026
  # Methods 'Population PK Model' and 'Population PK/PD Model'; no point
  # estimates are published for them).
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL/F and V/F and on the eGFR slope; not retained. Cohort median 59 kg (range 36.8-118, Table 1)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL/F and V/F and on the eGFR slope; not retained. Cohort median 56 years (range 22-69, Table 1)."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened on CL/F and V/F; not retained. 54.5 % female (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 44L,
    n_studies = 1L,
    n_observations = "3360 tacrolimus troughs (PK); 3504 eGFR measurements (PK/PD)",
    age_range = "22-69 years (median 56)",
    weight_range = "36.8-118 kg (median 59)",
    sex_female_pct = 54.5,
    race_ethnicity = "Japanese (single centre in Osaka, Japan); race not tabulated.",
    disease_state = "Adult liver transplant recipients (84.1 % living donor, 15.9 % deceased donor; 22.7 % ABO-incompatible) during the first post-transplant year. Observed eGFR 5.2-202 mL/min/1.73 m^2 (median 57.6) and hematocrit 13.8-50 % (median 28.7) over the study period (Table 1).",
    dose_range = "Once-daily modified-release oral tacrolimus (Graceptor, equivalent to Advagraf), physician-adjusted to whole-blood trough targets of 8-12 ng/mL in weeks 1-2, 6-10 ng/mL from week 3 to month 3 and 4-8 ng/mL thereafter; observed doses 0.1-32 mg (median 4 mg).",
    co_medication = "All received ganciclovir and sulfamethoxazole-trimethoprim, 93.2 % furosemide; regimens tacrolimus alone or with steroid and/or mycophenolate (Table 1).",
    regions = "University of Osaka Hospital, Japan; transplants 2010-2020.",
    notes = "Retrospective therapeutic-drug-monitoring data, troughs only, whole blood by chemiluminescent immunoassay (LLOQ 1.0 ng/mL; values below were set to 1.0). Of 89 recipients screened, 44 were analysed. Fitted in Phoenix NLME 8.3.4 with FOCE-ELS, sequentially: the PK/PD step used the empirical Bayes CL/F estimates of the PK model."
  )

  ini({
    # PK (Tang 2026 Table 2, 'Population PK model' block). Typical CL/F is at
    # POD >= 14 and HCT = 30 %.
    lka <- fixed(log(0.19)); label("Absorption rate constant ka (1/h)") # Table 2: ka = 0.19 h-1 (fixed); Methods: fixed from Lu 2019 (reference 8)
    lcl <- log(23.46); label("Apparent clearance CL/F at POD >= 14 days and HCT 30 percent (L/h)") # Table 2: CL/F = 23.46 L/h (RSE 9.21%)
    lvc <- log(1607.81); label("Apparent volume of distribution V/F (L)") # Table 2: V/F = 1607.81 L (RSE 12.36%)
    e_pod_cl <- 1.68; label("Power exponent of (min(POD, 14) / 14) on CL/F (unitless)") # Table 2: theta POD14X on CL = 1.68 (RSE 17.18%)
    e_hct_cl <- -0.59; label("Power exponent of (HCT / 30) on CL/F (unitless)") # Table 2: theta HCT on CL = -0.59 (RSE 42.15%)

    # IIV. Table 2 footnote: 'omega, standard deviation of the IIV', so the
    # omega(CL,V) row is a correlation (as a covariance of two SDs it would
    # imply r = 0.47 / (0.62 * 0.71) = 1.07). Variances are SD^2 and the
    # covariance is 0.47 * 0.62 * 0.71 = 0.206894.
    etalcl + etalvc ~ c(0.3844, 0.206894, 0.5041) # Table 2: omega CL/F = 0.62, omega V/F = 0.71 (SDs), omega CL,V = 0.47 (correlation)

    addSd <- 3.04; label("Additive residual error, tacrolimus whole-blood concentration (ng/mL)") # Table 2 PK block: sigma = 3.04 (RSE 7.45%)

    # PD (Table 2, 'Population PK/PD model' block). The table prints theta
    # Slope as 0.59; the Results text ('The population estimate of the slope
    # was -0.59') and the footnote ('Slope = -0.59 when eGFR baseline >= 90;
    # -0.19, when 60 <= eGFR baseline < 90; 0.08, when 30 <= eGFR baseline
    # < 60; -0.05, when eGFR baseline < 30') give the sign. The category
    # effects enter as Slope * (1 + theta_cat), which reproduces those values.
    slope <- -0.59; label("Change in eGFR per unit cAUC at pre-transplant eGFR >= 90 (mL/min/1.73 m^2 per mg*h/L)") # Table 2: theta Slope = 0.59 (RSE 23.73%), sign from Results and Table 2 footnote
    e_renalimp_mild_slope <- -0.68; label("Fractional change in slope for pre-transplant eGFR 60 to <90 (unitless)") # Table 2: 60 <= eGFR baseline < 90 = -0.68 (RSE 22.39%)
    e_renalimp_mod_slope <- -1.14; label("Fractional change in slope for pre-transplant eGFR 30 to <60 (unitless)") # Table 2: 30 <= eGFR baseline < 60 = -1.14 (RSE 17.26%)
    e_renalimp_sev_slope <- -0.92; label("Fractional change in slope for pre-transplant eGFR <30 (unitless)") # Table 2: eGFR baseline < 30 = -0.92 (RSE 13.25%)

    etaslope ~ 0.04 # Table 2: omega Slope = 0.20 (SD, additive IIV); variance 0.20^2

    addSd_egfrcfb <- 15.41; label("Additive residual error, change in eGFR (mL/min/1.73 m^2)") # Table 2 PK/PD block: sigma = 15.41 (RSE 12.06%)
  })
  model({
    # Post-operative day capped at 14 (Table 2 footnote, POD14X)
    pod14x <- POD
    if (POD > 14) pod14x <- 14

    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (pod14x / 14)^e_pod_cl * (HCT / 30)^e_hct_cl
    vc <- exp(lvc + etalvc)
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    # Integral of the whole-blood concentration in mg*h/L (central / vc is
    # mg/L).
    d/dt(auc_central) <- central / vc

    # The paper's cumulative AUC (Methods 'Population PK/PD Model'): each
    # day adds daily dose / CL/F of that day, in mg*h/L. Spread over the
    # 24 h after each once-daily dose as dose / (24 * CL/F); nothing is added
    # more than 24 h after the last dose. podo(depot) is NA before the first
    # dose.
    dose_last <- podo(depot)
    if (is.na(dose_last)) dose_last <- 0
    if (tad(depot) >= 24) dose_last <- 0
    d/dt(cauc) <- dose_last / (24 * cl)

    # mg/L to ng/mL
    Cc <- 1000 * central / vc

    # Additive IIV on the covariate-adjusted slope (Results: 'the direct
    # linear model with an additive IIV model').
    slope_i <- slope *
      (1 +
        e_renalimp_mild_slope * RENALIMP_MILD +
        e_renalimp_mod_slope * RENALIMP_MOD +
        e_renalimp_sev_slope * RENALIMP_SEV) +
      etaslope
    # Change in eGFR from the start of the cAUC accumulation (time 0)
    egfrcfb <- slope_i * cauc

    Cc ~ add(addSd)
    egfrcfb ~ add(addSd_egfrcfb)
  })
}
