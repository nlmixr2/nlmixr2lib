Bellon_2020_fibrinogen <- function() {
  description <- "One-compartment population PK model with non-centred estimated allometric body-weight scaling on CL and V for a triple-secured human plasma fibrinogen concentrate (CLOTTAFACT / FibCLOT; plasma fibrinogen activity by Clauss assay) in children, adolescents and adults with congenital afibrinogenaemia (Bellon 2020)"
  reference <- "Bellon A, Fuseau E, Roumanie O, Stevens W, Henriet C, Dahmane A, Lamazure J, Barthez-Toullec M, Golly D, Bridey F. Population pharmacokinetics of a triple-secured fibrinogen concentrate administered to afibrinogenaemic patients: Observed age- and body weight-related differences and consequences for dose adjustment in children. Br J Clin Pharmacol. 2020;86(2):329-337. doi:10.1111/bcp.14147"
  vignette <- "Bellon_2020_fibrinogen"
  units <- list(time = "h", dosing = "g", concentration = "g/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Actual body weight on the day of the PK infusion, entered WITHOUT centring (Bellon 2020 section 3.2: 'actual patient body weight on the day of PK assessment was used in the equations of CL without centring'). CL = 0.00288 * WT^0.556 and V = 0.0960 * WT^0.808 (Table 2; equations 4-5), so the typical values of lcl / lvc are for a 1 kg patient. Observed range 10.8-93.5 kg (Table 1).",
      source_name = "WT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened as a linear covariate in the stepwise forward search (Bellon 2020 section 3.1; supplementary Table S1) but not retained (delta OFV < 3.84)."
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened as a categorical covariate in the stepwise forward search (Bellon 2020 section 3.1; supplementary Table S1) but not retained (delta OFV < 3.84)."
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "fibrinogen (functional activity, Clauss assay)",
      units = "g",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 31L,
    n_studies = 3L,
    age_range = "1.5-48.7 years (17 months to 48 years)",
    age_median = "Study 1: 31.9 years; Study 2: 21.9 years; Study 3: 7.1 years",
    weight_range = "10.8-93.5 kg",
    weight_median = "Study 1: 72.4 kg; Study 2: 64.0 kg; Study 3: 20.0 kg",
    sex_female_pct = round(100 * 13 / 31, 1),
    race_ethnicity = "Not reported",
    disease_state = "Congenital afibrinogenaemia, nonbleeding state, no fibrinogen substitution for >= 2 weeks before the PK infusion; no hepatic or renal impairment, no pregnant women",
    dose_range = "Single intravenous dose of 0.06 g/kg fibrinogen concentrate (15 g/L solution) at a maximum rate of 4 mL/min; median infusion length 1.0 h (range 13 min to 2 h)",
    regions = "Not reported (multicentre studies sponsored by LFB, France)",
    weight_groups = "<40 kg: n = 12 (all Study 3, aged <= 12 years); >= 40 kg: n = 19 (Studies 1 and 2)",
    n_observations = "177 post-infusion fibrinogen activity concentrations (16 first-BLQ values imputed as LLOQ/2)",
    notes = "Pooled clinical-pharmacology parts of three open-label multicentre studies (Bellon 2020 Table 1): Study 1 (phase I/II, adults 18-65 y, n = 5, 8 post-infusion samples to 14 days), Study 2 (phase II/III, >= 40 kg, n = 14, 8 samples to 14 days), Study 3 (FGTW-1004, NCT02094430; phase II/III, <= 12 y, n = 12, 3 samples at 1 h, 3 d and 5 d). Clauss-assay LLOQ 0.09 g/L (Study 1 laboratory) or 0.12 / 0.3 g/L (second laboratory, by dilution). Fibrinogen antigen (nephelometry) was fitted with the same structural model, but the antigen parameter estimates are not reported, so only the activity model is encoded here."
  )

  ini({
    # Structural parameters: Bellon 2020 Table 2 and equations 4-5. Body
    # weight enters UNCENTRED, so these are the typical values for WT = 1 kg.
    lcl <- log(0.00288); label("Clearance for a 1 kg patient (L/h)")                          # Table 2 Theta1 = 0.00288 (CL (L/h) = Theta1 * WT^Theta3); eq. 4
    lvc <- log(0.0960);  label("Volume of distribution for a 1 kg patient (L)")               # Table 2 Theta2 = 0.0960 (V (L) = Theta2 * WT^Theta4); eq. 5

    # Estimated allometric exponents (Bellon 2020 Table 2; Discussion paragraph 1).
    e_wt_cl <- 0.556; label("Allometric exponent of body weight on CL (unitless)")            # Table 2 Theta3 = 0.556 (RSE 9.7%)
    e_wt_vc <- 0.808; label("Allometric exponent of body weight on V (unitless)")             # Table 2 Theta4 = 0.808 (RSE 5.7%)

    # Inter-individual variability: exponential (log-normal) model, variances
    # as printed in Table 2 (the printed CV column is sqrt(omega^2): sqrt(0.0387)
    # = 19.7%, sqrt(0.0250) = 15.8%). No CL-V covariance was estimated
    # (section 3.2).
    etalcl ~ 0.0387  # Table 2 Omega1^2 = 0.0387 (CV 19.7%)
    etalvc ~ 0.0250  # Table 2 Omega2^2 = 0.0250 (CV 15.8%)

    # Residual error: single proportional error shared by all Clauss assay
    # calibration curves (sections 2.3.1 and 3.2). Table 2 prints the variance
    # sigma^2 = 0.0120 with CV 11.0%; propSd = sqrt(0.0120) = 0.1095.
    propSd <- 0.1095; label("Proportional residual error (fraction)")                         # Table 2 sigma^2 = 0.0120 (CV 11.0%)
  })
  model({
    cl <- exp(lcl + etalcl) * WT^e_wt_cl
    vc <- exp(lvc + etalvc) * WT^e_wt_vc

    kel <- cl / vc

    # Intravenous infusion into the single (intravascular) compartment.
    # Patients are afibrinogenaemic, so no endogenous baseline is modelled.
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
