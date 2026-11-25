MorentinGutierrez_2019_AZD7687_human <- function() {
  description <- "Population PK/PD model for AZD7687, a diacylglycerol acyltransferase 1 (DGAT1) inhibitor, and postprandial plasma triglyceride (TAG) during an oral lipid tolerance test (standardized mixed meal) in healthy and overweight/obese men (Morentin Gutierrez 2019). Two-compartment first-order-absorption PK of the FREE (unbound) plasma AZD7687 concentration with body weight on CL, V1 and V2; meal fat enters a lagged first-order TAG depot whose input to plasma is inhibited by AZD7687 (1 - C/(C + IC50)), on top of an endogenous TAG turnover (baseline R0 with BMI covariate, loss rate kout)."
  reference <- "Morentin Gutierrez P, Yates J, Nilsson C, Birtles S. Evolving data analysis of an Oral Lipid Tolerance Test toward the standard for the Oral Glucose Tolerance Test: Cross species modeling effects of AZD7687 on plasma triacylglycerol. Pharmacol Res Perspect. 2019;7(2):e00465. doi:10.1002/prp2.465"
  vignette <- "MorentinGutierrez_2019_AZD7687"

  # `tg` (plasma triglyceride) and `depot_tg` (meal fat in the gut) are not yet
  # canonical compartments; `tg` is declared paper-specific in the same way by
  # Ahn_2023_vupanorsen.R and Ousey_2026_plozasiran.R.
  paper_specific_compartments <- c("depot_tg", "tg")

  # The meal is given into `depot_tg`, not into the drug depot.
  dosing <- c("depot", "depot_tg")

  units <- list(
    time = "h",
    dosing = "umol (AZD7687 into depot; convert mg with the molecular weight 367.4 g/mol) and mg (meal triglyceride into depot_tg; grams of fat = kcal from fat / 9)",
    concentration = "umol/L (free, unbound AZD7687, Cc) and mg/dL (plasma triglyceride, tg)"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power covariate on CL, V1 and V2, centred at the analysis mean of 85.3 kg (Supplementary section S1, eqs 1-3).",
      source_name = "BW"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power covariate on the baseline plasma TAG R0, centred at the analysis mean of 27.5 kg/m^2 (Supplementary section S2, eq 4).",
      source_name = "BMI"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "AZD7687",
      units = "umol",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "AZD7687",
      units = "umol",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "AZD7687",
      units = "umol",
      specimen = "plasma",
      verified = TRUE
    ),
    depot_tg = list(
      analyte = "triglycerides",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    tg = list(
      analyte = "triglycerides",
      units = "mg/dL",
      specimen = "serum",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 142,
    n_studies = 2,
    age_range = "not reported in this paper (adult men; see Denison 2013/2014)",
    weight_range = "not reported (analysis mean 85.3 kg)",
    weight_median = "mean 85.3 kg",
    bmi_mean = "27.5 kg/m^2",
    sex_female_pct = 0,
    disease_state = "Healthy men (SAD) and overweight or obese otherwise healthy men (MAD)",
    dose_range = "SAD: single oral doses 1-60 mg (10 cohorts, 6 active + 2 placebo each); MAD: 1, 2.5 (once or twice daily), 5, 10 or 20 mg/day orally for 1 week",
    meal_challenge = "Standardized mixed meal of 1100 kcal with 45% or 60% of energy from fat; SAD meal 4 h after dosing, MAD meal 1 h after dosing",
    regions = "USA (SAD) and UK (MAD)",
    notes = "80 healthy men in the single-ascending-dose study and 62 overweight/obese men (42 active, 20 placebo) in the multiple-ascending-dose study (Methods, Data: Human). PK modelled as free concentrations, fu = 0.055."
  )

  ini({
    # AZD7687 PK (free concentration) -- Supplementary Table 1 (also main Table 1)
    lka <- log(5.19); label("Absorption rate constant ka (1/h)") # Suppl Table 1 tvka = 5.19 1/hr
    lcl <- log(28.1); label("Apparent clearance of free AZD7687 CL (L/h)") # Suppl Table 1 tvCL = 28.1 L/hr
    lvc <- log(320.2); label("Apparent central volume of free AZD7687 V1 (L)") # Suppl Table 1 tvV1 = 320.2 L
    lq <- log(50.8); label("Apparent distributional clearance CLD (L/h)") # Suppl Table 1 tvCLD = 50.8 L/hr
    lvp <- log(136.8); label("Apparent peripheral volume V2 (L)") # Suppl Table 1 tvV2 = 136.8 L

    e_wt_cl <- 1.05; label("Power exponent of body weight on CL (unitless)") # Suppl Table 1 dCLdWeight = 1.05
    e_wt_vc <- 1.18; label("Power exponent of body weight on V1 (unitless)") # Suppl Table 1 dV1dWeight = 1.18
    e_wt_vp <- 1.89; label("Power exponent of body weight on V2 (unitless)") # Suppl Table 1 dV2dWeight = 1.89

    # Plasma TAG PK/PD -- Supplementary Table 4 (also main Table 1)
    ltlag_tg <- log(1.63); label("Lag time of meal TAG absorption TlagTAG (h)") # Suppl Table 4 tvTlagTAG = 1.63 hr
    lka_tg <- log(0.74); label("First-order absorption rate of meal TAG kaTAG (1/h)") # Suppl Table 4 tvKaTAG = 0.74 1/hr
    lvc_tg <- log(250.4); label("Volume of distribution of exogenous TAG VTAG (dL)") # Suppl Table 4 tvVTAG = 250.4 dL
    lrbase_tg <- log(97.8); label("Baseline plasma TAG R0 (mg/dL)") # Suppl Table 4 tvR0 = 97.8 mg/dL
    lkout_tg <- log(0.70); label("Fractional turnover rate of plasma TAG kout (1/h)") # Suppl Table 4 tvkout = 0.70 1/hr
    lic50_tg <- log(0.0078); label("Free AZD7687 concentration giving 50% inhibition of meal TAG input IC50 (umol/L)") # Suppl Table 4 tvIC50 = 0.0078 umol/L
    e_bmi_rbase_tg <- 1.28; label("Power exponent of BMI on R0 (unitless)") # Suppl Table 4 dR0dBMI = 1.28

    # IIV: log-normal; reported as %CV, converted with omega^2 = log(CV^2 + 1)
    etalcl + etalvc ~ c(
      0.046429,
      0.024436, 0.036945
    ) # Suppl Table 1 IIV CL 21.8%, IIV V1 19.4%, Cl-V1 correlation 0.59
    etalka ~ 0.32885 # Suppl Table 1 IIV ka 62.4%
    etalvp ~ 0.19959 # Suppl Table 1 IIV V2 47.0%
    etalka_tg ~ 0.040383 # Suppl Table 4 IIV KaTAG 20.3%
    etalvc_tg ~ 0.28306 # Suppl Table 4 IIV VTAG 57.2%
    etalrbase_tg ~ 0.13688 # Suppl Table 4 IIV R0 38.3%
    etalkout_tg ~ 0.35242 # Suppl Table 4 IIV kout 65.0%
    etalic50_tg ~ 1.7902 # Suppl Table 4 IIV IC50 223.4%

    # Residual error
    addSd <- 0.0005; label("Additive residual error, free AZD7687 (umol/L)") # Suppl Table 1 Conc Residual Additive = 0.0005 uM
    propSd <- 0.159; label("Proportional residual error, free AZD7687 (fraction)") # Suppl Table 1 Conc Residual Multiplicative = 0.159
    propSd_tg <- 0.192; label("Proportional residual error, plasma TAG (fraction)") # Suppl Table 4 TAG Residual Multiplicative = 0.192
  })

  model({
    # AZD7687 PK; body weight centred at the analysis mean (Suppl S1 eqs 1-3)
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 85.3)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 85.3)^e_wt_vc
    q <- exp(lq)
    vp <- exp(lvp + etalvp) * (WT / 85.3)^e_wt_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc

    # Plasma TAG; BMI centred at the analysis mean (Suppl S2 eq 4)
    tlag_tg <- exp(ltlag_tg)
    ka_tg <- exp(lka_tg + etalka_tg)
    vc_tg <- exp(lvc_tg + etalvc_tg)
    rbase_tg <- exp(lrbase_tg + etalrbase_tg) * (BMI / 27.5)^e_bmi_rbase_tg
    kout_tg <- exp(lkout_tg + etalkout_tg)
    ic50_tg <- exp(lic50_tg + etalic50_tg)

    # Eq 1: fraction of meal TAG input not blocked by AZD7687
    inh_tg <- 1 - Cc / (Cc + ic50_tg)

    # Eqs 2-4; the lag acts on the meal depot
    d/dt(depot_tg) <- -ka_tg * depot_tg
    alag(depot_tg) <- tlag_tg
    d/dt(tg) <- rbase_tg * kout_tg + ka_tg * depot_tg / vc_tg * inh_tg - kout_tg * tg
    tg(0) <- rbase_tg

    Cc ~ add(addSd) + prop(propSd)
    tg ~ prop(propSd_tg)
  })
}
