Thoueille_2023_lopinavir <- function() {
  description <- "One-compartment first-order-absorption population PK model for oral ritonavir-boosted lopinavir (LPV/r) in HIV-negative individuals receiving 5-day COVID-19 post-exposure prophylaxis (COPEP study) pooled with people living with HIV followed by routine therapeutic drug monitoring; body weight enters apparent oral clearance as a linear deviation from a 70 kg reference (Thoueille 2023)."
  reference <- "Thoueille P, Delfraysse M, Andre P, Buclin T, Decosterd LA, Fedeli C, Ustero P, Calmy A, Guidi M; Swiss HIV Cohort Study. Population pharmacokinetic analysis of lopinavir in HIV negative individuals exposed to SARS-CoV-2: a COPEP (COronavirus Post-Exposure Prophylaxis) sub-study. BMC Pharmacol Toxicol. 2023;24:47. doi:10.1186/s40360-023-00687-6"
  vignette <- "Thoueille_2023_lopinavir"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "The only covariate retained in the final model. Enters apparent oral clearance as a linear deviation from the 70 kg reference weight: TVCL_LPV = CL_LPV * (1 + theta_BW * (WT - 70) / 70) (Thoueille 2023, 'Final model' equation printed beneath Table 2; theta_BW reference weight of 70 kg stated in the Table 2 footnote). Cohort medians 73 kg (range 47-157) for the COPEP participants and 66 kg (range 40-147) for the PLWH on routine TDM (Thoueille 2023 Table 1). Time-fixed per subject in the source dataset.",
      source_name        = "BW"
    )
  )

  # Covariates screened by Thoueille 2023 but NOT retained in the final
  # model. Documented here for provenance only; none is referenced in
  # model(). Source: Thoueille 2023 Methods, "Models development,
  # evaluation and assessment" ("The following covariates were tested,
  # using linear or allometric functions as deemed appropriate: age, sex,
  # bodyweight, height, body mass index, smoking status and type of
  # population (i.e., COPEP vs PLWH). The latter was investigated on all
  # the PK parameters.").
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at sampling.",
      units       = "years",
      type        = "continuous",
      notes       = "Screened by forward inclusion / backward deletion on the PK parameters; not retained in the final model. Cohort medians 39 years (range 17-67) for COPEP and 41 years (range 19-78) for PLWH (Thoueille 2023 Table 1)."
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on the PK parameters; not retained in the final model. Cohort 60 female / 45 male (COPEP) and 69 female / 50 male (PLWH) (Thoueille 2023 Table 1). Note that sex nonetheless entered the analysis indirectly, outside the PK model: the assumed haematocrit used to convert the COPEP dried-blood-spot concentrations to plasma concentrations was set to 0.40 for women and 0.45 for men (Thoueille 2023 Methods, 'Population pharmacokinetic analysis'). That conversion is a data-preparation step applied before model fitting, not a covariate effect inside the PK model, so it is not encoded here."
    ),
    HT = list(
      description = "Body height.",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened on the PK parameters; not retained in the final model. Cohort medians 172 cm (range 149-192) for COPEP and 168 cm (range 148-190) for PLWH (Thoueille 2023 Table 1)."
    ),
    BMI = list(
      description = "Body mass index.",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened on the PK parameters; not retained in the final model (body weight was the retained size descriptor). Cohort median 24 kg/m^2 in both populations (ranges 17-49 and 16-51) (Thoueille 2023 Table 1)."
    ),
    SMOKE_CURRENT = list(
      description = "Current-smoker indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on the PK parameters; not retained in the final model. Thoueille 2023 does not tabulate the smoking-status distribution."
    ),
    HIV_POS = list(
      description = "Study-population indicator: 1 = person living with HIV enrolled in the Swiss HIV Cohort Study routine TDM programme, 0 = HIV-negative COPEP participant receiving LPV/r as COVID-19 post-exposure prophylaxis.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on ALL PK parameters and not retained: 'PopPK parameters between COPEP participants and PLWH did not differ significantly, thus supporting the use of a unique model for both COPEP and TDM data' (Thoueille 2023 Results, 'Structural, statistical and covariate models'). This null result is the central finding of the paper, so the single pooled parameterisation encoded in ini() applies to both populations."
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 224L,
    n_studies       = 2L,
    n_observations  = 275L,
    age_range       = "17-78 years (COPEP 17-67; PLWH 19-78)",
    age_median      = "39 years (COPEP) and 41 years (PLWH)",
    weight_range    = "40-157 kg (COPEP 47-157; PLWH 40-147)",
    weight_median   = "73 kg (COPEP) and 66 kg (PLWH)",
    height_range    = "148-192 cm (COPEP 149-192; PLWH 148-190)",
    height_median   = "172 cm (COPEP) and 168 cm (PLWH)",
    bmi_range       = "16-51 kg/m^2 (COPEP 17-49; PLWH 16-51)",
    bmi_median      = "24 kg/m^2 in both populations",
    sex_female_pct  = 57.6,
    disease_state   = "Two pooled populations. (1) 105 HIV-negative adults exposed to SARS-CoV-2 (more than 15 min at under 2 m, or a shared closed space for more than 2 h, with a person with confirmed infection) who received LPV/r as COVID-19 post-exposure prophylaxis in the COPEP trial. (2) 119 people living with HIV (PLWH) enrolled in the Swiss HIV Cohort Study and followed in the routine therapeutic drug monitoring programme of the Service of Clinical Pharmacology, Lausanne, between January 2010 and May 2022.",
    dose_range      = "COPEP: lopinavir/ritonavir 400/100 mg twice daily for 5 days (all 105 participants). PLWH on routine TDM: lopinavir 200 mg (9), 300 mg (1), 400 mg (85), 500 mg (5), 600 mg (7), 800 mg (10) and 1000 mg (2) per administration (Thoueille 2023 Table 1).",
    regions         = "Switzerland (COPEP recruitment in Geneva, Basel and Lugano; TDM programme in Lausanne).",
    notes           = "105 lopinavir concentrations from 105 COPEP participants (one dried-blood-spot sample each, taken on day 5 with the time of last intake documented; median time after dose 1.5 h, range 0.07-19 h) plus 170 sparse plasma concentrations from 119 PLWH (median 2 samples per patient, range 1-8; median time after dose 10.25 h, range 1.25-29.5 h). COPEP dried-blood-spot concentrations were converted to plasma concentrations as C_plasma = C_DBS * F_BP / (1 - HCT), with the lopinavir protein-binding ratio F_BP = 98.5% and haematocrit assumed at 0.40 for women and 0.45 for men (haematocrit was not measured in COPEP). Twelve COPEP measurements below 1000 ng/mL were excluded as indicating absolute non-adherence. ALL individuals were assumed to be at steady state (full adherence for COPEP, long treatment duration for the SHCS patients), so simulations reproducing this model should dose to steady state. Lopinavir was quantified by LC-MS/MS (LLOQ 10 ng/mL). NONMEM 7.4.3; final estimates from Thoueille 2023 Table 2, confirmed by a 2000-sample bootstrap and by prediction- and variability-corrected VPC (Fig. 1)."
  )

  ini({
    # ---- Structural parameters: Thoueille 2023 Table 2, "Final model /
    # Estimate" column. Reference subject is a 70 kg adult (Table 2
    # footnote: "theta_BW bodyweight effect on CL_LPV with reference
    # bodyweight of 70 kg"). Lopinavir was given orally only, so CL and V
    # are apparent (CL/F, V/F); the paper reports no bioavailability
    # parameter and none is introduced here.
    lka <- log(0.76)
    label("First-order absorption rate constant (ka, 1/h)")                  # Thoueille 2023 Table 2: ka = 0.76 1/h (RSE 3%; bootstrap median 0.75)
    lvc <- log(78.9)
    label("Apparent volume of distribution (V_LPV/F, L)")                    # Thoueille 2023 Table 2: V_LPV = 78.9 L (RSE 2%; bootstrap median 77.7, 95% CI 52.2-119.5)
    lcl <- log(4.02)
    label("Apparent oral clearance at the 70 kg reference weight (CL_LPV/F, L/h)")  # Thoueille 2023 Table 2: CL_LPV = 4.02 L/h (RSE 3%; bootstrap median 4.00, 95% CI 3.71-4.32)

    # ---- Covariate effect on CL/F. Thoueille 2023 prints the final-model
    # covariate equation directly beneath Table 2:
    #
    #   TVCL_LPV = CL_LPV * [1 + theta_BW * (BW - 70) / 70]
    #
    # The linear (rather than allometric) form was deliberately retained:
    # "Allometric scaling described equally well the effect of weight on
    # CL_LPV, but the linear function was retained for simplicity"
    # (Results, "Structural, statistical and covariate models").
    #
    # Cross-check on the coefficient: the same paragraph reports "a 19%
    # higher CL_LPV in a person of 100 kg vs 70 kg". The linear form gives
    # 0.447 * (100 - 70) / 70 = 0.1916, i.e. +19.2%, reproducing the
    # reported 19%. The allometric alternative would give
    # (100/70)^0.447 - 1 = +17.3%, which does not, confirming the linear
    # reading of the printed equation.
    e_wt_cl <- 0.447
    label("Proportionality coefficient relating CL/F to the relative deviation of body weight from 70 kg (unitless)")  # Thoueille 2023 Table 2: theta_BW = 0.447 (RSE 14%; bootstrap median 0.447, 95% CI 0.196-0.733)

    # ---- Between-subject variability. Thoueille 2023 Table 2 reports BSV
    # only on CL_LPV (28.5%, RSE 14%; bootstrap median 28.0%, 95% CI
    # 18.8-35.9%). BSV on V_LPV and ka was tested and rejected: "The
    # assignment of IIV on V_LPV and k_a did not improve data description
    # (delta-OFV = 0, p > 0.05)" (Results). The reported percentage is
    # read as the coefficient of variation of a log-normal (exponential)
    # random effect and back-transformed with omega^2 = log(1 + CV^2):
    #   omega^2 = log(1 + 0.285^2) = 0.078096
    etalcl ~ 0.078096  # Thoueille 2023 Table 2: BSV on CL_LPV = 28.5%

    # ---- Residual unexplained variability: combined proportional plus
    # additive ("A combined error model best described LPV residual
    # unexplained variability", Results). This model carries
    # concentrations in mg/L, so the paper's additive term of 1560 ng/mL
    # becomes 1.560 mg/L (1 mg/L = 1000 ng/mL).
    propSd <- 0.333
    label("Proportional residual error (fraction)")                          # Thoueille 2023 Table 2: sigma_prop = 33.3% (RSE 11%; bootstrap median 32.7%)
    addSd <- 1.560
    label("Additive residual error (mg/L)")                                  # Thoueille 2023 Table 2: sigma_add = 1560 ng/mL = 1.560 mg/L (RSE 1%; bootstrap median 1544 ng/mL)
  })

  model({
    # Apparent oral clearance with the linear body-weight effect
    # (Thoueille 2023 final-model equation beneath Table 2):
    #   TVCL_LPV = CL_LPV * [1 + theta_BW * (BW - 70) / 70]
    cl <- exp(lcl + etalcl) * (1 + e_wt_cl * (WT - 70) / 70)

    # Apparent volume of distribution: no covariates and no BSV retained.
    vc <- exp(lvc)

    # First-order absorption rate constant: no BSV retained.
    ka <- exp(lka)

    # Elimination micro-constant.
    kel <- cl / vc

    # One-compartment open model with first-order absorption from an oral
    # depot and first-order elimination ("a one-compartment model with
    # first-order absorption and elimination best described LPV plasma
    # concentrations", Results).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Observation. Doses in mg with vc in L give Cc in mg/L; the paper
    # tabulates concentrations in ng/mL (1 mg/L = 1000 ng/mL).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
