Zuo_2020_zanamivir <- function() {
  description <- "Two-compartment population PK model for intravenous zanamivir with linear elimination, a piecewise-linear (hinge) creatinine-clearance effect on CL, a hospitalized-patient (suspected or confirmed influenza) effect on CL and on the magnitude of CL IIV, estimated allometric weight exponents on V1/V2 and Q, and a study effect on V1/V2, in healthy adults and hospitalized adult and pediatric subjects with influenza (Zuo 2020)"
  reference <- paste(
    "Zuo P, Collins J, Okour M, Barth A, Shortino D, Yates P, Roberts G,",
    "Watson HA, Peppercorn A, Hossain M. Population pharmacokinetic/",
    "pharmacodynamic analysis of intravenous zanamivir in healthy adults and",
    "hospitalized adult and pediatric subjects with influenza.",
    "Clin Transl Sci. 2020;13(1):157-168. doi:10.1111/cts.12697.",
    "Parameter estimates are the final-model column of Zuo 2020 Table 2; the",
    "creatinine-clearance equation is printed in the PopPK analysis Results",
    "paragraph; validation targets are Supplementary Tables S4 and S6.",
    sep = " "
  )
  vignette <- "Zuo_2020_zanamivir"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on V1 and V2 with one shared estimated exponent (0.711) and on Q (exponent 0.658), reference weight 70 kg (Zuo 2020 Table 2 rows 'V1/V2 ~ WT' and 'Q ~ WT'; Figure 1 caption defines the reference individual as a 70 kg healthy subject). CL carries no weight effect in the final model; weight enters CL only indirectly through the Cockcroft-Gault / Schwartz creatinine clearance.",
      source_name = "WT"
    ),
    CRCL = list(
      description = "Creatinine clearance: Cockcroft-Gault in mL/min for subjects aged 13 years or older, bedside Schwartz in mL/min/1.73 m^2 for subjects younger than 13 years; for subjects on renal replacement therapy, computed from ultrafiltration rate and blood flow",
      units = "mL/min (>= 13 years) or mL/min/1.73 m^2 (< 13 years)",
      type = "continuous",
      reference_category = NULL,
      notes = "Piecewise-linear (hinge) effect on CL: factor = 1 for CRCL >= 97 mL/min and 1 + 0.00929 x (CRCL - 97) below it (Zuo 2020 PopPK analysis Results paragraph and Table 2 rows 'CL ~ CrCL (inflection point)' = 97.0 and 'CL ~ CrCL (slope)' = 0.00929). The Methods state CrCL was derived by Cockcroft-Gault for adults and by the Schwartz equation for infants and children under 13 years; Supplementary Tables S4-S6 footnotes give the units as mL/min for >= 13 years and mL/min/1.73 m^2 for < 13 years, and footnote a of Table S6 states that CrCL for subjects on renal replacement therapy was calculated from ultrafiltration rate and blood flow. The single column therefore mixes two estimating equations on two scales, exactly as the source dataset did.",
      source_name = "CrCL"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-participant indicator (1 = healthy adult volunteer from the six phase I studies, 0 = hospitalized subject with suspected or confirmed influenza from the phase II or phase III study)",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (healthy adult volunteer) is the structural reference of the paper; 0 (hospitalized patient) applies both the 0.756 CL multiplier and the 3.10 IIV scaling",
      notes = "The source flag is FLU (Table 2 rows 'CL ~ FLU' and 'IIV (CL) ~ FLU'), coded 1 for the 533 hospitalized subjects with suspected or confirmed influenza; model() computes flu = 1 - DIS_HEALTHY. The patient group includes the 11% of hospitalized subjects whose influenza test was negative or unknown (Table 1 'Influenza virus type unknown or negative'), which is why the healthy-participant canonical is used rather than an infection-specific indicator. Direction confirmed by the Results: CL for a typical healthy 70 kg subject is 6.82 L/h and 'decreased by 24% to 5.16 L/h if these subjects were infected with influenza' (6.82 x 0.756 = 5.16).",
      source_name = "FLU (inverted)"
    ),
    STUDY_NAI114346 = list(
      description = "Study NAI114346 indicator (1 = record from phase I study 4, the thorough-QT cardiac-conduction study NCT01353729 in healthy adults; 0 = any other study)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all other studies)",
      notes = "Multiplicative effect of 0.729 on both V1 and V2 (Zuo 2020 Table 2 row 'V1/V2 ~ Study (NAI114346)'). The Discussion notes the effect does not influence the patient simulations because NAI114346 enrolled only healthy volunteers; set STUDY_NAI114346 = 0 for any prospective simulation.",
      source_name = "STDY"
    )
  )

  covariatesDataExcluded <- list(
    RRT = list(
      description = "Renal replacement therapy indicator",
      units = "(binary)",
      type = "binary",
      notes = "Tested in the full covariate model on CL (0.900, 95% CI 0.645-1.155; Zuo 2020 Table 2 full-model column) but dropped from the final model because its effect range overlapped the null range (Figure 1)."
    ),
    ECMO = list(
      description = "Extracorporeal membrane oxygenation indicator",
      units = "(binary)",
      type = "binary",
      notes = "Tested in the full covariate model on CL (0.704, 95% CI 0.416-0.992; Zuo 2020 Table 2 full-model column) but dropped from the final model because its effect range overlapped the null range (Figure 1)."
    )
  )

  compartmentData <- list(
    central = list(analyte = "zanamivir", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "zanamivir", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 658L,
    n_studies = 8L,
    n_observations = 5273L,
    age_range = "0.6-101 years",
    age_median = "46 years",
    weight_range = "7.1-188.0 kg",
    weight_median = "72.0 kg",
    crcl_range = "12.1-274.9 mL/min (median 91.8)",
    sex_female_pct = 38.6,
    disease_state = "125 healthy adult volunteers (six phase I studies) and 533 hospitalized adult, adolescent and pediatric subjects with suspected or confirmed influenza (phase II NAI113678 and phase III NAI114373), including 57% of patients with renal impairment, 16 on renal replacement therapy and 10 on ECMO",
    dose_range = "Intravenous infusion: single doses of 100-1,200 mg and 300-600 mg twice daily; renally impaired and pediatric patients received an initial dose followed by a renal-function- and weight-adjusted twice-daily maintenance dose for 5-10 days (Zuo 2020 Methods and Table S1)",
    regions = "Multinational (Europe, North/South America, Asia, Africa); phase I studies in the UK/US, Thailand, Japan and China",
    notes = "Demographics from Zuo 2020 Table 1 (median (range)). 19% of subjects were healthy and 81% hospitalized; 58 pediatric or adolescent subjects contributed 289 observations. Serum zanamivir by LC-MS/MS; LLOQ 0.02-10 ng/mL across studies, and concentrations below the 0.01 ug/mL modelling LLOQ were imputed at the LLOQ. Estimation in NONMEM 7.2 with FOCE-I; full covariate model approach."
  )

  ini({
    # Structural parameters for the reference individual: healthy subject,
    # 70 kg, CrCL >= 97 mL/min, not in study NAI114346 (Figure 1 caption).
    lcl <- log(6.82); label("Clearance CL for a healthy 70 kg subject with CrCL >= 97 mL/min (L/h)") # Zuo 2020 Table 2 final model: CL = 6.82 L/h (RSE 2%)
    lvc <- log(12.3); label("Central volume of distribution V1 at 70 kg (L)") # Zuo 2020 Table 2 final model: V1 = 12.3 L (RSE 2%)
    lq <- log(4.82); label("Intercompartmental clearance Q at 70 kg (L/h)") # Zuo 2020 Table 2 final model: Q = 4.82 L/h (RSE 8%)
    lvp <- log(6.52); label("Peripheral volume of distribution V2 at 70 kg (L)") # Zuo 2020 Table 2 final model: V2 = 6.52 L (RSE 4%)

    # Creatinine-clearance hinge on CL: flat at 1 above the knot, linear below.
    lcrcl_hinge <- log(97.0); label("Creatinine clearance knot of the piecewise-linear CL relationship (mL/min)") # Zuo 2020 Table 2 final model: CL ~ CrCL (inflection point) = 97.0 mL/minute (RSE 2%)
    e_crcl_cl <- 0.00929; label("Slope of the fractional CL change per unit CrCL below the knot (min/mL)") # Zuo 2020 Table 2 final model: CL ~ CrCL (slope) = 0.00929 minute/mL (RSE 7%); Results prose quotes the full-model value 0.00923

    # Hospitalized-patient (FLU) effects.
    e_flu_cl <- 0.756; label("Multiplicative effect of being a hospitalized influenza patient on CL (unitless)") # Zuo 2020 Table 2 final model: CL ~ FLU = 0.756 (RSE 4%)
    e_flu_etalcl <- 3.10; label("Multiplicative scaling of the CL random effect in hospitalized influenza patients (unitless)") # Zuo 2020 Table 2 final model: IIV (CL) ~ FLU = 3.10 (RSE 11%)

    # Weight and study effects on distribution.
    e_wt_vc_vp <- 0.711; label("Power exponent of body weight on V1 and V2, reference 70 kg (unitless)") # Zuo 2020 Table 2 final model: V1/V2 ~ WT = 0.711 (RSE 4%)
    e_wt_q <- 0.658; label("Power exponent of body weight on Q, reference 70 kg (unitless)") # Zuo 2020 Table 2 final model: Q ~ WT = 0.658 (RSE 9%)
    e_study_nai114346_vc_vp <- 0.729; label("Multiplicative effect of study NAI114346 on V1 and V2 (unitless)") # Zuo 2020 Table 2 final model: V1/V2 ~ Study (NAI114346) = 0.729 (RSE 2%)

    # IIV: variances from CV% via omega^2 = log(CV^2 + 1).
    etalcl ~ 0.034372 # Zuo 2020 Table 2 final model: IIV CL = 18.7 CV% (healthy reference; scaled by e_flu_etalcl in patients)
    etalvc ~ 0.114313 # Zuo 2020 Table 2 final model: IIV V1 = 34.8 CV%

    # Combined residual error.
    propSd <- 0.262; label("Proportional residual error (fraction)") # Zuo 2020 Table 2 final model: proportional error = 26.2 CV%
    addSd <- 0.0269; label("Additive residual error (ug/mL)") # Zuo 2020 Table 2 final model: additive error = 0.0269 ug/mL SD
  })
  model({
    flu <- 1 - DIS_HEALTHY

    # Piecewise-linear CrCL effect (Results: 1 + (CrCL - 97) x slope below 97).
    crcl_hinge <- exp(lcrcl_hinge)
    crcl_eff <- min(CRCL, crcl_hinge)
    f_crcl <- 1 + e_crcl_cl * (crcl_eff - crcl_hinge)

    # CL random effect is scaled by 3.10 in hospitalized patients.
    eta_cl_scale <- 1 + (e_flu_etalcl - 1) * flu

    cl <- exp(lcl + etalcl * eta_cl_scale) * f_crcl * e_flu_cl^flu
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp * e_study_nai114346_vc_vp^STUDY_NAI114346
    vp <- exp(lvp) * (WT / 70)^e_wt_vc_vp * e_study_nai114346_vc_vp^STUDY_NAI114346
    q <- exp(lq) * (WT / 70)^e_wt_q

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d / dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
