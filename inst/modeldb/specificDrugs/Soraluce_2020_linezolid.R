Soraluce_2020_linezolid <- function() {
  description <- "Two-compartment IV population PK model for linezolid in 40 critically ill adults, 23 of them on continuous renal replacement therapy (Soraluce 2020). Total clearance is the sum of an estimated non-renal clearance (2.62 L/h), a renal clearance proportional to urine-measured creatinine clearance (4.35 L/h at 44 mL/min), and the individually measured extracorporeal clearance (sieving coefficient times effluent flow) supplied as the data column QEFF; a single exponential random effect scales the whole sum. Central volume 16.2 L, peripheral volume 29.0 L, intercompartmental clearance 71.7 L/h, with IIV on clearance and central volume and a combined additive + proportional residual error."
  reference <- "Soraluce A, Barrasa H, Asin-Prieto E, Sanchez-Izquierdo JA, Maynar J, Isla A, Rodriguez-Gascon A. Novel Population Pharmacokinetic Model for Linezolid in Critically Ill Patients and Evaluation of the Adequacy of the Current Dosing Recommendation. Pharmaceutics. 2020;12(1):54. doi:10.3390/pharmaceutics12010054"
  vignette <- "Soraluce_2020_linezolid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "linezolid", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "linezolid", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance MEASURED from a 10-hour urine collection, Clcr = (Cru * Vu) / (Crp * 600 min); raw mL/min, NOT BSA-normalized",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Soraluce 2020 Section 2.1: Clcr (mL/min) = (Cru x Vu) / (Crp x 600 min), with Cru and Crp",
        "the urine and plasma creatinine concentrations (mg/dL) and Vu the urine volume (mL)",
        "collected over 10 h. Not a Cockcroft-Gault estimate; the Discussion attributes the",
        "strength of the Clcr-clearance relationship in this cohort to the measured value.",
        "Enters the renal clearance arm as the through-origin linear ratio (CRCL / 44), so the",
        "4.35 L/h renal clearance is the value at 44 mL/min. Table 1: median 71.2 mL/min (range",
        "11.0-179.5) without CRRT and 6.0 mL/min (range 0.0-45.6) with CRRT; the external",
        "validation cohort on continuous infusion had a median of 111 mL/min (range 45-240)",
        "(Table 2). A value of 0 (anuric) is admissible and switches the renal arm off."
      ),
      source_name = "Clcr"
    ),
    QEFF = list(
      description = "Individually determined extracorporeal (CRRT) clearance of linezolid, CLEC = Sc * Qef; 0 for subjects not on CRRT",
      units = "L/h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Soraluce 2020 Section 2.1: the sieving coefficient Sc is the ratio of the effluent to the",
        "plasma linezolid AUC over the dosing interval, and CLEC = Sc x Qef with Qef the effluent",
        "flow. Section 2.3.2 and Table 3 footnote c: CLEC was 'considered as a fixed value per",
        "patient' and 'included as the value calculated for each patient ... only for those",
        "undergoing CRRT', i.e. it is a data column, not an estimated parameter, and it is ADDED",
        "to the non-renal and renal arms before the random effect is applied. Set QEFF = 0 for",
        "subjects not on CRRT (there is no separate on/off indicator in the model). Table 1:",
        "median 2.51 L/h (range 0.79-3.09) across the 23 CRRT patients, 2.61 (0.79-3.09) for 18",
        "on CVVHDF and 1.06 (1.00-2.73) for 5 on CVVHD; effluent flow 1.1-3.3 L/h and mean Sc",
        "about 0.8 (Results 3; Discussion)."
      ),
      source_name = "CLEC"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "year",
      type = "continuous",
      notes = "Soraluce 2020 Section 2.3.2: every Table 1 variable was screened by SCM (forward p < 0.05, backward p < 0.01); not retained. Table 1 median 72 years (22-85) without CRRT, 68 (37-79) with CRRT."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened (Section 2.3.2); not retained. Table 1 median 71 kg (60-95) without CRRT, 74 kg (55-110) with CRRT."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = "Screened (Section 2.3.2); not retained. Table 1 median 2.8 g/dL (1.9-4.0) without CRRT, 2.2 (1.7-3.6) with CRRT."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "mg/dL",
      type = "continuous",
      notes = "Section 2.3.2: grouped with GOT and GPT into three dichotomous out-of-normal-range liver-function covariates; Results 3.1.2: 'excluded from the final model, since a better fit was not obtained'."
    ),
    APACHE_II = list(
      description = "Acute Physiology and Chronic Health Evaluation (APACHE) II score",
      units = NA_character_,
      type = "continuous",
      notes = "Screened (Section 2.3.2); not retained. Table 1 median 16 (11-36) without CRRT, 22 (16-34) with CRRT."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 40L,
    n_studies = 1L,
    age_range = "22-85 years (median 72 without CRRT, 68 with CRRT)",
    weight_range = "55-110 kg (median 71 kg without CRRT, 74 kg with CRRT)",
    sex_female_pct = 27.5,
    race_ethnicity = "Not reported (Spanish multicentre cohort)",
    disease_state = "Critically ill adults in the intensive care unit treated with linezolid for a suspected Gram-positive infection (pulmonary 14, abdominal 10, neurological 9, biliary 2, other 5). 23 of 40 received continuous renal replacement therapy (18 CVVHDF, 5 CVVHD). APACHE II median 16 (11-36) without CRRT and 22 (16-34) with CRRT.",
    renal_function = "Urine-measured creatinine clearance median 71.2 mL/min (11.0-179.5) in the 17 patients without CRRT and 6.0 mL/min (0.0-45.6) in the 23 on CRRT. Extracorporeal clearance median 2.51 L/h (0.79-3.09); effluent flow 1.1-3.3 L/h.",
    dose_range = "Linezolid 600 mg IV every 12 h as a 30-min infusion (one subject 60 min); sampled at steady state after a mean of 8 doses.",
    regions = "Spain (Araba University Hospital, Vitoria-Gasteiz; Doce de Octubre University Hospital, Madrid; Joan XXIII University Hospital, Tarragona).",
    n_concentrations = 311L,
    notes = "Prospective open-label multicentre study; plasma (and effluent) sampled pre-dose, end of infusion and 1, 2, 3, 6, 8-10 and 12 h. HPLC-UV, plasma LLOQ 0.5 mg/L. NONMEM 7.3 FOCE-I on log-transformed concentrations; SCM, VPC and 2000-sample bootstrap in PsN 4.7.0. External validation in 11 further non-CRRT patients given a 600 mg loading dose then 50 mg/h continuous infusion (Table 2, Figure 3)."
  )

  ini({
    # Final model, Soraluce 2020 Table 3 footnote b:
    #   CL = (2.62 + 4.35 * (Clcr / 44) + Sc * Qef) * exp(eta1)
    #   V1 = 16.2 * exp(eta2)
    # Table 3, 'Final Model Estimate, RSE (%)' column for all values below.
    # Table 3 row CLNR: 2.62 L/h (RSE 18%; bootstrap median 2.65, 5th-95th 2.02-3.65)
    lcl_nonren <- log(2.62); label("Non-renal clearance CLNR (L/h)")
    # Table 3 row 'CLR = theta x (Clcr/44)': 4.35 L/h (RSE 19%; bootstrap 4.33, 2.99-5.84)
    lcl_renal <- log(4.35); label("Renal clearance CLR at CRCL = 44 mL/min (L/h)")
    # Table 3 row V1: 16.2 L (RSE 14%; bootstrap 16.6, 11.7-24.4)
    lvc <- log(16.2); label("Central volume of distribution V1 (L)")
    # Table 3 row Q: 71.7 L/h (RSE 14%; bootstrap 69.5, 40.4-92.0)
    lq <- log(71.7); label("Intercompartmental clearance Q (L/h)")
    # Table 3 row V2: 29.0 L (RSE 7%; bootstrap 28.6, 23.0-32.6)
    lvp <- log(29.0); label("Peripheral volume of distribution V2 (L)")

    # Exponential IIV (Results 3.1.1), no CL-V1 correlation. Table 3 prints
    # IIV as a percentage; converted with omega^2 = log(CV^2 + 1). The
    # Figure 3 inset (simulated Css median 3.30, 2.5th-97.5th 0.85-14.73 mg/L)
    # is reproduced on this scale and not on omega = CV^2; see the vignette.
    # Table 3 row IIV_CL: 61.5% (RSE 9%; shrinkage 1%)
    etalcl ~ 0.320546 # log(0.615^2 + 1)
    # Table 3 row IIV_V1: 65.9% (RSE 17%; shrinkage 18%)
    etalvc ~ 0.360602 # log(0.659^2 + 1)

    # Combined residual error (Results 3.1.1). NONMEM fitted log-transformed
    # concentrations; encoded as additive + proportional on the linear scale.
    # Table 3 row 'Residual error_additive (mg/L)': 0.266 (RSE 24%; bootstrap 0.267)
    addSd <- 0.266; label("Additive residual error (mg/L)")
    # Table 3 row 'Residual error_proportional (%)': 0.159 (RSE 19%; bootstrap 0.157)
    propSd <- 0.159; label("Proportional residual error (fraction)")
  })
  model({
    # Total clearance, Table 3 footnote b. The renal arm is the through-origin
    # linear ratio CRCL/44; QEFF (Sc * Qef, L/h) is the per-subject measured
    # extracorporeal clearance, 0 when the subject is not on CRRT. The single
    # random effect multiplies the whole sum, including QEFF.
    #
    # The total must be the variable named `cl`: rxode2 recognises a
    # cl/vc/q/vp set and may solve the system from those names, so any other
    # name for the sum would silently drop the renal and CRRT arms.
    cl_nonren <- exp(lcl_nonren)
    cl_renal <- exp(lcl_renal) * (CRCL / 44)
    cl <- (cl_nonren + cl_renal + QEFF) * exp(etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg and vc in L give mg/L (ug/mL), the units of the HPLC-UV plasma
    # assay the model was fitted to.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
