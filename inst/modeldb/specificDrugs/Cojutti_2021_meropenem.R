Cojutti_2021_meropenem <- function() {
  description <- paste(
    "One-compartment population PK model for meropenem given by continuous",
    "intravenous infusion to 74 critically ill adults (one adolescent) with",
    "Gram-negative infections at a tertiary hospital in Bologna, Italy, fitted",
    "non-parametrically with the NPAG algorithm in Pmetrics 1.5.0 to 183",
    "steady-state therapeutic-drug-monitoring concentrations. Clearance is",
    "ADDITIVE in an intercept and an arm linear in CKD-EPI creatinine clearance",
    "(CL = theta1 + theta2 * CLCR), and volume is a power function of absolute",
    "total body weight (V = theta3 * BW^theta4). All four thetas carry",
    "inter-individual variability taken from the CVs of the NPAG marginal",
    "distributions, including the body-weight exponent. Residual error is the",
    "Pmetrics assay-error polynomial (0.0798 + 0.0927 * C) multiplied by the",
    "gamma factor of 2.",
    sep = " "
  )
  reference <- paste(
    "Cojutti PG, Gatti M, Rinaldi M, Tonetti T, Laici C, Mega C, Siniscalchi A,",
    "Giannella M, Viale P, Pea F. Impact of Maximizing Css/MIC Ratio on Efficacy",
    "of Continuous Infusion Meropenem Against Documented Gram-Negative Infections",
    "in Critically Ill Patients and Population Pharmacokinetic/Pharmacodynamic",
    "Analysis to Support Treatment Optimization. Front Pharmacol.",
    "2021;12:781892. doi:10.3389/fphar.2021.781892. PMCID: PMC8694396.",
    "Structural and variability estimates are Table 4 ('Final model' columns);",
    "the covariate equations are printed in Section 3.3; the residual-error",
    "polynomial and gamma are in Section 2.4. No supplement was published and",
    "no erratum was found (Europe PMC search, 2026-09-30).",
    sep = " "
  )
  vignette <- "Cojutti_2021_meropenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Total (bound + unbound) plasma meropenem, measured by LC-MS/MS (Section
  # 2.1). The model is fitted to total concentrations; no protein-binding
  # conversion is applied (meropenem is about 2% bound).
  compartmentData <- list(
    central = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Creatinine clearance estimated with the CKD-EPI equation (Levey 2009),",
        "BSA-normalized to 1.73 m^2"
      ),
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Section 2.1: 'CLCR was estimated by means of the CKD-EPI formula'. Table",
        "1 reports it in 'ml/min/1.73 m2', median 91.5 (IQR 51.2-114.9); the",
        "Results give the range 7-192. It enters clearance UNCENTRED and",
        "LINEARLY: CL = theta1 + theta2 * CLCR, with theta1 the clearance at",
        "CLCR = 0 (Section 3.3). The mL/min/1.73 m^2 scale is confirmed by",
        "arithmetic: 1.040 + 0.103 * 91.5 = 10.5 L/h, against the Results median",
        "individual CL of 7.27 (IQR 4.53-10.41) L/h and the Table 1 median Css",
        "of 14.1 mg/L at a median 3 g/day (125 mg/h / 14.1 = 8.9 L/h). Patients",
        "on renal replacement therapy were excluded (Section 2.1). The Monte",
        "Carlo simulations drew CLCR uniformly within the classes 0-29, 30-79,",
        "80-129 and 130-200 (Section 3.4)."
      ),
      source_name = "CLCR"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters volume as an UNCENTRED power of absolute weight, V = theta3 *",
        "BW^theta4, 'theta3 is the distribution volume when BW = 1' (Section",
        "3.3; the PDF typesets the equation as 'Vi = theta3 * (BWi)^theta4').",
        "Table 1: median 79.0 kg (IQR 68.5-89.5); Results: range 50-160 kg.",
        "Baseline value."
      ),
      source_name = "BW"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 74L,
    n_studies = 1L,
    n_observations = 183L,
    age_range = "12-86 years (mean 60.1, SD 15.0)",
    weight_range = "50-160 kg (median 79.0, IQR 68.5-89.5)",
    sex_female_pct = 29.7,
    renal_function = "CKD-EPI CLCR 7-192 mL/min/1.73 m^2 (median 91.5, IQR 51.2-114.9); 13.5% augmented renal clearance (>= 130); renal replacement therapy excluded",
    disease_state = paste(
      "Critically ill patients (ICU) treated with continuous-infusion meropenem",
      "for Gram-negative infections: HAP/VAP 50.0%, bloodstream infection 25.7%,",
      "complicated intra-abdominal 14.9%, complicated urinary tract 5.4%,",
      "meningitis 2.7%, bone and joint 1.4%. SOFA median 7 (IQR 5-10)."
    ),
    dose_range = paste(
      "Loading dose 2 g over 2 h, then 1 g q6h over 6 h (CLCR >= 60) or 0.5 g",
      "q6h over 6 h (CLCR < 60), i.e. continuous infusion, adjusted by TDM;",
      "median 1 g q8h CI (IQR 0.5 g q6h - 1 g q6h)"
    ),
    regions = "Italy (IRCCS Azienda Ospedaliero-Universitaria di Bologna)",
    notes = paste(
      "Retrospective TDM cohort (Table 1). Median 2 (IQR 1-3) TDM samples per",
      "patient, drawn after at least 2 days of therapy; median Css 14.1 (IQR",
      "9.0-21.3) mg/L. Twenty-two of 74 (29.7%) were female."
    )
  )

  ini({
    # =====================================================================
    # Pmetrics / NPAG. Table 4 reports, for each theta, the MEAN of the
    # non-parametric marginal support-point distribution and its CV%. The
    # means are encoded as the typical (median) values of log-normal
    # marginals and the CVs converted with omega^2 = log(CV^2 + 1) -- the
    # convention of Sime_2019_ceftolozane.R and
    # Hughes_2024_vancomycin_nonparametric.R. The NPAG joint density is not
    # published, so the four marginals are INDEPENDENT here.
    # =====================================================================

    # --- Clearance: CL = theta1 + theta2 * CLCR (Section 3.3) ---------------
    lcl_nonren <- log(1.040)
    label("Clearance intercept at CLCR = 0, theta1 (L/h)")
    # Table 4, theta1: mean 1.040 (CV 77.016%); bootstrap median 0.908
    e_crcl_cl_renal <- 0.103
    label("Slope of CL on CKD-EPI CLCR, theta2 (L/h per mL/min/1.73 m^2)")
    # Table 4, theta2: mean 0.103 (CV 66.074%); bootstrap median 0.080

    # --- Volume: V = theta3 * BW^theta4 (Section 3.3) ------------------------
    lvc <- log(7.343)
    label("Volume of distribution at BW = 1 kg, theta3 (L)")
    # Table 4, theta3: mean 7.343 (CV 46.824%); bootstrap median 9.288
    e_wt_vc <- 0.612
    label("Power exponent of absolute body weight on V, theta4 (unitless)")
    # Table 4, theta4: mean 0.612 (CV 59.146%); bootstrap median 0.554
    # At 79 kg these means give V = 106 L, whereas Section 3.3 reports a
    # median individual (posterior) V of 20.0 L (IQR 17.16-23.59). Encoded as
    # printed: the Figure 5 PTA curves are reproduced only with the large V
    # (see the vignette, 'The loading dose, and why the volume is large').

    # --- Inter-individual variability: Table 4 'CV (%)' of the final model ---
    etalcl_nonren ~ 0.465711 # theta1 CV 77.016%: log(0.77016^2 + 1)
    etae_crcl_cl_renal ~ 0.362263 # theta2 CV 66.074%: log(0.66074^2 + 1)
    etalvc ~ 0.198235 # theta3 CV 46.824%: log(0.46824^2 + 1)
    etae_wt_vc ~ 0.299975 # theta4 CV 59.146%: log(0.59146^2 + 1)

    # --- Residual error (Section 2.4) -----------------------------------------
    # Assay-error polynomial 'coefficients of the four-term polynomial
    # functions were 0.0798, 0.0927, 0 and 0', i.e. SD = 0.0798 + 0.0927 * C,
    # multiplied by the extra-process-noise 'gamma model (gamma = 2)'. The
    # Pmetrics SD is the SUM gamma * (C0 + C1 * C), hence combined1 below.
    addSd <- 0.1596
    label("Additive residual SD, gamma * C0 = 2 * 0.0798 (mg/L)")
    # Section 2.4: C0 = 0.0798, gamma = 2
    propSd <- 0.1854
    label("Proportional residual SD, gamma * C1 = 2 * 0.0927 (fraction)")
    # Section 2.4: C1 = 0.0927, gamma = 2
  })

  model({
    # Clearance, additive in an intercept and a CLCR-proportional arm. Each
    # arm carries its own IIV because Table 4 reports a CV for each theta.
    cl_nonren <- exp(lcl_nonren + etalcl_nonren)
    cl_renal <- e_crcl_cl_renal * exp(etae_crcl_cl_renal) * CRCL
    cl <- cl_nonren + cl_renal

    # Volume, power of ABSOLUTE body weight (theta3 = V at BW = 1 kg). The
    # exponent is itself a random effect in the NPAG fit (Table 4 CV 59.146%).
    vc <- exp(lvc + etalvc) * WT^(e_wt_vc * exp(etae_wt_vc))

    kel <- cl / vc

    # One compartment, zero-order input (infusion given in the event table),
    # first-order elimination (Section 2.4).
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
