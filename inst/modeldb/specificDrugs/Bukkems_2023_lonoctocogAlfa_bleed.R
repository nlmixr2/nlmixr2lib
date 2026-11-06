Bukkems_2023_lonoctocogAlfa_bleed <- function() {
  description <- "Repeated time-to-event (RTTE) model for bleeding in children, adolescents and adults with severe hemophilia A on rFVIII-SingleChain (lonoctocog alfa, Afstyla) prophylaxis (Bukkems 2023). Weibull base hazard that declines with time since study start (scale 6.37 per year, shape 0.777) with an inhibitory FVIII effect, h(t) = lam * gam * (lam * t)^(gam - 1) * (1 - C / (C + IC50)), IC50 8.90 IU/dL, where lam carries a log-normal individual random effect (IIV 331% CV) and C is chromogenic-assay FVIII activity including the endogenous level. Successive bleeds are independent (no Markov term). FVIII activity is driven by the embedded two-compartment PK model of the same paper (weight on CL and V1, VWF antigen on CL, estimated endogenous FVIII). The cumulative hazard is an ODE state and the survivor function is exposed. Sister PopPK model: modellib('Bukkems_2023_lonoctocogAlfa')."
  reference <- paste(
    "Bukkems LH, Jonsson S, Cnossen MH, Karlsson MO, Mathot RAA; OPTI-CLOT",
    "studies and SYMPHONY consortium. Relationship between factor VIII levels",
    "and bleeding for rFVIII-SingleChain in severe hemophilia A: A repeated",
    "time-to-event analysis. CPT Pharmacometrics Syst Pharmacol.",
    "2023;12(5):706-718. doi:10.1002/psp4.12938. PMID:36965157.",
    "PMCID:PMC10196421. Open Access under CC BY-NC 4.0.",
    "RTTE estimates are in Table 3 and the hazard in Equation 1; the NONMEM",
    "control stream of the RTTE model is Appendix S3 of the supplement, and",
    "the embedded PK parameters are the Table 2 / Appendix S2 final estimates.",
    "Sister PopPK model from the same paper:",
    "modellib('Bukkems_2023_lonoctocogAlfa').",
    sep = " "
  )
  vignette <- "Bukkems_2023_lonoctocogAlfa"
  units <- list(time = "h", dosing = "IU", concentration = "IU/dL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Acts on the EMBEDDED PK model only (power exponents 0.702 on CL and 0.834 on V1, centred at 68 kg), and so on the hazard through FVIII activity. Body weight was also screened as a direct covariate on the bleeding hazard in the full fixed-effects covariate model and had no clinically important effect (Figure 3); it is not in the final RTTE model (Table 3; Appendix S3).",
      source_name = "WT"
    ),
    VWF = list(
      description = "Von Willebrand factor antigen (VWF:Ag)",
      units = "IU/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Acts on the EMBEDDED PK model only (power exponent -0.652 on CL, centred at 113 IU/dL). Screened as a direct covariate on the bleeding hazard with no clinically important effect (Figure 3); not in the final RTTE model. For the source's age-based imputation of missing VWF, see modellib('Bukkems_2023_lonoctocogAlfa').",
      source_name = "VWF"
    )
  )

  # Screened on the bleeding hazard in the full fixed-effects covariate model
  # (Figure 3), which reports only graphical effect sizes; none is in the
  # final RTTE model (Table 3; Appendix S3).
  covariatesDataExcluded <- list(
    ABR_PRESTUDY = list(
      description = "Annual bleeding rate in the period before study entry",
      units = "bleeds/year",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear model centred at the median of 8 bleeds/year, with separate slopes for values <= 100 and > 100 (n = 11 above 100, regarded as clinically improbable; Figure S6). A higher pre-study rate (<= 100) gave a clinically important (> 20%) increase in the bleeding hazard (Figure 3; Abstract 'the pre-study annual bleeding rate explained part of this variability'), but no coefficient is reported and the effect is not in the final model. Cohort median 8 (range 0-168; Table 1)."
    ),
    ARTHROPATHY = list(
      description = "Presence of hemophilic arthropathy at baseline: 1 = present, 0 = absent",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no arthropathy; the Figure 3 reference subject)",
      notes = "No clinically important effect on the bleeding hazard (Figure 3); correlated with the pre-study annual bleeding rate (Figure S7). Present in 112 of 241 subjects (46%; Table 1)."
    ),
    REGION_JAPAN = list(
      description = "Treated in Japan: 1 = yes, 0 = no",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Europe is the Figure 3 reference region)",
      notes = "No clinically important effect on the bleeding hazard (Figure 3). 9 subjects (4%; Table 1)."
    ),
    REGION_USA = list(
      description = "Treated in the United States: 1 = yes, 0 = no",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Europe is the Figure 3 reference region)",
      notes = "No clinically important effect on the bleeding hazard (Figure 3). 22 subjects (9%; Table 1)."
    ),
    REGION_ROW = list(
      description = "Treated in the rest of the world (Australia, Canada, Lebanon, Malaysia, Philippines, Russian Federation, South Africa, Ukraine): 1 = yes, 0 = no",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Europe is the Figure 3 reference region)",
      notes = "No clinically important effect on the bleeding hazard (Figure 3). 75 subjects (31%; Table 1)."
    )
  )

  compartmentData <- list(
    central = list(analyte = "lonoctocog alfa (FVIII activity)", units = "IU", specimen = "plasma", verified = TRUE),
    peripheral1 = list(
      analyte = "lonoctocog alfa (FVIII activity)",
      units = "IU",
      specimen = "plasma",
      verified = TRUE
    ),
    cumhaz = list(
      analyte = "Cumulative hazard of a bleed",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 241L,
    n_studies = 3L,
    age_range = "1-58 years",
    age_median = "19 years",
    weight_range = "10-120 kg",
    weight_median = "60 kg",
    sex_female_pct = NULL,
    race_ethnicity = NULL,
    disease_state = "Severe hemophilia A (endogenous FVIII < 1 IU/dL) on rFVIII-SingleChain prophylaxis; persons treated on demand were excluded.",
    dose_range = "Intravenous rFVIII-SingleChain prophylaxis, 15-50 IU/kg every second day or two to three times per week (median weekly dose 91.6 IU/kg, range 19.6-246).",
    regions = "Europe 56%, Japan 4%, United States 9%, rest of the world 31%.",
    biomarkers = "2080 bleeds (65% joint, 43% spontaneous, 44% trauma-related, 91% treated); median 4 bleeds per subject (range 0-153); 40 subjects (17%) had no bleed. Median annualized bleeding rate during the study 1.6 (range 0-48.5); pre-study median 8 (range 0-168).",
    notes = "Pooled from CSL627_1001, CSL627_3002 and the CSL627_3001 extension (Table 1); median follow-up 1090 days (range 29-1910). Sequential PK-RTTE analysis: individual PK parameters from the PopPK model drove the hazard. Bleeds with only the day known were interval-censored over that day; follow-up was right-censored. Estimated with IMPMAP and the Laplacian likelihood in NONMEM 7.4.1."
  )

  ini({
    # ------------------------------------------------------------------
    # Embedded PK -- identical to modellib('Bukkems_2023_lonoctocogAlfa'),
    # CSA scale. Appendix S3 drives the hazard with A(1)/V1 + BASE from the
    # individual PK parameters; the residual-previous-dose term and the PK
    # residual error are not part of the RTTE model.
    # ------------------------------------------------------------------
    lcl <- log(2.07866); label("Clearance at 68 kg and VWF 113 IU/dL (dL/h)") # Appendix S2 THETA(3) = 2.07866; Table 2 CL 2.08 dL/h
    lvc <- log(34.3587); label("Central volume of distribution at 68 kg (dL)") # Appendix S2 THETA(1) = 34.3587; Table 2 V1 34.4 dL
    lq <- log(1.52841); label("Intercompartmental clearance (dL/h)") # Appendix S2 THETA(4) = 1.52841; Table 2 Q 1.53 dL/h
    lvp <- log(2.40027); label("Peripheral volume of distribution (dL)") # Appendix S2 THETA(2) = 2.40027; Table 2 V2 2.40 dL
    lrbase <- log(0.265462); label("Endogenous FVIII activity (IU/dL)") # Appendix S2 THETA(7) = 0.265462; Table 2 0.27 IU/dL
    e_wt_vc <- 0.833885; label("Power exponent of body weight on V1, centred at 68 kg (unitless)") # Appendix S2 THETA(5) = 0.833885; Table 2 0.84
    e_wt_cl <- 0.702087; label("Power exponent of body weight on CL, centred at 68 kg (unitless)") # Appendix S2 THETA(6) = 0.702087; Table 2 0.70
    e_vwf_cl <- -0.652459; label("Power exponent of VWF antigen on CL, centred at 113 IU/dL (unitless)") # Appendix S2 THETA(10) = -0.652459; Table 2 -0.65

    etalvc ~ 0.0238443 # Appendix S2 $OMEGA 0.0238443 (IIV V1); Table 2 15.5 CV%
    etalcl + etalrbase ~ c(
      0.0611255,
      0.0427114, 1.31873
    ) # Appendix S2 $OMEGA BLOCK(2): IIV CL 0.0611255, covariance 0.0427114, IIV BASE 1.31873

    # ------------------------------------------------------------------
    # RTTE -- Appendix S3 $THETA / $OMEGA final estimates; Table 3 prints
    # them rounded. Appendix S3 converts the scale from per year to per hour
    # as THETA(1)/24/365.
    # ------------------------------------------------------------------
    llam_haz <- log(6.36754 / (24 * 365)); label("Scale of the Weibull bleeding hazard at zero FVIII, lambda (1/h)") # Appendix S3 THETA(1) = 6.36754 per year; Table 3 'Scale' 6.37 year^-1 (RSE 14.3%)
    lgam_haz <- log(0.776767); label("Shape of the Weibull bleeding hazard, gamma (unitless)") # Appendix S3 THETA(3) = 0.776767; Table 3 'Shape' 0.78 (RSE 2.2%)
    lic50 <- log(8.90353); label("FVIII activity halving the bleeding hazard (IU/dL)") # Appendix S3 THETA(2) = 8.90353; Table 3 'IC50' 8.90 IU/dL (RSE 11.7%)

    # The random effect sits on lambda (Appendix S3: SCALE = PSCALE *
    # EXP(ETA(1))), not on the whole hazard as printed in Equation 1.
    etallam_haz ~ 2.47759 # Appendix S3 $OMEGA 2.47759 (IIV_HAZ); Table 3 330.8 CV% = sqrt(exp(2.47759) - 1)

    # The source fits the RTTE likelihood directly, so there is no
    # observation-error model to translate. This placeholder additive
    # residual on the survivor function lets nlmixr2 accept the model for
    # forward simulation; it is not from the source.
    addSd <- 0.001; label("Placeholder additive residual error on the survivor-function output `sur` (unitless); not from the source")
  })

  model({
    # --- embedded PK (identical to modellib('Bukkems_2023_lonoctocogAlfa')) ---
    cl <- exp(lcl + etalcl) * (WT / 68)^e_wt_cl * (VWF / 113)^e_vwf_cl
    vc <- exp(lvc + etalvc) * (WT / 68)^e_wt_vc
    q <- exp(lq)
    vp <- exp(lvp)
    rbase <- exp(lrbase + etalrbase)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Chromogenic-assay FVIII activity including the endogenous level
    # (Appendix S3 $DES: FVIII = A(1)/V1 + BASE).
    Cc <- central / vc + rbase

    # --- bleeding hazard ------------------------------------------------
    lam <- exp(llam_haz + etallam_haz)
    gam <- exp(lgam_haz)
    ic50 <- exp(lic50)

    # Time is time since the subject's first record (study start); it is not
    # reset after a bleed. The 1e-6 offset is the source's DEL, which keeps
    # (lam * t)^(gam - 1) finite at t = 0 for gam < 1.
    hazard <- lam * gam * (lam * t + 1e-6)^(gam - 1) * (1 - Cc / (Cc + ic50))

    d/dt(cumhaz) <- hazard
    cumhaz(0) <- 0

    # Probability of no bleed since the start of the integration window.
    sur <- exp(-cumhaz)

    sur ~ add(addSd)
  })
}
