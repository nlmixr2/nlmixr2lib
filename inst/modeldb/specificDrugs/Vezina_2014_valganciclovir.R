Vezina_2014_valganciclovir <- function() {
  description <- "Two-compartment population PK model for ganciclovir after oral valganciclovir prophylaxis in paediatric and adult solid organ transplant recipients (Vezina 2014). First-order absorption with fixed lag time and rate, allometric (WT/70 kg) scaling on apparent CL/F and Q/F (exponent 0.75) and on V2/F and V3/F (exponent 1.0), and a power-form effect of body-weight-adjusted creatinine clearance on CL/F (reference 60 mL/min)."
  reference   <- "Vezina HE, Brundage RC, Balfour HH Jr. Population pharmacokinetics of valganciclovir prophylaxis in paediatric and adult solid organ transplant recipients. Br J Clin Pharmacol. 2014;78(2):343-352. doi:10.1111/bcp.12343"
  vignette    <- "Vezina_2014_valganciclovir"
  units       <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight (most-recent value at the sampling occasion; piecewise-constant within subject between blood draws, with linear interpolation when no same-day measurement was recorded).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Used for allometric scaling on CL/F and Q/F (exponent 0.75) and on V2/F and V3/F (exponent 1.0), with reference 70 kg fixed (Vezina 2014 Methods 'Population pharmacokinetic analysis'). Study median 71.6 kg in adults (range 8.05-115) and 33.0 kg in children (range 6.9-61.1) per Table 1.",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Creatinine clearance, in mL/min, NOT BSA-normalized to 1.73 m^2. Computed by Cockcroft-Gault for subjects 18 years of age or older and by the Schwartz equation for subjects under 18 years; the inherent BSA standardization in the Schwartz output was reverse-corrected so that all CRCL values are expressed in mL/min directly comparable to the Cockcroft-Gault output.",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying within subject. Power-form effect on CL/F as ((CRCL / 60))^0.492 with reference 60 mL/min (Vezina 2014 Methods: 'standardized to the approximate population median value (i.e. 60 ml min-1) and adjusted by body weight'). The body-weight adjustment described in the Methods is what produces the (WT/70)^0.75 allometric term on CL/F; it is not a second factor on CRCL itself. The CRCL canonical register entry is BSA-normalized; in this model the per-model unit override 'mL/min' is load-bearing -- to use a dataset whose creatinine clearance is BSA-normalized to 1.73 m^2, divide by 1.73 and multiply by the patient BSA before passing it to this model.",
      source_name        = "CrCL"
    )
  )

  population <- list(
    n_subjects        = 95L,
    n_studies         = 1L,
    n_observations    = 269L,
    age_range_overall = "6 months - 78 years",
    age_range_children   = "6 months - 17 years (median 7 years; n = 13: 3 aged 0-24 months, 4 aged 2-11 years, 6 aged 12-17 years)",
    age_range_adults     = "18 - 78 years (median 53 years; n = 82)",
    weight_range_children = "6.9 - 61.1 kg (median 33.0 kg)",
    weight_range_adults   = "8.05 - 115 kg (median 71.6 kg)",
    sex_female_pct    = 36.8,
    race_ethnicity    = c(`Caucasian/White` = 87.4, `African American or Black` = 7.4, `Asian` = 4.2, `American Indian or Alaska Native` = 1.0),
    disease_state     = "First-time solid organ transplant recipients (kidney 56.8%, liver 25.3%, lung 11.6%, kidney/pancreas 4.2%, pancreas 1.1%, kidney/liver 1.1%) on valganciclovir prophylaxis for cytomegalovirus (CMV) or Epstein-Barr virus (EBV) prevention. Maintenance immunosuppression with tacrolimus, ciclosporin, or sirolimus plus mycophenolate.",
    dose_range        = "Tablet: 900 mg every 24 h, or 450 mg every 12, 24, or 48 h (most subjects). Oral solution: 350, 300, 270, 225, 150, or 75 mg every 24 h (n = 8, primarily children). Doses adjusted by Cockcroft-Gault or Schwartz creatinine clearance per the Valcyte package insert.",
    crcl_range_children = "30.2 - 154 mL/min (median 72.1)",
    crcl_range_adults   = "29 - 108 mL/min (median 60.7)",
    donor_source      = c(Deceased = 52.6, `Living unrelated` = 21.1, `Living related` = 26.3),
    regions           = "Single centre, University of Minnesota Medical Center, Fairview (USA).",
    co_medications    = "All subjects also received mycophenolate as part of maintenance immunosuppression (not tested as a covariate). Induction immunosuppression with thymoglobulin, basiliximab, or methylprednisolone; maintenance with tacrolimus, ciclosporin, or sirolimus.",
    notes             = "Prospective natural-history study enrolling between February 2010 and June 2011. Sparse sampling at approximately weeks 2, 4, 8, and 12 post-transplant (and additionally at months 4, 6, 8, and 12 in subjects on prophylaxis longer than 3 months); samples drawn opportunistically at routine post-transplant clinic visits with self-reported dosing histories. 269 of 333 measurements were retained for the population PK analysis; 64 (19.2%) were excluded as below LOD (n = 21), CWRES > 3 SD (n = 6), or internally inconsistent (n = 37). Methods text reports 82 adults + 13 children = 95 analysed; Table 1 lists 83 adult and 13 children baseline rows (96), so the population figures here use the analysis-set N from the Methods narrative."
  )

  ini({
    # Structural PK parameters -- Vezina 2014 Table 2 final-model estimates.
    # Reference subject is a 70 kg individual with creatinine clearance 60
    # mL/min. Time in hours; apparent clearances (CL/F, Q/F) in L/h; apparent
    # volumes (V2/F, V3/F) in L; absorption-rate constant in 1/h; lag time in
    # h. CL/F refers to apparent oral clearance of ganciclovir based on the
    # valganciclovir dose; the F factor absorbs both absolute oral
    # bioavailability and the molar conversion from valganciclovir to
    # ganciclovir.
    lka   <- fixed(log(3.0))   ; label("Absorption rate constant ka (1/h; FIXED)")                                # Vezina 2014 Table 2: Ka = 3.0 (fixed)
    ltlag <- fixed(log(0.5))   ; label("Absorption lag time (h; FIXED)")                                          # Vezina 2014 Table 2: Lag time = 0.5 (fixed)
    lcl   <- log(14.5)         ; label("Apparent oral clearance CL/F at WT = 70 kg, CRCL = 60 mL/min (L/h)")      # Vezina 2014 Table 2: CL/F = 14.5 L/h
    lvc   <- log(87.5)         ; label("Apparent central volume of distribution V2/F at WT = 70 kg (L)")         # Vezina 2014 Table 2: V2/F = 87.5 L
    lq    <- log(4.80)         ; label("Apparent inter-compartmental clearance Q/F at WT = 70 kg (L/h)")         # Vezina 2014 Table 2: Q/F = 4.80 L/h
    lvp   <- log(42.6)         ; label("Apparent peripheral volume of distribution V3/F at WT = 70 kg (L)")      # Vezina 2014 Table 2: V3/F = 42.6 L

    # Allometric exponents -- Vezina 2014 Methods 'Population pharmacokinetic
    # analysis': "Weight was standardized to 70 kg and its effect on CL/F,
    # V2/F, Q/F and V3/F was described by a fixed exponent power function
    # using standard allometric scaling values of 0.75 and 1.0 on clearance
    # and volume of distribution terms, respectively."
    e_wt_cl_q  <- fixed(0.75) ; label("Shared allometric exponent of (WT/70 kg) on CL/F and Q/F (unitless; FIXED)")  # Vezina 2014 Methods: 0.75 on clearance terms
    e_wt_vc_vp <- fixed(1.0)  ; label("Shared allometric exponent of (WT/70 kg) on V2/F and V3/F (unitless; FIXED)") # Vezina 2014 Methods: 1.0 on volume terms

    # Covariate effect on CL/F -- Vezina 2014 Results final-model equation:
    # CL/F (L/h) = 14.5 * ((CRCL/60))^0.492 * (WT/70)^0.75. The CRCL
    # exponent is the only freely estimated covariate exponent in the final
    # model.
    e_crcl_cl <- 0.492 ; label("Power exponent of (CRCL / 60 mL/min) on CL/F (unitless)")  # Vezina 2014 Table 2: CRCL exponent = 0.492

    # Inter-individual variability -- Vezina 2014 Table 2 reports IIV on
    # CL/F as a CV%. Convert to log-scale variance via omega^2 = log(CV^2 + 1):
    #   CL/F  CV 33.5% -> log(0.335^2 + 1) = 0.1063
    # The paper's Methods section explicitly states: "The data were not
    # sufficient to support estimates of interindividual variability on V2/F,
    # Q/F and V3/F."
    etalcl ~ 0.1063  # Vezina 2014 Table 2 IIV CL/F = 33.5%

    # Residual unexplained variability -- Vezina 2014 Table 2 reports a
    # proportional error model. RUV CV = 32.7% -> propSd = 0.327 on the
    # linear concentration scale.
    propSd <- 0.327 ; label("Proportional residual error (fraction)")  # Vezina 2014 Table 2: Residual variability = 32.7%
  })

  model({
    # Individual PK parameters with Vezina 2014 covariate equations.
    # Reference subject: WT 70 kg, CRCL 60 mL/min.
    cl <- exp(lcl + etalcl) * (CRCL / 60)^e_crcl_cl * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc)          * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp)          * (WT / 70)^e_wt_vc_vp

    ka   <- exp(lka)
    tlag <- exp(ltlag)

    # Two-compartment oral PK with first-order absorption + lag time.
    # NONMEM ADVAN4/TRANS4 equivalent. Dose is valganciclovir mg landing in
    # `depot`; ganciclovir is the measured analyte and the apparent CL/F,
    # V2/F, Q/F, V3/F absorb the F factor (oral bioavailability x molar
    # conversion valganciclovir -> ganciclovir).
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # The ganciclovir assay reports concentrations in ng/mL. Internal mass /
    # volume is mg / L; multiply by 1000 to match the assay (and the units
    # field declared above).
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
