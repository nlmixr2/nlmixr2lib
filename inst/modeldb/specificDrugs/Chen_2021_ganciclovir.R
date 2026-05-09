Chen_2021_ganciclovir <- function() {
  description <- "Two-compartment population PK model for oral ganciclovir (the active metabolite of valganciclovir) in adult Chinese renal allograft recipients (Chen 2021), with first-order absorption after a lag time and a linear creatinine-clearance effect on apparent oral clearance (CL/F)."
  reference <- "Chen B, Hu SS, Rui WB, An HM, Zhai XH, Wang XH, Lu JQ, Shao K, Zhou PJ. Population Pharmacokinetics and Bayesian Estimation of the Area Under the Concentration-Time Curve for Ganciclovir in Adult Chinese Renal Allograft Recipients After Valganciclovir Administration. J Clin Pharmacol. 2021;61(3):328-338. doi:10.1002/jcph.1735"
  vignette <- "Chen_2021_ganciclovir"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source paper computes CrCl with the Cockcroft-Gault equation in raw mL/min and does NOT divide by BSA. Dataset values must therefore be raw Cockcroft-Gault CrCl, not the BSA-normalized form normally implied by the canonical CRCL register entry. Linear-deviation effect on CL/F as `(1 + e_crcl_cl * CRCL / 68.3)`; the centring value 68.3 mL/min is the value the source equation uses (close to but not equal to the modeling-group median CrCl of 64.5 mL/min reported in Table 1).",
      source_name        = "CLcr"
    )
  )

  population <- list(
    n_subjects        = 70L,
    n_studies         = 1L,
    age_range         = "17-61 years",
    age_mean          = "42.3 years (SD 9.95)",
    weight_range      = "40-85 kg",
    weight_mean       = "61.1 kg (SD 11.0)",
    sex_female_pct    = 34.3,
    race_ethnicity    = "Chinese (Han presumed; not specified in source).",
    disease_state     = "Adult kidney transplant recipients on triple immunosuppression (cyclosporin or tacrolimus + mycophenolate mofetil or sodium + prednisone) with valganciclovir CMV prophylaxis started 3 weeks post-transplant.",
    dose_range        = "Oral valganciclovir 450 mg or 900 mg once daily for 5-7 days (steady state).",
    regions           = "China (Ruijin Hospital, Shanghai Jiaotong University School of Medicine).",
    crcl_range        = "32.4-118.4 mL/min (modeling group, range); median 64.5 mL/min.",
    notes             = "Demographics from Chen 2021 Table 1. 70 patients split into a modeling group (n = 40, 23 received 450 mg, 17 received 900 mg) and a validation group (n = 30, 18 received 450 mg, 12 received 900 mg). Plasma ganciclovir was measured at predose and 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, and 24 h after dosing on day 5-7. Induction therapy was rabbit anti-thymocyte globulin (ATG, n = 39) or basiliximab (Simulect, n = 25); not modeled as a covariate."
  )

  ini({
    # Structural PK -- Chen 2021 Table 3 'Final' column (n = 70). Time in hours,
    # apparent clearances (CL/F, Q/F) in L/h, apparent volumes (V2/F, V3/F) in L,
    # absorption rate constant ka in 1/h, lag time Tlag in hours. The Table 3
    # reference subject for CL/F is a patient with CrCl = 0 mL/min (the
    # CL/F covariate equation is `theta1 * (1 + theta_CLcr * CrCl / 68.3)`, so
    # exp(lcl) = theta1 = 7.09 L/h is the apparent clearance at CrCl = 0).
    lcl   <- log(7.09)  ; label("Apparent oral clearance CL/F at CrCl = 0 mL/min (L/h)") # Chen 2021 Table 3 Final theta1 = 7.09 L/h
    lvc   <- log(10.8)  ; label("Apparent central volume V2/F (L)")                       # Chen 2021 Table 3 Final theta2 = 10.8 L
    lq    <- log(3.96)  ; label("Apparent inter-compartmental clearance Q/F (L/h)")       # Chen 2021 Table 3 Final theta3 = 3.96 L/h
    lvp   <- log(174)   ; label("Apparent peripheral volume V3/F (L)")                    # Chen 2021 Table 3 Final theta4 = 174 L
    lka   <- log(0.23)  ; label("First-order absorption rate constant ka (1/h)")          # Chen 2021 Table 3 Final theta5 = 0.23 1/h
    ltlag <- log(0.93)  ; label("Absorption lag time Tlag (h)")                            # Chen 2021 Table 3 Final theta6 = 0.93 h

    # Covariate effect on CL/F -- Chen 2021 Table 3 footnote covariate equation:
    # CL/F_i = theta1 * (1 + CLcr/68.3 * theta_CLcr) * exp(eta_CL). The
    # centring value 68.3 mL/min appears in the published equation; raw
    # Cockcroft-Gault CrCl in mL/min (NOT BSA-normalized).
    e_crcl_cl <- 1.08 ; label("Cockcroft-Gault CrCl linear coefficient on CL/F via (1 + e_crcl_cl * CRCL / 68.3) (unitless)") # Chen 2021 Table 3 Final theta7 = 1.08

    # Inter-individual variability -- Chen 2021 Table 3 Final 'omega (%)' is the
    # standard deviation of the log-scale random effect expressed as a
    # percentage, i.e. omega = 0.272 for CL means SD(eta_CL) = 0.272 (variance
    # 0.0740). Values for V2/F (153%), Q/F (63.1%), V3/F (107%) are large but
    # reproduce the published estimates verbatim.
    etalcl ~ 0.0740                                                                       # Chen 2021 Table 3 Final omega(CL/F)  = 27.2% -> 0.272^2
    etalvc ~ 2.3409                                                                       # Chen 2021 Table 3 Final omega(V2/F)  = 153%  -> 1.53^2
    etalq  ~ 0.3982                                                                       # Chen 2021 Table 3 Final omega(Q/F)   = 63.1% -> 0.631^2
    etalvp ~ 1.1449                                                                       # Chen 2021 Table 3 Final omega(V3/F)  = 107%  -> 1.07^2

    # Residual unexplained variability -- Chen 2021 Methods: log-transformed
    # observation, ln(Cobs) = ln(Cpred) + epsilon with var(epsilon) = sigma^2.
    # Per the NONMEM-to-nlmixr2 mapping, "additive on the log scale" is
    # proportional in linear space. Table 3 reports 'Residual variance delta =
    # 42.9' alongside %-formatted IIV; we read delta as a percentage (i.e. the
    # log-scale residual SD x 100, sigma = 0.429), consistent with the
    # IIV-row convention. The interpretation is documented in the vignette
    # Errata.
    propSd <- 0.429 ; label("Proportional residual error (fraction)")                     # Chen 2021 Table 3 Final delta = 42.9
  })

  model({
    # Individual PK parameters with the Chen 2021 covariate equation. The
    # Table 3 footnote gives CL/F as theta1 * (1 + CLcr/68.3 * theta_CLcr); a
    # patient at CLcr = 68.3 mL/min has typical CL/F = 7.09 * (1 + 1.08) =
    # 14.75 L/h, close to the base-model CL/F of 15.8 L/h.
    cl   <- exp(lcl + etalcl) * (1 + e_crcl_cl * CRCL / 68.3)
    vc   <- exp(lvc + etalvc)
    q    <- exp(lq  + etalq)
    vp   <- exp(lvp + etalvp)
    ka   <- exp(lka)
    tlag <- exp(ltlag)

    # Two-compartment oral PK with first-order absorption + lag time.
    # Bioavailability of GCV from oral VGCV is absorbed into the apparent
    # parameters (the model is parameterized in CL/F, V2/F, Q/F, V3/F), so no
    # explicit f(depot) is needed. Dose enters `depot` with a lag of `tlag`
    # hours, then transits to `central` with rate ka.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # Concentration in mg/L (dose mg, vc L). The source paper reports
    # GCV plasma concentrations in mg/L; the unit-rescaling factor is 1.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
