Chien_2022_imatinib <- function() {
  description <- "Two-compartment population PK model for oral imatinib in healthy adult volunteers (Chien 2022); first-order absorption preceded by a Savic 2007-style analytical transit-compartment chain (mean transit time and number of transit compartments estimated), first-order elimination, and an OMEGA BLOCK between the IIV on CL and V1 motivated by their estimated correlation r > 0.9. No covariates were retained in the final model."
  reference <- "Chien YH, Wuerthwein G, Zubiaur P, Posocco B, Pena MA, Borobia AM, Gagno S, Abad-Santos F, Hempel G. Population pharmacokinetic modelling of imatinib in healthy subjects receiving a single dose of 400 mg. Cancer Chemother Pharmacol. 2022;90(2):125-136. doi:10.1007/s00280-022-04454-y"
  vignette <- "Chien_2022_imatinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  population <- list(
    n_subjects     = 26L,
    n_studies      = 2L,
    n_observations = 472L,
    age_range      = "19.7-31.0 years",
    age_median     = "23.0 years",
    weight_range   = "52.0-96.0 kg",
    weight_median  = "69.5 kg",
    bmi_range      = "20.0-30.0 kg/m^2",
    bmi_median     = "22.5 kg/m^2",
    bsa_range      = "1.52-2.22 m^2",
    bsa_median     = "1.86 m^2",
    sex_female_pct = 30.8,
    race_ethnicity = "Healthy Caucasian volunteers (single-ancestry cohort enrolled in Spain)",
    disease_state  = "Healthy adult volunteers (no concurrent drugs, no organ dysfunction or inflammation, BMI 18-30 kg/m^2, all clinical labs within normal range or judged acceptable per investigators)",
    dose_range     = "Oral imatinib 400 mg single dose",
    administration = "Oral (single 400 mg dose; pooled from two randomised crossover bioequivalence studies)",
    regions        = "Spain (Hospital Universitario de La Paz, Madrid; Hospital General de Alicante)",
    notes          = "Demographics summarised in Chien 2022 page 4 paragraph 1 (median 23.0 y, 69.5 kg, 22.5 kg/m^2, 1.86 m^2 BSA) and Supplement Table S1 (not on disk). Sex distribution 8 female / 18 male = 30.8 % female. Genotype distribution for CYP3A4, CYP3A5, CYP2C9, CYP2C19, CYP2C8, CYP2B6, CYP2D6 and ABCB1 reported in Methods 'Sampling and analysis' but not retained as covariates in the final model. Sampling: 16-19 plasma samples per volunteer between 0.5 and 72 h post-dose."
  )

  ini({
    # Final estimates of population PK parameters (Chien 2022 Table 2,
    # 'PopPK model for healthy volunteers' column). Bootstrap medians and
    # 95 % CIs are also given in Table 2 and agree with the point estimates.
    lcl  <- log(13.2);  label("Apparent oral clearance CL/F (L/h)")                       # Chien 2022 Table 2: CL/F = 13.2 L/h (RSE 5.0 %)
    lq   <- log(3.75);  label("Apparent inter-compartmental clearance Q/F (L/h)")          # Chien 2022 Table 2: Q/F = 3.75 L/h (RSE 20.5 %)
    lvc  <- log(172);   label("Apparent volume of central compartment V1/F (L)")           # Chien 2022 Table 2: V1/F = 172 L (RSE 4.7 %)
    lvp  <- log(43.6);  label("Apparent volume of peripheral compartment V2/F (L)")        # Chien 2022 Table 2: V2/F = 43.6 L (RSE 10.4 %)
    lka  <- log(1.22);  label("First-order absorption rate constant Ka (1/h)")             # Chien 2022 Table 2: Ka = 1.22 1/h (RSE 18.2 %)
    lmtt <- log(0.537); label("Mean absorption transit time MTT (h)")                      # Chien 2022 Table 2: MTT = 0.537 h (RSE 14.9 %)
    lnn  <- log(3.62);  label("Number of absorption transit compartments NN (continuous)") # Chien 2022 Table 2: N = 3.62 (RSE 12.1 %)

    # Inter-individual variability. Chien 2022 Table 2 reports IIV as CV %
    # on CL = 24.8 %, V1 = 27.7 %, Ka = 88.3 %, MTT = 80.5 %.
    # Variances on the log scale: omega^2 = log(1 + CV^2) for log-normal.
    #   CL  : log(1 + 0.248^2) = 0.05969
    #   V1  : log(1 + 0.277^2) = 0.07393
    #   Ka  : log(1 + 0.883^2) = 0.57644
    #   MTT : log(1 + 0.805^2) = 0.49958
    # Chien 2022 page 5 paragraph 1 ('Structure model development') reports
    # an OMEGA BLOCK between IIV on CL and V1 (r > 0.9) added at the end
    # of model development (model D1 in Table 1). The numeric correlation
    # is not reported; r = 0.9 is used here as the stated lower bound (a
    # documented approximation; see vignette Errata).
    #   cov(CL, V1) = 0.9 * sqrt(0.05969 * 0.07393) = 0.05978
    # Note: Table 1 row D1's verbal label says 'OMEGA BLOCK between IIV on
    # CL and Ka', which conflicts with the body-text statement 'CL and V1
    # (r > 0.9)'. The body-text wording is more specific and is consistent
    # with the much smaller RSE on CL (12.8 %) and V1 (11.4 %) IIVs in
    # Table 2 vs Ka (15.7 %) and MTT (17.3 %). The body text is followed
    # here; see vignette Errata for the discrepancy note.
    etalcl + etalvc ~ c(
      0.05969,
      0.05978, 0.07393
    )                                                                                       # Chien 2022 Table 2 IIV CL = 24.8 %, V1 = 27.7 %; page 5 paragraph 1 OMEGA BLOCK r > 0.9 (using r = 0.9 as documented approximation)
    etalka  ~ 0.57644                                                                       # Chien 2022 Table 2 IIV Ka = 88.3 %; log(1 + 0.883^2)
    etalmtt ~ 0.49958                                                                       # Chien 2022 Table 2 IIV MTT = 80.5 %; log(1 + 0.805^2)

    # Residual error: proportional only (Chien 2022 page 4 paragraph 1
    # 'A proportional error model was selected' and Table 2 'Prop error
    # 13.6 %').
    propSd <- 0.136; label("Proportional residual error (fraction)")                       # Chien 2022 Table 2: Prop error = 13.6 % (RSE 10.0 %)
  })

  model({
    # Individual PK parameters. No covariates retained in the final model
    # (Chien 2022 page 5 paragraph 2 'the final popPK model is constructed
    # without any covariate').
    cl  <- exp(lcl + etalcl)
    vc  <- exp(lvc + etalvc)
    vp  <- exp(lvp)
    q   <- exp(lq)
    ka  <- exp(lka  + etalka)
    mtt <- exp(lmtt + etalmtt)
    nn  <- exp(lnn)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment disposition with a Savic 2007 analytical transit-
    # compartment chain (rxode2 transit() built-in) feeding the depot,
    # which then absorbs into central at rate ka. The transit() function
    # reads the dose amount from the depot input; f(depot) <- 0 suppresses
    # the bolus contribution so the chain delivers the full dose at the
    # gamma-density input rate.
    d/dt(depot)       <- transit(nn, mtt) - ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- 0  # suppress bolus into depot; transit() drives the input rate

    # Plasma imatinib concentration. Dose in mg, V in L, observed
    # concentrations in ng/mL (Chien 2022 reports Cmax and trough levels
    # in ng/mL and ug/L throughout); central (mg) / vc (L) gives mg/L =
    # ug/mL, multiplied by 1000 to recover ng/mL.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
