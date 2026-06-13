Wilhelm_2012_ciclosporin <- function() {
  description <- "Two-compartment population PK model for ciclosporin (CsA) in adults undergoing haematopoietic allogeneic stem cell transplantation, with first-order oral absorption + lag time and a 3 h intravenous infusion directly into the central compartment (Wilhelm 2012). Twenty subjects on routine fluconazole antimycotic prophylaxis (a CYP3A4 inhibitor) were included; ciclosporin was assayed in whole blood by FPIA (AxSYM, Abbott). Body weight, body surface area, co-medication with CYP3A4 inducers and co-medication with CYP3A4 inhibitors were tested but none reached statistical or clinical significance, so no covariates are retained in the final model. Inter-individual variability was reported on every PK parameter (CL, Vc, Q, Vp, ka, F, tlag); the paper estimated a full omega variance-covariance matrix but did not publish the off-diagonal elements, so the packaged model uses diagonal IIVs only (see vignette Assumptions and deviations)."
  reference <- paste(
    "Wilhelm AJ, de Graaf P, Veldkamp AI, Janssen JJWM, Huijgens PC, Swart EL.",
    "Population pharmacokinetics of ciclosporin in haematopoietic allogeneic",
    "stem cell transplantation with emphasis on limited sampling strategy.",
    "Br J Clin Pharmacol. 2012;73(4):553-563.",
    "doi:10.1111/j.1365-2125.2011.04116.x.",
    sep = " "
  )
  vignette <- "Wilhelm_2012_ciclosporin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight (kg). Tested but not retained.",
      units       = "kg",
      type        = "continuous",
      notes       = "Wilhelm 2012 Results / Covariate model building: weight (median 84 kg, range 53-110 kg) was entered into the basic model individually with CL and V1 as candidate parameters; the effect did not meet the statistical (delta-OFV > 6.6) and clinical (>= 20% change in typical value across observed range) inclusion criteria. Reported in Results: 'None of the covariates investigated for their influence on CL and V1 (WT, BSA, INH, IND) proved to be significant.'",
      source_name = "WT"
    ),
    BSA = list(
      description = "Body surface area (m^2). Tested but not retained.",
      units       = "m^2",
      type        = "continuous",
      notes       = "Wilhelm 2012 Results: tested on CL and V1, did not meet inclusion criteria. Cohort BSA median 2.02 m^2 (range 1.48-2.43, Table 1).",
      source_name = "BSA"
    ),
    INH = list(
      description = "Co-medication with CYP3A4 inhibitors (binary). Tested but not retained.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Wilhelm 2012 Methods/Results: dichotomous indicator for cotreatment with CYP3A4 inhibitors; all 20 subjects received fluconazole 50 mg once daily as antimycotic prophylaxis so within-cohort variability on this covariate was negligible. Did not meet inclusion criteria.",
      source_name = "INH"
    ),
    IND = list(
      description = "Co-medication with CYP3A4 inducers (binary). Tested but not retained.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Wilhelm 2012 Methods/Results: dichotomous indicator for cotreatment with CYP3A4 inducers; tested on CL via the multiplicative form CL = theta1 * (1 + theta2 * IND). Did not meet inclusion criteria.",
      source_name = "IND"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 20L,
    n_studies      = 1L,
    n_observations = "436 whole-blood ciclosporin concentrations across the 20 subjects following IV and oral CsA administration",
    age_range      = "37-66 years",
    age_median     = "54 years",
    weight_range   = "53-110 kg",
    weight_median  = "84 kg",
    bsa_range      = "1.48-2.43 m^2",
    bsa_median     = "2.02 m^2",
    sex_female_pct = 35.0,
    race_ethnicity = "Not reported in source paper (single-centre Dutch cohort).",
    disease_state  = "Adults undergoing allogeneic haematopoietic stem cell transplantation (HSCT) for haematological malignancies (acute myeloid leukaemia n=7, non-Hodgkin lymphoma n=6, chronic lymphoblastic leukaemia n=2, other n=5). Conditioning regimens were fludarabine/cyclophosphamide (n=10), fludarabine/total body irradiation (n=6), or cyclophosphamide/total body irradiation (n=4) per Table 1.",
    dose_range     = "2.5 mg/kg IV by 3 h infusion at study start; subsequently BID with dose adjusted to target trough 200-400 ug/L; oral phase used Neoral microemulsion at a mean of 2.77 +/- 0.81 mg/kg during the recorded oral profile.",
    co_medication  = "All subjects received fluconazole 50 mg once daily (CYP3A4 inhibitor) as antimycotic prophylaxis. Other routine prophylaxis: ciprofloxacin 500 mg BID, valacyclovir 500 mg BID, phenethicillin 250 mg QID, co-trimoxazole 960 mg BID twice weekly.",
    regions        = "The Netherlands (single centre: VU University Medical Center, Amsterdam).",
    sampling_window = "IV profile: predose and 0.25, 0.5, 1, 1.5, 2.5, 4, 4.5, 6.5, 9, 10, 12 h after start of infusion. Oral profile (post-tolerance): predose and 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 5, 8, 12 h after the morning dose.",
    assay          = "Validated specific fluorescence polarization immunoassay (FPIA) on AxSYM (Abbott Diagnostics). LLOQ 80 ug/L; intra-assay CV < 10%; accuracy 96% at 400 ug/L. Low cross-reactivity with CsA metabolites (AM1 5.5%, AM9 13.7%, AM4n 2.1%, AM19 2.5%).",
    notes          = "Study period January 2005 - February 2008. Approved by VU University Medical Center Ethics Committee, written informed consent obtained. Demographic details from Table 1; sampling schedule from Methods 'Study design'."
  )

  ini({
    # Final-model parameter estimates from Table 3 of Wilhelm 2012
    # ('Final parameters estimates of the pharmacokinetic model of
    # ciclosporin'). Inter-individual variability is reported as
    # percent CV; the internal log-scale variance is omega^2 =
    # log(1 + CV^2) for log-normal-distributed parameters (Methods
    # 'Basic pharmacokinetic model': 'The distribution of individual
    # clearance (CL), absorption rate constant (ka), oral
    # bioavailability (F), volumes of distribution (V1 central; V2
    # peripheral) and intercompartmental clearance (Q) was assumed to
    # be lognormal').

    # Structural PK parameters
    lcl   <- log(21.9);   label("Clearance CL (L/h)")                                           # Table 3: CL = 21.9 L/h (RSD 5.2%)
    lvc   <- log(16.6);   label("Central volume of distribution V1 (L)")                        # Table 3: V1 = 16.6 L (RSD 8.7%) -- the abstract reports 18.3 L (= bootstrap median ~ 18) but the formal final estimate from the table is 16.6 L
    lq    <- log(24.2);   label("Inter-compartmental clearance Q (L/h)")                        # Table 3: Q  = 24.2 L/h (RSD 9.3%)
    lvp   <- log(59.0);   label("Peripheral volume of distribution V2 (L)")                     # Table 3: V2 = 59.0 L (RSD 8.8%)

    # Absorption (oral route only; IV is a 3 h infusion direct to central)
    lka     <- log(0.280); label("First-order absorption rate constant ka (1/h)")               # Table 3: ka   = 0.280 1/h (RSD 14.6%)
    lfdepot <- log(0.710); label("Oral bioavailability F (fraction)")                            # Table 3: F    = 0.710 (RSD 9.9%)
    ltlag   <- log(0.440); label("Absorption lag time tlag (h)")                                 # Table 3: tlag = 0.440 h (RSD 5.5%)

    # Inter-individual variability (log-normal; omega^2 = log(1 + CV^2)).
    # Methods note the source estimated a full omega variance-covariance
    # matrix, but Table 3 reports only the diagonal CVs; off-diagonal
    # elements are not published, so the packaged model uses diagonal
    # IIVs only (see vignette Assumptions and deviations).
    etalcl     ~ 0.04812  # Table 3: IIV CL   = 22.2% CV -> log(1 + 0.222^2)
    etalvc     ~ 0.06988  # Table 3: IIV V1   = 26.9% CV -> log(1 + 0.269^2)
    etalq      ~ 0.07650  # Table 3: IIV Q    = 28.2% CV -> log(1 + 0.282^2)
    etalvp     ~ 0.08947  # Table 3: IIV V2   = 30.6% CV -> log(1 + 0.306^2)
    etalka     ~ 0.17559  # Table 3: IIV ka   = 43.8% CV -> log(1 + 0.438^2)
    etalfdepot ~ 0.06062  # Table 3: IIV F    = 25.0% CV -> log(1 + 0.250^2)
    etaltlag   ~ 0.03224  # Table 3: IIV tlag = 18.1% CV -> log(1 + 0.181^2)

    # Residual error: combined additive + proportional (Methods 'Basic
    # pharmacokinetic model'; Results 'combined additional and
    # proportional error of 65 ug/L and 9%').
    addSd  <- 65;    label("Additive residual SD (ug/L)")                                       # Table 3: additive error = 65 ug/L (RSD 86%)
    propSd <- 0.088; label("Proportional residual SD (fraction)")                                # Table 3: proportional error = 8.8% (RSD 84%)
  })

  model({
    # Individual PK parameters
    cl   <- exp(lcl   + etalcl)
    vc   <- exp(lvc   + etalvc)
    q    <- exp(lq    + etalq)
    vp   <- exp(lvp   + etalvp)
    ka   <- exp(lka   + etalka)
    bio  <- exp(lfdepot + etalfdepot)
    tlag <- exp(ltlag + etaltlag)

    # Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment disposition with first-order absorption from depot.
    # IV doses target the central compartment (3 h infusion in the source
    # protocol) and bypass both the depot and the bioavailability factor.
    # Oral doses target the depot, with bioavailability bio = F applied
    # via f(depot) and absorption lag via alag(depot).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot)    <- bio
    alag(depot) <- tlag

    # Whole-blood ciclosporin concentration (ug/L). Dose in mg, Vc in L,
    # so central (mg) / vc (L) -> mg/L; multiply by 1000 to obtain ug/L
    # to match the FPIA assay units used throughout the paper.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
