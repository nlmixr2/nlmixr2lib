Chelle_2019_factorviii_fanhdi <- function() {
  description <- "Two-compartment population PK model for Fanhdi/Alphanate (plasma-derived factor VIII concentrate, Grifols) in hemophilia A patients pooled from 12 hemophilia centers in the WAPPS-Hemo platform (Chelle 2019). Final model has fat-free mass (FFM) as a power-form covariate on CL, V1, and V2, and a piecewise-linear age effect on CL above the median age of 25 years; between-subject variability is a BLOCK(2) on CL and V1 with correlation 0.797; residual error is proportional only."
  reference <- "Chelle P, Yeung CHT, Bonanad S, Morales Munoz JC, Ozelo MC, Megias Vericat JE, Iorio A, Spears J, Mir R, Edginton A. Routine clinical care data for population pharmacokinetic modeling: the case for Fanhdi/Alphanate in hemophilia A patients. J Pharmacokinet Pharmacodyn. 2019 Oct;46(5):427-438. doi:10.1007/s10928-019-09637-4. PMID:31115793."
  vignette <- "Chelle_2019_factorviii_fanhdi"
  units <- list(time = "h", dosing = "IU", concentration = "IU/mL")

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL, V1, and V2 with reference FFM = 50.5 kg, the derivation-cohort median (Chelle 2019 Eq. 5; Table 1 derivation-population median = 50.5 kg). The paper computes FFM from a dual-energy x-ray-absorptiometry-validated formula that spans ages 3-82 years (Chelle 2019 Discussion; reference 22 in the paper). When body height is missing, the paper imputes HT from a multilinear regression on body weight and age before computing FFM. Treated as time-fixed at baseline.",
      source_name        = "FFM"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Piecewise-linear effect on CL only (Chelle 2019 Eq. 5): for AGE <= 25 years the typical CL is independent of age; for AGE > 25 years the typical CL scales as (1 + e_age_cl * (AGE - 25)/25) with e_age_cl = -0.302, so older patients have lower CL. Reference 25 years is the derivation-cohort median (Chelle 2019 Table 1). The paper used age as a surrogate for von Willebrand factor (vWF) protection of FVIII activity, since vWF data were unavailable for many subjects (Chelle 2019 Discussion).",
      source_name        = "AGE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 92L,
    n_studies      = 1L,
    age_range      = "1-72 years",
    age_median     = "25 years (mean 26.1, SD 18.3)",
    weight_range   = "9.68-119 kg",
    weight_median  = "63.5 kg (mean 59.9, SD 25.9)",
    height_range   = "73.8-188 cm (n = 87 reported; HT imputed for 5 subjects)",
    height_median  = "167 cm (mean 155.4, SD 26.6)",
    bmi_range      = "11.1-39.3 kg/m^2 (n = 87)",
    bmi_median     = "23.9 kg/m^2 (mean 23.4, SD 5.5)",
    ffm_range      = "7.5-73.0 kg (n = 87)",
    ffm_median     = "50.5 kg (mean 45.3, SD 18.0); used as the FFM covariate-centering reference",
    sex_female_pct = 0,
    race_ethnicity = "not reported in source",
    disease_state  = "Hemophilia A; 87.0% (80 / 92) had severe disease with endogenous FVIII activity < 0.01 IU/mL; remaining patients had endogenous FVIII up to 0.169 IU/mL. Patients with current FVIII inhibitors were excluded; patients with a history of inhibitors (now negative) were included.",
    dose_range     = "Single intravenous infusion of Fanhdi or Alphanate; one occasion per subject; 1-8 post-infusion FVIII activity samples per subject (median 5, mean 4.2, SD 1.5; total 386 observations) measured by one-stage clotting assay (LLOQ 0.01 IU/mL; 13 / 386 = 3.4% observations BLQ handled via M3 censoring during estimation)",
    regions        = "12 hemophilia centers worldwide; the three largest (67 / 92 subjects) were Campinas (Brazil), Valencia (Spain), and Santiago (Chile); the remaining 25 subjects from 9 other centers",
    notes          = "Data were extracted from the WAPPS-Hemo (Web-Accessible Population Pharmacokinetic Service - Hemophilia) database on 16 February 2018 under clinicaltrials.gov NCT02061072 / NCT03533504 (McMaster University HIREB). Hemophilia A is X-linked recessive so the cohort is all-male (92 / 92). The paper additionally reports an external-evaluation cohort of 49 patients (Chelle 2019 Table 1, 'Evaluation population'); only the 92-subject derivation cohort is encoded here. The model's intended use is as a prior for Bayesian forecasting on the WAPPS-Hemo platform; for that purpose the paper also published comparative results vs the McEneny-King 2019 generic plasma-derived-FVIII WAPPS model on the same 49-subject external cohort."
  )

  ini({
    # Structural parameters - typical values for the paper's reference patient
    # (FFM = 50.5 kg, AGE = 25 years; derivation-cohort medians). Volumes in L,
    # clearances in L/h match the paper. Source: Chelle 2019 Table 2.
    lcl <- log(0.195);  label("Clearance for the reference 50.5 kg FFM, 25 year-old patient (CL, L/h)") # Chelle 2019 Table 2: CLpop = 0.195 L/h (5.69% RSE; bootstrap 95% CI 0.176-0.217)
    lvc <- log(2.30);   label("Central volume for the reference 50.5 kg FFM patient (V1, L)")          # Chelle 2019 Table 2: V1pop = 2.30 L (7.45% RSE; bootstrap 95% CI 1.95-2.62)
    lq  <- log(0.078);  label("Inter-compartmental clearance (Q, L/h)")                                # Chelle 2019 Table 2: Qpop = 0.078 L/h (21.3% RSE; bootstrap 95% CI 0.047-0.120)
    lvp <- log(0.449);  label("Peripheral volume for the reference 50.5 kg FFM patient (V2, L)")       # Chelle 2019 Table 2: V2pop = 0.449 L (27.1% RSE; bootstrap 95% CI 0.279-0.776)

    # Covariate effects: power form for FFM on CL, V1, V2 (Chelle 2019 Eq. 5
    # and Eq. 2: TVPi = Ppop * (covi / covmed)^theta_cov); linear piecewise
    # form for AGE on CL above the median age of 25 years.
    e_ffm_cl <-  0.701; label("Power exponent of FFM on CL (unitless)")  # Chelle 2019 Table 2: FFM effect on CL = 0.701 (12.0% RSE; 95% CI 0.527-0.872)
    e_ffm_vc <-  0.726; label("Power exponent of FFM on V1 (unitless)")  # Chelle 2019 Table 2: FFM effect on V1 = 0.726 (13.0% RSE; 95% CI 0.542-0.903)
    e_ffm_vp <-  0.842; label("Power exponent of FFM on V2 (unitless)")  # Chelle 2019 Table 2: FFM effect on V2 = 0.842 (72.7% RSE; 95% CI 0.365-3.976; the wide CI reflects sparse posterior identification, kept on physiological grounds per Discussion)
    e_age_cl <- -0.302; label("Linear coefficient of (AGE - 25)/25 on CL above age 25 (unitless)") # Chelle 2019 Table 2: AGE effect on CL = -0.302 (19.1% RSE; 95% CI -0.407 to -0.167)

    # Inter-individual variability: BLOCK(2) on CL and V1 (Chelle 2019 Eq. 5
    # and Table 2). The paper reports BSV as "CV: coefficient of variation
    # (defined as standard deviation of eta)", i.e. omega = SD on the eta
    # scale. Therefore:
    #   omega^2(CL)         = 0.456^2 = 0.207936
    #   omega^2(V1)         = 0.542^2 = 0.293764
    #   cov(etalcl, etalvc) = 0.797 * 0.456 * 0.542 = 0.196980
    # BSV on Q and V2 was tested but rejected due to shrinkage > 44%
    # (Chelle 2019 Results, "Development of the PopPK model"). The paper
    # reports no IOV.
    etalcl + etalvc ~ c(0.207936,
                        0.196980,  0.293764) # Chelle 2019 Table 2: CV(CL) = 0.456 (9.22% RSE), Corr(CL, V1) = 0.797 (7.50% RSE), CV(V1) = 0.542 (11.3% RSE)

    # Residual error: proportional only. The paper tested additive,
    # proportional, and combined errors; addition of any additive component
    # did not significantly decrease the OFV (Chelle 2019 Results,
    # "Development of the PopPK model"). The "CV of proportional RUV" in
    # Table 2 is the SD of the proportional residual (eps_prop) so that
    # observation = prediction * (1 + eps_prop).
    propSd <- 0.205; label("Proportional residual error (fraction)") # Chelle 2019 Table 2: CV of proportional RUV = 0.205 (8.23% RSE; 95% CI 0.169-0.232)
  })
  model({
    # Piecewise-linear AGE term on CL: (1 + e_age_cl * max(0, AGE - 25)/25).
    # The factor equals 1 for AGE <= 25 (no effect) and (1 + e_age_cl *
    # (AGE - 25)/25) for AGE > 25. Negative e_age_cl drives CL down with age.
    age_dev_pos <- (AGE - 25) * (AGE > 25)
    age_eff_cl  <- 1 + e_age_cl * age_dev_pos / 25

    # Individual PK parameters with covariate effects (Chelle 2019 Eq. 5).
    cl <- exp(lcl + etalcl) * (FFM / 50.5)^e_ffm_cl * age_eff_cl
    vc <- exp(lvc + etalvc) * (FFM / 50.5)^e_ffm_vc
    q  <- exp(lq)
    vp <- exp(lvp) * (FFM / 50.5)^e_ffm_vp

    # Micro-constants for the explicit two-compartment ODEs.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment IV model: Fanhdi/Alphanate is administered as a single
    # IV infusion; doses enter the central compartment directly. Infusion
    # duration is supplied through the user's event table (rate / dur), not
    # the structural model.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Observation: FVIII activity in IU/mL. Dose in IU and V1 in L give
    # central / vc in IU/L; divide by 1000 to convert to IU/mL, the unit the
    # paper uses throughout (LLOQ 0.01 IU/mL; severe-disease threshold 0.01
    # IU/mL; baseline endogenous 0.005 IU/mL).
    Cc <- central / vc / 1000
    Cc ~ prop(propSd)
  })
}
