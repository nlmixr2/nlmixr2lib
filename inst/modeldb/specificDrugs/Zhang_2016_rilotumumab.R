Zhang_2016_rilotumumab <- function() {
  description <- "Two-compartment IV population PK model for rilotumumab (fully human anti-HGF IgG2 monoclonal antibody) in patients with MET-positive gastric or gastroesophageal-junction adenocarcinoma receiving rilotumumab in combination with epirubicin / cisplatin / capecitabine (ECX). The structural model and parameter values were inherited from the previously developed population PK analysis of rilotumumab (Zhu et al. 2014, J Pharm Sci 103:328-336); Zhang 2016 reports the typical-value point estimates and IIV %CV from that prior model and uses it as the reference for an external visual predictive check assessing whether ECX co-administration alters rilotumumab PK."
  reference   <- paste(
    "Zhang Y, Kondragunta V, Han T-H, et al.",
    "Assessment of pharmacokinetic interaction between rilotumumab and",
    "epirubicin, cisplatin and capecitabine (ECX) in a Phase 3 study in",
    "gastric cancer. Br J Clin Pharmacol. 2017;83(5):1048-1055.",
    "doi:10.1111/bcp.13179.",
    "Structural PK model and parameter values inherited from",
    "Zhu M, Doshi S, Gisleskog PO, et al.",
    "Population pharmacokinetics of rilotumumab, a fully human monoclonal",
    "antibody against hepatocyte growth factor, in cancer patients.",
    "J Pharm Sci. 2014;103(1):328-336. doi:10.1002/jps.23763."
  )
  vignette    <- "Zhang_2016_rilotumumab"
  units       <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 279L,
    n_observations = 1460L,
    n_studies      = 1L,
    age_range      = "19-85 years",
    age_median     = "59 years",
    weight_range   = "39-120 kg",
    weight_median  = "69 kg",
    disease_state  = "Unresectable, locally advanced or metastatic MET-positive gastric or gastroesophageal-junction (GEJ) adenocarcinoma (MET-positivity defined as tumour membrane staining >= 25% by central immunohistochemistry).",
    dose_range     = "15 mg/kg IV every 3 weeks (Q3W) in combination with epirubicin 50 mg/m^2 IV bolus Q3W, cisplatin 60 mg/m^2 IV infusion Q3W, and capecitabine 625 mg/m^2 orally twice daily.",
    regions        = "Multicenter Phase 3 (ClinicalTrials.gov NCT00719550).",
    notes          = "Population summary refers to the Phase 3 cohort that informed the rilotumumab serum concentrations used for the external VPC (Zhang 2016 Results: '279 subjects had measured PK concentrations for rilotumumab'). Of these, 53 had intensive PK sampling and 226 had sparse PK sampling. A total of 1460 serum concentrations were included in the population PK evaluation (34 records excluded as outliers or below-quantitation). The structural model and parameter values themselves were estimated by Zhu 2014 using data from seven Phase 1 and Phase 2 studies; the Phase 3 weight (39-120 kg) and age (19-85 years) ranges fall within the Zhu 2014 dataset (Zhang 2016 page 1050). Body weight and age were the significant covariates retained in the Zhu 2014 final model, but Zhang 2016 does not reproduce the covariate equations or coefficient estimates; sex, cancer type, ECX co-administration, baseline HGF and MET levels, and organ functions were tested and not retained (Zhang 2016 page 1053)."
  )

  ini({
    # Structural parameters - Zhu 2014 final-model typical values as reproduced
    # by Zhang 2016 in the Results section "Population PK results of rilotumumab"
    # (page 1052-1053). Volume notation in the source is inconsistent: the
    # central volume is called both V1 and Vc in adjacent sentences, and the
    # peripheral volume is called both Vp and V2 (see vignette Errata for the
    # exact quoted text). The four typical values map unambiguously to CL, Vc,
    # Q, Vp by position in the listing.
    lcl <- log(0.184); label("Systemic clearance CL (L/day)")               # Zhang 2016 p1053: 0.184 L/day (RSE 2.5%) - Zhu 2014 final model
    lvc <- log(3.56);  label("Central volume of distribution Vc (L)")       # Zhang 2016 p1053: 3.56 L (RSE 1.5%) - Zhu 2014 final model
    lq  <- log(0.833); label("Intercompartmental clearance Q (L/day)")      # Zhang 2016 p1053: 0.833 L/day (RSE 12.3%) - Zhu 2014 final model
    lvp <- log(2.50);  label("Peripheral volume of distribution Vp (L)")    # Zhang 2016 p1053: 2.50 L (RSE 6.8%) - Zhu 2014 final model

    # Inter-individual variability. Zhang 2016 reports omega as %CV on the
    # log-normal parameters; converted to log-normal variance via
    # omega^2 = log(CV^2 + 1). No off-diagonal correlations are reported in
    # Zhang 2016, so the etas are encoded as independent diagonal terms;
    # if the upstream Zhu 2014 model reports a CL-Vc correlation block,
    # the Zhu_2014_rilotumumab extraction (task 127) is the place to encode
    # it. CV values are taken verbatim from Zhang 2016 p1053:
    # "interindividual variabilities were 29.8%, 19.8%, 71.0% and 37.3% for
    # the model parameters CL, Vc, Q and V2, respectively" (V2 here is a
    # notation slip for Vp; see vignette Errata).
    etalcl ~ 0.0851   # CL CV 29.8% -> omega^2 = log(1 + 0.298^2) = 0.0851 - Zhang 2016 p1053
    etalvc ~ 0.0385   # Vc CV 19.8% -> omega^2 = log(1 + 0.198^2) = 0.0385 - Zhang 2016 p1053
    etalq  ~ 0.4082   # Q  CV 71.0% -> omega^2 = log(1 + 0.710^2) = 0.4082 - Zhang 2016 p1053
    etalvp ~ 0.1303   # Vp CV 37.3% -> omega^2 = log(1 + 0.373^2) = 0.1303 - Zhang 2016 p1053

    # Residual error: Zhang 2016 does not report the residual error model
    # (form or magnitude) for rilotumumab. Population predictions in
    # Figure 5 of the paper were generated with the Zhu 2014 residual
    # error model but Zhang 2016 does not reproduce those values. The
    # Zhu_2014_rilotumumab extraction (task 127) is the canonical source
    # for the residual error structure; this Zhang 2016 model file
    # therefore omits a residual error declaration so the gap is visible
    # rather than papered over with a guess. See vignette Errata.
  })

  model({
    # Two-compartment IV model with first-order linear elimination from
    # the central compartment (Zhang 2016 Methods, page 1050: "The base
    # model was a two-compartment model... parameterized by systemic
    # clearance (CL), volume of distribution for the central compartment
    # (Vc), volume of distribution for the peripheral compartment (Vp),
    # and intercompartmental clearance between the central and peripheral
    # compartment (Q)"). Zhang 2016 reports body weight and age as
    # retained covariates in the Zhu 2014 final model but does not
    # publish the covariate equations or coefficients, so no covariate
    # effects are encoded here; the Zhu_2014_rilotumumab extraction
    # (task 127) is the canonical source.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg and volumes in L -> central/vc has units mg/L = ug/mL.
    Cc <- central / vc
  })
}
