Cook_2016_paracetamol <- function() {
  description <- "Population PK model for paracetamol (APAP) and its glucuronide and sulphate conjugates with cumulative urinary excretion in term and preterm newborns (Cook 2016), as packaged in DDMORE Foundation Model Repository entry DDMODEL00000271."
  reference <- paste(
    "Cook SF, Stockmann C, Samiee-Zafarghandy S, King AD, Deutsch N,",
    "Williams EF, Wilkins DG, Sherwin CMT, van den Anker JN (2016).",
    "Neonatal maturation of paracetamol (acetaminophen) glucuronidation,",
    "sulfation, and oxidation based on a parent-metabolite population",
    "pharmacokinetic model.",
    "Clin Pharmacokinet 55(11):1395-1411.",
    "doi:10.1007/s40262-016-0408-1.",
    "DDMORE Foundation Model Repository: DDMODEL00000271.",
    sep = " "
  )
  vignette <- "Cook_2016_paracetamol"
  units <- list(time = "min", dosing = "mg", concentration = "mg/L")

  ddmore_id    <- "DDMODEL00000271"
  replicate_of <- NULL

  covariateData <- list(
    WT = list(
      description        = "Body weight at time of study",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying weight at the time of each PK sample.",
        "Source column 'BWS' (body weight at study) is mapped to the",
        "canonical 'WT' on input. Used as a multiplicative scaler with no",
        "reference weight: vc, cl_apapg, cl, and cl_apaps are all per-kg",
        "parameters multiplied by BWS (with a fitted exponent of 1.40 on",
        "BWS for cl_apaps only)."
      ),
      source_name        = "BWS"
    )
  )

  population <- list(
    n_subjects     = 54,
    n_studies      = "Not extractable from DDMORE bundle (Cook 2016 PDF not on disk).",
    age_range      = "Term and preterm newborns. Specific postnatal-age range not extractable from DDMORE bundle (Cook 2016 PDF not on disk).",
    weight_range   = "Newborn body weights. The bundle's simulated dataset (Simulated_ParacetamolPKnewborns.csv) contains BWS values 0.5-4 kg for 9 of 10 subjects plus one outlier at 6.5 kg; the simulated dataset is a smoke-test cohort and does not represent the publication's reported demographics.",
    sex_female_pct = "Not extractable from DDMORE bundle (Cook 2016 PDF not on disk).",
    race_ethnicity = "Not extractable from DDMORE bundle (Cook 2016 PDF not on disk).",
    disease_state  = "Term and preterm newborns receiving IV paracetamol. Specific clinical setting not extractable from DDMORE bundle alone; Cook 2016 reports a parent-metabolite population PK analysis describing the maturation of paracetamol glucuronidation, sulfation, and oxidation in newborns.",
    dose_range     = "IV paracetamol given as a short infusion. Bundle's simulated dataset uses ~10 mg/kg single doses (5, 10, 20, 35 mg paired with BWS 0.5, 1, 2, 3.5 kg respectively) infused over approximately 15 minutes (RATE=AMT/15 mg/min).",
    regions        = "Not extractable from DDMORE bundle (Cook 2016 PDF not on disk).",
    notes          = "N=54 subjects taken from the .lst FINAL ETABAR / shrinkage block ('N: 54 54 54 54'). Full demographics, study design, and inclusion criteria could not be cross-checked because the Cook 2016 publication PDF is not on disk. The DDMORE Model_Accommodations.txt note states only that the publication reported additive errors on urine recoveries while the deposited code uses (correct) proportional errors."
  )

  ini({
    # Structural typical values come from the Output_real_ParacetamolInNewborns.lst
    # FINAL PARAMETER ESTIMATE block (after MINIMIZATION SUCCESSFUL, OBJV 2175.556).
    # The .mod $PK applies the scale factor /1000 inline (TH2/TH3/TH4); the
    # log() arguments below carry that /1000 so the bare values reappear in
    # mechanistic units (L/min/kg). The per-kg parameterisation matches the
    # source: vc, cl_apapg, cl all scale linearly with BWS, while cl_apaps
    # scales as BWS^e_wt_cl_apaps with the exponent estimated.
    lvc       <- log(1.06)         ; label("V1 plasma APAP per kg (L/kg)")  # .lst FINAL TH1
    lcl_apapg <- log(0.266 / 1000) ; label("Formation clearance to APAP-G per kg (L/min/kg)")  # .lst FINAL TH2 with /1000 scale per .mod
    lcl_apaps <- log(1.46  / 1000) ; label("Formation clearance to APAP-S coefficient (L/min at BWS=1 kg)")  # .lst FINAL TH3 with /1000 scale per .mod
    lcl       <- log(0.285 / 1000) ; label("Renal clearance of unchanged APAP per kg (L/min/kg)")  # .lst FINAL TH4 with /1000 scale per .mod

    # Allometric exponent estimated only on cl_apaps (the only fitted exponent
    # on BWS in the .mod $PK; the other clearances and V1 use exponent fixed
    # at 1.0 implicitly via direct multiplication by BWS).
    e_wt_cl_apaps <- 1.40 ; label("Power exponent on BWS for cl_apaps (unitless)")  # .lst FINAL TH6

    # Ratio relating the metabolite-to-urine elimination rates K24 = K36 to
    # the parent unchanged-renal rate K15 = cl/vc. Paper-named, dimensionless;
    # documents the .mod construct K24 = TH5*K15 / K36 = K24 with TH5 fitted.
    kratio_urine <- 11.3 ; label("Metabolite/parent renal-excretion rate ratio K24 = K36 = kratio_urine * K15 (unitless)")  # .lst FINAL TH5

    # Inter-individual variability — diagonal $OMEGA. NONMEM ETA is on
    # the log scale (V1 = TH1*BWS*EXP(ETA(1)) etc.), so OMEGA values are
    # variances on log-clearance / log-volume.
    etalvc        ~ 0.0925  # .lst FINAL OMEGA(1,1) — IIV on log V1
    etalcl_apapg  ~ 0.599   # .lst FINAL OMEGA(2,2) — IIV on log CLG
    etalcl_apaps  ~ 0.312   # .lst FINAL OMEGA(3,3) — IIV on log CLS
    etalcl        ~ 0.0879  # .lst FINAL OMEGA(4,4) — IIV on log CLA

    # Residual error — NONMEM SIGMAs are variances; nlmixr2 propSd / addSd
    # are SDs, hence sqrt(.). The .mod $ERROR uses combined prop+add on
    # plasma APAP (Cc) and proportional only on each urine output. This
    # matches the deposited code (per Model_Accommodations.txt the
    # publication described urinary errors as additive but the deposited
    # code uses the correct proportional structure).
    propSd               <- sqrt(0.0198) ; label("Cc plasma APAP proportional residual (fraction)")     # .lst FINAL SIGMA(1,1)
    addSd                <- sqrt(0.354)  ; label("Cc plasma APAP additive residual (mg/L)")              # .lst FINAL SIGMA(5,5)
    propSd_Aurine_apapg  <- sqrt(0.223)  ; label("Cumulative urine APAP-G amount proportional residual (fraction)")  # .lst FINAL SIGMA(2,2)
    propSd_Aurine_apap   <- sqrt(0.188)  ; label("Cumulative urine APAP amount proportional residual (fraction)")    # .lst FINAL SIGMA(3,3)
    propSd_Aurine_apaps  <- sqrt(0.332)  ; label("Cumulative urine APAP-S amount proportional residual (fraction)")  # .lst FINAL SIGMA(4,4)
  })
  model({
    # Individual PK parameters (per .mod $PK block).
    # vc, cl_apapg, cl scale linearly with BWS; cl_apaps scales with BWS^1.40.
    vc       <- exp(lvc       + etalvc)       * WT
    cl_apapg <- exp(lcl_apapg + etalcl_apapg) * WT
    cl_apaps <- exp(lcl_apaps + etalcl_apaps) * WT^e_wt_cl_apaps
    cl       <- exp(lcl       + etalcl)       * WT

    # Volume of metabolite distribution: V2 = V3 = 0.18 * V1 per .mod $PK.
    v_metab <- 0.18 * vc

    # Micro-constants (per .mod $PK):
    #   K12 = CLG/V1                (formation rate APAP -> APAP-G)
    #   K13 = CLS/V1                (formation rate APAP -> APAP-S)
    #   K15 = CLA/V1                (urinary excretion rate of unchanged APAP)
    #   K24 = K36 = TH5 * K15       (urinary excretion rate of metabolites)
    k_apap_apapg  <- cl_apapg / vc
    k_apap_apaps  <- cl_apaps / vc
    k_apap_urine  <- cl       / vc
    k_apapg_urine <- kratio_urine * k_apap_urine
    k_apaps_urine <- kratio_urine * k_apap_urine

    # ODE system (per .mod $DES block, with named compartments):
    #   A(1) -> central        (plasma APAP)
    #   A(2) -> central_apapg  (plasma APAP-G)
    #   A(3) -> central_apaps  (plasma APAP-S)
    #   A(4) -> urine_apapg    (cumulative urine APAP-G amount, mg)
    #   A(5) -> urine_apap     (cumulative urine APAP amount, mg)
    #   A(6) -> urine_apaps    (cumulative urine APAP-S amount, mg)
    d/dt(central)        <- -k_apap_apapg * central - k_apap_apaps * central - k_apap_urine * central
    d/dt(central_apapg)  <-  k_apap_apapg * central - k_apapg_urine * central_apapg
    d/dt(central_apaps)  <-  k_apap_apaps * central - k_apaps_urine * central_apaps
    d/dt(urine_apapg)    <-  k_apapg_urine * central_apapg
    d/dt(urine_apap)     <-  k_apap_urine  * central
    d/dt(urine_apaps)    <-  k_apaps_urine * central_apaps

    # Observations:
    #   Cc = A1 / V1 — plasma APAP concentration in mg/L (.mod S1=V1)
    #   Aurine_* = A4/A5/A6 — cumulative urinary amounts in mg (.mod S4=S5=S6=1)
    Cc           <- central / vc
    Aurine_apap  <- urine_apap
    Aurine_apapg <- urine_apapg
    Aurine_apaps <- urine_apaps

    Cc           ~ add(addSd) + prop(propSd)
    Aurine_apap  ~ prop(propSd_Aurine_apap)
    Aurine_apapg ~ prop(propSd_Aurine_apapg)
    Aurine_apaps ~ prop(propSd_Aurine_apaps)
  })
}
