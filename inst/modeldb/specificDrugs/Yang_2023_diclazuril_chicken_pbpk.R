Yang_2023_diclazuril_chicken_pbpk <- function() {
  description <- "PBPK (whole-body, flow-limited; broiler chicken). Nine-compartment physiologically based pharmacokinetic model for the anticoccidial diclazuril in broiler chickens after continuous oral exposure via medicated feed or drinking water, comprising intestinal contents (absorption site), liver, kidney, lumped skin + fat, muscle, a lumped rest-of-body compartment, lung, arterial plasma and venous plasma; all tissues are perfusion (flow) limited with tissue:plasma partition coefficients, absorption from the gut lumen is first order (Ka) in competition with first-order fecal loss of unabsorbed drug (Kgut), and elimination is hepatic (Clhe) plus fecal excretion. Built to predict edible-tissue residues and withdrawal periods against Chinese and European maximum residue limits (Yang 2023)."
  reference   <- "Yang F, Zhang M, Jin Y-G, Chen J-C, Duan M-H, Liu Y, Li Z-E, Li X-P, Yang F. Development and Application of a Physiologically Based Pharmacokinetic Model for Diclazuril in Broiler Chickens. Animals (Basel). 2023;13(9):1512. doi:10.3390/ani13091512"
  vignette    <- "Yang_2023_diclazuril"
  units       <- list(time = "h", dosing = "ug", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Scales every compartment volume (Vxx = Vcxx * WT) and cardiac output (Qtot = CO * WT), and converts the hepatic clearance from L/h/kg to L/h (CLhe = Clhe * WT). Yang 2023 held body weight constant over the whole treatment and simulation period because breed-specific growth data were unavailable (Yang 2023 Section 2.3 and Table 3 footnote 1); the supplementary acslX code sets `constant bw = 1.5`, with 1.34 kg and 1.52 kg used for the two Wen/He datasets. Per-study values are 1.52 kg (refs [1,19]), 1.34 kg (ref [20]), 1.5 kg (ref [22]) and 1.23 kg (ref [21]) (Yang 2023 Table 1).",
      source_name        = "BW"
    )
  )

  population <- list(
    species        = "broiler chicken (Gallus gallus domesticus; Lingnan Yellow Chicken and Ross 308 broilers)",
    n_subjects     = NA_integer_,
    n_studies      = 4L,
    age_range      = "15-50 days (15, 21, 30 and 50 days across the four contributing studies; Yang 2023 Table 1)",
    weight_range   = "1.23-1.52 kg (study mean body weights; Yang 2023 Table 1)",
    sex_female_pct = NA_real_,
    disease_state  = "Healthy (coccidiosis-free) broilers; diclazuril given prophylactically",
    dose_range     = "Single oral gavage 80 ug/kg BW or 1 mg/kg BW; medicated feed 730 ug/kg or 1 mg/kg; medicated water 0.5-1 mg/L (recommended) or 3 mg/L (validation study), for 7-10 consecutive days (Yang 2023 Table 1 and Section 2.6)",
    regions        = "China (Lingnan Yellow Chicken studies) and Europe (Ross 308 study, Mortier et al.)",
    n_observations = NA_integer_,
    notes          = "Model built by digitising / tabulating published plasma and tissue concentration-time data from five datasets in four source studies (Yang 2023 Table 1): three sets for parameter optimisation (refs [1,19] single 1 mg/kg BW gavage plus 7-day medicated feed at 1 mg/kg; ref [21] single 80 ug/kg BW gavage) and two for external validation (ref [20] medicated water 3 mg/L for 9 days; ref [22] medicated feed 730 ug/kg for 10 days). Physiological parameters are population means for broilers from the Wang et al. 2020 review (Yang 2023 ref [24]); no individual-level data were fitted, so the model carries no inter-individual random effects. Yang 2023 propagated parameter uncertainty instead through a 1000-iteration Monte Carlo analysis over the normal distributions in Yang 2023 Table 4. Light regime governs daily exposure duration because chickens stop eating and drinking in the dark: 12/12 h (refs [1,19,20]), 21/3 h (ref [22]) and 18/6 h (ref [21]). Daily intake was assumed to be 1.1 kg feed and 0.5 L water per bird, evenly spread over the light hours (Yang 2023 Section 2.2)."
  )

  ini({
    # ---------------------------------------------------------------------
    # Tissue:plasma partition coefficients (Pxx, unitless), log-transformed
    # following the Zhang_2011_nutlin3a.R `lk_<organ>` convention. Muscle,
    # kidney and skin + fat were calculated by the area method (Gallo 1987;
    # Yang 2023 ref [25]); liver, lung and the rest compartment were
    # optimised against the published concentration-time data. All are
    # FIXED at the published point estimates.
    # ---------------------------------------------------------------------
    lk_mus <- fixed(log(0.1299)); label("Muscle:plasma partition coefficient (Pmu, unitless)")       # Yang 2023 Table 3 / Section 3.1 (area method)
    lk_skf <- fixed(log(0.0955)); label("Skin + fat:plasma partition coefficient (Psf, unitless)")   # Yang 2023 Table 3 / Section 3.1 (area method)
    lk_kid <- fixed(log(0.6813)); label("Kidney:plasma partition coefficient (Pki, unitless)")       # Yang 2023 Table 3 / Section 3.1 (area method)
    lk_liv <- fixed(log(0.9613)); label("Liver:plasma partition coefficient (Pli, unitless)")        # Yang 2023 Table 3 + footnote 5 / Section 3.1 (optimised)
    lk_lun <- fixed(log(0.5603)); label("Lung:plasma partition coefficient (Plu, unitless)")         # Yang 2023 Table 3 + footnote 5 / Section 3.1 (optimised)
    lk_res <- fixed(log(1.2965)); label("Rest-of-body:plasma partition coefficient (Pre, unitless)") # Yang 2023 Table 3 + footnote 5 / Section 3.1 (optimised)

    # ---------------------------------------------------------------------
    # Absorption and elimination, all optimised in acslXtreme with the
    # Nelder-Mead maximum-likelihood algorithm (Yang 2023 Section 2.3) and
    # reported as point estimates in Section 3.1. FIXED: Yang 2023 reports
    # no standard errors on these, and the Monte Carlo SDs in Table 4 are
    # uncertainty for sensitivity analysis, not estimated random effects.
    # ---------------------------------------------------------------------
    lka       <- fixed(log(0.1234));   label("First-order absorption rate constant from intestinal contents (Ka, 1/h)")        # Yang 2023 Section 3.1 / Table 4 / supplement `constant ka=0.1234`
    lkgut     <- fixed(log(0.3838));   label("First-order fecal excretion rate constant for unabsorbed drug (Kgut, 1/h)")      # Yang 2023 Section 3.1 / Table 4 / supplement `CONSTANT kgut=0.3838`
    lclhe     <- fixed(log(0.00344));  label("Hepatic clearance per unit body weight (Clhe, L/h/kg)")                          # Yang 2023 Section 3.1 / Table 4 / supplement `constant clhe=0.00344`

    # ---------------------------------------------------------------------
    # Residual error. Yang 2023 is a deterministic PBPK model fitted to
    # pooled/digitised literature means; it reports no residual-error model
    # (validation used MAPE and linear regression instead, Yang 2023
    # Table 5). The placeholders below exist only so the model is a
    # syntactically complete nlmixr2 object for forward simulation; they
    # are NOT paper-derived. Same convention as
    # Gaohua_2012_pregnancy_pbpk_midazolam.R and
    # An_2012_mitoxantrone_mouse_pbpk.R.
    # ---------------------------------------------------------------------
    propSd            <- fixed(0.10); label("Proportional residual error placeholder, venous plasma (fraction)") # not reported in Yang 2023; placeholder for syntactic completeness only
    propSd_Cmuscle    <- fixed(0.10); label("Proportional residual error placeholder, muscle (fraction)")        # not reported in Yang 2023; placeholder for syntactic completeness only
    propSd_Cliver     <- fixed(0.10); label("Proportional residual error placeholder, liver (fraction)")         # not reported in Yang 2023; placeholder for syntactic completeness only
    propSd_Ckidney    <- fixed(0.10); label("Proportional residual error placeholder, kidney (fraction)")        # not reported in Yang 2023; placeholder for syntactic completeness only
    propSd_Cskin_fat  <- fixed(0.10); label("Proportional residual error placeholder, skin + fat (fraction)")    # not reported in Yang 2023; placeholder for syntactic completeness only
    propSd_Clung      <- fixed(0.10); label("Proportional residual error placeholder, lung (fraction)")          # not reported in Yang 2023; placeholder for syntactic completeness only
  })

  model({
    # ==================================================================
    # Broiler physiology. Population means from the Wang et al. 2020
    # food-animal physiology review (Yang 2023 ref [24]), reproduced in
    # Yang 2023 Table 3 / Table 4 and in the supplementary acslX INITIAL
    # block. These are literature physiology rather than fitted
    # quantities, so they are carried as traceable literals here rather
    # than as ini() parameters - the same convention as
    # Gaohua_2012_pregnancy_pbpk_midazolam.R and
    # An_2012_mitoxantrone_mouse_pbpk.R.
    #
    # Skin and fat are carried SEPARATELY and summed, and blood is
    # carried WHOLE, because that is what the executable supplement
    # does. Yang 2023 Table 3 prints only the lumped forms, which agree
    # once decomposed: skin + fat 0.1338 + 0.134 = 0.2678 (weight) and
    # 0.1505 + 0.1 = 0.2505 (flow); arterial 0.0322 + venous 0.0161 =
    # 0.0483 = the whole-blood fraction Vcbl of Table 4.
    # ==================================================================
    co        <- 9.88                     # cardiac output, L/h/kg BW; Yang 2023 Table 3 footnote 2 / supplement `constant CO=9.88`
    pcv       <- 0.32                     # hematocrit, converts blood volume to plasma volume; Yang 2023 Section 2.3 (32 +/- 2.76%) / supplement `constant pcv=0.32`

    vc_muscle <- 0.5712                   # muscle weight / BW;  Yang 2023 Table 3 / supplement `constant vcmu=0.5712`
    vc_liver  <- 0.0214                   # liver weight / BW;   Yang 2023 Table 3 / supplement `constant vcli=0.0214`
    vc_kidney <- 0.0064                   # kidney weight / BW;  Yang 2023 Table 3 / supplement `constant vcki=0.0064`
    vc_lung   <- 0.0071                   # lung weight / BW;    Yang 2023 Table 3 / supplement `constant vclu=0.0071`
    vc_blood  <- 0.0483                   # whole-blood vol / BW; Yang 2023 Table 4 (Vcbl 4.83%) / supplement `constant vcbl=0.0483`
    vc_skin   <- 0.1338                   # skin weight / BW;    Yang 2023 Table 4 (Vcsk 13.38%) / supplement `constant vcsk=0.1338`
    vc_fat    <- 0.134                    # fat weight / BW;     Yang 2023 Table 4 (Vcfa 13.4%) / supplement `constant vcfa=0.134`

    qc_muscle <- 0.0764                   # muscle flow / Qtot;  Yang 2023 Table 3 / supplement `constant qcmu=0.0764`
    qc_liver  <- 0.2526                   # liver flow / Qtot (hepatic artery + portal vein); Yang 2023 Table 3 + footnote 4 / supplement `constant qcli=0.2526`
    qc_kidney <- 0.2012                   # kidney flow / Qtot;  Yang 2023 Table 3 / supplement `constant qcki=0.2012`
    qc_lung   <- 1.0                      # lung flow / Qtot = 1; Yang 2023 Table 3 + footnote 6 / supplement `constant qclu=1`
    qc_skin   <- 0.1505                   # skin flow / Qtot;    Yang 2023 supplement `constant qcsk=0.1505`
    qc_fat    <- 0.1                      # fat flow / Qtot;     Yang 2023 supplement `constant qcfa=0.1`

    # ------------------------------------------------------------------
    # Drug-specific parameters, back-transformed from the log scale.
    # ------------------------------------------------------------------
    ka         <- exp(lka)
    kgut       <- exp(lkgut)
    clhe       <- exp(lclhe)
    k_mus      <- exp(lk_mus)
    k_skf      <- exp(lk_skf)
    k_kid      <- exp(lk_kid)
    k_liv      <- exp(lk_liv)
    k_lun      <- exp(lk_lun)
    k_res      <- exp(lk_res)

    # ------------------------------------------------------------------
    # Compartment volumes (L, or kg taken as L at unit tissue density).
    # Yang 2023 supplement acslX INITIAL block.
    # ------------------------------------------------------------------
    v_muscle   <- vc_muscle * WT
    v_liver    <- vc_liver  * WT
    v_kidney   <- vc_kidney * WT
    v_lung     <- vc_lung   * WT
    v_blood    <- vc_blood  * WT
    v_skin_fat <- (vc_skin + vc_fat) * WT                                  # supplement: vsf = vfa + vsk
    # Rest of body is the balance of body weight after every named tissue
    # (blood is counted whole here, matching the supplement).
    v_other    <- (1 - (vc_muscle + vc_liver + vc_kidney + vc_blood +
                        vc_fat + vc_lung + vc_skin)) * WT                  # supplement: vcre = 1 - (...); = 0.0778 (Yang 2023 Table 3 footnote 7)

    # Arterial plasma is one third and venous plasma two thirds of blood
    # volume, each converted from blood to plasma with the hematocrit.
    # This is the supplement's assignment; Yang 2023 Table 3 prints the
    # arterial and venous rows the other way round AND without the
    # (1 - pcv) factor - see the vignette Errata.
    v_arterial <- v_blood / 3     * (1 - pcv)                              # supplement: vab = vbl/3;   vap = vab*(1-pcv)
    v_venous   <- v_blood * 2 / 3 * (1 - pcv)                              # supplement: vvb = vbl*2/3; vvp = vvb*(1-pcv)

    # ------------------------------------------------------------------
    # Plasma flows (L/h) and hepatic clearance (L/h).
    # ------------------------------------------------------------------
    q_total    <- co * WT                                                  # supplement: QTOT = CO*bw
    q_muscle   <- q_total * qc_muscle
    q_liver    <- q_total * qc_liver
    q_kidney   <- q_total * qc_kidney
    q_lung     <- q_total * qc_lung                                        # equals q_total (Qclu = 1)
    q_skin_fat <- q_total * (qc_skin + qc_fat)                             # supplement: qsf = qfa + qsk
    q_other    <- q_total * (1 - (qc_muscle + qc_liver + qc_kidney +
                                  qc_fat + qc_skin))                       # supplement: qcre = 1 - (...) excluding lung; = 0.2193 (Yang 2023 Table 3 footnote 8)
    cl_hepatic <- clhe * WT                                                # supplement: CCLhe = clhe*bw

    # ------------------------------------------------------------------
    # Compartment concentrations (ug/L, or ug/kg at unit tissue density).
    # States are drug AMOUNTS in ug, exactly as the supplement integrates.
    # ------------------------------------------------------------------
    Cc         <- venous   / v_venous                                      # venous plasma: the sampled matrix
    c_arterial <- arterial / v_arterial
    Clung      <- lung     / v_lung
    Cliver     <- liver    / v_liver
    Ckidney    <- kidney   / v_kidney
    Cskin_fat  <- skin_fat / v_skin_fat
    Cmuscle    <- muscle   / v_muscle
    c_other    <- other    / v_other

    # ------------------------------------------------------------------
    # Mass-balance ODEs (Yang 2023 Table 2; supplement DERIVATIVE block).
    #
    # Dosing: the daily dose (ug) is the diclazuril concentration in feed
    # or water times the daily intake, delivered to `a_gut` as a
    # zero-order input over the daily light period (tlen h) and repeated
    # every 24 h. In the acslX code this is `dose/tlen * PULSE(...)`;
    # here it is carried by the event table (rate / duration on the dose
    # record), so tlen is a property of the data, not of the model.
    # ------------------------------------------------------------------
    d/dt(a_gut)    <- -ka * a_gut - kgut * a_gut                           # supplement: ragicon = Rdose - rabsp - routgut
    d/dt(liver)    <- ka * a_gut +
                      q_liver * (c_arterial - Cliver / k_liv) -
                      cl_hepatic * Cliver / k_liv                       # absorption is into the liver (first pass)
    d/dt(kidney)   <- q_kidney   * (c_arterial - Ckidney   / k_kid)
    d/dt(skin_fat) <- q_skin_fat * (c_arterial - Cskin_fat / k_skf)
    d/dt(muscle)   <- q_muscle   * (c_arterial - Cmuscle   / k_mus)
    d/dt(other)    <- q_other    * (c_arterial - c_other   / k_res)
    d/dt(lung)     <- q_lung     * (Cc         - Clung     / k_lun)      # venous plasma perfuses the lung
    d/dt(arterial) <- q_lung     * (Clung / k_lun - c_arterial)
    d/dt(venous)   <- q_muscle   * Cmuscle   / k_mus +
                      q_skin_fat * Cskin_fat / k_skf +
                      q_kidney   * Ckidney   / k_kid +
                      q_liver    * Cliver    / k_liv +
                      q_other    * c_other   / k_res -
                      q_lung     * Cc

    # ------------------------------------------------------------------
    # Observations. Venous plasma plus the four regulated edible tissues
    # (muscle, skin + fat, kidney, liver) and lung.
    # ------------------------------------------------------------------
    Cc        ~ prop(propSd)
    Cmuscle   ~ prop(propSd_Cmuscle)
    Cliver    ~ prop(propSd_Cliver)
    Ckidney   ~ prop(propSd_Ckidney)
    Cskin_fat ~ prop(propSd_Cskin_fat)
    Clung     ~ prop(propSd_Clung)
  })
}
