Vasalou_2023_dendriticNanoparticle_mouse <- function() {
  description <- "Preclinical (mouse, SCID CB-17). PBPK (semi-mechanistic dendritic-nanoparticle biodistribution). Dual-species model carrying nanoparticle-conjugated API and nanoparticle-released (free) API through four compartments each -- blood, liver, spleen and a lumped 'other'/rest compartment (8 ODEs). Conjugated API extravasates into liver and spleen at fitted rates and releases free API by first-order, compartment-specific rate constants (t50 = 5.5 h blood, 43 h liver, 110 h spleen); free API distributes by blood flow with saturable liver and spleen partition coefficients and is cleared from the liver. No separate nanoparticle clearance term. This is the mouse parameterisation, the anchor that was fitted to data; the rat, dog and human files scale it per Eqs 22-25. Deterministic typical-value simulator: no IIV and no residual error are reported."
  reference <- "Vasalou C, Harding J, Jones RDO, Hariparsad N, McGinnity DF. Interspecies evaluation of a physiologically based pharmacokinetic model to predict the biodistribution dynamics of dendritic nanoparticles. PLoS ONE. 2023;18(5):e0285798. doi:10.1371/journal.pone.0285798. Mouse parameters from Tables 4 and 5; ODEs from Eqs 4-21; numeric values taken to full precision from the authors' MATLAB source (S2 File = ODE right-hand sides, S3 File = per-species parameter driver)."
  vignette <- "Vasalou_2023_dendriticNanoparticle"
  units <- list(
    time = "h",
    dosing = "mg/kg",
    concentration = "ng/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    blood_np  = list(analyte = "nanoparticle-conjugated API", units = NA_character_, specimen = "administration site", verified = FALSE),
    liver_np  = list(analyte = "nanoparticle-conjugated API", units = NA_character_, specimen = "tissue", verified = FALSE),
    spleen_np = list(analyte = "nanoparticle-conjugated API", units = NA_character_, specimen = "tissue", verified = FALSE),
    other_np  = list(analyte = "nanoparticle-conjugated API", units = NA_character_, specimen = "administration site", verified = FALSE),
    blood     = list(analyte = "free API", units = NA_character_, specimen = "plasma", verified = FALSE),
    liver     = list(analyte = "free API", units = NA_character_, specimen = "tissue", verified = FALSE),
    spleen    = list(analyte = "free API", units = NA_character_, specimen = "tissue", verified = FALSE),
    other     = list(analyte = "free API", units = NA_character_, specimen = "tissue", verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species        = "mouse (SCID CB-17)",
    n_subjects     = 21,
    n_studies      = 1,
    age_range      = "not reported",
    weight_range   = "0.02 kg (reference body weight, Tables 3 and 4)",
    sex_female_pct = NA_real_,
    disease_state  = "Healthy immunodeficient (SCID CB-17) mice; no tumour was implanted in the biodistribution study reported here.",
    dose_range     = "Single 10 mg/kg IV bolus of the nanoparticle (dose expressed as mg of API per kg body weight; dose volume 5 mL/kg).",
    regions        = "Single-centre preclinical (AstraZeneca R&D Boston; animals from Charles River Laboratories).",
    notes          = "Blood sampled from 3 mice at each of 20 min, 1, 6, 24, 48, 72 and 96 h post dose (21 mice total). Plasma, liver and spleen concentrations of total and released API were measured by LC-MS/MS; observed values are tabulated in Table S5 of S1 File. The mouse model was originally developed in reference [12] of the source paper and is re-used here without re-fitting, so every parameter is fixed().",
    scope_note     = "Mechanistic PBPK simulator intended for typical-value simulation: no inter-individual variability and no residual-error model are reported by the source paper. Dosing, volumes and flows are all body-weight-normalised (mg/kg, L/kg, L/h/kg), so simulate a nominal 1 kg subject and read concentrations directly."
  )

  ini({
    # ======================================================================
    # Nanoparticle-RELEASED (free) API disposition -- source paper Table 4.
    # Every value is fixed(): the mouse model was fitted in reference [12]
    # and is re-used here unchanged, so nothing is estimated in this paper.
    # Values are quoted to the precision of the authors' S3 File driver;
    # the Table 4 rounded value is given after each comment for cross-check.
    # ======================================================================
    # Vb -- "Volume of distribution in blood; apparent volume of central
    # compartment" (Table 4), hence the canonical central-volume name.
    lvc      <- fixed(log(0.33937))  ; label("Blood (central) volume of released API Vb (L/kg)")        # S3 VB_mouse; Table 4 = 0.34
    # VR -- "Rest volume; apparent volume of peripheral compartment" (Table 4).
    lvp      <- fixed(log(0.10923))  ; label("Rest (peripheral) volume of released API VR (L/kg)")      # S3 VR_mouse; Table 4 = 0.1
    # QBR -- described in the text as "intercompartmental clearance of the
    # API (QBR)", hence the canonical inter-compartmental-clearance name.
    lq       <- fixed(log(0.015265)) ; label("Blood/rest flow rate QBR (L/h/kg)")                       # S3 QPR_mouse; Table 4 = 0.015
    # CL -- released API clearance from the liver (Eq 5). Nominal mouse value
    # is the in-vivo blood CL from NCA of the unconjugated-API IV profile
    # (Table 2); the paper also simulates CL = 1 (ref [12]) and 2 (in-vitro
    # scaled) in Fig 4, and S3 records all three in a comment.
    lcl      <- fixed(log(2.3))      ; label("Released API clearance CL (L/h/kg)")                      # S3 CL_mouse = 2.3 (NCA); Table 2 in-vivo blood CL = 2.3
    v_liver  <- fixed(0.065)         ; label("Liver volume VL (L/kg)")                                  # S3 VL_mouse; Table 4 = 0.065 (physiological, Davies et al [5])
    v_spleen <- fixed(0.005)         ; label("Spleen volume VS (L/kg)")                                 # S3 VS_mouse = 0.1/1000/0.02; Table 4 = 0.005
    qbl_o    <- fixed(9.12)          ; label("Liver to blood flow rate QBL,o (L/h/kg)")                 # S3 QPLO_mouse; Table 4 = 9.1 (physiological)
    qbs      <- fixed(0.744)         ; label("Blood/spleen flow rate QBS (L/h/kg)")                     # S3 QPS_mouse = 14.88/1000/0.02; Table 4 = 0.74
    bmaxl    <- fixed(49.649)        ; label("Max binding capacity in liver BmaxL (ng/mL)")             # S3 BMAXL_mouse; Table 4 = 49.6
    kdl      <- fixed(0.59915)       ; label("Dissociation constant in liver KDL (ng/mL)")              # S3 KDL; Table 4 = 0.6
    pl       <- fixed(0.75973)       ; label("Non-specific binding in liver PL (unitless)")             # S3 PTL_mouse; Table 4 = 0.76
    bmaxs    <- fixed(125.91)        ; label("Max binding capacity in spleen BmaxS (ng/mL)")            # S3 BMAXS_mouse; Table 4 = 125.9
    kds      <- fixed(0.63238)       ; label("Dissociation constant in spleen KDS (ng/mL)")             # S3 KDS; Table 4 = 0.63
    ps       <- fixed(0.16919)       ; label("Non-specific binding in spleen PS (unitless)")            # S3 PTS_mouse; Table 4 = 0.17
    bpr      <- fixed(0.84)          ; label("Blood-to-plasma ratio of API BPR (unitless)")             # S3 BPR_mouse; Table 4 = 0.84
    # Hematocrit appears ONLY in the authors' S3 File, never in the paper
    # text or tables. It converts conjugated-API blood to plasma (Eq 19).
    hct      <- fixed(0.45)          ; label("Hematocrit H (unitless)")                                 # S3 Hem = 0.45 (not reported in the paper text/tables)

    # ======================================================================
    # Nanoparticle-CONJUGATED API disposition -- source paper Table 5.
    # ======================================================================
    vnb          <- fixed(0.085)      ; label("Blood volume for conjugated API VNb (L/kg)")             # S3 V1_mouse; Table 5 = 0.085 (physiological blood volume)
    nbl          <- fixed(0.00029801) ; label("Blood/liver extravasation rate NBL (L/h/kg)")            # S3 NPL_mouse; Table 5 = 0.0003
    nbs          <- fixed(1.9732e-5)  ; label("Blood/spleen extravasation rate NBS (L/h/kg)")           # S3 NPS_mouse; Table 5 = 2e-5
    # NBR fixed at exactly zero: Table 5 footnote "+ NBR was set to zero
    # based on previous assessments [12], since if it was allowed to float,
    # obtained value approximately equal to zero."
    nbr          <- fixed(0)          ; label("Blood/rest distribution rate NBR (L/h/kg)")              # S3 NPR_mouse = 0.0; Table 5 = 0
    knbl         <- fixed(0.89999)    ; label("Liver/blood partition coefficient KNBL (unitless)")      # S3 KNPL_mouse; Table 5 = 0.9
    # KNBS fixed at 1000: Table 5 footnote "++ KNBS was set to 1000 based on
    # previous model fits [12] ... Any distribution rate from spleen back to
    # blood is effectively removed".
    knbs         <- fixed(1000)       ; label("Spleen/blood partition coefficient KNBS (unitless)")     # S3 parameters(32); Table 5 = 1000
    krel_b       <- fixed(0.1246)     ; label("API release rate in blood krelb (1/h)")                  # S3 KREL; Table 5 = 0.125
    krel_l       <- fixed(0.016374)   ; label("API release rate in liver krelL (1/h)")                  # S3 KRLV; Table 5 = 0.016
    krel_s       <- fixed(0.0063168)  ; label("API release rate in spleen krelS (1/h)")                 # S3 KRLS; Table 5 = 0.0063
    # Vascular volume fractions used by Eqs 20-21 to add residual blood to
    # TOTAL tissue concentrations. Table 5 prints 0.125 / 0.016, but those
    # two cells are byte-identical to the krelb / krelL cells one row-group
    # above and are a transcription error; the authors' S3 File uses
    # 0.10739 / 0.055634, which are the values that generated Figs 4-9.
    # Operator-ratified 2026-07-29 (sidecar request-001 q2 = A). See vignette
    # "Assumptions and deviations" for the full comparison.
    fvasc_liver  <- fixed(0.10739)    ; label("Liver vascular volume fraction vliver (unitless)")       # S3 parameters(39); Table 5 prints 0.125 (transcription error)
    fvasc_spleen <- fixed(0.055634)   ; label("Spleen vascular volume fraction vspleen (unitless)")     # S3 parameters(41); Table 5 prints 0.016 (transcription error)
  })

  model({
    # ==================================================================
    # Canonical log-scale parameters back to the linear scale
    # ==================================================================
    vc <- exp(lvc)   # Vb, L/kg
    vp <- exp(lvp)   # VR, L/kg
    q  <- exp(lq)    # QBR, L/h/kg
    cl <- exp(lcl)   # CL, L/h/kg

    # ==================================================================
    # Derived structural quantities
    # ==================================================================
    # Blood-to-liver (arterial) flow. Table 4 defines QBL,i as QBL,o - QBS
    # because the spleen drains into the liver via the portal route.
    qbl_i <- qbl_o - qbs

    # Rest volume for the conjugated API. Table 5: VNR = 1 - sum(Vtissues).
    # S3 File: parameters(25) = 1 - V1 - V3 - V5 = 1 - VNb - VL - VS.
    vnr <- 1 - vnb - v_liver - v_spleen

    # Released API blood concentration in ng/mL (paper Eq 8, scaled by 1000
    # because amounts are mg/kg and volumes L/kg so amount/volume is ug/mL).
    # This drives the saturable partition coefficients, whose BmaxL / KDL
    # are reported in ng/mL, so the unit conversion is load-bearing.
    cb_rel <- blood / vc * 1000

    # Saturable tissue/blood partition coefficients of the released API
    # (paper Eqs 12 and 13; S3 File KPL / KPS).
    kbl <- bmaxl / (cb_rel + kdl) + pl
    kbs <- bmaxs / (cb_rel + kds) + ps

    # ==================================================================
    # Nanoparticle-CONJUGATED API amounts (paper Eqs 14-17)
    # The dose is administered into blood_np: only the conjugate is dosed.
    # krelR = krelb per Table 5 ("krelR = krel b").
    # ==================================================================
    d/dt(blood_np) <- -nbl * (blood_np / vnb - liver_np / (v_liver * knbl)) -
                       nbs * (blood_np / vnb - spleen_np / (v_spleen * knbs)) -
                       nbr * (blood_np / vnb - other_np / vnr) -
                       krel_b * blood_np
    d/dt(liver_np)  <-  nbl * (blood_np / vnb - liver_np / (v_liver * knbl)) -
                        krel_l * liver_np
    d/dt(spleen_np) <-  nbs * (blood_np / vnb - spleen_np / (v_spleen * knbs)) -
                        krel_s * spleen_np
    d/dt(other_np)  <-  nbr * (blood_np / vnb - other_np / vnr) -
                        krel_b * other_np

    # ==================================================================
    # Nanoparticle-RELEASED (free) API amounts (paper Eqs 4-7)
    # Each compartment gains free API from release of the co-located
    # conjugate. The liver receives arterial inflow QBL,i plus the spleen's
    # portal outflow QBS, and returns the whole QBL,o to blood.
    # ==================================================================
    d/dt(blood)  <- -qbl_i * (blood / vc) +
                     qbl_o * (liver / (v_liver * kbl)) -
                     qbs * (blood / vc) -
                     q * (blood / vc - other / vp) +
                     krel_b * blood_np
    d/dt(liver)  <-  qbl_i * (blood / vc) -
                     qbl_o * (liver / (v_liver * kbl)) -
                     cl * (liver / (v_liver * kbl)) +
                     qbs * (spleen / (v_spleen * kbs)) +
                     krel_l * liver_np
    d/dt(spleen) <-  qbs * (blood / vc - spleen / (v_spleen * kbs)) +
                     krel_s * spleen_np
    d/dt(other)  <-  q * (blood / vc - other / vp) +
                     krel_b * other_np

    # ==================================================================
    # Observation variables, all in ng/mL
    # ==================================================================
    # Released (free) API. Cc is the canonical central observation and holds
    # the released API in PLASMA, which is what the assay reports (Eq 9).
    Cc               <- cb_rel / bpr                        # Eq 9
    CbloodReleased   <- cb_rel                              # Eq 8
    CliverReleased   <- liver / v_liver * 1000              # Eq 10
    CspleenReleased  <- spleen / v_spleen * 1000            # Eq 11

    # Nanoparticle-conjugated API. The conjugate is confined to plasma, so
    # the blood-to-plasma conversion is 1/(1 - H) rather than 1/BPR.
    Cc_np            <- blood_np / (vnb * (1 - hct)) * 1000  # Eq 19, conjugated term
    Cblood_np        <- blood_np / vnb * 1000                # Eq 18, conjugated term
    Cliver_np        <- liver_np / v_liver * 1000
    Cspleen_np       <- spleen_np / v_spleen * 1000

    # TOTAL API = released + conjugated (paper Eqs 18-21). Tissue totals add
    # residual blood via the vascular volume fractions.
    CbloodTotal      <- CbloodReleased + Cblood_np                                  # Eq 18
    CplasmaTotal     <- Cc + Cc_np                                                  # Eq 19
    CliverTotal      <- CliverReleased + Cliver_np + fvasc_liver * CbloodTotal      # Eq 20
    CspleenTotal     <- CspleenReleased + Cspleen_np + fvasc_spleen * CbloodTotal   # Eq 21
  })
}
