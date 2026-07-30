Vasalou_2023_dendriticNanoparticle_rat <- function() {
  description <- "Preclinical (rat, male Han Wistar). PBPK (semi-mechanistic dendritic-nanoparticle biodistribution). Dual-species model carrying nanoparticle-conjugated API and nanoparticle-released (free) API through four compartments each -- blood, liver, spleen and a lumped 'other'/rest compartment (8 ODEs). Structure is identical to the mouse anchor; species-specific parameters were obtained PROSPECTIVELY by scaling the mouse fit, not by re-fitting rat data: blood and rest volumes and the saturable-binding parameters by the blood-protein-binding ratio (Eqs 22 and 24), blood/rest flow by a BW^0.7 power law (Eq 23), and the conjugated-API extravasation rates by organ blood-flow ratios (Eq 25). Deterministic typical-value simulator: no IIV and no residual error are reported."
  reference <- "Vasalou C, Harding J, Jones RDO, Hariparsad N, McGinnity DF. Interspecies evaluation of a physiologically based pharmacokinetic model to predict the biodistribution dynamics of dendritic nanoparticles. PLoS ONE. 2023;18(5):e0285798. doi:10.1371/journal.pone.0285798. Rat parameters from Tables 4 and 5 and scaling Eqs 22-25; ODEs from Eqs 4-21; numeric values taken to full precision from the authors' MATLAB source (S2 File = ODE right-hand sides, S3 File = per-species parameter driver, species = 2)."
  vignette <- "Vasalou_2023_dendriticNanoparticle"
  units <- list(
    time = "h",
    dosing = "mg/kg",
    concentration = "ng/mL"
  )

  covariateData <- list()

  population <- list(
    species        = "rat (male Han Wistar)",
    n_subjects     = 9,
    n_studies      = 1,
    age_range      = "not reported",
    weight_range   = "0.25 kg (reference body weight, Tables 3 and 4)",
    sex_female_pct = 0,
    disease_state  = "Healthy male Han Wistar rats; no tumour implanted.",
    dose_range     = "55, 110 or 505 mg/kg of the nanoparticle as a 30-minute IV infusion (dose volume 10 mL/kg), given on day 1 and day 8; dose expressed as mg of API per kg body weight.",
    regions        = "Single-centre preclinical (AstraZeneca R&D Boston; animals from Charles River Laboratories).",
    notes          = "Blood sampled at 0.5, 1, 8, 24 and 72 h post dose on day 1; tissue collected 24 h post dose on day 8. Plasma sampling was dense but tissues were sparse -- ONE liver time point and NO spleen samples were obtained in rat, so model performance in rat can only be judged from plasma. Observed 55 mg/kg values are tabulated in Table S6 of S1 File. The paper reports that total API in plasma was captured accurately while released API in plasma was underpredicted by roughly 3-5-fold.",
    scope_note     = "Mechanistic PBPK simulator intended for typical-value simulation: no inter-individual variability and no residual-error model are reported by the source paper. Parameters were scaled prospectively from mouse (Eqs 22-25), NOT fitted to rat data. Dosing, volumes and flows are all body-weight-normalised (mg/kg, L/kg, L/h/kg), so simulate a nominal 1 kg subject and read concentrations directly."
  )

  ini({
    # ======================================================================
    # Nanoparticle-RELEASED (free) API disposition -- source paper Table 4.
    # Every value is fixed(): rat parameters are prospective scalings of the
    # mouse fit (Eqs 22-25), so nothing is estimated from rat data.
    # Values are quoted to the precision of the authors' S3 File driver
    # (species = 2); the Table 4 rounded value follows for cross-check.
    # Rat blood-protein-binding ratio used by Eqs 22 and 24 is
    # fu,blood(rat)/fu,blood(mouse) = (0.0006/0.84)/(0.0003/0.84) = 2.000.
    # ======================================================================
    lvc      <- fixed(log(0.6787))   ; label("Blood (central) volume of released API Vb (L/kg)")        # S3 parameters(5); Eq 22: 0.33937 * 2.000; Table 4 = 0.67
    lvp      <- fixed(log(0.2185))   ; label("Rest (peripheral) volume of released API VR (L/kg)")      # S3 parameters(8); Eq 22: 0.10923 * 2.000; Table 4 = 0.2
    # Eq 23 (BW^0.7 power law on the L/h flow):
    # (0.015265 * 0.02) * (0.25/0.02)^0.7 / 0.25 = 0.00716 -> S3 uses 0.0072.
    lq       <- fixed(log(0.0072))   ; label("Blood/rest flow rate QBR (L/h/kg)")                       # S3 parameters(3); Eq 23; Table 4 = 0.007
    # Nominal rat CL is the in-vivo blood CL from NCA of the unconjugated-API
    # IV profile (Table 2 = 2.1); the paper also simulates 2.8 (in-vitro
    # scaled) in Fig 5, and Fig 7 states 2.1 was used for the rat GOF plot.
    lcl      <- fixed(log(2.1))      ; label("Released API clearance CL (L/h/kg)")                      # S3 CL_rat = 2.1 (NCA); Table 2 in-vivo blood CL = 2.1
    v_liver  <- fixed(0.078)         ; label("Liver volume VL (L/kg)")                                  # S3 VL_rat; Table 4 = 0.078 (physiological, Davies et al [5])
    v_spleen <- fixed(0.0052)        ; label("Spleen volume VS (L/kg)")                                 # S3 VS_rat; Table 4 = 0.0052
    qbl_o    <- fixed(5.4)           ; label("Liver to blood flow rate QBL,o (L/h/kg)")                 # S3 QPLO_rat; Table 4 = 5.4 (physiological)
    qbs      <- fixed(0.15)          ; label("Blood/spleen flow rate QBS (L/h/kg)")                     # S3 QPS_rat; Table 4 = 0.15
    bmaxl    <- fixed(99)            ; label("Max binding capacity in liver BmaxL (ng/mL)")             # S3 parameters(33); Eq 24: 49.649 * 2.000; Table 4 = 99
    # KDL and KDS are held at the mouse values: "the dissociation constants
    # in liver and spleen, KDL and KDS, were considered constant across
    # species" (text under Eq 24; Table 4 prints "= mouse").
    kdl      <- fixed(0.59915)       ; label("Dissociation constant in liver KDL (ng/mL)")              # S3 KDL (mouse value, constant across species); Table 4 = "= mouse"
    pl       <- fixed(1.51)          ; label("Non-specific binding in liver PL (unitless)")             # S3 parameters(35); Eq 24: 0.75973 * 2.000; Table 4 = 1.5
    bmaxs    <- fixed(251.8)         ; label("Max binding capacity in spleen BmaxS (ng/mL)")            # S3 parameters(36); Eq 24: 125.91 * 2.000; Table 4 = 251.8
    kds      <- fixed(0.63238)       ; label("Dissociation constant in spleen KDS (ng/mL)")             # S3 KDS (mouse value, constant across species); Table 4 = "= mouse"
    ps       <- fixed(0.3384)        ; label("Non-specific binding in spleen PS (unitless)")            # S3 parameters(38); Eq 24: 0.16919 * 2.000; Table 4 = 0.33
    bpr      <- fixed(0.84)          ; label("Blood-to-plasma ratio of API BPR (unitless)")             # S3 BPR_rat; Table 4 = 0.84
    hct      <- fixed(0.45)          ; label("Hematocrit H (unitless)")                                 # S3 Hem = 0.45 (not reported in the paper text/tables)

    # ======================================================================
    # Nanoparticle-CONJUGATED API disposition -- source paper Table 5.
    # ======================================================================
    vnb          <- fixed(0.054)      ; label("Blood volume for conjugated API VNb (L/kg)")             # S3 V1_rat; Table 5 = 0.054 (physiological blood volume)
    # Eq 25: extravasation rates scale by the organ blood-flow ratio.
    nbl          <- fixed(1.7645e-4)  ; label("Blood/liver extravasation rate NBL (L/h/kg)")            # S3 NPL_mouse*(QPLO_rat/QPLO_mouse) = 0.00029801*(5.4/9.12); Table 5 = 0.0002
    nbs          <- fixed(3.97823e-6) ; label("Blood/spleen extravasation rate NBS (L/h/kg)")           # S3 NPS_mouse*(QPS_rat/QPS_mouse) = 1.9732e-5*(0.15/0.744); Table 5 = 4e-6
    nbr          <- fixed(0)          ; label("Blood/rest distribution rate NBR (L/h/kg)")              # S3 parameters(19) = 0.0; Table 5 = 0 (footnote +)
    knbl         <- fixed(1.8)        ; label("Liver/blood partition coefficient KNBL (unitless)")      # S3 parameters(30); Eq 24: 0.89999 * 2.000; Table 5 = 1.8
    knbs         <- fixed(1000)       ; label("Spleen/blood partition coefficient KNBS (unitless)")     # S3 parameters(32); Table 5 = 1000 (footnote ++)
    # Release rate constants are held at the mouse values across all species:
    # "release rate constants across plasma, liver and spleen were assumed
    # identical across mouse, rat, dog and human" (Discussion).
    krel_b       <- fixed(0.1246)     ; label("API release rate in blood krelb (1/h)")                  # S3 KREL; Table 5 = 0.125, "= mouse"
    krel_l       <- fixed(0.016374)   ; label("API release rate in liver krelL (1/h)")                  # S3 KRLV; Table 5 = 0.016, "= mouse"
    krel_s       <- fixed(0.0063168)  ; label("API release rate in spleen krelS (1/h)")                 # S3 KRLS; Table 5 = 0.0063, "= mouse"
    # See the mouse file and the vignette Errata: Table 5 prints 0.125 /
    # 0.016 for these two, byte-identical to the krelb / krelL cells one
    # row-group above; the authors' S3 File values below are what generated
    # Figs 4-9. Operator-ratified 2026-07-29 (sidecar request-001 q2 = A).
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
