Vasalou_2023_dendriticNanoparticle_human <- function() {
  description <- "PBPK (semi-mechanistic dendritic-nanoparticle biodistribution), human first-in-human PROJECTION -- no human data were fitted. Dual-species model carrying nanoparticle-conjugated API and nanoparticle-released (free) API through four compartments each -- blood, liver, spleen and a lumped 'other'/rest compartment (8 ODEs). Structure is identical to the mouse anchor; human parameters were obtained entirely by prospective scaling of the mouse fit (Eqs 22-25), with the released-API clearance set to 0.32 L/h/kg from the in-vitro human hepatocyte Clint of 15.5 uL/min/1e6 cells via the well-stirred model. Intended as a translational / clinical-projection tool; the simulated dose-normalised profiles overlay the observed mouse, rat and dog data (Fig 9). Deterministic typical-value simulator: no IIV and no residual error are reported."
  reference <- "Vasalou C, Harding J, Jones RDO, Hariparsad N, McGinnity DF. Interspecies evaluation of a physiologically based pharmacokinetic model to predict the biodistribution dynamics of dendritic nanoparticles. PLoS ONE. 2023;18(5):e0285798. doi:10.1371/journal.pone.0285798. Human parameters from Tables 3, 4 and 5 and scaling Eqs 22-25; ODEs from Eqs 4-21; numeric values taken to full precision from the authors' MATLAB source (S2 File = ODE right-hand sides, S3 File = per-species parameter driver, species = 4)."
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
    blood     = list(analyte = "released-API (free API)", units = NA_character_, specimen = "blood cell", verified = FALSE),
    liver     = list(analyte = "released-API (free API)", units = NA_character_, specimen = "tissue", verified = FALSE),
    spleen    = list(analyte = "released-API (free API)", units = NA_character_, specimen = "tissue", verified = FALSE),
    other     = list(analyte = "released-API (free API)", units = NA_character_, specimen = "administration site", verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 0,
    n_studies      = 0,
    age_range      = "not applicable (simulation only)",
    weight_range   = "70 kg (reference body weight, Tables 3 and 4)",
    sex_female_pct = NA_real_,
    disease_state  = "Not applicable. No human subjects were studied; this is a prospective clinical projection built by scaling the preclinical model.",
    dose_range     = "Simulated 10 mg/kg of the nanoparticle; the authors' S3 File administers it as a 3-hour IV infusion. Model-derived profiles were dose-normalised for comparison against preclinical data (Fig 9).",
    regions        = "Not applicable (in silico).",
    notes          = "IMPORTANT: this file contains NO fitted human parameters and NO human observations exist to validate it. Physiological inputs (liver and spleen volumes, hepatic and splenic blood flows, blood volume, hepatocellularity, liver weight) come from Tables 3-5 and Davies et al [5]. The released-API clearance is the only species-specific non-physiological input and was predicted from in-vitro human hepatocyte Clint = 15.5 uL/min/1e6 cells via Eqs 1-2 (well-stirred model with a regression offset of 3), giving 0.32 L/h/kg -- roughly 5-7-fold lower than the preclinical in-vivo values. Because the paper showed released API in plasma was already underpredicted 3-5-fold in rat and dog at their nominal CL, treat simulated human released-API concentrations as a lower bound and run a CL sensitivity analysis as the paper recommends.",
    scope_note     = "Mechanistic PBPK simulator intended for typical-value simulation and clinical translation only: no inter-individual variability and no residual-error model are reported by the source paper, and no human data were fitted. Dosing, volumes and flows are all body-weight-normalised (mg/kg, L/kg, L/h/kg), so simulate a nominal 1 kg subject and read concentrations directly."
  )

  ini({
    # ======================================================================
    # Nanoparticle-RELEASED (free) API disposition -- source paper Table 4.
    # Every value is fixed(): human parameters are prospective scalings of
    # the mouse fit (Eqs 22-25) plus in-vitro-predicted clearance. Nothing
    # was estimated from human data because no human data exist.
    # Values are quoted to the precision of the authors' S3 File driver
    # (species = 4); the Table 4 rounded value follows for cross-check.
    # Human blood-protein-binding ratio used by Eqs 22 and 24 is
    # fu,blood(human)/fu,blood(mouse) = (0.0003/0.766)/(0.0003/0.84) = 1.0966.
    # ======================================================================
    lvc      <- fixed(log(0.3722))   ; label("Blood (central) volume of released API Vb (L/kg)")        # S3 parameters(5); Eq 22: 0.33937 * 1.0966; Table 4 = 0.37
    lvp      <- fixed(log(0.1198))   ; label("Rest (peripheral) volume of released API VR (L/kg)")      # S3 parameters(8); Eq 22: 0.10923 * 1.0966; Table 4 = 0.12
    # Eq 23 (BW^0.7 power law on the L/h flow):
    # (0.015265 * 0.02) * (70/0.02)^0.7 / 70 = 0.00132 -> S3 uses 0.0013.
    lq       <- fixed(log(0.0013))   ; label("Blood/rest flow rate QBR (L/h/kg)")                       # S3 parameters(3); Eq 23; Table 4 = 0.0013
    # Human CL is NOT measured in vivo. It is predicted from the in-vitro
    # human hepatocyte Clint of 15.5 uL/min/1e6 cells through Eqs 1-2, on
    # the basis of the interspecies IVIVC established in Table 2 (Results,
    # "the released API human CL was set to 0.32 L/h/Kg").
    lcl      <- fixed(log(0.324))    ; label("Released API clearance CL (L/h/kg)")                      # S3 CL_human = 5.4/1000*60 = 0.324; Table 4 = 0.32
    v_liver  <- fixed(0.0241)        ; label("Liver volume VL (L/kg)")                                  # S3 VL_human; Table 4 = 0.024 (physiological, Davies et al [5])
    v_spleen <- fixed(0.0027)        ; label("Spleen volume VS (L/kg)")                                 # S3 VS_human; Table 4 = 0.0027
    # S3 carries 1.2429 while Table 4 rounds to 1.3; the S3 value is used
    # because it is what generated Fig 9. Consequently the derived QBL,i is
    # 1.2429 - 0.066 = 1.1769 rather than Table 4's rounded 1.23.
    qbl_o    <- fixed(1.2429)        ; label("Liver to blood flow rate QBL,o (L/h/kg)")                 # S3 QPLO_human; Table 4 = 1.3 (physiological)
    qbs      <- fixed(0.066)         ; label("Blood/spleen flow rate QBS (L/h/kg)")                     # S3 QPS_human; Table 4 = 0.066
    bmaxl    <- fixed(54.4)          ; label("Max binding capacity in liver BmaxL (ng/mL)")             # S3 parameters(33); Eq 24: 49.649 * 1.0966; Table 4 = 54
    # KDL and KDS are held at the mouse values: "the dissociation constants
    # in liver and spleen, KDL and KDS, were considered constant across
    # species" (text under Eq 24; Table 4 prints "= mouse").
    kdl      <- fixed(0.59915)       ; label("Dissociation constant in liver KDL (ng/mL)")              # S3 KDL (mouse value, constant across species); Table 4 = "= mouse"
    pl       <- fixed(0.83)          ; label("Non-specific binding in liver PL (unitless)")             # S3 parameters(35); Eq 24: 0.75973 * 1.0966; Table 4 = 0.83
    bmaxs    <- fixed(138)           ; label("Max binding capacity in spleen BmaxS (ng/mL)")            # S3 parameters(36); Eq 24: 125.91 * 1.0966; Table 4 = 138
    kds      <- fixed(0.63238)       ; label("Dissociation constant in spleen KDS (ng/mL)")             # S3 KDS (mouse value, constant across species); Table 4 = "= mouse"
    ps       <- fixed(0.1855)        ; label("Non-specific binding in spleen PS (unitless)")            # S3 parameters(38); Eq 24: 0.16919 * 1.0966; Table 4 = 0.19
    bpr      <- fixed(0.766)         ; label("Blood-to-plasma ratio of API BPR (unitless)")             # S3 BPR_human; Table 4 = 0.77
    hct      <- fixed(0.45)          ; label("Hematocrit H (unitless)")                                 # S3 Hem = 0.45 (not reported in the paper text/tables)

    # ======================================================================
    # Nanoparticle-CONJUGATED API disposition -- source paper Table 5.
    # ======================================================================
    vnb          <- fixed(0.0743)     ; label("Blood volume for conjugated API VNb (L/kg)")             # S3 V1_human; Table 5 = 0.074 (physiological blood volume)
    # Eq 25: extravasation rates scale by the organ blood-flow ratio. The
    # Discussion notes this makes human extravasation about 10-fold slower
    # than mouse, consistent with published nanoparticle PBPK models.
    nbl          <- fixed(4.06136e-5) ; label("Blood/liver extravasation rate NBL (L/h/kg)")            # S3 NPL_mouse*(QPLO_human/QPLO_mouse) = 0.00029801*(1.2429/9.12); Table 5 = 0.00004
    nbs          <- fixed(1.75042e-6) ; label("Blood/spleen extravasation rate NBS (L/h/kg)")           # S3 NPS_mouse*(QPS_human/QPS_mouse) = 1.9732e-5*(0.066/0.744); Table 5 = 1.7e-6
    nbr          <- fixed(0)          ; label("Blood/rest distribution rate NBR (L/h/kg)")              # S3 parameters(19) = 0.0; Table 5 = 0 (footnote +)
    knbl         <- fixed(0.9869)     ; label("Liver/blood partition coefficient KNBL (unitless)")      # S3 parameters(30); Eq 24: 0.89999 * 1.0966; Table 5 = 1
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
