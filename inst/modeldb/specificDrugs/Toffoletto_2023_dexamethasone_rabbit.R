Toffoletto_2023_dexamethasone_rabbit <- function() {
  description <- paste(
    "Preclinical (rabbit). Physiology-based ocular PK model of dexamethasone",
    "released from a drug-eluting soft contact lens (Ross 2019 sandwich lens),",
    "with well-mixed tear film, lumped sclera-choroid, retina, vitreous humour",
    "and aqueous humour (plus stroma) compartments. Transport across each",
    "barrier is area x permeability x the difference of the partition-corrected",
    "concentrations; drug leaves by tear drainage, palpebral-conjunctival",
    "uptake, choroidal blood-flow clearance, uveoscleral outflow and",
    "aqueous-humour turnover. The lens is a first-order depot (Eqs. 12-13):",
    "a fraction F of the released drug crosses the cornea directly into the",
    "aqueous humour and the rest enters the tears (Eqs. 10-11). Deterministic",
    "(no IIV or residual error); every value is a literature or hand-tuned",
    "input, not a statistical estimate.",
    sep = " "
  )
  reference <- paste(
    "Toffoletto N, Saramago B, Serro AP, Chauhan A. A Physiology-Based",
    "Mathematical Model to Understand Drug Delivery from Contact Lenses to the",
    "Back of the Eye. Pharm Res. 2023;40(8):1939-1951.",
    "doi:10.1007/s11095-023-03560-7.",
    "Validation data: Ross AE et al. Topical sustained drug delivery to the",
    "retina with a drug-eluting contact lens. Biomaterials. 2019;217:119285.",
    "doi:10.1016/j.biomaterials.2019.119285.",
    sep = " "
  )
  vignette <- "Toffoletto_2023_dexamethasone_rabbit"
  units <- list(time = "day", dosing = "ug", concentration = "ng/mL")

  # The paper lumps sclera and choroid into one well-mixed state with a
  # weighted-average partition coefficient (Methods, Model Description); no
  # registered compartment covers the lumped tissue.
  paper_specific_compartments <- c("sclera_choroid")

  compartmentData <- list(
    depot = list(analyte = "dexamethasone", units = "ug", specimen = "administration site", verified = TRUE),
    tear = list(analyte = "dexamethasone", units = "ug", specimen = "tear fluid", verified = TRUE),
    sclera_choroid = list(analyte = "dexamethasone", units = "ug", specimen = "tissue", verified = TRUE),
    retina = list(analyte = "dexamethasone", units = "ug", specimen = "retina", verified = TRUE),
    vitreous = list(analyte = "dexamethasone", units = "ug", specimen = "vitreous", verified = TRUE),
    aqueous_humour = list(analyte = "dexamethasone", units = "ug", specimen = "aqueous humour", verified = TRUE)
  )

  population <- list(
    species = "rabbit",
    n_subjects = NA_integer_,
    n_studies = 1L,
    disease_state = "healthy",
    dose_range = paste(
      "One dexamethasone-eluting sandwich contact lens per eye, worn for 7",
      "days; about 1.05 mg released in vitro over the wear period",
      "(Supplementary Fig. S1)."
    ),
    notes = paste(
      "No animals were studied in Toffoletto 2023: the model is parameterised",
      "from rabbit physiology (Tables II and IV) and literature drug",
      "properties (Table III), and is compared against the published rabbit",
      "tissue concentrations of Ross 2019 (aqueous humour, choroid, sclera,",
      "retina and vitreous humour at 0.25, 1, 2, 5 and 7 days; 3-5 eyes per",
      "time in Fig. 3). Ross 2019 is not on disk, so its strain, sex and",
      "animal counts are not recorded here."
    )
  )

  ini({
    # Lens release (Eqs. 12-13): f = M_Released * (1 - exp(-t/T)), so the
    # release rate R = M_Released / T * exp(-t/T) is first-order release of a
    # dose M_Released from a depot with rate constant 1/T. The paper does not
    # print T or M_Released; both were digitised by the maintainers from the
    # dashed exponential fit in Supplementary Fig. S1 (right panel, 'Release
    # from CL'): M_Released = 1.047e6 ng, T = 2.55 days (fit RMSE 3e3 ng).
    lka <- log(1 / 2.55); label("Lens release rate constant, 1/T (1/day)") # Supplementary Fig. S1 right-panel fit, digitised: T = 2.55 day
    f_aqueous_humour <- 0.0208; label("Fraction of released drug entering the aqueous humour directly across the cornea, F (fraction)") # Table III 'F' = 2.08 (percent), calculated by Eq. 14

    # Drug-specific permeabilities (Table III, cm/s)
    p_conj <- fixed(2.5e-6); label("Conjunctival permeability P_Conj (cm/s)") # Table III 'PConj' = 2.50E-06 cm/s
    p_rpe <- fixed(10e-6); label("Retinal-pigment-epithelium permeability P_ScCh-Ret (cm/s)") # Table III 'PScCh-Ret' = 10E-06 cm/s
    p_ilm <- fixed(13e-6); label("Inner-limiting-membrane permeability P_Ret-Vit (cm/s)") # Table III 'PRet-Vit' = 13E-06 cm/s (30 percent above PScCh-Ret)

    # Partition coefficients, tissue with respect to buffer (unitless). KScCh
    # and KRet were adapted for dexamethasone and tuned by hand during the
    # validation (Table I footnote d); the others are assumed equal to buffer.
    kp_sclera_choroid <- 15; label("Sclera-choroid:buffer partition coefficient K_ScCh (unitless)") # Table III 'KScCh' = 15
    kp_retina <- 10; label("Retina:buffer partition coefficient K_Ret (unitless)") # Table III 'KRet' = 10
    kp_vitreous <- fixed(1); label("Vitreous:buffer partition coefficient K_Vit (unitless)") # Table IV 'KVit' = 1, assumed equal to buffer
    kp_tear <- fixed(1); label("Tear:buffer partition coefficient K_t (unitless)") # Table IV 'Kt' = 1, assumed equal to buffer

    # Rabbit-eye physiology (Table IV), flows in the paper's units
    q_uveoscleral <- fixed(0.176); label("Uveoscleral outflow Q_UvSc (uL/min)") # Table IV 'QUvSc' = 0.176 uL/min
    q_vitreous_aqueous_humour <- fixed(0.19); label("Vitreous-to-aqueous elimination flow Q_Vit-Aq (uL/min)") # Table IV 'QVit-Aq' = 0.19 uL/min
    q_aqueous_humour <- fixed(4.2); label("Aqueous-humour renovation rate Q_Aq (uL/min)") # Table IV 'QAq' = 4.2 uL/min
    q_tear_drainage <- fixed(0.5); label("Tear drainage Q_Drain (uL/min)") # Table IV 'QDrain' = 0.5 uL/min
    cl_sclera_choroid <- fixed(1); label("Sclera-choroid clearance by choroidal blood flow (mL/min)") # Table IV 'ClearanceScCh' = 1 mL/min

    # Rabbit-eye anatomy (Table II)
    v_tear <- fixed(7e-3); label("Tear volume V_t (mL)") # Table II 'Vt' = 7E-03 mL
    v_vitreous <- fixed(1.7); label("Vitreous-humour volume V_Vit (mL)") # Table II 'VVit' = 1.7 mL
    v_sclera_choroid <- fixed(0.361); label("Sclera-choroid volume V_ScCh (mL)") # Table II 'VScCh' = 0.361 mL
    v_retina <- fixed(0.086); label("Retina volume V_Ret (mL)") # Table II 'VRet' = 0.086 mL
    v_aqueous_humour <- fixed(0.325); label("Aqueous-humour volume V_Aq (mL)") # Table II 'VAq' = 0.325 mL
    a_globe <- fixed(8.6); label("Area of the eye globe A_Globe (cm^2)") # Table II 'AGlobe' = 8.6 cm2
    a_palpebral <- fixed(14); label("Exposed area of palpebral conjunctiva A_Palp (cm^2)") # Table II 'APalp' = 14 cm2
    a_bulbar <- fixed(3); label("Exposed area of bulbar conjunctiva A_Bulb (cm^2)") # Table II 'ABulb' = 3 cm2
  })

  model({
    ka <- exp(lka)

    # Barrier conductances: cm/s x cm^2 = mL/s; x 86400 s/day gives mL/day
    ps_palpebral <- p_conj * a_palpebral * 86400
    ps_bulbar <- p_conj * a_bulbar * 86400
    ps_rpe <- p_rpe * a_globe * 86400
    ps_ilm <- p_ilm * a_globe * 86400

    # Flows: uL/min x 1440 min/day / 1000 uL/mL = mL/day; mL/min x 1440
    q_uv <- q_uveoscleral * 1.44
    q_va <- q_vitreous_aqueous_humour * 1.44
    q_aq <- q_aqueous_humour * 1.44
    q_drain <- q_tear_drainage * 1.44
    cl_sc <- cl_sclera_choroid * 1440

    # Concentrations in ug/mL
    c_tear <- tear / v_tear
    c_sc <- sclera_choroid / v_sclera_choroid
    c_ret <- retina / v_retina
    c_vit <- vitreous / v_vitreous
    c_aq <- aqueous_humour / v_aqueous_humour

    # Eq. 13 release rate R = ka * depot
    release <- ka * depot

    # Transport terms (Fig. 2 arrows): d bulbar conjunctiva, h RPE, i ILM
    flux_bulbar <- ps_bulbar * (c_tear / kp_tear - c_sc / kp_sclera_choroid)
    flux_rpe <- ps_rpe * (c_sc / kp_sclera_choroid - c_ret / kp_retina)
    flux_ilm <- ps_ilm * (c_ret / kp_retina - c_vit / kp_vitreous)

    d/dt(depot) <- -release
    # Eq. 10 (tears)
    d/dt(tear) <- (1 - f_aqueous_humour) * release - q_drain * c_tear - ps_palpebral * c_tear - flux_bulbar
    # Eq. 3 (sclera-choroid); the choroidal clearance acts on C_ScCh as printed
    d/dt(sclera_choroid) <- q_uv * c_aq - q_uv * c_sc / kp_sclera_choroid + flux_bulbar - cl_sc * c_sc - flux_rpe
    # Eq. 4 (retina)
    d/dt(retina) <- flux_rpe - flux_ilm
    # Eq. 5 (vitreous humour)
    d/dt(vitreous) <- flux_ilm - q_va * c_vit
    # Eq. 11 (aqueous humour, stroma lumped in)
    d/dt(aqueous_humour) <- q_va * c_vit + f_aqueous_humour * release - q_aq * c_aq

    # Outputs in ng/mL (ug/mL x 1000), read as ng/g tissue at unit density
    Ctear <- 1000 * c_tear
    Csclera_choroid <- 1000 * c_sc
    Cretina <- 1000 * c_ret
    Cvitreous <- 1000 * c_vit
    Caqueous_humour <- 1000 * c_aq
  })
}
