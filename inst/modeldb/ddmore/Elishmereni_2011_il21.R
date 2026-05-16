# Pharmacokinetic model for recombinant murine interleukin-21 (rIL-21) in
# mice across three modes of administration (intraperitoneal, intravenous,
# subcutaneous), extracted from the DDMORE Foundation Model Repository
# entry DDMODEL00000230. The DDMORE deposit is a Monolix 3.2 re-fit of a
# simplified (one-compartment-removed) version of the model published in
# Elishmereni et al. 2011 (PLoS Comp Biol). The fit is in mice; no human
# scaling is included.

Elishmereni_2011_il21 <- function() {
  description <- "Preclinical (mouse). Seven-compartment population PK model for recombinant murine interleukin-21 (rIL-21) administered by intraperitoneal (ip), intravenous (iv), or subcutaneous (sc) routes, extracted from the DDMORE Foundation Model Repository (DDMODEL00000230). Linear elimination from the central compartment plus a direct first-pass loss from the ip depot, parallel ip and sc absorption routes each with a single intermediate transit compartment, and two peripheral distribution compartments. Saturable transfer terms present in the published model (sa0, sa1, sa2) and the sc-depot direct elimination rate (k30) are held at zero in the DDMORE deposit, leaving an effectively linear system. Per-subject random effects on CL, central volume, and the ip and sc absorption rates are correlated as a 4 x 4 block. Dose route is encoded by the cmt column: ip = depot, iv = central, sc = depot2. Units are not declared in the DDMORE bundle (no IL_21_PK.csv shipped); see the validation vignette Errata."
  reference <- paste(
    "Elishmereni M, Kheifetz Y, Sondergaard H, Overgaard RV, Agur Z (2011).",
    "An integrated disease/pharmacokinetic/pharmacodynamic model suggests",
    "improved interleukin-21 regimens validated prospectively for mouse",
    "solid cancers.",
    "PLoS Computational Biology 7(9):e1002206.",
    "doi:10.1371/journal.pcbi.1002206.",
    "DDMORE Foundation Model Repository: DDMODEL00000230",
    "(simplified Monolix 3.2 re-fit; one compartment removed and",
    "parameters re-estimated by mixed-effects, per Model_Accommodations.txt).",
    sep = " "
  )
  vignette <- "Elishmereni_2011_il21"
  units <- list(
    time          = "unspecified (assumed h)",
    dosing        = "unspecified",
    concentration = "unspecified"
  )

  ddmore_id    <- "DDMODEL00000230"
  replicate_of <- NULL

  covariateData <- list()

  population <- list(
    species        = "mouse (strain not declared in DDMORE bundle; Elishmereni 2011 used C57BL/6 for the in vivo IL-21 PK experiments)",
    n_subjects     = NA_integer_,
    n_studies      = NA_integer_,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    disease_state  = "Tumour-bearing and tumour-free mice receiving recombinant murine interleukin-21 (rIL-21) by ip, iv, or sc routes; the DDMORE deposit description states the model encompasses melanoma and renal cell carcinoma settings.",
    dose_range     = "Multiple administration routes (ip, iv, sc); specific dose levels are not reproduced in the DDMORE bundle (the IL_21_PK.csv source dataset is not shipped).",
    regions        = NA_character_,
    notes          = "The DDMORE deposit (DDMODEL00000230) is a simplified Monolix 3.2 re-fit of the original Elishmereni 2011 PK model: one compartment was removed and all parameters re-estimated by mixed-effects (the original paper used step-wise least-squares fitting). Per the deposit's RDF metadata: model-implementation-conforms-to-literature-controlled = 'No'; model-origin-of-code-in-literature-controlled = 'No'. The original Elishmereni 2011 publication is not on disk in this worktree; parameter values and equations were taken verbatim from the DDMORE bundle's IL_21_PK_model.mdl. See the validation vignette Errata for the full list of bundle-versus-publication caveats."
  )

  ini({
    # Final estimates from the DDMORE bundle's IL_21_PK_model.mdl PARAMETER
    # OBJECT (parObj), STRUCTURAL block. The DDMORE deposit ships the .mdl
    # and its rendered PharmML (.xml) as the authoritative encoding of the
    # Monolix 3.2 fit; no NONMEM .lst is shipped (the bundle is a Monolix
    # deposit, not a NONMEM run). The .mdl `value = ...` entries in a
    # DDMORE deposit are by convention the final fitted point estimates
    # rather than initial values.
    #
    # Several published-model parameters are held at zero in this DDMORE
    # version of the model: the saturable-transfer scalars sa0, sa1, sa2
    # and the sc-depot direct elimination rate k30. Holding them at zero
    # makes the corresponding terms drop out of the ODE system, so they
    # are omitted from ini() / model() rather than carried as fixed(0).
    # Two further structural parameters (POP_v1 = 0.001 ip-depot volume
    # and POP_v3 = 0.001 sc-depot volume) are present in the .mdl but
    # never appear in any rate equation; they are also omitted.

    # Central elimination and central volume (the central compartment is
    # the source's A2; v2 is the central volume).
    lcl  <- log(0.0229)  ; label("Central clearance CL (volume/time; units unspecified by DDMORE deposit)")  # IL_21_PK_model.mdl POP_cl   = 0.0229
    lvc  <- log(0.00551) ; label("Central volume V2 (volume; units unspecified)")                            # POP_v2   = 0.00551

    # Peripheral volumes V4 and V5 (the source's two peripheral
    # distribution compartments). Note that V5 (24.4) is several orders
    # of magnitude larger than V4 (0.0009), reflecting a deep tissue
    # distribution compartment in the source's parameterisation.
    lvp  <- log(0.0009)  ; label("Peripheral volume 1 V4 (volume; units unspecified)")                       # POP_v4   = 0.0009
    lvp2 <- log(24.4)    ; label("Peripheral volume 2 V5 (volume; units unspecified)")                       # POP_v5   = 24.4

    # Absorption / transit rate constants. The source labels these as
    # "Q" parameters in the parObj but treats them as first-order rate
    # constants in MODEL_PREDICTION (k12 = q12; k32 = q23; k24 = q24;
    # k25 = q25). Per-route mapping:
    #   q12 -> k_ip2c  (ip-depot/transit -> central; rate constant 1/time)
    #   q23 -> k_sc2c  (sc-depot/transit -> central)
    #   q24 -> k_c2p1  (central -> peripheral 1; the return rate uses
    #                   k_p12c = k_c2p1 * v2/v4 to preserve mass balance)
    #   q25 -> k_c2p2  (central -> peripheral 2; analogously k_p22c)
    lk_ip2c <- log(0.693) ; label("ip-route -> central rate constant k_ip2c = q12 (1/time)")                 # POP_q12  = 0.693
    lk_sc2c <- log(0.727) ; label("sc-route -> central rate constant k_sc2c = q23 (1/time)")                 # POP_q23  = 0.727
    lk_c2p1 <- log(0.48)  ; label("Central -> peripheral 1 rate constant k_c2p1 = q24 (1/time)")             # POP_q24  = 0.48
    lk_c2p2 <- log(6.38)  ; label("Central -> peripheral 2 rate constant k_c2p2 = q25 (1/time)")             # POP_q25  = 6.38

    # Direct first-pass elimination from the ip depot (no analogous loss
    # from the sc depot; the source's k30 is held at zero).
    lk_ipelim <- log(0.4) ; label("ip-depot first-pass elimination rate k10 (1/time)")                       # POP_k10  = 0.4

    # Inter-individual variability. The .mdl VARIABILITY block declares
    # all sd values as `type is sd`; non-zero IIV is reported on cl, v2,
    # q12 and q23 only. The OMEGA correlation matrix (`type is corr`,
    # parameter = [ETA_cl, ETA_v2, ETA_q12, ETA_q23]) gives the lower-
    # triangular off-diagonals; covariances below are computed as
    # cov_ij = corr_ij * sd_i * sd_j.
    #   sd_cl  = 0.407   var_cl  = 0.165649
    #   sd_v2  = 0.74    var_vc  = 0.5476
    #   sd_q12 = 0.235   var_q12 = 0.055225
    #   sd_q23 = 0.304   var_q23 = 0.092416
    # Lower-triangular correlations (.mdl row-major order):
    #   corr(v2,  cl ) = -0.27
    #   corr(q12, cl ) =  0.10
    #   corr(q12, v2 ) = -0.29
    #   corr(q23, cl ) =  0.12
    #   corr(q23, v2 ) = -0.28
    #   corr(q23, q12) =  0.10
    etalcl + etalvc + etalk_ip2c + etalk_sc2c ~ c(
      0.165649,
      -0.081303,  0.547600,
       0.009565, -0.050443,  0.055225,
       0.014847, -0.062989,  0.007144,  0.092416
    )

    # Residual error. The .mdl OBSERVATION block uses
    # proportionalError(proportional = b, eps = EPS_Y, prediction = output1)
    # with var(EPS_Y) = 1, so b is a proportional residual error scalar
    # that maps directly onto nlmixr2 propSd.
    propSd <- 0.00975 ; label("Proportional residual error (fraction)")                                       # b        = 0.00975
  })

  model({
    # Compartment ordering is chosen so that the implicit positional
    # numbering matches the DDMORE source's CMT codes for dosing:
    #   CMT = 1 -> depot   (ip dose)
    #   CMT = 2 -> central    (iv dose)
    #   CMT = 3 -> depot2   (sc dose)
    #   CMT = 4 -> peripheral1
    #   CMT = 5 -> peripheral2
    #   CMT = 6 -> transit1
    #   CMT = 7 -> transit2
    # The source's original A1..A7 mapping is preserved one-to-one.

    # Individual parameters
    cl      <- exp(lcl  + etalcl)
    vc      <- exp(lvc  + etalvc)
    vp     <- exp(lvp)
    vp2     <- exp(lvp2)
    k_ip2c  <- exp(lk_ip2c + etalk_ip2c)
    k_sc2c  <- exp(lk_sc2c + etalk_sc2c)
    k_c2p1  <- exp(lk_c2p1)
    k_c2p2  <- exp(lk_c2p2)
    k_ipelim <- exp(lk_ipelim)

    # Micro-rate constants. Forward rates k_c2p1, k_c2p2 are first-order
    # (paper uses Q-symbols but the equations treat them as rate
    # constants); return rates are scaled by volume ratios to preserve
    # k_ij * V_i = k_ji * V_j mass-balance.
    k_p12c <- k_c2p1 * vc / vp
    k_p22c <- k_c2p2 * vc / vp2
    kel    <- cl / vc

    # ODE system. Verbatim translation of MODEL_PREDICTION DEQ from
    # IL_21_PK_model.mdl (lines 160-168) with sa0 = sa1 = sa2 = 0 and
    # k30 = 0, so all (1 + A_i * sa) saturation factors collapse to 1.
    # The published 7-compartment topology is preserved: ip and sc
    # routes each pass through a single intermediate "transit" or
    # delay compartment (A6, A7 in the source) on the way to central
    # (A2), with a parallel direct elimination from the ip depot.
    d/dt(depot)    <- -k_ip2c * depot - k_ipelim * depot                              # source A1 dA1/dt
    d/dt(central)     <-  k_ip2c * transit1 - kel * central +                                # source A2 dA2/dt
                          k_sc2c * transit2 -
                          k_c2p1 * central + k_p12c * peripheral1 -
                          k_c2p2 * central + k_p22c * peripheral2
    d/dt(depot2)    <- -k_sc2c * depot2                                                    # source A3 dA3/dt (k30 = 0 dropped)
    d/dt(peripheral1) <-  k_c2p1 * central - k_p12c * peripheral1                              # source A4 dA4/dt
    d/dt(peripheral2) <-  k_c2p2 * central - k_p22c * peripheral2                              # source A5 dA5/dt
    d/dt(transit1)  <-  k_ip2c * depot - k_ip2c * transit1 - k_ipelim * transit1      # source A6 dA6/dt
    d/dt(transit2)  <-  k_sc2c * depot2 - k_sc2c * transit2                              # source A7 dA7/dt (k30 = 0 dropped)

    # Observation: free IL-21 in central, mapped onto the source's
    # output1 = A2 / v2 .mdl line 169.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
