Bram_2026_warfarin_node <- function() {
  description <- "Neural-ODE (NODE) warfarin PK/PD model. One-compartment oral PK with an absorption lag time drives an indirect-response model for prothrombin complex activity (PCA), where the PRODUCTION term is not a mechanistic function but a trained 5-neuron softplus neural network of the plasma concentration Cc. This is the fitted low-dimensional NODE of Braem 2026 Equation 16, the starting point of that paper's automated NODE-LASSO model-development workflow: the network is fitted first, its derivative-versus-state behaviour is then read out and distilled by LASSO regression into the explicit Imax-plus-exponential structural model of Equation 24. Only the NODE has published parameter estimates (deposited Monolix run Warf_node_mlx_file_ind), so only the NODE is packaged here; Equation 24 is documented in the vignette but its estimates are printed nowhere. All 21 parameters carry inter-individual variability, log-normal on the six mechanistic parameters and NORMAL (additive) on the 15 network weights and biases, which take negative values. Because the network replaces the production rate outright, the drug-free state is NOT a steady state: at Cc = 0 the network returns 2.31 per hour against a loss of kout * PCA0 = 3.88 per hour, so PCA declines from baseline even without drug. That asymmetry is precisely what the paper's Equation 24 rearrangement (kin = k + q) was introduced to repair."
  reference <- paste(
    "Braem DS, Steiert B, Steffens B, Pfister M, Koch G.",
    "Automated pharmacometric model development by leveraging",
    "low-dimensional neural ODEs and LASSO regression.",
    "CPT Pharmacometrics Syst Pharmacol. 2026.",
    "doi:10.1002/psp4.70285.",
    "Parameter estimates are from the article's Data S2 code deposit",
    "(Monolix project Warf_node_mlx_file_ind, populationParameters.txt);",
    "no parameter table appears in the article text.",
    "Warfarin data from O'Reilly RA, Aggeler PM. Circulation.",
    "1968;38(1):169-177. doi:10.1161/01.cir.38.1.169.",
    sep = " "
  )
  vignette <- "Bram_2026_node_lasso_model_development"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "warfarin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "warfarin", units = "mg", specimen = "plasma", verified = TRUE),
    PCA = list(analyte = "prothrombin complex", units = "% of normal activity", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 32,
    n_studies = 1,
    disease_state = "Adults initiating oral anticoagulation; the warfarin PK/PD dataset distributed as an example with Monolix and nlmixr, originating from O'Reilly and Aggeler 1968.",
    dose_range = "Single oral dose, 60-153 mg (mean 105 mg) per the deposited warfarin_pkpd_data.csv; approximately 1.5 mg/kg.",
    observations = "Median 6 (IQR 6-10) warfarin plasma concentrations and 7 (IQR 7-8) prothrombin complex activity observations per subject over 144 h (Braem 2026 Section 2.3.3).",
    notes = "Braem 2026 reports no baseline demographic table for this dataset; age, weight, sex and race are not available in the deposited data file. The cohort size and observation counts are taken from Section 2.3.3 of the article."
  )

  ini({
    # ---------------------------------------------------------------
    # Mechanistic PK/PD parameters. Every value is the final population
    # estimate from the deposited Monolix run
    # Warf_node_mlx_file_ind/populationParameters.txt (Data S2). This is
    # the second of the two runs in the deposited pipeline (population
    # fit first without IIV on the network, then this run with it), and
    # is the fit whose MARE appears as MARE(NODE) in Braem 2026 Table 1.
    # ---------------------------------------------------------------
    ltlag <- log(0.808327021863975); label("Absorption lag time (h)")      # Data S2 Tlag_pop
    lka <- log(1.24103750716188); label("Absorption rate constant (1/h)")  # Data S2 ka_pop
    lvc <- log(8.00330109670927); label("Central volume of distribution (L)")  # Data S2 V_pop
    lcl <- log(0.130798108892028); label("Clearance (L/h)")                # Data S2 Cl_pop
    lrbase <- log(95.1494895315056); label("Baseline prothrombin complex activity (% of normal)")  # Data S2 R0_pop
    lkout <- log(0.04074788557565); label("PCA turnover elimination rate constant (1/h)")  # Data S2 kout_pop

    # ---------------------------------------------------------------
    # Neural network on the PCA production rate (Braem 2026 Eq 16
    # fNN_Cc; expanded form in the deposit's Warf_node_converted.txt).
    # One hidden layer, 5 neurons, softplus activation:
    #   h_j    = (1 / beta) * log(1 + exp(beta * (-(w1_j^2) * Cc) + b1_j))
    #   NN(Cc) = sum_j w2_j * h_j
    # The input-to-hidden weight is SQUARED and NEGATED, which is how
    # pmxNODE enforces the monotone-decreasing restriction requested by
    # `time_nn = TRUE` in Warf_node.txt; the sign of w1_j is therefore
    # not identifiable and only its magnitude matters.
    # ---------------------------------------------------------------
    nn_w1_rc_1 <- -0.831182806282164; label("Input-to-hidden weight, neuron 1 (unitless)")  # Data S2 Wrc_11_pop
    nn_w1_rc_2 <- 0.847738875709482; label("Input-to-hidden weight, neuron 2 (unitless)")   # Data S2 Wrc_12_pop
    nn_w1_rc_3 <- -0.734778745543524; label("Input-to-hidden weight, neuron 3 (unitless)")  # Data S2 Wrc_13_pop
    nn_w1_rc_4 <- -0.0329204148187805; label("Input-to-hidden weight, neuron 4 (unitless)") # Data S2 Wrc_14_pop
    nn_w1_rc_5 <- 0.388797941794558; label("Input-to-hidden weight, neuron 5 (unitless)")   # Data S2 Wrc_15_pop

    nn_b1_rc_1 <- 3.49729824388274; label("Hidden-layer bias, neuron 1 (unitless)")   # Data S2 brc_11_pop
    nn_b1_rc_2 <- 1.89601613209086; label("Hidden-layer bias, neuron 2 (unitless)")   # Data S2 brc_12_pop
    nn_b1_rc_3 <- 3.00836503572554; label("Hidden-layer bias, neuron 3 (unitless)")   # Data S2 brc_13_pop
    nn_b1_rc_4 <- 0.674859444377844; label("Hidden-layer bias, neuron 4 (unitless)")  # Data S2 brc_14_pop
    nn_b1_rc_5 <- 24.3761239925392; label("Hidden-layer bias, neuron 5 (unitless)")   # Data S2 brc_15_pop

    nn_w2_rc_1 <- 1.20105606854911; label("Hidden-to-output weight, neuron 1 (% of normal per h)")     # Data S2 Wrc_21_pop
    nn_w2_rc_2 <- -0.443457148068339; label("Hidden-to-output weight, neuron 2 (% of normal per h)")   # Data S2 Wrc_22_pop
    nn_w2_rc_3 <- -0.0864056336436934; label("Hidden-to-output weight, neuron 3 (% of normal per h)")  # Data S2 Wrc_23_pop
    nn_w2_rc_4 <- 0.249689692646057; label("Hidden-to-output weight, neuron 4 (% of normal per h)")    # Data S2 Wrc_24_pop
    nn_w2_rc_5 <- 1.75547754212371; label("Hidden-to-output weight, neuron 5 (% of normal per h)")     # Data S2 Wrc_25_pop

    nn_beta_rc <- fixed(20); label("Softplus sharpness (unitless)")  # Braem 2026 Supporting Information Eq (2); pmxNODE default

    # ---------------------------------------------------------------
    # IIV. Monolix reports omega as a STANDARD DEVIATION; nlmixr2 takes a
    # VARIANCE, so every entry below is the deposited omega squared.
    # logNormal in Monolix (Cl, R0, Tlag, V, ka, kout) -> multiplicative
    # eta on the log scale, i.e. the usual exp(l<p> + eta) form.
    # ---------------------------------------------------------------
    etaltlag ~ 0.135884619  # Data S2 omega_Tlag = 0.368625309283149, squared
    etalka ~ 0.727167650  # Data S2 omega_ka = 0.852623975404321, squared
    etalvc ~ 0.046649532  # Data S2 omega_V = 0.215985025386417, squared
    etalcl ~ 0.073129713  # Data S2 omega_Cl = 0.270425063129657, squared
    etalrbase ~ 0.000795364  # Data S2 omega_R0 = 0.0282021949253777, squared
    etalkout ~ 0.005966186  # Data S2 omega_kout = 0.0772410892003664, squared

    # normal (not logNormal) in Monolix for the 15 network parameters,
    # because weights and biases take negative values -> ADDITIVE eta on
    # the natural scale.
    etann_w1_rc_1 ~ 0.399789455  # Data S2 omega_Wrc_11 = 0.632289068423692, squared
    etann_w1_rc_2 ~ 0.017284918  # Data S2 omega_Wrc_12 = 0.131472120409786, squared
    etann_w1_rc_3 ~ 0.078655468  # Data S2 omega_Wrc_13 = 0.280455833481093, squared
    etann_w1_rc_4 ~ 0.078470845  # Data S2 omega_Wrc_14 = 0.280126500848085, squared
    etann_w1_rc_5 ~ 0.001302042  # Data S2 omega_Wrc_15 = 0.0360838014002875, squared

    etann_b1_rc_1 ~ 0.150360302  # Data S2 omega_brc_11 = 0.38776370899177, squared
    etann_b1_rc_2 ~ 0.138665843  # Data S2 omega_brc_12 = 0.37237863689953, squared
    etann_b1_rc_3 ~ 2.324411637  # Data S2 omega_brc_13 = 1.52460212360593, squared
    etann_b1_rc_4 ~ 0.647774701  # Data S2 omega_brc_14 = 0.804844521853478, squared
    etann_b1_rc_5 ~ 0.267902460  # Data S2 omega_brc_15 = 0.517592942787506, squared

    etann_w2_rc_1 ~ 0.011190963  # Data S2 omega_Wrc_21 = 0.105787363578368, squared
    etann_w2_rc_2 ~ 0.197472715  # Data S2 omega_Wrc_22 = 0.444379036542218, squared
    etann_w2_rc_3 ~ 0.097629996  # Data S2 omega_Wrc_23 = 0.312457971500172, squared
    etann_w2_rc_4 ~ 0.008916540  # Data S2 omega_Wrc_24 = 0.0944274462726303, squared
    etann_w2_rc_5 ~ 0.395331612  # Data S2 omega_Wrc_25 = 0.62875402067328, squared

    # ---------------------------------------------------------------
    # Residual error. The deposited .mlxtran declares
    # errorModel = combined1(a, b) on BOTH outputs, i.e. SD = a + b * f,
    # which is nlmixr2's add() + prop() + combined1().
    # ---------------------------------------------------------------
    addSd <- 2.22044604925031e-16; label("Additive residual SD on Cc (mg/L)")  # Data S2 acp; collapsed onto the lower bound during estimation
    propSd <- 0.124023750458425; label("Proportional residual SD on Cc (fraction)")  # Data S2 bcp
    addSd_PCA <- 0.0453462873737628; label("Additive residual SD on PCA (% of normal)")  # Data S2 apca
    propSd_PCA <- 0.0973560412001303; label("Proportional residual SD on PCA (fraction)")  # Data S2 bpca
  })

  model({
    # 1. Individual PK parameters
    tlag <- exp(ltlag + etaltlag)
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl)
    kel <- cl / vc

    # 2. Individual PD parameters
    rbase <- exp(lrbase + etalrbase)
    kout <- exp(lkout + etalkout)

    # 3. Individual network parameters (additive etas; see ini()).
    w1_rc_1 <- nn_w1_rc_1 + etann_w1_rc_1
    w1_rc_2 <- nn_w1_rc_2 + etann_w1_rc_2
    w1_rc_3 <- nn_w1_rc_3 + etann_w1_rc_3
    w1_rc_4 <- nn_w1_rc_4 + etann_w1_rc_4
    w1_rc_5 <- nn_w1_rc_5 + etann_w1_rc_5

    b1_rc_1 <- nn_b1_rc_1 + etann_b1_rc_1
    b1_rc_2 <- nn_b1_rc_2 + etann_b1_rc_2
    b1_rc_3 <- nn_b1_rc_3 + etann_b1_rc_3
    b1_rc_4 <- nn_b1_rc_4 + etann_b1_rc_4
    b1_rc_5 <- nn_b1_rc_5 + etann_b1_rc_5

    w2_rc_1 <- nn_w2_rc_1 + etann_w2_rc_1
    w2_rc_2 <- nn_w2_rc_2 + etann_w2_rc_2
    w2_rc_3 <- nn_w2_rc_3 + etann_w2_rc_3
    w2_rc_4 <- nn_w2_rc_4 + etann_w2_rc_4
    w2_rc_5 <- nn_w2_rc_5 + etann_w2_rc_5

    # 4. ODE system. PK first (Braem 2026 Eq 16, first two lines), which
    #    is Monolix's pkmodel(Tlag, ka, V, Cl) written out.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag

    Cc <- central / vc

    # 5. Softplus hidden layer evaluated at the CURRENT concentration.
    #    Operator precedence follows Warf_node_converted.txt exactly:
    #    the bias is OUTSIDE the -(w^2 * Cc) product but INSIDE the
    #    beta-scaled exponent. Evaluating this at Cc = 0.5 with the
    #    population weights gives 2.02, which is the left-hand edge of
    #    Braem 2026 Figure 4E; the alternative reading
    #    beta * (-(w^2) * Cc + b) gives about 42 and is excluded.
    h_rc_1 <- log(1 + exp(nn_beta_rc * (-(w1_rc_1^2) * Cc) + b1_rc_1)) / nn_beta_rc
    h_rc_2 <- log(1 + exp(nn_beta_rc * (-(w1_rc_2^2) * Cc) + b1_rc_2)) / nn_beta_rc
    h_rc_3 <- log(1 + exp(nn_beta_rc * (-(w1_rc_3^2) * Cc) + b1_rc_3)) / nn_beta_rc
    h_rc_4 <- log(1 + exp(nn_beta_rc * (-(w1_rc_4^2) * Cc) + b1_rc_4)) / nn_beta_rc
    h_rc_5 <- log(1 + exp(nn_beta_rc * (-(w1_rc_5^2) * Cc) + b1_rc_5)) / nn_beta_rc

    nnrc <- w2_rc_1 * h_rc_1 + w2_rc_2 * h_rc_2 + w2_rc_3 * h_rc_3 +
      w2_rc_4 * h_rc_4 + w2_rc_5 * h_rc_5

    # 6. Indirect response on prothrombin complex activity. The network
    #    IS the production term; Warf_node.txt also computes
    #    kin = R0 * kout but never uses it, so no kin appears here.
    d/dt(PCA) <- nnrc - kout * PCA
    PCA(0) <- rbase

    # 7. Observations
    Cc ~ add(addSd) + prop(propSd) + combined1()
    PCA ~ add(addSd_PCA) + prop(propSd_PCA) + combined1()
  })
}
