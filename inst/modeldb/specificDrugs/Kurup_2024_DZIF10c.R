Kurup_2024_DZIF10c <- function() {
  description <- paste(
    "Pooled human + cynomolgus macaque semi-mechanistic population PK model for DZIF-10c",
    "(BI 767551), a SARS-CoV-2 neutralizing human IgG1 kappa monoclonal antibody developed",
    "for both intravenous infusion and inhaled (nebulized) administration (Kurup 2024).",
    "Five compartments in three regions: a lung epithelial lining fluid (ELF) space carried as",
    "three parallel pools that share one volume and one clearance and differ only by route of",
    "entry (elf_lrt = lower-respiratory-tract deposition site for an inhaled dose,",
    "elf_trachea = tracheal deposition site for an intratracheal dose, elf = the rest of the",
    "lung, the only pool that exchanges bidirectionally with the systemic circulation), plus a",
    "central and a peripheral compartment with linear clearance. The observed ELF concentration",
    "is the summed amount of all three lung pools divided by the shared ELF volume. Drug leaves",
    "the lung by two routes: distribution into the systemic circulation (Q1) and loss from the",
    "lung itself (CL_L, representing e.g. mucociliary clearance). Human and macaque data were",
    "bridged by body-weight allometry alone, with fixed exponents of 0.85 on every clearance",
    "and 1 on every volume, both referenced to 70 kg for BOTH species; the only",
    "species-dependent parameter is the fraction of an inhaled dose deposited in the lung",
    "(29.6% in humans vs 1.00% in macaques), estimated on the logit scale. The intratracheal",
    "dose is assumed to be fully deposited (F = 1). Fitted to 640 observations from 76 subjects",
    "across 3 preclinical macaque studies (serum/plasma plus urea-corrected bronchoalveolar",
    "lavage ELF) and 2 phase I/IIa clinical trials (serum only)."
  )
  reference <- paste(
    "Kurup S, Velez de Mendizabal N, Becker S, Bolella E, De Sousa D, Fatkenheuer G,",
    "Gruell H, Klein F, Malin JJ, Schmid U, Korell J.",
    "Semi-mechanistic population pharmacokinetic modeling of DZIF-10c, a neutralizing",
    "antibody against SARS-Cov-2: predicting systemic and lung exposure following inhaled",
    "and intravenous administration.",
    "J Pharmacokinet Pharmacodyn. 2024;52(1):3. doi:10.1007/s10928-024-09947-2"
  )
  vignette <- "Kurup_2024_DZIF10c"
  # Declared explicitly because the registry's automatic detection only
  # recognises `depot` and `central`, and would otherwise record this model as
  # central-dosed only -- which would hide the two routes the paper exists to
  # compare. An inhaled (nebulized) dose targets cmt = "elf_lrt", an
  # intratracheal dose cmt = "elf_trachea", and an IV infusion cmt = "central".
  dosing <- c("elf_lrt", "elf_trachea", "central")
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Every entry was checked against Kurup 2024 Fig. 1 and
  # supplement equations (1)-(7).
  compartmentData <- list(
    elf_lrt     = list(analyte = "DZIF-10c", units = "mg", specimen = "epithelial lining fluid", verified = TRUE),
    elf_trachea = list(analyte = "DZIF-10c", units = "mg", specimen = "epithelial lining fluid", verified = TRUE),
    elf         = list(analyte = "DZIF-10c", units = "mg", specimen = "epithelial lining fluid", verified = TRUE),
    # Systemic observations were serum in both clinical trials and in preclinical
    # study 1, and plasma in preclinical studies 2 and 3; the paper pooled them
    # after finding no significant Vc difference between the two matrices
    # (Discussion, "Systemic PK observations collected in the clinical trials...").
    central     = list(analyte = "DZIF-10c", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "DZIF-10c", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric size descriptor on every clearance and every volume, referenced to 70 kg",
        "(supplement eq. 11, P_ri = theta_r * (WT_i / 70 kg)^beta_s; control stream applies",
        "(WT/70)**ACL and (WT/70)**AV unconditionally). NOTE the deliberate departure from the",
        "SPECIES_MACAQUE register note that each species should carry its own allometric",
        "reference weight: this paper normalises BOTH species to a single 70 kg reference, and",
        "that single-reference normalisation IS the cross-species bridge the analysis relies on.",
        "Baseline (time-fixed) in this analysis. Macaque weights are far below 70 kg, so the",
        "macaque parameter values are extrapolations of the human-referenced typical values."
      ),
      source_name        = "WT"
    ),
    SPECIES_MACAQUE = list(
      description        = "Cynomolgus macaque indicator, 1 = cynomolgus macaque, 0 = human",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (human)",
      notes              = paste(
        "Selects which of the two logit-scale inhaled-deposition parameters applies.",
        "Every other PK parameter is species independent (Results: 'The rest of the PK",
        "parameters in this selected model were species independent'), so this indicator",
        "acts on the inhaled bioavailability only. Source column HUMAN in the control stream",
        "is the complement of this canonical: HUMAN = 1 - SPECIES_MACAQUE. The $PK block reads",
        "'PHI = PHIM; IF (HUMAN.EQ.1) PHI = PHIH', a stratum switch rather than a",
        "reference-plus-offset contrast, so no effect-coefficient sign is involved."
      ),
      source_name        = "HUMAN"
    )
  )

  population <- list(
    species        = "human + cynomolgus macaque (Macaca fascicularis)",
    n_subjects     = 76,
    n_studies      = 5,
    age_range      = "adults (clinical trials enrolled participants aged 18 years and older)",
    weight_range   = "not reported per subject; typical values referenced to 70 kg",
    sex_female_pct = NA_real_,
    disease_state  = paste(
      "SARS-CoV-2 infected and uninfected humans (healthy volunteers and COVID-19 patients),",
      "plus healthy and SARS-CoV-2 infected cynomolgus macaques."
    ),
    dose_range     = paste(
      "Human 2.5-80 mg/kg IV (~60 min infusion) and 50-250 mg inhaled (15-20 min nebulization);",
      "macaque 50 mg/kg IV bolus, 3.6 mg/kg intratracheal, and 500 or 1000 mg inhaled."
    ),
    notes          = paste(
      "Analysis dataset: 640 observations from 76 subjects, 104 (16%) below the LLOQ and",
      "excluded by IGNORE(BLQ>0). Preclinical: 14 macaques, 166 concentrations (150",
      "serum/plasma, 16 ELF). Human: 62 subjects, 474 serum concentrations (44 subjects /",
      "361 observations from the IV trial NCT04631666, 18 subjects / 113 observations from the",
      "inhaled trial NCT04631705). Bronchoalveolar-lavage concentrations from preclinical",
      "study 3 could not be urea-corrected to ELF and were excluded from the final dataset.",
      "Study-level detail in supplement Tables S1 and S2. LLOQ 20 ng/mL (macaque serum and",
      "BAL fluid) and 300 ng/mL (human serum). No demographic covariates were evaluated",
      "because of the small study sizes and the pooled cross-species design."
    )
  )

  ini({
    # Structural parameters. All typical values are the Bayesian posterior
    # medians of Kurup 2024 Table 1, back-transformed from the log domain by
    # the authors (Table 1 footnote). The supplement's $THETA block is
    # explicitly labelled ";INITIAL ESTIMATES" and is NOT used here.
    # Reference weight is 70 kg for both species.

    # Inhaled deposited fraction, estimated on the logit scale so that it is
    # bounded in (0, 1): F_INH = exp(phi) / (1 + exp(phi)), supplement
    # eqs. (8)-(9). Human is the reference species (SPECIES_MACAQUE = 0).
    logitfdepot         <- log(0.296 / (1 - 0.296))
    label("Logit of the fraction of an inhaled dose deposited in the lung, human (logit units)")   # Table 1 'F1h (%) ... 29.6 (21.9-40.0)'; $THETA 2
    logitfdepot_macaque <- log(0.0100 / (1 - 0.0100))
    label("Logit of the fraction of an inhaled dose deposited in the lung, cynomolgus macaque (logit units)")  # Table 1 'F1m (%) ... 1.00 (0.737-1.38)'; $THETA 1

    # Lung (ELF) disposition
    lq_elf  <- log(3.68e-5)
    label("Intercompartmental clearance between the lung ELF and the central compartment at WT = 70 kg (Q1, L/h)")  # Table 1 'Q1 (L/h) ... 3.68e-5 (2.58e-5-5.20e-5)'; $THETA 3
    lv_elf  <- log(0.0364)
    label("Volume of the lung ELF compartment at WT = 70 kg, shared by all three ELF pools (VL, L)")  # Table 1 'VL (L) ... 0.0364 (0.0255-0.0522)'; $THETA 4
    lcl_elf <- log(0.000412)
    label("Clearance out of the lung ELF, e.g. mucociliary clearance, at WT = 70 kg (CL_L, L/h)")  # Table 1 'CLL (L/h) ... 0.000412 (0.000291-0.000578)'; $THETA 5

    # Systemic disposition
    lcl <- log(0.0122)
    label("Clearance from the central compartment at WT = 70 kg (CL, L/h)")  # Table 1 'CL (L/h) ... 0.0122 (0.0116-0.0128)'; $THETA 6
    lvc <- log(3.08)
    label("Central compartment volume of distribution at WT = 70 kg (Vc, L)")  # Table 1 'Vc (L) ... 3.08 (2.84-3.33)'; $THETA 7
    lq  <- log(0.0359)
    label("Intercompartmental clearance between the central and peripheral compartments at WT = 70 kg (Q2, L/h)")  # Table 1 'Q2 (L/h) ... 0.0359 (0.0293-0.0435)'; $THETA 8
    lvp <- log(4.36)
    label("Peripheral compartment volume of distribution at WT = 70 kg (Vp, L)")  # Table 1 'VP (L) ... 4.36 (4.06-4.67)'; $THETA 9

    # Allometric exponents. Both were held constant rather than estimated
    # ("scaled using fixed exponents of 0.85 for clearance rates and 1 for the
    # different volumes of distribution"); the control stream hardcodes them as
    # ACL = 0.85 and AV = 1 outside the $THETA block. e_wt_cl_q applies to all
    # FOUR clearances (CL, Q1, Q2, CL_L) and e_wt_vc_vp to all THREE volumes
    # (Vc, Vp, VL), not only to the two parameters each name lists.
    e_wt_cl_q  <- fixed(0.85)
    label("Allometric (WT/70) exponent shared across every clearance: CL, Q1, Q2 and CL_L (unitless)")  # Methods 'Covariate analyses'; supplement eq. 11; $PK 'ACL = 0.85'
    e_wt_vc_vp <- fixed(1)
    label("Allometric (WT/70) exponent shared across every volume: Vc, Vp and VL (unitless)")  # Methods 'Covariate analyses'; supplement eq. 11; $PK 'AV = 1'

    # Inter-individual variability. The final model carries an eta on every
    # structural parameter, but only the two on Vc and Q2 were estimated; the
    # other seven were held at a variance of 0.0225 because their estimates
    # showed high uncertainty ("IIV was initially estimated for all PK
    # parameters, with estimates showing high uncertainty fixed to 0.0225, as
    # recommended for SAEM and BAYES estimation methods"), which is exactly the
    # seven '0.0225 FIX' rows of the supplement $OMEGA block.
    #
    # The two estimated IIVs are reported in Table 1 as percentages; they are
    # converted here with omega^2 = log(CV^2 + 1), the log-normal convention in
    # the skill's verification checklist. See vignette Errata for the
    # alternative reading (omega = CV/100), which the paper does not state.
    etalogitfdepot         ~ fixed(0.0225)  # $OMEGA 2 'F1H'   (additive on the logit scale, supplement eq. 13)
    etalogitfdepot_macaque ~ fixed(0.0225)  # $OMEGA 1 'F1M'   (additive on the logit scale, supplement eq. 13)
    etalq_elf              ~ fixed(0.0225)  # $OMEGA 3 'Q1'
    etalv_elf              ~ fixed(0.0225)  # $OMEGA 4 'VL'
    etalcl_elf             ~ fixed(0.0225)  # $OMEGA 5 'CLL'
    etalcl                 ~ fixed(0.0225)  # $OMEGA 6 'CL'
    etalvc                 ~ 0.0688626      # $OMEGA 7 estimated; Table 1 'IIV in VC (%) 26.7 (20.7-35.2)', shrinkage 28.3%; log(1 + 0.267^2)
    etalq                  ~ 0.3717297      # $OMEGA 8 estimated; Table 1 'IIV in Q2 (%) 67.1 (48.5-95.1)', shrinkage 35.0%; log(1 + 0.671^2)
    etalvp                 ~ fixed(0.0225)  # $OMEGA 9 'VP'

    # Residual error. The supplement's eq. (14) describes a combined additive
    # plus proportional model, but that is what was "initially modelled": the
    # FINAL control stream's $ERROR block is Y = IPRED * (1 + W * EPS(1)) with
    # $SIGMA 1 FIX, i.e. purely proportional, with a separate W for each
    # observed matrix. Encoded proportional-only to match the final model.
    propSd_Celf <- 0.547
    label("Proportional residual error for lung ELF concentrations (fraction)")  # Table 1 'Proportional error ELF compartment 0.547 (0.368-0.875)'; Results '54.7%'; $ERROR 'W = EXP(THETA(10))'
    propSd      <- 0.145
    label("Proportional residual error for serum/plasma concentrations (fraction)")  # Table 1 'Proportional error central compartment 0.145 (0.134-0.158)'; $ERROR 'W = EXP(THETA(11))'
  })

  model({
    # Species-specific inhaled deposition, on the logit scale. The control
    # stream selects one of two logit-scale fixed effects, each with its own
    # eta, rather than applying an offset to a common reference:
    #   PHI = PHIM ; IF (HUMAN.EQ.1) PHI = PHIH ; F1 = EXP(PHI)/(1+EXP(PHI))
    # Each fixed effect gets its own simple line so that both stay in a
    # mu-referenced position; combining them inline makes rxode2 fall back to
    # non-mu-referenced parsing for the second one.
    phih <- logitfdepot + etalogitfdepot
    phim <- logitfdepot_macaque + etalogitfdepot_macaque
    logitfinh <- phih * (1 - SPECIES_MACAQUE) + phim * SPECIES_MACAQUE
    finh <- exp(logitfinh) / (1 + exp(logitfinh))

    # Individual parameters. Every clearance carries the 0.85 exponent and
    # every volume the 1.0 exponent, both on WT/70, for humans and macaques
    # alike (supplement eq. 11).
    q_elf  <- exp(lq_elf + etalq_elf) * (WT / 70)^e_wt_cl_q
    v_elf  <- exp(lv_elf + etalv_elf) * (WT / 70)^e_wt_vc_vp
    cl_elf <- exp(lcl_elf + etalcl_elf) * (WT / 70)^e_wt_cl_q
    cl     <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q
    vc     <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q      <- exp(lq + etalq) * (WT / 70)^e_wt_cl_q
    vp     <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp

    # Micro-constants. Note that the lung-to-central transfer uses the ELF
    # volume while the central-to-lung transfer uses the central volume, so the
    # exchange is NOT symmetric in rate: Q1 is a clearance applied to whichever
    # compartment the drug is leaving.
    kel     <- cl / vc      # central elimination
    kelfc   <- q_elf / v_elf   # lung ELF -> central
    kcelf   <- q_elf / vc      # central -> lung ELF
    klosself <- cl_elf / v_elf # loss out of the lung ELF (e.g. mucociliary)
    k12     <- q / vc
    k21     <- q / vp

    # Supplement equations (1)-(5) / $DES, term for term. There is no transfer
    # from either deposition pool into `elf`: all three pools drain directly to
    # `central` (via Q1) and to lung loss (via CL_L), and only `elf` receives
    # back-transfer from `central`.
    d/dt(elf_lrt)     <- -klosself * elf_lrt - kelfc * elf_lrt
    d/dt(elf_trachea) <- -klosself * elf_trachea - kelfc * elf_trachea
    d/dt(elf)         <- kcelf * central - kelfc * elf
    d/dt(central)     <- kelfc * (elf_lrt + elf_trachea + elf) - kcelf * central +
      k21 * peripheral1 - k12 * central - kel * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Only a fraction of an inhaled dose reaches the lung; the whole
    # intratracheal dose was assumed to be deposited there ($PK 'F2 = 1',
    # supplement eq. 10). The IV dose likewise has F = 1 ($PK 'F4 = 1'), which
    # is rxode2's default for `central` and so is not restated.
    f(elf_lrt)     <- finh
    f(elf_trachea) <- 1

    # Observed lung ELF concentration is the SUM of all three lung pools over
    # the single shared ELF volume (supplement eq. 6); systemic concentration
    # is the central amount over Vc (supplement eq. 7). Amounts in mg over
    # volumes in L give mg/L = ug/mL.
    Celf <- (elf_lrt + elf_trachea + elf) / v_elf
    Cc   <- central / vc

    Celf ~ prop(propSd_Celf)
    Cc ~ prop(propSd)
  })
}
