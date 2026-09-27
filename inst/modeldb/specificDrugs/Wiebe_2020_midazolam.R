Wiebe_2020_midazolam <- function() {
  description <- paste(
    "Composite parent-metabolite population PK model for midazolam and",
    "1'-OH-midazolam in healthy adults, with CYP3A drug-drug-interaction",
    "effects for reversible inhibition (ketoconazole, voriconazole),",
    "irreversible inhibition (ritonavir) and induction (efavirenz) (Wiebe",
    "2020). Both analytes have three-compartment disposition. Oral",
    "midazolam is absorbed first-order from a depot that simultaneously",
    "feeds 1'-OH-midazolam through a pre-systemic formation rate constant;",
    "all systemic midazolam clearance forms 1'-OH-midazolam (fraction",
    "metabolised fixed to 1). Body weight is a power covariate on the first",
    "metabolite intercompartmental clearance. Treatment effects are",
    "additive shifts on the typical midazolam central volume, midazolam",
    "clearance, bioavailability, pre-systemic formation rate and",
    "metabolite clearance. Inter-occasion variability on midazolam",
    "clearance, metabolite clearance and bioavailability. Amounts are nmol",
    "and concentrations nM. Parameter values are from the deposited final",
    "NONMEM control stream (Online Resource 2)."
  )
  reference <- paste(
    "Wiebe ST, Meid AD, Mikus G. Composite midazolam and 1'-OH midazolam",
    "population pharmacokinetic model for constitutive, inhibited and",
    "induced CYP3A activity. J Pharmacokinet Pharmacodyn. 2020;47(6):527-542.",
    "doi:10.1007/s10928-020-09704-1.",
    "Parameter values from Online Resource 2 (final NONMEM control stream",
    "'Adopted Composite Model Control Stream with Interaction')."
  )
  vignette <- "Wiebe_2020_midazolam"
  units <- list(
    time = "h",
    dosing = "nmol",
    concentration = "nM"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on the first 1'-OH-midazolam intercompartmental",
        "clearance only: TVQMP = THETA(13)*(WT/70)**THETA(24) in the",
        "Online Resource 2 $PK block, with THETA(24) = 0.9863. Methods: the",
        "covariate is 'normalized to [its] approximate mean value (70 kg)'.",
        "Table 2: development-set weight mean 71.7 kg, range 47-111 kg."
      ),
      source_name = "WT"
    ),
    CONMED_KETOCONAZOLE = list(
      description = "Concomitant ketoconazole (reversible CYP3A inhibitor)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no ketoconazole)",
      notes = paste(
        "One of the two drugs pooled into the paper's reversible-inhibition",
        "category (control stream INH1 = 1 when TRT3 = 2). Development data:",
        "study K345, midazolam 1 mg or 3 mg oral on days 2 and 8 of",
        "ketoconazole 400 mg once daily (Online Resource 1 Table S1).",
        "Inside model() the reversible-inhibition indicator is",
        "CONMED_KETOCONAZOLE OR CONMED_VORICONAZOLE. The source dataset",
        "never set two treatment categories at once, so do not combine this",
        "with CONMED_RTV = 1 or CONMED_CYP3A4_IND = 1 on the same record.",
        "Reversible inhibition also selects the inhibition residual-error",
        "magnitudes (control stream TRT = 2)."
      ),
      source_name = "TRT3 (= 2, reversible inhibition)"
    ),
    CONMED_VORICONAZOLE = list(
      description = "Concomitant voriconazole (reversible CYP3A inhibitor)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no voriconazole)",
      notes = paste(
        "The second drug pooled into the paper's reversible-inhibition",
        "category (control stream INH1). Development data: study K257, days",
        "1, 2, 3 and 9 of voriconazole 400 mg bid then 200 mg bid oral, and",
        "study K380, 400 mg voriconazole oral or iv with 3 microgram oral",
        "midazolam (Online Resource 1 Table S1; the 50 mg voriconazole arms",
        "were not modelled). The Discussion notes that voriconazole also",
        "shows time-dependent inhibition, which the pooled reversible",
        "category does not capture. See CONMED_KETOCONAZOLE for the pooling",
        "rule."
      ),
      source_name = "TRT3 (= 2, reversible inhibition)"
    ),
    CONMED_RTV = list(
      description = "Concomitant ritonavir (irreversible CYP3A inhibitor)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no ritonavir)",
      notes = paste(
        "The paper's irreversible-inhibition category (control stream",
        "INH2 = 1 when TRT3 = 5). Development data: study K257, days 2 and 3",
        "of ritonavir 300 mg bid oral, and study K119 condition 3, after 14",
        "days of ritonavir 300 mg bid given together with St. John's wort",
        "300 mg tid (net inhibition; Online Resource 1 Table S1).",
        "Irreversible inhibition also selects the inhibition residual-error",
        "magnitudes (control stream TRT = 2)."
      ),
      source_name = "TRT3 (= 5, irreversible inhibition)"
    ),
    CONMED_CYP3A4_IND = list(
      description = "Concomitant CYP3A inducer",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no CYP3A inducer)",
      notes = paste(
        "The paper's induction category (control stream IND = 1 when",
        "TRT3 = 3). The only development data are study K155: oral",
        "midazolam 3 mg after 14 days of efavirenz 400 mg once daily",
        "(12 subjects; Online Resource 1 Table S1). The paper applied the",
        "same effect in its external validation to midazolam 1 day after",
        "stopping St. John's wort plus ritonavir (strong induction) and to",
        "1 and 6 days after a single 400 mg efavirenz dose (weak induction",
        "or activation). There the model under-predicted midazolam",
        "concentrations, so the effect is calibrated to potent induction."
      ),
      source_name = "TRT3 (= 3, induction)"
    ),
    OCC = list(
      description = "Occasion index for inter-occasion variability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Integer 1-5. The control stream uses $ABBREVIATED REPLACE",
        "ETA(OCC_QMET)=ETA(,10 to 14), ETA(OCC_CLM)=ETA(,15 to 19) and",
        "ETA(OCC_F1)=ETA(,20 to 24), so each of the three IOV-carrying",
        "parameters has five occasion etas that share one variance",
        "($OMEGA BLOCK(1) SAME). Decomposed inside model() into five",
        "mutually exclusive indicators. For a single-occasion simulation",
        "set OCC = 1 throughout."
      ),
      source_name = "OCC"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "midazolam",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "midazolam",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "midazolam",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral2 = list(
      analyte = "midazolam",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    central_1ohm = list(
      analyte = "1'-OH-midazolam",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_1ohm = list(
      analyte = "1'-OH-midazolam",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral2_1ohm = list(
      analyte = "1'-OH-midazolam",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 99,
    n_studies = 7,
    age_range = "19-52 years",
    age_mean = "26.9 years",
    weight_range = "47-111 kg",
    weight_mean = "71.7 kg",
    sex_female_pct = 40.4,
    disease_state = "healthy volunteers",
    dose_range = paste(
      "Midazolam 3 microgram to 4 mg oral, 1 microgram to 2 mg iv, and",
      "semi-simultaneous 4 mg oral followed by 2 mg iv 6 h later; alone or",
      "with ketoconazole, voriconazole, ritonavir or efavirenz."
    ),
    regions = "Germany (University Hospital Heidelberg)",
    notes = paste(
      "Model development set: 99 healthy adults from seven studies (K119,",
      "K155, K169, K257, K380, K194, K345; Table 2), 59 male and 40",
      "female. 2371 midazolam and 2197 1'-OH-midazolam concentrations for",
      "constitutive CYP3A activity, plus 1077 and 961 from the inhibition",
      "and induction arms (6606 in total). BLQ records were omitted (29",
      "midazolam, 190 1'-OH-midazolam). External validation used 46 more",
      "subjects from three limited-sampling studies (K292, K342, K363)",
      "and further arms of K119, K169 and K194. NONMEM 7.3, ADVAN6,",
      "FOCE-I. The interaction effects were estimated with every",
      "composite-model parameter fixed, so the composite model is the",
      "special case with all treatment indicators 0."
    )
  )

  ini({
    # =================================================================
    # SOURCING NOTE. Every value below is from the deposited final
    # control stream (Online Resource 2), whose $THETA and $OMEGA hold
    # the final estimates: the composite-model values are FIXed there
    # and the interaction values are the Table 3 'Interaction' estimates
    # to more digits. Table 3 of the article prints the $OMEGA rows as
    # STANDARD DEVIATIONS (and the Vc-Vmet covariance as a CORRELATION)
    # under an 'omega^2' header: sqrt(0.142583) = 0.378 = Table 3
    # 'omega2 Vc', and 0.109819 / sqrt(0.158449 * 0.142583) = 0.731 =
    # Table 3 'omega2 Vmet * omega2 Vc'. The stream variances are used.
    # =================================================================

    # ---- Midazolam absorption and disposition ----
    lka <- fixed(log(2.30635))
    label("Midazolam absorption rate constant (1/h)")                            # Online Resource 2 THETA(6) = 2.30635 FIX ';KA'; Table 3 'ka parent' 2.31 (RSE 4.44%)
    lfdepot <- fixed(log(0.275824))
    label("Midazolam oral bioavailability (fraction)")                           # Online Resource 2 THETA(7) = 0.275824 FIX ';F'; Table 3 'F Combined bioavailability' 0.276
    lvc <- fixed(log(19.4507))
    label("Midazolam central volume (L)")                                        # Online Resource 2 THETA(1) = 19.4507 FIX ';VC'; Table 3 'Vc' 19.5 L
    lvp <- fixed(log(41.0258))
    label("Midazolam peripheral 1 volume (L)")                                   # Online Resource 2 THETA(2) = 41.0258 FIX ';VP1'; Table 3 'Vp1' 41.0 L
    lvp2 <- fixed(log(23.8228))
    label("Midazolam peripheral 2 volume (L)")                                   # Online Resource 2 THETA(3) = 23.8228 FIX ';VP2'; Table 3 'Vp2' 23.8 L
    lq <- fixed(log(8.00001))
    label("Midazolam intercompartmental clearance 1 (L/h)")                      # Online Resource 2 THETA(4) = 8.00001 FIX ';QP1'; Table 3 'Qp1' 8.00 (unit printed as h-1; used as a clearance, K23 = QP1/VC)
    lq2 <- fixed(log(46.0527))
    label("Midazolam intercompartmental clearance 2 (L/h)")                      # Online Resource 2 THETA(5) = 46.0527 FIX ';QP2'; Table 3 'Qp2' 46.1 (unit printed as h-1; used as a clearance, K24 = QP2/VC)
    lcl <- fixed(log(24.1486))
    label("Midazolam clearance, all forming 1'-OH-midazolam (L/h)")              # Online Resource 2 THETA(10) = 24.1486 FIX ';QMET'; Table 3 'Qmet MDZ clearance via metabolism' 24.1 L/h

    # ---- 1'-OH-midazolam formation and disposition ----
    lk_1ohm_form <- fixed(log(5.31099))
    label("Pre-systemic 1'-OH-midazolam formation rate constant from the depot (1/h)")  # Online Resource 2 THETA(8) = 5.31099 FIX ';KMET'; Table 3 'kmet Pre-systemic metabolism rate' 5.31
    lvc_1ohm <- fixed(log(175.643))
    label("1'-OH-midazolam central volume (L)")                                  # Online Resource 2 THETA(9) = 175.643 FIX ';VMET'; Table 3 'Vmet' 175.6 L
    lvp_1ohm <- fixed(log(684.933))
    label("1'-OH-midazolam peripheral 1 volume (L)")                             # Online Resource 2 THETA(12) = 684.933 FIX ';VMP'; Table 3 'VMP' 684.9 L
    lvp2_1ohm <- fixed(log(67.145))
    label("1'-OH-midazolam peripheral 2 volume (L)")                             # Online Resource 2 THETA(14) = 67.145 FIX ';VMP2'; Table 3 'VMP2' 67.1 L
    lq_1ohm <- fixed(log(59.5749))
    label("1'-OH-midazolam intercompartmental clearance 1 at 70 kg (L/h)")       # Online Resource 2 THETA(13) = 59.5749 FIX ';QMP'; Table 3 'QMP' 59.6
    lq2_1ohm <- fixed(log(127.374))
    label("1'-OH-midazolam intercompartmental clearance 2 (L/h)")                # Online Resource 2 THETA(15) = 127.374 FIX ';QMP2'; Table 3 'QMP2' 127.4
    lcl_1ohm <- fixed(log(196.827))
    label("1'-OH-midazolam clearance (L/h)")                                     # Online Resource 2 THETA(11) = 196.827 FIX ';CLM'; Table 3 'CLmet' 196.8 L/h
    fm <- fixed(1)
    label("Fraction of systemic midazolam clearance forming 1'-OH-midazolam (fraction)")  # Online Resource 2 $PK 'FM = 1 ; systemic fraction metabolized, fixed to 1'; Methods 'systemic fraction of midazolam metabolized fixed to 1'

    e_wt_q_1ohm <- fixed(0.9863)
    label("Power exponent of body weight on 1'-OH-midazolam Q1 (unitless)")      # Online Resource 2 THETA(24) = 0.9863 FIX ';QMP~WT'; Table 3 'QMP,WT' 0.986

    # ---- Treatment effects: additive shifts on the typical values ----
    # Reversible inhibition (ketoconazole, voriconazole; control stream INH1)
    e_inh_rev_cl <- -16.336
    label("Additive shift in midazolam clearance, reversible inhibition (L/h)")  # Online Resource 2 THETA(26) = -16.336 ';QMET~INH1'; Table 3 'Qmet,INH1' -16.3
    e_inh_rev_fdepot <- 0.378806
    label("Additive shift in bioavailability, reversible inhibition (fraction)")  # Online Resource 2 THETA(29) = 0.378806 ';F~INH1'; Table 3 'F INH1 [%]' 37.9
    e_inh_rev_cl_1ohm <- -64.1022
    label("Additive shift in 1'-OH-midazolam clearance, reversible inhibition (L/h)")  # Online Resource 2 THETA(34) = -64.1022 ';CLM~INH1'; Table 3 'CLmet,INH1' -64.1
    # Irreversible inhibition (ritonavir; control stream INH2)
    e_conmed_rtv_vc <- 51.8404
    label("Additive shift in midazolam central volume, ritonavir (L)")           # Online Resource 2 THETA(25) = 51.8404 ';VC~INH2'; Table 3 'Vc,INH2' 51.8
    e_conmed_rtv_cl <- -12.67
    label("Additive shift in midazolam clearance, ritonavir (L/h)")              # Online Resource 2 THETA(27) = -12.67 ';QMET~INH2'; Table 3 'Qmet,INH2' -12.7
    e_conmed_rtv_fdepot <- 1.34248
    label("Additive shift in bioavailability, ritonavir (fraction)")             # Online Resource 2 THETA(30) = 1.34248 ';F~INH2'; Table 3 'F INH2 [%]' 134
    e_conmed_rtv_k_1ohm_form <- -5.23662
    label("Additive shift in pre-systemic formation rate constant, ritonavir (1/h)")  # Online Resource 2 THETA(32) = -5.23662 ';KMET~INH2'; Table 3 'kmet,INH2' -5.24
    e_conmed_rtv_cl_1ohm <- 1117.36
    label("Additive shift in 1'-OH-midazolam clearance, ritonavir (L/h)")        # Online Resource 2 THETA(35) = 1117.36 ';CLM~INH2'; Table 3 'CLmet,INH2' 1117
    # Induction (efavirenz; control stream IND)
    e_conmed_cyp3a4_ind_cl <- 37.9834
    label("Additive shift in midazolam clearance, induction (L/h)")              # Online Resource 2 THETA(28) = 37.9834 ';QMET~IND'; Table 3 'Qmet,IND' 38.0
    e_conmed_cyp3a4_ind_fdepot <- -0.200032
    label("Additive shift in bioavailability, induction (fraction)")             # Online Resource 2 THETA(31) = -0.200032 ';F~IND'; Table 3 'F IND [%]' -20.0
    e_conmed_cyp3a4_ind_k_1ohm_form <- 12.4264
    label("Additive shift in pre-systemic formation rate constant, induction (1/h)")  # Online Resource 2 THETA(33) = 12.4264 ';KMET~IND'; Table 3 'kmet,IND' 12.4

    # ---- Between-subject variability (control stream $OMEGA variances) ----
    etalvc_1ohm + etalvc ~ c(0.158449, 0.109819, 0.142583)  # Online Resource 2 $OMEGA BLOCK(2) FIX: IIV_VMET 0.158449, covariance 0.109819, IIV_VC 0.142583; Table 3 prints SDs 0.398 / 0.378 and correlation 0.731
    etalcl ~ fixed(0.0087203)          # Online Resource 2 $OMEGA IIV_QMET 0.0087203 FIX; Table 3 'omega2 Qmet' 0.0934 = SD
    etalcl_1ohm ~ fixed(0.0176367)     # Online Resource 2 $OMEGA IIV_CLM 0.0176367 FIX; Table 3 'omega2 CLmet' 0.133 = SD
    etalvp ~ fixed(0.167708)           # Online Resource 2 $OMEGA IIV_VP1 0.167708 FIX; Table 3 'omega2 Vp1' 0.410 = SD
    etalq ~ fixed(0.251747)            # Online Resource 2 $OMEGA IIV_QP1 0.251747 FIX; Table 3 'omega2 Qp1' 0.502 = SD
    etalfdepot ~ fixed(0.0545658)      # Online Resource 2 $OMEGA IIV_F 0.0545658 FIX; Table 3 'omega2 F' 0.234 = SD
    etalk_1ohm_form ~ fixed(0.135776)  # Online Resource 2 $OMEGA IIV_KMET 0.135776 FIX; Table 3 'omega2 kmet' 0.369 = SD
    etalq_1ohm ~ fixed(0.235254)       # Online Resource 2 $OMEGA IIV_QMP 0.235254 FIX; Table 3 'omega2 QMP' 0.485 = SD

    # ---- Inter-occasion variability: five occasions per parameter ----
    # $OMEGA BLOCK(1) FIX then four BLOCK(1) SAME; each occasion slot is a
    # separate eta held at the shared variance.
    etaiov_cl_1 ~ fixed(0.0231946)     # Online Resource 2 $OMEGA IOV_QMET 0.0231946 FIX; Table 3 IOV 'omega2 Qmet' 0.152 = SD
    etaiov_cl_2 ~ fixed(0.0231946)     # SAME as occasion 1
    etaiov_cl_3 ~ fixed(0.0231946)     # SAME as occasion 1
    etaiov_cl_4 ~ fixed(0.0231946)     # SAME as occasion 1
    etaiov_cl_5 ~ fixed(0.0231946)     # SAME as occasion 1
    etaiov_cl_1ohm_1 ~ fixed(0.0859661)  # Online Resource 2 $OMEGA IOV_CLM 0.0859661 FIX; Table 3 IOV 'omega2 CLmet' 0.293 = SD
    etaiov_cl_1ohm_2 ~ fixed(0.0859661)  # SAME as occasion 1
    etaiov_cl_1ohm_3 ~ fixed(0.0859661)  # SAME as occasion 1
    etaiov_cl_1ohm_4 ~ fixed(0.0859661)  # SAME as occasion 1
    etaiov_cl_1ohm_5 ~ fixed(0.0859661)  # SAME as occasion 1
    etaiov_fdepot_1 ~ fixed(0.0261532)   # Online Resource 2 $OMEGA IOV_F 0.0261532 FIX; Table 3 IOV 'omega2 F' 0.162 = SD
    etaiov_fdepot_2 ~ fixed(0.0261532)   # SAME as occasion 1
    etaiov_fdepot_3 ~ fixed(0.0261532)   # SAME as occasion 1
    etaiov_fdepot_4 ~ fixed(0.0261532)   # SAME as occasion 1
    etaiov_fdepot_5 ~ fixed(0.0261532)   # SAME as occasion 1

    # ---- Residual error: two-step (early/late) proportional ----
    # The control stream codes W = SQRT(FPROP**2 * IPRED**2 + FADD**2)
    # with $SIGMA 1 FIX, so each FPROP is an SD and its sign is
    # immaterial; the two negative estimates are stored as magnitudes.
    propSd_early <- fixed(0.503162)
    label("Midazolam proportional residual SD, time <= 0.5 h, no inhibition (fraction)")  # Online Resource 2 THETA(16) = 0.503162 FIX ';FPROP1'; Table 3 'Early MDZ,prop' 0.503
    propSd_late <- fixed(0.148938)
    label("Midazolam proportional residual SD, time > 0.5 h, no inhibition (fraction)")   # Online Resource 2 THETA(20) = 0.148938 FIX ';FPROP3'; Table 3 'Late MDZ,prop' 0.149
    propSd_1ohm_early <- fixed(0.555696)
    label("1'-OH-midazolam proportional residual SD, time <= 0.5 h, no inhibition (fraction)")  # Online Resource 2 THETA(18) = 0.555696 FIX ';FPROP2'; Table 3 'Early 1'-OHMDZ,prop' 0.556
    propSd_1ohm_late <- fixed(0.215375)
    label("1'-OH-midazolam proportional residual SD, time > 0.5 h, no inhibition (fraction)")   # Online Resource 2 THETA(22) = -0.215375 FIX ';FPROP4' (enters squared); Table 3 'Late 1'-OHMDZ,prop' -0.215
    addSd_1ohm <- fixed(0.00001)
    label("1'-OH-midazolam additive residual SD, no inhibition (nM)")            # Online Resource 2 THETA(19) = THETA(23) = 0.00001 FIX ';FADD2' / ';FADD4'; the midazolam FADD1 / FADD3 and inhibition FADD5 / FADD6 are 0 FIX
    propSd_inh_early <- 0.482949
    label("Proportional residual SD, both analytes, time <= 1.5 h, inhibition (fraction)")  # Online Resource 2 THETA(36) = 0.482949 ';FPROP5'; Table 3 'Early prop,INH' 0.483
    propSd_inh_late <- 0.267324
    label("Proportional residual SD, both analytes, time > 1.5 h, inhibition (fraction)")   # Online Resource 2 THETA(38) = -0.267324 ';FPROP6' (enters squared); Table 3 'Late prop,INH' -0.267
  })

  model({
    # ---------------------------------------------------------------
    # 1. Treatment indicators (control stream INH1 / INH2 / IND).
    # ---------------------------------------------------------------
    inh_rev <- CONMED_KETOCONAZOLE + CONMED_VORICONAZOLE
    if (inh_rev > 1) inh_rev <- 1
    inh_any <- inh_rev + CONMED_RTV
    if (inh_any > 1) inh_any <- 1

    # Occasion indicators for the IOV slots
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5
    iov_cl_1ohm <- oc1 * etaiov_cl_1ohm_1 + oc2 * etaiov_cl_1ohm_2 +
      oc3 * etaiov_cl_1ohm_3 + oc4 * etaiov_cl_1ohm_4 + oc5 * etaiov_cl_1ohm_5
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 +
      oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4 + oc5 * etaiov_fdepot_5

    # ---------------------------------------------------------------
    # 2. Individual parameters. Treatment effects are ADDITIVE on the
    #    linear-scale typical value, as in the control stream, e.g.
    #    TVQMET = THETA(10)+INH1*THETA(26)+INH2*THETA(27)+IND*THETA(28).
    # ---------------------------------------------------------------
    ka <- exp(lka)
    vc <- (exp(lvc) + e_conmed_rtv_vc * CONMED_RTV) * exp(etalvc)
    vp <- exp(lvp + etalvp)
    vp2 <- exp(lvp2)
    q <- exp(lq + etalq)
    q2 <- exp(lq2)
    cl <- (exp(lcl) + e_inh_rev_cl * inh_rev + e_conmed_rtv_cl * CONMED_RTV +
      e_conmed_cyp3a4_ind_cl * CONMED_CYP3A4_IND) * exp(etalcl + iov_cl)
    fdepot <- (exp(lfdepot) + e_inh_rev_fdepot * inh_rev +
      e_conmed_rtv_fdepot * CONMED_RTV +
      e_conmed_cyp3a4_ind_fdepot * CONMED_CYP3A4_IND) *
      exp(etalfdepot + iov_fdepot)

    k_1ohm_form <- (exp(lk_1ohm_form) +
      e_conmed_rtv_k_1ohm_form * CONMED_RTV +
      e_conmed_cyp3a4_ind_k_1ohm_form * CONMED_CYP3A4_IND) *
      exp(etalk_1ohm_form)
    vc_1ohm <- exp(lvc_1ohm + etalvc_1ohm)
    vp_1ohm <- exp(lvp_1ohm)
    vp2_1ohm <- exp(lvp2_1ohm)
    q_1ohm <- exp(lq_1ohm + etalq_1ohm) * (WT / 70)^e_wt_q_1ohm
    q2_1ohm <- exp(lq2_1ohm)
    cl_1ohm <- (exp(lcl_1ohm) + e_inh_rev_cl_1ohm * inh_rev +
      e_conmed_rtv_cl_1ohm * CONMED_RTV) * exp(etalcl_1ohm + iov_cl_1ohm)

    # ---------------------------------------------------------------
    # 3. Micro-constants (control stream K23 ... KM0).
    # ---------------------------------------------------------------
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2
    kmet <- fm * cl / vc
    k12_1ohm <- q_1ohm / vc_1ohm
    k21_1ohm <- q_1ohm / vp_1ohm
    k13_1ohm <- q2_1ohm / vc_1ohm
    k31_1ohm <- q2_1ohm / vp2_1ohm
    kel_1ohm <- cl_1ohm / vc_1ohm

    # ---------------------------------------------------------------
    # 4. ODEs (control stream $DES, COMP 1-7). The control stream sets
    #    KA = THETA(6) on midazolam records and KA = KMET on metabolite
    #    records, and its depot equation is DADT(1) = -KA*A(1) while the
    #    metabolite receives KMET*A(1). NONMEM advances each interval
    #    with the PK values of the record that closes it; with ka
    #    estimable (Table 3 RSE 4.44%) the intervals must close on
    #    midazolam records, so the depot empties at ka into the
    #    midazolam central compartment while k_1ohm_form * depot is
    #    added to 1'-OH-midazolam in parallel, without draining the
    #    depot. This as-run form reproduces the 1'-OH-midazolam
    #    concentrations of the paper's Figure 3 (see the vignette).
    # ---------------------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2 - kmet * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    d/dt(central_1ohm) <- k_1ohm_form * depot + kmet * central -
      k13_1ohm * central_1ohm + k31_1ohm * peripheral2_1ohm -
      k12_1ohm * central_1ohm + k21_1ohm * peripheral1_1ohm -
      kel_1ohm * central_1ohm
    d/dt(peripheral1_1ohm) <- k12_1ohm * central_1ohm - k21_1ohm * peripheral1_1ohm
    d/dt(peripheral2_1ohm) <- k13_1ohm * central_1ohm - k31_1ohm * peripheral2_1ohm

    f(depot) <- fdepot

    # ---------------------------------------------------------------
    # 5. Observations (nmol / L = nM) and the two-step residual error.
    #    The early/late split is on time since the occasion's first
    #    midazolam dose (control stream TIME), at 0.5 h without and
    #    1.5 h with inhibition; the inhibition magnitudes are shared by
    #    both analytes.
    # ---------------------------------------------------------------
    Cc <- central / vc
    Cc_1ohm <- central_1ohm / vc_1ohm

    early <- (t <= 0.5)
    early_inh <- (t <= 1.5)
    propSd_mdz <- (1 - inh_any) * (early * propSd_early + (1 - early) * propSd_late) +
      inh_any * (early_inh * propSd_inh_early + (1 - early_inh) * propSd_inh_late)
    propSd_met <- (1 - inh_any) * (early * propSd_1ohm_early + (1 - early) * propSd_1ohm_late) +
      inh_any * (early_inh * propSd_inh_early + (1 - early_inh) * propSd_inh_late)

    Cc ~ prop(propSd_mdz)
    Cc_1ohm ~ add(addSd_1ohm) + prop(propSd_met)
  })
}
