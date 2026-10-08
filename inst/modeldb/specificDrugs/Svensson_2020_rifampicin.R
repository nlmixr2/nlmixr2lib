Svensson_2020_rifampicin <- function() {
  description <- paste(
    "Semi-mechanistic two-compartment population PK model for rifampicin in",
    "plasma and lumbar cerebrospinal fluid (CSF) in adults with tuberculous",
    "meningitis given oral rifampicin 450-1350 mg or a 600 mg intravenous",
    "infusion once daily, pooled from three Indonesian phase 2 trials",
    "(Svensson 2020). Oral absorption is a Savic analytical transit chain",
    "(4.23 transit compartments, mean transit time 0.673 h) feeding a",
    "first-order absorption compartment (ka 1.41 1/h) that empties into a",
    "liver compartment, so first-pass extraction is structural; prehepatic",
    "oral bioavailability is 0.776 (logit scale) and intravenous doses enter",
    "the central compartment as a fixed 1.5 h infusion. Elimination is a",
    "well-stirred liver with saturable intrinsic clearance, CLint = CLint,max",
    "* Km / (CH + Km) and EH = CLint * fu / (CLint * fu + QH), with liver",
    "volume 1 L, hepatic plasma flow 50 L/h and fraction unbound 0.2 fixed.",
    "Autoinduction is a step change: CLint,max is 48% higher and the central",
    "and peripheral volumes 19.3% lower on the second PK occasion (from day 4",
    "on study) than on the first. All disposition parameters are",
    "allometrically scaled on fat-free mass (reference 44.55 kg) with fixed",
    "0.75 / 1 exponents. CSF is an effect compartment holding a",
    "concentration that equilibrates with plasma at a 2.06 h half-life",
    "toward a partition coefficient of 0.0545, which rises with log10 CSF",
    "total protein. Random effects are between-subject variability on V, ka,",
    "Q and the CSF partition coefficient and two-occasion between-occasion",
    "variability on hepatic clearance, mean transit time and logit",
    "bioavailability. Residual error is combined additive plus proportional",
    "in plasma and additive in CSF. Individual day-2 plasma AUC0-24 from",
    "this model drives the companion mortality model",
    "Svensson_2020_rifampicin_survival."
  )
  reference <- paste(
    "Svensson EM, Dian S, te Brake L, Ganiem AR, Yunivita V, van Laarhoven A,",
    "van Crevel R, Ruslami R, Aarnoutse RE (2020).",
    "Model-Based Meta-analysis of Rifampicin Exposure and Mortality in",
    "Indonesian Tuberculous Meningitis Trials.",
    "Clin Infect Dis 71(8):1817-1823. doi:10.1093/cid/ciz1071.",
    "Parameter estimates from Online Data Supplement Table E1 and the",
    "commented NONMEM control stream reproduced in the same supplement",
    "('NONMEM code pharmacokinetic model'), which also carries the model",
    "equations. The saturable-hepatic-extraction structure follows",
    "Chirehwa et al. (2016) Antimicrob Agents Chemother 60(1):487-494",
    "doi:10.1128/AAC.01830-15; the transit absorption follows Savic et al.",
    "(2007) J Pharmacokinet Pharmacodyn 34(5):711-726.",
    sep = " "
  )
  vignette <- "Svensson_2020_rifampicin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The `csf` state holds a CONCENTRATION, not an amount: the source $DES
  # integrates DADT(4) = LOG(2)/HLCSF*(PC*A(3)/V - A(4)) and $ERROR reads
  # IPRED2 = A(4) with no volume division ("The relation was coded without
  # actual mass transfer in the model", supplement, Pharmacokinetic model).
  # Same convention as Abdelgawad_2025_rifampicin.R.
  compartmentData <- list(
    depot = list(
      analyte = "rifampicin",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    liver = list(
      analyte = "rifampicin",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    central = list(
      analyte = "rifampicin",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "rifampicin",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    csf = list(
      analyte = "rifampicin",
      units = "mg/L",
      specimen = "CSF",
      verified = TRUE
    )
  )

  covariateData <- list(
    FFM = list(
      description = "Fat-free mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric size descriptor for every disposition parameter with",
        "exponents fixed a priori to 0.75 for clearances (CLint,max, Q, QH)",
        "and 1 for volumes (V, Vp, VH) (supplement, Pharmacokinetic model:",
        "'Allometric scaling with fat-free mass as body size descriptor and",
        "coefficients fixed to the theoretical values (0.75 for clearance and",
        "1 for volume) was included a priori'). The control stream normalises",
        "to 44.55 kg (ALLMV = (NFMV/44.55), ALLMCL = (NFMV/44.55)**0.75) for",
        "both the estimated and the fixed hepatic parameters, and imputes",
        "44.55 kg when FFM is missing (IF(FFM.EQ.-99) NFMV = 44.55). The",
        "paper does not state the FFM formula; body weight median 46 kg",
        "(range 34-78) in the pooled cohort (Table 1)."
      ),
      source_name = "FFM"
    ),
    OCC = list(
      description = "PK sampling occasion: 1 = first occasion (day 2 +/- 1 of study treatment), 2 = second occasion (day 12 +/- 4, i.e. from day 4 on study onward)",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Source column PERIOD (values 1 and 2). It plays two roles in the",
        "control stream. (1) It multiplexes the two-occasion between-occasion",
        "etas on hepatic clearance, mean transit time and logit",
        "bioavailability (OCC1 = PERIOD.EQ.1, OCC2 = PERIOD.EQ.2, each backed",
        "by $OMEGA BLOCK(1) followed by SAME). (2) It is the autoinduction",
        "switch: IND = 1 when PERIOD.EQ.2 raises CLint,max by 48%, and the",
        "factor (1 - (PERIOD - 1) * THETA(9)) lowers V and Vp by 19.3%. The",
        "paper describes the induction as 'a factor change in intrinsic",
        "clearance from day 4 on study and onwards' and explains that a",
        "gradual change could not be estimated because patients could have",
        "received up to 3 days of rifampicin before enrolment. For",
        "simulation, set OCC = 1 for doses and observations before day 4 of",
        "treatment and OCC = 2 from day 4 onward. Only values 1 and 2 are",
        "meaningful."
      ),
      source_name = "PERIOD"
    ),
    CSF_TPRO = list(
      description = "Cerebrospinal-fluid total protein concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source column CSFPROT, in mg/dL; converted inside model() as",
        "CSF_TPRO * 100. Enters the CSF partition coefficient as a linear",
        "function of the relative deviation of log10(protein) from",
        "log10(165 mg/dL) (control stream: PC = THETA*(1 + (LCSFPROT -",
        "LOG10(165))/LOG10(165)*THETA)); 165 mg/dL is the cohort median",
        "(Table 1; Table E1 footnote b, which misprints the unit as",
        "'cells/mL'). The control stream imputes 165 mg/dL when the value is",
        "missing. Cohort range 9-3869 mg/dL (0.09-38.7 g/L)."
      ),
      source_name = "CSFPROT"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Not a model covariate; fat-free mass was chosen a priori as the body-size descriptor. Cohort median 46 kg, range 34-78 (Table 1)."
    ),
    HIV_POS = list(
      description = "HIV infection (1 = positive)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (HIV-negative)",
      notes = "Present in the control stream $INPUT (HIV) but not referenced in $PK; 18 of 148 patients (12%) were HIV-infected (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 133L,
    n_studies = 3L,
    n_observations = "1150 rifampicin concentrations (980 plasma, 170 CSF)",
    age_range = "16-81 years (median 30; whole 148-patient cohort)",
    age_median = "30 years",
    weight_range = "34-78 kg (median 46; whole 148-patient cohort)",
    weight_median = "46 kg",
    sex_female_pct = 45,
    race_ethnicity = "Indonesian (all three trials were conducted in Bandung, Indonesia)",
    disease_state = paste(
      "Adults with definite (56%), probable or possible tuberculous",
      "meningitis; 95% had British Medical Research Council grade 2 or 3",
      "disease; baseline Glasgow Coma Scale median 13 (3-15); CSF total",
      "protein median 165 mg/dL (9-3869); 12% HIV-infected. All patients",
      "received other first-line TB drugs and adjunctive dexamethasone."
    ),
    dose_range = paste(
      "Rifampicin once daily: oral 450 mg (about 10 mg/kg, control arms),",
      "750 mg (17 mg/kg), 900 mg (20 mg/kg) or 1350 mg (30 mg/kg), or 600 mg",
      "(13 mg/kg) as a 1.5 h intravenous infusion, for the first 14 or 30",
      "days; standard treatment thereafter."
    ),
    regions = "Indonesia (Bandung)",
    notes = paste(
      "Individual-patient pooled analysis of Ruslami 2013 (Lancet Infect",
      "Dis), Yunivita 2016 (Int J Antimicrob Agents) and Dian 2018",
      "(Antimicrob Agents Chemother). 133 of the 148 patients had PK data.",
      "Rich plasma sampling (6 time points) on day 2 +/- 1 and, in two of",
      "the studies, again on day 12 +/- 4; single CSF samples 3-9 h after",
      "dose. Samples below the limit of quantification were excluded from",
      "estimation. Demographics are Table 1 (whole 148-patient cohort)."
    )
  )

  ini({
    # =======================================================================
    # Point estimates are the $THETA / $OMEGA / $SIGMA final values of the
    # commented NONMEM control stream in the Online Data Supplement
    # ('NONMEM code pharmacokinetic model'); every one reproduces the
    # rounded value of supplement Table E1. Time h, dose mg, conc. mg/L.
    #
    # THETA numbering: the printed $THETA block lists 14 values, but $PK
    # reads the CSF partition coefficient as THETA(14) and the protein
    # effect as THETA(15). The two values are matched here by their $THETA
    # comments and Table E1 labels (0.0545 'Penetration coefficient CSF',
    # 0.631 'Effect of CSFPROT on PC'), not by index.
    # =======================================================================

    lclint_max <- log(40.468)
    label("Maximal intrinsic hepatic clearance CLint,max on the first PK occasion at FFM 44.55 kg (L/h)")  # Table E1 intrinsic clearance 40.5 L/h (RSE 8.5%); $THETA (0,40.468) ; 1 CLint L/h
    lkm <- log(16.1378)
    label("Michaelis-Menten constant Km for intrinsic clearance, on liver concentration (mg/L)")           # Table E1 Km 16.1 mg/L (21%); $THETA (0,16.1378,100) ; 2 KM
    lvc <- log(7.37836)
    label("Central volume of distribution V on the first PK occasion at FFM 44.55 kg (L)")                  # Table E1 volume of distribution - central 7.38 L (29%); $THETA (0,7.37836) ; 3 V2 L
    lka <- log(1.40998)
    label("First-order absorption rate constant ka from the absorption compartment into the liver (1/h)")  # Table E1 rate of absorption 1.41 1/h (17%); $THETA (0,1.40998) ; 4 KA /h
    lmtt <- log(0.672668)
    label("Mean transit time MTT through the absorption transit chain (h)")                                # Table E1 mean transit time 0.67 h (12%); $THETA (0,0.672668) ; 5 MTT h
    lntr <- log(4.23047)
    label("Number of absorption transit compartments (unitless)")                                          # Table E1 number of transit compartments 4.23 (14%); $THETA (0,4.23047,100) ; 6 NN
    logitfdepot <- logit(0.77584)
    label("Prehepatic oral bioavailability (logit scale)")                                                 # Table E1 bioavailability 77.6% (4.4%); $THETA (0,0.77584,1) ; 7 Bioavailability; $PK PHI = LOG(TVBIO/(1-TVBIO)) + IOVF
    e_occ2_clint_max <- 0.480062
    label("Fractional increase in CLint,max on the second PK occasion (autoinduction) (fraction)")         # Table E1 induction 47.9% (22%), footnote a 'factor change after day 4'; $THETA (0,0.480062) ; 8 Induction
    # One estimated theta shared by V and Vp (the control stream applies
    # (1-(PERIOD-1)*THETA(9)) to both), so it is a single parameter named
    # for the central volume and reused on Vp in model().
    e_occ2_vc <- 0.19257
    label("Fractional decrease in central (and peripheral) volume on the second PK occasion (fraction)")    # Table E1 difference on volumes after day 4 -19.3% (19%); $THETA (0,0.19257,1) ; 9 Difference in V with PERIOD
    lq <- log(93.1391)
    label("Intercompartmental clearance Q at FFM 44.55 kg (L/h)")                                          # Table E1 intercompartmental clearance 93.1 L/h (25%); $THETA (0,93.1391) ; 10 Q
    lvp <- log(26.2304)
    label("Peripheral volume of distribution Vp on the first PK occasion at FFM 44.55 kg (L)")              # Table E1 volume of distribution - peripheral 26.2 L (5.0%); $THETA (0,26.2304) ; 11 Volume peripheral

    # CSF equilibration is parameterised in the source as a HALF-LIFE
    # (DADT(4) = LOG(2)/HLCSF*(...)); converted to the canonical rate
    # constant ke0 = log(2) / 2.06319 h = 0.33596 1/h.
    lke0 <- log(log(2) / 2.06319)
    label("Plasma-to-CSF equilibration rate constant ke0 (1/h); equivalent to a half-life of 2.06 h")      # Table E1 half-life distribution 2.07 h (20%); $THETA (0,2.06319) ; 12 HL equilibrium CSF
    lppc <- log(0.0544996)
    label("CSF partition coefficient at the median CSF protein of 165 mg/dL (fraction)")                   # Table E1 penetration coefficient 5.46% (9.2%); $THETA (0,0.0544996,1) ; 13 Penetration coefficent CSF (read as THETA(14) in $PK)
    e_csf_tpro_ppc <- 0.630943
    label("Slope of the CSF partition coefficient on the relative deviation of log10 CSF protein from log10(165 mg/dL) (unitless)")  # Table E1 effect of protein 63.1% (46%); $THETA (-1,0.630943,10) ; 14 Effect of CSFPROT on PC linear (read as THETA(15) in $PK)

    # ----- Fixed hepatic physiology and dosing constants (control stream) --
    lqh <- fixed(log(50))
    label("Hepatic plasma flow QH at FFM 44.55 kg (L/h)")                                                  # supplement 'Hepatic plasma flow and liver volume were fixed to 50 L/h and 1 L'; $PK QH = 50*ALLMCL
    lvh <- fixed(log(1))
    label("Liver volume VH at FFM 44.55 kg (L)")                                                           # supplement 'liver volume ... fixed to ... 1 L'; $PK VH = 1*ALLMV
    fub <- fixed(0.2)
    label("Fraction of rifampicin unbound in plasma, fu (unitless)")                                       # supplement 'Rifampin protein binding was assumed to be 20%'; $PK FU = 0.2 ; proportion free unbound
    ldur <- fixed(log(1.5))
    label("Duration of the intravenous infusion into the central compartment (h)")                         # $PK D3 = 1.5 ; Duration of infusion; main text '600-mg ... intravenous infusion (1.5 hours)'

    # Allometric exponents, hardcoded in the control stream as ALLMCL =
    # (NFMV/44.55)**0.75 (on CLint,max, Q and QH) and ALLMV = (NFMV/44.55)
    # (on V, Vp and VH).
    e_ffm_clint_max <- fixed(0.75)
    label("Allometric exponent of fat-free mass on CLint,max (unitless)")  # $PK TVCLINTB = THETA(1)*ALLMCL*...
    e_ffm_q <- fixed(0.75)
    label("Allometric exponent of fat-free mass on Q (unitless)")          # $PK TVQ = THETA(10)*ALLMCL
    e_ffm_qh <- fixed(0.75)
    label("Allometric exponent of fat-free mass on QH (unitless)")         # $PK QH = 50*ALLMCL
    e_ffm_vc <- fixed(1)
    label("Allometric exponent of fat-free mass on V (unitless)")          # $PK TVV = THETA(3)*ALLMV
    e_ffm_vp <- fixed(1)
    label("Allometric exponent of fat-free mass on Vp (unitless)")         # $PK TVVP = THETA(11)*ALLMV
    e_ffm_vh <- fixed(1)
    label("Allometric exponent of fat-free mass on VH (unitless)")         # $PK VH = 1*ALLMV

    # =======================================================================
    # Random effects: $OMEGA final values. Table E1 prints them as
    # sqrt(exp(omega^2) - 1), e.g. sqrt(exp(1.14544) - 1) = 147% for V.
    # $OMEGA BLOCK(1) ... SAME pairs become occasion-1 / occasion-2 etas with
    # the occasion-2 copy fixed to the same variance.
    # =======================================================================
    etaiov_cl_1 ~ 0.0597849           # Table E1 IOV intrinsic clearance 24.8% (11%); $OMEGA BLOCK(1) 0.0597849 ; 1 IOV in CL (applied to hepatic clearance CLH)
    etaiov_cl_2 ~ fixed(0.0597849)    # $OMEGA BLOCK(1) SAME
    etalvc ~ 1.14544                  # Table E1 IIV central volume 147% (17%); $OMEGA 1.14544 ; 3 IIV in V
    etalka ~ 0.575944                 # Table E1 IIV rate of absorption 88.3% (20%); $OMEGA 0.575944 ; 4 IIV in KA
    etaiov_mtt_1 ~ 0.379271           # Table E1 IOV mean transit time 67.7% (14%); $OMEGA BLOCK(1) 0.379271 ; 5 IOV in MTT
    etaiov_mtt_2 ~ fixed(0.379271)    # $OMEGA BLOCK(1) SAME
    etaiov_logitfdepot_1 ~ 1.0329     # Table E1 IOV bioavailability 134% (11%); $OMEGA BLOCK(1) 1.0329 ; 7 IOV in F (logit scale)
    etaiov_logitfdepot_2 ~ fixed(1.0329)  # $OMEGA BLOCK(1) SAME
    etalq ~ 0.94784                   # Table E1 IIV intercompartmental clearance 126% (23%); $OMEGA 0.94784 ; 9 IIV in Q
    etalppc ~ 0.124736                # Table E1 IIV penetration coefficient 36.3% (14%); $OMEGA 0.124736 ; 10 IIV in PC

    # =======================================================================
    # Residual error: $SIGMA variances converted to SDs.
    # Plasma Y = IPRED*(1+EPS(2)) + EPS(1) with a diagonal $SIGMA is the
    # combined additive + proportional model; CSF Y = IPRED + EPS(3).
    # =======================================================================
    addSd <- sqrt(0.0100628)
    label("Additive residual error for plasma rifampicin (mg/L)")         # Table E1 additive residual error 0.1 mg/L (33%); $SIGMA 0.0100628 ; ADD ERROR
    propSd <- sqrt(0.0572954)
    label("Proportional residual error for plasma rifampicin (fraction)") # Table E1 proportional residual error 24.3% (2.8%) = sqrt(exp(0.0572954)-1); $SIGMA 0.0572954 ; PROP ERROR
    addSd_Ccsf <- sqrt(0.0349855)
    label("Additive residual error for CSF rifampicin (mg/L)")            # Table E1 CSF additive residual error 0.19 mg/L (9%); $SIGMA 0.0349855 ; ADD ERROR CSF
  })

  model({
    # --- 1. Occasion indicators (control stream OCC1 / OCC2 from PERIOD). --
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2
    iov_logitfdepot <- oc1 * etaiov_logitfdepot_1 + oc2 * etaiov_logitfdepot_2


    # --- 3. CSF protein effect on the partition coefficient. ---------------
    #     $PK: PC = THETA(14)*(1 + (LCSFPROT - LOG10(165))/LOG10(165)*THETA(15))
    #     with CSFPROT in mg/dL; CSF_TPRO is carried in g/L.
    csf_tpro_mgdl <- CSF_TPRO * 100
    ppc_prot <- 1 + (log10(csf_tpro_mgdl) - log10(165)) / log10(165) * e_csf_tpro_ppc

    # --- 4. Individual parameters. -----------------------------------------
    #     Autoinduction and the volume change switch on at OCC == 2:
    #     TVCLINTB = THETA(1)*ALLMCL*(1+IND*THETA(8)),
    #     V = TVV*(1-(PERIOD-1)*THETA(9))*EXP(IIVV), VP likewise without eta.
    clint_max <- exp(lclint_max) * (FFM / 44.55)^e_ffm_clint_max * (1 + oc2 * e_occ2_clint_max)
    km <- exp(lkm)
    vc <- exp(lvc + etalvc) * (FFM / 44.55)^e_ffm_vc * (1 - oc2 * e_occ2_vc)
    vp <- exp(lvp) * (FFM / 44.55)^e_ffm_vp * (1 - oc2 * e_occ2_vc)
    q <- exp(lq + etalq) * (FFM / 44.55)^e_ffm_q
    ka <- exp(lka + etalka)
    mtt <- exp(lmtt + iov_mtt)
    ntr <- exp(lntr)
    fdepot <- expit(logitfdepot + iov_logitfdepot)
    qh <- exp(lqh) * (FFM / 44.55)^e_ffm_qh
    vh <- exp(lvh) * (FFM / 44.55)^e_ffm_vh
    ke0 <- exp(lke0)
    ppc <- exp(lppc + etalppc) * ppc_prot

    # --- 5. Well-stirred liver with saturable intrinsic clearance. ---------
    #     $DES: CH = A(2)/VH; CLINT = CLINTB*KM/(CH + KM);
    #           EH = CLINT*FU/(CLINT*FU + QH); CLH = EH*QH*EXP(IOVCL)
    #     The between-occasion eta scales only the eliminated flux CLH, not
    #     the (1 - EH) outflow to plasma, exactly as in the control stream.
    c_liver <- liver / vh
    clint <- clint_max * km / (c_liver + km)
    eh <- clint * fub / (clint * fub + qh)
    clh <- eh * qh * exp(iov_cl)

    # --- 6. ODE system, $DES DADT(1)-DADT(5). ------------------------------
    #     transit() is rxode2's implementation of the Savic 2007 analytical
    #     transit input used by the control stream:
    #       EXP(LOG(BIO*PD) + LOG(KTR) - L + NN*LOG(KTR*T) - KTR*T),
    #     KTR = (NN + 1)/MTT, L = Stirling approximation of log(NN!).
    #     Bioavailability enters as the third argument, exactly as BIO does.
    d/dt(depot) <- transit(ntr, mtt, fdepot) - ka * depot
    d/dt(liver) <- ka * depot - qh * (1 - eh) / vh * liver + qh / vc * central - clh / vh * liver
    d/dt(central) <- qh * (1 - eh) / vh * liver - qh / vc * central - q / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1

    Cc <- central / vc

    # CSF concentration state, no mass transfer: DADT(4) = LOG(2)/HLCSF*(PC*A(3)/V - A(4)).
    d/dt(csf) <- ke0 * (ppc * Cc - csf)

    # F1 = 0 in the control stream: the oral dose enters only through the
    # analytical transit input, never as a bolus into the depot.
    f(depot) <- 0

    # Intravenous doses go into central with F3 = 1 over D3 = 1.5 h. Dose
    # records targeting `central` must set rate = -2 so rxode2 uses dur().
    dur(central) <- exp(ldur)

    # --- 7. Observations. --------------------------------------------------
    Ccsf <- csf
    Cc ~ add(addSd) + prop(propSd)
    Ccsf ~ add(addSd_Ccsf)
  })
}
