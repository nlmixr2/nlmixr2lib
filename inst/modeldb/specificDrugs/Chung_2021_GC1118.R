Chung_2021_GC1118 <- function() {
  description <- paste(
    "Two-compartment target-mediated drug disposition (TMDD) population PK",
    "model for GC1118, a fully human anti-EGFR IgG1 monoclonal antibody, in",
    "adults with advanced solid tumours (Chung 2021; phase I, n = 32, 2-h IV",
    "infusions of 0.3-5 mg/kg once weekly or 8 mg/kg every two weeks).",
    "EGFR is confined to the peripheral compartment; drug-EGFR binding is at",
    "rapid equilibrium (quasi-equilibrium closed form for the complex), the",
    "total peripheral EGFR pool is constant at VP * RB / CLB, and the EGFR",
    "and complex clearances are assumed equal. Linear clearance of free",
    "drug scales with (WT/70)^0.8. The infusion duration is estimated (D1),",
    "a full 9 x 9 IIV block includes IIV on the residual-error magnitude,",
    "and inter-occasion variability (two occasions) sits on Q. The",
    "peripheral EGFR occupancy (RO, %) is returned as a derived output."
  )
  reference <- paste(
    "Chung TK, Lee HA, Park SI, Oh DY, Lee KW, Kim JW, Kim JH, Woo A,",
    "Lee SJ, Bang YJ, Lee H (2021). A target-mediated drug disposition",
    "population pharmacokinetic model of GC1118, a novel anti-EGFR",
    "antibody, in patients with solid tumors. Clinical and Translational",
    "Science 14(3):990-1001. doi:10.1111/cts.12963. Model code from the",
    "Supporting Information (NONMEM control stream CTS-14-990-s001.docx)."
  )
  vignette <- "Chung_2021_GC1118"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(
      analyte = "GC1118 (free)",
      units = "pmol",
      specimen = "serum",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "GC1118 (total: free + EGFR-bound)",
      units = "pmol",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on the linear clearance of free GC1118 (CLA) only,",
        "(WT/70)^0.8. Table 2 'Power term, effect of standardized body",
        "weight on CLA' = 0.8; the 70 kg reference is Methods Eq. 8 text",
        "('70 kg for body weight') and the control stream",
        "'CLA = ... *(WT/70)**THETA(11)'. The stream's comment on THETA(11)",
        "reads 'wt on clb', but the code applies it to CLA, as does Table 2",
        "and the Results. Baseline weight 63.1 +/- 10.4 kg (range",
        "42.5-89.2) per Table 1."
      ),
      source_name = "WT"
    ),
    OCC = list(
      description = "Occasion index for inter-occasion variability on Q (1 or 2)",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Two occasions: the control stream carries one '$OMEGA BLOCK(1)'",
        "IOV on Q followed by a single '$OMEGA BLOCK(1) SAME'",
        "('$ABBREVIATED REPLACE ETA(OCC_Q)=ETA(10,11)'). Chung 2021 does",
        "not state how an occasion was delimited; the study design had two",
        "rich-sampling dosing visits per subject (days 1 and 22 in the",
        "weekly cohorts, days 1 and 43 in the biweekly cohort). Records",
        "with OCC outside {1, 2} carry no IOV."
      ),
      source_name = "OCC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 32L,
    n_studies = 1L,
    n_observations = 793L,
    age_range = "34-72 years",
    age_mean = "56.7 years (SD 9.1)",
    weight_range = "42.5-89.2 kg",
    weight_mean = "63.1 kg (SD 10.4)",
    sex_female_pct = 37.5,
    race_ethnicity = "Korean (single-country phase I trial)",
    disease_state = paste(
      "Advanced solid tumours: colorectal cancer 56.2%; others (ampulla of",
      "Vater, appendix, breast, cholangiocarcinoma, oesophageal,",
      "gallbladder, gastric, nasal cavity, nasopharyngeal, pancreatic,",
      "tonsil) 43.8%. EGFR H-score 91.6 +/- 92.7 (range 0-290)."
    ),
    dose_range = paste(
      "2-h IV infusion of 0.3, 1, 3, 5 or 4 mg/kg once weekly on days 1, 8,",
      "15 and 22 (cohorts 1-5, n = 24) or 8 mg/kg every two weeks",
      "(cohort 6, n = 8)"
    ),
    regions = "South Korea",
    notes = paste(
      "Chung 2021 Table 1. All patients had normal renal and hepatic",
      "function; none developed anti-drug antibodies. Serum free GC1118 by",
      "ELISA (LLOQ 0.025 ug/mL)."
    )
  )

  ini({
    lcl <- log(16.2); label("Linear clearance of free GC1118, CLA, at 70 kg (mL/h)") # Table 2 CLA = 16.2 mL/h
    lvc <- log(3660); label("Central volume of distribution, VC (mL)") # Table 2 VC = 3660.0 mL
    lq <- log(627); label("Intercompartmental clearance of free GC1118, Q (mL/h)") # Table 2 Q = 627.0 mL/h
    lvp <- log(1180); label("Peripheral volume of distribution, VP (mL)") # Table 2 VP = 1180.0 mL
    ld1 <- log(2.2); label("Duration of the zero-order IV input into central, D1 (h)") # Table 2 D1 = 2.2 h
    lksyn <- log(2390); label("EGFR production rate in the peripheral compartment, RB (pmol/h)") # Table 2 RB = 2390.0 pmol/h
    lcl_complex <- log(63.4); label("Clearance of EGFR and of the GC1118-EGFR complex, CLB = CLC (mL/h)") # Table 2 'CLB and CLC' = 63.4 mL/h
    lkd <- fixed(log(0.16)); label("GC1118-EGFR equilibrium dissociation constant, KD (nM = pmol/mL)") # Table 2 KD = 0.16 nM, fixed (in vitro value, Methods ref 6)
    e_wt_cl <- 0.8; label("Power exponent of (WT/70) on CLA (unitless)") # Table 2 power term of standardized body weight on CLA = 0.8

    # IIV: full 9 x 9 block in the order of the control stream $OMEGA
    # BLOCK(9) and of the Table 2 correlation matrix (CLA, VC, Q, VP, D1, RB,
    # CLB, KD, RUV). Variances omega^2 = (CV%/100)^2 from the Table 2 IIV
    # column; covariances r * omega_i * omega_j from the Table 2 correlation
    # matrix after a nearest-positive-definite correction of at most 0.00075
    # per coefficient (the printed 3-decimal matrix has a smallest eigenvalue
    # of -0.00135). The variance convention and the correction are derived
    # and justified in the vignette.
    etalcl + etalvc + etalq + etalvp + etald1 + etalksyn + etalcl_complex + etalkd + etapropSd ~ c(
      0.131044, # CLA IIV 36.2% (Table 2)
      0.00169238, 0.045796, # VC IIV 21.4%; r(CLA,VC) = 0.022
      -0.0762606, -0.109048, 3.2906, # Q IIV 181.4%
      -0.140055, 0.033766, -0.446052, 0.3364, # VP IIV 58%
      0.0160029, 0.0028006, -0.027089, -0.00928948, 0.0036, # D1 IIV 6.0%
      -0.119006, 0.0689884, -0.0155626, 0.210083, -0.00800193, 0.247009, # RB IIV 49.7%
      -0.211703, 0.233042, -0.290194, 0.700378, -0.0147921, 0.752651, 3.53816, # CLB IIV 188.1%
      -0.452973, 0.162431, 1.43049, 0.770381, -0.0497125, 0.910539, 3.62864, 5.18928, # KD IIV 227.8%
      -0.00588526, -0.02131, 0.449898, -0.0269653, 0.000505381, -0.00508233, -0.00940309, 0.295891, 0.098596 # RUV IIV 31.4%
    )

    # IOV on Q, two occasions sharing one variance (control stream $OMEGA
    # BLOCK(1) followed by $OMEGA BLOCK(1) SAME).
    etaiov_q_1 ~ 6.838225 # Table 2 IOV on Q = 261.5% CV; omega^2 = 2.615^2
    etaiov_q_2 ~ fixed(6.838225) # Table 2 IOV on Q = 261.5% CV; shared variance (BLOCK(1) SAME)

    propSd <- 0.10; label("Proportional residual error (fraction)") # Table 2 'Proportional RUV' = 0.10
    addSd <- fixed(0.001); label("Additive residual error (ug/mL)") # Supplement control stream $THETA '(0.001) FIX ; Additive error' (not in Table 2)
  })

  model({
    # Molecular weight 150 kDa = 0.15 ug/pmol (control stream 'MWT = 0.15').
    # Dose in mg is converted to pmol (F1 = 1/MWT in the stream, which dosed
    # in ug); concentrations are returned in ug/mL.
    mw <- 0.15

    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_q <- oc1 * etaiov_q_1 + oc2 * etaiov_q_2

    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq + iov_q)
    vp <- exp(lvp + etalvp)
    d1 <- exp(ld1 + etald1)
    ksyn <- exp(lksyn + etalksyn)
    cl_complex <- exp(lcl_complex + etalcl_complex)
    kd <- exp(lkd + etalkd)

    # Eq. 5: constant total peripheral EGFR (pmol), with CLB = CLC
    rtot <- vp * ksyn / cl_complex
    # Eq. 4: quasi-equilibrium GC1118-EGFR complex in the peripheral
    # compartment (pmol); peripheral1 is the total peripheral drug TAP
    bterm <- kd * vp + peripheral1 + rtot
    complex <- (bterm - sqrt(bterm^2 - 4 * peripheral1 * rtot)) / 2
    # Eq. 3: free peripheral drug (pmol)
    free_p <- peripheral1 - complex

    # Eqs. 1-2
    d/dt(central) <- -cl * central / vc - q * central / vc + q * free_p / vp
    d/dt(peripheral1) <- q * central / vc - q * free_p / vp - complex * cl_complex / vp

    f(central) <- 1000 / mw
    dur(central) <- d1

    # Eq. 9: peripheral EGFR occupancy (%)
    RO <- 100 * complex / rtot

    # Residual error: W = SQRT(THETA(9)^2 + THETA(10)^2 * IPRED^2), Y = IPRED
    # + W * EPS(1) * EXP(ETA(9)); IIV on the residual magnitude scales both
    # components.
    propSdInd <- propSd * exp(etapropSd)
    addSdInd <- addSd * exp(etapropSd)

    Cc <- mw * central / vc
    Cc ~ add(addSdInd) + prop(propSdInd)
  })
}
