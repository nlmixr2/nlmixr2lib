Upton_2025_clofazimine <- function() {
  description <- paste(
    "Three-compartment population PK model for oral clofazimine in plasma",
    "and lumbar cerebrospinal fluid (CSF) in adults with pulmonary",
    "drug-resistant tuberculosis without tuberculous meningitis (Upton",
    "2025). The plasma model is the Abdelwahab 2020 model with every",
    "disposition, absorption and variability parameter fixed to its",
    "published value and re-expressed at this cohort's typical weight",
    "52 kg, fat-free mass 37.21 kg and fat mass 13.019 kg; only a relative",
    "bioavailability of 1.26 was re-estimated. Savic transit absorption",
    "(MTT 1.41 h, 4.75 transit compartments, ka 0.209 1/h), clearance",
    "11.1 L/h, central volume 240 L (fat-free mass), deep peripheral",
    "volume 12100 L scaled on fat mass, shallow peripheral volume 864 L,",
    "and intercompartmental clearances 54.7 and 83.5 L/h (total weight).",
    "CSF is a Sheiner-style effect compartment holding a concentration",
    "that equilibrates with plasma at an equilibration half-life of 55.4 h",
    "toward a CSF-to-plasma pseudo-partition coefficient of 0.0013 (total",
    "concentrations)."
  )
  reference <- paste(
    "Upton CM, Calderin JM, Diacon AH, Boeree MJ, Denti P, Wiesner L,",
    "Kellermann T, McCulloch M, Aarnoutse R, Wasserman S (2025).",
    "Cerebrospinal fluid penetration of cycloserine/terizidone and",
    "clofazimine in patients with pulmonary TB.",
    "Antimicrob Agents Chemother 69(12):e00931-25.",
    "doi:10.1128/aac.00931-25.",
    "Plasma model structure from: Abdelwahab MT, Wasserman S, Brust JCM,",
    "Gandhi NR, Meintjes G, Everitt D, Diacon A, Dawson R, Wiesner L,",
    "Svensson EM, Maartens G, Denti P (2020). Clofazimine pharmacokinetics",
    "in patients with TB: dosing implications. J Antimicrob Chemother",
    "75(11):3269. doi:10.1093/jac/dkaa310",
    sep = " "
  )
  vignette <- "Upton_2025_cycloserine_clofazimine"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The `csf` state holds a concentration, not an amount: the Supplement S1
  # clofazimine CSF control stream reads CE = A(5) directly and integrates
  # DADT(5) = KE0_CSF*(PPC_CSF*C2 - A(5)).
  compartmentData <- list(
    depot = list(
      analyte = "clofazimine",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "clofazimine",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "clofazimine",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    peripheral2 = list(
      analyte = "clofazimine",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    csf = list(
      analyte = "clofazimine",
      units = "mg/L",
      specimen = "CSF",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric descriptor for clearance and both intercompartmental",
        "clearances (fixed exponent 0.75) and the shallow peripheral",
        "volume Vp2 (fixed exponent 1), per the Abdelwahab 2020 control",
        "stream (TVCL, TVQ, TVQ2 x ALLMCL_WT; TVV4 x ALLMV_WT). Reference",
        "52 kg, the Upton 2025 cohort typical weight (Table S3 footnote d;",
        "TVWT = 52 in the same paper's cycloserine control stream)."
      ),
      source_name = "WT"
    ),
    FFM = list(
      description = "Fat-free mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric descriptor for the central volume (fixed exponent 1;",
        "Abdelwahab 2020 TVV = THETA(9)*ALLMV_FFM). Reference 37.21 kg,",
        "the Upton 2025 cohort typical value (Table S3 footnote e, 'fat-free",
        "mass of 37 kg'; TVFFM = 37.21 in the cycloserine control stream).",
        "Janmahasatian formula: FFM = WHSMAX * HT^2 * WT / (WHS50 * HT^2 +",
        "WT), HT in metres, WHSMAX / WHS50 = 37.99 / 35.98 for women and",
        "42.92 / 30.93 for men."
      ),
      source_name = "FFM"
    ),
    FM = list(
      description = "Fat mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Fat mass FM = WT - FFM (Janmahasatian FFM), the Abdelwahab 2020",
        "FAT variable. Allometric descriptor for the deep peripheral volume",
        "Vp1 (fixed exponent 1; TVV3 = THETA(1)*ALLMV_FAT), which carries",
        "the sex difference in clofazimine exposure. Reference 13.019 kg,",
        "the Upton 2025 cohort typical value (Table S3 footnote f, 'fat",
        "mass of 13 kg'; TVFAT = 13.019 in the cycloserine control",
        "stream). Must be positive: the deep volume is proportional to it."
      ),
      source_name = "FAT"
    ),
    OCC = list(
      description = "Dosing occasion for between-occasion variability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Abdelwahab 2020 carries between-occasion variability on",
        "bioavailability, ka and mean transit time over study-specific",
        "occasion codes (1, 2, 3, 10, 20, 30, ...). Upton 2025 does not",
        "print its clofazimine plasma control stream; two occasions are",
        "encoded here to match its sampling design (pre-dose sample after",
        "the self-reported dose, then the observed-dose profile), the same",
        "split as the cycloserine model. Any other OCC value switches the",
        "between-occasion variability off."
      ),
      source_name = "OCC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 23L,
    n_studies = 1L,
    age_range = "20-60 years (median 36, whole cohort of 27)",
    weight_range = "30-73 kg (median 52, whole cohort of 27)",
    sex_female_pct = 59,
    race_ethnicity = "Not reported (South African cohort, Cape Town).",
    disease_state = paste(
      "Pulmonary drug-resistant tuberculosis without tuberculous",
      "meningitis; intact blood-CSF barrier by CSF:serum albumin ratio.",
      "48% of the whole cohort living with HIV."
    ),
    dose_range = paste(
      "Clofazimine 100 mg once daily for at least 2 weeks (up to almost",
      "2 years) before the PK day, co-administered with terizidone."
    ),
    regions = "South Africa",
    notes = paste(
      "23 of the 27 participants received clofazimine and contributed",
      "184 plasma and 23 CSF samples (one lumbar puncture each at 2, 4, 6",
      "or 8 h post-dose). Plasma LLOQ 0.00781 mg/L, CSF LLOQ 0.00005 mg/L;",
      "all concentrations above the LLOQ. Plasma parameters were fixed to",
      "the Abdelwahab 2020 model (139 South African TB patients) apart",
      "from a re-estimated relative bioavailability; individual plasma",
      "parameters were then fixed and the CSF parameters estimated",
      "(NONMEM 7.5.1, FOCE-I, 95% CI by SIR). Final estimates are",
      "Table S3."
    )
  )

  ini({
    # --- Plasma structural parameters. Every value except the relative
    # bioavailability is Table S3 'Fixed', footnote b: 'fixed at the values
    # reported by the previously published model' (Abdelwahab 2020). Table
    # S3 prints them already re-expressed at this cohort's typical body
    # size (footnotes d, e, f). The Table S3 values are used as printed;
    # rescaling the rounded Abdelwahab 2020 Table 2 values reproduces them
    # to within 1-2% (Q1 56.3 x (52/54)^0.75 = 54.7; Vc 262 x 37/40.3 =
    # 240.5; CL 11.5 x (52/54)^0.75 = 11.2 against the printed 11.1).
    lcl <- fixed(log(11.1))
    label("Clearance at WT = 52 kg (L/h)")                                     # Table S3 CL 11.1 Fixed (footnotes b, d)
    lvc <- fixed(log(240))
    label("Central volume of distribution at FFM = 37.21 kg (L)")              # Table S3 Vc 240 Fixed (footnotes b, e)
    lq <- fixed(log(54.7))
    label("Intercompartmental clearance to the deep peripheral compartment at WT = 52 kg (L/h)")  # Table S3 Q1 54.7 Fixed (b, d)
    lvp <- fixed(log(12100))
    label("Deep peripheral volume Vp1 at fat mass = 13.019 kg (L)")            # Table S3 Vp1 12100 Fixed (b; footnote f gives the fat-mass reference)
    lq2 <- fixed(log(83.5))
    label("Intercompartmental clearance to the shallow peripheral compartment at WT = 52 kg (L/h)")  # Table S3 Q2 83.5 Fixed (b, d)
    lvp2 <- fixed(log(864))
    label("Shallow peripheral volume Vp2 at WT = 52 kg (L)")                   # Table S3 Vp2 864 Fixed (b, d)
    lka <- fixed(log(0.209))
    label("First-order absorption rate constant from depot to central (1/h)")  # Table S3 Ka 0.209 Fixed (b); Abdelwahab 2020 Table 2
    lmtt <- fixed(log(1.41))
    label("Mean absorption transit time (h)")                                  # Table S3 MTT 1.41 Fixed (b); Abdelwahab 2020 Table 2
    lntr <- fixed(log(4.75))
    label("Number of absorption transit compartments (unitless)")              # Table S3 NN 4.75 Fixed (b); Abdelwahab 2020 Table 2

    # Table S3 lists 'Bioavailability, F 1 Fixed' and a separately
    # estimated 'Scaling factor for F 1.26 (1.09-1.47)' (footnote g); the
    # product, the bioavailability in this cohort relative to the
    # Abdelwahab 2020 reference, is carried as the single lfdepot.
    lfdepot <- log(1.26)
    label("Relative bioavailability versus the Abdelwahab 2020 cohort (fraction)")  # Table S3 Scaling factor for F 1.26 (1.09-1.47) x F 1 Fixed; Results dOFV -9.05

    # --- Allometric exponents (Abdelwahab 2020 Supplementary data:
    # 'theta was fixed to 0.75 for clearance parameters (CL, Q1, Q2) and to
    # 1 for volume parameters (Vc, Vp1, Vp2)').
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of total body weight on clearance (unitless)")              # Abdelwahab 2020 ALLMCL_WT = (WT/TVWT)**0.75
    e_wt_q <- fixed(0.75)
    label("Allometric exponent of total body weight on Q1 (unitless)")                     # Abdelwahab 2020 TVQ = THETA(2)*ALLMCL_WT
    e_wt_q2 <- fixed(0.75)
    label("Allometric exponent of total body weight on Q2 (unitless)")                     # Abdelwahab 2020 TVQ2 = THETA(14)*ALLMCL_WT
    e_wt_vp2 <- fixed(1)
    label("Allometric exponent of total body weight on the shallow peripheral volume (unitless)")  # Abdelwahab 2020 TVV4 = THETA(15)*ALLMV_WT
    e_ffm_vc <- fixed(1)
    label("Allometric exponent of fat-free mass on central volume (unitless)")             # Abdelwahab 2020 TVV = THETA(9)*ALLMV_FFM
    e_fm_vp <- fixed(1)
    label("Allometric exponent of fat mass on the deep peripheral volume (unitless)")      # Abdelwahab 2020 TVV3 = THETA(1)*ALLMV_FAT

    # --- CSF effect-compartment parameters, estimated in Upton 2025.
    # ke0 = log(2) / 55.4 = 0.0125 1/h, the S1 CSF $THETA 4 value.
    lke0 <- log(log(2) / 55.4)
    label("Plasma-to-CSF equilibration rate constant (1/h)")                   # Table S3 t1/2 Plasma-CSF 55.4 h (22.3-321); S1 CSF $THETA 4 KE0_CSF = 0.0125
    lppc <- log(0.0013)
    label("CSF-to-plasma pseudo-partition coefficient (fraction)")             # Table S3 PPC CSF-Plasma 1.30x10-3 (1.13x10-3 to 1.47x10-3); S1 CSF $THETA 3 = 0.0013

    # --- Between-subject and between-occasion variability, all fixed
    # (Table S3 footnote b). The percentages are Abdelwahab 2020's
    # 'approximate %CV', i.e. the square root of the log-scale variance
    # (the same convention Table S3 follows for cycloserine, where
    # sqrt(0.0653) = 25.6%), so each variance is (percentage / 100)^2.
    etalcl ~ fixed(0.256^2)
    label("Between-subject variability in clearance (log-scale variance)")     # Table S3 BSV CL 25.6% Fixed; 0.256^2
    etalvc ~ fixed(0.235^2)
    label("Between-subject variability in central volume (log-scale variance)")  # Table S3 BSV Vc 23.5% Fixed; 0.235^2
    etalvp ~ fixed(0.296^2)
    label("Between-subject variability in the deep peripheral volume (log-scale variance)")  # Table S3 BSV Vp1 29.6% Fixed; 0.296^2
    etalvp2 ~ fixed(0.546^2)
    label("Between-subject variability in the shallow peripheral volume (log-scale variance)")  # Table S3 BSV Vp2 54.6% Fixed; 0.546^2
    etalfdepot ~ fixed(0.301^2)
    label("Between-subject variability in bioavailability (log-scale variance)")  # Table S3 BSV F 30.1% Fixed; 0.301^2

    etaiov_fdepot_1 ~ fixed(0.354^2)
    label("Between-occasion variability in bioavailability, occasion 1 (log-scale variance)")  # Table S3 BOV F 35.4% Fixed; 0.354^2
    etaiov_fdepot_2 ~ fixed(0.354^2)                                          # same variance on every occasion ($OMEGA BLOCK(1) SAME)
    etaiov_mtt_1 ~ fixed(0.466^2)
    label("Between-occasion variability in mean transit time, occasion 1 (log-scale variance)")  # Table S3 BOV MTT 46.6% Fixed; 0.466^2
    etaiov_mtt_2 ~ fixed(0.466^2)                                             # same variance on every occasion ($OMEGA BLOCK(1) SAME)
    etaiov_ka_1 ~ fixed(0.326^2)
    label("Between-occasion variability in ka, occasion 1 (log-scale variance)")  # Table S3 BOV Ka 32.6% Fixed; 0.326^2
    etaiov_ka_2 ~ fixed(0.326^2)                                              # same variance on every occasion ($OMEGA BLOCK(1) SAME)

    # --- Residual error, combined W = sqrt(ADD^2 + PROP^2). The plasma
    # additive SD is 20% of the 0.00781 mg/L plasma LLOQ and the CSF
    # additive SD 20% of the 0.00005 mg/L CSF LLOQ (Table S3 footnote h).
    propSd <- fixed(0.114)
    label("Proportional residual error for plasma clofazimine (fraction)")     # Table S3 plasma proportional 11.4% Fixed (b)
    addSd <- fixed(0.00156)
    label("Additive residual error for plasma clofazimine (mg/L)")             # Table S3 plasma additive 1.56x10-3 Fixed
    propSd_Ccsf <- 0.342
    label("Proportional residual error for CSF clofazimine (fraction)")        # Table S3 CSF proportional 34.2% (26.4-49.1); S1 CSF $THETA 5 = 0.342
    addSd_Ccsf <- fixed(0.00001)
    label("Additive residual error for CSF clofazimine (mg/L)")                # Table S3 CSF additive 1.00x10-5 Fixed, footnote h
  })

  model({
    # 1. Occasion indicators.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2

    # 2. Individual parameters (Abdelwahab 2020 control stream $PK, with
    # the reference sizes moved to the Upton 2025 cohort typical values).
    cl <- exp(lcl + etalcl) * (WT / 52)^e_wt_cl
    vc <- exp(lvc + etalvc) * (FFM / 37.21)^e_ffm_vc
    q <- exp(lq) * (WT / 52)^e_wt_q
    vp <- exp(lvp + etalvp) * (FM / 13.019)^e_fm_vp
    q2 <- exp(lq2) * (WT / 52)^e_wt_q2
    vp2 <- exp(lvp2 + etalvp2) * (WT / 52)^e_wt_vp2
    ka <- exp(lka + iov_ka)
    mtt <- exp(lmtt + iov_mtt)
    ntr <- exp(lntr)
    fdepot <- exp(lfdepot + etalfdepot + iov_fdepot)
    ke0 <- exp(lke0)
    ppc <- exp(lppc)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # 3. ODE system. Savic transit input as in Upton_2025_cycloserine
    # (F1 = 0 in the control streams, reproduced by f(depot) <- 0).
    d/dt(depot) <- transit(ntr, mtt, fdepot) - ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    Cc <- central / vc
    d/dt(csf) <- ke0 * (ppc * Cc - csf)
    f(depot) <- 0

    # 4. Observations.
    Ccsf <- csf
    Cc ~ add(addSd) + prop(propSd)
    Ccsf ~ add(addSd_Ccsf) + prop(propSd_Ccsf)
  })
}
