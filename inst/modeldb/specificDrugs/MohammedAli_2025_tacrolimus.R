MohammedAli_2025_tacrolimus <- function() {
  description <- paste(
    "Two-compartment population PK model jointly describing both oral",
    "tacrolimus formulations -- twice-daily immediate-release IR-Tac (Prograf)",
    "and once-daily extended-release LCP-Tac (Envarsus, MeltDose technology) --",
    "in the same stable adult renal transplant recipients studied in a",
    "within-patient conversion design (Mohammed Ali 2025). Parameterized in",
    "apparent elimination clearance CL/F, apparent distributional clearance",
    "CLD/F and apparent central and peripheral volumes Vc/F and Vp/F, with",
    "Vp/F fixed at 500 L. Delayed first-order absorption uses a separate",
    "absorption rate constant and a separate lag time for each formulation",
    "(ka 2.04 vs 0.111 per h; lag 0.465 vs 1.4 h). A 24 h cosinor rhythm rides",
    "on CL/F (acrophase 17:00, amplitude 3.42 L/h) and, for IR-Tac only, on its",
    "absorption rate constant (acrophase 03:08, amplitude 1.55 per h), so the",
    "model time variable must be clock time measured from midnight. CYP3A5",
    "expresser status crossed with formulation gives four relative",
    "bioavailabilities anchored at F = 1 for the LCP-Tac nonexpresser reference",
    "group, which is the pharmacogenetic basis for the paper's genotype-specific",
    "IR-Tac to LCP-Tac dose conversion ratios (1:0.6 for expressers, 1:0.7 for",
    "nonexpressers). Inter-individual variability is carried on CL/F and on a",
    "correlated block over Vc/F and the two formulation-specific absorption",
    "rate constants, with inter-occasion variability on CL/F and Vc/F and a",
    "proportional residual error on whole-blood concentrations.",
    sep = " "
  )
  reference <- paste(
    "Mohammed Ali Z, Fernandez-Alarcon B, Fontova P, Vidal-Alabro A,",
    "Rigo-Bonnin R, Melilli E, Montero N, Manonelles A, Coloma A, Fava A,",
    "Grinyo JM, Cruzado JM, Colom H, Lloberas N. Optimizing Dose Conversion",
    "from IR-Tac to LCP-Tac Formulations in Renal Transplant Recipients: A",
    "Population Pharmacokinetic Modeling Study. Pharmaceutics.",
    "2025;17(9):1185. doi:10.3390/pharmaceutics17091185.",
    sep = " "
  )
  vignette <- "MohammedAli_2025_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot       = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    CYP3A5_EXPR = list(
      description        = "Recipient CYP3A5 expresser status (rs776746)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5*3/*3 nonexpresser)",
      notes              = paste(
        "1 = at least one functional CYP3A5*1 allele (genotype *1/*3 or *1/*1);",
        "0 = CYP3A5*3/*3. Genotyped for CYP3A5*3 G>A (rs776746) by TaqMan SNP",
        "assay (Mohammed Ali 2025 Methods section 2.4). The paper labels",
        "expressers 'HM' and nonexpressers 'PM' in Table 3 and Table S1.",
        "Cohort distribution (Table 1): *1/*3 9 (30%), *1/*1 1 (3%),",
        "*3/*3 20 (67%), i.e. CYP3A5_EXPR = 1 in 10 of 30 subjects; the single",
        "*1/*1 patient was pooled with the *1/*3 heterozygotes because n = 1",
        "could not support its own stratum (Results section 3.1). Enters the",
        "final model only through the relative bioavailability F, crossed with",
        "formulation: CYP3A5 was significant on F rather than on CL/F",
        "(Results section 3.2.2, delta-MOFV = -51 units when added to F).",
        sep = " "
      ),
      source_name        = "CYP3A5 genotype"
    ),
    FORM_TAC_ENVARSUS = list(
      description        = "LCP-Tac (Envarsus, MeltDose) extended-release tacrolimus formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (IR-Tac, twice-daily immediate-release Prograf)",
      notes              = paste(
        "1 = record is on once-daily LCP-Tac (Envarsus, Chiesi Farmaceutici;",
        "life-cycle-pharma MeltDose melt-extrusion tablet); 0 = record is on",
        "twice-daily immediate-release IR-Tac (Prograf, Astellas). Per-occasion",
        "in this within-patient conversion design: every patient contributed",
        "one 24 h IR-Tac profile one week before conversion and one 24 h",
        "LCP-Tac profile four weeks after conversion (Methods section 2.2), so",
        "each subject supplies both levels. Switches three things at once in",
        "model(): the absorption rate constant, the absorption lag time, and",
        "the relative bioavailability F (the latter jointly with",
        "CYP3A5_EXPR). Registered as the sibling canonical that the",
        "FORM_TAC_IR entry of inst/references/covariate-columns.md directs",
        "future Envarsus / LCP-Tacro models to use rather than overloading",
        "FORM_TAC_IR, whose 0 level is the Advagraf / Astagraf XL",
        "prolonged-release product and not LCP-Tac.",
        sep = " "
      ),
      source_name        = "formulation"
    ),
    OCC = list(
      description        = "Sampling-occasion index for inter-occasion variability on CL/F and Vc/F",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Integer occasion index decomposed inside model() into",
        "mutually-exclusive binary indicators multiplexing per-occasion etas on",
        "log CL/F and log Vc/F. Mohammed Ali 2025 reports single IOV magnitudes",
        "(Table 3 'IOV CL 20.85' and 'IOV Vc 58.82') without stating the number",
        "of occasions. The design fixes it at two: each patient contributed",
        "exactly two full 24 h steady-state profiles, one per formulation",
        "(Methods section 2.2), so occasion 1 = the pre-conversion IR-Tac",
        "profile and occasion 2 = the post-conversion LCP-Tac profile. Occasion",
        "2 has its variance fixed equal to occasion 1 (the NONMEM $OMEGA",
        "BLOCK(1) SAME pattern; nlmixr2 has no SAME shortcut). Because OCC and",
        "FORM_TAC_ENVARSUS are one-to-one in this design, the IOV and the",
        "formulation effect on F are partially aliased -- see the vignette",
        "Assumptions and deviations section. Records with OCC outside 1..2",
        "carry no IOV.",
        sep = " "
      ),
      source_name        = "occasion"
    )
  )

  # Covariates screened by Mohammed Ali 2025 but NOT retained in the final
  # model. Documented for provenance only; none is referenced in model().
  # Results section 3.2.2: 'Graphical exploration of Bayesian estimates of the
  # pharmacokinetic parameters vs. demographic and biochemical covariates did
  # not show any significant trend. When covariates were entered univariately,
  # none of the size descriptors (body weight, body mass index) entered
  # allometrically or with any other relationship provided a significant drop
  # in the MOFV (p > 0.05) or improved the overall model. Similarly, this,
  # occurred with age and hematocrit.'
  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Tested univariately, allometrically and with other functional forms;",
        "no significant MOFV drop (p > 0.05) (Results section 3.2.2). The",
        "Discussion attributes the absence of any covariate effect other than",
        "CYP3A5 to the 30-patient sample size and the narrow covariate range",
        "in this stable cohort. Cohort: 72 kg mean (IQR 64-80) on IR-Tac,",
        "73 kg (IQR 64-80) on LCP-Tac (Table 1).",
        sep = " "
      )
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste(
        "Tested univariately as a size descriptor; no significant MOFV drop",
        "(Results section 3.2.2). Cohort: 26 kg/m^2 mean (IQR 21.5-29.3) on",
        "IR-Tac, 27 (IQR 21.5-29.3) on LCP-Tac (Table 1).",
        sep = " "
      )
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Tested univariately; no significant MOFV drop (Results section",
        "3.2.2). Cohort: 58 years mean (IQR 48-68) (Table 1).",
        sep = " "
      )
    ),
    HCT = list(
      description = "Hematocrit",
      units       = "%",
      type        = "continuous",
      notes       = paste(
        "Recorded at each sampling occasion (Methods section 2.2) and tested",
        "univariately; no significant MOFV drop (Results section 3.2.2).",
        "Cohort: 40.9% mean (IQR 37.6-44.8) on IR-Tac, 40.1% (IQR 37.1-43) on",
        "LCP-Tac (Table 1). Tacrolimus partitions strongly into erythrocytes,",
        "so hematocrit is a standard candidate covariate on apparent",
        "whole-blood clearance; it was not supported here.",
        sep = " "
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Listed among the demographic covariates tested for influence on",
        "tacrolimus pharmacokinetics (Methods section 2.6.2, 'gender'); not",
        "retained in the final model. Cohort: 8 of 30 female (27%) (Table 1).",
        sep = " "
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 30,
    n_studies      = 1,
    n_observations = 932,
    age_range      = "58 years mean (IQR 48-68)",
    weight_range   = "72 kg mean (IQR 64-80)",
    sex_female_pct = 26.7,
    disease_state  = paste(
      "stable adult renal transplant recipients at least 6 months",
      "post-transplant, on triple immunosuppression (tacrolimus +",
      "mycophenolate mofetil + prednisone), converted within-patient from",
      "twice-daily IR-Tac (Prograf) to once-daily LCP-Tac (Envarsus)",
      sep = " "
    ),
    dose_range     = paste(
      "IR-Tac median 5 mg/day (range 3-12) in CYP3A5 expressers and 3 mg/day",
      "(range 1.5-8) in nonexpressers; LCP-Tac median 3.75 mg/day (range",
      "2-8.5) and 2 mg/day (range 1-4.75) respectively (Table 2)",
      sep = " "
    ),
    regions        = "Spain (single centre, Hospital Universitari de Bellvitge, Barcelona)",
    renal_function = "eGFR (CKD-EPI) 49.6 mL/min mean (IQR 34-57) on IR-Tac; serum creatinine 141.9 umol/L mean (IQR 108-166)",
    hematocrit     = "40.9% mean (IQR 37.6-44.8) on IR-Tac, 40.1% (IQR 37.1-43) on LCP-Tac",
    genotypes      = paste(
      "CYP3A5 *1/*3 9 (30%), *1/*1 1 (3%), *3/*3 20 (67%); analysed as 10",
      "CYP3A5*1 expressers and 20 nonexpressers",
      sep = " "
    ),
    notes          = paste(
      "Mohammed Ali 2025 Tables 1 and 2. Open-label, prospective,",
      "non-randomized, investigator-initiated single-centre trial",
      "NCT02961608 (ethics ref. PR175/18). Each patient contributed two rich",
      "24 h steady-state profiles, one week before conversion on IR-Tac (481",
      "samples) and four weeks after conversion on LCP-Tac (451 samples), for",
      "932 whole-blood concentrations total; 10-18 samples per patient per",
      "profile at pre-dose and 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 12.5, 13, 13.5,",
      "14, 15, 20 and 24 h. Whole-blood tacrolimus by validated LC-MS/MS with",
      "a 1.0 ng/mL lower limit of quantitation. Exclusions included pregnancy,",
      "active infection, HIV, neoplasm, severe gastrointestinal disease,",
      "hepatitis B or C, and CYP3A-interacting comedication.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Disposition (Mohammed Ali 2025 Table 3, 'Disposition PK parameters').
    # CL/F is the MESOR of the 24 h cosinor rhythm defined below, and it is
    # the value for the reference formulation/genotype group (LCP-Tac in
    # CYP3A5*1 nonexpressers, whose F is fixed at 1). Because F is a pure
    # multiplier on dose, the apparent clearance of the other three groups is
    # 11.9 / F: the Discussion's quoted values of 15.97, 17.17 and 27.87 L/h
    # for IR-Tac nonexpressers, LCP-Tac expressers and IR-Tac expressers
    # reproduce exactly as 11.9/0.745, 11.9/0.693 and 11.9/0.427.
    # =====================================================================
    lcl <- log(11.9)
    label("Mesor of apparent elimination clearance CL/F (L/h)")
    # Table 3 row 'CL/F (L.h-1) = 11.9 (RSE 8.5%)'; bootstrap median 11.85 (90% CI 10.34-13.53)

    lvc <- log(78)
    label("Apparent central volume of distribution Vc/F (L)")
    # Table 3 row 'Vc/F (L) = 78 (RSE 14.7%)'; bootstrap median 81 (90% CI 63-100.22)

    lq <- log(25.8)
    label("Apparent distributional clearance CLd/F (L/h)")
    # Table 3 row 'CLd/F (L.h-1) = 25.8 (RSE 8.5%)'; bootstrap median 25.75 (90% CI 22.08-29.39)

    lvp <- fixed(log(500))
    label("Apparent peripheral volume of distribution Vp/F (L)")
    # Table 3 row 'Vp/F (L) = 500 FIX'. Results section 3.2.1: 'The peripheral
    # compartment distribution volume had to be fixed to the estimated amount
    # from the model, a value which is similar with our previous model ...
    # to increase the estimation precision of the remaining parameters of the
    # model and to avoid collinearities.'

    # =====================================================================
    # Absorption -- a separate first-order rate constant and a separate lag
    # time per formulation (Mohammed Ali 2025 Table 3, 'Absorption
    # parameters'). Results section 3.2.1: two distinct ka and two distinct
    # lag times were both statistically significant (delta-MOFV = -411 and
    # -196 units) and the two-ka structure cut IIV on ka by 47%. Transit
    # compartment models were tested and not supported by these data
    # (Discussion: 'probably due to overparameterization'), unlike the same
    # group's LCP-Tac-only model in MohammedAli_2023_tacrolimus.R.
    # =====================================================================
    lka_ir <- log(2.04)
    label("Mesor of the IR-Tac first-order absorption rate constant ka (1/h)")
    # Table 3 row 'Ka IR-Tac = 2.04 (RSE 40%)'; bootstrap median 2.17 (90% CI 1.23-3.72)

    lka_lcp <- log(0.111)
    label("LCP-Tac first-order absorption rate constant ka (1/h)")
    # Table 3 row 'Ka LCP_Tac = 0.111 (RSE 16.9%)'; bootstrap median 0.115 (90% CI 0.08-0.15)

    ltlag_ir <- log(0.465)
    label("IR-Tac absorption lag time (h)")
    # Table 3 row 'Lag-Time IR-Tac (h) = 0.465 (RSE 0.1%)'; bootstrap median 0.465 (90% CI 0.42-0.47)

    ltlag_lcp <- log(1.4)
    label("LCP-Tac absorption lag time (h)")
    # Table 3 row 'Lag-Time LCP-Tac (h) = 1.4 (RSE 2.4%)'; bootstrap median 1.39
    # (90% CI 1.32-1.57). The Discussion quotes this as 1.42 h; Table 3 and the
    # bootstrap CI both support 1.4, so the Table 3 value is used.

    # =====================================================================
    # Relative bioavailability F: formulation crossed with CYP3A5 expresser
    # status, four groups. Mohammed Ali 2025 Methods section 2.6.1: 'F value
    # was fixed to 1 for the combination of formulation and genetic variant
    # group taken as reference. In the other cases, the relative
    # bioavailability with respect to the reference group was estimated as
    # follows: F = 1 * theta_x' (equation 2). The reference group is LCP-Tac
    # in CYP3A5*1 nonexpressers. Within each genotype F is higher for LCP-Tac
    # than IR-Tac (the MeltDose bioavailability gain); within each formulation
    # F is lower in expressers (Results section 3.2.2).
    # =====================================================================
    lfdepot_lcp_nonexpr <- fixed(log(1))
    label("Relative bioavailability F, LCP-Tac in CYP3A5 nonexpressers (reference, fraction)")
    # Table 3 row 'F LCP-Tac_PM = 1 FIX'

    lfdepot_ir_nonexpr <- log(0.745)
    label("Relative bioavailability F, IR-Tac in CYP3A5 nonexpressers (fraction)")
    # Table 3 row 'F IR-Tac_PM = 0.745 (RSE 7.6%)'; bootstrap median 0.757 (90% CI 0.66-0.84)

    lfdepot_lcp_expr <- log(0.693)
    label("Relative bioavailability F, LCP-Tac in CYP3A5 expressers (fraction)")
    # Table 3 row 'F LCP-Tac_HM = 0.693 (RSE 13.7%)'; bootstrap median 0.695 (90% CI 0.52-0.85)

    lfdepot_ir_expr <- log(0.427)
    label("Relative bioavailability F, IR-Tac in CYP3A5 expressers (fraction)")
    # Table 3 row 'F IR-Tac_HM = 0.427 (RSE 13.4%)'; bootstrap median 0.428 (90% CI 0.34-0.52)

    # =====================================================================
    # 24 h circadian (cosinor) rhythm parameters, Mohammed Ali 2025 Table 3
    # 'Circadian rhythms parameters' and equation 3:
    #   P = theta_1 + theta_AMP * cos((2*pi/1440) * (TIME - ACROPHASE))
    # with a 24 h (1440 min) period, theta_1 the mesor, theta_AMP the
    # amplitude in the parameter's own units, and TIME measured in minutes
    # from midnight of the first PK profile. Table 3 reports both acrophases
    # in HOURS, so the implementation below uses a 24 h period with time in
    # hours, which is the same function. The amplitudes are absolute (an
    # amplitude of 3.42 read as a FRACTION of CL/F would be 342% and would
    # drive clearance negative), and both are smaller than their mesor
    # (3.42 < 11.9 and 1.55 < 2.04), so the rhythm never changes sign.
    # =====================================================================
    acrophase_cl <- 17
    label("Acrophase (clock time of peak) of the CL/F circadian rhythm (h)")
    # Table 3 row 'Acrophase CL/F (h) = 17 (RSE 3.6%)'; bootstrap median 16.94 (90% CI 15.94-17.98)

    amp_cl <- 3.42
    label("Amplitude of the CL/F circadian rhythm (L/h)")
    # Table 3 row 'Amp CL/F = 3.42 (RSE 17.1%)'; bootstrap median 3.41 (90% CI 2.33-4.39)

    acrophase_ka_ir <- 3.13
    label("Acrophase (clock time of peak) of the IR-Tac ka circadian rhythm (h)")
    # Table 3 row 'Acrophase ka (h) = 3.13 (RSE 18.3%)'; bootstrap median 3.17 (90% CI 1.82-4.52)

    amp_ka_ir <- 1.55
    label("Amplitude of the IR-Tac ka circadian rhythm (1/h)")
    # Table 3 row 'Amp ka = 1.55 (RSE 44.5%)'; bootstrap median 1.64 (90% CI 0.91-2.97)

    # =====================================================================
    # Inter-individual variability. Table 3 reports IIV as a CV percentage;
    # converted to the log-scale variance with omega^2 = log(CV^2 + 1).
    # Results section 3.2.1: 'A partial OMEGA block structure with an OMEGA
    # block on Vc/F, Ka IR-Tac, and Ka LCP-Tac was the most appropriate
    # structural model', so CL/F carries a diagonal eta and the three
    # remaining parameters share a 3x3 block. The block off-diagonals are
    # reconstructed from the three reported correlations as
    # cov(i,j) = r_ij * omega_i * omega_j; the resulting matrix is positive
    # definite (eigenvalues 1.4417, 0.2746, 0.0928).
    # =====================================================================
    etalcl ~ 0.067819
    # Table 3 'IIV CL/F = 26.49 (RSE 29.1%)' -> omega^2 = log(1 + 0.2649^2) = 0.067819

    # Lower-triangle order of the block below, with the Table 3 row each entry
    # comes from. rxode2 cannot parse a '#' comment inside the c() of an omega
    # block, so the provenance is recorded here instead of inline.
    #   var(etalvc)                = 0.251462  <- 'IIV Vc/F = 53.47 (RSE 42%)'; log(1 + 0.5347^2)
    #   cov(etalvc,   etalka_ir)   = 0.412802  <- 'Vc/F/Ka IR-Tac Correlation = 75.63 (RSE 16%)'; 0.7563*sqrt(0.251462*1.184742)
    #   var(etalka_ir)             = 1.184742  <- 'IIV Ka IR-Tac = 150.66 (RSE 25.6%)'; log(1 + 1.5066^2)
    #   cov(etalvc,   etalka_lcp)  = 0.135906  <- 'Vc/F/Ka LCP-Tac Correlation = 44.38 (RSE 10%)'; 0.4438*sqrt(0.251462*0.372933)
    #   cov(etalka_ir,etalka_lcp)  = 0.299116  <- 'Ka IR-Tac/Ka LCP-Tac Correlation = 45 (RSE 20.3%)'; 0.45*sqrt(1.184742*0.372933)
    #   var(etalka_lcp)            = 0.372933  <- 'IIV Ka LCP_Tac = 67.23 (RSE 46.5%)'; log(1 + 0.6723^2)
    etalvc + etalka_ir + etalka_lcp ~ c(
      0.251462,
      0.412802, 1.184742,
      0.135906, 0.299116, 0.372933
    )

    # =====================================================================
    # Inter-occasion variability on CL/F and Vc/F. A single magnitude is
    # reported for each (Table 3), and the design supplies exactly two
    # occasions (one 24 h profile per formulation), so occasion 2 is fixed
    # equal to occasion 1 (NONMEM $OMEGA BLOCK(1) SAME).
    # =====================================================================
    etaiov_cl_1 ~ 0.042554
    # Table 3 row 'IOV CL = 20.85 (RSE 23.9%)' -> omega^2 = log(1 + 0.2085^2) = 0.042554
    etaiov_cl_2 ~ fixed(0.042554)

    etaiov_vc_1 ~ 0.297122
    # Table 3 row 'IOV Vc = 58.82 (RSE 28.9%)' -> omega^2 = log(1 + 0.5882^2) = 0.297122
    etaiov_vc_2 ~ fixed(0.297122)

    # =====================================================================
    # Residual error -- proportional. Results section 3.2.1: 'The proportional
    # error model best described the residual error associated with
    # concentrations.' Results section 3.2.2: 'Residual error variability
    # associated with the final model was 13.3%'. Table 3 labels the row
    # 'Combined residual error', but no additive term is reported anywhere in
    # the paper and the prose states the proportional model was selected, so a
    # single proportional term is used.
    # =====================================================================
    propSd <- 0.133
    label("Proportional residual error (fraction)")
    # Table 3 row 'RE. (-) = 13.30 (RSE 8.2%)'; bootstrap median 13.11 (90% CI 11.83-14.14)
  })

  model({
    # ------------------------------------------------------------------
    # 1. Formulation and genotype indicators. The four formulation-by-
    #    genotype groups are mutually exclusive and sum to 1.
    # ------------------------------------------------------------------
    form_ir  <- 1 - FORM_TAC_ENVARSUS
    form_lcp <- FORM_TAC_ENVARSUS

    g_lcp_nonexpr <- form_lcp * (1 - CYP3A5_EXPR)
    g_ir_nonexpr  <- form_ir  * (1 - CYP3A5_EXPR)
    g_lcp_expr    <- form_lcp * CYP3A5_EXPR
    g_ir_expr     <- form_ir  * CYP3A5_EXPR

    # ------------------------------------------------------------------
    # 2. Occasion indicators multiplexing the per-occasion CL/F and Vc/F
    #    etas (occasion 1 = IR-Tac profile, occasion 2 = LCP-Tac profile).
    # ------------------------------------------------------------------
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2
    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2

    # ------------------------------------------------------------------
    # 3. Circadian (cosinor) modulation, Mohammed Ali 2025 equation 3.
    #    `t` MUST be clock time in hours measured from midnight, because the
    #    acrophases are absolute clock times (17:00 for CL/F, 03:08 for the
    #    IR-Tac ka) and not times after dose.
    #
    #    The rhythm multiplies the typical (mesor) value and the
    #    inter-individual random effect is applied to the product, i.e.
    #      P_i(t) = (theta_mesor + theta_AMP * cos(...)) * exp(eta)
    #    This is the standard NONMEM idiom and the only reading under which
    #    the paper's own 1000-subject simulations are reproducible: with the
    #    amplitude added to the INDIVIDUAL value instead, the 150.66 %CV IIV
    #    on the IR-Tac ka would put ka_i below the 1.55/h amplitude in about
    #    40% of subjects and drive the absorption rate negative for part of
    #    every day. See the vignette Assumptions and deviations section.
    #
    #    Both rhythms are 24 h; the CL/F rhythm applies on BOTH formulations
    #    (it is a property of the patient's diurnal metabolic capacity, and
    #    CL/F is a single shared parameter in this joint model), whereas the
    #    ka rhythm is IR-Tac-specific per Results section 3.2.1 and is gated
    #    by form_ir below. Table S1 of the paper confirms the CL/F rhythm
    #    acts on LCP-Tac: its simulated LCP-Tac AUC24 sits 2.5-3.7% above the
    #    time-invariant closed form dose*F/CL, which is the excess a
    #    full-cycle integral of a cosinor-modulated clearance produces.
    # ------------------------------------------------------------------
    #    These two symbols are the signed rhythm INCREMENTS at clock time t,
    #    carrying the units of the parameter they modulate (L/h and 1/h), not
    #    clearances or rate constants in their own right. They are named for
    #    the amplitude they scale rather than for the parameter they enter, per
    #    the parameter-names.md rule that a component must not be nameable as
    #    the total quantity the ODE consumes.
    amp_cl_t    <- amp_cl    * cos(2 * pi * (t - acrophase_cl)    / 24)
    amp_ka_ir_t <- amp_ka_ir * cos(2 * pi * (t - acrophase_ka_ir) / 24)

    # ------------------------------------------------------------------
    # 4. Individual PK parameters.
    # ------------------------------------------------------------------
    cl <- (exp(lcl) + amp_cl_t) * exp(etalcl + iov_cl)
    vc <- exp(lvc + etalvc + iov_vc)
    q  <- exp(lq)
    vp <- exp(lvp)

    ka_ir  <- (exp(lka_ir) + amp_ka_ir_t) * exp(etalka_ir)
    ka_lcp <- exp(lka_lcp + etalka_lcp)
    ka     <- ka_ir * form_ir + ka_lcp * form_lcp

    tlag <- exp(ltlag_ir) * form_ir + exp(ltlag_lcp) * form_lcp

    fdepot <- exp(
      lfdepot_lcp_nonexpr * g_lcp_nonexpr +
        lfdepot_ir_nonexpr  * g_ir_nonexpr +
        lfdepot_lcp_expr    * g_lcp_expr +
        lfdepot_ir_expr     * g_ir_expr
    )

    # ------------------------------------------------------------------
    # 5. Micro-constants.
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ------------------------------------------------------------------
    # 6. Two-compartment ODE system with delayed first-order absorption.
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central -
      k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    alag(depot) <- tlag
    f(depot)    <- fdepot

    # ------------------------------------------------------------------
    # 7. Observation. Doses are in mg and volumes in L, so central / vc is
    #    mg/L; the 1000 factor converts to the ng/mL units in which whole-
    #    blood tacrolimus concentrations are reported.
    # ------------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
