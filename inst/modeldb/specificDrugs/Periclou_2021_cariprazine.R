Periclou_2021_cariprazine <- function() {
  description <- "Population PK model (final model on the updated dataset) for oral cariprazine and its two active metabolites desmethyl-cariprazine (DCAR) and didesmethyl-cariprazine (DDCAR) in adults with schizophrenia or bipolar mania: three-compartment cariprazine with zero-order input into a depot followed by first-order absorption; cariprazine elimination forms DCAR (two-compartment), DCAR elimination forms DDCAR through a single delay transit compartment (two-compartment DDCAR); body weight, race (Black, Asian, Japanese) and sex covariates, and first-dose shifts on cariprazine Vc/Vp1/Q3 and DCAR Vc/Vp"
  reference <- "Periclou A, Phillips L, Ghahramani P, Kapas M, Carrothers T, Khariton T. Population Pharmacokinetics of Cariprazine and its Major Metabolites. Eur J Drug Metab Pharmacokinet. 2021;46(1):53-69. doi:10.1007/s13318-020-00650-4"
  vignette <- "Periclou_2021_cariprazine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects normalised to 79 kg on cariprazine CL/F and Vc/F, DCAR CL/F and Vc/F, and DDCAR CL/F and Vc/F (Supplemental Equations 10, 11, 14, 15, 17, 18). 79 kg is the rounded mean of the model-development cohort (78.9 kg, Supplemental Table 2).",
      source_name = "WTKG"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Proportional shift (1 - 0.160 * Female) on DCAR CL/F only (Supplemental Equation 14). Same orientation as the source 'Female' indicator.",
      source_name = "Female"
    ),
    RACE_BLACK = list(
      description = "Black or African-American race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White/Caucasian or Other; the paper's reference race group is 'white or other')",
      notes = "Proportional shifts on cariprazine CL/F (-0.0907), DCAR CL/F (+0.249), DDCAR CL/F (+0.547) and DDCAR Vc/F (+0.676) (Supplemental Equations 10, 14, 17, 18). Mutually exclusive with RACE_ASIAN and RACE_JAPANESE.",
      source_name = "Black"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator EXCLUDING Japanese (Periclou 2021: Asian patients were mainly from studies conducted in India)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White/Caucasian or Other)",
      notes = "Proportional shifts on cariprazine CL/F (-0.178), DCAR CL/F (-0.0861), DDCAR CL/F (-0.194) and DDCAR Vc/F (-0.240). In the updated analysis race was redefined as White, Black, Asian-Indian, Asian-Japanese and Other, so RACE_ASIAN = 1 selects the non-Japanese Asian group (314 patients, 2 of them in Phase 1) and must be 0 for a Japanese subject, who carries RACE_JAPANESE = 1 instead (same convention as Wade_2015_certolizumab).",
      source_name = "Asian"
    ),
    RACE_JAPANESE = list(
      description = "Japanese race indicator (Study A002-A11 only)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Japanese)",
      notes = "Proportional shifts on cariprazine CL/F (-0.111), DCAR CL/F (-0.145), DDCAR CL/F (-0.156) and DDCAR Vc/F (+0.0888). All 37 Japanese patients came from the single study A002-A11, so the paper cautions that the Japanese effect may be confounded with a study effect. Mutually exclusive with RACE_ASIAN.",
      source_name = "Japanese"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "cariprazine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "cariprazine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cariprazine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "cariprazine", units = "mg", specimen = "plasma", verified = TRUE),
    central_dcar = list(analyte = "desmethyl-cariprazine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_dcar = list(analyte = "desmethyl-cariprazine", units = "mg", specimen = "plasma", verified = TRUE),
    transit1 = list(analyte = "didesmethyl-cariprazine", units = "mg", specimen = "not applicable", verified = TRUE),
    central_ddcar = list(analyte = "didesmethyl-cariprazine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_ddcar = list(analyte = "didesmethyl-cariprazine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 2199L,
    n_studies = 11L,
    age_range = "18-65 years",
    age_mean = "39.2 years (SD 10.8)",
    weight_range = "33.1-155.1 kg",
    weight_mean = "78.9 kg (SD 18.7)",
    weight_median = "77.7 kg",
    sex_female_pct = 33.6,
    race_ethnicity = c(White = 45.6, Black = 34.9, Asian = 14.3, Japanese = 1.7, Other = 3.5),
    disease_state = "Adults with schizophrenia or with manic/mixed episodes of bipolar I disorder",
    dose_range = "Oral cariprazine 0.5-12.5 mg once daily in the model-development dataset (doses >= 15 mg/day excluded)",
    renal_function = "CrCL 31.4-360.5 mL/min; 81.9% normal, 17.1% mild, 0.9% moderate impairment (CrCL not a significant covariate)",
    regions = "Multinational (US, Europe, India, Japan)",
    notes = "Final (updated-dataset) model: 3 phase 1 and 8 phase 2/3 studies plus Japanese study A002-A11; the long-term open-label studies RGH-MD-11 and RGH-MD-17 were used only for external validation. 13,227 cariprazine, 12,462 DCAR and 12,092 DDCAR samples from 2199, 2180 and 2140 patients (Periclou 2021 Section 3.2.2; demographics Supplemental Table 2). Of 908 genotyped patients, 40 were CYP2D6 poor metabolisers; metaboliser status did not affect exposure."
  )

  ini({
    # Cariprazine (Periclou 2021 Table 2; Supplemental Equations 10-13)
    ld1 <- log(2.57)
    label("Duration of zero-order input of the dose into the depot, DUR (h)") # Table 2 'DUR 2.57 (2.35, 2.79)'
    lka <- log(0.352)
    label("First-order absorption rate constant, Ka (1/h)") # Table 2 'Ka 0.352 (0.32, 0.39)'
    lcl <- log(21.5)
    label("Cariprazine apparent clearance CL/F for a 79 kg White patient (L/h)") # Table 2 'CL/F 21.5 (21.1, 21.8)'; Suppl Eq 10
    e_wt_cl <- 0.0946
    label("Power exponent of WT/79 on cariprazine CL/F (unitless)") # Table 2 'Power of WTKG 0.0946'
    e_black_cl <- -0.0907
    label("Proportional shift on cariprazine CL/F for Black race (fraction)") # Table 2 'Proportional shift for race=black -0.0907'
    e_asian_cl <- -0.178
    label("Proportional shift on cariprazine CL/F for Asian (non-Japanese) race (fraction)") # Table 2 'race=Asian -0.178'
    e_japanese_cl <- -0.111
    label("Proportional shift on cariprazine CL/F for Japanese race (fraction)") # Table 2 'race=Japanese -0.111'
    lvc <- log(266)
    label("Cariprazine apparent central volume Vc/F for a 79 kg patient after the second and later doses (L)") # Table 2 'VC/F 266 (241, 288)'; Suppl Eq 11
    e_fd_vc <- fixed(2.84)
    label("Proportional shift on cariprazine Vc/F after the first dose (fraction)") # Table 2 'Proportional shift for first dose 2.84 FIXED'
    e_wt_vc <- 1.66
    label("Power exponent of WT/79 on cariprazine Vc/F (unitless)") # Table 2 'Power of WTKG 1.66'
    lq <- fixed(log(0.431))
    label("Cariprazine apparent first distribution clearance Q3/F (L/h)") # Table 2 'Q3/F 0.431 FIXED'; Suppl Eq 12
    e_fd_q <- fixed(39.4)
    label("Proportional shift on cariprazine Q3/F after the first dose (fraction)") # Table 2 'Proportional shift for first dose 39.4 FIXED'
    lvp <- fixed(log(149))
    label("Cariprazine apparent first peripheral volume VP1/F (L)") # Table 2 'VP1/F 149 FIXED'; Suppl Eq 13
    e_fd_vp <- fixed(2.61)
    label("Proportional shift on cariprazine VP1/F after the first dose (fraction)") # Table 2 'Proportional shift for first dose 2.61 FIXED'
    lq2 <- fixed(log(100))
    label("Cariprazine apparent second distribution clearance Q4/F (L/h)") # Table 2 'Q4/F 100 FIXED'
    lvp2 <- fixed(log(501))
    label("Cariprazine apparent second peripheral volume VP2/F (L)") # Table 2 'VP2/F 501 FIXED'

    # DCAR (Table 2; Supplemental Equations 14-16)
    lcl_dcar <- log(77.3)
    label("DCAR apparent clearance DCL/F for a 79 kg White male (L/h)") # Table 2 'DCL/F 77.3 (75.3, 79.4)'; Suppl Eq 14
    e_wt_cl_dcar <- 0.578
    label("Power exponent of WT/79 on DCAR CL/F (unitless)") # Table 2 'Power of WTKG 0.578'
    e_black_cl_dcar <- 0.249
    label("Proportional shift on DCAR CL/F for Black race (fraction)") # Table 2 'race=black 0.249'
    e_asian_cl_dcar <- -0.0861
    label("Proportional shift on DCAR CL/F for Asian (non-Japanese) race (fraction)") # Table 2 'race=Asian -0.0861'
    e_japanese_cl_dcar <- -0.145
    label("Proportional shift on DCAR CL/F for Japanese race (fraction)") # Table 2 'race=Japanese -0.145'
    e_sexf_cl_dcar <- -0.160
    label("Proportional shift on DCAR CL/F for female sex (fraction)") # Table 2 'sex=female -0.160'
    lvc_dcar <- log(128)
    label("DCAR apparent central volume DVC/F for a 79 kg patient after the second and later doses (L)") # Table 2 'DVC/F 128 (106, 150)'; Suppl Eq 15
    e_fd_vc_dcar <- fixed(1.27)
    label("Proportional shift on DCAR Vc/F after the first dose (fraction)") # Table 2 'Proportional shift for first dose 1.27 FIXED'
    e_wt_vc_dcar <- 1.18
    label("Power exponent of WT/79 on DCAR Vc/F (unitless)") # Table 2 'Power of WTKG 1.18'
    lq_dcar <- log(78.5)
    label("DCAR apparent distribution clearance DQ/F (L/h)") # Table 2 'DQ/F 78.5 (60.9, 105)'
    lvp_dcar <- log(347)
    label("DCAR apparent peripheral volume DVP/F after the second and later doses (L)") # Table 2 'DVP/F 347 (292, 411)'; Suppl Eq 16
    e_fd_vp_dcar <- fixed(0.535)
    label("Proportional shift on DCAR Vp/F after the first dose (fraction)") # Table 2 'Proportional shift for first dose 0.535 FIXED'

    # DDCAR (Table 2; Supplemental Equations 17-18)
    lcl_ddcar <- log(9.24)
    label("DDCAR apparent clearance DDCL/F for a 79 kg White patient (L/h)") # Table 2 'DDCL/F 9.24 (8.93, 9.57)'; Suppl Eq 17
    e_wt_cl_ddcar <- 0.427
    label("Power exponent of WT/79 on DDCAR CL/F (unitless)") # Table 2 'Power of WTKG 0.427'
    e_black_cl_ddcar <- 0.547
    label("Proportional shift on DDCAR CL/F for Black race (fraction)") # Table 2 'race=black 0.547'
    e_asian_cl_ddcar <- -0.194
    label("Proportional shift on DDCAR CL/F for Asian (non-Japanese) race (fraction)") # Table 2 'race=Asian -0.194'
    e_japanese_cl_ddcar <- -0.156
    label("Proportional shift on DDCAR CL/F for Japanese race (fraction)") # Table 2 'race=Japanese -0.156'
    lvc_ddcar <- log(1310)
    label("DDCAR apparent central volume DDVC/F for a 79 kg White patient (L)") # Table 2 'DDVC/F 1310 (1260, 1360)' (unit printed as 'l/h', a typo for L); Suppl Eq 18
    e_wt_vc_ddcar <- 0.881
    label("Power exponent of WT/79 on DDCAR Vc/F (unitless)") # Table 2 'Power of WTKG 0.881'
    e_black_vc_ddcar <- 0.676
    label("Proportional shift on DDCAR Vc/F for Black race (fraction)") # Table 2 'race=black 0.676'
    e_asian_vc_ddcar <- -0.240
    label("Proportional shift on DDCAR Vc/F for Asian (non-Japanese) race (fraction)") # Table 2 'race=Asian -0.240'
    e_japanese_vc_ddcar <- 0.0888
    label("Proportional shift on DDCAR Vc/F for Japanese race (fraction)") # Table 2 'race=Japanese 0.0888'
    lq_ddcar <- fixed(log(0.386))
    label("DDCAR apparent distribution clearance DDQ/F (L/h)") # Table 2 'DDQ/F 0.386 FIXED'
    lvp_ddcar <- fixed(log(258))
    label("DDCAR apparent peripheral volume DDVP/F (L)") # Table 2 'DDVP/F 258 FIXED'
    lktr <- fixed(log(0.0269))
    label("Rate constant of the transit compartment delaying DDCAR formation, DDKtr (1/h)") # Table 2 'DD Ktr 0.0269 FIXED'

    # IIV: Table 2 reports %CV; converted to log-normal variance omega^2 = log(1 + CV^2).
    # Off-diagonal elements are not reported, so the OMEGA is taken as diagonal.
    etalka ~ 0.87230 # Table 2 'Ka IIV 118% CV' (Phase 1 patients only, footnote b; shrinkage 76.5%)
    etalcl ~ 0.09982 # Table 2 'CL/F IIV 32.4% CV'
    etalvc ~ 0.78303 # Table 2 'VC/F IIV 109% CV' (Phase 1 patients only, footnote b)
    etalcl_dcar ~ 0.16532 # Table 2 'DCL/F IIV 42.4% CV'
    etalvc_dcar ~ 0.84264 # Table 2 'DVC/F IIV 115% CV' (Phase 1 patients only, footnote b)
    etalcl_ddcar ~ 0.28478 # Table 2 'DDCL/F IIV 57.4% CV'
    etalvc_ddcar ~ 0.47330 # Table 2 'DDVC/F IIV 77.8% CV'

    # Residual error: Table 2 footnote f states residual variability was estimated
    # separately for Phase 1 and Phase 2/3 studies, but neither the error model nor
    # its magnitudes are reported anywhere in the paper or supplement.
    propSd <- fixed(0)
    label("Proportional residual error on cariprazine (fraction; ZERO - not reported in source)") # not reported: Table 2 footnote f only
    propSd_dcar <- fixed(0)
    label("Proportional residual error on DCAR (fraction; ZERO - not reported in source)") # not reported: Table 2 footnote f only
    propSd_ddcar <- fixed(0)
    label("Proportional residual error on DDCAR (fraction; ZERO - not reported in source)") # not reported: Table 2 footnote f only
  })

  model({
    # First-dose indicator FD (Suppl Equation Set 3): 1 for records after the first
    # dose and before the second, 0 from the second dose onward. dosenum() is 0
    # before and 1 after the first dose record, so the indicator needs no data column.
    fd <- 1 * (dosenum() <= 1)
    wtn <- WT / 79

    # Cariprazine (Suppl Eq 10-13)
    d1 <- exp(ld1)
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) *
      (1 + e_black_cl * RACE_BLACK) *
      (1 + e_asian_cl * RACE_ASIAN) *
      (1 + e_japanese_cl * RACE_JAPANESE) *
      wtn^e_wt_cl
    vc <- exp(lvc + etalvc) * wtn^e_wt_vc * (1 + e_fd_vc * fd)
    q <- exp(lq) * (1 + e_fd_q * fd)
    vp <- exp(lvp) * (1 + e_fd_vp * fd)
    q2 <- exp(lq2)
    vp2 <- exp(lvp2)

    # DCAR (Suppl Eq 14-16)
    cl_dcar <- exp(lcl_dcar + etalcl_dcar) *
      (1 + e_black_cl_dcar * RACE_BLACK) *
      (1 + e_asian_cl_dcar * RACE_ASIAN) *
      (1 + e_japanese_cl_dcar * RACE_JAPANESE) *
      (1 + e_sexf_cl_dcar * SEXF) *
      wtn^e_wt_cl_dcar
    vc_dcar <- exp(lvc_dcar + etalvc_dcar) * wtn^e_wt_vc_dcar * (1 + e_fd_vc_dcar * fd)
    q_dcar <- exp(lq_dcar)
    vp_dcar <- exp(lvp_dcar) * (1 + e_fd_vp_dcar * fd)

    # DDCAR (Suppl Eq 17-18)
    cl_ddcar <- exp(lcl_ddcar + etalcl_ddcar) *
      (1 + e_black_cl_ddcar * RACE_BLACK) *
      (1 + e_asian_cl_ddcar * RACE_ASIAN) *
      (1 + e_japanese_cl_ddcar * RACE_JAPANESE) *
      wtn^e_wt_cl_ddcar
    vc_ddcar <- exp(lvc_ddcar + etalvc_ddcar) *
      (1 + e_black_vc_ddcar * RACE_BLACK) *
      (1 + e_asian_vc_ddcar * RACE_ASIAN) *
      (1 + e_japanese_vc_ddcar * RACE_JAPANESE) *
      wtn^e_wt_vc_ddcar
    q_ddcar <- exp(lq_ddcar)
    vp_ddcar <- exp(lvp_ddcar)
    ktr <- exp(lktr)

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2
    kel_dcar <- cl_dcar / vc_dcar
    k12_dcar <- q_dcar / vc_dcar
    k21_dcar <- q_dcar / vp_dcar
    kel_ddcar <- cl_ddcar / vc_ddcar
    k12_ddcar <- q_ddcar / vc_ddcar
    k21_ddcar <- q_ddcar / vp_ddcar

    # ODEs (Supplemental Equation Set 1, Eqs 1-9). All cariprazine eliminated is
    # converted to DCAR and all DCAR eliminated to DDCAR (Section 2.4.1), so the
    # DCAR and DDCAR parameters are apparent values. Amounts are carried in mg of
    # dose throughout, with no molecular-weight correction.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    d/dt(central_dcar) <- kel * central - (kel_dcar + k12_dcar) * central_dcar + k21_dcar * peripheral1_dcar
    d/dt(peripheral1_dcar) <- k12_dcar * central_dcar - k21_dcar * peripheral1_dcar
    d/dt(transit1) <- kel_dcar * central_dcar - ktr * transit1
    d/dt(central_ddcar) <- ktr * transit1 - (kel_ddcar + k12_ddcar) * central_ddcar + k21_ddcar * peripheral1_ddcar
    d/dt(peripheral1_ddcar) <- k12_ddcar * central_ddcar - k21_ddcar * peripheral1_ddcar

    # Zero-order input of the dose into the depot over DUR (Eq 1, Figure 1).
    # Dose records must carry rate = -2 for rxode2 to honour dur(depot).
    dur(depot) <- d1

    # mg/L -> ng/mL
    Cc <- 1000 * central / vc
    Cc_dcar <- 1000 * central_dcar / vc_dcar
    Cc_ddcar <- 1000 * central_ddcar / vc_ddcar

    Cc ~ prop(propSd)
    Cc_dcar ~ prop(propSd_dcar)
    Cc_ddcar ~ prop(propSd_ddcar)
  })
}
