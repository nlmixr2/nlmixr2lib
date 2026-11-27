Ismaeil_2026_cariprazine <- function() {
  description <- "Pediatric population PK model for oral cariprazine and its two active metabolites desmethyl-cariprazine (DCAR) and didesmethyl-cariprazine (DDCAR) in children (10-12 years), adolescents (13-17 years) and young adults with schizophrenia or bipolar I disorder: the Periclou 2021 adult structure (three-compartment cariprazine with zero-order input into a depot followed by first-order absorption; cariprazine elimination forms two-compartment DCAR, DCAR elimination forms two-compartment DDCAR through one delay transit compartment) with all structural, race, sex and first-dose parameters fixed to the adult values and the body-weight exponents on CL/F and Vc/F of all three analytes, the IIV and the full-profile-study residual errors re-estimated"
  reference <- "Ismaeil S, Siovitz L, Goyal RK, Eckert D, Stodtmann S, Xiong H. Extrapolating Efficacy of Cariprazine to Pediatric Schizophrenia and Bipolar Mania Through Population Pharmacokinetic Analysis. J Clin Pharmacol. 2026;66(7):e70243. doi:10.1002/jcph.70243"
  vignette <- "Ismaeil_2026_cariprazine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects normalised to 79 kg (the adult reference weight, Ismaeil 2026 Equation 1) on CL/F and Vc/F of cariprazine, DCAR and DDCAR; the six exponents were re-estimated on the pediatric data (Tables 3-5).",
      source_name = "BW"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Proportional shift (1 - 0.160 * Female) on DCAR CL/F only, fixed to the adult value (Table 4).",
      source_name = "sex = female"
    ),
    RACE_BLACK = list(
      description = "Black or African-American race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White or Other)",
      notes = "Proportional shifts fixed to the adult values on cariprazine CL/F (-0.0907), DCAR CL/F (+0.249), DDCAR CL/F (+0.547) and DDCAR Vc/F (+0.676) (Tables 3-5). Mutually exclusive with RACE_ASIAN and RACE_JAPANESE.",
      source_name = "race = Black"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator EXCLUDING Japanese",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White or Other)",
      notes = "Proportional shifts fixed to the adult values on cariprazine CL/F (-0.178), DCAR CL/F (-0.0861), DDCAR CL/F (-0.194) and DDCAR Vc/F (-0.240) (Tables 3-5). The adult model's five race levels (White, Black, Asian, Japanese, Other) are kept, so RACE_ASIAN = 1 is the non-Japanese Asian group and must be 0 for a Japanese subject.",
      source_name = "race = Asian"
    ),
    RACE_JAPANESE = list(
      description = "Japanese race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Japanese)",
      notes = "Proportional shifts fixed to the adult values on cariprazine CL/F (-0.111), DCAR CL/F (-0.145), DDCAR CL/F (-0.156) and DDCAR Vc/F (+0.0888) (Tables 3-5). No Japanese patient was in the pediatric dataset. Mutually exclusive with RACE_ASIAN.",
      source_name = "race = Japanese"
    ),
    STUDY_CARIPRAZINE_FULLPROFILE = list(
      description = "1 = observation from a cariprazine study with full-profile PK sampling; 0 = sparse-sampling study",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (sparse-sampling study)",
      notes = "Switches the proportional residual error between the re-estimated 'full profile studies' value and the fixed adult value (Tables 3-5). Both pediatric phase 1 studies (all observations of the model-development dataset) are full-profile studies, so set 1 to simulate the pediatric data.",
      source_name = "full profile studies"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Table 5 lists 'Power of Vc/F for age (fixed) 0' for DDCAR: the power is fixed to zero, so age has no effect and the reference age is not given."
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
    n_subjects = 101L,
    n_studies = 2L,
    age_range = "10-40 years (25 children 10-12 years, 60 adolescents 13-17 years, 16 adults 18-40 years)",
    age_mean = "16 years (SD 6)",
    weight_range = "26.9-95.4 kg",
    weight_mean = "59.6 kg (SD 14.7)",
    sex_female_pct = 39.6,
    race_ethnicity = c(White = 58, Black = 41, Asian = 1, Other = 0),
    disease_state = "Pediatric patients with schizophrenia (13-17 years) or bipolar I disorder (10-17 years), and adults with schizophrenia (18-40 years)",
    dose_range = "Oral cariprazine 1.5, 3, 4.5 or 6 mg once daily after titration, for 42 days (Study 1) or 28 days (Study 2)",
    regions = "Not stated (Study 2 has a EudraCT number)",
    notes = "Two phase 1 studies: Study 1 (N = 49; schizophrenia 13-17 years or bipolar I disorder 10-17 years; 1.5, 3 or 4.5 mg for 42 days) and Study 2 (N = 52; schizophrenia, 36 adolescents and 16 adults; 1.5, 3 or 6 mg for 28 days). 4237 concentrations (cariprazine 1570, DCAR 1269, DDCAR 1398). Nine patients weighed < 40 kg (Ismaeil 2026 Table 2)."
  )

  ini({
    # Cariprazine (Table 3). Structural and covariate values are fixed to the
    # Periclou 2021 adult model; the weight exponents are re-estimated.
    ld1 <- fixed(log(2.57))
    label("Duration of zero-order input of the dose into the depot, DUR (h)") # Table 3 'DUR/F, h (fixed) 2.57'
    lka <- fixed(log(0.352))
    label("First-order absorption rate constant, Ka (1/h)") # Table 3 'ka, 1/h (fixed) 0.352'
    lcl <- fixed(log(21.5))
    label("Cariprazine apparent clearance CL/F for a 79 kg White patient (L/h)") # Table 3 'CL/F, L/h (fixed) 21.5'
    e_wt_cl <- 0.0977
    label("Power exponent of WT/79 on cariprazine CL/F (unitless)") # Table 3 'Power of CL/F for weight 0.0977' (RSE 73.4%)
    e_black_cl <- fixed(-0.0907)
    label("Proportional shift on cariprazine CL/F for Black race (fraction)") # Table 3 'CL/F prop shift for race = Black (fixed) -0.0907'
    e_asian_cl <- fixed(-0.178)
    label("Proportional shift on cariprazine CL/F for Asian (non-Japanese) race (fraction)") # Table 3 'race = Asian (fixed) -0.178'
    e_japanese_cl <- fixed(-0.111)
    label("Proportional shift on cariprazine CL/F for Japanese race (fraction)") # Table 3 'race = Japanese (fixed) -0.111'
    lvc <- fixed(log(266))
    label("Cariprazine apparent central volume Vc/F for a 79 kg patient after the second and later doses (L)") # Table 3 'Vc/F, L (fixed) 266'
    e_fd_vc <- fixed(2.84)
    label("Proportional shift on cariprazine Vc/F after the first dose (fraction)") # Table 3 'Vc/F prop shift for 1st dose (fixed) 2.84'
    e_wt_vc <- 0.998
    label("Power exponent of WT/79 on cariprazine Vc/F (unitless)") # Table 3 'Power of Vc/F for weight 0.998' (RSE 28.6%)
    lq <- fixed(log(0.431))
    label("Cariprazine apparent first distribution clearance Q3/F (L/h)") # Table 3 'Q3/F, L (fixed) 0.431' (unit printed as L; L/h per Periclou 2021)
    e_fd_q <- fixed(39.4)
    label("Proportional shift on cariprazine Q3/F after the first dose (fraction)") # Table 3 'Q3/F prop shift for 1st dose (fixed) 39.4'
    lvp <- fixed(log(149))
    label("Cariprazine apparent first peripheral volume Vp1/F (L)") # Table 3 'Vp1/F, L (fixed) 149'
    e_fd_vp <- fixed(2.61)
    label("Proportional shift on cariprazine Vp1/F after the first dose (fraction)") # Table 3 'Vp1/F prop. shift for 1st dose (fixed) 2.61'
    lq2 <- fixed(log(100))
    label("Cariprazine apparent second distribution clearance Q4/F (L/h)") # Table 3 'Q4/F, L/h (fixed) 100'
    lvp2 <- fixed(log(501))
    label("Cariprazine apparent second peripheral volume Vp2/F (L)") # Table 3 'Vp2/F, L (fixed) 501'

    # DCAR (Table 4)
    lcl_dcar <- fixed(log(77.3))
    label("DCAR apparent clearance CL/F for a 79 kg White male (L/h)") # Table 4 'CL/F, L/h (fixed) 77.3'
    e_wt_cl_dcar <- 0.0697
    label("Power exponent of WT/79 on DCAR CL/F (unitless)") # Table 4 'Power of CL/F for weight (in kg) 0.0697' (RSE 121%)
    e_black_cl_dcar <- fixed(0.249)
    label("Proportional shift on DCAR CL/F for Black race (fraction)") # Table 4 'race = Black (fixed) 0.249'
    e_asian_cl_dcar <- fixed(-0.0861)
    label("Proportional shift on DCAR CL/F for Asian (non-Japanese) race (fraction)") # Table 4 'race = Asian (fixed) -0.0861'
    e_japanese_cl_dcar <- fixed(-0.145)
    label("Proportional shift on DCAR CL/F for Japanese race (fraction)") # Table 4 'race = Japanese (fixed) -0.145'
    e_sexf_cl_dcar <- fixed(-0.160)
    label("Proportional shift on DCAR CL/F for female sex (fraction)") # Table 4 'sex = female (fixed) -0.160'
    lvc_dcar <- fixed(log(128))
    label("DCAR apparent central volume Vc/F for a 79 kg patient after the second and later doses (L)") # Table 4 'Vc/F, L (fixed) 128'
    e_fd_vc_dcar <- fixed(1.27)
    label("Proportional shift on DCAR Vc/F after the first dose (fraction)") # Table 4 'Vc/F proportional shift for 1st dose (fixed) 1.27'
    e_wt_vc_dcar <- 0.284
    label("Power exponent of WT/79 on DCAR Vc/F (unitless)") # Table 4 'Power of Vc/F for weight (in kg) 0.284' (RSE 174%)
    lq_dcar <- fixed(log(78.5))
    label("DCAR apparent distribution clearance Q3/F (L/h)") # Table 4 'Q3/F, L/h (fixed) 78.5'
    lvp_dcar <- fixed(log(347))
    label("DCAR apparent peripheral volume Vp/F after the second and later doses (L)") # Table 4 'Vp/F, L (fixed) 347'
    e_fd_vp_dcar <- fixed(0.535)
    label("Proportional shift on DCAR Vp/F after the first dose (fraction)") # Table 4 'Vp/F proportional shift for 1st dose (fixed) 0.535'

    # DDCAR (Table 5)
    lcl_ddcar <- fixed(log(9.24))
    label("DDCAR apparent clearance CL/F for a 79 kg White patient (L/h)") # Table 5 'CL/F, L/h (fixed) 9.24'
    e_wt_cl_ddcar <- 0.511
    label("Power exponent of WT/79 on DDCAR CL/F (unitless)") # Table 5 'Power of CL/F for weight (in kg) 0.511' (RSE 20.1%)
    e_black_cl_ddcar <- fixed(0.547)
    label("Proportional shift on DDCAR CL/F for Black race (fraction)") # Table 5 'CL/F proportional shift for race = Black (fixed) 0.547'
    e_asian_cl_ddcar <- fixed(-0.194)
    label("Proportional shift on DDCAR CL/F for Asian (non-Japanese) race (fraction)") # Table 5 'race = Asian (fixed) -0.194'
    e_japanese_cl_ddcar <- fixed(-0.156)
    label("Proportional shift on DDCAR CL/F for Japanese race (fraction)") # Table 5 'race = Japanese (fixed) -0.156'
    lvc_ddcar <- fixed(log(1310))
    label("DDCAR apparent central volume Vc/F for a 79 kg White patient (L)") # Table 5 'Vc/F, L (fixed) 1310'
    e_wt_vc_ddcar <- 0.244
    label("Power exponent of WT/79 on DDCAR Vc/F (unitless)") # Table 5 'Power of Vc/F for weight (in kg) 0.244' (RSE 54.7%)
    e_black_vc_ddcar <- fixed(0.676)
    label("Proportional shift on DDCAR Vc/F for Black race (fraction)") # Table 5 'Vc/F proportional shift for race = Black (fixed) 0.676'
    e_asian_vc_ddcar <- fixed(-0.240)
    label("Proportional shift on DDCAR Vc/F for Asian (non-Japanese) race (fraction)") # Table 5 'race = Asian (fixed) -0.240'
    e_japanese_vc_ddcar <- fixed(0.0888)
    label("Proportional shift on DDCAR Vc/F for Japanese race (fraction)") # Table 5 'race = Japanese (fixed) 0.0888'
    lq_ddcar <- fixed(log(0.386))
    label("DDCAR apparent distribution clearance Q3/F (L/h)") # Table 5 'Q3/F, L/h (fixed) 0.386'
    lvp_ddcar <- fixed(log(258))
    label("DDCAR apparent peripheral volume Vp/F (L)") # Table 5 'Vp/F, L (fixed) 258'
    lktr <- fixed(log(0.0269))
    label("Rate constant of the transit compartment delaying DDCAR formation, DDKTR (1/h)") # Table 5 'DDKTR/F, 1/h (fixed) 0.0269'

    # IIV: Tables 3-5 report the variance omega^2 (the %CV column equals
    # sqrt(exp(omega^2) - 1) x 100). Off-diagonal elements are not reported.
    etalka ~ 0.906 # Table 3 IIV 'On k a 0.906' (121% CV, shrinkage 44.5%)
    etalcl ~ 0.0608 # Table 3 IIV 'On CL/F 0.0608' (25.0% CV)
    etalvc ~ 1.47 # Table 3 IIV 'On Vc/F 1.47' (184% CV)
    etalcl_dcar ~ 0.105 # Table 4 IIV 'On CL/F 0.105' (33.2% CV)
    etalvc_dcar ~ 1.56 # Table 4 IIV 'On Vc/F 1.56' (194% CV)
    etalcl_ddcar ~ 0.167 # Table 5 IIV 'On CL/F 0.167' (42.6% CV)
    etalvc_ddcar ~ 0.295 # Table 5 IIV 'On Vc/F 0.295' (58.6% CV)

    # Residual error: proportional, one magnitude for full-profile studies
    # (re-estimated) and one fixed to the adult value for the other studies.
    # The tabulated values are read as SDs (see the vignette).
    propSd_fullprofile <- 0.107
    label("Proportional residual SD of cariprazine, full-profile studies (fraction)") # Table 3 'Proportional error full profile studies 0.107'
    propSd_sparse <- fixed(0.154)
    label("Proportional residual SD of cariprazine, other studies (fraction)") # Table 3 'Proportional error (fixed) 0.154'
    propSd_fullprofile_dcar <- 0.0688
    label("Proportional residual SD of DCAR, full-profile studies (fraction)") # Table 4 'Proportional error full profile studies 0.0688'
    propSd_sparse_dcar <- fixed(0.141)
    label("Proportional residual SD of DCAR, other studies (fraction)") # Table 4 'Proportional error (fixed) 0.141'
    propSd_fullprofile_ddcar <- 0.0470
    label("Proportional residual SD of DDCAR, full-profile studies (fraction)") # Table 5 'Proportional error full profile studies 0.0470'
    propSd_sparse_ddcar <- fixed(0.0875)
    label("Proportional residual SD of DDCAR, other studies (fraction)") # Table 5 'Proportional error DCAR (fixed) 0.0875' (row label reads DCAR in a DDCAR table)
  })

  model({
    # First-dose indicator: 1 after the first dose and before the second, 0
    # from the second dose onward (Periclou 2021 Supplemental Equation Set 3).
    # dosenum() counts the doses given so far, so no data column is needed.
    fd <- 1 * (dosenum() <= 1)
    wtn <- WT / 79

    # Cariprazine (Equation 1 for the weight terms)
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

    # DCAR
    cl_dcar <- exp(lcl_dcar + etalcl_dcar) *
      (1 + e_black_cl_dcar * RACE_BLACK) *
      (1 + e_asian_cl_dcar * RACE_ASIAN) *
      (1 + e_japanese_cl_dcar * RACE_JAPANESE) *
      (1 + e_sexf_cl_dcar * SEXF) *
      wtn^e_wt_cl_dcar
    vc_dcar <- exp(lvc_dcar + etalvc_dcar) * wtn^e_wt_vc_dcar * (1 + e_fd_vc_dcar * fd)
    q_dcar <- exp(lq_dcar)
    vp_dcar <- exp(lvp_dcar) * (1 + e_fd_vp_dcar * fd)

    # DDCAR
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

    # All cariprazine eliminated forms DCAR, and all DCAR eliminated forms DDCAR
    # through one transit compartment, amount for amount in mg of dose (the
    # metabolite parameters are apparent values).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    d/dt(central_dcar) <- kel * central - (kel_dcar + k12_dcar) * central_dcar + k21_dcar * peripheral1_dcar
    d/dt(peripheral1_dcar) <- k12_dcar * central_dcar - k21_dcar * peripheral1_dcar
    d/dt(transit1) <- kel_dcar * central_dcar - ktr * transit1
    d/dt(central_ddcar) <- ktr * transit1 - (kel_ddcar + k12_ddcar) * central_ddcar + k21_ddcar * peripheral1_ddcar
    d/dt(peripheral1_ddcar) <- k12_ddcar * central_ddcar - k21_ddcar * peripheral1_ddcar

    # Zero-order input of the dose into the depot over DUR. Dose records must
    # carry rate = -2 for rxode2 to honour dur(depot).
    dur(depot) <- d1

    # Residual-error magnitude by study sampling design
    propSd_car <- propSd_fullprofile * STUDY_CARIPRAZINE_FULLPROFILE +
      propSd_sparse * (1 - STUDY_CARIPRAZINE_FULLPROFILE)
    propSd_car_dcar <- propSd_fullprofile_dcar * STUDY_CARIPRAZINE_FULLPROFILE +
      propSd_sparse_dcar * (1 - STUDY_CARIPRAZINE_FULLPROFILE)
    propSd_car_ddcar <- propSd_fullprofile_ddcar * STUDY_CARIPRAZINE_FULLPROFILE +
      propSd_sparse_ddcar * (1 - STUDY_CARIPRAZINE_FULLPROFILE)

    # mg/L -> ng/mL
    Cc <- 1000 * central / vc
    Cc_dcar <- 1000 * central_dcar / vc_dcar
    Cc_ddcar <- 1000 * central_ddcar / vc_ddcar

    Cc ~ prop(propSd_car)
    Cc_dcar ~ prop(propSd_car_dcar)
    Cc_ddcar ~ prop(propSd_car_ddcar)
  })
}
