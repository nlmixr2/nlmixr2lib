# Population pharmacokinetic model for oral clindamycin given with
# fosmidomycin and artesunate to Gabonese children and adults with
# uncomplicated Plasmodium falciparum malaria
# (Pfaffendorf 2026, Malaria Journal 25:152;
# doi:10.1186/s12936-026-05872-6).

Pfaffendorf_2026_clindamycin <- function() {
  description <- paste(
    "Population PK model for oral clindamycin co-administered with",
    "fosmidomycin and artesunate in Gabonese children and adults with",
    "uncomplicated Plasmodium falciparum malaria (Pfaffendorf 2026).",
    "One-compartment disposition with linear elimination, first-order",
    "absorption and an absorption lag time. Allometric body-weight",
    "scaling is applied with fixed exponents 0.75 on apparent clearance",
    "and 1 on apparent volume, centered at the cohort median weight of",
    "29.05 kg; no other covariate was retained. Correlated log-normal",
    "IIV is carried on the absorption rate constant and the lag time,",
    "and inter-occasion variability on apparent clearance across the",
    "six dosing occasions. Residual error is combined proportional plus",
    "additive.",
    sep = " "
  )
  reference <- paste(
    "Pfaffendorf C, Dejon-Agobe JC, Edoa JR, Maiga-Ascofare O,",
    "Ahenkan E, Adegnika AA, Ramharter M, Wicha SG, Mischlinger J (2026).",
    "Population pharmacokinetics of fosmidomycin and clindamycin in",
    "combination with artesunate for uncomplicated Plasmodium falciparum",
    "malaria in Gabonese children and adults.",
    "Malaria Journal 25:152. doi:10.1186/s12936-026-05872-6.",
    "Parameter values from Table 2 (clindamycin block); model structure",
    "from Additional file 1 section S3 (final NONMEM control stream).",
    sep = " "
  )
  vignette <- "Pfaffendorf_2026_fosmidomycin_clindamycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject at admission. Pfaffendorf 2026 Table 1",
        "reports mean (SD) 37.5 (20.2) kg and median [min, max]",
        "29.1 [12.0, 86.0] kg across the 40-patient cohort. Applied as",
        "allometric scaling with exponents fixed at the canonical 0.75",
        "on CL/F and 1 on V/F, centered at 29.05 kg:",
        "CL_i = CL_TV * (WGT/29.05)^0.75 and V_i = V_TV * (WGT/29.05)^1.",
        "Body weight was the only covariate retained by the stepwise",
        "covariate model for clindamycin (Results, Clindamycin model);",
        "in particular the albumin effect on volume reported by Smith",
        "et al. was tested and not found significant here (Discussion).",
        "The centering constant 29.05 kg is the study-population median",
        "weight used for the reference-corrected VPC (Methods, VPC",
        "paragraph) and appears verbatim in the final control stream",
        "(Additional file 1 S3). Note it differs slightly from the",
        "Table 1 median of 29.1 kg, which is the same quantity rounded",
        "to three significant figures.",
        sep = " "
      ),
      source_name        = "WGT"
    ),
    OCC = list(
      description        = "Dosing-occasion index for inter-occasion variability on CL/F",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Values 1-6 identify which of the six 12-hourly doses a record",
        "belongs to. The final control stream (Additional file 1 S3)",
        "decomposes OCC into six binary flags and multiplexes six IOV",
        "etas onto log-CL:",
        "IOV_CL = FLAG1*ETA(6) + ... + FLAG6*ETA(11), declared as",
        "$OMEGA BLOCK(1) 0.0798 followed by five BLOCK(1) SAME records,",
        "so all six occasions share one estimated variance. rxode2 has",
        "no `| occ` level, so the occasions are expanded here into six",
        "indicator-multiplexed etas; occasion 1 carries the estimated",
        "variance and occasions 2-6 are fixed to the same value",
        "(the registered idiom -- see Chu_2024_allopurinol.R,",
        "Blackman_2026_methotrexate.R, Jonsson_2011_ethambutol.R).",
        "Six occasions are encoded because the trial gave exactly six",
        "doses (10 mg/kg every 12 h for three days). For single-dose or",
        "single-occasion simulations pass OCC = 1 so the first IOV eta",
        "applies.",
        sep = " "
      ),
      source_name        = "OCC"
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "clindamycin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "clindamycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species         = "human",
    n_subjects      = 40L,
    n_studies       = 1L,
    age_range       = "3.50-57.2 years (Table 1)",
    age_median      = "10.8 years (Table 1)",
    weight_range    = "12.0-86.0 kg (Table 1)",
    weight_median   = "29.1 kg (Table 1)",
    sex_female_pct  = 42.5,
    disease_state   = paste(
      "Microscopically confirmed uncomplicated Plasmodium falciparum",
      "mono-infection with admission parasitaemia 1000-100,000 asexual",
      "parasites per microlitre and a history of fever within the",
      "previous 24 h. Patients with severe malaria (WHO definition),",
      "haemoglobin below 8 g/dL, significant medical disorders, or",
      "antimalarial treatment in the previous 6 weeks were excluded.",
      "Admission laboratory values (Table 1): GFR mean (SD)",
      "126 (34.9) mL/min/1.73m2, haematocrit 32.7 (3.96) %,",
      "haemoglobin 10.9 (1.39) g/dL, albumin 36.0 (5.49) g/L."
    ),
    dose_range      = paste(
      "Oral clindamycin 10 mg/kg every 12 h for three days (six doses",
      "total), co-administered with fosmidomycin 30 mg/kg and",
      "artesunate 2 mg/kg on the same schedule. Doses were rounded to",
      "the closest match achievable with the available 150 mg, 300 mg",
      "and 600 mg products (Ratiopharm, Ulm, Germany). The reference",
      "-corrected VPC used a 450 mg reference dose at 29.05 kg. The",
      "paper's dose-simulation recommendation (Results, Clindamycin",
      "dosing simulation, and Additional file 1 S8) is 12 mg/kg below",
      "35 kg and 10 mg/kg above, implemented as the band scheme",
      "150 mg below 18.7 kg, 300 mg for 18.7-31.2 kg, 450 mg for",
      "31.2-52.5 kg, 600 mg for 52.5-67.5 kg, 750 mg for 67.5-82.5 kg",
      "and 900 mg for 82.5-97.5 kg."
    ),
    regions         = "Gabon (Centre de Recherches Medicales de Lambarene)",
    trial_registration = "PACTR202008909968293 (pactr.samrc.ac.za)",
    notes           = paste(
      "Recruitment was stratified into three age groups: 20 patients",
      "aged 6 months-10 years, 10 aged 11-17 years, and 10 aged",
      "18-65 years. A total of 274 clindamycin plasma samples entered",
      "the PK analysis. Sampling was weight-banded to limit blood draws",
      "in young children: patients above 35 kg at 0.25, 0.5, 0.75, 1.5,",
      "3, 5, 8 h (day 0), 24 h, 48 h and day 7; patients 20-35 kg at 0,",
      "1, 5, 8 h (day 0), 24 h, 48 h and day 7; patients below 20 kg at",
      "0, 1, 8 h (day 0), 24 h, 48 h and day 7. Samples taken at the",
      "same time as a dose were drawn pre-dose. Bioanalysis by LC-MS/MS",
      "after acetonitrile protein precipitation, evaporation and",
      "reconstitution in 50/50 methanol / 20 mM ammonium formate;",
      "calibration range 0.005-0.5 mg/L with dilution of",
      "above-range samples. Data below the limit of quantification were",
      "handled by the M3 method. Estimation used NONMEM 7.5 with the",
      "Laplacian INTERACTION method; parameter uncertainty by sampling",
      "importance resampling (SIR), which is the source of the 95%",
      "confidence intervals in Table 2.",
      "This arm of the trial ran only at the Lambarene site.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters. Typical values are Table 2 (clindamycin
    # block, 'Typical value' column) and are identical to the $THETA
    # records of the final control stream in Additional file 1 S3, so
    # these are final estimates rather than initial values. They apply
    # to a patient of 29.05 kg. Values are reported on the linear
    # scale; log() is applied here for the nlmixr2 internal scale.
    # ------------------------------------------------------------------
    lcl   <- log(8.02)  ; label("Apparent clearance CL/F at 29.05 kg (L/h)")                 # Table 2 CL/F = 8.02 [7.2, 9.0]; S3 $THETA (0.001, 8.02) ;1_CL
    lvc   <- log(28.4)  ; label("Apparent central volume V/F at 29.05 kg (L)")               # Table 2 V/F  = 28.4 [25.4, 32.6]; S3 $THETA (0.001, 28.4) ;2_Vc
    lka   <- log(2.2)   ; label("Absorption rate constant ka (1/h)")                         # Table 2 ka   = 2.2  [1.3, 3.8]; S3 $THETA (0.001, 2.2) ;3_KA
    ltlag <- log(0.227) ; label("Absorption lag time (h)")                                   # Table 2 tlag = 0.227 [0.20, 0.24]; S3 $THETA (0, 0.227) ;5_LAG-Time

    # Relative bioavailability is a structural anchor fixed to 1. Table 2
    # tabulates no F row for clindamycin, but the final control stream
    # carries $THETA (1) FIX ;4_F1 with its eta also fixed to zero, so
    # F is identically 1 in this model.
    lfdepot <- fixed(log(1)) ; label("Relative bioavailability F (unitless)")                # S3 $THETA (1) FIX ;4_F1

    # ------------------------------------------------------------------
    # Allometric exponents. The paper tested 'allometric scaling with
    # fixed exponents' (Methods, Covariate model building) and the
    # control stream hardcodes 0.75 and 1 rather than estimating them,
    # so neither appears in Table 2 and both are fixed here.
    # ------------------------------------------------------------------
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent on CL/F (unitless)")                 # Methods 'Covariate model building'; S3 $PK (WGT/29.05)**0.75
    e_wt_vc <- fixed(1)    ; label("Allometric exponent on V/F (unitless)")                  # Methods 'Covariate model building'; S3 $PK (WGT/29.05)**1

    # ------------------------------------------------------------------
    # Inter-individual variability. S3 declares $OMEGA BLOCK(2) over the
    # absorption rate constant (ETA(4)) and the lag time (ETA(5)).
    # Table 2 reports the same quantities as %CV via
    # CV = sqrt(exp(omega^2) - 1) * 100 and the correlation as
    # cov / sqrt(var_ka * var_tlag), both reproduced below.
    # ------------------------------------------------------------------
    etalka + etaltlag ~ c(0.611,
                          -0.0524, 0.0114)
    # S3 $OMEGA BLOCK(2): 0.611 ;4_IIV_KA / -0.0524 0.0114 ;5_IIV_LAG-Time
    #   sqrt(exp(0.611)  - 1) = 91.8% = Table 2 omega_ka   91.8 [65.2, 129.1]
    #   sqrt(exp(0.0114) - 1) = 10.7% = Table 2 omega_tlag 10.7 [4.9, 20.6]
    #   -0.0524 / sqrt(0.611 * 0.0114) = -0.63 = Table 2 'Covariance between
    #   IIV on ka and tlag' -0.052 (-63%) [-0.11, -0.013]

    # S3 $OMEGA 1 (IIV on CL), 2 (IIV on Vc) and 3 (IIV on F1) are all
    # '0 FIX' and Table 2 reports no omega row for any of them: the
    # between-subject variability in clindamycin clearance is carried
    # entirely by the inter-occasion term below. The three zero-variance
    # etas are omitted rather than written as `~ fixed(0)` because a
    # zero-variance diagonal makes OMEGA singular and breaks the
    # Cholesky sampler used by rxSolve (same treatment as
    # Thoueille_2026_salmeterol.R and Wattanakul_2024_primaquine_motherinfant.R).

    # ------------------------------------------------------------------
    # Inter-occasion variability on apparent clearance across the six
    # 12-hourly doses. S3 declares one estimated BLOCK(1) variance
    # followed by five SAME records, so every occasion shares it:
    # occasion 1 carries the estimate and occasions 2-6 are fixed to it.
    # ------------------------------------------------------------------
    etaiov_cl_1 ~ 0.0798         # S3 $OMEGA BLOCK(1) 0.0798 ; IOV for CL; sqrt(exp(0.0798) - 1) = 28.8% = Table 2 kappa_CL/F 28.8 [23.7, 34.5]
    etaiov_cl_2 ~ fixed(0.0798)  # S3 $OMEGA BLOCK(1) SAME: shares the occasion-1 variance
    etaiov_cl_3 ~ fixed(0.0798)  # S3 $OMEGA BLOCK(1) SAME: shares the occasion-1 variance
    etaiov_cl_4 ~ fixed(0.0798)  # S3 $OMEGA BLOCK(1) SAME: shares the occasion-1 variance
    etaiov_cl_5 ~ fixed(0.0798)  # S3 $OMEGA BLOCK(1) SAME: shares the occasion-1 variance
    etaiov_cl_6 ~ fixed(0.0798)  # S3 $OMEGA BLOCK(1) SAME: shares the occasion-1 variance

    # ------------------------------------------------------------------
    # Residual error. S3 $ERROR is combined proportional plus additive
    # on the linear concentration scale:
    # Y = IPRED + IPRED*EPS(1) + EPS(2). $SIGMA holds variances, so each
    # SD below is the square root; Table 2 tabulates the same values as
    # a %CV and an mg/L SD respectively.
    # ------------------------------------------------------------------
    propSd <- 0.3286   ; label("Proportional residual error (fraction)")                     # S3 $SIGMA 0.108 ;prop Error; sqrt(0.108) = 0.3286 = Table 2 sigma_prop 32.9 %CV [27.5, 39.2]
    addSd  <- 0.004111 ; label("Additive residual error (mg/L)")                             # S3 $SIGMA 1.69E-05 ;add Error; sqrt(1.69e-05) = 0.004111 = Table 2 sigma_add 0.0041 mg/L [0.003, 0.006]
  })

  model({
    # 1. Decompose the integer occasion column into binary indicators to
    #    multiplex the six inter-occasion-variability etas onto log-CL.
    #    For single-occasion data pass OCC = 1 so the first eta applies.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)

    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5 + oc6 * etaiov_cl_6

    # 2. Individual parameters. Allometric weight scaling on both
    #    disposition parameters; clearance additionally carries the
    #    occasion-specific IOV eta. Neither CL/F nor V/F carries an IIV
    #    eta (see the ini() note on the '0 FIX' omegas).
    cl <- exp(lcl + iov_cl) * (WT / 29.05)^e_wt_cl
    vc <- exp(lvc)          * (WT / 29.05)^e_wt_vc
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag + etaltlag)

    # 3. Micro-constant.
    kel <- cl / vc

    # 4. ODE system: first-order absorption from an oral depot into a
    #    one-compartment disposition model (NONMEM ADVAN2 TRANS2).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 5. Relative bioavailability and absorption lag applied to the
    #    depot (S3 $PK: F1 = TVF1 * EXP(ETA(3)) with ETA(3) fixed to
    #    zero, ALAG1 = TVALAG1 * EXP(ETA(5))).
    f(depot)    <- exp(lfdepot)
    alag(depot) <- tlag

    # 6. Observation. Dose is in mg and V/F is in L (S3 sets S2 = V, so
    #    NONMEM's F is already a concentration), giving mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
