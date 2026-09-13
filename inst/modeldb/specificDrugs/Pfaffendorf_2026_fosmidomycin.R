# Population pharmacokinetic model for oral fosmidomycin given with
# clindamycin and artesunate to Gabonese children and adults with
# uncomplicated Plasmodium falciparum malaria
# (Pfaffendorf 2026, Malaria Journal 25:152;
# doi:10.1186/s12936-026-05872-6).

Pfaffendorf_2026_fosmidomycin <- function() {
  description <- paste(
    "Population PK model for oral fosmidomycin co-administered with",
    "clindamycin and artesunate in Gabonese children and adults with",
    "uncomplicated Plasmodium falciparum malaria (Pfaffendorf 2026).",
    "One-compartment disposition with linear elimination, first-order",
    "absorption and an absorption lag time. Allometric body-weight",
    "scaling is applied with fixed exponents 0.75 on apparent clearance",
    "and 1 on apparent volume, centered at the cohort median weight of",
    "29.05 kg. Body temperature increases apparent clearance as a power",
    "function centered at the cohort mean 37.1 degC. Relative",
    "bioavailability is fixed at 1 and carries log-normal IIV; further",
    "IIV is carried on the absorption rate constant. Residual error is",
    "combined proportional plus additive.",
    sep = " "
  )
  reference <- paste(
    "Pfaffendorf C, Dejon-Agobe JC, Edoa JR, Maiga-Ascofare O,",
    "Ahenkan E, Adegnika AA, Ramharter M, Wicha SG, Mischlinger J (2026).",
    "Population pharmacokinetics of fosmidomycin and clindamycin in",
    "combination with artesunate for uncomplicated Plasmodium falciparum",
    "malaria in Gabonese children and adults.",
    "Malaria Journal 25:152. doi:10.1186/s12936-026-05872-6.",
    "Parameter values from Table 2 (fosmidomycin block); model structure",
    "from Additional file 1 section S2 (final NONMEM control stream).",
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
        "CL_i = CL_TV * (WGT/29.05)^0.75 * (BT/37.1)^theta_TEMP and",
        "V_i  = V_TV  * (WGT/29.05)^1.",
        "The centering constant 29.05 kg is the study-population median",
        "weight used for the reference-corrected VPC (Methods, VPC",
        "paragraph) and appears verbatim in the final control stream",
        "(Additional file 1 S2). Note it differs slightly from the",
        "Table 1 median of 29.1 kg, which is the same quantity rounded",
        "to three significant figures.",
        sep = " "
      ),
      source_name        = "WGT"
    ),
    BODYTEMP = list(
      description        = "Body temperature at admission",
      units              = "degC",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Admission (screening) body temperature, treated as time-fixed",
        "per subject in this implementation. Pfaffendorf 2026 Table 1",
        "reports mean (SD) 37.1 (1.06) degC and median [min, max]",
        "36.8 [35.4, 39.0] degC. Applied as a POWER effect on apparent",
        "clearance centered at the cohort mean 37.1 degC:",
        "CL_i = CL_TV * (WGT/29.05)^0.75 * (BODYTEMP/37.1)^e_bodytemp_cl",
        "with e_bodytemp_cl = 6.68. This is a power form, unlike the",
        "exponential-deviation form used by Kloprogge_2014_quinine.R and",
        "Kloprogge_2013_lumefantrine.R for the same canonical covariate;",
        "the functional form follows the source control stream",
        "(Additional file 1 S2: TVCL = THETA(1) * ((WGT/29.05)**0.75) *",
        "((BT/37.1)**THETA(6))).",
        "The effect is steep and the authors caution against",
        "extrapolating outside the measured 35.4-39.0 degC range:",
        "temperature was recorded only twice daily, so the discrete",
        "measurements do not capture the full febrile time course",
        "(Discussion). Higher temperature INCREASES apparent clearance;",
        "the paper states the model predicts roughly a 60% increase at",
        "40 degC.",
        sep = " "
      ),
      source_name        = "BT"
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "fosmidomycin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "fosmidomycin", units = "mg", specimen = "plasma", verified = TRUE)
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
      "Oral fosmidomycin 30 mg/kg every 12 h for three days (six doses",
      "total), co-administered with clindamycin 10 mg/kg and artesunate",
      "2 mg/kg on the same schedule. Doses were rounded to the closest",
      "match achievable with the available 75 mg, 225 mg and 450 mg",
      "capsules (Nextpharma, Goettingen, Germany). The reference",
      "-corrected VPC used a 900 mg reference dose at 29.05 kg."
    ),
    regions         = "Gabon (Centre de Recherches Medicales de Lambarene)",
    trial_registration = "PACTR202008909968293 (pactr.samrc.ac.za)",
    notes           = paste(
      "Recruitment was stratified into three age groups: 20 patients",
      "aged 6 months-10 years, 10 aged 11-17 years, and 10 aged",
      "18-65 years. A total of 242 fosmidomycin plasma samples entered",
      "the PK analysis. Sampling was weight-banded to limit blood draws",
      "in young children: patients above 35 kg at 0.25, 0.5, 0.75, 1.5,",
      "3, 5, 8 h (day 0), 24 h, 48 h and day 7; patients 20-35 kg at 0,",
      "1, 5, 8 h (day 0), 24 h, 48 h and day 7; patients below 20 kg at",
      "0, 1, 8 h (day 0), 24 h, 48 h and day 7. Samples taken at the",
      "same time as a dose were drawn pre-dose. Bioanalysis by LC-MS/MS",
      "after trichloroacetic-acid protein precipitation, calibration",
      "range 0.25-15 mg/L; data below the limit of quantification were",
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
    # Structural parameters. Typical values are Table 2 (fosmidomycin
    # block, 'Typical value' column) and are identical to the $THETA
    # records of the final control stream in Additional file 1 S2, so
    # these are final estimates rather than initial values. They apply
    # to a patient of 29.05 kg with a body temperature of 37.1 degC.
    # Values are reported on the linear scale; log() is applied here for
    # the nlmixr2 internal scale.
    # ------------------------------------------------------------------
    lcl  <- log(58.3)  ; label("Apparent clearance CL/F at 29.05 kg and 37.1 degC (L/h)")    # Table 2 CL/F = 58.3 [51.0, 67.8]; S2 $THETA (0, 58.3) ;1_CL
    lvc  <- log(248)   ; label("Apparent central volume V/F at 29.05 kg (L)")                # Table 2 V/F  = 248  [201.9, 302.8]; S2 $THETA (0, 248) ;2_V
    lka  <- log(0.698) ; label("Absorption rate constant ka (1/h)")                          # Table 2 ka   = 0.698 [0.48, 1.0]; S2 $THETA (0, 0.698) ;3_KA
    ltlag <- log(0.105); label("Absorption lag time (h)")                                    # Table 2 tlag = 0.105 [0.03, 0.17]; S2 $THETA (0, 0.105) ;4_ALAG1

    # Relative bioavailability is a structural anchor: the paper reports
    # F = '1 FIXED' (Table 2) and the control stream carries
    # $THETA (1) FIX ;5_F1. All of the between-subject variability in
    # exposure that would otherwise sit on CL/F and V/F is carried by
    # the IIV on F below.
    lfdepot <- fixed(log(1)) ; label("Relative bioavailability F (unitless)")                # Table 2 F = 1 FIXED; S2 $THETA (1) FIX ;5_F1

    # ------------------------------------------------------------------
    # Allometric exponents. The paper tested 'allometric scaling with
    # fixed exponents' (Methods, Covariate model building) and the
    # control stream hardcodes 0.75 and 1 rather than estimating them,
    # so neither appears in Table 2 and both are fixed here.
    # ------------------------------------------------------------------
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent on CL/F (unitless)")                 # Methods 'Covariate model building'; S2 $PK ((WGT/29.05)**0.75)
    e_wt_vc <- fixed(1)    ; label("Allometric exponent on V/F (unitless)")                  # Methods 'Covariate model building'; S2 $PK (WGT/29.05)**1

    # ------------------------------------------------------------------
    # Covariate effect. Body temperature enters CL/F as a power term,
    # NOT as an exponential deviation; see covariateData[[BODYTEMP]].
    # ------------------------------------------------------------------
    e_bodytemp_cl <- 6.68 ; label("Power exponent for body temperature on CL/F (unitless)")  # Table 2 theta_TEMP = 6.68 [3.6, 9.7]; S2 $THETA (-100, 6.68, 100000) ;6_CLBT1

    # ------------------------------------------------------------------
    # Inter-individual variability. S2 $OMEGA holds the log-normal
    # variances; Table 2 reports the same quantities as %CV via
    # CV = sqrt(exp(omega^2) - 1) * 100, reproduced below.
    # ------------------------------------------------------------------
    etalka      ~ 0.383   # S2 $OMEGA 0.383 ;2_IIV_KA; sqrt(exp(0.383) - 1) = 68.3% = Table 2 omega_ka 68.3 [45.8, 102.1]
    etalfdepot  ~ 0.114   # S2 $OMEGA 0.114 ;3_IIV_F1; sqrt(exp(0.114) - 1) = 34.7% = Table 2 omega_F  34.7 [24.3, 48.8]

    # S2 $OMEGA 1 (IIV on V) is '0 FIX' and Table 2 reports no omega row
    # for V or for CL, so the fosmidomycin model carries no IIV on
    # either disposition parameter: CL is written as plain `CL = TVCL`
    # in the control stream with no eta at all. The zero-variance V eta
    # is omitted rather than written as `etalvc ~ fixed(0)` because a
    # zero-variance diagonal makes OMEGA singular and breaks the
    # Cholesky sampler used by rxSolve (same treatment as
    # Thoueille_2026_salmeterol.R and Wattanakul_2024_primaquine_motherinfant.R).

    # ------------------------------------------------------------------
    # Residual error. S2 $ERROR is combined proportional plus additive
    # on the linear concentration scale:
    # Y = IPRED + IPRED*EPS(1) + EPS(2). $SIGMA holds variances, so each
    # SD below is the square root; Table 2 tabulates the same values as
    # a %CV and an mg/L SD respectively.
    # ------------------------------------------------------------------
    propSd <- 0.3564 ; label("Proportional residual error (fraction)")                       # S2 $SIGMA 0.127 ;Prop.Err_PG; sqrt(0.127) = 0.3564 = Table 2 sigma_prop 35.6 %CV [29.7, 41.9]
    addSd  <- 0.1972 ; label("Additive residual error (mg/L)")                               # S2 $SIGMA 0.0389 ;add.Err_PG; sqrt(0.0389) = 0.1972 = Table 2 sigma_add 0.197 mg/L [0.15, 0.25]
  })

  model({
    # 1. Individual parameters. Allometric weight scaling on both
    #    disposition parameters; body temperature as a power term on
    #    apparent clearance. Neither CL/F nor V/F carries an eta (see
    #    the ini() note on the '0 FIX' omegas).
    cl <- exp(lcl) * (WT / 29.05)^e_wt_cl * (BODYTEMP / 37.1)^e_bodytemp_cl
    vc <- exp(lvc) * (WT / 29.05)^e_wt_vc
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag)

    # 2. Micro-constant.
    kel <- cl / vc

    # 3. ODE system: first-order absorption from an oral depot into a
    #    one-compartment disposition model (NONMEM ADVAN2 TRANS2).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 4. Relative bioavailability and absorption lag applied to the
    #    depot (S2 $PK: F1 = TVF1 * EXP(ETA(3)), ALAG1 = THETA(4)).
    f(depot)    <- exp(lfdepot + etalfdepot)
    alag(depot) <- tlag

    # 5. Observation. Dose is in mg and V/F is in L (S2 sets S2 = V, so
    #    NONMEM's F is already a concentration), giving mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
