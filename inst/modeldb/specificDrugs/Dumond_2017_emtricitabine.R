Dumond_2017_emtricitabine <- function() {
  description <- paste(
    "Two-compartment population PK model for plasma emtricitabine (FTC) after",
    "oral emtricitabine 200 mg once daily in HIV-infected adults (age 22-73",
    "years), with first-order absorption, linked to a",
    "peripheral-blood-mononuclear-cell (PBMC) emtricitabine-triphosphate",
    "(FTC-tp) compartment. Total plasma clearance CL is split into a fraction",
    "FRAC leaving the body and a fraction 1 - FRAC entering the PBMC",
    "compartment; FTC-tp leaves the PBMC compartment by first-order",
    "elimination (CLpbmc) that rises exponentially with log2 p16INK4a",
    "expression in peripheral blood CD3+ T cells, a marker of cellular",
    "senescence. Creatinine clearance is a power covariate on CL and age a",
    "power covariate on central volume. Doses and concentrations are molar",
    "(umol, umol/L); the PBMC volume is fixed to 1 L, so CLpbmc enters as the",
    "rate constant keff_pbmc."
  )
  reference <- paste(
    "Dumond JB, Collins JW, Cottrell ML, Trezza CR, Prince HMA, Sykes C,",
    "Torrice C, White N, Malone S, Wang R, Patterson KB, Sharpless NE,",
    "Forrest A. p16INK4a, a Senescence Marker, Influences",
    "Tenofovir/Emtricitabine Metabolite Disposition in HIV-Infected Subjects.",
    "CPT Pharmacometrics Syst Pharmacol. 2017;6(2):120-127.",
    "doi:10.1002/psp4.12150"
  )
  vignette <- "Dumond_2017_tenofovir_emtricitabine"

  # Methods, "Pharmacokinetic analysis": "Data were analyzed as micromoles/L."
  # Figure 3 labels both panels 'Concentration (micromoles/L)'. One 200 mg
  # emtricitabine capsule (MW 247.24 g/mol) is 200 / 247.24 * 1000 =
  # 808.93 umol.
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  # The PBMC state is an amount in a compartment whose volume the authors fixed
  # to 1 (Results, "Pharmacokinetic modeling"). The supplementary control stream
  # (ADVAN5) has no S3 scale, so the predicted FTC-tp observation is the state
  # itself, on the 'micromoles/L' axis of Figure 3.
  compartmentData <- list(
    depot = list(analyte = "emtricitabine", units = "umol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "emtricitabine", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "emtricitabine", units = "umol", specimen = "plasma", verified = TRUE),
    pbmc_ftctp = list(analyte = "emtricitabine triphosphate", units = "umol", specimen = "blood cell", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault equation using total body weight",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Raw Cockcroft-Gault creatinine clearance on total body weight, NOT",
        "BSA-normalized (Methods, 'Clinical study conduct'; Table 1 footnote).",
        "Power effect on CL centred on 107.24 mL/min, the value in the control",
        "stream ('CLCrCL = THETA(9)*LOG(CrCL/107.24)', added to log CL). Table 2",
        "footnote b mislabels the centring value as '2.04'; see the vignette."
      ),
      source_name = "CrCL"
    ),
    AGE = list(
      description = "Chronologic age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on central volume centred on the cohort median of 49",
        "years: Vc = Vc,pop * (AGE/49)^2.47 (Table 2 footnote c; control",
        "stream 'V2AGE = THETA(10)*LOG(AGE/49.00)'). Described as 'marginally",
        "significant' in Results."
      ),
      source_name = "AGE"
    ),
    P16INK4A_TCELL = list(
      description = "Relative p16INK4a mRNA expression in peripheral blood CD3+ T cells (linear scale)",
      units = "(relative expression, linear scale)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The paper and its control stream carry the log2-transformed",
        "expression (column LTSI; Table 1 'Log2 (p16INK4a)' median 2.0, range",
        "0.16-3.9). The register stores gene expression on the linear scale,",
        "so supply 2^LTSI here; the model takes log2() internally. Enters the",
        "PBMC elimination rate as exp(0.495 * (log2(P16INK4A_TCELL) - 2.04))",
        "(Table 2 footnote a; control stream 'KpLTSI = THETA(8)*(LTSI-2.04)').",
        "A missing value (LTSI = -99, all intensively sampled subjects) set",
        "the effect to 0, i.e. the median; supply 2^2.04 for such subjects."
      ),
      source_name = "LTSI"
    )
  )

  covariatesDataExcluded <- list(
    FRAILTY = list(
      description = "Fried frailty phenotype (count of positive components, 0-5)",
      units = "(count)",
      type = "categorical",
      notes = paste(
        "Tested as binary (0 vs 1+ positive components) and as three",
        "categories (non-frail, prefrail, frail) on CL and central volume; not",
        "significant (Results)."
      )
    ),
    REGIMEN_ATVR = list(
      description = "Background regimen (efavirenz vs atazanavir/ritonavir)",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Background antiretroviral regimen was tested and not retained, so the",
        "efavirenz and atazanavir/ritonavir arms were comodelled (Results)."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 91,
    n_studies = 1,
    age_range = "22-73 years",
    age_median = "49 years",
    bmi_range = "17.3-44.3 kg/m^2",
    bmi_median = "28.3 kg/m^2",
    sex_female_pct = 33,
    race_ethnicity = c(`African American` = 58, White = 35, Other = 7),
    disease_state = "HIV-1 infection on stable antiretroviral therapy",
    renal_function = "Cockcroft-Gault creatinine clearance median 107 (range 43-227) mL/min; >30 mL/min required",
    p16ink4a = "log2 p16INK4a expression in CD3+ T cells median 2.0 (range 0.16-3.9)",
    frailty = "65 (72%) zero, 14 (15%) one and 8 (9%) two or more positive Fried frailty components",
    dose_range = "tenofovir disoproxil fumarate 300 mg + emtricitabine 200 mg orally once daily with efavirenz 600 mg (n = 60) or atazanavir/ritonavir 300/100 mg (n = 31), at steady state (>= 2 weeks)",
    regions = "United States (Chapel Hill and Greensboro, North Carolina)",
    notes = paste(
      "NCT01180075. Table 1 demographics. 79 participants contributed four",
      "samples over a dosing interval (predose, 2 h, 4-6 h, 10-14 h) and 12",
      "(six per regimen, age >= 55, non-frail) contributed 11 samples; the",
      "intensively sampled participants had no p16INK4a measurement.",
      "Plasma emtricitabine and PBMC emtricitabine triphosphate were modelled",
      "jointly; two samples below the limit of quantification were treated as",
      "missing."
    )
  )

  ini({
    # Table 2, FTC/FTC-tp model, 'Final NONMEM estimate (RSE%)' column.
    lka <- log(0.298); label("Absorption rate constant ka (1/h)") # Table 2 'ka (1/hr)' = 0.298 (RSE 9%)
    lcl <- log(15.6); label("Apparent total clearance CL/F out of the central compartment at CRCL = 107.24 mL/min (L/h)") # Table 2 'CL (L/hr)' = 15.6 (RSE 2%)
    lvc <- log(19.7); label("Apparent central volume Vc/F at age 49 years (L)") # Table 2 'Vc (L)' = 19.7 (RSE 11%)
    lq <- log(7.77); label("Apparent intercompartmental clearance CLd/F (L/h)") # Table 2 'CLd (L/hr)' = 7.77 (RSE 13%)
    lvp <- log(98.5); label("Apparent peripheral volume Vp/F (L)") # Table 2 'Vp (L)' = 98.5 (RSE 5%)
    frac <- 0.986; label("Fraction of CL leaving the body; 1 - FRAC enters the PBMC FTC-tp compartment (fraction)") # Table 2 'FRAC' = 0.986 (RSE 0%; bootstrap 0.977-0.992)
    lkeff_pbmc <- log(2.70e-2); label("PBMC FTC-tp elimination CLpbmc / Vpbmc at log2 p16INK4a = 2.04 (1/h; L/h with Vpbmc = 1 L)") # Table 2 'CLpbmc (L/hr)' = 2.70*10^-2 (RSE 2%)

    e_crcl_cl <- 0.551; label("Power exponent of CRCL/107.24 on CL (unitless)") # Table 2 'CrCL on CLp' = 0.551 (RSE 22%), footnote b
    e_age_vc <- 2.47; label("Power exponent of AGE/49 on Vc (unitless)") # Table 2 'Age on Vc' = 2.47 (RSE 38%), footnote c
    e_p16ink4a_tcell_keff_pbmc <- 0.495; label("Exponential coefficient of (log2 p16INK4a - 2.04) on CLpbmc (unitless)") # Table 2 'p16INK4a on CLpmbc' = 0.495 (RSE 21%), footnote a

    # IIV: Table 2 reports %CV; omega^2 = (CV/100)^2, the reading the TFV block
    # of the same table proves and the deposited initial $OMEGA supports
    # (initial 0.446 / 0.117 / 0.288 against 0.460 / 0.120 / 0.289). Off-diagonal
    # elements are the printed covariances (final estimates); the printed
    # correlations (-0.397, 0.411, 0.248) are those of the control stream's
    # initial block and do not close on the final diagonal (see the vignette).
    etalcl + etalvc + etalkeff_pbmc ~ c(
      0.119716,
      -0.07, 0.459684,
      0.05, 0.13, 0.289444
    ) # Table 2 'IIV CL' 34.6%, 'IIV Vc' 67.8%, 'IIV CLpbmc' 53.8%; covariances Vc~CL -0.07, CL~CLpbmc 0.05, Vc~CLpbmc 0.13

    propSd <- 0.387; label("Proportional residual error, plasma FTC (fraction)") # Table 2 'Plasma (proportional)' 38.7% (RSE 10%)
    propSd_Cpbmc_ftctp <- 0.424; label("Proportional residual error, PBMC FTC-tp (fraction)") # Table 2 'PBMC (proportional)' 42.4% (RSE 11%)
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (CRCL / 107.24)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (AGE / 49)^e_age_vc
    q <- exp(lq)
    vp <- exp(lvp)
    keff_pbmc <- exp(lkeff_pbmc + etalkeff_pbmc) *
      exp(e_p16ink4a_tcell_keff_pbmc * (log2(P16INK4A_TCELL) - 2.04))

    # Figure 1 and the control stream: CLp = FRAC * CL leaves the body and
    # CLpbi = (1 - FRAC) * CL carries FTC from plasma into the PBMC
    # compartment ('CLt=X*CL', 'CLpb=(1-X)*CL', 'K20=CLt/V2', 'K23=CLpb/V2').
    # With Vpbmc = 1 L, CLpbi / Vpbmc is numerically the plasma-concentration-
    # driven influx kinf_pbmc (1/h), as in the companion tenofovir model.
    kel <- frac * cl / vc
    kinf_pbmc <- (1 - frac) * cl
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - kinf_pbmc * (central / vc) -
      k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(pbmc_ftctp) <- kinf_pbmc * (central / vc) - keff_pbmc * pbmc_ftctp

    Cc <- central / vc
    # PBMC volume fixed to 1 L: the amount is the reported concentration.
    Cpbmc_ftctp <- pbmc_ftctp

    Cc ~ prop(propSd)
    Cpbmc_ftctp ~ prop(propSd_Cpbmc_ftctp)
  })
}
