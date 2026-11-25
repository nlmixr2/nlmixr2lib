Dumond_2017_tenofovir <- function() {
  description <- paste(
    "Two-compartment population PK model for plasma tenofovir (TFV) after",
    "oral tenofovir disoproxil fumarate 300 mg once daily in HIV-infected",
    "adults (age 22-73 years), with first-order absorption after a lag time,",
    "linked to a peripheral-blood-mononuclear-cell (PBMC) tenofovir-diphosphate",
    "(TFV-dp) compartment whose input is driven by the plasma TFV",
    "concentration (CLpbi) and whose first-order elimination (CLpbmc) rises",
    "exponentially with log2 p16INK4a expression in peripheral blood CD3+",
    "T cells, a marker of cellular senescence. Creatinine clearance is a power",
    "covariate on plasma clearance. Doses and concentrations are molar",
    "(umol, umol/L); the PBMC volume is fixed to 1 L, so CLpbi and CLpbmc",
    "enter as the rate constants kinf_pbmc and keff_pbmc."
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
  # Figure 2 labels both panels 'Concentration (micromoles/L)'. One 300 mg
  # tenofovir disoproxil fumarate tablet (MW 635.52 g/mol) delivers
  # 300 / 635.52 * 1000 = 472.06 umol of tenofovir; the Figure 2 VPC median
  # (peak about 0.8 umol/L, trough about 0.4 umol/L) is reproduced at that
  # molar dose (see the vignette).
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  # The PBMC state is an amount in a compartment whose volume the authors fixed
  # to 1 (Results, "Pharmacokinetic modeling": "As volume of the PBMC
  # compartment is unidentifiable for both models, this parameter was fixed to a
  # value of one"). The supplementary control stream has no S4 scale, so the
  # predicted TFV-dp observation is the state itself, on the
  # 'micromoles/L' axis of Figure 2.
  compartmentData <- list(
    depot = list(analyte = "tenofovir disoproxil", units = "umol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tenofovir", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tenofovir", units = "umol", specimen = "plasma", verified = TRUE),
    pbmc_tfvdp = list(analyte = "tenofovir diphosphate", units = "umol", specimen = "blood cell", verified = TRUE)
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
        "Power effect on plasma clearance centred on 107.24 mL/min, the value",
        "in the control stream ('CLCRCL = ((CRCL/107.24)**THETA(9))'). Table 2",
        "footnote b prints the same equation but closes with '2.04 is the",
        "population median estimated creatinine clearance', a copy-over from",
        "footnote a; 107.24 is the creatinine-clearance median (Table 1:",
        "107 mL/min, range 43-227)."
      ),
      source_name = "CRCL"
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
        "so supply 2^LTSI here (median about 4.1, range about 1.1-14.9); the",
        "model takes log2() internally. Measured by TaqMan qRT-PCR in CD3+",
        "T cells isolated from 7 mL K2EDTA whole blood (Methods, 'Analytical",
        "methods'; method of the paper's reference 23). Enters the PBMC",
        "elimination rate as exp(0.460 * (log2(P16INK4A_TCELL) - 2.04)); 2.04",
        "is the dataset median (Table 2 footnote a). In the control stream a",
        "missing value (LTSI = -99, all intensively sampled subjects) set the",
        "covariate factor to 1, i.e. the median; supply 2^2.04 for such",
        "subjects."
      ),
      source_name = "LTSI"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Chronologic age",
      units = "years",
      type = "continuous",
      notes = paste(
        "Tested and not retained for tenofovir: age 'was not found to be a",
        "significant covariate on drug CL for either TFV or FTC' (Results).",
        "Age is retained on Vc in the companion emtricitabine model."
      )
    ),
    FRAILTY = list(
      description = "Fried frailty phenotype (count of positive components, 0-5)",
      units = "(count)",
      type = "categorical",
      notes = paste(
        "Tested as binary (0 vs 1+ positive components) and as three",
        "categories (non-frail, prefrail, frail) on CL and central volume; not",
        "significant (Results). Control-stream columns FRYE / FRAY / FRAP."
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
      "Plasma tenofovir and PBMC tenofovir diphosphate were modelled jointly;",
      "four samples below the limit of quantification were treated as missing."
    )
  )

  ini({
    # Table 2, TFV/TFV-dp model, 'Final NONMEM estimate (RSE%)' column.
    lka <- log(0.623); label("Absorption rate constant ka (1/h)") # Table 2 'ka (1/hr)' = 0.623 (RSE 40%)
    lcl <- log(43.2); label("Apparent plasma clearance CLp/F at CRCL = 107.24 mL/min (L/h)") # Table 2 'CLp (L/hr)' = 43.2 (RSE 4%)
    lvc <- log(186); label("Apparent central volume Vc/F (L)") # Table 2 'Vc (L)' = 186 (RSE 59%)
    lq <- log(176); label("Apparent intercompartmental clearance CLd/F (L/h)") # Table 2 'CLd (L/hr)' = 176 (RSE 20%)
    lvp <- log(800); label("Apparent peripheral volume Vp/F (L)") # Table 2 'Vp (L)' = 800 (RSE 17%)
    ltlag <- log(0.472); label("Absorption lag time (h)") # Table 2 'Lag time (hr)' = 0.472 (RSE 2%)
    lkinf_pbmc <- log(2.27e-3); label("Plasma-to-PBMC TFV-dp formation CLpbi / Vpbmc (1/h; L/h with Vpbmc = 1 L)") # Table 2 'CLpbi (L/hr)' = 2.27*10^-3 (RSE 116%)
    lkeff_pbmc <- log(2.20e-3); label("PBMC TFV-dp elimination CLpbmc / Vpbmc at log2 p16INK4a = 2.04 (1/h; L/h with Vpbmc = 1 L)") # Table 2 'CLpbmc (L/hr)' = 2.20*10^-3 (RSE 116%)

    e_crcl_cl <- 0.698; label("Power exponent of CRCL/107.24 on CLp (unitless)") # Table 2 'CrCL on CLp' = 0.698 (RSE 20%), footnote b
    e_p16ink4a_tcell_keff_pbmc <- 0.460; label("Exponential coefficient of (log2 p16INK4a - 2.04) on CLpbmc (unitless)") # Table 2 'p16INK4a on CLpmbc' = 0.460 (RSE 23%), footnote a

    # IIV: Table 2 reports %CV; omega^2 = (CV/100)^2. The printed covariance
    # 0.22 with correlation 0.434 for CLp~Vc closes only on this reading
    # (0.22 / sqrt(0.1399 * 1.8496) = 0.433); the lognormal reading
    # omega^2 = log(1 + CV^2) would give a correlation of 0.59.
    etalcl + etalvc ~ c(0.139876, 0.22, 1.8496) # Table 2 'IIV CLp' 37.4%, 'IIV Vc' 136%, 'Covariance CLp~Vc' 0.22
    etalvp ~ 0.119716 # Table 2 'IIV Vp' 34.6% (shrinkage 65%)
    etalkeff_pbmc ~ 0.110224 # Table 2 'IIV CLpbmc' 33.2% (shrinkage 46%)
    etalkinf_pbmc ~ 0.199809 # Table 2 'IIV CLpbi' 44.7% (shrinkage 25%)

    propSd <- 0.283; label("Proportional residual error, plasma TFV (fraction)") # Table 2 'Plasma (proportional)' 28.3% (RSE 10%)
    propSd_Cpbmc_tfvdp <- 0.387; label("Proportional residual error, PBMC TFV-dp (fraction)") # Table 2 'PBMC (proportional)' 38.7% (RSE 8%)
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (CRCL / 107.24)^e_crcl_cl
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp)
    kinf_pbmc <- exp(lkinf_pbmc + etalkinf_pbmc)
    keff_pbmc <- exp(lkeff_pbmc + etalkeff_pbmc) *
      exp(e_p16ink4a_tcell_keff_pbmc * (log2(P16INK4A_TCELL) - 2.04))

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Supplementary control stream ($DES): TFV-dp formation is driven by the
    # plasma concentration ('DADT(4)=KI*CP-K40*A(4)') and is NOT subtracted
    # from the plasma central compartment. Figure 1's printed equations
    # subtract CLpbi * Central from central as well; the fitted code is
    # followed here (the difference is 2.27e-3 against 43.2 L/h, 0.005% of
    # plasma clearance).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(pbmc_tfvdp) <- kinf_pbmc * (central / vc) - keff_pbmc * pbmc_tfvdp

    alag(depot) <- exp(ltlag)

    Cc <- central / vc
    # PBMC volume fixed to 1 L: the amount is the reported concentration.
    Cpbmc_tfvdp <- pbmc_tfvdp

    Cc ~ prop(propSd)
    Cpbmc_tfvdp ~ prop(propSd_Cpbmc_tfvdp)
  })
}
