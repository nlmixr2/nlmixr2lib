Dhoro_2015_efavirenz <- function() {
  description <- "One-compartment population PK model for oral efavirenz in HIV-positive and HIV/TB co-infected adults in Zimbabwe (Dhoro 2015), with apparent clearance CL/F stratified by CYP2B6 983T>C (CYP2B6*18, rs28399499) genotype and multiplicative fractional covariate effects of CYP2B6 516G>T (CYP2B6*6, rs3745274) genotype, body weight, and sex on CL/F. Absorption rate constant ka and apparent volume V/F are fixed from the upstream Nyakutira 2008 Zimbabwean cohort."
  reference <- paste(
    "Dhoro M, Zvada S, Ngara B, Nhachi C, Kadzirange G, Chonzi P, Masimirembwa C.",
    "CYP2B6*6, CYP2B6*18, Body weight and sex are predictors of efavirenz",
    "pharmacokinetics and treatment response: population pharmacokinetic modeling",
    "in an HIV/AIDS and TB cohort in Zimbabwe.",
    "BMC Pharmacol Toxicol. 2015;16:4. doi:10.1186/s40360-015-0004-2."
  )
  vignette <- "Dhoro_2015_efavirenz"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (baseline) per the Dhoro 2015 single-sample sparse PK design. Linear-additive fractional shift on CL/F centred on the cohort median weight; the paper does not report the median value explicitly (Table 1 reports mean 61.5 kg for males / 57.9 kg for females, weighted mean 59.1 kg across n = 60 males + 125 females), and the model file uses 60 kg as the round-figure median. Effect magnitude: +21.1 % per +10 kg above 60 kg (Dhoro 2015 Results and Table 3).",
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Time-fixed per subject. Dhoro 2015 Table 3 reports a +22.2 % fractional shift on CL/F for females vs the male reference. Cohort distribution (Table 1): 60 males (32 %) / 125 females (68 %).",
      source_name        = "Sex (Male / Female)"
    ),
    SNP_CYP2B6_RS3745274_T_COUNT = list(
      description        = "Count of CYP2B6 c.516G>T (rs3745274, p.Q172H, CYP2B6*6) T-alleles per subject (0 / 1 / 2). 0 = GG homozygous wild-type, 1 = GT heterozygous, 2 = TT homozygous variant.",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (germline genotype). Dhoro 2015 Methods 'DNA extraction and TaqMan Genotyping' and Table 1 report this SNP as CYP2B6*6 (also known as rs3745274 / 516G>T). The Dhoro 2015 final-model parameterization uses the heterozygote (GT, count = 1) as the reference category, with multiplicative fractional shifts on CL/F for the homozygous groups: +93.1 % for GG (count = 0, wild-type) and -63.4 % for TT (count = 2, variant). This non-monotonic (U-shaped) encoding is unusual; the canonical count column reconstructs the indicators in model() as (count == 0) for GG and (count == 2) for TT. Cohort allele-genotype frequencies (Table 1, n = 185): GG 30.8 % (57/185), GT 45.4 % (84/185), TT 21.1 % (39/185).",
      source_name        = "CYP2B6 G516T (rs3745274) / CYP2B6*6"
    ),
    SNP_CYP2B6_RS28399499_C_COUNT = list(
      description        = "Count of CYP2B6 c.983T>C (rs28399499, p.I328T, CYP2B6*18) C-alleles per subject (0 / 1 / 2). 0 = TT homozygous wild-type (extensive metaboliser), 1 = TC heterozygous (intermediate metaboliser), 2 = CC homozygous variant (poor metaboliser).",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (germline genotype). Dhoro 2015 Methods 'DNA extraction and TaqMan Genotyping' and Table 1 report this SNP as CYP2B6*18 (also known as rs28399499 / 983T>C). The Dhoro 2015 final-model stratifies the apparent oral clearance CL/F by CYP2B6*18 genotype with three distinct typical-value estimates (Table 3): CL/F = 7.01 L/h for TT (extensive metaboliser, count = 0, reference), 2.26 L/h for TC (intermediate metaboliser, count = 1), and 0.539 L/h for CC (poor metaboliser, count = 2). The reference category in the packaged model file is TT (extensive metaboliser); the non-reference shifts are encoded as log-ratio additive effects on log-CL/F to keep the single etalcl IIV applied uniformly. Cohort allele-genotype frequencies (Table 1, n = 185): TT 71.4 % (132/185), TC 25.4 % (47/185), CC 3.2 % (6/185).",
      source_name        = "CYP2B6*18 (rs28399499)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in Dhoro 2015 covariate analysis (Table 4: drop in OFV = 0.29, df = 1, p = 0.59, IIV explained = 0.01 %) but not retained in the final CL/F model. Cohort range (Table 1): mean 38.9 years (8.5 SD).",
      source_name        = "Age"
    ),
    SNP_CYP2A6_RS8192726 = list(
      description        = "CYP2A6 c.-48T>G (rs8192726, CYP2A6*9) genotype indicator",
      units              = "(genotype: GG / TT in source)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Screened in Dhoro 2015 covariate analysis (Table 4: drop in OFV = 3.836, df = 1, p = 0.0502, IIV explained = 2.7 %) but not retained in the final CL/F model at the 1 % backward-elimination threshold. Not registered as a canonical covariate column because no downstream nlmixr2lib model retains it; documented here for source provenance.",
      source_name        = "CYP2A6*9 (rs8192726)"
    ),
    SNP_CYP2A6_RS28399454 = list(
      description        = "CYP2A6 c.5065G>A (rs28399454, CYP2A6*17) genotype indicator",
      units              = "(genotype: GG / GA / AA in source)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Screened in Dhoro 2015 covariate analysis (Table 4: drop in OFV = 1.44, df = 2, p = 0.23, IIV explained = 1.8 %) but not retained in the final CL/F model. Not registered as a canonical covariate column because no downstream nlmixr2lib model retains it; documented here for source provenance.",
      source_name        = "CYP2A6*17 (rs28399454)"
    ),
    SNP_ABCB1_RS1128503 = list(
      description        = "ABCB1 c.1236C>T (rs1128503) genotype indicator",
      units              = "(genotype: CC / CT / TT in source)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Screened in Dhoro 2015 covariate analysis (Table 2 univariate ANOVA p = 0.841) and rejected; not in the final CL/F model. Not registered as a canonical covariate column because no downstream nlmixr2lib model retains it; documented here for source provenance.",
      source_name        = "ABCB1 1236C/T (rs1128503)"
    ),
    EFV_RIF = list(
      description        = "Concomitant rifampicin-containing anti-tuberculosis therapy indicator (1 = HIV/TB co-treatment, 0 = HIV treatment only)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (HIV treatment only)",
      notes              = "Screened in Dhoro 2015 covariate analysis (Table 4: drop in OFV = 1.932, df = 1, p = 0.165, IIV explained = 1.1 %) and rejected; not in the final CL/F model. The cohort comprised 95 HIV-only and 90 HIV/TB co-infected patients on EFV + rifampicin-containing anti-TB therapy.",
      source_name        = "EFV-RIF interaction"
    ),
    CNS_TOX = list(
      description        = "Central-nervous-system toxicity status indicator (1 = patient has CNS adverse effects, 0 = no CNS adverse effects)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no CNS adverse effects)",
      notes              = "Screened in Dhoro 2015 covariate analysis as an exploratory covariate (Table 4: drop in OFV = 1.011, df = 1, p = 0.315, IIV explained = 0.14 %) and not retained. Used in the paper's discussion of dose-response for CNS adverse effects (CNS-affected patients had ~27 % lower CL/F than unaffected patients in a stratified analysis) but not in the final population-PK structural model.",
      source_name        = "CNS Toxicity"
    ),
    HT = list(
      description        = "Body height at baseline",
      units              = "m",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in Dhoro 2015 univariate covariate analysis (Methods + Table 1) and not retained in the final CL/F model. Cohort range (Table 1): male mean 1.72 m, female mean 1.61 m.",
      source_name        = "Height"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 185,
    n_studies      = 1,
    age_range      = "mean 38.9 years (SD 8.5); male mean 40.2 years (SD 9.1), female mean 38.3 years (SD 8.1)",
    weight_range   = "male mean 61.5 kg (SD 10.1), female mean 57.9 kg (SD 11.3); weighted-mean 59.1 kg across the full cohort",
    weight_median  = "60 kg (estimated; the source paper does not report the exact median used to centre the WT covariate effect)",
    sex_female_pct = 67.6,
    race_ethnicity = c(Zimbabwean = 100),
    disease_state  = "HIV-positive adults receiving efavirenz-based combination antiretroviral therapy. 95 patients had HIV monoinfection on ART only; 90 patients were HIV/TB co-infected receiving ART plus a rifampicin-containing anti-tuberculosis regimen. All patients were evaluated for CNS adverse effects (sleep disorders, hallucinations, headaches) classified per WHO ADR grading.",
    dose_range     = "600 mg orally once daily as part of EFV-based combination antiretroviral therapy. Common backbone regimens (Table 2): TDF/3TC/EFV, AZT/3TC/EFV, D4T/3TC/EFV. The packaged model and validation vignette also simulate the 200 mg and 400 mg daily reduced-dose scenarios investigated by Monte-Carlo simulation in the source paper.",
    regions        = "Zimbabwe (Wilkins Hospital and Chitungwiza Hospital, Harare)",
    cyp2b6_freq    = "CYP2B6*6 (rs3745274 / 516G>T): GG 30.8 %, GT 45.4 %, TT 21.1 %. CYP2B6*18 (rs28399499 / 983T>C): TT 71.4 %, TC 25.4 %, CC 3.2 % (Dhoro 2015 Table 1).",
    notes          = "Single sparse plasma sample per subject collected 12-15 h post-dose by reverse-phase HPLC with UV detection (felodipine internal standard). 185 patients (60 male, 125 female) recruited from two Harare hospitals. The structural absorption-rate constant ka and apparent volume V/F were FIXED in NONMEM VI from the upstream Nyakutira et al. (2008) popPK in Zimbabwean patients (paper reference [37]). Clearance was the only structural parameter estimated; covariate screening retained CYP2B6*18, CYP2B6*6, body weight, and sex (Table 4) explaining 55 % of inter-individual variability in CL/F. Estimation method: FOCE INTER in NONMEM VI."
  )

  ini({
    # ---- Structural PK parameters ----
    # Dhoro 2015 Table 3: ka and V/F were FIXED from the upstream Nyakutira et al. (2008)
    # Zimbabwean popPK [Dhoro 2015 reference 37]. Only CL/F was estimated.
    lka <- fixed(log(0.18)) ; label("First-order absorption rate constant ka (1/h)") # Dhoro 2015 Table 3: ka = 0.18 1/h FIX (fixed from Nyakutira 2008)
    lvc <- fixed(log(150))  ; label("Apparent volume of distribution V/F (L)")        # Dhoro 2015 Table 3: V/F = 150 L FIX (fixed from Nyakutira 2008; the table header "L/hr" is a typo - V/F is a volume in L)

    # ---- CL/F stratified by CYP2B6*18 (rs28399499 / 983T>C) genotype ----
    # Dhoro 2015 Table 3 reports CL/F as three distinct typical values, one per *18
    # genotype, with a single shared IIV (etalcl). The TT (extensive metaboliser,
    # count = 0) group is the natural reference; the non-reference groups are encoded
    # as log-ratio additive effects on log-CL/F so the single etalcl applies uniformly.
    lcl       <- log(7.01)               ; label("CL/F (L/h) at CYP2B6*18 TT (extensive metaboliser, reference)") # Dhoro 2015 Table 3: CL/F TT = 7.01 L/h (RSE 10 %)
    e_18tc_cl <- log(2.26 / 7.01)        ; label("Log-ratio of CL/F at CYP2B6*18 TC vs TT reference (unitless)")   # Dhoro 2015 Table 3: CL/F TC = 2.26 L/h (RSE 12 %); log(2.26/7.01) = -1.1320
    e_18cc_cl <- log(0.539 / 7.01)       ; label("Log-ratio of CL/F at CYP2B6*18 CC vs TT reference (unitless)")   # Dhoro 2015 Table 3: CL/F CC = 0.539 L/h (RSE 24 %); log(0.539/7.01) = -2.5658

    # ---- Multiplicative fractional covariate effects on CL/F ----
    # Methods: "For categorical covariates ... the covariate model was expressed as a
    # fractional change (theta_cov) from the estimate for a typical value (theta_P)";
    # "The effect of continuous covariates was parameterized centred on the median value".
    # Applied as: CL = CL_*18 * (1 + e_6gg_cl * is_6gg + e_6tt_cl * is_6tt)
    #                         * (1 + e_wt_cl  * (WT - 60) / 10)
    #                         * (1 + e_sexf_cl * SEXF)
    #
    # CYP2B6*6 (rs3745274 / 516G>T) effects -- reference is GT heterozygote (count = 1)
    e_6gg_cl  <-  0.931                  ; label("Fractional change on CL/F for CYP2B6*6 GG (wild-type) vs GT reference (unitless)") # Dhoro 2015 Table 3: CYP2B6 GG = +93.1 % (RSE 24 %)
    e_6tt_cl  <- -0.634                  ; label("Fractional change on CL/F for CYP2B6*6 TT (variant homozygote) vs GT reference (unitless)") # Dhoro 2015 Table 3: CYP2B6 TT = -63.4 % (RSE 9 %)

    # Body-weight effect -- linear additive on the fractional-shift scale, centred on
    # the cohort median weight. The source paper states the centring was at the median
    # but does not report the exact median value; the model file uses 60 kg as the
    # round-figure approximation of the cohort weighted-mean weight (Table 1 + Methods).
    e_wt_cl   <-  0.211                  ; label("Fractional change on CL/F per +10 kg above the 60 kg reference (unitless / 10 kg)") # Dhoro 2015 Table 3 / Results: +21.1 % per +10 kg (RSE 21 %)

    # Sex effect -- reference is male (SEXF = 0)
    e_sexf_cl <-  0.222                  ; label("Fractional change on CL/F for female (SEXF = 1) vs male reference (unitless)") # Dhoro 2015 Table 3: Females = +22.2 % (RSE 67 %)

    # ---- IIV (only CL/F has IIV; ka and V/F are fixed structural parameters) ----
    # Dhoro 2015 Table 3: IIV CL/F = 70.3 % CV (RSE 7 %).
    # omega^2 = log(CV^2 + 1) = log(0.703^2 + 1) = 0.40114
    etalcl ~ 0.40114

    # ---- Residual error (proportional only) ----
    # Dhoro 2015 Table 3 reports PROP_ERR = 0.12 without explicit unit annotation. Per
    # NONMEM VI convention, $SIGMA reports the variance; the proportional SD on the
    # linear concentration scale is therefore sqrt(0.12) = 0.3464 (~34.6 % CV). This
    # follows the same reading the co-author group used in the Olagunju 2018 efavirenz
    # popPK (variance 0.085 -> propSd 0.292) packaged elsewhere in nlmixr2lib.
    propSd <- sqrt(0.12) ; label("Proportional residual error (fraction)") # Dhoro 2015 Table 3: PROP_ERR = 0.12 (NONMEM $SIGMA variance); propSd = sqrt(0.12) = 0.3464
  })

  model({
    # 1. Derive non-additive genotype indicators from the per-allele count columns.
    #    CYP2B6*18 (rs28399499): stratification reference is TT (count == 0).
    is_18_TC <- (SNP_CYP2B6_RS28399499_C_COUNT == 1)
    is_18_CC <- (SNP_CYP2B6_RS28399499_C_COUNT == 2)
    #    CYP2B6*6 (rs3745274): reference is GT heterozygote (count == 1).
    is_6_GG  <- (SNP_CYP2B6_RS3745274_T_COUNT  == 0)
    is_6_TT  <- (SNP_CYP2B6_RS3745274_T_COUNT  == 2)

    # 2. Typical log-CL/F: TT-reference + log-ratio shifts for the two non-reference
    #    *18 groups (mutually exclusive indicators).
    ltvcl <- lcl + e_18tc_cl * is_18_TC + e_18cc_cl * is_18_CC

    # 3. Individual CL/F with multiplicative fractional covariate effects of CYP2B6*6,
    #    body weight, and sex. IIV (etalcl) is applied on the log-CL scale; the
    #    fractional-effect multipliers act on linear-scale CL.
    cl <- exp(ltvcl + etalcl) *
          (1 + e_6gg_cl * is_6_GG + e_6tt_cl * is_6_TT) *
          (1 + e_wt_cl  * (WT - 60) / 10) *
          (1 + e_sexf_cl * SEXF)
    vc <- exp(lvc)
    ka <- exp(lka)

    # 4. Micro-constants
    kel <- cl / vc

    # 5. ODE system: one-compartment model with first-order absorption
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 6. Observation and error
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
