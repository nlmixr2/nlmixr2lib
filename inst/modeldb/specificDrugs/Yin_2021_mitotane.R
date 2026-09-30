Yin_2021_mitotane <- function() {
  description <- "Two-compartment population PK model with first-order absorption for oral mitotane in adults with adrenocortical carcinoma (Yin 2021), final pharmacogenetic model. Apparent clearance carries a power effect of lean body weight (Boer formula, derived from weight, height and sex) plus multiplicative genotype effects of CYP2C19*2 (rs4244285) carriage, SLCO1B3 699A>G (rs7311358) G-allele carriage and SLCO1B1 571T>C (rs4149057) CC or TT genotype (TC reference); apparent central volume carries a power effect of fat amount (total body weight minus lean body weight). Log-normal interindividual variability on CL/F, Vc/F, Vp/F and Q/F, interoccasion variability on CL/F with one occasion per 200 days of treatment, and combined additive plus proportional residual error. The absorption rate constant is fixed. Time in days. Yin_2021_mitotane_nogenotype is the authors' alternative model without genotype covariates, for patients whose genotype is unknown."
  reference <- paste(
    "Yin A, Ettaieb MHT, Swen JJ, van Deun L, Kerkhofs TMA,",
    "van der Straaten RJHM, Corssmit EPM, Gelderblom H, Kerstens MN,",
    "Feelders RA, Eekhoff M, Timmers HJLM, D'Avolio A, Cusato J,",
    "Guchelaar HJ, Haak HR, Moes DJAR. Population pharmacokinetic and",
    "pharmacogenetic analysis of mitotane in patients with adrenocortical",
    "carcinoma: towards individualized dosing. Clin Pharmacokinet.",
    "2021;60(1):89-102. doi:10.1007/s40262-020-00913-y.",
    "Parameter estimates from Table 2 (final model) and its footnotes a-c;",
    "IIV / IOV / residual-error equations (Eqs. S1-S2) and covariate forms",
    "(Eqs. S4-S5) from Online Resource 1. The Boer lean-body-weight",
    "coefficients, the standard-deviation scale of the printed CV% values and",
    "the ODE system are taken from the authors' published Shiny app script",
    "(github.com/AnyueYin/Shiny-app-script-for-model-simulation---Population-PK-and-PG-analysis-of-mitotane).",
    sep = " "
  )
  vignette <- "Yin_2021_mitotane"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "mitotane", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "mitotane", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "mitotane", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight at the start of mitotane treatment",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline value; Table 1 mean 80.0 kg (SD 15.9, range 52.5-120; 2 patients without a record were assigned the cohort median). Enters only through the derived lean body weight and fat amount; carries no effect of its own.",
      source_name = "WT"
    ),
    HT = list(
      description = "Body height at the start of mitotane treatment",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline value; Table 1 mean 172 cm (SD 10.0, range 154-193; 5 patients without a record were assigned the cohort median). Used only in the Boer lean-body-weight formula.",
      source_name = "HT"
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Selects the sex-specific Boer lean-body-weight formula: male LBW = 0.407*WT + 0.267*HT - 19.2, female LBW = 0.252*WT + 0.473*HT - 48.3 (Section 2.2 cites Boer; coefficients as coded in the authors' Shiny app). Sex was screened on the volumes but not retained. Cohort 27 of 48 female.",
      source_name = "SEX"
    ),
    OCC = list(
      description = "Occasion index for interoccasion variability on CL/F",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Every 200 days of treatment defines an occasion (Section 2.4, Table 2 footnote c). Decomposed inside model() into indicators occ1..occ8 selecting per-occasion IOV etas on CL. Eight occasions cover 0-1600 days; the cohort median treatment duration is 713.5 days (range 90-2856). Records outside 1-8 receive no IOV. See the vignette Errata.",
      source_name = "OCC"
    ),
    CYP2C19_S2_CARRIER = list(
      description = "CYP2C19*2 (rs4244285) loss-of-function allele carrier indicator: 1 = GA or AA, 0 = GG",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (GG, no *2 allele)",
      notes = "Dominant encoding per Table 2 row CL_SNP1 (GA/AA) = 0.551: carriers have 44.9% lower CL/F. Genotyped on the Affymetrix DMET Plus array. In this cohort rs4244285 was in 100% linkage disequilibrium with CYP2C18 1154C>T (rs2281891); the authors retained CYP2C19*2 as the functionally established variant.",
      source_name = "SNP1 CYP2C19*2 (rs4244285)"
    ),
    SNP_SLCO1B3_RS7311358_G_CARRIER = list(
      description = "SLCO1B3 699A>G (rs7311358, I233M) G-allele carrier indicator: 1 = AG or GG, 0 = AA",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (AA)",
      notes = "Dominant encoding per Table 2 row CL_SNP2 (AG/GG) = 0.601: G carriers have 39.9% lower CL/F. In this cohort rs7311358 was in 100% linkage disequilibrium with SLCO1B3 334G>T (rs4149117) and 1557G>A (rs2053098).",
      source_name = "SNP2 SLCO1B3 699A>G (rs7311358)"
    ),
    SNP_SLCO1B1_RS4149057_CC = list(
      description = "SLCO1B1 571T>C (rs4149057) homozygous CC genotype indicator: 1 = CC, 0 otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0 with SNP_SLCO1B1_RS4149057_TT = 0, i.e. the heterozygous TC group",
      notes = "Table 2 row CL_SNP3 (CC) = 0.753 relative to the TC reference. Eq. S5 sets the reference category multiplier to 1; the Shiny app assigns the reference to TC and estimates separate multipliers for CC and TT.",
      source_name = "SNP3 SLCO1B1 571T>C (rs4149057), CC genotype"
    ),
    SNP_SLCO1B1_RS4149057_TT = list(
      description = "SLCO1B1 571T>C (rs4149057) homozygous TT genotype indicator: 1 = TT, 0 otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0 with SNP_SLCO1B1_RS4149057_CC = 0, i.e. the heterozygous TC group",
      notes = "Table 2 row CL_SNP3 (TT) = 2.49 relative to the TC reference: TT patients have 2.49-fold higher CL/F than TC. The TC heterozygote is the reference (multiplier 1).",
      source_name = "SNP3 SLCO1B1 571T>C (rs4149057), TT genotype"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 48,
    n_studies = 1,
    disease_state = "adrenocortical carcinoma (ENSAT I-IV)",
    age_range = "22.6-76.8 years (mean 52.0)",
    weight_range = "52.5-120 kg (mean 80.0)",
    sex_female_pct = 56.3,
    dose_range = "0.5-16 g/day total daily dose, divided into 1-4 (occasionally 5, 6 or 8) administrations",
    regions = "Netherlands (Dutch Adrenal Network Registry)",
    notes = paste(
      "Retrospective TDM data, 914 mitotane plasma concentrations (33 below the",
      "2 mg/L LLOQ, omitted). Median 16.5 samples/patient (range 2-47); median",
      "treatment duration 713.5 days (range 90-2856). Multiple daily doses were",
      "simplified to a single daily dose equal to the total daily dose taken at",
      "08:00. Concentration target window 14-20 mg/L."
    )
  )

  ini({
    lka <- fixed(log(15.0))
    label("Absorption rate constant KA (1/day); Table 2, estimated on an absorption sub-dataset then held at 15.0")
    lcl <- log(298)
    label("Apparent clearance CL/F (L/day); Table 2 final model = 298")
    lvc <- log(6210)
    label("Apparent central volume Vc/F (L); Table 2 final model = 6210")
    lq <- log(883)
    label("Apparent intercompartmental clearance Q/F (L/day); Table 2 final model = 883")
    lvp <- log(18100)
    label("Apparent peripheral volume Vp/F (L); Table 2 final model = 18100")

    e_cyp2c19_s2_cl <- 0.551
    label("CL/F multiplier for CYP2C19*2 (GA/AA) carriers vs GG (unitless); Table 2 CL_SNP1 = 0.551")
    e_slco1b3_g_cl <- 0.601
    label("CL/F multiplier for SLCO1B3 699A>G (AG/GG) G carriers vs AA (unitless); Table 2 CL_SNP2 = 0.601")
    e_slco1b1_cc_cl <- 0.753
    label("CL/F multiplier for SLCO1B1 571T>C CC vs TC reference (unitless); Table 2 CL_SNP3 (CC) = 0.753")
    e_slco1b1_tt_cl <- 2.49
    label("CL/F multiplier for SLCO1B1 571T>C TT vs TC reference (unitless); Table 2 CL_SNP3 (TT) = 2.49")
    e_lbw_cl <- 1.10
    label("Power exponent of lean body weight on CL/F, (LBW/56.6)^e (unitless); Table 2 CL_LBW = 1.10, footnote a reference LBW 56.6 kg")
    e_fat_vc <- 1.22
    label("Power exponent of fat amount on Vc/F, (FAT/23.6)^e (unitless); Table 2 Vc_FAT = 1.22, footnote b reference FAT 23.6 kg")

    # IIV variances = (printed CV%/100)^2; Table 2 CV% are on the SD scale of
    # the log-normal eta (the authors' Shiny app uses CV%/100 as the rnorm SD
    # in exp(eta)), not the sqrt(exp(w^2)-1) log-normal CV.
    etalcl ~ 0.1849 # Table 2 IIV CL/F 43.0% -> SD 0.430
    etalvc ~ 0.222784 # Table 2 IIV Vc/F 47.2% -> SD 0.472
    etalq ~ 0.946729 # Table 2 IIV Q/F 97.3% -> SD 0.973
    etalvp ~ 0.788544 # Table 2 IIV Vp/F 88.8% -> SD 0.888

    # IOV on CL/F, one occasion per 200 days; Table 2 IOV 31.6% -> SD 0.316.
    # Occasions 2-8 share the occasion-1 variance (NONMEM $OMEGA BLOCK(1) SAME).
    etaiov_cl_1 ~ 0.099856
    etaiov_cl_2 ~ fixed(0.099856)
    etaiov_cl_3 ~ fixed(0.099856)
    etaiov_cl_4 ~ fixed(0.099856)
    etaiov_cl_5 ~ fixed(0.099856)
    etaiov_cl_6 ~ fixed(0.099856)
    etaiov_cl_7 ~ fixed(0.099856)
    etaiov_cl_8 ~ fixed(0.099856)

    propSd <- 0.166
    label("Proportional residual error (Table 2 PRO CV% = 16.6, SD scale)")
    addSd <- 0.920
    label("Additive residual error (mg/L) (Table 2 ADD = 0.920)")
  })

  model({
    # Boer lean body weight (kg) selected by sex, then fat amount (kg).
    # Coefficients from the authors' Shiny app (Section 2.2 cites Boer 1984).
    lbw <- SEXF * (0.252 * WT + 0.473 * HT - 48.3) +
      (1 - SEXF) * (0.407 * WT + 0.267 * HT - 19.2)
    fat <- WT - lbw

    # Multiplicative genotype effects on CL/F (Table 2 footnote a; the
    # reference genotypes GG, AA and TC each contribute a factor of 1).
    fgeno_cl <- ((1 - CYP2C19_S2_CARRIER) + CYP2C19_S2_CARRIER * e_cyp2c19_s2_cl) *
      ((1 - SNP_SLCO1B3_RS7311358_G_CARRIER) +
        SNP_SLCO1B3_RS7311358_G_CARRIER * e_slco1b3_g_cl) *
      ((1 - SNP_SLCO1B1_RS4149057_CC - SNP_SLCO1B1_RS4149057_TT) +
        SNP_SLCO1B1_RS4149057_CC * e_slco1b1_cc_cl +
        SNP_SLCO1B1_RS4149057_TT * e_slco1b1_tt_cl)

    # Interoccasion variability on CL/F: every 200 days is one occasion.
    occ1 <- (OCC == 1)
    occ2 <- (OCC == 2)
    occ3 <- (OCC == 3)
    occ4 <- (OCC == 4)
    occ5 <- (OCC == 5)
    occ6 <- (OCC == 6)
    occ7 <- (OCC == 7)
    occ8 <- (OCC == 8)
    iov_cl <- occ1 * etaiov_cl_1 + occ2 * etaiov_cl_2 +
      occ3 * etaiov_cl_3 + occ4 * etaiov_cl_4 +
      occ5 * etaiov_cl_5 + occ6 * etaiov_cl_6 +
      occ7 * etaiov_cl_7 + occ8 * etaiov_cl_8

    ka <- exp(lka)
    cl <- exp(lcl + etalcl + iov_cl) * fgeno_cl * (lbw / 56.6)^e_lbw_cl
    vc <- exp(lvc + etalvc) * (fat / 23.6)^e_fat_vc
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
