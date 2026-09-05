Liang_2024_lenalidomide <- function() {
  description <- "One-compartment oral population PK model with first-order absorption and elimination for lenalidomide in Chinese patients with multiple myeloma, lymphoma or myelodysplastic syndrome, with ABCB1 3435C>T (rs1045642) T-allele count and fed state as covariates on the apparent volume of distribution (Liang 2024)"
  reference <- "Liang X, Shi H, Bi K, Feng S, Chen S, Zhao W, Huang X. Population pharmacokinetics of lenalidomide in Chinese patients with influence of genetic polymorphisms of ABCB1. Sci Rep. 2024;14:2577. doi:10.1038/s41598-024-52460-2"
  vignette <- "Liang_2024_lenalidomide"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot   = list(analyte = "lenalidomide", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "lenalidomide", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    SNP_ABCB1_RS1045642_T_COUNT = list(
      description        = "ABCB1 3435C>T (rs1045642) T-allele count",
      units              = "(count, 0/1/2 alleles per subject)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "0 = CC homozygous wild type, 1 = CT heterozygous, 2 = TT homozygous variant.",
        "Time-fixed per subject (germline genotype), determined by first-generation DNA",
        "sequencing (Liang 2024 Methods, 'Measurements'); genotype detection rate was 100%",
        "and the distribution satisfied Hardy-Weinberg equilibrium (Liang 2024 Table 3:",
        "CC 19 (37.2%), CT 21 (41.2%), TT 11 (21.6%); T-allele frequency 41.1%, P = 0.27).",
        "Liang 2024 parameterises this as the paper's COVR_ABA column multiplied by a",
        "GENOTYPE-SPECIFIC coefficient theta_ABA (Table 4 rows 'ABCBA 3435 C > T = CT' and",
        "'= TT'), so the log-scale contribution to V/F is 0 for CC, 1 x 0.0151 for CT and",
        "2 x 0.335 for TT. The allele-count reading of COVR_ABA is not merely stated but is",
        "forced by the paper's own numbers: only 0/1/2 coding reproduces all six published",
        "genotype x diet V/F values (Liang 2024 'Safety study'), and in particular the TT",
        "fasted value of 98.96 L requires 2 x 0.335 rather than a single 0.335.",
        "The effect is markedly non-additive per allele (CT is essentially null at +1.5%",
        "while TT is +95%), which is why the model carries two coefficients rather than one",
        "per-allele slope."
      ),
      source_name        = "COVR_ABA (ABCB1 3435 C > T)"
    ),
    FED = list(
      description        = "Fed-state indicator at dosing",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted medication)",
      notes              = paste(
        "1 = lenalidomide taken within 1 h after a meal ('postprandial medication'),",
        "0 = fasted medication. This is Liang 2024's 'diet' covariate, whose operational",
        "definition is given in the Results final-model equation gloss: 'CHFY represent[s]",
        "whether or not taking lenalidomide within 1 h after a meal'. It is a general",
        "fed-vs-fasted flag rather than a high-fat-meal challenge, so the canonical FED",
        "applies rather than FED_HIGHFAT (same rationale as Chen_2023_nemonoxacin.R).",
        "Cohort split (Liang 2024 Table 2): postprandial 28 (54.9%), fasted 23 (45.1%).",
        "Liang 2024 recorded diet as a per-patient habit over the sampled treatment cycle",
        "rather than per dose record, but the covariate is encoded on the canonical",
        "per-dose-record FED column so that a mixed-habit subject can be simulated."
      ),
      source_name        = "COVR_CHFY (diet)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened by forward inclusion / backward elimination; not retained (Liang 2024 Methods 'PPK model' and Results 'PPK model establishment'). Cohort 65.2 +/- 10.8 years, median 67.0 (34.0-85.0) (Table 1)."
    ),
    SEXF = list(
      description = "Sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened; not retained. Cohort 24 women (47.1%) and 27 men (52.9%) (Liang 2024 Table 2)."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened; not retained. Cohort 166 +/- 8.48 cm, median 167 (150-181) (Liang 2024 Table 1)."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened; not retained. Cohort 67.5 +/- 11.0 kg, median 70.0 (40.0-90.0) (Liang 2024 Table 1). Note that Connarn 2018 (Liang 2024 reference 24) did retain body weight for lenalidomide, so the null result here is a genuine disagreement with the prior literature that Liang 2024 discusses explicitly."
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = "Screened; not retained. Cohort 1.79 +/- 0.180 m^2, median 1.79 (1.35-2.17) (Liang 2024 Table 1). Guglieri-Lopez 2017 (Liang 2024 reference 23) retained BSA on V/F for lenalidomide; Liang 2024 did not reproduce it."
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened; not retained. Cohort 91.1 +/- 33.2 mL/min, median 87.4 (33.3-191) (Liang 2024 Table 1). This is the most consequential null result in the paper: lenalidomide is predominantly renally excreted and every prior lenalidomide popPK study cited by Liang 2024 (references 22-25) retained creatinine clearance on CL/F. Liang 2024 found NO covariate for CL/F at all and attributes this to its opportunistic sparse sampling and small sample size."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened; not retained. Cohort 67.9 +/- 24.2 umol/L, median 63.0 (30.0-140) (Liang 2024 Table 1)."
    ),
    CYSC = list(
      description = "Cystatin C",
      units       = "mg/L",
      type        = "continuous",
      notes       = "Screened; not retained. Cohort 1.17 +/- 0.449 mg/L, median 1.07 (0.470-3.09) (Liang 2024 Table 1)."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened; not retained. Cohort 5.44 +/- 2.14 mmol/L, median 5.00 (2.50-13.6) (Liang 2024 Table 1)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened; not retained. Cohort 23.6 +/- 19.9 U/L, median 17.7 (4.80-98.5) (Liang 2024 Table 1)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened; not retained. Cohort 27.8 +/- 28.8 U/L, median 17.9 (4.50-155) (Liang 2024 Table 1)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened; not retained. Cohort 11.1 +/- 5.62 umol/L, median 9.90 (4.60-31.2) (Liang 2024 Table 1)."
    ),
    TPRO = list(
      description = "Total serum protein",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened; not retained. Cohort 62.9 +/- 11.2 g/L, median 60.1 (43.7-97.6) (Liang 2024 Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened; not retained. Cohort 35.9 +/- 5.32 g/L, median 36.9 (23.4-47.8) (Liang 2024 Table 1)."
    ),
    HGB = list(
      description = "Hemoglobin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened; not retained. Cohort 105 +/- 19.6 g/L, median 106 (57.0-153) (Liang 2024 Table 1)."
    ),
    B2M = list(
      description = "Serum beta-2-microglobulin",
      units       = "mg/L",
      type        = "continuous",
      notes       = "Screened; not retained. Liang 2024 Table 1 reports 2.14 +/- 0.183 mg/L with median 2.96 (1.07-10.8); the printed mean lies below the printed median and the SD is implausibly small for that range, so the Table 1 beta-2-microglobulin row appears to be mis-typeset. This does not affect the model, which does not use the covariate."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 51,
    n_observations   = 87,
    n_studies        = 1,
    age_range        = "34-85 years",
    age_median       = "67 years",
    weight_range     = "40-90 kg",
    weight_median    = "70 kg",
    sex_female_pct   = 47.1,
    race_ethnicity   = c(Asian = 100),
    disease_state    = "Hematologic malignancy: multiple myeloma 33 (64.7%), lymphoma 17 (33.3%), myelodysplastic syndrome 1 (2.0%)",
    dose_range       = "10 mg or 25 mg orally once daily",
    regions          = "China (Jinan, Shandong)",
    co_medication    = "Lenalidomide-containing regimens: VRD (bortezomib-lenalidomide-dexamethasone) 27 (52.9%), RR (rituximab-lenalidomide) 12 (23.5%), other 12 (23.5%). Frequent concomitant drugs: aspirin 28 (54.9%), omeprazole 27 (52.9%), acyclovir 23 (45.1%), glutathione 16 (31.4%), tiopronin 14 (27.4%), mecobalamin 11 (21.6%), P-gp substrates 7 (13.7%), other renally excreted drugs 38 (74.5%)",
    notes            = paste(
      "Prospective open-label opportunistic-sampling PPK study (ClinicalTrials.gov NCT06069024)",
      "at the First Affiliated Hospital of Shandong First Medical University, October 2021 to",
      "June 2023. All 51 patients were of Han ethnicity. Samples were drawn after at least five",
      "consecutive daily doses, i.e. at steady state, giving 87 concentrations across 51",
      "subjects (about 1.7 samples per subject) -- the sparseness that Liang 2024 cites as the",
      "reason Ka could not be estimated and no CL/F covariate could be identified.",
      "Baseline demographics: Liang 2024 Table 1 (continuous), Table 2 (categorical),",
      "Table 3 (genotypes). A safety analysis in the 39 patients with available safety data",
      "(Liang 2024 Tables 5 and 6) is reported alongside the PK model but is a non-model",
      "statistical comparison (Mann-Whitney U / Fisher exact tests in SPSS 26.0), not a",
      "fitted exposure-response model, so it is not encoded here."
    )
  )

  ini({
    # Structural parameters -- Liang 2024 Table 4 (final model estimates).
    lka <- fixed(log(6.55)); label("Absorption rate constant (1/h)")            # Liang 2024 Table 4 'theta Ka (h-1) = 6.55 (fixed)'; the Results state Ka "could not be robustly estimated" from the sparse opportunistic data and was fixed to 6.55 1/h from the earlier lenalidomide PPK studies (Liang 2024 references 23-25)
    lcl <- log(7.25);        label("Apparent oral clearance CL/F (L/h)")        # Liang 2024 Table 4 'theta CL/F (L/h) = 7.25' (RSE 7.10%; bootstrap median 7.25, 5-95% CI 6.23-8.56)
    lvc <- log(29.1);        label("Apparent oral volume of distribution V/F, scale term (L)")  # Liang 2024 Table 4 'theta V/F (L) = 29.1' (RSE 9.90%; bootstrap median 28.4, 5-95% CI 17.5-57.2). NOTE this is only the SCALE of the V/F equation, not the reference-subject V/F: the printed final-model equation multiplies it by exp(0.554) (see vcLogOffset in model()), so the reference CC / fasted subject has V/F = 50.64 L

    # Covariate effects on V/F -- all on the natural-log scale, entering the
    # printed final-model equation
    #   V/F_i (L) = 29.1 * Exp(0.554 + theta_ABA * COVR_ABA_i + theta_CHFY * COVR_CHFY_i)
    # (Liang 2024 Results, "PPK model establishment").
    e_snp_abcb1_rs1045642_ct_vc <- 0.0151; label("Per-allele log-scale effect of ABCB1 3435C>T on V/F in CT heterozygotes")   # Liang 2024 Table 4 'ABCBA 3435 C > T = CT = 0.0151' (RSE 42.2%; bootstrap median 0.0151, 5-95% CI 0.0151-0.140). Applied as coefficient x allele count = 0.0151 x 1
    e_snp_abcb1_rs1045642_tt_vc <- 0.335;  label("Per-allele log-scale effect of ABCB1 3435C>T on V/F in TT homozygotes")     # Liang 2024 Table 4 'ABCBA 3435 C > T = TT = 0.335' (RSE 63.9%; bootstrap median 0.350, 5-95% CI 0.0151-1.50). Applied as coefficient x allele count = 0.335 x 2 = 0.670, which is what reproduces the paper's published TT V/F values (98.96 L fasted, 172.9 L fed)
    e_fed_vc                    <- 0.558;  label("Log-scale effect of postprandial dosing on V/F")                            # Liang 2024 Table 4 'theta CHFY = 0.558' (RSE 24.4%; bootstrap median 0.569, 5-95% CI -0.980-1.18); exp(0.558) = 1.75, i.e. V/F is 75% higher when lenalidomide is taken within 1 h after a meal

    # Between-subject variability. Liang 2024 Methods specify the exponential
    # form theta_i = theta_pop * Exp(eta_i) with eta ~ N(0, omega^2), and
    # Table 4 reports eta as a percentage without defining the percentage
    # convention. The percentages are read as omega itself (omega^2 = pct^2),
    # the usual NONMEM reporting for an exponential BSV and the convention used
    # elsewhere in this package (e.g. Gu_2025_rivaroxaban.R). Under the
    # alternative exact log-normal CV reading, omega would be 0.345 and 0.517
    # instead -- 3% and 7% smaller.
    etalcl ~ 0.126025  # Liang 2024 Table 4 'eta CL/F (%) = 35.5' (RSE 17.1%; bootstrap median 31.0, 5-95% CI 0.300-43.5); omega^2 = 0.355^2
    etalvc ~ 0.306916  # Liang 2024 Table 4 'eta V/F (%) = 55.4' (RSE 17.0%; bootstrap median 53.4, 5-95% CI 28.3-96.8); omega^2 = 0.554^2

    # Residual error -- additive. Liang 2024 tested additive, proportional,
    # power-exponential and combined residual models and selected the additive
    # one (Results, "PPK model establishment").
    addSd <- 38.5; label("Additive residual error (ng/mL)")  # Liang 2024 Table 4 'epsilon (ng/mL) = 38.5' (RSE 21.6%; bootstrap median 39.9, 5-95% CI 20.8-62.9)
  })

  model({
    # Constant printed in the Liang 2024 final-model V/F equation
    #   V/F_i (L) = 29.1 * Exp(0.554 + theta_ABA * COVR_ABA_i + theta_CHFY * COVR_CHFY_i)
    # but absent from Table 4. It is load-bearing and is retained verbatim:
    # 29.1 * exp(0.554) = 50.64 L, which is exactly the V/F the paper reports
    # for its reference subject (CC genotype, fasted) in the "Safety study"
    # section, and all six of the paper's published genotype x diet V/F values
    # reproduce to four significant figures only when it is included. It also
    # gives a physiologically correct terminal half-life of
    # log(2) * 50.64 / 7.25 = 4.8 h, inside the 3-5 h reported for
    # lenalidomide (Liang 2024 Introduction); dropping it would give 2.8 h.
    vcLogOffset <- 0.554

    # ABCB1 3435C>T effect: the paper's genotype-specific coefficient
    # theta_ABA multiplied by the T-allele count COVR_ABA (0 / 1 / 2), so CC
    # contributes 0, CT contributes 1 x 0.0151 and TT contributes 2 x 0.335.
    abcb1LogEffect <-
      e_snp_abcb1_rs1045642_ct_vc * SNP_ABCB1_RS1045642_T_COUNT * (SNP_ABCB1_RS1045642_T_COUNT == 1) +
      e_snp_abcb1_rs1045642_tt_vc * SNP_ABCB1_RS1045642_T_COUNT * (SNP_ABCB1_RS1045642_T_COUNT == 2)

    # Individual parameters. Reference subject: CC genotype
    # (SNP_ABCB1_RS1045642_T_COUNT = 0), fasted (FED = 0), for whom
    # V/F = 50.64 L and CL/F = 7.25 L/h.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) * exp(vcLogOffset + abcb1LogEffect + e_fed_vc * FED)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Dose is in mg and vc is in L, so central / vc is mg/L; x 1000 converts to
    # the ng/mL used for the concentrations and for the additive residual error
    # throughout Liang 2024.
    Cc <- central / vc * 1000
    Cc ~ add(addSd)
  })
}
