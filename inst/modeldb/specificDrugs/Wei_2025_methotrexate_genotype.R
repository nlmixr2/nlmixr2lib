Wei_2025_methotrexate_genotype <- function() {
  description <- "Three-compartment population PK model for high-dose intravenous methotrexate in Chinese adults with primary central nervous system lymphoma (Wei 2025), GENE variant. Identical in structure to the companion Wei_2025_methotrexate.R but adds a composite ABCC4-ABCG2-ADORA2A genotype effect on clearance, derived inside the model from the variant-allele counts of ABCC4 rs2274407, ABCG2 rs2231142 and ADORA2A rs2298383. Clearance also carries power effects of BSA-normalized eGFR, blood urea nitrogen and alanine aminotransferase; the inter-compartmental clearance to the first peripheral compartment carries a power effect of total serum protein. Parameter values are taken from the publication's Table 5 ('Final gene-model' column) and the covariate equations 14 to 19. The authors fitted both models because genotyping is not routinely available in clinical practice; use the companion nongene file when genotype data are absent."
  reference <- paste(
    "Wei S, Zhang S, Wang D, Zhang D, Lu Q, Mo J, Yang Z, Guan L, He Y,",
    "Zhao Z, Mei S. (2025). Population pharmacokinetics of high-dose",
    "methotrexate in patients with primary central nervous system lymphoma.",
    "Front Pharmacol 16:1578033.",
    "doi:10.3389/fphar.2025.1578033.",
    sep = " "
  )
  vignette <- "Wei_2025_methotrexate"
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  covariateData <- list(
    CRCL = list(
      description        = "Estimated glomerular filtration rate, BSA-normalized",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Computed by the 2021 race-free CKD-EPI creatinine equation (Methods 'Study design', Equation 3, citing Inker 2021), NOT by Cockcroft-Gault. Normalized to 101.8 mL/min/1.73 m^2 in Equation 14, matching the Table 3 cohort median exactly. Cohort range 5.4-162.9 mL/min/1.73 m^2. Time-varying, monitored daily for at least three days after each administration. Enters clearance only, as a power term with the estimated exponent `e_crcl_cl` = 0.67 -- numerically identical to the nongene model.",
      source_name        = "eGFR"
    ),
    BUN = list(
      description        = "Blood urea nitrogen concentration",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "SI units (mmol/L), NOT mg/dL -- the Table 3 median of 4.6 with range 0.5-19 is the mmol/L scale. Normalized to 4.6 mmol/L in Equation 14, matching the Table 3 median exactly. Enters clearance only, as a power term with the estimated exponent `e_bun_cl` = -0.08, identical to the nongene model.",
      source_name        = "BUN"
    ),
    ALT = list(
      description        = "Serum alanine aminotransferase activity",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Normalized to 25 U/L in Equation 14, matching the Table 3 cohort median exactly. Cohort range 2.2-1141.7 U/L. Enters clearance only, as a power term with the estimated exponent `e_alt_cl` = +0.03, identical to the nongene model. The sign is POSITIVE: liver injury raises methotrexate clearance in this cohort, which the authors flag as counterintuitive (Discussion).",
      source_name        = "ALT"
    ),
    TPRO = list(
      description        = "Total serum protein concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters the inter-compartmental clearance Q1 only, as a power term with the estimated exponent `e_tpro_q` = -1.72. This is the ONE covariate exponent that differs between the two final models (-1.72 here against -1.68 in the nongene model). NOTE the same paper-internal mismatch recorded in the companion file, transcribed as printed: Equation 15 normalizes to 58 g/L while the Table 3 cohort median is 61.8 g/L, despite Methods 'Covariate model' stating that continuous covariates were standardized to their medians. The printed equation constant 58 is used.",
      source_name        = "TP"
    ),
    SNP_ABCC4_RS2274407_G_COUNT = list(
      description        = "ABCC4 rs2274407 (T>G) variant G-allele count",
      units              = "(count, 0/1/2 alleles per subject)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Germline genotype, time-invariant. 0 = TT wild-type homozygote, 1 = TG heterozygote, 2 = GG variant homozygote. ABCC4 encodes multidrug-resistance protein 4 (MRP4). One of the three variants combined into the composite genotype covariate; the count form is required because the paper's combination rule sums per-variant group scores rather than using carrier indicators. Direction: the Discussion states 'the T allele was associated with increased MTX clearance', so the G variant allele DECREASES clearance, placing this variant on Table 1's 'Decreased' scoring row where the group score is 1 + (variant allele count). Note the paper labels the SNP '(T > G)' in Results 3.2 but discusses it as 'G912T' with 'T allele carriers' in the Discussion, following the cited Mesrian Tanha 2017 nomenclature; the two are opposite strand/orientation conventions for the same variant. This column follows the dbSNP-consistent orientation the paper uses when it defines the variant.",
      source_name        = "ABCC4 rs2274407"
    ),
    SNP_ABCG2_RS2231142_T_COUNT = list(
      description        = "ABCG2 rs2231142 (G>T, Q141K) variant T-allele count",
      units              = "(count, 0/1/2 alleles per subject)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Germline genotype, time-invariant. 0 = GG wild-type homozygote, 1 = GT heterozygote, 2 = TT variant homozygote. ABCG2 encodes breast cancer resistance protein (BCRP). This is the variant the paper uses as its worked example of the grouping rule (Methods 'Grouping and combination of variants': 'the ABCG2 rs2231142 G > T variant, with the T allele, was associated with decreased MTX clearance. As a result, the genotypes GG, GT, and TT were categorized into groups 1, 2, and 3'), which fixes the scoring convention for the whole composite. Direction: the T variant allele DECREASES clearance (Discussion: 'our research also links the T allele to reduced MTX CL').",
      source_name        = "ABCG2 rs2231142"
    ),
    SNP_ADORA2A_RS2298383_T_COUNT = list(
      description        = "ADORA2A rs2298383 (C>T) variant T-allele count",
      units              = "(count, 0/1/2 alleles per subject)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Germline genotype, time-invariant. 0 = CC wild-type homozygote, 1 = CT heterozygote, 2 = TT variant homozygote. ADORA2A encodes the adenosine A2A receptor; rs2298383 sits in a putative promoter region and is associated with transcriptional regulation. Direction: the Discussion states 'the C allele was associated with increased MTX clearance', so the T variant allele DECREASES clearance, placing this variant on Table 1's 'Decreased' scoring row alongside the other two.",
      source_name        = "ADORA2A rs2298383"
    )
  )

  compartmentData <- list(
    central     = list(analyte = "methotrexate", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "methotrexate", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "methotrexate", units = "umol", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 752L,
    n_studies      = 1L,
    age_range      = "18.12-86.65 years (median 57.445)",
    age_median     = "57.445 years",
    weight_range   = "30-115 kg (median 68)",
    weight_median  = "68 kg",
    bsa_range      = "1.16-2.39 m^2 (median 1.73)",
    bsa_median     = "1.73 m^2",
    sex_female_pct = 44.1,
    renal_function = "eGFR 5.4-162.9 mL/min/1.73 m^2 (median 101.8) by the 2021 CKD-EPI equation; serum creatinine 24.8-641.7 umol/L (median 64.6); Cockcroft-Gault CLcr 5.9-361.8 mL/min (median 100.1). At least 17.4% of the cohort met the label definition of delayed elimination.",
    hepatic_function = "ALT 2.2-1141.7 U/L (median 25); AST 5-1915.2 U/L (median 20.4); total protein 27.4-95.7 g/L (median 61.8); albumin 19.9-51.8 g/L (median 37.5).",
    disease_state  = "Adults with primary central nervous system lymphoma (PCNSL) receiving high-dose methotrexate, most commonly combined with rituximab or cytarabine.",
    dose_range     = "Intravenous methotrexate 3.5 g/m^2, median infusion duration 3.1 h. Median of four infusions per patient (range 1-34).",
    genotyping     = "Twenty-nine single nucleotide polymorphisms with a minor allele frequency above 0.05 in the Chinese population were genotyped by MassARRAY, spanning MTHFR, MTR, ATIC, ABCG2, MTRR, ABCB1, ABCC2, ABCC4, MTHFD1, SLCO1B1, SLC28A2, TYMS and SLC19A1 (listed in Supplementary Appendix SA1). All but rs10760502, rs11045879, rs2413775 and rs3758149 were in Hardy-Weinberg equilibrium. Only the three-variant ABCC4-ABCG2-ADORA2A composite reached significance; SLCO1B1, ABCC2, ABCB1 and MTHFR did not.",
    regions        = "China (single center: Beijing Tiantan Hospital, Capital Medical University), September 2016 through August 2023.",
    notes          = "Retrospective therapeutic-drug-monitoring cohort of 752 adults contributing 6074 methotrexate plasma concentrations. Each methotrexate administration was treated as an INDEPENDENT event in the dataset because dosing intervals exceeded five elimination half-lives (Methods 'Base model'). Concentrations were total drug by UHPLC-MS/MS, lower limit of quantification 0.002 umol/L. Estimation was by first-order conditional estimation extended least squares in Phoenix NLME 8.3, with 200 bootstrap replicates and a 1000-replicate visual predictive check. Adding the composite genotype to the nongene model dropped the objective function by 9.95 units (Table 4, models 5 to 6), the smallest of the five covariate steps. Demographics from Table 3; parameter estimates from Table 5 ('Final gene-model' column)."
  )

  ini({
    # Structural PK parameters -- Wei 2025 Table 5, 'Final gene-model' column,
    # cross-checked against the closed forms printed as Equations 14 to 19.
    # Typical values at the reference covariates eGFR = 101.8 mL/min/1.73 m^2,
    # BUN = 4.6 mmol/L, ALT = 25 U/L, TP = 58 g/L and NON-carrier genotype.
    # Only CL, Vc, Vp1 and theta_TP differ from the nongene model, and only in
    # the third significant figure apart from CL.
    lcl  <- log(8.45)  ; label("Clearance CL at the reference covariates, non-carrier (L/h)")                     # Table 5 gene CL = 8.45 (%RSE 2.98) [7.95, 8.94]; also Equation 14
    lvc  <- log(33.29) ; label("Central volume of distribution Vc (L)")                                           # Table 5 gene Vc = 33.29 (%RSE 3.50) [31.00, 35.57]; Equation 17 rounds this to 33.3
    lq   <- log(0.04)  ; label("Inter-compartmental clearance Q1 to peripheral1 at TP = 58 g/L (L/h)")            # Table 5 gene Q1 = 0.04 (%RSE 8.14) [0.03, 0.05]; also Equation 15
    lvp  <- log(17.85) ; label("First peripheral volume of distribution Vp1 (L)")                                 # Table 5 gene Vp1 = 17.85 (%RSE 11.48) [13.84, 21.87]; Equation 18 rounds this to 17.9
    lq2  <- log(0.09)  ; label("Inter-compartmental clearance Q2 to peripheral2 (L/h)")                           # Table 5 gene Q2 = 0.09 (%RSE 5.14) [0.08, 0.10]; also Equation 16
    lvp2 <- log(1.14)  ; label("Second peripheral volume of distribution Vp2 (L)")                                # Table 5 gene Vp2 = 1.14 (%RSE 4.36) [1.04, 1.24]; also Equation 19

    # Covariate effects. Equation 14 is
    #   CL (L/h) = 8.45 * (eGFR/101.8)^0.67 * (BUN/4.6)^-0.08 * (ALT/25)^0.03
    #              * 0.91 (if ABCC-ABCG-ADORA2A mutation)
    # and Equation 15 is Q1 (L/h) = 0.04 * (TP/58)^-1.72.
    e_crcl_cl <-  0.67 ; label("Power exponent on (CRCL / 101.8 mL/min/1.73 m^2) for CL (unitless)")              # Table 5 gene theta_eGFR = 0.67 (%RSE 2.10) [0.64, 0.70]
    e_bun_cl  <- -0.08 ; label("Power exponent on (BUN / 4.6 mmol/L) for CL (unitless)")                          # Table 5 gene theta_BUN = -0.08 (%RSE 10.20) [-0.09, -0.06]
    e_alt_cl  <-  0.03 ; label("Power exponent on (ALT / 25 U/L) for CL (unitless)")                              # Table 5 gene theta_ALT = 0.03 (%RSE 14.23) [0.02, 0.03]
    e_tpro_q  <- -1.72 ; label("Power exponent on (TPRO / 58 g/L) for Q1 (unitless)")                             # Table 5 gene theta_TP = -1.72 (%RSE 8.12) [-1.99, -1.44]

    # Composite-genotype effect on clearance, encoded as the DIRECT
    # MULTIPLIER 0.91 that Equation 14 and the Abstract both print
    # ("a = 0.91 for gene-model if ABCC-ABCG-ADORA2A mutation, otherwise
    # a = 1"). Table 5 reports the underlying coefficient as
    # theta_ABCC4-ABCG2-ADORA2A = -0.09 with CI [-0.14, -0.03]. The printed
    # 0.91 equals 1 + (-0.09) EXACTLY and equals exp(-0.09) = 0.9139 only
    # after rounding, and the Discussion quotes 'an approximately 9%
    # reduction in MTX clearance', so the multiplier is taken as printed
    # rather than back-transformed. The two readings differ by 0.4% in CL.
    e_gene_cl <- 0.91 ; label("Multiplicative effect of the ABCC4-ABCG2-ADORA2A composite mutation on CL (unitless, non-carrier = reference)") # Table 5 gene theta_ABCC4-ABCG2-ADORA2A = -0.09 (%RSE -30.63) [-0.14, -0.03]; Equation 14 prints the multiplier as 0.91

    # Inter-individual variability -- Equation 10, theta_i = theta_TV *
    # exp(eta), eta ~ N(0, omega^2). Table 5 heads these rows 'IIV<param>
    # (CV%)' and the tabulated numbers are read as 100 * omega. See the
    # companion file's ini() comment and the vignette Errata for the
    # supporting %RSE argument and the one residual ambiguity. There is NO
    # inter-individual variability on Q2.
    etalcl  ~ 0.073441  # Table 5 gene IIV_CL  = 27.1 CV%  (%RSE 1.88);  variance = 0.271^2
    etalvc  ~ 0.041087  # Table 5 gene IIV_Vc  = 20.27 CV% (%RSE 5.91);  variance = 0.2027^2
    etalq   ~ 0.980298  # Table 5 gene IIV_Q1  = 99.01 CV% (%RSE 9.47);  variance = 0.9901^2
    etalvp  ~ 0.609961  # Table 5 gene IIV_Vp1 = 78.1 CV%  (%RSE 25.46); variance = 0.781^2
    etalvp2 ~ 0.073984  # Table 5 gene IIV_Vp2 = 27.2 CV%  (%RSE 5.23);  variance = 0.272^2

    # Residual error -- Equation 11, Cobs = Cpred * (1 + epsilon). Table 5
    # reports the identical sigma for both final models.
    propSd <- 0.7379 ; label("Proportional residual error (fraction)")                                            # Table 5 sigma (proportional) = 73.79 (%RSE 1.07) [72.24, 75.34]
  })

  model({
    # 1. Derive the composite ABCC4-ABCG2-ADORA2A genotype indicator from the
    #    three variant-allele counts. The paper builds it in two documented
    #    steps. Table 1 scores each variant into three groups running WITH
    #    decreasing clearance; all three of these variants have the variant
    #    allele associated with DECREASED clearance (see each covariateData
    #    note), so each variant's group score is 1 + (variant allele count)
    #    and the three-variant sum runs 3 to 9. Table 2 then bins that sum
    #    into two groups. Results 3.2 fixes which of Table 2's three binning
    #    rules applies: 'Patients were identified as mutation carriers ... if
    #    they exhibited more than three nucleotide mutations among these
    #    variants', i.e. more than 3 of the 6 possible variant alleles, i.e.
    #    a group sum of at least 7 -- which is exactly Table 2's 'Two groups,
    #    rule 2' row for 3 combined variants (3-6 -> group 1, 7-9 -> group 2).
    #    Working directly in variant-allele counts, the carrier condition is
    #    therefore simply a total count above 3.
    gene_allele_count <- SNP_ABCC4_RS2274407_G_COUNT +
      SNP_ABCG2_RS2231142_T_COUNT +
      SNP_ADORA2A_RS2298383_T_COUNT
    gene_carrier <- (gene_allele_count > 3)

    # 2. Individual PK parameters. Clearance carries the three power-scaled
    #    covariate terms plus the composite-genotype multiplier of
    #    Equation 14; the first inter-compartmental clearance carries the
    #    total-protein term of Equation 15. The multiplier is applied as
    #    e_gene_cl^gene_carrier so that non-carriers (indicator 0) get a
    #    factor of exactly 1, matching the 'otherwise a = 1' of the Abstract.
    cl  <- exp(lcl + etalcl) * (CRCL / 101.8)^e_crcl_cl *
      (BUN / 4.6)^e_bun_cl * (ALT / 25)^e_alt_cl *
      e_gene_cl^gene_carrier
    vc  <- exp(lvc + etalvc)
    q   <- exp(lq + etalq) * (TPRO / 58)^e_tpro_q
    vp  <- exp(lvp + etalvp)
    q2  <- exp(lq2)
    vp2 <- exp(lvp2 + etalvp2)

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # 4. Three-compartment intravenous disposition, Equations 4 to 9.
    d/dt(central)     <- -kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # 5. Observation. Dose units umol and vc units L give Cc in umol/L.
    #    Methotrexate has a molar mass of 454.44 g/mol, so 1 mg = 2.2005 umol.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
