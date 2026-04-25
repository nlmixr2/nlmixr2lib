Papachristos_2020_bevacizumab_pk <- function() {
  description <- "Two-compartment population PK model for IV bevacizumab in adults with metastatic colorectal cancer, with allometric weight scaling and ICAM-1 / VEGF-A genotype covariates (Papachristos 2020, Table 1)"
  reference <- "Papachristos A, Karatza E, Kalofonos H, Karalis V. Pharmacogenetics in Model-Based Optimization of Bevacizumab Therapy for Metastatic Colorectal Cancer. Int J Mol Sci. 2020;21(11):3753. doi:10.3390/ijms21113753"
  vignette <- "Papachristos_2020_bevacizumab_pk"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed in source analysis (single baseline value used per subject); applied as power-form covariate on CL with reference weight 70 kg.",
      source_name        = "weight"
    ),
    SNP_ICAM1_RS1799969 = list(
      description        = "ICAM-1 rs1799969 mutant-allele-presence indicator (1 = heterozygous or homozygous mutant; 0 = homozygous wild-type)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (homozygous wild-type)",
      notes              = "Time-fixed per subject. Multiplicative effect on CL. Mutant carriers (20% of the cohort) have lower bevacizumab clearance and higher trough concentrations.",
      source_name        = "cat (Papachristos 2020 Table 1 narrative; the paper does not name the data column)"
    ),
    SNP_VEGFA_RS1570360 = list(
      description        = "VEGF-A rs1570360 mutant-allele-presence indicator (1 = heterozygous or homozygous mutant; 0 = homozygous wild-type)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (homozygous wild-type)",
      notes              = "Time-fixed per subject. Multiplicative effect on inter-compartmental clearance Q. Mutant carriers (33% of the cohort) have higher Q.",
      source_name        = "cat1 (Papachristos 2020 Table 1 narrative)"
    ),
    SNP_VEGFA_RS699947 = list(
      description        = "VEGF-A rs699947 mutant-allele-presence indicator (1 = heterozygous or homozygous mutant; 0 = homozygous wild-type)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (homozygous wild-type)",
      notes              = "Time-fixed per subject. Multiplicative effect on Q. Mutant carriers (52% of the cohort) have lower Q.",
      source_name        = "cat2 (Papachristos 2020 Table 1 narrative)"
    )
  )

  population <- list(
    n_subjects        = 46,
    n_studies         = 1,
    age_range         = "IQR 53-72 years (full range not reported)",
    age_median        = "63 years",
    weight_range      = "IQR 64.15-81.75 kg (full range not reported)",
    weight_median     = "74.5 kg",
    sex_female_pct    = 39,
    race_ethnicity    = "Greek (single-centre study at the University Hospital of Patras, Greece); race not formally reported.",
    disease_state     = "Adults with histologically confirmed metastatic colorectal cancer (mCRC); ECOG performance status 0-2.",
    dose_range        = "5 mg/kg IV every 2 weeks (76% of patients) or 7.5 mg/kg IV every 3 weeks (24%) in combination with oxaliplatin/fluoropyrimidine or irinotecan/fluoropyrimidine chemotherapy.",
    regions           = "Greece",
    n_observations    = 156,
    co_medication     = "BEV-FOLFIRI, BEV-FOLFOX, BEV-CapIRI, or BEV-CapOX (chemotherapy partner did not have a statistically significant effect on PK parameters).",
    snp_carrier_rates = c(
      ICAM1_rs1799969_mutant = 20,
      VEGFA_rs1570360_mutant = 33,
      VEGFA_rs699947_mutant  = 52,
      ICAM1_rs5498_mutant    = 70,
      VEGFA_rs2010963_mutant = 41
    ),
    notes             = "Single-centre prospective observational study (Papachristos 2020 section 2.1 and 4.1). The non-significant SNPs (ICAM-1 rs5498 and VEGF-A rs2010963) are listed for completeness; only the three SNPs in covariateData entered the final PK model."
  )

  ini({
    # Structural PK parameters - reference values for a 70 kg adult, all SNPs wild-type
    lcl <- log(0.200); label("Clearance for a 70 kg wild-type adult (CL, L/day)")  # Table 1 row CLpop
    lv1 <- log(3.09);  label("Central volume of distribution (V1, L)")              # Table 1 row V1pop
    lq  <- log(0.35);  label("Inter-compartmental clearance for a wild-type adult (Q, L/day)")  # Table 1 row Qpop
    lv2 <- log(2.39);  label("Peripheral volume of distribution (V2, L)")           # Table 1 row V2pop

    # Allometric / covariate effects on log-CL (paper formula: CL = CLpop * (WT/70)^allo * exp(e * cat))
    allo_cl <- 1.04;  label("Allometric exponent on CL for log(WT/70) (unitless)")  # Table 1 row "log(weight/70) on CL"
    e_icam1_rs1799969_cl <- -0.423; label("ICAM-1 rs1799969 mutant effect on log-CL (unitless)")  # Table 1 row "ICAM-1 rs1799969 mutant on CL"

    # Covariate effects on log-Q
    e_vegfa_rs1570360_q <-  0.378; label("VEGF-A rs1570360 mutant effect on log-Q (unitless)")  # Table 1 row "VEGF-A rs1570360 mutant on Q"
    e_vegfa_rs699947_q  <- -0.429; label("VEGF-A rs699947 mutant effect on log-Q (unitless)")   # Table 1 row "VEGF-A rs699947 mutant on Q"

    # IIV (exponential model: P = Ppop * exp(eta), eta ~ N(0, omega^2))
    # Variances are squares of the SDs reported in Table 1; correlation rho(CL,Q) = -0.999 maps to
    # cov(etalcl, etalq) = rho * omega_cl * omega_q = -0.999 * 0.319 * 0.160 = -0.05099
    etalcl + etalq ~ c(0.10176,
                       -0.05099, 0.02560)  # Table 1 omega_CL=0.319, omega_Q=0.160, p(Q,CL)=-0.999
    etalv1 ~ 0.03028  # Table 1 omega_V1 = 0.174
    etalv2 ~ 0.45698  # Table 1 omega_V2 = 0.676

    # Residual error (proportional)
    propSd <- 0.246; label("Proportional residual error (fraction)")  # Table 1 row sigma_prop
  })
  model({
    # Individual PK parameters
    # log(CL_i) = lcl + allo_cl * log(WT/70) + e_icam * SNP_ICAM1_RS1799969 + etalcl  (Papachristos 2020 sec. 2.2 equation)
    cl <- exp(lcl + allo_cl * log(WT / 70) + e_icam1_rs1799969_cl * SNP_ICAM1_RS1799969 + etalcl)
    v1 <- exp(lv1 + etalv1)
    q  <- exp(lq + e_vegfa_rs1570360_q * SNP_VEGFA_RS1570360 + e_vegfa_rs699947_q * SNP_VEGFA_RS699947 + etalq)
    v2 <- exp(lv2 + etalv2)

    # Two-compartment linear PK with IV-infusion dosing into central
    kel <- cl / v1
    k12 <- q  / v1
    k21 <- q  / v2

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Concentration: dose in mg, V1 in L -> mg/L (= ug/mL)
    Cc <- central / v1
    Cc ~ prop(propSd)
  })
}
