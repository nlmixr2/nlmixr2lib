Papachristos_2020_bevacizumab_pkpd <- function() {
  description <- "Two-compartment population PK plus immediate-response Imax PK/PD model for IV bevacizumab and free VEGF-A in adults with metastatic colorectal cancer, with allometric weight scaling and ICAM-1 / VEGF-A genotype covariates (Papachristos 2020, Table 3)"
  reference <- "Papachristos A, Karatza E, Kalofonos H, Karalis V. Pharmacogenetics in Model-Based Optimization of Bevacizumab Therapy for Metastatic Colorectal Cancer. Int J Mol Sci. 2020;21(11):3753. doi:10.3390/ijms21113753"
  vignette <- "Papachristos_2020_bevacizumab_pkpd"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L (bevacizumab) and ng/L (free VEGF-A)")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed in source analysis; applied as power-form covariate on CL with reference weight 70 kg.",
      source_name        = "weight"
    ),
    SNP_ICAM1_RS1799969 = list(
      description        = "ICAM-1 rs1799969 mutant-allele-presence indicator (1 = heterozygous or homozygous mutant; 0 = homozygous wild-type)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (homozygous wild-type)",
      notes              = "Time-fixed per subject. Multiplicative effect on CL.",
      source_name        = "cat (Papachristos 2020 Table 3 narrative)"
    ),
    SNP_VEGFA_RS699947 = list(
      description        = "VEGF-A rs699947 mutant-allele-presence indicator (1 = heterozygous or homozygous mutant; 0 = homozygous wild-type)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (homozygous wild-type)",
      notes              = "Time-fixed per subject. Multiplicative effect on Q.",
      source_name        = "cat (Papachristos 2020 Table 3 narrative)"
    )
  )

  population <- list(
    n_subjects        = 46,
    n_studies         = 1,
    age_range         = "IQR 53-72 years",
    age_median        = "63 years",
    weight_range      = "IQR 64.15-81.75 kg",
    weight_median     = "74.5 kg",
    sex_female_pct    = 39,
    race_ethnicity    = "Greek (single-centre study); race not formally reported.",
    disease_state     = "Adults with metastatic colorectal cancer (mCRC); ECOG performance status 0-2.",
    dose_range        = "5 mg/kg IV every 2 weeks (76% of patients) or 7.5 mg/kg IV every 3 weeks (24%) with chemotherapy.",
    regions           = "Greece",
    n_observations    = "156 bevacizumab + 169 free VEGF-A serum concentrations",
    co_medication     = "BEV-FOLFIRI, BEV-FOLFOX, BEV-CapIRI, or BEV-CapOX.",
    notes             = "Single-centre prospective observational study (Papachristos 2020 sec. 2.1, 4.1)."
  )

  ini({
    # Structural PK parameters (reference 70 kg adult, all SNPs wild-type)
    lcl <- log(0.388); label("Clearance for a 70 kg wild-type adult (CL, L/day)")  # Table 3 row CLpop
    lv1 <- log(5.48);  label("Central volume of distribution (V1, L)")              # Table 3 row V1pop
    lq  <- log(0.315); label("Inter-compartmental clearance (Q, L/day)")            # Table 3 row Qpop
    lv2 <- log(8.81);  label("Peripheral volume of distribution (V2, L)")           # Table 3 row V2(L)

    # PD parameters (immediate-response Imax model on free VEGF-A)
    le0    <- log(684);   label("Baseline free VEGF-A concentration for a wild-type adult (E0, ng/L)")  # Table 3 row E0pop
    imaxv  <- 0.951;      label("Maximal fractional inhibition of free VEGF-A (Imax, fraction)")        # Table 3 row Imaxpop
    lic50  <- log(29.1);  label("Bevacizumab concentration producing 50% Imax (IC50, mg/L)")           # Table 3 row IC50pop

    # Allometric / SNP covariate effects on log-CL
    allo_cl <- 0.78;  label("Allometric exponent on CL for log(WT/70) (unitless)")  # Table 3 row "log(weight/70) on CL"
    e_icam1_rs1799969_cl <- -0.423; label("ICAM-1 rs1799969 mutant effect on log-CL (unitless)")  # Table 3 row "ICAM-1 rs1799969 mutant on CL"

    # SNP effect on Q (Table 3 value; the section-2.4 text gives "exp(0.378)" but that is a
    # copy-paste typo from the PK model — Table 3 reports -0.414 and the narrative confirms
    # the direction: "Inter-compartmental clearance (Q) was lower in patients with mutant
    # type gene VEGF-A rs699947".)
    e_vegfa_rs699947_q <- -0.414; label("VEGF-A rs699947 mutant effect on log-Q (unitless)")  # Table 3 row "VEGF-A rs699947 mutant on Q"

    # IIV (exponential model; correlation rho(CL,Q) = -0.979 -> cov = -0.979 * 0.338 * 0.601 = -0.19891)
    etalcl + etalq ~ c(0.11424,
                       -0.19891, 0.36120)  # Table 3 omega_CL = 0.338, omega_Q = 0.601, p(Q,CL) = -0.979
    etalv1 ~ 0.03098  # Table 3 omega_V1 = 0.176
    etalv2 ~ 0.33524  # Table 3 omega_V2 = 0.579
    etale0 ~ 0.02789  # Table 3 omega_E0 = 0.167

    # Residual error (proportional, two outputs)
    CcpropSd <- 0.238; label("Proportional residual error for bevacizumab (fraction)")   # Table 3 row sigma_BEVA
    EpropSd  <- 0.264; label("Proportional residual error for free VEGF-A (fraction)")   # Table 3 row sigma_VEGF
  })
  model({
    # Individual PK parameters (paper formulas, sec. 2.4)
    cl  <- exp(lcl + allo_cl * log(WT / 70) + e_icam1_rs1799969_cl * SNP_ICAM1_RS1799969 + etalcl)
    v1  <- exp(lv1 + etalv1)
    q   <- exp(lq + e_vegfa_rs699947_q * SNP_VEGFA_RS699947 + etalq)
    v2  <- exp(lv2 + etalv2)

    # Individual PD parameters
    e0  <- exp(le0 + etale0)
    ic50 <- exp(lic50)
    imax <- imaxv

    # Two-compartment linear PK with IV-infusion dosing into central
    kel <- cl / v1
    k12 <- q  / v1
    k21 <- q  / v2

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Bevacizumab observation: dose in mg, V1 in L -> mg/L
    Cc <- central / v1
    Cc ~ prop(CcpropSd)

    # Free VEGF-A observation: immediate-response Imax inhibition (no sigmoidicity)
    # E = E0 * (1 - Imax * Cc / (IC50 + Cc))
    E <- e0 * (1 - imax * Cc / (ic50 + Cc))
    E ~ prop(EpropSd)
  })
}
