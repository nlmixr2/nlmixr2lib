Duwal_2019_efavirenz <- function() {
  description <- paste(
    "QSP (systems-pharmacology PK-PD with HIV-1 viral dynamics). Efavirenz",
    "(EFV) for HIV-1 pre- and post-exposure prophylaxis, from Duwal et al.",
    "2019. One-compartment oral popPK model of the ENCORE1 trial (ka fixed",
    "0.6 1/h; CL/F 11.9 L/h and V/F 282 L for a 70 kg subject of reference",
    "CYP2B6/CYP2A6 genotype; allometric weight; a 16-level composite",
    "CYP2B6 516G>T / 983T>C / CYP2A6 *9B / *17 genotype multiplier on",
    "CL/F; IIV and IOV on CL/F; proportional residual error) extended by",
    "hyperbolic metabolic autoinduction from the first dose (day-1 CL/F",
    "58% of the steady-state value, half-time 245 h). Plasma EFV drives a",
    "sigmoid Emax inhibition of reverse transcription (protein-binding",
    "corrected ex vivo IC50 and Hill coefficient) for wild-type HIV-1 and",
    "the Y181C, K103N and G190S mutants, which also differ in fitness. For",
    "each strain the model carries the expected-value (mean-field) ODEs of",
    "the authors' three-state viral-replication reactions (free virus,",
    "early and productively infected T cells; uninfected T cells constant)",
    "and the closed-form prophylactic efficacy against a single",
    "transmitted virion that the current EFV concentration would give if",
    "held constant (Figure 2A). The stochastic EXTRANDE simulations behind",
    "Figures 3-4 are not reproduced."
  )
  reference <- paste(
    "Duwal S, Seeler D, Dickinson L, Khoo S, von Kleist M. The utility of",
    "efavirenz-based prophylaxis against HIV infection. A systems",
    "pharmacological analysis. Front Pharmacol. 2019;10:199.",
    "doi:10.3389/fphar.2019.00199. PMCID PMC6424904. PK structure from",
    "Eqs. (1)-(3) and Table 1; autoinduction from Eq. (4) (Zhu et al.,",
    "Antimicrob Agents Chemother. 2009;53:2346-2353,",
    "doi:10.1128/AAC.01120-08); genotype multipliers from Supplementary",
    "Table 2, reproduced from Dickinson et al., Clin Pharmacol Ther.",
    "2015;98:406-416, doi:10.1002/cpt.156; direct effect from Eq. (5),",
    "Table 2 and Supplementary Text 1 Eq. (S1.10); viral dynamics from",
    "Eqs. (6)-(11) and Table 3; prophylactic efficacy from Eq. (13)."
  )
  vignette <- "Duwal_2019_efavirenz"

  # Dose in mg and volume in L, so Cc is in mg/L (the unit of Figure 1). The
  # viral states are counts (virions, cells); an inoculum is dosed into a
  # virus state as a number of virions.
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  dosing <- c("depot", "virus", "virus_y181c", "virus_k103n", "virus_g190s")

  paper_specific_compartments <- c(
    "tcell_early",
    "tcell_infected",
    "virus_y181c",
    "tcell_early_y181c",
    "tcell_infected_y181c",
    "virus_k103n",
    "tcell_early_k103n",
    "tcell_infected_k103n",
    "virus_g190s",
    "tcell_early_g190s",
    "tcell_infected_g190s"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling centred on 70 kg: exponent 0.75 on CL/F and 1 on",
        "V/F (Table 1 legend). ENCORE1 median 65 kg (range 39-148)."
      ),
      source_name = "weight"
    ),
    SNP_CYP2B6_RS3745274_T_COUNT = list(
      description = paste(
        "CYP2B6 516G>T (rs3745274) T-allele count: 0 = GG, 1 = GT, 2 = TT."
      ),
      units = "(count, 0/1/2)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "First position of the composite CYP2B6 516G>T / 983T>C / CYP2A6",
        "*9B / *17 genotype of Supplementary Table 2; each level of the",
        "composite carries its own CL/F multiplier. ENCORE1 frequencies",
        "(Dickinson 2015 Table 1): GG 44%, GT 46%, TT 10%."
      ),
      source_name = "CYP2B6 516G>T"
    ),
    SNP_CYP2B6_RS28399499_C_COUNT = list(
      description = paste(
        "CYP2B6 983T>C (rs28399499) C-allele count: 0 = TT, 1 = TC, 2 = CC."
      ),
      units = "(count, 0/1/2)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Second position of the composite genotype. Supplementary Table 2",
        "pools TC with CC, so the model uses only whether the count is 0 or",
        "at least 1. ENCORE1 frequencies (Dickinson 2015 Table 1): TT 93%,",
        "TC 6.0%, CC 0.5%."
      ),
      source_name = "CYP2B6 983T>C"
    ),
    SNP_CYP2A6_RS8192726 = list(
      description = paste(
        "CYP2A6*9B (rs8192726) variant-carrier indicator: 1 = carries at",
        "least one variant allele (CA or AA in the source's notation),",
        "0 = CC."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CC)",
      notes = paste(
        "Third position of the composite genotype; Supplementary Table 2",
        "pools CA with AA. rs number from Dickinson 2015 Methods. ENCORE1",
        "frequencies (Dickinson 2015 Table 1): CC 82%, CA 16%, AA 3%. Has no",
        "effect in subjects with CYP2B6 516GG and 983TT (the reference",
        "group spans all CYP2A6 combinations)."
      ),
      source_name = "CYP2A6*9B"
    ),
    SNP_CYP2A6_RS28399454 = list(
      description = paste(
        "CYP2A6*17 (rs28399454) variant-carrier indicator: 1 = carries at",
        "least one variant allele (CT or TT in the source's notation),",
        "0 = CC."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CC)",
      notes = paste(
        "Fourth position of the composite genotype; Supplementary Table 2",
        "pools CT with TT. rs number from Dickinson 2015 Methods. ENCORE1",
        "frequencies (Dickinson 2015 Table 1): CC 90%, CT 10%, TT 0.5%."
      ),
      source_name = "CYP2A6*17"
    ),
    CYP2B6_MISSING = list(
      description = paste(
        "Composite CYP2B6 / CYP2A6 genotype unavailable: 1 = missing,",
        "0 = genotyped."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (genotyped)",
      notes = paste(
        "Supplementary Table 2 'Missing genotype' level, with its own CL/F",
        "multiplier (0.73). When 1, the four SNP columns are ignored. In",
        "ENCORE1, 32 of 606 subjects had no genotyping sample, and a few",
        "more failed amplification at one CYP2A6 position."
      ),
      source_name = "Missing genotype"
    ),
    OCC = list(
      description = "Occasion index (1, 2 or 3) for inter-occasion variability on CL/F",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "ENCORE1 contributed 1-3 occasions per patient (Dickinson 2015",
        "Results), so three occasion etas are carried; the second and third",
        "repeat the first's variance (fixed). Set OCC = 1 throughout for a",
        "single draw of IOV, which is how the paper's virtual patients",
        "(Supplementary Table 1) were generated."
      ),
      source_name = NA_character_
    )
  )

  compartmentData <- list(
    depot = list(analyte = "efavirenz", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "efavirenz", units = "mg", specimen = "plasma", verified = TRUE),
    virus = list(
      analyte = "free infectious HIV-1 virions (wild type)",
      units = "virions",
      specimen = "not applicable",
      verified = TRUE
    ),
    tcell_early = list(
      analyte = "early infected T cells, pre-integration (wild type)",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    tcell_infected = list(
      analyte = "productively infected T cells (wild type)",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    virus_y181c = list(
      analyte = "free infectious HIV-1 virions (Y181C mutant)",
      units = "virions",
      specimen = "not applicable",
      verified = TRUE
    ),
    tcell_early_y181c = list(
      analyte = "early infected T cells, pre-integration (Y181C mutant)",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    tcell_infected_y181c = list(
      analyte = "productively infected T cells (Y181C mutant)",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    virus_k103n = list(
      analyte = "free infectious HIV-1 virions (K103N mutant)",
      units = "virions",
      specimen = "not applicable",
      verified = TRUE
    ),
    tcell_early_k103n = list(
      analyte = "early infected T cells, pre-integration (K103N mutant)",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    tcell_infected_k103n = list(
      analyte = "productively infected T cells (K103N mutant)",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    virus_g190s = list(
      analyte = "free infectious HIV-1 virions (G190S mutant)",
      units = "virions",
      specimen = "not applicable",
      verified = TRUE
    ),
    tcell_early_g190s = list(
      analyte = "early infected T cells, pre-integration (G190S mutant)",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    tcell_infected_g190s = list(
      analyte = "productively infected T cells (G190S mutant)",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 606,
    n_studies = 1,
    age_range = "18-69 years",
    age_median = "35 years",
    weight_range = "39-148 kg",
    weight_median = "65 kg",
    sex_female_pct = 32,
    race_ethnicity = c(
      African = 37,
      Asian = 33,
      Hispanic = 17,
      Caucasian = 13,
      `Aboriginal and Torres Strait Islander` = 0.2
    ),
    disease_state = paste(
      "Treatment-naive adults with HIV-1 (baseline viral load",
      "162-10,000,000 copies/mL) starting efavirenz with tenofovir",
      "disoproxil fumarate / emtricitabine"
    ),
    dose_range = "Efavirenz 400 mg or 600 mg orally once daily",
    regions = "Africa, Asia, South America, Europe and Oceania",
    notes = paste(
      "PK model developed on 1491 plasma samples (median 2, range 1-9 per",
      "patient) from the ENCORE1 trial: single random samples 8-16 h",
      "post-dose at weeks 4 and 12 plus intensive sampling (0-24 h) in a",
      "subgroup between weeks 4 and 8 (Section 2 'Patients'). 574 patients",
      "were genotyped. The prophylaxis analysis applies the model to",
      "HIV-negative people; the viral-dynamics and PD constants come from",
      "the literature (Tables 2-3), not from ENCORE1."
    )
  )

  ini({
    # ---- PK (Section 3.1.1, Table 1; 70 kg, reference genotype) -----------
    lka <- fixed(log(0.6)); label("Absorption rate constant ka (1/h)") # Section 3.1.1, ka fixed to 0.6 h-1 (Arab-Alameddine 2009)
    lcl <- log(11.9); label("Steady-state (post-autoinduction) apparent clearance CLss/F (L/h)") # Table 1, CLss/Fbio = 11.9 L/h (RSE 2.4%)
    lvc <- log(282); label("Apparent volume of distribution V/F (L)") # Table 1, V/Fbio = 282 L (RSE 5.2%)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of WT/70 on CL/F (unitless)") # Table 1 legend, (weight/70)^0.75
    e_wt_vc <- fixed(1); label("Allometric exponent of WT/70 on V/F (unitless)") # Table 1 legend, V(i)/Fbio = V/Fbio * (weight(i)/70)

    # ---- Metabolic autoinduction (Section 3.1.2, Eq. 4; Zhu 2009) ----------
    # alpha = E(CL day 1) / E(CL steady state) = 5.76 / 9.86; the deficit
    # below the steady-state asymptote at the first dose is 1 - alpha.
    cl_time_fdef <- fixed(1 - 5.76 / 9.86); label("Fractional deficit of CL/F below steady state at the first dose, 1 - alpha (unitless)") # Section 3.1.2, alpha = 5.76/9.86 = 0.58 (Table 1 alpha = 0.58)
    lcl_t50 <- fixed(log(245)); label("Half-time of the autoinduction clearance rise T50 (h)") # Section 3.1.2, T50 = 245 h (Zhu 2009)

    # ---- Composite genotype multipliers on CL/F (Supplementary Table 2) ----
    # Code e_geno<a><b><c><d>_cl: a = CYP2B6 516 T-allele count (0 GG,
    # 1 GT, 2 TT); b = CYP2B6 983 C carrier (TC or CC); c = CYP2A6*9B
    # carrier (CA or AA); d = CYP2A6*17 carrier (CT or TT). The reference
    # (516GG / 983TT, any CYP2A6) is 1.
    e_geno0100_cl <- 0.73; label("CL/F multiplier, CYP2B6 GG/TC-CC, CYP2A6 CC/CC (unitless)") # Supplementary Table 2, 0.73 (RSE 18.4%)
    e_geno0101_cl <- 0.30; label("CL/F multiplier, CYP2B6 GG/TC-CC, CYP2A6 CC/CT-TT (unitless)") # Supplementary Table 2, 0.30 (RSE 20.6%)
    e_geno0110_cl <- 0.75; label("CL/F multiplier, CYP2B6 GG/TC-CC, CYP2A6 CA-AA/CC (unitless)") # Supplementary Table 2, 0.75 (RSE 32.8%)
    e_geno0111_cl <- 0.96; label("CL/F multiplier, CYP2B6 GG/TC-CC, CYP2A6 CA-AA/CT (unitless)") # Supplementary Table 2, 0.96 (RSE 2.4%)
    e_geno1000_cl <- 0.75; label("CL/F multiplier, CYP2B6 GT/TT, CYP2A6 CC/CC (unitless)") # Supplementary Table 2, 0.75 (RSE 3.7%)
    e_geno1001_cl <- 0.68; label("CL/F multiplier, CYP2B6 GT/TT, CYP2A6 CC/CT-TT (unitless)") # Supplementary Table 2, 0.68 (RSE 8.8%)
    e_geno1010_cl <- 0.69; label("CL/F multiplier, CYP2B6 GT/TT, CYP2A6 CA-AA/CC (unitless)") # Supplementary Table 2, 0.69 (RSE 5.7%)
    e_geno1011_cl <- 0.71; label("CL/F multiplier, CYP2B6 GT/TT, CYP2A6 CA-AA/CT (unitless)") # Supplementary Table 2, 0.71 (RSE 5.5%)
    e_geno1100_cl <- 0.42; label("CL/F multiplier, CYP2B6 GT/TC-CC, CYP2A6 CC/CC (unitless)") # Supplementary Table 2, 0.42 (RSE 24.5%)
    e_geno1101_cl <- 0.39; label("CL/F multiplier, CYP2B6 GT/TC-CC, CYP2A6 CC/CT-TT (unitless)") # Supplementary Table 2, 0.39 (RSE 30.1%)
    e_geno1110_cl <- 0.18; label("CL/F multiplier, CYP2B6 GT/TC-CC, CYP2A6 CA-AA/CC (unitless)") # Supplementary Table 2, 0.18 (RSE 18.1%)
    e_geno2000_cl <- 0.33; label("CL/F multiplier, CYP2B6 TT/TT, CYP2A6 CC/CC (unitless)") # Supplementary Table 2, 0.33 (RSE 6.6%)
    e_geno2001_cl <- 0.34; label("CL/F multiplier, CYP2B6 TT/TT, CYP2A6 CC/CT-TT (unitless)") # Supplementary Table 2, 0.34 (RSE 19.6%)
    e_geno2010_cl <- 0.44; label("CL/F multiplier, CYP2B6 TT/TT, CYP2A6 CA-AA/CC (unitless)") # Supplementary Table 2, 0.44 (RSE 25.0%)
    e_geno2011_cl <- 0.75; label("CL/F multiplier, CYP2B6 TT/TT, CYP2A6 CA-AA/CT (unitless)") # Supplementary Table 2, 0.75 (RSE 74.1%)
    e_genomiss_cl <- 0.73; label("CL/F multiplier, missing genotype (unitless)") # Supplementary Table 2, 'Missing genotype' 0.73 (RSE 8.3%)

    # ---- Direct effect (Section 3.2, Eq. 5, Table 2, Text S1) --------------
    lic50 <- fixed(log(5.4)); label("Ex vivo IC50 of wild-type HIV-1 in 50% human serum, before protein-binding correction (nM)") # Table 2, IC50 wild type = 5.4 (sd 0.9) nmol/L (Shen 2008)
    hill <- fixed(1.69); label("Hill coefficient m, wild-type HIV-1 (unitless)") # Table 2, m wild type = 1.69 (sd 0.08)
    fu <- fixed(0.002); label("Fraction of efavirenz unbound in plasma (unitless)") # Figure 2 legend and Text S1, fu,plasma = 0.2%
    ic50_fold_y181c <- fixed(2.8); label("IC50 fold change, Y181C vs wild type (unitless)") # Table 2, Y181C IC50 = 2.8 * IC50(wt) (Sampah 2011)
    ic50_fold_k103n <- fixed(89.1); label("IC50 fold change, K103N vs wild type (unitless)") # Table 2, K103N IC50 = 89.1 * IC50(wt)
    ic50_fold_g190s <- fixed(72.1); label("IC50 fold change, G190S vs wild type (unitless)") # Table 2, G190S IC50 = 72.1 * IC50(wt)
    hill_fold_y181c <- fixed(0.9); label("Hill coefficient fold change, Y181C vs wild type (unitless)") # Table 2, Y181C m = 0.9 * m(wt)
    hill_fold_k103n <- fixed(0.83); label("Hill coefficient fold change, K103N vs wild type (unitless)") # Table 2, K103N m = 0.83 * m(wt)
    hill_fold_g190s <- fixed(0.6); label("Hill coefficient fold change, G190S vs wild type (unitless)") # Table 2, G190S m = 0.6 * m(wt)
    fitness_y181c <- fixed(0.78); label("Relative fitness f of the Y181C mutant (unitless)") # Table 2, f = 0.78
    fitness_k103n <- fixed(0.74); label("Relative fitness f of the K103N mutant (unitless)") # Table 2, f = 0.74
    fitness_g190s <- fixed(0.24); label("Relative fitness f of the G190S mutant (unitless)") # Table 2, f = 0.24

    # ---- Viral dynamics (Section 3.3, Table 3; rates in 1/day) -------------
    lambda_t <- fixed(2e9); label("Production of uninfected T cells lambdaT (cells/day)") # Table 3, lambdaT = 2e9 (Wei 1995)
    delta_t <- fixed(0.02); label("Death rate of uninfected T cells deltaT (1/day)") # Table 3, deltaT = 0.02 (Sedaghat 2009)
    delta_t1 <- fixed(0.02); label("Death rate of early infected T cells deltaT1 (1/day)") # Table 3, deltaT1 = 0.02 (Sedaghat 2009)
    delta_t2 <- fixed(1); label("Death rate of productively infected T cells deltaT2 (1/day)") # Table 3, deltaT2 = 1 (Markowitz 2003)
    delta_pic_t <- fixed(0.35); label("Intracellular destruction of the pre-integration complex deltaPIC (1/day)") # Table 3, deltaPIC = 0.35 (Zhou 2005)
    k_t <- fixed(0.35); label("Proviral integration (early to productive) rate k (1/day)") # Table 3, k = 0.35 (Zhou 2005)
    beta_t <- fixed(8e-12); label("Infection rate of T cells beta (1/(day virion))") # Table 3, beta = 8e-12 (Sedaghat 2008)
    n_t <- fixed(670); label("Infectious virions released per productively infected T cell NT (1/day)") # Table 3, NT = 670 (Sedaghat 2009)
    cl_virus_naive <- fixed(2.3); label("Clearance rate of free virus in a virus-naive individual CLV (1/day)") # Table 3, CLV = 2.3 (Tan and Wu 1998)
    rho_rev <- fixed(0.5); label("Probability of completing reverse transcription without drug rho_rev (unitless)") # Section 3.3, rho_rev = 0.5 (Pierson 2002)

    # ---- Between-subject and inter-occasion variability --------------------
    # Log-normal, CV = 100 * sqrt(exp(omega^2) - 1) (Table 1 legend), so
    # omega^2 = log(1 + CV^2).
    etalcl ~ 0.1256 # Table 1, CV_IIV(CLss/Fbio) = 36.6%: log(1 + 0.366^2) = 0.1256
    etaiov_lcl_1 ~ 0.04315 # Table 1, CV_IOV(CLss/Fbio) = 21.0%: log(1 + 0.210^2) = 0.04315
    etaiov_lcl_2 ~ fixed(0.04315) # same IOV variance, occasion 2
    etaiov_lcl_3 ~ fixed(0.04315) # same IOV variance, occasion 3
    # PD parameter variability (Table 2 legend): IC50 log-normal with the
    # wild-type CV (0.9/5.4 = 16.7%); m normal with sd 0.08. Mutants share
    # the wild-type CV, so one eta of each scales all four strains.
    etalic50 ~ 0.02740 # Table 2, IC50 sd 0.9 on 5.4: log(1 + (0.9/5.4)^2) = 0.02740
    etahill ~ 0.0064 # Table 2, m sd 0.08: 0.08^2 = 0.0064 (additive, normal)

    # ---- Residual error -----------------------------------------------------
    propSd <- 0.20; label("Proportional residual error (fraction)") # Table 1, sigma = 0.20 (RSE 8.6%)
  })

  model({
    # ---- 1. Composite genotype multiplier on CL/F (Supplementary Table 2) --
    g516 <- SNP_CYP2B6_RS3745274_T_COUNT
    c983 <- 0
    if (SNP_CYP2B6_RS28399499_C_COUNT >= 1) c983 <- 1
    c9b <- SNP_CYP2A6_RS8192726
    c17 <- SNP_CYP2A6_RS28399454
    b0 <- (g516 == 0) * (1 - c983)
    b1 <- (g516 == 0) * c983
    b2 <- (g516 == 1) * (1 - c983)
    b3 <- (g516 == 1) * c983
    b4 <- (g516 == 2) * (1 - c983)
    a00 <- (1 - c9b) * (1 - c17)
    a01 <- (1 - c9b) * c17
    a10 <- c9b * (1 - c17)
    a11 <- c9b * c17
    geno_fac <- b0 +
      b1 * (a00 * e_geno0100_cl + a01 * e_geno0101_cl + a10 * e_geno0110_cl + a11 * e_geno0111_cl) +
      b2 * (a00 * e_geno1000_cl + a01 * e_geno1001_cl + a10 * e_geno1010_cl + a11 * e_geno1011_cl) +
      b3 * (a00 * e_geno1100_cl + a01 * e_geno1101_cl + a10 * e_geno1110_cl) +
      b4 * (a00 * e_geno2000_cl + a01 * e_geno2001_cl + a10 * e_geno2010_cl + a11 * e_geno2011_cl)
    # Combinations absent from Supplementary Table 2 (516TT with a 983 C
    # allele; 516GT / 983 C carrier / CYP2A6 *9B and *17 carrier) have no
    # estimate, so they return NA rather than borrowing another group's.
    geno_listed <- b0 + b1 + b2 + b4 + b3 * (1 - a11)
    if (geno_listed == 0) geno_fac <- NA
    if (CYP2B6_MISSING == 1) geno_fac <- e_genomiss_cl

    # ---- 2. Individual PK parameters ----------------------------------------
    iov_cl <- (OCC == 1) * etaiov_lcl_1 + (OCC == 2) * etaiov_lcl_2 +
      (OCC == 3) * etaiov_lcl_3
    ka <- exp(lka)
    cl_ss <- exp(lcl + etalcl + iov_cl) * geno_fac * (WT / 70)^e_wt_cl
    vc <- exp(lvc) * (WT / 70)^e_wt_vc

    # Eq. 4: CL(t) = CL_t0 + (CL_ss - CL_t0) * t / (t + T50), with
    # CL_t0 = alpha * CL_ss and t the time since the first efavirenz dose;
    # equivalently CL_ss * (1 - (1 - alpha) * T50 / (t + T50)).
    cl_t50 <- exp(lcl_t50)
    t_efv <- tafd(depot)
    if (is.na(t_efv)) t_efv <- 0
    cl <- cl_ss * (1 - cl_time_fdef * cl_t50 / (t_efv + cl_t50))
    kel <- cl / vc

    # ---- 3. PK ODEs (Eqs. 1-3) ------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    Cc <- central / vc

    # ---- 4. Direct effect per strain (Eq. 5, Text S1 Eq. S1.10) --------------
    # The ex vivo IC50 (50% human serum) is converted to a total-plasma IC50
    # with IC50_p = IC50_assay * 2 / (fu + 1) (Eq. S1.10), then from nM to
    # mg/L with the molar mass 315.6 g/mol (Section 3.1).
    mw_efv <- 315.6
    ic50_wt <- exp(lic50 + etalic50) * 2 / (fu + 1) * mw_efv * 1e-6
    hill_wt <- hill + etahill
    cpos <- max(Cc, 0)
    eta_d_wt <- cpos^hill_wt / (ic50_wt^hill_wt + cpos^hill_wt)
    hill_y181c <- hill_fold_y181c * hill_wt
    ic50_y181c <- ic50_fold_y181c * ic50_wt
    eta_d_y181c <- cpos^hill_y181c / (ic50_y181c^hill_y181c + cpos^hill_y181c)
    hill_k103n <- hill_fold_k103n * hill_wt
    ic50_k103n <- ic50_fold_k103n * ic50_wt
    eta_d_k103n <- cpos^hill_k103n / (ic50_k103n^hill_k103n + cpos^hill_k103n)
    hill_g190s <- hill_fold_g190s * hill_wt
    ic50_g190s <- ic50_fold_g190s * ic50_wt
    eta_d_g190s <- cpos^hill_g190s / (ic50_g190s^hill_g190s + cpos^hill_g190s)

    # ---- 5. Viral dynamics (Eqs. 6-11; rates in 1/day, divided by 24) -------
    # Expected values of the stochastic reactions, which are linear in the
    # viral states because uninfected T cells are held at Tu = lambdaT /
    # deltaT. a1 clears free virus at CLV + CLT * Tu with
    # CLT = (1/rho_rev - (1 - eta_D)) * beta * f; a4 infects at
    # (1 - eta_D) * beta * f * Tu; a6 produces NT virions per T2 cell.
    tu <- lambda_t / delta_t
    a_t1 <- delta_pic_t + delta_t1 + k_t

    inf_wt <- (1 - eta_d_wt) * beta_t * tu
    clr_wt <- cl_virus_naive + (1 / rho_rev - (1 - eta_d_wt)) * beta_t * tu
    d/dt(virus) <- (n_t * tcell_infected - (clr_wt + inf_wt) * virus) / 24
    d/dt(tcell_early) <- (inf_wt * virus - a_t1 * tcell_early) / 24
    d/dt(tcell_infected) <- (k_t * tcell_early - delta_t2 * tcell_infected) / 24

    inf_y181c <- (1 - eta_d_y181c) * beta_t * fitness_y181c * tu
    clr_y181c <- cl_virus_naive + (1 / rho_rev - (1 - eta_d_y181c)) * beta_t * fitness_y181c * tu
    d/dt(virus_y181c) <- (n_t * tcell_infected_y181c - (clr_y181c + inf_y181c) * virus_y181c) / 24
    d/dt(tcell_early_y181c) <- (inf_y181c * virus_y181c - a_t1 * tcell_early_y181c) / 24
    d/dt(tcell_infected_y181c) <- (k_t * tcell_early_y181c - delta_t2 * tcell_infected_y181c) / 24

    inf_k103n <- (1 - eta_d_k103n) * beta_t * fitness_k103n * tu
    clr_k103n <- cl_virus_naive + (1 / rho_rev - (1 - eta_d_k103n)) * beta_t * fitness_k103n * tu
    d/dt(virus_k103n) <- (n_t * tcell_infected_k103n - (clr_k103n + inf_k103n) * virus_k103n) / 24
    d/dt(tcell_early_k103n) <- (inf_k103n * virus_k103n - a_t1 * tcell_early_k103n) / 24
    d/dt(tcell_infected_k103n) <- (k_t * tcell_early_k103n - delta_t2 * tcell_infected_k103n) / 24

    inf_g190s <- (1 - eta_d_g190s) * beta_t * fitness_g190s * tu
    clr_g190s <- cl_virus_naive + (1 / rho_rev - (1 - eta_d_g190s)) * beta_t * fitness_g190s * tu
    d/dt(virus_g190s) <- (n_t * tcell_infected_g190s - (clr_g190s + inf_g190s) * virus_g190s) / 24
    d/dt(tcell_early_g190s) <- (inf_g190s * virus_g190s - a_t1 * tcell_early_g190s) / 24
    d/dt(tcell_infected_g190s) <- (k_t * tcell_early_g190s - delta_t2 * tcell_infected_g190s) / 24

    # ---- 6. Prophylactic efficacy at a constant concentration (Eq. 13) ------
    # Extinction probabilities of the branching process defined by
    # reactions (6)-(11) when eta_D is held constant: from one virion
    # P_V = 1 - A + A * (b1 + b3 * P_T2), with A = a4 / (a1 + a4) the chance
    # that the virion infects, b3 = k / (deltaPIC + deltaT1 + k) the chance
    # that an early infected cell integrates, b1 = 1 - b3, and P_T2 the
    # smallest root of g*A*b3*x^2 + (g*(1 - A*b3) - 1)*x + (1 - g) = 0, where
    # g = NT / (deltaT2 + NT). The root is written in the form that stays
    # finite as A -> 0. phi = 1 - (1 - P_V(drug)) / (1 - P_V(no drug)) for a
    # single transmitted virion (Y0 = [1, 0, 0], as in Figure 2A).
    pe_b3 <- k_t / a_t1
    pe_g <- n_t / (delta_t2 + n_t)
    pe_ec <- 1 - pe_g

    # Drug-free infection chance per virion: a4 / (a1 + a4) at eta_D = 0.
    pa_wt0 <- beta_t * tu / (cl_virus_naive + (1 / rho_rev) * beta_t * tu)
    pe_ea <- pe_g * pa_wt0 * pe_b3
    pe_eb <- pe_g - 1 - pe_ea
    pt2_wt0 <- 2 * pe_ec / (-pe_eb + sqrt(pe_eb^2 - 4 * pe_ea * pe_ec))
    pinf_wt0 <- pa_wt0 * pe_b3 * (1 - pt2_wt0)
    pa_wt <- inf_wt / (clr_wt + inf_wt)
    pe_ea <- pe_g * pa_wt * pe_b3
    pe_eb <- pe_g - 1 - pe_ea
    pt2_wt <- 2 * pe_ec / (-pe_eb + sqrt(pe_eb^2 - 4 * pe_ea * pe_ec))
    phi_wt <- 1 - pa_wt * pe_b3 * (1 - pt2_wt) / pinf_wt0

    pa_y181c0 <- beta_t * fitness_y181c * tu / (cl_virus_naive + (1 / rho_rev) * beta_t * fitness_y181c * tu)
    pe_ea <- pe_g * pa_y181c0 * pe_b3
    pe_eb <- pe_g - 1 - pe_ea
    pt2_y181c0 <- 2 * pe_ec / (-pe_eb + sqrt(pe_eb^2 - 4 * pe_ea * pe_ec))
    pinf_y181c0 <- pa_y181c0 * pe_b3 * (1 - pt2_y181c0)
    pa_y181c <- inf_y181c / (clr_y181c + inf_y181c)
    pe_ea <- pe_g * pa_y181c * pe_b3
    pe_eb <- pe_g - 1 - pe_ea
    pt2_y181c <- 2 * pe_ec / (-pe_eb + sqrt(pe_eb^2 - 4 * pe_ea * pe_ec))
    phi_y181c <- 1 - pa_y181c * pe_b3 * (1 - pt2_y181c) / pinf_y181c0

    pa_k103n0 <- beta_t * fitness_k103n * tu / (cl_virus_naive + (1 / rho_rev) * beta_t * fitness_k103n * tu)
    pe_ea <- pe_g * pa_k103n0 * pe_b3
    pe_eb <- pe_g - 1 - pe_ea
    pt2_k103n0 <- 2 * pe_ec / (-pe_eb + sqrt(pe_eb^2 - 4 * pe_ea * pe_ec))
    pinf_k103n0 <- pa_k103n0 * pe_b3 * (1 - pt2_k103n0)
    pa_k103n <- inf_k103n / (clr_k103n + inf_k103n)
    pe_ea <- pe_g * pa_k103n * pe_b3
    pe_eb <- pe_g - 1 - pe_ea
    pt2_k103n <- 2 * pe_ec / (-pe_eb + sqrt(pe_eb^2 - 4 * pe_ea * pe_ec))
    phi_k103n <- 1 - pa_k103n * pe_b3 * (1 - pt2_k103n) / pinf_k103n0

    pa_g190s0 <- beta_t * fitness_g190s * tu / (cl_virus_naive + (1 / rho_rev) * beta_t * fitness_g190s * tu)
    pe_ea <- pe_g * pa_g190s0 * pe_b3
    pe_eb <- pe_g - 1 - pe_ea
    pt2_g190s0 <- 2 * pe_ec / (-pe_eb + sqrt(pe_eb^2 - 4 * pe_ea * pe_ec))
    pinf_g190s0 <- pa_g190s0 * pe_b3 * (1 - pt2_g190s0)
    pa_g190s <- inf_g190s / (clr_g190s + inf_g190s)
    pe_ea <- pe_g * pa_g190s * pe_b3
    pe_eb <- pe_g - 1 - pe_ea
    pt2_g190s <- 2 * pe_ec / (-pe_eb + sqrt(pe_eb^2 - 4 * pe_ea * pe_ec))
    phi_g190s <- 1 - pa_g190s * pe_b3 * (1 - pt2_g190s) / pinf_g190s0

    # ---- 7. Observation -------------------------------------------------------
    Cc ~ prop(propSd)
  })
}
