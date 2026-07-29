Jeong_2022_torsemide <- function() {
  description <- "Two-compartment population PK model for oral torsemide in healthy Korean adult males (Jeong 2022), with first-order absorption after a lag time, proportional residual error, and categorical genotype covariates: OATP1B1 *15 haplotype (intermediate / poor transporter) reduces apparent central volume, and CYP2C9 extensive-metabolizer phenotype increases apparent oral clearance and apparent inter-compartmental clearance."
  reference <- "Jeong S-H, Jang J-H, Cho H-Y, Lee Y-B. Population Pharmacokinetic (Pop-PK) Analysis of Torsemide in Healthy Korean Males Considering CYP2C9 and OATP1B1 Genetic Polymorphisms. Pharmaceutics. 2022;14(4):771. doi:10.3390/pharmaceutics14040771"
  vignette <- "Jeong_2022_torsemide"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    SLCO1B1_HAP15_HET = list(
      description        = "SLCO1B1 *15 haplotype heterozygote indicator (intermediate transporter, IT)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no *15 allele; ET phenotype: *1a/*1a, *1a/*1b, or *1b/*1b)",
      notes              = "Time-fixed per subject. In Jeong 2022 (Table 1) the IT phenotype pools *1a/*15 (n = 4) and *1b/*15 (n = 19); together 23 of 112 subjects (20.5%) are IT carriers. Paired with SLCO1B1_HAP15_HOM to encode the three-level OATP1B1 phenotype (ET / IT / PT) with ET as the implicit reference (both indicators = 0).",
      source_name        = "OATP1B1 phenotype"
    ),
    SLCO1B1_HAP15_HOM = list(
      description        = "SLCO1B1 *15 haplotype homozygote indicator (poor transporter, PT)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no *15 allele; ET phenotype: *1a/*1a, *1a/*1b, or *1b/*1b)",
      notes              = "Time-fixed per subject. In Jeong 2022 (Table 1) only 3 of 112 subjects (2.68%) carry the *15/*15 PT genotype, so the dV/F coefficient for PT carries a wide 95% bootstrap CI (-1.00 to -0.10; Table 5). Paired with SLCO1B1_HAP15_HET.",
      source_name        = "OATP1B1 phenotype"
    ),
    CYP2C9_EM = list(
      description        = "CYP2C9 extensive-metabolizer phenotype indicator: 1 if *1/*1 (wild-type homozygote, EM), 0 if *1/*3 or *1/*13 (reduced-function-allele carrier, IM)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP2C9 intermediate metabolizer; carrier of *3 or *13)",
      notes              = "Time-fixed per subject. In Jeong 2022 (Table 1) the EM phenotype is *1/*1 (n = 97, 86.6%) and the IM phenotype pools *1/*3 (n = 12, 10.7%) and *1/*13 (n = 3, 2.7%); together 15 of 112 subjects (13.4%) are IM. No CYP2C9 poor metabolizers (*3/*3, *3/*13, or *13/*13) were observed in the cohort. The Jeong 2022 final model uses IM as the reference and a linear-deviation positive coefficient on EM for both CL/F and CL2/F (apparent inter-compartmental clearance), reflecting that EMs metabolize torsemide ~51% faster than IMs (Table 4).",
      source_name        = "CYP2C9 phenotype"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 112L,
    n_studies        = 4L,
    n_observations   = 1344L,
    age_range        = "19-29 years",
    age_mean         = "23.38 years",
    weight_range     = "44-88.4 kg",
    weight_mean      = "67.84 kg",
    height_range     = "159.9-188.1 cm",
    height_mean      = "173.78 cm",
    sex_female_pct   = 0,
    race_ethnicity   = "Korean (single-ancestry East Asian cohort recruited at Chonnam National University, Gwangju).",
    disease_state    = "Healthy adult males with no history of congenital or chronic disease; all clinical chemistry, hematology, and urine parameters within normal ranges prior to dosing.",
    dose_range       = "Oral torsemide 5 mg (n = 28), 10 mg (n = 28), or 20 mg (n = 56) single dose; no dose-proportionality deviations across 5-20 mg.",
    regions          = "Republic of Korea (Chonnam National University Bioequivalence and Bridging Study Centre, Gwangju).",
    oatp1b1_distribution = "ET (*1a/*1a, *1a/*1b, *1b/*1b) n = 86 (76.8%); IT (*1a/*15, *1b/*15) n = 23 (20.5%); PT (*15/*15) n = 3 (2.7%). Jeong 2022 Table 1.",
    cyp2c9_distribution  = "EM (*1/*1) n = 97 (86.6%); IM (*1/*3) n = 12 (10.7%); IM (*1/*13) n = 3 (2.7%). No PM (*3/*3, *3/*13, *13/*13) observed. Jeong 2022 Table 1.",
    notes            = "Pooled retrospective analysis of four open-label single-dose two-period crossover bioequivalence studies (2004, 2005, 2006, 2007) using a single reference torsemide tablet formulation. Only samples following the reference formulation were included. Serum torsemide measured by HPLC-UV (assay range 0.02-10 ug/mL). Blood sampled at 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, and 12 h post-dose (12 timepoints per subject including pre-dose). Biochemical covariates (albumin, total proteins, creatinine, CrCl, BMI, BSA, AST, ALT, ALP, GFR, BUN, cholesterol, total bilirubin) were tested but not retained in the final model because of the narrow physiological range in healthy young males. Jeong 2022 Methods Section 2.2 and Results Section 3.1."
  )

  ini({
    # Structural PK -- Jeong 2022 Table 4 final-model estimates (Phoenix NLME
    # FOCE-ELS with eta-eps interaction; ADVAN4/TRANS4 with ALAG1 in NONMEM
    # confirmation). Time in hours; apparent clearances (CL/F, CL2/F = Q/F)
    # in L/h; apparent volumes (V/F = Vc/F, V2/F = Vp/F) in L. The Table 4
    # reference subject is a CYP2C9 intermediate metabolizer (IM; *1/*3 or
    # *1/*13) who is an OATP1B1 extensive transporter (ET; no *15 allele);
    # the typical-value parameters reproduce that subject directly and the
    # genotype-effect coefficients below adjust to other phenotypes.
    lka   <- log(0.861) ; label("Absorption rate constant ka (1/h)")                            # Jeong 2022 Table 4 final tvKa = 0.861 1/h
    ltlag <- log(0.274) ; label("Absorption lag time (h)")                                       # Jeong 2022 Table 4 final tvTlag = 0.274 h
    lvc   <- log(1.027) ; label("Apparent central volume V/F (L) at OATP1B1 ET reference")       # Jeong 2022 Table 4 final tvV/F = 1.027 L
    lcl   <- log(1.740) ; label("Apparent oral clearance CL/F (L/h) at CYP2C9 IM reference")     # Jeong 2022 Table 4 final tvCL/F = 1.740 L/h
    lvp   <- log(3.914) ; label("Apparent peripheral volume V2/F (L)")                           # Jeong 2022 Table 4 final tvV2/F = 3.914 L
    lq    <- log(0.828) ; label("Apparent inter-compartmental clearance Q/F (L/h) at CYP2C9 IM reference") # Jeong 2022 Table 4 final tvCL2/F = 0.828 L/h

    # Genotype covariate effects -- Jeong 2022 Table 4 and Section 3.5
    # narrative. Linear-deviation parameterization:
    #   V/F  = tvV/F  * (1 + e_slco1b1_hap15_het_vc * SLCO1B1_HAP15_HET
    #                       + e_slco1b1_hap15_hom_vc * SLCO1B1_HAP15_HOM)
    #   CL/F = tvCL/F * (1 + e_cyp2c9_em_cl * CYP2C9_EM)
    #   Q/F  = tvQ/F  * (1 + e_cyp2c9_em_q  * CYP2C9_EM)
    # Reference subject (all three indicators = 0) is OATP1B1 ET / CYP2C9 IM,
    # for which the typical values above apply directly.
    e_slco1b1_hap15_het_vc <- -0.410 ; label("OATP1B1 *15 heterozygote (IT) linear-deviation coefficient on V/F (fraction)")  # Jeong 2022 Table 4 dV/FdOATP1B1 IT = -0.410
    e_slco1b1_hap15_hom_vc <- -0.646 ; label("OATP1B1 *15 homozygote (PT) linear-deviation coefficient on V/F (fraction)")    # Jeong 2022 Table 4 dV/FdOATP1B1 PT = -0.646
    e_cyp2c9_em_cl         <-  0.510 ; label("CYP2C9 extensive-metabolizer (EM) linear-deviation coefficient on CL/F (fraction)")  # Jeong 2022 Table 4 dCL/FdCYP2C9 EM = 0.510
    e_cyp2c9_em_q          <-  0.365 ; label("CYP2C9 extensive-metabolizer (EM) linear-deviation coefficient on Q/F (fraction)")   # Jeong 2022 Table 4 dCL2/FdCYP2C9 EM = 0.365

    # Inter-individual variability -- Jeong 2022 Table 4 omega^2 final-model
    # estimates (variances on the log scale; the paper's "IIV (%)" column is
    # 100 * sqrt(omega^2), e.g., sqrt(0.562) = 0.7497 ~= 75% on V/F). The
    # final model does NOT estimate IIV on ka (Model 08 in Table 2 dropped
    # IIV on ka because removal decreased -2LL by 8.8 while reducing nParam
    # by 1; the paper retains ka as a typical value only).
    etalvc   ~ 0.562  # Jeong 2022 Table 4 omega^2(V/F)   = 0.562 (IIV 74.95%)
    etalcl   ~ 0.029  # Jeong 2022 Table 4 omega^2(CL/F)  = 0.029 (IIV 17.10%)
    etalvp   ~ 0.036  # Jeong 2022 Table 4 omega^2(V2/F)  = 0.036 (IIV 18.94%)
    etalq    ~ 0.022  # Jeong 2022 Table 4 omega^2(CL2/F) = 0.022 (IIV 14.90%)
    etaltlag ~ 0.341  # Jeong 2022 Table 4 omega^2(Tlag)  = 0.341 (IIV 58.37%)

    # Residual unexplained variability -- Jeong 2022 Table 4 reports the
    # proportional error model selected in Table 2 step 03 vs additive /
    # log-additive / mixed / power alternatives. Phoenix NLME convention:
    # epsilon (Table 4 'epsilon' row, no squared superscript) is the
    # residual standard deviation on the proportional-error scale.
    propSd <- 0.136 ; label("Proportional residual error (fraction)")  # Jeong 2022 Table 4 final epsilon = 0.136 (RSE 7.15%)
  })

  model({
    # Individual PK parameters with Jeong 2022 covariate equations.
    # Reference subject: OATP1B1 ET (SLCO1B1_HAP15_HET = 0 AND
    # SLCO1B1_HAP15_HOM = 0) and CYP2C9 IM (CYP2C9_EM = 0).
    ka   <- exp(lka)
    tlag <- exp(ltlag + etaltlag)
    vc   <- exp(lvc + etalvc) *
            (1 + e_slco1b1_hap15_het_vc * SLCO1B1_HAP15_HET +
                 e_slco1b1_hap15_hom_vc * SLCO1B1_HAP15_HOM)
    cl   <- exp(lcl + etalcl) * (1 + e_cyp2c9_em_cl * CYP2C9_EM)
    vp   <- exp(lvp + etalvp)
    q    <- exp(lq  + etalq)  * (1 + e_cyp2c9_em_q  * CYP2C9_EM)

    # Micro-constants for the two-compartment disposition.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment oral PK with first-order absorption + absorption lag
    # time. Dose lands in `depot`; bioavailability F is implicit in the
    # apparent-clearance / apparent-volume scaling (Jeong 2022 reports V/F
    # and CL/F, not the absolute V and CL, because the analysis is oral-only
    # and F is not separately identifiable).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    lag(depot) <- tlag

    # Concentrations in central compartment. Dose units = mg, vc = L, so
    # central / vc is in mg/L = ug/mL (Jeong 2022 Figure 1 reports the
    # concentration axis in ug/mL).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
