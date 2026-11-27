Scholl_2018_egc <- function() {
  description <- "Two-compartment population PK model for epigallocatechin (EGC) after oral green tea extract in healthy adults (Scholl 2018), with lagged first-order absorption; oral-contraceptive use and ABCC2 c.-24C>T (rs717620) T alleles lower the apparent central volume, ABCC2 c.3972C>T (rs3740066) C/C raises the intercompartmental clearance, SLCO1B1 c.388A>G (rs2306283) A/A and UGT1A1*28 (rs8175347) 6/6 wild type lower the apparent clearance, COMT rs4680 A/A lowers relative bioavailability, and residual error is proportional."
  reference <- paste(
    "Scholl C, Lepper A, Lehr T, Hanke N, Schneider KL, Brockmoller J, Seufferlein T, Stingl JC. (2018).",
    "Population nutrikinetics of green tea extract.",
    "PLoS One 13(2):e0193074. doi:10.1371/journal.pone.0193074"
  )
  vignette <- "Scholl_2018_green_tea_catechins"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  compartmentData <- list(
    depot = list(analyte = "epigallocatechin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "epigallocatechin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "epigallocatechin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CONMED_BIRTHCONTROL = list(
      description = "Oral (hormonal) contraceptive use",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Scholl 2018 CONTRACEPTIVES column (No = 0, Yes = 1). Multiplies Vcentral/F by 0.450 when 1 (-55 percent; Table 5). Reference group pools men and women not using contraceptives. 37 of 54 women used contraceptives (Table 1).",
      source_name = "CONTRACEPTIVES"
    ),
    SNP_ABCC2_RS717620_HET = list(
      description = "ABCC2 (MRP2) c.-24C>T (rs717620) heterozygote (C/T) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Scholl 2018 coded MRP71T as the T-allele count (C/C = 0, C/T = 1, T/T = 2) and entered it linearly on Vcentral/F: V2 * (1 + theta * MRP71T), theta = -0.370 (-37 percent C/T, -74 percent T/T; Table 5). MRP71T = HET + 2 * HOM.",
      source_name = "MRP71T"
    ),
    SNP_ABCC2_RS717620_HOM = list(
      description = "ABCC2 (MRP2) c.-24C>T (rs717620) homozygous-variant (T/T) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Companion to SNP_ABCC2_RS717620_HET; derive from MRP71T == 2.",
      source_name = "MRP71T"
    ),
    SNP_ABCC2_RS3740066_HET = list(
      description = "ABCC2 (MRP2) c.3972C>T (rs3740066) heterozygote (C/T) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Scholl 2018 coded MRP37T as the T-allele count (C/C = 0, C/T = 1, T/T = 2). The effect applies to the C/C wild-type stratum only (HET = 0 and HOM = 0): Q/F multiplied by 2.33 (+133 percent; Table 5). Genotype frequencies C/C 40, C/T 32, T/T 12 (Table 4).",
      source_name = "MRP37T"
    ),
    SNP_ABCC2_RS3740066_HOM = list(
      description = "ABCC2 (MRP2) c.3972C>T (rs3740066) homozygous-variant (T/T) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Companion to SNP_ABCC2_RS3740066_HET; derive from MRP37T == 2.",
      source_name = "MRP37T"
    ),
    SNP_SLCO1B1_RS2306283_HET = list(
      description = "SLCO1B1 (OATP1B1) c.388A>G (rs2306283) heterozygote (A/G) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Scholl 2018 coded OATP23T as the G-allele count (A/A = 0, G/A = 1, G/G = 2). The effect applies to the A/A wild-type stratum only (HET = 0 and HOM = 0): CL/F multiplied by 0.653 (-35 percent; Table 5). Genotype frequencies A/A 29, G/A 40, G/G 15 (Table 4).",
      source_name = "OATP23T"
    ),
    SNP_SLCO1B1_RS2306283_HOM = list(
      description = "SLCO1B1 (OATP1B1) c.388A>G (rs2306283) homozygous-variant (G/G) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Companion to SNP_SLCO1B1_RS2306283_HET; derive from OATP23T == 2.",
      source_name = "OATP23T"
    ),
    SNP_UGT1A1_RS8175347_HET = list(
      description = "UGT1A1*28 (rs8175347) heterozygote ((TA)6/(TA)7) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Scholl 2018 coded UGT as (-) = 0, (-/TA) = 1, TA = 2. The effect applies to the (-) = (TA)6/(TA)6 wild-type stratum only (HET = 0 and HOM = 0): CL/F multiplied by 0.744 (-26 percent; Table 5). Genotype frequencies 37 / 37 / 10 (Table 4).",
      source_name = "UGT"
    ),
    SNP_UGT1A1_RS8175347_HOM = list(
      description = "UGT1A1*28 (rs8175347) homozygous-variant ((TA)7/(TA)7) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Companion to SNP_UGT1A1_RS8175347_HET; derive from UGT == 2.",
      source_name = "UGT"
    ),
    SNP_COMT_RS4680_HOM = list(
      description = "COMT rs4680 (Val158Met; c.1947G>A in Scholl 2018) homozygous A/A (Met/Met) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Scholl 2018 coded COMT as GG = 0, GA = 1, AA = 2; the effect applies when COMT == 2: F1 multiplied by 0.760 (-24 percent; Table 5). Reference pools G/G and G/A. Genotype frequencies G/G 14, G/A 47, A/A 23 (Table 4).",
      source_name = "COMT"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 84L,
    n_studies = 1L,
    age_range = "19-49 years",
    age_mean = "24.5 years (SD 3.9)",
    weight_mean = "68.0 kg (SD 11.7)",
    sex_female_pct = 64.3,
    race_ethnicity = "Western European ancestry (all participants)",
    disease_state = "Healthy volunteers (BMI 18-30 kg/m2; no acute or chronic disease; no regular medication except oral contraceptives).",
    dose_range = "Green tea extract capsules containing 18 mg EGC (standardised to 150 mg EGCG) taken orally twice daily with a meal for 4 days, then one capsule after an overnight fast on day 5; PK sampled on day 5 pre-dose and 0.5-9 h post-dose.",
    regions = "Germany (University of Ulm)",
    notes = "Scholl 2018 Table 1 (84 analysed of 89 enrolled; 30 male, 54 female, 37 of the women on oral contraceptives). NCT01360320. The three catechin models (EGCG, EGC, ECG) were developed independently on the same subjects."
  )

  ini({
    # Structural parameters -- Scholl 2018 Table 3, EGC column (final
    # estimates). Structure from the S6 File (Modeling Supplement) EGC NONMEM
    # control stream: ADVAN6 with explicit $DES, bolus dose into the depot
    # (RATE dropped), first-order K12 to central, lag time and F1 on the depot.
    lka <- log(0.830); label("First-order absorption rate constant ka, depot to central (1/h)") # Table 3 EGC ka = 0.830 1/h (S6 File K12 = THETA(1), no ETA)
    lcl <- log(560); label("Apparent clearance CL/F, SLCO1B1 rs2306283 G carriers and UGT1A1*28 carriers (L/h)") # Table 3 EGC CL/F = 560 L/h
    lvc <- log(1180); label("Apparent central volume Vcentral/F, no contraceptive and ABCC2 rs717620 C/C (L)") # Table 3 EGC Vcentral/F = 1180 L
    lq <- log(324); label("Apparent intercompartmental clearance Q/F, ABCC2 rs3740066 T carriers (L/h)") # Table 3 EGC Q/F = 324 L/h
    lvp <- log(10100); label("Apparent peripheral volume Vperipheral/F (L)") # Table 3 EGC Vperipheral/F = 10100 L
    ltlag <- log(0.413); label("Absorption lag time (h)") # Table 3 EGC lag time = 0.413 h
    lfdepot <- fixed(log(1)); label("Relative bioavailability F1, COMT rs4680 G carriers (unitless)") # Table 3 EGC F = 1.00 fixed; S6 File $THETA (1 FIX) ;F

    # Covariate effects -- Table 3 EGC rows; implementation from the S6 File
    # EGC $PK block.
    e_conmed_birthcontrol_vc <- 0.450; label("Vcentral/F multiplier with oral contraceptive use (unitless)") # Table 3 CONCEP = 0.450; S6 File IF (CONTRACEPTIVES.EQ.1) CONCEP=THETA(8)
    e_snp_abcc2_rs717620_vc <- -0.370; label("Fractional change in Vcentral/F per ABCC2 rs717620 T allele (unitless)") # Table 3 MRP2_rs717620_EGC = -0.370; S6 File V2 = ...*(1+MRP71EFF*MRP71T)
    e_snp_abcc2_rs3740066_cc_q <- 2.33; label("Q/F multiplier for ABCC2 rs3740066 C/C wild-type carriers (unitless)") # Table 3 MRP2_rs3740066 = 2.33; S6 File IF (MRP37T.EQ.0) MRP37EFF=THETA(10)
    e_snp_slco1b1_rs2306283_aa_cl <- 0.653; label("CL/F multiplier for SLCO1B1 rs2306283 A/A wild-type carriers (unitless)") # Table 3 OATP1B1_rs2306283 = 0.653; S6 File IF (OATP23T.EQ.0) OATP23EFF=THETA(9)
    e_snp_ugt1a1_rs8175347_wt_cl <- 0.744; label("CL/F multiplier for UGT1A1 (TA)6/(TA)6 wild-type carriers (unitless)") # Table 3 UGT1A1*28 = 0.744; S6 File IF (UGT.EQ.0) UGTEFF=THETA(12)
    e_snp_comt_rs4680_hom_fdepot <- 0.760; label("F1 multiplier for COMT rs4680 A/A carriers (unitless)") # Table 3 COMT_rs4680 = 0.760; S6 File IF (COMT.EQ.2) COMTEFF=THETA(13)

    # IIV -- Table 3 EGC IIV rows (%CV), omega^2 = (CV/100)^2 (see the EGCG
    # model and the vignette for the convention). Diagonal $OMEGA in the S6
    # File EGC stream.
    etalcl ~ 0.1444 # Table 3 EGC IIV CL/F 38 %CV -> 0.38^2
    etalvc ~ 0.64 # Table 3 EGC IIV Vcentral/F 80 %CV -> 0.80^2
    etalq ~ 0.8464 # Table 3 EGC IIV Q/F 92 %CV -> 0.92^2
    etalfdepot ~ 0.1369 # Table 3 EGC IIV F 37 %CV -> 0.37^2
    etaltlag ~ 0.0121 # Table 3 EGC IIV lag time 11 %CV -> 0.11^2

    # Residual error -- S6 File $ERROR Y = IPRED*(1+EPS(1)).
    propSd <- 0.28; label("Proportional residual error (fraction)") # Table 3 EGC prop. residual variability 28 %
  })

  model({
    # Genotype strata as coded in the S6 File EGC stream.
    mrp71t <- SNP_ABCC2_RS717620_HET + 2 * SNP_ABCC2_RS717620_HOM
    mrp37_cc <- 1 - SNP_ABCC2_RS3740066_HET - SNP_ABCC2_RS3740066_HOM
    oatp23_aa <- 1 - SNP_SLCO1B1_RS2306283_HET - SNP_SLCO1B1_RS2306283_HOM
    ugt_wt <- 1 - SNP_UGT1A1_RS8175347_HET - SNP_UGT1A1_RS8175347_HOM

    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * e_snp_slco1b1_rs2306283_aa_cl^oatp23_aa * e_snp_ugt1a1_rs8175347_wt_cl^ugt_wt
    vc <- exp(lvc + etalvc) * e_conmed_birthcontrol_vc^CONMED_BIRTHCONTROL * (1 + e_snp_abcc2_rs717620_vc * mrp71t)
    q <- exp(lq + etalq) * e_snp_abcc2_rs3740066_cc_q^mrp37_cc
    vp <- exp(lvp)
    tlag <- exp(ltlag + etaltlag)
    fdepot <- exp(lfdepot + etalfdepot) * e_snp_comt_rs4680_hom_fdepot^SNP_COMT_RS4680_HOM

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    alag(depot) <- tlag
    f(depot) <- fdepot

    # Amounts in mg, volumes in L; S2 = V2/1000 gives ug/L.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
