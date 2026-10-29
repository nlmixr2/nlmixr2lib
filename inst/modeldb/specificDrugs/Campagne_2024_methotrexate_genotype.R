Campagne_2024_methotrexate_genotype <- function() {
  description <- "Joint parent-metabolite population PK model for high-dose intravenous methotrexate (MTX) and 7-hydroxymethotrexate (7OHMTX) in genotyped infants and young children (0.09-4.94 years) with primary brain tumors (Campagne 2024, final BSA-eGFR-SNP covariate model). Same structure as Campagne_2024_methotrexate (two-compartment MTX with non-metabolic and 7OHMTX-forming metabolic clearances; two-compartment 7OHMTX whose central volume equals the MTX central volume; BSA-normalized dosing and parameters; BSA and eGFR power effects; inter-individual and inter-occasion variability), re-estimated in the 117 genotyped patients and adding dominant-model carrier effects of MTHFR rs1801131 on MTX non-metabolic clearance and central volume, ABCC2 rs8187710 on MTX central volume, and SLC19A1 rs4818789 on MTX metabolic clearance."
  reference <- "Campagne O, Huang J, Lin T, Reddick WE, Selvo NS, Onar-Thomas A, Ward D, Robinson G, Gajjar A, Stewart CF. Population pharmacokinetics of methotrexate and 7-hydroxymethotrexate and delayed excretion in infants and young children with brain tumors. Eur J Pharm Sci. 2024;193:106669. doi:10.1016/j.ejps.2023.106669"
  vignette <- "Campagne_2024_methotrexate"
  units <- list(time = "h", dosing = "mg/m^2", concentration = "umol/L")

  covariateData <- list(
    BSA = list(
      description = "Body surface area, baseline value before each high-dose MTX cycle",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect centered on 0.55 m^2, the median of the 117 genotyped patients (Table 1, pharmacogenetic analysis column; Methods: 'Continuous covariates were evaluated as a power model centered on the median covariate value'). Enters MTX non-metabolic clearance, MTX central volume, MTX metabolic clearance and 7OHMTX clearance, on top of the BSA normalization already carried by the mg/m^2 dose input. Baseline value of each occasion (induction cycle). Cohort range 0.25-0.83 m^2.",
      source_name = "BSA"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate, baseline value before each high-dose MTX cycle",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Estimated with the full equation cited by the paper (Millisor et al. 2017). Power effect centered on 114.3 mL/min/1.73 m^2, the median of the genotyped patients (Table 1, pharmacogenetic analysis column). Enters MTX non-metabolic clearance and 7OHMTX clearance. Baseline value of each occasion (induction cycle). Cohort range 56.4-244.8 mL/min/1.73 m^2.",
      source_name = "eGFR"
    ),
    OCC = list(
      description = "Occasion index = induction therapy cycle (1-4) for inter-occasion variability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "SJYC07 induction therapy delivered high-dose MTX on day 1 of each of four 28-day cycles, and pharmacokinetic sampling was performed in each cycle (Table 1 'Induction therapy cycle'). Decomposed inside model() into the indicators occ1..occ4 that select the per-cycle IOV etas on MTX non-metabolic clearance, MTX central volume, MTX intercompartmental clearance, MTX metabolic clearance and 7OHMTX clearance. Records with OCC outside 1-4 receive no IOV.",
      source_name = "cycle"
    ),
    SNP_MTHFR_RS1801131 = list(
      description = "MTHFR rs1801131 variant-carrier indicator (1 = heterozygous or homozygous variant, 0 = wild-type)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (wild-type)",
      notes = "Dominant inheritance model (Table 2 footnote: 'wild-type (reference) versus combined heterozygous and homozygous mutated (variant)'). Supplementary Table S1: 61 WT (52.1%), 48 HE (41.0%), 8 HOM (6.8%) of 117, so 56 carriers. Lowers MTX non-metabolic clearance and central volume.",
      source_name = "MTHFR"
    ),
    SNP_ABCC2_RS8187710 = list(
      description = "ABCC2 rs8187710 variant-carrier indicator (1 = heterozygous or homozygous variant, 0 = wild-type)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (wild-type)",
      notes = "Dominant inheritance model (Table 2 footnote). Supplementary Table S1: 101 WT (86.3%), 15 HE (12.8%), 1 HOM (0.9%) of 117, so 16 carriers. Lowers MTX central volume.",
      source_name = "ABCC2"
    ),
    SNP_SLC19A1_RS4818789 = list(
      description = "SLC19A1 rs4818789 variant-carrier indicator (1 = heterozygous or homozygous variant, 0 = wild-type)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (wild-type)",
      notes = "Dominant inheritance model (Table 2 footnote). Supplementary Table S1: 67 WT (57.3%), 40 HE (34.2%), 10 HOM (8.5%) of 117, so 50 carriers. Raises MTX metabolic (7OHMTX-forming) clearance.",
      source_name = "SLC19A1"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened in the first covariate analysis; highly correlated with BSA and explained no further variability once BSA was in the model (Results; Supplementary Fig. S5). Genotyped-cohort median 1.97 years, range 0.09-4.94 (Table 1)."
    ),
    WT = list(
      description = "Actual body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened; not retained. Genotyped-cohort median 11.8 kg, range 3.7-21.2 (Table 1)."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not retained. 56 of 117 genotyped patients female (Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = "Screened; not retained. Genotyped-cohort median 4.1, range 2.8-4.8; Table 1 prints the unit as U/L, which is a typesetting error for g/dL given the range."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened; not retained (eGFR was retained instead). Genotyped-cohort median 0.2 mg/dL, range 0.1-0.45 (Table 1)."
    )
  )

  compartmentData <- list(
    central = list(analyte = "methotrexate", units = "mg/m^2", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "methotrexate", units = "mg/m^2", specimen = "tissue", verified = TRUE),
    central_7ohmtx = list(
      analyte = "7-hydroxymethotrexate (as methotrexate mass equivalents)",
      units = "mg/m^2",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_7ohmtx = list(
      analyte = "7-hydroxymethotrexate (as methotrexate mass equivalents)",
      units = "mg/m^2",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 117,
    n_studies = 1,
    n_courses = 424,
    age_range = "0.09-4.94 years",
    age_median = "1.97 years",
    weight_range = "3.7-21.2 kg",
    weight_median = "11.8 kg",
    bsa_range = "0.25-0.83 m^2",
    bsa_median = "0.55 m^2",
    sex_female_pct = 47.9,
    renal_function = "eGFR median 114.3 mL/min/1.73 m^2 (range 56.4-244.8); serum creatinine median 0.2 mg/dL (range 0.1-0.45)",
    genotypes = "Illumina Omni2.5Exome-8 genome-wide genotyping; 22 SNPs in 9 MTX/folate transport and metabolism genes screened (Supplementary Table S1). Variant carriers: MTHFR rs1801131 56/117, ABCC2 rs8187710 16/117, SLC19A1 rs4818789 50/117.",
    disease_state = "Genotyped subset (117 of 142) of infants and young children (<= 5 years) with primary brain tumors after tumor resection, enrolled in the St. Jude SJYC07 trial (NCT00602667).",
    dose_range = "High-dose MTX 5 g/m^2 (2.5 g/m^2 for infants <= 31 days at enrollment) as a 24-h infusion on day 1 of each of four 28-day induction cycles: 10 percent of the dose over 1 h, the remaining 90 percent over 23 h. Hydration, urine alkalinization and leucovorin rescue from 42 h.",
    regions = "United States (St. Jude Children's Research Hospital, Memphis, TN)",
    notes = "Baseline demographics in Campagne 2024 Table 1 (pharmacogenetic analysis column); number of courses is the sum of the per-cycle counts in that column (109 + 110 + 103 + 102). MTX and 7OHMTX measured by LC-MS/MS. Fitted in Monolix 2020R1 (SAEM) with M3 handling of BLQ data."
  )

  ini({
    # Structural parameters at the median BSA (0.55 m^2) and median eGFR
    # (114.3 mL/min/1.73 m^2) of the genotyped cohort, for wild-type patients
    # at all three loci. Campagne 2024 Table 2 'BSA, eGFR, SNP model (N = 117)'
    # column. Every clearance and volume is per m^2 of BSA because the model
    # input was the BSA-normalized MTX dose (mg/m^2).
    lcl_nonmet <- log(4.9); label("MTX non-metabolic clearance CL (L/h/m^2)") # Table 2 'MTX non-metabolic clearance CL (L/h/m2)' = 4.9 (RSE 2.5%)
    lvc <- log(20.8); label("MTX central volume V (L/m^2); also the 7OHMTX central volume") # Table 2 'MTX volume, V (L/m2)' = 20.8 (RSE 1.8%)
    lq <- log(0.048); label("MTX intercompartmental clearance Q (L/h/m^2)") # Table 2 'MTX peripheral clearance, Q (L/h/m2)' = 0.048 (RSE 5.6%)
    lvp <- log(1.7); label("MTX peripheral volume Vp (L/m^2)") # Table 2 'MTX peripheral volume, Vp (L/m2)' = 1.7 (RSE 6.8%)
    lcl_met <- log(1.1); label("MTX metabolic clearance forming 7OHMTX, CLMET (L/h/m^2)") # Table 2 'MTX metabolic clearance, CLMET' = 1.1 (RSE 3.4%)
    lcl_7ohmtx <- log(3.0); label("7OHMTX clearance CL7OH (L/h/m^2)") # Table 2 '7OHMTX clearance, CL7OH (L/h/m2)' = 3.0 (RSE 1.8%)
    lq_7ohmtx <- log(0.37); label("7OHMTX intercompartmental clearance Q7OH (L/h/m^2)") # Table 2 '7OHMTX peripheral clearance, Q7OH (L/h/m2)' = 0.37 (RSE 5.9%)
    lvp_7ohmtx <- log(9.2); label("7OHMTX peripheral volume Vp,7OH (L/m^2)") # Table 2 '7OHMTX peripheral volume, Vp,7OH (L/m2)' = 9.2 (RSE 7.4%)

    # Covariate exponents (power model centered on the median covariate value;
    # Results final covariate equations).
    e_bsa_cl_nonmet <- 0.43; label("Power exponent of BSA/0.55 on CL (unitless)") # Table 2 'BSA effect, thetaBSA,CL' = 0.43 (RSE 19.5%)
    e_crcl_cl_nonmet <- 0.13; label("Power exponent of eGFR/114.3 on CL (unitless)") # Table 2 'eGFR effect, thetaeGFR,CL' = 0.13 (RSE 36.2%)
    e_bsa_vc <- 0.46; label("Power exponent of BSA/0.55 on V (unitless)") # Table 2 'BSA effect on V, thetaBSA,V' = 0.46 (RSE 12.7%)
    e_bsa_cl_met <- 1.2; label("Power exponent of BSA/0.55 on CLMET (unitless)") # Table 2 'BSA effect, thetaBSA,CLmet' = 1.2 (RSE 9.9%)
    e_bsa_cl_7ohmtx <- 0.61; label("Power exponent of BSA/0.55 on CL7OH (unitless)") # Table 2 'BSA effect, thetaBSA,CL7oh' = 0.61 (RSE 13.4%)
    e_crcl_cl_7ohmtx <- 0.18; label("Power exponent of eGFR/114.3 on CL7OH (unitless)") # Table 2 'eGFR effect, thetaeGFR,CL7oh' = 0.18 (RSE 33.1%)

    # Genotype effects (dominant model, carrier vs wild-type). Methods:
    # 'Categorical covariates were modeled using an exponential model', so a
    # carrier's parameter is multiplied by exp(theta). The printed equations
    # write (theta)^SNP, which with a negative theta would give a negative
    # clearance for carriers; see the vignette's Assumptions and deviations.
    e_snp_mthfr_rs1801131_cl_nonmet <- -0.083; label("Log-scale effect of MTHFR rs1801131 carrier status on CL (unitless)") # Table 2 'MTHFR (rs1801131) variant, thetaMTHFR,CL' = -0.083 (RSE 43.1%)
    e_snp_mthfr_rs1801131_vc <- -0.076; label("Log-scale effect of MTHFR rs1801131 carrier status on V (unitless)") # Table 2 'MTHFR (rs1801131) variant, thetaMTHFR,V' = -0.076 (RSE 32.3%)
    e_snp_abcc2_rs8187710_vc <- -0.10; label("Log-scale effect of ABCC2 rs8187710 carrier status on V (unitless)") # Table 2 'ABCC2 (rs8187710) variant, thetaABCC2' = -0.10 (RSE 34.9%)
    e_snp_slc19a1_rs4818789_cl_met <- 0.099; label("Log-scale effect of SLC19A1 rs4818789 carrier status on CLMET (unitless)") # Table 2 'SLC19A1 (rs4818789) variant, thetaSLC' = 0.099 (RSE 49.1%)

    # Inter-individual variability. Table 2 footnote: 'Inter-individual and
    # inter-occasion variabilities are reported as standard deviation', so
    # each variance below is the reported SD squared. No correlations were
    # estimated (Results).
    etalcl_nonmet ~ 0.0289 # Table 2 'IIV-CL' SD 0.17 -> 0.17^2
    etalvc ~ 0.0121 # Table 2 'IIV-V' SD 0.11 -> 0.11^2
    etalq ~ 0.2209 # Table 2 'IIV-Q' SD 0.47 -> 0.47^2
    etalvp ~ 0.4624 # Table 2 'IIV-V2' (MTX peripheral volume) SD 0.68 -> 0.68^2
    etalcl_met ~ 0.04 # Table 2 'IIV-CLMET' SD 0.20 -> 0.20^2
    etalcl_7ohmtx ~ 0.0121 # Table 2 'IIV-CL7OH' SD 0.11 -> 0.11^2
    etalq_7ohmtx ~ 0.2116 # Table 2 'IIV-Q7OH' SD 0.46 -> 0.46^2
    etalvp_7ohmtx ~ 0.4761 # Table 2 'IIV-Vp,7OH' SD 0.69 -> 0.69^2

    # Inter-occasion variability (occasion = induction cycle 1-4), reported as
    # SD in Table 2 and squared here. One variance per parameter is estimated;
    # occasions 2-4 share occasion 1's variance. Implemented by the
    # occasion-indicator expansion rather than rxode2's `eta ~ var | OCC`
    # syntax, which parses but cannot be simulated from an rxUi.
    etaiov_cl_nonmet_1 ~ 0.0144 # Table 2 'IOV-CL' SD 0.12 -> 0.12^2
    etaiov_cl_nonmet_2 ~ fixed(0.0144) # same variance as occasion 1
    etaiov_cl_nonmet_3 ~ fixed(0.0144) # same variance as occasion 1
    etaiov_cl_nonmet_4 ~ fixed(0.0144) # same variance as occasion 1
    etaiov_vc_1 ~ 0.002704 # Table 2 'IOV-V' SD 0.052 -> 0.052^2
    etaiov_vc_2 ~ fixed(0.002704) # same variance as occasion 1
    etaiov_vc_3 ~ fixed(0.002704) # same variance as occasion 1
    etaiov_vc_4 ~ fixed(0.002704) # same variance as occasion 1
    etaiov_q_1 ~ 0.1296 # Table 2 'IOV-Q' SD 0.36 -> 0.36^2
    etaiov_q_2 ~ fixed(0.1296) # same variance as occasion 1
    etaiov_q_3 ~ fixed(0.1296) # same variance as occasion 1
    etaiov_q_4 ~ fixed(0.1296) # same variance as occasion 1
    etaiov_cl_met_1 ~ 0.04 # Table 2 'IOV-CLMET' SD 0.20 -> 0.20^2
    etaiov_cl_met_2 ~ fixed(0.04) # same variance as occasion 1
    etaiov_cl_met_3 ~ fixed(0.04) # same variance as occasion 1
    etaiov_cl_met_4 ~ fixed(0.04) # same variance as occasion 1
    etaiov_cl_7ohmtx_1 ~ 0.0289 # Table 2 'IOV-CL7OH' SD 0.17 -> 0.17^2
    etaiov_cl_7ohmtx_2 ~ fixed(0.0289) # same variance as occasion 1
    etaiov_cl_7ohmtx_3 ~ fixed(0.0289) # same variance as occasion 1
    etaiov_cl_7ohmtx_4 ~ fixed(0.0289) # same variance as occasion 1

    # Residual error: proportional, with the additive component fixed to
    # 0.001 uM (Table 2 'Proportional residual error' rows and footnote 'In
    # both final models, the additive residual error was fixed to 0.001 uM').
    propSd <- 0.21; label("MTX proportional residual error (fraction)") # Table 2 'Proportional residual error, MTX' = 0.21 (RSE 2.3%)
    addSd <- fixed(0.001); label("MTX additive residual error (umol/L)") # Table 2 footnote: additive residual error fixed to 0.001 uM
    propSd_7ohmtx <- 0.24; label("7OHMTX proportional residual error (fraction)") # Table 2 'Proportional residual error, 7OHMTX' = 0.24 (RSE 2.2%)
    addSd_7ohmtx <- fixed(0.001); label("7OHMTX additive residual error (umol/L)") # Table 2 footnote: additive residual error fixed to 0.001 uM
  })

  model({
    # Inter-occasion variability, selected by the induction cycle OCC.
    occ1 <- (OCC == 1)
    occ2 <- (OCC == 2)
    occ3 <- (OCC == 3)
    occ4 <- (OCC == 4)
    iov_cl_nonmet <- occ1 * etaiov_cl_nonmet_1 + occ2 * etaiov_cl_nonmet_2 +
      occ3 * etaiov_cl_nonmet_3 + occ4 * etaiov_cl_nonmet_4
    iov_vc <- occ1 * etaiov_vc_1 + occ2 * etaiov_vc_2 +
      occ3 * etaiov_vc_3 + occ4 * etaiov_vc_4
    iov_q <- occ1 * etaiov_q_1 + occ2 * etaiov_q_2 +
      occ3 * etaiov_q_3 + occ4 * etaiov_q_4
    iov_cl_met <- occ1 * etaiov_cl_met_1 + occ2 * etaiov_cl_met_2 +
      occ3 * etaiov_cl_met_3 + occ4 * etaiov_cl_met_4
    iov_cl_7ohmtx <- occ1 * etaiov_cl_7ohmtx_1 + occ2 * etaiov_cl_7ohmtx_2 +
      occ3 * etaiov_cl_7ohmtx_3 + occ4 * etaiov_cl_7ohmtx_4

    # Covariate ratios, centered on the Table 1 genotyped-cohort medians.
    bsa_ratio <- BSA / 0.55
    egfr_ratio <- CRCL / 114.3

    # Individual parameters (per m^2 of BSA). Results, final covariate
    # equations for CL, V, CLMET and CL7OH, with the genotype terms on the
    # exponential (log-additive) scale.
    cl_nonmet <- exp(
      lcl_nonmet + etalcl_nonmet + iov_cl_nonmet +
        e_snp_mthfr_rs1801131_cl_nonmet * SNP_MTHFR_RS1801131
    ) *
      bsa_ratio^e_bsa_cl_nonmet * egfr_ratio^e_crcl_cl_nonmet
    vc <- exp(
      lvc + etalvc + iov_vc +
        e_snp_mthfr_rs1801131_vc * SNP_MTHFR_RS1801131 +
        e_snp_abcc2_rs8187710_vc * SNP_ABCC2_RS8187710
    ) *
      bsa_ratio^e_bsa_vc
    q <- exp(lq + etalq + iov_q)
    vp <- exp(lvp + etalvp)
    cl_met <- exp(
      lcl_met + etalcl_met + iov_cl_met +
        e_snp_slc19a1_rs4818789_cl_met * SNP_SLC19A1_RS4818789
    ) *
      bsa_ratio^e_bsa_cl_met
    cl_7ohmtx <- exp(lcl_7ohmtx + etalcl_7ohmtx + iov_cl_7ohmtx) *
      bsa_ratio^e_bsa_cl_7ohmtx * egfr_ratio^e_crcl_cl_7ohmtx
    q_7ohmtx <- exp(lq_7ohmtx + etalq_7ohmtx)
    vp_7ohmtx <- exp(lvp_7ohmtx + etalvp_7ohmtx)
    # The 7OHMTX central volume was set equal to the MTX central volume
    # for identifiability (Methods; Supplementary Figure S3 'V7OH = V').
    vc_7ohmtx <- vc

    # Micro-constants
    kel <- cl_nonmet / vc
    kmet <- cl_met / vc
    k12 <- q / vc
    k21 <- q / vp
    kel_7ohmtx <- cl_7ohmtx / vc_7ohmtx
    k12_7ohmtx <- q_7ohmtx / vc_7ohmtx
    k21_7ohmtx <- q_7ohmtx / vp_7ohmtx

    # IV infusion of the BSA-normalized dose (mg/m^2) into central. The
    # metabolic flux enters the 7OHMTX central compartment mole for mole, so
    # the 7OHMTX states carry MTX mass equivalents.
    d/dt(central) <- -(kel + kmet) * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(central_7ohmtx) <- kmet * central - kel_7ohmtx * central_7ohmtx -
      k12_7ohmtx * central_7ohmtx + k21_7ohmtx * peripheral1_7ohmtx
    d/dt(peripheral1_7ohmtx) <- k12_7ohmtx * central_7ohmtx -
      k21_7ohmtx * peripheral1_7ohmtx

    # Concentrations in umol/L: (mg/m^2) / (L/m^2) = mg/L, converted with the
    # MTX molecular weight 454.44 g/mol (C20H22N8O5). The 7OHMTX states hold
    # MTX mass equivalents, so the same factor gives 7OHMTX in umol/L.
    Cc <- central / vc * 1000 / 454.44
    Cc_7ohmtx <- central_7ohmtx / vc_7ohmtx * 1000 / 454.44

    Cc ~ add(addSd) + prop(propSd) + combined1()
    Cc_7ohmtx ~ add(addSd_7ohmtx) + prop(propSd_7ohmtx) + combined1()
  })
}
