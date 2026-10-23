Agergaard_2025_tacrolimus <- function() {
  description <- paste(
    "Two-compartment population PK model for oral immediate-release tacrolimus in stable adult",
    "kidney transplant recipients (Agergaard 2025), fitted jointly to whole-blood and peripheral",
    "blood mononuclear cell (PBMC) concentrations. First-order absorption (ka = 1/MTT) after a lag",
    "time and linear elimination. The PBMC concentration is an algebraic observable, not an ODE",
    "state: distribution into PBMCs was found to be instantaneous, so Cpbmc = R_C:PBMC * Cc, where",
    "the PBMC:whole-blood ratio carries a Box-Cox-transformed IIV, fractional-shift effects of",
    "CYP3A5 non-expressor status and NR1I2 -25385C>T (rs3814055) homozygous-mutant genotype, and a",
    "0.709 dilution multiplier for PBMC samples visibly contaminated with erythrocytes. CYP3A5",
    "non-expressor status lowers CL/F, female sex lowers Vc/F, and relative bioavailability (fixed",
    "to 1 at age 59 years, with IIV) increases exponentially with age. Separate proportional",
    "residual errors for whole blood and PBMC.",
    sep = " "
  )
  reference <- paste(
    "Agergaard K, Thiesson HC, Carstens J, Staatz CE, Jarvinen E, Nielsen F, Christensen HD,",
    "Juhl-Sandberg R, Brosen K, Stage TB, Andersen DT, Kjellsson MC, Bergmann TK.",
    "Population pharmacokinetics of tacrolimus whole blood and peripheral blood mononuclear cell",
    "concentrations in stable kidney-transplanted patients. Br J Clin Pharmacol.",
    "2025;91(3):761-773. doi:10.1111/bcp.16277",
    sep = " "
  )
  vignette <- "Agergaard_2025_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "tacrolimus", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Exponential effect on relative bioavailability centred by subtraction at the cohort",
        "median of 59 years: F = exp(0.0123 * (AGE - 59)) (Table 2 footnote d; Table 1 median).",
        "The paper states only that continuous covariates were normalised to the median and tested",
        "as linear, exponential and power functions; the exponential form is the one that",
        "reproduces the typical-individual line of Figure 4A (about 0.63 at 21 years and 1.30 at",
        "81 years). See the vignette for the comparison against the linear and power readings."
      ),
      source_name = "Age"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Fractional shift on Vc/F: Vc/F = 180 * (1 - 0.345 * SEXF). Table 2 footnote b states the",
        "typical Vc/F is that of a male, and the Discussion states women have a 34.5% lower Vc/F."
      ),
      source_name = "Sex"
    ),
    CYP3A5_EXPR = list(
      description = "CYP3A5 expresser status: 1 = at least one CYP3A5*1 allele, 0 = CYP3A5*3/*3 non-expresser",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (expresser) is the Table 2 typical subject; see notes",
      notes = paste(
        "The paper parameterises its effects on the NON-expresser (Table 2 rows 'CYP3A5",
        "non-expressor'; footnotes a and f put expressers at the typical value), so the model uses",
        "(1 - CYP3A5_EXPR) as the non-expresser indicator, the value inversion the register entry",
        "prescribes. Effects: CL/F * (1 - 0.466) and R_C:PBMC * (1 + 0.424) in non-expressers.",
        "Recipient genotype (rs776746); 11 of 63 (17.5%) expressers (Table 1)."
      ),
      source_name = "CYP3A5*1/*3"
    ),
    SNP_NR1I2_RS3814055_HOM = list(
      description = "NR1I2 (PXR) -25385C>T (rs3814055) homozygous-mutant indicator: 1 = TT, 0 = CC or CT",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CC or CT)",
      notes = paste(
        "Recessive encoding: the authors pooled CC and CT because pooling heterozygotes with the",
        "homozygous mutants fitted worse (dOFV 9.73; Results 3.3). Fractional shift on R_C:PBMC",
        "only: R * (1 + 0.607 * SNP_NR1I2_RS3814055_HOM). Figure 5C shows no whole-blood effect.",
        "Genotype counts CC 31, CT 21, TT 11 (Table 1)."
      ),
      source_name = "NR1I2-25385C>T homozygous mutant"
    ),
    SAMPLE_RBC_CONTAM = list(
      description = "Per-observation flag for a PBMC sample whose pellet had a considerable red colour (erythrocyte contamination)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no significant erythrocyte contamination)",
      notes = paste(
        "The PBMC isolation omitted an erythrocyte-lysis step; pellets with a considerable red",
        "colour were marked (Methods 2.2). Applied as a direct multiplier on R_C:PBMC:",
        "R * 0.709 when SAMPLE_RBC_CONTAM = 1 (Results 3.2: red samples had '70.9% of R_C:PBMC",
        "without redness'; 23 such samples). Meaningful on PBMC observation records only; it does",
        "not change the whole-blood prediction."
      ),
      source_name = "Erythrocyte dilution"
    )
  )

  covariatesDataExcluded <- list(
    HCT = list(
      description = "Haematocrit",
      units = "L/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Significant on CL/F in forward inclusion but removed at backward elimination (Results 3.3). Cohort median 0.40 L/L (Table 1)."
    ),
    ALB = list(
      description = "Plasma albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL/F, Vc/F and R_C:PBMC; not retained (Methods 2.6). Cohort median 44 g/L (Table 1)."
    ),
    LBW = list(
      description = "Lean body weight (Boers formula)",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on all parameters and tested as a replacement for sex on Vc/F, which fitted worse (Discussion). Cohort median 62.4 kg (Table 1)."
    ),
    SNP_ABCB1_RS1128503 = list(
      description = "ABCB1 1236C>T (rs1128503)",
      units = "(binary)",
      type = "binary",
      reference_category = NULL,
      notes = "Screened on CL/F, Vc/F, R_C:PBMC and F; not retained (Methods 2.6)."
    ),
    SNP_ABCB1_RS2032582 = list(
      description = "ABCB1 2677G>T (rs2032582)",
      units = "(binary)",
      type = "binary",
      reference_category = NULL,
      notes = "Screened on CL/F, Vc/F, R_C:PBMC and F; not retained (Methods 2.6)."
    ),
    SNP_ABCB1_RS1045642 = list(
      description = "ABCB1 3435C>T (rs1045642)",
      units = "(binary)",
      type = "binary",
      reference_category = NULL,
      notes = "Screened on CL/F, Vc/F and R_C:PBMC; not retained (Methods 2.6)."
    ),
    SNP_CYP3A4_RS35599367 = list(
      description = "CYP3A4*22 (rs35599367) carrier",
      units = "(binary)",
      type = "binary",
      reference_category = NULL,
      notes = "Screened on CL/F, Vc/F and R_C:PBMC; not retained (Methods 2.6). 6 of 63 carriers (Table 1)."
    ),
    SNP_NR1I2_RS1523127 = list(
      description = "NR1I2 -24381A>T (rs1523127)",
      units = "(binary)",
      type = "binary",
      reference_category = NULL,
      notes = "Screened on CL/F, Vc/F and R_C:PBMC; not retained (Methods 2.6)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL/F; not retained (Methods 2.6). Cohort median 18 U/L (Table 1)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL/F; not retained (Methods 2.6). Cohort median 7 umol/L (Table 1)."
    ),
    CONMED_STEROID = list(
      description = "Concomitant glucocorticoid",
      units = "(binary)",
      type = "binary",
      reference_category = NULL,
      notes = "Screened on CL/F; not retained (Methods 2.6). 12 of 63 patients (Table 1)."
    ),
    CONMED_CCB = list(
      description = "Concomitant calcium-channel blocker",
      units = "(binary)",
      type = "binary",
      reference_category = NULL,
      notes = "Screened on CL/F; not retained (Methods 2.6). 29 of 63 patients (Table 1)."
    ),
    CONMED_PPI = list(
      description = "Concomitant proton-pump inhibitor",
      units = "(binary)",
      type = "binary",
      reference_category = NULL,
      notes = "Screened on CL/F; not retained (Methods 2.6). 8 of 63 patients (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 63L,
    n_studies = 1L,
    n_observations = "366 whole-blood and 348 PBMC concentrations",
    age_range = "IQR 48-68 years",
    age_median = "59 years",
    weight_range = "IQR 73.4-98.3 kg",
    weight_median = "86.8 kg",
    sex_female_pct = 28.6,
    race_ethnicity = c(Caucasian = 95.2, Other = 4.8),
    disease_state = "Stable adult kidney transplant recipients (>= 3 months post-transplant; median 38 months, IQR 15-79) on immediate-release tacrolimus b.i.d. plus mycophenolate b.i.d., at steady state with no dose change in the previous 5 days.",
    dose_range = "Oral immediate-release tacrolimus twice daily; median dose 2 mg (IQR 1.2-2.0 mg as printed in Table 1); whole-blood trough target generally 5 ng/mL.",
    regions = "Denmark (Odense, Esbjerg and Kolding nephrology departments), June 2021 - October 2023.",
    cyp3a5_distribution = "11 expressers (17.5%), 52 non-expressers (82.5%) (Table 1).",
    nr1i2_rs3814055_distribution = "CC 31 (49.2%), CT 21 (33.3%), TT 11 (17.5%) (Table 1).",
    sampling = "One occasion per patient: dense (n = 18; predose and 0.5, 1, 1.5, 2, 3, 4, 6, 8, 10, 12 h) or sparse (n = 45; predose and 0.5, 1.5, 4 h), with paired whole-blood and PBMC samples.",
    notes = "Agergaard 2025 Methods 2.1, Results 3.1 and Table 1. Tacrolimus measured by LC-MS/MS in both matrices (LOQ 0.20 ng/mL)."
  )

  ini({
    # Whole-blood disposition -- Table 2, 'Final model with covariates'. Apparent
    # parameters: F is fixed to 1 for a 59-year-old (footnotes c, d) so CL, Vc, Q
    # and Vp are all relative to the bioavailable fraction.
    lcl <- log(28.2); label("Apparent clearance CL/F of a CYP3A5 expresser (L/h)") # Table 2: CL/F = 28.2 L/h (RSE 5.73%), footnote a: CYP3A5 expresser
    lvc <- log(180); label("Apparent central volume Vc/F of a male (L)") # Table 2: Vc/F = 180 L (RSE 5.60%), footnote b: male
    lq <- log(51.0); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2: Q/F = 51.0 L/h (RSE 8.38%)
    lvp <- log(1688); label("Apparent peripheral volume Vp/F (L)") # Table 2: Vp/F = 1688 L (RSE 19.3%)
    lfdepot <- fixed(log(1)); label("Relative bioavailability at age 59 years (fraction)") # Table 2: F = 1 fixed (footnotes c, d), fixed to allow IIV on F
    lmtt <- log(0.676); label("Mean absorption transit time MTT; ka = 1/MTT (h)") # Table 2: MTT = 0.676 h (RSE 15.7%)
    ltlag <- log(0.369); label("Absorption lag time (h)") # Table 2: Lag time = 0.369 h (RSE 5.00%)

    # Covariate effects on whole-blood PK -- Table 2 with the Discussion's
    # percentages confirming the fractional-shift form.
    e_cyp3a5_nonexpr_cl <- -0.466; label("Fractional change in CL/F for a CYP3A5 non-expresser (unitless)") # Table 2: CYP3A5 non-expressor on CL/F = -0.466 (RSE 10.6%); Discussion '46.6% lower CL/F'
    e_sexf_vc <- -0.345; label("Fractional change in Vc/F for a female (unitless)") # Table 2: Female on Vc/F = -0.345 (RSE 20.0%); Discussion '34.5% lower apparent oral Vc'
    e_age_fdepot <- 0.0123; label("Exponential age coefficient on F, centred at 59 years (1/year)") # Table 2: Age on F = 0.0123 (RSE 25.8%); form from Figure 4A

    # PBMC layer (Equation 2): Cpbmc = R_C:PBMC * A_C / (Vc/F).
    lratio_pbmc_central <- log(1.96); label("PBMC:whole-blood concentration ratio R_C:PBMC, CYP3A5 expresser with >= 1 NR1I2-25385 C allele (unitless)") # Table 2: R_C:PBMC = 1.96 (RSE 4.97%), footnote f
    e_cyp3a5_nonexpr_ratio_pbmc_central <- 0.424; label("Fractional change in R_C:PBMC for a CYP3A5 non-expresser (unitless)") # Table 2: CYP3A5 non-expressor on R_C:PBMC = 0.424 (RSE 20.8%); Discussion '42.4% higher'
    e_nr1i2_hom_ratio_pbmc_central <- 0.607; label("Fractional change in R_C:PBMC for NR1I2-25385 TT (unitless)") # Table 2: NR1I2-25385C>T homozygous mutant = 0.607 (RSE 23.8%); Discussion '60.7% higher'
    e_sample_rbc_contam_ratio_pbmc_central <- 0.709; label("Multiplier on R_C:PBMC for an erythrocyte-contaminated PBMC sample (unitless)") # Table 2: Erythrocyte dilution = 0.709 (RSE 8.58%); Results 3.2 '70.9% of R_C:PBMC without redness'
    boxcox_lratio_pbmc_central <- -2.50; label("Box-Cox shape parameter for the IIV on R_C:PBMC (unitless)") # Table 2: lambda = -2.50 (RSE 26.7%)

    # IIV -- Table 2 'IIV (%)' read as 100 * the eta standard deviation, so the
    # variance is (pct/100)^2. The vignette compares the alternative
    # arithmetic-CV reading omega^2 = log(1 + (pct/100)^2).
    etalcl ~ 0.0441 # Table 2: IIV CL/F = 21.0% (RSE 39.9%); 0.210^2
    etalfdepot ~ 0.1681 # Table 2: IIV F = 41.0% (RSE 15.7%); 0.410^2
    etalmtt ~ 1.6641 # Table 2: IIV MTT = 129% (RSE 13.7%); 1.29^2
    etalratio_pbmc_central ~ 0.142884 # Table 2: IIV R_C:PBMC = 37.8% (RSE 25.3%); 0.378^2, Box-Cox transformed

    # Residual error -- proportional in each matrix.
    propSd <- 0.137; label("Proportional residual SD, whole blood (fraction)") # Table 2: RUV_C = 13.7% (RSE 9.65%)
    propSd_Cpbmc <- 0.354; label("Proportional residual SD, PBMC (fraction)") # Table 2: RUV_PBMC = 35.4% (RSE 6.28%)
  })

  model({
    # Indicators. The paper parameterises its CYP3A5 effects on the
    # non-expresser, so invert the canonical expresser column.
    cyp3a5_nonexpr <- 1 - CYP3A5_EXPR

    # Whole-blood PK
    cl <- exp(lcl + etalcl) * (1 + e_cyp3a5_nonexpr_cl * cyp3a5_nonexpr)
    vc <- exp(lvc) * (1 + e_sexf_vc * SEXF)
    q <- exp(lq)
    vp <- exp(lvp)
    mtt <- exp(lmtt + etalmtt)
    ka <- 1 / mtt
    tlag <- exp(ltlag)
    fdepot <- exp(lfdepot + etalfdepot) * exp(e_age_fdepot * (AGE - 59))

    # PBMC:whole-blood ratio with a Box-Cox-transformed eta (Petersson 2009
    # form); lambda < 0 gives a long lower tail of R.
    bcLambda <- boxcox_lratio_pbmc_central
    phi_lratio_pbmc_central <- ((exp(etalratio_pbmc_central))^bcLambda - 1) / bcLambda
    ratio_pbmc_central <- exp(lratio_pbmc_central + phi_lratio_pbmc_central) *
      (1 + e_cyp3a5_nonexpr_ratio_pbmc_central * cyp3a5_nonexpr) *
      (1 + e_nr1i2_hom_ratio_pbmc_central * SNP_NR1I2_RS3814055_HOM) *
      (1 - SAMPLE_RBC_CONTAM + e_sample_rbc_contam_ratio_pbmc_central * SAMPLE_RBC_CONTAM)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    alag(depot) <- tlag

    # mg / L -> ng/mL
    Cc <- central / vc * 1000
    # Equation 2: distribution into PBMCs is instantaneous.
    Cpbmc <- ratio_pbmc_central * Cc

    Cc ~ prop(propSd)
    Cpbmc ~ prop(propSd_Cpbmc)
  })
}
