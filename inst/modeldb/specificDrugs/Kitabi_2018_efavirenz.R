Kitabi_2018_efavirenz <- function() {
  description <- "Two-compartment population PK model with first-order absorption for oral efavirenz 600 mg once daily at steady state in Tanzanian adults with HIV (efavirenz alone, week 16) and HIV/tuberculosis (week 20 during rifampicin-based anti-TB co-treatment and week 32, 8 weeks after it ended) (Kitabi 2018). Apparent clearance CL/F scales allometrically with weight (exponent 0.75, 70 kg), carries CYP2B6*6 (c.516G>T) genotype fractions, genotype-stratified fractions for the HIV/TB arm during and after co-treatment, and a sex factor; V/F and VP/F scale linearly and Q/F with exponent 0.75 on weight. Correlated IIV on CL/F and V/F, IIV on ka and Q/F, and between-occasion variability on CL/F in the HIV/TB arm."
  reference <- paste(
    "Kitabi EN, Minzi OMS, Mugusi S, Sasi P, Janabi M, Mugusi F,",
    "Bertilsson L, Burhenne J, Aklillu E. (2018). Long-term efavirenz",
    "pharmacokinetics is comparable between Tanzanian HIV and",
    "HIV/Tuberculosis patients with the same CYP2B6*6 genotype.",
    "Sci Rep 8:16316. doi:10.1038/s41598-018-34674-3.",
    sep = " "
  )
  vignette <- "Kitabi_2018_efavirenz"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling on a 70 kg reference with exponents hard-coded in the supplement control stream: 0.75 on CL/F and Q/F, 1 on V/F and VP/F (`0.75*(LOG(WT)-LOG(70))` in LTVCL and MU_4, `(LOG(WT)-LOG(70))` in LTVV and MU_5). Methods: 'The influence of weight on efavirenz clearance was modelled using a previous reported allometric scaling power model, with the allometric exponent fixed to 0.75'. Cohort weight median 48 kg (IQR 43-60); Table 1 mean (SD) 45.5 (13.15) kg in arm 1 and 55.14 (13.45) kg in arm 2.",
      source_name = "WT"
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (female)",
      notes = "The control stream codes `MALE = 1 IF(SEX.EQ.2)` and multiplies CL/F by THETA(18)^MALE, so the reference is female and male = 1 - SEXF. Sex was kept in the final model by the full-covariate approach for confounding control, not because it was significant (Results: 'male subjects had comparable clearance (approximately 1.02 folds) to females'). Cohort: 28 of 56 female (Table 1: 15/20 in arm 1, 13/36 in arm 2).",
      source_name = "SEX"
    ),
    SNP_CYP2B6_RS3745274_T_COUNT = list(
      description = "Count of CYP2B6 c.516G>T (rs3745274, CYP2B6*6) T-alleles per subject (0/1/2). 0 = CYP2B6*1/*1, 1 = CYP2B6*1/*6, 2 = CYP2B6*6/*6.",
      units = "(count, 0/1/2 alleles per subject)",
      type = "count",
      reference_category = "0 (CYP2B6*1/*1)",
      notes = "Kitabi 2018 genotypes 'CYP2B6 c.516 G > T(*6)' (Table 1) and codes it CP2B = 1/2/3 for *1/*1, *1/*6, *6/*6 in the control stream; decomposed in model() into the B1/B2/B3 indicators of Methods equation (3). Genotype was missing for 19 of 56 subjects; the authors imputed it with a NONMEM $MIX mixture model (P(*1/*1) = THETA(16) = 0.784, P(*1/*6) = (1 - P(*1/*1)) * THETA(17) with THETA(17) fixed to 0, P(*6/*6) = 0.216; Table 2 rows 'pi *1/*1' 0.784 and 'delta *1/*6' 0 fixed). That imputation is not encoded here: supply a genotype for every subject, or draw a missing one from those probabilities. Cohort (Table 1): *1/*1 18, *1/*6 13, *6/*6 6, missing 19.",
      source_name = "CP2B"
    ),
    TB_POS = list(
      description = "HIV/tuberculosis co-infected patient indicator (1 = study arm 2, HIV/TB; 0 = arm 1, HIV only)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (HIV only, arm 1)",
      notes = "The A2 indicator of Methods equation (3) and the control stream (`A2 = 1 IF(ARM.EQ.2)`). Time-fixed per subject. Arm 2 patients started 6 months of rifampicin/isoniazid/pyrazinamide/ethambutol anti-TB therapy and began efavirenz 4 weeks later; at OCC = 1 (week 20) they were on anti-TB co-treatment, at OCC = 2 (week 32) they had stopped it 8 weeks earlier. The arm-2 CL/F fractions therefore combine TB disease and anti-TB co-medication (rifampicin and isoniazid) and cannot be separated in this design.",
      source_name = "ARM"
    ),
    OCC = list(
      description = "Integer-valued occasion indicator (1 or 2).",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "OCC = 1 for every arm-1 record (week 16 of efavirenz) and for arm-2 records at week 20 (during anti-TB co-treatment); OCC = 2 for arm-2 records at week 32 (8 weeks after anti-TB co-treatment ended). OCC = 2 together with TB_POS = 1 forms the O2 indicator of Methods equation (3); TB_POS = 1 with OCC = 1 or 2 selects the between-occasion eta (control stream GROUP 21 / 22). Arm-1 subjects have a single occasion and no between-occasion variability.",
      source_name = "OCC"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "efavirenz", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "efavirenz", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "efavirenz", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 56L,
    n_studies = 1L,
    n_observations = 556L,
    age_range = "IQR 36-50 years (median 42)",
    age_median = "42 years",
    weight_range = "IQR 43-60 kg (median 48)",
    weight_median = "48 kg",
    sex_female_pct = 50,
    race_ethnicity = "Black African (Tanzanian)",
    disease_state = "Treatment-naive HIV-1 infection with CD4 < 200 cells/uL (arm 1, n = 20), or HIV with active tuberculosis on rifampicin-based anti-TB treatment (arm 2, n = 36)",
    dose_range = "600 mg efavirenz orally once daily (evening), with lamivudine plus zidovudine, stavudine or tenofovir",
    regions = "Tanzania (Muhimbili National Hospital, Dar es Salaam)",
    notes = "Intensive steady-state sampling at about 0, 1, 3, 6, 12, 16 and 24 h after an evening dose: arm 1 at week 16 of efavirenz; arm 2 at week 20 of rifampicin-based anti-TB therapy (about 16 weeks of efavirenz) and again 8 weeks after anti-TB therapy ended (week 32). No concentrations were below the 10 ng/mL LLOQ. Demographics in Table 1. NONMEM 7.3, FOCE-I, ADVAN6."
  )

  ini({
    # Structural parameters for a 70 kg female, CYP2B6*1/*1, HIV-only (arm 1)
    # subject. The supplement control stream $THETA initial values equal the
    # Table 2 final estimates to rounding, so Table 2 is used throughout.
    lcl <- log(23.4); label("Apparent clearance CL/F for CYP2B6*1/*1, arm 1, 70 kg female (L/h)") # Table 2 'CL11' 23.4; stream THETA(1) 23.3545
    lvc <- log(209); label("Apparent central volume V/F at 70 kg (L)") # Table 2 'V' 209; stream THETA(2) 208.817
    lka <- log(0.344); label("First-order absorption rate constant ka (1/h)") # Table 2 'KA' 0.344; stream THETA(3) 0.344479
    lq <- log(59.2); label("Apparent intercompartmental clearance Q/F at 70 kg (L/h)") # Table 2 'Q' 59.2; stream THETA(4) 59.1631
    lvp <- log(912); label("Apparent peripheral volume VP/F at 70 kg (L)") # Table 2 'VP' 912; stream THETA(5) 912.068

    # Allometric exponents hard-coded in the control stream.
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent of weight on CL/F and Q/F (unitless)") # Methods 'allometric exponent fixed to 0.75'; stream LTVCL and MU_4
    e_wt_vc_vp <- fixed(1); label("Allometric exponent of weight on V/F and VP/F (unitless)") # stream LTVV and MU_5 '(LOG(WT)-LOG(70))'

    # CYP2B6*6 genotype fractions of CL/F in arm 1 (equation 3, Theta_12, Theta_13).
    e_cyp2b6_6het_cl <- 0.618; label("Fraction of CL/F for CYP2B6*1/*6 vs *1/*1 (unitless)") # Table 2 'Theta12' 0.618; stream THETA(8)
    e_cyp2b6_6hom_cl <- 0.311; label("Fraction of CL/F for CYP2B6*6/*6 vs *1/*1 (unitless)") # Table 2 'Theta13' 0.311; stream THETA(9)

    # HIV/TB arm during anti-TB co-treatment (week 20) vs arm 1, per genotype
    # (equation 3, theta_21_g, raised to A2).
    e_tb_cyp2b6_6wt_cl <- 0.915; label("Fraction of CL/F, HIV/TB arm vs arm 1, CYP2B6*1/*1 (unitless)") # Table 2 'theta211' 0.915; stream THETA(10)
    e_tb_cyp2b6_6het_cl <- 0.852; label("Fraction of CL/F, HIV/TB arm vs arm 1, CYP2B6*1/*6 (unitless)") # Table 2 'theta212' 0.852; stream THETA(11)
    e_tb_cyp2b6_6hom_cl <- 1.02; label("Fraction of CL/F, HIV/TB arm vs arm 1, CYP2B6*6/*6 (unitless)") # Table 2 'theta213' 1.02 (row text says *1/*6; it is *6/*6 per equation 3); stream THETA(12) 1.01737

    # HIV/TB arm week 32 (8 weeks after anti-TB treatment) vs week 20, per
    # genotype (equation 3, phi_22_g, raised to O2).
    e_occ2_cyp2b6_6wt_cl <- 1.19; label("Fraction of CL/F, HIV/TB arm week 32 vs week 20, CYP2B6*1/*1 (unitless)") # Table 2 'phi221' 1.19; stream THETA(13) 1.18644
    e_occ2_cyp2b6_6het_cl <- 1.19; label("Fraction of CL/F, HIV/TB arm week 32 vs week 20, CYP2B6*1/*6 (unitless)") # Table 2 'phi222' 1.19; stream THETA(14) 1.19174
    e_occ2_cyp2b6_6hom_cl <- 1.05; label("Fraction of CL/F, HIV/TB arm week 32 vs week 20, CYP2B6*6/*6 (unitless)") # Table 2 'phi223' 1.05 (row text says *1/*6; it is *6/*6 per equation 3); stream THETA(15) 1.04729

    e_male_cl <- 1.02; label("Fraction of CL/F for male vs female (unitless)") # Table 2 'thetaSEX' 1.02; stream THETA(18) 1.01959

    # Between-subject variability (variances). Table 2 omits the CL-V
    # covariance; the supplement control stream carries $OMEGA BLOCK(2) with
    # diagonal values equal to Table 2 and off-diagonal 0.300906.
    # VP/F IIV was fixed to 0 (Table 2 'VP' 0a; stream '$OMEGA BLOCK(1) FIX 0')
    # and is omitted.
    etalcl + etalvc ~ c(0.104, 0.300906, 0.91) # Table 2 BSV 'CL11' 0.104 and 'V' 0.91; covariance from stream $OMEGA BLOCK(2) 0.300906
    etalka ~ 0.412 # Table 2 BSV 'KA' 0.412; stream 0.411878
    etalq ~ 0.556 # Table 2 BSV 'Q' 0.556; stream 0.556338

    # Between-occasion variability on CL/F, HIV/TB arm only (stream GROUP 21
    # and 22; the second $OMEGA is SAME).
    etaiov_cl_1 ~ 0.0934925 # Table 2 'Between Occasion Variance' CL11 0.0934925; stream BOVCL21 0.093522
    etaiov_cl_2 ~ fixed(0.0934925) # stream BOVCL22 '$OMEGA BLOCK(1) SAME'

    # Combined residual error; $SIGMA 1 FIX so THETA(6) and THETA(7) are SDs.
    addSd <- 0.0873; label("Additive residual error SD (mg/L)") # Table 2 'ADD' 0.0873 (printed unit ug/L; the data are in ug/mL, Figure 1); stream THETA(6) 0.0872526
    propSd <- 0.242; label("Proportional residual error SD (fraction)") # Table 2 'PROP' 0.242; stream THETA(7) 0.241611
  })

  model({
    # CYP2B6*6 genotype indicators (B1, B2, B3 of equation 3).
    b1 <- (SNP_CYP2B6_RS3745274_T_COUNT == 0)
    b2 <- (SNP_CYP2B6_RS3745274_T_COUNT == 1)
    b3 <- (SNP_CYP2B6_RS3745274_T_COUNT == 2)

    # A2 = HIV/TB arm; O2 = HIV/TB arm at occasion 2 (week 32).
    a2 <- TB_POS
    o2 <- TB_POS * (OCC == 2)

    # Between-occasion eta on CL/F, applied only in the HIV/TB arm.
    iov_cl <- TB_POS * ((OCC == 1) * etaiov_cl_1 + (OCC == 2) * etaiov_cl_2)

    male <- 1 - SEXF

    # Equation 3 / control stream LTVCL.
    cl_geno <- e_cyp2b6_6het_cl^b2 * e_cyp2b6_6hom_cl^b3
    cl_arm <- (e_tb_cyp2b6_6wt_cl^b1 * e_tb_cyp2b6_6het_cl^b2 * e_tb_cyp2b6_6hom_cl^b3)^a2
    cl_occ <- (e_occ2_cyp2b6_6wt_cl^b1 * e_occ2_cyp2b6_6het_cl^b2 * e_occ2_cyp2b6_6hom_cl^b3)^o2

    cl <- exp(lcl + etalcl + iov_cl) * cl_geno * cl_arm * cl_occ * e_male_cl^male * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    ka <- exp(lka + etalka)
    q <- exp(lq + etalq) * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp) * (WT / 70)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose mg, volume L -> mg/L (= ug/mL, the units of Figure 1).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
