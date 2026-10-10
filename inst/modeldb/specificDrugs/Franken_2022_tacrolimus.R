Franken_2022_tacrolimus <- function() {
  description <- "Two-compartment population PK model with first-order absorption and an absorption lag time for twice-daily oral immediate-release tacrolimus (Prograft) in adult kidney transplant recipients, extended with an intracellular peripheral blood mononuclear cell (PBMC) effect compartment (Franken 2022 final model). The whole-blood layer is the Andrews 2019 model, carried unchanged: apparent oral clearance CL/F depends on CYP3A5 expresser status, CYP3A4*22 carriage, haematocrit, serum creatinine, serum albumin, age and body surface area, and apparent central volume V1/F on lean body mass. The PBMC compartment has no mass transfer from the central compartment: the intracellular concentration equilibrates with the whole-blood concentration at a fixed rate constant (0.9 1/h) towards a steady state 14.1-fold higher than whole blood, and this ratio rises with lean body mass (power exponent 1.01, centred at 59.5 kg) and falls with haematocrit (power exponent -1.22, centred at 34 percent), with 38.9 percent inter-individual variability. Residual error is proportional for whole blood and combined additive plus proportional for the PBMC concentration."
  reference <- "Franken LG, Francke MI, Andrews LM, van Schaik RHN, Li Y, de Wit LEA, Baan CC, Hesselink DA, de Winter BCM. A Population Pharmacokinetic Model of Whole-Blood and Intracellular Tacrolimus in Kidney Transplant Recipients. Eur J Drug Metab Pharmacokinet. 2022;47(4):523-535. doi:10.1007/s13318-022-00767-8. Whole-blood layer from: Andrews LM, Hesselink DA, van Schaik RHN, et al. A population pharmacokinetic model to predict the individual starting dose of tacrolimus in adult renal transplant recipients. Br J Clin Pharmacol. 2019;85(3):601-615. doi:10.1111/bcp.13838"
  vignette <- "Franken_2022_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    AGE = list(
      description = "Recipient age at transplantation",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Whole-blood layer only (Andrews 2019). Enters CL/F as (AGE/55.72)^-0.43. Franken 2022 Table 2 reports the exponent; the centring value 55.72 years is the Andrews 2019 Supporting Information Data S1 control-stream constant, carried with the rest of the fixed whole-blood model (see Andrews_2019_tacrolimus.R). Franken 2022 cohort median 57 years (IQR 46-64; Table 1).",
      source_name = "AGE"
    ),
    ALB = list(
      description = "Serum albumin concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Whole-blood layer only (Andrews 2019). Enters CL/F as (ALB/42)^0.43; centring value from Andrews 2019 Equation (1) and Data S1. Franken 2022 cohort median 43 g/L (IQR 40-46; Table 1). Albumin was also screened on the whole-blood:intracellular ratio and not retained (Table 3, delta OFV -0.001).",
      source_name = "ALB"
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Whole-blood layer only (Andrews 2019). Enters CL/F as (BSA/1.93)^0.88; centring value from Andrews 2019 Equation (1). Franken 2022 Methods computes BSA by the Mosteller formula, sqrt(height_cm * weight_kg / 3600); cohort median 1.97 m^2 (IQR 1.80-2.14; Table 1). BSA was also screened on the whole-blood:intracellular ratio and dropped at backward elimination (Table 3).",
      source_name = "BSA"
    ),
    CREAT = list(
      description = "Serum creatinine concentration",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Whole-blood layer only (Andrews 2019). Enters CL/F as (CREAT/134.98)^-0.14; centring value from the Andrews 2019 Data S1 control stream. Franken 2022 cohort median 135 umol/L (IQR 109-171; Table 1).",
      source_name = "CREA"
    ),
    HCT = list(
      description = "Haematocrit, packed red blood cell volume fraction",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Used twice. (1) Whole-blood layer (Andrews 2019): CL/F scales as (HCT/34)^-0.76. (2) Intracellular layer (Franken 2022 Equation 4 and Supplementary Data S1 'RPIC = THETA(5)*(LBW/59.5)**THETA(6)*(HCT/0.34)**THETA(7)'): the PBMC:whole-blood ratio scales as (HCT/34)^-1.22. Franken 2022 records haematocrit as a fraction in L/L (cohort median 0.34, IQR 0.31-0.38; Table 1), but the canonical HCT column is on the percent scale, so both centring values are multiplied by 100 (0.34 L/L -> 34 percent), following Andrews_2019_tacrolimus.R. Supplying the column in L/L would misscale CL/F by 100^-0.76 and the PBMC ratio by 100^-1.22. Higher haematocrit lowers the intracellular share of tacrolimus, which the authors attribute to binding to erythrocytes (Discussion).",
      source_name = "HCT"
    ),
    LBM = list(
      description = "Lean body mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Reported as 'lean body weight (LBW)', a registered synonym of LBM. Franken 2022 Methods computes it by the James formula: female 1.07 * weight_kg - 148 * (weight_kg / height_cm)^2, male 1.1 * weight_kg - 128 * (weight_kg / height_cm)^2; cohort median 60.9 kg (IQR 53.8-66.9; Table 1). Used twice: on V1/F in the whole-blood layer as (LBM/58.94)^1.52 (Andrews 2019 Data S1 centring value), and on the PBMC:whole-blood ratio as (LBM/59.5)^1.01 (Franken 2022 Equation 4 and Supplementary Data S1). The two centring values differ because they are the medians of the two model-building cohorts.",
      source_name = "LBW"
    ),
    CYP3A5_EXPR = list(
      description = "CYP3A5 expresser indicator: 1 if the patient carries at least one functional CYP3A5*1 allele (genotype *1/*1 or *1/*3), 0 if a non-expresser.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP3A5 non-expresser, *3/*3)",
      notes = "Whole-blood layer only (Andrews 2019); multiplies CL/F by 1.63 (Franken 2022 Table 2). Franken 2022 cohort: 41 (22.3%) expressers, 143 (77.7%) non-expressers (Table 1).",
      source_name = "NEXP"
    ),
    SNP_CYP3A4_RS35599367 = list(
      description = "CYP3A4*22 carrier indicator (rs35599367): 1 if the patient carries at least one *22 allele, 0 if CYP3A4*1 homozygous wild-type or genotype unknown.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP3A4*1 wild-type, or genotype unknown -- Andrews 2019 pools these two groups)",
      notes = "Whole-blood layer only (Andrews 2019); multiplies CL/F by 0.80 (Franken 2022 Table 2). Franken 2022 cohort: 20 (10.9%) carriers, 155 (84.2%) non-carriers, 9 (4.9%) unknown (Table 1).",
      source_name = "CYP"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened on the whole-blood:intracellular ratio; significant in the univariate step (delta OFV -6.362, exponent 0.514) but dropped at backward elimination (Franken 2022 Table 3). Enters this model only indirectly, through LBM and BSA. Cohort median 80.0 kg (IQR 69.2-92.0; Table 1)."
    ),
    IBW = list(
      description = "Ideal body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened on the whole-blood:intracellular ratio; significant in the univariate step (delta OFV -8.405, exponent 1.07) but dropped at backward elimination (Franken 2022 Table 3). Computed by the Robinson formula (Methods)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Listed among the screened body-size covariates (Franken 2022 Methods 2.3.2) but absent from Table 3; not retained. Cohort median 25.9 kg/m^2 (IQR 23.7-29.5; Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Screened on the whole-blood:intracellular ratio; not significant (delta OFV -0.537; Franken 2022 Table 3)."
    ),
    SNP_ABCB1_RS2229109 = list(
      description = "ABCB1 1199G>A genotype (rs2229109)",
      units = "(binary)",
      type = "binary",
      notes = "Screened on the whole-blood:intracellular ratio by genotype (AA, GA, GG); no level significant (Franken 2022 Table 3)."
    ),
    SNP_ABCB1_RS1045642 = list(
      description = "ABCB1 3435C>T genotype (rs1045642)",
      units = "(binary)",
      type = "binary",
      notes = "Screened on the whole-blood:intracellular ratio by genotype (CC, CT, TT); no level significant (Franken 2022 Table 3)."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE),
    # The `effect` state holds a CONCENTRATION (ng/mL), not an amount:
    # Supplementary Data S1 integrates DADT(4) = K24*(RPIC*(A(2)/V2) - A(4))
    # and reports A(4) directly as the intracellular concentration (no S4).
    effect = list(analyte = "tacrolimus", units = "ng/mL", specimen = "blood cell", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 184L,
    n_studies = 1L,
    n_observations = 590L,
    age_median = "57 years (IQR 46-64)",
    weight_median = "80.0 kg (IQR 69.2-92.0)",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported",
    disease_state = "Adult kidney transplant recipients (living related 39.1 percent, living unrelated 60.9 percent donors) in the first 3 months after transplantation, treated with oral twice-daily immediate-release tacrolimus (Prograft) with mycophenolic acid, prednisolone and basiliximab induction.",
    dose_range = "Oral twice-daily tacrolimus titrated by therapeutic drug monitoring to whole-blood pre-dose targets of 10.0-15.0 ng/mL in weeks 1-2, 8.0-12.0 ng/mL in weeks 3-4 and 5.0-10.0 ng/mL after week 4.",
    regions = "The Netherlands (Erasmus MC, Rotterdam)",
    sampling_window = "406 whole-blood concentrations (median 2 per patient, range 1-10) and 184 intracellular PBMC concentrations (exactly one per patient, a pre-dose sample at month 3 post-transplantation).",
    assay = "Whole blood: antibody-conjugated magnetic immunoassay (ACMIA) and enzyme-multiplied immunoassay technique (EMIT), r = 0.97 between them. PBMC: Ficoll isolation, then LC-MS/MS; intracellular concentration computed from the PBMC cell count and mean cell volume. Samples with a visual red level above 2 (erythrocyte contamination) were excluded.",
    notes = "Post-hoc analysis of the randomized controlled trial of CYP3A5 genotype-based versus body-weight-based tacrolimus dosing (Shuker 2016, NTR2226); 184 of its 237 patients had a usable month-3 intracellular sample (Supplementary Figure S2). Sex and race distributions are not reported (Table 1). Baseline characteristics in Table 1; median haematocrit 0.34 L/L (IQR 0.31-0.38), albumin 43 g/L, creatinine 135 umol/L, lean body weight 60.9 kg."
  )

  ini({
    # Whole-blood layer. Franken 2022 fixed every whole-blood parameter to the
    # Andrews 2019 model (Table 2 footnote a: 'Whole-blood model parameter
    # values were fixed at individual values using the previously reported
    # model of Andrews et al.'), so each is fixed() here at the population
    # value printed in Franken 2022 Table 2 'Final model'. Bioavailability is
    # fixed at 1, so all clearances and volumes are apparent (CL/F, V/F).
    # Supplementary Data S1 hardcodes KA = 3.6 and ALAG1 = 0.382, the same
    # values at different rounding.
    ltlag <- fixed(log(0.38)); label("Absorption lag time tlag (h)") # Table 2 'T lag (h)' = 0.38, fixed; Data S1 ALAG1 = 0.382
    lka <- fixed(log(3.58)); label("Absorption rate constant ka (1/h)") # Table 2 'k a (h-1)' = 3.58, fixed; Data S1 KA = 3.6
    lcl <- fixed(log(23.0)); label("Apparent oral clearance CL/F at the reference subject (L/h)") # Table 2 'CL/F (l h-1)' = 23.0, fixed
    lvc <- fixed(log(692)); label("Apparent central volume V1/F at the reference subject (L)") # Table 2 'V 1 /F (l)' = 692, fixed
    lq <- fixed(log(11.6)); label("Apparent inter-compartmental clearance Q/F (L/h)") # Table 2 'Q 1 /F (l h-1)' = 11.6, fixed
    lvp <- fixed(log(5340)); label("Apparent peripheral volume V2/F (L)") # Table 2 'V 2 /F (l)' = 5340, fixed

    # Whole-blood covariate effects, Table 2 'Covariate effect on CL' and
    # 'Covariate effect on V 1', all fixed. The centring values in model()
    # are the Andrews 2019 Data S1 control-stream constants (not printed in
    # Franken 2022); see Andrews_2019_tacrolimus.R for their provenance.
    e_cyp3a5_expr_cl <- fixed(1.63); label("CYP3A5 expresser (*1/*1 or *1/*3) multiplier on CL/F") # Table 2 'CYP3A5*1' = 1.63, fixed
    e_cyp3a4_22_cl <- fixed(0.80); label("CYP3A4*22 carrier multiplier on CL/F") # Table 2 'CYP3A4*22' = 0.80, fixed
    e_hct_cl <- fixed(-0.76); label("Haematocrit power exponent on CL/F, centred at 34 percent (unitless)") # Table 2 'Haematocrit (l l-1)' = -0.76, fixed
    e_creat_cl <- fixed(-0.14); label("Serum creatinine power exponent on CL/F, centred at 134.98 umol/L (unitless)") # Table 2 'Creatinine' = -0.14, fixed
    e_alb_cl <- fixed(0.43); label("Serum albumin power exponent on CL/F, centred at 42 g/L (unitless)") # Table 2 'Albumin (g l-1)' = 0.43, fixed
    e_age_cl <- fixed(-0.43); label("Age power exponent on CL/F, centred at 55.72 years (unitless)") # Table 2 'Age (years)' = -0.43, fixed
    e_bsa_cl <- fixed(0.88); label("Body surface area power exponent on CL/F, centred at 1.93 m^2 (unitless)") # Table 2 'BSA (m2)' = 0.88, fixed
    e_lbm_vc <- fixed(1.52); label("Lean body mass power exponent on V1/F, centred at 58.94 kg (unitless)") # Table 2 'Covariate effect on V 1, Lean body weight' = 1.52, fixed

    # Intracellular (PBMC) layer, Franken 2022 Equations 3-4 and
    # Supplementary Data S1:
    #   DADT(4) = K24*(RPIC*(A(2)/V2) - A(4))
    #   RPIC = THETA(5)*(LBW/59.5)**THETA(6)*(HCT/0.34)**THETA(7)*EXP(ETA(1))
    # K24 (the paper's K WB-IC) is THETA(4) = 0.9 FIX.
    lke0 <- fixed(log(0.9)); label("Whole-blood to PBMC equilibration rate constant K WB-IC (1/h)") # Table 2 'K WB-IC' = 0.9, fixed; Results 3.2 'K WB-IC was fixed at 0.9'; Data S1 THETA(4) (0.9) FIX
    # The control stream multiplies RPIC by A(2)/V2, which is in mg/L (dose
    # mg, volume L), while A(4) is the intracellular concentration in ug/L
    # (Supplementary Table S1 'IC (ug/L)'). The estimate 14100 therefore
    # carries a mg/L -> ug/L factor of 1000. Here the effect compartment is
    # driven by Cc, already in ng/mL (= ug/L), so the dimensionless
    # PBMC:whole-blood ratio is 14100 / 1000 = 14.1 ('on average, there was a
    # 14-fold higher concentration in the PBMCs', Results 3.2).
    lppc <- log(14100 / 1000); label("PBMC to whole-blood steady-state concentration ratio at the reference subject (unitless)") # Table 2 'R WB:IC' final model = 14100 (bootstrap 14023, 12595-15633); Equation 4; divided by 1000 for the mg/L -> ug/L scale
    e_lbm_ppc <- 1.01; label("Lean body mass power exponent on the PBMC:whole-blood ratio, centred at 59.5 kg (unitless)") # Table 2 'Covariate effect on R WB:IC, Lean body weight' = 1.01 (bootstrap 1.002, 0.56-1.54); Equation 4
    e_hct_ppc <- -1.22; label("Haematocrit power exponent on the PBMC:whole-blood ratio, centred at 34 percent (unitless)") # Table 2 'Haematocrit' = -1.22 (bootstrap -1.21, -1.96 to -0.42); Equation 4

    # Inter-individual variability. Table 2 reports IIV as %CV; variances on
    # the log scale are omega^2 = log(1 + CV^2), the convention used for the
    # same Andrews 2019 values in Andrews_2019_tacrolimus.R:
    #   CL/F 38.6% -> 0.138889; V1/F 49.2% -> 0.216775; V2/F 53.0% -> 0.247563;
    #   Q/F 78.7% -> 0.482037; ratio 38.9% -> log(1 + 0.389^2) = 0.140910.
    # The four whole-blood variances are fixed (Table 2 footnote a).
    etalcl ~ fixed(0.138889) # Table 2 IIV CL 38.6%, held at the Andrews 2019 value (footnote a)
    etalvc ~ fixed(0.216775) # Table 2 IIV V1 49.2%, held at the Andrews 2019 value (footnote a)
    etalvp ~ fixed(0.247563) # Table 2 IIV V2 53.0%, held at the Andrews 2019 value (footnote a)
    etalq ~ fixed(0.482037) # Table 2 IIV Q 78.7%, held at the Andrews 2019 value (footnote a)
    etalppc ~ 0.140910 # Table 2 IIV R WB-IC final model 38.9% [shrinkage 33%]; Data S1 $OMEGA ETA(1) on RPIC

    # Residual error. Data S1 $ERROR: CMT 2 (whole blood) Y = F + F*EPS(1)*THETA(1);
    # CMT 4 (PBMC) Y = F + F*EPS(2)*THETA(2) + EPS(3)*THETA(3), with every
    # $SIGMA 1 FIX, so each THETA is a standard deviation. The whole-blood
    # proportional error was re-estimated in Franken 2022 (it carries a
    # bootstrap interval and no footnote a), so it differs from Andrews 2019.
    propSd <- 0.611; label("Proportional residual SD for whole-blood tacrolimus (fraction)") # Table 2 'Proportional WB' final model = 0.611 (bootstrap 0.584, 0.292-0.874); Data S1 THETA(1)
    propSd_Cpbmc <- 0.184; label("Proportional residual SD for intracellular PBMC tacrolimus (fraction)") # Table 2 'Proportional IC' final model = 0.184 (bootstrap 0.179, 0.139-0.210); Data S1 THETA(2)
    addSd_Cpbmc <- 36.8; label("Additive residual SD for intracellular PBMC tacrolimus (ng/mL)") # Table 2 'Additive IC' final model = 36.8 (bootstrap 36.1, 20.5-46.4), in ug/L = ng/mL; Data S1 THETA(3)
  })

  model({
    # Whole-blood layer (Andrews 2019, fixed). CYP3A5 and CYP3A4*22 enter as
    # multipliers whose reference level is exactly 1.
    f_cyp3a5 <- 1 + (e_cyp3a5_expr_cl - 1) * CYP3A5_EXPR
    f_cyp3a4 <- 1 + (e_cyp3a4_22_cl - 1) * SNP_CYP3A4_RS35599367
    f_age <- (AGE / 55.72)^e_age_cl
    f_alb <- (ALB / 42)^e_alb_cl
    f_bsa <- (BSA / 1.93)^e_bsa_cl
    f_creat <- (CREAT / 134.98)^e_creat_cl
    f_hct <- (HCT / 34)^e_hct_cl
    f_lbm <- (LBM / 58.94)^e_lbm_vc

    tlag <- exp(ltlag)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * f_cyp3a5 * f_cyp3a4 * f_age * f_alb * f_bsa * f_creat * f_hct
    vc <- exp(lvc + etalvc) * f_lbm
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    # Intracellular layer, Franken 2022 Equation 4 (power covariate model,
    # Equation 1, centred at the cohort medians LBW 59.5 kg and haematocrit
    # 0.34 L/L = 34 percent).
    ke0 <- exp(lke0)
    ppc <- exp(lppc + etalppc) * (LBM / 59.5)^e_lbm_ppc * (HCT / 34)^e_hct_ppc

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # Whole-blood concentration: dose mg / volume L = mg/L, times 1000 = ng/mL
    # (Data S1 'S2 = V2 / 1000').
    Cc <- central / vc * 1000

    # PBMC effect compartment without mass transfer (Franken 2022 Equation 3,
    # Figure 1): the intracellular concentration relaxes towards ppc * Cc at
    # rate ke0, and removes no drug from the central compartment. ppc sits
    # inside the derivative, as in Data S1, so a time-varying haematocrit or
    # lean body mass acts through the same lag.
    d/dt(effect) <- ke0 * (ppc * Cc - effect)
    Cpbmc <- effect

    Cc ~ prop(propSd)
    Cpbmc ~ add(addSd_Cpbmc) + prop(propSd_Cpbmc)
  })
}
