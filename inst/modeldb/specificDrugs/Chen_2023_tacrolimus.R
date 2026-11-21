Chen_2023_tacrolimus <- function() {
  description <- paste(
    "Population PK/PD model of oral tacrolimus in Chinese adults with",
    "myasthenia gravis (MG) (Chen 2023). PK: one-compartment model with",
    "first-order absorption and an absorption lag time, ka and tlag fixed to",
    "the values of the earlier Chen 2017 MG analysis because only trough",
    "whole-blood concentrations were collected; apparent clearance is 2.05-fold",
    "higher in CYP3A5 expressers (carriers of the functional star-1 allele) and",
    "apparent volume falls",
    "exponentially with serum total protein (centred on 65 g/L). PD: the",
    "Chinese quantitative MG total score (sum of nine item scores) falls from",
    "an Osserman-class-specific baseline through an Emax function of the",
    "cumulative tacrolimus AUC from the first dose, with Emax fixed at 1 and",
    "a cumulative AUC50 of 4110 ng h/mL. The relative total score",
    "(percentage improvement from baseline, responder threshold 25%) is",
    "provided as a derived output."
  )
  reference <- paste(
    "Chen D, Yao Q, Chen W, Yin J, Hou S, Tian X, Zhao M, Zhang H, Yang L,",
    "Zhou T, Jin P. Population PK/PD model of tacrolimus for exploring the",
    "relationship between accumulated exposure and quantitative scores in",
    "myasthenia gravis patients. CPT Pharmacometrics Syst Pharmacol.",
    "2023;12(7):963-976. doi:10.1002/psp4.12966"
  )
  vignette <- "Chen_2023_tacrolimus"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL",
    mgts = "quantitative MG score points",
    mgrts = "percent"
  )

  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE),
    auc_central = list(analyte = "tacrolimus", units = "ng*h/mL", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    CYP3A5_EXPR = list(
      description = "CYP3A5 expresser status: 1 = *1/*1 or *1/*3, 0 = *3/*3",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP3A5 *3/*3 nonexpresser)",
      notes = paste(
        "Source column CYP3A5 in the supplementary NONMEM code is coded the",
        "other way round (0 = *1/*1 and *1/*3, 1 = *3/*3), so",
        "CYP3A5_EXPR = 1 - CYP3A5. *1/*1 and *1/*3 were pooled because only 7",
        "patients were *1/*1. The 12 patients of unknown genotype were split",
        "by their base-model individual CL/F: the six with the higher CL/F",
        "were assigned to the expresser group and the other six to *3/*3",
        "(Methods 'Covariate effect in PK model'; Figure S1)."
      ),
      source_name = "CYP3A5"
    ),
    TPRO = list(
      description = "Serum total protein",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters V/F as exp(-0.0817 * (TP - 65)), centred on the cohort",
        "median of 65 g/L (Eq. 9; supplementary NONMEM $PK). Table 1 reports",
        "65.1 +/- 15.9 g/L (range 48-85). Lower total protein gives a larger",
        "V/F, which the authors attribute to the high plasma-protein binding",
        "of tacrolimus."
      ),
      source_name = "TP"
    ),
    DIS_MG_OSSERMAN2 = list(
      description = "Osserman class IIa or IIb myasthenia gravis indicator (1 = class IIa or IIb, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Osserman class I when DIS_MG_OSSERMAN2 and DIS_MG_OSSERMAN34 are both 0)",
      notes = paste(
        "Derived from the source OSM column (1 = I, 2 = IIa, 3 = IIb, 4 = III,",
        "5 = IV): DIS_MG_OSSERMAN2 = OSM in (2, 3). Selects the baseline score",
        "theta for IIa/IIb (10.4) in place of the class-I value (4.08)",
        "(supplementary $PRED 'OSSERMANB')."
      ),
      source_name = "OSM"
    ),
    DIS_MG_OSSERMAN34 = list(
      description = "Osserman class III or IV myasthenia gravis indicator (1 = class III or IV, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Osserman class I when DIS_MG_OSSERMAN2 and DIS_MG_OSSERMAN34 are both 0)",
      notes = paste(
        "Derived from the source OSM column: DIS_MG_OSSERMAN34 = OSM in (4, 5).",
        "Selects the baseline score theta for III/IV (19.0) in place of the",
        "class-I value (supplementary $PRED 'OSSERMANC'). No class-V patient",
        "was enrolled, so the model is not informed for class V."
      ),
      source_name = "OSM"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = paste(
        "Added in forward inclusion on the PK model (Table S1 model 7, dOFV",
        "-6.809) and on PD baseline TS0 (Table S2 model 2, dOFV -6.019) but",
        "removed in backward elimination at p < 0.001 in both. Uric acid and",
        "thymus status (normal / hyperplasia / thymoma; unknown assumed",
        "hyperplasia) were screened the same way and also not retained."
      )
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Added in PK forward inclusion (Table S1 model 4, dOFV -8.895) but removed in backward elimination."
    ),
    WBC = list(
      description = "White blood cell count",
      units = "10^9/L",
      type = "continuous",
      notes = "Added in PK forward inclusion (Table S1 model 6, dOFV -6.659) but removed in backward elimination."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 63L,
    n_studies = 1L,
    age_range = "21-87 years",
    age_mean = "52.4 years (SD 16.4)",
    weight_range = "40.5-97 kg",
    weight_mean = "69.2 kg (SD 13.4)",
    sex_female_pct = 46.0,
    race_ethnicity = "Chinese (Han)",
    disease_state = "Myasthenia gravis, Osserman class I (11.1%), IIa (20.6%), IIb (36.5%), III (6.3%), IV (25.4%); no class V",
    dose_range = "Oral tacrolimus starting at 1 mg/day given every 12 h, then titrated to a 5-10 ng/mL trough target with dosing intervals of 8, 12, 24 or 48 h",
    regions = "China (Beijing Hospital)",
    genotype = "CYP3A5 *3/*3 44.4%, *1/*3 25.4%, *1/*1 11.1%, unknown 19.0%",
    co_medication = "Pyridostigmine bromide in all patients; methylprednisolone 23.8%, amlodipine 14.3%, nifedipine 12.7%",
    n_observations = "107 trough whole-blood tacrolimus concentrations from 63 patients; 255 quantitative MG total scores from 55 of them",
    notes = "Retrospective therapeutic-drug-monitoring cohort, January 2013 to June 2018 (Chen 2023 Table 1). Trough samples drawn 1 h before the next dose and assayed by chemiluminescent microparticle immunoassay. The longest PD observation was on day 64."
  )

  ini({
    # PK (Table 2 final model; supplementary NONMEM 'PK model' $THETA)
    lcl <- log(15.7); label("Apparent clearance CL/F in CYP3A5 *3/*3 (L/h)") # Table 2 CL/F = 15.7 L/h (RSE 9.3%)
    lvc <- log(1410); label("Apparent volume V/F at total protein 65 g/L (L)") # Table 2 V/F = 1410 L (RSE 16.2%)
    lka <- fixed(log(0.502)); label("Absorption rate constant ka (1/h)") # Table 2 ka = 0.502 h-1 FIX (from Chen 2017)
    ltlag <- fixed(log(0.346)); label("Absorption lag time (h)") # Table 2 Tlag = 0.346 h FIX (from Chen 2017)
    e_cyp3a5_expr_cl <- 2.05; label("Fold change in CL/F for CYP3A5 expressers (unitless)") # Table 2 'CYP3A5*1/*1 and *1/*3 on CL/F' = 2.05 (RSE 14.1%)
    e_tpro_vc <- -0.0817; label("Exponential coefficient of total protein on V/F (per g/L)") # Table 2 'TP on V/F' = -0.0817 (RSE 28.9%)

    # PK IIV: $OMEGA variances; Table 2 IIV percent = sqrt(omega)
    etalcl ~ 0.233 # $OMEGA IIV CL = 0.233; Table 2 IIV CL/F 48.3%
    etalvc ~ 0.561 # $OMEGA IIV V = 0.561; Table 2 IIV V/F 74.9%

    propSd <- sqrt(0.0718); label("Proportional residual error, tacrolimus (fraction)") # $SIGMA 0.0718; Table 2 proportional error 26.8%

    # PD (Table 2 final PK/PD model; supplementary NONMEM 'PK/PD model' $THETA)
    lrbase <- log(4.08); label("Baseline quantitative MG total score, Osserman class I (score points)") # Table 2 TS0 for Osserman I = 4.08 (RSE 27.0%)
    lrbase_osserman2 <- log(10.4); label("Baseline quantitative MG total score, Osserman class IIa/IIb (score points)") # Table 2 TS0 for Osserman IIa/IIb = 10.4 (RSE 11.5%)
    lrbase_osserman34 <- log(19.0); label("Baseline quantitative MG total score, Osserman class III/IV (score points)") # Table 2 TS0 for Osserman III/IV = 19.0 (RSE 13.5%)
    lemax <- fixed(log(1)); label("Maximum fractional reduction of the total score (fraction)") # Table 2 Emax = 1 FIX
    lauc50 <- log(4110); label("Cumulative AUC giving half of Emax, EcAUC50 (ng*h/mL)") # Table 2 EcAUC50 = 4110 ng*h/mL (RSE 33.8%)

    # PD IIV: $OMEGA variances (IIV on Emax is 0 FIX in the code and is omitted)
    etalrbase ~ 0.355 # $OMEGA IIV TS0 = 0.355; Table 2 IIV TS0 59.6%
    etalauc50 ~ 2.75 # $OMEGA IIV EcAUC50 = 2.75; Table 2 IIV EcAUC50 165.8%

    propSd_mgts <- sqrt(0.0619); label("Proportional residual error, quantitative MG total score (fraction)") # $SIGMA 0.0619; Table 2 proportional error 24.9%
  })
  model({
    # Eq. 7-9: CL/F by CYP3A5 genotype, V/F exponential in total protein
    cl <- exp(lcl + etalcl) * e_cyp3a5_expr_cl^CYP3A5_EXPR
    vc <- exp(lvc + etalvc) * exp(e_tpro_vc * (TPRO - 65))
    ka <- exp(lka)
    tlag <- exp(ltlag)
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    lag(depot) <- tlag

    # mg/L -> ng/mL (NONMEM S2 = V/1000)
    Cc <- 1000 * central / vc

    # Cumulative AUC from the first dose (ng*h/mL), the paper's cAUC
    d/dt(auc_central) <- Cc

    # Eq. 10-12: Osserman-class-specific baseline score
    osm1 <- 1 - DIS_MG_OSSERMAN2 - DIS_MG_OSSERMAN34
    rbase <- (osm1 * exp(lrbase) + DIS_MG_OSSERMAN2 * exp(lrbase_osserman2) + DIS_MG_OSSERMAN34 * exp(lrbase_osserman34)) * exp(etalrbase)
    emax <- exp(lemax)
    auc50 <- exp(lauc50 + etalauc50)

    # Eq. 5: total score as an Emax function of cumulative AUC
    drug_eff <- emax * auc_central / (auc50 + auc_central)
    mgts <- rbase * (1 - drug_eff)

    # Eq. 6: relative total score (percent improvement from baseline)
    mgrts <- 100 * (rbase - mgts) / rbase

    Cc ~ prop(propSd)
    mgts ~ prop(propSd_mgts)
  })
}
