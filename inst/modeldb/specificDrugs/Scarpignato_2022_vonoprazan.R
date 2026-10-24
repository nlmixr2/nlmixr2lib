Scarpignato_2022_vonoprazan <- function() {
  description <- paste0(
    "Two-compartment population PK model for the potassium-competitive ",
    "acid blocker vonoprazan in adult healthy volunteers and patients with ",
    "gastroesophageal reflux disease (GERD) from Japan and Europe ",
    "(Scarpignato 2022: 1156 subjects, 15 studies). First-order absorption ",
    "after a lag time into the central compartment, with a dose-dependent ",
    "(over-proportional) relative bioavailability centred at 20 mg. The lag ",
    "time carries a two-class mixture (a delayed-absorption subpopulation ",
    "with a 107 percent longer lag; mixing fraction not reported). ",
    "Absorption rate carries high-fat-meal, non-Asian, female and GERD ",
    "patient effects; clearance carries body-weight and serum-creatinine ",
    "power effects and a patient effect; central volume carries ",
    "serum-creatinine and age power effects and a patient effect; ",
    "peripheral volume carries a non-Asian effect. Between-subject ",
    "variability on absorption rate, on clearance and central volume ",
    "(correlated) and on peripheral volume; combined proportional plus a ",
    "fixed, negligible additive residual error."
  )
  reference <- paste(
    "Scarpignato C, Leifke E, Smith N, Mulford DJ, Lahu G, Facius A,",
    "Howden CW. A Population Pharmacokinetic Model of Vonoprazan:",
    "Evaluating the Effects of Race, Disease Status, and Other Covariates",
    "on Exposure. J Clin Pharmacol. 2022;62(6):801-811.",
    "doi:10.1002/jcph.2019.",
    sep = " "
  )
  vignette <- "Scarpignato_2022_vonoprazan"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Doses are in mg and volumes in L, so central / vc is mg/L (= ug/mL);
  # the observation multiplies by 1000 to report ng/mL, the unit of the
  # Table 3 additive error and of Figures 2-4.
  compartmentData <- list(
    depot = list(analyte = "vonoprazan", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "vonoprazan", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vonoprazan", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    DOSE_VONOPRAZAN_MG = list(
      description = "Administered vonoprazan dose level carried on every record of the dosing interval, in mg.",
      units = "mg",
      type = "continuous",
      reference_category = "20 mg (the dose at which relative bioavailability equals 1)",
      notes = paste(
        "Power model on relative bioavailability, Frel = (DOSE/20)^0.389",
        "(Table 3 'Frel, Dose effect, exp'). The paper does not print the",
        "dose centring value, but 20 mg is the reference-subject dose of",
        "Figures 3-4, and the Figure 4 reference AUC of about 161 ng*h/mL",
        "equals 20 mg / 119 L/h only with Frel(20 mg) = 1. The positive",
        "exponent means higher doses have a higher relative bioavailability:",
        "doubling or halving the dose from 20 mg changes exposure by",
        "2^1.389, which reproduces the Figure 4 dose rows (+161.5 / -61.8",
        "percent). Studied doses span 1 to 120 mg (Table 1)."
      ),
      source_name = "DOSE"
    ),
    MIX_LONG_TLAG = list(
      description = "Latent mixture-model class indicator: 1 = delayed-absorption subpopulation with a longer lag time, 0 = reference lag time.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (reference lag time; the class the paper's reference-subject simulations represent)",
      notes = paste(
        "Methods: 'A mixture-model was used to describe the apparent bimodal",
        "distribution for the BSV on absorption rate (Ka), which was not",
        "related to any tested patient or treatment characteristic.' Table 3",
        "carries the class difference as the tlag row 'Delay, % = 107'",
        "(RSE 2.8 percent, 95 percent CI 101-113), applied as",
        "tlag * (1 + 107/100) for the respective subpopulation (Table 3",
        "footnote), so the delayed class has a lag of 0.225 * 2.07 =",
        "0.466 h. The mixing probability is NOT reported anywhere in the",
        "paper or its supplement, so no Bernoulli draw can be reproduced.",
        "Class 0 is treated as the reference: Figure 3 typical steady-state",
        "profiles rise from about 0.2 h, not 0.47 h, and the Figure 4",
        "reference tmax (about 1.41 h) matches 1.38 h simulated with",
        "tlag = 0.225 h against 1.62 h with the delay. Set",
        "MIX_LONG_TLAG = 0 for typical-value simulation."
      ),
      source_name = "mixture subpopulation (NONMEM $MIXTURE class; not named in the paper)"
    ),
    FED_HIGHFAT = list(
      description = "High-fat meal at dosing indicator, 1 = dosed after a high-fat meal, 0 = fasted or non-high-fat meal.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted or non-high-fat meal)",
      notes = paste(
        "Table 3 'Ka, Food effect, % = -58.8'. The reference subject of",
        "Figures 3-4 is dosed 'under fasted/non-high-fat meal conditions',",
        "so the food effect contrasts a high-fat meal against fasted and",
        "non-high-fat-meal dosing pooled. Methods lists the screened",
        "covariate as 'food (fed/fasted)'; the meal composition of the",
        "food-effect study is not given in the paper."
      ),
      source_name = "food"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator, 1 = Asian, 0 = non-Asian.",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (Asian; the paper's reference subject)",
      notes = paste(
        "Table 3 reports 'Non-Asian effect' on Ka (-20.2 percent) and on",
        "V3 (-11.7 percent), with the Asian reference subject of Figures",
        "3-4. The model applies each effect as",
        "(1 + effect * (1 - RACE_ASIAN)), keeping the coefficients verbatim",
        "and the Asian typical values unchanged. Non-Asian subjects in the",
        "dataset were White (32.1 percent), Black (1.1 percent) and Other",
        "(0.3 percent) (Table 2); Asian subjects were from the Japanese",
        "studies."
      ),
      source_name = "race (Asian vs non-Asian)"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Table 3 'Ka, Female effect, % = -17.2'; reference subject is male (Figures 3-4).",
      source_name = "sex"
    ),
    DIS_GERD = list(
      description = "GERD patient indicator, 1 = patient with gastroesophageal reflux disease, 0 = healthy volunteer.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy volunteer, including the renal-impairment volunteers of study 113)",
      notes = paste(
        "Table 3 'Patient effect' on Ka (-61.5 percent), CL (-33.8",
        "percent) and V2 (-32.9 percent). Methods defines the screened",
        "covariate as 'patient (healthy volunteer or patient with",
        "gastroesophageal reflux disease)'. The patient studies are",
        "CCT-001 (Japanese erosive esophagitis) and Vonoprazan-2001",
        "(European non-erosive reflux disease or mild erosive esophagitis)",
        "(Table 1). The 24 renal-impairment participants of study 113 were",
        "not GERD patients, so they carry DIS_GERD = 0."
      ),
      source_name = "patient"
    ),
    WT = list(
      description = "Baseline body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power model on CL, (WT/70)^0.409, centred at 70 kg (Table 3 and its footnote). Cohort mean 69 kg (SD 12.9; Table 2).",
      source_name = "baseline body weight"
    ),
    CREAT = list(
      description = "Baseline serum creatinine.",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power models on CL, (CREAT/72)^-0.36, and on V2,",
        "(CREAT/72)^-0.294, centred at 72 umol/L (Table 3 and its",
        "footnote). Cohort mean 73 umol/L (SD 21.6; Table 2). Participants",
        "with GFR < 30 mL/min were excluded from the pooled studies",
        "(Discussion), so very high creatinine values are extrapolation."
      ),
      source_name = "baseline serum creatinine"
    ),
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power model on V2, (AGE/48)^0.466, centred at 48 years (Table 3 and its footnote). Cohort mean 48 years (SD 17.7; Table 2).",
      source_name = "age"
    )
  )

  covariatesDataExcluded <- list(
    CYP2C19_PM = list(
      description = "CYP2C19 poor-metabolizer indicator, 1 = PM, 0 = other.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = paste(
        "CYP2C19 status (EM, IM, PM, other, not screened) was screened in",
        "the stepwise covariate search (Methods 'Covariate Screening') but",
        "not retained in the final model (Table 3). Distribution in Table",
        "2: EM 56.2, IM 3.2, PM 9.4, other 7.3, not screened 23.9 percent.",
        "The Discussion attributes the lack of effect to metabolism through",
        "several CYP enzymes (mainly CYP3A4/5)."
      ),
      source_name = "CYP2C19 status"
    ),
    CYP2C19_IM = list(
      description = "CYP2C19 intermediate-metabolizer indicator, 1 = IM, 0 = other.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Screened with CYP2C19 status but not retained; see CYP2C19_PM.",
      source_name = "CYP2C19 status"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1156L,
    n_studies = 15L,
    n_observations = 15168L,
    age_mean = "48 years (SD 17.7)",
    weight_mean = "69 kg (SD 12.9)",
    creatinine_mean = "73 umol/L (SD 21.6)",
    sex_female_pct = 26.0,
    race_ethnicity = "Asian 66.4 percent, White 32.1 percent, Black 1.1 percent, Other 0.3 percent",
    disease_state = "healthy volunteers (387), patients with GERD (745: Japanese erosive esophagitis and European non-erosive reflux disease or mild erosive esophagitis), and renal-impairment volunteers (24)",
    dose_range = "1-120 mg oral film-coated tablet; single doses and once- or twice-daily repeated dosing up to 8 weeks",
    regions = "Japan, Europe",
    notes = paste(
      "Demographics from Scarpignato 2022 Table 2; study designs, doses and",
      "PK sampling from Table 1 (13 phase 1 studies in healthy volunteers,",
      "including renal-impairment volunteers in study 113, and 2 phase 2",
      "studies in GERD patients). CYP2C19 status: EM 56.2, IM 3.2, PM 9.4,",
      "other 7.3, not screened 23.9 percent. Four film-coated tablet",
      "formulations with minor excipient differences were pooled. BSV on",
      "Ka was not estimated for study CCT-001, which had only trough",
      "samples. Fitted with NONMEM 7.4.1, FOCE with interaction."
    )
  )

  ini({
    # All values from Scarpignato 2022 Table 3 'Final Model Parameter
    # Estimates'. Categorical effects are printed in percent and applied as
    # tv * (1 + eff/100) (Table 3 footnote); they are stored here as
    # fractions. Continuous effects are power exponents on (cov/ref) with
    # references 48 years, 70 kg and 72 umol/L (Table 3 footnote).

    # Relative bioavailability is 1 at the 20 mg reference dose (no
    # intravenous data; the paper estimates only the dose exponent).
    lfdepot <- fixed(log(1)); label("Relative oral bioavailability at the 20 mg reference dose (unitless)")
    e_dose_fdepot <- 0.389; label("Power exponent on DOSE_VONOPRAZAN_MG/20 for relative bioavailability (unitless)") # Table 3 'Frel, Dose effect, exp' = 0.389 (RSE 4.1)

    ltlag <- log(0.225); label("Typical absorption lag time (h)") # Table 3 'tlag, TV, h' = 0.225 (RSE 1.4)
    e_mix_long_tlag_tlag <- 1.07; label("Fractional increase in lag time in the delayed-absorption mixture class (unitless)") # Table 3 'tlag, Delay, %' = 107 (RSE 2.8)

    lka <- log(2.51); label("Typical first-order absorption rate constant (1/h)") # Table 3 'Ka, TV' = 2.51 (RSE 12.3; printed unit 'L/h' is a typo for 1/h)
    e_fed_highfat_ka <- -0.588; label("Fractional change in Ka after a high-fat meal (unitless)") # Table 3 'Ka, Food effect, %' = -58.8 (RSE 10.3)
    e_race_nonasian_ka <- -0.202; label("Fractional change in Ka in non-Asian subjects (unitless)") # Table 3 'Ka, Non-Asian effect, %' = -20.2 (RSE 49.6)
    e_sexf_ka <- -0.172; label("Fractional change in Ka in females (unitless)") # Table 3 'Ka, Female effect, %' = -17.2 (RSE 57.8)
    e_dis_gerd_ka <- -0.615; label("Fractional change in Ka in GERD patients (unitless)") # Table 3 'Ka, Patient effect, %' = -61.5 (RSE 14.3)

    lcl <- log(119); label("Typical apparent clearance (L/h)") # Table 3 'CL, TV, L/h' = 119 (RSE 19.9)
    e_wt_cl <- 0.409; label("Power exponent on WT/70 for clearance (unitless)") # Table 3 'CL, Weight effect, exp' = 0.409 (RSE 19.9)
    e_creat_cl <- -0.36; label("Power exponent on CREAT/72 for clearance (unitless)") # Table 3 'CL, Creatinine effect, exp' = -0.36 (RSE 24.9)
    e_dis_gerd_cl <- -0.338; label("Fractional change in clearance in GERD patients (unitless)") # Table 3 'CL, Patient effect, %' = -33.8 (RSE 15.6)

    lvc <- log(1037); label("Typical apparent central volume of distribution (L)") # Table 3 'V2, TV, L' = 1037 (RSE 14.8)
    e_creat_vc <- -0.294; label("Power exponent on CREAT/72 for central volume (unitless)") # Table 3 'V2, Creatinine effect, exp' = -0.294 (RSE 27.6)
    e_age_vc <- 0.466; label("Power exponent on AGE/48 for central volume (unitless)") # Table 3 'V2, Age effect, exp' = 0.466 (RSE 7.4)
    e_dis_gerd_vc <- -0.329; label("Fractional change in central volume in GERD patients (unitless)") # Table 3 'V2, Patient effect (%)' = -32.9 (RSE 21.8)

    lq <- log(50.4); label("Typical apparent intercompartmental clearance (L/h)") # Table 3 'Q, TV, L/h' = 50.4 (RSE 7.2)

    lvp <- log(324); label("Typical apparent peripheral volume of distribution (L)") # Table 3 'V3, TV, L' = 324 (RSE 4.7)
    e_race_nonasian_vp <- -0.117; label("Fractional change in peripheral volume in non-Asian subjects (unitless)") # Table 3 'V3, Non-Asian effect, %' = -11.7 (RSE 25.9)

    # Between-subject variability: Table 3 footnote a gives the BSV rows on
    # the standard-deviation scale (exponential model), so each variance is
    # the square of the printed SD; the 95 percent CIs (e.g. Ka 0.8, 1.0)
    # are printed on that same SD scale. The CL-V2 covariance is
    # 0.915 * 0.41 * 0.361 = 0.135426.
    etalka ~ 0.876096 # Table 3 'Ka, BSV, %' = 93.6 as SD (0.936^2)
    etalcl + etalvc ~ c(0.1681,
                        0.135426, 0.130321) # Table 3 'CL, BSV, %' = 41 SD; 'BSV (CL,V2)' correlation = 0.915; 'V2, BSV (%)' = 36.1 SD
    etalvp ~ 0.037636 # Table 3 'V3, BSV, %' = 19.4 as SD (0.194^2)

    # Residual error: combined additive and proportional error model
    # (Results 'Final Model'); the additive term is fixed (Table 3 footnote b).
    propSd <- 0.287; label("Proportional residual error (fraction)") # Table 3 'Residual error, prop.err., %' = 28.7 (RSE 1.7)
    addSd <- fixed(0.00625); label("Additive residual error (ng/mL)") # Table 3 'Residual error, add.err., ng/mL' = 0.00625, fixed
  })

  model({
    # 1. Power-model covariate terms (Table 3 footnote references).
    dose_fdepot <- (DOSE_VONOPRAZAN_MG / 20)^e_dose_fdepot
    wt_cl <- (WT / 70)^e_wt_cl
    creat_cl <- (CREAT / 72)^e_creat_cl
    creat_vc <- (CREAT / 72)^e_creat_vc
    age_vc <- (AGE / 48)^e_age_vc

    # 2. Individual parameters. Categorical effects apply as (1 + eff) for
    #    the indicated subpopulation; the non-Asian effects use the
    #    complement of RACE_ASIAN because the reference subject is Asian.
    frel <- exp(lfdepot) * dose_fdepot
    tlag <- exp(ltlag) * (1 + e_mix_long_tlag_tlag * MIX_LONG_TLAG)
    ka <- exp(lka + etalka) * (1 + e_fed_highfat_ka * FED_HIGHFAT) *
      (1 + e_race_nonasian_ka * (1 - RACE_ASIAN)) * (1 + e_sexf_ka * SEXF) *
      (1 + e_dis_gerd_ka * DIS_GERD)
    cl <- exp(lcl + etalcl) * wt_cl * creat_cl * (1 + e_dis_gerd_cl * DIS_GERD)
    vc <- exp(lvc + etalvc) * creat_vc * age_vc * (1 + e_dis_gerd_vc * DIS_GERD)
    q <- exp(lq)
    vp <- exp(lvp + etalvp) * (1 + e_race_nonasian_vp * (1 - RACE_ASIAN))

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system: oral depot -> central (V2) <-> peripheral (V3)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Absorption lag and dose-dependent relative bioavailability on the
    #    oral dose.
    alag(depot) <- tlag
    f(depot) <- frel

    # 6. Observation: mg / L * 1000 = ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
