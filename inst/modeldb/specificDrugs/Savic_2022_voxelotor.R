Savic_2022_voxelotor <- function() {
  description <- paste(
    "Joint plasma and whole-blood population PK model for oral voxelotor in",
    "adults and adolescents (12-59 years) with sickle cell disease (Savic",
    "2022): two-compartment model with first-order absorption and",
    "elimination, linked to whole blood through a site-of-action effect",
    "compartment with a plasma-to-whole-blood transfer rate constant (Kbp)",
    "and a whole-blood-to-plasma concentration ratio (Rbp). Baseline blood",
    "volume scales the apparent central volume, time-varying hematocrit and",
    "nominal dose scale Rbp, and a concomitant weak CYP3A4 inducer raises",
    "apparent clearance; between-occasion variability on clearance over 13",
    "sampling occasions."
  )
  reference <- paste(
    "Savic RM, Green ML, Jorga K, Zager M, Washington CB. Model-informed",
    "drug development of voxelotor in sickle cell disease: Population",
    "pharmacokinetics in whole blood and plasma. CPT Pharmacometrics Syst",
    "Pharmacol. 2022;11(6):687-697. doi:10.1002/psp4.12731"
  )
  vignette <- "Savic_2022_voxelotor"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    BLOOD_VOLUME = list(
      description = "Baseline total blood volume",
      units = "L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline (time-fixed) estimated total blood volume. Savic 2022",
        "states only that blood volume was 'calculated based on body weight",
        "and sex' and does not print the formula, so supply the column",
        "directly. Enters as a power function on Vc/F with reference 3.89 L",
        "(Table 2 footnote c); cohort median 3.9 L, 10th-90th percentile",
        "2.9-5.2 L (Results), range 2-7 L (Table 1)."
      ),
      source_name = "BLV"
    ),
    HCT = list(
      description = "Hematocrit",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying hematocrit on the percent scale. Enters as a power",
        "function on Rbp with reference 27.8 % (Table 2). Baseline median",
        "27.0 % (range 17-40, Table 1); median time-varying value 30.5 %,",
        "10th-90th percentile 24.3-37.7 % (Results)."
      ),
      source_name = "HCT"
    ),
    DOSE_VOXELOTOR_MG = list(
      description = "Nominal voxelotor dose of the subject's regimen",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Nominal dose level in mg (500, 600, 700, 900, 1000 or 1500 in the",
        "source studies; Figure 2c). Enters as a power function on Rbp with",
        "reference 900 mg (Table 2 'Nominal dose on Rbp, (dose/900)^TH').",
        "For the once-daily regimens that make up almost all of the data",
        "this is the amount per administration; the paper does not say how",
        "the FIH 500 mg twice-daily arm was coded."
      ),
      source_name = "dose"
    ),
    CONMED_CYP3A4_IND_WEAK = list(
      description = "Concomitant weak CYP3A4 inducer (time-varying)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant weak CYP3A4 inducer)",
      notes = paste(
        "1 = weak CYP3A4 inducer coadministered at the record. Strong",
        "inducers were prohibited in the studies and no subject received a",
        "moderate inducer (Table 1: weak 14 of 264 PK-evaluable patients,",
        "moderate 0; Discussion cites nine concomitant users). The paper",
        "does not name the agents. Do not apply the coefficient to strong",
        "inducers, which the Discussion expects to reduce exposure much more."
      ),
      source_name = "time-varying weak CYP3A4 inducer"
    ),
    OCC = list(
      description = "Sampling occasion (1-13) for between-occasion variability on CL/F",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Integer occasion index 1..13; occasions were 'defined based on",
        "sampling visits' (Results) and Table 2 footnote d counts 13 of them",
        "(from day 25 in the FIH study to week 72 in HOPE). The visit-to-",
        "occasion mapping is not printed. Decomposed inside model() into",
        "indicators that multiplex 13 IOV etas sharing one variance; a",
        "record with OCC outside 1..13 receives no IOV."
      ),
      source_name = "OCC"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "voxelotor", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "voxelotor", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "voxelotor", units = "mg", specimen = "plasma", verified = TRUE),
    effect = list(analyte = "voxelotor", units = "ug/mL", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 279L,
    n_studies = 3L,
    age_range = "12-59 years",
    age_median = "22 years",
    weight_range = "28-135 kg",
    weight_median = "61 kg",
    sex_female_pct = 58,
    race_ethnicity = c(Black = 72, White = 10, `Arab/Middle Eastern` = 12, `Other/multiple/missing` = 6),
    disease_state = "Sickle cell disease (HbSS 79%, HbSbeta0 13%, HbSC 4%, HbSbeta+ 3%), Hb about 6-10.5 g/dL at entry, 62% on hydroxyurea",
    dose_range = "500-1500 mg oral once daily (also 600 or 1000 mg single doses and 500 mg twice daily) for up to 72 weeks",
    regions = "Multinational (FIH study, HOPE Kids 1, HOPE)",
    notes = paste(
      "Savic 2022 Table 1 (N = 279; 76 adolescents 12 to <18 years, 203",
      "adults). Final PK-evaluable dataset 264 patients with 2155 plasma and",
      "2168 whole-blood observations. Studies: FIH phase I/II in adults",
      "(NCT02285088), HOPE Kids 1 phase IIa in 12-17-year-olds",
      "(NCT02850406), HOPE phase III (NCT03036813). Baseline medians: Hb",
      "9 g/dL, HCT 27.0 %, blood volume 4 L, albumin 43 g/L."
    )
  )

  ini({
    # Structural parameters (Table 2; reported as exp(TH) for MU-referenced thetas)
    lka <- fixed(log(2.38)); label("Absorption rate constant Ka (1/h)") # Table 2 'Ka (1/h)' 2.38 (FIXED)
    lcl <- log(6.14); label("Apparent clearance CL/F (L/h)") # Table 2 'CL/F (L/h)' 6.14 (RSE 2.8%)
    lvc <- log(333); label("Apparent central volume Vc/F at blood volume 3.89 L (L)") # Table 2 'Vc/F (L)' 333 (RSE 0.5%)
    lq <- log(0.39); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2 'Q/F (L/h)' 0.39 (RSE 1.9%)
    lvp <- log(72.3); label("Apparent peripheral volume Vp/F (L)") # Table 2 'Vp/F (L)' 72.3 (RSE 0.8%)
    lke0 <- log(0.43); label("Plasma-to-whole-blood transfer rate constant Kbp (1/h)") # Table 2 'Kbp (1/h)' 0.43 (RSE 6.6%)
    lbpr <- log(16.6); label("Whole-blood-to-plasma concentration ratio Rbp at HCT 27.8 and 900 mg (unitless)") # Table 2 'Rbp' 16.6 (RSE 0.5%)

    # Covariate effects (Table 2)
    e_blood_volume_vc <- 0.74; label("Power exponent of baseline blood volume (BLOOD_VOLUME/3.89) on Vc/F (unitless)") # Table 2 'Blood volume on Vc/F, (BLV/3.89)^TH' 0.74 (RSE 14.1%)
    e_hct_bpr <- 0.77; label("Power exponent of time-varying hematocrit (HCT/27.8) on Rbp (unitless)") # Table 2 'Hematocrit on Rbp, (HCT/27.8)^TH' 0.77 (RSE 8.6%)
    e_cyp3a4_ind_weak_cl <- 0.39; label("Log-scale effect of a concomitant weak CYP3A4 inducer on CL/F, exp(TH) (unitless)") # Table 2 'CYP3A4 inducer on CL/F, exp TH' 0.39 (RSE 4.9%)
    e_dose_bpr <- -0.37; label("Power exponent of nominal dose (DOSE_VOXELOTOR_MG/900) on Rbp (unitless)") # Table 2 'Nominal dose on Rbp, (dose/900)^TH' -0.37 (RSE 10.8%)

    # Between-subject variability. Table 2 footnote a: %CV = 100 * sqrt(omega),
    # so omega = (CV/100)^2; covariance = correlation * SD1 * SD2.
    # CL/F 34.1% -> 0.116281; Vc/F 21.7% -> 0.047089; corr 0.1 -> 0.1*0.341*0.217 = 0.0073997
    etalcl + etalvc ~ c(0.116281, 0.0073997, 0.047089) # Table 2 'BSV CL/F' 34.1, 'BSV Vc/F' 21.7, correlation 0.1
    # Kbp 43.8% -> 0.191844; Rbp 15.2% -> 0.023104; corr -0.102 -> -0.102*0.438*0.152 = -0.0067908
    etalke0 + etalbpr ~ c(0.191844, -0.0067908, 0.023104) # Table 2 'BSV Kbp' 43.8, 'BSV Rbp' 15.2, correlation -0.102

    # Between-occasion variability on CL/F: 56.3% CV -> 0.563^2 = 0.316969,
    # 13 occasions (Table 2 footnote d) sharing one variance. Occasion 1
    # carries the estimate; occasions 2-13 are fixed to the same value.
    etaiov_cl_1 ~ 0.316969 # Table 2 'BOV on CL/F, %CV' 56.3 (estimated)
    etaiov_cl_2 ~ fixed(0.316969) # same variance as occasion 1
    etaiov_cl_3 ~ fixed(0.316969) # same variance as occasion 1
    etaiov_cl_4 ~ fixed(0.316969) # same variance as occasion 1
    etaiov_cl_5 ~ fixed(0.316969) # same variance as occasion 1
    etaiov_cl_6 ~ fixed(0.316969) # same variance as occasion 1
    etaiov_cl_7 ~ fixed(0.316969) # same variance as occasion 1
    etaiov_cl_8 ~ fixed(0.316969) # same variance as occasion 1
    etaiov_cl_9 ~ fixed(0.316969) # same variance as occasion 1
    etaiov_cl_10 ~ fixed(0.316969) # same variance as occasion 1
    etaiov_cl_11 ~ fixed(0.316969) # same variance as occasion 1
    etaiov_cl_12 ~ fixed(0.316969) # same variance as occasion 1
    etaiov_cl_13 ~ fixed(0.316969) # same variance as occasion 1

    # Residual error (Table 2)
    propSd <- 0.24; label("Proportional residual error, plasma (fraction)") # Table 2 'Proportional error, plasma (%)' 24.0 (RSE 2.1%)
    propSd_Cblood <- 0.17; label("Proportional residual error, whole blood (fraction)") # Table 2 'Proportional error, whole blood (%)' 17.0 (RSE 3.2%)
    addSd_Cblood <- 0.88; label("Additive residual error, whole blood (ug/mL)") # Table 2 'Additive error, whole blood (ng/ml)' 880 (RSE 12.5%); 880 ng/mL = 0.88 ug/mL
  })

  model({
    # Occasion indicators for the between-occasion variability on CL/F
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)
    oc7 <- (OCC == 7)
    oc8 <- (OCC == 8)
    oc9 <- (OCC == 9)
    oc10 <- (OCC == 10)
    oc11 <- (OCC == 11)
    oc12 <- (OCC == 12)
    oc13 <- (OCC == 13)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5 + oc6 * etaiov_cl_6 +
      oc7 * etaiov_cl_7 + oc8 * etaiov_cl_8 + oc9 * etaiov_cl_9 +
      oc10 * etaiov_cl_10 + oc11 * etaiov_cl_11 + oc12 * etaiov_cl_12 +
      oc13 * etaiov_cl_13

    # Individual parameters. Covariate forms per the Table 2 row labels and
    # the Table 2 note (continuous covariates as power functions).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl + iov_cl + e_cyp3a4_ind_weak_cl * CONMED_CYP3A4_IND_WEAK)
    vc <- exp(lvc + etalvc) * (BLOOD_VOLUME / 3.89)^e_blood_volume_vc
    q <- exp(lq)
    vp <- exp(lvp)
    ke0 <- exp(lke0 + etalke0)
    bpr <- exp(lbpr + etalbpr) * (HCT / 27.8)^e_hct_bpr * (DOSE_VOXELOTOR_MG / 900)^e_dose_bpr

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment disposition with first-order absorption (Figure 1)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Plasma concentration: dose in mg, volume in L, so mg/L = ug/mL
    Cc <- central / vc

    # Site-of-action (effect) compartment linking plasma to whole blood
    # (Figure 1): the effect concentration equilibrates with plasma at rate
    # Kbp and the whole-blood concentration is Rbp times that concentration.
    d/dt(effect) <- ke0 * (Cc - effect)
    Cblood <- bpr * effect

    Cc ~ prop(propSd)
    Cblood ~ add(addSd_Cblood) + prop(propSd_Cblood)
  })
}
