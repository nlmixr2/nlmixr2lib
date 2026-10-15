Han_2022_cabotegravir <- function() {
  description <- "Two-compartment population PK model for cabotegravir after once-daily oral tablets and long-acting (LA) intramuscular injection in adults with and without HIV-1 (Han 2022, 16 pooled phase 1-3 studies). Oral and LA doses are absorbed first-order from separate depots, with oral bioavailability expressed relative to the LA injection (LA absolute bioavailability is the apparent-parameter reference). Baseline body weight scales CL/F and Q/F with one shared estimated exponent and V2/F and V3/F with another, at a 76.6 kg reference. Current smoking raises CL/F by 17.4%. The LA absorption rate constant is 50.9% lower in women, 47.8% higher when the dose is given as two split injections, falls with body mass index (power -0.766 at 25.4 kg/m^2) and rises with needle length (power 0.478 at 1.5 inch); because LA absorption is rate-limiting (flip-flop), a slower LA absorption lowers the first trough but raises the steady-state trough."
  reference <- "Han K, Baker M, Lovern M, Paul P, Xiong Y, Patel P, Moore KP, Seal CS, Cutrell AG, D'Amico RD, Benn PD, Landovitz RJ, Marzinke MA, Spreen WR, Ford SL. Population pharmacokinetics of cabotegravir following administration of oral tablet and long-acting intramuscular injection in adult HIV-1-infected and uninfected subjects. Br J Clin Pharmacol. 2022;88(10):4607-4622. doi:10.1111/bcp.15439"
  vignette <- "Han_2022_cabotegravir"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Han 2022 used BASELINE body weight (BWT) in the final model; supply a per-subject constant. Power function centred at the model-building median 76.6 kg (Table 2; the reference is printed in the final-model equations). One estimated exponent is shared by CL/F and Q/F (0.618) and another by V2/F and V3/F (0.702), as Methods 2.2 states. Model-building cohort 41.2-168.3 kg.",
      source_name = "BWT"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Fractional change on the LA absorption rate constant, KALA x (1 - 50.9%) in women (Table 3 'Female on KALA 50.9%'; final-model equation '(1 - 50.9% if female)'). No transgender subjects were enrolled (Results 3.1). 25.7% of the model-building cohort was female.",
      source_name = "female"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power function on KALA centred at the population median 25.4 kg/m^2 (final-model equation; Results 3.2 'population median BMI (25.4 kg/m2)'). The exponent is NEGATIVE (-0.766): Table 3 and the bootstrap CI print without minus signs in the PDF text layer, but Discussion states KALA is '28.6% lower in subjects with BMI of 32.5 kg/m2 (90th percentile) than those with BMI of 20.9 kg/m2 (10th percentile)', and (32.5/20.9)^-0.766 = 0.713 reproduces that 28.7% reduction. Model-building cohort 15.3-69.51 kg/m^2 (Table 2).",
      source_name = "BMI"
    ),
    SMOKE = list(
      description = "Current smoker, 1 = current smoker, 0 = never, former, 'not current' or unknown",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not a current smoker)",
      notes = "Fractional change on CL/F, CL/F x (1 + 17.4%) in current smokers (Table 3 'Smoke on CL/F 17.4%'; final-model equation '(1 + 17.4% if current smoker)'). Han 2022 Table 2 records five smoking categories (never, former, current, 'not current', unknown); the final model contrasts current smokers with everyone else, which is exactly the two-level SMOKE encoding. 27.1% of the model-building cohort were current smokers.",
      source_name = "Smoke"
    ),
    INJ_SPLIT = list(
      description = "Long-acting dose given as two split intramuscular injections, 1 = split, 0 = single (unsplit) injection",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (unsplit injection)",
      notes = "Per-injection covariate: carry it on the LA dose record and forward in time (Han 2022 used last observation carried forward for time-varying covariates, Methods 2.3). Fractional change on KALA, KALA x (1 + 47.8%) for a split dose (Table 3 'Split on KALA 47.8%'). 726 of 12 169 injections (6.0%) were split in the model-building dataset and none in the ATLAS-2M validation set (Table 2). Value has no effect on oral-dosing intervals because KALA governs only the LA depot.",
      source_name = "Split"
    ),
    NEEDLE_LEN_IN = list(
      description = "Length of the needle used for the long-acting intramuscular injection",
      units = "inch",
      type = "continuous",
      reference_category = NULL,
      notes = "Per-injection covariate, carried forward from the LA dose record. Supply INCHES, as Han 2022 records it (Table 2 'Injection needle length (inch)': 1, 1.5, 2, 2.5); the power function is centred at 1.5 inch (Results 3.2 equation '(NDL/1.5)^0.478'), so a value in millimetres would give a wildly wrong absorption rate. Exponent 0.478 (Table 3 'NDL exponent on KALA'). 88.5% of injections used the standard 1.5-inch needle; 4.0% had unknown length, which Han 2022 imputed (Methods 2.3). The Figure 5 covariate simulations assigned 1.5 inch to every injection (Discussion).",
      source_name = "NDL"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened and not retained; Han 2022 Abstract and Conclusions: 'Race and age were not significant covariates'. Model-building cohort 18-74 years (Table 2).",
      source_name = "Age"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Significant in forward addition and backward elimination on CL/F and V2/F, then removed in the final refinement step because the effect on a typical subject was below the 15% threshold between the 10th and 90th percentiles (Discussion and Methods 2.3). No coefficient is reported.",
      source_name = "albumin"
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Significant on CL/F in the stepwise search but removed during model refinement under the 15% clinical-relevance rule (Discussion). No coefficient is reported.",
      source_name = "total bilirubin"
    ),
    HIV_POS = list(
      description = "HIV-1 infection status, 1 = HIV-1-infected",
      units = "(binary)",
      type = "binary",
      notes = "Significant on F1 in the stepwise search but removed during model refinement under the 15% clinical-relevance rule (Discussion). 72.2% of the model-building cohort was HIV-1-infected (Table 2).",
      source_name = "HIV-1 status"
    )
  )

  compartmentData <- list(
    depot_oral = list(analyte = "cabotegravir", units = "mg", specimen = "administration site", verified = TRUE),
    depot_im = list(analyte = "cabotegravir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "cabotegravir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cabotegravir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 1647L,
    n_studies = 16L,
    age_range = "18-74 years",
    age_median = "36 years",
    weight_range = "41.2-168.3 kg",
    weight_median = "76.6 kg",
    bmi_range = "15.3-69.51 kg/m^2",
    bmi_median = "25.37 kg/m^2",
    sex_female_pct = 25.7,
    race_ethnicity = c(White = 66.9, Black = 23.9, Asian = 3.4, Other = 5.6, Unknown = 0.12),
    disease_state = "Adults with HIV-1 infection (72.2%; treatment-naive, induction and virologically suppressed switch populations from LATTE, LATTE-2, FLAIR and ATLAS) and HIV-1-uninfected adults (27.8%; phase 1 volunteers and the HPTN 077 / ECLAIR pre-exposure-prophylaxis studies).",
    dose_range = "Oral cabotegravir sodium tablet 10, 30 or 60 mg once daily or as a single dose; long-acting intramuscular injection 200, 400, 600 or 800 mg, given every 4, 8 or 12 weeks after a 4-week (or longer) oral lead-in.",
    regions = "Multinational phase 1-3 programme (ViiV Healthcare / GSK, HPTN).",
    smoking = "Current smoker 27.1%, never 35.8%, former 12.8%, not current 7.3%, unknown 17.0%",
    n_concentrations = "23 926 plasma concentrations (537 BQL after LA injection, modelled with the M3 method)",
    external_validation = "ATLAS-2M (study 207966): 647 subjects, 5097 concentrations not used for model building, predicted adequately by pcVPC and NPC (Results 3.4)",
    notes = "Han 2022 Table 1 (studies) and Table 2 (baseline characteristics). Estimation in NONMEM 7.3, FOCE with interaction. The external-validation cohort is not part of n_subjects."
  )

  ini({
    # All point estimates are the 'Final model Estimate' column of Han 2022
    # Table 3 (bootstrap 90% CI from the same table quoted alongside). No
    # control stream was published; the supplement holds Figures S1-S4 only.
    #
    # IIV SCALE. Table 3 prints IIV as a percentage without saying how it was
    # computed. The register default for a log-normal IIV reported as a CV is
    # omega^2 = log(1 + CV^2), used here; the alternative omega = CV/100
    # reading gives a variance at most 36% larger (KAoral: 0.799 versus
    # 0.587) and is set out in the vignette Errata.

    # --- Absorption ------------------------------------------------------
    lka_oral <- log(1.41); label("Oral tablet absorption rate constant KAoral (1/h)") # Table 3 'KAoral (h-1) 1.41', 90% CI 1.23-1.52
    lka_im <- log(0.000733); label("LA intramuscular absorption rate constant KALA for a 25.4 kg/m^2 male, unsplit dose, 1.5-inch needle (1/h)") # Table 3 'KALA (h-1) 0.000733', 90% CI 0.000705-0.000761 (the equation in Results 3.2 rounds it to 0.00073)
    lfdepot_oral <- log(0.756); label("Oral tablet bioavailability relative to the LA intramuscular injection F1 (fraction)") # Table 3 'F1 75.6%', 90% CI 74.4%-76.8%

    # --- Disposition (typical subject 76.6 kg) ---------------------------
    lcl <- log(0.151); label("Apparent clearance CL/F at 76.6 kg (L/h)") # Table 3 'CL/F (L/h) 0.151', 90% CI 0.148-0.153
    lvc <- log(5.27); label("Apparent central volume V2/F at 76.6 kg (L)") # Table 3 'V2/F (L) 5.27', 90% CI 5.07-5.45
    lq <- log(0.507); label("Apparent intercompartmental clearance Q/F at 76.6 kg (L/h)") # Table 3 'Q/F (L/h) 0.507', 90% CI 0.455-0.579
    lvp <- log(2.43); label("Apparent peripheral volume V3/F at 76.6 kg (L)") # Table 3 'V3/F (L) 2.43', 90% CI 2.26-2.66

    # --- Covariate effects -----------------------------------------------
    e_wt_cl_q <- 0.618; label("Power exponent of baseline weight on CL/F and Q/F (unitless)") # Table 3 'BWT on CL/F and Q/F 0.618', 90% CI 0.562-0.68
    e_wt_vc_vp <- 0.702; label("Power exponent of baseline weight on V2/F and V3/F (unitless)") # Table 3 'BWT on V2/F and V3/F 0.702', 90% CI 0.622-0.79
    e_smoke_cl <- 0.174; label("Fractional change in CL/F for current smokers (unitless)") # Table 3 'Smoke on CL/F 17.4%', 90% CI 14.7%-20%
    e_sexf_ka_im <- -0.509; label("Fractional change in KALA for women (unitless)") # Table 3 'Female on KALA 50.9%', a reduction per the Results 3.2 equation '(1 - 50.9% if female)'
    e_bmi_ka_im <- -0.766; label("Power exponent of BMI on KALA, reference 25.4 kg/m^2 (unitless)") # Table 3 'BMI exponent on KALA 0.766'; negative sign from Discussion (KALA 28.6% lower at BMI 32.5 than 20.9)
    e_inj_split_ka_im <- 0.478; label("Fractional change in KALA for a split LA dose (unitless)") # Table 3 'Split on KALA 47.8%', 90% CI 36.9%-60.2%
    e_needle_len_in_ka_im <- 0.478; label("Power exponent of needle length on KALA, reference 1.5 inch (unitless)") # Table 3 'NDL exponent on KALA 0.478', 90% CI 0.184-0.747

    # --- Between-subject variability -------------------------------------
    # Table 3 reports IIV on every parameter except Q/F and V3/F. Results 3.2
    # states that off-diagonal covariances among KALA, CL/F and V2/F were
    # estimated, but their values are not published, so the block is encoded
    # diagonal (see vignette Errata).
    etalka_oral ~ 0.58736 # Table 3 'KAoral IIV 89.4%': log(1 + 0.894^2) = 0.58736; shrinkage 69.1%
    etalka_im ~ 0.28911 # Table 3 'KALA IIV 57.9%': log(1 + 0.579^2) = 0.28911; shrinkage 17.9%
    etalcl ~ 0.052867 # Table 3 'CL/F IIV 23.3%': log(1 + 0.233^2) = 0.052867; shrinkage 10.3%
    etalvc ~ 0.040383 # Table 3 'V2/F IIV 20.3%': log(1 + 0.203^2) = 0.040383; shrinkage 31.1%
    etalfdepot_oral ~ 0.029827 # Table 3 'F1 IIV 17.4%': log(1 + 0.174^2) = 0.029827; shrinkage 40%

    # --- Residual error (shared by oral and LA data, Results 3.2) --------
    addSd <- 0.0319; label("Additive residual error (ug/mL)") # Table 3 'Add err (ug/mL) 0.0319', 90% CI 0.0202-0.0496
    propSd <- 0.273; label("Proportional residual error (fraction)") # Table 3 'Prop err 27.3%', 90% CI 26.7%-27.8%
  })

  model({
    # Individual parameters. Covariate forms follow Han 2022 Methods 2.3:
    # power functions for continuous covariates and (1 + theta * IND) for
    # categorical ones; the final-model equations are in Results 3.2.
    ka_oral <- exp(lka_oral + etalka_oral)
    ka_im <- exp(lka_im + etalka_im) *
      (1 + e_sexf_ka_im * SEXF) *
      (1 + e_inj_split_ka_im * INJ_SPLIT) *
      (BMI / 25.4)^e_bmi_ka_im *
      (NEEDLE_LEN_IN / 1.5)^e_needle_len_in_ka_im
    cl <- exp(lcl + etalcl) * (WT / 76.6)^e_wt_cl_q * (1 + e_smoke_cl * SMOKE)
    vc <- exp(lvc + etalvc) * (WT / 76.6)^e_wt_vc_vp
    q <- exp(lq) * (WT / 76.6)^e_wt_cl_q
    vp <- exp(lvp) * (WT / 76.6)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # depot_oral = oral tablet; depot_im = long-acting intramuscular injection
    d/dt(depot_oral) <- -ka_oral * depot_oral
    d/dt(depot_im) <- -ka_im * depot_im
    d/dt(central) <- ka_oral * depot_oral + ka_im * depot_im - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Oral bioavailability relative to the LA injection (F1); the LA depot is
    # the bioavailability reference of the apparent parameters.
    f(depot_oral) <- exp(lfdepot_oral + etalfdepot_oral)

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
