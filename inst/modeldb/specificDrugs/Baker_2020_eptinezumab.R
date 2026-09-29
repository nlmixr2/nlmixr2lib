Baker_2020_eptinezumab <- function() {
  description <- "Two-compartment population PK model for intravenous eptinezumab (humanized IgG1 anti-CGRP mAb) in healthy adults and adults with episodic or chronic migraine, with linear elimination; body weight on all four disposition parameters (shared exponent on CL/Q and on Vc/Vp), disease state (healthy/episodic/chronic migraine), capped creatinine clearance and baseline monthly migraine days on CL, and disease state and sex on Vc (Baker 2020)"
  reference <- paste(
    "Baker B, Schaeffler B, Pederson S, Trinh M, Smith J, Latham J, Beliveau M, Rubets I.",
    "Population pharmacokinetic and exposure-response analysis of eptinezumab in the treatment of episodic and chronic migraine.",
    "Pharmacol Res Perspect. 2020;8(2):e00567. doi:10.1002/prp2.567.",
    "Covariate coefficients, Q, Vp and the CLp/Vp BSV are from the FDA Clinical Pharmacology Review of BLA 761119 (Vyepti), Table 12",
    "(reproducing sponsor report ALD403-088-PK Table 7), and the shared weight exponents from the EMA Vyepti assessment report (EMA/9446/2022, procedure EMEA/H/C/005287/0000) Table 6."
  )
  vignette <- "Baker_2020_eptinezumab"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power model (WT/70)^theta. One exponent (0.709) is shared by CL and the intercompartmental clearance Q, and one (0.544) by Vc and Vp: EMA EPAR Table 6 adds each pair as a single fitted parameter (steps 1 and 2, #Param 11 -> 12 -> 13).",
      source_name = "WT"
    ),
    CRCL = list(
      description = "Creatinine clearance (Cockcroft-Gault), absolute, capped at 150 mL/min",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters CL as (min(CRCL, 150)/118)^0.162. The cap is applied inside model(), so supply the uncapped Cockcroft-Gault value. Not BSA-normalized (the FDA review names the Cockcroft-Gault and MDRD markers screened; the cap and the 118 mL/min reference are in mL/min).",
      source_name = "CLcr_cap"
    ),
    DIS_MIGRAINE_EPISODIC = list(
      description = "Episodic migraine patient indicator (1 = episodic migraine, including the frequent-episodic-migraine studies CLIN-002 and CLIN-006; 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy participant, when DIS_MIGRAINE_CHRONIC is also 0)",
      notes = "Mutually exclusive with DIS_MIGRAINE_CHRONIC; healthy participants have both indicators 0. Exponentiated factor exp(theta) on CL and on Vc.",
      source_name = "DS = EM"
    ),
    DIS_MIGRAINE_CHRONIC = list(
      description = "Chronic migraine patient indicator (1 = chronic migraine; 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy participant, when DIS_MIGRAINE_EPISODIC is also 0)",
      notes = "Mutually exclusive with DIS_MIGRAINE_EPISODIC. Exponentiated factor exp(theta) on CL and on Vc.",
      source_name = "DS = CM"
    ),
    MIGRAINE_DAYS_BL = list(
      description = "Baseline number of monthly migraine days (screening period)",
      units = "days/month",
      type = "continuous",
      reference_category = NULL,
      notes = "Power model (MIGRAINE_DAYS_BL/13)^0.044 on CL; 13 days is the population median. Healthy participants have no migraine days; the paper does not say what value they were given, so the typical healthy reference subject of the forest plot (Figure 4) is simulated at the 13-day reference.",
      source_name = "MDBASE"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (female)",
      notes = "Source coefficient is for male sex relative to female; applied as exp(e_male_vc * (1 - SEXF)) so females (the reference) have factor 1.",
      source_name = "SEX (Male)"
    )
  )

  compartmentData <- list(
    central = list(analyte = "eptinezumab", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "eptinezumab", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 2123,
    n_studies = 8,
    age_range = "18-71 years",
    age_median = "39.0 years",
    weight_range = "39.2-190 kg",
    weight_median = "74.2 kg",
    sex_female_pct = 83.8,
    race_ethnicity = c(White = 88.5),
    disease_state = "Healthy adults (n = 83, including healthy overweight/obese and type 1 diabetes participants) and adults with episodic (n = 727) or chronic (n = 1313) migraine.",
    dose_range = "1-1000 mg IV infused over ~30 min to 1 h (5 subjects ~2 h); single doses, or every 12 weeks in the phase 3 studies (100 and 300 mg).",
    renal_function = "Normal 55.3%, mild decrease in eGFR 41.9%, moderate decrease 2.7%.",
    ada_negative_pct = 83.8,
    regions = "Pooled analysis of studies CLIN-001, -002, -005, -006, -010, -011, -012 and -013; regions not reported.",
    analysis_dataset = "15135 quantifiable plasma concentrations of free eptinezumab (FDA review section 4.3).",
    notes = "Demographics from Baker 2020 Results section 3.1 and Table 1; disease-group counts from the FDA Clinical Pharmacology Review of BLA 761119."
  )

  ini({
    # Structural parameters, typical healthy female subject, WT 70 kg,
    # capped CrCl 118 mL/min, baseline MMD 13 days.
    lcl <- log(0.00620); label("Clearance (L/h)") # Baker 2020 Results 3.2 = 0.00620 L/h; FDA review Table 12 'CL (L/h)' = 0.00620
    lvc <- log(3.636); label("Central volume of distribution (L)") # FDA review Table 12 'Vc (L)' = 3.636; Baker 2020 Results 3.2 = 3.64 L
    lq <- log(0.039); label("Intercompartmental clearance (L/h)") # FDA review Table 12 'CLp (L/h)' = 0.039
    lvp <- log(2.012); label("Peripheral volume of distribution (L)") # FDA review Table 12 'Vd (L)' = 2.012 (peripheral volume per table footnote)

    # Body weight, power model (WT/70)^theta; exponents shared per EMA EPAR Table 6
    e_wt_cl_q <- 0.709; label("Body-weight exponent shared by CL and Q (unitless)") # FDA review Table 12 'x (WT/70)^theta' on CL = 0.709
    e_wt_vc_vp <- 0.544; label("Body-weight exponent shared by Vc and Vp (unitless)") # FDA review Table 12 'x (WT/70)^theta' on Vc = 0.544

    # Other covariate effects on CL
    e_dis_migraine_episodic_cl <- -0.231; label("Episodic migraine effect on CL, exponentiated (unitless)") # FDA review Table 12 'x theta if EM' on CL = -0.231
    e_dis_migraine_chronic_cl <- -0.272; label("Chronic migraine effect on CL, exponentiated (unitless)") # FDA review Table 12 'x theta if CM' on CL = -0.272
    e_migraine_days_bl_cl <- 0.044; label("Baseline monthly migraine days exponent on CL (unitless)") # FDA review Table 12 'x (MDBASE/13.0)^theta' = 0.044
    e_crcl_cl <- 0.162; label("Capped creatinine clearance exponent on CL (unitless)") # FDA review Table 12 'x (CLcr_cap/118)^theta' = 0.162

    # Other covariate effects on Vc
    e_dis_migraine_episodic_vc <- -0.311; label("Episodic migraine effect on Vc, exponentiated (unitless)") # FDA review Table 12 'x theta if EM' on Vc = -0.311
    e_dis_migraine_chronic_vc <- -0.422; label("Chronic migraine effect on Vc, exponentiated (unitless)") # FDA review Table 12 'x theta if CM' on Vc = -0.422
    e_male_vc <- 0.091; label("Male sex effect on Vc, exponentiated (unitless)") # FDA review Table 12 'x theta if Male' on Vc = 0.091

    # Between-subject variability. Table 12 Note 1: CV% = sqrt(exp(omega) - 1),
    # so omega = log(CV^2 + 1). EPAR Table 6 base model has a CLp-Vp eta block;
    # its covariance is not printed in any source and is set to 0 here.
    etalcl ~ 0.08075 # Baker 2020 Results 3.2 BSV CL 29%; FDA review Table 12 = 29.0% -> log(0.290^2 + 1)
    etalvc ~ 0.09176 # Baker 2020 Results 3.2 BSV Vc 31%; FDA review Table 12 = 31.0% -> log(0.310^2 + 1)
    etalq ~ 0.80594 # FDA review Table 12 BSV CLp = 111.3% -> log(1.113^2 + 1)
    etalvp ~ 0.11431 # FDA review Table 12 BSV Vd = 34.8% -> log(0.348^2 + 1)

    # Residual error: C = f * (1 + eps_p) + eps_a (Baker 2020 section 2.2.1 equation)
    addSd <- 0.0375; label("Additive residual error (ug/mL)") # Baker 2020 Results 3.2 = 37.5 ng/mL; FDA review Table 12 'Additive (ng/mL)' = 37.5
    propSd <- 0.252; label("Proportional residual error (fraction)") # FDA review Table 12 'Prop Error (%)' = 25.2
  })
  model({
    # Covariate terms (FDA review Table 12 footnote: power models standardized
    # by the median for continuous covariates; exponentiated factors for
    # categorical covariates relative to the reference category)
    crcl_cap <- min(CRCL, 150)
    wt_cl <- (WT / 70)^e_wt_cl_q
    wt_v <- (WT / 70)^e_wt_vc_vp
    dis_cl <- exp(e_dis_migraine_episodic_cl * DIS_MIGRAINE_EPISODIC + e_dis_migraine_chronic_cl * DIS_MIGRAINE_CHRONIC)
    dis_vc <- exp(e_dis_migraine_episodic_vc * DIS_MIGRAINE_EPISODIC + e_dis_migraine_chronic_vc * DIS_MIGRAINE_CHRONIC)

    cl <- exp(lcl + etalcl) * wt_cl * dis_cl * (MIGRAINE_DAYS_BL / 13)^e_migraine_days_bl_cl * (crcl_cap / 118)^e_crcl_cl
    vc <- exp(lvc + etalvc) * wt_v * dis_vc * exp(e_male_vc * (1 - SEXF))
    q <- exp(lq + etalq) * wt_cl
    vp <- exp(lvp + etalvp) * wt_v

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volume in L -> mg/L = ug/mL
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
