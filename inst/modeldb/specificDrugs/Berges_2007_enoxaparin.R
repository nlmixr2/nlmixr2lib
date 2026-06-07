Berges_2007_enoxaparin <- function() {
  description <- "Two-compartment first-order absorption population PK model of anti-factor Xa activity in elderly patients (>75 years) receiving prophylactic subcutaneous enoxaparin 4000 IU once daily (Berges 2007 PROPHRE.75 study)"
  reference <- "Berges A, Laporte S, Epinat M, Zufferey P, Alamartine E, Tranchand B, Decousus H, Mismetti P. Anti-factor Xa activity of enoxaparin administered at prophylactic dosage to patients over 75 years old. Br J Clin Pharmacol. 2007;64(4):428-438. doi:10.1111/j.1365-2125.2007.02920.x"
  vignette <- "Berges_2007_enoxaparin"
  units <- list(time = "hour", dosing = "IU", concentration = "IU/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight (baseline)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference weight 65 kg (PROPHRE.75 population median). Used with power scaling on CL (exponent 0.78) and Vc (exponent 1.25). Mean 66 kg, range 38-108 kg per Berges 2007 Table 1.",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Creatinine clearance estimated by the simplified Modification of Diet in Renal Disease (MDRD / simplified Levey) formula",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference value 69 mL/min (PROPHRE.75 population median for simplified-MDRD CrCl). Stored under the canonical CRCL column per inst/references/covariate-columns.md; the simplified-MDRD output is conceptually in mL/min/1.73 m^2 but Berges 2007 reports it in mL/min (Table 1: mean 69, range 27-127). Used with power scaling on CL (exponent 0.25). Note that the paper retained simplified-MDRD CrCl in preference to Cockcroft-Gault CrCl (also reported in Table 1) because it produced the larger reduction in objective function during forward selection (Table 2). The Cockcroft-Gault alternative is documented for reference only.",
      source_name        = "CrCl_MDRD"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 189L,
    n_studies      = 1L,
    age_range      = "75-95 years",
    age_median     = "81 years",
    age_mean_sd    = "82 +/- 5 years",
    weight_range   = "38-108 kg",
    weight_median  = "65 kg",
    weight_mean_sd = "66 +/- 14 kg",
    sex_female_pct = 62,
    race_ethnicity = "Not reported (French university-hospital cohort)",
    disease_state  = "Elderly (>75 years) hospitalized patients receiving venous thromboembolism (VTE) prophylaxis; 63% immobile due to acute medical disease, 15% orthopaedic surgery, 22% stroke; 50% with moderate or severe renal failure (CrCl < 50 mL/min by Cockcroft-Gault) or 18% (by simplified MDRD); 24% diabetes, 19% cancer, 57% hypertension; 36% concomitant antiplatelet agents",
    dose_range     = "4000 IU subcutaneous once daily",
    regions        = "France (University Hospital of Saint-Etienne)",
    renal_function = "Simplified-MDRD CrCl median 69 mL/min (range 27-127); Cockcroft-Gault CrCl median 52 mL/min (range 24-93)",
    notes          = "Baseline demographics per Berges 2007 Table 1 (PROPHRE.75 study, prospective cohort, mean treatment duration 7 days). A total of 451 anti-factor Xa activity samples (mean 2.4 per patient, very sparse) were analysed; 56 below-LOQ values (LOQ 0.05 IU/mL) were excluded from estimation."
  )

  ini({
    # Structural parameters at the population median (WT = 65 kg, CRCL = 69 mL/min);
    # Berges 2007 Table 3 final-parameter-estimates column.
    lka <- log(0.63); label("Absorption rate constant (Ka, 1/h)")                            # Berges 2007 Table 3: KA = 0.63 (95% CI 0.44, 0.81)
    lcl <- log(0.70); label("Apparent clearance at WT = 65 kg and CRCL = 69 mL/min (CL/F, L/h)") # Berges 2007 Table 3: theta1 = 0.70 (95% CI 0.66, 0.75)
    lvc <- log(6.43); label("Apparent central volume at WT = 65 kg (V2/F, L)")               # Berges 2007 Table 3: theta2 = 6.43 (95% CI 5.47, 7.39)
    lq  <- log(0.34); label("Apparent intercompartmental clearance (Q/F, L/h)")              # Berges 2007 Table 3: Q = 0.34 (95% CI 0.17, 0.49)
    lvp <- log(8.18); label("Apparent peripheral volume (V3/F, L)")                          # Berges 2007 Table 3: V3 = 8.18 (95% CI 1.97, 14.36)

    # Allometric / power covariate effects (Berges 2007 Table 3).
    e_wt_cl   <- 0.78; label("Power exponent of (WT / 65) on CL/F (unitless)")                # Berges 2007 Table 3: theta6 = 0.78 (95% CI 0.47, 1.08)
    e_crcl_cl <- 0.25; label("Power exponent of (CRCL / 69) on CL/F (unitless)")              # Berges 2007 Table 3: theta7 = 0.25 (95% CI 0.05, 0.45)
    e_wt_vc   <- 1.25; label("Power exponent of (WT / 65) on V2/F (unitless)")                # Berges 2007 Table 3: theta8 = 1.25 (95% CI 0.72, 1.78)

    # Inter-individual variability. The paper reports IIV as %CV from a log-normal
    # exponential parameterisation; omega^2 = log(CV^2 + 1).
    # The paper states (Results > Model building) that a block matrix was added to
    # take into account the correlation between CL and V2, but the covariance value
    # itself is not reported in Table 3 (only the diagonal CV%s). The model is
    # encoded with diagonal etas (no off-diagonal); see vignette Assumptions and
    # deviations for the deviation from the published block structure.
    etalcl ~ 0.065361 # log(0.26^2 + 1); 26% CV on CL/F (95% CI 20, 31), Berges 2007 Table 3
    etalvc ~ 0.022275 # log(0.15^2 + 1); 15% CV on V2/F (95% CI  0, 23), Berges 2007 Table 3
    etalvp ~ 0.622924 # log(0.93^2 + 1); 93% CV on V3/F (95% CI 22, 130), Berges 2007 Table 3
    # IIV on KA and Q were fixed to 0 (Berges 2007 Table 3: "0 FIXED").

    # Residual error (Berges 2007 Table 3: proportional model, CV 30% (95% CI 26, 33)).
    propSd <- 0.30; label("Proportional residual error on anti-Xa activity (fraction)")     # Berges 2007 Table 3: sigma = 30% (95% CI 26, 33)
  })
  model({
    # Individual PK parameters. CL has power-style covariate effects of WT and
    # CRCL (each ratioed to the population median); V2 has a power-style WT effect.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 65)^e_wt_cl * (CRCL / 69)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 65)^e_wt_vc
    q  <- exp(lq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # Dose in IU, volumes in L. central/vc has units IU/L; divide by 1000 to
    # convert to IU/mL (the paper's reporting unit for anti-Xa activity).
    Cc <- central / (vc * 1000)
    Cc ~ prop(propSd)
  })
}
