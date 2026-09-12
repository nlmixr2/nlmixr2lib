Bertin_2026_propofol <- function() {
  description <- "Two-compartment population PK model for propofol in critically ill adults receiving a continuous intravenous infusion, half with and half without extracorporeal membrane oxygenation (ECMO), with a linear body-weight effect on clearance centred on 70 kg; ECMO was tested as a covariate on CL and V1 and was not retained" # nolint: line_length_linter.
  reference <- "Bertin S, Haefliger D, Mercier T, Decosterd LA, Giraud R, Assouline B, Schneider A, Buclin T, Guidi M, Livio F. Population Pharmacokinetics of Propofol in Critically Ill Patients with and Without Extracorporeal Membrane Oxygenation. Clin Pharmacokinet. 2026. doi:10.1007/s40262-025-01585-2. PMCID: PMC12881008" # nolint: line_length_linter.
  vignette <- "Bertin_2026_propofol"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  # CL and Q are published in L/h and V1/V2 in L, so time is hours. Doses are
  # continuous infusions in mg/kg/h plus occasional mg boluses. central holds
  # mg and vc is in L, so central/vc is mg/L -- the unit the HPLC-MS/MS assay
  # reports (calibration range 0.05-30 mg/L, LLOQ 0.1 mg/L, Sect. 2.2).

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The only covariate retained in the final model. It enters clearance as a",
        "LINEAR fractional deviation from a 70 kg reference, not as an allometric",
        "power: the Table 2 final-model equation is",
        "CL_i = CL * (1 + (BW - 70)/70 * theta_BW) * exp(eta_i^CL).",
        "The paper compared both forms and reports the linear model as slightly",
        "stronger (dOFV = -6.64, p = 0.01) than the allometric one",
        "(dOFV = -6.39, p < 0.05); the effect explained 15% of the CL",
        "between-subject variability (Sect. 3.1).",
        "The normalising constant is 70 kg, stated in the Table 2 equation and its",
        "legend ('CL, clearance of the reference patient of 70 kg') -- it is NOT the",
        "cohort median body weight, which was 78 kg on ECMO and 88 kg in controls",
        "(Table 1). Baseline, not time-varying, in the source analysis.",
        "Cohort range 51-120 kg (Table 1)."
      ),
      source_name = "BW"
    )
  )

  covariatesDataExcluded <- list(
    ECMO_STATUS = list(
      description = "Extracorporeal membrane oxygenation support, 1 = on ECMO, 0 = no ECMO",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no ECMO)",
      notes = paste(
        "The primary covariate of interest, and the reason the study was run.",
        "Screening initially associated ECMO with reduced CL",
        "(theta_ECMO = -0.2, dOFV = -4.43, p < 0.05) and reduced V1",
        "(theta_ECMO = -0.7, dOFV = -4.89, p < 0.05), but a sensitivity analysis",
        "excluding the two non-ECMO patients with V1 of 631 and 1040 L removed the",
        "significance of both (dOFV = -3.78, p > 0.05 for CL; dOFV = -1.0,",
        "p > 0.05 for V1). Those two control patients were judged to exert a",
        "leverage effect and ECMO was therefore NOT retained (Sect. 3.1).",
        "Documentation only -- deliberately not referenced in model(), because the",
        "published final model carries no ECMO term. 20 of 40 patients were on ECMO",
        "(16 veno-arterial, 3 veno-venous, 1 veno-arteriovenous; Table 1)."
      ),
      source_name = "ECMO"
    ),
    T_ECMO = list(
      description = "Time since ECMO initiation",
      units = "h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Tested for significance on the base-model parameters and not retained",
        "(Sect. 2.4.1 covariate list; no effect reported in Table 2 or in",
        "Supplementary Tables 2a/2b as summarised in Sect. 3.1). ECMO had been",
        "running for a median of 50 h (range 12-135 h) at the time of the first",
        "sample (Table 1), so the design could not resolve the transient",
        "post-cannulation clearance change reported by Morales et al. Documentation",
        "only -- not referenced in model()."
      ),
      source_name = "time since ECMO initiation"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Tested and not retained (Sect. 2.4.1). Cohort medians 28 g/L on ECMO and",
        "27 g/L in controls, range 19-37 g/L (Table 1). Documentation only -- not",
        "referenced in model()."
      ),
      source_name = "albumin"
    ),
    TBILI = list(
      description = "Total serum bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Tested and not retained (Sect. 2.4.1). Cohort medians 17 umol/L on ECMO",
        "and 16 umol/L in controls, range 4-238 umol/L (Table 1). Documentation",
        "only -- not referenced in model()."
      ),
      source_name = "bilirubin"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Tested and not retained (Sect. 2.4.1). No height summary is tabulated in",
        "Table 1. Documentation only -- not referenced in model()."
      ),
      source_name = "HT"
    ),
    AGE = list(
      description = "Chronological age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Tested and not retained (Sect. 2.4.1). Cohort medians 57 years on ECMO and",
        "56 years in controls, overall range 18-75 years (Table 1). Documentation",
        "only -- not referenced in model()."
      ),
      source_name = "age"
    ),
    SEXF = list(
      description = "Sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Tested as a dichotomous covariate and not retained (Sect. 2.4.1). The",
        "cohort was 87.5% male (35 of 40; Table 1), so the female stratum was too",
        "small to support an effect. Documentation only -- not referenced in",
        "model()."
      ),
      source_name = "sex"
    ),
    CIRRHOSIS = list(
      description = "Hepatic cirrhosis, 1 = present, 0 = absent",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no cirrhosis)",
      notes = paste(
        "Tested as a dichotomous covariate and not retained (Sect. 2.4.1). Only 4 of",
        "40 patients were cirrhotic (1 on ECMO, 3 controls; Table 1, graded by",
        "Child-Pugh score). Documentation only -- not referenced in model(), and",
        "deliberately NOT proposed as a canonical covariate column, because no",
        "retained model in the library uses it. The register's nearest existing",
        "entries are disease-specific (for example DIS_PBC, primary biliary",
        "cirrhosis) and do not match a general all-cause-cirrhosis indicator."
      ),
      source_name = "hepatic cirrhosis"
    )
  )

  compartmentData <- list(
    central = list(analyte = "propofol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "propofol", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 39,
    n_studies = 1,
    age_range = "18-75 years",
    age_median = "57 years (ECMO), 56 years (controls)",
    weight_range = "51-120 kg",
    weight_median = "78 kg (ECMO), 88 kg (controls)",
    sex_female_pct = 12.5,
    disease_state = "critically ill adults in intensive care receiving a continuous propofol infusion for sedation; 20 supported by ECMO (16 veno-arterial, 3 veno-venous, 1 veno-arteriovenous) and 20 matched non-ECMO controls", # nolint: line_length_linter.
    dose_range = "continuous intravenous infusion, median 2.2 mg/kg/h (range 0.6-4.6) on ECMO and 3.0 mg/kg/h (range 0.7-4.7) in controls, with intermittent boluses in 45% of ECMO and 70% of control patients", # nolint: line_length_linter.
    regions = "Switzerland (Lausanne CHUV and Geneva HUG)",
    renal_function = "6 of 40 on continuous veno-venous haemodialysis (3 ECMO, 3 control) and 3 controls on intermittent haemodialysis; median creatinine in the non-RRT patients 84 umol/L (ECMO) and 72 umol/L (controls)", # nolint: line_length_linter.
    hepatic_function = "4 of 40 with hepatic cirrhosis (1 ECMO, 3 control); median total bilirubin 17 umol/L (ECMO) and 16 umol/L (controls)",
    notes = paste(
      "Prospective bicentric observational study, January-December 2023",
      "(ethics project-ID 2022-01262). 40 patients (20 ECMO, 20 matched controls)",
      "contributed 300 plasma samples, up to 8 per patient over 9 h. Controls were",
      "matched to ECMO patients on sex, age band, body-weight band, renal",
      "replacement therapy, cirrhosis, total bilirubin band and cardiac function.",
      "The final popPK analysis used 289 concentrations from 39 patients: one",
      "control patient's 8 concentrations were removed for unexplained bias, two",
      "early concentrations from one ECMO patient were removed on CWRES grounds,",
      "and one sample collected after propofol discontinuation was removed",
      "(Sect. 3.1). The population metadata above therefore reports n_subjects = 39",
      "(the analysis dataset) while Table 1 describes all 40 enrolled patients.",
      "Baseline demographics and clinical data are in Table 1."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters, Bertin 2026 Table 2 (final-model column).
    # The reference individual weighs 70 kg -- stated both in the Table 2
    # equation, which divides (BW - 70) by 70, and in the Table 2 legend
    # ('CL, clearance of the reference patient of 70 kg'). Cross-checks
    # against the paper's own prose: V1 = 82 L is quoted as 1.2 L/kg and
    # V2 = 100 L as 1.4 L/kg (Sect. 4), both of which are the L values
    # divided by 70 kg rather than by the cohort median weight.
    # ------------------------------------------------------------------
    lcl <- log(68)
    label("Clearance for a 70 kg reference patient (L/h)") # Table 2 CL: 68 L/h (RSE 6.5%; bootstrap median 67, 95% CI 58-76)
    lvc <- log(82)
    label("Central volume of distribution (L)") # Table 2 V1: 82 L (RSE 33%; bootstrap median 83, 95% CI 28-169)
    lq <- log(26)
    label("Intercompartmental clearance (L/h)") # Table 2 Q: 26 L/h (RSE 16%; bootstrap median 26, 95% CI 5-124)
    lvp <- log(100)
    label("Peripheral volume of distribution (L)") # Table 2 V2: 100 L (RSE 27%; bootstrap median 101, 95% CI 14-377)

    # Covariate effect. LINEAR fractional slope, not an allometric
    # exponent: the Methods 'Linear equation' form is
    # Param_cov = Param * (1 + (COV - COVmed)/COVmed * theta), and Table 2
    # instantiates it with COVmed = 70 kg. Estimated, so not fixed().
    e_wt_cl <- 0.66
    label("Linear fractional effect of body weight on CL, applied as 1 + (WT - 70)/70 * e_wt_cl (unitless)") # Table 2 theta BW: 0.66 (RSE 38%; bootstrap median 0.69, 95% CI 0.22-1.40)

    # ------------------------------------------------------------------
    # Between-subject variability, on CL and V1 only. Table 2 reports both
    # as CV%, and footnote b gives the conversion the authors used:
    # CV% = sqrt(exp(omega^2) - 1). Inverting it, omega^2 = log(CV^2 + 1).
    # BSV on Q was tested and had no statistical impact (dOFV = 0,
    # p > 0.05); BSV on V2 gave only a modest improvement (dOFV = -3.93,
    # p = 0.05) with unsuccessful RSE estimation and was not retained
    # (Sect. 3.1). Neither is therefore an unreported variance -- the
    # published model deliberately has none, so no fixed(0) placeholder is
    # used for them.
    # ------------------------------------------------------------------
    etalcl ~ 0.109392 # Table 2 'omega CL' = 34 CV% (RSE 12%, shrinkage 2%); omega^2 = log(0.34^2 + 1) = 0.109392
    etalvc ~ 1.838961 # Table 2 'omega V1 (CV%)' = 230 (RSE 48%, shrinkage 39%); omega^2 = log(2.30^2 + 1) = 1.838961

    # Residual unexplained variability. A proportional model was selected
    # over additive and mixed alternatives (Sect. 2.4.1, 3.1).
    propSd <- 0.13
    label("Proportional residual error (fraction)") # Table 2 'sigma prop' = 13% (RSE 8.8%, shrinkage 9%; bootstrap median 13, 95% CI 11-15)
  })

  model({
    # 1. Individual parameters. Table 2 prints the final-model equation as
    #      CL_i = CL * (1 + (BW - 70)/70 * theta_BW) * exp(eta_i^CL)
    #    Numerical check against the paper's own worked example
    #    (Sect. 3.1): at 100 kg, 68 * (1 + 30/70 * 0.66) = 87.2 L/h, which
    #    the paper rounds to 87 L/h and describes as making the 70 kg
    #    clearance '22% lower'. 68/87.2 = 0.78, i.e. 22% lower. An
    #    un-normalised linear slope (68 + 0.66 * 30 = 87.8) or a
    #    cohort-median normalisation (86.2 at COVmed = 82 kg) would both
    #    round to a different number.
    cl <- exp(lcl + etalcl) * (1 + (WT - 70) / 70 * e_wt_cl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)

    # 2. Micro-constants for the two-compartment mammillary system
    #    parameterised in CL, V1, Q and V2 (Sect. 3.1).
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. ODE system. First-order elimination from the central compartment;
    #    a two-compartment model was preferred over one-compartment
    #    (dOFV = -13.6, p < 0.01) and gave no improvement over
    #    three-compartment (dOFV = 0, p > 0.05).
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 4. Observation. Doses in mg and volumes in L, so central/vc is mg/L,
    #    the unit the assay and every concentration in the paper use.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
