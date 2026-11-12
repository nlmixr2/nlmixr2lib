Wu_2025_idecabtageneVicleucel_auc <- function() {
  description <- paste0(
    "Multivariable linear regression of the natural log of the ",
    "noncompartmental day-0-28 transgene-level AUC (AUC0-28D, ",
    "day*copies/ug, droplet digital PCR) of the BCMA-directed CAR T-cell ",
    "therapy idecabtagene vicleucel (ide-cel) on the actual infused dose ",
    "and baseline covariates, in adults with triple-class-exposed ",
    "relapsed/refractory multiple myeloma (Wu 2025; ide-cel arm of the ",
    "phase 3 KarMMa-3 trial, N = 220). The dose slope is 0.0022 per ",
    "10^6 CAR+ T cells (95% CI -0.001 to 0.0055), i.e. about +0.4 log ",
    "units across the 300-510 x 10^6 range: the dose-exposure relationship ",
    "is essentially flat. Exposure falls with body weight, extramedullary ",
    "plasmacytoma and the Formulation B product and rises with baseline ",
    "IL-15 and with post-infusion ADA positivity. Continuous covariates are ",
    "centred on the reference values the paper reports (body weight 75 kg; ",
    "medians of ferritin, soluble BCMA and IL-15), the dose is not centred. ",
    "No ODE and no cellular kinetic layer; the model returns the predicted ",
    "log(AUC0-28D) (lauc0_28) and AUC0-28D (auc0_28) as derived variables. ",
    "The residual standard deviation of the regression is not reported, ",
    "so no endpoint is declared."
  )
  reference <- paste(
    "Wu F, Zheng X, Burnett J, Masilamani M, Zhang W, Zhong X, Caia A,",
    "Cook M, Piasecki J, Kondic A, Lamba M, Zhou J.",
    "Characterizing cellular expansion of idecabtagene vicleucel and",
    "association with clinical efficacy and safety in patients with",
    "triple-class-exposed relapsed/refractory multiple myeloma.",
    "J Clin Pharmacol. 2025;65(11):1568-1576.",
    "doi:10.1002/jcph.70075.",
    "Parameter estimates in Supporting Information Table S4 (full = final",
    "model) and covariate reference values in Table S5.",
    sep = " "
  )
  vignette <- "Wu_2025_idecabtageneVicleucel"
  units <- list(
    time = "n/a (static regression of a per-subject exposure metric; no time dimension)",
    dosing = "n/a (no dose events; the actual infused dose enters as the DOSE_IDECEL_MCELLS covariate column)",
    concentration = "lauc0_28 (natural log of AUC0-28D in day*copies/ug); auc0_28 (AUC0-28D, day*copies/ug)"
  )

  covariateData <- list(
    DOSE_IDECEL_MCELLS = list(
      description = "Actual infused dose of idecabtagene vicleucel, in millions (10^6) of CAR+ T cells.",
      units = "10^6 CAR+ T cells",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters linearly and UNCENTRED: the Table S4 intercept is the predicted log(AUC0-28D) of the reference patient at zero dose (Figure 1: the reference-patient line is 13.1186 + 0.0022 * dose, about 13.5 at 170 and 14.3 at 530 x 10^6 cells). KarMMa-3 target range 150-450, up to 540 x 10^6 permitted.",
      source_name = "Actual dose (x10^6 cells)"
    ),
    FERRITIN_BL = list(
      description = "Baseline serum ferritin.",
      units = "ug/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Natural log, centred on the analysis-set median 163 ug/L (Table S5 note; Table S3 prints the N = 220 median as 162). Missing values (3 subjects) were imputed to the population median.",
      source_name = "Log(baseline ferritin) (ug/L)"
    ),
    WT = list(
      description = "Baseline body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear, centred on 75 kg, the reference value Table S5 states for body weight (the analysis-set median is 83 kg, Table S3). Reproduces Table S5 exactly: exp(-0.0111 * (54.97 - 75)) - 1 = +24.9% at the 5th percentile.",
      source_name = "Baseline body weight (kg)"
    ),
    SBCMA = list(
      description = "Baseline serum soluble B-cell maturation antigen, measured on the day of the ide-cel infusion.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Natural log, centred on the analysis-set median 144.5 ng/mL (Table S5 note).",
      source_name = "Log(baseline soluble BCMA) (ng/mL)"
    ),
    FORM_IDECEL_B = list(
      description = "1 = ide-cel manufactured as Formulation B, 0 = Formulation A (reference).",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Product attribute; the paper states only that the two formulations differ in the transduction media. 72% of the analysis set received Formulation A, 26% Formulation B (2% missing, imputed to the most frequent category).",
      source_name = "Product Attributes (Formulation B vs A)"
    ),
    LINE_4L = list(
      description = "1 = exactly three prior multiple-myeloma regimens (ide-cel is the fourth line), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "With LINE_5L decomposes the three-level 'number of prior MM regimens' (2 / 3 / 4). Reference level = 2 prior regimens (both indicators 0).",
      source_name = "No. of prior MM regimens (3 vs 2)"
    ),
    LINE_5L = list(
      description = "1 = exactly four prior multiple-myeloma regimens (ide-cel is the fifth line), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "KarMMa-3 enrolled patients with two to four prior regimens, so four is the top level.",
      source_name = "No. of prior MM regimens (4 vs 2)"
    ),
    IL15 = list(
      description = "Baseline serum interleukin-15, measured on the day of the ide-cel infusion.",
      units = "pg/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Natural log, centred on the analysis-set median 43 pg/mL (Table S5 note). IL-15 is a homeostatic T-cell cytokine raised by the lymphodepleting chemotherapy.",
      source_name = "Log(baseline IL-15) (pg/mL)"
    ),
    DIS_EMD = list(
      description = "1 = extramedullary plasmacytoma present at baseline, 0 = absent.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "23.6% of the analysis set (Table S3).",
      source_name = "Baseline extramedullary plasmacytoma (Yes vs No)"
    ),
    ADA_POS = list(
      description = "1 = anti-CAR antidrug antibody positive at at least one post-infusion visit, 0 = negative at every visit.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "A POST-BASELINE covariate ('ever positive' after the infusion; 61.4% of the analysis set). Most positivity appeared from about month 3, after the day-0-28 exposure window, so the authors caution it may mark longer survival and follow-up rather than cause higher exposure.",
      source_name = "Post-infusion ADA positive (Yes vs No)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Baseline age",
      units = "years",
      type = "continuous",
      notes = "Screened; not found to influence AUC0-28D (Results)."
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not found to influence AUC0-28D (Results)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 220L,
    n_studies = 1L,
    n_observations = "220 per-subject noncompartmental AUC0-28D values",
    age_range = "30-81 years (median 63; 43% aged 65 or older)",
    weight_range = "41.1-144 kg (median 83)",
    sex_female_pct = 37.7,
    race_ethnicity = c(White = 70.0, Black = 7.3, Asian = 2.3, Other = 0.9, Unknown = 19.5),
    disease_state = "Triple-class-exposed relapsed/refractory multiple myeloma after two to four prior regimens (including daratumumab, an immunomodulatory agent and a proteasome inhibitor) with progression within 60 days of the last regimen; 61.8% triple-class refractory",
    dose_range = "Single IV infusion; target 150-450 x 10^6 CAR+ T cells (up to 540 x 10^6 permitted); the approved range was extended to 300-510 x 10^6 on these data",
    regions = "International (KarMMa-3, NCT03651128)",
    notes = "Ide-cel arm of KarMMa-3, CK data cutoff 1 March 2022; 225 infused, 220 with evaluable exposure parameters. Transgene levels by droplet digital PCR (LLOQ 40 copies/ug); NCA in Phoenix WinNonlin 8.3.5. The paper notes that the KarMMa (phase 2) and KarMMa-3 exposures used different bioanalytical methods and were not compared directly. Baseline covariates: Table S3."
  )

  ini({
    # Wu 2025 Table S4, 'Parameter Estimates of Dose-log(AUC0-28D)
    # Full/Final Model' (natural log; the final model equals the full
    # model, Results). Continuous covariates are centred on the Table S5
    # reference values (ferritin 163 ug/L, soluble BCMA 144.5 ng/mL,
    # IL-15 43 pg/mL, body weight 75 kg); with that centring every Table S5
    # percent change is reproduced (vignette). The dose is not centred:
    # Figure 1's reference-patient line runs 13.1186 + 0.0022 * dose.
    lauc_ref <- 13.1186; label("Natural log of AUC0-28D (day*copies/ug) for the reference patient at zero dose (unitless)") # Table S4 Intercept 13.1186 (SE 0.7565)
    e_dose_idecel_lauc <- 0.0022; label("Slope of log(AUC0-28D) per 10^6 CAR+ T cells of actual dose (1/(10^6 cells))") # Table S4 'Actual dose (x10^6 cells)' 0.0022 (SE 0.0016)
    e_ferritin_bl_lauc <- 0.0458; label("Slope of log(AUC0-28D) per unit log(ferritin / 163 ug/L) (unitless)") # Table S4 'Log(baseline ferritin)' 0.0458 (SE 0.0793)
    e_wt_lauc <- -0.0111; label("Slope of log(AUC0-28D) per kg of body weight above 75 kg (1/kg)") # Table S4 'Baseline body weight (kg)' -0.0111 (SE 0.0049)
    e_sbcma_lauc <- 0.1615; label("Slope of log(AUC0-28D) per unit log(soluble BCMA / 144.5 ng/mL) (unitless)") # Table S4 'Log(baseline soluble BCMA)' 0.1615 (SE 0.0854)
    e_form_idecel_b_lauc <- -0.3939; label("Shift in log(AUC0-28D) for Formulation B vs A (unitless)") # Table S4 'Product Attributes (Formulation B vs A)' -0.3939 (SE 0.1985)
    e_line_4l_lauc <- -0.1576; label("Shift in log(AUC0-28D) for 3 vs 2 prior MM regimens (unitless)") # Table S4 'No. of prior MM regimens (3 vs 2)' -0.1576 (SE 0.2093)
    e_line_5l_lauc <- -0.2029; label("Shift in log(AUC0-28D) for 4 vs 2 prior MM regimens (unitless)") # Table S4 'No. of prior MM regimens (4 vs 2)' -0.2029 (SE 0.2203)
    e_il15_lauc <- 0.5157; label("Slope of log(AUC0-28D) per unit log(IL-15 / 43 pg/mL) (unitless)") # Table S4 'Log(baseline IL-15)' 0.5157 (SE 0.1596)
    e_emd_lauc <- -0.5036; label("Shift in log(AUC0-28D) for baseline extramedullary plasmacytoma (unitless)") # Table S4 'Baseline extramedullary plasmacytoma (Yes vs No)' -0.5036 (SE 0.2089)
    e_ada_pos_lauc <- 0.4621; label("Shift in log(AUC0-28D) for post-infusion ADA positivity (unitless)") # Table S4 'Post-infusion ADA positive (Yes vs No)' 0.4621 (SE 0.1847)
  })

  model({
    lauc0_28 <- lauc_ref +
      e_dose_idecel_lauc * DOSE_IDECEL_MCELLS +
      e_ferritin_bl_lauc * log(FERRITIN_BL / 163) +
      e_wt_lauc * (WT - 75) +
      e_sbcma_lauc * log(SBCMA / 144.5) +
      e_form_idecel_b_lauc * FORM_IDECEL_B +
      e_line_4l_lauc * LINE_4L +
      e_line_5l_lauc * LINE_5L +
      e_il15_lauc * log(IL15 / 43) +
      e_emd_lauc * DIS_EMD +
      e_ada_pos_lauc * ADA_POS

    auc0_28 <- exp(lauc0_28)
  })
}
