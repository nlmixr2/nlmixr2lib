Yu_2022_colistinSulfate <- function() {
  description <- paste(
    "One-compartment population PK model for intravenous colistin sulfate in",
    "critically ill adults with carbapenem-resistant organism infections",
    "(Yu 2022; n = 42 Chinese patients, 112 sparse steady-state therapeutic-",
    "drug-monitoring plasma concentrations spanning 0.28-6.20 mg/L). Linear",
    "elimination with intravenous-infusion input. Cockcroft-Gault creatinine",
    "clearance enters clearance as an additive linear term,",
    "CL = 0.994 + 0.525 * CrCL / 66.47 L/h, with exponential inter-individual",
    "variability on CL only and a proportional residual error. Colistin",
    "sulfate is the active drug and must not be confused with colistimethate",
    "sodium (CMS), the inactive prodrug modelled in Plachouras 2009,",
    "Mohamed 2012, Jacobs 2016 and Karaiskos 2015. Doses are in mg; the",
    "paper reports doses only in international units and its stated",
    "conversion is contradicted by its own simulations (see the vignette)."
  )
  reference <- paste(
    "Yu XB, Zhang XS, Wang YX, Wang YZ, Zhou HM, Xu FM, Yu JH, Zhang LW,",
    "Dai Y, Zhou ZY, Zhang CH, Lin GY, Pan JY (2022).",
    "Population pharmacokinetics of colistin sulfate in critically ill",
    "patients: exposure and clinical efficacy.",
    "Front Pharmacol 13:915958. doi:10.3389/fphar.2022.915958.",
    "See also the published commentary: Chen H, Li P (2022). Commentary:",
    "Population pharmacokinetics of colistin sulfate in critically ill",
    "patients: exposure and clinical efficacy. Front Pharmacol 13:992085.",
    "doi:10.3389/fphar.2022.992085.",
    sep = " "
  )
  vignette <- "Yu_2022_colistinSulfate"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    # Methods "Colistin Sulfate Administration and Sample Collection": plasma
    # was separated by centrifugation and colistin was quantified by a
    # validated HPLC-MS/MS assay (calibration range 0.1-20 ug/mL). The
    # analyte is colistin itself; colistin sulfate is not a prodrug.
    central = list(analyte = "colistin sulfate", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Creatinine clearance estimated with the Cockcroft-Gault equation,",
        "reported as RAW mL/min and NOT normalised to 1.73 m^2 body surface",
        "area."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table 1 footnote b: 'Creatinine clearance calculated using the",
        "Cockcroft-Gault equation'; cohort 79.54 +/- 53.99 mL/min (mean +/-",
        "SD). Additive linear effect on CL per Yu 2022 Eq. 1,",
        "CL (L/h) = 0.994 + 0.525 * CrCL / 66.47, so 0.994 L/h is the",
        "clearance extrapolated to zero creatinine clearance and 0.525 L/h is",
        "the increment per 66.47 mL/min. The divisor 66.47 appears only in",
        "Eq. 1 and in Table 4; it is not the tabulated cohort mean (79.54) and",
        "the paper does not report a median, so it is most plausibly the",
        "unreported cohort median. Exponential IIV multiplies the whole",
        "bracket. Concentrations drawn during renal replacement therapy or",
        "ECMO were excluded from the modelling dataset, so the model carries",
        "no information about extracorporeal clearance. The paper's Monte",
        "Carlo simulations span 10-120 mL/min."
      ),
      source_name = "CrCL"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 42L,
    n_studies = 1L,
    n_observations = 112L,
    age_range = "adults >= 18 years (inclusion criterion); range not reported",
    age_mean = "67.90 +/- 13.74 years",
    weight_mean = "63.47 +/- 9.64 kg",
    sex_female_pct = 100 * 5 / 42,
    race_ethnicity = c(Asian = 100),
    disease_state = paste(
      "Critically ill adults treated with intravenous colistin sulfate for",
      "at least 3 days for confirmed carbapenem-resistant Gram-negative",
      "infections. Infection sites: respiratory tract 69.05%, blood 26.19%,",
      "abdomen 23.81%, intracranial 9.52%; 33.33% had multiple sites.",
      "Isolates: Acinetobacter baumannii 73.81%, Pseudomonas aeruginosa",
      "16.67%, Klebsiella pneumoniae 7.14%, Enterobacter cloacae. 76.19%",
      "were mechanically ventilated and 35.71% received vasoactive agents;",
      "APACHE II 17 [14, 26]. Every patient received a concomitant",
      "antibacterial (carbapenems 52.38%, cefoperazone-sulbactam 47.62%,",
      "tigecycline 16.67%)."
    ),
    dose_range = paste(
      "Intravenous drip of colistin sulfate (Shanghai New Asia",
      "Pharmaceuticals; 500,000 IU per vial). Label regimen: 1 million IU",
      "loading dose then 1.5 million IU/day in 2-3 divided doses. Table 1",
      "daily dose 150 [150, 200] x 10^4 IU; treatment 11.63 +/- 5.96 days.",
      "22 patients (52.38%) also received inhaled colistin sulfate (250,000",
      "IU q12h) and 4 received intraventricular/intrathecal colistin",
      "sulfate (50,000 IU q24h); neither route is represented in the model.",
      "IMPORTANT: doses are stated only in IU; this model takes mg (see",
      "`notes`)."
    ),
    sampling = paste(
      "Sparse therapeutic drug monitoring at steady state (after at least",
      "six doses), with dosing and sampling clock times indexed from the",
      "medical records. 112 concentrations, 0.28-6.20 mg/L (Figure 1)."
    ),
    renal_function = paste(
      "Cockcroft-Gault creatinine clearance 79.54 +/- 53.99 mL/min; serum",
      "creatinine 116.64 +/- 105.49 umol/L (Table 1). Samples taken during",
      "renal replacement therapy or ECMO were excluded from the PK analysis."
    ),
    regions = "People's Republic of China (single centre; First Affiliated Hospital of Wenzhou Medical University).",
    notes = paste(
      "Baseline demographics from Yu 2022 Table 1. Retrospective cohort,",
      "January 2020 to December 2021. Fitted in NONMEM 7.4 with Pirana",
      "2.9.7; evaluated with goodness-of-fit plots, a 1000-sample",
      "nonparametric bootstrap (Supplementary Table S2) and a prediction-",
      "and variability-corrected VPC (Supplementary Figure S2).",
      "",
      "DOSE UNITS. The Introduction states '1 mg of pure colistin base =",
      "17,000 IU of colistin', i.e. 58.8 mg per million IU. The published",
      "commentary (Chen and Li 2022) points out that 17,000 IU/mg is only",
      "the Chinese Pharmacopoeia lower potency limit for colistin SULFATE,",
      "and that the Chinese consensus conversion is about 22,300 IU/mg",
      "(44.8 mg per million IU). Neither reproduces the paper's own",
      "model-based output: fitting the paper's Figure 3 probability-of-",
      "target-attainment curves with this model recovers about 32.5 mg per",
      "million IU (close to the colistin-base-activity convention of",
      "30,000 IU/mg, 33.3 mg per million IU) with the IIV reported in",
      "Table 3, and the paper's MAP AUC0-24,ss of 39.39 +/- 14.47 mg*h/L",
      "(Table 2) also points below 44.8. The conversion is a property of",
      "the product and dataset, not a model parameter, so it is not",
      "encoded; the vignette derives it and shows the evidence.",
      "",
      "The unbound fraction of 0.5 and the fAUC/MIC >= 15 target used for",
      "the paper's probability-of-target-attainment simulations are",
      "literature assumptions (Cheah 2015), not fitted parameters, and are",
      "not encoded."
    )
  )

  ini({
    # Structural parameters. Yu 2022 Eq. 1: CL (L/h) = 0.994 + 0.525 x
    # CrCL / 66.47; Eq. 2: V (L) = 20.7. Table 3 'Population pharmacokinetic
    # parameter estimates from the final model'.
    lcl <- log(0.994); label("Intercept of additive linear CL ~ CrCL (CL at CrCL = 0, L/h)") # Table 3 'TVCL (L/h)' 0.994 (RSE 16%); Eq. 1 intercept
    lvc <- log(20.7); label("Volume of distribution (V, L)") # Table 3 'TVV (L)' 20.7 (RSE 10%); Eq. 2

    # Covariate effect: additive linear CrCL term on CL, divided by 66.47
    # mL/min (Eq. 1 and Table 4 'CL = 0.994 + 0.525 x CrCL/ 66.47 L/h').
    e_crcl_cl <- 0.525; label("Renal CL slope per (CrCL / 66.47 mL/min) (L/h)") # Table 3 'CrCL on CL (theta 1)' 0.525 (RSE 22%); Eq. 1

    # IIV. Methods: 'Between-subject variability (BSV) was assessed using an
    # exponential function.' Table 3 reports BSV on CL only, as a %CV, so
    # omega^2 = log(CV^2 + 1) = log(0.304^2 + 1) = 0.08835. (BSV on V was
    # estimated in the base model but is not in the final-model table.)
    etalcl ~ 0.08835 # Table 3 'BSV_CL [%CV]' 30.40% (RSE 15%, shrinkage 9%); omega^2 = log(0.304^2 + 1)

    # Residual error. Results: 'The proportional error model was selected to
    # evaluate the residual variability.'
    propSd <- 0.251; label("Proportional residual error (fraction)") # Table 3 'Proportional error [%CV]' 25.10% (RSE 16%, shrinkage 14%)
  })

  model({
    # Individual parameters. Additive linear CrCL term (Eq. 1) with
    # exponential IIV on the whole clearance; V has no covariate and no IIV
    # in the final model (Eq. 2).
    cl <- (exp(lcl) + e_crcl_cl * (CRCL / 66.47)) * exp(etalcl)
    vc <- exp(lvc)

    kel <- cl / vc

    # Intravenous drip: dose records target `central` with a rate or
    # duration; there is no depot.
    d/dt(central) <- -kel * central

    # Plasma colistin concentration in mg/L (= ug/mL).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
