Mehta_2025_burosumab_fracture <- function() {
  description <- paste(
    "Poisson exposure-response model for the count of active fractures plus",
    "pseudofractures on skeletal survey in adults with X-linked",
    "hypophosphatemia treated with burosumab or placebo (UX023-CL303). The",
    "Poisson mean is a baseline count with log-normal IIV, multiplied by a",
    "log-linear time effect, an additional log-linear time effect in",
    "patients randomised to burosumab, and a log-linear effect of the",
    "model-predicted average serum phosphorus per 4-week dosing cycle",
    "centred at 2.5 mg/dL. The serum phosphorus input comes from the",
    "companion PK/PD model Mehta_2025_burosumab."
  )
  reference <- paste(
    "Mehta K, Storopoli J, Ramwani N, Quattrocchi E, Gobburu J, Weber T,",
    "Hruska MW, Marsteller D. Pharmacodynamic Exposure-Response Analysis of",
    "Fracture Count Data Following Treatment with Burosumab in Patients with",
    "XLH. J Clin Pharmacol. 2025;65(2):253. doi:10.1002/jcph.6140"
  )
  vignette <- "Mehta_2025_burosumab"
  units <- list(
    time = "week",
    dosing = "mg",
    concentration = "count (the model output is the mean number of fractures plus pseudofractures, not a drug concentration)"
  )

  covariateData <- list(
    CAV_PHOS = list(
      description = "Model-predicted average serum phosphorus over the 4-week dosing cycle matched to the fracture assessment (time-varying)",
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The paper's 'Serum phosphate Cavg' (Table 4), computed per dosing",
        "cycle from the individual predictions of the PK/PD model",
        "(Mehta_2025_burosumab) and merged time-matched with the counts",
        "(Methods 'Population Pharmacokinetic-Pharmacodynamic Model'). The",
        "paper works in mg/dL; the model converts inline with x 3.0975 and",
        "centres at the lower limit of normal, 2.5 mg/dL (Abstract and",
        "Discussion). Without burosumab the cycle average equals the",
        "patient's baseline serum phosphorus."
      ),
      source_name = "Cavg"
    ),
    TRT_BUROSUMAB = list(
      description = "Randomised treatment arm: 1 = burosumab from week 0, 0 = placebo for 24 weeks then open-label burosumab",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (randomised to placebo)",
      notes = paste(
        "The paper's treatment indicator, 'introduced ... to allow",
        "differentiation between participants who received the drug under",
        "the study and those who received placebo for 24 weeks' (Results",
        "'Pharmacodynamic Exposure-Fracture Response Modeling'). Coded here as the",
        "static randomised arm: it stays 0 for placebo-randomised patients",
        "after their week-24 crossover, whose burosumab benefit then enters",
        "through CAV_PHOS only. The paper does not state the post-crossover",
        "coding."
      ),
      source_name = "treatment/placebo indicator"
    )
  )

  covariatesDataExcluded <- list(
    PHOS = list(
      description = "Baseline serum phosphorus",
      units = "mmol/L",
      type = "continuous",
      notes = "Supplement S2 run 4 (baseline phosphate levels) lowered the objective function by only 0.13 and was not carried forward."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Not correlated with the baseline fracture count in the exploratory analysis (Results, Supplement S3); not tested in the model."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Not correlated with the baseline fracture count in the exploratory analysis (Results, Supplement S3); not tested in the model."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Not correlated with the baseline fracture count in the exploratory analysis (Results, Supplement S3); not tested in the model."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 70,
    n_studies = 1,
    age_range = "18-65 years (UX023-CL303 eligibility)",
    weight_range = "not reported",
    disease_state = "X-linked hypophosphatemia (XLH) with a confirmed PHEX mutation and at least one active fracture or pseudofracture at baseline",
    dose_range = "burosumab 1.0 mg/kg SC (rounded to the nearest 10 mg, maximum 90 mg) every 4 weeks; placebo arm crossed over to burosumab at week 24",
    regions = "United States, European Union, Japan and Korea",
    notes = paste(
      "UX023-CL303 (NCT02526160), a randomised double-blind placebo-controlled",
      "phase 3 study with open-label extension (Table 1): 134 adults",
      "randomised 1:1, of whom 70 (32 burosumab, 38 placebo) had at least one",
      "fracture or pseudofracture at baseline and entered this analysis.",
      "Counts were assessed by skeletal survey at weeks 0, 12, 24, 36, 48,",
      "84, 96, 108 and 116 (Table 3). Fitted in Pumas 2.3.1."
    )
  )

  ini({
    # ---- Poisson mean, Mehta 2025 Table 4 ----
    lbase <- log(1.87)
    label("Baseline mean count of fractures plus pseudofractures (count)") # Table 4: tvbase = 1.87 (RSE 11.74%)
    e_time_base <- -0.0058
    label("Log-scale time effect on the mean count (1/week)") # Table 4: tv-alpha 'Time effect (week)' = -0.0058 (RSE 16.57%)
    e_trt_time_base <- -0.0058
    label("Additional log-scale time effect in patients randomised to burosumab (1/week)") # Table 4: tv-beta 'Treatment effect over each unit of time (week)' = -0.0058 (RSE 16.57%)
    e_cavphos_base <- -0.26
    label("Log-scale effect of cycle-average serum phosphorus on the mean count (per mg/dL)") # Table 4: tvdPhos 'Serum phosphate Cavg effect (mg/dL)' = -0.26 (RSE 40.21%)

    # ---- IIV ----
    etalbase ~ log(1 + 0.4515^2) # Table 4: IIV-base = 45.15%
  })

  model({
    # ---- Serum phosphorus input in the paper's units ----
    # CAV_PHOS is carried in mmol/L; 1 mmol/L = 3.0975 mg/dL elemental
    # phosphorus. The effect is centred at 2.5 mg/dL, the lower limit of
    # normal (Abstract: 'by 23% with each unit increase (1.0 mg/dL) in
    # average serum phosphate from lower limit of normal (2.5 mg/dL)').
    cav_phos_mgdl <- CAV_PHOS * 3.0975

    # ---- Poisson mean (Equation 2: lambda = base * e^time * e^covariate) ----
    # t is weeks since the first randomised dose (week 0 = baseline survey).
    base <- exp(lbase + etalbase)
    lambda <- base *
      exp(
        e_time_base * t +
          e_trt_time_base * TRT_BUROSUMAB * t +
          e_cavphos_base * (cav_phos_mgdl - 2.5)
      )

    # ---- Observation (Equation 1: Y ~ Poisson(lambda)) ----
    fracture_count <- lambda
    fracture_count ~ pois(lambda)
  })
}
