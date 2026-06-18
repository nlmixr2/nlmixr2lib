elDesoky_1997_theophylline_pediatric_asthma <- function() {
  description <- "One-compartment IV PK model for theophylline in 15 Egyptian pediatric patients (age 2-12 yr, weight 12-30 kg) treated for an acute asthma attack (elDesoky 1997). Aminophylline given as a 30-min loading infusion (6 mg/kg) followed by 12 hr of continuous maintenance infusion (1 mg/kg/hr); theophylline concentrations measured at 0.75, 7, and 13.25 hr. Parameter values taken from the Standard Calculations (SC) column of Table 2, which is independent of the Bayesian-prior population data and is treated by the authors as the reference (true) values."
  reference <- "El Desoky E, Ghazal MH, Mohamed MA, Klotz U. Disposition of intravenous theophylline in asthmatic children: Bayesian approach vs direct pharmacokinetic calculations. Japanese Journal of Pharmacology. 1997;75(1):13-20. doi:10.1254/jjp.75.13"
  vignette <- "elDesoky_1997_theophylline_pediatric_asthma"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear scaling on CL and Vc; reference weight 19.7 kg (Table 1 cohort mean, n=15, range 12.0-30.0 kg). The source paper reports CL and Vd per kg of body weight (Table 2), implying a linear (exponent = 1) weight effect, reproduced here as (WT/19.7)^1. Simulations at weights outside the studied 12-30 kg range are an extrapolation; the linearity assumption may break down in larger or smaller subjects.",
      source_name        = "WT"
    )
  )

  population <- list(
    species         = "human (Egyptian pediatric patients with acute bronchial asthma)",
    n_subjects      = 15L,
    n_studies       = 1L,
    age_range       = "2-12 years",
    age_mean        = "6.4 years (SD 3.4)",
    weight_range    = "12.0-30.0 kg",
    weight_mean     = "19.7 kg (SD 5.9)",
    height_range    = "80-126 cm",
    height_mean     = "104.6 cm (SD 15.1)",
    sex_female_pct  = 40,
    race_ethnicity  = "Egyptian (single-centre cohort at Assiut University Hospital, Egypt)",
    disease_state   = "Acute attack of allergic bronchial asthma (1-4 years' disease duration) with dyspnea, cough, wheezy hyperinflated chest; precipitated by an upper respiratory infection. Normal liver and renal function. Prior maintenance on oral albuterol (salbutamol) 0.1 mg/kg every 8 hr; theophylline-free for at least 1 week prior to admission.",
    dose_range      = "Aminophylline 6 mg/kg IV loading dose over 30 min, then 1 mg/kg/hr continuous IV maintenance for 12 hr (with a 15-min sampling pause at 6.25 hr). Aminophylline is 86 percent anhydrous theophylline, so the equivalent theophylline-base doses are LD 5.16 mg/kg and MD 0.86 mg/kg/hr.",
    regions         = "Egypt (Pediatric Department, Faculty of Medicine, Assiut University Hospital).",
    notes           = "Three blood samples per patient: C1 = 15 min after end of loading infusion (t = 0.75 hr); C2 = 6 hr after start of maintenance (t = 6.25 hr, after a 15-min infusion pause); C3 = end of second maintenance period (t = 13.25 hr). Theophylline measured by EMIT (Syva); assay precision CV 4.1-5.8 percent across 7.5-25 ug/mL. The paper compares 3 estimation methods on a single one-compartment IV model: SC (direct calculation from the 3 measured levels using equations a-d of Tozer), Bay 1 (Bayesian using only C1 and the Abbott Pharmacokinetics System with population prior CL = 0.09 L/hr/kg CV 30 percent and Vd = 0.5 L/kg CV 15 percent), Bay 3 (Bayesian using C1+C2+C3 and the same prior). SC is treated as the reference because it is independent of the population prior. Mean steady-state plasma level achieved ~12-13 ug/mL (within the therapeutic 8-20 ug/mL window)."
  )

  ini({
    # Structural parameters at reference body weight 19.7 kg (cohort mean, Table 1).
    # CL and Vd taken from the SC (Standard Calculations) column of Table 2 (mean +/- SD across 15 patients).
    # SC mean CL = 1.10 mL/min/kg = 0.066 L/hr/kg; at 19.7 kg, CL = 0.066 * 19.7 = 1.30 L/hr.
    # SC mean Vd = 0.50 L/kg; at 19.7 kg, Vd = 0.50 * 19.7 = 9.85 L.
    lcl <- log(1.30); label("Clearance at reference 19.7 kg (L/hr)")                     # elDesoky 1997 Table 2 SC column: mean CL = 1.10 mL/min/kg, equivalent to 0.066 L/hr/kg
    lvc <- log(9.85); label("Apparent central volume at reference 19.7 kg (L)")          # elDesoky 1997 Table 2 SC column: mean Vd = 0.50 L/kg

    # IIV approximated from inter-subject SD/mean of the SC parameter estimates across the 15
    # individual fits (Table 2 SC column, last row); these are descriptive variances of the
    # individual-fit estimates, not formal popPK omega estimates from a NONMEM-style fit. The
    # caveat is documented in the vignette Assumptions and deviations section.
    # CV(CL) = 0.20 / 1.10 = 0.1818;  omega^2 = log(1 + 0.1818^2) = 0.03257
    # CV(Vd) = 0.08 / 0.50 = 0.1600;  omega^2 = log(1 + 0.1600^2) = 0.02529
    etalcl ~ 0.03257  # elDesoky 1997 Table 2 SC column SD/mean: CV(CL) ~ 18 percent
    etalvc ~ 0.02529  # elDesoky 1997 Table 2 SC column SD/mean: CV(Vd) ~ 16 percent

    # Residual error not reported by the source paper. The EMIT theophylline assay has
    # within-run CV 4.1-5.8 percent and between-run CV 5.3 percent (Methods, page 14);
    # propSd = 0.10 is a conservative simulation-only default that envelopes the assay
    # precision plus model misspecification. Documented in the vignette Assumptions and
    # deviations section.
    propSd <- 0.10; label("Proportional residual error (fraction)")                     # not reported by elDesoky 1997; conservative simulation-only default (assay CV 4-6 percent per Methods page 14)
  })

  model({
    cl <- exp(lcl + etalcl) * (WT / 19.7)
    vc <- exp(lvc + etalvc) * (WT / 19.7)

    kel <- cl / vc

    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
