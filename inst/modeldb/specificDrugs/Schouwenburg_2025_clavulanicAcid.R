Schouwenburg_2025_clavulanicAcid <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and",
    "estimated oral bioavailability for clavulanic acid in preterm and term",
    "neonates and infants up to 1 year of age, pooled from four datasets",
    "(RAIN, Staph Trio/Durham, De Cock, Dhont) covering both oral and",
    "intravenous amoxicillin/clavulanic acid and ticarcillin/clavulanic acid.",
    "Clearance carries a priori allometric weight scaling (fixed exponent",
    "0.75, reference 3.9 kg) plus one estimated postnatal-age power term",
    "(reference 55.5 days); central volume carries fixed linear weight",
    "scaling (exponent 1, reference 3.9 kg). Oral bioavailability is 24.4%,",
    "estimated from neonates up to 10 days of age, and absorption is rapid",
    "(Ka 0.781 1/h). Interindividual variability is on clearance only;",
    "residual error is proportional (log-transform-both-sides, 77.8%)",
    "(Schouwenburg 2025)."
  )
  reference <- paste(
    "Schouwenburg S, Keij FM, Tramper-Stranders GA, Kornelisse RF,",
    "Reiss IKM, De Cock PAJG, Dhont E, Watt KM, Muller AE, Flint RB,",
    "Koch BCP, Allegaert K, Preijers T. A Pooled Population",
    "Pharmacokinetic Study of Oral and Intravenous Administration of",
    "Clavulanic Acid in Neonates and Infants: Targeting Effective",
    "Beta-Lactamase Inhibition. Clin Pharmacol Ther. 2025;117(1):193-202.",
    "doi:10.1002/cpt.3423.",
    "The companion amoxicillin model this paper simulates alongside",
    "clavulanic acid (its reference 18) is available as",
    "modellib('Keij_2023_amoxicillin')."
  )
  vignette <- "Schouwenburg_2025_clavulanicAcid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Both entries confirmed against Schouwenburg 2025:
  # the model is a one-compartment model with a first-order absorption
  # depot, and all reported concentrations are plasma clavulanic acid
  # (Results: "403 clavulanic acid plasma concentrations").
  compartmentData <- list(
    depot   = list(analyte = "clavulanic acid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "clavulanic acid", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Current (not birth) body weight of the neonate or infant.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column BW ('current bodyweight'). Enters both clearance and",
        "central volume as an allometric power term normalized to 3.9 kg,",
        "the pooled dataset median (Schouwenburg 2025 Table 2 equations and",
        "Table 2 footnote: 'Current bodyweight and postnatal age are scaled",
        "to the dataset median (3.9 kg and 55.5 days)'). Both exponents were",
        "FIXED a priori, not estimated: Methods ('Covariate relationship",
        "analysis') states 'a priori allometric scaling with a fixed exponent",
        "(i.e., 0.75 on clearances, 1 on distribution volume)'. Unlike the",
        "companion Keij 2023 amoxicillin model, which anchors allometry at",
        "70 kg, this model anchors at the 3.9 kg cohort median, so",
        "exp(lcl) = 0.675 L/h and exp(lvc) = 1.88 L are directly",
        "interpretable as neonatal typical values rather than extrapolated",
        "adult anchors. Time-varying in principle; the source datasets",
        "recorded weight during therapy. Observed range 0.55-9.0 kg",
        "(Table 1), so simulations outside that range are extrapolation."
      ),
      source_name        = "BW"
    ),
    PNA = list(
      description        = "Postnatal (chronological) age since birth.",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Schouwenburg 2025 reports postnatal age in DAYS (the Table 2",
        "clearance equation uses (PNA / 55.5) with PNA in days; the",
        "55.5-day reference is the pooled dataset median per the Table 2",
        "footnote 'Current bodyweight and postnatal age are scaled to the",
        "dataset median (3.9 kg and 55.5 days)'). The canonical nlmixr2lib",
        "PNA column carries MONTHS, so model() reparameterises the reference",
        "as 55.5 / 30.4375 = 1.8234 months; numerator and denominator carry",
        "the same units factor, so the ratio and hence the exponent are",
        "unchanged. This follows the Zhao_2018_omeprazole and",
        "Keij_2023_amoxicillin precedents. Users supply PNA in months.",
        "NOTE 1: Table 1 reports a pooled median PNA of 54.5 days, one day",
        "less than the 55.5-day centering constant printed in the Table 2",
        "equation and footnote. The equation is authoritative and 55.5 is",
        "used here (the same Table-1-vs-equation discrepancy occurs in",
        "Keij 2023 for gestational age).",
        "NOTE 2: the Results text describes the retained PNA relationship",
        "as 'an exponential function', but the Table 2 equation prints a",
        "POWER form, (PNA / 55.5)^theta_PNA, and the Methods state that",
        "continuous covariates were 'centered on the median and were",
        "evaluated as exponential or power relationships'. The printed",
        "equation is authoritative, so the power form is encoded.",
        "NOTE 3: the power form is singular at PNA = 0 -- clearance goes to",
        "zero on the day of birth, and concentration diverges. The authors'",
        "own simulations nevertheless include PNA = 0 patients (Figure 3",
        "'youngest (PNA 0 days, weight 3.2 kg)' and Figure 4 panels A-B),",
        "so the published model cannot be evaluated at exactly PNA = 0 as",
        "written. Use PNA >= 1 day (0.0329 months) when simulating.",
        "Observed range was PNA 0-365 days (Table 1); the ORAL sub-cohort",
        "that supports the bioavailability estimate spans only PNA 0-12",
        "days (RAIN, median 2.9 days), and the authors capped their own oral",
        "simulations at PNA 10 days 'to prevent extrapolation of oral data'."
      ),
      source_name        = "PNA"
    )
  )

  # Covariates that Schouwenburg 2025 screened on clearance during the
  # stepwise covariate analysis but did NOT retain in the final model
  # (Methods: "The following covariates were tested on CLclav: sex, PNA,
  # GA, PMA, and study center"; Conclusions: "No additional covariate
  # relationships with clearance were found"). Documentation only -- no
  # published point estimates exist for these, so they carry no effect in
  # model().
  covariatesDataExcluded <- list(
    SEXF = list(
      description        = "Female sex indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Tested on clearance as a categorical covariate (Methods:",
        "'Categorical variables were modeled using a proportional model')",
        "via stepwise forward inclusion (P < 0.05) and backward elimination",
        "(P < 0.01), and not retained. Pooled cohort was 32.2% female",
        "(Table 1). No point estimate is published, so no effect is",
        "encoded."
      ),
      source_name        = "SEX"
    ),
    GA = list(
      description        = "Gestational age at birth. Time-fixed per subject.",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Tested on clearance and not retained; only postnatal age survived",
        "(Results: 'PNA was found to best describe maturation of CLclav",
        "(P < 0.05; -14.0 difference in objective function value (dOFV))').",
        "Pooled median GA 37.4 weeks, range 23.0-41.7 (Table 1). Note that",
        "GA was NOT recorded in the De Cock et al. dataset and was imputed:",
        "Methods ('Covariate relationship analysis') states 'The De Cock",
        "et al. study did not register GA. However, as no ex-premature",
        "children were included the GA was set at 40 weeks.' No point",
        "estimate is published, so no effect is encoded. Contrast with the",
        "companion Keij_2023_amoxicillin model, where GA WAS retained on",
        "clearance with a steep exponent of 2.37."
      ),
      source_name        = "GA"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age plus postnatal age).",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Tested on clearance as a continuous covariate and not retained;",
        "postnatal age entered alone described maturation of clavulanic",
        "acid clearance better (Results: 'PNA was found to best describe",
        "maturation of CLclav'). Pooled median PMA 45.1 weeks, range",
        "24.9-92.1 (Table 1). No point estimate is published, so no effect",
        "is encoded."
      ),
      source_name        = "PMA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 89,
    n_studies      = 4,
    n_observations = 403,
    age_range      = "postnatal age 0-365 days (pooled median 54.5 days)",
    ga_range       = "gestational age 23.0-41.7 weeks (pooled median 37.4 weeks)",
    pma_range      = "postmenstrual age 24.9-92.1 weeks (pooled median 45.1 weeks)",
    weight_range   = "0.55-9.0 kg",
    weight_median  = "3.9 kg",
    sex_female_pct = 32.2,
    disease_state  = paste(
      "preterm and term neonates and infants treated for probable,",
      "suspected or confirmed bacterial infection; includes critically ill",
      "children in paediatric and cardiac intensive care"
    ),
    dose_range     = paste(
      "clavulanic acid 1.4-80 mg per administration (median 22 mg);",
      "2.3-10.8 mg/kg/day (median 5.2); median 5 mg/kg/dose (range 1.0-12.8)"
    ),
    routes         = paste(
      "oral (42 subjects, 82 concentrations) and intravenous (47 subjects,",
      "321 concentrations); 6 of the 42 oral subjects switched from",
      "intravenous to oral during sampling"
    ),
    regions        = "The Netherlands (RAIN), Belgium (Ghent: De Cock and Dhont), USA (Durham/Staph Trio)",
    notes          = paste(
      "Pooled analysis of four datasets (Schouwenburg 2025 Table 1):",
      "(1) RAIN (Reduction of intravenous Antibiotics In Neonates,",
      "The Netherlands), n = 42, oral and intravenous, PMA >= 35 weeks,",
      "PNA 0-28 days, weight >= 2 kg, amoxicillin/clavulanic acid 4:1;",
      "(2) Staph Trio / Durham (USA), n = 15, intravenous, PNA < 91 days,",
      "GA < 30 weeks, normal renal function, ticarcillin/clavulanic acid",
      "30:1; (3) De Cock et al. (Ghent University Hospital paediatric ICU,",
      "Belgium), n = 6, intravenous, amoxicillin/clavulanic acid 5:1,",
      "restricted here to patients younger than 1 year; (4) Dhont et al.",
      "(Ghent cardiac ICU, Belgium, not yet published at time of writing),",
      "n = 26, intravenous, amoxicillin/clavulanic acid 5:1 and 10:1,",
      "critically ill children after cardiac surgery.",
      "Race/ethnicity was not reported. Modelling was in NONMEM 7.4.3 with",
      "PsN 4.2.0, Pirana 3.0.0 and R 4.1.2. Log-transform-both-sides was",
      "applied. 50 of 403 samples (12%) were below the lower limit of",
      "quantification (430 ng/mL = 0.43 mg/L per Figure 2) and were",
      "RETAINED and handled with Beal's M3 method rather than excluded;",
      "M3 is an estimation-time likelihood method and has no",
      "simulation-time counterpart, so it is not encoded here.",
      "Free (unbound) concentrations were NOT measured: only total",
      "clavulanic acid concentrations were available, and the model was",
      "NOT corrected for protein binding. The authors note adult",
      "clavulanic acid protein binding is only moderate (up to 30%) and",
      "that neonates start with lower protein binding, so the paper's",
      "%fT>CT target attainment is computed on total concentrations.",
      "Serum creatinine was available for every dataset EXCEPT the RAIN",
      "oral subjects and was therefore omitted from the covariate analysis",
      "entirely (Methods), which the authors flag as a possible source of",
      "unexplained clearance variability. Study centre was screened on",
      "clearance and not retained. Covariates screened and NOT retained:",
      "sex, gestational age, postmenstrual age, and study centre; only",
      "postnatal age was retained (dOFV -14.0, P < 0.05).",
      "Model evaluation was internal only (nonparametric bootstrap,",
      "n = 1000, of which 82.8% of runs were successful; GOF plots; VPC).",
      "No external validation was performed.",
      "Limited sampling in the oral subjects precluded estimation of",
      "interindividual variability on either bioavailability or Ka, so",
      "neither carries an eta."
    )
  )

  ini({
    # --- Structural parameters (reference: WT 3.9 kg, PNA 55.5 d) ---
    lka <- log(0.781)
    label("Absorption rate constant (1/h)")
    # Schouwenburg 2025 Table 2, 'Ka (hour-1)': 0.781 (RSE 7%);
    # bootstrap 0.801 (95% CI 0.499-1.610). Discussion restates it as
    # "moderate absorption (0.781 hour-1)" and contrasts it with
    # amoxicillin's 0.091 1/h.

    lcl <- log(0.675)
    label("Typical clearance TVCL at WT 3.9 kg and PNA 55.5 d (L/h)")
    # Schouwenburg 2025 Table 2, 'TVCL (L/h)': 0.675 (RSE 9%);
    # bootstrap 0.685 (95% CI 0.581-0.782). Discussion restates it as
    # "clavulanic acid clearance is ~20% of that of amoxicillin
    # (0.68 vs. 3.22 L hour-1)".

    lvc <- log(1.88)
    label("Typical central volume TVV at WT 3.9 kg (L)")
    # Schouwenburg 2025 Table 2, 'TVV (L)': 1.88 (RSE 11%);
    # bootstrap 1.907 (95% CI 1.558-2.274).

    lfdepot <- log(0.244)
    label("Oral bioavailability (fraction)")
    # Schouwenburg 2025 Table 2, 'F (%)': 24.4 (RSE 14%);
    # bootstrap 24.7 (95% CI 19-31.4). The abstract, Results and
    # Discussion all quote 24.4%; only the Study Highlights box says
    # "24.5%", which is a rounding slip in a summary box and is not used.
    # Estimated from the oral (RAIN) sub-cohort, PNA 0-12 days, so this is
    # a neonatal value; no covariate on F was estimated.

    # --- Allometric exponents: FIXED a priori, not estimated ---
    # Methods, 'Covariate relationship analysis': "Refinements applied to
    # the model to provide better parameter estimations were a priori
    # allometric scaling with a fixed exponent (i.e., 0.75 on clearances,
    # 1 on distribution volume)". Table 2 prints 0.75 inline in the CL
    # equation with no RSE and no bootstrap CI; the volume equation is
    # printed as V = TVV * (BW/3.9) with the exponent of 1 implicit,
    # made explicit by the Methods sentence above.
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on clearance (unitless)")
    # Schouwenburg 2025 Table 2 equation:
    #   CLclav (L/hour) = TVCL * (BW/3.9)^0.75 * (PNA/55.5)^theta_PNA

    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on central volume (unitless)")
    # Schouwenburg 2025 Table 2 equation: Vclav (L) = TVV * (BW/3.9)

    # --- Estimated maturation covariate effect on clearance ---
    e_pna_cl <- 0.207
    label("Postnatal-age power exponent on clearance (unitless; reference 55.5 days)")
    # Schouwenburg 2025 Table 2, 'theta_PNA' (under 'Covariate
    # relationships'): 0.207 (RSE 24%); bootstrap 0.211 (95% CI
    # 0.124-0.290). Retained with dOFV -14.0 (P < 0.05).
    # The Results text calls this "an exponential function" while Table 2
    # prints the POWER form (PNA/55.5)^theta_PNA; the printed equation is
    # authoritative (see covariateData$PNA$notes NOTE 2).

    # --- Interindividual variability (clearance only) ---
    etalcl ~ 0.142129
    # Schouwenburg 2025 Table 2, 'Inter-individual variability (IIV) /
    # Clearance (%CV)': 37.7 (RSE 15%) [shrinkage 27%]; bootstrap
    # 35.9 (95% CI 24.5-46.0). The bootstrap CI is reported on the %CV
    # scale here (24.5-46.0 brackets the 35.9 %CV point estimate), so the
    # row is unambiguously a CV percentage, not a variance.
    # omega^2 = 0.377^2 = 0.142129, matching the direct-square convention
    # used by the same author group in Keij_2023_amoxicillin. The strict
    # log-normal conversion log(CV^2 + 1) = 0.13290 (36.5 %CV) is an
    # equally defensible reading; the difference is 3% on the SD scale and
    # immaterial. Methods state IIV was placed on clearance only; the
    # Discussion confirms "Limited sampling in oral patients precluded the
    # estimation of the inter-individual variability from bioavailability
    # or Ka", so no eta is encoded on lka, lvc or lfdepot.

    # --- Residual unexplained variability ---
    propSd <- 0.778
    label("Proportional residual error (fraction)")
    # Schouwenburg 2025 Table 2, 'Residual variability / Proportional
    # error (%)': 77.8 (RSE 7%); bootstrap 77.5 (95% CI 68.0-87.0).
    # Methods: "Log transformation both sides (model-based predictions and
    # observed clavulanic acid concentrations) was applied and a
    # proportional residual error model (additive on the log-scale) was
    # used"; Table 2 footnote repeats "Log transformation on both sides,
    # log-scale additive error model applied". A log-scale additive SD is
    # proportional error in nlmixr2's linear space, so it maps to
    # propSd = 0.778 (the row is a percentage of the prediction, matching
    # the fractional-scale 0.132 'Proportional error' row in the same
    # group's Keij 2023 Table 2).
    # This RUV is unusually large; the Discussion independently accounts
    # for it: clavulanic acid is chemically unstable ("Unpublished
    # research from our group demonstrated clavulanic acid stability to be
    # poor when samples were not immediately stored ... frosting/defrosting
    # severely affects the stability"), and a group of peak samples from
    # the Dhont cohort showed CWRES <= -2.5, "which might indicate
    # inaccurate documentation of sampling and dose administration".
  })

  model({
    # 1. Derived covariate terms
    #    Schouwenburg 2025 Table 2:
    #      CLclav (L/hour) = TVCL * (BW/3.9)^0.75 * (PNA/55.5)^theta_PNA
    #      Vclav  (L)      = TVV  * (BW/3.9)
    #    The PNA reference is 55.5 DAYS in the paper; the canonical PNA
    #    column is in MONTHS, so the reference is converted once here.
    #    The ratio is unit-invariant, so the exponent transfers unchanged.
    pna_ref_months <- 55.5 / 30.4375
    f_pna <- (PNA / pna_ref_months)^e_pna_cl

    # 2. Individual parameters
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 3.9)^e_wt_cl * f_pna
    vc <- exp(lvc) * (WT / 3.9)^e_wt_vc

    # 3. Micro-constants
    kel <- cl / vc

    # 4. ODE system (one compartment with a first-order absorption depot)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # 5. Bioavailability -- applies to oral (depot) doses only. Intravenous
    #    doses are administered directly into `central` and therefore carry
    #    F = 1, matching the study design (Table 1: RAIN oral and
    #    intravenous; Durham, De Cock and Dhont intravenous only).
    f(depot) <- exp(lfdepot)

    # 6. Observation and error
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
