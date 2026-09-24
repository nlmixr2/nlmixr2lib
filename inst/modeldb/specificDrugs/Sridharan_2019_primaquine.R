Sridharan_2019_primaquine <- function() {
  description <- "One-compartment population PK model for primaquine after a single oral 15-mg dose in 53 Indian adults: 13 healthy volunteers, 12 with mild and 6 with moderate hepatic dysfunction (Child-Pugh), and 22 with renal dysfunction. First-order absorption and first-order elimination; during estimation the absorption rate constant was constrained above the elimination rate constant to avoid flip-flop. Apparent volume of distribution is normalized linearly to a 70-kg person and rises 3.86-fold in moderate hepatic dysfunction, the only covariate retained; mild hepatic dysfunction, renal dysfunction, age and sex were screened and dropped, and neither hepatic nor renal dysfunction affected clearance. Exponential between-subject variability on CL/F, V/F and Ka, with combined proportional-plus-additive residual error."
  reference <- paste(
    "Sridharan K, Sannala CKR, Mallayasamy S, Chaturvedula A, Kadam P, Hase N,",
    "Shukla A, Gogtay N, Thatte U.",
    "Population pharmacokinetics of primaquine and the effect of hepatic and",
    "renal dysfunction: An exploratory approach.",
    "Indian J Pharmacol. 2019;51(1):17-23. doi:10.4103/ijp.ijp_230_16.",
    "Structural model from Methods, 'Population pharmacokinetic modeling', and",
    "Results, 'Model development and evaluation'; parameter values from Table 2",
    "('Population estimates' column); covariate function from the Methods",
    "equation 'TVP = P x (1 + theta_mild x FLAG) x (1 + theta_mod x FLAG1)'.",
    sep = " "
  )
  vignette <- "Sridharan_2019_primaquine"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline. Methods, 'Population pharmacokinetic modeling': 'Volume of distribution was modeled, which was normalized to 70-kg person ... in all models tested.' The paper states a normalization to 70 kg but never prints an exponent, so the scaling is the linear one that phrase denotes and e_wt_vc is held at 1 rather than estimated. No weight scaling is applied to CL/F: body weight was among the covariates screened by stepwise forward inclusion and backward elimination and was not retained anywhere except as this a priori volume normalization. Group median (range) weights, Table 1: healthy 66 kg (62-95), mild hepatic dysfunction 65.5 kg (49-75), moderate hepatic dysfunction 56 kg (47-67), renal dysfunction 55 kg (40-74).",
      source_name = "WT"
    ),
    HEPIMP_MOD = list(
      description = "Moderate hepatic dysfunction indicator",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = normal hepatic function, mild hepatic dysfunction, or renal dysfunction (the paper's pooled 'All Other subjects' stratum)",
      notes = "Classification scheme is Child-Pugh, not NCI ODWG. Methods, 'Ethics and study participants': 'the individuals with hepatic dysfunction were classified into mild or moderate degree based on Child-Pugh's criteria'. The paper does not state which Child-Pugh class letters map to 'mild' and 'moderate'. Carried in the source NONMEM dataset as FLAG1 (Methods: 'FLAG1 = 1 and FLAG2 = 0, for moderate hepatic function'), which the covariate equation writes as the second factor. 6 of 53 participants were moderate. This is the only covariate retained in the final model, and it acts only on V/F: Results, 'Only the moderate hepatic dysfunction showed to be significant on the volume of distribution', and 'There was no significant effect on absorption rate constant in the moderate hepatic failure group.' It explained 25% of the between-subject variability in V/F.",
      source_name = "FLAG1"
    )
  )

  # Covariates the authors screened by stepwise forward inclusion (P = 0.05)
  # and backward elimination (P = 0.01) but did NOT retain in the final model
  # (Methods, 'Population pharmacokinetic modeling': "Covariates tested were
  # body weight, age, gender, hepatic dysfunction, and renal dysfunction").
  # Table 2 publishes no point estimate for any of them, so none can be
  # encoded; they are recorded here so the paper's covariate screen survives.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on Ka, V/F and CL/F, not retained. Results: 'Covariate model was developed to study the impact of age and hepatic disorder on absorption (Ka), volume of distribution, and CL. Of these, the hepatic function had a significant effect on volume of distribution.' Group median (range), Table 1: healthy 25.5 years (19-34), mild hepatic dysfunction 45 (22-61), moderate hepatic dysfunction 50.5 (26-61), renal dysfunction 43 (20-60)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "categorical",
      notes = "Screened as 'gender', not retained. Methods gives the power form used for binary covariates, 'TVP = P x theta_COV^Gender', but no estimate is published. Male:female ratios, Table 1: healthy 5:1, mild hepatic dysfunction 5:1, moderate hepatic dysfunction all male, renal dysfunction 2:1."
    ),
    HEPIMP_MILD = list(
      description = "Mild hepatic dysfunction indicator (Child-Pugh)",
      units = "(binary)",
      type = "categorical",
      notes = "Screened as the first factor of the hepatic covariate equation (theta_mild, carried in the source dataset as FLAG2), not retained: Results, 'Only the moderate hepatic dysfunction showed to be significant on the volume of distribution.' Table 2 has no 'Mild HD on V' row, so theta_mild has no published value. 12 of 53 participants were mild. Discussion offers the mechanism: 'there was no effect of mild hepatic dysfunction on the volume of distribution and could possibly that the protein binding differences may not be apparent until moderate dysfunction develops.'"
    ),
    RENALIMP = list(
      description = "Renal dysfunction indicator (any degree)",
      units = "(binary)",
      type = "categorical",
      notes = "Screened with the same proportional FLAG-variable function as hepatic dysfunction, not retained on any parameter. Methods, 'Ethics and study participants': 'patients with renal dysfunction were diagnosed according to the National Kidney Foundation Kidney Disease Outcomes Quality Initiative based on their serum creatinine levels.' Discussion: 'we did not see a difference in CL for renal and hepatic dysfunction subjects ... the drug may not have significant impact of renal dysfunction as the major elimination pathway was metabolic CL.' 22 of 53 participants had renal dysfunction; they sit in the reference stratum of HEPIMP_MOD."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "primaquine", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "primaquine", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 53,
    n_studies = 3,
    age_range = "19-61 years",
    age_median = "25.5 years (healthy), 45 years (mild hepatic dysfunction), 50.5 years (moderate hepatic dysfunction), 43 years (renal dysfunction)",
    weight_range = "40-95 kg",
    weight_median = "66 kg (healthy), 65.5 kg (mild hepatic dysfunction), 56 kg (moderate hepatic dysfunction), 55 kg (renal dysfunction)",
    sex_female_pct = NA_real_,
    disease_state = "13 normal healthy individuals, 12 patients with mild and 6 with moderate hepatic dysfunction graded by Child-Pugh criteria, and 22 patients with renal dysfunction diagnosed by National Kidney Foundation KDOQI criteria from serum creatinine.",
    dose_range = "Single oral 15-mg primaquine phosphate tablet (Bharat Parenterals, India) given post-breakfast with 200 mL water after an overnight fast; liquids restricted 2 h and food 4 h post-dose.",
    regions = "India (Seth GS Medical College and KEM Hospital, Mumbai); retrospective pooling of three single-centre studies registered as CTRI/2011/06/001803 (healthy), CTRI/2011/06/001794 (hepatic dysfunction) and CTRI/2010/091/000356 (renal dysfunction), conducted April-December 2013.",
    notes = "Baseline demographics from Table 1; 458 concentration records across the 53 participants. Sampling: 0 h (pre-dose) and 0.5, 1.0, 1.5, 2, 3, 4, 6, 8, 12 and 24 h post-dose, assayed by reversed-phase HPLC. Sex is reported only as per-group male:female ratios (5:1, 5:1, all male, 2:1), which do not resolve to exact counts for the 12- and 22-subject groups, so sex_female_pct is left NA. Model qualification used a hepatic-dysfunction-stratified VPC (n = 1000 simulations) and a nonparametric bootstrap (n = 2000 resamples, 98% minimizing successfully); condition number 10.83."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters, Table 2, 'Population estimates / Population
    # mean' column, with the relative standard error from the adjacent
    # 'RSE (%)' column. Values are apparent (per-bioavailability): the
    # paper never estimates F, so CL is CL/F and V is V/F throughout, as
    # Table 1 labels them.
    #
    # Reported in linear units, not on the log scale, so each is wrapped
    # in log() here to match the repository's log-parameterization.
    # ------------------------------------------------------------------
    lka <- log(0.9462)
    label("Apparent first-order absorption rate constant (1/h)") # Table 2 'Ka theta 0.9462', RSE 18.08%; bootstrap median 0.9461, 95% CI 0.65-1.34
    lcl <- log(39.14)
    label("Apparent clearance (L/h)") # Table 2 'CL theta 39.14', RSE 10.34%; bootstrap median 39.175, 95% CI 31.36-48.30
    lvc <- log(438)
    label("Apparent central volume of distribution at WT = 70 kg and normal-to-mild hepatic function (L)") # Table 2 'V theta 438', RSE 11.9%; bootstrap median 439.646, 95% CI 345.37-550.53

    # Body-weight normalization of V/F. Methods: "Volume of distribution
    # was modeled, which was normalized to 70-kg person ... in all models
    # tested." The exponent is never printed; a plain normalization to a
    # reference weight is the linear one, so it is held at 1 rather than
    # estimated. See covariateData$WT$notes.
    e_wt_vc <- fixed(1)
    label("Body-weight exponent on the apparent central volume of distribution (unitless)") # Methods, 'Population pharmacokinetic modeling'; implied by 'normalized to 70-kg person', not separately tabulated

    # Moderate hepatic dysfunction on V/F. The Methods covariate function
    # is proportional, TVP = P x (1 + theta_mild x FLAG) x (1 + theta_mod
    # x FLAG1), so this coefficient is the FRACTIONAL increase and the
    # moderate-dysfunction multiplier is 1 + 2.859 = 3.859. The paper's
    # own Figure 3 pins that reading: at 0.75 mg/kg single dose it shows
    # a peak near 95 ng/mL for 'All Other Subjects' and near 30 ng/mL for
    # 'Moderate HD Subjects', which the 3.859-fold volume reproduces
    # (93.7 and 28.3 ng/mL) and a bare 2.859-fold volume does not
    # (93.7 and 37.3 ng/mL). theta_mild has no published estimate because
    # mild dysfunction was not retained; see covariatesDataExcluded.
    e_hepimp_mod_vc <- 2.859
    label("Fractional increase in the apparent central volume of distribution in moderate hepatic dysfunction (unitless)") # Table 2 'Moderate HD on V theta 2.859', RSE 30.17%; bootstrap median 2.8589, 95% CI 1.26-5.08

    # ------------------------------------------------------------------
    # Between-subject variability, Table 2 'Between subject variability'
    # block, published as % CV. The same '(% CV)' column header is used
    # for the proportional residual error, where it can only mean
    # 100 x sqrt(variance); the etas are read on that same convention, so
    # each variance below is (%CV/100)^2. The alternative log-normal
    # reading, log(1 + CV^2), would give 0.353 / 0.367 / 0.470.
    # ------------------------------------------------------------------
    etalka ~ 0.599695 # Table 2 'BSV on Ka (per cent CV) 77.44', RSE 27.44; bootstrap median 74.5582; 0.7744^2
    etalcl ~ 0.423541 # Table 2 'BSV on CL (per cent CV) 65.08', RSE 19.982; bootstrap median 64.1992; 0.6508^2
    etalvc ~ 0.443290 # Table 2 'BSV on V (per cent CV) 66.58', RSE 19.642; bootstrap median 64.8353; 0.6658^2

    # ------------------------------------------------------------------
    # Residual error. Methods equation 2 is a combined additive and
    # proportional model; as printed the two epsilon subscripts are
    # transposed (the multiplicative term is labelled 'add' and the
    # standalone term 'prop'). Table 2's units settle the orientation:
    # the additive term carries ng/mL and the proportional term % CV.
    # Variances add, which is nlmixr2's default add() + prop() form and
    # matches NONMEM's Y = F * (1 + EPS_prop) + EPS_add.
    # ------------------------------------------------------------------
    propSd <- 0.3206
    label("Proportional residual error (fraction of the predicted concentration)") # Table 2 'Proportional (per cent CV) 32.06', RSE 12.526; bootstrap median 31.9786
    addSd <- 1.5
    label("Additive residual error standard deviation (ng/mL)") # Table 2 'Additive (ng/ml) 1.5', RSE 102.5 (Results text says 116%); bootstrap median 1.4913, 95% CI 0.015-2.88
  })

  model({
    # ---- Individual parameters -----------------------------------------
    # Methods equation 1, P_i = TVP * exp(eta_i), log-normal on CL and V;
    # Ka carries the same exponential form (Table 2 publishes a BSV on Ka).
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    # Weight normalization to 70 kg, then the proportional moderate-
    # hepatic-dysfunction factor from the Methods covariate equation.
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc * (1 + e_hepimp_mod_vc * HEPIMP_MOD)

    kel <- cl / vc

    # ---- ODE system (amounts in ug, volumes in L) -----------------------
    # One compartment with first-order absorption and elimination:
    # "One compartment model with first-order absorption best described
    # the observed data."
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # ---- Observation ----------------------------------------------------
    # Amounts in ug over volumes in L give ug/L, which is ng/mL, the unit
    # Table 1 and Table 2 report concentrations and the additive residual
    # error in.
    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
