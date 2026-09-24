Mian_2019_acetaminophen <- function() {
  description <- "Parent-and-metabolites population PK model for intravenous acetaminophen (paracetamol) and its sulfate, glucuronide, and combined oxidative (cysteine + mercapturate) metabolites in infants and young children after cardiac surgery with cardiopulmonary bypass (Mian 2019). One-compartment plasma disposition for the parent and for each of the three metabolite pools. Formation clearance of each metabolite is a literature-assumed fixed fraction of the total acetaminophen elimination clearance (sulfation 0.49, glucuronidation 0.36, oxidation 0.10, unchanged 0.05). Total body weight enters the parent and all three metabolite elimination clearances as a linear (not allometric power) function centred on the population median weight of 6.1 kg. Down syndrome, age, sex, cardiopulmonary bypass time and RACHS-1 category were screened and none was retained."
  reference <- paste(
    "Mian P, Valkenburg AJ, Allegaert K, Koch BCP, Breatnach CV,",
    "Knibbe CAJ, Tibboel D, Krekels EHJ (2019).",
    "Population pharmacokinetic modeling of acetaminophen and metabolites",
    "in children after cardiac surgery with cardiopulmonary bypass.",
    "J Clin Pharmacol 59(6):847-855.",
    "doi:10.1002/jcph.1373.",
    sep = " "
  )
  vignette <- "Mian_2019_acetaminophen"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. All concentrations in Mian 2019 were expressed in
  # ACETAMINOPHEN EQUIVALENTS via molecular-weight conversion (Methods,
  # 'Population PK Analysis and Internal Validation'), so every metabolite
  # state holds acetaminophen-equivalent mass and the parent-to-metabolite
  # transfer is 1:1 on a mass basis.
  compartmentData <- list(
    central = list(analyte = "acetaminophen", units = "mg", specimen = "plasma", verified = TRUE),
    central_sulf = list(
      analyte = "acetaminophen sulfate (acetaminophen equivalents)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_gluc = list(
      analyte = "acetaminophen glucuronide (acetaminophen equivalents)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_cysmer = list(
      analyte = "acetaminophen cysteine + acetaminophen mercapturate (acetaminophen equivalents)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight at surgery",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The ONLY covariate retained in the final model (Results, 'Covariate Model').",
        "Enters the parent and all three metabolite elimination clearances as a LINEAR",
        "function of WT / 6.1 (6.1 kg = population median weight), not as an allometric",
        "power function: Table 2 prints CLE = theta_slope * (BW / 6.1) + theta_intercept.",
        "The paper explicitly tested an exponential (power) form and it did not improve",
        "the fit. Time-fixed at the weight recorded at surgery. Extrapolation outside the",
        "studied 4.0-12.9 kg range is NOT justified (Discussion): the linear relationship",
        "crosses zero clearance at about 2.2 kg and is negative below it.",
        sep = " "
      ),
      source_name = "BW"
    )
  )

  # Covariates screened in the covariate analysis but NOT retained in the final
  # model (Methods, 'Covariate Model'; Results, 'Covariate Model'). Documented
  # here so the paper's covariate screen is preserved without declaring
  # covariates that model() never references.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at surgery",
      units = "years",
      type = "continuous",
      notes = paste(
        "Tested as a descriptor of growth and maturation alongside body weight.",
        "Body weight was the better descriptor and age was not retained",
        "(Results, 'Covariate Model'; Discussion). Source column was age in days",
        "(median 177 days, range 92-944).",
        sep = " "
      ),
      source_name = "Age at surgery, days"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "male (SEXF = 0)",
      notes = "Screened as a categorical covariate; not statistically significant on any parameter.",
      source_name = "Male"
    ),
    T_CPB = list(
      description = "Cardiopulmonary bypass duration",
      units = "minutes",
      type = "continuous",
      notes = "Screened; not statistically significant. Median 111 min (Down syndrome) / 106 min (no Down syndrome), Table 1.",
      source_name = "CPB time, min"
    ),
    RACHS1 = list(
      description = "Risk Adjustment for Congenital Heart Surgery (RACHS-1) surgical-risk category",
      units = "(category; integer 1-6)",
      type = "categorical",
      reference_category = "paper-specific; only categories 2 and 3 occur in this cohort",
      notes = "Screened; not statistically significant. Category 2 in 16/30 and category 3 in 14/30 children (Table 1).",
      source_name = "RACHS score"
    ),
    DIS_DOWN = list(
      description = "Down syndrome (trisomy 21) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "no Down syndrome (DIS_DOWN = 0)",
      notes = paste(
        "The paper's primary covariate of interest: 17 of 30 children had Down syndrome.",
        "Tested as a fractional change on every model parameter at each stage of the",
        "covariate analysis and never statistically significant, so the two groups were",
        "pooled (Results, 'Covariate Model'; Conclusions). Documentation only -- the",
        "final model carries no Down-syndrome term.",
        sep = " "
      ),
      source_name = "Down syndrome"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 30L,
    n_studies = 1L,
    age_range = "92-944 days (inclusion criterion 3-36 months)",
    age_median = "177 days",
    weight_range = "4.0-12.9 kg",
    weight_median = "6.1 kg",
    sex_female_pct = 56.7,
    race_ethnicity = NA_character_,
    disease_state = paste(
      "Infants and young children in the immediate postoperative period after cardiac",
      "surgery with cardiopulmonary bypass for atrial septal defect, ventricular septal",
      "defect, atrioventricular septal defect or tetralogy of Fallot repair.",
      "17 of 30 children had Down syndrome (trisomy 21).",
      sep = " "
    ),
    dose_range = paste(
      "Three intravenous acetaminophen doses at 8-hour intervals, each infused over",
      "15 minutes: 7.5 mg/kg for children < 10 kg and 15 mg/kg for children >= 10 kg.",
      sep = " "
    ),
    regions = "Ireland (Our Lady's Children's Hospital, Dublin)",
    n_observations = paste(
      "161 acetaminophen, 161 acetaminophen sulfate, 161 acetaminophen glucuronide,",
      "161 acetaminophen cysteine and 153 acetaminophen mercapturate concentrations",
      "(3-9 samples per patient). All acetaminophen glutathione concentrations were",
      "below the limit of quantification and that metabolite is not in the model.",
      sep = " "
    ),
    notes = paste(
      "Demographics from Table 1 of Mian 2019. Acetaminophen concentrations were measured",
      "in scavenged blood samples from a previously published morphine / midazolam study",
      "in the same cohort. Sex is reported as 7/17 male in the Down-syndrome group and",
      "6/13 male in the group without Down syndrome, i.e. 13/30 male and 17/30 (56.7%)",
      "female. Estimation used NONMEM 7.2 FOCE-I with ADVAN13.",
      sep = " "
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # All final parameter estimates are from Table 2 of Mian 2019 (page 851).
    #
    # Body-weight model. Table 2 prints each elimination clearance as
    #     CLE = theta_slope * (BW / 6.1) + theta_intercept
    # (Equation 1 with the exponent fixed to 1, i.e. the linear form). That is
    # algebraically identical to the centred form used in model():
    #     CLE = CLE_6.1 + theta_slope * (BW / 6.1 - 1),  CLE_6.1 = slope + intercept
    # so the reference clearance below is `theta_slope + theta_intercept` and the
    # covariate-effect parameter is `theta_slope` itself (units L/h, a LINEAR
    # slope -- not an allometric exponent). The parent reference value
    # 1.50 - 0.54 = 0.96 L/h is confirmed verbatim by the paper's own text
    # ('total elimination clearance and distribution volume of acetaminophen are
    # 0.96 L/h and 7.96 L' for the typical 6.1 kg child, Results, 'Comparison
    # With Noncardiac Surgery').
    # ---------------------------------------------------------------------

    # Parent acetaminophen disposition.
    lvc <- log(7.96); label("Acetaminophen central volume of distribution (L)") # Table 2: V_APAP = 7.96 L (RSE 10%); no covariate retained on V
    lcl <- log(1.50 - 0.54); label("Acetaminophen total elimination clearance at WT = 6.1 kg (L/h)") # Table 2: CLE_APAP = theta1 * (BW/6.1) + theta2 with theta1 = 1.50 (RSE 27%), theta2 = -0.54 (RSE 61%); 1.50 - 0.54 = 0.96 L/h, matching the paper's stated typical value
    e_wt_cl <- 1.50; label("Linear slope of WT / 6.1 on acetaminophen elimination clearance (L/h)") # Table 2: theta1 = 1.50 L/h (RSE 27%)

    # Acetaminophen sulfate.
    lvc_sulf <- log(0.68); label("Acetaminophen sulfate volume of distribution (L)") # Table 2: V_sulf = 0.68 L (RSE 29%)
    lcle_sulf <- log(0.65 - 0.24); label("Acetaminophen sulfate elimination clearance at WT = 6.1 kg (L/h)") # Table 2: CLE_sulf = theta3 * (BW/6.1) + theta4 with theta3 = 0.65 (RSE 19%), theta4 = -0.24 (RSE 35%)
    e_wt_cle_sulf <- 0.65; label("Linear slope of WT / 6.1 on acetaminophen sulfate elimination clearance (L/h)") # Table 2: theta3 = 0.65 L/h (RSE 19%)

    # Acetaminophen glucuronide.
    lvc_gluc <- log(1.69); label("Acetaminophen glucuronide volume of distribution (L)") # Table 2: V_gluc = 1.69 L (RSE 29%)
    lcle_gluc <- log(1.41 - 0.53); label("Acetaminophen glucuronide elimination clearance at WT = 6.1 kg (L/h)") # Table 2: CLE_gluc = theta5 * (BW/6.1) + theta6 with theta5 = 1.41 (RSE 21%), theta6 = -0.53 (RSE 50%)
    e_wt_cle_gluc <- 1.41; label("Linear slope of WT / 6.1 on acetaminophen glucuronide elimination clearance (L/h)") # Table 2: theta5 = 1.41 L/h (RSE 21%)

    # Combined oxidative metabolites (acetaminophen cysteine + acetaminophen
    # mercapturate), modelled in a single compartment on the assumption that
    # both NAPQI-derived species share a distribution volume (Methods,
    # 'Structural and Statistical Model').
    lvc_cysmer <- log(0.042); label("Combined oxidative (cysteine + mercapturate) volume of distribution (L)") # Table 2: V_ox = 0.042 L (RSE 18%)
    lcle_cysmer <- log(40.86 - 1.26); label("Combined oxidative (cysteine + mercapturate) elimination clearance at WT = 6.1 kg (L/h)") # Table 2: CLE_ox = theta7 * (BW/6.1) + theta8 with theta7 = 40.86 (RSE 25%), theta8 = -1.26 (RSE 28%)
    e_wt_cle_cysmer <- 40.86; label("Linear slope of WT / 6.1 on combined oxidative elimination clearance (L/h)") # Table 2: theta7 = 40.86 L/h (RSE 25%)

    # Fractions of the total acetaminophen elimination clearance routed to each
    # pathway. Not estimated: the metabolite sub-model is unidentifiable without
    # them, so the authors held them at published values (Methods, 'Structural
    # and Statistical Model': sulfation 0.49 and glucuronidation 0.36 from
    # children aged 3-9 years, oxidation 0.10 and unchanged 0.05 from healthy
    # adults). The four shares sum to exactly 1.
    fm_sulf <- fixed(0.49); label("Fraction of acetaminophen elimination clearance forming sulfate (unitless; literature value assumed by the authors)") # Methods: 'a fraction of 0.49 ... of total elimination acetaminophen clearance'
    fm_gluc <- fixed(0.36); label("Fraction of acetaminophen elimination clearance forming glucuronide (unitless; literature value assumed by the authors)") # Methods: 'a fraction of ... 0.36 of total elimination acetaminophen clearance'
    fm_cysmer <- fixed(0.10); label("Fraction of acetaminophen elimination clearance forming oxidative metabolites (unitless; literature value assumed by the authors)") # Methods: 'the oxidative metabolites and the unchanged clearance of acetaminophen ... namely, 0.10 and 0.05 of the total acetaminophen clearance'
    fm_other <- fixed(0.05); label("Fraction of acetaminophen elimination clearance eliminated unchanged (unitless; literature value assumed by the authors)") # Methods: '... 0.10 and 0.05 of the total acetaminophen clearance'

    # Inter-individual variability. Table 2 reports these directly as omega^2
    # (the 'Interindividual variability [omega^2]' block), log-normally
    # distributed on every parameter (Methods, 'Structural and Statistical
    # Model'), so the printed value is used as the variance with no CV%
    # back-transformation. Shrinkage is given in square brackets in Table 2.
    etalvc ~ 0.189 # Table 2 row 'V APAP' = 0.189 (RSE 27%) [shrinkage 11%]
    etalcl ~ 0.185 # Table 2 row 'CLE APAP' = 0.185 (RSE 27%) [shrinkage 6%]
    etalvc_sulf ~ 0.726 # Table 2 row 'V sulf' = 0.726 (RSE 52%) [shrinkage 12%]
    etalcle_sulf ~ 0.189 # Table 2 row 'CLE sulf' = 0.189 (RSE 32%) [shrinkage 6%]
    etalvc_gluc ~ 0.927 # Table 2 row 'V gluc' = 0.927 (RSE 50%) [shrinkage 15%]
    etalcle_gluc ~ 0.129 # Table 2 row 'CLE gluc' = 0.129 (RSE 39%) [shrinkage 13%]
    etalvc_cysmer ~ 0.600 # Table 2 row 'V ox' = 0.600 (RSE 49%) [shrinkage 9%]
    etalcle_cysmer ~ 0.552 # Table 2 row 'CLE ox' = 0.552 (RSE 32%) [shrinkage 9%]

    # Residual error: a proportional model best described each compound
    # (Results, 'Structural and Statistical Model'). Table 2 reports the
    # 'Residual variability [sigma^2]' block as variances, so the nlmixr2
    # proportional SD is the square root of the printed value.
    propSd <- sqrt(0.146); label("Acetaminophen proportional residual SD (fraction)") # Table 2: sigma^2 = 0.146 (RSE 29%) -> SD = 0.382
    propSd_sulf <- sqrt(0.0507); label("Acetaminophen sulfate proportional residual SD (fraction)") # Table 2: sigma^2 = 0.0507 (RSE 15%) -> SD = 0.225
    propSd_gluc <- sqrt(0.0813); label("Acetaminophen glucuronide proportional residual SD (fraction)") # Table 2: sigma^2 = 0.0813 (RSE 14%) -> SD = 0.285
    propSd_cysmer <- sqrt(0.0494); label("Combined oxidative metabolites proportional residual SD (fraction)") # Table 2: sigma^2 = 0.0494 (RSE 12%) -> SD = 0.222
  })

  model({
    # 1. Derived covariate term. Body weight centred on the population median of
    #    6.1 kg. `wt_cen` is zero at WT = 6.1 kg, so each clearance below equals
    #    its Table 2 reference value there.
    wt_cen <- WT / 6.1 - 1

    # 2. Individual parameters. The linear weight relationship acts on the
    #    typical value and the log-normal IIV multiplies the covariate-adjusted
    #    typical value, matching NONMEM's TVCL * EXP(ETA) parameterisation.
    vc <- exp(lvc + etalvc)
    cl <- (exp(lcl) + e_wt_cl * wt_cen) * exp(etalcl)
    vc_sulf <- exp(lvc_sulf + etalvc_sulf)
    cle_sulf <- (exp(lcle_sulf) + e_wt_cle_sulf * wt_cen) * exp(etalcle_sulf)
    vc_gluc <- exp(lvc_gluc + etalvc_gluc)
    cle_gluc <- (exp(lcle_gluc) + e_wt_cle_gluc * wt_cen) * exp(etalcle_gluc)
    vc_cysmer <- exp(lvc_cysmer + etalvc_cysmer)
    cle_cysmer <- (exp(lcle_cysmer) + e_wt_cle_cysmer * wt_cen) * exp(etalcle_cysmer)

    # 3. Micro-constants. `kel` is the TOTAL acetaminophen elimination rate
    #    constant; each pathway takes its fixed share of it.
    kel <- cl / vc
    kel_sulf <- cle_sulf / vc_sulf
    kel_gluc <- cle_gluc / vc_gluc
    kel_cysmer <- cle_cysmer / vc_cysmer

    # 4. ODE system, Figure 1 of Mian 2019. One compartment for the parent with
    #    four parallel elimination arms (sulfation, glucuronidation, oxidation,
    #    unchanged) and one compartment for each measured metabolite pool. The
    #    four fractions are written out in the parent equation rather than
    #    collapsed to `-kel * central` so that the mass balance
    #    fm_sulf + fm_gluc + fm_cysmer + fm_other = 1 is visible and auditable.
    #    Because every concentration is expressed in acetaminophen equivalents,
    #    the parent-to-metabolite transfer is 1:1 in mass.
    d/dt(central) <- -(fm_sulf + fm_gluc + fm_cysmer + fm_other) * kel * central
    d/dt(central_sulf) <- fm_sulf * kel * central - kel_sulf * central_sulf
    d/dt(central_gluc) <- fm_gluc * kel * central - kel_gluc * central_gluc
    d/dt(central_cysmer) <- fm_cysmer * kel * central - kel_cysmer * central_cysmer

    # 5. Observations. Acetaminophen-equivalent plasma concentrations (mg/L);
    #    doses are in mg and volumes in L.
    Cc <- central / vc
    Cc_sulf <- central_sulf / vc_sulf
    Cc_gluc <- central_gluc / vc_gluc
    Cc_cysmer <- central_cysmer / vc_cysmer

    Cc ~ prop(propSd)
    Cc_sulf ~ prop(propSd_sulf)
    Cc_gluc ~ prop(propSd_gluc)
    Cc_cysmer ~ prop(propSd_cysmer)
  })
}
