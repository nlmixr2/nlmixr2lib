Galluppi_2016_BIIB023 <- function() {
  description <- "Two-compartment population PK model for intravenous BIIB023, an anti-TWEAK monoclonal antibody, in healthy Chinese, Japanese and Caucasian volunteers and adults with rheumatoid arthritis (Galluppi 2016); parallel first-order linear and Michaelis-Menten elimination from the central compartment with body weight on CL and V, and sex on V."
  reference <- "Galluppi GR, Wisniacki N, Stebbins C. Population pharmacokinetic and pharmacodynamic analysis of BIIB023, an anti-TNF-like weak inducer of apoptosis (anti-TWEAK) monoclonal antibody. Br J Clin Pharmacol. 2016;82(1):118-128. doi:10.1111/bcp.12914"
  vignette <- "Galluppi_2016_BIIB023"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight (baseline)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on linear CL and on central volume V per Galluppi 2016 Table 4 (dCLdBodyWeight = 0.830; dVdBodyWeight = 0.459). The paper does not report the centering / reference weight used in the model; per the canonical undefined-reference policy this implementation uses 70 kg (close to the combined-cohort mean of ~67 kg across studies 211HV102 [58.6-69.1 kg per ethnic group, Table 1] and 211RA101 [74.1 kg, Table 2]).",
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Proportional effect on central volume V per Galluppi 2016 Table 4 dVdFemale = -0.105 (95% CI -0.168, -0.040), i.e., V_female = V_male * (1 - 0.105). The paper's reference category is male, matching the canonical SEXF = 0 (male) convention.",
      source_name        = "SEX"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Age (years)",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a covariate of V, V2, CL and CL2 but had no significant effect on the objective function value and was not retained in the final model (Galluppi 2016 Methods; Results, p124)."
    ),
    RACE_ASIAN = list(
      description        = "Asian ethnicity indicator (Chinese or Japanese vs Caucasian)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Caucasian)",
      notes              = "Screened (as 'ethnicity': Chinese, Japanese, Caucasian) but no clinically meaningful difference between groups; not retained in the final model (Galluppi 2016 Results, p124-125)."
    ),
    STUDY_RA = list(
      description        = "Rheumatoid arthritis study indicator (study 211RA101 vs healthy-volunteer study 211HV102)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteers, study 211HV102)",
      notes              = "Disease status (RA patients vs healthy volunteers) was screened as a covariate but had no significant effect on the objective function value and was not retained in the final model (Galluppi 2016 Methods; Results, p124)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 96L,
    n_observations = 1583L,
    n_studies      = 2L,
    age_range      = "18-55 years (healthy volunteers, study 211HV102); rheumatoid arthritis adults (study 211RA101, mean 54.0 +/- 9.0 years)",
    weight_range   = "58.6-69.1 kg per ethnic group (healthy volunteers, study 211HV102, Table 1 means); 74.1 +/- 11.0 kg overall (study 211RA101 adults, Table 2)",
    sex_female_pct = 56.3,
    race_ethnicity = c("Chinese (healthy, n=16)" = 16.7, "Japanese (healthy, n=16)" = 16.7, "Caucasian (healthy, n=16)" = 16.7, "Caucasian (RA, n=48)" = 50.0),
    disease_state  = "Pooled healthy volunteers (study 211HV102, Chinese / Japanese / Caucasian, n=48) and adults with rheumatoid arthritis (study 211RA101, Caucasian, n=48). Disease status was screened as a covariate and was not retained.",
    dose_range     = "Single 1-hour intravenous infusion: 3 or 20 mg/kg in study 211HV102; ascending 0.03-20 mg/kg in study 211RA101.",
    regions        = "Hong Kong (Chinese healthy volunteers) and Australia (Japanese and Caucasian healthy volunteers) in study 211HV102; study 211RA101 ascending-dose phase 1 in subjects with rheumatoid arthritis (Caucasian, see reference [14] of the paper).",
    notes          = "Combined dataset of 1583 BIIB023 serum concentrations from 96 subjects (48 per study). Predose BLQ samples set to zero; 120 postdose BLQ (~8% of dataset) treated as missing per Galluppi 2016 Methods. Concentrations measured by validated ELISA (calibration range 0.800-20.0 ug/mL in neat serum). Soluble TWEAK and TWEAK:BIIB023 complex were analysed descriptively but no structural PD model was fit (Galluppi 2016 Discussion: 'Whether the observed nonlinearity in clearance of BIIB023 is due to target-mediated drug disposition could not be confirmed due to insufficient information')."
  )

  ini({
    # Structural PK parameters - Galluppi 2016 Table 4 final-model estimates (p126).
    # Reference covariate values for the typical subject: total body weight 70 kg
    # (paper does not report the centering weight; canonical undefined-reference
    # default), male (SEXF = 0). The paper reports time in hours; central and
    # peripheral volumes in mL; clearances in mL/h; Km in ug/mL; Vmax in ug/h.
    # Volumes are converted to L (1 L = 1000 mL) and Vmax to mg/h
    # (1 mg = 1000 ug) so that all ODE terms share the unit system
    # central [mg] / vc [L] = Cc [mg/L = ug/mL] (since 1 mg/L = 1 ug/mL).
    lvc   <- log(3050 / 1000);  label("Central volume of distribution V (L)")             # Galluppi 2016 Table 4, tvV = 3050 mL
    lvp   <- log(2480 / 1000);  label("Peripheral volume of distribution V2 (L)")         # Galluppi 2016 Table 4, tvV2 = 2480 mL
    lcl   <- log(7.42 / 1000);  label("Linear clearance CL (L/h)")                        # Galluppi 2016 Table 4, tvCL = 7.42 mL/h
    lq    <- log(23.3 / 1000);  label("Intercompartmental clearance Q (L/h)")             # Galluppi 2016 Table 4, tvCL2 = 23.3 mL/h
    lkm   <- log(0.792);        label("Michaelis-Menten constant Km (ug/mL)")             # Galluppi 2016 Table 4, tvKm = 0.792 ug/mL
    lvmax <- log(29.2 / 1000);  label("Maximum Michaelis-Menten elimination rate Vmax (mg/h)")  # Galluppi 2016 Table 4, tvVmax = 29.2 ug/h

    # Covariate effects on CL and V - Galluppi 2016 Table 4.
    # Body weight enters as a power-of-ratio on both CL and V with reference 70 kg
    # (undefined-reference default); sex enters as a proportional shift on V only.
    # Final-model form (paper does not print the closed-form equation but states
    # "body weight effect on CL and V, and sex effect on V"):
    #   CL = tvCL * (WT/70)^dCLdWT
    #   V  = tvV  * (WT/70)^dVdWT * (1 + dVdSEXF * SEXF)
    dCLdWT  <-  0.830; label("Power exponent of WT/70 on linear CL (unitless)")       # Galluppi 2016 Table 4, dCLdBodyWeight
    dVdWT   <-  0.459; label("Power exponent of WT/70 on central volume V (unitless)") # Galluppi 2016 Table 4, dVdBodyWeight
    dVdSEXF <- -0.105; label("Proportional sex effect on V for SEXF = 1 (unitless)")   # Galluppi 2016 Table 4, dVdFemale

    # Inter-individual variability - Galluppi 2016 Table 4 ("Omega" rows).
    # The paper footnote explicitly defines Omega as "standard deviation of
    # the interindividual variability distribution" and the shrinkage
    # formula 1-SD(eta)/Omega uses Omega as the model-predicted SD; the
    # reported magnitudes are also consistent with the SD interpretation when
    # compared to the noncompartmental between-subject variability in Table 3
    # (NCA CL CV ~14-30% per ethnic group; SD interpretation gives 11.8% on
    # CL after weight + sex covariates absorb part of the variability).
    # nlmixr2 ini() expects variances, so each reported Omega-SD is squared:
    #   Omega-V  = 0.0345 SD -> variance 0.0345^2 = 0.001190
    #   Omega-V2 = 0.0734 SD -> variance 0.0734^2 = 0.005388
    #   Omega-CL = 0.118  SD -> variance 0.118^2  = 0.013924
    #   Omega-CL2= 0.461  SD -> variance 0.461^2  = 0.212521
    # The Correlation V/CL = 0.663 is between eta_V and eta_CL on the
    # log scale; covariance = rho * Omega-V * Omega-CL = 0.663 * 0.0345 *
    # 0.118 = 0.002699. Km and Vmax IIV were removed from the final model
    # (high shrinkage; Galluppi 2016 Methods).
    etalvc + etalcl ~ c(0.001190,
                        0.002699, 0.013924)   # Galluppi 2016 Table 4: Omega-V^2, rho*Omega-V*Omega-CL, Omega-CL^2
    etalvp ~ 0.005388                         # Galluppi 2016 Table 4 Omega-V2^2 (SD 0.0734)
    etalq  ~ 0.212521                         # Galluppi 2016 Table 4 Omega-CL2^2 (SD 0.461)

    # Residual error - Galluppi 2016 Methods ("proportional error model
    # (multiplicative) with FOCE-extended least squares") and Table 4
    # (stdev0 = 0.0953). Phoenix multiplicative form Cobs = Cpred * (1 + eps)
    # with eps ~ N(0, propSd^2) maps to nlmixr2 prop(propSd).
    propSd <- 0.0953; label("Proportional residual error (fraction)")                  # Galluppi 2016 Table 4, stdev0
  })
  model({
    # Individual PK parameters with Galluppi 2016 final-model covariate equations.
    # Reference body weight 70 kg (canonical undefined-reference default); sex
    # reference = male (SEXF = 0). The sex effect enters V as a proportional
    # shift (1 + dVdSEXF * SEXF); since dVdSEXF = -0.105, V_female = V_male * 0.895.
    cl   <- exp(lcl + etalcl) * (WT / 70)^dCLdWT
    vc   <- exp(lvc + etalvc) * (WT / 70)^dVdWT * (1 + dVdSEXF * SEXF)
    vp   <- exp(lvp + etalvp)
    q    <- exp(lq  + etalq)
    km   <- exp(lkm)
    vmax <- exp(lvmax)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment IV PK with parallel linear (CL) and Michaelis-Menten
    # (Vmax / (Km + Cc)) elimination from the central compartment. Doses enter
    # the central compartment directly via the 1-hour IV infusion. There is no
    # depot or absorption compartment. Concentrations: central in ug, vc in mL,
    # Cc in ug/mL (matches the paper's concentration units).
    Cc <- central / vc

    d/dt(central)     <- -kel * central -
                          k12 * central + k21 * peripheral1 -
                          vmax * Cc / (km + Cc)
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    Cc ~ prop(propSd)
  })
}
