Sun_2025_teicoplanin <- function() {
  description <- "One-compartment IV infusion population PK model for teicoplanin in 86 critically ill adults with sepsis, 20 of whom were receiving continuous renal replacement therapy (Sun 2025). Clearance is 0.98 L/h with no retained covariate; the central volume of distribution is 108.69 L in a male not receiving CRRT and carries two exponential categorical covariate effects, V = 108.69 * exp(-0.71 * RRT_CRRT_STATUS) * exp(-1.07 * SEXF). Both effects reduce V: CRRT halves it (to 53.4 L) and female sex reduces it to about a third (37.3 L), so a female receiving CRRT has V = 18.3 L. The CRRT direction is opposite to the usual literature finding of an expanded volume during renal replacement; the authors attribute it to fluid-overload correction being the indication for starting CRRT in this cohort. Interindividual variability is exponential on both CL (omega^2 = 0.31) and V (omega^2 = 0.09), and residual variability is additive at 0.23 mg/L. Age, body weight, serum creatinine, serum albumin and ECMO were screened but not retained, and no covariate improved the fit on clearance."
  reference   <- "Sun Q, Jian J, Zhou X, Hong Z, Yang S, Zheng Y, Wang S, Zhao M. Population pharmacokinetics of teicoplanin and dosage optimization in sepsis patients based on continuous renal replacement therapy. Front Pharmacol. 2025;16:1621959. doi:10.3389/fphar.2025.1621959"
  vignette    <- "Sun_2025_teicoplanin"
  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Sun 2025 Section 2.3 (total teicoplanin
  # in plasma by HPLC-UV, calibration range 5.63-125.00 mg/L).
  compartmentData <- list(
    central = list(analyte = "teicoplanin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    RRT_CRRT_STATUS = list(
      description        = "Subject-level binary indicator for continuous renal replacement therapy during the teicoplanin sampling period",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Source column CRRT. 1 = subject was receiving CRRT; 0 = no CRRT (Sun 2025 Table 1: 20/86, 23.26%). Modality mix within the CRRT subgroup: 14 continuous venovenous hemofiltration (CVVH) only, 5 continuous venovenous hemodiafiltration (CVVHD) only, and 1 subject who received both (Sun 2025 Results 3.1); the model treats all three as a single binary indicator and does not distinguish modality. Time-fixed at the subject level, matching the RRT_CRRT_STATUS canonical's stated convention. Enters V as exp(-0.71 * RRT_CRRT_STATUS) per Sun 2025 Equation 5, i.e. CRRT REDUCES the volume of distribution by 50.8%. This direction is the opposite of most published CRRT covariate effects; Sun 2025 Discussion attributes it to volume overload being the indication for initiating CRRT, so that CRRT-treated subjects had their expanded interstitial volume corrected. Adding CRRT to clearance did NOT improve the fit (Sun 2025 Discussion), which is also atypical for a renally eliminated drug.",
      source_name        = "CRRT"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Source column reported as 'Gender' (Sun 2025 Table 1: 51 males, 59.30%; 35 females, 40.70%). Sun 2025 Equation 5 writes the effect as exp(-1.07 * (if is Female)), so the source indicator is already female = 1 and maps onto the canonical SEXF orientation with no inversion; the reference subject is male. Sun 2025 Section 3.3 confirms the coding by naming the female simulation cohort 'Group Sex = 1'. Female sex reduces V by 65.7% relative to male. Median body weight differed by sex (males 66.0 kg, IQR 54.0-73.0; females 57.0 kg, IQR 47.0-63.5; Sun 2025 Discussion), so the sex effect on V is partly confounded with body size, which was screened but not retained.",
      source_name        = "Gender"
    )
  )

  # Screened during covariate model building but NOT retained in the final model
  # (Sun 2025 Section 2.5 lists the tested set; Section 3.2 reports that only
  # gender and CRRT on V survived forward inclusion and backward elimination, and
  # that "No covariate was found to significantly influence the CL of
  # teicoplanin"). Documentation only -- none of these appears in model().
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Sun 2025 Table 1 median (IQR) 62.00 (51.88, 70.00) kg. Screened as a covariate but not retained on either CL or V. Note that no allometric scaling was applied at all, so the reported CL and V are unnormalized whole-body values."
    ),
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = "Sun 2025 Table 1 median (IQR) 62.00 (53.00, 71.25) years. Screened but not retained."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Sun 2025 Table 1 median (IQR) 109.00 (74.00, 184.80). Table 1's column header reads 'Serum creatinine concentration (mg/dL)', which is a units error in the source: 109 mg/dL is not physiologically possible, whereas 109 umol/L (about 1.23 mg/dL) is an unremarkable ICU value consistent with the cohort's renal impairment and CRRT use. Recorded here in umol/L. Screened but not retained; the paper notes serum creatinine is a poor renal-function marker in CRRT-treated patients."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Sun 2025 Table 1 median (IQR) 31.60 (29.48, 36.90). Table 1's column header reads 'Serum albumin concentration (mg/L)', a units error in the source; the Discussion quotes the same median as '31.6 (29.5, 39.4) g/L', confirming g/L. Recorded here in g/L. (The Discussion's upper quartile of 39.4 also disagrees with Table 1's 36.90; Table 1 is taken as authoritative for the IQR.) Screened but not retained, despite teicoplanin being over 90% albumin-bound."
    ),
    ECMO_STATUS = list(
      description = "Extracorporeal membrane oxygenation treatment-status indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Sun 2025 Table 1: 4/86 subjects (4.65%) received ECMO. Collected (Section 2.1) but not retained in the final model; with only 4 ECMO subjects the effect was not estimable."
    ),
    HT = list(
      description = "Body height at baseline",
      units       = "cm",
      type        = "continuous",
      notes       = "Sun 2025 Table 1 median (IQR) 165.00 (156.30, 172.00) cm. Collected but not retained."
    ),
    WBC = list(
      description = "White blood cell count",
      units       = "10^9/L",
      type        = "continuous",
      notes       = "Listed in Sun 2025 Section 2.1 among the collected physiological and biochemical parameters. No summary statistics are reported in Table 1 and it was not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Listed in Sun 2025 Section 2.1 among the collected parameters. No summary statistics reported; not retained."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Listed in Sun 2025 Section 2.1 among the collected parameters. No summary statistics reported; not retained."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 86L,
    n_studies      = 1L,
    age_range      = "Median (IQR) 62.00 (53.00, 71.25) years; adults aged 18 years and over (Sun 2025 Table 1 and inclusion criteria)",
    weight_range   = "Median (IQR) 62.00 (51.88, 70.00) kg; by sex, males 66.0 (54.0, 73.0) and females 57.0 (47.0, 63.5) (Sun 2025 Table 1 and Discussion)",
    sex_female_pct = 40.7,
    race_ethnicity = "Not reported (single-centre Chinese ICU cohort, presumed predominantly Han Chinese)",
    disease_state  = "Adults with sepsis by Sepsis 3.0 criteria admitted to the intensive care unit with confirmed or suspected Gram-positive infection, treated with teicoplanin for at least 4 days. 20/86 (23.26%) received CRRT (14 CVVH only, 5 CVVHD only, 1 both); 4/86 (4.65%) received ECMO. Children, pregnant women, and patients with joint, bone or endocardial infection were excluded.",
    dose_range     = "Per-protocol study regimen: teicoplanin 400 mg intravenously q12h for the first three doses, then 400 mg q24h maintenance, each given as a 1-hour infusion (Sun 2025 Section 2.2). Actual loading doses received ranged 200-800 mg (median 400) and maintenance doses 200-1,000 mg (median 400). The Monte Carlo dose-optimization simulations explored loading doses of 600-1,200 mg q12h for 3 or 5 doses with 200-1,000 mg q24h maintenance, plus continuous regimens of 400-1,000 mg q12h or 1,000-1,800 mg q24h.",
    regions        = "China (Beijing Jishuitan Hospital Guizhou Hospital, Guiyang, Guizhou; single-centre ICU)",
    renal_function = "Serum creatinine median (IQR) 109.00 (74.00, 184.80) umol/L (Table 1 header mislabels the unit as mg/dL). 20/86 subjects required CRRT for sepsis-associated acute kidney injury. CRRT effluent flow rate, filter adsorption capacity, dialysis timing relative to dosing, and modality were NOT captured in the retrospective dataset, which the authors list as a source of unexplained residual variability.",
    co_medication  = "Not reported beyond the study drug.",
    notes          = "Retrospective single-centre study, 1 June 2022 to 1 June 2024; IRB No. KT2022102101. IMPORTANT for interpreting the variance estimates: only 86 teicoplanin concentrations were available from 86 patients, i.e. essentially ONE trough sample per subject, drawn within 30 minutes before a dose at steady state (Sections 2.2 and 3.1). Observed trough concentrations had median (IQR) 13.40 (10.48, 19.83) mg/L. With a single observation per subject the residual error and the interindividual variability are only weakly separable, which is the likely explanation for the unusually small additive residual (0.23 mg/L) and for the wide bootstrap CI on omega^2 for V, which includes zero. Assay: HPLC-UV at 220 nm after protein precipitation, piperacillin internal standard, calibration range 5.63-125.00 mg/L (r^2 > 0.99), accuracy 2.98-10.36% and precision 7.33-11.25% at QC concentrations of 7.81, 31.25 and 90.00 mg/L. Estimation: Phoenix NLME 8.1, FOCE with interaction. Model evaluation: goodness-of-fit plots, 1,000-replicate bootstrap, and prediction-corrected VPC."
  )

  ini({
    # Structural fixed-effect parameters from Sun 2025 Table 2 ("Final model
    # (n = 86)" Estimate column) and confirmed by the printed final-model
    # Equations 4 and 5.
    lcl       <- log(0.98);    label("Clearance (L/h)")                      # Sun 2025 Table 2: CL = 0.98 L/h (RSE 6.92%; bootstrap median 0.97, 95% CI 0.83-1.12). Equation 4: CL = 0.98 * EXP(eta)
    lvc       <- log(108.69);  label("Central volume of distribution (L)")   # Sun 2025 Table 2: V = 108.69 L (RSE 9.89%; bootstrap median 107.89, 95% CI 84.30-151.56). Equation 5 reference subject is male without CRRT

    # Categorical covariate effects on V, entered on the log scale. Sun 2025
    # Equation 3 gives the categorical covariate form as
    #   tvP' = tv(P) * EXP(theta_cov * Cov_cat),
    # and Equation 5 instantiates it as
    #   V(L) = 108.69 * EXP(-0.71 * (if with CRRT)) * EXP(-1.07 * (if is Female)).
    # Both coefficients are therefore log-scale multipliers, NOT fractional
    # changes, and both are negative (each covariate reduces V).
    e_crrt_vc <- -0.71;        label("Log-scale CRRT effect on V (unitless)")   # Sun 2025 Table 2: theta_CRRT,V = -0.71 (RSE 22.55%; bootstrap median -0.70, 95% CI -1.15 to -0.19). exp(-0.71) = 0.492, a 50.8% reduction. Cross-check: the Results state V was "102.48% lower" with CRRT, and exp(+0.7055) - 1 = 102.48%, reproducing the printed ratio from the unrounded coefficient
    e_sexf_vc <- -1.07;        label("Log-scale female-sex effect on V (unitless)") # Sun 2025 Table 2: theta_Sex,V = -1.07 (RSE 28.53%; bootstrap median -1.07, 95% CI -1.61 to -0.25). exp(-1.07) = 0.343, a 65.7% reduction. Cross-check: the Results state males had "1.90-fold higher" V, and exp(1.07) - 1 = 1.92, i.e. a 1.9-fold increment over the female value

    # Between-subject variability. Sun 2025 Table 2 reports these under the
    # heading "Between-subject variation" as omega^2 CL and omega^2 V, and
    # Equations 4 and 5 print the SAME two numbers inside the exponential IIV
    # term ("* EXP(0.31)" and "* EXP(0.09)"). Read literally the printed
    # equations would be fixed multipliers, which is nonsense; read as the
    # authors intended they pin the table column to the VARIANCE slot of an
    # exponential IIV term, so the values below go into ini() as-is. Two
    # independent cross-checks confirm the variance reading:
    #   (1) the omega^2 CL row's RSE of 17.58% sits just above sqrt(2/86) =
    #       15.25%, the Cramer-Rao floor for a variance, and far above the
    #       sqrt(1/(2*86)) = 7.62% floor that an SD-scale parameter would have;
    #   (2) Section 2.4's generic Equation 1 prints the IIV as ADDITIVE,
    #       Pj = tv(P) + eta_j, but that reading is untenable against these
    #       magnitudes -- an additive eta with variance 0.09 on a V of 108.69 L
    #       is an SD of 0.30 L (0.3% of V), which no modeller would retain with
    #       an estimated RSE of 32%, while an additive eta with variance 0.31 on
    #       a CL of 0.98 L/h would make 3.9% of subjects have negative clearance.
    # Equations 4 and 5 are the final-model equations and are specific, so they
    # govern over the generic Equation 1; this is also Phoenix NLME's default
    # structural-parameter parameterization. See the vignette Errata.
    etalcl ~ 0.31  # Sun 2025 Table 2: omega^2 = 0.31 on CL (RSE 17.58%; bootstrap median 0.31, 95% CI 0.20-0.41). Equivalent to sqrt(exp(0.31) - 1) = 60.2% CV on the linear scale
    etalvc ~ 0.09  # Sun 2025 Table 2: omega^2 = 0.09 on V (RSE 32.11%; bootstrap median 0.08, 95% CI -0.05 to 0.21 -- the bootstrap CI includes zero, so this variance is only weakly identified from single-trough data). Equivalent to sqrt(exp(0.09) - 1) = 30.7% CV

    # Residual variability. Sun 2025 Section 3.2 states the final model
    # incorporated "an additive residual variability", and Table 2 reports it
    # under "Within-subject variation" as sigma_additive with explicit units of
    # mg/L, which fixes it on the SD scale in concentration units.
    addSd  <- 0.23;  label("Additive residual error (mg/L)")  # Sun 2025 Table 2: sigma_additive = 0.23 mg/L (RSE 10.93%; bootstrap median 0.22, 95% CI 0.03-0.26). Small relative to the observed troughs (median 13.40 mg/L) because each subject contributed only one concentration, so the etas absorb nearly all of the variability
  })
  model({
    # Individual PK parameters. Clearance carries no covariate (Sun 2025
    # Section 3.2: "No covariate was found to significantly influence the CL of
    # teicoplanin"). The two categorical effects on V are additive on the log
    # scale, which is algebraically identical to the product of exponentials
    # printed in Equation 5.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + e_crrt_vc * RRT_CRRT_STATUS + e_sexf_vc * SEXF + etalvc)

    kel <- cl / vc

    # One-compartment model with first-order elimination (Sun 2025 Section 3.2).
    # Teicoplanin was given as a 1-hour intravenous infusion (Section 2.2), so
    # doses enter `central` with a rate or duration set in the event table.
    # Dose in mg and vc in L give central / vc in mg/L.
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
