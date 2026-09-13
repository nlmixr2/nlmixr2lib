Helfer_2026_pentobarbital <- function() {
  description <- "Two-compartment population PK model for intravenous pentobarbital in children (birth to 21 years) given the drug as standard of care for preoperative sedation, deep sedation during mechanical ventilation, and seizure control. All four disposition parameters are allometrically scaled to total body weight against a 70 kg reference with exponents fixed at 0.75 for clearance and intercompartmental clearance and 1 for both volumes; no other covariate survived backward elimination. Inter-individual variability uses a shared-eta construction in which the single estimated eta on clearance is reused on the central volume after multiplication by an estimated scaling factor of 1.13, so the two random effects are perfectly correlated. Residual error is proportional."
  reference <- paste(
    "Helfer VE, Medina-Aymerich L, Muller WJ, Meyer M, Al-Uzri A, McCulloh R,",
    "Hornik CD, Balevic SJ, Greenberg RG, Benjamin DK Jr, Anderson SG, Gonzalez D;",
    "on behalf of the Best Pharmaceuticals for Children Act-Pediatric Trials",
    "Network Steering Committee.",
    "Population Pharmacokinetics and Dosing Simulations of Pentobarbital in the",
    "Pediatric Population.",
    "J Clin Pharmacol. 2026;66(5):e70204.",
    "doi:10.1002/jcph.70204.",
    sep = " "
  )
  vignette <- "Helfer_2026_pentobarbital"

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Two ODE states. Pentobarbital was administered intravenously only (IV bolus
  # and/or continuous infusion), so both states hold parent drug in mg and
  # Cc = central/vc is in mg/L. Bioanalytical Methods confirms the matrix:
  # blood was "collected in EDTAK2 Microtainers and processed into plasma",
  # and "plasma pentobarbital concentrations were quantified using a validated
  # HPLC/MS-MS assay" over 50-50,000 ng/mL (0.05-50 mg/L).
  compartmentData <- list(
    central     = list(analyte = "pentobarbital", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "pentobarbital", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "The only covariate in the final model. Allometric size descriptor normalised to a 70 kg reference (Equation 4: PAR_ij = theta_Pop,j * (WT_i / W_standard)^beta, with W_standard = 70 kg for total body weight). Exponents are fixed at 0.75 for CL and Q and 1 for V1 and V2. Table S1 records that estimating the exponents instead gave 0.758 for clearance - essentially the theoretical value - and about 1.3 for the volumes, but improved the fit by only dOFV = -2.427, so the fixed exponents were retained. Total body weight beat lean body mass (dOFV -52.642) and fat-free mass (dOFV -50.569) as the size descriptor, reaching the lowest OFV of 1195.373 (dOFV -52.984 versus the unscaled base model); the bodyweight-dependent allometric exponent of Wang 2012 was also tried but the covariance step failed and no RSEs were computed. Cohort weights span 3.14-65.0 kg with an overall median of 17.1 kg (Table 1), so the 70 kg reference itself sits above every subject in the analysis dataset. Body weight was recorded at or closest to the time of first sample collection, i.e. treated as a baseline value rather than a time-varying one.",
      source_name        = "WT"
    )
  )

  # Covariates the paper screened but did not retain in the final model. None of
  # these is referenced in model(); they are recorded so the provenance of the
  # covariate screen is not lost. Every dOFV below is from Table S1 (forward
  # inclusion, relative to the total-body-weight base model at OFV 1195.373);
  # inclusion required dOFV <= -3.84 and backward deletion required a rise of
  # at least 6.64.
  covariatesDataExcluded <- list(
    BMI = list(
      description = "Body mass index, and the derived paediatric obesity indicator",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened on both CL (dOFV -0.278) and V1 (dOFV -1.211) as a linear-deviation effect centred at the cohort median of 17.11 kg/m^2; neither reached the 3.84 forward-inclusion threshold. The binary obesity indicator derived from it - BMI at or above the 95th percentile for age and sex on the CDC 2000 growth charts, defined only for participants over 2 years, with missing treated as non-obese - was the one covariate that DID pass forward inclusion, on V1 only (dOFV -4.088), giving V1 = 94.9 L/70 kg in children with obesity versus 33.2 L/70 kg without. It was not retained during backward elimination and is absent from the final model. Only 5 of 39 participants (12.8%, or 18.5% of those aged 2 years and over) were classified as obese. The Discussion links this to the estimated volume exponent of about 1.3 and to pentobarbital's lipophilicity (log P 2.1)."
    ),
    PNA = list(
      description = "Postnatal age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened on CL as a sigmoidal Emax maturation function, both with an estimated Hill coefficient (dOFV -2.53) and with the Hill coefficient fixed to 1 (dOFV -2.517); neither reached the forward-inclusion threshold. Figure S1B shows no trend of the shared eta against postnatal age. Results states that lower absolute CL and V1 in children under 6 years disappeared after adjusting for weight, 'demonstrating that age had no additional effect beyond weight'. Cohort ages span 2 days to 20.8 years, median 4.18 years (Table 1)."
    ),
    SEXF = list(
      description = "Sex, female",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a proportional shift relative to male on CL (dOFV -0.179) and V1 (dOFV -1.283); neither was significant. 17 of 39 participants (43.6%) were female (Table 1)."
    ),
    RACE_BLACK = list(
      description = "Black or African American race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Part of the race screen, tested in three encodings, all relative to White (missing coded as White): a three-level shift for Black or African American / Asian / Multiple races on CL (dOFV -1.516) and V1 (dOFV -0.266); a single non-White indicator on CL (dOFV -0.348) and V1 (dOFV -0.006); and a two-level Black-or-African-American / other-non-White shift on CL (dOFV -1.44) and V1 (dOFV -0.008). No encoding reached the forward-inclusion threshold. 7 of 39 participants (17.9%) were Black or African American (Table 1)."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested as one level of the three-level race shift on CL and V1 described under RACE_BLACK; not significant. 2 of 39 participants (5.1%) were Asian (Table 1)."
    ),
    RACE_HISPANIC = list(
      description = "Hispanic or Latino ethnicity indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a proportional shift relative to not-Hispanic-or-Latino (missing coded as not Hispanic) on CL (dOFV -1.625) and V1 (dOFV -0.064); neither was significant. 9 of 39 participants (23.1%) were Hispanic or Latino (Table 1)."
    ),
    ECMO_STATUS = list(
      description = "Receiving extracorporeal membrane oxygenation at the time of sampling",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL (dOFV -0.408) and V1 (dOFV -0.04); neither was significant. ECMO was nominally an exclusion criterion but 4 participants (10.3%) were on ECMO at data collection and were retained because they contributed 5 PK samples. The Discussion cautions that this sample size limits any conclusion that ECMO has no effect."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 39L,
    n_studies      = 1L,
    age_range      = "0.01-20.8 years (youngest 2 days; all other participants at least 6 months)",
    age_median     = "4.18 years",
    weight_range   = "3.14-65.0 kg",
    weight_median  = "17.1 kg",
    sex_female_pct = 43.6,
    race_ethnicity = c(White = 74.4, Black = 17.9, Asian = 5.1, Multiple = 2.6),
    disease_state  = "Children under 21 years receiving intravenous pentobarbital as part of standard of care. Sedation was the most common indication, followed by seizure control. 4 participants (10.3%) received vasopressors and 1 (2.6%) received valproic acid on the day of PK sampling; 4 (10.3%) were on ECMO. 5 participants (12.8%) met the paediatric obesity definition.",
    dose_range     = "Intravenous bolus and/or continuous infusion. Median IV bolus dose 2 mg/kg (range 0.39-8.9), absolute 4.7-290 mg (median 45 mg); median IV infusion rate 2.4 mg/kg/h (range 0.00175-16.1). Median 5 doses per participant (range 1-45), median bolus interval 2.34 h (range 0-73.7), median infusion duration 7.2 h (range 0.3-146.42).",
    regions        = "United States (multicentre; Pediatric Trials Network sites)",
    notes          = "Opportunistically collected standard-of-care data from the POP01 study, 'Pharmacokinetics of Understudied Drugs Administered to Children Per Standard of Care' (NICHD-2011-POP01, ClinicalTrials.gov NCT01431326). 42 participants were enrolled; 3 were excluded (one who also received a single oral dose, one whose only sample was below the limit of quantification, and one with a negative infusion duration from a dosing-entry error), leaving 39 participants and 70 plasma samples (median 2 samples per participant, range 1-5). The data are sparse, which the authors identify as the main limitation on covariate detection. Metabolic panel values (direct and total bilirubin, serum creatinine, AST, ALT, albumin) were missing for more than 48.7% of participants and so were never evaluated - covariates with more than 10% missingness were excluded from testing by protocol. C-reactive protein, a significant covariate on CL in Ketharanathan 2023, was not collected."
  )

  ini({
    # Structural parameters - Table 2, "Estimate" column. Every value is the
    # typical value at the 70 kg allometric reference, as the row labels state.
    lcl <- log(5.21); label("Clearance at the 70 kg reference (L/h)")                              # Table 2: CL = 5.21 L/h/70 kg (RSE 13.2%); bootstrap median 5.13, 95% CI 3.80-6.36
    lvc <- log(37.4); label("Central volume of distribution at the 70 kg reference (L)")           # Table 2: V1 = 37.4 L/70 kg (RSE 19.3%); bootstrap median 35.43, 95% CI 24.43-53.79
    lq  <- log(18.1); label("Intercompartmental clearance at the 70 kg reference (L/h)")           # Table 2: Q = 18.1 L/h/70 kg (RSE 20.2%); bootstrap median 17.81, 95% CI 13.15-53.34
    lvp <- log(63.9); label("Peripheral volume of distribution at the 70 kg reference (L)")        # Table 2: V2 = 63.9 L/70 kg (RSE 14.2%); bootstrap median 62.71, 95% CI 44.39-109.65

    # Allometric exponents - fixed, not estimated, and shared within each pair.
    # Equation 4 and the Methods text: beta "was either estimated or fixed to
    # 0.75 for clearance parameters (CL and intercompartmental CL [Q]) and 1 for
    # volume of distribution parameters (central and peripheral)". Results
    # confirms the fixed values were carried into the final model: "fixed
    # allometric exponents (0.75 for CL and Q; 1 for V1 and V2) were retained in
    # the final model". Table S1's total-body-weight row writes all four
    # exponents out as literals. They carry no RSE in Table 2 because they are
    # not estimated parameters.
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent of body weight shared by CL and Q (unitless)")   # Equation 4 / Results; Table S1 total body weight row: (WT/70)^0.75 on both CL and Q
    e_wt_vc_vp <- fixed(1);    label("Allometric exponent of body weight shared by V1 and V2 (unitless)")  # Equation 4 / Results; Table S1 total body weight row: (WT/70)^1 on both V1 and V2

    # Shared-variability scaler. Methods: IIV on CL and V1 was supported, but
    # "upon estimating the correlation between these random effects, a high
    # degree of correlation was observed. To account for this, a shared
    # variability approach was adopted." Equations 2 and 3 give the construction
    # verbatim:
    #     PAR_ij = theta_Pop,j * exp(eta_ij)                (2)
    #     PAR_ik = theta_Pop,k * exp(theta_var * eta_ij)    (3)
    # with "theta_var denotes the ratio of the standard deviation between
    # parameters j and k, where the variance of PAR_ik can be computed as
    # Var(PAR_ik) = theta_var^2 * omega^2_ij". The Table 2 footnote fixes which
    # parameter is which: "Estimated shared eta_1 approach for CL and V1, with
    # Var(V) = theta^2_shared variability * omega^2_CL" - so CL carries the base
    # eta and V1 is the scaled one. This is the registered vc_eta_scale pattern
    # (Tang_2023_tenecteplase.R, Hirt_2009_efavirenz.R, Prytula_2016_tacrolimus.R);
    # the two random effects are perfectly correlated by construction.
    vc_eta_scale <- 1.13; label("Scaling factor relating the V1 random effect to the CL random effect, correlation fixed to 1 (unitless)")  # Table 2 row 'Shared variability CL and V1' = 1.13 (RSE 36.9%); bootstrap median 1.08, 95% CI 0.15-2.33

    # Inter-individual variability - a single eta, on CL, reused on V1 through
    # vc_eta_scale above. Table 2 tabulates it as a coefficient of variation,
    # not a variance: the table footnote defines "CV(%) = sqrt(exp(eta) - 1) x
    # 100", i.e. CV = sqrt(exp(omega^2) - 1), so omega^2 = log(CV^2 + 1) =
    # log(0.619^2 + 1) = 0.3244.
    #
    # The paper's own arithmetic confirms this scale. Applying the footnote
    # formula to the scaled variance gives the V1 CV that Table 2 quotes:
    # theta_var^2 x omega^2_CL = 1.13^2 x 0.3244 = 0.4142, and
    # sqrt(exp(0.4142) - 1) = 71.6%, matching the footnote's "approximately
    # 71.5% CV for V1". Reading the 61.9 as a variance instead would not
    # reproduce that number.
    etalcl ~ 0.3244  # Table 2 row 'eta 1' = 61.9 %CV (RSE 23.1%) [7% shrinkage]; bootstrap median 57.1%, 95% CI 32.9-106.4%. omega^2 = log(0.619^2 + 1) = 0.3244

    # Residual error - proportional only. Results: "Pentobarbital concentration
    # data were best described by a two-compartment PK model with proportional
    # residual error." Additive and combined models were also explored (Methods)
    # and Table 2 reports a single residual row, in percent.
    propSd <- 0.229; label("Proportional residual error (fraction)")  # Table 2 row 'Proportional error (%)' = 22.9% (RSE 31.4%) [24% shrinkage]; bootstrap median 21.8%, 95% CI 9.7-31.9%
  })

  model({
    # 1. Allometric size scaling on a 70 kg reference (Equation 4). Both
    #    clearance terms share one exponent and both volumes share the other,
    #    exactly as the paper constrains them.
    sizeCl <- (WT / 70)^e_wt_cl_q
    sizeV  <- (WT / 70)^e_wt_vc_vp

    # 2. Individual parameters. Exponential IIV (Equation 1) on CL, and the same
    #    eta scaled by vc_eta_scale on V1 (Equations 2-3). Q and V2 carry no
    #    random effect - Table 2 reports a single eta.
    cl <- exp(lcl + etalcl) * sizeCl
    vc <- exp(lvc + vc_eta_scale * etalcl) * sizeV
    q  <- exp(lq) * sizeCl
    vp <- exp(lvp) * sizeV

    # 3. Micro-constants for the two-compartment system.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system. Intravenous administration only - bolus loading doses and
    #    continuous maintenance infusions both enter `central` directly, so
    #    there is no depot and no bioavailability term.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # 5. Observation and error.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
