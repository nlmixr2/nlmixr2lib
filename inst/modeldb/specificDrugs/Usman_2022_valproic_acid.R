Usman_2022_valproic_acid <- function() {
  description <- "One-compartment population PK model with first-order elimination for intravenous valproic acid in adult Pakistani and South Korean patients, fitted to pooled routine therapeutic-drug-monitoring data (Usman 2022; 191 patients, 553 serum concentrations). Clearance carries a linear body-weight effect centred on the 67 kg pooled median and a multiplicative linear effect of the Pakistani centre (South Korean patients are the reference); volume of distribution carries a linear body-weight effect centred on 67 kg. Exponential IIV on CL and V, proportional residual error. Fitted in NONMEM 7.4.4 (ADVAN1 TRANS2, FOCE-I) with the covariate model built by PsN stepwise covariate modelling."
  reference <- "Usman M, Shaukat Q-u-A, Khokhar MI, Bilal R, Khan RR, Saeed HA, Ali M, Khan HM. Comparative pharmacokinetics of valproic acid among Pakistani and South Korean patients: A population pharmacokinetic study. PLoS One. 2022;17(8):e0272622. doi:10.1371/journal.pone.0272622. PMCID PMC9401156. Parameters from Table 2 and Eqs 1-4; covariate-model building from S1 Table (PsN scm output)."
  vignette <- "Usman_2022_valproic_acid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "valproic acid", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters CL and V as linear terms centred on 67 kg, the pooled-cohort median (Usman 2022 Eqs 3-4 and Table 1, 'Weight (kg)' row: pooled 67 (40-101) kg; Pakistani 71 (47-101); South Korean 60 (40-91)). The linear form (1 + slope * (WT - 67)) is the PsN scm 'state 2' linear relation (S1 Table rows 'CLWT-2' and 'VWT-2'). It turns negative below WT = 67 - 1/0.0143 = -3 kg for CL and 67 - 1/0.009 = -44 kg for V, so it is positive for every physiological weight.",
      source_name = "WT"
    ),
    RACE_KOREAN = list(
      description = "South Korean patient indicator; 1 = South Korean (Park 2002 dataset), 0 = Pakistani (Aziz Fatima Hospital Faisalabad dataset)",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (South Korean). The paper's covariate CENT is 0 for South Korean and 1 for Pakistani patients (text after Eq 2), so CENT = 1 - RACE_KOREAN and the effect enters on (1 - RACE_KOREAN).",
      notes = "The paper calls CENT both 'center' and 'ethnicity': in this pooled dataset the two coincide, because every patient from the Korean centre is Korean and every patient from the Pakistani centre is Pakistani. The fitted cohort holds only these two groups, so for any patient who is neither Korean nor Pakistani the model is an extrapolation. RACE_KOREAN = 0 then applies the Pakistani clearance. 99 South Korean and 92 Pakistani patients (Table 1). South Korean is the reference because PsN scm sets the most common category as the reference (Eq 1, 'IF(CENT.EQ.0) ... Most common'). This follows the inverted-indicator precedent of Pohl_2022_linzagolix.R and Yao_2018_guselkumab.R.",
      source_name = "CENT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested on CL and V in the scm. V-AGE was added in the fifth forward step (dOFV -6.19) and removed in the backward step (dOFV +6.19 < 6.63, p = 0.0128; S1 Table). It is not in the final model. The prose under Eq 4 still names age, but Eq 4 itself and Table 2 carry no age term."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Listed among the candidate covariates in the Covariate analysis methods. It does not appear in the scm output (S1 Table) and was not retained."
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      notes = "Listed among the candidate covariates in the Covariate analysis methods. It does not appear in the scm output (S1 Table) and was not retained."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Listed among the candidate covariates in the Covariate analysis methods. It does not appear in the scm output (S1 Table) and was not retained."
    ),
    SEXF = list(
      description = "Sex indicator; 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      notes = "Gender is listed as a candidate categorical covariate in the Covariate analysis methods. It does not appear in the scm output (S1 Table) and was not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 191,
    n_studies = 2,
    n_observations = 553,
    age_range = "18-90 years (pooled; Pakistani 19-90, South Korean 18-81)",
    age_median = "48 years (pooled; Pakistani 54, South Korean 44)",
    weight_range = "40-101 kg (pooled; Pakistani 47-101, South Korean 40-91)",
    weight_median = "67 kg (pooled; Pakistani 71, South Korean 60)",
    height_range = "144-190 cm (pooled median 167 cm)",
    bmi_range = "15.6-37.6 kg/m^2 (pooled median 23.2 kg/m^2)",
    sex_female_pct = 33.5,
    race_ethnicity = c(Pakistani = 48.2, Korean = 51.8),
    disease_state = "Adult patients receiving intravenous valproic acid under routine therapeutic drug monitoring. The paper does not state the indication or the clinical setting.",
    dose_range = "Intravenous valproic acid. Single dose median 1000 mg (range 500-1800 mg) in the pooled data (Table 1). The Methods give the Pakistani daily dose as 500-1600 mg. The infusion duration and dosing interval are not reported.",
    regions = "Pakistan (Aziz Fatima Hospital, Faisalabad; 92 patients, 218 samples) and South Korea (99 patients, 335 samples from Park 2002, doi:10.1046/j.1365-2710.2002.00440.x, supplied by the corresponding author of that study).",
    notes = "Retrospective TDM data. Serum valproic acid was assayed by ELISA. Samples were drawn at peak and trough. Observed concentrations ranged from 3.38 to 106.4 mg/L (Table 1). Demographics from Table 1."
  )

  ini({
    # Structural parameters: final-model estimates from Usman 2022 Table 2.
    # They are the typical values for a South Korean patient (RACE_KOREAN = 1)
    # weighing 67 kg.
    lcl <- log(0.931); label("Clearance for a 67 kg South Korean patient (L/h)") # Table 2: CL = 0.931 L/h (RSE 5%)
    lvc <- log(16.6); label("Volume of distribution for a 67 kg patient (L)") # Table 2: Vd = 16.6 L (RSE 2%)

    # Covariate effects (PsN scm linear relations, S1 Table 'CLCENT-2',
    # 'CLWT-2', 'VWT-2').
    e_nonkorean_cl <- 0.386; label("Fractional change in CL for Pakistani vs South Korean patients (unitless)") # Table 2: CL-CENT = 0.386 (RSE 28%); Eqs 1-2
    e_wt_cl <- 0.0143; label("Linear effect of body weight on CL, centred at 67 kg (1/kg)") # Table 2: CL-WT = 0.0143 (RSE 17%); Eq 3
    e_wt_vc <- 0.009; label("Linear effect of body weight on V, centred at 67 kg (1/kg)") # Table 2: Vd-WT = 0.009 (RSE 19%); Eq 4

    # IIV. Table 2 gives IIV as %. Its RSE is on the variance scale: the
    # bootstrap half-width / (1.96 * estimate * RSE) is 0.51 for CL and 0.50
    # for V, which is the 0.5 expected when % = 100 * sqrt(omega). So
    # omega = (IIV/100)^2: 0.434^2 = 0.188356 and 0.223^2 = 0.049729.
    etalcl ~ 0.188356 # Table 2: IIV CL = 43.4% (RSE 13%)
    etalvc ~ 0.049729 # Table 2: IIV Vd = 22.3% (RSE 16%)

    # Proportional residual error. The Table 2 rows follow NONMEM THETA order
    # (CL, Vd, Proportional Error, then the scm covariate thetas), so the error
    # was estimated as a THETA scaling EPS with SIGMA fixed at 1: 0.148 is the SD.
    propSd <- 0.148; label("Proportional residual error (fraction)") # Table 2: Proportional Error = 0.148 (RSE 7%)
  })

  model({
    # Eqs 1-3. PsN scm writes the covariate factors as separate multiplicative
    # terms, TVCL = THETA(CL) * CLCENT * CLWT, with CLCENT = 1 for CENT = 0 and
    # (1 + THETA) for CENT = 1. CENT = 1 is Pakistani = 1 - RACE_KOREAN. The
    # paper's Eq 2 prints the Pakistani CL as 0.931 + 0.386 = 1.317 L/h; the
    # fitted scm code gives 0.931 * 1.386 = 1.290 L/h (see vignette Errata).
    cl <- exp(lcl + etalcl) *
      (1 + e_nonkorean_cl * (1 - RACE_KOREAN)) *
      (1 + e_wt_cl * (WT - 67))

    # Eq 4. The paper writes the V random effect as eta1; Table 2 reports a
    # separate IIV on Vd, so it is its own eta here.
    vc <- exp(lvc + etalvc) * (1 + e_wt_vc * (WT - 67))

    kel <- cl / vc

    # One compartment, intravenous input, first-order elimination (ADVAN1 TRANS2).
    d/dt(central) <- -kel * central

    # Dose in mg, V in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
