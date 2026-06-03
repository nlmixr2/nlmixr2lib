Zhao_2014_vancomycin <- function() {
  description <- "One-compartment IV-infusion population PK model for vancomycin in 70 children with malignant hematological disease (Zhao 2014). Clearance scales with body weight by power exponent (reference 20.2 kg, exponent 0.677) and with Schwartz-formula creatinine clearance by power exponent (reference 191 mL/min/1.73 m^2, exponent 1.03); central volume scales with body weight by power exponent (reference 20.2 kg, exponent 0.838). Vancomycin clearance was substantially higher than in pediatric populations without cancer; the published patient-tailored daily dose is target AUC * CL_i."
  reference <- "Zhao W, Zhang D, Fakhoury M, Fahd M, Duquesne F, Storme T, Baruchel A, Jacqz-Aigrain E. Population pharmacokinetics and dosing optimization of vancomycin in children with malignant hematological disease. Antimicrob Agents Chemother. 2014;58(6):3191-3199. doi:10.1128/AAC.02564-13"
  vignette <- "Zhao_2014_vancomycin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Zhao 2014 Table 1: mean 25.7 kg (SD 15.5), median 20.2 kg (range 5.6-71.0). Reference value 20.2 kg (population median) used in the power-scaling terms in Table 3 footnote for both V (V = theta1 * (WT/20.2)^theta2) and CL (CL = theta3 * (WT/20.2)^theta4 * RF).",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Schwartz-formula creatinine clearance (BSA-normalized eGFR)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Zhao 2014 Table 1: mean 199.8 (SD 63.3), median 191, range 48.7-457 (units reported as mL/min in Table 1 but produced by the Schwartz formula, which yields eGFR in mL/min/1.73 m^2 by construction; Zhao 2014 Methods footnote b 'Creatinine clearance was calculated by the Schwartz formula'). Reference value 191 mL/min/1.73 m^2 (population median) used in the renal-function factor RF = (CLCR/191)^theta5 multiplying CL (Table 3 footnote and Results equation 7). Stored under canonical CRCL per inst/references/covariate-columns.md, which accepts the Schwartz-formula eGFR with its native BSA-normalized units.",
      source_name        = "CLCR"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 70L,
    n_studies        = 1L,
    age_range        = "0.3-17.7 years",
    age_median       = "5.6 years (mean 6.8, SD 4.8)",
    weight_range     = "5.6-71.0 kg",
    weight_median    = "20.2 kg (mean 25.7, SD 15.5)",
    sex_female_pct   = 41.4,
    race_ethnicity   = "Not reported (single-center cohort at Robert Debre Hospital, Paris)",
    disease_state    = "Children with malignant hematological disease (acute lymphoblastic leukemia 40/70, acute myeloblastic leukemia 17/70, juvenile myelomonocytic leukemia 5/70, lymphoma 5/70, other 3/70); 25/70 had bone marrow transplantation. Vancomycin given as empirical antibiotic therapy with therapeutic drug monitoring.",
    dose_range       = "Vancomycin IV infusion over 60 min, empirical initial dose 40-60 mg/kg/day in four divided doses; observed individual doses 50-950 mg (median 250; mean 13.0 mg/kg, range 7.0-31.5 mg/kg). Monitoring target was steady-state trough concentration 10-20 mg/L.",
    regions          = "France (single center, Department of Pediatric Hematology-Oncology, Robert Debre Hospital, APHP, Paris)",
    renal_function   = "Schwartz-formula creatinine clearance median 191 mL/min/1.73 m^2 (range 48.7-457); serum creatinine median 30 umol/L (range 10-141).",
    n_concentrations = 98L,
    notes            = "Baseline demographics from Zhao 2014 Table 1 (cohort enrolled 2010-2011). 98 vancomycin concentrations from 70 children analyzed; blood samples drawn at median 54 h after initiation of treatment, concentrations 1.8-27.3 mg/L. Vancomycin assay: fluorescence polarization immunoassay on Cobas Integra 400 plus (calibration 0.74-80 mg/L, LLOQ 0.74 mg/L). Model fit using NONMEM 7.2.0 FOCE-INTER. Final-model parameter values were confirmed by 500-replicate nonparametric bootstrap (Table 3) and externally validated in an independent group of 20 children (Bayesian estimation, r^2 = 0.99, mean PE 1.0%, mean APE 4.7%). The paper observes that vancomycin CL in this cohort (mean 0.22 L/h/kg) was substantially higher than in pediatric populations without cancer (Discussion Table 4)."
  )

  ini({
    # Structural parameters (Zhao 2014 Table 3 final-model "PK parameter value"
    # column). Reference subject: WT = 20.2 kg, CRCL = 191 mL/min/1.73 m^2.
    lvc <- log(119);  label("Central volume at WT = 20.2 kg (L)")                                       # Zhao 2014 Table 3: theta1 = 119 L (RSE 13.4%)
    lcl <- log(4.37); label("Clearance at WT = 20.2 kg, CRCL = 191 mL/min/1.73 m^2 (L/h)")              # Zhao 2014 Table 3: theta3 = 4.37 L/h (RSE 4.8%)

    # Covariate effects (Zhao 2014 Table 3 footnote and Results Eq. 5-7):
    #   V  = theta1 * (WT / 20.2)^theta2
    #   CL = theta3 * (WT / 20.2)^theta4 * (CRCL / 191)^theta5
    e_wt_vc   <- 0.838; label("Power exponent on (WT/20.2) for V")                                      # Zhao 2014 Table 3: theta2 = 0.838 (RSE 25.1%)
    e_wt_cl   <- 0.677; label("Power exponent on (WT/20.2) for CL")                                     # Zhao 2014 Table 3: theta4 = 0.677 (RSE 12.2%)
    e_crcl_cl <- 1.03;  label("Power exponent on (CRCL/191) for CL")                                    # Zhao 2014 Table 3: theta5 = 1.03  (RSE 21.2%)

    # Inter-individual variability (Zhao 2014 Table 3 "Interindividual variability
    # (%)" rows, exponential model theta_i = theta_TV * exp(eta_i)); for
    # log-normal etas omega^2 = log(CV^2 + 1).
    etalvc ~ 0.46586  # log(0.770^2 + 1); 77.0% CV on V  (RSE 35.9%)
    etalcl ~ 0.11437  # log(0.348^2 + 1); 34.8% CV on CL (RSE 22.9%)

    # Combined proportional + additive residual error (Zhao 2014 Table 3).
    propSd <- 0.053; label("Proportional residual error (fraction)")                                    # Zhao 2014 Table 3: proportional = 5.3% (RSE 79.2%)
    addSd  <- 1.17;  label("Additive residual error (mg/L)")                                            # Zhao 2014 Table 3: additive    = 1.17 mg/L (RSE 26.6%)
  })
  model({
    # Individual PK parameters. V scales by (WT/20.2)^0.838; CL scales by
    # (WT/20.2)^0.677 and by the renal-function factor (CRCL/191)^1.03.
    vc <- exp(lvc + etalvc) * (WT / 20.2)^e_wt_vc
    cl <- exp(lcl + etalcl) * (WT / 20.2)^e_wt_cl * (CRCL / 191)^e_crcl_cl

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, vc in L -> central/vc has units mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
