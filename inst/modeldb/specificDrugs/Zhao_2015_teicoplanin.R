Zhao_2015_teicoplanin <- function() {
  description <- "Two-compartment IV-injection population PK model for teicoplanin in 85 children with malignant hematological disease (Zhao 2015). Body weight enters Vc and Vp with the fixed allometric exponent 1 and enters CL and Q with the fixed allometric exponent 0.75; Schwartz-formula creatinine clearance enters CL via a power exponent estimated at 0.606. Reference subject: WT = 27.1 kg, CRCL = 179 mL/min. The published model was used to derive age-band mg/kg dosing (18 mg/kg for infants, 14 mg/kg for children, 12 mg/kg for adolescents) and a patient-tailored daily dose (target AUC * CL_i) to attain the AUC(0,24 h) target of 750 mg.L/h."
  reference <- "Zhao W, Zhang D, Storme T, Baruchel A, Decleves X, Jacqz-Aigrain E. Population pharmacokinetics and dosing optimization of teicoplanin in children with malignant haematological disease. Br J Clin Pharmacol. 2015;80(5):1197-1207. doi:10.1111/bcp.12710"
  vignette <- "Zhao_2015_teicoplanin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Zhao 2015 Table 1: mean 32.3 kg (SD 17.8), median 27.1 kg (range 7.7-90.6). Reference value 27.1 kg (population median) used in the allometric size terms in Table 2 footnote for Vc, Vp, Q, and CL: Vc = theta1 * (WT/27.1)^1, Vp = theta2 * (WT/27.1)^1, Q = theta3 * (WT/27.1)^0.75, and CL = theta4 * (WT/27.1)^0.75 * RF. The two volume exponents (1) and the two clearance exponents (0.75) are fixed a priori per Results: 'allometric coefficients of 0.75 for CL and Q, 1 for V1 and V2'.",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Schwartz-formula creatinine clearance (BSA-normalized eGFR)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Zhao 2015 Table 1: mean 191.2 (SD 76.2), median 178.9, range 48.6-464.1 (units reported as mL/min in Table 1 but produced by the Schwartz formula, which yields eGFR in mL/min/1.73 m^2 by construction; Zhao 2015 Methods 'creatinine clearance (Schwartz formula)'). Reference value 179 mL/min/1.73 m^2 (population median; the Results text uses '179 mL/min' and Table 2 footnote writes 'RF = (CLcr/179)^theta5'). Stored under canonical CRCL per inst/references/covariate-columns.md, which accepts the Schwartz-formula eGFR with its native BSA-normalized units.",
      source_name        = "CLCR"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 85L,
    n_studies        = 1L,
    age_range        = "0.5-16.9 years",
    age_median       = "8.1 years (mean 8.4, SD 4.6)",
    weight_range     = "7.7-90.6 kg",
    weight_median    = "27.1 kg (mean 32.3, SD 17.8)",
    sex_female_pct   = 37.6,
    race_ethnicity   = "Not reported (single-center cohort at Robert Debre Hospital, Paris)",
    disease_state    = "Children with malignant hematological disease (acute lymphoblastic leukemia 41/85, acute myeloblastic leukemia 27/85, biphenotypic acute leukemia 4/85, juvenile myelomonocytic leukemia 3/85, lymphoma 5/85, other 5/85); 34/85 had bone marrow transplantation. Teicoplanin given as empirical antibiotic therapy with therapeutic drug monitoring.",
    dose_range       = "Teicoplanin (Targocid, Sanofi-Aventis) IV injection over 3-5 min. Initial regimen 10 mg/kg every 12 h for three loading doses followed by 10 mg/kg once daily maintenance. Observed individual doses 90-600 mg (median 270; mean 9.4 mg/kg, range 4.9-13.2 mg/kg). Monitoring target was steady-state trough concentration Css,min >= 10 mg/L.",
    regions          = "France (single center, Department of Paediatric Haemato-Oncology, Robert Debre Hospital, APHP, Paris)",
    renal_function   = "Schwartz-formula creatinine clearance median 178.9 mL/min/1.73 m^2 (range 48.6-464.1); serum creatinine median 33 umol/L (range 12-121).",
    n_concentrations = 143L,
    age_groups       = "Infants (1 month to 2 years) 10/85; children (2-12 years) 49/85; adolescents (12-18 years) 26/85.",
    notes            = "Baseline demographics from Zhao 2015 Table 1 (cohort enrolled 2012-2013). 143 teicoplanin concentrations from 85 children analyzed (TDM n=123 + opportunistic n=20), concentrations <LLOQ to 35.1 mg/L. Teicoplanin assay: quantitative microsphere system (QMS) on CDX automate (Thermo Fisher), calibration 0-100 mg/L, LLOQ 3 mg/L (CV<5%). Model fit using NONMEM 7.2.0 FOCE-INTER. Final-model parameter values were confirmed by 500-replicate nonparametric bootstrap (Table 2) and externally validated in an independent group of 15 children (Bayesian estimation, r^2 = 0.99, mean PE 0.7%, mean APE 5.3%)."
  )

  ini({
    # Structural parameters (Zhao 2015 Table 2 final-model "PK parameter value"
    # column). Reference subject: WT = 27.1 kg, CRCL = 179 mL/min/1.73 m^2.
    lvc <- log(12.9);  label("Central volume at WT = 27.1 kg (L)")                                       # Zhao 2015 Table 2: theta1 = 12.9 L  (RSE 24.7%)
    lvp <- log(25.2);  label("Peripheral volume at WT = 27.1 kg (L)")                                    # Zhao 2015 Table 2: theta2 = 25.2 L  (RSE 19.2%)
    lq  <- log(0.341); label("Inter-compartmental clearance at WT = 27.1 kg (L/h)")                      # Zhao 2015 Table 2: theta3 = 0.341 L/h (RSE 25.8%)
    lcl <- log(0.491); label("Clearance at WT = 27.1 kg, CRCL = 179 mL/min/1.73 m^2 (L/h)")              # Zhao 2015 Table 2: theta4 = 0.491 L/h (RSE 10.1%)

    # Covariate effects (Zhao 2015 Table 2 footnote and Results equations):
    #   Vc = theta1 * (WT/27.1)^1
    #   Vp = theta2 * (WT/27.1)^1
    #   Q  = theta3 * (WT/27.1)^0.75
    #   CL = theta4 * (WT/27.1)^0.75 * (CRCL/179)^theta5
    # The two volume exponents (1) and the two clearance exponents (0.75) are
    # fixed a priori (Results: "allometric coefficients of 0.75 for CL and Q,
    # 1 for V1 and V2"). The CRCL exponent is estimated.
    e_wt_vc   <- fixed(1);     label("Power exponent on (WT/27.1) for Vc (fixed)")                       # Zhao 2015 Results: V1 allometric coefficient fixed at 1
    e_wt_vp   <- fixed(1);     label("Power exponent on (WT/27.1) for Vp (fixed)")                       # Zhao 2015 Results: V2 allometric coefficient fixed at 1
    e_wt_q    <- fixed(0.75);  label("Power exponent on (WT/27.1) for Q (fixed)")                        # Zhao 2015 Results: Q  allometric coefficient fixed at 0.75
    e_wt_cl   <- fixed(0.75);  label("Power exponent on (WT/27.1) for CL (fixed)")                       # Zhao 2015 Results: CL allometric coefficient fixed at 0.75
    e_crcl_cl <- 0.606;        label("Power exponent on (CRCL/179) for CL")                              # Zhao 2015 Table 2: theta5 = 0.606 (RSE 25.7%)

    # Inter-individual variability (Zhao 2015 Table 2 "Inter-individual
    # variability (%)" rows, exponential model theta_i = theta_TV * exp(eta_i));
    # for log-normal etas omega^2 = log(CV^2 + 1). The paper does not report
    # IIV on Vp (V2), so no etalvp is included.
    etalvc ~ 0.04811  # log(0.222^2 + 1); 22.2% CV on Vc (RSE 110.5%)
    etalq  ~ 1.00402  # log(1.315^2 + 1); 131.5% CV on Q  (RSE 46.8%)
    etalcl ~ 0.09633  # log(0.318^2 + 1); 31.8% CV on CL (RSE 31.9%)

    # Combined proportional + additive residual error (Zhao 2015 Table 2).
    propSd <- 0.148; label("Proportional residual error (fraction)")                                     # Zhao 2015 Table 2: proportional = 14.8% (RSE 33.8%)
    addSd  <- 1.1;   label("Additive residual error (mg/L)")                                             # Zhao 2015 Table 2: additive     = 1.1 mg/L (RSE 68.3%)
  })
  model({
    # Individual PK parameters. Volumes scale linearly with body weight;
    # clearances scale by (WT/27.1)^0.75; CL is additionally scaled by the
    # renal-function factor (CRCL/179)^0.606.
    vc <- exp(lvc + etalvc) * (WT / 27.1)^e_wt_vc
    vp <- exp(lvp)          * (WT / 27.1)^e_wt_vp
    q  <- exp(lq  + etalq)  * (WT / 27.1)^e_wt_q
    cl <- exp(lcl + etalcl) * (WT / 27.1)^e_wt_cl * (CRCL / 179)^e_crcl_cl

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg, vc in L -> central/vc has units mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
