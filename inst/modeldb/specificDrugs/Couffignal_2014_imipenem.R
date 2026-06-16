Couffignal_2014_imipenem <- function() {
  description <- paste(
    "Two-compartment IV population PK model for imipenem in 51 critically",
    "ill adult ICU patients with suspected ventilator-associated pneumonia",
    "due to Gram-negative bacilli (Couffignal 2014). All patients received",
    "imipenem as a 0.5 h IV infusion every 8 hours; the protocol dose",
    "(500, 750 or 1000 mg) was chosen by Cockcroft-Gault creatinine",
    "clearance per the European Medicine Agency renal-adjustment table.",
    "Central clearance scales as a power of measured 4-hour creatinine",
    "clearance (reference 86.4 mL/min, the cohort median); central volume",
    "scales jointly with total bodyweight (reference 77 kg) and serum",
    "albumin (reference 18 g/L). The model was fitted in Monolix 4.1.2",
    "using the SAEM algorithm with M3-equivalent BQL handling."
  )
  reference <- paste(
    "Couffignal C, Pajot O, Laouenan C, Burdet C, Foucrier A, Wolff M,",
    "Armand-Lefevre L, Mentre F, Massias L.",
    "Population pharmacokinetics of imipenem in critically ill patients",
    "with suspected ventilator-associated pneumonia and evaluation of",
    "dosage regimens.",
    "Br J Clin Pharmacol. 2014;78(5):1022-1034.",
    "doi:10.1111/bcp.12435.",
    sep = " "
  )
  vignette <- "Couffignal_2014_imipenem"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total bodyweight at inclusion",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Total bodyweight at inclusion (paper Table 1; cohort median 77 kg,",
        "range 45-126 kg). Used as the reference value in the multiplicative",
        "power model V1 = 20.4 * (WT/77)^1.3 * (ALB/18)^(-1.1) (paper Results",
        "paragraph below Table 2). The 77 kg reference is the cohort median",
        "per Methods 'Covariate analysis' (COVmedian is the median value of",
        "covariates). Supplementary Table S1 reports body-weight percentiles",
        "10th = 53 kg, 50th = 77 kg, 90th = 111 kg, which bracket the cohort",
        "covariate sensitivity analysis."
      ),
      source_name        = "Total bodyweight (kg)"
    ),
    CRCL = list(
      description        = "Measured 4-hour urinary creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Measured creatinine clearance over a 4-hour urine collection",
        "starting at the fourth imipenem infusion (paper Methods 'Covariate",
        "analysis'; the paper cites internal references 30-32 establishing",
        "this measurement as a true reflection of renal function during the",
        "fourth infusion). Reported in raw mL/min and NOT BSA-normalised to",
        "mL/min/1.73 m^2; cohort median 86.4 mL/min, range 9.1-571.4 mL/min",
        "(paper Table 1). Used as the reference value in the multiplicative",
        "power model CL = 13.2 * (CRCL/86.4)^0.2 (paper Results paragraph",
        "below Table 2). Supplementary Table S1 reports CrCL4h percentiles",
        "10th = 17 mL/min, 50th = 86.4 mL/min, 90th = 258 mL/min, which",
        "bracket the cohort covariate sensitivity analysis. Stored under the",
        "canonical CRCL column per inst/references/covariate-columns.md",
        "(CRCL accepts raw measured CrCL in mL/min when the source paper",
        "does not BSA-normalise; the per-model description records the",
        "assay form -- precedent: Delattre 2010 amikacin, Lamoth 2009",
        "imipenem, AbdulAziz 2016 doripenem)."
      ),
      source_name        = "CrCL4h"
    ),
    ALB = list(
      description        = "Serum albumin measured at the fourth imipenem dose",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Serum albumin (g/L) measured at the time of the fourth imipenem",
        "dose (paper Table 1; cohort median 18 g/L, range 10-28 g/L, with",
        "the median reported as a footnote derived from 9 patients with",
        "complete data). Missing values for the 8 remaining patients were",
        "imputed to the cohort median 18 g/L before model estimation",
        "(paper Methods 'Covariate analysis' and Results 'Model with",
        "covariates'). Used as the reference value in the multiplicative",
        "power model V1 = 20.4 * (WT/77)^1.3 * (ALB/18)^(-1.1) (paper",
        "Results paragraph below Table 2). Supplementary Table S1 reports",
        "serum-albumin percentiles 10th = 11 g/L, 50th = 18 g/L, 90th =",
        "23 g/L, which bracket the cohort covariate sensitivity analysis."
      ),
      source_name        = "Serum albumin (g/L)"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 51L,
    n_studies        = 1L,
    age_range        = "28-84 years",
    age_median       = "60 years",
    weight_range     = "45-126 kg",
    weight_median    = "77 kg",
    sex_female_pct   = 19.6,
    race_ethnicity   = NULL,
    disease_state    = paste(
      "Critically ill adult ICU patients (age >= 18) mechanically",
      "ventilated for more than 48 h with clinical suspicion of",
      "ventilator-associated pneumonia (VAP) due to Gram-negative",
      "bacilli, at high risk of multiresistant bacteria (at least 6",
      "days of mechanical ventilation or antibiotic treatment within",
      "15 days; paper Methods 'Study design and population'). Septic",
      "shock present in 35% at the fourth dose; SAPS II at admission",
      "median 40 (range 19-74), SOFA score at the fourth dose median",
      "6 (range 2-14), oedema score median 7 (range 0-18). Reasons for",
      "ICU admission: medical 78%, surgical 22%. 94% had received",
      "antibiotic therapy in the 3 months before admission (30% of",
      "those had received imipenem)."
    ),
    dose_range       = paste(
      "0.5 h IV infusion of imipenem every 8 hours; dose chosen per",
      "the European Medicine Agency renal-adjustment table by",
      "Cockcroft-Gault CrCL (ClCG > 70 mL/min/1.73 m^2: 1000 mg q8h;",
      "ClCG > 30 and <= 70 mL/min/1.73 m^2: 750 mg q8h; ClCG <= 30",
      "mL/min/1.73 m^2: 500 mg q8h). Observed cohort distribution: 4",
      "patients (9%) on 500 mg, 15 (29%) on 750 mg, 32 (62%) on 1000",
      "mg, all q8h (paper Results 'Patients')."
    ),
    regions          = "France (multicentre: 3 ICUs across 2 hospitals -- Hopital V Dupouy polyvalent ICU in Argenteuil, AP-HP Hopital Bichat medical and surgical ICUs in Paris)",
    renal_function   = paste(
      "Measured 4-hour creatinine clearance median 86.4 mL/min (range",
      "9.1-571.4 mL/min, paper Table 1; raw mL/min, not BSA-normalised).",
      "Patients with Cockcroft-Gault CrCL < 10 mL/min or on renal",
      "replacement therapy were excluded."
    ),
    n_concentrations = 297L,
    notes            = paste(
      "Prospective open-label multicentre IMPACT study",
      "(ClinicalTrials.gov NCT00950222; paper Methods 'Study design",
      "and population'). 63 patients screened, 12 excluded (3 lacked",
      "a kinetic profile, 9 did not receive 4 doses); 51 included in",
      "the PK analysis. Six plasma imipenem samples per patient drawn",
      "around the fourth dose (trough immediately before the dose,",
      "then 0.5, 1, 2, 5 and 8 h after the start of the 0.5 h",
      "infusion). One patient received the fourth dose 5 h late.",
      "Plasma stabilised within 0.5 h of collection by",
      "4-morpholine propane sulphonic acid (MOPS) in ethylene glycol",
      "and immediately frozen at -80 degrees C. Imipenem measured by",
      "HPLC-UV at 302 nm after ultrafiltration (Interchrome YP5C18",
      "25QS reverse phase column); LLOQ 0.5 mg/L; 9% of",
      "concentrations were BQL. Model fitted with Monolix 4.1.2 by",
      "stochastic-approximation expectation-maximisation (SAEM); BQL",
      "data handled by left-censoring (Monolix M3-equivalent).",
      "Final-model 95% confidence intervals from 1000 nonparametric",
      "bootstrap resamples."
    )
  )

  ini({
    # ===== Structural PK (Couffignal 2014 Table 2 'Final model' column) =====
    # Reference subject: CRCL = 86.4 mL/min, WT = 77 kg, ALB = 18 g/L
    # (cohort medians; paper Methods 'Covariate analysis' centres each
    # covariate on its median).
    lcl <- log(13.2); label("Typical CL at CRCL=86.4 mL/min (L/h)")  # Couffignal 2014 Table 2 final model: CL = 13.2 L/h (RSE 5%; bootstrap median 13.2, 95% CI 11.4-15.3)
    lvc <- log(20.4); label("Typical central volume V1 at WT=77 kg, ALB=18 g/L (L)")  # Couffignal 2014 Table 2 final model: V1 = 20.4 L (RSE 7%; bootstrap median 19.8, 95% CI 14.9-25.4)
    lq  <- log(12.2); label("Typical intercompartmental clearance Q (L/h)")  # Couffignal 2014 Table 2 final model: Q = 12.2 L/h (RSE 25%; bootstrap median 12.3, 95% CI 4.7-20.3)
    lvp <- log(9.8);  label("Typical peripheral volume V2 (L)")  # Couffignal 2014 Table 2 final model: V2 = 9.8 L (RSE 13%; bootstrap median 10.5, 95% CI 6.9-13.7)

    # ===== Covariate effects (estimated power exponents) =====
    # Paper Methods 'Covariate analysis' defines the multiplicative power
    # form: parameter_i = theta_pop * (COV_i / COV_median)^beta. The paper
    # writes "the unit of beta is the logarithm of the unit of the
    # associated parameter"; numerically each beta is the unitless
    # exponent on the (COV_i / COV_median) ratio. All three coefficients
    # are estimated by the SAEM algorithm and retained after backward
    # likelihood-ratio selection at p < 0.05 (paper Table 2 final model).
    e_crcl_cl <- 0.2;  label("Power exponent on (CRCL/86.4) for CL (unitless)")  # Couffignal 2014 Table 2 final model: beta_CrCL4h = 0.2 (RSE 19%; p = 6.4e-5; bootstrap median 0.25, 95% CI 0.1-0.4)
    e_wt_vc   <- 1.3;  label("Power exponent on (WT/77) for V1 (unitless)")  # Couffignal 2014 Table 2 final model: beta_Weight = 1.3 (RSE 17%; p = 1.3e-4; bootstrap median 1.2, 95% CI 0.6-2.2)
    e_alb_vc  <- -1.1; label("Power exponent on (ALB/18) for V1 (unitless)")  # Couffignal 2014 Table 2 final model: beta_Serum_albumin = -1.1 (RSE 18%; p = 1.8e-4; bootstrap median -1.0, 95% CI -1.8 to -0.5)

    # ===== Inter-individual variability (Couffignal 2014 Table 2 final model) =====
    # Monolix 4.1.2 exponential random-effects parameterisation: for
    # subject i, CL_i = CL_pop * exp(eta_CL,i) with eta ~ N(0, omega^2).
    # The "omega CL (%)" column in Monolix's parameter table reports
    # omega (the SD of eta) as a percentage; hence variance = (omega/100)^2.
    # IIV reported only on CL and V1 (paper Results 'Population PK
    # analysis': "Since the variability of intercompartmental clearance
    # (Q) and the volume of distribution of the peripheral compartment
    # (V2) were very low, the between-subject variability was not
    # estimated and was taken as zero."). A correlation between eta_CL
    # and eta_V1 was retained throughout model development.
    # Variance diagonals: omega_CL = 0.38 -> 0.38^2 = 0.1444 ; omega_V1
    # = 0.31 -> 0.31^2 = 0.0961. Off-diagonal covariance: cov(eta_CL,
    # eta_V1) = correlation * omega_CL * omega_V1 = 0.51 * 0.38 * 0.31
    # = 0.0601.
    etalcl + etalvc ~ c(0.1444,
                        0.0601, 0.0961)  # Couffignal 2014 Table 2 final model: omega_CL = 0.38 (RSE 13%; bootstrap median 0.36, 95% CI 0.26-0.49), omega_V1 = 0.31 (RSE 18%; bootstrap median 0.22, 95% CI 0.01-0.45), correlation eta_CL,eta_V1 = 0.51 (RSE 28%; bootstrap median 0.79, 95% CI -1 to 1)

    # ===== Residual error (Couffignal 2014 Table 2 final model) =====
    # Proportional model: y_obs_ij = y_pred_ij * (1 + eps_ij), eps_ij ~
    # N(0, sigma^2). Paper Results 'Population PK analysis': "A
    # proportional model was used to describe the residual variability."
    propSd <- 0.33; label("Proportional residual error (fraction)")  # Couffignal 2014 Table 2 final model: sigma = 0.33 (RSE 3%; bootstrap median 0.34, 95% CI 0.26-0.41)
  })

  model({
    # ----- Individual PK parameters -----
    # Multiplicative power covariate model (paper Results paragraph below
    # Table 2): "CL = 13.2 * (CrCL/86.4)^0.2 and V1 = 20.4 * (Weight/77)^1.3
    # * (Albumin/18)^(-1.1)". Q and V2 have no covariate effects and no
    # IIV.
    cl <- exp(lcl + etalcl) * (CRCL / 86.4)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 77)^e_wt_vc * (ALB / 18)^e_alb_vc
    q  <- exp(lq)
    vp <- exp(lvp)

    # ----- Micro-constants -----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ----- ODE system -----
    # IV imipenem given as a 0.5 h infusion every 8 h. Dose is delivered
    # directly into the central compartment via the event-table rate
    # column; no absorption compartment (paper Methods 'Sampling
    # procedure': "0.5 h infusions").
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # ----- Output -----
    # Plasma imipenem concentration: dose in mg, vc in L -> mg/L (matches
    # the HPLC-UV assay reporting units in paper Methods 'Sampling
    # procedure and analytical methods').
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
