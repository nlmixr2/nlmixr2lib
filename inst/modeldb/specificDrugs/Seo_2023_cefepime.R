Seo_2023_cefepime <- function() {
  description <- "Two-compartment intravenous population PK model for cefepime in 21 critically ill Korean adults with hospital-acquired or ventilator-associated pneumonia (Seo 2023); total clearance scales as a power function of raw Cockcroft-Gault creatinine clearance referenced to 77.21 mL/min, with exponential inter-individual variability on CL, Vc, Q and Vp and a proportional residual error. The plasma unbound fraction fixed at 81% for the paper's fT>MIC Monte Carlo target-attainment simulations is exposed as fu, with the unbound concentration returned as Cu."
  reference <- paste(
    "Seo H, Kim YK, Park S, Kim HI, Lee DH. Population Pharmacokinetics",
    "and Monte Carlo Simulation of Cefepime in Critically Ill Patients",
    "with Hospital-Acquired/Ventilator-Associated Pneumonia.",
    "Infect Chemother. 2023 Mar;55(1):29-41.",
    "doi:10.3947/ic.2022.0087",
    sep = " "
  )
  vignette <- "Seo_2023_cefepime"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central     = list(analyte = "cefepime", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cefepime", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance estimated by the Cockcroft-Gault equation (raw, NOT BSA-normalized).",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate retained in the Seo 2023 final model, entering",
        "total clearance as the power term (CG / 77.21)^theta_2",
        "(Table 2, 'Structural model' row 'CL = theta_1 x (CG/77.21)^theta_2').",
        "Carried as the RAW Cockcroft-Gault value in mL/min, NOT the canonical",
        "BSA-normalized mL/min/1.73 m^2 form -- the same deviation documented",
        "in Jonckheere_2019_cefepime.R and Delattre_2010_amikacin.R for the",
        "same paper-reported quantity. Seo 2023 Table 1 reports CLCR by",
        "Cockcroft-Gault in mL/min (mean 79.02, SD 47.67; median 77.2, IQR",
        "45.3-125); the model's reference constant 77.21 mL/min is that cohort",
        "median. The paper computed MDRD, CKD-EPI, modified MDRD and modified",
        "CKD-EPI renal-function estimates as well (Methods 'Population PK",
        "analysis'; Table 1), but only the Cockcroft-Gault value was retained",
        "as the covariate on CL. The paper reports no min-max for the",
        "model-development CLCR range (Table 1 gives only mean, SD, median and",
        "IQR), so the fitted range is bounded only by that IQR of 45.3-125",
        "mL/min; the Monte Carlo simulations extrapolate the term over 0-130",
        "mL/min (simulation 1) and 0-170 mL/min (simulations 2-3)."
      ),
      source_name        = "CG"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 21L,
    n_studies      = 1L,
    n_observations = 144L,
    age_range      = "IQR 62-76 years (mean 67.10, SD 10.47; median 67)",
    age_median     = "67 years",
    weight_range   = "IQR 54.4-65 kg (mean 59.68, SD 10.53; median 60)",
    weight_median  = "60 kg",
    height_range   = "IQR 155-172 cm (mean 164.05, SD 8.64; median 165)",
    bsa_range      = "IQR 1.56-1.75 m^2 (mean 1.64, SD 0.17; median 1.60)",
    sex_female_pct = 42.9,
    disease_state  = "Critically ill adults with hospital-acquired pneumonia (HAP) or ventilator-associated pneumonia (VAP); cefepime also indicated for empirical management of septic shock of unknown source. Severity at enrollment: SOFA mean 5.05 (SD 3.04, median 4, IQR 4-5); APACHE II mean 17.33 (SD 4.63, median 17, IQR 14-22). Patients with a history of beta-lactam allergy were excluded.",
    renal_function = "Cockcroft-Gault CLCR mean 79.02 mL/min (SD 47.67), median 77.2 (IQR 45.3-125). Serum creatinine mean 1.07 mg/dL (SD 0.80), median 0.80 (IQR 0.53-1.12); only two patients exceeded 2 mg/dL (3.44 and 2.71 mg/dL). BUN mean 22.48 mg/dL (SD 9.60).",
    dose_range     = "2 g cefepime as a 30-minute IV infusion every 8, 12 or 24 h; seven blood samples per patient over the first dose interval",
    regions        = "Republic of Korea (single 840-bed university-affiliated tertiary referral centre: Hallym University Sacred Heart Hospital, Anyang). Prospective enrollment September-November 2019; IRB 2019-05-033.",
    notes          = paste(
      "Baseline demographics per Seo 2023 Table 1 (n = 21; 12 men, 9 women).",
      "Serum protein mean 5.42 g/dL (SD 0.67), albumin mean 2.87 g/dL (SD 0.43).",
      "144 plasma concentrations entered the population PK analysis. Fit in",
      "NONMEM 7.5 with FOCE-I; one-, two- and three-compartment structures were",
      "compared (OFV 307.124 / 162.452 / 147.837 without covariates) and the",
      "three-compartment model failed to converge, so the two-compartment model",
      "was selected. Adding Cockcroft-Gault CLCR on CL reduced the OFV to",
      "140.143 and the IIV on CL from 59.0% to 33.7%. Evaluated by",
      "prediction- and variability-corrected VPC (1,000 datasets,",
      "Supplementary Figure 1) and a 2,000-sample nonparametric bootstrap",
      "(Table 2 bootstrap column). Sex, age, weight, height, BSA, serum",
      "protein, serum albumin, serum creatinine and five renal-function",
      "estimates were screened; only Cockcroft-Gault CLCR on CL was retained.",
      "Cefepime assayed by HPLC-MS/MS with a 0.5 mg/L lower limit of",
      "quantitation."
    )
  )

  ini({
    # ----- Structural PK: Seo 2023 Table 2, 'Structural model'. -----
    # The table prints total clearance as CL = theta_1 * (CG / 77.21)^theta_2,
    # where CG is the Cockcroft-Gault creatinine clearance in mL/min.
    lcl <- log(6.60);  label("Typical CL at CRCL = 77.21 mL/min (L/h)")            # Seo 2023 Table 2: theta_1 = 6.60 L/h, RSE 7.91% (bootstrap 6.63 [5.55-7.64])
    e_crcl_cl <- 0.656; label("Power exponent on (CRCL / 77.21 mL/min) for CL (unitless)") # Seo 2023 Table 2: theta_2 = 0.656, RSE 10.7% (bootstrap 0.650 [0.438-0.796])
    lvc <- log(13.3);  label("Central volume of distribution, V1 (L)")             # Seo 2023 Table 2: V1 = 13.3 L, RSE 9.79% (bootstrap 13.4 [10.8-16.5])
    lq  <- log(16.5);  label("Intercompartmental clearance between V1 and V2, Q (L/h)") # Seo 2023 Table 2: Q = 16.5 L/h, RSE 18.7% (bootstrap 16.2 [9.59-23.7])
    lvp <- log(13.0);  label("Peripheral volume of distribution, V2 (L)")          # Seo 2023 Table 2: V2 = 13.0 L, RSE 10.6% (bootstrap 12.8 [9.87-16.3])

    # Reference (centring) covariate value -- a structural constant of the
    # published clearance equation, not an estimated quantity. 77.21 mL/min is
    # the cohort median Cockcroft-Gault CLCR (Table 1 reports it rounded to
    # 77.2 mL/min).
    crcl_ref <- fixed(77.21); label("Reference Cockcroft-Gault creatinine clearance (mL/min)") # Seo 2023 Table 2 denominator (CG/77.21); Table 1 cohort median 77.2

    # Plasma unbound fraction, fixed a priori (not estimated) for the Monte
    # Carlo PK/PD target-attainment analysis rather than for the PK fit itself.
    # Exposed here so a simulation can compute the paper's fT>MIC index
    # directly from the model output.
    fu <- fixed(0.81); label("Fraction of cefepime unbound in plasma (unitless)")  # Seo 2023 Methods 'PD target attainment': "The f was fixed at 81%", citing Okamoto 1993 (reference 14)

    # ----- Inter-individual variability -- Seo 2023 Table 2, 'Interindividual
    # variability', reported as percentages under the exponential IIV model
    # declared in Methods ("theta_i = theta x exp(eta_i)", eta ~ N(0, omega^2)).
    # Encoded with the NONMEM convention that the reported percentage is
    # omega x 100 (the approximate CV), so omega^2 = (pct / 100)^2:
    #   CL: 0.337^2 = 0.113569
    #   V1: 0.341^2 = 0.116281
    #   Q : 0.508^2 = 0.258064
    #   V2: 0.401^2 = 0.160801
    # The paper gives no footnote defining its percentage, but the SAME "(%)"
    # formatting is used one block lower for the proportional residual error
    # (7.62%), where sqrt(sigma^2) x 100 is the only reading that makes sense --
    # a proportional error's CV IS its standard deviation. Reading the IIV
    # column the same way is therefore the internally consistent choice.
    # The alternative log-normal reading omega^2 = log(1 + CV^2) would give
    # 0.107588 / 0.110037 / 0.230287 / 0.148840 -- 5-11% smaller in variance,
    # immaterial against these parameters' own uncertainty (RSE 18-25%, and
    # bootstrap 95% CIs spanning e.g. 17.0-46.8% for CL). No correlations
    # between etas are reported, so each is independent.
    etalcl ~ 0.113569 # Seo 2023 Table 2: IIV CL 33.7%, RSE 25.2%, shrinkage 0.000% (bootstrap 32.0 [17.0-46.8])
    etalvc ~ 0.116281 # Seo 2023 Table 2: IIV V1 34.1%, RSE 24.9%, shrinkage 11.2% (bootstrap 32.9 [7.52-49.5])
    etalq  ~ 0.258064 # Seo 2023 Table 2: IIV Q 50.8%, RSE 18.4%, shrinkage 19.6% (bootstrap 46.5 [0.000-65.8])
    etalvp ~ 0.160801 # Seo 2023 Table 2: IIV V2 40.1%, RSE 21.2%, shrinkage 8.84% (bootstrap 38.6 [16.7-54.3])

    # ----- Residual error -- proportional only. Seo 2023 Table 2 'Residual
    # variability / Proportional error (%) 7.62'; Methods list additive,
    # proportional and combined models as candidates.
    propSd <- 0.0762; label("Proportional residual error (fraction)") # Seo 2023 Table 2: 7.62%, RSE 9.43% (bootstrap 7.61 [6.10-9.23])
  })

  model({
    # ----- Individual PK parameters -----
    # CL = theta_1 * (CG / 77.21)^theta_2 (Seo 2023 Table 2). V1, Q and V2
    # carry no covariate; only CL does.
    cl <- exp(lcl + etalcl) * (CRCL / crcl_ref)^e_crcl_cl
    vc <- exp(lvc + etalvc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    # ----- Micro-constants -----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ----- ODE system (two-compartment, IV infusion into central) -----
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ----- Observation -----
    # Dose in mg, vc in L -> central / vc has units mg/L (= ug/mL).
    Cc <- central / vc
    # Unbound plasma concentration driving the paper's fT>MIC PK/PD index
    # (Methods 'PD target attainment'). Returned alongside Cc; it carries no
    # residual error of its own.
    Cu <- fu * Cc
    Cc ~ prop(propSd)
  })
}
