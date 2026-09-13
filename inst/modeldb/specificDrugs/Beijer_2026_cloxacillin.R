Beijer_2026_cloxacillin <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous cloxacillin surgical prophylaxis in 200",
    "adults undergoing primary elective total hip or total knee arthroplasty (Beijer 2026). The",
    "disposition is parameterised entirely on UNBOUND cloxacillin -- Figure S1 shows clearance and",
    "intercompartmental clearance both acting on the unbound pool, with V1 and V2 the unbound",
    "distribution volumes -- so central / vc is the unbound plasma concentration Cc that the study",
    "measured directly by ultrafiltration. Unbound clearance carries body weight and relative",
    "eGFR (Lund-Malmo Revised 2018) as power covariates. The observed TOTAL plasma concentration",
    "Ctot is then recovered algebraically from Cc through a one-site saturable plasma-protein",
    "binding model, which reproduces the concentration-dependent rise in unbound fraction the",
    "paper reports (median plasma protein binding 91%, range 69-98%). Both outputs carry their own",
    "log-additive residual error. The paper uses the model to show that 18-22% of patients fall",
    "below an unbound 2 mg/L target within the recommended 2 h interval between the first two 2 g",
    "doses, and that a 1 g/h continuous infusion after a 1 g loading dose holds >99% of patients",
    "above target."
  )
  reference <- paste(
    "Beijer G, Wallander K, Soderquist B, Giske CG, Breuer O, Eriksen J, Eliasson E.",
    "Optimizing cloxacillin prophylaxis in hip and knee arthroplasty based on population",
    "pharmacokinetics of unbound plasma concentrations.",
    "J Antimicrob Chemother. 2026. doi:10.1093/jac/dkag116.",
    "All parameter values are from Table S2 of the Supplementary Material; the protein-binding",
    "equation is Supplementary Material Eq. 1 and the model schematic is Figure S1."
  )
  vignette <- "Beijer_2026_cloxacillin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power covariate on unbound clearance, exponent 0.67, centred on 84 kg (Table S2 CL row:",
        "'CL (L/h) = theta_CL x (BW/84)^beta1 x (eGFR/67)^beta2 x exp(eta_CL)'). The 84 kg",
        "reference is the typical patient of the Figure 1 caption ('a body weight of 84 kg and",
        "eGFR 67 mL/min/1.73 m^2'), not the Table 1 cohort median of 83 kg. Cohort median 83 kg,",
        "IQR 73-95, range 53-185 (Table 1). Total body weight was chosen over ideal and adjusted",
        "body weight, BMI and BSA, all of which were screened (Methods, Pharmacokinetic",
        "modelling). It was retained alongside relative eGFR specifically to break the",
        "collinearity that absolute eGFR had with every anthropometric covariate (Supplementary",
        "Material, Renal function covariate selection); adding it cut unexplained CL variability",
        "from 50% to 47%."
      ),
      source_name        = "BW"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate, Lund-Malmo Revised 2018 equation, BSA-normalized",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power covariate on unbound clearance, exponent 0.51, centred on 67 mL/min/1.73 m^2",
        "(Table S2 CL row). The estimating equation is the Lund-Malmo Revised 2018 formula",
        "(LMR18), which Table S1 shows correlated more strongly with eta_CL than CKD-EPI 2021 or",
        "Cockcroft-Gault in either relative or absolute form (r = 0.32 for LMR18rel against 0.15",
        "for CKDEPI-2021rel). This column is the RELATIVE (BSA-normalized) LMR18 estimate; the",
        "paper's absolute LMR18 estimate correlated even better (r = 0.43) but was rejected",
        "because it left body size independently correlated with clearance. The 67",
        "mL/min/1.73 m^2 reference is the Figure 1 typical patient, not the Table 1 cohort median",
        "of 72. Cohort median relative eGFR 72 mL/min/1.73 m^2, IQR 61-84, range 18-142 (Table 1).",
        "Renal function was generally well preserved with few patients at very low levels, so the",
        "paper states the model carries little information about severe renal impairment",
        "(Strengths and limitations)."
      ),
      source_name        = "eGFR-LMR18rel"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "year",
      type        = "continuous",
      notes       = "Screened but not retained (Methods, Pharmacokinetic modelling; Results: 'No other covariates were found to significantly improve the model fit')."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained. Cohort 105/200 (53%) female (Table 1)."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened but not retained; used only to compute BMI and Du Bois BSA."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened but not retained; total body weight was the retained body-size covariate. Cohort median 28 kg/m^2, IQR 25-32, range 20-61 (Table 1)."
    ),
    BSA = list(
      description = "Body surface area (Du Bois formula)",
      units       = "m^2",
      type        = "continuous",
      notes       = "Screened but not retained as a covariate; used to convert between relative and absolute eGFR estimates (Methods, Pharmacokinetic modelling)."
    ),
    CREAT = list(
      description = "Plasma creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened but not retained; entered the model only through the eGFR equations. Cohort median 70 umol/L, IQR 58-89, range 11-236 (Table 1)."
    ),
    ALB = list(
      description = "Plasma albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = paste(
        "Screened and significantly correlated with protein binding, but NOT retained: the",
        "Supplementary Material states that adding albumin as a covariate on Bmax (and Kd) 'led to",
        "over-fitting and did not improve overall model performance'. Cohort median 34 g/L,",
        "IQR 32-36, range 15-42 (Table 1)."
      )
    ),
    ASA_CLASS = list(
      description = "American Society of Anesthesiologists preoperative physical-status class",
      units       = "(ordinal, I-V)",
      type        = "categorical",
      notes       = "Screened but not retained. Cohort I 29 (15%), II 84 (42%), III 87 (43%) (Table 1)."
    ),
    SURGSITE_KNEE = list(
      description = "Knee (rather than hip) arthroplasty indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained; the type of arthroplasty and the type of prosthesis were both evaluated as categorical covariates. Cohort 105 knee (53%), 95 hip (47%) (Table 1)."
    )
  )

  compartmentData <- list(
    central = list(
      analyte  = "cloxacillin, unbound",
      units    = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte  = "cloxacillin, unbound",
      units    = "mg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 200L,
    n_studies      = 1L,
    n_observations = "496 paired total and unbound plasma cloxacillin samples; median 3 per patient (range 1-3); 193 (97%) evaluable at start of surgery, 178 (89%) at end of surgery, 125 (63%) at 90 min postoperatively",
    age_range      = "36-90 years",
    age_median     = "73 years (IQR 65-78)",
    weight_range   = "53-185 kg",
    weight_median  = "83 kg (IQR 73-95)",
    sex_female_pct = 52.5,
    disease_state  = "Adults (>18 years) scheduled for primary elective total hip arthroplasty (95, 47%) or total knee arthroplasty (105, 53%), non-allergic to penicillin; ASA class I 15%, II 42%, III 43%",
    renal_function = "Relative eGFR (LMR18) median 72 mL/min/1.73 m^2 (IQR 61-84, range 18-142); absolute eGFR median 83 mL/min (IQR 68-93, range 21-158); plasma creatinine median 70 umol/L (range 11-236). Generally well preserved, with few patients at very low levels",
    dose_range     = "Intravenous cloxacillin 2 g at three time points: 30-45 min before surgical incision, then 2 h and 6 h after the start of the first dose. Recommended infusion duration 20-30 min; actual durations were 20-30 min in 100 (50%), <20 min in 90 (45%) and >30 min in 10 (5%) of patients, and only 42 (21%) of preoperative doses complied with the guideline in full",
    regions        = "Sweden (two centres)",
    notes          = paste(
      "Prospective two-centre study, enrolment 2022-2024 (Swedish Ethical Review Authority",
      "2021-02358 and 2024-06520-02). Baseline characteristics from Table 1. Total and unbound",
      "plasma cloxacillin were both measured by reversed-phase HPLC-MS/MS, the unbound fraction",
      "separated by ultrafiltration through a 10 kDa membrane; analytical ranges 0.1-100 mg/L",
      "total and 0.01-15 mg/L unbound. Sampling times were dictated by surgical events (start of",
      "surgery, end of surgery, 90 min postoperatively) rather than by dosing, which the paper",
      "notes gave good coverage of the whole post-dose time course. Observed median plasma protein",
      "binding 91% (range 69-98%), with 98/200 patients (49%) below 90% in at least one sample.",
      "Estimation used Monolix 2024R1; a 1000-sample non-parametric bootstrap supports every",
      "estimate in Table S2."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- all from Table S2 of the Supplementary
    # Material, which also reports the RSE and the 1000-sample bootstrap
    # median with 95% CI for each estimate.
    #
    # The disposition is on UNBOUND cloxacillin: Figure S1 draws CL and Q
    # leaving and connecting the unbound pools, with V1 and V2 the
    # unbound distribution volumes and plasma protein binding hanging off
    # V1 as an equilibrium (Kd, Bmax) rather than as a mass-transfer arm.
    #
    # Monolix reports omega as the STANDARD DEVIATION of the (log-normal)
    # random effect, whereas nlmixr2's ini() takes the VARIANCE, so every
    # omega below is squared. That the Table S2 omegas are log-scale SDs
    # is confirmed by its own CV column: 0.45 -> sqrt(exp(0.45^2) - 1) =
    # 47%, 0.64 -> 71%, 0.79 -> 93%, matching the printed 47 / 72 / 94%.
    # ------------------------------------------------------------------
    lcl <- log(76.4); label("Unbound cloxacillin clearance at 84 kg and eGFR 67 mL/min/1.73 m^2 (L/h)")  # Table S2: theta_CL = 76.4 (RSE 4.4%), bootstrap median 74.3 (61.4-83.3)
    lvc <- log(57.1); label("Central volume of distribution V1 for unbound cloxacillin (L)")             # Table S2: theta_V1 = 57.1 (RSE 7.7%), bootstrap median 56.5 (43.8-70.8)
    lq  <- log(57.2); label("Intercompartmental clearance Q for unbound cloxacillin (L/h)")              # Table S2: theta_Q = 57.2 (RSE 12%), bootstrap median 56.1 (42.6-96.3)
    lvp <- log(69.0); label("Peripheral volume of distribution V2 for unbound cloxacillin (L)")          # Table S2: theta_V2 = 69.0 (RSE 10%), bootstrap median 72.3 (53.4-153)

    # Covariate effects on unbound clearance. Table S2 prints the CL row
    # as an equation, 'CL (L/h) = theta_CL x (BW/84)^beta1 x
    # (eGFR/67)^beta2 x exp(eta_CL)'; the two centring constants are the
    # typical patient of the Figure 1 caption (84 kg, eGFR 67
    # mL/min/1.73 m^2) and NOT the Table 1 cohort medians of 83 kg and
    # 72 mL/min/1.73 m^2. Neither exponent is fixed -- both carry an RSE
    # and a bootstrap CI that excludes zero.
    e_wt_cl   <- 0.67; label("Power exponent on (WT/84) for unbound clearance (unitless)")    # Table S2: beta1 = 0.67 (RSE 22%), bootstrap median 0.71 (0.39-1.07)
    e_crcl_cl <- 0.51; label("Power exponent on (CRCL/67) for unbound clearance (unitless)")  # Table S2: beta2 = 0.51 (RSE 18%), bootstrap median 0.53 (0.33-0.75)

    # ------------------------------------------------------------------
    # Plasma protein binding. Supplementary Material Eq. 1 gives the
    # one-site specific binding model
    #
    #   Ctot = Cu + Bmax * Cu / (Kd + Cu)
    #
    # with Ctot and Cu in mg/L, Bmax the maximum binding capacity in
    # mg/L and Kd the dissociation constant in mg/L. Both were estimated
    # inside the population model (Table S2), not carried in from an
    # external binding study, and Bmax carries IIV. There is no linear
    # non-saturable arm in this paper's form, i.e. the register's `kns`
    # term is absent rather than zero.
    # ------------------------------------------------------------------
    lbmax <- log(559);  label("Maximum plasma protein binding capacity Bmax (mg/L)")           # Table S2: theta_Bmax = 559 (RSE 1.6%), bootstrap median 564 (511-625)
    lkd   <- log(47.8); label("Plasma protein binding dissociation constant Kd (mg/L)")        # Table S2: theta_Kd = 47.8 (RSE 2.6%), bootstrap median 48.7 (43.0-55.0)

    # ------------------------------------------------------------------
    # Between-subject variability. Table S2 reports omega on the SD scale
    # (see the header comment); squared here. IIV was estimable on CL,
    # V1, Q and Bmax only -- V2 and Kd carry none.
    #
    # CL and V1 are positively correlated. Monolix estimated the
    # correlation itself as a parameter rather than the covariance
    # (Supplementary Material, Correlation between parameters:
    # rho = omega_CL,V1 / (omega_CL x omega_V1)), so the off-diagonal is
    # reconstructed as rho x omega_CL x omega_V1 =
    # 0.68 x 0.45 x 0.64 = 0.19584.
    # ------------------------------------------------------------------
    etalcl + etalvc ~ c(0.2025,
                        0.19584, 0.4096)  # Table S2: omega_CL 0.45 (CV 47%, RSE 7.1%), omega_V1 0.64 (CV 72%, RSE 9.2%), rho_CL~V1 0.68 (RSE 9.8%); variances 0.45^2 and 0.64^2, covariance 0.68 x 0.45 x 0.64
    etalq   ~ 0.6241                      # Table S2: omega_Q = 0.79 (CV 94%, RSE 13%), bootstrap median 0.77 (0.49-1.19); variance = 0.79^2
    etalbmax ~ 0.0025                     # Table S2: omega_Bmax = 0.05 (CV 4.7%, RSE 25%), bootstrap median 0.05 (0.03-0.10); variance = 0.05^2. EBE eta-shrinkage was 77% for this parameter

    # ------------------------------------------------------------------
    # Residual error. Table S2 states the error model as
    # 'log(Y) = log(f) + a x eps', i.e. additive on the LOG scale, which
    # is nlmixr2's lnorm() exponential error. One term per output; the
    # study measured total and unbound concentrations with two separate
    # HPLC-MS/MS assays over different analytical ranges, so the two
    # magnitudes differ.
    # ------------------------------------------------------------------
    expSd      <- 0.24; label("Log-scale additive residual SD for unbound cloxacillin, the Cc output (log units)")  # Table S2: a(unbound) = 0.24 (RSE 4.4%), bootstrap median 0.24 (0.20-0.27); epsilon-shrinkage 19.5%
    expSd_Ctot <- 0.19; label("Log-scale additive residual SD for total cloxacillin (log units)")                    # Table S2: a(total) = 0.19 (RSE 4.6%), bootstrap median 0.19 (0.16-0.22)
  })

  model({
    # 1. Individual parameters. Unbound clearance carries body weight and
    #    relative eGFR as power terms centred on the Figure 1 typical
    #    patient (84 kg, eGFR 67 mL/min/1.73 m^2).
    cl   <- exp(lcl + etalcl) * (WT / 84)^e_wt_cl * (CRCL / 67)^e_crcl_cl
    vc   <- exp(lvc + etalvc)
    q    <- exp(lq + etalq)
    vp   <- exp(lvp)
    bmax <- exp(lbmax + etalbmax)
    kd   <- exp(lkd)

    # 2. Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 3. ODE system. Both states hold UNBOUND cloxacillin: per Figure S1
    #    the elimination arm CL and the distribution arm Q both act on the
    #    unbound pools, and V1 / V2 are unbound distribution volumes. The
    #    administered amount enters `central` unchanged, with no
    #    bioavailability term -- the paper doses 2 g intravenously and its
    #    typical-patient profile is reproduced by putting the whole 2000 mg
    #    into this compartment (see the vignette's mass-balance and
    #    typical-patient gates).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 4. Observations. `central` holds the unbound pool, so central / vc is
    #    the UNBOUND plasma concentration; by the library convention that
    #    quantity is `Cc` whatever the pool contains, and the label records
    #    that it is unbound (same treatment as Hennig_2015_phenytoin). The
    #    total plasma concentration follows algebraically from
    #    Supplementary Material Eq. 1. The binding term is an output
    #    transformation, not a reservoir: no mass moves into it, so `Ctot`
    #    does not feed back on the disposition.
    Cc   <- central / vc
    Ctot <- Cc + bmax * Cc / (kd + Cc)

    Cc   ~ lnorm(expSd)
    Ctot ~ lnorm(expSd_Ctot)
  })
}
