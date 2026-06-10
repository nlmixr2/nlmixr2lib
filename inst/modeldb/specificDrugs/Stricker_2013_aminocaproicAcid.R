Stricker_2013_aminocaproicAcid <- function() {
  description <- paste(
    "Two-compartment IV population PK model for epsilon-aminocaproic acid",
    "(EACA) in infants (2-24 months) undergoing craniofacial reconstruction",
    "surgery. Allometric scaling on body weight (reference 8.82 kg; fixed",
    "exponents 0.75 on CL and Q, 1.0 on V1 and V2), an asymptotically",
    "increasing post-natal age maturation effect on clearance (age50 = 7.36",
    "weeks), and binary intra-operative-period multipliers on CL (0.89) and",
    "V1 (0.80) capturing the composite effect of anaesthesia, blood loss,",
    "and surgical fluid management. Parameter values from Stricker 2013",
    "Table 4."
  )
  reference <- paste(
    "Stricker PA, Zuppa AF, Fiadjoe JE, Maxwell LG, Sussman EM, Pruitt EY,",
    "Goebel TK, Gastonguay MR, Taylor JA, Bartlett SP, Schreiner MS.",
    "Population pharmacokinetics of epsilon-aminocaproic acid in infants",
    "undergoing craniofacial reconstruction surgery.",
    "Br J Anaesth 2013;112(5):788-795. doi:10.1093/bja/aes507.",
    sep = " "
  )
  vignette <- "Stricker_2013_aminocaproicAcid"
  units <- list(
    time          = "min",
    dosing        = "mg",
    concentration = "mg/L"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight at study inclusion",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed within the surgical window. Allometric scaling per",
        "Stricker 2013 Table 4 footnote: CL and Q are scaled by",
        "(WT/8.82)^0.75 and V1 and V2 by (WT/8.82)^1.0. Reference weight",
        "8.82 kg is the median weight of the 18 infants in the cohort",
        "(Methods 'Base model'). The allometric exponents 0.75 and 1 are",
        "fixed per physiologic-size convention (Anderson and Holford 2008)."
      ),
      source_name        = "Weight"
    ),
    PNA = list(
      description        = "Postnatal age",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed within the surgical window. The paper reports postnatal",
        "age in weeks and models an asymptotically increasing maturation",
        "effect on CL: CL_age = AGE / (7.36 + AGE) with AGE in weeks (Table",
        "4 footnote). The canonical PNA carries months, so inside model()",
        "the maturation formula is reparameterised in months using the",
        "conversion 1 month = 4.34524 weeks: the published age50 of",
        "7.36 weeks becomes 7.36 / 4.34524 = 1.6938 months. The Zhao 2018",
        "precedent established this PNA-unit-reparameterisation pattern.",
        "Source column 'age (weeks)' -> canonical PNA in months."
      ),
      source_name        = "age (weeks)"
    ),
    INTRAOP = list(
      description        = "Intra-operative period indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (pre- or post-operative)",
      notes              = paste(
        "Time-varying within subject. 1 = the observation falls within the",
        "intra-operative period (between the post-loading-dose PK sample",
        "and skin closure / end of surgery, per Stricker 2013 Methods 'PK",
        "modelling and simulation' and Discussion); 0 = pre- or",
        "post-operative period. CL and V1 carry a multiplicative shift",
        "during the intra-operative window: CL_intraop = 0.89 * CL_postop,",
        "V1_intraop = 0.80 * V1_postop (Table 4)."
      ),
      source_name        = "intra-operative period"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 18L,
    n_studies      = 1L,
    age_range      = "27-107 weeks postnatal (2-24 months)",
    age_median     = "39 weeks",
    weight_range   = "6.7-11.8 kg",
    weight_median  = "8.8 kg",
    sex_female_pct = 55.6,
    disease_state  = paste(
      "Healthy infants 2-24 months of age undergoing craniofacial",
      "reconstruction surgery (most commonly fronto-orbital advancement",
      "for unicoronal, metopic, lambdoid or sagittal synostosis, and a",
      "smaller number with Pfeiffer or Saethre-Chotzen syndromes).",
      "Patients with abnormal renal function or known coagulation",
      "disorders were excluded."
    ),
    dose_range     = paste(
      "Three dose-escalation cohorts of 6 infants each. Each subject",
      "received an IV loading dose over 10 min followed by a continuous IV",
      "infusion that ran until skin closure. Cohort 1: 25 mg/kg loading +",
      "10 mg/kg/h CIVI. Cohort 2: 50 mg/kg loading + 20 mg/kg/h CIVI.",
      "Cohort 3: 100 mg/kg loading + 40 mg/kg/h CIVI. Median (range)",
      "duration of CIVI was 230 (111-342) min in Cohort 1, 227 (169-365)",
      "in Cohort 2, and 254 (219-280) in Cohort 3."
    ),
    regions        = "USA (single-centre, Children's Hospital of Philadelphia)",
    notes          = paste(
      "Open-label, non-randomised, dose-escalation PK trial conducted",
      "under FDA IND 105,301. 18 infants enrolled in three sequential",
      "cohorts; up to 12 PK samples per subject collected via arterial or",
      "central venous catheter at pre-dose, post-loading-dose, during the",
      "continuous infusion (0.5, 2, 4-6 h after CIVI start, plus end of",
      "infusion), and post-infusion (0.5, 3, 6, 9, 12, 15 h). Plasma EACA",
      "measured by validated HPLC-MS/MS (LLOQ 1 mg/L; intraday precision",
      "0.3-2%, accuracy 89-102%; validated range 1-250 mg/L). Baseline",
      "demographics in Table 1; intra-operative blood loss and fluid",
      "balance in Table 2. NONMEM VI Level 2.0, ADVAN3 TRANS4, FOCE-I.",
      "Sex distribution: 10 female / 8 male per Table 1."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural PK parameters -- Stricker 2013 Table 4.
    # Reference subject: WT = 8.82 kg and PNA -> infinity (typical-value
    # postoperative clearance for a fully-mature subject -- the intra-op
    # multiplier and age maturation are applied inside model()).
    # NONMEM parameterisation: ADVAN3 TRANS4 (two-compartment, CL/V1/Q/V2).
    # The paper reports clearances in mL/min; here they are stored on a
    # log(L/min) scale so kel = cl/vc is in 1/min when vc is in L. The
    # source-trace comments preserve the paper's mL/min values.
    # ------------------------------------------------------------------
    lcl   <- log(37.6 / 1000); label("Postoperative clearance (L/min)")               # Table 4: 37.6 mL/min = 0.0376 L/min
    lvc   <- log(1.27);        label("Postoperative central volume of distribution (L)")  # Table 4: V1 = 1.27 L
    lq    <- log(42.3 / 1000); label("Inter-compartmental clearance (L/min)")          # Table 4: 42.3 mL/min = 0.0423 L/min
    lvp   <- log(2.53);        label("Peripheral volume of distribution (L)")           # Table 4: V2 = 2.53 L

    # ------------------------------------------------------------------
    # Allometric exponents (fixed per Table 4 footnote: 0.75 on
    # clearances, 1 on volumes; physiologic-size convention).
    # ------------------------------------------------------------------
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent of (WT/8.82) on CL and Q (unitless)")   # Table 4 footnote (fixed)
    e_wt_vc_vp <- fixed(1.0);  label("Allometric exponent of (WT/8.82) on V1 and V2 (unitless)")  # Table 4 footnote (fixed)

    # ------------------------------------------------------------------
    # Age maturation on CL: asymptotic Emax form
    #     Mat_CL = AGE / (age50 + AGE)
    # where AGE is postnatal age in weeks and age50 = 7.36 weeks per
    # Table 4 footnote. Re-expressed in canonical PNA (months) inside
    # model() using the conversion 1 month = 4.34524 weeks:
    #     age50_months = 7.36 / 4.34524 = 1.6938 months.
    # Stored on log scale (positive-constrained parameter).
    # ------------------------------------------------------------------
    lage50_pna <- log(7.36 / 4.34524); label("Postnatal age at 50% mature CL (log months)")  # Table 4 footnote (7.36 weeks)

    # ------------------------------------------------------------------
    # Intra-operative multiplicative shifts on CL and V1 (Table 4).
    # Stored on the log scale so they enter the model as
    #     cl = exp(lcl + lr_intraop_cl * INTRAOP)
    # giving cl_typical = 37.6 * 0.89 = 33.5 mL/min during the
    # intra-operative window (INTRAOP = 1) and cl_typical = 37.6
    # mL/min pre/postoperatively (INTRAOP = 0). Table 4 reports the
    # ratios CL_intra/CL_post = 0.89 and V1_intra/V1_post = 0.80.
    # ------------------------------------------------------------------
    lr_intraop_cl <- log(0.89); label("Intra-op CL ratio (log of CL_intra/CL_post)")   # Table 4
    lr_intraop_vc <- log(0.80); label("Intra-op V1 ratio (log of V1_intra/V1_post)")   # Table 4

    # ------------------------------------------------------------------
    # IIV -- correlated block on CL and V1 (Table 4 + caption).
    # The paper reports "between-subject variability = sqrt(variance) *
    # 100", so the table values are %CV-on-the-eta scale:
    #     omega^2_CL = (16.79/100)^2 = 0.0281904
    #     omega^2_V1 = (47.01/100)^2 = 0.2209940
    # Covariance between eta_CL and eta_V1 is 0.03 (Table 4 caption).
    # That implies a correlation of 0.03 / sqrt(0.02819 * 0.22099) = 0.382.
    # Q and V2 had no IIV estimated (over-parameterisation; the covariance
    # step failed when BSV was added to Q or V2 -- Table 3 runs 3 and 4).
    # ------------------------------------------------------------------
    etalcl + etalvc ~ c(0.02819,
                        0.03, 0.22099)  # Table 4 caption (sqrt(var)*100 convention) and cov(eta_CL,eta_V1) = 0.03

    # ------------------------------------------------------------------
    # Residual error -- combined additive + proportional (Table 4
    # "Residual variability", reported as variance):
    #     propSd = sqrt(0.03) = 0.1732 (fraction)
    #     addSd  = sqrt(0.6)  = 0.7746 (mg/L; same units as Cc)
    # Paper Methods: "C_obs = C_pred * (1 + eps_P) + eps_A"; matches
    # the add() + prop() error model in nlmixr2.
    # ------------------------------------------------------------------
    propSd <- sqrt(0.03); label("Proportional residual SD (fraction)")    # Table 4 (variance 0.03)
    addSd  <- sqrt(0.6);  label("Additive residual SD (mg/L)")            # Table 4 (variance 0.6)
  })

  model({
    # 1. Derived covariate terms.
    # Convert canonical PNA (months) into the maturation factor on CL.
    # The Emax-style ratio is dimensionless: with PNA and exp(lage50_pna)
    # both in months, mat_cl is identical to the paper's
    # AGE_weeks / (7.36 + AGE_weeks) reparameterised in months.
    mat_cl <- PNA / (exp(lage50_pna) + PNA)

    # 2. Individual PK parameters.
    # Postoperative typical-value PK parameters, scaled allometrically by
    # weight (reference 8.82 kg) with WT exponents 0.75 (CL/Q) and 1
    # (V1/V2). Age maturation is applied multiplicatively to CL only; V1
    # has no age term in Stricker 2013. The INTRAOP binary indicator
    # multiplies the postoperative typical values by the published
    # intra-op ratios.
    cl <- exp(lcl + etalcl + lr_intraop_cl * INTRAOP) *
          (WT / 8.82)^e_wt_cl_q * mat_cl
    vc <- exp(lvc + etalvc + lr_intraop_vc * INTRAOP) *
          (WT / 8.82)^e_wt_vc_vp
    q  <- exp(lq)  * (WT / 8.82)^e_wt_cl_q
    vp <- exp(lvp) * (WT / 8.82)^e_wt_vc_vp

    # 3. Micro-constants (1/min with cl, q in L/min and vc, vp in L).
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 4. ODE system -- two-compartment disposition with IV dosing
    # (loading-dose bolus over 10 min + CIVI directly to the central
    # compartment; no extravascular absorption).
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-           k12 * central - k21 * peripheral1

    # 5. Observation and error.
    # central amount in mg, vc in L -> Cc in mg/L (= ug/mL), matching
    # the paper's concentration units in Figure 1 and the 130 mg/L
    # therapeutic target.
    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
