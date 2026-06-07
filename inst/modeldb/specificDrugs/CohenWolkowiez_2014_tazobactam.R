CohenWolkowiez_2014_tazobactam <- function() {
  description <- "One-compartment population PK model for tazobactam in premature and term infants under 61 days postnatal age (Cohen-Wolkowiez 2014); linear body-weight scaling on CL and V (fixed exponent = 1), and PMA, serum creatinine and concomitant gentamicin coadministration as covariates on CL."
  reference <- paste(
    "Cohen-Wolkowiez M, Watt KM, Zhou C, Bloom BT, Poindexter B,",
    "Castro L, Gao J, Capparelli EV, Benjamin DK Jr, Smith PB.",
    "Developmental pharmacokinetics of piperacillin and tazobactam",
    "using plasma and dried blood spots from infants.",
    "Antimicrob Agents Chemother. 2014;58(5):2856-2865.",
    "doi:10.1128/AAC.02139-13.",
    sep = " "
  )
  vignette <- "CohenWolkowiez_2014_piperacillin_tazobactam"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Linear allometric scaling on both CL and V with the WT",
        "exponent fixed at 1 (Cohen-Wolkowiez 2014 Table 2 and Population PK",
        "model building paragraph: estimation of a free body-size exponent and",
        "3/4-power scaling were both rejected for lack of improvement in fit).",
        "The paper reports CL and V in L/kg/h and L/kg respectively, so the",
        "per-kg typical values multiply directly by WT with no centring (the",
        "cohort spanned 0.473-3.990 kg, median 1.439 kg, Table 1)."
      ),
      source_name        = "WT"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age in weeks + postnatal age in weeks)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Drives the power-form maturation effect on CL with",
        "exponent 1.35 (Cohen-Wolkowiez 2014 Table 4) and reference 33 weeks",
        "(close to the cohort median PMA of 32 weeks, Table 1). The source",
        "paper reports PMA in weeks; this model converts the canonical PAGE",
        "(months) back to weeks via pma_wk = PAGE * 4.35 so the published",
        "exponent and 33-week reference apply unchanged."
      ),
      source_name        = "PMA (weeks)"
    ),
    CREAT = list(
      description        = "Serum creatinine concentration",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Power-form effect on CL: (CREAT / 0.5)^theta_SCR with",
        "theta_SCR = -0.35 (Cohen-Wolkowiez 2014 Table 4). Reference 0.5 mg/dL",
        "is below the cohort median (0.8 mg/dL, Table 1) but matches the value",
        "embedded in the Table 2 model-building equation (SCR/0.5)^COVSCR.",
        "Source paper column SCR; higher SCR (worse renal function) lowers CL.",
        "Missing baseline SCR values were imputed by the source authors using",
        "an exponential SCR-vs-PMA relationship derived from the study cohort",
        "(Population PK analysis paragraph)."
      ),
      source_name        = "SCR"
    ),
    CONMED_GEN = list(
      description        = "Concomitant gentamicin coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concurrent gentamicin)",
      notes              = paste(
        "Time-varying indicator: 1 if gentamicin was administered concurrently",
        "with piperacillin-tazobactam, 0 if not. Multiplicative power-form",
        "effect on CL: COVGENT^GENT with COVGENT = 1.52 (Cohen-Wolkowiez 2014",
        "Table 4). Multiplier is 1 when GENT = 0 and 1.52 when GENT = 1.",
        "63% of cohort received concurrent gentamicin (Table 1). Note from the",
        "Population PK model building paragraph: the gentamicin covariate was",
        "confounded by PMA / PNA differences between infants who received",
        "gentamicin and those who did not, so the effect should be interpreted",
        "as a marker of an older / sicker subset rather than a mechanistic",
        "drug-drug interaction."
      ),
      source_name        = "GENT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 32,
    n_studies      = 1,
    age_range      = "Postnatal age 1-60 days; postmenstrual age 25-48 weeks; gestational age 23-40 weeks",
    age_median     = "Postnatal age 8 days; postmenstrual age 32 weeks; gestational age 30 weeks",
    weight_range   = "0.473-3.990 kg",
    weight_median  = "1.439 kg",
    sex_female_pct = 37,
    race_ethnicity = "Not reported in the source paper",
    disease_state  = "Critically ill premature and term infants under 61 days postnatal age with suspected systemic infection (e.g. bacteraemia, necrotising enterocolitis, complicated intra-abdominal infection)",
    dose_range     = paste(
      "Tazobactam component of piperacillin-tazobactam dosed at the fixed 8:1",
      "piperacillin:tazobactam ratio: 10 mg/kg q8h IV (30-min infusion) for",
      "infants who received 80 mg/kg piperacillin (cohorts 1-3) and",
      "12.5 mg/kg q8h for infants who received 100 mg/kg piperacillin (cohort 4)."
    ),
    regions        = "United States (four neonatal intensive care units)",
    ga_range       = "23-40 weeks at birth (Table 1)",
    pma_range      = "25-48 weeks postmenstrual age (Table 1)",
    pna_range      = "1-60 days postnatal age (Table 1)",
    creat_range    = "0.3-2.0 mg/dL serum creatinine (Table 1; cohort median 0.8 mg/dL)",
    conmed_pct     = c(gentamicin = 63, dopamine = 13, epinephrine = 6),
    notes          = paste(
      "32 critically ill infants enrolled 2010-2011 at four US neonatal",
      "intensive care units. 128 timed plasma piperacillin samples and",
      "matched tazobactam samples were retained (median 4 samples per infant);",
      "concurrent dried-blood-spot samples were also collected and analysed",
      "but the structural PK model on which this entry is based was fit to",
      "plasma alone. The source authors selected the SCR + GENT + PMA",
      "tazobactam model as the final model on goodness-of-fit grounds but",
      "explicitly did not use it for dosing simulations or recommendations",
      "(Population PK model building paragraph and Results). Cohen-Wolkowiez",
      "2014 Table 1, Table 4, and PK specimens paragraph."
    )
  )

  ini({
    # ===== Structural PK (Cohen-Wolkowiez 2014 Table 4 final tazobactam model) =====
    # Typical CL and V are reported per kg of body weight (L/kg/h and L/kg);
    # the linear WT scaling (exponent fixed at 1) multiplies these directly.
    lcl <- log(0.088); label("Typical CL per kg body weight at PMA = 33 weeks, SCR = 0.5 mg/dL, no concurrent gentamicin (L/kg/h)")  # Cohen-Wolkowiez 2014 Table 4: typical CL value 0.088 L/kg/h (RSE 8.9%)
    lvc <- log(0.57);  label("Typical V per kg body weight (L/kg)")                                                                  # Cohen-Wolkowiez 2014 Table 4: typical V  value 0.57  L/kg (RSE 5.4%)

    # Allometric exponents (Cohen-Wolkowiez 2014 Table 2 base model and Population PK model building paragraph)
    e_wt_cl <- fixed(1.0); label("Linear WT exponent on CL (unitless, fixed)")  # Population PK model building paragraph: WT exponent fixed at 1
    e_wt_vc <- fixed(1.0); label("Linear WT exponent on V (unitless, fixed)")   # Population PK model building paragraph: WT exponent fixed at 1

    # Covariate effects on CL (Cohen-Wolkowiez 2014 Table 4 final tazobactam model)
    e_page_cl       <- 1.35;  label("Power exponent on (PMA / 33 weeks) for CL (unitless)")                                    # Cohen-Wolkowiez 2014 Table 4: PMA covariate 1.35  (RSE 27.3%)
    e_creat_cl      <- -0.35; label("Power exponent on (CREAT / 0.5 mg/dL) for CL (unitless)")                                 # Cohen-Wolkowiez 2014 Table 4: SCR covariate -0.35 (RSE 37.0%)
    e_conmed_gen_cl <- 1.52;  label("Multiplicative factor on CL for concurrent gentamicin coadministration (GENT = 1, unitless)")  # Cohen-Wolkowiez 2014 Table 4: GENT covariate 1.52  (RSE 10.1%)

    # ===== IIV (Cohen-Wolkowiez 2014 Table 4 final tazobactam model) =====
    # Exponential (log-normal) IIV; CV%-to-variance: omega^2 = log(1 + CV^2)
    # 23% CV -> log(1 + 0.23^2) = 0.05154
    etalcl ~ 0.05154  # Table 4: CL interindividual variability 23% CV (RSE 34.0%)
    # No etalvc: between-subject variability on V was fixed to 0 (boundary
    # value); see Population PK model building paragraph: "Inclusion of weight
    # as a covariate on V resulted in a between-subject variability of V
    # estimate close to 0 (boundary value). Consequently, between-subject
    # variability of V was fixed to 0 for subsequent model-building steps."

    # ===== Residual error (Cohen-Wolkowiez 2014 Table 4: combined proportional + additive) =====
    propSd <- 0.24; label("Proportional residual error (fraction)")     # Table 4: proportional residual error 24% CV (RSE 19.8%)
    addSd  <- 1.43; label("Additive residual error (mg/L)")             # Table 4: plasma additive residual error 1.43 mg/L (RSE 51.2%)
  })

  model({
    # ----- Derived covariate terms -----
    # Convert canonical PAGE (months) back to source-paper PMA (weeks) so the
    # published exponent and 33-week reference apply unchanged.
    pma_wk <- PAGE * 4.35

    # Cohen-Wolkowiez 2014 Table 2 / Table 4 final tazobactam CL covariate equation:
    #   CL = theta_CL * WT^1 * (PMA / 33)^COVPMA * (SCR / 0.5)^COVSCR * COVGENT^GENT
    fpma         <- (pma_wk / 33)^e_page_cl
    creat_factor <- (CREAT / 0.5)^e_creat_cl
    gent_factor  <- e_conmed_gen_cl^CONMED_GEN

    # ----- Individual PK parameters (Cohen-Wolkowiez 2014 Table 4 final tazobactam model) -----
    cl <- exp(lcl + etalcl) * WT^e_wt_cl * fpma * creat_factor * gent_factor
    vc <- exp(lvc)          * WT^e_wt_vc

    # ----- Micro-constants -----
    kel <- cl / vc

    # ----- ODE system -----
    # IV piperacillin-tazobactam (8:1 ratio) dosed into the central compartment;
    # this entry models the tazobactam component.
    d/dt(central) <- -kel * central

    # ----- Output -----
    # Plasma tazobactam concentration: dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
