CohenWolkowiez_2014_piperacillin <- function() {
  description <- "One-compartment population PK model for piperacillin in premature and term infants under 61 days postnatal age (Cohen-Wolkowiez 2014); linear body-weight scaling on CL and V (fixed exponent = 1) and a power effect of postmenstrual age on CL."
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
        "Time-varying. Linear allometric scaling on both CL and V with the",
        "WT exponent fixed at 1 (Cohen-Wolkowiez 2014 Table 2 and Population",
        "PK model building paragraph: estimation of a free body-size exponent",
        "and 3/4-power scaling were both rejected for lack of improvement in",
        "fit). The paper reports CL and V in L/kg/h and L/kg respectively, so",
        "the per-kg typical values multiply directly by WT with no centring",
        "(the cohort spanned 0.473-3.990 kg, median 1.439 kg, Table 1)."
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
        "exponent 1.76 (Cohen-Wolkowiez 2014 Table 3) and reference 33 weeks",
        "(close to the cohort median PMA of 32 weeks, Table 1). The source",
        "paper reports PMA in weeks; this model converts the canonical",
        "PAGE (months) back to weeks via pma_wk = PAGE * 4.35 so the",
        "published exponent and 33-week reference apply unchanged."
      ),
      source_name        = "PMA (weeks)"
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
      "Piperacillin 80 mg/kg every 8 h IV (30-min infusion) in cohorts 1-3",
      "(GA < 32 wk PNA < 14 d; GA < 32 wk PNA >= 14 d; GA >= 32 wk PNA < 14 d);",
      "100 mg/kg every 8 h in cohort 4 (GA >= 32 wk PNA >= 14 d). Piperacillin",
      "and tazobactam were coadministered at the fixed 8:1 piperacillin:tazobactam ratio."
    ),
    regions        = "United States (four neonatal intensive care units)",
    ga_range       = "23-40 weeks at birth (Table 1)",
    pma_range      = "25-48 weeks postmenstrual age (Table 1)",
    pna_range      = "1-60 days postnatal age (Table 1)",
    creat_range    = "0.3-2.0 mg/dL serum creatinine (Table 1; cohort median 0.8 mg/dL)",
    conmed_pct     = c(gentamicin = 63, dopamine = 13, epinephrine = 6),
    notes          = paste(
      "32 critically ill infants enrolled 2010-2011 at four US neonatal",
      "intensive care units. 128 timed plasma samples retained (median 4",
      "samples per infant); concurrent dried-blood-spot samples were also",
      "collected and analysed but the structural PK model on which this entry",
      "is based was fit to plasma alone. Cohen-Wolkowiez 2014 Table 1 and",
      "PK specimens paragraph."
    )
  )

  ini({
    # ===== Structural PK (Cohen-Wolkowiez 2014 Table 3 final piperacillin PMA-based model) =====
    # Typical CL and V are reported per kg of body weight (L/kg/h and L/kg);
    # the linear WT scaling (exponent fixed at 1) multiplies these directly.
    lcl <- log(0.080); label("Typical CL per kg body weight at PMA = 33 weeks (L/kg/h)")  # Cohen-Wolkowiez 2014 Table 3: typical CL value 0.080 L/kg/h (RSE 7.9%)
    lvc <- log(0.42);  label("Typical V per kg body weight (L/kg)")                       # Cohen-Wolkowiez 2014 Table 3: typical V  value 0.42  L/kg (RSE 9.6%)

    # Allometric exponents (Cohen-Wolkowiez 2014 Table 2 base model and Population PK model building paragraph)
    e_wt_cl <- fixed(1.0); label("Linear WT exponent on CL (unitless, fixed)")  # Population PK model building paragraph: WT exponent fixed at 1
    e_wt_vc <- fixed(1.0); label("Linear WT exponent on V (unitless, fixed)")   # Population PK model building paragraph: WT exponent fixed at 1

    # PMA maturation exponent on CL (Cohen-Wolkowiez 2014 Table 3 final piperacillin model)
    e_page_cl <- 1.76; label("Power exponent on (PMA / 33 weeks) for CL (unitless)")  # Cohen-Wolkowiez 2014 Table 3: PMA covariate 1.76 (RSE 33.6%)

    # ===== IIV (Cohen-Wolkowiez 2014 Table 3 final piperacillin model) =====
    # Exponential (log-normal) IIV; CV%-to-variance: omega^2 = log(1 + CV^2)
    # 37% CV -> log(1 + 0.37^2) = 0.12839
    etalcl ~ 0.12839  # Table 3: CL interindividual variability 37% CV (RSE 27.5%)
    # No etalvc: between-subject variability on V was fixed to 0 (boundary
    # value); see Population PK model building paragraph: "Inclusion of weight
    # as a covariate on V resulted in a between-subject variability of V
    # estimate close to 0 (boundary value). Consequently, between-subject
    # variability of V was fixed to 0 for subsequent model-building steps."

    # ===== Residual error (Cohen-Wolkowiez 2014 Table 3: combined proportional + additive) =====
    propSd <- 0.33; label("Proportional residual error (fraction)")     # Table 3: proportional residual error 33% CV (RSE 9.9%)
    addSd  <- 6.90; label("Additive residual error (mg/L)")             # Table 3: plasma additive residual error 6.90 mg/L (RSE 42.6%)
  })

  model({
    # ----- Derived covariate terms -----
    # Convert canonical PAGE (months) back to source-paper PMA (weeks) so the
    # published exponent and 33-week reference apply unchanged.
    pma_wk <- PAGE * 4.35

    # Cohen-Wolkowiez 2014 Table 2 / Table 3 final piperacillin CL covariate equation:
    #   CL = theta_CL * WT^1 * (PMA / 33)^COVPMA
    fpma <- (pma_wk / 33)^e_page_cl

    # ----- Individual PK parameters (Cohen-Wolkowiez 2014 Table 3 final piperacillin model) -----
    cl <- exp(lcl + etalcl) * WT^e_wt_cl * fpma
    vc <- exp(lvc)          * WT^e_wt_vc

    # ----- Micro-constants -----
    kel <- cl / vc

    # ----- ODE system -----
    # IV piperacillin-tazobactam (8:1 ratio) dosed into the central compartment.
    d/dt(central) <- -kel * central

    # ----- Output -----
    # Plasma piperacillin concentration: dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
