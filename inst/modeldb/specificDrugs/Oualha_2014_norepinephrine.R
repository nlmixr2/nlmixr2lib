Oualha_2014_norepinephrine <- function() {
  description <- paste(
    "Population PK/PD model for continuous IV norepinephrine in",
    "hypotensive critically ill children (Oualha 2014). One-compartment",
    "open PK with first-order elimination plus an endogenous zero-order",
    "production rate q0 and circulating-volume-anchored Vc = 0.08 * WT;",
    "allometric scaling of CL and q0 on body weight (exponents fixed to",
    "3/4). Emax PD sub-model on mean arterial pressure (MAP) with a",
    "power-of-postmenstrual-age effect on basal MAP0 and a categorical",
    "organ-dysfunction-count effect on the maximal drug-induced MAP",
    "increase dMAP (32 mmHg for <=3 dysfunctions vs 12 mmHg for >=4).",
    sep = " "
  )
  reference <- paste(
    "Oualha M, Treluyer JM, Lesage F, de Saint Blanquat L, Dupic L,",
    "Hubert P, Spreux-Varoquaux O, Urien S (2014).",
    "Population pharmacokinetics and haemodynamic effects of",
    "norepinephrine in hypotensive critically ill children.",
    "British Journal of Clinical Pharmacology 78(4):886-897.",
    "doi:10.1111/bcp.12412.",
    sep = " "
  )
  vignette <- "Oualha_2014_norepinephrine"
  units <- list(time = "hour", dosing = "ug", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling on CL and q0 with paper-FIXED exponent 0.75,",
        "and linear scaling on Vc (Vc = 0.08 * WT L per the circulating-",
        "volume rule cited by Oualha 2014, Methods PK/PD modelling). Cohort",
        "WT range 2-85 kg (median 6.7 kg), Table 1.",
        sep = " "
      ),
      source_name        = "BW"
    ),
    PAGE = list(
      description        = "Postmenstrual age",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Oualha 2014 uses 'post-conceptional age (PCA)', which in",
        "paediatric clinical-pharmacology literature of that era is used",
        "synonymously with postmenstrual age (PNA + GA from LMP). The",
        "paper-reported PCA-in-months values are stored in the canonical",
        "PAGE column. Enters the MAP0 equation as a power-of-PCA effect:",
        "MAP0_i = MAP0 * (PCA / 9)^0.166. Reference value 9 months",
        "corresponds to a term newborn (GA 40 weeks ~= 9 months). Cohort",
        "range encompasses 7 premature neonates (GA 32-36 weeks) and older",
        "children up to chronological age 182 months (Results / Patient",
        "data).",
        sep = " "
      ),
      source_name        = "PCA"
    ),
    ORG_FAIL_COUNT = list(
      description        = "Number of failing organs in a critically ill patient.",
      units              = "(count)",
      type               = "categorical",
      reference_category = "<=3 organ dysfunctions (typical-value dMAP = 32 mmHg).",
      notes              = paste(
        "Oualha 2014 pools the integer count into two strata: <=3",
        "(reference) versus >=4. dMAP = 32 mmHg for the reference stratum",
        "and 12 mmHg for the >=4 stratum (Table 3). Decomposed inside",
        "model() as a binary indicator orgf_high = (ORG_FAIL_COUNT >= 4)",
        "that selects an additive log-shift on log(dMAP). Cohort: 14 of 38",
        "children with >3 dysfunctions (Table 1 / Patient data).",
        sep = " "
      ),
      source_name        = "number of organ dysfunctions"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 38L,
    n_studies      = 1L,
    age_range      = "0-182 months (paediatric; including 7 premature neonates at GA 32-36 weeks)",
    age_median     = "7.6 months",
    weight_range   = "2-85 kg",
    weight_median  = "6.7 kg",
    sex_female_pct = 29.0,
    race_ethnicity = "Not reported (single-centre French paediatric ICU).",
    disease_state  = paste(
      "Critically ill children with systemic arterial hypotension",
      "(septic shock n=16, non-traumatic cerebral injury n=6, heavy",
      "sedation n=8, severe congenital diaphragmatic defect n=8)",
      "requiring continuous IV norepinephrine for haemodynamic support.",
      sep = " "
    ),
    dose_range     = paste(
      "Continuous IV infusion 0.05-2 ug/kg/min (median 0.5) via central",
      "venous catheter; infusion duration median 1.5 days (range 1-13).",
      sep = " "
    ),
    regions        = "France (Hopital Necker Enfants-Malades, Paris).",
    notes          = paste(
      "Per Table 1 of Oualha 2014. 27 of 38 patients were male (71%); 14",
      "were malnourished (< -2 SD); 17 (45%) died during ICU stay; PELOD",
      "score median 31. ORG_FAIL_COUNT distribution: 14 of 38 with >3",
      "dysfunctions. Norepinephrine concentrations were measured at",
      "baseline before infusion (C0, n=22), at >=60 min into the infusion",
      "(C1, n=38), and 40 min after a flow-rate change or before the",
      "6-24 h cut-off (C2, n=16), for 76 PK observations; MAP was",
      "recorded at infusion onset and hourly thereafter (431 MAP",
      "observations).",
      sep = " "
    )
  )

  ini({
    # ---- PK structural parameters (Table 2 of Oualha 2014) -------------
    # theta_CL and theta_q0 are per-kg^0.75 typical values; both BW
    # exponents were estimated and then FIXED to 3/4 per the allometric
    # rule (Results paragraph 1 of "Norepinephrine pharmacokinetics";
    # Table 2 "fixed" flag).
    lcl <- log(6.6)        ; label("Typical unit clearance (L/h/kg^0.75)")              # Table 2: theta_CL = 6.6 L/h/kg, RSE 11%
    lq0 <- log(3.12)       ; label("Typical unit endogenous NorEp production rate (ug/h/kg^0.75)") # Table 2: theta_q0 = 3.12 ug/h/kg, RSE 23%
    e_wt_cl <- fixed(0.75) ; label("Body weight allometric exponent on CL (unitless, FIXED to 3/4)") # Table 2: theta_BW(CL) = 0.75 FIX
    e_wt_q0 <- fixed(0.75) ; label("Body weight allometric exponent on q0 (unitless, FIXED to 3/4)") # Table 2: theta_BW(q0) = 0.75 FIX

    # Vc was not separately estimated; ascribed to the circulating volume
    # V_circ = 0.08 * BW (L) per the rule cited from Linderkamp et al. in
    # Methods/PK section. Encoded as FIXED log(Vc/kg) plus a FIXED linear
    # weight exponent so the parameter is provenance-traceable. V (10 kg)
    # = 0.8 L is the paper's worked example (Table 2 footnote).
    lvc     <- fixed(log(0.08)) ; label("Central volume per kg (Vc / WT, L/kg, FIXED to circulating volume)") # Methods: V_Circ = 0.08 * BW; Table 2 footnote
    e_wt_vc <- fixed(1)         ; label("Body weight allometric exponent on Vc (unitless, FIXED to 1)")      # Same V_Circ rule

    # IIV on log(CL) and log(q0). Paper reports eta_CL = 0.6 and eta_q0 =
    # 1.1 as sqrt(omega^2) (Table 2); variances are 0.6^2 = 0.36 and
    # 1.1^2 = 1.21. No correlation between eta_CL and eta_q0 reported.
    etalcl ~ 0.36           # SD = 0.6  ; Table 2: eta_CL, RSE 14%
    etalq0 ~ 1.21           # SD = 1.1  ; Table 2: eta_q0, RSE 17%

    # ---- Hemodynamic Emax sub-model (Table 3 of Oualha 2014) -----------
    # Paper PK/PD modelling Eq.: MAP(t) = MAP0 + dMAP * Cc / (Cc + EC50)
    # where dMAP = MAPmax - MAP0.
    lmap0 <- log(34)        ; label("Typical basal MAP at PCA = 9 months (mmHg)")           # Table 3: MAP0 = 34 mmHg, RSE 5%
    ldmap <- log(32)        ; label("Typical drug-induced max MAP increase, org dysfn <=3 (mmHg)") # Table 3: dMAP = 32 mmHg (org dysfn 1-3), RSE 24%
    lec50 <- log(4.11)      ; label("NorEp concentration producing 50% of dMAP (ug/L)")     # Table 3: EC50_MAP = 4.11 ug/L, RSE 43%

    # Power-of-PCA covariate effect on basal MAP0:
    # MAP0_i = MAP0 * (PCA / 9)^theta_PCA, theta_PCA = 0.166 (Table 3).
    e_page_map0 <- 0.166    ; label("Postmenstrual age (months) power exponent on MAP0 (unitless)") # Table 3 / Results: theta_PCA = 0.166, RSE 19%

    # Categorical (binary) ORG_FAIL_COUNT effect on dMAP: 32 mmHg for org
    # dysfn <=3 (reference); 12 mmHg for org dysfn >=4. Encoded as an
    # additive log-shift on log(dMAP): e = log(12 / 32) = -0.981 applied
    # only when orgf_high = 1.
    e_orgf_high_dmap <- log(12 / 32) ; label("Additive log-shift on dMAP for ORG_FAIL_COUNT >= 4 (unitless)") # Table 3: dMAP 32 -> 12 across strata

    # IIV on log(MAP0) and log(dMAP). Table 3 reports sqrt(omega^2) =
    # 0.17 and 0.3, respectively; the Results paragraph 2 of
    # "Norepinephrine pharmacodynamics" gives the final-model BSVs as
    # 0.17 and 0.32 (Table 3 prints 0.3 = rounded 0.32). The Table 3
    # parenthetical "sqrt(omega^2 C50MAP)" on the dMAP IIV row is an
    # apparent transcription typo: no EC50 IIV is reported anywhere
    # else in the paper, while the dMAP IIV value matches the Results
    # paragraph 2 number.
    etalmap0 ~ 0.0289       # SD = 0.17 ; Table 3: eta_MAP0, RSE 15%
    etaldmap ~ 0.09         # SD = 0.30 ; Table 3: eta_dMAP, RSE 36%

    # ---- Residual error ------------------------------------------------
    propSd     <- 0.25      ; label("Proportional residual SD on NorEp concentration (fraction)") # Table 2: residual proportional = 0.25, RSE 17%
    propSd_MAP <- 0.14      ; label("Proportional residual SD on MAP (fraction)")                  # Table 3: residual proportional = 0.14, RSE 4%
  })

  model({
    # ---- Derived covariate terms ----------------------------------------
    # Binary high-organ-dysfunction indicator (>=4 vs reference <=3).
    orgf_high <- (ORG_FAIL_COUNT >= 4)

    # ---- Individual PK parameters ---------------------------------------
    # CL and q0 carry FIXED 0.75 exponent on WT (kg); Vc = 0.08 * WT (L)
    # with FIXED linear exponent. Time unit = hour.
    cl <- exp(lcl + etalcl) * WT^e_wt_cl
    q0 <- exp(lq0 + etalq0) * WT^e_wt_q0
    vc <- exp(lvc) * WT^e_wt_vc        # Vc = 0.08 * WT (L)

    # ---- Individual hemodynamic parameters ------------------------------
    # Power-of-PAGE effect on basal MAP0; additive log-shift on dMAP for
    # the high organ-dysfunction stratum.
    map0 <- exp(lmap0 + etalmap0) * (PAGE / 9)^e_page_map0
    dmap <- exp(ldmap + etaldmap + e_orgf_high_dmap * orgf_high)
    ec50 <- exp(lec50)

    # ---- Initial conditions ---------------------------------------------
    # central holds total (endogenous + exogenous) norepinephrine amount.
    # At baseline the endogenous production q0 balances elimination
    # cl * Cc, giving Cc(0) = q0 / cl and amount A(0) = q0 * vc / cl.
    central(0) <- q0 * vc / cl

    # ---- ODE system -----------------------------------------------------
    # Methods Eq. for total NorEp: dA/dt = q0 + Rate - cl * Cc. The
    # exogenous infusion Rate enters via the event table.
    d/dt(central) <- q0 - cl * central / vc

    # ---- Observations and residual error --------------------------------
    Cc  <- central / vc
    MAP <- map0 + dmap * Cc / (Cc + ec50)

    Cc  ~ prop(propSd)
    MAP ~ prop(propSd_MAP)
  })
}
