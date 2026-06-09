Oualha_2018_enoxaparin <- function() {
  description <- "Population PK model for subcutaneous enoxaparin in 22 children during the first post-operative week after paediatric liver transplantation (Oualha 2018). One-compartment open model with first-order absorption (ka fixed at 1/h) and first-order elimination, measured as anti-Xa activity (target 0.2-0.4 IU/mL). Apparent clearance CL/F is allometrically scaled by pre-operative bodyweight BWPREOP (fixed exponent 0.75); apparent central volume V/F is allometrically scaled (fixed exponent 1) by a time-varying post-operative bodyweight BW(t) that captures peri-operative fluid resuscitation followed by post-operative diuresis: BW(t) = (BWPREOP + PFA/1000) * (1 - (1 - fbw) * t^hill_bw / (tbw50^hill_bw + t^hill_bw)). Bodyweight-evolution parameters fbw / hill_bw / tbw50 are jointly estimated with the enoxaparin PK and carry their own between-subject variability."
  reference <- paste(
    "Oualha M, Chardot C, Debray D, Lesage F, Harroche A,",
    "Renolleau S, Treluyer J-M, Urien S (2018).",
    "Population pharmacokinetics of enoxaparin in early stage",
    "of paediatric liver transplantation.",
    "British Journal of Clinical Pharmacology 84(8):1736-1745.",
    "doi:10.1111/bcp.13543.",
    sep = " "
  )
  vignette <- "Oualha_2018_enoxaparin"
  units <- list(time = "hour", dosing = "IU", concentration = "IU/mL")

  covariateData <- list(
    WT = list(
      description        = "Pre-operative body weight (BWPREOP)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject at the pre-operative measurement. Drives the allometric scaling of CL/F with fixed exponent 0.75 against the 70 kg adult reference, and (in combination with PFA) sets the initial value of the time-varying post-operative bodyweight BW(t) computed inside model() which then drives the allometric scaling of V/F with fixed exponent 1. Cohort median 10.6 kg, range 6.7-34 kg (Table 1 of Oualha 2018).",
      source_name        = "BWPREOP"
    ),
    PFA = list(
      description        = "Perioperative intra-operative fluid administration volume",
      units              = "mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Total volume of saline / albumin / fresh-frozen-plasma / blood / platelet replacement administered intra-operatively, summed over the entire liver-transplant procedure. Time-fixed per subject. Enters the time-varying post-operative bodyweight curve BW(t) algebraically as (BWPREOP + PFA/1000), with the /1000 conversion applying the conventional 1 mL ~ 1 g fluid-density-1 mapping so that PFA in mL maps onto a kg increment of body weight. Cohort median 2634 mL, range 1008-6520 mL (Table 1 of Oualha 2018).",
      source_name        = "PFA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 22L,
    n_studies      = 1L,
    age_range      = "5-154 months (paediatric; 0.4-12.8 years)",
    age_median     = "21.5 months",
    weight_range   = "6.7-34 kg (BWPREOP)",
    weight_median  = "10.6 kg (BWPREOP)",
    sex_female_pct = 63.6,
    race_ethnicity = "Not reported (single-centre French paediatric intensive-care cohort).",
    disease_state  = "Paediatric liver-transplant recipients in the first post-operative week. Indication for transplantation: biliary cirrhosis (n=20), metabolic disease (n=1), tumour (n=1). Left split liver graft n=18 (81.8%); median percent of graft weight to BWPREOP 3.4% (range 1.2-8.7).",
    dose_range     = "Initial subcutaneous 50 IU/kg every 12 hours (range 43-59 IU/kg at initiation), with empirical dose adjustments up to 108 IU/kg every 12 hours based on anti-Xa activity therapeutic-drug-monitoring.",
    regions        = "France (Hopital Necker Enfants-Malades, Paris, paediatric intensive care unit).",
    notes          = "Per Table 1 of Oualha 2018. 19 of 22 patients required perioperative vasopressor support (norepinephrine). Acute renal dysfunction (transient) in 7 patients (31.8%); baseline creatinine 25 umol/L (range 13-100); baseline creatinine clearance 143.3 mL/min/1.73 m^2 (range 79-308). Median perioperative fluid administration 2634 mL (range 1008-6520 mL). 136 anti-Xa activity samples collected within 10-156 h after enoxaparin initiation; 37 samples below the 0.1 IU/mL limit of quantification (Monolix M3 censored handling). Target therapeutic anti-Xa activity 0.2-0.4 IU/mL was reached in only 23 of 74 steady-state samples (31%); all 22 children experienced at least one sub-therapeutic exposure during the first post-operative week."
  )

  ini({
    # ---- Structural PK parameters (Table 3 of Oualha 2018) -----------------
    # Apparent ka, CL/F and V/F at the 70 kg adult allometric reference.
    # Anti-Xa activity units throughout (IU dose, IU/mL concentration); V
    # in L, CL in L/h.
    lka <- fixed(log(1))     ; label("First-order absorption rate (1/h, FIXED)")             # Table 3: ka = 1 (fixed)
    lcl <- log(1.23)         ; label("Apparent clearance CL/F at 70 kg BWPREOP (L/h)")       # Table 3: CL_TYP = 1.23 L/h for 70 kg BWPREOP (RSE 15%)
    lvc <- log(14.6)         ; label("Apparent central volume V/F at 70 kg and fbw=0 (L)")   # Table 3: V_TYP = 14.6 L for 70 kg BWPREOP and f_BW = 0 (RSE 27%)

    # Allometric exponents -- FIXED per Table 3 (theta_BW(CL) = 3/4, theta_BW(V) = 1).
    e_wt_cl <- fixed(0.75)   ; label("BWPREOP allometric exponent on CL/F (unitless, FIXED to 3/4)") # Table 3: theta_BW(CL) = 0.75 (fixed)
    e_wt_vc <- fixed(1)      ; label("BW(t) allometric exponent on V/F (unitless, FIXED to 1)")     # Table 3: theta_BW(V) = 1 (fixed)

    # ---- BW(t) evolution sub-model (Table 3 of Oualha 2018) ----------------
    # Post-operative bodyweight curve: BW(t) = (BWPREOP + PFA/1000) *
    # (1 - (1 - fbw) * t^hill_bw / (tbw50^hill_bw + t^hill_bw)).
    # Encoded under the canonical lfbw / lhill_bw / ltbw50 family registered
    # alongside this paper in parameter-names.md section "Body-weight evolution".
    lfbw     <- log(0.89)    ; label("Asymptotic fractional body-weight retention after fluid loss (unitless)") # Table 3: theta_fBW = 0.89 (RSE 1%)
    lhill_bw <- log(4.43)    ; label("Hill steepness of the BW(t) decline curve (unitless)")                    # Table 3: theta_Hill = 4.43 (RSE 7%)
    ltbw50   <- log(61.3)    ; label("Time to 50% of the post-operative BW loss (h)")                           # Table 3: theta_tBW50 = 61.3 h (RSE 18%)

    # ---- Inter-individual variability (Table 3) ----------------------------
    # Monolix exponential BSV model with omega reported as sqrt(omega^2).
    # No correlations reported; encoded as independent etas.
    etalcl     ~ 0.3969      # eta_CL  : omega = 0.63 (shrinkage 10.0%); variance = 0.63^2
    etalvc     ~ 1.5129      # eta_V   : omega = 1.23 (shrinkage 20.6%); variance = 1.23^2
    etalfbw    ~ 0.0036      # eta_FBW : omega = 0.06 (shrinkage 4.08%); variance = 0.06^2
    etaltbw50  ~ 0.64        # eta_tBW50: omega = 0.80 (shrinkage 4.60%); variance = 0.80^2

    # ---- Residual error (Table 3) ------------------------------------------
    # Proportional residual on anti-Xa activity; the additive 0.054 kg and
    # proportional 0.012 BW residuals reported in Table 3 belong to the
    # joint estimation of the BW(t) curve against measured per-day body
    # weights and are not relevant when this model is used purely for
    # forward simulation of anti-Xa exposure (BW(t) is then a deterministic
    # output of the structural curve).
    propSd <- 0.42           ; label("Proportional residual SD on anti-Xa activity (fraction)") # Table 3: proportional on anti-Xa = 0.42 (RSE 3%)
  })

  model({
    # ---- Time-varying post-operative bodyweight BW(t) ----------------------
    # Hill-type saturation curve combining the immediate fluid loading
    # (BWPREOP + PFA/1000, with PFA in mL converted to kg) with a
    # post-operative diuresis governed by the fractional retention fbw,
    # the Hill steepness hill_bw, and the half-time tbw50 (all in hours).
    fbw     <- exp(lfbw     + etalfbw)
    hill_bw <- exp(lhill_bw)
    tbw50   <- exp(ltbw50   + etaltbw50)
    bwt     <- (WT + PFA / 1000) * (1 - (1 - fbw) * t^hill_bw / (tbw50^hill_bw + t^hill_bw))

    # ---- Individual PK parameters ------------------------------------------
    # ka is fixed; CL scaled allometrically by static BWPREOP (=WT);
    # V scaled allometrically by the time-varying BW(t).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT  / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (bwt / 70)^e_wt_vc

    # ---- ODE system --------------------------------------------------------
    # depot holds absorbed dose (IU); central holds circulating drug (IU);
    # Cc converts IU/L (= central / vc) to the published IU/mL unit by
    # dividing by 1000.
    kel <- cl / vc
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ---- Observation and residual error ------------------------------------
    Cc <- central / vc / 1000
    Cc ~ prop(propSd)
  })
}
