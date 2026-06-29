Su_2016_dexmedetomidine <- function() {
  description <- "Two-compartment IV population PK model for dexmedetomidine in neonates (1 day-1 month) and infants (1-24 months) after open heart surgery, with a priori allometric weight scaling on CL, Q (exponent 0.75) and V1, V2 (exponent 1) at a 70 kg reference weight; an Emax-form postnatal-age maturation on CL with TM50 = 0.032 months; a power-form effect of total cardiopulmonary bypass duration on CL centred at 60 min; and a 1.24-fold multiplicative increase in CL for patients with right-to-left intracardiac shunt (Qp:Qs < 1) (Su 2016 Table 4, allometric weight-normalized full covariate model)"
  reference <- paste(
    "Su F, Gastonguay MR, Nicolson SC, DiLiberto M,",
    "Ocampo-Pelland A, Zuppa AF.",
    "Dexmedetomidine pharmacology in neonates and infants after",
    "open heart surgery.",
    "Anesth Analg. 2016;122(5):1556-1566.",
    "doi:10.1213/ANE.0000000000000869",
    sep = " "
  )
  vignette <- "Su_2016_dexmedetomidine"
  units <- list(time = "h", dosing = "ug", concentration = "pg/mL") # Methods: dose ug (loading dose 0.25-1 ug/kg over 10 min then CIVI 0.2-0.75 ug/kg/h); validated HPLC-MS/MS assay reported in pg/mL with LLOQ 5 pg/mL. Final model parameterised with CL, Q in mL/min and V1, V2 in L (Table 4); CL and Q converted to L/h here (657 mL/min = 39.42 L/h; 6780 mL/min = 406.8 L/h). With central in ug and vc in L, central/vc gives ug/L = ng/mL; the observation Cc multiplies by 1000 to convert to pg/mL.

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject (per-subject weight at study entry, ranging 2.3-11.9 kg across the 59 evaluable subjects; median 5.97 kg per Table 3). A priori allometric scaling per Su 2016 Methods 'Full Covariate Model': CL and Q scale as (WT/70)^0.75 and V1 and V2 scale as (WT/70)^1, with reference weight 70 kg chosen for comparison with adult populations.",
      source_name        = "WT"
    ),
    PNA = list(
      description        = "Postnatal age",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject (per-subject postnatal age at study entry, ranging 0.03 months / 1 day-20.4 months across the 59 evaluable subjects; median 4.3 months per Table 3). Su 2016 Methods 'Full Covariate Model' enters age as an Emax-form maturation on CL: CL_maturation = Age / (TM50 + Age) with TM50 = 0.032 months (~1 day). The model is silent on whether the covariate is gestational, postmenstrual, or postnatal -- the cohort is full-term neonates and infants so postmenstrual age equals (40 weeks + PNA) and the maturation effect is anchored at chronological birth (PNA = 0). The covariate column is therefore mapped to the canonical PNA. At PNA = 0 the maturation multiplier evaluates to 0; for prospective simulation the youngest subject is set at PNA = 0.03 months (the cohort minimum) so the multiplier is well-defined.",
      source_name        = "Age"
    ),
    T_CPB = list(
      description        = "Total cardiopulmonary bypass time during the immediately preceding open heart surgery",
      units              = "minutes",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject (single per-subject scalar set at surgery close; cohort median 63 min, range 16-169 min per Table 3). Reference 60 minutes (the cohort median, used in the Abstract and Discussion as the typical-patient anchor). Su 2016 Methods 'Full Covariate Model': power-form effect on CL, (T_CPB / 60)^(-0.31), so longer bypass time reduces post-operative CL (negative exponent).",
      source_name        = "TBYP"
    ),
    ICSHUNT_R2L = list(
      description        = "Right-to-left intracardiac shunt indicator (Qp:Qs < 1)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no right-to-left shunt; Qp:Qs >= 1)",
      notes              = "Time-fixed per subject. 1 = patient cardiac anatomy produces right-to-left intracardiac shunting with pulmonary-blood-flow-to-systemic-blood-flow ratio Qp:Qs < 1; 0 = otherwise. Most common shunt anatomy in this cohort was single-ventricle physiology after stage 2 palliation (Glenn or hemi-Fontan procedure). 19 of 59 subjects (32%) had a right-to-left shunt per Table 3. Su 2016 Methods 'Full Covariate Model': multiplicative effect on CL of 1.24 when ICSHUNT_R2L = 1, hypothesised mechanism is shunt-induced reduction in first-pass lung extraction combined with increased hepatic perfusion via increased systemic blood flow.",
      source_name        = "intracardiac shunt"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 59L,
    n_studies       = 1L,
    age_range       = "0.03-20.4 months (1 day-20.4 months postnatal)",
    age_median      = "4.3 months",
    weight_range    = "2.3-11.9 kg",
    weight_median   = "5.97 kg",
    sex_female_pct  = 45.8,                                # 27/59 female per Table 3
    race_ethnicity  = NULL,                                # Not separately tabulated in Table 3
    disease_state   = "Full-term neonates (1 day-1 month, n = 23) and infants (1-24 months, n = 36) with congenital heart disease requiring mechanical ventilation after open heart surgery; adequate hepatic and renal function; no evidence of heart block. The most common shunt anatomy in shunt-positive subjects was single-ventricle physiology after stage 2 palliation (Glenn or hemi-Fontan procedure).",
    dose_range      = "Loading dose 0.25-1 ug/kg administered IV over 10 minutes followed by a continuous IV infusion (CIVI) of 0.2-0.75 ug/kg/h for up to 24 hours (Table 1; cohort dependent). Median duration of dexmedetomidine administration 10.1 hours (range 3.2-24.3 h).",
    regions         = "USA (single-centre dose-escalation trial conducted at The Children's Hospital of Philadelphia under FDA IND #69,758).",
    cpb_time        = "Median 63 min (range 16-169 min; Table 3)",
    intracardiac_shunt_pct = "32% (19/59 with Qp:Qs < 1)",
    notes           = "Dose-escalation cohorts: infants 1 (0.35 ug/kg + 0.25 ug/kg/h, n = 12), 2 (0.7 + 0.5, n = 12), 3 (1 + 0.75, n = 12); neonates 4 (0.25 + 0.2, n = 9), 4A (0.35 + 0.3, n = 9), 5 (0.5 + 0.4, n = 5; stopped early at MTD). NONMEM v6 level 2.0 with ADVAN 3 TRANS 4 and FOCE-I. PK samples below the LLOQ of 5 pg/mL were excluded (3% of neonatal, 8% of infant samples)."
  )

  ini({
    # ============================================================
    # Structural PK parameters - reference 70 kg adult-equivalent
    # (Su 2016 Table 4, allometric weight-normalized model)
    # CL and Q reported in mL/min; converted to L/h (* 60 / 1000)
    # for natural pairing with V in L and time in hours.
    # ============================================================
    lcl <- log(39.42)  ; label("Clearance at 70 kg (L/h)")                            # Table 4 (allometric): theta_CL = 657 mL/min/70kg (LLP 95% CI 600-750; SE% 6.24); 657 mL/min * 60 / 1000 = 39.42 L/h
    lvc <- log(88)     ; label("Central volume V1 at 70 kg (L)")                      # Table 4: theta_V1 = 88 L/70kg (LLP 95% CI 70-110; SE% 11.92)
    lq  <- log(406.8)  ; label("Intercompartmental clearance Q at 70 kg (L/h)")       # Table 4: theta_Q = 6780 mL/min/70kg (LLP 95% CI 4000-20000; SE% 40.71); 6780 mL/min * 60 / 1000 = 406.8 L/h
    lvp <- log(112)    ; label("Peripheral volume V2 at 70 kg (L)")                   # Table 4: theta_V2 = 112 L/70kg (LLP 95% CI 90-130; SE% 9.73)

    # Allometric exponents - fixed a priori per Su 2016 Methods 'Full Covariate Model'
    # (theoretical allometric exponents, Anderson and Holford)
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent on CL (unitless, FIXED)")     # Su 2016 Methods: theta_allometric = 0.75 for clearances, fixed
    e_wt_vc <- fixed(1)    ; label("Allometric exponent on V1 (unitless, FIXED)")     # Su 2016 Methods: exponent = 1 for volumes, fixed
    e_wt_q  <- fixed(0.75) ; label("Allometric exponent on Q  (unitless, FIXED)")     # Su 2016 Methods: theta_allometric = 0.75 for clearances, fixed
    e_wt_vp <- fixed(1)    ; label("Allometric exponent on V2 (unitless, FIXED)")     # Su 2016 Methods: exponent = 1 for volumes, fixed

    # Covariate effects on CL (Su 2016 Table 4 footnote, allometric model)
    tm50_cl     <- log(0.032) ; label("Postnatal age at which CL reaches 50% of mature value (months)")  # Table 4: Age 50% CL = 0.032 months (LLP 95% CI 0.015-0.055; SE% 34.16); log-transformed to keep estimates positive
    e_tcpb_cl   <- -0.31      ; label("Power exponent of total CPB time on CL (unitless)")               # Table 4: TBYP effect on CL = -0.31 (LLP 95% CI -0.45 to -0.15; SE% -23.48)
    e_icshunt_cl <- 1.24      ; label("Multiplicative factor on CL when right-to-left intracardiac shunt is present (unitless)")  # Table 4: Intracardiac shunting CL = 1.24 (LLP 95% CI 1.1-1.5; SE% 7.49)

    # ============================================================
    # Inter-individual variability (Su 2016 Table 4, allometric model)
    # Table footer: "Interindividual variability ... point estimates
    # are presented as percent coefficient of variation (square root
    # of variance) x 100"; i.e., the variance omega^2 has sqrt(omega^2)
    # * 100 = CV%, so omega^2 = (CV% / 100)^2 directly. Covariances
    # are reported as point estimates (with correlations in parentheses
    # for cross-checking); the correlations verify against the diagonal
    # variances to within rounding (e.g. 0.116 / sqrt(0.0817 * 0.387)
    # = 0.652, matching the reported 0.65).
    # Block structure: 3x3 block on CL, V1, Q (Table 4 covariance rows
    # Cov CL V1, Cov CL Q, Cov V1 Q); diagonal IIV on V2.
    # ============================================================
    etalcl + etalvc + etalq ~ c(0.08168,
                                0.116, 0.38701,
                                0.165, 0.775, 2.46993)        # Table 4: CV% 28.58 (CL) -> var 0.08168; CV% 62.21 (V1) -> var 0.38701; CV% 157.16 (Q) -> var 2.46993; Cov CL V1 = 0.116 (corr 0.65); Cov CL Q = 0.165 (corr 0.37); Cov V1 Q = 0.775 (corr 0.79)
    etalvp                  ~ 0.08202                          # Table 4: CV% 28.64 (V2) -> var 0.08202 (SE% 39.88); no covariance with the CL/V1/Q block reported

    # ============================================================
    # Residual error - combined additive + proportional (Su 2016
    # Methods 'Base Model' and Table 4)
    # ============================================================
    propSd <- 0.1975 ; label("Proportional residual SD on Cc (fraction)")              # Table 4: sigma^2_proportional reported as CV% = 19.75 -> SD = 0.1975 (SE% 12.67); footer "proportional residual variability point estimates are presented as percent coefficient of variation (square root of variance) x 100"
    addSd  <- 3.30   ; label("Additive residual SD on Cc (pg/mL)")                     # Table 4: sigma^2_additive = 3.30 (SE% 29.91); footer "sigma^2_additive point estimate is expressed as a standard deviation"
  })

  model({
    # 1. Derived terms for covariate effects on CL
    #    Maturation:                 Age / (TM50 + Age), Emax-form
    #    CPB-time power effect:      (T_CPB / 60)^(-0.31)
    #    Right-to-left shunt factor: 1.24 when ICSHUNT_R2L = 1, else 1
    tm50          <- exp(tm50_cl)
    maturation_cl <- PNA / (tm50 + PNA)
    tcpb_cl       <- (T_CPB / 60)^e_tcpb_cl
    icshunt_cl    <- e_icshunt_cl^ICSHUNT_R2L

    # 2. Individual PK parameters with a priori allometric weight scaling (reference 70 kg)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * maturation_cl * tcpb_cl * icshunt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q  <- exp(lq  + etalq)  * (WT / 70)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vp

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 4. Two-compartment IV disposition (dose into central as IV infusion)
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-   k12 * central        - k21 * peripheral1

    # 5. Plasma DEX concentration: central in ug, vc in L gives ug/L = ng/mL;
    #    multiply by 1000 to convert to pg/mL (Su 2016 assay reporting unit).
    Cc <- central / vc * 1000

    # 6. Observation model (combined additive + proportional residual error)
    Cc ~ add(addSd) + prop(propSd)
  })
}
