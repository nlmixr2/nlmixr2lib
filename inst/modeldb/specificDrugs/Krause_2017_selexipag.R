Krause_2017_selexipag <- function() {
  description <- "Joint two-compartment parent + two-compartment metabolite population PK model for oral selexipag and its active metabolite ACT-333679 in adults with pulmonary arterial hypertension (Krause 2017, GRIPHON study). First-order absorption with a fixed 0.668 h absorption lag delivers selexipag into a two-compartment disposition with linear total clearance CL/F (apparent total clearance, of which the rate constant kmet describes the fraction converted to ACT-333679); the metabolite has its own two-compartment disposition with first-order elimination via km. Body weight (allometric on V_p/F and CL/F; on V_m/F), total bilirubin (power on CL/F), sex (multiplicative on km), and a four-level PAH-comedication categorical (naive / ERA only / PDE5 inhibitor only / ERA + PDE5 combined; multiplicative on km) were retained as statistically significant covariates."
  reference   <- "Krause A, Machacek M, Lott D, Hurst N, Bruderer S, Dingemanse J. Population modeling of selexipag pharmacokinetics and clinical response parameters in patients with pulmonary arterial hypertension. CPT Pharmacometrics Syst Pharmacol. 2017;6(7):477-485. doi:10.1002/psp4.12202"
  vignette    <- "Krause_2017_selexipag"
  units       <- list(time = "hour", dosing = "ug", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight at baseline.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effects on selexipag V_p/F (exponent 1.20), selexipag CL/F (exponent 0.61), and ACT-333679 V_m/F (exponent 0.88); reference body weight 70 kg per Krause 2017 Methods ('typical covariate values were derived from the data using a round value close to mean and median... 70 kg body weight').",
      source_name        = "WT"
    ),
    TBILI = list(
      description        = "Total serum bilirubin at baseline (umol/L).",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form negative effect on selexipag CL/F (exponent -0.40); reference total bilirubin 10 umol/L (Krause 2017 Methods). Hepatic-function marker; higher bilirubin reduces selexipag CL.",
      source_name        = "BILI"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female) in the source paper -- Krause 2017 takes females as the reference category for the categorical sex effect on the ACT-333679 elimination rate constant k_m.",
      notes              = "Krause 2017 reports a +15% effect of male sex on k_m (Table 1, b = 0.15). To preserve the paper's female-reference parameterisation while using the canonical SEXF column (1 = female), the model() block applies the effect via the male indicator (1 - SEXF): k_m *= (1 + 0.15 * (1 - SEXF)). Males therefore have ~15% higher k_m (faster ACT-333679 elimination) than females, which translates to ~13% lower ACT-333679 AUCss (Krause 2017 Table 3).",
      source_name        = "SEX"
    ),
    CONMED_ERA = list(
      description        = "Concomitant endothelin-receptor-antagonist (ERA) monotherapy indicator (1 = on an ERA but not on a PDE5 inhibitor, 0 = otherwise). One of three orthogonal mutually-exclusive indicators decomposing a four-level PAH-comedication categorical with PAH-comedication-naive as the reference (all three indicators = 0).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the PAH-comedication-naive stratum; mutually exclusive with CONMED_PDE5I and CONMED_ERA_PDE5I).",
      notes              = "Multiplicative effect (1 + 0.15) on the ACT-333679 elimination rate constant k_m (Krause 2017 Table 1; +15% relative to naive, translating to ~13% lower ACT-333679 AUCss). ERA agents in the GRIPHON cohort were bosentan, macitentan, and ambrisentan.",
      source_name        = "PAHCOMED category 'ERA' (decomposed from a four-level categorical {naive, ERA, PDE5, ERA+PDE5})"
    ),
    CONMED_PDE5I = list(
      description        = "Concomitant phosphodiesterase type 5 inhibitor (PDE5I) monotherapy indicator (1 = on a PDE5I but not on an ERA, 0 = otherwise). One of three orthogonal mutually-exclusive indicators decomposing a four-level PAH-comedication categorical with PAH-comedication-naive as the reference.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the PAH-comedication-naive stratum; mutually exclusive with CONMED_ERA and CONMED_ERA_PDE5I).",
      notes              = "Multiplicative effect (1 + 0.07) on the ACT-333679 elimination rate constant k_m (Krause 2017 Table 1; +7% relative to naive). The PDE5I-only coefficient is statistically not significant (p = 0.19) but is retained in the final model because the four-level categorical encoding is kept intact. PDE5I agents in the GRIPHON cohort were sildenafil and tadalafil.",
      source_name        = "PAHCOMED category 'PDE5 inh.' (decomposed from a four-level categorical {naive, ERA, PDE5, ERA+PDE5})"
    ),
    CONMED_ERA_PDE5I = list(
      description        = "Concomitant ERA + PDE5-inhibitor combination indicator (1 = on both an ERA and a PDE5 inhibitor, 0 = otherwise). One of three orthogonal mutually-exclusive indicators decomposing a four-level PAH-comedication categorical with PAH-comedication-naive as the reference.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the PAH-comedication-naive stratum; mutually exclusive with CONMED_ERA and CONMED_PDE5I).",
      notes              = "Multiplicative effect (1 + 0.37) on the ACT-333679 elimination rate constant k_m (Krause 2017 Table 1; +37% relative to naive, translating to a 30% reduction in ACT-333679 AUCss for this stratum vs naive per Table 3). The combined-stratum coefficient is its own categorical level rather than the multiplicative product of the ERA-only and PDE5I-only effects -- the four-level PAH-comedication taxonomy was kept intact for parameterisation.",
      source_name        = "PAHCOMED category 'ERA and PDE5 inh.' (decomposed from a four-level categorical {naive, ERA, PDE5, ERA+PDE5})"
    )
  )

  population <- list(
    n_subjects     = 512,
    n_studies      = 1,
    age_range      = "Adult PAH patients; reference 40 years (Krause 2017 Methods 'typical covariate values'). Full age summary in GRIPHON main paper (Sitbon 2015 NEJM).",
    weight_range   = "Reference 70 kg; covariate exploration spanned 51-96 kg (Krause 2017 Table 3 reference vs comparison patients).",
    sex_female_pct = "Female-majority cohort consistent with the PAH epidemiology; per-subject sex was retained as a covariate on k_m. Exact percentage not reported in the popPK paper (source: Sitbon 2015 NEJM GRIPHON main paper).",
    race_ethnicity = "Tested as a covariate (reference category Asian) but not retained as statistically significant.",
    disease_state  = "Pulmonary arterial hypertension (idiopathic, drug- or toxin-induced, associated with connective tissue disease, congenital heart disease, or other aetiologies) on stable background PAH comedication (endothelin-receptor antagonist and/or phosphodiesterase type 5 inhibitor) or PAH-comedication-naive at study entry.",
    dose_range     = "Selexipag 200 ug b.i.d. up-titrated weekly in 200 ug increments to the individual maximum tolerated dose (IMD), with maximum allowed 1,600 ug b.i.d. Doses 200-1,600 ug b.i.d. were observed in the population-PK cohort; 68/559 patients (12%) had IMD = 200 ug and 163/559 (29%) had IMD = 1,600 ug.",
    regions        = "Multinational phase III GRIPHON study, 181 centres globally (NCT01106014).",
    n_observations = "5,068 quantifiable post-baseline plasma concentrations (after exclusion of 2,034 records with missing values for concentration/sampling/dosing/dose). Below-quantification samples imputed by sampling from the conditional distribution within Monolix SAEM (351 LLOQ selexipag, 58 LLOQ ACT-333679).",
    co_medication  = "PAH-specific comedication was a covariate (four-level taxonomy {naive, ERA only, PDE5 inhibitor only, ERA + PDE5}; reference = naive). Other concomitant medications including digoxin and CYP3A4/CYP2C8 inhibitors were tested but not retained (Krause 2017 Tables 2-3).",
    notes          = "Population PK fit by Monolix SAEM with diagonal random-effects covariance matrix and Monolix 'auto' iteration count. Five doses prior to each PK sample were carried in the dataset to reach steady state. The population-typical absorption lag time t_lag was fixed to the healthy-subject popPK estimate (Bruderer 2014 Pharmacology / Hoch 2015 thorough-QT analysis) because GRIPHON's sparse sampling could not estimate it robustly; individual lag-time estimates were allowed to vary between 0-2 h via a Monolix logit-bounded transformation (Krause 2017 Methods)."
  )

  ini({
    # ---------------- Parent (selexipag) structural parameters ----------------
    # Apparent parameters because absolute bioavailability of selexipag is
    # unknown; all volumes and clearances carry an implicit "/F". Reference
    # covariate values: 70 kg body weight, 10 umol/L total bilirubin
    # (Krause 2017 Methods 'typical covariate values'); reference sex = female
    # and reference PAH comedication = naive (Krause 2017 Methods 'Covariate
    # selection (PK)').
    lka   <- log(0.71);  label("First-order absorption rate constant ka (1/h)")                                # Krause 2017 Table 1: ka = 0.71 1/h, RSE 5%
    lcl   <- log(19.10); label("Apparent selexipag clearance CL/F at 70 kg / 10 umol-L bilirubin (L/h)")       # Krause 2017 Table 1: CL/F = 19.10 L/h, RSE 8%
    lvc   <- log(12.90); label("Apparent selexipag central volume V_p/F at 70 kg (L)")                         # Krause 2017 Table 1: V_p/F = 12.90 L, RSE 16%
    lk12  <- log(0.09);  label("Selexipag central-to-peripheral rate constant k_12 (1/h)")                     # Krause 2017 Table 1: k_12 = 0.09 1/h, RSE 18%
    lk21  <- log(0.06);  label("Selexipag peripheral-to-central rate constant k_21 (1/h)")                     # Krause 2017 Table 1: k_21 = 0.06 1/h, RSE 17%
    ltlag <- fixed(log(0.668)); label("Selexipag absorption lag time t_lag, fixed (h)")                        # Krause 2017 Table 1 and Methods 'Population PK': t_lag = 0.668 h, fixed to the estimate from the healthy subject popPK model

    # ---------------- Metabolite (ACT-333679) structural parameters -----------
    # Apparent V_m/F is the metabolite central volume; k_m, k_34, k_43 are the
    # metabolite elimination and inter-compartmental rate constants; k_met is
    # the parent->metabolite formation rate constant (independent of CL_parent;
    # fraction-metabolised is implicit at k_met * V_c / CL_parent = ~45%).
    lvc_act  <- log(4.65); label("Apparent ACT-333679 central volume V_m/F at 70 kg (L)")                      # Krause 2017 Table 1: V_m/F = 4.65 L, RSE 17%
    lkm_act  <- log(0.49); label("ACT-333679 elimination rate constant k_m at female / PAH-naive reference (1/h)") # Krause 2017 Table 1: k_m = 0.49 1/h, RSE 16% (reference female on no PAH comedication)
    lk34_act <- log(1.04); label("ACT-333679 central-to-peripheral rate constant k_34 (1/h)")                  # Krause 2017 Table 1: k_34 = 1.04 1/h, RSE 22%
    lk43_act <- log(0.18); label("ACT-333679 peripheral-to-central rate constant k_43 (1/h)")                  # Krause 2017 Table 1: k_43 = 0.18 1/h, RSE 14%
    lkmet    <- log(0.67); label("Selexipag-to-ACT-333679 metabolite-formation rate constant k_met (1/h)")     # Krause 2017 Table 1: k_met = 0.67 1/h, RSE 18%

    # ---------------- Covariate effects ---------------------------------------
    # Power-form on continuous covariates and multiplicative-from-1 categorical
    # form (Krause 2017 Methods 'Covariate selection (PK)'). The four-level
    # PAH-comedication categorical is decomposed into three orthogonal
    # mutually-exclusive binary indicators (CONMED_ERA, CONMED_PDE5I,
    # CONMED_ERA_PDE5I); the PAH-naive stratum is the reference (all three
    # indicators = 0).
    e_wt_vc        <-  1.20;  label("Allometric exponent on selexipag V_p/F (unitless)")                       # Krause 2017 Table 1: body weight on V_p/F, b = 1.20, RSE 25%
    e_wt_cl        <-  0.61;  label("Allometric exponent on selexipag CL/F (unitless)")                        # Krause 2017 Table 1: body weight on CL/F, b = 0.61, RSE 25%
    e_tbili_cl     <- -0.40;  label("Power exponent of total bilirubin on selexipag CL/F (unitless)")          # Krause 2017 Table 1: total bilirubin on CL/F, b = -0.40, RSE 18%
    e_wt_vc_act    <-  0.88;  label("Allometric exponent on ACT-333679 V_m/F (unitless)")                      # Krause 2017 Table 1: body weight on V_m/F, b = 0.88, RSE 21%
    e_sex_km_act   <-  0.15;  label("Multiplicative effect of male sex on ACT-333679 k_m (unitless)")          # Krause 2017 Table 1: sex on k_m, b = 0.15, RSE 31% (males have +15% k_m vs female reference)
    e_era_km_act   <-  0.15;  label("Multiplicative effect of PAH comedication ERA only on ACT-333679 k_m")    # Krause 2017 Table 1: PAH comedication (ERA) on k_m, b = 0.15, RSE 38%
    e_pde5_km_act  <-  0.07;  label("Multiplicative effect of PAH comedication PDE5I only on ACT-333679 k_m")  # Krause 2017 Table 1: PAH comedication (PDE5 inh.) on k_m, b = 0.07, RSE 77%, p = 0.19 (NS but retained because the four-level categorical was kept intact)
    e_combo_km_act <-  0.37;  label("Multiplicative effect of PAH comedication ERA + PDE5I on ACT-333679 k_m") # Krause 2017 Table 1: PAH comedication (ERA and PDE5 inh.) on k_m, b = 0.37, RSE 14%

    # ---------------- IIV (variances on log scale) ----------------------------
    # Krause 2017 reports the IIV standard deviation (omega) on the log scale
    # for each parameter; ini() takes variances, so each entry below is
    # omega^2 with the source-table omega quoted in the trailing comment.
    # Rate-constant IIVs (etalk12, etalk21, etalk34_act, etalk43_act) are
    # carried in their paper-faithful form rather than re-parameterised into
    # Q / V_p / Q_act / V_p_act so the source IIV structure is preserved
    # without cross-parameter covariance translation.
    etalka      ~ 0.1521  # 0.39^2; Krause 2017 Table 1: omega(ka) = 0.39, RSE 12%
    etalcl      ~ 0.5329  # 0.73^2; Krause 2017 Table 1: omega(CL/F) = 0.73, RSE 4%
    etalvc      ~ 0.0961  # 0.31^2; Krause 2017 Table 1: omega(V_p/F) = 0.31, RSE 38%
    etalk12     ~ 0.0625  # 0.25^2; Krause 2017 Table 1: omega(k_12) = 0.25, RSE 73%
    etalk21     ~ 1.1236  # 1.06^2; Krause 2017 Table 1: omega(k_21) = 1.06, RSE 12%
    etalvc_act  ~ 0.0100  # 0.10^2; Krause 2017 Table 1: omega(V_m/F) = 0.10, RSE 241% (small IIV; kept estimated to support SAEM convergence per the paper)
    etalkm_act  ~ 0.0729  # 0.27^2; Krause 2017 Table 1: omega(k_m) = 0.27, RSE 20%
    etalk34_act ~ 0.2209  # 0.47^2; Krause 2017 Table 1: omega(k_34) = 0.47, RSE 39%
    etalk43_act ~ 0.7921  # 0.89^2; Krause 2017 Table 1: omega(k_43) = 0.89, RSE 16%
    etalkmet    ~ 0.0025  # 0.05^2; Krause 2017 Table 1: omega(k_met) = 0.05, RSE 788% (small IIV; kept estimated to support SAEM convergence per the paper)

    # ---------------- Residual error (proportional, by output) ----------------
    propSd     <- 0.75; label("Selexipag proportional residual SD (fraction)")                                 # Krause 2017 Table 1: b_1 (selexipag) = 0.75, RSE 2%
    propSd_act <- 0.49; label("ACT-333679 proportional residual SD (fraction)")                                # Krause 2017 Table 1: b_2 (ACT-333679) = 0.49, RSE 2%
  })

  model({
    # Reference covariate values for the typical patient (Krause 2017 Methods
    # 'Covariate selection (PK)' and Methods 'typical covariate values').
    ref_wt    <- 70    # kg
    ref_tbili <- 10    # umol/L

    # ----- Individual parameters: selexipag (parent) ---------------------------
    ka  <- exp(lka  + etalka)
    cl  <- exp(lcl  + etalcl) * (WT / ref_wt)^e_wt_cl * (TBILI / ref_tbili)^e_tbili_cl
    vc  <- exp(lvc  + etalvc) * (WT / ref_wt)^e_wt_vc
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    tlag <- exp(ltlag)

    # ----- Individual parameters: ACT-333679 (metabolite) ----------------------
    # Sex effect is on the male indicator (1 - SEXF) because the source paper
    # defines females as the reference category for the categorical sex
    # contrast on k_m. PAH-comedication effects are mutually exclusive on the
    # three non-reference strata (Krause 2017 Methods 'Covariate selection
    # (PK)' and Table 2 footnote: each category has its own b coefficient).
    vc_act  <- exp(lvc_act + etalvc_act) * (WT / ref_wt)^e_wt_vc_act
    km_act  <- exp(lkm_act + etalkm_act) *
               (1 + e_sex_km_act * (1 - SEXF)) *
               (1 + e_era_km_act   * CONMED_ERA +
                    e_pde5_km_act  * CONMED_PDE5I +
                    e_combo_km_act * CONMED_ERA_PDE5I)
    k34_act <- exp(lk34_act + etalk34_act)
    k43_act <- exp(lk43_act + etalk43_act)
    kmet    <- exp(lkmet + etalkmet)

    # Derived parent elimination rate (CL / V_c). The parent loses material at
    # rate kel * central; only a fraction kmet * V_c / CL of that loss
    # generates ACT-333679 in the metabolite central compartment (the
    # remainder is eliminated via routes not represented in the model).
    kel <- cl / vc

    # ----- ODE system -----------------------------------------------------------
    # Two-compartment selexipag (depot -> central -> peripheral1) feeding a
    # two-compartment ACT-333679 metabolite (central_act <-> peripheral1_act).
    # The metabolite formation arrow uses k_met (independent first-order rate)
    # so the parent's CL captures total parent clearance while the metabolite
    # appears at rate k_met * A_parent_central (Krause 2017 Figure 1).
    d/dt(depot)           <- -ka * depot
    d/dt(central)         <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1)     <-  k12 * central - k21 * peripheral1
    d/dt(central_act)     <-  kmet * central - km_act * central_act - k34_act * central_act + k43_act * peripheral1_act
    d/dt(peripheral1_act) <-  k34_act * central_act - k43_act * peripheral1_act

    # Absorption lag time on the oral depot (fixed at 0.668 h).
    alag(depot) <- tlag

    # ----- Observations ---------------------------------------------------------
    Cc     <- central     / vc
    Cc_act <- central_act / vc_act

    Cc     ~ prop(propSd)
    Cc_act ~ prop(propSd_act)
  })
}
