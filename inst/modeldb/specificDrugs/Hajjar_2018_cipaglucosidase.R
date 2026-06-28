Hajjar_2018_cipaglucosidase <- function() {
  description <- paste(
    "Two-compartment population PK model for IV cipaglucosidase alfa",
    "(ATB200; recombinant human acid alpha-glucosidase / rhGAA) in adult",
    "patients with Pompe disease (Hajjar 2018 ACCP poster, phase 1/2",
    "study ATB200-02 / NCT02675465). Disposition has parallel linear",
    "clearance (CL 0.569 L/h per 70 kg) and Michaelis-Menten saturable",
    "elimination (Vmax 98.6 mg/h per 70 kg, Km 62.4 mg/L) from the",
    "central compartment. Clearance and volume parameters are",
    "allometrically scaled by total body weight (exponent 0.75 fixed on",
    "all clearances, exponent 1 fixed on all volumes; reference 70 kg).",
    "Co-administration of the pharmacological chaperone miglustat",
    "(AT2221) reduces ATB200 linear CL: 130 mg AT2221 multiplies CL by",
    "0.738 (26.2% reduction), 260 mg AT2221 multiplies CL by 0.595 (40.5%",
    "reduction); both effects are estimated as separate categorical",
    "covariate multipliers (paper Methods 'A categorical covariate effect",
    "model was implemented'). Residual error is proportional",
    "(variance 0.0317, SD 0.178 on the linear-scale concentration)."
  )
  reference <- paste(
    "Hajjar JL, Mondick JT, Gastonguay MR, Mulberg AE, Johnson FK.",
    "Population pharmacokinetic modeling of enzyme replacement therapy",
    "ATB200 and pharmacological chaperone AT2221 in adult patients with",
    "Pompe disease and simulation to predict adolescent exposures. ACCP",
    "Annual Meeting Poster; September 2018; Seattle, WA.",
    "https://metrumrg.com/wp-content/uploads/Pubs/2018-ACCP-Population-PK-of-ATB200-AT221-in-Pompe-Patients_2018-09-18-Poster_L1e.pdf"
  )
  vignette <- "Hajjar_2018_pompe_disease"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline body weight (kg). Drives the allometric scaling of all",
        "clearances (exponent 0.75 fixed) and all volumes (exponent 1.0",
        "fixed) with reference 70 kg per Hajjar 2018 Methods 'Modeling'",
        "bullet 'Clearance and volume parameters were allometrically",
        "scaled by individual weights normalized to 70 kg body weight,",
        "with exponents fixed to values of 0.75 and 1, respectively'.",
        "Adult cohort baseline mean age 49.4 years; sex 10 M / 5 F across",
        "the 15 adults included in the modeling (Hajjar 2018 Table 1)."
      ),
      source_name        = "WT"
    ),
    DOSE_130MG = list(
      description        = "Indicator that 130 mg AT2221 (miglustat) is co-administered with the ATB200 dose",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ATB200 administered without AT2221, or with the 260 mg AT2221 dose; reference linear CL applies when both DOSE_130MG and DOSE_260MG are 0)",
      notes              = paste(
        "Per-dose-record indicator. 1 = the ATB200 dose was administered",
        "with 130 mg AT2221 co-dosing on that occasion. Enters the",
        "ATB200 linear-CL equation as a categorical multiplier:",
        "CL = CL_ref * e_dose130mg_cl ^ DOSE_130MG (Hajjar 2018",
        "Methods 'Modeling' bullet 'Fractional changes in CL were modeled",
        "as CL * Eff, where Eff is the fractional change in CL due to",
        "AT2221 dosing' and Table 2 'Fractional change in ATB200 CL for",
        "130 mg AT2221 = 0.738'). Member of the canonical DOSE_<N>MG",
        "family (siblings DOSE_50MG, DOSE_70MG, DOSE_400MG). Mutually",
        "exclusive with DOSE_260MG (a single dose occasion uses at most",
        "one AT2221 dose level)."
      ),
      source_name        = NA_character_
    ),
    DOSE_260MG = list(
      description        = "Indicator that 260 mg AT2221 (miglustat) is co-administered with the ATB200 dose",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ATB200 administered without AT2221, or with the 130 mg AT2221 dose; reference linear CL applies when both DOSE_130MG and DOSE_260MG are 0)",
      notes              = paste(
        "Per-dose-record indicator. 1 = the ATB200 dose was administered",
        "with 260 mg AT2221 co-dosing on that occasion. Enters the",
        "ATB200 linear-CL equation as a categorical multiplier:",
        "CL = CL_ref * e_dose260mg_cl ^ DOSE_260MG (Hajjar 2018",
        "Methods 'Modeling' and Table 2 'Fractional change in ATB200 CL",
        "for 260 mg AT2221 = 0.595'). Member of the canonical",
        "DOSE_<N>MG family. Mutually exclusive with DOSE_130MG."
      ),
      source_name        = NA_character_
    )
  )

  covariatesDataExcluded <- list(
    BACT = list(
      description        = "Prior enzyme-replacement-therapy experience indicator (1 = ERT-experienced with alglucosidase alfa; 0 = ERT-naive)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ERT-naive)",
      notes              = paste(
        "Screened as a categorical covariate on linear CL via the model",
        "CL * theta_ERT (Hajjar 2018 Methods 'Modeling' last two",
        "bullets: 'An exploratory covariate analysis was performed for",
        "both ATB200 and AT2221 data to investigate the effects of ERT",
        "experience on CL ... ERT = 0 for naive and ERT = 1 for",
        "experienced'). NOT retained in the final model: the mean ERT",
        "covariate effect was 1.16 with 95% CI [0, 2.86] for ATB200, and",
        "the wide CI was attributed to the small ERT-naive sample size",
        "(n = 5 vs n = 10 ERT-experienced) (Results 'Population",
        "Pharmacokinetic Models' bullets 4-6). Documented here for",
        "covariate-screen provenance; not referenced inside model()."
      ),
      source_name        = "ERT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 15L,
    n_studies      = 1L,
    age_range      = "24-66 years",
    age_mean       = "49.4 years",
    weight_range   = "not reported (paper states only that allometric scaling normalised to 70 kg)",
    sex_female_pct = 100 * 6 / 15,
    disease_state  = paste(
      "Adults with Pompe disease (genetic deficiency of acid",
      "alpha-glucosidase / GAA). 10 ERT-experienced adults previously",
      "treated with alglucosidase alfa (mean 4.8 years on prior ERT,",
      "SD 1.42) and 5 ERT-naive adults. One additional ERT-experienced",
      "subject was excluded from the modeling for missing PK values",
      "(Table 1 footnote a). Baseline 6-Minute Walk Test 392 m",
      "(SD 93) ERT-experienced vs 400 m (SD 84) ERT-naive; baseline",
      "forced vital capacity (upright) 52% predicted (SD 13) vs 53%",
      "predicted (SD 20)."
    ),
    dose_range     = paste(
      "ERT-experienced (n = 10): successive single doses of 5, 10, and",
      "20 mg/kg ATB200 alone, then 20 mg/kg ATB200 + 130 mg AT2221, then",
      "20 mg/kg ATB200 + 260 mg AT2221. ERT-naive (n = 5): single dose",
      "of 20 mg/kg ATB200 + 260 mg AT2221. Plasma samples collected over",
      "24 h periods for both ATB200 and AT2221 concentration assays."
    ),
    regions        = "United States; ATB200-02 (NCT02675465), a phase 1/2 study sponsored by Amicus Therapeutics.",
    notes          = paste(
      "Adult cohort only; the same poster reports Monte-Carlo",
      "simulations forward-projecting to adolescents (12 to <18 years",
      "old) with CDC age- and body-weight distributions, but the",
      "structural model is fit to adult data only. Parameter estimation",
      "used NONMEM 7.3 with First Order Conditional Estimation (FOCE)",
      "and Maximum Likelihood Estimation. 500-trial Monte-Carlo",
      "predictive checks (quantile-quantile plots and visual predictive",
      "checks for the 5/50/95 percentiles) supported the structural",
      "model. The 2018 ACCP poster source carries no peer-reviewed DOI;",
      "the on-disk metadata DOI 10.1038/mt.2009.53 is a vendor-index",
      "defect (it points to an unrelated 2009 Molecular Therapy paper)",
      "and is not used as the canonical reference for this extraction."
    )
  )

  ini({
    # =========================================================================
    # ATB200 (cipaglucosidase alfa, IV) structural parameters at the reference
    # subject WT = 70 kg, no AT2221 co-dosing (DOSE_130MG = 0, DOSE_260MG = 0).
    # Estimated point values from Hajjar 2018 Table 2 'ATB200 model parameter
    # estimates' column. All values entered as point estimates (no fixed()
    # wrapper) because the poster reports a non-zero RSE for each.
    # =========================================================================
    lcl   <- log(0.569);  label("Linear clearance CL (L/h, at WT 70 kg, no AT2221)")                              # Hajjar 2018 Table 2: CL = 0.569 L/h (RSE 27.2%)
    lvc   <- log(2.63);   label("Central volume of distribution Vc (L, at WT 70 kg)")                             # Hajjar 2018 Table 2: Vc = 2.63 L (RSE 19.1%)
    lq    <- log(0.151);  label("Inter-compartmental clearance Q (L/h, at WT 70 kg)")                             # Hajjar 2018 Table 2: Q = 0.151 L/h (RSE 13.3%)
    lvp   <- log(0.85);   label("Peripheral volume of distribution Vp (L, at WT 70 kg)")                          # Hajjar 2018 Table 2: Vp = 0.85 L (RSE 8.51%)

    # Michaelis-Menten saturable elimination from the central compartment. The
    # poster reports Vmax = 98.6 with units 'ug/L*h' and Km = 62.4 ug/mL. The
    # 'ug/L*h' label is dimensionally inconsistent with the standard NONMEM
    # parameterisation `dA(central)/dt = -Vmax * Cc / (Km + Cc)` (Vmax has
    # units of mass/time, Km has units of mass/volume). It is also implausible
    # under any conc/time interpretation: with Cmax ~ 530 ug/mL at a 20 mg/kg
    # IV dose, saturated elimination at 98.6 ug/L/h would require ~5400 h to
    # deplete the dose, ~1000x slower than the reported terminal half-life.
    # The only physically consistent interpretation is Vmax in mg/h (matching
    # the reported NCA AUC of ~2040 ug.h/mL at 20 mg/kg in 70 kg adults, which
    # requires the MM arm to remove ~700 mg over a few hours of saturated
    # elimination -- consistent with Vmax = 98.6 mg/h). Encoded accordingly
    # and documented in the vignette Errata.
    lvmax <- log(98.6);   label("Michaelis-Menten Vmax (mg/h, at WT 70 kg; see vignette Errata for unit reconciliation)")  # Hajjar 2018 Table 2: Vmax = 98.6 (RSE 16.5%); source unit 'ug/L*h' reinterpreted as mg/h on dimensional grounds
    lkm   <- log(62.4);   label("Michaelis-Menten Km (mg/L = ug/mL)")                                             # Hajjar 2018 Table 2: Km = 62.4 ug/mL (RSE 18.2%); 62.4 ug/mL identically 62.4 mg/L

    # =========================================================================
    # AT2221 (miglustat) co-administration effects on ATB200 linear CL.
    # Categorical multipliers applied as CL = CL_ref * e^DOSE indicator
    # (Hajjar 2018 Methods 'Fractional changes in CL were modeled as CL *
    # Eff'). Stored on the natural (multiplicative) scale -- a value < 1
    # represents reduced CL with chaperone co-dosing.
    # =========================================================================
    e_dose130mg_cl <- 0.738; label("Fractional change in ATB200 linear CL with 130 mg AT2221 co-dosing (unitless)")  # Hajjar 2018 Table 2: 0.738 (RSE 2.80%)
    e_dose260mg_cl <- 0.595; label("Fractional change in ATB200 linear CL with 260 mg AT2221 co-dosing (unitless)")  # Hajjar 2018 Table 2: 0.595 (RSE 4.94%)

    # =========================================================================
    # Allometric exponents (Hajjar 2018 Methods 'Clearance and volume
    # parameters were allometrically scaled by individual weights normalized
    # to 70 kg body weight, with exponents fixed to values of 0.75 and 1,
    # respectively'). Applied to ALL clearances (CL, Q, Vmax-mass-flow)
    # and ALL volumes (Vc, Vp) via the canonical Anderson-Holford scaling.
    # =========================================================================
    allo_cl <- fixed(0.75);  label("Allometric exponent on clearances (unitless; fixed at theoretical 0.75)")     # Hajjar 2018 Methods 'Modeling' bullet: exponents fixed to 0.75 and 1
    allo_v  <- fixed(1.00);  label("Allometric exponent on volumes (unitless; fixed at theoretical 1.00)")        # Hajjar 2018 Methods 'Modeling' bullet: exponents fixed to 0.75 and 1

    # =========================================================================
    # Between-subject variability (Hajjar 2018 Table 2 'BSV' column reported
    # as a percent CV alongside the fixed-effect estimate). Lognormal
    # variances computed as omega^2 = log(1 + CV^2):
    #   CL 28.35% -> log(1 + 0.2835^2) = 0.07730566 (RSE of BSV 93.1%)
    #   Vc 15.48% -> log(1 + 0.1548^2) = 0.02368043 (RSE of BSV 68.1%)
    # Q and Vp do not have a BSV reported in Table 2 (single-value cells), so
    # no eta is included for them. Single etas; no correlation block reported.
    # =========================================================================
    etalcl ~ 0.07730566  # Hajjar 2018 Table 2: BSV CL 28.35% (RSE 93.1%)
    etalvc ~ 0.02368043  # Hajjar 2018 Table 2: BSV Vc 15.48% (RSE 68.1%)

    # =========================================================================
    # Residual error. Table 2 reports 'Variance of residual error = 0.0317
    # (RSE 5.12%)'. With concentrations in the 10s-100s of ug/mL, a variance
    # this small is consistent with a proportional / log-additive residual
    # rather than a linear-scale additive variance. Encoded as propSd =
    # sqrt(0.0317) = 0.178 on the linear-scale concentration.
    # =========================================================================
    propSd <- 0.1780449; label("Proportional residual error (fraction CV) on linear-scale concentration")          # Hajjar 2018 Table 2: variance = 0.0317; SD = sqrt(0.0317)
  })

  model({
    # Reference body weight for allometric scaling.
    ref_wt <- 70

    # ---- Individual PK parameters ----------------------------------------
    # Linear CL: log-normal eta + allometric WT^0.75 + categorical AT2221
    # multiplier. The two AT2221 dose effects are applied as ratio^indicator
    # (the no-AT2221 reference state has both indicators = 0 so neither
    # multiplier contributes).
    cl <- exp(lcl + etalcl) * (WT / ref_wt) ^ allo_cl *
          e_dose130mg_cl ^ DOSE_130MG * e_dose260mg_cl ^ DOSE_260MG

    # Central volume: log-normal eta + allometric WT^1.
    vc <- exp(lvc + etalvc) * (WT / ref_wt) ^ allo_v

    # Inter-compartmental clearance and peripheral volume: allometric only
    # (no BSV reported by paper).
    q  <- exp(lq) * (WT / ref_wt) ^ allo_cl
    vp <- exp(lvp) * (WT / ref_wt) ^ allo_v

    # MM Vmax allometrically scales as a clearance (mass/time, same exponent
    # as CL). Km is a concentration and does not scale with body weight.
    vmax <- exp(lvmax) * (WT / ref_wt) ^ allo_cl
    km   <- exp(lkm)

    # ---- ODE system ------------------------------------------------------
    # dA(central)/dt = - CL/Vc * A_central                 (linear elimination)
    #                  - Vmax * (A/Vc) / (Km + A/Vc)       (Michaelis-Menten)
    #                  - Q/Vc * A_central + Q/Vp * A_p     (distribution)
    # The MM term carries units of mass/time directly because vmax is mg/h
    # and (central/vc)/(km + central/vc) is dimensionless.
    # NOTE on parameterisation: the MM expression is written inline in the
    # d/dt() line (NOT as an intermediate `mm_rate` LHS). This is required
    # to defeat rxode2's `rxSolve.rxUi` `useLinCmt = TRUE` heuristic: when
    # the linear-backbone params (cl/vc, q/vc, q/vp) match the 2-cpt linear
    # template AND the MM term is a free symbol from rxode2's linear-RHS
    # parser perspective (which does NOT substitute intermediate LHSs),
    # rxode2 silently converts the model to a linCmtA() analytic solver
    # that drops the MM contribution. Inlining `vmax*(central/vc)/(km +
    # central/vc)` makes the d/dt() RHS visibly nonlinear in the state and
    # disables the conversion -- so the saturable elimination is honoured
    # in both the rxUi solve path and the rxSolve.default path.
    d/dt(central)     <- -(cl / vc) * central -
                          vmax * (central / vc) / (km + central / vc) -
                          (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1

    # Concentration in the central compartment (dose in mg, vc in L => mg/L,
    # which equals ug/mL for matching the paper's plasma-concentration axes).
    Cc <- central / vc

    # Residual error: proportional on the linear-scale concentration. Cc in
    # mg/L = ug/mL matches Hajjar 2018 Figure 1 axes ('Concentrations, ug/mL').
    Cc ~ prop(propSd)
  })
}
