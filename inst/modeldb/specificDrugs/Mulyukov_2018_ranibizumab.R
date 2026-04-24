Mulyukov_2018_ranibizumab <- function() {
  description <- "Indirect-response PK/PD model of intravitreal ranibizumab on best-corrected visual acuity (BCVA, ETDRS letters) in anti-VEGF-naive adults with neovascular age-related macular degeneration (Mulyukov 2018). BCVA is driven by an indirect-response ODE in which drug concentration stimulates the BCVA production rate (kin) through a Michaelis-Menten-like term with a time-dependent maximum effect Emax(t) = Emax_ss + dEmax_0 * exp(-kEmax * t). The PK is a fixed first-order vitreous-elimination placeholder (kel = 0.077/day, vitreous volume = 4 mL, no IIV) borrowed from a previous population PK analysis (reference 20 of the paper) because vitreous PK data were not collected in the development studies."
  reference <- "Mulyukov Z, Weber S, Pigeolet E, Clemens A, Lehr T, Racine A. Neovascular Age-Related Macular Degeneration: A Visual Acuity Model of Natural Disease Progression and Ranibizumab Treatment Effect. CPT Pharmacometrics Syst Pharmacol. 2018;7(10):660-669. doi:10.1002/psp4.12322. PMID: 30043524."
  vignette <- "Mulyukov_2018_ranibizumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L (equivalent to ug/mL)", response = "BCVA (ETDRS letters, 0-100)")

  covariateData <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Power-form effect on steady-state drug effect Emax_ss normalized as (AGE/77)^beta_Emax_ss,AGE (Mulyukov 2018 Eq. 3, Table 2). Reference 77 years is the study-population mean baseline age. Paper narrative: a 4 ETDRS letter reduction in 12-month BCVA improvement is expected for an 85-year-old vs a 65-year-old patient.",
      source_name        = "AGE"
    ),
    BCVA = list(
      description        = "Observed baseline best-corrected visual acuity",
      units              = "ETDRS letters (0-100)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Used as the per-subject center for the baseline BCVA initial condition g0_i = BVA_i + eta_{1,i} (Mulyukov 2018 Methods, paragraph above Eq. 3: 'We modeled baseline VA as normally distributed around observed value BVA'). Canonical register alias: source column 'BVA' in Mulyukov 2018. Study-population baseline mean (SD) = 54 (13) letters across the four development datasets.",
      source_name        = "BVA"
    )
  )

  population <- list(
    n_subjects          = 1524L,
    n_observations      = 29754L,
    n_studies           = 4L,
    age_range           = "50-96 years; mean (SD) 77 (7.5) years",
    age_median          = "77 years (mean)",
    weight_range        = NULL,
    sex_female_pct      = 60,
    race_ethnicity      = NULL,
    disease_state       = "Anti-VEGF treatment-naive neovascular (wet) age-related macular degeneration with best-corrected visual acuity 25-70 ETDRS letters at entry.",
    dose_range          = "Ranibizumab 0.3 mg or 0.5 mg intravitreal injection, monthly (q4w) or quarterly (q12w) after three monthly loading doses; sham comparator arms also modelled simultaneously.",
    regions             = "Multinational (ANCHOR, MARINA, PIER, and EXCITE phase III programmes).",
    baseline_BCVA       = "Mean (SD) 54 (13) ETDRS letters; range 3-84 letters.",
    external_validation = "HARBOR study (ranibizumab 0.5 mg and 2.0 mg q4w, summary data only) used for predictive check only; not part of model fitting.",
    notes               = "Pooled from the ANCHOR (2 y), MARINA (2 y), PIER (1 y), and EXCITE (1 y) phase III ranibizumab trials (Mulyukov 2018 Table 1). The PDT (verteporfin) arm of ANCHOR and the year-2 PIER data were excluded from the modelled dataset per the paper's Methods. Approximately 40% of patients were male. Weight / race-ethnicity distributions were not reported in Table 1."
  )

  ini({
    # ------------------------------------------------------------------------
    # PK backbone (fixed, no IIV) borrowed from a previous population PK
    # analysis cited as reference 20 in Mulyukov 2018. Paper (Methods - Model
    # development, paragraph 1): "for all modeling we used results of a previous
    # population PK analysis that concluded the vitreous concentration of
    # ranibizumab follows first-order elimination PKs, with a half-life of 9
    # days ... Due to the absence of vitreous PK data, intersubject variability
    # of PK parameters was disregarded." Vitreous volume set to 4 mL (Methods
    # paragraph above Eq. 1).
    # ------------------------------------------------------------------------
    lkel  <- fixed(log(0.077));   label("Vitreous elimination rate constant kel (1/day)")               # Mulyukov 2018 Table 2: k_elim 0.077/day (t1/2 = 9 days), fixed
    lvvit <- fixed(log(0.004));   label("Vitreous volume of distribution (L)")                          # Mulyukov 2018 Methods (Model development): vitreous volume 4 mL = 0.004 L, fixed

    # ------------------------------------------------------------------------
    # BCVA indirect-response PD model (Mulyukov 2018 Eq. 1, Table 2).
    # Natural-progression steady state g_ss = kin/kout (paper Methods, paragraph
    # after Eq. 1). With no treatment, the BCVA state decays from the observed
    # baseline g0 toward g_ss at rate kout.
    # ------------------------------------------------------------------------
    lgss       <- log(11);          label("Equilibrium BCVA under natural disease progression g_ss (ETDRS letters)")  # Mulyukov 2018 Table 2: g_ss 11 letters (RSE 5%)
    lkout      <- log(0.19 / 365.25); label("BCVA deterioration rate constant kout (1/day)")                          # Mulyukov 2018 Table 2: k_out 0.19/year (t1/2 = 3.6 years); converted to 1/day by /365.25
    lEmaxss    <- log(6.1);         label("Steady-state maximum drug effect Emax_ss on BCVA kin (unitless)")          # Mulyukov 2018 Table 2: Emax_ss 6.1 (RSE 7%)
    ldEmax0    <- log(41);          label("Additional drug effect at treatment onset dEmax_0 (unitless)")             # Mulyukov 2018 Table 2: dEmax_0 41 (RSE 12%)
    lkEmax     <- fixed(log(0.046));label("Rate of Emax decay kEmax (1/day)")                                         # Mulyukov 2018 Table 2: k_Emax 0.046/day (t1/2 = 15 days), fixed
    lEC50      <- log(2.1);         label("Drug concentration at half-maximum effect EC50 (mg/L = ug/mL)")            # Mulyukov 2018 Table 2: EC50 2.1 ug/mL (= mg/L) (RSE 35%)

    # ------------------------------------------------------------------------
    # Covariate effect (Mulyukov 2018 Table 2 / Eq. 3).
    # Age / baseline-VA / sex effects on k_out, Emax_ss, and dEmax_0 were tested.
    # Only age on Emax_ss was significant and retained in the final model.
    # ------------------------------------------------------------------------
    e_age_Emaxss <- -1.4;           label("Power exponent of AGE/77 on Emax_ss (unitless)")                           # Mulyukov 2018 Table 2: beta_Emax_ss,AGE -1.4 (RSE 18%)

    # ------------------------------------------------------------------------
    # Typical-value placeholder for the additive random effect on baseline BCVA.
    # Paper's parameterization g_{i,0} = BVA_i + eta_{1,i} has no typical-value
    # theta for the baseline (BVA is subject-level data, not a population
    # parameter). g0res is a fixed(0) theta kept only to pair naming-wise with
    # etag0res per the nlmixr2lib convention (eta<x> must pair with theta <x>).
    # ------------------------------------------------------------------------
    g0res <- fixed(0); label("Typical additive residual on baseline BCVA around the observed value (ETDRS letters); paper Eq. g_{i,0} = BVA_i + eta_{1,i}, so typical value is 0")  # Mulyukov 2018 Methods (paragraph above Eq. 3): baseline modelled as normally distributed around observed value with additive random effect

    # ------------------------------------------------------------------------
    # Inter-individual variability.
    # Baseline BCVA uses an ADDITIVE per-subject random effect in letters (paper
    # Eq.: g_{i,0} = BVA_i + eta_{1,i}; Table 2 reports the SD, not a CV).
    # The other three PD parameters (k_out, Emax_ss, dEmax_0) use log-normal
    # IIV with large CV%; convert to the log-scale variance via
    #   omega^2 = log(CV^2 + 1)
    #     k_out   CV 730%   -> log(7.30^2 + 1)  = 3.994
    #     Emax_ss CV 110%   -> log(1.10^2 + 1)  = 0.7929
    #     dEmax_0 CV 1100%  -> log(11.0^2 + 1)  = 4.804
    # ------------------------------------------------------------------------
    etag0res   ~ 16.81   # Mulyukov 2018 Table 2: IIV g0 SD = 4.1 letters (additive); omega^2 = 4.1^2 = 16.81
    etalkout   ~ 3.994   # Mulyukov 2018 Table 2: IIV k_out CV 730%; log(7.30^2 + 1) = 3.994
    etalEmaxss ~ 0.7929  # Mulyukov 2018 Table 2: IIV Emax_ss CV 110%; log(1.10^2 + 1) = 0.7929
    etaldEmax0 ~ 4.804   # Mulyukov 2018 Table 2: IIV dEmax_0 CV 1100%; log(11.0^2 + 1) = 4.804

    # ------------------------------------------------------------------------
    # Residual error. Paper Table 2 reports two additive residual SDs:
    #   sigma_treatment = 5 letters (treated arms)
    #   sigma_sham      = 7 letters (untreated arms)
    # The library exposes the treated-arm value (primary use case for
    # ranibizumab simulation); the sham-arm value is documented in the
    # vignette's Assumptions and deviations.
    # ------------------------------------------------------------------------
    bcvaaddSd  <- 5; label("Additive residual error on BCVA for treated subjects (ETDRS letters)")  # Mulyukov 2018 Table 2: sigma_treatment 5 letters (RSE 0.5%)
  })

  model({
    # ------------------------------------------------------------------
    # PK parameters (fixed, no IIV).
    # ------------------------------------------------------------------
    kel  <- exp(lkel)
    vvit <- exp(lvvit)

    # ------------------------------------------------------------------
    # Individual PD parameters (Mulyukov 2018 Eq. 3 / Table 2).
    # Only age-on-Emax_ss survived covariate selection; other covariate
    # effects (baseline VA, sex, age on k_out / dEmax_0) were dropped.
    # ------------------------------------------------------------------
    gss    <- exp(lgss)
    kout   <- exp(lkout + etalkout)
    Emaxss <- exp(lEmaxss + etalEmaxss) * (AGE / 77)^e_age_Emaxss
    dEmax0 <- exp(ldEmax0 + etaldEmax0)
    kEmax  <- exp(lkEmax)
    EC50   <- exp(lEC50)

    # Baseline BCVA: additive IIV around the observed per-subject value.
    # Paper: g_{i,0} = BVA_i + eta_{1,i}; typical-value theta g0res fixed(0).
    g0  <- BCVA + g0res + etag0res

    # Natural-progression mass balance: kin chosen so untreated steady state
    # equals gss (paper: g_ss = kin / kout).
    kin <- kout * gss

    # Time-dependent maximum drug effect (Mulyukov 2018 Eq. 2). `time` is the
    # integration time in days, equal to time from first dose for the event
    # datasets this model expects.
    Emax_t <- Emaxss + dEmax0 * exp(-kEmax * time)

    # ------------------------------------------------------------------
    # PK ODE: first-order vitreous elimination. The intravitreal dose is
    # administered directly into the central compartment (F = 1 by default);
    # there is no absorption compartment. Concentration units: central (mg) /
    # vvit (L) = mg/L, which equals ug/mL in the paper's tables.
    # ------------------------------------------------------------------
    Cc <- central / vvit

    d/dt(central) <- -kel * central

    # ------------------------------------------------------------------
    # PD ODE: indirect response on BCVA with stimulation of kin.
    # Paper (Results - Model development): "Stimulation of kin rather than
    # suppression of kout was selected to describe the drug effect as it led
    # to better stability of model fit with no other significant differences."
    # ------------------------------------------------------------------
    bcva(0)    <- g0
    d/dt(bcva) <- kin * (1 + Emax_t * Cc / (EC50 + Cc)) - kout * bcva

    # ------------------------------------------------------------------
    # Observation and error.
    # ------------------------------------------------------------------
    bcva ~ add(bcvaaddSd)
  })
}
