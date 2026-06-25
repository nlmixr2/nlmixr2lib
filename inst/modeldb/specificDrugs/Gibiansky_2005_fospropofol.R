Gibiansky_2005_fospropofol <- function() {
  description <- "Joint two-compartment fospropofol (GPI 15715, AQUAVAN) prodrug + intermediate delay compartment + two-compartment propofol active-metabolite population PK model in adults receiving IV bolus AQUAVAN for procedural sedation (Gibiansky 2005, ASCPT poster, colonoscopy sedation Phase II study). The model assumes complete metabolism of GPI 15715 to propofol via systemic alkaline-phosphatase hydrolysis; the intermediate compartment captures the appearance delay between GPI 15715 elimination from plasma and the corresponding rise in propofol concentration. Lean body mass (LBM, reference 55 kg) was retained as a linear-fractional covariate on GPI 15715 central volume Vc_GPI, GPI 15715 metabolic clearance CL_GPI, and propofol central volume Vc_PR; fentanyl premedication exposure, age, sex, and other demographics/laboratory covariates were tested but not retained. Propofol Vc_PR (6.91 L) was fixed as the data were insufficient for joint estimation with CL_PR (the model is identifiable on CL_PR = K10_PR * Vc_PR = 4.53 L/min)."
  reference   <- paste(
    "Gibiansky E, Gibiansky L, Enriquez J.",
    "Population pharmacokinetic model of sedative doses of GPI 15715 and propofol",
    "liberated from GPI 15715. Clin Pharmacol Ther. 2005;77(2):P32 (PII-87, ASCPT",
    "Annual Meeting poster). doi:10.1016/j.clpt.2004.12.076.",
    "Poster PDF hosted at https://metrumrg.com/wp-content/uploads/2018/08/ascpt_2005_ppkmodelgpi.pdf",
    sep = " "
  )
  vignette    <- "Gibiansky_2005_fospropofol"
  paper_specific_compartments <- "delay"

  units       <- list(time = "min", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    LBM = list(
      description        = "Lean body mass (lean body weight) at study entry.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-fractional effect on GPI 15715 central volume Vc_GPI, GPI 15715 metabolic clearance CL_GPI, and propofol central volume Vc_PR, centered at LBM = 55 kg per the poster Results section ('for each kilogram of LBW (from 55 kg)'). LBM range in the study population was 37-81 kg.",
      source_name        = "LBW"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight at study entry (range 45-140 kg).",
      units       = "kg",
      type        = "continuous",
      notes       = "Tested in covariate screening; not retained in the final model because LBM was the better predictor and WT, LBM, and BMI were highly co-linear (poster Modeling Steps step 5)."
    ),
    BMI = list(
      description = "Body mass index at study entry.",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Tested in covariate screening; not retained (co-linear with WT and LBM)."
    ),
    AGE = list(
      description = "Age at study entry (range 20-85 years; 18 patients > 65 y).",
      units       = "years",
      type        = "continuous",
      notes       = "Tested; no association with GPI 15715 or propofol PK detected (poster Covariate Effects section: 'Older age (>65 y) was not associated with changes in the PK of either GPI 15715 or propofol')."
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested; effect on PK was fully explained by LBM (poster: 'gender was strongly correlated with WT and LBW. After accounting for LBW, no additional dependencies on gender were evident')."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 158L,
    n_studies      = 1L,
    age_range      = "20-85 years; 18 patients (about 11%) older than 65 y",
    weight_range   = "WT 45-140 kg; LBM 37-81 kg",
    sex_female_pct = 100 * 89 / 158,
    race_ethnicity = "Tested in covariate screening; not retained in the final model. Per-stratum percentages not reported in the poster.",
    disease_state  = "Adults undergoing elective colonoscopy receiving procedural sedation with AQUAVAN (GPI 15715, water-soluble prodrug of propofol) following premedication with IV fentanyl citrate.",
    dose_range     = "AQUAVAN initial IV bolus 7.5 to 12.5 mg/kg, with up to 4 supplemental doses (each ~25% of the initial bolus, occasional 50%) at intervals > 3 min; total cumulative dose 495 to 1,680 mg. Fentanyl citrate premedication 0.5 to 1.5 ug/kg given 5 minutes before the AQUAVAN bolus (cumulative fentanyl 11-201 ug).",
    regions        = "Phase II colonoscopy sedation study; specific countries / centres not enumerated in the poster.",
    n_observations = "597 GPI 15715 and 599 propofol plasma concentrations; 282 AQUAVAN doses administered.",
    notes          = "Randomized, open-label, dose-ranging, adaptive-dose Phase II trial. Source: ASCPT 2005 Annual Meeting poster (PII-87) authored by Gibiansky, Gibiansky, and Enriquez of Guilford Pharmaceuticals and Metrum Research Group. The poster is the only published account of the model; no peer-reviewed journal article expands on the results. Detailed demographic breakdowns (race/ethnicity stratification, weight medians, regional distribution) are not in the poster."
  )

  ini({
    # ---------------- GPI 15715 (parent prodrug) structural parameters --------
    # Final estimates from Gibiansky 2005 Table 1 (Fixed-Effect Parameters of
    # the Final PK Model). All rate constants are in 1/min; volumes in L;
    # clearances in L/min. The reference covariate value for LBM is 55 kg
    # (poster Results section). The K12 / K21 rates are carried as primary
    # ini() parameters rather than re-parameterised into Q / V_p so the
    # paper's IIV variance/covariance structure (Table 2) is preserved
    # without cross-parameter covariance translation, following the same
    # convention as Krause_2017_selexipag.R.
    lvc      <- log(6.08);    label("GPI 15715 central volume Vc_GPI at LBM = 55 kg (L)")                                  # Table 1: Vc_GPI = 6.08 L, RSE 4.6%
    lcl      <- log(0.298);   label("GPI 15715 metabolic clearance CL_GPI to propofol at LBM = 55 kg (L/min)")              # Table 1: CL_GPI = 0.298 L/min, RSE 8.0%
    lk12     <- log(0.0198);  label("GPI 15715 central-to-peripheral rate constant K12_GPI (1/min)")                        # Table 1: K12_GPI = 0.0198 1/min, RSE 24%
    lk21     <- log(0.00617); label("GPI 15715 peripheral-to-central rate constant K21_GPI (1/min)")                        # Table 1: K21_GPI = 0.00617 1/min, RSE 21%

    # ---------------- Delay -> propofol appearance rate -----------------------
    # K_GPI-PR is the first-order rate at which propofol appears in central
    # circulation from the intermediate delay state that follows GPI 15715
    # metabolism. Mechanistically analogous to a first-order absorption
    # step, hence the canonical 'ka' parameter prefix with the propofol
    # metabolite suffix.
    lka_ppf  <- log(0.982);   label("Delay-to-propofol appearance rate K_GPI-PR (1/min)")                                   # Table 1: K_GPI-PR = 0.982 1/min, RSE 14%

    # ---------------- Propofol (active metabolite) structural parameters ------
    # Vc_PR was held fixed at 6.91 L because the data were insufficient to
    # estimate it jointly with CL_PR (poster Table 1 footnote). The model is
    # identifiable on CL_PR = K10_PR * Vc_PR = 4.53 L/min; the K10_PR
    # estimate moves to maintain CL_PR as Vc_PR is varied. The IIV on Vc_PR
    # was still estimated (Table 2: omega^2 = 0.0699, RSE 29%).
    lvc_ppf  <- fixed(log(6.91)); label("Propofol central volume Vc_PR (L; fixed, see poster Table 1 footnote)")            # Table 1: Vc_PR = 6.91 L FIXED
    lk10_ppf <- log(0.655);   label("Propofol elimination rate K10_PR (1/min)")                                              # Table 1: K10_PR = 0.655 1/min, RSE 8.0%
    lk12_ppf <- log(0.732);   label("Propofol central-to-peripheral rate constant K12_PR (1/min)")                          # Table 1: K12_PR = 0.732 1/min, RSE 11%
    lk21_ppf <- log(0.0383);  label("Propofol peripheral-to-central rate constant K21_PR (1/min)")                          # Table 1: K21_PR = 0.0383 1/min, RSE 20%

    # ---------------- Covariate effects (LBM, linear fractional from 55 kg) ---
    # Poster Results: 'GPI 15715 volume (Vc_GPI) and clearance (CL_GPI), and
    # propofol volume (Vc_PR) increased with increasing LBW, by 1.8%, 2.5%,
    # and 1.4%, respectively, for each kilogram of LBW (from 55 kg).' The
    # per-kg slope is the unambiguous quantity in the paper; the printed
    # Table 1 values 0.194, 0.270, 0.155 reconcile with this prose under
    # the linear-fractional parameterisation TVP = TV * (1 + slope * (LBM -
    # 55)) only when read as slopes per 10 kg of LBM deviation (equivalent
    # to per-kg slopes 0.0194, 0.0270, 0.0155). Both readings are
    # mathematically equivalent; this model carries the per-kg form (the
    # form the abstract and Results section state explicitly) and
    # documents the interpretation in the vignette Errata.
    e_lbm_vc      <- 0.0194;  label("LBM slope on GPI 15715 Vc_GPI (per kg from 55)")    # Table 1 / Results: Vc_GPI LBW = 1.8%/kg (Table 1 prints 0.194, read as per-10-kg)
    e_lbm_cl      <- 0.0270;  label("LBM slope on GPI 15715 CL_GPI (per kg from 55)")    # Table 1 / Results: CL_GPI LBW = 2.5%/kg (Table 1 prints 0.270, read as per-10-kg)
    e_lbm_vc_ppf  <- 0.0155;  label("LBM slope on propofol Vc_PR (per kg from 55)")      # Table 1 / Results: Vc_PR LBW = 1.4%/kg (Table 1 prints 0.155, read as per-10-kg)

    # ---------------- IIV (Table 2; correlated 4x4 block) ---------------------
    # Variance / covariance estimates on the log scale of the four
    # primary parameters carrying inter-individual variability:
    # Vc_GPI, K_GPI-PR (== lka_ppf here), K12_GPI, Vc_PR. The off-diagonal
    # zeros in row 1 / column 3, row 1 / column 4, and row 2 / column 4
    # are the structural zeros of the banded matrix the paper fit.
    #
    # Row 1 (Vc_GPI):
    #   omega^2(Vc_GPI)          = 0.0727 (Table 2; RSE 29%, CV 27.5%)
    # Row 2 (K_GPI-PR):
    #   omega(Vc_GPI, K_GPI-PR)  = -0.239 (Table 2; RSE 15%, R = -0.886)
    #   omega^2(K_GPI-PR)        = 1      (Table 2; RSE 23%, CV 131%; near upper boundary)
    # Row 3 (K12_GPI):
    #   omega(Vc_GPI, K12_GPI)   = 0      (Table 2; structural zero in the banded matrix)
    #   omega(K_GPI-PR, K12_GPI) = 0.271  (Table 2; RSE 46%, R = 0.266)
    #   omega^2(K12_GPI)         = 1.04   (Table 2; RSE 49%, CV 135%; near upper boundary)
    # Row 4 (Vc_PR):
    #   omega(Vc_GPI, Vc_PR)     = 0      (Table 2; structural zero in the banded matrix)
    #   omega(K_GPI-PR, Vc_PR)   = 0      (Table 2; structural zero in the banded matrix)
    #   omega(K12_GPI, Vc_PR)    = -0.136 (Table 2; RSE 48%, R = -0.504)
    #   omega^2(Vc_PR)           = 0.0699 (Table 2; RSE 29%, CV 26.9%; estimated although the
    #                                      typical value of Vc_PR was fixed -- see lvc_ppf above)
    etalvc + etalka_ppf + etalk12 + etalvc_ppf ~ c(
      0.0727,
      -0.239, 1,
      0,      0.271,  1.04,
      0,      0,     -0.136, 0.0699
    )

    # ---------------- Residual error (combined add+prop, per output) ----------
    # Poster Error Models section: Ln(Y) = Ln(F) + W * eps with W^2 =
    # theta1 / F^2 + theta2, which expands in linear space to Var(Y - F)
    # = theta1 + theta2 * F^2 (standard add+prop). Table 2 reports the
    # variance parameters theta1 (..._ADD) and theta2 (..._PROP); the SDs
    # carried below are the square roots and match the SD column of
    # Table 2. Per-output residual-error parameters follow the
    # nlmixr2lib multi-output convention (parent suffix-free, metabolite
    # suffixed by _ppf).
    propSd     <- 0.4;    label("GPI 15715 proportional residual SD (fraction)")          # Table 2: sigma^2_GPI_PROP = 0.16, RSE 20%, CV 40.0% -> SD = sqrt(0.16) = 0.4
    addSd      <- 2.74;   label("GPI 15715 additive residual SD (ug/mL)")                 # Table 2: sigma^2_GPI_ADD = 7.52, RSE 44%, SD = 2.74 ug/mL
    propSd_ppf <- 0.378;  label("Propofol proportional residual SD (fraction)")           # Table 2: sigma^2_PR_PROP = 0.143, RSE 12%, CV 37.8% -> SD = sqrt(0.143) = 0.378
    addSd_ppf  <- 0.104;  label("Propofol additive residual SD (ug/mL)")                  # Table 2: sigma^2_PR_ADD = 0.0109, RSE 33%, SD = 0.104 ug/mL
  })

  model({
    # Reference LBM for the linear-fractional covariate term (poster
    # Results: 'from 55 kg').
    ref_lbm <- 55

    # ----- Individual structural parameters -----------------------------------
    # GPI 15715 (parent) carries IIV on Vc only (per Table 2 OMEGA block);
    # CL_GPI, K21_GPI have no IIV component in the final model.
    vc     <- exp(lvc + etalvc)       * (1 + e_lbm_vc     * (LBM - ref_lbm))
    cl     <- exp(lcl)                * (1 + e_lbm_cl     * (LBM - ref_lbm))
    k12    <- exp(lk12 + etalk12)
    k21    <- exp(lk21)
    kel    <- cl / vc

    # Delay -> propofol appearance rate (paper-faithful K_GPI-PR with IIV).
    ka_ppf <- exp(lka_ppf + etalka_ppf)

    # Propofol (active metabolite) -- Vc_PR fixed typical value with
    # estimated IIV, K10_PR / K12_PR / K21_PR estimated without IIV.
    vc_ppf  <- exp(lvc_ppf + etalvc_ppf) * (1 + e_lbm_vc_ppf * (LBM - ref_lbm))
    k10_ppf <- exp(lk10_ppf)
    k12_ppf <- exp(lk12_ppf)
    k21_ppf <- exp(lk21_ppf)

    # ----- ODE system ---------------------------------------------------------
    # Five-compartment structure (poster Figure 2):
    #   AQUAVAN bolus -> central (GPI 15715) <-> peripheral1 (GPI 15715)
    #   central (GPI 15715) -- kel --> delay (intermediate)
    #   delay -- ka_ppf --> central_ppf (propofol) <-> peripheral1_ppf
    #   central_ppf -- k10_ppf --> eliminated
    # The model assumes complete (100%) metabolism of GPI 15715 to propofol
    # via systemic alkaline-phosphatase hydrolysis; the delay state has no
    # explicit volume because the published parameterisation collapses any
    # mass-conversion factor into the apparent Vc_PR.
    d/dt(central)         <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1)     <-  k12 * central - k21 * peripheral1
    d/dt(delay)           <-  kel * central - ka_ppf * delay
    d/dt(central_ppf)     <-  ka_ppf * delay - k10_ppf * central_ppf - k12_ppf * central_ppf + k21_ppf * peripheral1_ppf
    d/dt(peripheral1_ppf) <-  k12_ppf * central_ppf - k21_ppf * peripheral1_ppf

    # ----- Observations and residual error ------------------------------------
    Cc     <- central     / vc
    Cc_ppf <- central_ppf / vc_ppf

    Cc     ~ add(addSd)     + prop(propSd)
    Cc_ppf ~ add(addSd_ppf) + prop(propSd_ppf)
  })
}
