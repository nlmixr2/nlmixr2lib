Mulder_2025_ribavirin <- function() {
  description <- paste(
    "Integrated population PK/PD model for oral ribavirin (RBV) in solid",
    "organ transplant (SOT) recipients with chronic hepatitis E virus",
    "(HEV) infection (Mulder 2025). PK: two-compartment model with first-",
    "order absorption; ka, Vc, Q, and Vp fixed at the Wu 2015 (HCV-",
    "patient) starting-model estimates and CL re-estimated. Covariates",
    "carried over from Wu 2015: allometric weight on Vc (exponent 1.29)",
    "and Vp (exponent 0.725) with reference weight 79 kg, and a female-",
    "factor 0.732 on Vp. New covariate estimated on CL: MDRD eGFR with",
    "a capped power effect (exponent 1.32) above an estimated threshold",
    "of 57 mL/min/1.73 m^2. Haemoglobin: indirect-response (kin/kout) on",
    "an endogenous Hb pool, with RBV producing a linear concentration-",
    "proportional acceleration of the haemoglobin loss rate (1 + slope *",
    "Cc) so that haemoglobin declines with increasing RBV exposure; kin",
    "is set per subject from the baseline-Hb covariate HGB_BL so the Hb",
    "state is at steady state pre-treatment. Viral load: target-cell-",
    "limited (Baccam / Dahari) model with three paper-specific",
    "compartments (healthy hepatocytes, infected hepatocytes, virions);",
    "RBV inhibits viral replication via an Imax/IC50 sigmoidal Emax",
    "form (Imax = 0.999, IC50 = 1000 ng/L) so production of virions from",
    "infected cells is essentially fully suppressed throughout the",
    "observed RBV concentration range; healthy hepatocyte half-life and",
    "the infected:healthy hepatocyte decay-rate ratio are fixed to the",
    "literature values (Dahari 2007) and the viral elimination rate is",
    "estimated. Initial conditions for the viral compartments are",
    "computed at baseline steady state: H(0) = 1 (arbitrary unit), I(0)",
    "= rho = 0.001 (fixed), V(0) = HEV_VLOAD per subject; the synthesis",
    "rates ksyn (healthy), beta (infection), and p (virion production)",
    "are derived inside model() from these initial conditions and the",
    "rate constants so that all three viral compartments start at",
    "baseline steady state. Hb and the three viral compartments are",
    "declared paper_specific_compartments."
  )
  reference <- paste(
    "Mulder MB, van Noort M, de Man RA, Kamar N, de Bruijne J, Knoester",
    "M, Blokzijl H, Vanwolleghem T, Roosens L, Izopet J, Gandia P, van",
    "der Eijk AA, Metselaar HJ, Ahsman MJ, van Steeg TJ, Hesselink DA,",
    "de Winter BCM. (2025). Development of a ribavirin dosing regimen",
    "in transplant recipients with chronic hepatitis E virus infection:",
    "a population pharmacokinetic and -dynamic model.",
    "J Antimicrob Chemother 80(8):2158-2168.",
    "doi:10.1093/jac/dkaf183. PMID 40581807.",
    "PK structural parameters Ka, Vc, Q, and Vp and the WT and SEXF",
    "effects on Vc and Vp are fixed from Wu LS, Rower JE, Burton JR Jr",
    "et al. (2015) Antimicrob Agents Chemother 59:2179-2188",
    "(doi:10.1128/AAC.04618-14; the reference weight 79 kg is taken",
    "from Wu 2015 Table 2 / covariate-equation form).",
    "Infected and healthy hepatocyte half-lives derive from Dahari H,",
    "Lo A, Ribeiro RM et al. (2007) J Theor Biol 247:371-381",
    "(doi:10.1016/j.jtbi.2007.03.006).",
    sep = " "
  )
  vignette <- "Mulder_2025_ribavirin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  paper_specific_compartments <- c("hb_state", "healthy", "infected", "virus")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used for allometric scaling on Vc (exponent 1.29 fixed) and Vp",
        "(exponent 0.725 fixed). Reference weight 79 kg taken from Wu",
        "2015 (upstream PK source); Mulder 2025 reports cohort median",
        "74 kg (Table 1), and the per-sex medians used in simulations",
        "are 70 kg (male) and 75 kg (female). All structural body-weight",
        "effects are carried over from Wu 2015 verbatim."
      ),
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Used as multiplicative factor 0.732^SEXF on Vp (so Vp is 26.8%",
        "smaller in females vs males). Fixed from Wu 2015 (the source",
        "paper reports the same effect as 0.732^SEX with SEX = 1",
        "indicating female)."
      ),
      source_name        = "SEX"
    ),
    CRCL = list(
      description        = "MDRD-estimated glomerular filtration rate (BSA-normalized)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Cohort median 50 mL/min/1.73 m^2 (Mulder 2025",
        "Table 1); simulations in the source paper use 57 as the",
        "median surrogate. Enters CL via a power effect with exponent",
        "1.32 (estimated, RSE 14%); the effect is centred on (and",
        "capped at) an estimated threshold of 57 mL/min/1.73 m^2 above",
        "which CL is held at its typical value. Renamed from source",
        "column eGFR to the canonical CRCL per covariate-columns.md",
        "(MDRD eGFR is a registered CRCL alias)."
      ),
      source_name        = "eGFR"
    ),
    HGB_BL = list(
      description        = "Per-subject baseline (pre-treatment) haemoglobin concentration",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used as initial condition for hb_state and as the baseline",
        "anchor for kin = kout * HGB_BL so the Hb state is at steady",
        "state pre-RBV. Cohort median 8.3 mmol/L (Mulder 2025 Table 1).",
        "Per the source paper, when an individual baseline was not",
        "available (n = 1) the cohort typical (median) was used; encode",
        "the same fallback in the data-preparation step. Renamed from",
        "source column HBBASE to the canonical HGB_BL per covariate-",
        "columns.md."
      ),
      source_name        = "HBBASE"
    ),
    HEV_VLOAD = list(
      description        = "Per-subject baseline (pre-treatment) HEV viral load",
      units              = "IU/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used as the initial condition virus(0) = HEV_VLOAD and enters",
        "the derivation of the infection-rate constant beta = kloss *",
        "rho / HEV_VLOAD and the virion production rate p = elim *",
        "HEV_VLOAD / rho so all viral compartments start at baseline",
        "steady state. Cohort median 1.886e6 IU/mL (Mulder 2025 Table",
        "1; range 527-1.68e8). Renamed from source column VLBASE to",
        "the canonical HEV_VLOAD per covariate-columns.md."
      ),
      source_name        = "VLBASE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 107L,
    n_studies      = 1L,
    age_range      = "22-84 years",
    age_median     = "56.9 years",
    weight_range   = "43.5-140 kg",
    weight_median  = "74 kg (cohort); 70 kg (male median) / 75 kg (female median) used in simulations",
    sex_female_pct = 32.7,
    disease_state  = paste(
      "Adult solid organ transplant (SOT) recipients with chronic",
      "hepatitis E virus (HEV) infection treated with oral ribavirin",
      "between September 2009 and November 2019. Organ-transplant",
      "distribution: kidney 43.9%, liver 17.8%, heart 14.9%, lung 14%,",
      "kidney+pancreas 3.7%, kidney+heart 2.8%, pancreas 0.9%,",
      "lung+liver 0.9%, lung+heart 0.9%. Most subjects on tacrolimus-",
      "based immunosuppression (84.9% tacrolimus, 71.7% glucocorticoids,",
      "58.5% mycophenolic acid, 13.2% everolimus, 6.6% sirolimus)."
    ),
    dose_range     = paste(
      "Oral ribavirin median 600 mg/day (range 100-2400) for median 90",
      "days (range 21-1333). No loading dose used in the modelled",
      "regimens; the source paper's simulation arm evaluating a 33",
      "mg/kg loading dose (Burrows 2015 protocol) is illustrative",
      "only and not part of the fitted dataset."
    ),
    regions        = "Multicentre: Netherlands, France, Belgium (5 hospitals)",
    notes          = paste(
      "Retrospective study (MEC-2018-1326). 92 subjects overlap with",
      "Mulder 2021 J Viral Hepat 28:431-435 (PMID 33135238) but no",
      "popPK/PD analysis was performed in that prior paper. 305 RBV",
      "plasma levels (range 0.1-6.2 mg/L), 443 haemoglobin",
      "observations, and 592 viral loads (38% BLQ, included via the",
      "Beal M3 likelihood method). 68.2% (73/107) achieved sustained",
      "virologic response (SVR); 28% (30/107) did not; 3.7% (4/107)",
      "unknown. Bootstrap analysis on 1000 datasets (Mulder 2025 Table",
      "2 columns 'Bootstrap of the final model')."
    )
  )

  ini({
    # ================================================================
    # Population PK model (Mulder 2025 Table 2 - "Population
    # pharmacokinetic model" block). Time unit is hours throughout.
    # ka, V2 (= Vc), Q1 (= Q), V3 (= Vp), and all body-weight and sex
    # covariate effects on Vc and Vp are fixed at the Wu 2015 starting-
    # model values; clearance was re-estimated because of deviations
    # observed in the visual predictive check (Mulder 2025 Methods,
    # "Pharmacokinetic, haemoglobin and viral load analysis").
    # ================================================================
    lka  <- fixed(log(2.91)); label("Absorption rate constant ka (1/h)")        # Table 2 row 'K a (h-1)': 2.91 (fixed) - Wu 2015
    lcl  <- log(26.4);         label("Typical apparent clearance CL/F (L/h)")    # Table 2 row 'CL (L/h)': 26.4 (RSE 15%, bootstrap median 24.3, 95% CI 15.8-36.1)
    lvc  <- fixed(log(769));   label("Central compartment volume Vc/F (L)")     # Table 2 row 'V2 (L)': 769 (fixed) - Wu 2015
    lq   <- fixed(log(104));   label("Intercompartmental clearance Q/F (L/h)")  # Table 2 row 'Q1 (L/h)': 104 (fixed) - Wu 2015
    lvp  <- fixed(log(3570));  label("Peripheral compartment volume Vp/F (L)")  # Table 2 row 'V3 (L)': 3570 (fixed) - Wu 2015

    # eGFR (CRCL) effect on CL: capped power model centred on the
    # estimated cut-off value. CL = TVCL * (min(CRCL, 57) / 57)^1.32.
    # Below 57 CL scales sub-typically; at or above 57 CL stays at
    # TVCL. The cut-off 57 was estimated (RSE 0.1%, bootstrap median
    # 57.4, 95% CI 52-180) but is encoded as fixed here because the
    # operational role is structural (the centring/cap reference) and
    # the very low RSE makes the value effectively deterministic; the
    # discrepancy with the wide bootstrap CI is documented in the
    # validation vignette Errata.
    e_crcl_cl  <- 1.32;        label("Power exponent of capped CRCL on CL (unitless)")  # Table 2 row 'eGFR on Cl': 1.32 (RSE 14%, bootstrap median 1.22, 95% CI 0.89-1.7)
    crcl_ref   <- fixed(57);   label("CRCL reference / capping threshold (mL/min/1.73 m^2)")  # Table 2 row 'Cut-off value on kidney function (mL/min)': 57 (RSE 0.1%, bootstrap median 57.4)

    # Body-weight allometric exponents and female sex factor on Vp:
    # fixed at the Wu 2015 published values. The reference weight is
    # 79 kg per Wu 2015 Table 2 (Vc/F = 756 * (WT/79)^1.29; Mulder
    # 2025 uses a slightly different fixed Vc point estimate of 769
    # but the same covariate equation form).
    e_wt_vc    <- fixed(1.29);   label("Allometric exponent of WT on Vc (unitless)")   # Table 2 row 'WGT on V2': 1.29 (fixed) - Wu 2015
    e_wt_vp    <- fixed(0.725);  label("Allometric exponent of WT on Vp (unitless)")   # Table 2 row 'WGT on V3': 0.725 (fixed) - Wu 2015
    e_sexf_vp  <- fixed(0.732);  label("Female (SEXF = 1) multiplicative factor on Vp (unitless)")  # Table 2 row 'Sex on V3': 0.732 (fixed) - Wu 2015 form 0.732^SEX with SEX = 1 for female

    # IIV on CL only (Mulder 2025 Table 2 "IIV CL (%CV)" = 50.5%).
    # Variance on the log-CL eta: omega^2 = log(CV^2 + 1).
    etalcl ~ log(0.505^2 + 1)    # Table 2 row 'IIV CL (%CV)': 50.5 -> omega^2 = log(0.255 + 1) = 0.227

    # Proportional residual error on plasma RBV concentration.
    propSd <- 0.377; label("Proportional residual error on Cc (fraction)")  # Table 2 row 'Standard deviation proportional error': 0.377 (RSE 11%)

    # ================================================================
    # Haemoglobin indirect-response model (Mulder 2025 Table 2
    # "Haemoglobin population model" block; Methods "Pharmacokinetic,
    # haemoglobin and viral load analysis" subsection). Indirect-
    # response on a haemoglobin turnover pool with RBV acting via a
    # linear-in-Cc factor (1 + slope_hb * Cc) on the haemoglobin loss
    # rate constant so haemoglobin declines with increasing RBV
    # exposure. The paper text says "linear inhibitory effect of RBV
    # on the degradation rate (kout)"; the only form that reproduces
    # the observed haemoglobin decline with increasing RBV (the entire
    # point of the simulations in Figures 3-5) is RBV ENHANCING the
    # loss rate, i.e. kout * (1 + slope * Cc) rather than kout * (1 -
    # slope * Cc). This sign-of-effect interpretation is documented
    # in the validation vignette Assumptions and deviations section.
    # ================================================================
    lkout      <- log(0.556);    label("Haemoglobin loss rate constant kout (1/h)")       # Table 2 row 'k out (h-1)': 0.556 (RSE 26%, bootstrap median 0.565, 95% CI 0.308-0.911)
    lslope_hb  <- log(0.102);    label("Linear slope of the RBV effect on kout (L/mg)")    # Table 2 row 'Slope (slope of RBV effect)': 0.102 (RSE 11%, bootstrap median 0.102, 95% CI 0.078-0.129)

    etalkout    ~ log(4.63^2 + 1)   # Table 2 row 'IIV on k out (%CV)': 463 -> omega^2 = log(21.4369 + 1) = log(22.4369) = 3.111 (paper reported %CV very large; shrinkage 31%)
    etalslope_hb ~ log(0.521^2 + 1) # Table 2 row 'IIV on slope (%CV)': 52.1 -> omega^2 = log(0.2714 + 1) = 0.240 (shrinkage 36%)

    addSd_hb   <- 0.406;         label("Additive residual error on Hb (mmol/L)")           # Table 2 row 'Standard deviation additive error': 0.406 (RSE 8%)

    # ================================================================
    # Viral load target-cell-limited model (Mulder 2025 Table 2
    # "Viral load population model" block; Methods "Pharmacokinetic,
    # haemoglobin and viral load analysis" subsection). Three paper-
    # specific compartments: healthy hepatocytes (healthy), infected
    # hepatocytes (infected), virions (virus). The fraction of
    # infected hepatocytes at baseline rho is fixed at 0.001 to
    # prevent unrealistic liver growth on clearance of the virus
    # (Mulder 2025 Methods). The maximum inhibition Imax = 0.999 and
    # IC50 = 1000 ng/L are both fixed; IC50 was tested in a
    # sensitivity analysis and the viral load was found to be
    # insensitive to it because observed RBV concentrations (range
    # 0.1-6.2 mg/L = 100-6200 ug/L) are several orders of magnitude
    # above the fixed IC50. The healthy-hepatocyte half-life TDEG =
    # 6398 h (~ 266 days) and the infected:healthy hepatocyte death-
    # rate ratio Factor = 100 (so the infected-hepatocyte half-life
    # is 6398 / 100 = 63.98 h, kloss = 100 * kdeg) are fixed to the
    # Dahari 2007 literature values (Mulder 2025 reference 20). The
    # only estimated structural parameter in this block is the viral
    # elimination rate elim. The paper does not report a Hill
    # coefficient for the sigmoidal Emax inhibition function; this
    # extraction uses Hill = 1 (basic Emax / Imax form) - see the
    # vignette Assumptions and deviations section.
    # ================================================================
    ltdeg    <- fixed(log(6398));            label("Half-life of healthy hepatocytes TDEG (h)")  # Table 2 row 'TDEG (h)': 6398 (fixed) - Dahari 2007
    lfactor  <- fixed(log(100));             label("Ratio kloss/kdeg (factor for infected hepatocyte half-life, unitless)")  # Table 2 row 'Factor (factor for half-life of infected hepatocytes)': 100 (fixed) - Dahari 2007
    lelim    <- log(0.0123);                  label("Viral elimination rate elim (1/h)")           # Table 2 row 'Elimination rate of virions (h-1)': 0.0123 (RSE 7%, bootstrap median 0.0126, 95% CI 0.0090-0.0171)
    lic50    <- fixed(log(1000 / 1e6));      label("IC50 of RBV on viral replication (mg/L; = 1000 ng/L)")  # Table 2 row 'IC50 (ng/L)': 1000 (fixed); inline conversion 1000 ng/L = 1 ug/L = 1e-3 mg/L (Cc units)
    imax_vl  <- fixed(0.999);                label("Maximum inhibition of viral replication by RBV Imax (fraction)")  # Table 2 row 'Imax': 0.999 (fixed)
    rho      <- fixed(0.001);                label("Fraction of infected hepatocytes at baseline rho (unitless)")  # Table 2 row 'Rho': 0.001 (fixed)

    etalelim ~ log(0.717^2 + 1)              # Table 2 row 'IIV on elimination rate (%CV)': 71.7 -> omega^2 = log(0.5141 + 1) = 0.415 (shrinkage < 30%)

    addSd_vload <- 2.01;        label("Additive residual error on log(viral load) (log-IU/mL)")  # Table 2 row 'Standard deviation additive error': 2.01 (RSE 16.8%); applied to log-transformed viral load via the lvload <- log(virus) output
  })

  model({
    # ================================================================
    # 1. Derived covariate effects (eGFR cap + power, body-weight
    #    allometry, female sex factor on Vp).
    # ================================================================
    crcl_eff <- min(CRCL, crcl_ref) / crcl_ref

    # ================================================================
    # 2. Individual PK parameters. Reference weight 79 kg comes from
    #    Wu 2015 (the upstream PK source); the WT effects on Vc and
    #    Vp and the female factor on Vp are inherited verbatim.
    # ================================================================
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * crcl_eff^e_crcl_cl
    vc <- exp(lvc) * (WT / 79)^e_wt_vc
    vp <- exp(lvp) * (WT / 79)^e_wt_vp * e_sexf_vp^SEXF
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ================================================================
    # 3. PK ODE system (oral, 2-compartment with first-order
    #    absorption). Dose in mg, volumes in L -> central / vc in mg/L.
    # ================================================================
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    Cc <- central / vc

    # ================================================================
    # 4. Haemoglobin indirect-response model. RBV increases the
    #    haemoglobin loss rate linearly with concentration so that
    #    Hb declines with RBV exposure. kin is set to kout * HGB_BL
    #    so that the hb_state is at steady state pre-treatment
    #    (i.e. at t = 0 with Cc = 0, d/dt(hb_state) = 0).
    # ================================================================
    kout     <- exp(lkout + etalkout)
    slope_hb <- exp(lslope_hb + etalslope_hb)

    kin_hb <- kout * HGB_BL

    d/dt(hb_state) <- kin_hb - kout * (1 + slope_hb * Cc) * hb_state
    hb_state(0)    <- HGB_BL
    hb <- hb_state

    # ================================================================
    # 5. Viral load target-cell-limited model. Three paper-specific
    #    compartments: healthy hepatocytes (healthy, arbitrary unit
    #    with H(0) = 1), infected hepatocytes (infected, I(0) = rho =
    #    0.001), virions (virus, V(0) = HEV_VLOAD in IU/mL). The
    #    synthesis, infection-rate, and virion-production parameters
    #    ksyn_h, beta, and p are derived from baseline steady-state
    #    constraints so that the three viral states start at
    #    pre-treatment steady state and only diverge once RBV begins
    #    inhibiting viral replication.
    # ================================================================
    elim   <- exp(lelim + etalelim)
    tdeg   <- exp(ltdeg)
    factor <- exp(lfactor)
    ic50   <- exp(lic50)

    kdeg   <- log(2) / tdeg          # healthy-hepatocyte first-order death rate (1/h)
    kloss  <- factor * kdeg          # infected-hepatocyte first-order death rate (1/h)

    # Baseline steady-state derivations (Mulder 2025 Methods,
    # "Pharmacokinetic, haemoglobin and viral load analysis"):
    #   H(0) = 1 (arbitrary)
    #   I(0) = rho   (fixed at 0.001)
    #   V(0) = HEV_VLOAD
    #   d/dt(I)|t=0 = 0 -> beta * H(0) * V(0) = kloss * I(0)
    #     -> beta = kloss * rho / HEV_VLOAD
    #   d/dt(V)|t=0 = 0 (pre-treatment) -> p * I(0) = elim * V(0)
    #     -> p = elim * HEV_VLOAD / rho
    #   d/dt(H)|t=0 = 0 -> ksyn_h = kdeg * H(0) + beta * H(0) * V(0)
    #                            = kdeg + kloss * rho
    beta_h <- kloss * rho / HEV_VLOAD
    p_v    <- elim  * HEV_VLOAD / rho
    ksyn_h <- kdeg + kloss * rho

    # Inhibitory effect of RBV on viral replication (Imax / IC50
    # sigmoidal Emax form, Hill = 1; the paper reports IC50 in ng/L
    # and Cc is in mg/L, so the inline conversion in ini() (lic50 =
    # log(1000 / 1e6)) puts IC50 in mg/L matching Cc).
    eff_vl <- imax_vl * Cc / (Cc + ic50)

    d/dt(healthy)  <- ksyn_h - kdeg * healthy - beta_h * healthy * virus
    d/dt(infected) <- beta_h * healthy * virus - kloss * infected
    d/dt(virus)    <- p_v * (1 - eff_vl) * infected - elim * virus

    healthy(0)  <- 1
    infected(0) <- rho
    virus(0)    <- HEV_VLOAD

    # Log-transformed viral load observation (1e-30 floor to keep
    # the log finite if numerical noise pushes virus below 0).
    vload <- log(virus + 1e-30)

    # ================================================================
    # 6. Observation models and residual errors.
    # ================================================================
    Cc    ~ prop(propSd)
    hb    ~ add(addSd_hb)
    vload ~ add(addSd_vload)   # additive residual on log-scale -> log-normal residual on linear viral load
  })
}
