Brekkan_2018_pegfilgrastim <- function() {
  description <- "Bidirectional population PK/PD model for pegfilgrastim (PG) in healthy volunteers after single 6 mg subcutaneous doses. PK is one-compartment with sequential zero- and first-order absorption (zero-order input rate R1 into depot followed by first-order Ka into central), and parallel elimination via a linear ANC-dependent pathway (cl_anc * ANC) and a saturable non-specific Michaelis-Menten pathway (Vmax / Km). PD is a Friberg/Quartino-style maturation cascade (4 transit compartments + circulating compartment) with the production rate set by baseline ANC and a fixed 7-hour circulating neutrophil half-life, plus three Emax drug effects: proliferation (scaling production), maturation (scaling transit rate ktr), and a margination effect on circulating-pool clearance kcirc that is parameterised by scaled Emax,prol and EC50 (Emax,Scale = 0.0622, EC50,Scale = 0.477)."
  reference <- paste(
    "Brekkan A, Lopez-Lazaro L, Yngman G, Plan EL, Acharya C, Hooker AC, Kankanwadi S, Karlsson MO. (2018).",
    "A Population Pharmacokinetic-Pharmacodynamic Model of Pegfilgrastim.",
    "AAPS J 20(5):91. doi:10.1208/s12248-018-0249-y."
  )
  vignette <- "Brekkan_2018_pegfilgrastim"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL", ANC = "10^9 cells/L")

  covariateData <- list(
    # Covariates were evaluated using a full random effects model (FREM)
    # approach with sex, age, lean body weight (LBW = LBM), body mass index,
    # weight, and race (time-constant) plus formulation and period
    # (time-varying). Per the Brekkan 2018 Discussion, no individual
    # covariate explained more than 2% of variability in the derived
    # metrics AUC and Cmax (Appendix 2). The structural Table I parameters
    # used in this model file do not depend on covariates; the covariate
    # tests are documented here for transparency but the effects from the
    # Figs 6-9 forest plots are not encoded in model() because the
    # forest-plot values are graphical only and the paper concludes the
    # effects are not clinically meaningful.
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-constant covariate tested in the FREM analysis. Forest plots (Brekkan 2018 Figs 8-9) show small effects on Vc/F and R1; not retained in the final structural model. Reference value not tabulated in the paper.",
      source_name        = "WT"
    ),
    LBM = list(
      description        = "Lean body mass (Janmahasatian et al.)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Brekkan 2018 uses the alias 'LBW' (lean body weight); canonical column in nlmixr2lib is LBM. Computed via the Janmahasatian formula (Brekkan 2018 Methods cites Janmahasatian et al. 2005). Tested via FREM as a time-constant covariate; second most influential covariate after sex per the Discussion but still explains less than 2% of derived-metric variability. Reference value not tabulated.",
      source_name        = "LBW"
    ),
    BMI = list(
      description        = "Body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-constant covariate tested in the FREM analysis (Brekkan 2018 Methods). Not retained in the final structural model.",
      source_name        = "BMI"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-constant covariate tested in the FREM analysis. Not retained.",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Brekkan 2018 reports 43% female (Results, Data). Tested via FREM as a time-constant covariate; most influential covariate by FREM forest plots but still explains a small fraction of overall variability.",
      source_name        = "SEX"
    ),
    RACE_WHITE = list(
      description        = "Race indicator (1 = White / Caucasian, 0 = other)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Caucasian)",
      notes              = "Brekkan 2018 reports 86% Caucasian (Results, Data). Tested via FREM as a time-constant covariate; effect modest in forest plots.",
      source_name        = "RACE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 174L,
    n_studies      = 1L,
    n_occasions    = 445L,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = 43,
    race_ethnicity = c(Caucasian = 86, Other = 14),
    disease_state  = "Healthy adult volunteers, drug-naive (no measurable baseline PG, no anti-PEG or anti-drug antibody positivity).",
    dose_range     = "Single 6 mg subcutaneous pegfilgrastim per period in a three-way crossover comparing one biosimilar candidate (BIOS_PG) and two reference Neulasta(R) products (US_PG and EU_PG), separated by at least 5 weeks washout. Median delivered doses by syringe batch: BIOS_PG 6.6 and 6.7 mg (two batches), US_PG 6.3 mg, EU_PG 6.2 mg (Brekkan 2018 Results, Data).",
    regions        = "USA",
    notes          = "Study PG-01-003 (Dr. Reddy's Laboratories). 192 enrolled; 174 included after excluding 15 subjects with measurable baseline PG and ADA-positive occasions (Brekkan 2018 Results, Data). 445 dosing occasions analysed. The paper does not tabulate age or weight ranges; baseline demographics are referenced via Fig 1 and prose only."
  )

  ini({
    # -----------------------------------------------------------------
    # Structural PK -- Brekkan 2018 Table I (PKPD model parameter estimates).
    # Reference subject is the typical healthy adult (no covariates retained
    # in the final structural model). Bioavailability is anchored at Frel=1
    # because no intravenous data were available (Methods, Pharmacokinetic/
    # Pharmacodynamic Model Development).
    # -----------------------------------------------------------------
    lka     <- log(0.0114) ; label("First-order absorption rate Ka (1/h)")                            # Brekkan 2018 Table I (Ka 0.0114, RSE 2.78%)
    lr1     <- log(2.84)   ; label("Zero-order input rate R1 into depot (mg/h)")                      # Brekkan 2018 Table I (R1 2.84 mg/h, RSE 23.8%)
    lfdepot <- fixed(log(1)); label("Relative bioavailability Frel into depot (fixed)")                # Brekkan 2018 Table I (Frel 1, Fixed; no IV data)
    lvc     <- log(1.81)   ; label("Apparent central volume Vc/F (L)")                                 # Brekkan 2018 Table I (Vc/F 1.81 L, RSE 7.31%)
    lvmax   <- log(0.0467) ; label("Maximum saturable elimination rate Vmax/F (mg/h)")                 # Brekkan 2018 Table I (Vmax/F 0.0467 mg/h, RSE 0.726%)
    lkm     <- log(2.05)   ; label("Michaelis-Menten constant Km (ng/mL)")                             # Brekkan 2018 Table I (Km 2.05 ng/mL, RSE 3.38%)
    lcl_anc <- log(2.80)   ; label("Linear ANC-dependent clearance coefficient (L/h per 10^9 cells/L)") # Brekkan 2018 Table I (CL/F 2.80, RSE 13.0%)

    # -----------------------------------------------------------------
    # Structural PD -- Brekkan 2018 Table I.
    # Mean maturation time (MMT) and circulating-neutrophil half-life
    # (T1/2) are fixed at literature values (Methods, Pharmacokinetic/
    # Pharmacodynamic Model Development): MMT = 120 h (~5 days for
    # cells to mature and migrate from bone marrow to bloodstream);
    # T1/2 = 7 h (frequently reported half-life of circulating
    # neutrophils). The only estimated PD turnover parameter is BASE.
    # -----------------------------------------------------------------
    lbase  <- log(2.70)        ; label("Baseline neutrophil count BASE (10^9 cells/L)")               # Brekkan 2018 Table I (BASE 2.70, RSE 3.75%)
    lthalf <- fixed(log(7))    ; label("Circulating neutrophil half-life T1/2 (h, fixed)")             # Brekkan 2018 Table I (T1/2 7 h, Fixed)
    lmmt   <- fixed(log(120))  ; label("Mean maturation time MMT (h, fixed)")                          # Brekkan 2018 Table I (MMT 120 h, Fixed)

    # Drug-effect parameters -- Brekkan 2018 Table I and Eqs. 2-4 (page 5).
    # EC50 is shared between the proliferation and the maturation effects;
    # the margination effect uses (Emax,prol * Emax,Scale) and
    # (EC50 * EC50,Scale) so that a single Emax,prol / EC50 pair is
    # estimated and the margination magnitude is scaled relative to it
    # (page 5: "to achieve a parsimonious model, a scaling parameter was
    # estimated to scale both Emax,prol and EC50").
    lec50       <- log(9.24)   ; label("EC50 for proliferation and maturation effects (ng/mL)")          # Brekkan 2018 Table I (EC50 9.24 ng/mL, RSE 7.57%)
    lemax_mat   <- log(102)    ; label("Emax,mat: maximum proportional effect on maturation rate (unitless)")    # Brekkan 2018 Table I (Emax,mat 102, RSE 8.78%)
    lemax_prol  <- log(109)    ; label("Emax,prol: maximum proportional effect on proliferation rate (unitless)") # Brekkan 2018 Table I (Emax,prol 109, RSE 7.64%)
    lec50_scale <- log(0.477)  ; label("EC50 scaling factor for margination effect (unitless)")          # Brekkan 2018 Table I (EC50Scale 0.477, RSE 5.78%)
    lemax_scale <- log(0.0622) ; label("Emax,prol scaling factor for margination effect (unitless)")     # Brekkan 2018 Table I (Emax,Scale 0.0622, RSE 1.23%)

    # -----------------------------------------------------------------
    # Inter-individual variability -- Brekkan 2018 Table I.
    # Footnote a: "RSE for IIV and RUV parameters are reported on the
    # approximate SD scale"; the IIV column itself reports SD. Internal
    # variance is therefore SD^2. Brekkan 2018 estimated a full
    # variance-covariance block for PK (with off-diagonal correlations
    # tabulated in Table III) and a 3x3 Omega block for (EC50,
    # Emax,mat, Emax,prol); PK off-diagonal correlations are NOT
    # captured here for simplicity (see vignette Assumptions and
    # deviations). IIV on Vmax was reported as "unsupported" by the
    # data and is excluded.
    # Inter-occasion variability (IOV) on Ka, R1, Frel and T1/2 is also
    # reported in Table I but is NOT separated from IIV in this library
    # implementation; see vignette Assumptions and deviations.
    # -----------------------------------------------------------------
    etalka     ~ 0.05336 ; label("IIV variance on lka")        # Brekkan 2018 Table I (IIV Ka SD 0.231; var = 0.231^2)
    etalvc     ~ 0.28196 ; label("IIV variance on lvc")        # Brekkan 2018 Table I (IIV Vc SD 0.531; var = 0.531^2)
    etalcl_anc ~ 0.25503 ; label("IIV variance on lcl_anc")    # Brekkan 2018 Table I (IIV CL SD 0.505; var = 0.505^2)
    etalkm     ~ 0.27144 ; label("IIV variance on lkm")        # Brekkan 2018 Table I (IIV Km SD 0.521; var = 0.521^2)
    etalr1     ~ 0.10433 ; label("IIV variance on lr1")        # Brekkan 2018 Table I (IIV R1 SD 0.323; var = 0.323^2)
    etalfdepot ~ 0.04580 ; label("IIV variance on lfdepot")    # Brekkan 2018 Table I (IIV Frel SD 0.214; var = 0.214^2)

    etalbase   ~ 0.06350 ; label("IIV variance on lbase")      # Brekkan 2018 Table I (IIV BASE SD 0.252; var = 0.252^2)
    etalthalf  ~ 0.07508 ; label("IIV variance on lthalf")     # Brekkan 2018 Table I (IIV T1/2 SD 0.274; var = 0.274^2)

    # 3x3 Omega block on (EC50, Emax,mat, Emax,prol) per Brekkan 2018
    # Pharmacokinetic/Pharmacodynamic Model Development:
    # "IIV in the PD model was included on BASE, neutrophil half-life,
    # EC50, Emax,mat, and Emax,prol, the latter three in an Omega block."
    # Diagonal variances from Table I IIV SDs (0.817, 0.417, 0.892);
    # off-diagonals from Table III FREM IIV correlation coefficients:
    #   r(EC50, Emax,mat)   = 0.27  -> cov = 0.27 * sqrt(0.66745 * 0.17389) = 0.09195
    #   r(EC50, Emax,prol)  = 0.58  -> cov = 0.58 * sqrt(0.66745 * 0.79566) = 0.42245
    #   r(Emax,mat, prol)   = 0.24  -> cov = 0.24 * sqrt(0.17389 * 0.79566) = 0.08926
    etalec50 + etalemax_mat + etalemax_prol ~
      c(0.66745,
        0.09195, 0.17389,
        0.42245, 0.08926, 0.79566)                              # Brekkan 2018 Table I (IIV SDs 0.817, 0.417, 0.892) + Table III (correlations 0.27, 0.58, 0.24)

    # -----------------------------------------------------------------
    # Residual unexplained variability -- Brekkan 2018 Table I.
    # Methods (Pharmacokinetic/Pharmacodynamic Model Development) states
    # "Residual errors in the PK and PD models were additive on the log
    # scale, resulting in proportional residual errors on a normal
    # scale." RUV values reported as approximate SD.
    # -----------------------------------------------------------------
    propSd     <- 0.272 ; label("Proportional residual error on PG concentration (fraction)")  # Brekkan 2018 Table I (RUV PK SD 0.272, RSE 3.98%)
    propSd_ANC <- 0.255 ; label("Proportional residual error on ANC (fraction)")               # Brekkan 2018 Table I (RUV PD SD 0.255, RSE 2.75%)
  })

  model({
    # Individual PK parameters
    ka     <- exp(lka     + etalka)
    vc     <- exp(lvc     + etalvc)
    vmax   <- exp(lvmax)                        # no IIV on Vmax (Brekkan 2018: "IIV ... included as log-normal distributions on all parameters in the PK model apart from Vmax (unsupported)")
    km     <- exp(lkm     + etalkm)
    cl_anc <- exp(lcl_anc + etalcl_anc)
    r1     <- exp(lr1     + etalr1)
    frel   <- exp(lfdepot + etalfdepot)

    # Individual PD parameters
    base_anc   <- exp(lbase      + etalbase)
    thalf_anc  <- exp(lthalf     + etalthalf)
    mmt        <- exp(lmmt)
    ec50       <- exp(lec50      + etalec50)
    emax_mat   <- exp(lemax_mat  + etalemax_mat)
    emax_prol  <- exp(lemax_prol + etalemax_prol)
    ec50_scale <- exp(lec50_scale)
    emax_scale <- exp(lemax_scale)

    # Maturation cascade kinetics. Brekkan 2018 has four transit
    # compartments + a circulating compartment with no separate
    # proliferation compartment (production enters the first transit
    # compartment at a constant rate). The 4 transits each have mean
    # residence time 1 / ktr, so the total chain transit time is
    # 4 / ktr. Equating this to MMT = 120 h (~5 days) gives
    # ktr = 4 / MMT, mirroring the Friberg 2002 / Quartino 2014
    # ktr = (N+1)/MTT convention applied to N = 3 transits + 1 prol.
    # The paper's loose phrasing "the transit rate ... (the inverse of
    # the MMT)" (page 4) is interpreted here as referring to the chain's
    # overall maturation completion rate; the per-compartment ktr in
    # nlmixr2 must be N/MMT to give a 5-day total transit time.
    ktr   <- 4 / mmt
    kcirc <- log(2) / thalf_anc

    # Baseline production rate kin, set so the chain is at steady state
    # with each pool equal to BASE.
    kin <- kcirc * base_anc

    # PG concentration in ng/mL: central[mg] / vc[L] = mg/L; multiply by
    # 1000 to convert to ng/mL (the unit of Km and EC50 in Table I).
    Cc <- central / vc * 1000

    # Drug effects on the PD system -- Brekkan 2018 Eqs 2-4 (page 5).
    #   Maturation effect:   Emax,mat * Cc / (EC50 + Cc)
    #   Proliferation effect: Emax,prol * Cc / (EC50 + Cc)
    #   Margination effect:   (Emax,prol * Emax,Scale * Cc) /
    #                         (EC50 * EC50,Scale + Cc)
    prol_eff <- emax_prol * Cc / (ec50 + Cc)
    mat_eff  <- emax_mat  * Cc / (ec50 + Cc)
    marg_eff <- (emax_prol * emax_scale * Cc) /
                (ec50 * ec50_scale + Cc)

    # PK ODEs. The depot receives a zero-order input at rate R1 (mg/h)
    # followed by first-order absorption Ka into central. Elimination
    # is the sum of (i) a linear ANC-dependent clearance
    # cl_anc * circ * (central / vc) -- units check: cl_anc
    # [L/h per 10^9 cells/L] * circ [10^9 cells/L] * (central / vc)
    # [mg/L] = mg/h -- and (ii) a saturable non-specific Michaelis-
    # Menten elimination vmax * Cc / (km + Cc) -- units check:
    # vmax [mg/h] * Cc [ng/mL] / (km + Cc) [ng/mL] = mg/h (the ratio is
    # dimensionless). The circulating-neutrophil state circ couples PK
    # and PD bidirectionally: PG increases ANC via prol_eff / mat_eff /
    # marg_eff, and ANC increases PG elimination via cl_anc * circ.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot -
                      cl_anc * circ * (central / vc) -
                      vmax * Cc / (km + Cc)

    # PD ODEs: 4 transit compartments + circulating ANC. Drug effects
    # scale the production rate (prol_eff), the transit rate ktr
    # (mat_eff), and the circulating-pool first-order clearance kcirc
    # (marg_eff). The margination effect represents an "apparent volume
    # expansion" / rapid sequestration of circulating neutrophils per
    # the Discussion (page 7); encoding it as an increase in apparent
    # kcirc reproduces the transient drop in circulating ANC observed
    # in the early post-dose period (see Fig 1).
    d/dt(precursor1) <- kin * (1 + prol_eff) - ktr * (1 + mat_eff) * precursor1
    d/dt(precursor2) <- ktr * (1 + mat_eff) * (precursor1 - precursor2)
    d/dt(precursor3) <- ktr * (1 + mat_eff) * (precursor2 - precursor3)
    d/dt(precursor4) <- ktr * (1 + mat_eff) * (precursor3 - precursor4)
    d/dt(circ)       <- ktr * (1 + mat_eff) * precursor4 -
                        kcirc * (1 + marg_eff) * circ

    # Steady-state initial conditions. At flux balance with constant
    # production kin = kcirc * BASE flowing through the transit chain
    # and out through the circulating-pool first-order clearance kcirc,
    # the chain's steady-state amounts are NOT uniform when ktr != kcirc.
    # Mass-balance gives:
    #   precursor_i(0) = kin / ktr     (all four transit pools)
    #   circ(0)        = kin / kcirc   = base_anc
    # With ktr = 4/MMT = 0.0333 /h and kcirc = ln(2)/7 = 0.099 /h, the
    # transit-pool baseline is roughly 3x the circulating baseline.
    precursor1(0) <- kin / ktr
    precursor2(0) <- kin / ktr
    precursor3(0) <- kin / ktr
    precursor4(0) <- kin / ktr
    circ(0)       <- base_anc

    # Bioavailability and zero-order absorption rate into depot.
    # rxode2 reads rate() at solve time and converts a bolus dose to a
    # constant-rate input until the dose amount is depleted; the
    # duration is therefore (amt * frel) / r1.
    f(depot)    <- frel
    rate(depot) <- r1

    # Observations: PG plasma concentration Cc (ng/mL) and circulating
    # absolute neutrophil count ANC (10^9 cells/L), each with a
    # proportional residual error.
    Cc  ~ prop(propSd)
    ANC <- circ
    ANC ~ prop(propSd_ANC)
  })
}
