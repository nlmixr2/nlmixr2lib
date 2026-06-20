Borghardt_2016_olodaterol <- function() {
  description <- paste(
    "Population PK model for inhaled and intravenous olodaterol (long-acting",
    "beta-2-adrenergic receptor agonist) in 148 healthy adult volunteers from",
    "three Phase I trials (Borghardt 2016). Four-compartment systemic",
    "disposition (central + 3 peripheral) fitted to IV plasma + urine data,",
    "with two parallel first-order elimination processes from the central",
    "compartment: renal (cl_renal) and nonrenal (cl_nonren). For inhaled",
    "administration via the Respimat inhaler, three parallel first-order",
    "absorption depots (slow, intermediate, fast) feed the central",
    "compartment, with absorption half-lives of 21.8 h, 2.00 h, and 0.268 h",
    "respectively. The pulmonary bioavailable fraction (49.4% of the nominal",
    "ex-mouthpiece dose) is split across the three depots by two",
    "logit-transformed proportionality parameters. Smoking is a covariate on",
    "the slow and fast absorption rate constants (active smokers vs",
    "ex-smokers and never-smokers pooled). Systemic disposition parameters",
    "were estimated from IV data and fixed when fitting the inhalation data."
  )
  reference <- paste(
    "Borghardt JM, Weber B, Staab A, Kunz C, Formella S, Kloft C.",
    "Investigating pulmonary and systemic pharmacokinetics of inhaled",
    "olodaterol in healthy volunteers using a population pharmacokinetic",
    "approach.",
    "Br J Clin Pharmacol. 2016;81(3):538-552.",
    "doi:10.1111/bcp.12780.",
    sep = " "
  )
  vignette <- "Borghardt_2016_olodaterol"
  units <- list(time = "hour", dosing = "ug", concentration = "pg/mL")

  covariateData <- list(
    SMOKE = list(
      description        = "Active-smoker binary indicator at trial entry (1 = current smoker, 0 = ex-smoker or never-smoker).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ex-smoker or never-smoker)",
      notes              = paste(
        "Paper grouped ex-smokers with never-smokers based on the observed",
        "lack of difference between those subgroups (Borghardt 2016 Results,",
        "Covariate model paragraph). Encoded here using the canonical SMOKE",
        "register entry (1 = current smoker, 0 otherwise)."
      ),
      source_name        = "Smoking status (active smoker vs ex-smoker/never-smoker)"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 148L,
    n_studies       = 3L,
    age_range       = "22-49 years (median 30.5-34 across trials)",
    weight_range    = "53-105 kg (median 78-81 across trials)",
    sex_female_pct  = 10.1,
    race_ethnicity  = NA,
    disease_state   = "healthy adult volunteers",
    dose_range      = paste(
      "Intravenous infusion (trial 1): 0.5-25 ug single dose over 30 min or",
      "3 h. Oral inhalation via Respimat (trial 2 single rising-dose):",
      "2.5-70 ug. Oral inhalation via Respimat (trial 3 multiple-dose):",
      "2.5, 10, or 30 ug once daily for 14 days. Pulmonary bioavailable",
      "fraction is 49.4% of the nominal ex-mouthpiece dose."
    ),
    regions         = "Germany (trials 1 and 2), Netherlands (trial 3)",
    notes           = paste(
      "Borghardt 2016 Table 1 (trial-by-trial demographics). Pooled cohort:",
      "133 male / 15 female (about 10% female; trial 1 was 100% male).",
      "Smoking status across pooled trials (nonsmoker/ex-smoker/smoker):",
      "approximately 65/16/19 percent. Creatinine clearance medians were",
      "116-121 mL/min across trials. ClinicalTrials.gov: NCT02172131",
      "(trial 1, IV), NCT02171780 (trial 2, single inhaled rising dose),",
      "NCT02171806 (trial 3, multiple-dose inhalation)."
    )
  )

  ini({
    # ============================================================
    # Systemic disposition model parameters (Borghardt 2016 Table 2).
    # Estimated from IV data in trial 1 (Methods, Structural model
    # paragraph; "Parameters were estimated based on intravenous data").
    # Four-compartment model: central + 3 peripheral compartments;
    # two parallel first-order elimination routes (renal, nonrenal)
    # from the central compartment.
    # NONMEM-to-canonical mapping for inter-compartmental terms:
    #   paper V2 (CMT 2) -> canonical peripheral1 / lvp / lq
    #   paper V3 (CMT 3) -> canonical peripheral2 / lvp2 / lq2
    #   paper V4 (CMT 4) -> canonical peripheral3 / lvp3 / lq3
    # ============================================================
    lvc        <- log(23.5);  label("Central volume of distribution VC (L)")                       # Borghardt 2016 Table 2 (VC = 23.5 L, 4.35% RSE)
    lvp        <- log(2590);  label("First peripheral volume V2 (L)")                              # Borghardt 2016 Table 2 (V2 = 2590 L, 35.7% RSE)
    lvp2       <- log(473);   label("Second peripheral volume V3 (L)")                             # Borghardt 2016 Table 2 (V3 = 473 L, 10.7% RSE)
    lvp3       <- log(16.1);  label("Third peripheral volume V4 (L)")                              # Borghardt 2016 Table 2 (V4 = 16.1 L, 19.7% RSE)
    lq         <- log(31.7);  label("Inter-compartmental clearance Q2 (L/h)")                      # Borghardt 2016 Table 2 (Q2 = 31.7 L/h, 12.3% RSE)
    lq2        <- log(65.7);  label("Inter-compartmental clearance Q3 (L/h)")                      # Borghardt 2016 Table 2 (Q3 = 65.7 L/h, 5.28% RSE)
    lq3        <- log(22.5);  label("Inter-compartmental clearance Q4 (L/h)")                      # Borghardt 2016 Table 2 (Q4 = 22.5 L/h, 8.06% RSE)
    lcl_renal  <- log(10.5);  label("Renal clearance CL_R (L/h)")                                  # Borghardt 2016 Table 2 (CL_R = 10.5 L/h, 4.55% RSE)
    lcl_nonren <- log(63.7);  label("Nonrenal clearance CL_NR (L/h)")                              # Borghardt 2016 Table 2 (CL_NR = 63.7 L/h, 8.49% RSE)

    # ============================================================
    # Pulmonary absorption model parameters (Borghardt 2016 Table 2;
    # Results, Structural model paragraph). Three parallel first-order
    # absorption processes feed the systemic central compartment from
    # three depot compartments. Each depot receives a fraction of the
    # nominal inhaled dose multiplied by the pulmonary bioavailable
    # fraction; the fractions are parameterised through two
    # proportionality parameters (FF1 and FF2; Methods, Structural
    # model paragraph: "the number of proportionality parameters was
    # equal to the number of depot compartments representing pulmonary
    # drug absorption minus one").
    #   ka_slow corresponds to ka_slow in the paper
    #   ka_int  corresponds to ka_int  in the paper
    #   ka_fast corresponds to ka_fast in the paper
    # ============================================================
    lka_slow <- log(0.0318);  label("Slow absorption rate constant ka_slow (1/h)")                 # Borghardt 2016 Table 2 (ka_slow = 0.0318 1/h, 5.23% RSE)
    lka_int  <- log(0.347);   label("Intermediate absorption rate constant ka_int (1/h)")          # Borghardt 2016 Table 2 (ka_int  = 0.347 1/h,  7.04% RSE)
    lka_fast <- log(2.59);    label("Fast absorption rate constant ka_fast (1/h)")                 # Borghardt 2016 Table 2 (ka_fast = 2.59 1/h,   9.97% RSE)

    # Logit-transformed pulmonary bioavailable fraction and proportionality
    # parameters (Borghardt 2016 Methods, Random-effects model paragraph:
    # "BSV was assumed to be normally distributed and added to the
    # logit-transformed parameter estimate. In a second step, the parameter
    # associated with BSV was transformed using an inverse logit
    # transformation, so that the resulting individual parameter estimate
    # was again constrained between zero and one"). The same logit
    # transformation is applied to the typical values used here, so that
    # the model parameterisation is on the same scale as the paper.
    #   pbio = 0.494  -> logit = log(0.494 / 0.506) = -0.024
    #   ff1  = 0.701  -> logit = log(0.701 / 0.299) =  0.852
    #   ff2  = 0.889  -> logit = log(0.889 / 0.111) =  2.080
    # Per-depot fractions of the bioavailable dose:
    #   slow         = ff1            = 70.1%  (Borghardt 2016 Table 2/3)
    #   intermediate = (1 - ff1)*ff2  = 26.6%
    #   fast         = (1 - ff1)*(1 - ff2) = 3.32%
    logitpbio <- log(0.494 / (1 - 0.494));   label("Logit of pulmonary bioavailable fraction (PBIO; fraction of nominal ex-mouthpiece dose; unitless)")  # Borghardt 2016 Table 2 (PBIO = 49.4% of ND, 3.32% RSE)
    logitff1  <- log(0.701 / (1 - 0.701));   label("Logit of first proportionality parameter FF1 (slow-absorption fraction of pulmonary bioavailable dose; unitless)")  # Borghardt 2016 Table 2 (FF1 = 70.1%, 2.33% RSE)
    logitff2  <- log(0.889 / (1 - 0.889));   label("Logit of second proportionality parameter FF2 (intermediate fraction within the non-slow remainder; unitless)")  # Borghardt 2016 Table 2 (FF2 = 88.9%, 1.89% RSE)

    # ============================================================
    # Smoking covariate effects on the slow and fast absorption rate
    # constants (Borghardt 2016 Table 2 and Results, Covariate model
    # paragraph). Modelled multiplicatively on the linear (back-
    # transformed) rate constant:
    #   ka_slow_smoker = ka_slow * (1 + e_smoke_ka_slow * SMOKE)
    #   ka_fast_smoker = ka_fast * (1 + e_smoke_ka_fast * SMOKE)
    # The reported smoker-modified values in Table 3 (ka_slow_smoker
    # = 0.0197 1/h; ka_fast_smoker = 5.39 1/h) reproduce the paper's
    # numbers under this parameterisation:
    #   0.0318 * (1 + (-0.380)) = 0.01972
    #   2.59   * (1 + ( 1.08 )) = 5.387
    # Ex-smokers and never-smokers are pooled (SMOKE = 0) per Results.
    # ============================================================
    e_smoke_ka_slow <- -0.380;  label("Smoking effect on slow absorption rate constant (relative change; unitless)")  # Borghardt 2016 Table 2 (Impact of active smoking on ka_slow = -0.380, 10.2% RSE)
    e_smoke_ka_fast <-  1.08;   label("Smoking effect on fast absorption rate constant (relative change; unitless)")  # Borghardt 2016 Table 2 (Impact of active smoking on ka_fast =  1.08,  18.0% RSE)

    # ============================================================
    # Between-subject variability (BSV; Borghardt 2016 Table 2,
    # BSV [% CV] column). Reported as %CV on the linear scale; log-
    # scale variance computed as omega^2 = log(1 + (CV/100)^2). For
    # logit-transformed parameters (logitff1, logitff2, logitpbio)
    # the same %CV is applied directly on the logit scale per the
    # paper's parameterisation; the back-transformed CV via the
    # delta method is approximate but follows the published Random-
    # effects model description and the Weber 2015 fluticasone
    # precedent in this package. The IIV terms are diagonal (no
    # inter-parameter correlation; Methods, Random-effects model
    # paragraph: "Additional covariances did not significantly
    # improve the pulmonary absorption model").
    # NOTE on PBIO: Borghardt 2016 reports BOV (32.2% CV) rather than
    # BSV on the pulmonary bioavailable fraction; the rxode2 forward-
    # simulation pipeline samples a single random effect per subject
    # and so cannot distinguish BSV from BOV without a per-occasion
    # event grouping. The 32.2% magnitude is therefore encoded as
    # BSV (etalogitpbio) here for simulation purposes and the
    # deviation is documented in the vignette's Assumptions and
    # deviations section.
    # ============================================================
    etalvc          ~ 0.0665   # Borghardt 2016 Table 2 (BSV VC      26.2% CV -> log(1 + 0.262^2) = 0.0665)
    etalq           ~ 0.0639   # Borghardt 2016 Table 2 (BSV Q2      25.7% CV -> log(1 + 0.257^2) = 0.0639)
    etalq2          ~ 0.0279   # Borghardt 2016 Table 2 (BSV Q3      16.8% CV -> log(1 + 0.168^2) = 0.0279)
    etalcl_nonren   ~ 0.0693   # Borghardt 2016 Table 2 (BSV CL_NR   26.8% CV -> log(1 + 0.268^2) = 0.0693)
    etalogitff1     ~ 0.0129   # Borghardt 2016 Table 2 (BSV FF1     11.4% CV -> log(1 + 0.114^2)  = 0.0129; applied on the logit scale)
    etalogitff2     ~ 0.00866  # Borghardt 2016 Table 2 (BSV FF2      9.33% CV -> log(1 + 0.0933^2) = 0.00866; applied on the logit scale)
    etalogitpbio    ~ 0.0987   # Borghardt 2016 Table 2 (BOV PBIO    32.2% CV -> log(1 + 0.322^2)  = 0.0987; encoded as IIV here; see Errata)

    # ============================================================
    # Residual error (Borghardt 2016 Table 2, plasma rows). The
    # plasma additive component (0.00001 pg/mL fixed) is structurally
    # present to stabilise the NONMEM estimation at very low
    # concentrations; numerically it is negligible relative to the
    # 15.8% proportional component. The paper's urine residual error
    # (37.7% CV proportional plus 0.00001 pg/mL additive) is not
    # encoded here because this model exposes only the plasma
    # observation Cc; the urine ODE state tracks cumulative renal
    # excretion for users who want to inspect it but is not declared
    # as a model output. See vignette Assumptions and deviations.
    # ============================================================
    propSd <- 0.158;          label("Proportional residual error on plasma Cc (fraction)")  # Borghardt 2016 Table 2 (plasma proportional RSV = 15.8% CV, 1.66% RSE)
    addSd  <- fixed(0.00001); label("Additive residual error on plasma Cc (pg/mL; fixed)")  # Borghardt 2016 Table 2 (plasma additive RSV = 0.00001 pg/mL, fixed)
  })

  model({
    # ------------------------------------------------------------
    # 1. Individual structural parameters. Log-transforms back-
    #    transform via exp(); logit-transformed parameters back-
    #    transform via the inverse-logit 1 / (1 + exp(-x)).
    # ------------------------------------------------------------
    vc        <- exp(lvc        + etalvc)
    vp        <- exp(lvp)
    vp2       <- exp(lvp2)
    vp3       <- exp(lvp3)
    q         <- exp(lq         + etalq)
    q2        <- exp(lq2        + etalq2)
    q3        <- exp(lq3)
    cl_renal  <- exp(lcl_renal)
    cl_nonren <- exp(lcl_nonren + etalcl_nonren)

    # Smoker-modified absorption rate constants. SMOKE = 1 selects the
    # smoker rate; SMOKE = 0 leaves the rate at the typical value.
    ka_slow <- exp(lka_slow) * (1 + e_smoke_ka_slow * SMOKE)
    ka_int  <- exp(lka_int)
    ka_fast <- exp(lka_fast) * (1 + e_smoke_ka_fast * SMOKE)

    # Logit back-transform for pulmonary bioavailable fraction and the
    # two proportionality parameters. Per-subject etas are added on the
    # logit scale before the inverse-logit transformation, matching the
    # paper's parameterisation.
    logit_pbio_i <- logitpbio + etalogitpbio
    logit_ff1_i  <- logitff1  + etalogitff1
    logit_ff2_i  <- logitff2  + etalogitff2
    pbio <- 1 / (1 + exp(-logit_pbio_i))
    ff1  <- 1 / (1 + exp(-logit_ff1_i))
    ff2  <- 1 / (1 + exp(-logit_ff2_i))

    # Per-depot fractions of the pulmonary bioavailable dose.
    frac_slow <- ff1
    frac_int  <- (1 - ff1) * ff2
    frac_fast <- (1 - ff1) * (1 - ff2)

    # ------------------------------------------------------------
    # 2. Micro-constants for the linear 4-compartment systemic
    #    disposition with two parallel elimination routes (renal,
    #    nonrenal) from the central compartment.
    # ------------------------------------------------------------
    kel_renal  <- cl_renal  / vc
    kel_nonren <- cl_nonren / vc
    kel        <- kel_renal + kel_nonren
    k12        <- q  / vc
    k21        <- q  / vp
    k13        <- q2 / vc
    k31        <- q2 / vp2
    k14        <- q3 / vc
    k41        <- q3 / vp3

    # ------------------------------------------------------------
    # 3. ODE system. Three parallel absorption depots feed the
    #    central compartment (combined PK model; Borghardt 2016
    #    Figure 2 top box). For IV administration the dose enters
    #    `central` directly via the event table; the depots stay
    #    at zero.
    # ------------------------------------------------------------
    d/dt(depot)       <- -ka_slow * depot
    d/dt(depot2)      <- -ka_int  * depot2
    d/dt(depot3)      <- -ka_fast * depot3
    d/dt(central)     <-  ka_slow * depot + ka_int * depot2 + ka_fast * depot3 -
                          (kel + k12 + k13 + k14) * central +
                          k21 * peripheral1 +
                          k31 * peripheral2 +
                          k41 * peripheral3
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2
    d/dt(peripheral3) <-  k14 * central - k41 * peripheral3
    d/dt(urine)       <-  kel_renal * central

    # ------------------------------------------------------------
    # 4. Bioavailability. For inhaled administration via Respimat,
    #    each dose event splits across the three depots with
    #    fractions (pbio * frac_slow, pbio * frac_int, pbio * frac_fast).
    #    For IV administration the depots receive dose = 0 (or are
    #    skipped entirely in the event table) and `central` takes
    #    the full dose with default bioavailability of 1.
    # ------------------------------------------------------------
    f(depot)  <- pbio * frac_slow
    f(depot2) <- pbio * frac_int
    f(depot3) <- pbio * frac_fast

    # ------------------------------------------------------------
    # 5. Observation. Plasma olodaterol concentration in pg/mL.
    #    central [ug] / vc [L] = ug/L = 1000 pg/mL.
    # ------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
