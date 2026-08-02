Roepcke_2018_tak_079 <- function() {
  description <- paste(
    "Preclinical (cynomolgus monkey). Two-compartment QSS-TMDD population PK",
    "model plus three PK-PD lymphocyte-depletion models (NK-cell turnover with",
    "Emax on elimination; B-cell four-transit-compartment chain with Emax on",
    "the circulating depletion rate; T-cell direct-response with Emax) for the",
    "cytolytic anti-CD38 human IgG1-kappa monoclonal antibody TAK-079 pooled",
    "from eight monkey studies (dose 0.03-100 mg/kg IV / SC). PK covariate:",
    "route of administration on Vc (SC vs IV).",
    sep = " "
  )
  reference <- paste(
    "Roepcke S, Plock N, Yuan J, Fedyk ER, Lahu G, Zhao L, Smithson G.",
    "Pharmacokinetics and pharmacodynamics of the cytolytic anti-CD38",
    "human monoclonal antibody TAK-079 in monkey - model assisted",
    "preparation for the first in human trial.",
    "Pharmacol Res Perspect. 2018;6(3):e00402.",
    "doi:10.1002/prp2.402. PMID 29864242.",
    sep = " "
  )
  vignette <- "Roepcke_2018_tak_079"
  paper_specific_compartments <- c("nkcell")

  units <- list(
    time          = "day",
    dosing        = "mg",
    concentration = "ug/mL"
  )

  covariateData <- list(
    ROUTE_IV = list(
      description        = "Route of administration indicator (1 = IV, 0 = SC).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (IV)",
      notes              = paste(
        "Per-subject route covariate carrying the paper's route-of-administration",
        "effect on central volume of distribution. IV is the reference category",
        "(Vc = exp(lvc) = 0.141 L); SC scales Vc by (1 - e_route_iv_vc) = 0.303",
        "-> 0.043 L, ca. 70% smaller (Roepcke 2018 Results Section 3.2 and",
        "Table 2 'ROUT on V_C'). SC data come from the four lower single-dose",
        "groups (<= 1 mg/kg) in studies 7 and 8. For simulation, set",
        "ROUTE_IV = 1 for IV cohorts (dose into central) and 0 for SC cohorts",
        "(dose into depot). Distinct from the dosing-event `cmt` column, which",
        "names the target compartment.",
        sep = " "
      ),
      source_name        = "ROUT"
    )
  )

  population <- list(
    species        = "cynomolgus monkey (Macaca fascicularis)",
    n_subjects     = 140L,
    n_studies      = 8L,
    n_male         = 58L,
    n_female       = 82L,
    weight_range   = "2.1-4.7 kg (pooled PK dataset)",
    weight_typical = "2.6 kg (reference weight used for human-scaling; Berkeley Madonna code and Supplemental 'Scaling of monkey PK parameters')",
    disease_state  = "Healthy cynomolgus monkeys (placebo groups excluded from PK dataset).",
    dose_range     = paste(
      "0.03-100 mg/kg TAK-079. Single-dose IV bolus / 30-min IV infusion",
      "(studies 1-4 infusion; studies 5-8 bolus) and repeat-dose QW / Q2W IV",
      "(studies 1-6). SC single-dose 0.03, 0.1, 0.3, 1 mg/kg in one group of",
      "study 7 and three groups of study 8 (15 animals total).",
      sep = " "
    ),
    n_observations_pk = 2199L,
    ada_status     = paste(
      "Anti-drug antibodies (ADA) developed over time in repeat-dose studies;",
      "229 ADA-affected PK observations were flagged (ADAF=1) and excluded",
      "during model development (Roepcke 2018 Section 3.1 and Supplemental",
      "'Data set preparation'). Not a covariate in the final model.",
      sep = " "
    ),
    notes          = paste(
      "Data pooled from eight preclinical monkey studies (Roepcke 2018 Table 1):",
      "studies 1-6 were toxicology / repeat-dose (weekly QW or every-other-week",
      "Q2W); studies 2, 7, 8 were single-dose PK / PD. Placebo animals",
      "contributed cell-count data only (no drug exposure) and are included in",
      "the pooled PD datasets. In study 5 a dosing error assigned 0.01 mg/kg",
      "(instead of 0.1 mg/kg) at the second dose in the lowest-dose group;",
      "actual dosing amounts were used. Body weights and per-subject data by sex",
      "are summarised in Section 3.1 of the paper; sex was not retained as a",
      "covariate in the final PK model.",
      sep = " "
    )
  )

  ini({
    # =========================================================================
    # PK model -- Roepcke 2018 Table 2 and Figure 2 (final 2-compartment
    # QSS-TMDD structure). Point estimates are the paper-reported typical
    # values; %SE columns are captured as inline comments.
    # =========================================================================

    # Bioavailability: paper models F via logit transform F = exp(PAR)/(1+exp(PAR))
    # so that estimates lie in (0, 1). Table 2 reports PAR = 0.227 (%SE 121).
    # Berkeley Madonna F1 = 0.5565979 = exp(0.227366)/(1 + exp(0.227366)) is the
    # transformed bioavailability actually used in the simulations.
    logitfdepot <- 0.227;    label("Logit-transformed SC bioavailability (F = exp(PAR)/(1+exp(PAR)) ~ 0.557)")  # Table 2: F=0.227 (%SE 121)

    lka         <- log(0.399);   label("Absorption rate constant Ka after SC dose (1/day)")            # Table 2: Ka=0.399 /day (%SE 20.5); Berkeley Madonna KA=0.398734
    lcl         <- log(0.0187);  label("Systemic clearance CL (L/day) for IV reference monkey")        # Table 2: CL=0.0187 L/day (%SE 5.17); Berkeley Madonna CL=0.0187399
    lvc         <- log(0.141);   label("Central volume of distribution Vc (L) for IV reference monkey")  # Table 2: Vc=0.141 L (%SE 3.23); Berkeley Madonna V1=0.140556
    lq          <- log(0.127);   label("Inter-compartmental clearance Q (L/day)")                       # Table 2: Q=0.127 L/day (%SE 14.7); Berkeley Madonna Q=0.127358
    lvp         <- log(0.127);   label("Peripheral volume of distribution Vp (L)")                      # Table 2: Vp=0.127 L (%SE 6.45); Berkeley Madonna V2=0.127326

    # TMDD-QSS receptor parameters (paper Figure 2 blue box + Table 2).
    # KINT is fixed at 0.1/day (Table 2 %SE 'FIXED'); KSYN is fixed at
    # 0.04 units/L/day. Kss and Kdeg were estimated from the low-dose
    # single-dose studies 7 and 8 and then held fixed while fitting on
    # the entire dataset (Section 3.2 'Model development').
    lkint       <- fixed(log(0.1));    label("Drug-target complex internalisation rate KINT (1/day)")  # Table 2: KINT=0.1 /day FIXED; Berkeley Madonna KINT=0.1
    lkss        <- log(5.68);          label("QSS binding constant KSS = (KOFF + KINT)/KON (ug/mL)")           # Table 2: KSS=5.68 (%SE 38.7); Berkeley Madonna KSS=5.68
    lksyn       <- fixed(log(0.04));   label("CD38 receptor synthesis rate KSYN (u/L per day)")         # Table 2: KSYN=0.04 u/L/day FIXED; Berkeley Madonna KSYN=0.04
    lkdeg       <- log(0.00452);       label("CD38 receptor degradation rate KDEG (1/day)")                    # Table 2: KDEG=0.00452 /day (%SE 30.1); Berkeley Madonna KDEG=0.00452

    # Route-of-administration effect on Vc: source paper models Vc_SC =
    # Vc_IV * (1 - 0.697). With canonical ROUTE_IV=1 (IV, reference), rewrite
    # as Vc = Vc_typical * (1 - e_route_iv_vc * (1 - ROUTE_IV)).
    e_route_iv_vc <- 0.697;  label("Fractional reduction of Vc when SC administered (SC vs IV reference)")  # Table 2: ROUT on V_C = 0.697 (%SE 6.51); Berkeley Madonna 0.696753

    # =========================================================================
    # PK IIV (Table 2 %CV column). Reported as exponential-BSV %CV;
    # converted to omega^2 via log(1 + CV^2) so that the on-run %CV back-
    # reconstructs to the reported value: %CV = 100 * sqrt(exp(omega^2) - 1).
    # =========================================================================
    etalka    ~ 0.16332  # Table 2: Ka BSV 42.1% CV (%SE 53.8) -> log(1 + 0.421^2) = 0.16332
    etalcl    ~ 0.16909  # Table 2: CL BSV 42.9% CV (%SE 20.1) -> log(1 + 0.429^2) = 0.16909
    etalvc    ~ 0.03849  # Table 2: Vc BSV 19.8% CV (%SE 22.5) -> log(1 + 0.198^2) = 0.03849
    etalvp    ~ 0.14435  # Table 2: Vp BSV 39.4% CV (%SE 20.8) -> log(1 + 0.394^2) = 0.14435
    etalkint  ~ 0.21785  # Table 2: KINT BSV 49.3% CV (%SE 36.7; typical value FIXED; BSV estimated) -> log(1 + 0.493^2) = 0.21785

    # =========================================================================
    # PK residual error (Table 2). Combined additive + proportional. Paper
    # reports 'Magnitude' as variance; nlmixr2 takes SDs, so values are the
    # square root of the reported variances.
    # =========================================================================
    addSd  <- sqrt(3.17e-04);  label("Additive residual SD on TAK-079 serum concentration (ug/mL)")   # Table 2: RV_add var = 3.17e-04 (%SE 21.2); SD = 0.01780 ug/mL
    propSd <- sqrt(0.0677);    label("Proportional residual SD on TAK-079 serum concentration (fraction)")  # Table 2: RV_prop var = 0.0677 (%SE 1.58); SD = 0.2602 fraction

    # =========================================================================
    # NK-cell PK-PD -- Roepcke 2018 Table 3 (NK Cells block) and the
    # equations in Section 3.4:
    #   d/dt(NK) = K_IN - K_OUT * NK - NK * EFF
    #   EFF     = E_MAX * c / (C50 + c)    with K_OUT = K_IN / BL
    # BL for NK is taken as the paper's median baseline (Section 3.3):
    # 685 cells/uL; individual BL enters as a random effect on lrbase_nk.
    # =========================================================================
    lkin_nk    <- log(13957);   label("NK-cell zero-order production rate KIN (cells / uL / day)")  # Table 3: KIN=13 957 (%SE 4.52)
    lec50_nk   <- log(27.5);    label("NK-cell concentration for half-maximal additional depletion C50 (ug/mL)")  # Table 3: C50=27.5 ug/mL (%SE 20.8)
    lemax_nk   <- log(414.6);   label("NK-cell maximal drug-driven depletion rate EMAX (1/day)")     # Table 3: EMAX=414.6 (%SE 10.4)
    lrbase_nk  <- log(685);     label("NK-cell typical baseline count (cells / uL)")                 # Section 3.3: NK median 685 cells/uL (IQR 482.8-970.1)

    etalkin_nk    ~ 0.80333    # Table 3: NK KIN BSV 111% CV (%SE 21.0) -> log(1 + 1.11^2) = 0.80333
    etalec50_nk   ~ 1.14164    # Table 3: NK C50 BSV 146% CV (%SE 24.5) -> log(1 + 1.46^2) = 1.14164
    etalrbase_nk  ~ 0.07866    # Table 3: NK baseline BSV 28.6% CV (%SE 20.6) -> log(1 + 0.286^2) = 0.07866

    propSd_nkcell <- sqrt(0.2905);  label("NK-cell proportional residual SD on cell count (fraction)")  # Table 3: NK residual var = 0.2905 (%SE 2.49); SD = 0.5390 fraction

    # =========================================================================
    # B-cell PK-PD -- Roepcke 2018 Table 3 (B Cells block) and the five
    # equations in Section 3.4. Four transit compartments (precursor1..
    # precursor4) followed by circulating B cells (bcell). All three rate
    # constants share the same value K_PROL = K_TR = K_CIRC = 4 / MTT. No
    # feedback loop on progenitor cells. Drug effect (E_MAX * c/(C50+c))
    # acts as an extra depletion term on the circulating pool only.
    # =========================================================================
    lmtt_b     <- log(8.19);    label("B-cell mean transit time MTT (day)")                        # Table 3: MTT=8.19 day (%SE 15.3)
    lec50_b    <- log(19.8);    label("B-cell concentration for half-maximal additional depletion C50 (ug/mL)")  # Table 3: C50=19.8 ug/mL (%SE 8.95)
    lemax_b    <- log(2.43);    label("B-cell maximal drug-driven depletion rate EMAX (1/day)")     # Table 3: EMAX=2.43 (%SE 4.96)
    lrbase_b   <- log(1279);    label("B-cell typical baseline count (cells / uL)")                # Section 3.3: B median 1279 cells/uL (IQR 860.8-1890)

    etalmtt_b     ~ 1.03775    # Table 3: B MTT BSV 135% CV (%SE 17.6) -> log(1 + 1.35^2) = 1.03775
    etalrbase_b   ~ 0.05617    # Table 3: B baseline BSV 24.03% CV (%SE 10.7) -> log(1 + 0.2403^2) = 0.05617

    propSd_bcell  <- sqrt(0.136);   label("B-cell proportional residual SD on cell count (fraction)")   # Table 3: B residual var = 0.136 (%SE 2.26); SD = 0.3688 fraction

    # =========================================================================
    # T-cell PK-PD -- Roepcke 2018 Table 3 (T Cells block) and Section 3.4:
    #   T(c) = BL * (1 - EMAX * c / (c + C50))
    # Direct-response, algebraic (no ODE). EMAX is fractional (dimensionless)
    # and bounded above by ~0.47, so only ~half of T-cell counts can be
    # depleted by TAK-079 at any concentration.
    # =========================================================================
    lec50_t    <- log(11.86);   label("T-cell concentration for half-maximal depletion C50 (ug/mL)")  # Table 3: C50=11.86 ug/mL (%SE 7.267)
    lemax_t    <- log(0.4656);  label("T-cell maximal fractional depletion EMAX (dimensionless, <=1)")  # Table 3: EMAX=0.4656 (%SE 6.578)
    lrbase_t   <- log(3732);    label("T-cell typical baseline count (cells / uL)")                # Section 3.3: T median 3732 cells/uL (IQR 2881-5176)

    etalemax_t    ~ 0.39366    # Table 3: T EMAX BSV 69.46% CV (%SE 29.50) -> log(1 + 0.6946^2) = 0.39366
    etalrbase_t   ~ 0.08122    # Table 3: T baseline BSV 29.08% CV (%SE 15.50) -> log(1 + 0.2908^2) = 0.08122

    propSd_tcell  <- sqrt(0.1343);  label("T-cell proportional residual SD on cell count (fraction)")   # Table 3: T residual var = 0.1343 (%SE 2.406); SD = 0.3665 fraction
  })

  model({
    # =========================================================================
    # PK layer -- 2-compartment linear disposition with QSS-TMDD elimination
    # and first-order SC absorption. Equations from Roepcke 2018 Figure 2
    # and the Berkeley Madonna reference implementation in Supplement
    # (lines 155-213 of PMID_29864242_trimmed.md).
    # =========================================================================

    # Bioavailability (logit -> fraction). Only applies to depot dosing (SC).
    fdepot <- exp(logitfdepot) / (1 + exp(logitfdepot))

    # Individual PK parameters. Route effect: SC reduces Vc to ~30% of the IV
    # typical value; IV is the reference (ROUTE_IV = 1).
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) * (1 - e_route_iv_vc * (1 - ROUTE_IV))
    q  <- exp(lq)
    vp <- exp(lvp + etalvp)

    # QSS-TMDD receptor pool (individual only on KINT).
    kint <- exp(lkint + etalkint)
    kss  <- exp(lkss)
    ksyn <- exp(lksyn)
    kdeg <- exp(lkdeg)

    # Micro-rate constants for the linear 2-compartment disposition.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Working concentration in the central compartment (ug/mL).
    conc <- central / vc

    # ODE system -- linear 2-compartment PK with QSS-TMDD elimination on
    # the central drug and a receptor turnover ODE. `total_target` is the
    # receptor concentration Rtot (u/L equivalent) with steady-state
    # Rtot_ss = KSYN / KDEG.
    d/dt(depot)        <- -ka * depot
    d/dt(central)      <-  ka * depot -
                          k12 * central + k21 * peripheral1 -
                          kel * central -
                          kint * total_target * central / (kss + conc)
    d/dt(peripheral1)  <-  k12 * central - k21 * peripheral1
    d/dt(total_target) <-  ksyn - kdeg * total_target -
                          (kint - kdeg) * total_target * conc / (kss + conc)

    # SC bioavailability applied to the depot. IV dosing goes directly to
    # central (F is not applied). No lag time was estimated in the paper.
    f(depot) <- fdepot

    # Receptor steady-state initial condition (Berkeley Madonna
    # `init rtot = RBASE = KSYN/KDEG`).
    total_target(0) <- ksyn / kdeg

    # PK observation.
    Cc <- conc
    Cc ~ add(addSd) + prop(propSd)

    # =========================================================================
    # NK-cell PD layer -- indirect-response / turnover with drug-driven
    # additional depletion. Baseline is individualised through the
    # exponential eta on lrbase_nk; K_OUT is derived from the individual
    # baseline so that at zero drug concentration NK stays at BL.
    # =========================================================================
    kin_nk    <- exp(lkin_nk + etalkin_nk)
    ec50_nk   <- exp(lec50_nk + etalec50_nk)
    emax_nk   <- exp(lemax_nk)
    rbase_nk  <- exp(lrbase_nk + etalrbase_nk)
    kout_nk   <- kin_nk / rbase_nk
    eff_nk    <- emax_nk * conc / (ec50_nk + conc)

    d/dt(nkcell) <- kin_nk - kout_nk * nkcell - nkcell * eff_nk
    nkcell(0)    <- rbase_nk

    nkcell ~ prop(propSd_nkcell)

    # =========================================================================
    # B-cell PD layer -- four-transit-compartment chain with drug-driven
    # additional depletion on the circulating pool. K_PROL = K_TR = K_CIRC =
    # 4 / MTT, so the transit chain and the circulating loss share the same
    # first-order rate. Compartments are named canonically as
    # precursor1..precursor4 (transit chain) and bcell (circulating).
    # No feedback loop from bcell to the precursor pool (Section 3.4:
    # 'No additional effects and no feedback loop on progenitor cells were
    # necessary to describe the available monkey B-cell data').
    # =========================================================================
    mtt_b     <- exp(lmtt_b + etalmtt_b)
    ec50_b    <- exp(lec50_b)
    emax_b    <- exp(lemax_b)
    rbase_b   <- exp(lrbase_b + etalrbase_b)
    ktr_b     <- 4 / mtt_b
    eff_b     <- emax_b * conc / (ec50_b + conc)

    # Chain of maturation compartments plus circulating B pool. TR1 uses a
    # constant-source parameterisation K_PROL * rbase_b so that in absence
    # of drug the steady-state (precursor1 = precursor2 = ... = bcell =
    # rbase_b) is preserved regardless of the compartment's initial
    # condition. TR1..TR4 = paper's TR1..TR4 (Section 3.4); bcell = paper's
    # B (circulating).
    d/dt(precursor1) <- ktr_b * rbase_b - ktr_b * precursor1
    d/dt(precursor2) <- ktr_b * (precursor1 - precursor2)
    d/dt(precursor3) <- ktr_b * (precursor2 - precursor3)
    d/dt(precursor4) <- ktr_b * (precursor3 - precursor4)
    d/dt(bcell)      <- ktr_b * precursor4 - ktr_b * bcell - bcell * eff_b

    precursor1(0) <- rbase_b
    precursor2(0) <- rbase_b
    precursor3(0) <- rbase_b
    precursor4(0) <- rbase_b
    bcell(0)      <- rbase_b

    bcell ~ prop(propSd_bcell)

    # =========================================================================
    # T-cell PD layer -- direct-response algebraic model. No ODE state.
    # Note (Section 3.4): the direct-response model fits the low-dose and
    # later-time-point T-cell data adequately but under-predicts the sharp
    # initial depletion after a first high dose (e.g. the 3 mg/kg group).
    # This is a known limitation carried into the vignette as a caveat.
    # =========================================================================
    ec50_t    <- exp(lec50_t)
    emax_t    <- exp(lemax_t + etalemax_t)
    rbase_t   <- exp(lrbase_t + etalrbase_t)

    tcell <- rbase_t * (1 - emax_t * conc / (ec50_t + conc))
    tcell ~ prop(propSd_tcell)
  })
}
