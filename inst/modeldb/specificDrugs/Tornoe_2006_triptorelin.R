Tornoe_2006_triptorelin <- function() {
  description <- paste(
    "Population PK/PD model of the hypothalamic-pituitary-gonadal (HPG)",
    "axis after a single 3.75 mg subcutaneous (s.c.) depot injection of",
    "the GnRH agonist triptorelin in healthy adult males. PK is a",
    "two-compartment disposition model with a combined zero-order burst",
    "(fraction Fr of dose released over duration t into central) and a",
    "two-step first-order s.c. absorption (lymphatic delay) for the",
    "remaining (1 - Fr) fraction. PD is a four-state HPG-axis feedback",
    "model (feedback compartment F, LH pool P, LH, testosterone Te) with",
    "sigmoidal Emax stimulation of LH pool release by triptorelin and a",
    "negative interaction (F^-1) from the feedback compartment on LH",
    "synthesis and release; testosterone secretion is stimulated by LH",
    "via a sigmoidal Emax model. ke_LH, ke_F, lambda, LH_base, and",
    "Te_base are the triptorelin-study-specific values from Table 4."
  )
  reference <- paste(
    "Tornoe CW, Agerso H, Senderovitz T, Nielsen HA, Madsen H,",
    "Karlsson MO, Jonsson EN.",
    "Population pharmacokinetic/pharmacodynamic (PK/PD) modelling of the",
    "hypothalamic-pituitary-gonadal axis following treatment with GnRH",
    "analogues.",
    "Br J Clin Pharmacol 2007;63(6):648-664.",
    "doi:10.1111/j.1365-2125.2006.02820.x.",
    sep = " "
  )
  vignette <- "Tornoe_2006_HPG_axis"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  paper_specific_compartments <- c("feedback", "lhpool", "lh", "testosterone")
  paper_specific_etas <- c(
    "etalksc1", "etalksc2", "etalogitfr",
    "etalkrel", "etalec50", "etalkef", "etallmax", "etall50"
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 30L,
    n_studies      = 1L,
    age_range      = "20-74 years (median 41)",
    weight_range   = "60-111 kg (median 80)",
    sex_female_pct = 0,
    disease_state  = "Healthy adult males",
    dose_range     = "Single subcutaneous (s.c.) injection of 3.75 mg Decapeptyl Depot",
    regions        = "Not reported",
    notes          = paste(
      "Triptorelin sub-study: 30 healthy males (s.c. arm) of a 58-subject",
      "single-dose, unblinded, randomized, parallel-group study",
      "investigating PK/PD/tolerability of Decapeptyl Depot after single",
      "s.c. or i.m. injections of 3.75 mg. Blood samples were collected",
      "pre-dose and at 0.25, 1, 2, 3, 4, 8, 12, 24 h and 2, 4, 7, 14, 21,",
      "28, 35, 42, 49, 56, 70 days post-dose with additional samples on",
      "days 84, 98 and at 2-week intervals until testosterone returned to",
      "the physiological range. Triptorelin LLOQ 0.01 ng/mL; LH LLOQ 0.1",
      "IU/L; testosterone LLOQ 0.05 ng/mL. The packaged model covers the",
      "s.c. arm only; the i.m. arm (28 subjects, t_1/2,im = 17.0 days,",
      "first-order absorption with no lymphatic delay) is not implemented",
      "in this file. See Tornoe 2007 Tables 1-2 and Materials and methods."
    )
  )

  ini({
    # ----------------------------------------------------------------
    # Triptorelin PK -- Table 3 (Br J Clin Pharmacol 2007;63:648-664).
    # Two-compartment disposition (CL/F, Vc/F, Q/F, Vp/F) with a
    # combined zero-order burst (fraction Fr over duration trel into
    # central) and a two-step first-order s.c. absorption pathway
    # (depot -> transit1 -> central) with rates ksc,1 and ksc,2 derived
    # from the published absorption half-lives.
    # ----------------------------------------------------------------
    lcl    <- log(63.2)    ; label("Apparent clearance CL/F (L/h)")                                # Table 3: 63.2 (RSE 4.13%)
    lvc    <- log(640)     ; label("Apparent central volume Vc/F (L)")                             # Table 3: 640  (RSE 5.77%)
    lq     <- log(76.3)    ; label("Apparent inter-compartmental clearance Q/F (L/h)")             # Table 3: 76.3 (RSE 10.5%)
    lvp    <- log(698)     ; label("Apparent peripheral volume Vp/F (L)")                          # Table 3: 698  (RSE 7.95%)

    # SC absorption: ksc,1 = ln(2) / t_1/2,sc,1, ksc,2 = ln(2) / t_1/2,sc,2.
    # Table 3 t_1/2,sc,1 = 11.3 d -> ksc,1 = 0.002556 1/h
    # Table 3 t_1/2,sc,2 = 7.92 d -> ksc,2 = 0.003647 1/h
    lksc1  <- log(0.002556); label("SC absorption rate ksc,1 (1/h)")                                # Table 3: t_1/2,sc,1 = 11.3 d (RSE 9.78%)
    lksc2  <- log(0.003647); label("SC absorption rate ksc,2 (1/h)")                                # Table 3: t_1/2,sc,2 = 7.92 d (RSE 17.9%)

    # Burst input: fraction Fr of dose released as zero-order infusion
    # of duration trel into central. Fr is logit-transformed per the
    # Table 3 IIV footnote (CV(q) = (1-q) * w_q).
    # logit(0.605) = log(0.605 / (1 - 0.605)) = 0.4263.
    logitfr <- 0.4263      ; label("Logit-transformed burst fraction Fr (logit units)")             # Table 3: Fr = 0.605 (RSE 2.08%)
    ltrel   <- log(1.77)   ; label("Duration of zero-order burst infusion t (h)")                    # Table 3: t = 1.77 h (RSE 4.45%)

    base_trip <- fixed(0.0107) ; label("Triptorelin baseline plasma concentration (ng/mL)")          # Table 3: Base = 0.0107 (RSE 5.02%); used as additive offset on Cc

    # ----------------------------------------------------------------
    # HPG-axis PD -- Table 4. Triptorelin-study-specific values are
    # taken from the rows marked * (asterisk). Sigmoidal Emax
    # stimulation of LH pool release by triptorelin (H1: Emax, EC50,
    # gamma). LH stimulates testosterone secretion via H3 (Lmax, L50,
    # kappa). Testosterone stimulates the feedback compartment via H4
    # ((Te/Tebase)^lambda). Pulsatile LH synthesis and release is
    # *negatively* modulated by the feedback compartment for the GnRH
    # agonist via H5(F) = H6(F) = F^-1.
    # ----------------------------------------------------------------
    lkelh   <- log(0.0082) ; label("LH elimination rate constant ke,LH (1/h; triptorelin-specific)") # Table 4: ke,LH* = 0.0082 (RSE 1.56%)
    lkrel   <- log(0.00241); label("LH pool release rate constant krel,LH (1/h)")                    # Table 4: krel,LH = 0.00241 (RSE 6.36%)
    lemax   <- log(1330)   ; label("Maximal triptorelin stimulation of LH pool release Emax (unitless)") # Table 4: Emax = 1330 (RSE 8.58%)
    lec50   <- log(0.047)  ; label("Potency of triptorelin stimulation of LH pool release EC50 (ng/mL)") # Table 4: EC50 = 0.047 (RSE 5.53%)
    lgamma  <- log(4.87)   ; label("Sigmoidicity of triptorelin stimulation gamma (unitless)")       # Table 4: gamma = 4.87 (RSE 3.69%)

    lkef    <- log(0.00107); label("Feedback compartment elimination ke,F (1/h; triptorelin-specific)") # Table 4: ke,F* = 0.00107 (RSE 8.26%)
    llambda <- log(8.26)   ; label("Feedback-stimulation exponent lambda (unitless; triptorelin-specific)") # Table 4: lambda* = 8.26 (RSE 3.65%)

    lketet  <- log(0.0901) ; label("Testosterone elimination rate constant ke,Te (1/h)")             # Table 4: ke,Te = 0.0901 (RSE 2.72%)
    llmax   <- log(77.5)   ; label("Maximal LH stimulation of testosterone secretion Lmax (unitless)") # Table 4: Lmax = 77.5 (RSE 3.51%)
    ll50    <- log(5.18)   ; label("LH concentration producing half-maximal testosterone stimulation L50 (IU/L)") # Table 4: L50 = 5.18 (RSE 2.98%)
    lkappa  <- log(1.9)    ; label("Sigmoidicity of LH-to-testosterone stimulation kappa (unitless)") # Table 4: kappa = 1.9 (RSE 0.836%)

    lhbase  <- log(4.76)   ; label("Baseline LH concentration LH_base (IU/L; triptorelin-specific)") # Table 4: LH_base* = 4.76 (RSE 1.80%)
    ltbase  <- log(4.85)   ; label("Baseline testosterone concentration Te_base (ng/mL; triptorelin-specific)") # Table 4: Te_base* = 4.85 (RSE 1.64%)

    # ----------------------------------------------------------------
    # Inter-individual variability -- omega^2 = log(CV^2 + 1) for
    # exponential / log-normal IIV. For logit-transformed Fr the
    # paper's CV is converted to an underlying logit-scale variance via
    # w_q = CV(q) / (1 - q); variance is w_q^2.
    # ----------------------------------------------------------------
    etalcl     ~ log(1 + 0.216^2)  # Table 3 IIV CL/F CV 21.6% (RSE 10.4)
    etalvc     ~ log(1 + 0.329^2)  # Table 3 IIV Vc/F CV 32.9% (RSE 14.1)
    etalksc1   ~ log(1 + 0.340^2)  # Table 3 IIV t_1/2,sc,1 CV 34.0% (RSE 21.3); same on log-rate
    etalksc2   ~ log(1 + 0.616^2)  # Table 3 IIV t_1/2,sc,2 CV 61.6% (RSE 25.2); same on log-rate
    etalogitfr ~ (0.0660 / (1 - 0.605))^2  # Table 3 IIV Fr CV 6.60% (RSE 27.6) on logit scale

    etalkrel   ~ log(1 + 0.834^2)  # Table 4 IIV krel,LH CV 83.4% (RSE 3.25)
    etalec50   ~ log(1 + 0.327^2)  # Table 4 IIV EC50    CV 32.7% (RSE 13.4)
    etalkef    ~ log(1 + 0.591^2)  # Table 4 IIV ke,F    CV 59.1% (RSE 4.69)
    etallmax   ~ log(1 + 0.517^2)  # Table 4 IIV Lmax    CV 51.7% (RSE 6.24)
    etall50    ~ log(1 + 0.460^2)  # Table 4 IIV L50     CV 46.0% (RSE 5.42)

    # ----------------------------------------------------------------
    # Residual error: proportional ("constant CV") on the untransformed
    # scale for each of triptorelin Cc, LH and Te (Methods, Residual
    # variability).
    # ----------------------------------------------------------------
    propSd    <- 0.278 ; label("Proportional residual error for triptorelin Cc (fraction)") # Table 3: s_prop = 27.8% (RSE 4.78)
    propSd_LH <- 0.419 ; label("Proportional residual error for LH (fraction)")             # Table 4: s_LH  = 41.9% (RSE 1.00)
    propSd_Te <- 0.494 ; label("Proportional residual error for testosterone (fraction)")   # Table 4: s_Te  = 49.4% (RSE 0.983)
  })

  model({
    # ----------------------------------------------------------------
    # 1. Individual structural PK parameters
    # ----------------------------------------------------------------
    cl    <- exp(lcl + etalcl)
    vc    <- exp(lvc + etalvc)
    q     <- exp(lq)
    vp    <- exp(lvp)
    ksc1  <- exp(lksc1 + etalksc1)
    ksc2  <- exp(lksc2 + etalksc2)
    fr    <- expit(logitfr + etalogitfr)
    trel  <- exp(ltrel)

    # 2-compartment disposition micro-constants
    kel   <- cl / vc
    k12   <- q  / vc
    k21   <- q  / vp

    # ----------------------------------------------------------------
    # 2. Individual PD parameters (triptorelin-study-specific values
    #    plus shared structural PD parameters)
    # ----------------------------------------------------------------
    kelh    <- exp(lkelh)
    krel    <- exp(lkrel + etalkrel)
    emax    <- exp(lemax)
    ec50    <- exp(lec50 + etalec50)
    gamma   <- exp(lgamma)
    kef     <- exp(lkef  + etalkef)
    lambda  <- exp(llambda)
    ketet   <- exp(lketet)
    lmax    <- exp(llmax + etallmax)
    l50     <- exp(ll50  + etall50)
    kappa   <- exp(lkappa)
    lhbasei <- exp(lhbase)
    tbasei  <- exp(ltbase)

    # Basal secretion rates from steady-state balance (Tornoe 2007
    # Results section equations for beta_F, beta_LH, beta_Te):
    #   beta_F  = ke,F                            (since F0 = 1)
    #   beta_LH = ke,LH * LH_base                  (since H6(F=1) = 1)
    #   beta_Te = ke,Te * Te_base
    #             / (1 + Lmax * LH_base^kappa
    #                   / (L50^kappa + LH_base^kappa))
    bf    <- kef
    blh   <- kelh * lhbasei
    h3_ss <- lmax * lhbasei^kappa / (l50^kappa + lhbasei^kappa)
    bte   <- ketet * tbasei / (1 + h3_ss)

    # ----------------------------------------------------------------
    # 3. PK ODE system. Dose-event convention:
    #    - bolus to depot of amount Dose (with f(depot) = 1 - fr giving
    #      (1-Fr)*Dose entering the s.c. slow-release pathway), and
    #    - infusion to central of amount Dose (with f(central) = fr and
    #      dur(central) = trel giving a zero-order Fr*Dose burst over
    #      trel hours).
    # ----------------------------------------------------------------
    d/dt(depot)       <- -ksc1 * depot
    d/dt(transit1)    <-  ksc1 * depot - ksc2 * transit1
    d/dt(central)     <-  ksc2 * transit1 - kel * central -
                          k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(central)   <- fr
    dur(central) <- trel
    f(depot)     <- 1 - fr

    # ----------------------------------------------------------------
    # 4. HPG-axis PD ODE system (Tornoe 2007 Results, system of ODEs
    #    for dF/dt, dP/dt, dLH/dt, dTe/dt). For the GnRH agonist
    #    triptorelin, H5(F) = H6(F) = F^-1.
    # ----------------------------------------------------------------
    # central (mg) / vc (L) gives mg/L = 1000 ng/mL; multiply by 1000 to
    # convert to ng/mL (the paper's reporting unit for triptorelin and
    # the unit of EC50, Te_base, base_trip).
    cp_trip <- central / vc * 1000                # triptorelin Cp (ng/mL)
    h1      <- emax * cp_trip^gamma / (ec50^gamma + cp_trip^gamma)
    h3      <- lmax * lh^kappa / (l50^kappa + lh^kappa)
    h4      <- (testosterone / tbasei)^lambda
    h5_f    <- 1 / feedback
    h6_f    <- 1 / feedback

    d/dt(feedback)     <-  bf  * h4 - kef * feedback
    d/dt(lhpool)       <-  blh * h5_f -
                           krel * lhpool * (1 + h1) * h6_f
    d/dt(lh)           <-  krel * lhpool * (1 + h1) * h6_f - kelh * lh
    d/dt(testosterone) <-  bte * (1 + h3) - ketet * testosterone

    # Steady-state initial conditions (Tornoe 2007 Results):
    #   F0 = 1; P0 = beta_LH / krel,LH; LH0 = LH_base; Te0 = Te_base.
    feedback(0)     <- 1
    lhpool(0)       <- blh / krel
    lh(0)           <- lhbasei
    testosterone(0) <- tbasei

    # ----------------------------------------------------------------
    # 5. Observations and residual error.
    # Cc = triptorelin plasma concentration with additive baseline
    # offset Base; LH and Te are paper-named outputs.
    # ----------------------------------------------------------------
    Cc <- cp_trip + base_trip
    LH <- lh
    Te <- testosterone

    Cc ~ prop(propSd)
    LH ~ prop(propSd_LH)
    Te ~ prop(propSd_Te)
  })
}
