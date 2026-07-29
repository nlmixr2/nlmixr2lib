Tornoe_2006_degarelix <- function() {
  description <- paste(
    "Population PK/PD model of the hypothalamic-pituitary-gonadal (HPG)",
    "axis after repeated subcutaneous (s.c.) injections of the GnRH",
    "receptor blocker degarelix in prostate-cancer patients. PK is a",
    "two-compartment disposition model with two parallel first-order",
    "absorption routes from a self-forming s.c. depot: a rapid release",
    "(fraction Fr via ka,fast) and a prolonged slow release",
    "((1 - Fr) via ka,slow). The packaged values correspond to the",
    "40 mg/mL dose-concentration arm of Tornoe 2007 Table 3 (Fr_40,",
    "F_40 and t_1/2,slow,40); the 20 and 60 mg/mL alternatives are",
    "tabulated in the validation vignette. PD is a four-state HPG-axis",
    "feedback model (feedback compartment F, LH pool P, LH, testosterone",
    "Te) with sigmoidal Imax inhibition of LH pool release by degarelix",
    "and a positive interaction (F) from the feedback compartment on LH",
    "synthesis and release; testosterone secretion is stimulated by LH",
    "via a sigmoidal Emax model. ke_LH, ke_F, lambda, LH_base and Te_base",
    "are the degarelix-study-specific values from Table 4."
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
    "etalkaslow", "etalogitfr", "etalogitfdeg",
    "etalkrel", "etalic50", "etaldeltapd", "etalkef", "etallmax", "etall50"
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 170L,
    n_studies      = 1L,
    age_range      = "19-89 years (median 73)",
    weight_range   = "45-117 kg (median 78)",
    sex_female_pct = 0,
    disease_state  = "Prostate cancer patients receiving ascending repeated-dose degarelix",
    dose_range     = paste(
      "Initial s.c. loading doses of 120-320 mg in 20-60 mg/mL injection",
      "solutions, followed by maintenance doses identical to the loading",
      "dose when testosterone > 0.5 ng/mL. Eight dose groups: 120 mg @ 20",
      "mg/mL, 120 mg @ 40 mg/mL, 160 mg @ 40 mg/mL, 200 mg @ 40 mg/mL,",
      "200 mg @ 60 mg/mL, 240 mg @ 40 mg/mL, 240 mg @ 60 mg/mL, 320 mg @",
      "60 mg/mL (Table 1)."
    ),
    regions        = "Not reported",
    notes          = paste(
      "Degarelix study: 170 prostate-cancer patients in an open-label,",
      "multicentre, parallel and sequential, ascending repeated-dose",
      "trial. Blood samples were collected pre-dose and at 3 h and 1, 2,",
      "3, 7, 14, 21, 28, 35, 42, 49, 56 days post-dose and once weekly",
      "until insufficient testosterone suppression (testosterone > 0.5",
      "ng/mL) prompted a maintenance dose. Degarelix LLOQ 0.5 ng/mL; LH",
      "LLOQ 0.1 IU/L; testosterone LLOQ 0.05 ng/mL. Patients withdrew at",
      "testosterone > 1.0 ng/mL >= 2 weeks post-dose or 0.5-1.0 ng/mL on",
      "two consecutive >= 4-week visits. The packaged model uses the",
      "40 mg/mL Table 3 values (Fr_40 = 0.0573, F_40 = 0.240, t_1/2,slow,40",
      "= 73.7 d); the 20 and 60 mg/mL parameter sets are documented in",
      "the validation vignette. See Tornoe 2007 Tables 1-2 and Materials",
      "and methods."
    )
  )

  ini({
    # ----------------------------------------------------------------
    # Degarelix PK -- Table 3 (Br J Clin Pharmacol 2007;63:648-664).
    # Two-compartment disposition (CL, Vc, Q, Vp) with two parallel
    # first-order absorption routes from the s.c. depot:
    #   - Rapid release: fraction Fr at ka,fast (t_1/2,fast = 1.98 d)
    #   - Slow release: (1 - Fr) at ka,slow (t_1/2,slow varies with
    #                                       dose concentration)
    # The 40 mg/mL dose-concentration values are packaged as the
    # typical-value reference; the 20 and 60 mg/mL alternatives are
    # listed in the vignette (Tornoe 2007 Table 3 rows Fr_20/Fr_60,
    # F_20/F_60 and the 53.3 / 95.4 d slow-release half-lives).
    # ----------------------------------------------------------------
    lcl       <- log(2.54)                       ; label("Clearance CL (L/h)")                                          # Table 3: 2.54 (RSE 5.43%)
    lvc       <- log(13.2)                       ; label("Central volume Vc (L)")                                       # Table 3: 13.2 (RSE 9.24%)
    lq        <- log(6.59)                       ; label("Inter-compartmental clearance Q (L/h)")                       # Table 3: 6.59 (RSE 7.36%)
    lvp       <- log(36.1)                       ; label("Peripheral volume Vp (L)")                                    # Table 3: 36.1 (RSE 4.99%)

    # ka,fast = ln(2) / (1.98 d * 24 h/d) = 0.01459 1/h
    # ka,slow,40 = ln(2) / (73.7 d * 24 h/d) = 3.919e-4 1/h
    lkafast   <- log(0.01459)                    ; label("Fast absorption rate ka,fast (1/h)")                          # Table 3: t_1/2,fast = 1.98 d (RSE 6.17%)
    lkaslow   <- log(3.919e-4)                   ; label("Slow absorption rate ka,slow at 40 mg/mL (1/h)")              # Table 3: t_1/2,slow,40 = 73.7 d (RSE 4.74%)

    # Fr and F are reported on the natural (0, 1) scale with IIVs
    # given as logit-transformed approximations (CV(q) = (1-q) * w_q).
    # logit(0.0573) = log(0.0573 / (1 - 0.0573)) = -2.798.
    # logit(0.240)  = log(0.240  / (1 - 0.240))  = -1.153.
    logitfr   <- log(0.0573 / (1 - 0.0573))      ; label("Logit-transformed rapid-release fraction Fr at 40 mg/mL (logit units)") # Table 3: Fr_40 = 0.0573 (RSE 6.30%)
    logitfdeg <- log(0.240  / (1 - 0.240))       ; label("Logit-transformed bioavailability F at 40 mg/mL (logit units)")         # Table 3: F_40  = 0.240  (RSE 6.83%)

    # ----------------------------------------------------------------
    # HPG-axis PD -- Table 4. Degarelix-study-specific values are
    # taken from the rows marked dagger (sup-script #). Sigmoidal Imax
    # inhibition of LH pool release by degarelix (H2: Imax, IC50,
    # delta). Pulsatile LH synthesis and release is *positively*
    # modulated by the feedback compartment for the GnRH receptor
    # blocker via H5(F) = H6(F) = F.
    # ----------------------------------------------------------------
    lkelh     <- log(0.535)                      ; label("LH elimination rate ke,LH (1/h; degarelix-specific)")         # Table 4: ke,LH^# = 0.535 (RSE 4.14%)
    lkrel     <- log(0.00241)                    ; label("LH pool release rate krel,LH (1/h)")                          # Table 4: krel,LH = 0.00241 (RSE 6.36%)
    limax     <- log(0.942)                      ; label("Maximal degarelix inhibition of LH pool release Imax (unitless)") # Table 4: Imax = 0.942 (RSE 0.155%)
    lic50     <- log(1.49)                       ; label("Potency of degarelix inhibition IC50 (ng/mL)")                # Table 4: IC50 = 1.49 (RSE 5.04%)
    ldeltapd  <- log(1.97)                       ; label("Sigmoidicity of degarelix inhibition delta (unitless)")       # Table 4: delta = 1.97 (RSE 3.45%)

    lkef      <- log(0.00497)                    ; label("Feedback elimination ke,F (1/h; degarelix-specific)")          # Table 4: ke,F^# = 0.00497 (RSE 4.87%)
    llambda   <- log(0.56)                       ; label("Feedback-stimulation exponent lambda (unitless; degarelix-specific)") # Table 4: lambda^# = 0.56 (RSE 1.00%)

    lketet    <- log(0.0901)                     ; label("Testosterone elimination rate ke,Te (1/h)")                   # Table 4: ke,Te = 0.0901 (RSE 2.72%)
    llmax     <- log(77.5)                       ; label("Maximal LH stimulation of Te secretion Lmax (unitless)")      # Table 4: Lmax = 77.5 (RSE 3.51%)
    ll50      <- log(5.18)                       ; label("LH producing half-maximal Te stimulation L50 (IU/L)")         # Table 4: L50 = 5.18 (RSE 2.98%)
    lkappa    <- log(1.9)                        ; label("Sigmoidicity of LH-to-Te stimulation kappa (unitless)")       # Table 4: kappa = 1.9 (RSE 0.836%)

    lhbase    <- log(6.98)                       ; label("Baseline LH LH_base (IU/L; degarelix-specific)")              # Table 4: LH_base^# = 6.98 (RSE 1.34%)
    ltbase    <- log(3.21)                       ; label("Baseline testosterone Te_base (ng/mL; degarelix-specific)")   # Table 4: Te_base^# = 3.21 (RSE 1.47%)

    # ----------------------------------------------------------------
    # Inter-individual variability -- omega^2 = log(CV^2 + 1) for
    # exponential / log-normal IIV. Logit-scale variances for Fr and
    # F are derived from the paper-approximation CV(q) = (1-q) * w_q.
    # ----------------------------------------------------------------
    etalcl       ~ log(1 + 0.281^2)              # Table 3 IIV CL    CV 28.1% (RSE 9.61)
    etalvc       ~ log(1 + 0.246^2)              # Table 3 IIV Vc    CV 24.6% (RSE 32.3)
    etalkaslow   ~ log(1 + 0.444^2)              # Table 3 IIV t_1/2,slow,40 CV 44.4% (RSE 7.30); same on log-rate
    etalogitfr   ~ (0.378 / (1 - 0.0573))^2      # Table 3 IIV Fr_40 CV 37.8% (RSE 10.5) on logit scale
    etalogitfdeg ~ (0.232 / (1 - 0.240))^2       # Table 3 IIV F_40  CV 23.2% (RSE 13.4) on logit scale

    etalkrel     ~ log(1 + 0.834^2)              # Table 4 IIV krel,LH CV 83.4% (RSE 3.25)
    etalic50     ~ log(1 + 0.595^2)              # Table 4 IIV IC50    CV 59.5% (RSE 6.44)
    etaldeltapd  ~ log(1 + 0.389^2)              # Table 4 IIV delta   CV 38.9% (RSE 6.73)
    etalkef      ~ log(1 + 0.591^2)              # Table 4 IIV ke,F    CV 59.1% (RSE 4.69)
    etallmax     ~ log(1 + 0.517^2)              # Table 4 IIV Lmax    CV 51.7% (RSE 6.24)
    etall50      ~ log(1 + 0.460^2)              # Table 4 IIV L50     CV 46.0% (RSE 5.42)

    # ----------------------------------------------------------------
    # Residual error: proportional ("constant CV") on the untransformed
    # scale for each of degarelix Cc, LH and Te (Methods, Residual
    # variability).
    # ----------------------------------------------------------------
    propSd    <- 0.287 ; label("Proportional residual error for degarelix Cc (fraction)") # Table 3: s_prop = 28.7% (RSE 2.53)
    propSd_LH <- 0.419 ; label("Proportional residual error for LH (fraction)")           # Table 4: s_LH   = 41.9% (RSE 1.00)
    propSd_Te <- 0.494 ; label("Proportional residual error for testosterone (fraction)") # Table 4: s_Te   = 49.4% (RSE 0.983)
  })

  model({
    # ----------------------------------------------------------------
    # 1. Individual structural PK parameters
    # ----------------------------------------------------------------
    cl     <- exp(lcl + etalcl)
    vc     <- exp(lvc + etalvc)
    q      <- exp(lq)
    vp     <- exp(lvp)
    kafast <- exp(lkafast)
    kaslow <- exp(lkaslow + etalkaslow)
    fr     <- expit(logitfr   + etalogitfr)
    fdeg   <- expit(logitfdeg + etalogitfdeg)

    # 2-compartment disposition micro-constants
    kel    <- cl / vc
    k12    <- q  / vc
    k21    <- q  / vp

    # ----------------------------------------------------------------
    # 2. Individual PD parameters (degarelix-study-specific values
    #    plus shared structural PD parameters)
    # ----------------------------------------------------------------
    kelh    <- exp(lkelh)
    krel    <- exp(lkrel + etalkrel)
    imax    <- exp(limax)
    ic50    <- exp(lic50 + etalic50)
    deltapd <- exp(ldeltapd + etaldeltapd)
    kef     <- exp(lkef  + etalkef)
    lambda  <- exp(llambda)
    ketet   <- exp(lketet)
    lmax    <- exp(llmax + etallmax)
    l50     <- exp(ll50  + etall50)
    kappa   <- exp(lkappa)
    lhbasei <- exp(lhbase)
    tbasei  <- exp(ltbase)

    # Basal secretion rates from steady-state balance (Tornoe 2007
    # Results section equations):
    #   beta_F  = ke,F
    #   beta_LH = ke,LH * LH_base
    #   beta_Te = ke,Te * Te_base
    #             / (1 + Lmax * LH_base^kappa
    #                   / (L50^kappa + LH_base^kappa))
    bf    <- kef
    blh   <- kelh * lhbasei
    h3_ss <- lmax * lhbasei^kappa / (l50^kappa + lhbasei^kappa)
    bte   <- ketet * tbasei / (1 + h3_ss)

    # ----------------------------------------------------------------
    # 3. PK ODE system. Two parallel first-order absorption routes
    #    from the s.c. depot:
    #    - depot  (rapid): receives bolus Dose with f(depot)  = Fr * F
    #    - depot2 (slow):  receives bolus Dose with f(depot2) = (1 - Fr) * F
    # ----------------------------------------------------------------
    d/dt(depot)       <- -kafast * depot
    d/dt(depot2)      <- -kaslow * depot2
    d/dt(central)     <-  kafast * depot + kaslow * depot2 -
                          kel * central -
                          k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot)  <- fr   * fdeg
    f(depot2) <- (1 - fr) * fdeg

    # ----------------------------------------------------------------
    # 4. HPG-axis PD ODE system (Tornoe 2007 Results, system of ODEs).
    #    For the GnRH receptor blocker degarelix, H5(F) = H6(F) = F.
    # ----------------------------------------------------------------
    # central (mg) / vc (L) gives mg/L = 1000 ng/mL; multiply by 1000 to
    # convert to ng/mL (the paper's reporting unit for degarelix and
    # the unit of IC50, Te_base).
    cp_deg <- central / vc * 1000                 # degarelix Cp (ng/mL)
    h2     <- -imax * cp_deg^deltapd / (ic50^deltapd + cp_deg^deltapd)
    h3     <- lmax * lh^kappa / (l50^kappa + lh^kappa)
    h4     <- (testosterone / tbasei)^lambda
    h5_f   <- feedback
    h6_f   <- feedback

    d/dt(feedback)     <-  bf  * h4 - kef * feedback
    d/dt(lhpool)       <-  blh * h5_f -
                           krel * lhpool * (1 + h2) * h6_f
    d/dt(lh)           <-  krel * lhpool * (1 + h2) * h6_f - kelh * lh
    d/dt(testosterone) <-  bte * (1 + h3) - ketet * testosterone

    # Steady-state initial conditions (Tornoe 2007 Results):
    #   F0 = 1; P0 = beta_LH / krel,LH; LH0 = LH_base; Te0 = Te_base.
    feedback(0)     <- 1
    lhpool(0)       <- blh / krel
    lh(0)           <- lhbasei
    testosterone(0) <- tbasei

    # ----------------------------------------------------------------
    # 5. Observations and residual error.
    # ----------------------------------------------------------------
    Cc <- cp_deg
    LH <- lh
    Te <- testosterone

    Cc ~ prop(propSd)
    LH ~ prop(propSd_LH)
    Te ~ prop(propSd_Te)
  })
}
