Savic_2010_warfarin <- function() {
  description <- paste(
    "Population PKPD model for orally dosed warfarin in adult subjects,",
    "presented as the worked illustration of MONOLIX 3.1's SAEM",
    "algorithm for ordered-categorical PD data. PK: one-compartment",
    "with first-order absorption and a lag time. PD link: effect",
    "compartment driven by central amount via rate constant ke0. PD",
    "endpoint: a three-category recoding of percent prothrombin complex",
    "activity (PCA) with cutoffs 50% and 33% (Y=0 if PCA > 50%, Y=1 if",
    "33% <= PCA <= 50%, Y=2 if PCA < 33%), described by a proportional-",
    "odds (cumulative-logit) model with random intercept driven by",
    "effect-site warfarin concentration. The PD categorisation is",
    "acknowledged by the authors (Page 6) as 'done for illustration",
    "purpose only ... not recommended in the real analysis'; this",
    "extraction is the registry's founding example of an ordered-",
    "categorical PD likelihood and the authors' caveat applies. All",
    "parameter values are from the MONOLIX output in Fig. 4."
  )
  reference <- paste(
    "Savic RM, Mentre F, Lavielle M. (2011).",
    "Implementation and Evaluation of the SAEM Algorithm for",
    "Longitudinal Ordered Categorical Data with an Illustration in",
    "Pharmacokinetics-Pharmacodynamics.",
    "The AAPS Journal 13(1):44-53.",
    "doi:10.1208/s12248-010-9238-5.",
    sep = " "
  )
  vignette <- "Savic_2010_warfarin"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/L"
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 33L,
    n_studies      = 1L,
    age_range      = "(not reported in Savic 2010)",
    weight_range   = "(not reported in Savic 2010; parameters reported in absolute units, e.g., V = 7.96 L)",
    sex_female_pct = NA_real_,
    disease_state  = "adult subjects after a single oral dose of warfarin (Savic 2010 Page 6 wording: '33 patients'; underlying O'Reilly 1968 / 1963 cohort was healthy adult volunteers)",
    dose_range     = "single oral dose of warfarin (Savic 2010 reports parameter point estimates in absolute units; in the original O'Reilly 1968 study the standard dose was 1.5 mg/kg orally as a single dose)",
    regions        = "United States (San Francisco VA Medical Center; underlying O'Reilly studies)",
    notes          = paste(
      "33 subjects from the historical O'Reilly & Aggeler (1968)",
      "Circulation 38:169 and O'Reilly et al. (1963) J Clin Invest",
      "42:1542 warfarin PKPD dataset, re-used by Savic 2010 as the",
      "canonical PKPD illustration dataset for the SAEM-on-",
      "categorical-data demonstration (Page 6: '251 PK observations",
      "and 232 PD observations ... for 140 h post dose'). PD endpoint",
      "in Savic 2010 is a three-category recoding of the originally",
      "continuous PCA (prothrombin complex activity, 0-100%): Y = 0",
      "if PCA > 50%, Y = 1 if 33% <= PCA <= 50%, Y = 2 if PCA < 33%.",
      "The authors note (Page 6) that this categorisation is 'done for",
      "illustration purpose only, and it is not recommended to be done",
      "in the real analysis'; the cutoffs were chosen to approximate",
      "INR clinical-target ranges (low INR <2 ~ category 0; targeted",
      "INR 2-3 ~ category 1; high INR >3 ~ category 2). See vignette",
      "Assumptions and deviations for the model-scope rationale."
    )
  )

  ini({
    # ============================================================
    # PK parameters -- Savic 2010 Fig. 4 (MONOLIX output, page 8 of
    # PDF). All values are MONOLIX SAEM final estimates; RSE in
    # parentheses on each line trails the value.
    # ============================================================
    ltlag <- log(0.9)     ; label("Absorption lag time (h)")               # Fig. 4: Tlag = 0.9 h (RSE 21%)
    lka   <- log(1.45)    ; label("Absorption rate constant (1/h)")        # Fig. 4: ka = 1.45 1/h (RSE 37%)
    lvc   <- log(7.96)    ; label("Apparent central volume V (L)")         # Fig. 4: V = 7.96 L (RSE 4%)
    lcl   <- log(0.132)   ; label("Apparent clearance Cl (L/h)")           # Fig. 4: Cl = 0.132 L/h (RSE 5%)
    lke0  <- log(0.0179)  ; label("Effect-compartment equilibration rate ke0 (1/h)")  # Fig. 4: ke0 = 0.0179 1/h (RSE 6%)

    # ============================================================
    # PD parameters -- proportional-odds (cumulative-logit) model
    # with random intercept (Savic 2010 Fig. 3 MONOLIX
    # $CATEGORICAL block):
    #   LOGIT1(Y >= 2) = alpha1 + beta * Ce
    #   LOGIT1(Y >= 1) = alpha1 + alpha2 + beta * Ce
    # alpha1 is the cumulative-logit intercept for P(Y >= 2);
    # alpha2 is the additive increment that makes the intercept
    # for P(Y >= 1) less negative (alpha2 > 0 ensures P(Y >= 1)
    # >= P(Y >= 2) at every Ce); beta is the PD slope of the
    # cumulative logit on effect-site concentration. alpha1 is
    # estimated on the linear scale and is negative; alpha2 and
    # beta are estimated on the linear scale and are positive.
    # They are kept on the linear scale here because the MONOLIX
    # estimation used a linear-scale random intercept on alpha1
    # (omega2_alpha1 = 8.74; linear-scale variance).
    # ============================================================
    alpha1 <- -10.5       ; label("Cumulative-logit intercept for P(Y >= 2) (unitless)")                                       # Fig. 4: alpha1 = -10.5 (RSE 15%)
    alpha2 <- 5.41        ; label("Increment from intercept for P(Y >= 2) to intercept for P(Y >= 1) (unitless, positive)")    # Fig. 4: alpha2 = 5.41 (RSE 17%)
    beta   <- 4.5         ; label("PD slope of cumulative logit on Ce (1/(mg/L))")                                              # Fig. 4: beta = 4.5 (RSE 12%)

    # ============================================================
    # IIV. Savic 2010 Fig. 4 reports omega^2 for Tlag, ka, V, Cl,
    # ke0 (log-scale variances since the PK parameters were
    # estimated as log-normal), and for alpha1 (linear-scale
    # variance since alpha1 was estimated on the linear scale).
    # ============================================================
    etaltlag  ~ 0.252     # Fig. 4: omega2_Tlag = 0.252 (RSE 59%)
    etalka    ~ 0.689     # Fig. 4: omega2_ka = 0.689 (RSE 66%)
    etalvc    ~ 0.0478    # Fig. 4: omega2_V = 0.0478 (RSE 28%)
    etalcl    ~ 0.0797    # Fig. 4: omega2_Cl = 0.0797 (RSE 26%)
    etalke0   ~ 0.0229    # Fig. 4: omega2_ke0 = 0.0229 (RSE 87%)
    etaalpha1 ~ 8.74      # Fig. 4: omega2_alpha1 = 8.74 (RSE 49%); linear-scale eta on alpha1 (cumulative-logit intercept)

    # ============================================================
    # PK residual error. Savic 2010 Fig. 4 reports both an
    # additive (a = 0.231) and a proportional (b = 0.0632) term
    # for the PK observation. The PD observation is categorical
    # and has no residual-error parameter in Fig. 4.
    # ============================================================
    addSd  <- 0.231       ; label("Additive PK residual SD (mg/L)")          # Fig. 4: a = 0.231 (RSE 20%)
    propSd <- 0.0632      ; label("Proportional PK residual SD (fraction)")  # Fig. 4: b = 0.0632 (RSE 15%)
  })

  model({
    # ============================================================
    # 1. Individual structural parameters. Log-normal IIV on each
    # PK parameter (Tlag, ka, V, Cl, ke0). Linear-scale random
    # intercept on alpha1 (Savic 2010 used a linear-scale eta in
    # MONOLIX; see ini() comment on etaalpha1). alpha2 and beta
    # carry no IIV in Fig. 4.
    # ============================================================
    tlag <- exp(ltlag + etaltlag)
    ka   <- exp(lka   + etalka)
    vc   <- exp(lvc   + etalvc)
    cl   <- exp(lcl   + etalcl)
    ke0  <- exp(lke0  + etalke0)

    alpha1_ind <- alpha1 + etaalpha1

    # ============================================================
    # 2. Micro-constants and ODEs. One-compartment oral PK with a
    # depot, first-order absorption ka, and elimination rate
    # constant kel = cl/vc; absorption lag tlag applied to the
    # depot. The effect compartment is driven by the central
    # amount via ke0 (Savic 2010 Fig. 3 $ODE block: DDT_Qc = -k*Qc;
    # DDT_Qe = ke0*Qc - ke0*Qe). Both `central` and `effect` are
    # amount states; Cc and Ce are obtained by dividing by vc.
    # The effect compartment is mass-balanced rather than the
    # Sheiner-style notional-volume formulation because Savic
    # 2010's MONOLIX code drives the effect-compartment ODE
    # with the central AMOUNT (Qc), not the central concentration
    # (Cc); after dividing through by vc this reduces to the
    # standard d(Ce)/dt = ke0 * (Cc - Ce) form.
    # ============================================================
    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central
    d/dt(effect)  <-  ke0 * central - ke0 * effect

    alag(depot) <- tlag

    # ============================================================
    # 3. Observation quantities. Cc is plasma concentration (mg/L)
    # and is the PK observation. Ce is the effect-site
    # concentration (mg/L; same vc as central per Savic 2010
    # Fig. 3) that drives the PD logit. P_ge1 and P_ge2 are the
    # cumulative-probability outputs P(Y >= 1) and P(Y >= 2);
    # P_eq0, P_eq1, P_eq2 are the discrete category probabilities.
    # These are exposed as derived states so downstream
    # simulations can either plot the typical-value probability
    # trajectory directly (replicating Savic 2010 Fig. 5) or
    # sample a categorical observation post-hoc using
    # category index k with probability P_eq_k.
    # ============================================================
    Cc <- central / vc
    Ce <- effect  / vc

    P_ge2 <- 1.0 / (1.0 + exp(-(alpha1_ind + beta * Ce)))
    P_ge1 <- 1.0 / (1.0 + exp(-(alpha1_ind + alpha2 + beta * Ce)))
    P_eq0 <- 1.0 - P_ge1
    P_eq1 <- P_ge1 - P_ge2
    P_eq2 <- P_ge2

    # ============================================================
    # 4. Observation likelihood. The PK observation Cc carries
    # the combined additive + proportional residual error from
    # Savic 2010 Fig. 4 (a and b). The ordered-categorical PD
    # observation has no native rxode2 likelihood form; the
    # typical-value probabilities P_eq{0,1,2} are exposed as
    # derived states for downstream simulation, and a
    # categorical observation can be drawn post-hoc by
    # sampling category index k with probability P_eq_k.
    # See vignette Assumptions and deviations.
    # ============================================================
    Cc ~ add(addSd) + prop(propSd)
  })
}
