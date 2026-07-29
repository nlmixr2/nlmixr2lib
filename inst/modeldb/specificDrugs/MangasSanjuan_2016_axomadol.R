MangasSanjuan_2016_axomadol <- function() {
  description <- paste(
    "Semi-physiological population pharmacokinetic and joint",
    "pharmacodynamic model of axomadol (a racemic analgesic with",
    "opioid agonistic and monoamine-reuptake-inhibitor activity)",
    "and its O-demethyl (ODM) metabolite in healthy adult",
    "volunteers. The PK structure carries two parallel enantiomer",
    "chains (RR-suffix r and SS-suffix s), each consisting of a",
    "first-order absorption depot, a liver compartment mimicking",
    "first-pass conversion, a parent central compartment, and a",
    "metabolite central compartment. Within each enantiomer the",
    "parent and metabolite share a single apparent volume of",
    "distribution (VP = VM) and share a single typical",
    "first-order elimination rate constant (kP0 = kM0), although",
    "between-subject variability is estimated separately for the",
    "two elimination pathways. The PD layer is shared across the",
    "enantiomer chains and is driven by the SS parent in plasma",
    "(mydriatic Emax) and by the RR metabolite at a hysteresis",
    "effect site (linearly miotic). Pupil diameter is the sum of",
    "those two opposing effects; cold-pressor analgesic AUC is a",
    "linear function of the parent and metabolite contributions",
    "to pupil diameter. Parameter values are from Mangas-Sanjuan",
    "et al. 2016 Tables 2, 4, and 5."
  )
  reference <- paste(
    "Mangas-Sanjuan V, Pastor JM, Rengelshausen J, Bursi R,",
    "Troconiz IF. (2016). Population",
    "pharmacokinetic/pharmacodynamic modelling of the effects of",
    "axomadol and its O-demethyl metabolite on pupil diameter and",
    "nociception in healthy subjects. Br J Clin Pharmacol",
    "82(1):112-128. doi:10.1111/bcp.12921.",
    sep = " "
  )
  vignette <- "MangasSanjuan_2016_axomadol"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  # Per-paper named states. The metabolite plasma compartments
  # central_odm_r / central_odm_s and the RR-metabolite effect-site
  # compartment effect_odm_r are paper-mechanistic states (the
  # "_odm" token is the O-demethyl-axomadol metabolite designator
  # used in the source paper). They sit alongside the canonical
  # parent chains depot_r / liver_r / central_r (and _s).
  paper_specific_compartments <- c(
    "central_odm_r",
    "central_odm_s",
    "effect_odm_r"
  )

  # Paper-mechanistic etas whose typical-value parameter is shared
  # with another fixed effect (etalkm0_* share the lkp0_* typical
  # value because Mangas-Sanjuan 2016 Table 2 footnote * constrains
  # kP0 = kM0 typical, but lists separate IIV variances for kP0 and
  # kM0) or whose name encodes a shared multiplicative IIV
  # (etaltotal multiplies the SLP3 and SLP4 contributions to cold-
  # pressor AUC simultaneously per the Table 5 "IIV Total effect"
  # row).
  paper_specific_etas <- c(
    "etalkm0_r",
    "etalkm0_s",
    "etaltotal"
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 74L,
    n_studies      = 2L,
    age_range      = "40-65 years (Study A 46-64, Study B 40-65)",
    weight_range   = "48-100 kg (Study A 61.0-92.2, Study B 48-100)",
    sex_female_pct = 51,
    race_ethnicity = c(White = 100),
    disease_state  = paste(
      "Healthy volunteers, CYP2D6 extensive metabolizers. Subjects",
      "using CYP2D6 inhibitors within the previous 4 weeks were",
      "excluded; Study A genotyped via dextromethorphan urinary",
      "metabolic ratio and Study B via TaqMan CYP2D6 alleles *3, *4,",
      "*5, *6, *7, and *8."
    ),
    dose_range     = paste(
      "Oral axomadol in Study A: 66 mg or 111 mg b.i.d. with a",
      "single-dose day-1 and day-8 PK profile around five days of",
      "twice-daily maintenance. Study B: dose-escalation from b.i.d.",
      "day 1 to t.i.d. days 2-4 at 100 / 125 / 150 mg axomadol,",
      "with a separate group at 150 b.i.d. day 1, 225 b.i.d. days",
      "2-4, and 225 morning of day 5."
    ),
    regions        = "(not reported in the source publication)",
    notes          = paste(
      "Pooled analysis of two Caucasian-volunteer Phase I trials.",
      "Study A: n = 24 (12 male + 14 female), 2-period crossover.",
      "Study B: n = 48 (24 male + 24 female), parallel-group dose",
      "escalation with a 2-way crossover. Methods reports the",
      "study total as n = 74; Table 1 dose-group counts sum to 72",
      "and may differ from n = 74 because of placebo-only or",
      "non-evaluable subjects. The current analysis is in healthy",
      "volunteers, so per Methods 'covariate selection was limited",
      "to dose level' and neither dose level nor time after start",
      "of dosing was retained on any PK or PD parameter (Mangas-",
      "Sanjuan 2016 Population PK modelling, Pupil-diameter PD",
      "modelling, and Cold-pressor PD modelling)."
    )
  )

  ini({
    # -----------------------------------------------------------------
    # PK structural parameters (RR enantiomer) -- Mangas-Sanjuan 2016
    # Table 2 (page numeric estimates with RSE% in parentheses).
    # All rate constants are first-order (1/h); the volume of the
    # parent central compartment VP is also the metabolite central
    # volume VM (Methods: 'VP and VM were assumed to be equal').
    # Bioavailability F1 is fixed at 1 (Methods: 'the typical
    # absolute bioavailability F was assumed to be complete') and
    # carries no IIV in the final model.
    # -----------------------------------------------------------------
    lkgl_r  <- log(1.68)          ; label("RR axomadol first-order absorption rate kGL (1/h)")          # Table 2 RR: 1.68 (RSE 22.4%)
    lklp_r  <- log(18.8)          ; label("RR axomadol liver-to-parent-central rate kLP (1/h)")          # Table 2 RR: 18.8 (RSE 35.3%)
    lkpl_r  <- log(0.0469)        ; label("RR axomadol parent-central-to-liver rate kPL (1/h)")         # Table 2 RR: 4.69e-2 (RSE 11.4%)
    lkp0_r  <- log(0.0793)        ; label("RR axomadol parent + metabolite elimination rate kP0 = kM0 (1/h)")  # Table 2 RR: 7.93e-2 (RSE 3.9%); shared typical kP0|kM0 per footnote *
    lklm_r  <- log(0.0115)        ; label("RR liver-to-metabolite formation rate kLM (1/h)")             # Table 2 RR: 1.15e-2 (RSE 37.0%)
    lvc_r   <- log(424)           ; label("RR apparent parent + metabolite central volume VP/F = VM/F (L)") # Table 2 RR: VP/F = 424 (RSE 3.2%); F is fixed at 1
    ltlag_r <- log(0.26)          ; label("RR absorption lag time tlag (h)")                              # Table 2 RR: 0.26 (RSE 5.7%)

    # -----------------------------------------------------------------
    # PK structural parameters (SS enantiomer) -- Mangas-Sanjuan 2016
    # Table 2.
    # -----------------------------------------------------------------
    lkgl_s  <- log(2.81)          ; label("SS axomadol first-order absorption rate kGL (1/h)")          # Table 2 SS: 2.81 (RSE 27.4%)
    lklp_s  <- log(18.0)          ; label("SS axomadol liver-to-parent-central rate kLP (1/h)")          # Table 2 SS: 18.0 (RSE 29.2%)
    lkpl_s  <- log(0.0997)        ; label("SS axomadol parent-central-to-liver rate kPL (1/h)")         # Table 2 SS: 9.97e-2 (RSE 10.5%)
    lkp0_s  <- log(0.0894)        ; label("SS axomadol parent + metabolite elimination rate kP0 = kM0 (1/h)")  # Table 2 SS: 8.94e-2 (RSE 4.1%); shared typical kP0|kM0 per footnote *
    lklm_s  <- log(0.0157)        ; label("SS liver-to-metabolite formation rate kLM (1/h)")             # Table 2 SS: 1.57e-2 (RSE 30.6%)
    lvc_s   <- log(528)           ; label("SS apparent parent + metabolite central volume VP/F = VM/F (L)") # Table 2 SS: VP/F = 528 (RSE 4.3%); F is fixed at 1
    ltlag_s <- log(0.27)          ; label("SS absorption lag time tlag (h)")                              # Table 2 SS: 0.27 (RSE 4.8%)

    # -----------------------------------------------------------------
    # Pupil-diameter PD parameters (joint across enantiomers) --
    # Mangas-Sanjuan 2016 Table 4 and Equation [1]. The active
    # components on the pupil diameter are SS parent (Emax-type
    # mydriatic action) and RR metabolite at a hysteresis effect
    # site (linear miotic action). Equation [1]:
    #   PDiameter(t) = PDiameter(0) +
    #                  Emax * Cssp / (EC50 + Cssp) -
    #                  SLP2 * CeRRm
    # -----------------------------------------------------------------
    lpd0    <- log(5.45)          ; label("Baseline pupil diameter PDiameter(0) (mm)")                   # Table 4: 5.45 (RSE 1.9%)
    lemax   <- log(0.79)          ; label("Maximum SS-parent-driven pupil-diameter increase Emax (mm)")   # Table 4: 0.79 (RSE 17.4%)
    lec50   <- log(90.7)          ; label("SS-parent plasma EC50 on pupil diameter (ng/mL)")              # Table 4: 90.7 (RSE 27.0%)
    lke0    <- log(0.0124)        ; label("RR-metabolite effect-compartment equilibration rate ke0 (1/h)") # Table 4: ke0 = 1.24e-2 (RSE 2.7%); corresponding t1/2 ~ 55.9 h
    lslp2   <- log(0.00967)       ; label("RR-metabolite effect-site slope on pupil diameter SLP2 (mm/(ng/mL))")  # Table 4: SLP2 = 9.67e-3 (RSE 18.7%)

    # -----------------------------------------------------------------
    # Cold-pressor AUC PD parameters (joint across enantiomers) --
    # Mangas-Sanjuan 2016 Table 5 and Equation [3]:
    #   coldPressorAUC(t) = coldPressorAUC(0) -
    #                        SLP3 * ERRm(t) -
    #                        SLP4 * ESSp(t)
    # where ERRm(t) = SLP2 * CeRRm and ESSp(t) = Emax*Cssp/(EC50+Cssp)
    # are the absolute pupil-diameter contributions of the RR
    # metabolite and SS parent, respectively. The interaction term
    # alpha was not significant (Table 3).
    # -----------------------------------------------------------------
    laucbase <- log(3760)         ; label("Baseline cold-pressor AUC coldPressorAUC(0) (pain unit*s)")    # Table 5: 3760 (RSE 8.2%)
    lslp3    <- log(0.0939)       ; label("Slope of RR-metabolite pupil effect on cold-pressor AUC SLP3 (pain unit*s / mm)")  # Table 5: SLP3 = 9.39e-2 (RSE 43.8%)
    lslp4    <- log(0.23)         ; label("Slope of SS-parent pupil effect on cold-pressor AUC SLP4 (pain unit*s / mm)")      # Table 5: SLP4 = 0.23 (RSE 23.3%)

    # -----------------------------------------------------------------
    # PK between-subject variability (Mangas-Sanjuan 2016 Table 2
    # 'IIV ...' rows). Reported as %CV of log-normally distributed
    # individual parameters; the on-internal-variance scale is
    # omega^2 = log(CV^2 + 1). Per Methods: 'Interindividual
    # variability was supported in all model parameters except the
    # first-order distribution rate constants between the liver and
    # central compartments of the parent compounds, and Relative
    # Bioavailability parameter (F1)' -- so kLP, kPL, and F1 carry
    # no eta. Table 2 lists separate IIV variances for kP0 (parent
    # elimination application) and kM0 (metabolite elimination
    # application), even though the typical values share lkp0_*;
    # encoded as etalkp0_* + etalkm0_* with paper-specific-eta
    # declarations.
    # -----------------------------------------------------------------
    etalkgl_r ~ log(0.77^2 + 1)   # Table 2 RR: IIV kGL = 77% (RSE 13.8%)
    etalkgl_s ~ log(0.66^2 + 1)   # Table 2 SS: IIV kGL = 66% (RSE 46.4%)
    etalkp0_r ~ log(0.37^2 + 1)   # Table 2 RR: IIV kP0 = 37% (RSE 19.7%) on parent elimination
    etalkp0_s ~ log(0.40^2 + 1)   # Table 2 SS: IIV kP0 = 40% (RSE 17.5%)
    etalkm0_r ~ log(0.20^2 + 1)   # Table 2 RR: IIV kM0 = 20% (RSE 21.2%) on metabolite elimination
    etalkm0_s ~ log(0.18^2 + 1)   # Table 2 SS: IIV kM0 = 18% (RSE 19.2%)
    etalklm_r ~ log(0.48^2 + 1)   # Table 2 RR: IIV kLM = 48% (RSE 15.6%)
    etalklm_s ~ log(0.35^2 + 1)   # Table 2 SS: IIV kLM = 35% (RSE 12.0%)
    etalvc_r  ~ log(0.24^2 + 1)   # Table 2 RR: IIV VP  = 24% (RSE 9.7%)
    etalvc_s  ~ log(0.31^2 + 1)   # Table 2 SS: IIV VP  = 31% (RSE 12.6%)
    etaltlag_r ~ log(0.32^2 + 1)  # Table 2 RR: IIV tlag = 32% (RSE 34.8%)
    etaltlag_s ~ log(0.27^2 + 1)  # Table 2 SS: IIV tlag = 27% (RSE 33.6%)

    # -----------------------------------------------------------------
    # Pupil-PD between-subject variability (Mangas-Sanjuan 2016
    # Table 4 'IIV ...' rows). 'The contribution of the nondiagonal
    # elements of the Omega matrix was found to be negligible'.
    # -----------------------------------------------------------------
    etalpd0   ~ log(0.16^2 + 1)   # Table 4: IIV PDiameter(0) = 16% (RSE 9.0%)
    etalslp2  ~ log(0.84^2 + 1)   # Table 4: IIV SLP2 = 84% (RSE 13.1%)
    etalec50  ~ log(2.57^2 + 1)   # Table 4: IIV C50 = 257% (RSE 22.4%); estimated very high

    # -----------------------------------------------------------------
    # Cold-pressor-PD between-subject variability (Mangas-Sanjuan
    # 2016 Table 5). Two etas: one on coldPressorAUC(0), and one
    # shared multiplicative eta on the SLP3+SLP4 combined effect
    # term per Table 5 'IIV Total effect'.
    # -----------------------------------------------------------------
    etalaucbase ~ log(0.68^2 + 1) # Table 5: IIV coldPressorAUC(0) = 68% (RSE 10%)
    etaltotal   ~ log(0.93^2 + 1) # Table 5: IIV Total effect = 93% (RSE 21%); applied multiplicatively to both SLP3 and SLP4 contributions

    # -----------------------------------------------------------------
    # Residual variability (Mangas-Sanjuan 2016 Table 2, 4, 5). The
    # paper transforms parent / metabolite plasma concentrations and
    # cold-pressor AUC logarithmically before fitting and applies an
    # additive residual on the log scale -- expSd (lnorm error in
    # nlmixr2). Pupil diameter is fit on the linear scale with an
    # additive residual.
    #
    # Mangas-Sanjuan 2016 estimated study-specific residuals for the
    # plasma analytes (Study A vs Study B). To preserve a single set
    # of parameters per analyte, the model encodes the Study A
    # residuals here (the smaller, two-period crossover study) and
    # documents the Study B residuals in the validation vignette.
    # -----------------------------------------------------------------
    expSd_r              <- 0.20  ; label("RR-parent plasma residual SD (Study A, log(ng/mL))")          # Table 2 RR: 0.20 (RSE 9.2%); Study B value 0.33 documented in vignette
    expSd_s              <- 0.26  ; label("SS-parent plasma residual SD (Study A, log(ng/mL))")          # Table 2 SS: 0.26 (RSE 13.5%); Study B value 0.49 documented in vignette
    expSd_Cm_r           <- 0.17  ; label("RR-metabolite plasma residual SD (Study A, log(ng/mL))")     # Table 2 RR: 0.17 (RSE 8.3%); Study B value 0.24 documented in vignette
    expSd_Cm_s           <- 0.14  ; label("SS-metabolite plasma residual SD (Study A, log(ng/mL))")     # Table 2 SS: 0.14 (RSE 11.4%); Study B value 0.27 documented in vignette
    addSd_pdiameter      <- 0.46  ; label("Pupil-diameter additive residual SD (mm)")                    # Table 4: 0.46 (RSE 5.1%)
    expSd_coldPressorAUC <- 0.42  ; label("Cold-pressor AUC residual SD (log(pain unit*s))")             # Table 5: 0.42 (RSE 9.4%)
  })

  model({
    # =================================================================
    # Individual PK parameters (RR enantiomer). kLP and kPL carry no
    # IIV per Methods. kP0 (parent elimination) and kM0 (metabolite
    # elimination) share the same typical value lkp0_r but carry
    # separate IIV variances per Table 2.
    # =================================================================
    kgl_r  <- exp(lkgl_r  + etalkgl_r)
    klp_r  <- exp(lklp_r)
    kpl_r  <- exp(lkpl_r)
    kp0_r  <- exp(lkp0_r  + etalkp0_r)
    km0_r  <- exp(lkp0_r  + etalkm0_r)
    klm_r  <- exp(lklm_r  + etalklm_r)
    vc_r   <- exp(lvc_r   + etalvc_r)
    tlag_r <- exp(ltlag_r + etaltlag_r)

    # =================================================================
    # Individual PK parameters (SS enantiomer).
    # =================================================================
    kgl_s  <- exp(lkgl_s  + etalkgl_s)
    klp_s  <- exp(lklp_s)
    kpl_s  <- exp(lkpl_s)
    kp0_s  <- exp(lkp0_s  + etalkp0_s)
    km0_s  <- exp(lkp0_s  + etalkm0_s)
    klm_s  <- exp(lklm_s  + etalklm_s)
    vc_s   <- exp(lvc_s   + etalvc_s)
    tlag_s <- exp(ltlag_s + etaltlag_s)

    # =================================================================
    # Individual PD parameters. The Table-5 'IIV Total effect'
    # variance is applied multiplicatively as exp(etaltotal) to BOTH
    # the SLP3 and SLP4 contributions to cold-pressor AUC (Methods:
    # 'Intersubject variability was estimated for the cold pressor
    # AUC(0) and for [SLP3 * EM-RR(t) + SLP4 * EP-SS(t)]').
    # =================================================================
    pd0     <- exp(lpd0   + etalpd0)
    emax    <- exp(lemax)
    ec50    <- exp(lec50  + etalec50)
    ke0     <- exp(lke0)
    slp2    <- exp(lslp2  + etalslp2)
    aucbase <- exp(laucbase + etalaucbase)
    total_pain_iiv <- exp(etaltotal)
    slp3    <- exp(lslp3) * total_pain_iiv
    slp4    <- exp(lslp4) * total_pain_iiv

    # =================================================================
    # PK ODEs (RR enantiomer chain). Depot dose enters the liver
    # compartment at rate kGL, with reversible exchange between
    # liver and parent central (kLP forward, kPL reverse).
    # Metabolite formation occurs in the liver at rate kLM, with
    # the metabolite then distributing into a central compartment
    # sharing VC with the parent (per Methods: VP = VM).
    # =================================================================
    d/dt(depot_r)       <- -kgl_r * depot_r
    d/dt(liver_r)       <-  kgl_r * depot_r +
                            kpl_r * central_r -
                            klp_r * liver_r -
                            klm_r * liver_r
    d/dt(central_r)     <-  klp_r * liver_r -
                            kpl_r * central_r -
                            kp0_r * central_r
    d/dt(central_odm_r) <-  klm_r * liver_r -
                            km0_r * central_odm_r

    # =================================================================
    # PK ODEs (SS enantiomer chain).
    # =================================================================
    d/dt(depot_s)       <- -kgl_s * depot_s
    d/dt(liver_s)       <-  kgl_s * depot_s +
                            kpl_s * central_s -
                            klp_s * liver_s -
                            klm_s * liver_s
    d/dt(central_s)     <-  klp_s * liver_s -
                            kpl_s * central_s -
                            kp0_s * central_s
    d/dt(central_odm_s) <-  klm_s * liver_s -
                            km0_s * central_odm_s

    # =================================================================
    # Absorption lag times. Applied at the depot for each enantiomer
    # so the rate-constant chain begins kGL * (depot evaluated after
    # tlag).
    # =================================================================
    alag(depot_r) <- tlag_r
    alag(depot_s) <- tlag_s

    # =================================================================
    # Plasma concentrations. Parent and metabolite of the same
    # enantiomer share VC; the two enantiomer chains use enantiomer-
    # specific VCs (per Table 2 vp/F differs RR vs SS). The factor of
    # 1000 converts amount/volume from mg/L into ng/mL so the
    # downstream PD equations (EC50 = 90.7 ng/mL, SLP2 = 9.67e-3
    # mm/(ng/mL)) are dimensionally consistent. Output names:
    #   Cc_r  = RR axomadol parent plasma concentration (ng/mL)
    #   Cc_s  = SS axomadol parent plasma concentration (ng/mL)
    #   Cm_r  = RR O-demethyl-axomadol metabolite plasma (ng/mL)
    #   Cm_s  = SS O-demethyl-axomadol metabolite plasma (ng/mL)
    # The Cm_<r/s> outputs are paper-named (the model file's
    # convention for the O-demethyl metabolite, which has no
    # registered metabolite-suffix token in the nlmixr2lib register).
    # =================================================================
    Cc_r <- central_r     / vc_r * 1000
    Cc_s <- central_s     / vc_s * 1000
    Cm_r <- central_odm_r / vc_r * 1000
    Cm_s <- central_odm_s / vc_s * 1000

    # =================================================================
    # RR-metabolite effect-compartment hysteresis. Sheiner 1979
    # one-compartment effect model with concentration-equality
    # parameterisation: d Ce / dt = ke0 * (Cm - Ce). The state
    # carries the effect-site concentration directly so Cem_r is
    # the effect-site value used by the pupil-diameter equation.
    # =================================================================
    d/dt(effect_odm_r) <- ke0 * (Cm_r - effect_odm_r)
    Cem_r <- effect_odm_r

    # =================================================================
    # Pupil diameter PD output (Mangas-Sanjuan 2016 Equation [1]).
    # The SS parent in plasma contributes a sigmoid-Emax mydriasis;
    # the RR metabolite at the effect site contributes a linear
    # meiosis. The two opposing effects are combined additively to
    # give the observed pupil diameter.
    # =================================================================
    e_parent_s <- emax * Cc_s / (ec50 + Cc_s)
    e_metab_r  <- slp2 * Cem_r
    pdiameter  <- pd0 + e_parent_s - e_metab_r

    # =================================================================
    # Cold-pressor AUC PD output (Mangas-Sanjuan 2016 Equation [3]).
    # Both the RR-metabolite and SS-parent absolute pupil-diameter
    # contributions linearly reduce cold-pressor AUC. The shared
    # multiplicative IIV (total_pain_iiv) has already been applied
    # to slp3 and slp4 above.
    # =================================================================
    coldPressorAUC <- aucbase - slp3 * e_metab_r - slp4 * e_parent_s

    # =================================================================
    # Residual error. Parent / metabolite plasma and cold-pressor
    # AUC carry log-additive (lnorm) residual error per the paper's
    # log-transformation prior to fitting; pupil diameter carries
    # additive residual error on the linear scale.
    # =================================================================
    Cc_r ~ lnorm(expSd_r)
    Cc_s ~ lnorm(expSd_s)
    Cm_r ~ lnorm(expSd_Cm_r)
    Cm_s ~ lnorm(expSd_Cm_s)
    pdiameter ~ add(addSd_pdiameter)
    coldPressorAUC ~ lnorm(expSd_coldPressorAUC)
  })
}
