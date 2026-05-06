Schindler_2016_sunitinib <- function() {
  description <- "Joint pharmacodynamic model for sunitinib in advanced GIST coupling a five-lesion indirect-response model of [18F]FDG-PET SUVmax with a per-subject sum-of-longest-diameters (SLD) tumor-growth-inhibition module and constant-baseline-hazard Weibull time-to-event sub-models for overall survival and study dropout (Schindler 2016 / DDMODEL00000221). Sunitinib exposure enters via an effect compartment driven by a per-day AUC = DAILY_DOSE / CL_SUNITINIB; the OS hazard depends on the per-subject week-1 maximum across-lesion relative SUVmax change from baseline, RCFB1MAX."
  reference <- paste(
    "Schindler E, Amantea MA, Karlsson MO, Friberg LE. (2016).",
    "PK-PD modeling of individual lesion FDG-PET response to predict overall",
    "survival in patients with sunitinib-treated gastrointestinal stromal tumor.",
    "CPT Pharmacometrics Syst Pharmacol 5(4):173-181.",
    "doi:10.1002/psp4.12057.",
    "DDMORE Foundation Model Repository: DDMODEL00000221.",
    sep = " "
  )
  vignette     <- "Schindler_2016_sunitinib"
  units        <- list(time = "hour", dosing = "mg", concentration = "n/a (non-PK outputs only: SUVmax unitless and SLD in mm)")
  ddmore_id    <- "DDMODEL00000221"
  replicate_of <- NULL

  covariateData <- list(
    CL_SUNITINIB = list(
      description        = "Per-subject post-hoc apparent oral clearance of sunitinib carried in from a previously-developed sunitinib popPK model.",
      units              = "L/h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Used together with `DAILY_DOSE` to drive the effect-compartment AUC = DAILY_DOSE / CL_SUNITINIB. The original Schindler 2016 NONMEM dataset supplies this column as `CL` (`$INPUT` comment: 'post-hoc clearance from previously-developed PK model'); the upstream sunitinib popPK model is not part of the DDMODEL00000221 bundle and is not currently in nlmixr2lib. The vignette virtual cohort uses a single literature-typical adult sunitinib CL of 50 L/h, consistent with Houk et al. (2010) J Clin Pharmacol 50:843-858, and references Houk 2010 narratively rather than reproducing its popPK structure inline.",
      source_name        = "CL"
    ),
    DAILY_DOSE = list(
      description        = "Time-varying daily sunitinib dose in mg, switching between the prescribed dose level on dosing-cycle records and 0 on off-cycle / dose-holiday records.",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. The Schindler 2016 NONMEM dataset supplies this column as `DOS` (`$INPUT` comment: 'daily dose in mg'); for a sunitinib 4-weeks-on / 2-weeks-off schedule the column toggles between the prescribed daily-dose level (e.g., 50 mg/day) and 0 at each cycle boundary. Drives the per-day AUC fed to the effect compartment.",
      source_name        = "DOS"
    ),
    RCFB1MAX = list(
      description        = "Per-subject scalar predictor for the overall-survival Weibull hazard, defined as the maximum (across the up-to-five tracked target lesions) of the relative change in SUVmax at one week of sunitinib therapy: max((SUVmax(t = 168 h) - SUVmax(0)) / SUVmax(0)).",
      units              = "(unitless)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Sign convention: more negative `RCFB1MAX` (greater week-1 SUVmax suppression) reduces the OS hazard. In the source NONMEM `.mod` `RCFB1MAX` is computed inline at FLAG = 1 / TIME = 168 h from the running SUVmax compartment values and reused on subsequent records — a NONMEM record-loop construct without an idiomatic rxode2 / nlmixr2 equivalent. The model file therefore consumes `RCFB1MAX` as a per-subject input covariate; reproducing the source's behavior requires a two-stage simulation (run the SUVmax + SLD ODEs first, compute `RCFB1MAX` per subject from the t = 168 h SUVmax values, then run the OS / dropout TTE arms with `RCFB1MAX` bound). The vignette virtual cohort follows this pattern.",
      source_name        = "RCFB1MAX"
    )
  )

  population <- list(
    n_subjects     = 66L,
    n_studies      = 1L,
    age_range      = NA_character_,
    age_median     = NA_character_,
    weight_range   = NA_character_,
    weight_median  = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = "Imatinib-resistant or imatinib-intolerant advanced gastrointestinal stromal tumor (GIST) on second-line oral sunitinib therapy. The Schindler 2016 dataset comprises 66 GIST patients pooled from a Phase I/II program (Demetri et al. 2009; Casali et al. 2008-style cohort) who underwent serial [18F]FDG-PET SUVmax assessments and target-lesion sum-of-longest-diameters (SLD) measurements over the sunitinib 4-weeks-on / 2-weeks-off cycle, with overall survival followed to event or right-censoring.",
    dose_range     = "Oral sunitinib 50 mg/day on a 4-weeks-on / 2-weeks-off schedule (the standard GIST regimen). The DDMODEL00000221 simulated dataset uses 50 mg/day during dose-on weeks and 0 mg/day during the off-cycle.",
    regions        = NA_character_,
    notes          = "Detailed baseline demographics (age, weight, sex, race / ethnicity distribution, prior-line distribution) for the 66-patient cohort are reported in Schindler 2016 Table 1 of the linked publication; this publication is not on disk in the worktree, so the per-field values are not reproduced here. The DDMODEL00000221 bundle's `Simulated_SLD_SUV_OS_GIST.csv` ships a single virtual subject with the standard 50 mg/day 4-on/2-off schedule as a regression-style smoke test rather than a representative cohort. See the validation vignette's Errata section for the full list of bundle-versus-publication caveats."
  )

  ini({
    # ----- Structural SUVmax indirect-response parameters -----
    # Final estimates from Output_real_SLD_SUV_OS_GIST.lst, FINAL PARAMETER
    # ESTIMATE block (lines 868-882). Unit conversions follow the source `.mod`
    # `$PK` block (lines 54, 63, 72, 81): the listed THETA values are scaled
    # (typically /1000/24/7 or /24/7) to give 1/hour rate constants matching
    # the dataset's TIME column in hours. Each `log()` wraps the back-
    # transformed 1/hour rate so that `exp(l...)` inside `model()` recovers it.
    lkout    <- log(555.634 / 1000 / 24 / 7) ; label("KOUT: SUVmax loss-of-response rate constant (1/h)")  # THETA(1) = 555.634 FIX (1/week * 1000); .mod line 54 KOUT = THETA(1)/1000/24/7 -> 0.003307 1/h
    ldrug    <- log(0.94569)                 ; label("DRUG_SUV: sunitinib drug-effect coefficient on SUVmax (per mg*h/L of effect-compartment AUC)")  # THETA(3) = 0.94569 FIX; .mod line 72 TVDRUG = THETA(3)
    lbase_suv <- log(7.58866)                ; label("Typical baseline SUVmax (unitless)")                  # THETA(4) = 7.58866 FIX; .mod line 81 TVBASE = THETA(4)

    # Disease-progression slope. The .mod fixes ALPHA (and all its lesion-
    # specific etas) at zero because no SUVmax disease progression was
    # identified during model building. Kept here as `fixed(0)` rather than
    # log-transformed so the structural form `KIN_i*(1 + alpha*t)` is
    # preserved if a future user wants to re-estimate ALPHA on a different
    # dataset.
    alpha    <- fixed(0)                     ; label("ALPHA: SUVmax disease-progression slope (1/h^2); FIX = 0 in source model")  # THETA(2) = 0 FIX; .mod line 63 TVALPHA = THETA(2)/1000/24/7

    # ----- Structural SLD tumor-growth parameters -----
    lkg      <- log(10.4796 / 1000 / 24 / 7) ; label("K_GROW: tumor growth rate constant on SLD (1/h)")    # THETA(7) = 10.4796 FIX (KG*1000 in 1/week); .mod line 111 TVKG = THETA(7)/24/7/1000
    llambda  <- log(20.0529 / 1000 / 24 / 7) ; label("LAMBDA: drug-resistance decay rate constant (1/h)")  # THETA(8) = 20.0529 FIX (LAMBDA*1000 in 1/week); .mod line 113 TVLAM = THETA(8)/24/7/1000
    lkdrug   <- log(16.5971 / 1000 / 24 / 7) ; label("DRUG_SLD: sunitinib drug-effect rate constant on SLD (1/h per mg*h/L of effect-compartment AUC)")  # THETA(9) = 16.5971 FIX (KDRUG*1000 in 1/week); .mod line 115 TVDRU = THETA(9)/24/7/1000
    lbase_sld <- log(263.064)                ; label("SLD_0: baseline sum of longest diameters (mm)")      # THETA(10) = 263.064 FIX; .mod line 117 IBASE = THETA(10)

    # ----- Time-to-event parameters (overall survival + dropout) -----
    # ALPHH and ALPHD are FIX = 1, degenerating the Weibull pdf into a
    # constant-baseline-hazard form. THETA(12), THETA(14), THETA(16) are
    # the only THETAs estimated during minimization; final estimates from
    # iteration 8 of the .lst (lines 788-790, NPARAMETR row): 1.9030E-02,
    # 1.9924E-02, 5.3558E+00.
    llambh     <- log(0.019030 / 24 / 7)     ; label("LAMBH: overall-survival Weibull scale parameter (1/h)")  # THETA(12) final = 0.019030 (1/week per .mod line 123 LAMBH = THETA(12)/24/7); converted to 1/h
    alphh      <- fixed(1)                   ; label("ALPHH: overall-survival Weibull shape (unitless); FIX = 1 yields constant baseline hazard")  # THETA(13) = 1 FIX; .mod line 124
    llambd     <- log(0.019924 / 24 / 7)     ; label("LAMBD: dropout Weibull scale parameter (1/h)")            # THETA(14) final = 0.019924 (1/week per .mod line 125 LAMBD = THETA(14)/24/7); converted to 1/h
    alphd      <- fixed(1)                   ; label("ALPHD: dropout Weibull shape (unitless); FIX = 1 yields constant baseline hazard")  # THETA(15) = 1 FIX; .mod line 126
    theta_pred <- 5.3558                     ; label("Coefficient of RCFB1MAX in the OS hazard: hazard *= exp(theta_pred * RCFB1MAX) (unitless)")  # THETA(16) final = 5.3558; .mod line 175 DADT(7) ... * EXP(THETA(16)*RCFB1MAX)

    # ----- Inter-individual variability -----
    # Final estimates from Output_real_SLD_SUV_OS_GIST.lst FINAL PARAMETER
    # ESTIMATE OMEGA block (lines 928-1000). All omegas are declared FIX
    # in the source `.mod` and were not modified during minimization.
    #
    # ETAs 1-12 (KOUT, ALPHA + 5 SAME-block lesion deviations each) are
    # FIX = 0 in the .mod and are therefore omitted from the nlmixr2 ini
    # (no IIV / no inter-lesion variability on KOUT or ALPHA was identified
    # during model building, per the .mod comments at lines 53 and 62).
    #
    # BLOCK(2) on (ETA13 = KDRUG SLD typical IIV, ETA14 = DRUG SUV typical
    # IIV): variances 0.398557 and 0.548112 with covariance 0.397941
    # (.lst lines 928-934; .mod $OMEGA BLOCK(2) FIX block at lines 375-377).
    etalkdrug + etaldrug ~ fixed(c(0.398557, 0.397941, 0.548112))

    # BLOCK(1) SAME on the 5 lesion-specific drug-effect deviations
    # (ETA15-ETA19 in the .mod): variance 0.324481 each, no cross-lesion
    # covariance (.lst lines 942-960; .mod $OMEGA BLOCK(1) ... SAME at
    # lines 379-384). nlmixr2 has no `SAME` shortcut, so the five etas
    # are spelled out individually with the same FIXed variance.
    etaldrug_les1 ~ fixed(0.324481)
    etaldrug_les2 ~ fixed(0.324481)
    etaldrug_les3 ~ fixed(0.324481)
    etaldrug_les4 ~ fixed(0.324481)
    etaldrug_les5 ~ fixed(0.324481)

    # Independent diagonal omegas on (ETA20 = SLD baseline IIV) and
    # (ETA21 = typical-subject SUV baseline IIV). .lst lines 962-968;
    # .mod $OMEGA at lines 386-387.
    etalbase_sld ~ fixed(0.288669)
    etalbase     ~ fixed(0.105177)

    # BLOCK(1) SAME on the 5 lesion-specific baseline deviations
    # (ETA22-ETA26 in the .mod): variance 0.0549036 each, no cross-lesion
    # covariance (.lst lines 970-990; .mod $OMEGA BLOCK(1) ... SAME at
    # lines 389-394).
    etalbase_les1 ~ fixed(0.0549036)
    etalbase_les2 ~ fixed(0.0549036)
    etalbase_les3 ~ fixed(0.0549036)
    etalbase_les4 ~ fixed(0.0549036)
    etalbase_les5 ~ fixed(0.0549036)

    # Diagonal omega on (ETA27 = K_GROW SLD growth-rate IIV); ETA28
    # (LAMBDA IIV) is FIX = 0 and is omitted. .lst lines 992-1000;
    # .mod $OMEGA at lines 396-397.
    etalkg ~ fixed(0.361258)

    # ----- Residual error -----
    # SUVmax outputs (5 lesions): the source `.mod` `$ERROR` block uses a
    # log-transform-both-sides residual model, `IPRED = LOG(A(i))` and
    # `Y = IPRED + EPSN`, where EPSN is one of five EPS components
    # generated by Cholesky-decomposition of a 5x5 residual-error matrix
    # with diagonal variance THETA(5) = 0.173794 and across-lesion
    # correlation THETA(6) = 0.466827 (.lst lines 880, 1010-1026; .mod
    # lines 90-106 and 188-225). NONMEM "additive on log scale" maps to
    # nlmixr2's `prop()` (proportional in linear space); the per-lesion
    # standard deviation is sqrt(THETA(5)) ≈ 0.4169. The Cholesky-
    # induced cross-lesion correlation (THETA(6) ≈ 0.467) is not
    # represented natively because nlmixr2 does not support cross-output
    # residual-error correlation in its standard `~ prop()` syntax;
    # marginal per-lesion residual variance is preserved. See the
    # validation vignette's "Assumptions and deviations" section.
    propSd_suv1 <- fixed(sqrt(0.173794)) ; label("SUVmax lesion 1 proportional residual error (sqrt of variance THETA(5))")  # THETA(5) = 0.173794 FIX; .lst line 880
    propSd_suv2 <- fixed(sqrt(0.173794)) ; label("SUVmax lesion 2 proportional residual error (sqrt of variance THETA(5))")  # THETA(5) = 0.173794 FIX
    propSd_suv3 <- fixed(sqrt(0.173794)) ; label("SUVmax lesion 3 proportional residual error (sqrt of variance THETA(5))")  # THETA(5) = 0.173794 FIX
    propSd_suv4 <- fixed(sqrt(0.173794)) ; label("SUVmax lesion 4 proportional residual error (sqrt of variance THETA(5))")  # THETA(5) = 0.173794 FIX
    propSd_suv5 <- fixed(sqrt(0.173794)) ; label("SUVmax lesion 5 proportional residual error (sqrt of variance THETA(5))")  # THETA(5) = 0.173794 FIX

    # SLD output: linear-space proportional residual error from the
    # `Y = IPRED + W1 * EPS(6)` block with W1 = IPRED * THETA(11). .mod
    # lines 253-263; .lst lines 880-881.
    propSd_sld  <- fixed(0.0667619)      ; label("SLD proportional residual error (THETA(11), fraction)")  # THETA(11) = 0.0667619 FIX
  })

  model({
    # ============================================================
    # 1. Derived inputs (per-subject covariates -> internal scalars)
    # ============================================================
    # Daily AUC of sunitinib over a 24-hour dosing interval. The source
    # .mod $PK block (line 22) computes AUC = DOS / CL on every record
    # and feeds it to the effect compartment. With DAILY_DOSE in mg/day
    # and CL_SUNITINIB in L/h the product has units mg*h/L, which is the
    # standard concentration-time-area unit for AUC_24,ss; A(9) (effect
    # compartment) carries the same units.
    auc_daily <- DAILY_DOSE / CL_SUNITINIB

    # ============================================================
    # 2. Individual structural parameters
    # ============================================================
    # SUVmax indirect-response sub-model. `kout` and `alpha` are
    # population-typical (no IIV / inter-lesion variability identified
    # in the source model). Lesion-specific drug effect and baseline
    # combine a typical-subject eta (etaldrug, etalbase) with a
    # lesion-specific eta drawn independently from the BLOCK(1) SAME
    # variance (etaldrug_lesK, etalbase_lesK).
    #
    # The composite `theta + eta_typical + eta_lesion` form in the
    # source `.mod` cannot be expressed on a single mu-referenced
    # line in nlmixr2 (mu-reference allows at most one eta per
    # theta-anchored line). The typical-subject contribution is
    # built first as a single mu-referenced product; each lesion
    # then multiplies by its own `exp(lesion-specific eta)`.
    kout       <- exp(lkout)
    drug_typ   <- exp(ldrug + etaldrug)
    base_typ   <- exp(lbase_suv + etalbase)
    drug1  <- drug_typ * exp(etaldrug_les1)
    drug2  <- drug_typ * exp(etaldrug_les2)
    drug3  <- drug_typ * exp(etaldrug_les3)
    drug4  <- drug_typ * exp(etaldrug_les4)
    drug5  <- drug_typ * exp(etaldrug_les5)
    ibase1 <- base_typ * exp(etalbase_les1)
    ibase2 <- base_typ * exp(etalbase_les2)
    ibase3 <- base_typ * exp(etalbase_les3)
    ibase4 <- base_typ * exp(etalbase_les4)
    ibase5 <- base_typ * exp(etalbase_les5)

    # SLD tumor-growth sub-model parameters.
    kg        <- exp(lkg + etalkg)
    lambda    <- exp(llambda)
    kdrug     <- exp(lkdrug + etalkdrug)
    base_sld  <- exp(lbase_sld + etalbase_sld)

    # Time-to-event sub-model parameters.
    lambh <- exp(llambh)
    lambd <- exp(llambd)

    # Effect-compartment equilibration rate (KEO = ln(2)/50, with the
    # 50-hour half-life from the source .mod line 118; in 1/h).
    keo <- log(2) / 50

    # ============================================================
    # 3. Compartment initial conditions
    # ============================================================
    # SUVmax compartments start at the per-subject per-lesion baseline.
    suv1(0) <- ibase1
    suv2(0) <- ibase2
    suv3(0) <- ibase3
    suv4(0) <- ibase4
    suv5(0) <- ibase5
    # SLD compartment starts at the per-subject SLD baseline.
    sld(0)  <- base_sld
    # Cumulative-hazard compartments and the effect compartment start
    # at zero by default; left implicit.

    # ============================================================
    # 4. ODE system
    # ============================================================
    # Effect compartment: equilibrates to the daily AUC at rate keo.
    d/dt(effect) <- keo * (auc_daily - effect)

    # SUVmax indirect-response model for each lesion: KIN = KOUT * IBASE
    # (so that A(t = 0) = IBASE in steady state with no drug); EFF = drug
    # * effect-compartment value drives loss-of-response. The (1 + alpha
    # * t) progression term is preserved with alpha = 0 so the structural
    # form remains visible.
    d/dt(suv1) <- kout * ibase1 * (1 + alpha * t) - kout * (1 + drug1 * effect) * suv1
    d/dt(suv2) <- kout * ibase2 * (1 + alpha * t) - kout * (1 + drug2 * effect) * suv2
    d/dt(suv3) <- kout * ibase3 * (1 + alpha * t) - kout * (1 + drug3 * effect) * suv3
    d/dt(suv4) <- kout * ibase4 * (1 + alpha * t) - kout * (1 + drug4 * effect) * suv4
    d/dt(suv5) <- kout * ibase5 * (1 + alpha * t) - kout * (1 + drug5 * effect) * suv5

    # SLD tumor-growth-inhibition model. Resistance term exp(-lambda * t)
    # progressively reduces drug efficacy, so the tumor regrows late in
    # treatment even at constant exposure.
    d/dt(sld)  <- kg * sld - kdrug * effect * exp(-lambda * t) * sld

    # Cumulative-hazard compartments for overall survival and dropout.
    # With ALPHH = ALPHD = 1 the Weibull instantaneous hazard reduces to
    # a constant baseline (* exp(theta_pred * RCFB1MAX) for OS), so the
    # cumulative hazards are linear in time.
    d/dt(cumHaz_os)   <- lambh * alphh * (t + 1e-16)^(alphh - 1) * exp(theta_pred * RCFB1MAX)
    d/dt(cumHaz_drop) <- lambd * alphd * (t + 1e-16)^(alphd - 1)

    # ============================================================
    # 5. Observation model
    # ============================================================
    # SUVmax outputs: log-additive residual error in the source model
    # maps to nlmixr2's prop() (proportional in linear space) per
    # lesion. Cross-lesion residual correlation is dropped; see
    # vignette Errata.
    suv1 ~ prop(propSd_suv1)
    suv2 ~ prop(propSd_suv2)
    suv3 ~ prop(propSd_suv3)
    suv4 ~ prop(propSd_suv4)
    suv5 ~ prop(propSd_suv5)
    # SLD output: linear-space proportional residual error.
    sld  ~ prop(propSd_sld)
  })
}
