Clewe_2018_rifampicin <- function() {
  description <- "Multistate Tuberculosis Pharmacometric (MTP) model coupled with the General Pharmacodynamic Interaction (GPDI) model for the triple combination of rifampicin (RIF), isoniazid (INH), and ethambutol (EMB) against in vitro Mycobacterium tuberculosis B1585 (Clewe 2018, scenario = 4). Three bacterial subpopulations (fast-multiplying Fbugs, slow-multiplying Sbugs, non-replicating Nbugs) exchange via first-order rates and a time-dependent F-to-S transfer; INH adaptive resistance is captured by a two-state ARON / AROFF system that dynamically shifts the INH EC50 on the F and S subpopulations. Each drug acts on each subpopulation through a Hill or hyperbolic exposure-response, combined across drugs via Bliss independence on Fbugs and linear addition on Sbugs and Nbugs; pairwise GPDI interaction parameters shift the Emax / EC50 of each affected drug-effect term. Drug exposures (RIF, INH, EMB) are time-fixed in vitro concentrations supplied as data covariates."
  reference <- paste(
    "Clewe O, Wicha SG, de Vogel CP, de Steenwinkel JEM, Simonsson USH. (2018).",
    "A model-informed preclinical approach for prediction of clinical pharmacodynamic interactions of anti-tuberculosis drug combinations.",
    "J Antimicrob Chemother 73(2):437-447.",
    "doi:10.1093/jac/dkx380.",
    "DDMORE Foundation Model Repository: DDMODEL00000259 (scenario = 4, triple-combination MTP-GPDI).",
    sep = " "
  )
  vignette <- "Clewe_2018_rifampicin"
  units <- list(time = "day", dosing = "CFU/mL", concentration = "mg/L")
  ddmore_id    <- "DDMODEL00000259"
  replicate_of <- NULL

  covariateData <- list(
    RIF = list(
      description        = "Rifampicin in vitro exposure concentration (time-fixed per replicate).",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source dataset column `RIF`. Constant per experimental replicate (the .mod initialises a static drug compartment from this column and sets DADT = 0). Supply 0 for RIF-free arms; the GPDI interaction terms multiplying RIF reduce to zero in that case.",
      source_name        = "RIF"
    ),
    INH = list(
      description        = "Isoniazid in vitro exposure concentration (time-fixed per replicate).",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source dataset column `INH`. Constant per experimental replicate; drives the Hill exposure-response on Fbugs / Sbugs and the kinetics of the ARON / AROFF adaptive-resistance switch (kon * AROFF * INH).",
      source_name        = "INH"
    ),
    EMB = list(
      description        = "Ethambutol in vitro exposure concentration (time-fixed per replicate).",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source dataset column `EMB`. Constant per experimental replicate. The .mod source guards five GPDI parameters (sdieh, sdihe, fdier, sdier, sdierh) behind `IF(A(8) > 0)` blocks so that EMB-mediated interaction shifts collapse to zero in EMB-free arms; this implementation reproduces that behaviour with `(EMB > 0) * <param>` factors.",
      source_name        = "EMB"
    )
  )

  population <- list(
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    disease_state  = "in vitro Mycobacterium tuberculosis B1585 culture (no human subjects).",
    dose_range     = "Initial inoculum f0 + s0 = 209 + 324 = 533 (THETA(7..8) units; the .mod scales these by 1000 in $PK so the bacterial-state initial conditions are 209,000 + 324,000 = 533,000 CFU/mL with a non-replicating seed of 1e-5 CFU/mL). RIF, INH, EMB exposures are static in vitro concentrations supplied per replicate; the source dataset spans single-drug, two-drug, and three-drug exposures across NATG=1 normal-growth experiments (EXPR codes 1,2,3,7,8,11,13).",
    regions        = NA_character_,
    notes          = "In vitro time-kill model fit to M. tuberculosis B1585 cultures under combinations of rifampicin, isoniazid, and ethambutol. The DDMORE bundle reports `Scenario = 4`, which the .mod and the publication identify as the triple-combination MTP-GPDI fit (MTP block fixed from an earlier mono-data MTP fit; INH adaptive-resistance and per-drug exposure-response parameters fixed from the corresponding mono-data INH / RIF / EMB fits; only the GPDI interaction parameters were estimated jointly on the combination data). The companion `Clewe_2016_rifampicin` task covers the underlying single-drug RIF / MTP model (DDMODEL00000220 lineage); the present model adds the multidrug GPDI layer."
  )

  ini({
    # ============================================================
    # Final estimates from DDMODEL00000259's `Executable_MTP-GPDI.mod`
    # $THETA / $OMEGA / $SIGMA blocks (the .mod's initial values are the
    # scenario-4 final estimates carried forward from the publication;
    # the bundled `Output_simulated_MTP-GPDI.lst` re-prints these as the
    # FINAL PARAMETER ESTIMATE of an `$ESTIMATION MAXEVAL=0` evaluation
    # run. The shipped `Output_real_MTP-GPDI.lst` is a $TABLE export
    # rather than a NONMEM listing, so a direct `MINIMIZATION
    # SUCCESSFUL` cross-check is not available; see the validation
    # vignette's Errata section.) Comments cite the .mod line and slot.
    # ============================================================

    # ---- MTP model: bacterial sub-population kinetics ----
    # All MTP parameters are FIXED in scenario 4 because they were
    # estimated upstream on the single-drug data and held fixed when the
    # GPDI interaction parameters were fit on the combination data.
    kg      <- fixed(0.796)        ; label("Exponential growth-rate constant kG of the F (fast-multiplying) subpopulation (1/day)")          # .mod $THETA(1) FIX
    kfslin  <- fixed(0.166 / 100)  ; label("Time-dependent F-to-S transfer-rate slope kFSLIN (1/day^2; KFS = kfslin * t)")                   # .mod $THETA(2)/100 FIX
    kfn     <- fixed(0.897 / 1e6)  ; label("First-order F-to-N transfer rate kFN (1/day)")                                                    # .mod $THETA(3)/1e6 FIX
    ksf     <- fixed(0.145 / 10)   ; label("First-order S-to-F transfer rate kSF (1/day)")                                                    # .mod $THETA(4)/10 FIX
    ksn     <- fixed(0.186)        ; label("First-order S-to-N transfer rate kSN (1/day)")                                                    # .mod $THETA(5) FIX
    kns     <- fixed(0.123 / 100)  ; label("First-order N-to-S transfer rate kNS (1/day)")                                                    # .mod $THETA(6)/100 FIX
    f0      <- fixed(209 * 1000)   ; label("Initial F (fast-multiplying) bacterial number F0 (CFU/mL)")                                     # .mod $THETA(7)*1000 FIX
    s0      <- fixed(324 * 1000)   ; label("Initial S (slow-multiplying) bacterial number S0 (CFU/mL)")                                     # .mod $THETA(8)*1000 FIX

    # ---- INH (H) exposure-response (FIXED from mono-data fit) ----
    hfdemax <- fixed(22.2209)      ; label("INH Emax on Fbugs (FD; 1/h-equivalent killing-rate scalar inside the Bliss combination)")        # .mod $THETA(9) FIX
    hfdec50 <- fixed(0.168431)     ; label("INH EC50 on Fbugs (FD; mg/L) - reference, modulated by adaptive resistance")                     # .mod $THETA(10) FIX
    hfdgam  <- fixed(1.90157)      ; label("INH Hill exponent on Fbugs (FD)")                                                                # .mod $THETA(11) FIX
    hsdemax <- fixed(8.55316)      ; label("INH Emax on Sbugs (SD; 1/h)")                                                                    # .mod $THETA(12) FIX
    hsdec50 <- fixed(0.0328672)    ; label("INH EC50 on Sbugs (SD; mg/L) - reference, modulated by adaptive resistance")                     # .mod $THETA(13) FIX
    hsdgam  <- fixed(1.74098)      ; label("INH Hill exponent on Sbugs (SD)")                                                                # .mod $THETA(14) FIX

    # ---- INH adaptive resistance (FIXED) ----
    kon      <- fixed(0.0205994)   ; label("INH adaptive-resistance kON: rate of AROFF -> ARON conversion ((mg/L)^-1 day^-1, multiplied by INH)") # .mod $THETA(15) FIX
    koff     <- fixed(0)           ; label("INH adaptive-resistance kOFF: rate of ARON -> AROFF reversion (1/day)")                            # .mod $THETA(16) FIX (= 0)
    arlinfd  <- fixed(522.42)      ; label("Linear-resistance slope on the INH EC50 for Fbugs ARLINFD (unitless)")                            # .mod $THETA(17) FIX
    arlinsd  <- fixed(2352.28)     ; label("Linear-resistance slope on the INH EC50 for Sbugs ARLINSD (unitless)")                            # .mod $THETA(18) FIX

    # ---- RIF (R) exposure-response (FIXED from mono-data fit) ----
    rfdemax <- fixed(1.96922)      ; label("RIF Emax on Fbugs (FD; ratio relative to hfdemax: the Bliss equation pre-multiplies by hfdemax)") # .mod $THETA(19) FIX
    rfdec50 <- fixed(0.0030258)    ; label("RIF EC50 on Fbugs (FD; mg/L)")                                                                   # .mod $THETA(20) FIX
    rfgemax <- fixed(1)            ; label("RIF growth-inhibition Emax on Fbugs growth-rate (FG; unitless 1 - sigmoid)")                     # .mod $THETA(21) FIX
    rfgec50 <- fixed(0.388407)     ; label("RIF growth-inhibition EC50 on Fbugs growth-rate (FG; mg/L)")                                     # .mod $THETA(22) FIX
    rfggam  <- fixed(2.80234)      ; label("RIF growth-inhibition Hill exponent on Fbugs growth-rate (FG)")                                  # .mod $THETA(23) FIX
    rsdemax <- fixed(1.79211)      ; label("RIF Emax on Sbugs (SD; 1/h)")                                                                    # .mod $THETA(24) FIX
    rsdec50 <- fixed(0.0112798)    ; label("RIF EC50 on Sbugs (SD; mg/L)")                                                                   # .mod $THETA(25) FIX
    rndk    <- fixed(3.28587)      ; label("RIF first-order kill-rate constant on Nbugs ((mg/L)^-1 day^-1)")                                   # .mod $THETA(26) FIX

    # ---- EMB (E) exposure-response (FIXED from mono-data fit) ----
    efdemax <- fixed(2.2073)       ; label("EMB Emax on Fbugs (FD; ratio relative to hfdemax)")                                              # .mod $THETA(27) FIX
    efdec50 <- fixed(0.860332)     ; label("EMB EC50 on Fbugs (FD; mg/L)")                                                                   # .mod $THETA(28) FIX
    efdgam  <- fixed(2.45751)      ; label("EMB Hill exponent on Fbugs (FD)")                                                                # .mod $THETA(29) FIX
    esdk    <- fixed(4.38781)      ; label("EMB first-order kill-rate constant on Sbugs ((mg/L)^-1 day^-1, scaled by interaction terms)")      # .mod $THETA(30) FIX

    # ---- GPDI: INH(H) - RIF(R) PD interaction ----
    fdirh   <- -0.683351           ; label("GPDI: RIF -> INH-FD interaction shift on EC50_H_FD (unitless; lower bound -1)")                  # .mod $THETA(31) (estimable, lower bound -1)
    fdihr   <- fixed(0)            ; label("GPDI: INH -> RIF-FD interaction shift on EC50_R_FD (unitless; FIX = 0)")                         # .mod $THETA(32) FIX (= 0)
    sdirh   <- 1.52908             ; label("GPDI: RIF -> INH-SD interaction shift on EC50_H_SD (unitless)")                                  # .mod $THETA(33) (estimable)
    sdihr   <- 10.7494             ; label("GPDI: INH -> RIF-SD interaction shift on EC50_R_SD (unitless)")                                  # .mod $THETA(34) (estimable)

    # ---- GPDI: INH(H) - EMB(E) PD interaction ----
    fdieh   <- 1.80936             ; label("GPDI: EMB -> INH-FD interaction shift on EC50_H_FD (unitless)")                                  # .mod $THETA(35) (estimable)
    fdihe   <- fixed(0)            ; label("GPDI: INH -> EMB-FD interaction shift on EC50_E_FD (unitless; FIX = 0)")                         # .mod $THETA(36) FIX (= 0)
    sdieh   <- 0.0854746           ; label("GPDI: EMB -> INH-SD interaction shift on EC50_H_SD (unitless; conditional on EMB > 0)")          # .mod $THETA(37) (estimable, IF EMB > 0)
    sdihe   <- 91.4222             ; label("GPDI: INH -> EMB-SD interaction shift on EMB-SD term (unitless; conditional on INH > 0)")        # .mod $THETA(38) (estimable, IF INH > 0)

    # ---- GPDI: RIF(R) - EMB(E) PD interaction ----
    fdier   <- -0.662812           ; label("GPDI: EMB -> RIF-FD interaction shift on EC50_R_FD (unitless; conditional on EMB > 0)")          # .mod $THETA(39) (estimable, IF EMB > 0)
    fdire   <- fixed(-0.99999)     ; label("GPDI: RIF -> EMB-FD interaction shift on EC50_E_FD (unitless; FIX at -0.99999, near-floor)")     # .mod $THETA(40) FIX
    sdier   <- 1.70602             ; label("GPDI: EMB -> RIF-SD interaction shift on EC50_R_SD (unitless; conditional on EMB > 0)")          # .mod $THETA(41) (estimable, IF EMB > 0)
    sdire   <- 479.458             ; label("GPDI: RIF -> EMB-SD interaction shift on EMB-SD term (unitless)")                                # .mod $THETA(42) (estimable)

    # ---- GPDI: EMB(E) modulating the RIF-INH interaction ----
    sdierh  <- -0.677036           ; label("GPDI: EMB -> (RIF-INH SDIRH) interaction shift on EC50_H_SD (unitless; conditional on EMB > 0)") # .mod $THETA(43) (estimable, IF EMB > 0)

    # ---- Residual error: additive on the natural-log bacterial-density scale (M3 method
    # in the source .mod for below-LOQ data; the censoring is handled at fit time via
    # the data's `cens` column rather than inside the model file). ----
    addSd   <- fixed(sqrt(0.936573))  ; label("Additive residual SD on log(Fbugs + Sbugs); SD = sqrt(SIGMA(1)) on the natural-log scale") # .mod $SIGMA(1) FIX = 0.936573 (variance)
  })

  model({
    # ====================================================================
    # 1. Time-dependent F-to-S transfer rate (.mod $DES `KFS = KFSLIN*T`).
    # `t` is the simulation time in rxode2.
    # ====================================================================
    kfs <- kfslin * t

    # ====================================================================
    # 2. INH adaptive-resistance modulation of the INH EC50s.
    # The .mod $DES uses A(5) = aron in `AREFD = 1 + ARLINFD*A(5)` and
    # `ARESD = 1 + ARLINSD*A(5)` to inflate the INH EC50 as adaptive
    # resistance accumulates.
    # ====================================================================
    arefd     <- 1 + arlinfd * aron
    aresd     <- 1 + arlinsd * aron
    hfdec50a  <- hfdec50 * arefd       # adaptive-resistance-shifted INH EC50 on Fbugs
    hsdec50a  <- hsdec50 * aresd       # adaptive-resistance-shifted INH EC50 on Sbugs

    # ====================================================================
    # 3. EMB-presence guards on the GPDI parameters that the source .mod
    # gates with `IF(A(8) > 0) THEN ... ELSE ... = 0` (and an INH-presence
    # guard on sdihe). Since RIF / INH / EMB are time-fixed covariates
    # supplied per replicate, the guard reduces to a multiplicative factor
    # `(EMB > 0)` / `(INH > 0)`.
    # ====================================================================
    emb_on    <- (EMB > 0)
    inh_on    <- (INH > 0)
    sdieh_eff  <- emb_on * sdieh
    sdihe_eff  <- inh_on * sdihe
    fdier_eff  <- emb_on * fdier
    sdier_eff  <- emb_on * sdier
    sdierh_eff <- emb_on * sdierh

    # ====================================================================
    # 4. Drug-effect terms - reproduce the .mod $DES verbatim.
    # The Bliss-independence combination on FD (.mod L156) requires INHFD,
    # RIFFD, EMBFD on a 0..1 scale (each is a Hill ratio with Emax = 1).
    # The final FD multiplies the Bliss combination by hfdemax (so the
    # absolute kill-rate scale is set by hfdemax). Compare:
    #   .mod L132 INHFD = 1 * A(4)^HFDGAM / ((HFDEC50A * (...) * (...))^HFDGAM + A(4)^HFDGAM)
    #   .mod L138 RIFFD = (RFDEMAX/HFDEMAX) * A(7) / ((RFDEC50 * (...) * (...)) + A(7))
    #   .mod L151 EMBFD = (EFDEMAX/HFDEMAX) * A(8)^EFDGAM / ((EFDEC50 * (...) * (...))^EFDGAM + A(8)^EFDGAM)
    # ====================================================================

    # INH effect on Fbugs (saturating; Emax = 1 internally, scaled to hfdemax via Bliss).
    inh_ec50_fd <- hfdec50a *
      (1 + (fdirh * RIF / (rfdec50 + RIF))) *
      (1 + (fdieh * EMB / (efdec50 + EMB)))
    inhfd <- INH^hfdgam / (inh_ec50_fd^hfdgam + INH^hfdgam)

    # INH effect on Sbugs (saturating with absolute Emax = hsdemax).
    inh_ec50_sd <- hsdec50a *
      (1 + (sdirh * (1 + sdierh_eff) * RIF / (rsdec50 + RIF))) *
      (1 + sdieh_eff)
    inhsd <- hsdemax * INH^hsdgam / (inh_ec50_sd^hsdgam + INH^hsdgam)

    # RIF effect on Fbugs (Hill exponent = 1 in source; saturating).
    rif_ec50_fd <- rfdec50 *
      (1 + (fdihr * INH / (hfdec50a + INH))) *
      (1 + fdier_eff)
    riffd <- (rfdemax / hfdemax) * RIF / (rif_ec50_fd + RIF)

    # RIF growth-inhibition effect on Fbugs growth-rate (1 - sigmoidal).
    # The source .mod has a stale `IF(RIFEFG.LT.0) RIFEFG=0` guard that
    # references an undeclared variable (`RIFEFG` vs the declared `RIFFG`)
    # and so never fires; with rfgemax = 1 the expression is mathematically
    # bounded to [0, 1] and the guard is unnecessary. See vignette Errata.
    rif_growth_inhib <- rfgemax * RIF^rfggam / (rfgec50^rfggam + RIF^rfggam)
    riffg <- 1 - rif_growth_inhib

    # RIF effect on Sbugs (saturating).
    rif_ec50_sd <- rsdec50 *
      (1 + (sdihr * INH / (hsdec50a + INH))) *
      (1 + sdier_eff)
    rifsd <- rsdemax * RIF / (rif_ec50_sd + RIF)

    # RIF first-order kill on Nbugs (linear in RIF).
    rifnd <- rndk * RIF

    # EMB effect on Fbugs.
    emb_ec50_fd <- efdec50 *
      (1 + (fdihe * INH / (hfdec50a + INH))) *
      (1 + (fdire * RIF / (rfdec50 + RIF)))
    embfd <- (efdemax / hfdemax) * EMB^efdgam / (emb_ec50_fd^efdgam + EMB^efdgam)

    # EMB effect on Sbugs - reproduced verbatim from .mod L154; the
    # parenthesisation `(1 + sdire*RIF/rsdec50 + RIF)` was carried over
    # from the source even though it parses as a non-standard sum rather
    # than the canonical Hill-shift `RIF/(rsdec50 + RIF)`. Flagged in the
    # vignette's Errata section as a likely upstream typo that
    # nonetheless ships with the published scenario-4 fit.
    embsd <- EMB * (esdk / ((1 + sdihe_eff * INH) * (1 + sdire * RIF / rsdec50 + RIF)))

    # ====================================================================
    # 5. Bliss-independence combination on F, linear summation on S and N.
    # FG only sees RIF (.mod L160 `FG = RIFFG`). ND only sees RIF.
    # ====================================================================
    fd <- (inhfd + riffd + embfd
           - inhfd * riffd - inhfd * embfd - riffd * embfd
           + inhfd * riffd * embfd) * hfdemax
    fg <- riffg
    sd <- inhsd + rifsd + embsd
    nd <- rifnd

    # ====================================================================
    # 6. Growth function on Fbugs: GROWTHFUNC = kg with a non-negativity
    # floor at 0 (.mod L120). With kg = 0.796 fixed positive and constant
    # in time, the floor is decorative; reproduced verbatim for
    # source-trace fidelity.
    # ====================================================================
    growthfunc <- kg
    if (growthfunc < 0) growthfunc <- 0

    # ====================================================================
    # 7. Bacterial-state ODEs (.mod $DES L168-170) and adaptive-resistance
    # ODEs (.mod $DES L172-173). The .mod's INH / RIF / EMB compartments
    # (DADT(4) = DADT(7) = DADT(8) = 0) are dropped here in favour of the
    # scalar covariates RIF / INH / EMB; the observable bacterial outputs
    # are unchanged.
    # ====================================================================
    Fbugs(0) <- f0
    Sbugs(0) <- s0
    Nbugs(0) <- 1e-5         # .mod $PK A_0(3) = 0.00001
    aron(0)  <- 0            # .mod $PK A_0(5) = 0
    aroff(0) <- 1            # .mod $PK A_0(6) = 1

    d/dt(Fbugs) <- Fbugs * fg * growthfunc + ksf * Sbugs - kfs * Fbugs - kfn * Fbugs - fd * Fbugs
    d/dt(Sbugs) <- kfs * Fbugs + kns * Nbugs - ksn * Sbugs - ksf * Sbugs - sd * Sbugs
    d/dt(Nbugs) <- ksn * Sbugs + kfn * Fbugs - kns * Nbugs - nd * Nbugs
    d/dt(aron)  <- kon * aroff * INH - koff * aron
    d/dt(aroff) <- koff * aron - kon * aroff * INH

    # ====================================================================
    # 8. Observation: log of the replicating bacterial mass (Fbugs + Sbugs),
    # matching .mod $ERROR L189 `IPRED = LOG(A(1) + A(2))`. Residual error
    # is additive on the natural-log scale (.mod $SIGMA additive on log).
    # `totBugs` and `nBugs` are exposed as auxiliary outputs for users who
    # want to inspect Fbugs + Sbugs + Nbugs and the persister sub-population
    # alongside the fitted observation.
    # ====================================================================
    totBugs   <- Fbugs + Sbugs + Nbugs
    nBugs     <- Nbugs
    logFSbugs <- log(Fbugs + Sbugs)
    logFSbugs ~ add(addSd)
  })
}
