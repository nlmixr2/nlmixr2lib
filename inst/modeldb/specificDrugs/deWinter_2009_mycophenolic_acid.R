deWinter_2009_mycophenolic_acid <- function() {
  description <- "Semi-mechanistic competitive-protein-binding population PK model for mycophenolic acid (MPA, the active moiety of mycophenolate mofetil MMF) and its glucuronide metabolite MPAG in adult renal transplant recipients (de Winter 2009). Free MPA (fMPA) follows a two-compartment disposition with first-order oral absorption (lag-time TLAG; fixed ka = 4.00 1/h); free MPAG (fMPAG) follows a one-compartment disposition. Both species bind competitively to a saturable plasma protein binding pool with capacity BMAX and species-specific association / dissociation rate constants k24 / k42 (MPA) and k56 / k65 (MPAG). The fMPAG-to-gallbladder transport rate constant k57 drives enterohepatic recirculation: fMPAG accumulates in a gallbladder compartment and empties into the fMPA central compartment during a fixed window (TGB to TGB+DGB post-dose) at rate constant k72, completing the EHC loop. Three covariates: a power effect of creatinine clearance (CRCL) on CL fMPAG (exponent 1.36; CRCL reference 45 mL/min); a power effect of plasma albumin (ALB) on BMAX (exponent 1.39; ALB reference 0.5 mmol/L); and a multiplicative power-form effect of cyclosporine cotreatment (CONMED_CSA) on k57 (multiplier 0.002, reducing EHC by ~99.8% under cyclosporine vs the tacrolimus reference). Total MPA (tMPA) and total MPAG (tMPAG) plasma concentrations are reported as the sum of the unbound and bound concentrations of each species. Dosing is BID by default (tau = 12 h hardcoded in model() for the gallbladder-emptying window). Concentrations are in molar units (umol/L) per the source paper's choice to analyse MMF / MPA / MPAG on a molar basis (MMF MW 433.5; MPA MW 320.3; MPAG MW 496.5)."
  reference <- paste(
    "de Winter BCM, van Gelder T, Sombogaard F, Shaw LM, van Hest RM, Mathot RAA.",
    "Pharmacokinetic role of protein binding of mycophenolic acid and its",
    "glucuronide metabolite in renal transplant recipients.",
    "J Pharmacokinet Pharmacodyn. 2009;36(6):541-564.",
    "doi:10.1007/s10928-009-9136-6.",
    sep = " "
  )
  vignette <- "deWinter_2009_mycophenolic_acid"
  units <- list(time = "hour", dosing = "umol", concentration = "umol/L")

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance computed by the Cockcroft-Gault formula in mL/min (NOT BSA-normalized). Baseline value carried forward across the modeled occasions.",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Cockcroft-Gault raw mL/min, per Methods 'Renal function was tested by calculation of the creatinine clearance (CrCL) according to Cockcroft and Gault [32]'. Used as a power-covariate on CL fMPAG with reference 45 mL/min (cohort median across cyclosporine and tacrolimus arms; Table 1 medians 44 and 45 mL/min). Effect equation derived from text 'A decrease in CrCL from 45 to 25 ml/min resulted in a decrease from 4.75 to 2.14 l/h in clearance of fMPAG' (consistent with power exponent 1.36 reported in Table 2 covariate effects).",
      source_name        = "CrCL"
    ),
    ALB = list(
      description        = "Plasma albumin concentration (mass concentration; SI canonical g/L per the 2026-06-19 register standardization audit).",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source paper reports plasma albumin in mmol/L (Table 1 medians 0.51 mmol/L for both cyclosporine and tacrolimus cohorts; range 0.35-0.68 mmol/L). Canonical column is now SI g/L per the 2026-06-19 register standardization audit; convert inline in model() via `alb_mmolL <- ALB / 66.5` (albumin MW 66.5 g/mmol) to recover the mmol/L value used in the de Winter 2009 calibration. Used as a power-covariate on BMAX with reference 0.5 mmol/L (= 33.25 g/L in SI). The Eq. 8 text of the paper writes the form as P_i = P_pop * exp(theta * (Alb - 0.5)) but the numerical predictions in the Results ('A decrease in albumin from 0.6 to 0.4 mmol/l resulted in a decrease in the number of binding sites from 45200 to 25700 l mol') match a power form P_i = P_pop * (Alb / 0.5)^theta with theta = 1.39 to three significant figures (45100 vs 45200 and 25800 vs 25700); the power form is the one implemented here. See vignette Errata.",
      source_name        = "Alb"
    ),
    CONMED_CSA = list(
      description        = "Concomitant cyclosporine indicator: 1 = patient cotreated with cyclosporine as the calcineurin inhibitor (CNI), 0 = patient cotreated with tacrolimus.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (tacrolimus cotreatment; in the source cohort, n = 45 profiles)",
      notes              = "Time-fixed per patient in the de Winter 2009 cohort because subjects were randomised to one CNI regimen for the entire post-transplant observation window. Used as a multiplicative power-form effect on the fMPAG-to-gallbladder rate constant k57: cyclosporine inhibits MRP2-mediated biliary efflux of MPAG, suppressing enterohepatic recirculation. Tacrolimus does not have this MRP2 effect, so MPA exposure differs between the two regimens for the same MMF dose. Encoded as `k57 <- exp(lk57) * e_csa_k57^CONMED_CSA` with `e_csa_k57 = 0.002` so that k57 drops from 0.0796 1/h under tacrolimus to 0.000159 1/h under cyclosporine (Table 2 / Eq. 9 of de Winter 2009).",
      source_name        = "CsA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 75L,
    n_studies      = 2L,
    n_profiles     = 93L,
    n_observations = "489 tMPA + 489 fMPA + 488 tMPAG + 210 fMPAG plasma concentrations (Results: Patients).",
    age_range      = "19-76 years (cyclosporine cohort 21-70 years; tacrolimus cohort 19-76 years)",
    age_median     = "51 years (cyclosporine cohort) / 53 years (tacrolimus cohort), Table 1",
    weight_range   = "42-113 kg (cyclosporine cohort 42-99 kg; tacrolimus cohort 44-113 kg)",
    weight_median  = "67 kg (cyclosporine cohort) / 78 kg (tacrolimus cohort), Table 1",
    sex_female_pct = 36.0,
    race_ethnicity = "Not reported in source paper.",
    disease_state  = "Adult de novo renal transplant recipients. The cyclosporine cohort received MMF + cyclosporine + corticosteroids; the tacrolimus cohort received MMF + tacrolimus + corticosteroids. Both cohorts pooled from two prior randomised trials (a randomised concentration-controlled trial of MPA AUC target groups + the IMPDH-activity study).",
    dose_range     = "MMF dose adjusted by clinicians; medians 1350 mg BID (range 400-2200, cyclosporine cohort) and 1000 mg BID (500-1500, tacrolimus cohort), Table 1. Default packaged dosing in the vignette is 1000 mg MMF BID = 2306 umol MPA-equivalent (molar dose) per 12 h.",
    regions        = "The Netherlands (Erasmus University Medical Center, Rotterdam) and USA (University of Pennsylvania, Philadelphia).",
    sampling_window= "Up to nine pharmacokinetic occasions post-transplantation per protocol. On day-3, day-7, day-11 (RCCT cohort) and day-6 (IMPDH cohort) occasions, full PK profiles were collected at predose and 0.33, 0.66, 1.25, 2, 6, 8, 12 h post-dose; later occasions used sparser sampling.",
    cni_distribution = "Cyclosporine cotreatment n = 48 profiles (47 patients reported by male/female; 1 missing); tacrolimus cotreatment n = 45 profiles (28 patients reported by male/female; some missing).",
    notes          = "Free-MPA and free-MPAG concentrations were measured via ultrafiltration on a subset of occasions. The final NONMEM run did not minimize successfully (rounding errors) so the source paper does not report parameter standard errors; the visual predictive check confirmed adequate fit. Patient characteristics from Table 1; full cohort summary in Methods 'Patients'."
  )

  ini({
    # Final-model parameter estimates from Table 2 ('Parameter estimates of
    # the pharmacokinetic model'). The source paper's NONMEM final run did
    # not minimize successfully so no standard errors / RSEs are reported in
    # the publication; this packaged model carries only the point estimates.
    #
    # Concentration / amount units throughout this file are molar (umol for
    # amounts; umol/L for concentrations; rate constants involving binding
    # carry units of L/(h*umol) for k_on and 1/h for k_off; BMAX is in umol
    # of total binding sites treated numerically as the saturation value of
    # bound species in the V_binding = 1 L convention used by the source).
    # See the vignette Assumptions and deviations section for the
    # dimensional analysis.

    # Absorption parameters
    ltlag   <- log(0.231);        label("Absorption lag time TLAG (h)")                                                     # Table 2 TLAG = 0.231 h
    lka     <- fixed(log(4.00));  label("First-order absorption rate constant ka (1/h), fixed at 4.00")                     # Table 2 ka = 4.00 1/h (fixed)

    # fMPA disposition (2-compartment, central + peripheral). Apparent CL / V
    # / Q values per Methods 'Because bioavailability (F) could not be
    # quantified, V, CL and Q values correspond to the ratios V/F, CL/F, and
    # Q/F'. CL fMPA represents the rate of conversion of MPA to MPAG (the
    # dominant elimination pathway).
    lvc     <- log(189);          label("Apparent central volume of distribution for fMPA Vc/F (L)")                        # Table 2 Vc fMPA = 189 L
    lcl     <- log(747);          label("Apparent clearance of fMPA to fMPAG CL/F (L/h)")                                   # Table 2 CL fMPA = 747 L/h
    lvp     <- log(34300);        label("Apparent peripheral volume of distribution for fMPA Vp/F (L)")                     # Table 2 Vp fMPA = 34300 L
    lq      <- log(2010);         label("Apparent inter-compartmental clearance for fMPA Q/F (L/h)")                        # Table 2 Q fMPA = 2010 L/h

    # Competitive binding to plasma proteins (Eqs. 1-7 of de Winter 2009).
    # k24 and k56 are bimolecular association rate constants in L/(h*umol);
    # k42 and k65 are first-order dissociation rate constants in 1/h. KD =
    # k_off / k_on -> KD_MPA = 169 / 0.153 = 1104 umol/L (paper reports 1100),
    # KD_MPAG = 93.1 / 0.0133 = 6997 umol/L (paper reports 7000). BMAX is the
    # saturating amount of bound species achievable in the V_binding = 1 L
    # convention; numerically 35100 (treated as umol or equivalently umol/L
    # under the V_binding = 1 L convention).
    lk24    <- log(0.153);        label("Plasma protein binding association rate constant for MPA k24 (L/(h*umol))")       # Table 2 k24 = 0.153 1/(h*umol/L)
    lbmax   <- log(35100);        label("Maximum protein binding capacity BMAX (umol; saturation level at ALB = 0.5 mmol/L)") # Table 2 BMAX = 35100 umol
    lk42    <- log(169);          label("Plasma protein binding dissociation rate constant for MPA k42 (1/h)")              # Table 2 k42 = 169 1/h

    # fMPAG disposition (1-compartment, central only). The fMPAG central
    # compartment receives elimination flux from fMPA (1:1 molar conversion
    # via glucuronidation).
    lvc_mpag  <- log(8.56);       label("Apparent central volume of distribution for fMPAG Vc_mpag/F (L)")                  # Table 2 Vc fMPAG = 8.56 L
    lk56      <- log(0.0133);     label("Plasma protein binding association rate constant for MPAG k56 (L/(h*umol))")       # Table 2 k56 = 0.0133 1/(h*umol/L)
    lk65      <- log(93.1);       label("Plasma protein binding dissociation rate constant for MPAG k65 (1/h)")             # Table 2 k65 = 93.1 1/h
    lcl_mpag  <- log(4.75);       label("Apparent clearance of fMPAG (renal) CL_mpag/F (L/h)")                              # Table 2 CL fMPAG = 4.75 L/h

    # Gallbladder enterohepatic recirculation. fMPAG is continuously
    # transported to the gallbladder via rate constant k57; the gallbladder
    # empties into the fMPA central compartment during a fixed window starting
    # at TGB post-dose with duration DGB at rate constant k72. DGB and k72
    # were fixed during NONMEM estimation due to insufficient data between
    # 4-10 h postdose to identify the gallbladder-emptying kinetics.
    ltgb      <- log(7.90);       label("Time of gallbladder emptying TGB (h post-dose)")                                   # Table 2 TGB = 7.90 h
    ldgb      <- fixed(log(1.00)); label("Duration of gallbladder emptying DGB (h), fixed")                                 # Table 2 DGB = 1.00 h (fixed)
    lk72      <- fixed(log(10.0)); label("Gallbladder-to-fMPA-central rate constant during emptying k72 (1/h), fixed")      # Table 2 k72 = 10.0 1/h (fixed)
    lk57      <- log(0.0796);     label("fMPAG-to-gallbladder transport rate constant k57 (1/h, tacrolimus reference)")     # Table 2 k57 = 0.0796 1/h

    # Covariate effects
    # CRCL on CL fMPAG: power form with reference 45 mL/min. Verified against
    # the text 'A decrease in CrCL from 45 to 25 ml/min resulted in a decrease
    # from 4.75 to 2.14 l/h in clearance of fMPAG': (25/45)^1.36 = 0.451 so
    # CL = 4.75 * 0.451 = 2.14 L/h.
    e_crcl_cl_mpag <- 1.36;       label("Power exponent of CRCL on CL fMPAG (unitless; reference 45 mL/min)")               # Table 2 CrCL on CL fMPAG = 1.36

    # ALB on BMAX: power form with reference 0.5 mmol/L. Numerical
    # predictions in the Results section match a power form rather than the
    # exponential form written in Eq. 8 of the paper (see vignette Errata):
    # 35100 * (0.6/0.5)^1.39 = 45209 (paper reports 45200);
    # 35100 * (0.4/0.5)^1.39 = 25801 (paper reports 25700).
    e_alb_bmax     <- 1.39;       label("Power exponent of plasma ALB on BMAX (unitless; reference 0.5 mmol/L)")            # Table 2 Albumine on BMAX = 1.39

    # CONMED_CSA on k57: multiplicative power-form effect. Table 2 reports
    # the value 0.002 as the multiplier (= 1 + theta_CsA per Eq. 9) so that
    # k57 under cyclosporine = 0.0796 * 0.002 = 0.000159 1/h, matching the
    # Results text 'In patients cotreated with cyclosporine k56 [sic; refers
    # to k57] is very small with a value of 0.000159 h-1 compared to 0.0796
    # h-1 in patients cotreated with tacrolimus'.
    e_csa_k57      <- 0.002;      label("Multiplicative factor on k57 under cyclosporine cotreatment (unitless)")           # Table 2 CsA on k57 = 0.002

    # Inter-patient variability. Source paper reports IPV (interpatient
    # variability) as %CV for exponential error models, translated to
    # log-normal log-scale variance via omega^2 = log(1 + CV^2). TGB IPV in
    # the source paper was described with an 'additional error model'
    # (additive eta) per Results 'IPV was described with an additional error
    # model for TGB and with an exponential error model for k57'; this
    # packaged model uses an exponential / log-normal IIV form on ltgb for
    # uniformity (the 141% CV is encoded as log(1 + 1.41^2) = 1.095). The
    # difference between additive and log-normal at the population mean is
    # negligible for typical-value simulations; see vignette Assumptions and
    # deviations.
    #
    # k57 IPV was fixed at 71% in the source paper 'to prevent variability on
    # EHC to take on extreme values' (per Results), so etalk57 is wrapped in
    # fixed() with variance log(1 + 0.71^2) = 0.4084.
    etaltlag    ~ 1.279                                                                                                       # Table 2 IPV TLAG    = 161% CV -> log(1 + 1.61^2)
    etalvc      ~ 0.853                                                                                                       # Table 2 IPV Vc fMPA = 116% CV -> log(1 + 1.16^2)
    etalcl      ~ 0.663                                                                                                       # Table 2 IPV CL fMPA =  97% CV -> log(1 + 0.97^2)
    etalbmax    ~ 0.207                                                                                                       # Table 2 IPV BMAX    =  48% CV -> log(1 + 0.48^2)
    etalcl_mpag ~ 0.753                                                                                                       # Table 2 IPV CL fMPAG= 106% CV -> log(1 + 1.06^2)
    etaltgb     ~ 1.095                                                                                                       # Table 2 IPV TGB     = 141% CV -> log(1 + 1.41^2); source used additive eta, packaged as log-normal (see vignette)
    etalk57     ~ fixed(0.408)                                                                                                # Table 2 IPV k57     =  71% CV (FIXED) -> fixed(log(1 + 0.71^2))

    # Residual variability. Source paper reports 'Residual variability
    # between observed and predicted MPA plasma concentrations was described
    # using an additional error model for logarithmically transformed data',
    # which is the NONMEM LTBS pattern: an additive error on log-scale
    # corresponds to a proportional error on the linear scale with CV ~ the
    # log-scale SD. The Table 2 values labelled 'Additive error ... (mmol/l)'
    # for the four observed analytes are therefore interpreted here as
    # proportional residual SDs (unitless fractions); the unit annotation in
    # Table 2 is preserved verbatim in the source-trace comment but does not
    # apply to the linear-space proportional residual encoded here. The
    # ordering of magnitudes (tMPAG smallest, tMPA next, fMPAG next, fMPA
    # largest) is consistent with the relative difficulty of measuring the
    # four analytes (the highly-bound species are easier to measure than the
    # free fractions). See vignette Assumptions and deviations.
    propSd          <- 0.520;     label("Proportional residual error on tMPA (fraction)")                                    # Table 2 Additive error tMPA  = 0.52 (source units mmol/l; interpreted as log-scale additive ~ linear-space proportional)
    propSd_fMPA     <- 0.993;     label("Proportional residual error on fMPA (fraction)")                                    # Table 2 Additive error fMPA  = 0.993
    propSd_mpag     <- 0.186;     label("Proportional residual error on tMPAG (fraction)")                                   # Table 2 Additive error tMPAG = 0.186
    propSd_fMPAG    <- 0.551;     label("Proportional residual error on fMPAG (fraction)")                                   # Table 2 Additive error fMPAG = 0.551
  })

  model({
    # ------------------------------------------------------------------
    # SI -> source-paper unit conversion. The canonical ALB column is in
    # SI g/L per the 2026-06-19 register standardization audit; de Winter
    # 2009 calibrated the BMAX power exponent against ALB in mmol/L
    # (reference 0.5 mmol/L = 33.25 g/L). Convert inline so the Table 2
    # exponent (1.39) and reference (0.5) stay load-bearing.
    # ------------------------------------------------------------------
    alb_mmolL <- ALB / 66.5  # SI g/L -> mmol/L (albumin MW 66.5 g/mmol)

    # ------------------------------------------------------------------
    # Individual parameters. Inter-patient variability is log-normal
    # (exponential IIV); the dosing-interval-related parameters (tlag,
    # tgb) carry their etas through exp() / log-normal accordingly.
    # ------------------------------------------------------------------
    tlag    <- exp(ltlag    + etaltlag)
    ka      <- exp(lka)
    vc      <- exp(lvc      + etalvc)
    cl      <- exp(lcl      + etalcl)
    vp      <- exp(lvp)
    q       <- exp(lq)
    k24     <- exp(lk24)
    bmax    <- exp(lbmax    + etalbmax) * (alb_mmolL / 0.5)^e_alb_bmax
    k42     <- exp(lk42)
    vc_mpag <- exp(lvc_mpag)
    k56     <- exp(lk56)
    k65     <- exp(lk65)
    cl_mpag <- exp(lcl_mpag + etalcl_mpag) * (CRCL / 45)^e_crcl_cl_mpag
    tgb     <- exp(ltgb     + etaltgb)
    dgb     <- exp(ldgb)
    k72     <- exp(lk72)
    k57     <- exp(lk57     + etalk57) * e_csa_k57^CONMED_CSA

    # ------------------------------------------------------------------
    # Gallbladder emptying window. The source paper used BID dosing
    # (tau = 12 h); the emptying occurs at TGB post-dose for duration
    # DGB and recurs each dosing interval. The packaged model hardcodes
    # tau = 12 to match the source paper's dose regimen. Users wanting
    # a different dosing interval can edit this line.
    # ------------------------------------------------------------------
    tau           <- 12
    tpost         <- t - floor(t / tau) * tau
    in_emptying   <- (tpost >= tgb) * (tpost <= (tgb + dgb))
    empty_rate    <- k72 * in_emptying

    # ------------------------------------------------------------------
    # Free-species plasma concentrations (umol/L). Used for the binding
    # kinetics and for the unbound-concentration observations.
    # ------------------------------------------------------------------
    fMPA_conc    <- central / vc
    fMPAG_conc   <- central_mpag / vc_mpag

    # ------------------------------------------------------------------
    # Competitive protein binding. Both species bind to the same
    # plasma-protein pool with saturation BMAX (umol; treated as the
    # numerical saturation value of the bound state under the
    # V_binding = 1 L convention used by the source). Binding rate
    # equations:
    #   d(bMPA)/dt  = k24 * [fMPA]  * (BMAX - bMPA - bMPAG) - k42 * bMPA
    #   d(bMPAG)/dt = k56 * [fMPAG] * (BMAX - bMPA - bMPAG) - k65 * bMPAG
    # ------------------------------------------------------------------
    free_sites          <- bmax - complex - complex_mpag
    rate_mpa_bind       <- k24 * fMPA_conc  * free_sites
    rate_mpa_release    <- k42 * complex
    rate_mpag_bind      <- k56 * fMPAG_conc * free_sites
    rate_mpag_release   <- k65 * complex_mpag

    # ------------------------------------------------------------------
    # Inter-compartmental and elimination rates for the fMPA central
    # compartment. CL fMPA represents the conversion of MPA to MPAG; the
    # corresponding flux enters the fMPAG central compartment 1:1 in
    # molar terms (MPA + UDPGA -> MPAG by UGT1A9 / UGT2B7).
    # ------------------------------------------------------------------
    kel        <- cl / vc
    k12_fmpa   <- q  / vc
    k21_fmpa   <- q  / vp
    kel_mpag   <- cl_mpag / vc_mpag

    # ------------------------------------------------------------------
    # ODE system. State variables and units:
    #   depot              (umol)   -- oral MMF / MPA-equivalent depot
    #   central            (umol)   -- fMPA central
    #   peripheral1        (umol)   -- fMPA peripheral
    #   complex            (umol)   -- bMPA (protein-bound MPA)
    #   central_mpag       (umol)   -- fMPAG central
    #   complex_mpag       (umol)   -- bMPAG (protein-bound MPAG)
    #   gallbladder_mpag   (umol)   -- fMPAG accumulated in gallbladder
    # MMF is hydrolyzed to MPA in the gut so the depot dose is
    # expressed in umol MPA-equivalent (1 mol MMF -> 1 mol MPA via ester
    # cleavage in the intestinal wall and plasma).
    # ------------------------------------------------------------------
    # MPA central: oral absorption from depot, elimination to fMPAG, two-
    # compartment distribution, competitive protein binding, and EHC
    # reconversion from gallbladder (MPAG -> MPA upon emptying).
    d/dt(depot)            <- -ka * depot
    d/dt(central)          <- ka * depot - kel * central - k12_fmpa * central + k21_fmpa * peripheral1 - rate_mpa_bind + rate_mpa_release + empty_rate * gallbladder_mpag
    d/dt(peripheral1)      <- k12_fmpa * central - k21_fmpa * peripheral1
    d/dt(complex)          <- rate_mpa_bind - rate_mpa_release
    # MPAG central: 1:1 molar conversion from fMPA, renal elimination,
    # competitive protein binding, and transport to gallbladder.
    d/dt(central_mpag)     <- kel * central - kel_mpag * central_mpag - rate_mpag_bind + rate_mpag_release - k57 * central_mpag
    d/dt(complex_mpag)     <- rate_mpag_bind - rate_mpag_release
    d/dt(gallbladder_mpag) <- k57 * central_mpag - empty_rate * gallbladder_mpag

    # Absorption lag for the oral depot.
    lag(depot) <- tlag

    # ------------------------------------------------------------------
    # Observation outputs (all in umol/L). The bound-species "amounts"
    # complex / complex_mpag are interpreted as concentrations directly
    # under the V_binding = 1 L convention; tMPA / tMPAG are computed as
    # the sum of free and bound concentrations (paper Results: 'The
    # concentrations of tMPA and tMPAG were modeled as the sum of the
    # unbound and bound concentrations').
    # ------------------------------------------------------------------
    Cc        <- fMPA_conc  + complex            # tMPA  (total MPA  in umol/L)
    fMPA      <- fMPA_conc                       # fMPA  (free  MPA  in umol/L)
    Cc_mpag   <- fMPAG_conc + complex_mpag       # tMPAG (total MPAG in umol/L)
    fMPAG     <- fMPAG_conc                      # fMPAG (free  MPAG in umol/L)

    Cc        ~ prop(propSd)
    fMPA      ~ prop(propSd_fMPA)
    Cc_mpag   ~ prop(propSd_mpag)
    fMPAG     ~ prop(propSd_fMPAG)
  })
}
