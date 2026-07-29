Zhang_2012_lopinavir_ritonavir <- function() {
  description <- paste(
    "Simultaneous integrated population pharmacokinetic model of oral",
    "lopinavir (LPV, parent) and ritonavir (RTV, sibling-drug suffix _rtv)",
    "in 21 HIV-infected South African adults with and without concomitant",
    "antitubercular rifampicin (Zhang 2012). Structure: LPV one-compartment",
    "with first-order absorption (ka 0.991 1/h) and LPV CL/F dynamically",
    "inhibited by RTV plasma concentration via a sigmoid Imax (Imax = 0.953,",
    "IC50 = 0.0351 mg/L); RTV two-compartment with a Savic transit-",
    "compartment absorption chain (NN = 2.03, MTT = 1.44 h) feeding RTV",
    "depot at rate ktr = (NN+1)/MTT and absorbed to RTV central at ka_rtv",
    "= 3.28 1/h. Allometric scaling fixed at the Holford / Anderson",
    "literature values: fat-free mass (Janmahasatian) drives CL/F (exponent",
    "0.75) and total body weight drives Vc/F and Vp/F (exponent 1.0).",
    "Rifampicin (CONMED_RIF) increases LPV CL/F by 71.0% and RTV CL/F by",
    "36.0%, reduces LPV F by 20.0% and RTV F by 45.0% (at the 100 mg",
    "reference RTV dose), and the RTV F when on rifampicin scales upward",
    "with RTV dose at 8.1% per 10 mg above the 100 mg reference (saturation",
    "of first-pass metabolism / P-gp self-inhibition; identifiable only",
    "within the RIF-coadministered arm of the source study). Diurnal",
    "variation is encoded via the simulation convention t = clock-hours-",
    "from-midnight: doses given during the overnight window (clock 20:00 to",
    "08:00) carry +42.0% (LPV) and +45.0% (RTV) relative bioavailability",
    "vs morning doses, and oral CL/F of both drugs is reduced by 32.7%",
    "overnight."
  )
  reference <- paste(
    "Zhang C, Denti P, Decloedt E, Maartens G, Karlsson MO, Simonsson USH,",
    "McIlleron H. Model-based approach to dose optimization of",
    "lopinavir/ritonavir when co-administered with rifampicin.",
    "Br J Clin Pharmacol. 2012;73(5):758-767.",
    "doi:10.1111/j.1365-2125.2011.04154.x."
  )
  vignette <- "Zhang_2012_lopinavir_ritonavir"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description       = "Total body weight (baseline or per-record).",
      units             = "kg",
      type              = "continuous",
      notes             = "Allometric scaling of LPV / RTV central and peripheral volumes (V/F) referenced to 70 kg with exponent fixed at 1.0 per Holford / Anderson convention (Zhang 2012 Methods: 'allometric scaling was applied to ... volume of distribution (V/F) ... as described by Holford et al.'). The paper does not explicitly state the reference weight; 70 kg is the Holford-school convention and is documented in the vignette Assumptions and deviations section.",
      source_name       = "WT"
    ),
    FFM = list(
      description       = "Fat-free mass derived from total body weight, height, and sex via the Janmahasatian formula.",
      units             = "kg",
      type              = "continuous",
      notes             = "Allometric scaling of LPV / RTV apparent clearance (CL/F) referenced to 50 kg with exponent fixed at 0.75 per Holford / Anderson convention (Zhang 2012 Methods: 'fat free mass ... as described by Holford et al.'; Janmahasatian 2005 formula referenced for FFM derivation). The paper does not explicitly state the reference FFM; 50 kg is the Holford-school convention (approximate FFM of a 70 kg adult male) and is documented in the vignette Assumptions and deviations section. Users supply FFM directly; the Janmahasatian formula is FFM_M = 9.27e3 * WT / (6.68e3 + 216 * BMI) for males and FFM_F = 9.27e3 * WT / (8.78e3 + 244 * BMI) for females, with BMI = WT / (HT/100)^2.",
      source_name       = "FFM"
    ),
    CONMED_RIF = list(
      description       = "Binary indicator of concomitant rifampicin co-administration (chronic 600 mg once daily for >=7 days, post-induction steady state).",
      units             = "(binary)",
      type              = "binary",
      reference_category = "0 (no rifampicin)",
      notes             = "Set to 1 for the rifampicin co-administration arms (PK2, PK3, PK4 in Zhang 2012); set to 0 for the baseline arm (PK1: LPV/r 400/100 mg twice daily without rifampicin). The paper assumes full CYP3A4 induction is reached 1 week after starting daily 600 mg rifampicin (Methods Discussion). The indicator drives the multiplicative effects on LPV CL/F (+71.0%), RTV CL/F (+36.0%), LPV F (-20.0%), RTV F (-45.0% at the 100 mg reference RTV dose), and gates the RTV dose-dependent F effect (active only when CONMED_RIF = 1).",
      source_name       = "RIF"
    ),
    DOSE = list(
      description       = "Per-record ritonavir dose level in mg used by the RTV-bioavailability dose-effect term.",
      units             = "mg",
      type              = "continuous",
      notes             = "Anchored at 100 mg (reference). The 8.1% multiplicative increment per 10 mg of RTV dose is applied only when CONMED_RIF = 1 (the dose-by-bioavailability interaction was identified only within the rifampicin-coadministered arm of the source study; Zhang 2012 Results 'Model description' paragraph 5). Set DOSE = 100 mg for the standard 400/100 mg LPV/r regimen, 150 mg for 600/150 mg LPV/r, 200 mg for 800/200 mg LPV/r. Outside the 100-200 mg range the linear extrapolation is unvalidated and may produce non-physiological F values.",
      source_name       = "(per-protocol RTV dose level; not explicitly tabulated in the source)"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 21L,
    n_studies        = 1L,
    n_observations   = 800L,
    age_range        = "26-58 years",
    age_median       = "36 years",
    weight_range     = "43.0-110.0 kg",
    weight_median    = "64.5 kg",
    height_range     = "148.0-186.5 cm",
    height_median    = "160.5 cm",
    bmi_range        = "17.4-41.4 kg/m^2",
    bmi_median       = "26.7 kg/m^2",
    ffm_range        = "30.6-65.9 kg",
    ffm_median       = "39.5 kg",
    sex_female_pct   = 85.7,
    race_ethnicity   = "South African adults; specific race / ethnicity composition not separately reported by Zhang 2012",
    disease_state    = "HIV-1 infection without active tuberculosis; protease-inhibitor-naive and virologically suppressed on LPV/r plus two NRTIs at study entry (Zhang 2012 Methods 'Study design and drug analysis').",
    dose_range       = paste(
      "Four sequential treatment conditions tested:",
      "(PK1) LPV/r 400/100 mg twice daily without rifampicin (reference);",
      "(PK2) LPV/r 400/100 mg twice daily + rifampicin 600 mg once daily;",
      "(PK3) LPV/r 600/150 mg twice daily + rifampicin 600 mg once daily;",
      "(PK4) LPV/r 800/200 mg twice daily + rifampicin 600 mg once daily.",
      "Morning doses given after overnight fast; evening doses given with a meal.",
      "Intensive sampling at 0 (pre-dose), 1.5, 2, 2.5, 3, 4, 5, 6, 8, 12 h after morning dose,",
      "one week after each dose adjustment."
    ),
    regions          = "South Africa (Cape Town, University of Cape Town)",
    notes            = paste(
      "Three patients withdrew before completing all four occasions (two transaminitis, one nausea);",
      "partial data retained. All patients had 100% pill-count adherence. Lopinavir LLOQ 0.05 mg/L,",
      "ritonavir LLOQ 0.025 mg/L; values below LLOQ (1% LPV, 2% RTV) set to LLOQ/2 in the model build.",
      "Bootstrap confidence intervals based on 10 samples only owing to model complexity; CIs reported",
      "in Table 2 are illustrative rather than precise. The model has no IIV on the lopinavir-ritonavir",
      "interaction parameters (Imax, IC50) because the dynamic perpetrator-substrate relationship is",
      "shared at the population typical-value level. The reference body weight (70 kg) and reference",
      "fat-free mass (50 kg) are Holford / Anderson convention defaults; the source paper does not",
      "explicitly state either reference."
    )
  )

  ini({
    # ==========================================================================
    # LOPINAVIR (LPV, parent) structural parameters
    # ----- Zhang 2012 Table 2 'Lopinavir Estimates' column -----
    # LPV: one-compartment open model with first-order absorption.
    # ==========================================================================
    lka <- log(0.991);   label("Lopinavir absorption rate constant ka (1/h)")                                 # Table 2 LPV ka = 0.991 h^-1 (95% CI 0.63-1.43)
    lcl <- log(37.9);    label("Lopinavir baseline apparent clearance CL/F when no ritonavir effect (L/h)")   # Table 2 LPV CL/F = 37.9 L/h (95% CI 28.5-52.1); footnote dagger: CL/F of lopinavir without ritonavir
    lvc <- log(54.7);    label("Lopinavir apparent central volume of distribution Vc/F (L)")                  # Table 2 LPV Vc/F = 54.7 L (95% CI 50.5-64.7)
    lfdepot <- fixed(log(1));  label("Lopinavir reference relative bioavailability F (unitless, fixed at 1)") # Table 2 footnote: relative bioavailability of standard dose LPV/r without rifampicin assumed 100% reference

    # ==========================================================================
    # RITONAVIR (RTV, sibling-drug suffix _rtv) structural parameters
    # ----- Zhang 2012 Table 2 'Ritonavir Estimates' column -----
    # RTV: two-compartment open model with Savic-style transit-compartment
    # absorption (NN = 2.03 transit compartments, continuous-valued) feeding a
    # central first-order absorption ka_rtv = 3.28 1/h. NN is non-integer; the
    # rxode2 built-in transit() function accepts the continuous parameterisation
    # via the closed-form gamma-PDF input rate (Savic 2007).
    # ==========================================================================
    lka_rtv <- log(3.28);    label("Ritonavir absorption rate constant ka (1/h)")                                # Table 2 RTV ka = 3.28 h^-1 (95% CI 2.90-3.38)
    lcl_rtv <- log(19.2);    label("Ritonavir apparent oral clearance CL/F (L/h)")                               # Table 2 RTV CL/F = 19.2 L/h (95% CI 18.4-22.2)
    lvc_rtv <- log(22.6);    label("Ritonavir apparent central volume of distribution Vc/F (L)")                 # Table 2 RTV Vc/F = 22.6 L (95% CI 21.9-24.6)
    lq_rtv  <- log(31.0);    label("Ritonavir intercompartmental clearance Q/F (L/h)")                           # Table 2 RTV Q/F = 31.0 L/h (95% CI 25.7-34.7)
    lvp_rtv <- log(56.6);    label("Ritonavir peripheral volume of distribution Vp/F (L)")                       # Table 2 RTV Vp/F = 56.6 L (95% CI 50.8-66.0)
    lmtt_rtv <- log(1.44);   label("Ritonavir mean transit time MTT through the absorption chain (h)")           # Table 2 RTV MTT = 1.44 h (95% CI 1.39-1.53)
    lnn_rtv  <- log(2.03);   label("Ritonavir number of Savic-style transit compartments NN (unitless)")         # Table 2 RTV NN = 2.03 (95% CI 1.83-2.37); continuous parameter, supports non-integer values
    lfdepot_rtv <- fixed(log(1));  label("Ritonavir reference relative bioavailability F (unitless, fixed at 1)") # Table 2 footnote: relative bioavailability of standard dose LPV/r without rifampicin assumed 100% reference

    # ==========================================================================
    # Drug-drug interaction: dynamic sigmoidal Imax inhibition of LPV CL/F by
    # RTV plasma concentration (Zhang 2012 Methods 'Population pharmacokinetic
    # analysis' equation; Table 2 'Lopinavir-ritonavir interaction' block).
    #     CL/F_LPV(t) = CL0/F_LPV * (1 - Imax * C_RTV(t) / (IC50 + C_RTV(t)))
    # Bare names imax / ic50 (linear scale) follow the Schipani 2013
    # atazanavir + ritonavir precedent (Schipani_2013_atazanavir_ritonavir.R).
    # ==========================================================================
    imax <- 0.953;       label("Maximum fractional inhibition of LPV CL/F by RTV (unitless, bounded [0,1])")    # Table 2 LPV-RTV interaction Emax = 95.3% (95% CI 94.5%-96.3%)
    ic50 <- 0.0351;      label("Ritonavir plasma concentration producing 50% of Imax on LPV CL/F (mg/L)")       # Table 2 LPV-RTV interaction EC50 = 0.0351 mg/L (95% CI 0.0194-0.0438)

    # ==========================================================================
    # Allometric scaling exponents -- fixed at the Holford / Anderson convention
    # (Zhang 2012 Methods 'Population pharmacokinetic analysis' paragraph 4:
    # 'allometric scaling was applied to apparent clearance (CL/F) and volume
    # of distribution (V/F) ... as described by Holford et al.'). The paper
    # reports that 'fat free mass was more appropriate for the scaling of
    # clearances' and 'total body weight was found to be appropriate for
    # allometric scaling of central and peripheral volumes' (Results 'Model
    # description' paragraph 4).
    # ==========================================================================
    e_ffm_cl <- fixed(0.75);  label("Allometric exponent of FFM on CL/F for both LPV and RTV (unitless, fixed at Holford convention)")
    e_wt_vc  <- fixed(1.0);   label("Allometric exponent of WT on Vc/F for both LPV and RTV (unitless, fixed at Holford convention)")
    e_wt_vp  <- fixed(1.0);   label("Allometric exponent of WT on Vp/F for RTV (unitless, fixed at Holford convention)")

    # ==========================================================================
    # Concomitant rifampicin (CONMED_RIF) effects (Zhang 2012 Table 2 rows
    # 'RIF on CL/F (+)' and 'Relative bioavailability when given with RIF').
    # Implemented as multiplicative perturbations gated by the binary
    # CONMED_RIF indicator.
    # ==========================================================================
    e_rif_cl     <- 0.710;  label("Multiplicative effect of CONMED_RIF on LPV CL/F (fraction increase)")                         # Table 2 RIF on LPV CL/F (+) = 71.0% (95% CI 65.7%-75.4%)
    e_rif_cl_rtv <- 0.360;  label("Multiplicative effect of CONMED_RIF on RTV CL/F (fraction increase)")                         # Table 2 RIF on RTV CL/F (+) = 36.0% (95% CI 35.2%-40.0%)
    e_rif_f      <- -0.20;  label("Multiplicative effect of CONMED_RIF on LPV F (fraction change)")                              # Table 2 LPV Relative F with RIF = 0.80 (95% CI 0.76-0.85); -20% vs no-RIF reference
    e_rif_f_rtv  <- -0.45;  label("Multiplicative effect of CONMED_RIF on RTV F at the 100 mg reference RTV dose (fraction change)") # Table 2 RTV Relative F with RIF = 0.55 (95% CI 0.54-0.59) at standard 100 mg RTV dose; -45% vs no-RIF reference; the dose-dependent boost (e_dose_f_rtv) is applied multiplicatively on top of this RIF reduction
    e_dose_f_rtv <- 0.0081; label("Multiplicative effect on RTV F per mg of RTV dose above the 100 mg reference (fraction per mg, active only when CONMED_RIF = 1)") # Table 2 RTV 'Bioavailability/10 mg ritonavir (+)' = 8.1% (95% CI 5.7%-11.2%); per-10-mg coefficient -> 0.081/10 = 0.0081 per mg; reproduces Results 'Model description' paragraph 5: RTV F with RIF at 100 mg = 0.55, at 150 mg = 0.77, at 200 mg = 0.996

    # ==========================================================================
    # Diurnal variation (Zhang 2012 Results 'Model description' paragraph 6 and
    # Table 2 rows 'Evening effect on bioavailability (+)' and 'Evening effect
    # on CL/F (-)'). Encoded via the simulation convention t = clock-hours-
    # from-midnight: the overnight indicator is 1 during the 20:00-08:00 window
    # (post-evening dose + pre-morning dose) and 0 during the daytime 08:00-
    # 20:00 window. The same evening / overnight indicator drives the dose-time
    # bioavailability boost (+42.0% LPV, +45.0% RTV) and the time-varying CL/F
    # reduction (-32.7% for both drugs).
    # ==========================================================================
    e_eve_f        <-  0.420;  label("Multiplicative effect of evening-window dose on LPV F (fraction increase)")        # Table 2 Evening on LPV F (+) = 42.0% (95% CI 38.0%-48.2%)
    e_eve_f_rtv    <-  0.450;  label("Multiplicative effect of evening-window dose on RTV F (fraction increase)")        # Table 2 Evening on RTV F (+) = 45.0% (95% CI 41.4%-53.6%)
    e_overnight_cl <- -0.327;  label("Multiplicative effect of the overnight window on LPV and RTV CL/F (fraction change)") # Table 2 Evening on CL/F (-) = 32.7% (95% CI 29.6%-38.4%); shared between LPV and RTV per paper

    # ==========================================================================
    # Inter-individual variability (Zhang 2012 Table 2 rows 'IIV CL/F (%CV)',
    # 'IIV V/F (%CV)', 'IIV F (%CV)' and 'IOV ...' for ka / F / CL/F / MTT /
    # RUV). The paper reports both IIV (between-subject) and IOV (inter-
    # occasion) variability separately. Per the nlmixr2lib convention used in
    # Bienczak_2016_nevirapine.R, Chirehwa_2017_pyrazinamide.R, and
    # Svensson_2014_bedaquiline.R: BOV is dropped where a BSV term is reported
    # on the same parameter, and BOV is folded in as a BSV-equivalent where
    # only BOV is reported. The IOV on residual unexplained variability (17.1%
    # for LPV) cannot be straightforwardly encoded as a BSV-equivalent on a
    # residual-error magnitude and is dropped with a note in the vignette.
    # omega^2 = log(1 + CV^2) converts the CV%-on-SD scale to the log-scale
    # variance.
    # ==========================================================================
    etalcl     ~ 0.04004    # Table 2 LPV IIV CL/F = 20.2% (95% CI 12.7%-25.1%); IOV CL/F 11.8% dropped per convention; omega^2 = log(1 + 0.202^2) = 0.04004
    etalvc     ~ 0.07113    # Table 2 LPV IIV Vc/F = 27.2% (95% CI 10.3%-41.4%); omega^2 = log(1 + 0.272^2) = 0.07113
    etalka     ~ 0.63100    # Table 2 LPV IOV ka  = 94.2% (95% CI 46.5%-150.1%) folded in as BSV-equivalent (no separate BSV reported); omega^2 = log(1 + 0.942^2) = 0.63100
    etalfdepot ~ 0.04692    # Table 2 LPV IOV F   = 21.9% (95% CI 17.1%-24.0%) folded in as BSV-equivalent (no separate BSV reported); omega^2 = log(1 + 0.219^2) = 0.04692

    etalcl_rtv     ~ 0.04525  # Table 2 RTV IIV CL/F = 21.5% (95% CI 11.5%-31.7%); IOV CL/F 20.4% dropped per convention; omega^2 = log(1 + 0.215^2) = 0.04525
    etalvc_rtv     ~ 0.01035  # Table 2 RTV IIV Vc/F = 10.2% (95% CI 9.85%-10.5%); omega^2 = log(1 + 0.102^2) = 0.01035
    etalfdepot_rtv ~ 0.08731  # Table 2 RTV IIV F   = 30.3% (95% CI 17.4%-49.6%); IOV F 30.3% dropped per convention; omega^2 = log(1 + 0.303^2) = 0.08731
    etalmtt_rtv    ~ 0.07480  # Table 2 RTV IOV MTT = 27.9% (95% CI 19.6%-38.2%) folded in as BSV-equivalent (no separate BSV reported); omega^2 = log(1 + 0.279^2) = 0.07480

    # ==========================================================================
    # Residual error (Zhang 2012 Table 2 rows 'Residual variability
    # (proportional %)'). Proportional-only structure for both LPV and RTV;
    # the IOV on LPV residual variability (17.1% CV) is dropped with a note
    # in the vignette Assumptions and deviations section.
    # ==========================================================================
    propSd     <- 0.127;  label("Lopinavir proportional residual error (fraction)")  # Table 2 LPV Residual variability proportional = 12.7% (95% CI 11.6%-13.6%)
    propSd_rtv <- 0.188;  label("Ritonavir proportional residual error (fraction)")  # Table 2 RTV Residual variability proportional = 18.8% (95% CI 17.1%-20.3%)
  })

  model({
    # ------------------------------------------------------------------------
    # 1. Time-of-day diurnal indicator. Simulation convention: t is in hours,
    #    anchored to clock midnight (t = 0 -> 00:00). The morning dose is
    #    administered at clock 08:00 (t = 8 modulo 24) and the evening dose
    #    at clock 20:00 (t = 20 modulo 24). The overnight window covers
    #    clock 20:00-08:00 (= hod in [20, 24) U [0, 8)).
    #    Same convention as Bienczak_2016_nevirapine.R; documented in the
    #    vignette Assumptions and deviations section.
    # ------------------------------------------------------------------------
    hod       <- t - floor(t / 24) * 24
    overnight <- (hod >= 20) + (hod < 8)

    # ------------------------------------------------------------------------
    # 2. Allometric scaling factors (Holford / Anderson convention; reference
    #    WT = 70 kg, reference FFM = 50 kg). The reference values are not
    #    explicitly stated by Zhang 2012; the 70 / 50 kg defaults follow the
    #    most common Holford-school convention.
    # ------------------------------------------------------------------------
    all_cl <- (FFM / 50)^e_ffm_cl
    all_vc <- (WT  / 70)^e_wt_vc
    all_vp <- (WT  / 70)^e_wt_vp

    # ------------------------------------------------------------------------
    # 3. Ritonavir individual PK parameters. Computed first so that the RTV
    #    plasma concentration C_RTV is available for the dynamic LPV CL/F
    #    inhibition term below.
    # ------------------------------------------------------------------------
    ka_rtv  <- exp(lka_rtv)
    mtt_rtv <- exp(lmtt_rtv + etalmtt_rtv)
    nn_rtv  <- exp(lnn_rtv)
    vc_rtv  <- exp(lvc_rtv + etalvc_rtv) * all_vc
    vp_rtv  <- exp(lvp_rtv)              * all_vp
    q_rtv   <- exp(lq_rtv)               * all_cl
    cl_rtv  <- exp(lcl_rtv + etalcl_rtv) * all_cl *
               (1 + e_rif_cl_rtv * CONMED_RIF) *
               (1 + e_overnight_cl * overnight)

    # ------------------------------------------------------------------------
    # 4. Ritonavir bioavailability. Composition (gated by the binary
    #    CONMED_RIF indicator on the rifampicin co-administration arm):
    #      Without RIF: F = 1 (reference).
    #      With RIF:    F = (1 + e_rif_f_rtv) * (1 + e_dose_f_rtv * (DOSE-100))
    #                     = 0.55 at DOSE = 100 mg, 0.773 at 150 mg, 0.996 at 200 mg
    #                     (Zhang 2012 Results 'Model description' paragraph 5).
    #    On top of this, the evening-dose window applies a +45% multiplier
    #    and the IIV / fixed F-anchor enters via exp(lfdepot_rtv +
    #    etalfdepot_rtv).
    # ------------------------------------------------------------------------
    fdepot_rtv <- ((1 - CONMED_RIF) +
                   CONMED_RIF * (1 + e_rif_f_rtv) * (1 + e_dose_f_rtv * (DOSE - 100))) *
                  (1 + e_eve_f_rtv * overnight) *
                  exp(lfdepot_rtv + etalfdepot_rtv)

    # ------------------------------------------------------------------------
    # 5. Dynamic sigmoidal Imax inhibition of LPV CL/F by ritonavir plasma
    #    concentration (Zhang 2012 Methods 'Population pharmacokinetic
    #    analysis' equation):
    #      CL/F_LPV(t) = CL0/F_LPV * (1 - Imax * C_RTV / (IC50 + C_RTV))
    # ------------------------------------------------------------------------
    crtv  <- central_rtv / vc_rtv
    inhib <- imax * crtv / (ic50 + crtv)

    # ------------------------------------------------------------------------
    # 6. Lopinavir individual PK parameters.
    # ------------------------------------------------------------------------
    ka_lpv <- exp(lka + etalka)
    vc_lpv <- exp(lvc + etalvc) * all_vc
    cl_lpv <- exp(lcl + etalcl) * all_cl *
              (1 + e_rif_cl * CONMED_RIF) *
              (1 - inhib) *
              (1 + e_overnight_cl * overnight)

    # ------------------------------------------------------------------------
    # 7. Lopinavir bioavailability. Composition: CONMED_RIF reduces F by 20%,
    #    evening-window dose adds +42%, and the fixed F = 1 anchor plus IIV
    #    enter via exp(lfdepot + etalfdepot).
    # ------------------------------------------------------------------------
    fdepot_lpv <- (1 + e_rif_f * CONMED_RIF) *
                  (1 + e_eve_f * overnight) *
                  exp(lfdepot + etalfdepot)

    # ------------------------------------------------------------------------
    # 8. ODE system. LPV occupies a single depot + central pair with classical
    #    first-order absorption. RTV uses Savic's analytical transit-chain
    #    input via rxode2's transit(nn, mtt, bio) built-in (bio = fdepot_rtv),
    #    with the bolus content into depot_rtv suppressed (f(depot_rtv) <- 0)
    #    so the transit chain alone delivers the dose smoothly. RTV central
    #    has two-compartment disposition with linear elimination.
    # ------------------------------------------------------------------------
    d/dt(depot)            <- -ka_lpv * depot
    d/dt(central)          <-  ka_lpv * depot - cl_lpv / vc_lpv * central

    d/dt(depot_rtv)        <-  transit(nn_rtv, mtt_rtv, fdepot_rtv) - ka_rtv * depot_rtv
    d/dt(central_rtv)      <-  ka_rtv * depot_rtv -
                               cl_rtv / vc_rtv * central_rtv -
                               q_rtv  / vc_rtv * central_rtv +
                               q_rtv  / vp_rtv * peripheral1_rtv
    d/dt(peripheral1_rtv)  <-  q_rtv  / vc_rtv * central_rtv -
                               q_rtv  / vp_rtv * peripheral1_rtv

    # ------------------------------------------------------------------------
    # 9. Bioavailability hooks. LPV uses the standard f(depot) <- fdepot_lpv
    #    bolus-scaling hook; RTV suppresses the bolus contribution so the
    #    transit() chain alone provides the input.
    # ------------------------------------------------------------------------
    f(depot)     <- fdepot_lpv
    f(depot_rtv) <- 0

    # ------------------------------------------------------------------------
    # 10. Observation variables and residual error. Cc = LPV plasma
    #     concentration (mg/L); Cc_rtv = RTV plasma concentration (mg/L).
    # ------------------------------------------------------------------------
    Cc     <- central     / vc_lpv
    Cc_rtv <- central_rtv / vc_rtv

    Cc     ~ prop(propSd)
    Cc_rtv ~ prop(propSd_rtv)
  })
}
