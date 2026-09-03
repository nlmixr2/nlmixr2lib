Snelder_2013_cardiovascular_rat <- function() {
  description <- "QSP (rat). Drug-independent cardiovascular systems (CVS) model linking mean arterial blood pressure (MAP), cardiac output (CO) and total peripheral resistance (TPR) in conscious male spontaneously hypertensive rats (SHR), developed by Snelder 2013 from the simultaneous analysis of six antihypertensive compounds with diverse mechanisms of action (enalapril, fasudil, amlodipine, prazosin, propranolol, hydrochlorothiazide). CO and TPR are two coupled indirect-response (turnover) states whose zero-order production rates are each suppressed by a linear negative feedback of MAP (feedback constants FB1 on CO and FB2 on TPR); MAP is the product CO * TPR plus a five-harmonic cosine circadian rhythm. Production rate constants Kin_CO and Kin_TPR are not estimated but derived at baseline steady state from BSL_MAP, BSL_TPR, the two k_out constants and the two feedback constants (Snelder 2013 Equation 3), and BSL_CO is derived as BSL_MAP / BSL_TPR. Each of the six compounds inhibits the production rate of either CO (propranolol, hydrochlorothiazide) or TPR (enalapril, fasudil, amlodipine, prazosin) through its own Emax function of plasma concentration; all twelve drug-specific Emax / IC50 values of Snelder 2013 Table 6 are carried in the file. The model has NO internal PK: the six plasma concentrations enter as time-varying CP_<drug>_NGML covariate columns, because Snelder 2013 took every concentration-time profile from separate literature PK models and published no PK parameter values. The five system parameters are drug-independent -- omitting any one compound from the fit left them unchanged except FB1, which depends on the amlodipine data."

  reference <- paste(
    "Snelder N, Ploeger BA, Luttringer O, Rigel DF, Webb RL, Feldman D,",
    "Fu F, Beil M, Jin L, Stanski DR, Danhof M.",
    "PKPD modelling of the interrelationship between mean arterial BP,",
    "cardiac output and total peripheral resistance in conscious rats.",
    "Br J Pharmacol. 2013 Aug;169(7):1510-1524.",
    "doi:10.1111/bph.12190. PMCID:PMC3724108.",
    "A later re-use of the Snelder CVS framework (with the heart-rate /",
    "stroke-volume split of Snelder 2014 and its own parameter set) is",
    "packaged as modellib('Fu_2023_cardiovascular_qsp').",
    sep = " "
  )

  vignette <- "Snelder_2013_cardiovascular_rat"

  # Hemodynamic turnover states of the CVS system (Snelder 2013 Equation 1:
  # dCO/dt and dTPR/dt). Declared here so checkModelConventions() accepts them
  # as legitimate paper-mechanistic compartments, following the precedent set
  # by Fu_2023_cardiovascular_qsp.R (which declares "hr", "svt", "tpr" for the
  # Snelder 2014 extension of this same framework).
  paper_specific_compartments <- c("co", "tpr")

  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL",
    dosing_notes  = "This model has NO dosing compartment and no internal PK. Snelder 2013 derived every plasma-concentration-time profile from separate literature PK models (Table 3: Lin 1988 enalapril, Ikegaki 2001 fasudil, Stopher 1988 amlodipine, Hamilton 1985 prazosin, van Steeg 2010 propranolol, Asdaq & Inamdar 2009 hydrochlorothiazide) and published no PK parameter values for any of them. Drug exposure therefore enters only through the six time-varying CP_<drug>_NGML covariate columns, in ng/mL to match the IC50 units of Snelder 2013 Table 6. The `dosing` entry is nominal and exists only to satisfy the metadata schema; supply concentrations, not doses. The doses actually administered in the two source studies are recorded in population$dose_range."
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Both states are hemodynamic quantities rather than
  # drug amounts in a specimen, so specimen is "not applicable".
  compartmentData <- list(
    co  = list(analyte = "cardiac output (CO)",              units = "mL/min",           specimen = "not applicable", verified = TRUE),
    tpr = list(analyte = "total peripheral resistance (TPR)", units = "mmHg/(mL/min)",   specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    CP_ENALAPRIL_NGML = list(
      description        = "Instantaneous enalapril plasma concentration driving the Emax inhibition of TPR production.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Set to 0 outside the enalapril exposure window (and for every rat not receiving enalapril). Snelder 2013 Table 3 obtained the enalapril profile from a two-compartment model with Michaelis-Menten elimination refitted in NONMEM to data read off Lin et al. 1988 (Sprague-Dawley rats); no PK parameter values are published in Snelder 2013, so this column must be supplied externally. Study 1 dosed 30 mg/kg p.o. once daily for 6 days. Scale reference: the fitted IC50 is 2410 ng/mL (Table 6).",
      source_name        = "C (plasma concentration, Equation 6)"
    ),
    CP_FASUDIL_NGML = list(
      description        = "Instantaneous fasudil plasma concentration driving the Emax inhibition of TPR production.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Set to 0 outside the fasudil exposure window. Snelder 2013 Table 3 used a one-compartment model whose Ka and lag time were back-derived in Berkeley Madonna from the half-life, AUC and Cmax reported by Ikegaki et al. 2001 (Wistar-Kyoto rats); no PK parameter values are published in Snelder 2013. Study 1 dosed 30 mg/kg p.o. once daily for 6 days. Scale reference: the fitted IC50 is 321 ng/mL (Table 6).",
      source_name        = "C (plasma concentration, Equation 6)"
    ),
    CP_AMLODIPINE_NGML = list(
      description        = "Instantaneous amlodipine plasma concentration driving the Emax inhibition of TPR production.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Set to 0 outside the amlodipine exposure window. Snelder 2013 Table 3 used a one-compartment model whose Ka was back-derived in Berkeley Madonna from the half-life, Vd, F and Tmax reported by Stopher et al. 1988 (Sprague-Dawley rats); no PK parameter values are published in Snelder 2013. Study 1 dosed 10 mg/kg p.o. once daily for 6 days; Study 2 dosed 0.3, 1, 3 and 10 mg/kg p.o. on separate days. Scale reference: the fitted IC50 is 185 ng/mL (Table 6). Amlodipine is the one compound whose omission changes a system parameter (FB1 moves from 0.00379 to 0.00454 1/mmHg; Snelder 2013 Figure 5).",
      source_name        = "C (plasma concentration, Equation 6)"
    ),
    CP_PRAZOSIN_NGML = list(
      description        = "Instantaneous prazosin plasma concentration driving the Emax inhibition of TPR production.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Set to 0 outside the prazosin exposure window. Snelder 2013 Table 3 used a one-compartment model whose CL and Vd were allometrically scaled to the rat from the New Zealand white rabbit model of Hamilton et al. 1985, with Ka estimated simultaneously with the PD; no PK parameter values are published in Snelder 2013. Study 2 dosed 0.04, 0.2, 1 and 5 mg/kg p.o. on separate days. Scale reference: the fitted IC50 is 0.133 ng/mL, the single imprecise drug parameter in the paper (CV 110%, 95% CI -0.15 to 0.4 ng/mL; Table 6).",
      source_name        = "C (plasma concentration, Equation 6)"
    ),
    CP_PROPRANOLOL_NGML = list(
      description        = "Instantaneous propranolol plasma concentration driving the Emax inhibition of CO production.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Set to 0 outside the propranolol exposure window. Snelder 2013 Table 3 used the three-compartment model of van Steeg et al. 2010 (Wistar-Kyoto rats) with absorption described as an infusion of fixed 24 h duration and Ka estimated; no PK parameter values are published in Snelder 2013. Study 1 supplied propranolol continuously at 1 mg/mL in drinking water, hence the 24 h infusion representation. Scale reference: the fitted IC50 is 9.82 ng/mL (Table 6). Propranolol is one of only two compounds acting primarily on CO.",
      source_name        = "C (plasma concentration, Equation 6)"
    ),
    CP_HCTZ_NGML = list(
      description        = "Instantaneous hydrochlorothiazide (HCTZ) plasma concentration driving the Emax inhibition of CO production.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Set to 0 outside the hydrochlorothiazide exposure window. Snelder 2013 Table 3 used a one-compartment model taking Ke, Ka and Vd from Asdaq & Inamdar 2009 (Wistar-Kyoto rats) with F back-calculated from the reported AUC; no PK parameter values are published in Snelder 2013. Study 2 dosed 0.1, 0.3, 1 and 3 mg/kg p.o. on separate days. Scale reference: the fitted IC50 is 12300 ng/mL (Table 6). Snelder 2013 notes that the 1 mg/kg HCTZ response was under-predicted, which they attribute to the Emax = 1 constraint rather than to the system model.",
      source_name        = "C (plasma concentration, Equation 6)"
    )
  )

  population <- list(
    species        = "rat (male spontaneously hypertensive rat, SHR; Taconic Farms)",
    n_subjects     = 12L,
    n_studies      = 2L,
    study_names    = c(
      "Study 1 (multiple-dose telemetry; MAP only; enalapril, fasudil, amlodipine or propranolol once daily for 6 days, n = 5 SHR per drug drawn from 10 animals)",
      "Study 2 (single ascending doses on separate days; MAP and CO, with TPR derived; amlodipine, prazosin or hydrochlorothiazide, n = 2 SHR)"
    ),
    age_range      = "21-45 weeks at time of study.",
    weight_range   = "269 to 490 g.",
    sex_female_pct = 0,
    sex_notes      = "All animals were male (Methods, 'Animals').",
    disease_state  = "Spontaneous (genetic) hypertension; the SHR strain is hypertensive at baseline (BSL_MAP = 147 mmHg) with no induced disease model. Snelder 2013 states explicitly (Discussion) that the identified system-parameter set is specific to the SHR strain and that applications using these values are limited to that strain.",
    dose_range     = "Study 1: enalapril 30 mg/kg p.o. QD x 6, fasudil 30 mg/kg p.o. QD x 6, amlodipine 10 mg/kg p.o. QD x 6, propranolol 1 mg/mL in drinking water. Study 2: amlodipine 0.3, 1, 3, 10 mg/kg p.o.; prazosin 0.04, 0.2, 1, 5 mg/kg p.o.; hydrochlorothiazide 0.1, 0.3, 1, 3 mg/kg p.o., one dose per day on separate days after a vehicle day.",
    regions        = "Preclinical; Novartis Institutes for BioMedical Research, East Hanover, NJ, USA (in-life) with modelling at Leiden Academic Centre for Drug Research / LAP&P Consultants BV, The Netherlands.",
    instrumentation = "Study 1: femoral-artery catheter / radiotransmitter (DSI PA-C40), BP recorded 15 s every 10 min. Study 2: ascending-aortic transit-time flow probe (Transonic 2.5PS or 3PS) plus femoral catheter / radiotransmitter; MAP, heart rate, stroke volume, CO and TPR derived for all beats and averaged over consecutive 2 min intervals. Study 2 animals were reused for up to 6 months with washout between experiments.",
    n_ode_states   = 2L,
    notes          = "The total number of animals was 12 (10 in Study 1, 2 in Study 2); the per-drug n = 5 of Table 2 therefore reflects reuse of the Study 1 animals across compounds. TPR was never measured directly: in the experiments it was derived from measured MAP and CO, and in the modelling BSL_MAP and BSL_TPR were the estimated parameters with BSL_CO derived from them 'for reasons of model stability' (Data analysis). CO was measured only in Study 2, and among the two CO-acting compounds only hydrochlorothiazide had CO measured -- which is why kout_TPR and FB2 remain strongly correlated (-0.984) in the final fit."
  )

  ini({
    # ==================================================================
    # SYSTEM-SPECIFIC PARAMETERS -- Snelder 2013 Table 5.
    #
    # These five (plus the fixed kout_CO) are the paper's primary
    # result: the drug-independent description of the SHR
    # cardiovascular system. Omitting any single compound from the
    # fit left all of them unchanged except FB1 (Figure 5).
    # ==================================================================
    lrbase_MAP <- log(147)    ; label("Baseline mean arterial pressure BSL_MAP (mmHg)")                          # Snelder 2013 Table 5 (BSL_MAP = 147 mmHg, SE 1.38, CV 0.939%, 95% CI 144-150)
    lrbase_TPR <- log(2.32)   ; label("Baseline total peripheral resistance BSL_TPR (mmHg/(mL/min))")            # Snelder 2013 Table 5 (BSL_TPR = 2.32 mmHg*mL^-1*min, SE 0.132, CV 5.69%, 95% CI 2.06-2.58)
    # BSL_CO is NOT a parameter: the paper derives it as BSL_MAP / BSL_TPR
    # ("in the modelling, the baseline values of MAP (BSL_MAP) and BSL_TPR
    # were estimated and BSL_CO was derived from these parameters for
    # reasons of model stability", Data analysis). 147 / 2.32 = 63.4 mL/min.

    lkout_CO   <- fixed(log(99))    ; label("First-order dissipation rate constant of CO, kout_CO (1/h)")        # Snelder 2013 Table 5 (kout_CO = 99 1/h FIXED). Results: "The dissipation rate of CO (kout_CO) was found to be very high and could not be estimated accurately. Therefore, this parameter was fixed to a high value (99 1/h) prior the other model parameters were estimated."
    lkout_TPR  <- log(0.260)        ; label("First-order dissipation rate constant of TPR, kout_TPR (1/h)")      # Snelder 2013 Table 5 (kout_TPR = 0.260 1/h, SE 0.129, CV 49.6%, 95% CI 0.00716-0.513)

    # Linear negative feedback of MAP on the two production rates. Named
    # lfb_CO / lfb_TPR after the affected state, extending the single
    # `lfb` of Fu_2023_cardiovascular_qsp.R (which has one shared
    # feedback constant) with the same _<STATE> suffix that file already
    # uses for lkout_<STATE> and lrbase_<STATE>.
    lfb_CO     <- log(0.00378)      ; label("Magnitude of the negative feedback of MAP on CO production, FB1 (1/mmHg)")   # Snelder 2013 Table 5 (FB1 = 0.00378 1/mmHg, SE 0.000148, CV 3.92%, 95% CI 0.00349-0.00407). The Results text quotes 0.00379 for the same estimate when contrasting it with the amlodipine-omitted value of 0.00454; Table 5 is the final-model number and is used here.
    lfb_TPR    <- log(0.00492)      ; label("Magnitude of the negative feedback of MAP on TPR production, FB2 (1/mmHg)")  # Snelder 2013 Table 5 (FB2 = 0.00492 1/mmHg, SE 0.00101, CV 20.5%, 95% CI 0.00294-0.00690)

    # ==================================================================
    # CIRCADIAN RHYTHM ON MAP -- Snelder 2013 Equation 2 and Results.
    #
    # Equation 2 writes the rhythm as the sum of up to 10 cosine
    # harmonics of a 24 h fundamental, amplitude amp_n (mmHg) and a
    # single shared horizontal displacement HOR (h):
    #   MAP = CO*TPR + sum_{n=1}^{10} amp_n * cos(n*2*pi*(t + HOR)/24)
    # "The number of cosine functions was systematically reduced
    # following the criteria for statistical significance": five
    # harmonics survived, the other five were fixed to zero. All are
    # carried here so the reduction is explicit and auditable.
    #
    # Kept on the linear scale (not log) because amp_3 is negative.
    # ==================================================================
    amp_MAP1   <- 3.17          ; label("Amplitude of the 1st circadian harmonic on MAP, amp_1 (mmHg; period 24 h)")     # Snelder 2013 Results, Model development ("Amp1, amp3, amp4, amp5 and amp7 were estimated to be 3.17, -2.03, 1.15, 1.63 and 1.28 mmHg, respectively")
    amp_MAP2   <- fixed(0)      ; label("Amplitude of the 2nd circadian harmonic on MAP, amp_2 (mmHg; period 12 h)")     # Snelder 2013 Results ("Amp2, amp6, amp8, amp9 and amp10 were fixed to 0, implying that these harmonics did not contribute to the circadian rhythm")
    amp_MAP3   <- -2.03         ; label("Amplitude of the 3rd circadian harmonic on MAP, amp_3 (mmHg; period 8 h)")      # Snelder 2013 Results, Model development (amp3 = -2.03 mmHg)
    amp_MAP4   <- 1.15          ; label("Amplitude of the 4th circadian harmonic on MAP, amp_4 (mmHg; period 6 h)")      # Snelder 2013 Results, Model development (amp4 = 1.15 mmHg)
    amp_MAP5   <- 1.63          ; label("Amplitude of the 5th circadian harmonic on MAP, amp_5 (mmHg; period 4.8 h)")    # Snelder 2013 Results, Model development (amp5 = 1.63 mmHg)
    amp_MAP6   <- fixed(0)      ; label("Amplitude of the 6th circadian harmonic on MAP, amp_6 (mmHg; period 4 h)")      # Snelder 2013 Results (fixed to 0)
    amp_MAP7   <- 1.28          ; label("Amplitude of the 7th circadian harmonic on MAP, amp_7 (mmHg; period 24/7 h)")   # Snelder 2013 Results, Model development (amp7 = 1.28 mmHg)
    amp_MAP8   <- fixed(0)      ; label("Amplitude of the 8th circadian harmonic on MAP, amp_8 (mmHg; period 3 h)")      # Snelder 2013 Results (fixed to 0)
    amp_MAP9   <- fixed(0)      ; label("Amplitude of the 9th circadian harmonic on MAP, amp_9 (mmHg; period 24/9 h)")   # Snelder 2013 Results (fixed to 0)
    amp_MAP10  <- fixed(0)      ; label("Amplitude of the 10th circadian harmonic on MAP, amp_10 (mmHg; period 2.4 h)")  # Snelder 2013 Results (fixed to 0)

    # NOT REPORTED IN THE SOURCE. Snelder 2013 defines HOR in Equation 2
    # and in the Abbreviations list but never publishes its estimate:
    # it is absent from Table 5 (which carries the other five system
    # parameters), from the Results narrative that lists the five
    # amplitudes, and from every figure panel (Figures 3, 4, 5, 6, A1
    # were inspected at 200-600 dpi). Anchored at 0 h so the harmonic
    # sum peaks at t = 0 of the model's own time axis. The circadian
    # SHAPE and AMPLITUDES are the paper's; only the PHASE relative to
    # study clock time is arbitrary. Set hor_MAP to align the rhythm
    # with your own light/dark cycle. See the vignette's "Assumptions
    # and deviations" section.
    hor_MAP    <- fixed(0)      ; label("Horizontal displacement of the MAP circadian rhythm, HOR (h) -- NOT REPORTED, anchored at 0")  # Snelder 2013 Equation 2 defines HOR; no estimate is published anywhere in the paper

    # ==================================================================
    # DRUG-SPECIFIC PARAMETERS -- Snelder 2013 Table 6.
    #
    # Every compound acts through the same Emax function of plasma
    # concentration (Equation 6) on the production rate of its primary
    # hemodynamic state (Table 1, "Primary effect" column), and the
    # effect is SUBTRACTED inside the production term of Equation 4, so
    # a positive Emax is an inhibition. Each compound carries a stratum
    # suffix; none keeps a bare canonical name.
    #
    # Assumption 3 of Table 4: for amlodipine, fasudil, enalapril and
    # hydrochlorothiazide the maximum effect was never observed, so
    # Emax was FIXED to 1 (complete inhibition of Kin at infinite
    # concentration) to make IC50 identifiable. Fixing amlodipine's
    # Emax to 0.8 instead of 1 moved no system parameter materially
    # (Table 5, last column).
    #
    # Table 6 labels the potency parameter IC50 while Equation 6 writes
    # EC50; the file follows Table 6 (the effect is an inhibition).
    # ------------------------------------------------------------------
    # Compounds acting on TPR (Table 1)
    # ------------------------------------------------------------------
    emax_enalapril  <- fixed(1)       ; label("Maximum fractional inhibition of TPR production by enalapril, Emax []")        # Snelder 2013 Table 6 (Emax = 1 fixed; Table 4 Assumption 3)
    lic50_enalapril <- log(2410)      ; label("Enalapril concentration giving half-maximal inhibition of TPR production, IC50 (ng/mL)")   # Snelder 2013 Table 6 (IC50 = 2410 ng/mL, SE 373, CV 15.5%, 95% CI 1679-3141)
    emax_fasudil    <- fixed(1)       ; label("Maximum fractional inhibition of TPR production by fasudil, Emax []")          # Snelder 2013 Table 6 (Emax = 1 fixed; Table 4 Assumption 3)
    lic50_fasudil   <- log(321)       ; label("Fasudil concentration giving half-maximal inhibition of TPR production, IC50 (ng/mL)")     # Snelder 2013 Table 6 (IC50 = 321 ng/mL, SE 60.3, CV 18.8%, 95% CI 203-439)
    emax_amlodipine <- fixed(1)       ; label("Maximum fractional inhibition of TPR production by amlodipine, Emax []")       # Snelder 2013 Table 6 (Emax = 1 fixed; Table 4 Assumption 3). Table 5 also reports the whole system-parameter set re-estimated with this value fixed to 0.8 instead, as a sensitivity analysis.
    lic50_amlodipine <- log(185)      ; label("Amlodipine concentration giving half-maximal inhibition of TPR production, IC50 (ng/mL)")  # Snelder 2013 Table 6 (IC50 = 185 ng/mL, SE 26.2, CV 14.2%, 95% CI 134-236)
    emax_prazosin   <- 0.213          ; label("Maximum fractional inhibition of TPR production by prazosin, Emax []")         # Snelder 2013 Table 6 (Emax = 0.213, SE 0.0158, CV 7.42%, 95% CI 0.182-0.244); estimated, not fixed
    lic50_prazosin  <- log(0.133)     ; label("Prazosin concentration giving half-maximal inhibition of TPR production, IC50 (ng/mL)")    # Snelder 2013 Table 6 (IC50 = 0.133 ng/mL, SE 0.146, CV 109.8%, 95% CI -0.15 to 0.4). The single imprecise drug parameter in the paper: Results notes that fixing Emax to 1 sharpened the IC50 but significantly worsened the objective function, so the imprecise joint estimate was retained.
    # ------------------------------------------------------------------
    # Compounds acting on CO (Table 1)
    # ------------------------------------------------------------------
    emax_propranolol  <- 0.335        ; label("Maximum fractional inhibition of CO production by propranolol, Emax []")       # Snelder 2013 Table 6 (Emax = 0.335, SE 0.0624, CV 18.6%, 95% CI 0.213-0.457); estimated, not fixed
    lic50_propranolol <- log(9.82)    ; label("Propranolol concentration giving half-maximal inhibition of CO production, IC50 (ng/mL)")  # Snelder 2013 Table 6 (IC50 = 9.82 ng/mL, SE 3.8, CV 38.7%, 95% CI 2.37-17.3)
    emax_hctz         <- fixed(1)     ; label("Maximum fractional inhibition of CO production by hydrochlorothiazide, Emax []")           # Snelder 2013 Table 6 (Emax = 1 fixed; Table 4 Assumption 3)
    lic50_hctz        <- log(12300)   ; label("Hydrochlorothiazide concentration giving half-maximal inhibition of CO production, IC50 (ng/mL)")  # Snelder 2013 Table 6 (IC50 = 12 300 ng/mL, SE 780, CV 6.34%, 95% CI 10 771-13 829)

    # ==================================================================
    # INTER-INDIVIDUAL VARIABILITY -- structure reported, magnitudes NOT.
    #
    # Snelder 2013 Results: "In Study 1, BSL_MAP was allowed to vary
    # between individual rats [inter-individual variability (IIV)].
    # Study 2 provided information to estimate IIV on both BSL_MAP and
    # BSL_TPR." Computation: "Random effects were included as
    # exponential terms reflecting log normal distributions of model
    # parameters" -- i.e. log-normal IIV on BSL_MAP and BSL_TPR, exactly
    # the etalrbase_MAP / etalrbase_TPR this file would carry.
    #
    # No omega value is published: Table 5 lists only the six structural
    # system parameters and Table 6 only the drug parameters. The etas
    # are therefore OMITTED rather than written as `~ fixed(0)`, because
    # a zero-variance diagonal makes OMEGA singular and breaks the
    # Cholesky sampler used by rxSolve (the Thoueille_2026_salmeterol.R
    # precedent). The model is a typical-value mechanism as packaged;
    # add `etalrbase_MAP ~ <var>` / `etalrbase_TPR ~ <var>` and the
    # matching `+ etalrbase_*` terms in model() if you have magnitudes.
    # ==================================================================

    # ==================================================================
    # RESIDUAL ERROR -- structure reported, magnitudes NOT.
    #
    # Snelder 2013 Results: "The residual errors of MAP and TPR were
    # best described by additive residual error models, whereas the
    # residual error of CO was best described by a proportional error
    # model." The published SIGMA values appear in neither Table 5 nor
    # Table 6. The Discussion reports derived detection limits instead
    # ("the model is qualified to distinguish changes in MAP, CO and TPR
    # larger than 7.6 mmHg, 4.3 mL/min and 0.5 mmHg/(mL/min) ... from
    # noise"); those are noise thresholds on a stated but unspecified
    # scale, not the sigma estimates, and the 4.3 mL/min figure is in
    # additive units while the CO error model is proportional -- so they
    # cannot be back-converted. The reported STRUCTURE is preserved with
    # magnitudes at 0 (Machavaram 2013 / Takada 2025 precedent).
    # ==================================================================
    addSd_MAP  <- fixed(0)  ; label("Additive residual SD on MAP (mmHg; ZERO - structure reported, magnitude not published)")                  # Snelder 2013 Results, Model development (additive residual error model on MAP; no sigma published)
    propSd_CO  <- fixed(0)  ; label("Proportional residual SD on CO (fraction; ZERO - structure reported, magnitude not published)")           # Snelder 2013 Results, Model development (proportional residual error model on CO; no sigma published)
    addSd_TPR  <- fixed(0)  ; label("Additive residual SD on TPR (mmHg/(mL/min); ZERO - structure reported, magnitude not published)")         # Snelder 2013 Results, Model development (additive residual error model on TPR; no sigma published)
  })

  model({
    # ==================================================================
    # 1. Linear-scale system parameters and derived baseline.
    #
    # BSL_CO is derived, not estimated (Snelder 2013 Data analysis):
    #   BSL_MAP = BSL_CO * BSL_TPR  =>  BSL_CO = BSL_MAP / BSL_TPR
    # With the published values, 147 / 2.32 = 63.36 mL/min.
    # ==================================================================
    bsl_MAP  <- exp(lrbase_MAP)
    bsl_TPR  <- exp(lrbase_TPR)
    bsl_CO   <- bsl_MAP / bsl_TPR
    kout_CO  <- exp(lkout_CO)
    kout_TPR <- exp(lkout_TPR)
    fb1      <- exp(lfb_CO)
    fb2      <- exp(lfb_TPR)

    # ==================================================================
    # 2. Steady-state production rate constants -- Snelder 2013
    #    Equation 3, transcribed in the paper's published algebraic
    #    form rather than the simplified equivalent:
    #
    #      Kin_CO  = -kout_CO * BSL_CO / (-1 + FB1*BSL_CO*BSL_TPR)
    #      Kin_TPR = kout_TPR * (Kin_CO*FB1*BSL_TPR + kout_CO) * BSL_TPR
    #                / (Kin_CO*FB1*BSL_TPR + kout_CO - FB2*Kin_CO*BSL_TPR)
    #
    #    Both reduce to the familiar turnover form
    #    Kin_X = kout_X * BSL_X / (1 - FB_X * BSL_MAP); the vignette
    #    proves that equivalence numerically. Units: Kin_CO in
    #    (mL/min)/h, Kin_TPR in (mmHg/(mL/min))/h.
    # ==================================================================
    kin_CO  <- -kout_CO * bsl_CO / (-1 + fb1 * bsl_CO * bsl_TPR)
    kin_TPR <- kout_TPR * (kin_CO * fb1 * bsl_TPR + kout_CO) * bsl_TPR /
      (kin_CO * fb1 * bsl_TPR + kout_CO - fb2 * kin_CO * bsl_TPR)

    # ==================================================================
    # 3. Circadian rhythm on MAP -- Snelder 2013 Equation 2. Additive,
    #    in mmHg, on top of the CO * TPR product. Harmonic n has period
    #    24/n hours. Five amplitudes were estimated and five fixed to 0.
    # ==================================================================
    circ_MAP <-
      amp_MAP1  * cos(1  * 2 * pi * (t + hor_MAP) / 24) +
      amp_MAP2  * cos(2  * 2 * pi * (t + hor_MAP) / 24) +
      amp_MAP3  * cos(3  * 2 * pi * (t + hor_MAP) / 24) +
      amp_MAP4  * cos(4  * 2 * pi * (t + hor_MAP) / 24) +
      amp_MAP5  * cos(5  * 2 * pi * (t + hor_MAP) / 24) +
      amp_MAP6  * cos(6  * 2 * pi * (t + hor_MAP) / 24) +
      amp_MAP7  * cos(7  * 2 * pi * (t + hor_MAP) / 24) +
      amp_MAP8  * cos(8  * 2 * pi * (t + hor_MAP) / 24) +
      amp_MAP9  * cos(9  * 2 * pi * (t + hor_MAP) / 24) +
      amp_MAP10 * cos(10 * 2 * pi * (t + hor_MAP) / 24)

    # ==================================================================
    # 4. Drug effects -- Snelder 2013 Equation 6, one Emax term per
    #    compound, routed to the compound's primary hemodynamic state
    #    (Table 1, "Primary effect"). Table 4 Assumption 1: each
    #    compound acts on ONE of CO or TPR; any change in the other is
    #    a consequence of the feedback, not a second direct effect.
    #
    #    The terms are summed so that a dataset containing several
    #    compounds (as the paper's own simultaneous fit did) works
    #    unchanged: in practice only one CP_ column is non-zero on any
    #    given record, and every column defaults to 0.
    # ==================================================================
    eff_TPR <-
      emax_enalapril  * CP_ENALAPRIL_NGML  / (exp(lic50_enalapril)  + CP_ENALAPRIL_NGML)  +
      emax_fasudil    * CP_FASUDIL_NGML    / (exp(lic50_fasudil)    + CP_FASUDIL_NGML)    +
      emax_amlodipine * CP_AMLODIPINE_NGML / (exp(lic50_amlodipine) + CP_AMLODIPINE_NGML) +
      emax_prazosin   * CP_PRAZOSIN_NGML   / (exp(lic50_prazosin)   + CP_PRAZOSIN_NGML)

    eff_CO <-
      emax_propranolol * CP_PROPRANOLOL_NGML / (exp(lic50_propranolol) + CP_PROPRANOLOL_NGML) +
      emax_hctz        * CP_HCTZ_NGML        / (exp(lic50_hctz)        + CP_HCTZ_NGML)

    # ==================================================================
    # 5. ODE system -- Snelder 2013 Equations 1, 2 and 4.
    #
    #    MAP is the model's single blood-pressure quantity: the product
    #    CO * TPR plus the circadian term (Equation 2). It is that MAP
    #    which drives the negative feedback, per the paper's statement
    #    that "at baseline, MAP oscillates around its baseline value,
    #    which equals the product of the baseline values of CO and TPR"
    #    -- the rhythm is part of MAP, and Equation 3 sets Kin so that
    #    the MEAN of MAP is BSL_MAP.
    #
    #    Drug effect enters as a SUBTRACTION inside the production term
    #    (Equation 4), so a positive EFF inhibits production.
    # ==================================================================
    map <- co * tpr + circ_MAP

    d/dt(co)  <- kin_CO  * (1 - fb1 * map - eff_CO)  - kout_CO  * co
    d/dt(tpr) <- kin_TPR * (1 - fb2 * map - eff_TPR) - kout_TPR * tpr

    co(0)  <- bsl_CO
    tpr(0) <- bsl_TPR

    # ==================================================================
    # 6. Observations. MAP and CO were measured directly (MAP in both
    #    studies, CO in Study 2 only); TPR was derived from them
    #    experimentally but carries its own residual error model in the
    #    fit, so it is exposed as a third endpoint.
    # ==================================================================
    MAP <- map
    CO  <- co
    TPR <- tpr

    MAP ~ add(addSd_MAP)
    CO  ~ prop(propSd_CO)
    TPR ~ add(addSd_TPR)
  })
}
