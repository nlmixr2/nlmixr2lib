Snelder_2014_cardiovascular_rat <- function() {
  description <- "QSP (rat). Extended drug-independent cardiovascular systems (CVS) model linking mean arterial pressure (MAP), cardiac output (CO), heart rate (HR), stroke volume (SV) and total peripheral resistance (TPR) in conscious spontaneously hypertensive (SHR) and normotensive Wistar-Kyoto (WKY) rats. Snelder 2014 extends the Snelder 2013 two-state (CO, TPR) model by parsing CO into HR and SV, giving THREE coupled indirect-response states -- HR, an SV-turnover state SVT, and TPR -- whose zero-order production rates are each suppressed by one shared linear negative feedback of MAP. Actual stroke volume adds a direct inverse log-linear coupling to heart rate, SV = SVT * (1 - HR_SV * ln(HR / BSL_HR)), representing the shortening of left-ventricular filling time as HR rises; CO = HR * SV and MAP = CO * TPR. Two cosine circadian rhythms multiply the HR and TPR production rates, and an exponentially decaying handling effect (brief manual restraint at each oral dosing) transiently raises both. Production rate constants are derived at baseline steady state rather than estimated (Equation 5), and BSL_SV and BSL_TPR are derived from the estimated BSL_HR, BSL_MAP and BSL_CO. Baselines are strain-specific and the feedback constant declines with the individual baseline MAP through a power relationship (Equation 9), making feedback about twofold stronger in normotensive WKY than in hypertensive SHR rats. Seven of the eight challenge compounds carry drug-specific parameters, each routed to its site of action per Table 4: atropine stimulates HR (linear); amiloride and hydrochlorothiazide inhibit SV, amlodipine and fasudil inhibit TPR, and prazosin inhibits TPR through a power model; enalapril inhibits BOTH TPR and SV with one shared EC50, delayed through an effect compartment. Propranolol's effect was too small to quantify and carries no parameters. The model has NO internal PK: the seven plasma concentrations enter as time-varying CP_<drug>_NGML covariate columns, because Snelder 2014 took every concentration-time profile from a separate literature PK model (Table 2) and published no PK parameter values."

  reference <- paste(
    "Snelder N, Ploeger BA, Luttringer O, Rigel DF, Fu F, Beil M,",
    "Stanski DR, Danhof M.",
    "Drug effects on the CVS in conscious rats: separating cardiac output",
    "into heart rate and stroke volume using PKPD modelling.",
    "Br J Pharmacol. 2014 Nov;171(22):5076-5092.",
    "doi:10.1111/bph.12824. PMCID PMC4253457.",
    sep = " "
  )
  vignette <- "Snelder_2014_cardiovascular_rat"

  # Paper-mechanistic states of the extended CVS system (Snelder 2014
  # Equations 4 and 6): the three linked turnover states hr / svt / tpr,
  # plus a first-order decay state carrying the handling artefact. The
  # hr / svt / tpr names are reused verbatim from
  # Fu_2023_cardiovascular_qsp.R, which encodes this same Snelder 2014
  # system with a hypothetical drug, and follow the two-state co / tpr
  # naming of the Snelder 2013 predecessor. `handling` is the empirical
  # HD function of Equation 4 realised as an ODE state so that repeated
  # handling events superpose (Study 2 handled the rats twice, at 0 h and
  # 3 h); `effect` is the canonical effect compartment of Equation 8,
  # used only by enalapril.
  paper_specific_compartments <- c("hr", "svt", "tpr", "handling")

  units <- list(
    time          = "h",
    dosing        = "unitless impulse into the `handling` state",
    concentration = "ng/mL (all seven CP_<drug>_NGML driver columns and every EC50 / slope)",
    dosing_notes  = "This model has NO drug dosing compartment and no internal PK. Snelder 2014 Table 2 derived every plasma-concentration-time profile from a separate literature PK model (Segre 1998 amiloride, Stopher 1988 amlodipine, Perlstein 2002 atropine, Lin 1988 + Li 2007 enalapril, Ikegaki 2001 fasudil, Asdaq & Inamdar 2009 hydrochlorothiazide, Hamilton 1985 prazosin, van Steeg 2010 + Belpaire 1990 propranolol) and published no CL / V / ka / F value for any of them. Drug exposure therefore enters only through the seven time-varying CP_<drug>_NGML covariate columns, in ng/mL to match the EC50 and slope units of Table 5. The only dosing records this model consumes are unit impulses (amt = 1) into the `handling` compartment, one per manual-restraint / oral-gavage event, which drive the Equation 4 handling artefact. The doses actually administered in the two source studies are recorded in population$dose_range.",
    output_units  = "HR in beats/min, SV in mL/beat, CO in mL/min, TPR in mmHg/(mL/min), MAP in mmHg."
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Four of the five states are hemodynamic quantities
  # or an empirical artefact rather than drug amounts in a specimen, so
  # specimen is "not applicable" for those.
  compartmentData <- list(
    hr       = list(analyte = "heart rate (HR)",                                          units = "beats/min",       specimen = "not applicable", verified = TRUE),
    svt      = list(analyte = "stroke-volume turnover state (SV*, the SV driven by MAP feedback before the direct HR coupling)", units = "mL/beat", specimen = "not applicable", verified = TRUE),
    tpr      = list(analyte = "total peripheral resistance (TPR)",                        units = "mmHg/(mL/min)",   specimen = "not applicable", verified = TRUE),
    effect   = list(analyte = "enalapril effect-site concentration (Ce, Equation 8)",     units = "ng/mL",           specimen = "not applicable", verified = TRUE),
    handling = list(analyte = "handling-artefact decay state (unit impulse per restraint / gavage event, Equation 4)", units = "unitless", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    STRAIN_SHR = list(
      description        = "1 = spontaneously hypertensive rat (SHR); 0 = normotensive Wistar-Kyoto (WKY) rat.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (WKY, normotensive)",
      notes              = "Time-fixed per animal. Selects the strain-specific baseline triple BSL_HR / BSL_MAP / BSL_CO of Snelder 2014 Table 5; every other system parameter is shared across the two strains. The strain does NOT enter the feedback constant directly -- FB is a power function of the individual baseline MAP (Equation 9), so the twofold-stronger feedback the paper reports in WKY rats emerges from their lower BSL_MAP rather than from a strain term. Study 2 (atropine / propranolol) enrolled SHR only; among the Study 1 compounds only amlodipine, hydrochlorothiazide and prazosin were given to WKY rats.",
      source_name        = "strain (Methods, 'Animals'; Table 1 'Strain' column; the _SHR / _WKY suffixes on the Table 5 baseline rows)"
    ),
    CP_AMILORIDE_NGML = list(
      description        = "Instantaneous amiloride plasma concentration driving the Emax inhibition of SV production.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Set to 0 outside the amiloride exposure window and for every rat not receiving amiloride. Snelder 2014 Table 2 used the two-compartment model with a liver compartment of Segre et al. 1998 (Wistar rats); no PK parameter values are published in Snelder 2014, so this column must be supplied externally. Study 1 dosed 10 mg/kg p.o. in 3 SHR. Scale reference: the fitted EC50 is 245 ng/mL (Table 5), the least precise EC50 in the paper (RSE 25.1%, 95% CI 125-365 ng/mL).",
      source_name        = "C (plasma concentration, Equations 6 and 7)"
    ),
    CP_AMLODIPINE_NGML = list(
      description        = "Instantaneous amlodipine plasma concentration driving the Emax inhibition of TPR production.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Set to 0 outside the amlodipine exposure window. Snelder 2014 Table 2 used a one-compartment model whose Ka was back-derived in Berkeley Madonna from the half-life, Vd, F and Tmax reported by Stopher et al. 1988 (Sprague-Dawley rats); no PK parameter values are published in Snelder 2014. Study 1 dosed 0.3, 1, 3 and 10 mg/kg p.o. on separate days in 2 SHR and 2 WKY rats. Scale reference: the fitted EC50 is 82.8 ng/mL (Table 5). Amlodipine is the paradigm compound of the HR-and-MAP-only analysis: re-fitting it to HR and MAP alone with the system parameters fixed recovered EC50 = 84.9 ng/mL (95% CI 75.4-94.4), statistically indistinguishable from the full-data estimate.",
      source_name        = "C (plasma concentration, Equations 6 and 7)"
    ),
    CP_ATROPINE_NGML = list(
      description        = "Instantaneous atropine plasma concentration driving the linear stimulation of HR production.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Set to 0 outside the atropine exposure window. Snelder 2014 Table 2 used the two-compartment model of Perlstein et al. 2002 (Sabra rats), given i.v. there but p.o. here, so the absorption rate constant was estimated simultaneously with the PD (Ka = 1.17 1/h, RSE 59.9%; Table 5). The distribution and elimination parameters are not published in Snelder 2014, so the profile must be supplied externally -- Ka alone does not reconstruct it. Study 2 dosed 10 mg/kg p.o., alone or 3 h before / after propranolol, in 8 SHR. Scale reference: the fitted linear slope is 0.00149 per ng/mL (Table 5); atropine is the only compound in the paper with a stimulating effect.",
      source_name        = "C (plasma concentration, Equations 6 and 7)"
    ),
    CP_ENALAPRIL_NGML = list(
      description        = "Instantaneous enalapril plasma concentration; equilibrates into the effect compartment that drives the Emax inhibition of BOTH TPR and SV production.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Set to 0 outside the enalapril exposure window. Snelder 2014 Table 2 refitted a two-compartment model with Michaelis-Menten elimination in NONMEM to data read off Lin et al. 1988 and Li et al. 2007 (Sprague-Dawley rats); no PK parameter values are published in Snelder 2014. Study 1 dosed 3, 10 and 30 mg/kg p.o. on separate days in 4 SHR. Enalapril is the only compound whose effect is delayed through an effect compartment (ke0 = 0.163 1/h, half-life 4.3 h; Equation 8), and the only one acting on two sites at once, with a single shared EC50 = 1200 ng/mL because the separately estimated TPR and SV values had overlapping confidence intervals.",
      source_name        = "C (plasma concentration, Equations 6, 7 and 8)"
    ),
    CP_FASUDIL_NGML = list(
      description        = "Instantaneous fasudil plasma concentration driving the Emax inhibition of TPR production.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Set to 0 outside the fasudil exposure window. Snelder 2014 Table 2 used a one-compartment model whose Ka and lag time were back-derived in Berkeley Madonna from the half-life, AUC and Cmax reported by Ikegaki et al. 2001 (Wistar-Kyoto rats); no PK parameter values are published in Snelder 2014. Study 1 dosed 3, 10 and 30 mg/kg p.o. on separate days in 4 SHR. Scale reference: the fitted EC50 is 0.172 ng/mL (Table 5) -- roughly 1900-fold more potent than the 321 ng/mL reported for fasudil by the Snelder 2013 predecessor on the same experimental platform. Both values are internally consistent with their own reported RSE and confidence interval, so this is a genuine discrepancy between the two publications and not a transcription artefact; it is discussed in the vignette's Errata. Because the driver is a covariate column supplied by the user, the practical consequence is that a fasudil profile built for one paper's EC50 scale must not be reused with the other's.",
      source_name        = "C (plasma concentration, Equations 6 and 7)"
    ),
    CP_HCTZ_NGML = list(
      description        = "Instantaneous hydrochlorothiazide (HCTZ) plasma concentration driving the Emax inhibition of SV production.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Set to 0 outside the hydrochlorothiazide exposure window. Snelder 2014 Table 2 used a one-compartment model taking Ke, Ka and Vd from Asdaq & Inamdar 2009 (Wistar-Kyoto rats) with AUC/F calculated from them; no PK parameter values are published in Snelder 2014. Study 1 dosed 0.1, 0.3, 1 and 3 mg/kg p.o. on the first occasion (2 SHR, 2 WKY) and 10 and 30 mg/kg on a second occasion (4 SHR), the higher doses added precisely because the Snelder 2013 dose range had not reached the maximum effect. Scale reference: the fitted EC50 is 28900 ng/mL (Table 5).",
      source_name        = "C (plasma concentration, Equations 6 and 7)"
    ),
    CP_PRAZOSIN_NGML = list(
      description        = "Instantaneous prazosin plasma concentration driving the power-model inhibition of TPR production.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Set to 0 outside the prazosin exposure window. Snelder 2014 Table 2 scaled the CL and Vd of the New Zealand white rabbit model of Hamilton et al. 1985 allometrically to the rat; the absorption rate constant could not be estimated with precision and was FIXED to 99 1/h (Results, 'Drug effects'), which is a near-instantaneous-absorption idealisation rather than a measurement. No other PK parameter values are published in Snelder 2014. Study 1 dosed 0.04, 0.2, 1 and 5 mg/kg p.o. on separate days in 2 SHR and 1 WKY rat. Prazosin is the only compound described by a power rather than an Emax model; the low exponent (0.0910) means the maximum effect was not reached at the highest dose, so the effect is nearly flat in concentration and this driver column is the one whose absolute scale matters least.",
      source_name        = "C (plasma concentration, Equations 6 and 7)"
    )
  )

  population <- list(
    species        = "rat (male spontaneously hypertensive rat, SHR, and male normotensive Wistar-Kyoto rat, WKY; both Taconic Farms)",
    n_subjects     = 12L,
    n_studies      = 2L,
    study_names    = c(
      "Study 1 (single administrations of different doses on separate days, one vehicle day first; MAP, HR and SV measured with CO and TPR derived; amiloride, amlodipine, enalapril, fasudil, hydrochlorothiazide or prazosin; SHR and WKY rats)",
      "Study 2 (single, sequential or combined administration of atropine 10 mg/kg and/or propranolol 30 mg/kg with a 3 h interval; SHR only, 8 rats)"
    ),
    age_range      = "41-54 weeks (SHR) and 35-38 weeks (WKY) at time of study.",
    weight_range   = "367-504 g (SHR) and 499-600 g (WKY).",
    sex_female_pct = 0,
    sex_notes      = "All animals were male (Methods, 'Animals').",
    disease_state  = "Spontaneous (genetic) hypertension in the SHR arm (BSL_MAP = 155 mmHg) versus normotension in the WKY arm (BSL_MAP = 102 mmHg); no induced disease model. Snelder 2014 states in the Conclusions that applications of the identified system-parameter set are limited to SHR and WKY rats.",
    dose_range     = "Study 1 (p.o., one dose per day on separate days after a vehicle day): amiloride 10 mg/kg; amlodipine 0.3, 1, 3, 10 mg/kg; enalapril 3, 10, 30 mg/kg; fasudil 3, 10, 30 mg/kg; hydrochlorothiazide 0.1, 0.3, 1, 3 mg/kg (first occasion) and 10, 30 mg/kg (second occasion); prazosin 0.04, 0.2, 1, 5 mg/kg. Study 2 (p.o.): atropine 10 mg/kg and propranolol 30 mg/kg, alone, sequentially 3 h apart, or combined. All compounds were given by oral gavage at 2 mL/kg.",
    regions        = "Preclinical; in-life work at Novartis Institutes for BioMedical Research, East Hanover, NJ, USA, with modelling at Leiden Academic Centre for Drug Research, The Netherlands.",
    instrumentation = "Rats were surgically instrumented with BOTH an ascending-aortic transit-time flow probe and a femoral-artery catheter / radiotransmitter (as in Snelder et al. 2013a), giving continuous MAP, HR and CO. Flow cables were disconnected between 1700 h and 0700 h, so overnight only MAP and HR were captured. On experiment days baseline data were collected 0700-1000 h, drug was given at 1000 h (and at 1300 h in Study 2), and collection continued to 1700 h. Rats were housed on a 12 h light/dark cycle with lights on 0600-1800 h.",
    n_ode_states   = 5L,
    notes          = "Ten SHR and two WKY rats were used across both studies, with repeated experiments in the same animals over periods of up to 6 months and sufficient washout between them, so the per-compound counts in Table 1 reflect reuse rather than distinct animals. Data from one SHR in Study 2 were excluded because it learned to disconnect its flow cable and responded far more strongly than the others. SV and TPR were never estimated directly: experimentally they were derived from the measured MAP, CO and HR, and in the modelling BSL_HR, BSL_MAP and BSL_CO were the estimated parameters with BSL_SV and BSL_TPR derived from them. Only HR, MAP and CO carry residual-error models for the same reason. The system was initialised at t = 0 h and pharmacological intervention started at t = 336 h (two weeks, determined empirically) so that the circadian oscillation had settled into its oscillating steady state before dosing; since dosing occurred at 1000 h clock time, model time t = 0 corresponds to 1000 h. Propranolol was administered and modelled but its effect was too small to quantify, so it contributes no drug-specific parameters and no covariate column."
  )

  ini({
    # ==================================================================
    # SYSTEM-SPECIFIC PARAMETERS -- Snelder 2014 Table 5.
    #
    # These are the paper's primary result: the drug-independent
    # description of the rat cardiovascular system, shared by SHR and
    # WKY rats except for the three baselines. Every Table 5 estimate
    # below was cross-checked against its own reported RSE and 95%
    # confidence interval (estimate +/- 1.96 * RSE% * estimate / 100
    # reproduces the printed LLCI / ULCI for all 11 system parameters
    # and all 8 drug parameters); the vignette source-trace table
    # records that check.
    # ==================================================================

    # ---- Strain-specific hemodynamic baselines -----------------------
    # Snelder 2014 estimated BSL_HR, BSL_MAP and BSL_CO separately in
    # each strain within one joint fit, so each carries an explicit
    # stratum suffix (parameter-names.md, "Stratum-suffixed
    # parameters"). BSL_SV and BSL_TPR are NOT parameters: the paper
    # derives them from these three ("in the modelling, the BSL_MAP and
    # BSL_CO and BSL_HR were estimated and BSL_SV and BSL_TPR were
    # derived from these parameters", Data analysis).
    lrbase_HR_shr  <- log(310)  ; label("Baseline heart rate in SHR, BSL_HR_SHR (beats/min)")                        # Snelder 2014 Table 5 (310, RSE 1.12%, 95% CI 303-317)
    lrbase_HR_wky  <- log(323)  ; label("Baseline heart rate in WKY rats, BSL_HR_WKY (beats/min)")                   # Snelder 2014 Table 5 (323, RSE 1.61%, 95% CI 313-333)
    lrbase_MAP_shr <- log(155)  ; label("Baseline mean arterial pressure in SHR, BSL_MAP_SHR (mmHg)")                # Snelder 2014 Table 5 (155, RSE 0.684%, 95% CI 153-157)
    lrbase_MAP_wky <- log(102)  ; label("Baseline mean arterial pressure in WKY rats, BSL_MAP_WKY (mmHg)")           # Snelder 2014 Table 5 (102, RSE 0.884%, 95% CI 100-104)
    lrbase_CO_shr  <- log(69.0) ; label("Baseline cardiac output in SHR, BSL_CO_SHR (mL/min)")                       # Snelder 2014 Table 5 (69.0, RSE 4.17%, 95% CI 63.4-74.6)
    lrbase_CO_wky  <- log(129)  ; label("Baseline cardiac output in WKY rats, BSL_CO_WKY (mL/min)")                  # Snelder 2014 Table 5 (129, RSE 1.47%, 95% CI 125-133)

    # ---- First-order dissipation rate constants ----------------------
    # One per turnover state. kout_HR is ~90-fold faster than kout_SV,
    # which is what separates the fast HR-mediated and slow
    # SV-mediated routes to a MAP change (Figure 4).
    lkout_HR  <- log(11.6)  ; label("First-order dissipation rate constant of HR, kout_HR (1/h)")                    # Snelder 2014 Table 5 (11.6, RSE 19.1%, 95% CI 7.27-15.9)
    lkout_SV  <- log(0.126) ; label("First-order dissipation rate constant of SV, kout_SV (1/h)")                    # Snelder 2014 Table 5 (0.126, RSE 30.7%, 95% CI 0.0501-0.202)
    lkout_TPR <- log(3.58)  ; label("First-order dissipation rate constant of TPR, kout_TPR (1/h)")                  # Snelder 2014 Table 5 (3.58, RSE 29.1%, 95% CI 1.54-5.62)

    # ---- MAP negative feedback ---------------------------------------
    # ONE feedback constant is shared by all three production rates
    # (Snelder 2014 Equation 2: "FB is a constant representing the
    # magnitude of the negative feedback of MAP on HR, SV and TPR"),
    # unlike the Snelder 2013 predecessor which had separate FB1 on CO
    # and FB2 on TPR. FB itself is not constant across animals: it is a
    # power function of the individual baseline MAP normalised to the
    # SHR typical value of 155 mmHg (Equation 9), so a low-BSL_MAP
    # animal has stronger feedback. The exponent is negative.
    lfb0        <- log(0.00290) ; label("Negative feedback of MAP on HR, SV and TPR production for a typical SHR, FB0 (1/mmHg)")   # Snelder 2014 Table 5 (0.00290, RSE 5.93%, 95% CI 0.00256-0.00324)
    e_bslmap_fb <- -1.98        ; label("Power exponent of the individual baseline MAP (normalised to 155 mmHg) on FB, FB0_MAP (unitless)")  # Snelder 2014 Table 5 (-1.98, RSE 10.6%, 95% CI -2.39 to -1.57) and Equation 9

    # ---- Direct HR-on-SV coupling ------------------------------------
    # Log-linear: SV = SVT * (1 - HR_SV * ln(HR / BSL_HR)) (Equation 2).
    # Represents the shortening of left-ventricular filling time when
    # the cardiac interval shortens. At HR = BSL_HR the term is 1.
    lhrsv <- log(0.312) ; label("Magnitude of the direct inverse log-linear HR-on-SV coupling, HR_SV (unitless)")    # Snelder 2014 Table 5 (0.312, RSE 15.6%, 95% CI 0.216-0.408)

    # ---- Handling effect (Equation 4) --------------------------------
    # Brief manual restraint and oral gavage -- experienced directly or
    # sensed from a bystander rat in the same room -- transiently raise
    # HR, TPR, CO and MAP and lower SV, independently of drug. The
    # empirical HD function of Visser 2006 is P_X * exp(-kHD*(t - tHD))
    # for t > tHD, applied multiplicatively alongside the drug effect
    # on the HR and TPR production rates only. kHD = 4.70 1/h gives a
    # 8.8 min half-life, so the artefact is essentially gone within an
    # hour -- the paper shows it at 3 h post-vehicle in Figure 2.
    lkhd     <- log(4.70)  ; label("First-order rate of disappearance of the handling effect, kHD (1/h)")            # Snelder 2014 Table 5 (4.70, RSE 8.19%, 95% CI 3.95-5.45)
    lphd_HR  <- log(0.632) ; label("Magnitude of the handling effect on HR production, P_HR (unitless)")             # Snelder 2014 Table 5 (0.632, RSE 9.67%, 95% CI 0.512-0.752)
    lphd_TPR <- log(0.331) ; label("Magnitude of the handling effect on TPR production, P_TPR (unitless)")           # Snelder 2014 Table 5 (0.331, RSE 12.9%, 95% CI 0.247-0.415)

    # ---- Circadian rhythm (Equation 3) -------------------------------
    # Two 24 h cosines, one multiplying Kin_HR and one multiplying
    # Kin_TPR; the rhythms in SV, CO and MAP follow from the feedback
    # rather than being modelled directly. The amplitudes could not be
    # distinguished, so ampTPR is FIXED to ampHR (Table 5, "ampTPR:
    # Fixed to ampHR", encoded here as a fixed ratio of 1); the
    # horizontal displacements ARE significantly different and both are
    # estimated. Amplitudes and displacements are kept on the linear
    # scale because a phase offset can in principle be signed.
    hor_HR        <- 8.73   ; label("Horizontal displacement of the HR circadian cosine, hor_HR (h)")                # Snelder 2014 Table 5 (8.73, RSE 3.10%, 95% CI 8.20-9.26)
    amp_HR        <- 0.0918 ; label("Amplitude of the HR circadian cosine as a fraction of Kin_HR, amp_HR")          # Snelder 2014 Table 5 (0.0918, RSE 5.15%, 95% CI 0.0825-0.101); Results: "estimated to be 0.09 indicating that the variation in Kin_HR and Kin_TPR is maximally 9% during the day"
    hor_TPR       <- 19.3   ; label("Horizontal displacement of the TPR circadian cosine, hor_TPR (h)")              # Snelder 2014 Table 5 (19.3, RSE 1.92%, 95% CI 18.6-20.0)
    amp_TPR_ratio <- fixed(1) ; label("Ratio of the TPR circadian amplitude to the HR circadian amplitude (unitless)")  # Snelder 2014 Table 5 row "ampTPR: Fixed to ampHR"; corroborated by Fu 2023 Suppl. S1 $THETA TH15 = 1 FIX, which re-encodes this same model

    # ==================================================================
    # DRUG-SPECIFIC PARAMETERS -- Snelder 2014 Table 5.
    #
    # SITE and DIRECTION of each effect come from Table 4 ("Effect"
    # column) and are restated in the figure captions (Figure 3,
    # Figure S1, Figure S2). Table 5 reports only the MAGNITUDE of each
    # effect, so the direction is applied as an explicit sign in
    # model() rather than being folded into these values -- that keeps
    # every number here byte-identical to the printed table.
    #
    # Emax is FIXED to 1 for the five Emax compounds: Table 3
    # Assumption 2 states that "for compounds for which the maximum
    # effect was not observed, complete inhibition (i.e. Emax = 1) was
    # assumed at infinite concentrations to ensure identification of
    # the EC50 parameter", and each Table 5 sub-heading repeats "Emax
    # model with Emax fixed to 1".
    #
    # Propranolol was dosed (Study 2, 30 mg/kg p.o.) but "the effect of
    # propranolol was too small to be quantified" (Results, 'Drug
    # effects'), so it has no row in Table 5 and no parameters here.
    # ==================================================================

    # ---- Amiloride: Emax inhibition of SV production ------------------
    emax_amiloride  <- fixed(1) ; label("Maximum fractional inhibition of SV production by amiloride, Emax (unitless)")  # Snelder 2014 Table 5 sub-heading "Amiloride: Emax model with Emax fixed to 1"; Table 3 Assumption 2
    lec50_amiloride <- log(245) ; label("Amiloride concentration giving half-maximal inhibition of SV production, EC50 (ng/mL)")  # Snelder 2014 Table 5 (245, RSE 25.1%, 95% CI 125-365)

    # ---- Amlodipine: Emax inhibition of TPR production ----------------
    emax_amlodipine  <- fixed(1)  ; label("Maximum fractional inhibition of TPR production by amlodipine, Emax (unitless)")  # Snelder 2014 Table 5 sub-heading "Amlodipine: Emax model with Emax fixed to 1"
    lec50_amlodipine <- log(82.8) ; label("Amlodipine concentration giving half-maximal inhibition of TPR production, EC50 (ng/mL)")  # Snelder 2014 Table 5 (82.8, RSE 4.99%, 95% CI 74.7-90.9)

    # ---- Atropine: linear stimulation of HR production ----------------
    # The only stimulating effect in the paper. Results, 'Drug effects':
    # "As atropine had a stimulating effect on Kin_HR, applying a linear
    # concentration-effect relationship did not result in problems with
    # parameter optimization."
    slope_atropine <- 0.00149 ; label("Linear slope of the atropine stimulation of HR production, SL (per ng/mL)")   # Snelder 2014 Table 5 (0.00149, RSE 32.3%, 95% CI 0.000547-0.00243)

    # ---- Enalapril: Emax inhibition of BOTH TPR and SV production -----
    # One shared EC50: "Initially, different EC50 values were estimated.
    # However, confidence intervals overlapped indicating that the EC50
    # values for the two effects could not be distinguished." The effect
    # is driven by an effect-compartment concentration rather than
    # plasma (Equation 8), the only compound for which that was needed.
    emax_enalapril  <- fixed(1)   ; label("Maximum fractional inhibition of TPR and SV production by enalapril, Emax (unitless)")  # Snelder 2014 Table 5 sub-heading "Enalapril: Emax model with Emax fixed to 1"
    lec50_enalapril <- log(1200)  ; label("Enalapril effect-site concentration giving half-maximal inhibition of TPR and SV production, EC50 (ng/mL)")  # Snelder 2014 Table 5 (1200, RSE 4.03%, 95% CI 1110-1290)
    lke0_enalapril  <- log(0.163) ; label("Enalapril effect-compartment equilibration rate constant, ke0 (1/h)")     # Snelder 2014 Table 5 (0.163, RSE 5.07%, 95% CI 0.147-0.179); Results: "The half-life of this additional delay was 4.3 h" (ln(2)/0.163 = 4.25 h)

    # ---- Fasudil: Emax inhibition of TPR production -------------------
    emax_fasudil  <- fixed(1)    ; label("Maximum fractional inhibition of TPR production by fasudil, Emax (unitless)")  # Snelder 2014 Table 5 sub-heading "Fasudil: Emax model with Emax fixed to 1"
    lec50_fasudil <- log(0.172)  ; label("Fasudil concentration giving half-maximal inhibition of TPR production, EC50 (ng/mL)")  # Snelder 2014 Table 5 (0.172, RSE 18.4%, 95% CI 0.110-0.234). NOTE: ~1900-fold more potent than the 321 ng/mL of the Snelder 2013 predecessor; both are internally consistent with their own RSE and CI, so this is a between-publication discrepancy, not a transcription error. See vignette Errata.

    # ---- Hydrochlorothiazide: Emax inhibition of SV production --------
    emax_hctz  <- fixed(1)     ; label("Maximum fractional inhibition of SV production by hydrochlorothiazide, Emax (unitless)")  # Snelder 2014 Table 5 sub-heading "HCTZ: Emax model with Emax fixed to 1"
    lec50_hctz <- log(28900)   ; label("Hydrochlorothiazide concentration giving half-maximal inhibition of SV production, EC50 (ng/mL)")  # Snelder 2014 Table 5 (28 900, RSE 7.65%, 95% CI 24 600-33 200)

    # ---- Prazosin: power-model inhibition of TPR production -----------
    # EFF = SL * C^POW (Equation 7). Results: "The exponent of this
    # relationship was low (0.0910) indicating that the maximum effect
    # was not reached for the highest dose evaluated." Unlike an Emax
    # form this is unbounded above, so it must not be extrapolated far
    # beyond the studied 0.04-5 mg/kg range.
    slope_prazosin <- 0.328  ; label("Magnitude coefficient of the prazosin power model on TPR production, SL (per ng/mL raised to POW)")  # Snelder 2014 Table 5 (0.328, RSE 5.58%, 95% CI 0.292-0.364); Table 5 prints the unit as (ng/mL)^-1, which is exact only when POW = 1
    pow_prazosin   <- 0.0910 ; label("Exponent of the prazosin power model on TPR production, POW (unitless)")       # Snelder 2014 Table 5 (0.0910, RSE 6.05%, 95% CI 0.0802-0.102)

    # ==================================================================
    # INTER-INDIVIDUAL VARIABILITY -- Snelder 2014 Table 5.
    #
    # Log-normal IIV on the three estimated baselines only (Results:
    # "interindividual variations in the baseline values of the
    # parameters, BSL_HR, BSL_MAP and BSL_CO, were allowed"). Table 5
    # reports these as CV%, and the CV% is the approximate form
    # sqrt(omega) * 100 rather than sqrt(exp(omega) - 1) * 100: Fu 2023
    # Supplemental S1, which re-encodes this same model as a NONMEM
    # control stream, carries $OMEGA 0.00372 / 0.00137 / 0.0515 whose
    # square roots are exactly 6.1% / 3.7% / 22.7%. The variances below
    # are therefore (CV / 100)^2 from the Table 5 CV% values.
    # ==================================================================
    etalrbase_HR  ~ 0.003721 # Snelder 2014 Table 5 BSL_HR (CV 6.1%, 95% CI 4.36-7.47); variance = (6.1/100)^2
    etalrbase_MAP ~ 0.001369 # Snelder 2014 Table 5 BSL_MAP (CV 3.7%, 95% CI 2.67-4.49); variance = (3.7/100)^2
    etalrbase_CO  ~ 0.051529 # Snelder 2014 Table 5 BSL_CO (CV 22.7%, 95% CI 18.09-26.57); variance = (22.7/100)^2

    # ==================================================================
    # RESIDUAL VARIABILITY -- Snelder 2014 Table 5.
    #
    # Proportional on each of the three MEASURED readouts (Results:
    # "The residual errors of HR, MAP and CO were best described by
    # proportional residual error models. The residual errors of TPR
    # and SV were derived from these parameters."), so SV and TPR are
    # exposed as outputs but carry no error model of their own.
    # nlmixr2 propSd is an SD, and the Table 5 CV% for a proportional
    # model is that SD as a percentage, so propSd = CV / 100 directly.
    # ==================================================================
    propSd_HR  <- 0.078 ; label("Proportional residual SD on heart rate (fraction)")            # Snelder 2014 Table 5 Prop. Res. Error HR (CV 7.8%, 95% CI 7.26-8.22)
    propSd_MAP <- 0.060 ; label("Proportional residual SD on mean arterial pressure (fraction)") # Snelder 2014 Table 5 Prop. Res. Error MAP (CV 6.0%, 95% CI 5.44-6.57)
    propSd_CO  <- 0.069 ; label("Proportional residual SD on cardiac output (fraction)")         # Snelder 2014 Table 5 Prop. Res. Error CO (CV 6.9%, 95% CI 5.72-7.83)
  })

  model({
    # ==================================================================
    # 1. Strain-specific baselines.
    #
    # STRAIN_SHR selects between the two Table 5 baseline triples. For
    # a binary indicator the interpolation below is exact: it returns
    # the WKY value at STRAIN_SHR = 0 and the SHR value at 1. The IIV
    # etas act on the selected log baseline, so an animal's baseline
    # distribution is centred on its own strain's typical value.
    #
    # BSL_SV and BSL_TPR are derived, not estimated (Data analysis):
    #   CO = HR * SV      =>  BSL_SV  = BSL_CO  / BSL_HR
    #   MAP = CO * TPR    =>  BSL_TPR = BSL_MAP / BSL_CO
    # With the SHR values: 69.0 / 310 = 0.2226 mL/beat and
    # 155 / 69.0 = 2.246 mmHg/(mL/min).
    # ==================================================================
    lrbase_HR  <- lrbase_HR_wky  + STRAIN_SHR * (lrbase_HR_shr  - lrbase_HR_wky)
    lrbase_MAP <- lrbase_MAP_wky + STRAIN_SHR * (lrbase_MAP_shr - lrbase_MAP_wky)
    lrbase_CO  <- lrbase_CO_wky  + STRAIN_SHR * (lrbase_CO_shr  - lrbase_CO_wky)

    bsl_HR  <- exp(lrbase_HR  + etalrbase_HR)
    bsl_MAP <- exp(lrbase_MAP + etalrbase_MAP)
    bsl_CO  <- exp(lrbase_CO  + etalrbase_CO)
    bsl_SV  <- bsl_CO  / bsl_HR
    bsl_TPR <- bsl_MAP / bsl_CO

    # ==================================================================
    # 2. Linear-scale system parameters.
    #
    # The feedback constant is a power function of the INDIVIDUAL
    # baseline MAP normalised to the SHR typical value of 155 mmHg
    # (Snelder 2014 Equation 9):
    #   FB = FB0 * (IBSL_MAP / TVBSL_MAP_SHR)^FB0_MAP
    # The reference 155 mmHg is the SHR typical BSL_MAP of Table 5 and
    # is used for BOTH strains -- with the WKY typical BSL_MAP of 102
    # mmHg this gives (102/155)^-1.98 = 2.29, reproducing the paper's
    # statement that "the feedback is about twofold higher in WKY rats
    # as compared with SHR".
    # ==================================================================
    kout_HR  <- exp(lkout_HR)
    kout_SV  <- exp(lkout_SV)
    kout_TPR <- exp(lkout_TPR)
    hrsv     <- exp(lhrsv)
    khd      <- exp(lkhd)
    phd_HR   <- exp(lphd_HR)
    phd_TPR  <- exp(lphd_TPR)
    ke0_enalapril <- exp(lke0_enalapril)

    fb <- exp(lfb0) * (bsl_MAP / 155)^e_bslmap_fb

    # ==================================================================
    # 3. Steady-state production rate constants -- Snelder 2014
    #    Equation 5. Kin is expressed in terms of BSL and kout WITHOUT
    #    accounting for the circadian rhythm, because the oscillating
    #    steady state has no analytical solution; the paper compensates
    #    by initialising the system two weeks (336 h) before the first
    #    pharmacological intervention so the oscillation has settled.
    #
    #    Units: Kin_HR in (beats/min)/h, Kin_SV in (mL/beat)/h,
    #    Kin_TPR in (mmHg/(mL/min))/h.
    # ==================================================================
    fb_bsl  <- 1 - fb * bsl_MAP
    kin_HR  <- kout_HR  * bsl_HR  / fb_bsl
    kin_SV  <- kout_SV  * bsl_SV  / fb_bsl
    kin_TPR <- kout_TPR * bsl_TPR / fb_bsl

    # ==================================================================
    # 4. Circadian rhythm -- Snelder 2014 Equation 3. Two 24 h cosines,
    #    multiplicative on the HR and TPR production rates. ampTPR is
    #    fixed to ampHR (Table 5), encoded as the fixed ratio
    #    amp_TPR_ratio = 1.
    #
    #    Model time t = 0 corresponds to 1000 h clock time: dosing
    #    happened at 1000 h (Experimental design) and pharmacological
    #    intervention started at t = 336 h = exactly 14 x 24 h. On that
    #    clock the HR cosine peaks at t = 15.3 h (0116 h) and the TPR
    #    cosine at t = 4.7 h (1442 h), i.e. HR peaks in the dark phase
    #    and TPR in the light phase, as expected for a nocturnal
    #    species on the paper's 0600-1800 h light cycle.
    # ==================================================================
    cr_HR  <- amp_HR                 * cos(2 * pi * (t + hor_HR)  / 24)
    cr_TPR <- amp_HR * amp_TPR_ratio * cos(2 * pi * (t + hor_TPR) / 24)

    # ==================================================================
    # 5. Handling effect -- Snelder 2014 Equation 4:
    #      HD_X = P_X * exp(-kHD * (t - tHD))  for t > tHD
    #    realised here as a first-order decay state that receives a unit
    #    impulse at each handling event, which is the exact same
    #    function for a single event and superposes correctly when a rat
    #    is handled more than once (Study 2 dosed at 1000 h and 1300 h).
    #    Handling raises the HR and TPR production rates only; the fall
    #    in SV that Figure 2 shows is a consequence of the direct
    #    HR-on-SV coupling, not a separate handling term.
    # ==================================================================
    hd_HR  <- phd_HR  * handling
    hd_TPR <- phd_TPR * handling

    # ==================================================================
    # 6. Drug effects -- Snelder 2014 Equations 6 and 7, one term per
    #    compound, routed to that compound's site of action per Table 4
    #    and the Figure 3 / S1 / S2 captions.
    #
    #    Equation 6 applies the drug through the factor (1 + EFF), so
    #    EFF is NEGATIVE for the six inhibitory compounds and POSITIVE
    #    for atropine, the only stimulator. Table 5 reports magnitudes
    #    only, so the sign is written explicitly here.
    #
    #    Terms are summed so a dataset containing several compounds --
    #    as the paper's own simultaneous fit did, and as Study 2's
    #    atropine + propranolol combination arm required -- works
    #    unchanged; in practice only one CP_ column is non-zero on any
    #    given record and every column defaults to 0.
    #
    #    Enalapril is driven by the effect-site concentration rather
    #    than plasma, and appears in BOTH eff_TPR and eff_SV with the
    #    same EC50.
    # ==================================================================
    eff_HR <-
      slope_atropine * CP_ATROPINE_NGML

    eff_SV <-
      -emax_amiloride * CP_AMILORIDE_NGML / (exp(lec50_amiloride) + CP_AMILORIDE_NGML) +
      -emax_hctz      * CP_HCTZ_NGML      / (exp(lec50_hctz)      + CP_HCTZ_NGML)      +
      -emax_enalapril * effect            / (exp(lec50_enalapril) + effect)

    eff_TPR <-
      -emax_amlodipine * CP_AMLODIPINE_NGML / (exp(lec50_amlodipine) + CP_AMLODIPINE_NGML) +
      -emax_fasudil    * CP_FASUDIL_NGML    / (exp(lec50_fasudil)    + CP_FASUDIL_NGML)    +
      -emax_enalapril  * effect             / (exp(lec50_enalapril)  + effect)             +
      -slope_prazosin  * CP_PRAZOSIN_NGML^pow_prazosin

    # ==================================================================
    # 7. Derived hemodynamic quantities -- Snelder 2014 Equation 2.
    #    Actual SV adds the direct inverse log-linear HR coupling on top
    #    of the feedback-driven turnover state SVT; CO and MAP are pure
    #    products, with no circadian or error term of their own.
    # ==================================================================
    sv  <- svt * (1 - hrsv * log(hr / bsl_HR))
    co  <- hr * sv
    map <- co * tpr

    # ==================================================================
    # 8. ODE system -- Snelder 2014 Equation 6 (which is Equation 4 with
    #    the drug effect added). Note the asymmetry, which is the
    #    paper's: the circadian rhythm and the handling effect act on HR
    #    and TPR only, while all three states receive the MAP feedback
    #    and can receive a drug effect.
    # ==================================================================
    d/dt(hr)  <- kin_HR  * (1 + cr_HR)  * (1 - fb * map) * (1 + eff_HR  + hd_HR)  - kout_HR  * hr
    d/dt(svt) <- kin_SV                 * (1 - fb * map) * (1 + eff_SV)           - kout_SV  * svt
    d/dt(tpr) <- kin_TPR * (1 + cr_TPR) * (1 - fb * map) * (1 + eff_TPR + hd_TPR) - kout_TPR * tpr

    # Enalapril effect compartment -- Snelder 2014 Equation 8,
    # dCe/dt = ke0 * (C - Ce), driven by the plasma covariate column.
    d/dt(effect) <- ke0_enalapril * (CP_ENALAPRIL_NGML - effect)

    # Handling artefact -- unit impulses dosed into this state; see
    # section 5. First-order decay at kHD.
    d/dt(handling) <- -khd * handling

    # Initial conditions: the system starts at its (non-oscillating)
    # baseline steady state, which section 3 constructs Kin to satisfy.
    hr(0)       <- bsl_HR
    svt(0)      <- bsl_SV
    tpr(0)      <- bsl_TPR
    effect(0)   <- 0
    handling(0) <- 0

    # ==================================================================
    # 9. Observations. HR, MAP and CO were measured directly and each
    #    carries its own proportional residual error. SV and TPR were
    #    derived experimentally from those three and are exposed as
    #    outputs without an error model of their own, matching the
    #    paper ("The residual errors of TPR and SV were derived from
    #    these parameters").
    # ==================================================================
    HR  <- hr
    MAP <- map
    CO  <- co
    SV  <- sv
    TPR <- tpr

    HR  ~ prop(propSd_HR)
    MAP ~ prop(propSd_MAP)
    CO  ~ prop(propSd_CO)
  })
}
