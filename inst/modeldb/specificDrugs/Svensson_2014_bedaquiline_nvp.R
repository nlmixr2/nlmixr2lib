Svensson_2014_bedaquiline_nvp <- function() {
  description <- "Three-compartment population PK model for bedaquiline (BDQ) and a two-compartment N-desmethyl metabolite M2 in HIV-1-infected ART-naive adult volunteers following single 400 mg oral doses, with Savic 2007 analytical transit-compartment absorption (non-integer NN feeding a first-order depot at rate ka), fixed allometric scaling on disposition (0.75 on CL/Q at 70 kg, 1 on Vc/Vp), and multiplicative nevirapine (NVP) drug-drug-interaction factors of 0.915 on bedaquiline and 1.05 on M2 apparent clearances during steady-state NVP co-administration (study C117). The factors are fixed-effects only because BSV on the NVP interaction effects was not estimated."
  reference <- paste(
    "Svensson E. M., Dooley K. E., Karlsson M. O. (2014).",
    "Impact of lopinavir-ritonavir or nevirapine on bedaquiline exposures",
    "and potential implications for patients with tuberculosis-HIV",
    "coinfection.",
    "Antimicrobial Agents and Chemotherapy 58(11):6406-6412.",
    "doi:10.1128/AAC.03246-14.",
    "Structural model inherited from Svensson 2013 (efavirenz DDI;",
    "doi:10.1128/AAC.00191-13); see modellib('Svensson_2013_bedaquiline').",
    sep = " "
  )
  vignette <- "Svensson_2014_bedaquiline_arv"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight (used for allometric scaling around 70 kg)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline body weight. Allometric scaling applied with fixed exponents 0.75 on apparent clearances (CL/F, Q1/F, Q2/F, CL_M2, Q_M2) and 1 on apparent volumes (V/F, VP1/F, VP2/F, V_M2, VP_M2) around a 70 kg reference adult. Svensson 2014 Materials and Methods 'Nonlinear mixed-effects modeling' states 'Allometric scaling of disposition parameters with body weight as the size descriptor and fixed coefficients (0.75 for clearance [CL] and 1 for volume of distribution) was applied.' Reference weight of 70 kg is confirmed by Supplementary Table S1b footnote b: 'Disposition parameters for a typical individual of 70 kg, allometric scaling with body weight and fixed coefficients 0.75 for CL and 1 for V applied.'",
      source_name        = "WT"
    ),
    CONMED_NVP = list(
      description        = "Concomitant nevirapine (NVP) co-administration at steady state with full CYP3A4 induction (1 = on twice-daily 200 mg NVP at steady state, 0 = not on NVP or pre-induction lag).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not on NVP or pre-induction window)",
      notes              = "Subject- and time-varying indicator of concomitant NVP co-administration at full induction. Study C117 (n = 16 HIV-1-positive ART-naive adult volunteers) was a single-sequence study with two 400 mg bedaquiline doses; NVP at standard doses (200 mg once daily for 2 weeks, followed by 200 mg twice daily) was started before the second bedaquiline dose, and participants received at least 4 weeks of the twice-daily NVP dosing before the second bedaquiline dose was administered. Svensson 2014 Materials and Methods state 'the impacts of NVP (induction) were assumed to start after 2 weeks of twice-daily administration.' Multiplicative factor on apparent CL during co-administration: cl_bdq_eff = cl_bdq_base * e_nvp_cl ^ CONMED_NVP with e_nvp_cl = 0.915 (RSE 5.9%) for bedaquiline (BDQ CL falls to 91.5% of the no-NVP value) and cl_m2_eff = cl_m2_base * e_nvp_cl_m2 ^ CONMED_NVP with e_nvp_cl_m2 = 1.05 (RSE 10.3%) for M2 (M2 CL rises to 105% of the no-NVP value). Svensson 2014 Supplementary Table S1b 'EFF NVP on BDQ CL = 0.915' and 'EFF NVP on M2 CL = 1.05'. For simulation, set CONMED_NVP = 1 on observation rows >= 2 weeks after the start of NVP 200 mg twice-daily dosing and 0 otherwise.",
      source_name        = "NVP"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 16L,
    n_studies      = 1L,
    age_range      = "22-51 years",
    age_median     = "30 years",
    weight_range   = "48-71 kg",
    weight_median  = "55 kg",
    sex_female_pct = 37.5,
    race_ethnicity = c(Black = 50.0, White = 0.0, `Mixed race` = 50.0),
    disease_state  = "HIV-1-infected ART-naive adult volunteers (18-65 years; BMI 18-30 kg/m^2) with a medical indication to start ART with NVP plus two NRTIs. Documented HIV-1 infection but no history of ART. Women of childbearing potential, individuals with active AIDS-defining illnesses or TB, and individuals with a history of substance abuse were excluded. Subjects previously enrolled in trials involving bedaquiline were ineligible.",
    dose_range     = "Two single 400 mg oral doses of bedaquiline given several weeks apart in a single-sequence design; NVP at standard doses (200 mg once daily for 2 weeks, then 200 mg twice daily) was started before the second bedaquiline dose with at least 4 weeks of twice-daily NVP prior to the second bedaquiline dose. PK samples were collected pre-dose and at 1, 2, 3, 4, 5, 6, 8, 12, 24, 48, 72, 120, 168, 216, 264, and 336 h after each bedaquiline dose (17 samples per dose per analyte).",
    regions        = "Not stated in the publication beyond study sponsor (Tibotec / Janssen) and ClinicalTrials.gov registration NCT00910806.",
    study_id       = "C117 (NVP DDI; NCT00910806).",
    notes          = "Baseline demographics from Svensson 2014 Table 1. 528 BDQ + 528 M2 PK samples available; 1 BDQ and 33 M2 samples were below limit of quantification (LLOQ 1.00 ng/mL) and were omitted from modelling. Bedaquiline and M2 concentrations were determined by LC-MS/MS validated to FDA guidelines."
  )

  ini({
    # Final parameter estimates from Svensson 2014 Supplementary Table S1b
    # 'Parameter estimates including precision (RSE) for study C117 (bedaquiline
    # drug-drug interaction with nevirapine)'. The structural model
    # (3-compartment BDQ + 2-compartment M2 with NN transit absorption +
    # first-order ka into the depot) was inherited from Svensson 2013
    # (doi:10.1128/AAC.00191-13); see modellib('Svensson_2013_bedaquiline'). All
    # final parameter values were re-estimated on the C117 data in this paper.

    # Structural absorption: Savic 2007 transit-compartment model (non-integer
    # number of transit compartments NN feeding the depot at rate
    # KTR = (NN+1)/MTT, which then absorbs into central at rate ka).
    lmtt    <- log(3.37)    ; label("Mean transit time MTT (h, on log scale)")                                            # Svensson 2014 Supplementary Table S1b 'MTT = 3.37 h' (RSE 7.7%)
    lka     <- log(0.131)   ; label("First-order absorption rate constant ka (1/h, on log scale)")                         # Svensson 2014 Supplementary Table S1b 'KA = 0.131 1/h' (RSE 9.8%)
    lnn     <- log(4.48)    ; label("Number of transit compartments NN (unitless, on log scale)")                          # Svensson 2014 Supplementary Table S1b 'NN = 4.48' (RSE 9.3%)

    # Structural disposition: 3-compartment bedaquiline (parent) +
    # 2-compartment M2 metabolite. Apparent volumes and clearances are CL/F,
    # V/F for the parent and CL/(F*fm), V/(F*fm) for the metabolite (fm not
    # separately identifiable from F-relative parameters).
    lcl     <- log(3.34)    ; label("Apparent bedaquiline clearance CL/F at 70 kg (L/h)")                                   # Svensson 2014 Supplementary Table S1b 'CL/F = 3.34 L/h' (RSE 9.5%)
    lvc     <- log(11.2)    ; label("Apparent bedaquiline central volume V/F at 70 kg (L)")                                 # Svensson 2014 Supplementary Table S1b 'V/F = 11.2 L' (RSE 42.0%)
    lq      <- log(7.03)    ; label("Apparent bedaquiline first-peripheral inter-compartmental clearance Q1/F at 70 kg (L/h)") # Svensson 2014 Supplementary Table S1b 'Q1/F = 7.03 L/h' (RSE 8.3%)
    lvp     <- log(4000)    ; label("Apparent bedaquiline first-peripheral volume VP1/F at 70 kg (L)")                       # Svensson 2014 Supplementary Table S1b 'VP1/F = 4000 L' (RSE 11.4%)
    lq2     <- log(5.26)    ; label("Apparent bedaquiline second-peripheral inter-compartmental clearance Q2/F at 70 kg (L/h)") # Svensson 2014 Supplementary Table S1b 'Q2/F = 5.26 L/h' (RSE 19.2%)
    lvp2    <- log(164)     ; label("Apparent bedaquiline second-peripheral volume VP2/F at 70 kg (L)")                      # Svensson 2014 Supplementary Table S1b 'VP2/F = 164 L' (RSE 8.2%)
    lcl_m2  <- log(16.0)    ; label("Apparent M2 metabolite clearance CL_M2/(F*fm) at 70 kg (L/h)")                          # Svensson 2014 Supplementary Table S1b 'CLM2/F/fm = 16.0 L/h' (RSE 9.9%)
    lvc_m2  <- log(824)     ; label("Apparent M2 metabolite central volume V_M2/(F*fm) at 70 kg (L)")                        # Svensson 2014 Supplementary Table S1b 'VM2/F/fm = 824 L' (RSE 10.5%)
    lq_m2   <- log(131)     ; label("Apparent M2 metabolite inter-compartmental clearance Q_M2/(F*fm) at 70 kg (L/h)")       # Svensson 2014 Supplementary Table S1b 'Q1M2/F/fm = 131 L/h' (RSE 11.4%)
    lvp_m2  <- log(3090)    ; label("Apparent M2 metabolite peripheral volume VP_M2/(F*fm) at 70 kg (L)")                    # Svensson 2014 Supplementary Table S1b 'VP1M2/F/fm = 3090 L' (RSE 8.3%)

    # Bioavailability is implicitly F = 1 because CL and V are reported as
    # apparent F-relative values. Svensson 2014 Materials and Methods
    # 'All disposition parameters were estimated as relative to the bioavailability.'

    # Allometric scaling (theory-based, fixed): exponent 0.75 on apparent
    # clearances and 1 on apparent volumes around 70 kg (Svensson 2014
    # Supplementary Table S1b footnote b).
    e_wt_cl_q  <- fixed(0.75) ; label("Allometric exponent on apparent CL and Q (CL/F, Q1/F, Q2/F, CL_M2, Q_M2; unitless)")  # Svensson 2014 Methods (fixed)
    e_wt_vc_vp <- fixed(1)    ; label("Allometric exponent on apparent Vc and Vp (V/F, VP1/F, VP2/F, V_M2, VP_M2; unitless)") # Svensson 2014 Methods (fixed)

    # Drug-drug-interaction multiplicative factors on the bedaquiline and M2
    # apparent clearances during nevirapine (NVP) co-administration at full
    # induction. Svensson 2014 Results 'NVP was not found to have a
    # significant effect on BDQ bioavailability, and the effects on CL
    # were small, i.e., changes to 91.5% (relative standard error [RSE],
    # 5.9%) and 105% (RSE, 10.3%) of BDQ and M2 CL values without
    # comedication, respectively.' The paper reports two SEPARATE factors
    # (one for BDQ CL, one for M2 CL) but BSV on the interaction effects
    # was NOT estimated (Supplementary Table S1b leaves BSV columns blank
    # for EFF NVP on BDQ CL and EFF NVP on M2 CL), so the factors apply
    # as fixed-effects-only structural multipliers.
    e_nvp_cl    <- 0.915 ; label("Multiplicative factor on bedaquiline apparent CL during NVP co-administration (unitless)") # Svensson 2014 Supplementary Table S1b 'EFF NVP on BDQ CL = 0.915' (RSE 5.9%)
    e_nvp_cl_m2 <- 1.05  ; label("Multiplicative factor on M2 apparent CL during NVP co-administration (unitless)")          # Svensson 2014 Supplementary Table S1b 'EFF NVP on M2 CL = 1.05' (RSE 10.3%)

    # Inter-individual variability (between-subject, BSV). Final estimates
    # converted from the paper's CV% notation via omega^2 = log(1 + CV^2)
    # for log-normal IIV. Svensson 2014 Supplementary Table S1b reports
    # diagonal BSV terms on CL, CLM2, V, Q1, VM2, VP1M2 only (the C117
    # study did NOT estimate a correlated CL ~ CL_M2 block in the same
    # way the C110 study did; Supplementary Table S1b row 'BSV CL~BSV
    # CLM2 = 24.2%' reports a correlation but the off-diagonal entry of
    # the BLOCK(2) is the cross-correlation of the diagonal BSVs as
    # reported, not a separately-estimated covariance term). The
    # encoding here keeps the diagonal BSVs on CL, V, Q1, CL_M2, V_M2,
    # VP1_M2 and DROPS the 24.2% cross-correlation between BSV CL and
    # BSV CL_M2, BOV F = 22.6% CV, BSV F = 20.5% CV, and BOV MTT = 32.9%
    # CV (between-occasion variabilities and between-subject
    # variability on bioavailability have no idiomatic encoding in
    # nlmixr2lib). See vignette Assumptions and deviations.
    etalcl     ~ 0.040039  # var = log(1 + 0.204^2) = 0.040039 ; Svensson 2014 Supplementary Table S1b BSV CL = 20.4% (RSE 18.0%)
    etalcl_m2  ~ 0.048966  # var = log(1 + 0.226^2) = 0.048966 ; Svensson 2014 Supplementary Table S1b BSV CLM2 = 22.6% (RSE 19.2%)
    etalvc     ~ 0.625537  # var = log(1 + 0.980^2) = 0.625537 ; Svensson 2014 Supplementary Table S1b BSV V = 98.0% (RSE 35.8%)
    etalq      ~ 0.030745  # var = log(1 + 0.177^2) = 0.030745 ; Svensson 2014 Supplementary Table S1b BSV Q1 = 17.7% (RSE 20.9%)
    etalvc_m2  ~ 0.008820  # var = log(1 + 0.0943^2) = 0.008820 ; Svensson 2014 Supplementary Table S1b BSV VM2 = 9.43% (RSE 152%)
    etalvp_m2  ~ 0.009704  # var = log(1 + 0.0990^2) = 0.009704 ; Svensson 2014 Supplementary Table S1b BSV VP1M2 = 9.90% (RSE 75.8%)

    # Residual error. The paper reports proportional residual errors on
    # bedaquiline (22.7% CV) and M2 (16.2% CV) plus a 54.2% cross-output
    # correlation between the two residuals and a 2.42-fold amplification
    # of the residual SD for samples drawn during the absorption phase
    # (TAD < 6 h). The cross-output residual correlation and the
    # absorption-phase weighting are dropped here because nlmixr2lib has
    # no idiomatic encoding for either; see vignette Assumptions and
    # deviations for the original paper values.
    propSd     <- 0.227 ; label("Bedaquiline residual error (proportional SD, fraction)")   # Svensson 2014 Supplementary Table S1b 'Prop error TMC = 22.7%' (RSE 8.4%) (TMC = bedaquiline, formerly TMC207)
    propSd_m2  <- 0.162 ; label("M2 metabolite residual error (proportional SD, fraction)") # Svensson 2014 Supplementary Table S1b 'Prop error M2 = 16.2%' (RSE 9.0%)
  })

  model({
    # 1. Allometric scaling on apparent clearances (CL/F, Q1/F, Q2/F,
    #    CL_M2, Q_M2) and apparent volumes (V/F, VP1/F, VP2/F, V_M2,
    #    VP_M2) around a 70 kg reference adult.
    allcl <- (WT / 70)^e_wt_cl_q
    allv  <- (WT / 70)^e_wt_vc_vp

    # 2. Drug-drug-interaction multiplicative factors on bedaquiline and
    #    M2 apparent CL. CONMED_NVP is a time- and subject-varying binary
    #    indicator that the subject has reached full NVP induction. When
    #    CONMED_NVP = 0 (no co-administration), both factors collapse to
    #    1.
    ddi_bdq <- e_nvp_cl    ^ CONMED_NVP
    ddi_m2  <- e_nvp_cl_m2 ^ CONMED_NVP

    # 3. Individual PK parameters (etas applied multiplicatively on the
    #    log scale around the typical-value log-transformed structural
    #    parameter).
    mtt   <- exp(lmtt)
    ka    <- exp(lka)
    nn    <- exp(lnn)

    cl    <- exp(lcl    + etalcl)    * allcl * ddi_bdq
    vc    <- exp(lvc    + etalvc)    * allv
    q     <- exp(lq     + etalq)     * allcl
    vp    <- exp(lvp)                * allv
    q2    <- exp(lq2)                * allcl
    vp2   <- exp(lvp2)               * allv
    cl_m2 <- exp(lcl_m2 + etalcl_m2) * allcl * ddi_m2
    vc_m2 <- exp(lvc_m2 + etalvc_m2) * allv
    q_m2  <- exp(lq_m2)              * allcl
    vp_m2 <- exp(lvp_m2 + etalvp_m2) * allv

    # 4. ODE system. Same parent + metabolite structure as the
    #    Svensson_2014_bedaquiline_lpvr.R sibling (and the upstream
    #    Svensson_2013_bedaquiline.R BDQ + M2 + M3 structure with M3
    #    removed -- C117 did not collect M3 samples).
    d/dt(depot)          <- transit(nn, mtt, 1) - ka * depot
    d/dt(central)        <-  ka * depot -
                              cl  * central / vc -
                              q   * central / vc + q  * peripheral1 / vp -
                              q2  * central / vc + q2 * peripheral2 / vp2
    d/dt(peripheral1)    <-  q   * central / vc - q  * peripheral1 / vp
    d/dt(peripheral2)    <-  q2  * central / vc - q2 * peripheral2 / vp2
    d/dt(central_m2)     <-  cl    * central    / vc -
                              cl_m2 * central_m2 / vc_m2 -
                              q_m2  * central_m2 / vc_m2 + q_m2 * peripheral1_m2 / vp_m2
    d/dt(peripheral1_m2) <-  q_m2  * central_m2 / vc_m2 - q_m2 * peripheral1_m2 / vp_m2

    # 5. Bioavailability. F fixed at 1 implicitly.
    f(depot) <- 0

    # 6. Observed plasma concentrations (mg/L = ug/mL).
    Cc    <- central    / vc
    Cc_m2 <- central_m2 / vc_m2

    Cc    ~ prop(propSd)
    Cc_m2 ~ prop(propSd_m2)
  })
}
