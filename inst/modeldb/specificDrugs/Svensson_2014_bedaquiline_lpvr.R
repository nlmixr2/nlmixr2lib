Svensson_2014_bedaquiline_lpvr <- function() {
  description <- "Three-compartment population PK model for bedaquiline (BDQ) and a two-compartment N-desmethyl metabolite M2 in healthy adult volunteers following single 400 mg oral doses, with Savic 2007 analytical transit-compartment absorption (non-integer NN feeding a first-order depot at rate ka), fixed allometric scaling on disposition (0.75 on CL/Q at 70 kg, 1 on Vc/Vp), and multiplicative ritonavir-boosted lopinavir (LPV/r) drug-drug-interaction factors of 0.347 on bedaquiline apparent clearance and 0.578 on M2 apparent clearance during LPV/r co-administration (study C110)."
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
      notes              = "Time-fixed baseline body weight. Allometric scaling applied with fixed exponents 0.75 on apparent clearances (CL/F, Q1/F, Q2/F, CL_M2, Q_M2) and 1 on apparent volumes (V/F, VP1/F, VP2/F, V_M2, VP_M2) around a 70 kg reference adult. Svensson 2014 Materials and Methods 'Nonlinear mixed-effects modeling' states 'Allometric scaling of disposition parameters with body weight as the size descriptor and fixed coefficients (0.75 for clearance [CL] and 1 for volume of distribution) was applied.' Reference weight of 70 kg is confirmed by Supplementary Table S1a footnote b: 'Disposition parameters for a typical individual of 70 kg, allometric scaling with body weight and fixed coefficients 0.75 for CL and 1 for V applied.'",
      source_name        = "WT"
    ),
    CONMED_LPV = list(
      description        = "Concomitant ritonavir-boosted lopinavir (LPV/r) co-administration (1 = on twice-daily 400/100 mg LPV/r, 0 = not on LPV/r).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not on LPV/r)",
      notes              = "Subject- and time-varying indicator of concomitant LPV/r co-administration. Study C110 was a crossover study (n = 16 healthy seronegative volunteers) in which a single 400 mg bedaquiline dose was given alone in one period and after 10 days of LPV/r 400/100 mg twice daily in the other period; the LPV/r dosing continued throughout the second-period 14-day PK sampling window. The DDI is mediated by ritonavir's potent CYP3A4 inhibition. Svensson 2014 Materials and Methods state 'The impacts of LPV/r (inhibition) were assumed to start immediately upon initiation of administration and to vanish within 1 day after the last LPV/r dose.' Multiplicative factor on apparent CL during co-administration: cl_bdq_eff = cl_bdq_base * e_lpv_cl ^ CONMED_LPV with e_lpv_cl = 0.347 (RSE 9.3%) for bedaquiline (BDQ CL falls to 35% of the no-LPV/r value) and cl_m2_eff = cl_m2_base * e_lpv_cl_m2 ^ CONMED_LPV with e_lpv_cl_m2 = 0.578 (RSE 8.7%) for M2 (M2 CL falls to 58% of the no-LPV/r value). Svensson 2014 Supplementary Table S1a 'EFF1 LPV/r BDQ CL = 0.347' and 'EFF2 LPV/r M2 CL = 0.578'. For simulation, set CONMED_LPV = 1 on observation rows while LPV/r is being given and 0 otherwise.",
      source_name        = "LPV/r"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 16L,
    n_studies      = 1L,
    age_range      = "20-54 years",
    age_median     = "25.5 years",
    weight_range   = "65-103 kg",
    weight_median  = "75 kg",
    sex_female_pct = 6.25,
    race_ethnicity = c(Black = 37.5, White = 62.5, `Mixed race` = 0.0),
    disease_state  = "HIV-seronegative healthy adult volunteers (18-55 years; BMI 18-32 kg/m^2; less than 10 cigarettes/day). Women of childbearing potential, individuals with TB / HIV-1 / HIV-2, and individuals with a history of substance abuse were excluded. Subjects previously enrolled in trials involving bedaquiline were ineligible.",
    dose_range     = "Two single 400 mg oral doses of bedaquiline given 4 weeks apart in a crossover sequence; LPV/r at 400/100 mg twice daily started 10 days before either the first or second bedaquiline dose and continued through the 14-day PK sampling window. PK samples were collected pre-dose and at 1, 2, 3, 4, 5, 6, 8, 12, 24, 48, 72, 120, 168, 216, 264, and 336 h after each bedaquiline dose (17 samples per dose per analyte).",
    regions        = "Not stated in the publication beyond study sponsor (Tibotec / Janssen) and ClinicalTrials.gov registration NCT00828529.",
    study_id       = "C110 (LPV/r DDI; NCT00828529).",
    notes          = "Baseline demographics from Svensson 2014 Table 1. 532 BDQ + 532 M2 PK samples available; 20 M2 samples were below limit of quantification (LLOQ 1.00 ng/mL) and were omitted from modelling. Bedaquiline and M2 concentrations were determined by LC-MS/MS validated to FDA guidelines."
  )

  ini({
    # Final parameter estimates from Svensson 2014 Supplementary Table S1a
    # 'Parameter estimates including precision (RSE) for study C110 (bedaquiline
    # drug-drug interaction with ritonavir-boosted lopinavir)'. The structural
    # model (3-compartment BDQ + 2-compartment M2 with NN transit absorption
    # + first-order ka into the depot) was inherited from Svensson 2013
    # (doi:10.1128/AAC.00191-13); see modellib('Svensson_2013_bedaquiline'). All
    # final parameter values were re-estimated on the C110 data in this paper.

    # Structural absorption: Savic 2007 transit-compartment model (non-integer
    # number of transit compartments NN feeding the depot at rate
    # KTR = (NN+1)/MTT, which then absorbs into central at rate ka).
    lmtt    <- log(1.02)    ; label("Mean transit time MTT (h, on log scale)")                                            # Svensson 2014 Supplementary Table S1a 'MTT = 1.02 h' (RSE 18.6%)
    lka     <- log(0.0983)  ; label("First-order absorption rate constant ka (1/h, on log scale)")                         # Svensson 2014 Supplementary Table S1a 'KA = 0.0983 1/h' (RSE 11.2%)
    lnn     <- log(5.77)    ; label("Number of transit compartments NN (unitless, on log scale)")                          # Svensson 2014 Supplementary Table S1a 'NN = 5.77' (RSE 41.8%)

    # Structural disposition: 3-compartment bedaquiline (parent) +
    # 2-compartment M2 metabolite (formed from BDQ via CYP3A4 N-demethylation).
    # Apparent volumes and clearances are CL/F, V/F for the parent and
    # CL/(F*fm), V/(F*fm) for the metabolite (fm = fraction metabolised to M2,
    # not separately identifiable from F-relative parameters).
    lcl     <- log(3.09)    ; label("Apparent bedaquiline clearance CL/F at 70 kg (L/h)")                                   # Svensson 2014 Supplementary Table S1a 'CL/F = 3.09 L/h' (RSE 17.3%)
    lvc     <- log(16.1)    ; label("Apparent bedaquiline central volume V/F at 70 kg (L)")                                 # Svensson 2014 Supplementary Table S1a 'V/F = 16.1 L' (RSE 22.8%)
    lq      <- log(5.97)    ; label("Apparent bedaquiline first-peripheral inter-compartmental clearance Q1/F at 70 kg (L/h)") # Svensson 2014 Supplementary Table S1a 'Q1/F = 5.97 L/h' (RSE 4.6%)
    lvp     <- log(4890)    ; label("Apparent bedaquiline first-peripheral volume VP1/F at 70 kg (L)")                       # Svensson 2014 Supplementary Table S1a 'VP1/F = 4890 L' (RSE 17%)
    lq2     <- log(3.52)    ; label("Apparent bedaquiline second-peripheral inter-compartmental clearance Q2/F at 70 kg (L/h)") # Svensson 2014 Supplementary Table S1a 'Q2/F = 3.52 L/h' (RSE 9.3%)
    lvp2    <- log(174)     ; label("Apparent bedaquiline second-peripheral volume VP2/F at 70 kg (L)")                      # Svensson 2014 Supplementary Table S1a 'VP2/F = 174 L' (RSE 27.2%)
    lcl_m2  <- log(14.6)    ; label("Apparent M2 metabolite clearance CL_M2/(F*fm) at 70 kg (L/h)")                          # Svensson 2014 Supplementary Table S1a 'CLM2/F/fm = 14.6 L/h' (RSE 16.6%)
    lvc_m2  <- log(746)     ; label("Apparent M2 metabolite central volume V_M2/(F*fm) at 70 kg (L)")                        # Svensson 2014 Supplementary Table S1a 'VM2/F/fm = 746 L' (RSE 25.1%)
    lq_m2   <- log(75.5)    ; label("Apparent M2 metabolite inter-compartmental clearance Q_M2/(F*fm) at 70 kg (L/h)")       # Svensson 2014 Supplementary Table S1a 'Q1M2/F/fm = 75.5 L/h' (RSE 17.7%)
    lvp_m2  <- log(3140)    ; label("Apparent M2 metabolite peripheral volume VP_M2/(F*fm) at 70 kg (L)")                    # Svensson 2014 Supplementary Table S1a 'VP1M2/F/fm = 3140 L' (RSE 19.5%)

    # Bioavailability is implicitly F = 1 because CL and V are reported as
    # apparent F-relative values. Svensson 2014 Materials and Methods
    # 'All disposition parameters were estimated as relative to the bioavailability.'

    # Allometric scaling (theory-based, fixed): exponent 0.75 on apparent
    # clearances and 1 on apparent volumes around 70 kg (Svensson 2014
    # Supplementary Table S1a footnote b).
    e_wt_cl_q  <- fixed(0.75) ; label("Allometric exponent on apparent CL and Q (CL/F, Q1/F, Q2/F, CL_M2, Q_M2; unitless)")  # Svensson 2014 Methods (fixed)
    e_wt_vc_vp <- fixed(1)    ; label("Allometric exponent on apparent Vc and Vp (V/F, VP1/F, VP2/F, V_M2, VP_M2; unitless)") # Svensson 2014 Methods (fixed)

    # Drug-drug-interaction multiplicative factors on the bedaquiline and M2
    # apparent clearances during ritonavir-boosted lopinavir (LPV/r)
    # co-administration. Svensson 2014 Results 'LPV/r was found to decrease
    # the CL of both BDQ and M2 substantially, to 35% (RSE, 9.2%) and 58%
    # (RSE, 8.4%) of CL values without comedication, respectively.' The
    # paper reports two SEPARATE factors (one for BDQ CL, one for M2 CL),
    # unlike the rifampicin paper which used a single shared factor.
    e_lpv_cl    <- 0.347 ; label("Multiplicative factor on bedaquiline apparent CL during LPV/r co-administration (unitless)") # Svensson 2014 Supplementary Table S1a 'EFF1 LPV/r BDQ CL = 0.347' (RSE 9.3%)
    e_lpv_cl_m2 <- 0.578 ; label("Multiplicative factor on M2 apparent CL during LPV/r co-administration (unitless)")          # Svensson 2014 Supplementary Table S1a 'EFF2 LPV/r M2 CL = 0.578' (RSE 8.7%)

    # Inter-individual variability (between-subject, BSV). Final estimates
    # converted from the paper's CV% notation via omega^2 = log(1 + CV^2)
    # for log-normal IIV. Svensson 2014 Supplementary Table S1a reports a
    # 2x2 correlated BSV block on (CL, CL_M2) at +77.4% correlation, plus
    # diagonal BSV terms on V, Q1, V_M2, VP1_M2, and a paired BSV block on
    # the LPV/r interaction effects (BSV EFF1 + scaled BSV EFF2; the paper
    # text confirms 'The BSV values for the interaction effects for BDQ
    # and M2 were 100% correlated but were not significantly correlated
    # with the BSV in the CL of BDQ or M2', so EFF2's eta is encoded as
    # the EFF1 eta scaled by 0.335). nlmixr2lib has no idiomatic encoding
    # for a 100% correlation between two etas where one is a scaled
    # version of the other, nor for between-occasion variability (BOV) on
    # F or MTT, nor for an eta gated by a binary covariate (the LPV/r
    # interaction etas only apply when CONMED_LPV = 1). The encoding here
    # therefore captures the 2x2 BSV(CL, CL_M2) block and the diagonal
    # BSV terms on V, Q1, V_M2, VP1_M2, and DROPS:
    #   - BSV EFF1 = 34.6% (CV) and the scaled BSV EFF2 = 0.335 * EFF1 eta
    #     (Svensson 2014 Supplementary Table S1a 'BSV EFF1, scaled BSV EFF2 = 34.6%'
    #     and 'Scale ETA EFF M2 CL = 0.335')
    #   - BOV F = 13.4% CV (between-occasion variability on bioavailability)
    #   - BSV F = 9.40% CV (between-subject variability on bioavailability)
    #   - BOV MTT = 71.1% CV (between-occasion variability on mean transit time)
    # See the vignette Assumptions and deviations section for the full
    # list of dropped random effects with the original paper values.
    #
    # BLOCK(2) on (CL, CL_M2): var_CL = log(1 + 0.392^2) = 0.142111,
    # var_CL_M2 = log(1 + 0.405^2) = 0.150661,
    # cov = 0.774 * sqrt(0.142111 * 0.150661) = +0.113304
    etalcl + etalcl_m2 ~ c(0.142111,
                            0.113304, 0.150661)  # Svensson 2014 Supplementary Table S1a BSV CL = 39.2% (RSE 18.3%), BSV CLM2 = 40.5% (RSE 20.5%), correlation BSV CL ~ BSV CLM2 = 77.4% (RSE 23.6%)
    etalvc     ~ 0.201115  # var = log(1 + 0.475^2) = 0.201115 ; Svensson 2014 Supplementary Table S1a BSV V = 47.5% (RSE 18.9%)
    etalq      ~ 0.019408  # var = log(1 + 0.140^2) = 0.019408 ; Svensson 2014 Supplementary Table S1a BSV Q1 = 14.0% (RSE 30.7%)
    etalvc_m2  ~ 0.235566  # var = log(1 + 0.520^2) = 0.235566 ; Svensson 2014 Supplementary Table S1a BSV VM2 = 52.0% (RSE 20.6%)
    etalvp_m2  ~ 0.137692  # var = log(1 + 0.385^2) = 0.137692 ; Svensson 2014 Supplementary Table S1a BSV VP1M2 = 38.5% (RSE 40.5%)

    # Residual error. The paper reports proportional residual errors on
    # bedaquiline (17.1% CV) and M2 (14.5% CV) plus a 54.3% cross-output
    # correlation between the two residuals and a 1.89-fold amplification
    # of the residual SD for samples drawn during the absorption phase
    # (TAD < 6 h). The cross-output residual correlation and the
    # absorption-phase weighting are dropped here because nlmixr2lib has
    # no idiomatic encoding for either; see vignette Assumptions and
    # deviations for the original paper values.
    propSd     <- 0.171 ; label("Bedaquiline residual error (proportional SD, fraction)")   # Svensson 2014 Supplementary Table S1a 'Prop error TMC = 17.1%' (RSE 6.9%) (TMC = bedaquiline, formerly TMC207)
    propSd_m2  <- 0.145 ; label("M2 metabolite residual error (proportional SD, fraction)") # Svensson 2014 Supplementary Table S1a 'Prop error M2 = 14.5%' (RSE 6.7%)
  })

  model({
    # 1. Allometric scaling on apparent clearances (CL/F, Q1/F, Q2/F,
    #    CL_M2, Q_M2) and apparent volumes (V/F, VP1/F, VP2/F, V_M2,
    #    VP_M2) around a 70 kg reference adult.
    allcl <- (WT / 70)^e_wt_cl_q
    allv  <- (WT / 70)^e_wt_vc_vp

    # 2. Drug-drug-interaction multiplicative factors on bedaquiline and
    #    M2 apparent CL. CONMED_LPV is a time- and subject-varying binary
    #    indicator that the subject is currently on ritonavir-boosted
    #    lopinavir. When CONMED_LPV = 0 (no co-administration), both
    #    factors collapse to 1.
    ddi_bdq <- e_lpv_cl    ^ CONMED_LPV
    ddi_m2  <- e_lpv_cl_m2 ^ CONMED_LPV

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

    # 4. ODE system. Savic 2007 analytical transit-compartment input
    #    feeds the depot at rate (NN+1)/MTT (via the rxode2 transit()
    #    built-in) and the depot then absorbs into bedaquiline central at
    #    rate ka. Three-compartment bedaquiline disposition (central +
    #    peripheral1 + peripheral2) with linear elimination feeding the
    #    M2 central compartment; M2 has its own two-compartment
    #    disposition (central_m2 + peripheral1_m2) with linear
    #    elimination from central_m2. The transit() function reads dose
    #    amount from podo(depot); f(depot) <- 0 suppresses the bolus
    #    contribution so the chain delivers the full dose smoothly
    #    through the transit compartments.
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

    # 5. Bioavailability. F is fixed at 1 (CL and V reported as
    #    apparent F-relative values). f(depot) <- 0 suppresses the bolus
    #    contribution; transit() drives the input rate over the
    #    absorption window.
    f(depot) <- 0

    # 6. Observed plasma concentrations. Compartmental amount (mg) /
    #    volume (L) = concentration in mg/L (= ug/mL).
    Cc    <- central    / vc
    Cc_m2 <- central_m2 / vc_m2

    Cc    ~ prop(propSd)
    Cc_m2 ~ prop(propSd_m2)
  })
}
