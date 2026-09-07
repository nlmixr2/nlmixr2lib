Bausch_2024_cefazolin <- function() {
  description <- "One-compartment joint total/unbound population PK model for intravenous cefazolin in adults with invasive Staphylococcus aureus infection (Bausch 2024). The disposition is carried on the UNBOUND concentration Cu = central/V with linear unbound clearance; the measured TOTAL plasma concentration is reconstructed as Cc = Cu + Cb, where the bound concentration Cb = Bmax * Cu / (kd + Cu) + NS * Cu combines a saturable albumin-binding site with a linear non-saturable component. Total and unbound concentrations are BOTH observed endpoints, each with its own proportional residual error. Clearance carries an allometric Cockcroft-Gault creatinine-clearance term, volume an allometric weight term with the exponent FIXED at 1, and the non-saturable binding constant an allometric serum-albumin term. Inter-individual variability on V and Bmax; inter-occasion variability on CL and NS across the four therapeutic-drug-monitoring sampling days."
  reference <- "Bausch S, Draeger S, Charitos-Fragkakis P, Egli A, Moser S, Hinic V, Kuehl R, Bassetti S, Siegemund M, Rentsch KM, Hermann L, Schoening V, Hammann F, Sendi P, Osthoff M. Target Attainment and Population Pharmacokinetics of Cefazolin in Patients with Invasive Staphylococcus aureus Infections: A Prospective Cohort Study. Antibiotics (Basel). 2024;13(10):928. doi:10.3390/antibiotics13100928. PMCID PMC11504871. Binding equations from the Results Sect. 2.5 display equation; all parameter estimates and covariate relationships from Supplementary Table S2; model-building sequence from Supplementary Table S3."
  vignette <- "Bausch_2024_cefazolin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central = list(
      analyte = "cefazolin", units = "mg", specimen = "plasma", verified = TRUE,
      notes = "Holds the administered cefazolin amount. Dividing by V yields the UNBOUND plasma concentration Cu, not the total: V is the unbound-referenced volume, and CL the unbound-referenced clearance. The measured total concentration is the algebraic sum Cu + Cb. This reading is what makes the parameter magnitudes physiological -- a total-referenced cefazolin volume is roughly 10 L and a total-referenced clearance roughly 4 L/h, whereas Table S2 reports 67.2 L and 15.6 L/h, i.e. both inflated by about the reciprocal of the 27% unbound fraction the paper reports. It also reproduces all four of the paper's own Table 2 concentration medians (see the vignette)."
    )
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance estimated by the Cockcroft-Gault equation (eGFR-CG), raw mL/min and NOT BSA-normalized",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 80.1634 mL/min, the weighted mean eGFR-CG of the cohort, given verbatim in the Bausch 2024 Supplementary Table S2 legend. Enters as the power term (CRCL / 80.1634)^1.2 on unbound clearance; the exponent 1.2 is estimated (5.0% RSE, bootstrap 95% CI 0.9-1.4), so the renal effect is slightly MORE than proportional. Cockcroft-Gault returns mL/min and is not indexed to body surface area, so this column is a raw clearance and is NOT interchangeable with the BSA-normalized eGFR reported in Bausch 2024 Table 1 (median 76 mL/min/1.73 m^2, IQR 41-91), which was computed by CKD-EPI. The paper tested both equations and Supplementary Table S3 records that Cockcroft-Gault won on objective function value (2833.13 vs 2845.04 for CKD-EPI in the joint fu-parameterized model), which is why eGFR-CG and not eGFR-CKD-EPI is the covariate carried here. Patients on renal replacement therapy were excluded by design and none required it during the study, so the model carries no information about dialysis. The paper's own Monte Carlo simulations (Figure 2) extrapolate this term across 20-100 mL/min.",
      source_name        = "eGFR-CG"
    ),
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 70 kg, printed as the literal denominator of the Bausch 2024 Supplementary Table S2 volume covariate equation (it is a rounded standard, not the cohort median of 73 kg from Table 1). Enters as the power term (WT / 70)^1.0 on unbound volume of distribution with the exponent FIXED at 1.0 -- Table S2 marks the row 'Weight_V (fixed)' and reports no RSE and no bootstrap interval for it -- so the effect is simple linear proportionality, not the 0.75/1.0 allometric pair. Cohort median 73 kg (IQR 67-93); BMI median 24.5 kg/m^2 (IQR 21.8-29.6), so the model was not fitted in an obese cohort.",
      source_name        = "Weight"
    ),
    ALB = list(
      description        = "Serum albumin concentration measured in the same blood sample as the cefazolin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 24.0844 g/L, the weighted mean albumin, given verbatim in the Bausch 2024 Supplementary Table S2 legend and re-used as the fixed albumin of the paper's own Monte Carlo virtual patients (Methods Sect. 4.7). Enters as the power term (ALB / 24.0844)^2.4 on the NON-SATURABLE binding constant NS. The exponent is steep (2.4, 12.0% RSE, bootstrap 95% CI 1.8-3.9), so albumin drives the unbound fraction strongly, which is the paper's central clinical message. Note the reference is markedly below a healthy 35-50 g/L: the cohort median albumin at onset of infection was 29.0 g/L (Table 1) and 25 g/L at first drug measurement (Table S1), i.e. this is a hypoalbuminaemic acute-infection population, and the weighted mean is pulled further down by the repeated sampling of the sickest patients. Albumin was measured at every visit but the authors explicitly chose NOT to treat it as time-varying (Methods Sect. 4.7: the parameters 'did not change substantially across occasions', so only the first measurement per patient and occasion was used).",
      source_name        = "Albumin"
    ),
    OCC = list(
      description        = "Therapeutic-drug-monitoring sampling occasion",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Values 1-4 identify the four protocol sampling days on which cefazolin concentrations were drawn: study day 1, day 3 (+/- 1), day 7 (+/- 2) and day 14 (+/- 5) (Bausch 2024 Methods Sect. 4.4 and Supplementary Figure S6). The paper introduced inter-occasion variability precisely because 'blood concentrations [were] measured in different inter-dose intervals' (Methods Sect. 4.7) and fits ONE shared IOV magnitude per parameter across occasions rather than a per-occasion set, so the encoding below carries the estimated variance on occasion 1 and fixes occasions 2-4 to the same value (the registered `$OMEGA BLOCK(1) SAME`-equivalent idiom). Unlike Stoschus 2025 and Ding 2026, where the occasion count had to be inferred, the count of four is set by the paper's own protocol schedule. For single-occasion records pass OCC = 1 so the first IOV eta applies.",
      source_name        = "occasion"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 51,
    n_studies      = 1,
    age_range      = "median 74.1 years (IQR 56.6-81.8)",
    weight_range   = "median 73 kg (IQR 67-93)",
    sex_female_pct = 25.5,
    race_ethnicity = "Not reported (single-centre Swiss tertiary-care cohort)",
    disease_state  = "Adults with invasive methicillin-susceptible Staphylococcus aureus infection (76.5% bloodstream infection; foci included endocarditis 19.6%, osteomyelitis or septic arthritis 19.6%, catheter or prosthetic material 25.5%). Median Charlson comorbidity score 5, median SOFA score 1, 19.6% ICU admission, 30-day mortality 5.9%.",
    renal_function = "Median eGFR (CKD-EPI) 76 mL/min/1.73 m^2 (IQR 41-91) at onset of infection; chronic kidney disease stage G3 in 9.8% and G4 in 9.8%; acute kidney injury in 11.8%. Haemodialysis was an exclusion criterion and no patient required renal replacement therapy during the study.",
    dose_range     = "Cefazolin 2 g as a 30-minute intermittent bolus infusion every 6 or 8 h (2 g q12h if eGFR 10-30 and 2 g q24h if eGFR < 10 mL/min/1.73 m^2); 92% of patients received 6 g/day and 8% 4 g/day",
    regions        = "Switzerland (University Hospital Basel, single centre, January 2020 - December 2021; ClinicalTrials.gov NCT04503252)",
    notes          = "226 paired total and unbound cefazolin plasma concentrations (mean 4.4 per patient) measured by HPLC-MS/MS with ultracentrifugation for the unbound fraction; LLOQ 0.5 mg/L for both, assay linear 0.5-100.0 mg/L, and no observation fell below the LLOQ. Median serum albumin 29.0 g/L (IQR 24.0-32.8) at onset of infection. The mean unbound fraction was 27.0% (SD 13.4, range 8.9-79.7%), far above and far more variable than the 20% conventionally assumed. Fitted in Monolix 2023R1 and validated by a 1000-replicate non-parametric bootstrap (Supplementary Table S2). Note the population is elderly, hypoalbuminaemic and predominantly male, and the observed S. aureus MIC distribution was narrow (median 1.0 mg/L, IQR 0.75-1.0, maximum 1.5 mg/L) and from a single geographical area."
  )

  ini({
    # ----------------------------------------------------------------
    # Structural disposition parameters (Bausch 2024 Supplementary
    # Table S2, "Final joint model" column; estimate [%RSE]). Both are
    # UNBOUND-referenced: elimination acts on Cu, and Cu = central / V.
    # See compartmentData$central for why the unbound reading is the
    # only one consistent with the magnitudes.
    # ----------------------------------------------------------------
    lcl <- log(15.6); label("Unbound plasma clearance CL (L/h)")               # Table S2 CL_pop = 15.6 l h-1 [3.8% RSE]; bootstrap median 16.2, 95% CI 13.4-19.1
    lvc <- log(67.2); label("Unbound volume of distribution V (L)")            # Table S2 V_pop = 67.2 l [7.3% RSE]; bootstrap median 71.2, 95% CI 55.7-90.1

    # ----------------------------------------------------------------
    # Albumin-binding isotherm constants (Bausch 2024 Results Sect. 2.5
    # display equation). Bmax and kd are in concentration units (mg/L,
    # the units of every concentration in the paper); NS is the
    # dimensionless slope of the linear non-saturable component.
    # Table S2 tabulates all three without unit labels; the mg/L
    # reading is confirmed by the reproduction of the paper's own
    # Table 2 concentration medians documented in the vignette.
    # ----------------------------------------------------------------
    lbmax_pb <- log(55.2); label("Maximum saturable binding capacity Bmax (mg/L)")            # Table S2 Bmax_pop = 55.2 [10.7% RSE]; bootstrap median 54.1, 95% CI 33.5-80.1
    lkd_pb   <- log(12.7); label("Dissociation constant kd for the saturable site (mg/L)")    # Table S2 kd_pop = 12.7 [10.1% RSE]; bootstrap median 12.9, 95% CI 7.1-21.2
    lns_pb   <- log(0.5);  label("Slope of the linear non-saturable binding component (unitless)") # Table S2 NS_pop = 0.5 [18.6% RSE]; bootstrap median 0.6, 95% CI 0.4-1.0

    # ----------------------------------------------------------------
    # Covariate effects (Bausch 2024 Supplementary Table S2, "Covariate
    # Relationships" block). All three are allometric power terms on a
    # normalized covariate:
    #   CL = CL_pop * (GFR      / GFR_mean)^GFR_CL
    #   V  = V_pop  * (Weight   / 70)^Weight_V
    #   NS = NS_pop * (Albumin  / Albumin_mean)^Albumin_NS
    # with GFR_mean = 80.1634 and Albumin_mean = 24.0844 given in the
    # Table S2 legend.
    #
    # Albumin acts on NS, NOT on kd. The Table S2 legend glosses
    # Albumin_NS as "exponent for the allometrically scaled albumin on
    # kd", which contradicts both (a) the Covariate Relationships block
    # of the same table, whose third row is headed "NS" and prints
    # NS_pop * (Albumin/Albumin_mean)^Albumin_NS, and (b) Supplementary
    # Table S3, whose winning row in the final "Protein binding
    # non-linear (kd), saturable (Bmax) and unsaturable (NS)" block is
    # "GFR-CG on Cl / Albumin on NS / Weight on Vd" at OFV 2816.93 --
    # better than the "Albumin on kd" alternative at 2854.78. The
    # legend gloss is a leftover from the earlier kd-only variant that
    # Table S3 shows the authors tried and rejected.
    # ----------------------------------------------------------------
    e_crcl_cl <- 1.2;        label("Exponent on (CRCL / 80.1634) for unbound clearance")            # Table S2 GFR_CL = 1.2 [5.0% RSE]; bootstrap median 1.2, 95% CI 0.9-1.4
    e_wt_vc   <- fixed(1.0); label("Exponent on (WT / 70) for unbound volume of distribution")      # Table S2 "Weight_V (fixed)" = 1.0; no RSE, bootstrap n/a
    e_alb_ns  <- 2.4;        label("Exponent on (ALB / 24.0844) for the non-saturable constant NS") # Table S2 Albumin_NS = 2.4 [12.0% RSE]; bootstrap median 2.8, 95% CI 1.8-3.9

    # ----------------------------------------------------------------
    # Inter-individual variability (Bausch 2024 Supplementary Table S2,
    # "Inter-individual variability (IIV)" block). Only V and Bmax
    # carry IIV; CL and NS carry inter-occasion variability instead
    # (next block), and kd carries neither.
    #
    # The fit was run in Monolix (Methods Sect. 4.7), whose population-
    # parameter table reports random effects as the STANDARD DEVIATION
    # omega of the log-scale random effect, so the nlmixr2 variances
    # below are the squares of the tabulated numbers:
    #   V_IIV    = 0.4 -> variance 0.16 (CV 41.6%)
    #   Bmax_IIV = 0.2 -> variance 0.04 (CV 20.2%)
    # ----------------------------------------------------------------
    etalvc      ~ 0.16  # Table S2 V_IIV = 0.4 [15.6% RSE], bootstrap 95% CI 0.2-0.6; variance = 0.4^2
    etalbmax_pb ~ 0.04  # Table S2 Bmax_IIV = 0.2 [21.0% RSE], bootstrap 95% CI 0.1-0.3; variance = 0.2^2

    # ----------------------------------------------------------------
    # Inter-occasion variability (Bausch 2024 Supplementary Table S2,
    # "Inter-occasional variability (IOV)" block), on CL and NS. As
    # with IIV the tabulated numbers are Monolix standard deviations,
    # so the variances are their squares:
    #   CL_IOV = 0.4 -> variance 0.16
    #   NS_IOV = 0.2 -> variance 0.04
    #
    # The paper reports ONE shared magnitude per parameter rather than
    # a per-occasion set, so occasion 1 carries the estimated variance
    # and occasions 2-4 fix it to the same value -- the registered
    # `$OMEGA BLOCK(1) SAME`-equivalent idiom (see
    # Jonsson_2011_ethambutol.R, Stoschus_2025_phenobarbital.R,
    # Ding_2026_vancomycin.R).
    # ----------------------------------------------------------------
    etaiov_cl_1 ~ 0.16        # Table S2 CL_IOV = 0.4 [7.12% RSE], bootstrap 95% CI 0.3-0.5; variance = 0.4^2 (estimated)
    etaiov_cl_2 ~ fixed(0.16) # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_cl_3 ~ fixed(0.16) # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_cl_4 ~ fixed(0.16) # SAME-equivalent: equal to the occasion-1 IOV variance

    etaiov_ns_1 ~ 0.04        # Table S2 NS_IOV = 0.2 [31.0% RSE], bootstrap 95% CI 0.1-0.3; variance = 0.2^2 (estimated)
    etaiov_ns_2 ~ fixed(0.04) # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_ns_3 ~ fixed(0.04) # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_ns_4 ~ fixed(0.04) # SAME-equivalent: equal to the occasion-1 IOV variance

    # ----------------------------------------------------------------
    # Residual error (Bausch 2024 Supplementary Table S2, "Residual
    # error" block). The joint model observes BOTH endpoints, so there
    # are two error terms, each purely proportional. Table S2 heads the
    # block "Proportional (b)" and lists no additive term for either
    # endpoint. Cc (total) is the parent observation and takes the
    # suffix-free name; Cu (unbound) takes the per-output suffix.
    # ----------------------------------------------------------------
    propSd    <- 0.2; label("Proportional residual SD on total concentration Cc (fraction)")    # Table S2 Residual error, Proportional (b), Total = 0.2 [8.33% RSE]; bootstrap 95% CI 0.1-0.2
    propSd_Cu <- 0.2; label("Proportional residual SD on unbound concentration Cu (fraction)")  # Table S2 Residual error, Proportional (b), Unbound = 0.2 [7.79% RSE]; bootstrap 95% CI 0.2-0.3
  })

  model({
    # 1. Decompose the integer occasion column into binary indicators to
    #    multiplex the per-occasion IOV etas. For single-occasion data
    #    pass OCC = 1 so the first IOV eta applies.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)

    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 + oc4 * etaiov_cl_4
    iov_ns <- oc1 * etaiov_ns_1 + oc2 * etaiov_ns_2 + oc3 * etaiov_ns_3 + oc4 * etaiov_ns_4

    # 2. Individual PK parameters. CL carries the Cockcroft-Gault
    #    creatinine-clearance power term and its occasion-specific IOV
    #    eta but no IIV eta; V carries the weight power term (exponent
    #    fixed at 1) and its IIV eta.
    cl <- exp(lcl + iov_cl) * (CRCL / 80.1634)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    # 3. Albumin-binding constants. Bmax carries IIV; NS carries the
    #    albumin power term and its occasion-specific IOV eta; kd is a
    #    typical value with no random effect and no covariate.
    bmax_pb <- exp(lbmax_pb + etalbmax_pb)
    kd_pb   <- exp(lkd_pb)
    ns_pb   <- exp(lns_pb + iov_ns) * (ALB / 24.0844)^e_alb_ns

    # 4. Micro-constant. Elimination acts on the unbound concentration,
    #    CL * Cu = (CL / V) * central, so the amount kinetics stay
    #    first-order; the binding non-linearity sits entirely in the
    #    observation equations below.
    kel <- cl / vc

    # 5. One-compartment ODE system. Cefazolin is given intravenously,
    #    so doses go straight to `central` (as a 30-minute infusion for
    #    the intermittent regimens, or as a constant-rate infusion for
    #    the continuous ones).
    d/dt(central) <- -kel * central

    # 6. Observations (Bausch 2024 Results Sect. 2.5 display equation):
    #       Cbound = Bmax * Cunbound / (kd + Cunbound) + NS * Cunbound
    #       Ctotal = Cunbound + Cbound
    #    Both Cc and Cu are measured endpoints in the joint model, each
    #    with its own proportional residual error.
    Cu <- central / vc
    Cb <- bmax_pb * Cu / (kd_pb + Cu) + ns_pb * Cu
    Cc <- Cu + Cb

    Cc ~ prop(propSd)
    Cu ~ prop(propSd_Cu)
  })
}
