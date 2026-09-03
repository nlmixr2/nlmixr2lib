Gafar_2026_rifampicin <- function() {
  description <- "Semi-mechanistic one-compartment population PK model for standard- and high-dose oral rifampicin given as tuberculosis preventive therapy (10, 20 or 30 mg/kg once daily) to adolescents and adults with tuberculosis infection in the 2R2 randomized trial (Canada, Indonesia, Vietnam). Absorption is a Savic analytical transit chain (17.9 transit compartments, mean transit time 0.768 h) feeding a first-order absorption compartment that empties into a liver compartment, so the drug is subject to first-pass metabolism before reaching the systemic circulation and recirculates to the liver with hepatic blood flow. Elimination is a well-stirred hepatic model with saturable intrinsic clearance: CLint = CLint,max * Km / (CH + Km), EH = CLint * fu / (CLint * fu + QH), which produces the greater-than-proportional rise in exposure observed at 20 and 30 mg/kg. All four disposition parameters (CLint,max, V, QH, VH) are allometrically scaled by fat-free mass, with exponents fixed at 0.75 for the clearance terms and 1 for the volume terms; CLint,max and V are referenced to the cohort median fat-free mass of 41 kg and QH and VH to the 56 kg fat-free mass of a 70 kg reference adult. Prehepatic bioavailability is fixed to 1 in Indonesia (Indofarma capsules) and reduced by 21.8% in Canada (Rofact) and 12.3% in Vietnam (Svizera / Mekophar). Between-subject variability is carried on CLint,max, V and the Michaelis-Menten constant; between-occasion variability on prehepatic bioavailability, the absorption rate constant and the mean transit time; residual error is additive plus a proportional term whose magnitude switches between sparse and intensive PK sampling."
  reference <- paste(
    "Gafar F., Svensson E. M., Yunivita V., Fregonese F., Fisher D.,",
    "Fox G. J., Nguyen T. A., Nguyen B. H., Johnston J., Long R.,",
    "Valiquette C., Aarnoutse R. E., Ruslami R., Menzies D. (2026).",
    "Population Pharmacokinetic Modeling of Standard- and High-Dose",
    "Rifampicin for Tuberculosis Preventive Therapy in the 2R2",
    "Randomized Controlled Trial.",
    "The Journal of Infectious Diseases 233(4):e999-e1010.",
    "doi:10.1093/infdis/jiag052.",
    "Parameter estimates from Table 2; model equations from the Figure 1",
    "caption and from the NONMEM control stream reproduced verbatim in",
    "Supplementary Appendix 1.",
    "The saturable-hepatic-extraction structure was adapted from",
    "Chirehwa et al. (2016) Antimicrob Agents Chemother 60(1):487-494",
    "doi:10.1128/AAC.01084-15, with autoinduction omitted because",
    "sampling was performed only once at week 4.",
    "Fat-free mass follows Janmahasatian et al. (2005)",
    "Clin Pharmacokinet 44(10):1051-1065 doi:10.2165/00003088-200544100-00004.",
    sep = " "
  )
  vignette <- "Gafar_2026_rifampicin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass, computed from total body weight, height and sex with the Janmahasatian formula.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Janmahasatian formula, reproduced in the Table 1 footnote and the",
        "Table S4 footnote of the source paper:",
        "males FFM = 9270 * WT / (6680 + 216 * BMI);",
        "females FFM = 9270 * WT / (8780 + 244 * BMI); BMI = WT / HT^2 with",
        "WT in kg and HT in m. Allometric scaling by FFM was preferred over",
        "total body weight (199-point OFV drop versus no allometry, and a",
        "further 100-point drop versus scaling by total body weight;",
        "Results 'Covariate Model' paragraph 1). Two different reference",
        "values are used in the same model, exactly as in the control stream:",
        "CLint,max and V are normalised to the cohort median FFM of 41 kg",
        "(ALLMCL = (FFM/41)**0.75, ALLMV = (FFM/41)**1) whereas hepatic blood",
        "flow QH and liver volume VH are normalised to 56 kg, the FFM of the",
        "70 kg reference adult from which the fixed 90 L/h and 1 L values are",
        "taken (QH = THETA(7)*(FFM/56)**0.75, VH = THETA(9)*(FFM/56)).",
        "Cohort median 41 kg, IQR 35.3-49.2 kg (Table 1); per-country medians",
        "51.3 kg Canada, 40.1 kg Indonesia, 38.0 kg Vietnam (Table S7).",
        "Source column FFM."
      ),
      source_name        = "FFM"
    ),
    REGION_CANADA = list(
      description        = "1 = participant enrolled at a Canadian study site (Calgary, Edmonton, Montreal, Vancouver), 0 = otherwise. Reference country is Indonesia.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Indonesia, the reference country; Vietnam is carried by REGION_VIETNAM)",
      notes              = paste(
        "Country and rifampicin formulation are completely confounded in the",
        "2R2 trial, so the source paper labels this effect a",
        "'country-specific formulation effect' (Table 2). Every Canadian",
        "participant received re-compounded 300/450 mg Rofact capsules",
        "(Bausch Health, re-compounded by Pharmacie Linda Frayne, Montreal),",
        "every Indonesian participant received Indofarma capsules and every",
        "Vietnamese participant received Svizera or Mekophar capsules",
        "(Methods 'Study Procedures and Interventions'). In the control",
        "stream the driver column is FRM (formulation), decomposed as",
        "IF(FRM.EQ.1) FRM1 = 1 (Canada / Rofact) and",
        "IF(FRM.EQ.2.OR.FRM.EQ.3) FRM2 = 1 (Vietnam / Svizera or Mekophar),",
        "with the remaining level (Indonesia / Indofarma) as the reference;",
        "the two Vietnamese formulations were also fitted separately and did",
        "not change the result (dOFV 1 point, 3 df, P > .05), so the retained",
        "effect is a country-level contrast. Enters multiplicatively on",
        "prehepatic bioavailability as THETA(11)**FRM1. Source column FRM",
        "(value 1)."
      ),
      source_name        = "FRM"
    ),
    REGION_VIETNAM = list(
      description        = "1 = participant enrolled at a Vietnamese study site (Hanoi, Ho Chi Minh City), 0 = otherwise. Reference country is Indonesia.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Indonesia, the reference country; Canada is carried by REGION_CANADA)",
      notes              = paste(
        "See REGION_CANADA for the country / formulation confounding. Enters",
        "multiplicatively on prehepatic bioavailability as THETA(12)**FRM2,",
        "where FRM2 = 1 when the source column FRM is 2 (Svizera) or 3",
        "(Mekophar). Source column FRM (values 2 and 3 pooled)."
      ),
      source_name        = "FRM"
    ),
    OCC = list(
      description        = "Integer occasion indicator for between-occasion variability. Two occasions in the source analysis.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Between-occasion variability was estimated over 2 consecutive doses",
        "with at least 1 postdose observation in the same individual",
        "(Methods 'Stochastic Model'). The control stream $PK block builds",
        "each occasion effect as IF (OCC.EQ.1) THEN ... ELSE ... over a",
        "single $OMEGA BLOCK(1) plus SAME pair, i.e. exactly two occasions",
        "sharing one variance. Decomposed inside model() into oc1 and oc2 and",
        "multiplexed onto the bioavailability, absorption-rate and",
        "mean-transit-time IOV etas. For a single-occasion simulation pass",
        "OCC = 1 so the first IOV eta applies. Source column OCC."
      ),
      source_name        = "OCC"
    ),
    SAMPLE_INTENSIVE = list(
      description        = "1 = the observation belongs to the intensive PK substudy (serial sampling at 0, 1, 2, 4, 8, 12 h after directly observed dosing); 0 = sparse PK substudy (2 and 4 h after self-administered dosing at home).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (sparse sampling)",
      notes              = paste(
        "Record-level switch between the two proportional residual-error",
        "magnitudes. Stratifying by sampling strategy dropped the objective",
        "function value by 23 points (1 df, P < .001) versus a single",
        "proportional error (Results 'Stochastic Model'). Of the 440 modelled",
        "participants, 51 (11.6%) contributed intensive profiles (all from",
        "Bandung, Indonesia) and 389 (88.4%) contributed sparse samples",
        "(Table 1). Control stream: IF(ITS.EQ.1) PROP = IPRED*THETA(14).",
        "Source column ITS."
      ),
      source_name        = "ITS"
    )
  )

  # Covariates the source paper screened but did not retain in the final
  # model. Documented here so the provenance of the covariate search survives
  # without declaring covariates that model() never references.
  covariatesDataExcluded <- list(
    WT = list(
      description        = "Total body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling by total body weight was tested and rejected in favour of fat-free mass (100-point higher OFV; Results 'Covariate Model'). WT is still needed upstream of the model as an input to the Janmahasatian FFM formula. Cohort median 60.0 kg, IQR 51.0-71.0 (Table 1)."
    ),
    HT = list(
      description        = "Body height at baseline.",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Not a model covariate; an input to the Janmahasatian FFM formula via BMI. Reported in metres in the source (median 1.59 m, IQR 1.52-1.66; Table 1); the canonical column is in cm."
    ),
    SEXF = list(
      description        = "1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Screened and not retained: 'Sex was not an independent predictor of exposure, as its effects were already accounted for by fat-free mass' (Results 'Covariate Model'). Still needed upstream as an input to the sex-specific Janmahasatian FFM formula. 255/440 (57.9%) of the modelled cohort were female (Table 1)."
    ),
    AGE = list(
      description        = "Age at baseline.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on absorption and disposition parameters and not retained (Results 'Covariate Model'). Cohort median 40 years, IQR 26-50; 32/440 (7.3%) were adolescents aged 10-17 years (Table 1)."
    ),
    BMI = list(
      description        = "Body mass index at baseline.",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters the Janmahasatian FFM formula but is not itself a model covariate. Nutritional status, the categorical BMI / BMI-for-age-Z-score derivative defined in the Table 1 footnote, was screened and not retained. Cohort median 24.1 kg/m^2, IQR 20.7-27.8 (Table 1)."
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "Rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    liver   = list(analyte = "Rifampicin", units = "mg", specimen = "tissue", verified = TRUE),
    central = list(analyte = "Rifampicin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 440L,
    n_studies      = 1L,
    age_range      = "10 years and older (eligibility); modelled cohort median 40 years, IQR 26-50; 32 (7.3%) adolescents aged 10-17 years and 408 (92.7%) adults aged 18 years or older",
    weight_range   = "at least 25 kg (eligibility); modelled cohort median 60.0 kg, IQR 51.0-71.0 kg",
    ffm_range      = "median 41 kg, IQR 35.3-49.2 kg (Table 1); per country 51.3 kg Canada, 40.1 kg Indonesia, 38.0 kg Vietnam (Table S7)",
    sex_female_pct = 57.9,
    race_ethnicity = "not reported; participants were enrolled in Canada, Indonesia and Vietnam",
    disease_state  = "Tuberculosis infection (positive tuberculin skin test or interferon-gamma release assay with tuberculosis disease excluded) with an indication for tuberculosis preventive therapy. Generally healthy, in contrast with the tuberculosis-disease cohorts of most published rifampicin popPK models. 14 (3.2%) were HIV-1 positive and received lamivudine/tenofovir/efavirenz with efavirenz 800 mg once daily; 12 (2.7%) had diabetes mellitus; 63 (14.3%) took comedications with a potential drug-drug interaction.",
    dose_range     = "Oral rifampicin once daily, dispensed by the prespecified weight bands of Table S3 (25-35, >35-55 and >55 kg) targeting 10 mg/kg/day for 120 days (4R10; 300/450/600 mg) or 20 or 30 mg/kg/day for 60 days (2R20 600/900/1200 mg; 2R30 900/1350/1800 mg). Treatment arms for the popPK analysis were redefined by the actual mg/kg dose received: 5.1-15.0 mg/kg (n = 191), 15.1-25.0 mg/kg (n = 159) and 25.1-35.0 mg/kg (n = 90). Cohort median actual dose 17.1 mg/kg, IQR 9.4-23.2.",
    regions        = "Canada (Calgary, Edmonton, Montreal, Vancouver; n = 87, 19.8%), Indonesia (Bandung; n = 265, 60.2%) and Vietnam (Hanoi, Ho Chi Minh City; n = 88, 20.0%)",
    notes          = paste(
      "Two pharmacokinetic substudies nested in the 2R2 phase 2b randomized",
      "trial (NCT03988933). 51 participants from Bandung, Indonesia",
      "contributed intensive profiles (0, 1, 2, 4, 8 and 12 h after directly",
      "observed dosing at the research clinic) and 389 participants from all",
      "sites contributed sparse samples (2 and 4 h after self-administered",
      "dosing at home). Sampling was performed at approximately 4 weeks of",
      "treatment (median 30 days, IQR 27-36), assumed to be steady state and",
      "past the completion of rifampicin autoinduction, which is why the",
      "reference model's autoinduction component was omitted. Participants",
      "fasted 4-8 h before dosing and until 2 h postdose; 17 (3.9%) had a",
      "light meal. Of 460 participants enrolled in the substudies, 18 with",
      "undetectable concentrations at both 2 and 4 h (likely missed dosing)",
      "were excluded and 11 outlier concentrations censored, leaving 1041",
      "concentrations from 440 participants, of which 998 (96%) were above",
      "the 0.125 mg/L lower limit of quantification. Baseline characteristics",
      "are Table 1; Table S4 gives the same table by randomized rather than",
      "actual dose group, and Table S7 the per-country breakdown."
    )
  )

  ini({
    # =======================================================================
    # Structural parameters. Point estimates are Table 2 of Gafar 2026,
    # cross-checked against the $THETA block of the NONMEM control stream in
    # Supplementary Appendix 1 (whose "initial estimates" are the converged
    # final values -- every one of THETA(1), (2), (3), (5), (6), (10), (11)
    # and (12) reproduces the corresponding Table 2 row). Where the control
    # stream carries an extra significant figure it is used here.
    # Time unit h, dose mg, concentration mg/L.
    # =======================================================================
    lclint_max <- log(46.7)   ; label("Maximum intrinsic hepatic clearance CLint,max (L/h) at the median fat-free mass of 41 kg")     # Table 2 CLint,max = 46.7 L/h (RSE 4.4%, 90% CI 42.1-48.9); control stream $THETA (0, 46.7) ;[CL]
    lvc        <- log(22.8)   ; label("Central volume of distribution V (L) at the median fat-free mass of 41 kg")                    # Table 2 V = 22.8 L (RSE 2.6%, 90% CI 21.4-23.4); control stream $THETA (0, 22.8) ;[V]
    lka        <- log(4.38)   ; label("First-order absorption rate constant Ka (1/h) out of the absorption compartment")              # Table 2 Ka = 4.4 /h (RSE 19.8%, 90% CI 3.7-5.3); control stream $THETA (0, 4.38) ;[KA] carries the extra digit
    lmtt       <- log(0.768)  ; label("Mean transit time MTT (h) of the absorption transit chain")                                    # Table 2 MTT = 0.8 h (RSE 9.7%, 90% CI 0.7-0.9); control stream $THETA (0, 0.768) ;[MTT] carries the extra digit
    lnn        <- log(17.9)   ; label("Number of transit compartments NN (unitless, non-integer)")                                    # Table 2 NN = 17.9 (RSE 28.4%, 90% CI 12.9-23.1); control stream $THETA (0, 17.9) ;[NN]
    lfdepot    <- fixed(log(1))  ; label("Typical prehepatic bioavailability F in the reference country Indonesia (unitless)")        # Table 2 'Prehepatic bioavailability, F: 1 fixed'; control stream $THETA (1) FIX ;[BIO]

    # Km was estimated on the log scale to stabilise the model (Table 2
    # footnote d), so THETA(10) IS the log of Km and needs no log() wrapper
    # here: exp(4.16) = 64.07 mg/L, the 64.1 mg/L of Table 2.
    lkm        <- 4.16        ; label("Michaelis-Menten constant Km (mg/L, on the log scale)")                                        # Table 2 Km = 64.1 mg/L (RSE 15.8%, 90% CI 52.3-98.2); control stream $THETA (0, 4.16) ;[LOGKM], exp(4.16) = 64.07

    # ----- Fixed physiological constants (Results 'Parameter Estimates') ---
    lqh        <- fixed(log(90)) ; label("Hepatic blood flow QH (L/h) for a 70 kg reference adult of 56 kg fat-free mass")            # Table 2 QH = 90 L/h fixed, citing Abdelgawad 2025 and Yang 2007; control stream $THETA (90) FIX ;[QH]
    lvh        <- fixed(log(1))  ; label("Liver volume VH (L) for a 70 kg reference adult of 56 kg fat-free mass")                    # Table 2 VH = 1 L fixed; control stream $THETA (1) FIX ;[VH]
    fub        <- fixed(0.2)     ; label("Fraction of rifampicin unbound in blood fu (unitless)")                                     # Table 2 fu = 0.2 fixed, citing Acocella 1978; control stream $THETA (0.2) FIX ;[FU]

    # ----- Allometric exponents, fixed a priori ---------------------------
    # Methods 'Covariate Model': "with power exponents fixed at 0.75 for
    # clearance and 1 for volume parameters". Hardcoded in the control
    # stream's $PK block rather than estimated as THETAs.
    e_ffm_clint_max <- fixed(0.75) ; label("Allometric exponent of fat-free mass on CLint,max")                                       # control stream ALLMCL = (FFM/41)**0.75
    e_ffm_vc        <- fixed(1)    ; label("Allometric exponent of fat-free mass on V")                                               # control stream ALLMV  = (FFM/41)**1
    e_ffm_qh        <- fixed(0.75) ; label("Allometric exponent of fat-free mass on QH")                                              # control stream QH = THETA(7)*(FFM/56)**0.75
    e_ffm_vh        <- fixed(1)    ; label("Allometric exponent of fat-free mass on VH")                                              # control stream VH = THETA(9)*(FFM/56)

    # ----- Country / formulation effects on prehepatic bioavailability ----
    # Reported in Table 2 as percentage changes; the control stream carries
    # them as multiplicative factors raised to the country indicator.
    e_region_canada_fdepot  <- 0.782 ; label("Multiplicative effect of Canadian enrolment (Rofact capsules) on prehepatic bioavailability")                    # Table 2 -21.8% (RSE 3.3%, 90% CI -27.9 to -18.0); 1 - 0.782 = 0.218; control stream $THETA (0, 0.782) ;[BIOFRM1]
    e_region_vietnam_fdepot <- 0.877 ; label("Multiplicative effect of Vietnamese enrolment (Svizera / Mekophar capsules) on prehepatic bioavailability")      # Table 2 -12.3% (RSE 3.6%, 90% CI -17.7 to -7.9); 1 - 0.877 = 0.123; control stream $THETA (0, 0.877) ;[BIOFRM2]

    # =======================================================================
    # Between-subject variability. Table 2 reports CV%, computed with the
    # footnote formula CV(%) = sqrt(exp(omega^2) - 1) * 100 (the radical is
    # dropped by the PDF's symbol font in the running text but printed in
    # Figure 2's caption area); the OMEGA values below are the control
    # stream's $OMEGA block and reproduce every printed CV to 3 figures.
    # BSV on Ka, MTT and NN was tested and fixed to zero in the final model.
    # =======================================================================
    etalclint_max ~ 0.0137   # control stream $OMEGA 0.0137 ;[IIV in CL];    sqrt(exp(0.0137) - 1) = 11.7% CV = Table 2 'BSV in CL 11.7'
    etalvc        ~ 0.0221   # control stream $OMEGA 0.0221 ;[IIV in V];     sqrt(exp(0.0221) - 1) = 14.9% CV = Table 2 'BSV in V 14.9'
    etalkm        ~ 0.0213   # control stream $OMEGA 0.0213 ;[IIV in LOGKM]; sqrt(exp(0.0213) - 1) = 14.7% CV = Table 2 'BSV in Km 14.7'

    # =======================================================================
    # Between-occasion variability, two occasions sharing one variance
    # ($OMEGA BLOCK(1) <value> followed by $OMEGA BLOCK(1) SAME). nlmixr2 has
    # no SAME keyword, so the second occasion's slot repeats the value with
    # fix(). Ka and MTT carry IOV only (their BSV was fixed to zero); F
    # carries IOV only because its typical value is a fixed anchor of 1.
    # =======================================================================
    etaiov_fdepot_1 ~ 0.0421       # control stream $OMEGA BLOCK(1) 0.0421 ;[IOV in BIO]; sqrt(exp(0.0421) - 1) = 20.7% CV = Table 2 'BOV in F 20.7'
    etaiov_fdepot_2 ~ fix(0.0421)  # control stream $OMEGA BLOCK(1) SAME
    etaiov_ka_1     ~ 1.18         # control stream $OMEGA BLOCK(1) 1.18 ;[IOV in KA];    sqrt(exp(1.18) - 1)   = 150.1% CV = Table 2 'BOV in Ka 150.1'
    etaiov_ka_2     ~ fix(1.18)    # control stream $OMEGA BLOCK(1) SAME
    etaiov_mtt_1    ~ 0.176        # control stream $OMEGA BLOCK(1) 0.176 ;[IOV in MTT];  sqrt(exp(0.176) - 1)  = 43.9% CV = Table 2 'BOV in MTT 43.9'
    etaiov_mtt_2    ~ fix(0.176)   # control stream $OMEGA BLOCK(1) SAME

    # =======================================================================
    # Residual error. The control stream's $ERROR block is
    #   PROP = IPRED*THETA(13); IF(ITS.EQ.1) PROP = IPRED*THETA(14)
    #   ADD  = LLOQ/5 + THETA(15)          with LLOQ = 0.125, THETA(15) = 0 FIX
    #   SD   = SQRT(ADD**2 + PROP**2);  Y = IPRED + SD*ERR(1);  $SIGMA 1 FIX
    # so THETA(13) and THETA(14) are proportional standard deviations on the
    # linear scale, matching nlmixr2's root-sum-square prop() + add().
    #
    # ERRATUM (see the vignette's Assumptions and deviations): Table 2 prints
    # these two rows as "35.1" and "45.3" CV%, which are the log-normal
    # back-transformations sqrt(exp(0.116) - 1) = 0.3507 and
    # sqrt(exp(0.187) - 1) = 0.4535 of the control stream's THETA(13) = 0.116
    # and THETA(14) = 0.187. That transformation belongs to the OMEGA rows
    # (Table 2 footnote a) and does not apply to a linear proportional-error
    # coefficient, so the table's residual-error CV% column overstates the
    # residual error roughly threefold. Both reproduce to three significant
    # figures, so the table was generated by applying the OMEGA formula to
    # these thetas. The control stream values are used here.
    # =======================================================================
    propSdSparse    <- 0.116        ; label("Proportional residual error, sparse PK sampling (fraction of the individual prediction)")     # control stream $THETA (0, 0.116) ;[PROP-ITS0]; printed in Table 2 as 35.1 CV% after an inapplicable log-normal back-transformation
    propSdIntensive <- 0.187        ; label("Proportional residual error, intensive PK sampling (fraction of the individual prediction)")  # control stream $THETA (0, 0.187) ;[PROP-ITS1]; printed in Table 2 as 45.3 CV% after an inapplicable log-normal back-transformation
    addSd           <- fixed(0.025) ; label("Additive residual error (mg/L)")                                                              # Table 2 'Additive, mg/L: 0.025 fixed' with footnote e; control stream ADD = LLOQ/5 + THETA(15) = 0.125/5 + 0
  })

  model({
    # --- 1. Occasion indicators for the between-occasion variability. -----
    #     The control stream's IF (OCC.EQ.1) THEN ... ELSE ... construct maps
    #     every occasion other than the first onto the second (SAME) eta.
    oc1 <- (OCC == 1)
    oc2 <- (OCC != 1)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2
    iov_ka     <- oc1 * etaiov_ka_1     + oc2 * etaiov_ka_2
    iov_mtt    <- oc1 * etaiov_mtt_1    + oc2 * etaiov_mtt_2

    # --- 2. Allometric scaling by fat-free mass. Two reference values are --
    #     in play: 41 kg (the cohort median) for the estimated disposition
    #     parameters and 56 kg (the fat-free mass of a 70 kg reference adult)
    #     for the fixed hepatic physiology, per Table 2 footnotes b and c.
    allm_cl <- (FFM / 41)^e_ffm_clint_max
    allm_v  <- (FFM / 41)^e_ffm_vc
    allm_qh <- (FFM / 56)^e_ffm_qh
    allm_vh <- (FFM / 56)^e_ffm_vh

    # --- 3. Country / formulation effect on prehepatic bioavailability. ----
    #     control stream: BIOFRM = (THETA(11)**FRM1) * (THETA(12)**FRM2)
    bio_form <- e_region_canada_fdepot^REGION_CANADA *
                e_region_vietnam_fdepot^REGION_VIETNAM

    # --- 4. Individual parameters. ----------------------------------------
    clint_max <- exp(lclint_max + etalclint_max) * allm_cl
    vc        <- exp(lvc + etalvc) * allm_v
    ka        <- exp(lka + iov_ka)
    mtt       <- exp(lmtt + iov_mtt)
    nn        <- exp(lnn)
    fdepot    <- exp(lfdepot + iov_fdepot) * bio_form
    qh        <- exp(lqh) * allm_qh
    vh        <- exp(lvh) * allm_vh

    # The BSV on Km is multiplicative on the LOG-scale value, not additive on
    # it: the control stream codes LOGKM = THETA(10)*EXP(IIVKM) and then uses
    # EXP(LOGKM) wherever Km is needed. That is why the eta is written
    # multiplying lkm here rather than being added to it.
    km   <- exp(lkm * exp(etalkm))
    # VMAX = CLINT*EXP(LOGKM), i.e. the product of the INDIVIDUAL CLint,max
    # (which already carries etalclint_max and the allometric term) and the
    # individual Km. It is the numerator of the saturable intrinsic clearance
    # and has units of mg/h.
    vmax <- clint_max * km

    # --- 5. Well-stirred liver with saturable intrinsic clearance. --------
    #     Figure 1 caption: CLint = (CLint,max * Km) / (CH + Km),
    #     EH = CLint * fu / (CLint * fu + QH), CLH = QH * EH.
    #     At CH = 0 the intrinsic clearance equals CLint,max and the
    #     extraction ratio is at its maximum; as the liver concentration
    #     rises, CLint falls, EH falls and exposure grows faster than
    #     proportionally with dose. The control stream guards the division
    #     with IF (CH>0), which is a no-op because every term it feeds is
    #     multiplied by the liver amount; the smooth form is used here so the
    #     solver sees no discontinuity.
    c_liver <- liver / vh
    clint   <- vmax / (c_liver + km)
    eh      <- (clint * fub) / (clint * fub + qh)
    fh      <- 1 - eh

    k20 <- qh * eh / vh   # liver -> eliminated
    k23 <- qh * fh / vh   # liver -> central (the fraction escaping extraction)
    k32 <- qh / vc        # central -> liver (recirculation with hepatic blood flow)

    # --- 6. ODE system, matching $DES DADT(1)-DADT(3) exactly. ------------
    #     transit() is rxode2's implementation of the Savic 2007 analytical
    #     transit-chain input rate and reproduces the control stream's
    #     PIZZA / TRANSIT construction term for term:
    #       TRANSIT = BIO*PD*KTR * (KTR*TAD)^NN * exp(-KTR*TAD) / gamma(NN+1)
    #     with KTR = (NN+1)/MTT. Prehepatic bioavailability enters as the
    #     bio argument, exactly as BIO enters PIZZA. The absorbed drug is
    #     delivered to the liver, so first-pass extraction is structural
    #     rather than a separate F term.
    d/dt(depot)   <- transit(nn, mtt, fdepot) - ka * depot
    d/dt(liver)   <- ka * depot - k23 * liver + k32 * central - k20 * liver
    d/dt(central) <- k23 * liver - k32 * central

    # F1 = 0 in the control stream: the dose event must not deposit a bolus
    # into the absorption compartment because transit() already feeds the
    # whole dose in through the analytical chain. transit() reads the raw
    # amount from podo(depot) regardless of f(depot).
    f(depot) <- 0

    # --- 7. Observation. ---------------------------------------------------
    Cc <- central / vc

    # Residual-error magnitude switched per record by the sampling strategy:
    #   SAMPLE_INTENSIVE = 1 -> intensive substudy, 18.7% proportional
    #   SAMPLE_INTENSIVE = 0 -> sparse substudy,    11.6% proportional
    propSd <- propSdIntensive * SAMPLE_INTENSIVE +
              propSdSparse * (1 - SAMPLE_INTENSIVE)
    Cc ~ prop(propSd) + add(addSd)
  })
}
