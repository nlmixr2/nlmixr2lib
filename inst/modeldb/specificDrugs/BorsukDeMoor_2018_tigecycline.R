BorsukDeMoor_2018_tigecycline <- function() {
  description <- "Two-compartment linear IV population PK model for high-dose tigecycline in adult ICU patients with sepsis or septic shock. Carries interindividual variability on clearance and on both volumes, plus interoccasion variability on clearance and peripheral volume across up to eight 12-hourly dosing occasions, and a combined additive-plus-proportional residual error. No covariate was retained: a visual search over age, weight, height, sex, ECMO, CRRT, dialysis volume, ultrafiltration rate, extravascular lung water index, cardiac output, SOFA score and procalcitonin found no systematic relationship with the individual PK parameters. Fitted by NONMEM 7.3 FOCE-I with eta-epsilon interaction (ADVAN3 TRANS4) to 940 plasma concentrations from 37 patients receiving a 200 mg loading dose followed by 100 mg every 12 h as 30-minute infusions."
  reference <- paste(
    "Borsuk-De Moor A, Rypulak E, Potrec B, Piwowarczyk P, Borys M, Sysiak J,",
    "Onichimowski D, Raszewski G, Czuczwar M, Wiczling P. Population",
    "pharmacokinetics of high-dose tigecycline in patients with sepsis or",
    "septic shock. Antimicrob Agents Chemother. 2018;62(4):e02273-17.",
    "doi:10.1128/AAC.02273-17.",
    "Structural model, parameter estimates and the random-effects structure are",
    "taken from Borsuk-De Moor 2018 Table 2 and the 'Pharmacokinetic modeling'",
    "subsection of Materials and Methods. The article's supplemental material",
    "(AAC.02273-17_zac004187006s1.pdf) contains only the chromatographic method,",
    "its validation, and Figures S1-S6 (time-dependent covariate summaries,",
    "goodness-of-fit, individual fits and eta-versus-covariate scatter plots);",
    "it carries no NONMEM control stream and no additional parameter values.",
    sep = " "
  )
  vignette <- "BorsukDeMoor_2018_tigecycline"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "tigecycline", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tigecycline", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    OCC = list(
      description = "Integer-valued dosing-occasion indicator for interoccasion variability on clearance and peripheral volume",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "An occasion is one 12-hourly tigecycline administration. Patients received 2 to 8 doses",
        "over 1 to 4 consecutive days, with samples drawn 0.5, 2, 4, 8 and 12 h after each dose",
        "(Borsuk-De Moor 2018, Patients and study design), so the model carries eight occasion slots.",
        "Borsuk-De Moor 2018 Table 2 confirms the count directly: the interoccasion rows for both",
        "pi^2 CL and pi^2 V2 print a separate shrinkage value for Occasion 1 through Occasion 8.",
        "All eight slots of a parameter share a single estimated magnitude -- the paper states the",
        "IOV variances 'were assumed to be constant across occasions' -- so occasions 2-8 are encoded",
        "as fixed() at the occasion-1 value, the nlmixr2 spelling of NONMEM $OMEGA BLOCK(1) SAME.",
        "Set OCC = 0 on every record to switch interoccasion variability off entirely, or OCC = 1",
        "to simulate a single typical occasion."
      ),
      source_name = "occasion (k in the P[i,k] = theta_P exp(eta_P,i) exp(kappa_P,i,k) parameterisation)"
    )
  )

  # Covariates that Borsuk-De Moor 2018 screened but did NOT retain. The paper
  # performed a visual covariate search -- plotting eta estimates against
  # time-independent covariates and kappa estimates against time-dependent
  # covariates (Figs. 3, 4 and S4-S6) -- found no regular trend, and therefore
  # never proceeded to formal statistical testing. No effect size is published
  # for any of them, so none can be encoded; they are recorded here so the
  # provenance of the covariate screen survives.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      notes = "Time-independent. Median 61 years, range 25-79 (Table 1); the protocol enrolled 18-75 year olds. Screened against eta_CL, eta_V1, eta_V2 and individual CL / Vss (Figs. 3 and S4-S6); no trend."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Time-independent. Median 80 kg, range 50-129 (Table 1). Screened directly and, per the Results, additionally through the derived body surface area and body mass index; no trend was discovered for any of the three. The paper contrasts this with Van Wart 2006 and Rubino 2010, both of which retained weight-related terms on CL."
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      notes = "Time-independent. Median 175 cm, range 158-190 (Table 1). Screened alone and as an input to body surface area and body mass index; no trend."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Time-independent. 26 male / 11 female, i.e. 29.7% female (Table 1). Screened as a categorical covariate in Fig. 3; no trend. Van Wart 2006 retained gender on CL, Borsuk-De Moor 2018 did not."
    ),
    ECMO_STATUS = list(
      description = "Extracorporeal membrane oxygenation treatment-status indicator (1 = on ECMO, 0 = not)",
      units = "(binary)",
      type = "binary",
      notes = "Time-independent as screened. 35 no / 2 yes (Table 1). Only two patients received ECMO, so the categorical panel of Fig. 3 cannot support an effect estimate even had one been apparent."
    ),
    RRT_CRRT_STATUS = list(
      description = "Continuous renal replacement therapy status (no / yes / started during therapy)",
      units = "(category)",
      type = "categorical",
      notes = "Three levels in Table 1: 6 no / 30 yes / 1 started during therapy. Screened in the Fig. 3 categorical panel; no trend. The Discussion notes explicitly that the non-CRRT group was too small to assess the impact of CRRT on individual PK parameters, against the suggestion of Honore 2013 that CRRT influences tigecycline PK."
    ),
    DIALYSIS_VOLUME = list(
      description = "Dialysis volume delivered over the occasion, normalised to body weight",
      units = "mL/kg",
      type = "continuous",
      notes = "Time-dependent; median 23.8 mL/kg, range 14.2-40.0 (Table 1), with its day-by-day distribution in Fig. S1. Screened against the kappa estimates for CL and V2 in Fig. 4; no trend. No canonical register entry is minted because the covariate is documentation-only here -- it is never referenced in model()."
    ),
    ULTRAFILTRATION_RATE = list(
      description = "Ultrafiltration rate of the renal-replacement circuit, normalised to body weight",
      units = "mL/kg/h",
      type = "continuous",
      notes = "Time-dependent; median 1.54 mL/kg/h, range 0.34-6.6 (Table 1), day-by-day distribution in Fig. S1. Screened against the kappa estimates for CL and V2 in Fig. 4; no trend. Documentation-only, never referenced in model()."
    ),
    ELWI = list(
      description = "Extravascular lung water index",
      units = "mL/kg",
      type = "continuous",
      notes = "Time-dependent; median 9 mL/kg, range 5-41 (Table 1), day-by-day distribution in Fig. S1. Screened against the kappa estimates for CL and V2 in Fig. 4 as a marker of the capillary-leak fluid shift the Discussion invokes to explain the raised V1; no trend. Documentation-only, never referenced in model()."
    ),
    CARDIAC_OUTPUT = list(
      description = "Cardiac output",
      units = "L/min",
      type = "continuous",
      notes = "Time-dependent; median 7.49, range 2.55-15.8 (Table 1, where the unit is printed as 'liters' -- cardiac output is conventionally a flow, L/min, and the quoted range matches the hyperdynamic septic values expected on that scale), day-by-day distribution in Fig. S1. Screened against the kappa estimates for CL and V2 in Fig. 4; no trend. Documentation-only, never referenced in model()."
    ),
    SOFA = list(
      description = "Sequential Organ Failure Assessment score",
      units = "points",
      type = "continuous",
      notes = "Time-dependent; median 13 points, range 2.0-21 (Table 1), day-by-day distribution in Fig. S1. Screened against the kappa estimates for CL and V2 in Fig. 4 as the summary organ-failure severity measure; no trend. Documentation-only, never referenced in model()."
    ),
    PROCALCITONIN = list(
      description = "Serum procalcitonin concentration",
      units = "ug/L",
      type = "continuous",
      notes = "Time-dependent; median 8.22, range 0.16-122 (Table 1). Table 1 prints the unit as 'umol/liter', which is not a plausible procalcitonin scale -- procalcitonin is a 13 kDa peptide reported clinically in ng/mL (= ug/L), and 8.22 ng/mL is a typical septic value whereas 8.22 umol/L would be roughly 10^8-fold above any reported concentration. The canonical ug/L reading is used here. Screened against the kappa estimates for CL and V2 in Fig. 4; no trend. Documentation-only, never referenced in model()."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 37L,
    n_studies = 1L,
    age_range = "25-79 years",
    age_median = "61 years",
    weight_range = "50-129 kg",
    weight_median = "80 kg",
    height_range = "158-190 cm",
    height_median = "175 cm",
    sex_female_pct = 29.7,
    race_ethnicity = NULL,
    disease_state = "Adults admitted to a tertiary medical/surgical ICU with sepsis or septic shock of medical or surgical origin, with suspected nosocomial infection by multidrug-resistant or extensively drug-resistant strains requiring empirical broad-spectrum antibiotics. Median SOFA score 13 points (range 2.0-21), median albumin 2.2 g/dL (range 1.5-3.6), 23 of 37 patients died. 30 of 37 were receiving continuous renal replacement therapy and 2 were on extracorporeal membrane oxygenation.",
    dose_range = "High-dose regimen throughout: a 200 mg loading dose as a 30-minute intravenous infusion, followed by 100 mg every 12 h as 30-minute infusions. Patients received 2 to 8 doses over 1 to 4 consecutive days. This is twice the licensed 100 mg / 50 mg q12h regimen.",
    regions = "Poland (two tertiary ICUs, Lublin and Olsztyn)",
    notes = "Prospective observational cohort study. 942 tigecycline concentrations were collected; 2 were removed as outliers during model building (CWRES > 5), leaving 940 observations from 37 patients. Arterial blood was sampled 0.5, 2, 4, 8 and 12 h after every dose. Exclusion criteria: no life-threatening condition at symptom onset, HIV infection, terminal cancer, prior tigecycline intolerance or allergy, high probability of infection by a tigecycline-resistant strain such as Pseudomonas aeruginosa, and tigecycline exposure within the preceding 3 months. Plasma was assayed by HPLC-UV at 225 nm, linear over 0.078-2.5 ug/mL with LOQ below 0.078 ug/mL and LOD below 0.02 ug/mL; intra- and interday accuracy 98.4-103.1% and CV 0.6-3.8% (Borsuk-De Moor 2018 Assay and supplemental material). Estimation used NONMEM 7.3 with FOCE-I (eta-epsilon interaction) and the ADVAN3 TRANS4 routine, driven by Wings for NONMEM WFN741; uncertainty came from a 1000-replicate nonparametric bootstrap with 1 unsuccessful run. The paper's headline derived quantity is Vss = V1 + V2 = 250 L, and it reports an accumulation ratio of 1.57 on this regimen."
  )

  ini({
    # Structural parameters. Two-compartment disposition parameterised in
    # clearances and volumes (NONMEM ADVAN3 TRANS4), fitted to intravenous
    # infusion data only, so CL and V1 are absolute rather than apparent.
    # Borsuk-De Moor 2018 Table 2, 'Mean population parameters (theta)'.
    lcl <- log(22.1); label("Elimination clearance CL (L/h)")                     # Borsuk-De Moor 2018 Table 2: theta CL = 22.1 L/h (RSE 3.16%; bootstrap median 22.1, 90% CI 20.9-23.2)
    lvc <- log(162);  label("Volume of distribution of the central compartment V1 (L)")    # Borsuk-De Moor 2018 Table 2: theta V1 = 162 L (RSE 5.3%; bootstrap median 163, 90% CI 150-176)
    lq  <- log(69.4); label("Intercompartmental clearance Q (L/h)")               # Borsuk-De Moor 2018 Table 2: theta Q = 69.4 L/h (RSE 32.6%; bootstrap median 67.3, 90% CI 41.9-98.4). The Discussion attributes the low precision to the sampling design carrying little information about Q.
    lvp <- log(87.9); label("Volume of distribution of the peripheral compartment V2 (L)") # Borsuk-De Moor 2018 Table 2: theta V2 = 87.9 L (RSE 8.67%; bootstrap median 87.6, 90% CI 76.1-101). V1 + V2 = 250 L = the Vss quoted in the Discussion.

    # Interindividual variability. Table 2 labels these rows 'omega^2 CL (% CV)'
    # etc., i.e. the estimated NONMEM variance displayed on the %CV scale, so
    # each variance below is the printed percentage squared. No omega^2 for Q was
    # estimable -- Table 2 prints '0 FIX' and the Results state it 'was not
    # possible to estimate the interindividual variability (IIV) for the
    # intercompartmental clearance (Q2)' -- so q carries no eta in model().
    # Table 2 gives no CV formula in a footnote. The conventional NONMEM reading
    # CV% = 100*sqrt(omega^2) is used; the exact log-normal alternative
    # CV% = 100*sqrt(exp(omega^2)-1) would give 0.02949 / 0.03625 / 0.13966,
    # a difference of at most 3.6% on the omega scale. See the vignette
    # 'Assumptions and deviations' section.
    etalcl ~ 0.029929 # Borsuk-De Moor 2018 Table 2, row 'omega^2 CL (% CV)' = 17.3 (RSE 19%, shrinkage 7.3%; bootstrap median 17.1, 90% CI 14.2-19.7); 0.173^2 = 0.029929
    etalvc ~ 0.036864 # Borsuk-De Moor 2018 Table 2, row 'omega^2 V1 (% CV)' = 19.2 (RSE 29.2%, shrinkage 6.7%; bootstrap median 19.1, 90% CI 14.3-23.7); 0.192^2 = 0.036864
    etalvp ~ 0.149769 # Borsuk-De Moor 2018 Table 2, row 'omega^2 V2 (% CV)' = 38.7 (RSE 40.8%, shrinkage 22.4%; bootstrap median 37.4, 90% CI 20.4-48.5); 0.387^2 = 0.149769

    # Interoccasion variability on CL and V2 over eight 12-hourly dosing
    # occasions, one shared magnitude per parameter. Table 2 prints a single
    # pi^2 estimate per parameter and then a separate shrinkage value for each
    # of Occasions 1-8, and Materials and Methods states the IOV variances
    # 'were assumed to be constant across occasions', so occasions 2-8 are
    # fixed at the occasion-1 value (NONMEM $OMEGA BLOCK(1) SAME).
    etaiov_lcl_1 ~ 0.020736        # Borsuk-De Moor 2018 Table 2, row 'pi^2 CL (% CV)' = 14.4 (RSE 35%; bootstrap median 14.2, 90% CI 9.4-18.1); 0.144^2 = 0.020736. Occasion 1, shrinkage 27.4%.
    etaiov_lcl_2 ~ fixed(0.020736) # Occasion 2, shrinkage 38.5%; same magnitude per the constant-across-occasions assumption
    etaiov_lcl_3 ~ fixed(0.020736) # Occasion 3, shrinkage 52.4%
    etaiov_lcl_4 ~ fixed(0.020736) # Occasion 4, shrinkage 42.6%
    etaiov_lcl_5 ~ fixed(0.020736) # Occasion 5, shrinkage 35.4%
    etaiov_lcl_6 ~ fixed(0.020736) # Occasion 6, shrinkage 48.4%
    etaiov_lcl_7 ~ fixed(0.020736) # Occasion 7, shrinkage 94.9%; the Results note only a few observations were available on occasions 7 and 8
    etaiov_lcl_8 ~ fixed(0.020736) # Occasion 8, shrinkage 99.8%
    etaiov_lvp_1 ~ 0.043264        # Borsuk-De Moor 2018 Table 2, row 'pi^2 V2 (% CV)' = 20.8 (RSE 66.4%; bootstrap median 21.7, 90% CI 0.200-30.9); 0.208^2 = 0.043264. Occasion 1, shrinkage 40.2%.
    etaiov_lvp_2 ~ fixed(0.043264) # Occasion 2, shrinkage 50.2%
    etaiov_lvp_3 ~ fixed(0.043264) # Occasion 3, shrinkage 58.5%
    etaiov_lvp_4 ~ fixed(0.043264) # Occasion 4, shrinkage 56.4%
    etaiov_lvp_5 ~ fixed(0.043264) # Occasion 5, shrinkage 52.9%
    etaiov_lvp_6 ~ fixed(0.043264) # Occasion 6, shrinkage 57.2%
    etaiov_lvp_7 ~ fixed(0.043264) # Occasion 7, shrinkage 90.6%
    etaiov_lvp_8 ~ fixed(0.043264) # Occasion 8, shrinkage 99.6%

    # Residual error: combined additive plus proportional. The two rows are on
    # DIFFERENT scales, as Table 2's own footnote a makes explicit -- 'sigma_add,
    # additive residual random error; sigma^2_prop, variance of proportional
    # residual random error'. The additive row is therefore already a standard
    # deviation in ug/mL (= mg/L) and is used as printed, while the proportional
    # row follows the same 'variance displayed as % CV' convention as the omegas.
    addSd  <- 0.0210; label("Additive residual error (mg/L)")        # Borsuk-De Moor 2018 Table 2: sigma_add = 0.0210 ug/mL = 0.0210 mg/L (bootstrap median 0.0224, 90% CI 0.000209-0.0357). The printed RSE of 0.41% is inconsistent with that bootstrap interval, whose width implies an RSE near 100%; the point estimate is used regardless.
    propSd <- 0.130;  label("Proportional residual error (fraction)") # Borsuk-De Moor 2018 Table 2: sigma^2_prop = 13.0 % CV (RSE 17.7%; bootstrap median 12.7, 90% CI 8.89-16.1); 13.0% -> 0.130 on the standard-deviation scale nlmixr2 expects
  })

  model({
    # Interoccasion variability, multiplexed by the occasion indicator. OCC = 0
    # zeroes every indicator and switches interoccasion variability off.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)
    oc7 <- (OCC == 7)
    oc8 <- (OCC == 8)
    iov_cl <- oc1 * etaiov_lcl_1 + oc2 * etaiov_lcl_2 +
      oc3 * etaiov_lcl_3 + oc4 * etaiov_lcl_4 +
      oc5 * etaiov_lcl_5 + oc6 * etaiov_lcl_6 +
      oc7 * etaiov_lcl_7 + oc8 * etaiov_lcl_8
    iov_vp <- oc1 * etaiov_lvp_1 + oc2 * etaiov_lvp_2 +
      oc3 * etaiov_lvp_3 + oc4 * etaiov_lvp_4 +
      oc5 * etaiov_lvp_5 + oc6 * etaiov_lvp_6 +
      oc7 * etaiov_lvp_7 + oc8 * etaiov_lvp_8

    # 1. Individual PK parameters. Borsuk-De Moor 2018 Materials and Methods:
    #    P[i,k] = theta_P * exp(eta_P,i) * exp(kappa_P,i,k), i.e. log-normal IIV
    #    and log-normal IOV multiplying the typical value. No covariate was
    #    retained in the final model. Q carries neither eta nor kappa.
    cl <- exp(lcl + etalcl + iov_cl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp + iov_vp)

    # 2. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. ODE system. Tigecycline is administered as a 30-minute intravenous
    #    infusion into the central compartment; the infusion is supplied through
    #    the event table (rate = -2 with dur(central), or an explicit rate).
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 4. Observation and error. Dose in mg, volumes in L -> central / vc has
    #    units mg/L, which equals the ug/mL of the paper's assay and Table 2.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
