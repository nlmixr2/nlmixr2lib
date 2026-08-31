Storgaard_2026_tetrahydrocannabinol <- function() {
  description <- paste(
    "Parent-plus-metabolite population pharmacokinetic model for",
    "delta-9-tetrahydrocannabinol (THC) and its active metabolite",
    "11-hydroxy-delta-9-tetrahydrocannabinol (11-OH-THC, written THC-OH by",
    "the authors) after OROMUCOSAL administration of the THC/CBD spray",
    "Sativex in acutely hospitalised older medical patients with poor",
    "appetite. THC absorption is a chain of three sequential first-order",
    "transit steps sharing one rate constant ktr, and THC disposition is",
    "one-compartment; 11-OH-THC is two-compartment. Total apparent THC",
    "clearance is decomposed into an apparent FORMATION clearance into",
    "11-OH-THC (CLF/F = 765 L/h) and an apparent clearance through other",
    "pathways (CLO/F = 162 L/h), summing to the CL/F of 927 L/h the paper",
    "reports. Because THC was given only oromucosally, every clearance and",
    "volume is apparent (X/F). Variability is unusually large and is placed",
    "as inter-occasion variability on MTT and CLF (the two doses 4 h apart",
    "are the two occasions) and as inter-individual variability on V, CLO",
    "and CLMET; the authors attribute it to variability in oromucosal",
    "bioavailability. No covariate entered the final model: age, body",
    "composition and renal function were all screened and none was",
    "retained. The terminal elimination phase was not captured (sampling",
    "stopped 4 h after the last dose), so the apparent clearances are not",
    "directly comparable with models built on longer profiles."
  )
  reference <- paste(
    "Storgaard IK, Nielsen RL, Houlind MB, Bornaes O, Christensen LWS,",
    "Andersen AL, Juul-Larsen HG, Jorgensen LM, Breindahl T, Jawad BN,",
    "Altintas I, Andersen O, Lund TM. Population pharmacokinetic modelling",
    "revealed large variability in oromucosal absorption of",
    "delta-9-tetrahydrocannabinol in older patients with poor appetite.",
    "Br J Clin Pharmacol. 2026;92(2):504-514. doi:10.1002/bcp.70284.",
    "PMCID: PMC12850555."
  )
  vignette <- "Storgaard_2026_tetrahydrocannabinol"

  # Methods 2.1 and Table 2: doses are quoted in mg of THC per
  # administration (5.2 mg for two sprays, 8.1 mg for three), clearances in
  # L/h and volumes in L, while every plasma concentration in Figure 1,
  # Figure 5 and the Results text is in ug/L. amount(mg)/volume(L) is
  # therefore mg/L and needs the explicit 1000 ug/mg factor applied in
  # model() to land on the paper's reporting scale.
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    OCC = list(
      description        = "Dosing occasion index (1 = first dose at t = 0 h, 2 = second dose at t = 4 h)",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Methods 2.2: 'Variability could be implemented on parameters as",
        "interindividual variability (IIV, or between-subject variability)",
        "and/or interoccasion variability (IOV, or within-subject",
        "variability), where each dosing of Sativex was considered the",
        "beginning of a new occasion.' Methods 2.1 gives exactly two doses,",
        "4 h apart, so the model has exactly two occasions. OCC is decomposed",
        "in model() into the binary indicators oc1 / oc2 that select the",
        "per-occasion IOV slots on MTT and CLF. Records before the first dose",
        "should carry OCC = 1."
      ),
      source_name        = "OCC"
    )
  )

  # Screened but NOT retained in the final model. Results 3.2: 'None of the
  # potential covariates appeared to be correlated with any of the model
  # parameters based on visual screening of eta-covariate plots (see
  # Supporting Information Figures S5, S6 and S7). Allometric weight scaling
  # was implemented on clearance and volume parameters (as a power model with
  # an exponent of 0.75 and 1 for clearances and volumes, respectively), but
  # this did not improve the model significantly (dOFV = -3.31). Ultimately,
  # no covariates were included in the final model.' No point estimate is
  # published for any of these, so none can be encoded; they are recorded
  # here so the covariate screen is not lost.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened. Allometric scaling on clearances (exponent 0.75) and",
        "volumes (exponent 1) was explicitly tested and REJECTED: Results 3.2",
        "reports dOFV = -3.31, short of the 3.84 forward-addition threshold",
        "set in Methods 2.2. Table 1 cohort weight was 56 (48, 70) kg, median",
        "(IQR), range 38-82 kg."
      )
    ),
    HT = list(
      description = "Body height",
      units       = "cm",
      type        = "continuous",
      notes       = "Measured on each trial day (Methods 2.1); Table 1 median 165 (157, 172) cm, range 149-179. Screened via eta-covariate plots, not retained."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Table 1 median 21.8 (20.0, 25.0) kg/m^2, range 14.7-27.3. BMI <= 30 was an eligibility criterion (Methods 2.1). Screened, not retained."
    ),
    BSA = list(
      description = "Body surface area (Du Bois formula)",
      units       = "m^2",
      type        = "continuous",
      notes       = "Table 1 median 1.58 (1.45, 1.85) m^2, range 1.29-1.98; Table 1 footnote cites the Du Bois formula. Screened, not retained."
    ),
    BODYFAT_PCT = list(
      description = "Percent total body fat by bioelectrical impedance",
      units       = "%",
      type        = "continuous",
      notes       = paste(
        "Measured with a portable bedside direct segmental multifrequency",
        "bioelectrical impedance analyser (InBody S10; Biospace), Methods 2.1.",
        "Table 1 median 29.6 (25.6, 33.5) %, range 17.3-41.1. Of particular",
        "interest a priori because THC is highly lipophilic and the cohort was",
        "selected for low or declining body weight and fat percentage",
        "(Introduction), but the Discussion records the negative result: 'The",
        "large inter- and intra-individual variability could not be explained",
        "by any of the measured physiological/biological characteristics such",
        "as body composition indicators (eg, fat mass and fat percentage)'.",
        "The same impedance measurement also produced absolute body fat mass",
        "(16.4 [13.8, 22.2] kg) and muscle mass (20.2 [18.2, 27.0] kg), which",
        "were screened alongside it; neither has a separate register entry",
        "because neither is used by any model."
      )
    ),
    FFM = list(
      description = "Fat-free mass by bioelectrical impedance",
      units       = "kg",
      type        = "continuous",
      notes       = "Table 1 median 39 (35, 50) kg, range 27-55. Screened, not retained."
    ),
    CRCL = list(
      description = "BSA-normalised glomerular filtration rate (both the CKD-EPI 2009 creatinine estimate eGFR and the 99m-Tc-DTPA tracer measurement mGFR)",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = paste(
        "BOTH renal measures were screened and neither was retained. Table 1:",
        "eGFR 62 (44, 83), range 39-98 (2009 CKD-EPI without race factor,",
        "Methods 2.1); mGFR 62 (47, 83), range 34-121, available in only 14 of",
        "20 patients (Table 1 footnote '*n = 14'). Discussion: 'While we",
        "screened eGFR and mGFR as potential covariates, we did not see any",
        "correlation with pharmacokinetic parameters. As renal excretion of",
        "unchanged THC is low, these results were not unexpected.'"
      )
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Table 1 median 77 (71, 85) years, range 66-94; age >= 65 years was an eligibility criterion. Screened, not retained."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Table 1: 14 of 20 female (70%), 6 male (30%). Screened, not retained."
    )
  )

  compartmentData <- list(
    # Absorption chain. Figure 3 draws THREE boxes each labelled 'Transit'
    # in series between Dose and Central THC, joined by three ktr arrows,
    # and annotates 'MTT = 3/ktr'. nlmixr2lib names the DOSED box of such a
    # chain `depot` and the remainder `transit<n>` (the Savic 2007 layout
    # used by Wattanakul_2024_primaquine_motherinfant.R and
    # Svensson_2012_nevirapine.R), so the paper's three boxes map to
    # depot -> transit1 -> transit2 -> central, which is the same three
    # first-order ktr steps.
    depot = list(
      analyte = "delta-9-tetrahydrocannabinol", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    transit1 = list(
      analyte = "delta-9-tetrahydrocannabinol", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    transit2 = list(
      analyte = "delta-9-tetrahydrocannabinol", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "delta-9-tetrahydrocannabinol", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    # Methods 2.2: 'THC-OH concentrations were normalized using the molecular
    # weight ratio with THC.' The metabolite states are therefore carried in
    # mg of THC EQUIVALENTS and the formation flux transfers mass 1:1 with no
    # conversion factor; Cc_11oh is likewise a THC-equivalent concentration,
    # which is the scale Table 2, Figure 1 and Figure 5 all report.
    central_11oh = list(
      analyte = "11-hydroxy-delta-9-tetrahydrocannabinol", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1_11oh = list(
      analyte = "11-hydroxy-delta-9-tetrahydrocannabinol", units = "mg",
      specimen = "plasma", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 20L,
    n_studies      = 1L,
    n_observations = paste(
      "370 observations from the Sativex trial days entered model building:",
      "185 for THC (one excluded as an outlier, leaving 184; 38 [20.5%] below",
      "the LLOQ) and 185 for THC-OH (no outliers; 34 [18.4%] below the LLOQ).",
      "19.5% of all observations were below the LLOQ and were handled with",
      "Beal's M3 method (Methods 2.2). Blood was drawn at nine to 14 time",
      "points under one of three sampling schemes (Supporting Information",
      "Table S1); 16 of 20 patients (80%) followed scheme C and two followed",
      "each of the others."
    ),
    age_median     = "77 years (IQR 71-85)",
    age_range      = "66-94 years",
    weight_median  = "56 kg (IQR 48-70)",
    weight_range   = "38-82 kg",
    height_median  = "165 cm (IQR 157-172)",
    bmi_median     = "21.8 kg/m^2 (IQR 20.0-25.0), range 14.7-27.3",
    sex_female_pct = 70,
    body_fat       = "29.6% (IQR 25.6-33.5), range 17.3-41.1; fat mass 16.4 kg (IQR 13.8-22.2); fat-free mass 39 kg (IQR 35-50)",
    renal_function = paste(
      "eGFR 62 mL/min/1.73 m^2 (IQR 44-83), range 39-98, by the 2009 CKD-EPI",
      "creatinine equation without race factor; mGFR 62 (IQR 47-83), range",
      "34-121, by 99m-Tc-DTPA in the 14 patients who joined that sub-study."
    ),
    disease_state  = paste(
      "Acutely hospitalised older medical patients with poor appetite,",
      "enrolled at an emergency department and studied after discharge.",
      "Eligibility (Methods 2.1) required age >= 65 years, BMI <= 30 and a",
      "Simplified Nutritional Appetite Questionnaire (SNAQ) score <= 14;",
      "cohort SNAQ was 12 (IQR 10-13), range 5-14. Thirteen of 20 (65%) had",
      "mild to severe xerostomia. Regular cannabis use before or within two",
      "weeks prior to baseline was not permitted and a positive morning-urine",
      "cannabinoid immunoassay (cut-off 50 ng/mL) was exclusionary."
    ),
    dose_range     = paste(
      "Two oromucosal doses 4 h apart on a single trial day. Scheme A gave two",
      "sprays per dose (5.2 mg THC and 5 mg CBD); schemes B and C gave three",
      "sprays per dose (8.1 mg THC and 7.5 mg CBD). 16 of 20 patients (80%)",
      "followed scheme C."
    ),
    co_medication  = paste(
      "Cannabidiol is co-administered in every Sativex spray (5 or 7.5 mg per",
      "dose) but is not modelled here. The Discussion argues the CYP",
      "inhibition is negligible at these doses: 'the maximum dose used in this",
      "study of 7.5 mg of CBD (with 8.1 mg of THC) is in this context",
      "considered low.'"
    ),
    regions        = "Denmark (single centre: Copenhagen University Hospital, patients enrolled at Hvidovre and dosed at the Zelo Phase 1 Unit, Bispebjerg)",
    notes          = paste(
      "Protocolised sub-study of a triple-blinded, randomised,",
      "placebo-controlled cross-over controlled-feeding trial run February",
      "2022 to January 2024; only the Sativex arm is used here (no cannabinoid",
      "was detected on placebo days). Patients fasted 12 h beforehand, then",
      "received a standardised 125 mL nutritional drink 25 min after the first",
      "dose and a standardised lunch 25 min after the second. They were asked",
      "not to talk, swallow or move their tongue for 15 min after each spray.",
      "Assay was UHPLC-MS/MS with an LLOQ of 0.25 ng/mL for every analyte and a",
      "CV of at most 15%. Demographics are Table 1. EudraCT 2021-002318-15,",
      "ClinicalTrials.gov NCT05503147.",
      "MODEL RELIABILITY CAVEAT, stated by the authors: only 338 of 1000",
      "bootstrap runs succeeded (33.8%) and 'in general, stability was an",
      "issue during model development ... When applying this model to",
      "different data sets, this low reliability must be considered.' CLO/F is",
      "singled out as poorly estimated (bootstrap mean 39.6 vs the estimate of",
      "162 L/h, 95% CI 1.24-647)."
    )
  )

  ini({
    # ==================================================================
    # Structural parameters -- Storgaard 2026 Table 2, 'Parameter
    # estimates for the final population pharmacokinetic
    # parent-metabolite model'. The 'Estimate' column is used throughout;
    # the bootstrap mean and 95% CI are recorded per line but are NOT the
    # model values.
    #
    # Every clearance and volume is APPARENT (X/F). Methods 2.2: 'As THC
    # was administered exclusively via the oromucosal route, the
    # estimated clearances and volumes are apparent, ie, the volume of
    # the central THC compartment will be the apparent volume of
    # distribution, defined as V/F, and the same will be the case for
    # other clearance and volume parameters.' No bioavailability
    # parameter is estimated, so none is encoded and f(depot) is left at
    # its default of 1.
    # ==================================================================

    lmtt <- log(1.35)
    label("Mean transit time of the oromucosal absorption chain, MTT (h)")            # Table 2: MTT 1.35 (bootstrap mean 1.39, 95% CI 1.05-1.87). Discussion repeats it: 'Absorption through the oromucosal route after administration of Sativex had an estimated mean transit time of 1.35 h.'

    lcl_met <- log(765)
    label("Apparent formation clearance of THC to 11-OH-THC, CLF/F (L/h)")            # Table 2: CLF/F 765 (bootstrap mean 696, 95% CI 378-1158)

    lcl_nonmet <- log(162)
    label("Apparent clearance of THC through other pathways, CLO/F (L/h)")            # Table 2: CLO/F 162 (bootstrap mean 39.6, 95% CI 1.24-647). Poorly estimated -- see the Discussion caveat recorded in population$notes. Results 3.2 confirms the sum: 'Based on the estimated CLF/F and CLO/F, the total apparent clearance (CL/F) of THC could be calculated as the sum of the two and had a value of 927 L h-1' (765 + 162 = 927).

    lvc <- log(1669)
    label("Apparent central volume of distribution for THC, V/F (L)")                 # Table 2: V/F 1669 (bootstrap mean 1532, 95% CI 925-2588)

    lcl_11oh <- log(281)
    label("Apparent clearance of 11-OH-THC, CLMET/F (L/h)")                           # Table 2: CLMET/F 281 (bootstrap mean 204, 95% CI 3.52-516). Results 3.2: CLMET 'describes the apparent clearance of THC-OH (mainly to THC-COOH)'.

    lvc_11oh <- log(77.5)
    label("Apparent central volume of distribution for 11-OH-THC, VMET/F (L)")        # Table 2: VMET/F 77.5 (bootstrap mean 66.4, 95% CI 36.9-110)

    lq_11oh <- log(183)
    label("Apparent intercompartmental clearance for 11-OH-THC, QMET/F (L/h)")        # Table 2: QMET/F 183 (bootstrap mean 185, 95% CI 82.8-354)

    lvp_11oh <- log(539)
    label("Apparent peripheral volume of distribution for 11-OH-THC, V2MET/F (L)")    # Table 2: V2MET/F 539 (bootstrap mean 679, 95% CI 287-2727)

    # ------------------------------------------------------------------
    # Inter-individual variability. Table 2 reports IIV and IOV as CV%,
    # and Methods 2.2 equation (3) gives the conversion the authors used:
    # CV% = sqrt(exp(omega^2) - 1) * 100%. Inverting it, omega^2 =
    # log(1 + CV^2). Each line below records the published CV% and the
    # variance it implies; every one round-trips back to the printed CV%
    # to the tabulated precision.
    #
    # Results 3.2 fixes which parameter carries which: 'The final
    # parent-metabolite model had IOV on MTT and CLF/F, and IIV on V/F,
    # CLO/F and CLMET/F.' No parameter carries both, and VMET, QMET and
    # V2MET carry neither. No correlation between random effects is
    # reported, so all slots are diagonal.
    # ------------------------------------------------------------------
    etalvc         ~ 0.169658   # Table 2 IIVV     43.0% CV (bootstrap 42.1, 95% CI 11.0-60.0; shrinkage 32%); omega^2 = log(1 + 0.430^2)
    etalcl_11oh    ~ 0.329753   # Table 2 IIVCLMET 62.5% CV (bootstrap 101,  95% CI 50.8-1383; shrinkage 10%); omega^2 = log(1 + 0.625^2)
    etalcl_nonmet  ~ 1.197069   # Table 2 IIVCLO   152%  CV (bootstrap 206,  95% CI 80.9-1682; shrinkage 12%); omega^2 = log(1 + 1.52^2). The Discussion quotes the same quantity to one more digit as 152.3% ('With IOV and IIV parameters estimated at 40.2-152.3% (CV%)'); Table 2 is used, and the difference is immaterial (omega^2 1.19707 vs 1.19919).

    # ------------------------------------------------------------------
    # Inter-occasion variability. Methods 2.2: 'each dosing of Sativex
    # was considered the beginning of a new occasion', and Methods 2.1
    # gives exactly two doses 4 h apart, so there are exactly two
    # occasions. The paper reports ONE IOV magnitude per parameter, i.e.
    # the occasions share a variance -- the NONMEM '$OMEGA BLOCK(1) SAME'
    # pattern -- so occasion 2 is fixed() to occasion 1's value, as in
    # Wattanakul_2024_primaquine_motherinfant.R and Goggin_2004_emfilermin.R.
    # ------------------------------------------------------------------
    etaiov_mtt_1     ~ 0.299572        # Table 2 IOVMTT 59.1% CV (bootstrap 58.4, 95% CI 42.1-73.3; shrinkage 4%);  omega^2 = log(1 + 0.591^2)
    etaiov_mtt_2     ~ fixed(0.299572) # same variance on occasion 2 ($OMEGA BLOCK(1) SAME); only one IOVMTT row is published
    etaiov_cl_met_1  ~ 0.149802        # Table 2 IOVCLF 40.2% CV (bootstrap 40.5, 95% CI 18.4-58.6; shrinkage 13%); omega^2 = log(1 + 0.402^2)
    etaiov_cl_met_2  ~ fixed(0.149802) # same variance on occasion 2 ($OMEGA BLOCK(1) SAME); only one IOVCLF row is published

    # ------------------------------------------------------------------
    # Residual unexplained variability. Methods 2.2: 'For model
    # development, plasma concentrations were log-transformed and an
    # exponential error model was used to handle residual unexplained
    # variability (RUV).' An additive error on log-transformed
    # concentrations is a log-normal (exponential) error on the linear
    # scale, so the observation model is lnorm(SD-on-log-scale).
    #
    # SCALE OF THE TABLE 2 VALUES. The two residual rows are printed as
    # bare numbers (0.447, 0.315) while every IIV/IOV row is printed as a
    # CV%, which is the signature of untransformed NONMEM $SIGMA output:
    # NONMEM reports $OMEGA and $SIGMA as VARIANCES, and the authors put
    # only $OMEGA through their equation (3) CV% conversion. Applying
    # equation (3) to 0.447 would have printed 75.1, not 0.447.
    #
    # Independently confirmed by the bootstrap CI widths, which have
    # different sampling theory for a variance than for an SD. Relative
    # 95%-CI width for THC is (0.509 - 0.333)/0.431 = 0.408. A variance
    # estimated from N observations predicts 3.92*sqrt(2/N) = 0.409 at
    # the N = 184 THC observations actually modelled; an SD predicts
    # 3.92*sqrt(1/(2N)) = 0.204, off by a factor of two. THC-OH agrees
    # (observed 0.438 against 0.408 predicted at N = 185).
    #
    # So these are variances and the log-scale SD is their square root,
    # following the idiom of Wu_2024_gfr_maturation.R.
    # ------------------------------------------------------------------
    expSd <- sqrt(0.447)
    label("Log-normal (exponential) residual error for THC, SD on the log scale")     # Table 2: ERRTHC 0.447 (bootstrap mean 0.431, 95% CI 0.333-0.509); sqrt(0.447) = 0.6686, i.e. 75.1% CV by the paper's equation (3)

    expSd_11oh <- sqrt(0.315)
    label("Log-normal (exponential) residual error for 11-OH-THC, SD on the log scale") # Table 2: ERRTHC-OH 0.315 (bootstrap mean 0.299, 95% CI 0.224-0.355); sqrt(0.315) = 0.5612, i.e. 60.8% CV by the paper's equation (3)
  })

  model({
    # ---- 1. Occasion indicators and IOV ------------------------------
    # Two dosing occasions (Methods 2.1: two doses 4 h apart).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)

    iov_mtt    <- oc1 * etaiov_mtt_1    + oc2 * etaiov_mtt_2
    iov_cl_met <- oc1 * etaiov_cl_met_1 + oc2 * etaiov_cl_met_2

    # ---- 2. Individual parameters ------------------------------------
    # Methods 2.2 equation (2): theta_i = theta_pop * exp(eta_i), i.e. all
    # random effects are log-normal, which is what exp(l<param> + eta)
    # encodes. MTT and CLF carry IOV only; V, CLO and CLMET carry IIV
    # only; VMET, QMET and V2MET carry neither.
    mtt        <- exp(lmtt       + iov_mtt)
    cl_met     <- exp(lcl_met    + iov_cl_met)
    cl_nonmet  <- exp(lcl_nonmet + etalcl_nonmet)
    vc         <- exp(lvc        + etalvc)

    cl_11oh    <- exp(lcl_11oh   + etalcl_11oh)
    vc_11oh    <- exp(lvc_11oh)
    q_11oh     <- exp(lq_11oh)
    vp_11oh    <- exp(lvp_11oh)

    # ---- 3. Micro-constants ------------------------------------------
    # Absorption. Figure 3 annotates 'MTT = 3/ktr' for its three-box
    # chain, so ktr = 3/MTT. Methods 2.2 equation (1) writes the same
    # thing as MTT = (n + 1)/ktr 'where n is the number of transit
    # compartments in a model with an absorption compartment between the
    # transit compartments and the central compartment' -- i.e. n = 2
    # transit compartments plus one absorption compartment = the three
    # boxes of Figure 3, which Results 3.2 in turn calls 'absorption via
    # three transit compartments'. All three statements describe the same
    # system: three sequential first-order steps at rate ktr between the
    # dose site and central. In the nlmixr2lib depot + transit<n> layout
    # that is depot -> transit1 -> transit2 -> central, and the general
    # relation ktr = (number of boxes)/MTT gives ktr = 3/MTT, matching
    # Figure 3 exactly. At the typical MTT of 1.35 h, ktr = 2.222 /h.
    ktr <- 3 / mtt

    # THC leaves central by both elimination arms; their sum is the CL/F
    # of 927 L/h that Results 3.2 reports.
    kel <- (cl_met + cl_nonmet) / vc

    kel_11oh <- cl_11oh / vc_11oh
    k12_11oh <- q_11oh  / vc_11oh
    k21_11oh <- q_11oh  / vp_11oh

    # ---- 4. ODE system -----------------------------------------------
    # Absorption chain: three equal first-order ktr steps (Figure 3).
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot    - ktr * transit1
    d/dt(transit2) <-  ktr * transit1 - ktr * transit2

    # THC: one compartment (Results 3.2, 'A one-compartment model with
    # absorption via three transit compartments was the best fit for
    # THC').
    d/dt(central)  <-  ktr * transit2 - kel * central

    # 11-OH-THC: two compartments, fed ONLY by the formation arm CLF/F.
    # The CLO/F arm leaves the system without producing metabolite --
    # this is exactly why the authors kept it. Discussion: removing CLO/F
    # 'resulted in a worse model fit and would require the assumption
    # that all THC is converted to THC-OH, which is not the case.'
    # Concentrations were molecular-weight-normalised to THC equivalents
    # (Methods 2.2), so the formation flux transfers mass 1:1.
    d/dt(central_11oh)     <-  cl_met * central / vc -
                               kel_11oh * central_11oh -
                               k12_11oh * central_11oh +
                               k21_11oh * peripheral1_11oh
    d/dt(peripheral1_11oh) <-  k12_11oh * central_11oh -
                               k21_11oh * peripheral1_11oh

    # ---- 5. Observations and error -----------------------------------
    # Amounts are in mg and volumes in L, so the 1000 converts mg/L to
    # the ug/L of Figure 1, Figure 5 and the Results text.
    Cc      <- 1000 * central      / vc
    Cc_11oh <- 1000 * central_11oh / vc_11oh

    Cc      ~ lnorm(expSd)
    Cc_11oh ~ lnorm(expSd_11oh)
  })
}
