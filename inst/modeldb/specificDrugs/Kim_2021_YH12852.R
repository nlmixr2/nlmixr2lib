Kim_2021_YH12852 <- function() {
  description <- "Joint population PK-PD model of the 5-HT4 receptor agonist YH12852 in healthy adults and adults with functional constipation (Kim 2021). The PK is a two-compartment model with first-order absorption, body weight on the peripheral volume, and between-occasion variability on CL/F, V2/F and Ka. The prokinetic effect is a semi-mechanistic three-compartment model of the 13C-Spirulina gastric emptying breath test (GEBT). 13C in the test meal leaves the stomach at rate (K45 + SLP * Cc), and the fraction FC13 is absorbed into the systemic circulation. It then moves to the lung (K56) and is exhaled (Kout = K56). The observed kPCD is the per-minute percent of the 13C dose exhaled, multiplied by 1000. Baseline gastric emptying t10 scales the drug slope SLP."
  reference <- "Kim S, Lee HA, Jang SB, Lee H. A population pharmacokinetic-pharmacodynamic model of YH12852, a highly selective 5-hydroxytryptamine 4 receptor agonist, in healthy subjects and patients with functional constipation. CPT Pharmacometrics Syst Pharmacol. 2021;10(8):902-913. doi:10.1002/psp4.12664"
  vignette <- "Kim_2021_YH12852"
  units <- list(time = "h", dosing = "mg", concentration = "pg/mL")

  # Drug doses go to `depot` (mg); each GEBT test meal is a 100,000 kPCD-unit
  # bolus of 13C into `stomach` given at the time the meal is finished
  # (Methods: "the initial amount of 13C in compartment 4 was set to 100,000").
  dosing <- c("depot", "stomach")

  # The two downstream 13C states have no canonical compartment name: they hold
  # the absorbed 13C label (not drug) in the systemic circulation and in the
  # lung (Kim 2021 Figure 1b; control stream GEINTER / GECENT).
  paper_specific_compartments <- c("c13_blood", "c13_lung")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the peripheral volume V3/F, (WT / 59.25)^0.93. The 59.25 kg reference is from the supplementary NONMEM control stream (Supplementary Method 2, V3 line); the Methods text rounds the median to 59.2 kg. Table 1 labels the row 'BWT effect on V2/F' but the Results text ('Body weight ... significant covariates on V3/F') and the control stream both put it on V3/F.",
      source_name = "WTKG"
    ),
    GE_T10_BL = list(
      description = "Baseline (pre-treatment) gastric emptying t10: time for 10% of the 13C test meal to leave the stomach, from the drug-free GEBT.",
      units = "min",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the drug slope SLP, (GE_T10_BL / 30)^3.57, reference 30 min (Methods; control stream HT10/30). Kim 2021 derived it by interpolating the gastric emptying fractions that the Szarka 2008 regressions (Table S1) predict from the baseline kPCD values, sex and BMI. MLD cohort mean +/- SD 30.3 +/- 15.5 min, range 11.1-59.4 min (Table 2). Time-fixed per subject.",
      source_name = "HT10"
    ),
    OCC = list(
      description = "Occasion index for between-occasion variability on CL/F, V2/F and Ka. 1 = Day 1; 2 = any later sampling (the Day 7 GEBT and the Day 14 profile).",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Kim 2021 Methods: 'Occasion was defined as a set of sampling times clearly separated between two adjacent occasions (i.e., 1 for Day 1 and 2 the other)'. Control stream: $OMEGA BLOCK(3) for occasion 1 followed by $OMEGA BLOCK(3) SAME for occasion 2. OCC = 0 (or any other value) switches the between-occasion variability off.",
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened (typical value 27.6 years) and not retained; Kim 2021 Results."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained; Kim 2021 Results. Sex still enters the Szarka 2008 regressions that turn simulated kPCD values into gastric emptying fractions (Table S1), which are applied outside the model."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened (typical value 22.0 kg/m^2) and not retained; Kim 2021 Results."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened (typical value 15 U/L) and not retained; Kim 2021 Methods and Results."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened (typical value 10 U/L) and not retained; Kim 2021 Methods and Results."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units = "mmol/L",
      type = "continuous",
      notes = "Screened (typical value 10 mmol/L) and not retained; Kim 2021 Methods and Results."
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "YH12852",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "YH12852",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "YH12852",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    stomach = list(
      analyte = "13C (GEBT test-meal label)",
      units = "percent of the 13C dose x 1000",
      specimen = "administration site",
      verified = TRUE
    ),
    c13_blood = list(
      analyte = "13C (absorbed GEBT label)",
      units = "percent of the 13C dose x 1000",
      specimen = "not applicable",
      verified = TRUE
    ),
    c13_lung = list(
      analyte = "13C (absorbed GEBT label)",
      units = "percent of the 13C dose x 1000",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 49,
    n_studies = 1,
    n_observations = "1287 YH12852 plasma concentrations (49 subjects) and 196 kPCD values (14 MLD-cohort subjects)",
    age_range = "19-53 years",
    age_mean = "27.3 years overall (MD cohort 28.6 +/- 7.7; MLD cohort 24.2 +/- 3.6)",
    weight_range = "45.9-78.8 kg",
    weight_mean = "MD cohort 60.4 +/- 8.2 kg; MLD cohort 58.2 +/- 8.1 kg",
    bmi_range = "18.2-25.0 kg/m^2",
    sex_female_pct = 71.4,
    race_ethnicity = "Korean adults (single-centre study in Seoul); race is not tabulated.",
    disease_state = "Healthy adults reporting 3 or fewer spontaneous bowel movements per week, and adults with functional constipation by the Rome III criteria (17 of the 35 MD-cohort subjects). The PD (GEBT) data come only from the 14 healthy MLD-cohort subjects.",
    dose_range = "YH12852 0.05, 0.1, 0.3, 0.5, 1, 2 or 3 mg orally once daily after breakfast for 14 days",
    regions = "South Korea",
    notes = "Randomized, double-blind, placebo-controlled phase I/IIa trial NCT02538367. The multiple-dose (MD) cohort received 0.3-3 mg and the multiple-low-dose (MLD) cohort 0.05 or 0.1 mg. The GEBT was done in the MLD cohort at baseline and on Day 7, with breath samples 45-240 min after the meal. Baseline demographics are in Kim 2021 Table 2. The plasma LLOQ was 30 pg/mL."
  )

  ini({
    # ---- PK (Kim 2021 Table 1, final PK-PD model) ---------------------------
    lcl <- log(88.8)
    label("Apparent clearance CL/F (L/h)") # Table 1: CL/F 88.8 L/hr
    lvc <- log(1380.2)
    label("Apparent central volume V2/F (L)") # Table 1: V2/F 1380.2 L
    lvp <- log(989.7)
    label("Apparent peripheral volume V3/F for a 59.25 kg subject (L)") # Table 1: V3/F 989.7 L
    lq <- log(137.4)
    label("Apparent intercompartmental clearance Q/F (L/h)") # Table 1: Q 137.4 L/hr
    lka <- log(0.48)
    label("First-order absorption rate constant Ka (1/h)") # Table 1: Ka 0.48 1/hr
    e_wt_vp <- 0.93
    label("Power exponent of body weight on V3/F (unitless)") # Table 1: 'BWT effect on V2/F' 0.93; applied to V3 per Results text and control stream

    # ---- GEBT 13C PD (Kim 2021 Table 1; Equations 1-4) ----------------------
    lf13c <- log(0.22)
    label("Fraction of the test-meal 13C that is absorbed, FC13 (fraction)") # Table 1: FC13 0.22
    lkge <- log(0.39)
    label("Drug-free rate constant of 13C leaving the stomach, K45 (1/h)") # Table 1: K45 0.39 1/hr
    lk13c <- log(0.78)
    label("Rate constant of 13C transfer blood->lung and lung->air, K56 = Kout (1/h)") # Table 1: K56, Kout 0.78 1/hr
    lslp_kge <- log(0.0009)
    label("Linear slope of plasma YH12852 on the stomach emptying rate, SLP (1/h per pg/mL)") # Table 1: SLP 0.0009
    e_t10_slp <- 3.57
    label("Power exponent of baseline gastric emptying t10 on SLP (unitless)") # Table 1: t10 effect on SLP 3.57

    # ---- Between-subject variability ----------------------------------------
    # Table 1 reports each IIV/IOV as a percentage only. Converted as
    # omega^2 = (pct/100)^2, i.e. the percentage read as 100*sqrt(omega^2).
    # The supplementary control stream's initial $OMEGA for IOV CL, 0.0822,
    # gives 100*sqrt(0.0822) = 28.7%, the Table 1 value exactly. The
    # log(1 + CV^2) reading gives 29.3%. See the vignette for the full reasoning.
    # The PK etas were fitted as a BLOCK(4) (CL, V2, Q, V3), but the final
    # covariances are not reported, so they are encoded as independent.
    etalcl ~ 0.152881 # Table 1 IIV CL/F 39.1%
    etalvc ~ 0.037636 # Table 1 IIV V2/F 19.4%
    etalvp ~ 0.104976 # Table 1 IIV V3/F 32.4%
    etalq ~ 0.086436 # Table 1 IIV Q 29.4%
    etalka ~ 0.0025 # Table 1 IIV Ka 5.0% (shrinkage 90.1%)
    etalf13c ~ 0.010609 # Table 1 IIV FC13 10.3%
    etalk13c ~ 0.051984 # Table 1 IIV K56, Kout 22.8%
    etalslp_kge ~ 1.221025 # Table 1 IIV SLP 110.5%

    # ---- Between-occasion variability (2 occasions, shared magnitude) -------
    # rxode2 has no NONMEM-style occasion level, so each occasion gets its own
    # eta. The second occasion's variance is fixed to the first, which is the
    # $OMEGA BLOCK(3) SAME idiom. The BLOCK(3) covariances are not reported.
    etaiov_cl_1 ~ 0.082369 # Table 1 IOV CL/F 28.7%
    etaiov_cl_2 ~ fixed(0.082369) # shared magnitude ($OMEGA BLOCK(3) SAME)
    etaiov_vc_1 ~ 0.200704 # Table 1 IOV V2/F 44.8%
    etaiov_vc_2 ~ fixed(0.200704) # shared magnitude ($OMEGA BLOCK(3) SAME)
    etaiov_ka_1 ~ 0.259081 # Table 1 IOV Ka 50.9%
    etaiov_ka_2 ~ fixed(0.259081) # shared magnitude ($OMEGA BLOCK(3) SAME)

    # ---- Residual error -----------------------------------------------------
    # Control stream $ERROR: W = sqrt(THETA_prop^2 * IPRED^2 + THETA_add^2)
    # with $SIGMA 1 FIX, on the untransformed scale. The additive SDs are
    # THETA(14) and THETA(15), both (0, 0.001) FIX.
    propSd <- 0.177
    label("Proportional residual error, YH12852 plasma concentration (fraction)") # Table 1: 17.7%
    addSd <- fixed(0.001)
    label("Additive residual error, YH12852 plasma concentration (pg/mL)") # control stream THETA(14) 0.001 FIX
    propSd_kpcd <- 0.125
    label("Proportional residual error, kPCD (fraction)") # Table 1: 12.5%
    addSd_kpcd <- fixed(0.001)
    label("Additive residual error, kPCD (percent of dose x 1000 per min)") # control stream THETA(15) 0.001 FIX
  })

  model({
    # ---- 1. Between-occasion variability -----------------------------------
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2
    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2

    # ---- 2. Individual parameters ------------------------------------------
    cl <- exp(lcl + etalcl + iov_cl)
    vc <- exp(lvc + etalvc + iov_vc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp) * (WT / 59.25)^e_wt_vp
    ka <- exp(lka + etalka + iov_ka)

    f13c <- exp(lf13c + etalf13c)
    kge <- exp(lkge)
    k13c <- exp(lk13c + etalk13c)
    slp_kge <- exp(lslp_kge + etalslp_kge) * (GE_T10_BL / 30)^e_t10_slp

    # ---- 3. Micro-constants ------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- 4. Plasma concentration driving the PD ----------------------------
    # Doses in mg and volumes in L give mg/L; 1 mg/L = 1e6 pg/mL. The control
    # stream uses CONC = A(2)/V2 directly with pg/mL data, so its AMT was in ng.
    Cc <- 1e6 * central / vc

    # ---- 5. ODE system -----------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Kim 2021 Equations 1-3. All 13C that leaves the stomach is lost from it,
    # but only the fraction FC13 reaches the systemic circulation. The control
    # stream's $DES defines EFF1 = SLP * A(2)/V2 and then writes DADT(4) and
    # DADT(5) with EFF; Equations 1-2 confirm that the term is SLP * CONC.
    d/dt(stomach) <- -(kge + slp_kge * Cc) * stomach
    d/dt(c13_blood) <- (kge + slp_kge * Cc) * f13c * stomach - k13c * c13_blood
    d/dt(c13_lung) <- k13c * c13_blood - k13c * c13_lung

    # ---- 6. Observations ---------------------------------------------------
    # Equation 4: rate constants are per hour and kPCD is per minute, hence /60.
    kpcd <- k13c * c13_lung / 60

    Cc ~ add(addSd) + prop(propSd)
    kpcd ~ add(addSd_kpcd) + prop(propSd_kpcd)
  })
}
