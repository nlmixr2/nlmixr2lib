Kloprogge_2015_artesunate <- function() {
  description <- paste(
    "Simultaneous parent-metabolite population PK model for intravenous and oral",
    "artesunate (AS) and its active metabolite dihydroartemisinin (DHA) in 20",
    "pregnant women with acute uncomplicated Plasmodium falciparum malaria on the",
    "Thailand-Myanmar border, 15 of whom were restudied as healthy volunteers 3",
    "months post-partum (Kloprogge 2015). Both species have two-compartment",
    "disposition and AS is assumed to be completely converted to DHA, so all AS",
    "elimination clearance is metabolic formation of DHA. Oral absorption uses a",
    "six-transit-compartment chain (ktr = (n + 1) / MTT) followed by a first-order",
    "absorption step, at which a first-pass effect splits the absorbed flux: a",
    "fraction fp = 17.1% enters the systemic circulation as AS and the remaining",
    "82.9% is converted pre-systemically and enters the DHA central compartment",
    "directly. Malaria and pregnancy act in opposite directions on the absolute",
    "oral bioavailability of artesunate and on nothing else: acute malaria raises",
    "F by 86.6% and pregnancy lowers it by 23.3%, so neither affects the",
    "disposition of intravenous artesunate. Body weight is applied allometrically",
    "to all clearance (exponent 0.75) and all volume (exponent 1) parameters,",
    "referenced to 46 kg. Residual error is additive on the natural-log scale and",
    "is stratified by route for each analyte. The paper's inter-occasion",
    "variability on the first-pass fraction is omitted from this packaged model",
    "per library convention; see the validation vignette."
  )
  reference <- paste(
    "Kloprogge F, McGready R, Phyo AP, Rijken MJ, Hanpithakpon W, Than HH,",
    "Hlaing N, Zin NT, Day NP, White NJ, Nosten F, Tarning J (2015).",
    "Opposite malaria and pregnancy effect on oral bioavailability of artesunate",
    "- a population pharmacokinetic evaluation.",
    "Br J Clin Pharmacol 80(3):642-653. doi:10.1111/bcp.12660."
  )
  vignette <- "Kloprogge_2015_artesunate"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  # Kloprogge 2015 Table 2 reports four residual-error magnitudes: one per
  # analyte per route of administration (sigma ARS i.v., sigma ARS oral,
  # sigma DHA i.v., sigma DHA oral). The canonical expSd / expSd_<output>
  # matcher has no slot for a route suffix, so the four names are declared
  # here as paper-specific, following the Ahmed_2015_topiramate.R
  # (propSdOral / propSdIv) and AitOudhia_2024_sotatercept.R precedents.
  paper_specific_residual_sds <- c(
    "expSdIv", "expSdOral",
    "expSd_dihydroartIv", "expSd_dihydroartOral"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Amounts are molar (nmol) because Kloprogge 2015
  # Methods state "Plasma concentrations were converted into molar units and
  # modelled simultaneously as their natural logarithms"; the 1:1 molar
  # conversion of AS to DHA therefore needs no molecular-weight factor.
  compartmentData <- list(
    depot                  = list(analyte = "artesunate", units = "nmol", specimen = "administration site", verified = TRUE),
    transit1               = list(analyte = "artesunate", units = "nmol", specimen = "administration site", verified = TRUE),
    transit2               = list(analyte = "artesunate", units = "nmol", specimen = "administration site", verified = TRUE),
    transit3               = list(analyte = "artesunate", units = "nmol", specimen = "administration site", verified = TRUE),
    transit4               = list(analyte = "artesunate", units = "nmol", specimen = "administration site", verified = TRUE),
    transit5               = list(analyte = "artesunate", units = "nmol", specimen = "administration site", verified = TRUE),
    transit6               = list(analyte = "artesunate", units = "nmol", specimen = "administration site", verified = TRUE),
    transit7               = list(analyte = "artesunate", units = "nmol", specimen = "administration site", verified = TRUE),
    central                = list(analyte = "artesunate", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1            = list(analyte = "artesunate", units = "nmol", specimen = "plasma", verified = TRUE),
    central_dihydroart     = list(analyte = "dihydroartemisinin", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1_dihydroart = list(analyte = "dihydroartemisinin", units = "nmol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling with exponents fixed to 0.75 on every clearance",
        "parameter (CL and Q of both artesunate and dihydroartemisinin) and to 1",
        "on every volume parameter (Vc and Vp of both species), per Kloprogge",
        "2015 Results: 'Bodyweight was implemented using allometry on all",
        "clearance (power fixed to 3/4) and volume parameters (power fixed to",
        "1)'. The reference weight is 46 kg, stated in the Table 2 footnote:",
        "'Population parameter estimates are given for a typical non-pregnant",
        "patient with a body weight of 46 kg'. Weight is time-varying across the",
        "two visits in the source design -- the same women were heavier when",
        "pregnant (median 48.0 kg, range 40.0-64.0) than post-partum (median 46.0",
        "kg, range 37.0-52.0; Table 1) -- which is precisely why the authors",
        "included allometry, 'to correct for differences in bodyweight between",
        "the pregnancy and post-partum visit' (Methods)."
      ),
      source_name        = "WT"
    ),
    PREG = list(
      description        = "Pregnancy status",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = second or third trimester of pregnancy; 0 = not pregnant (the same",
        "women restudied 3 months post-partum). Time-varying within subject: the",
        "source design re-enrolled 15 of the 20 pregnant patients as healthy",
        "post-partum volunteers, so a subject carries PREG = 1 on the pregnancy",
        "visit and PREG = 0 on the post-partum visit. Applied as a fractional",
        "multiplicative effect on the absolute oral bioavailability of artesunate",
        "only: fdepot = fdepot_typ * (1 + e_preg_fdepot * PREG) with",
        "e_preg_fdepot = -0.233, i.e. 23.3% lower bioavailability during",
        "pregnancy (Kloprogge 2015 Table 2 'Pregnancy effect on F (%)' = -23.3,",
        "%RSE 18.7). Pregnancy was also significant on dihydroartemisinin",
        "elimination clearance (+27%, dOFV = -45.1, a larger drop than the",
        "-22.5 obtained on bioavailability) but the authors carried the",
        "bioavailability parameterisation forward because covariate modelling on",
        "the intravenous data alone could not detect a pregnancy effect on",
        "elimination clearance, in agreement with the non-compartmental analysis",
        "(Results and Discussion). Trimester as a categorical covariate and",
        "estimated gestational age as a continuous covariate were both tested and",
        "were unstable and not superior to the binary indicator."
      ),
      source_name        = "PREG"
    ),
    DIS_MALARIA_ACUTE = list(
      description        = "Acute symptomatic phase of a Plasmodium falciparum malaria episode",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = the record falls in the acute phase of the malaria episode (study",
        "days 1 and 2); 0 = convalescent or healthy (study day 7 of the treated",
        "episode, and every record from the post-partum healthy-volunteer visit).",
        "Time-varying WITHIN subject, not a cohort-membership indicator: each",
        "pregnant patient contributes acute records on days 1-2 and convalescent",
        "records on day 7 of the same admission. This encoding follows the",
        "central modelling assumption Kloprogge 2015 state in Methods -- 'that",
        "pregnant malaria patients during their convalescent phase (day 7 during",
        "the first visit) have a similar disease state (i.e. healthy) to",
        "post-partum healthy volunteers' -- which is what allows the malaria and",
        "pregnancy effects to be dissected. Applied as a fractional multiplicative",
        "effect on the absolute oral bioavailability of artesunate only:",
        "fdepot = fdepot_typ * (1 + e_dis_malaria_acute_fdepot * DIS_MALARIA_ACUTE)",
        "with e_dis_malaria_acute_fdepot = 0.866, i.e. 86.6% higher bioavailability",
        "during acute malaria (Kloprogge 2015 Table 2 'Disease effect on F (%)' =",
        "86.6, %RSE 3.50; dOFV = -66.4). The authors attribute the increase to a",
        "reduced first-pass effect during acute infection. Alternative disease",
        "descriptors were tested and not retained: a mildly/moderately-unwell",
        "dichotomy based on raised creatinine, blood urea nitrogen or liver",
        "function tests (or fever > 37.5 degC with tachycardia > 100 beats/min),",
        "and a disease-severity count of fulfilled criteria (Methods)."
      ),
      source_name        = "DISEASE"
    ),
    ROUTE_IV = list(
      description        = "Intravenous (vs oral) route of administration for the record",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = intravenous artesunate, 0 = ORAL artesunate. The reference category",
        "here is oral, not subcutaneous. Per-dose-record / per-observation",
        "indicator: every patient received both routes in the source design",
        "(group 1 had intravenous artesunate on admission then oral for 6 days,",
        "group 2 had oral on admission, intravenous on day 2, then oral for 5",
        "days). Used only to select the residual-error magnitude -- Kloprogge",
        "2015 Table 2 reports four log-scale additive residual SDs, one per",
        "analyte per route (sigma ARS i.v. 0.856, sigma ARS oral 1.16, sigma DHA",
        "i.v. 0.333, sigma DHA oral 0.793) -- exactly the residual-only role of",
        "ROUTE_IV in Zierhut_2008_osteoprotegerin.R, Wang_2021_pertuzumab.R and",
        "Kuroda_2024_quinidine_horse.R. No structural parameter differs by route.",
        "When simulating, set ROUTE_IV = 1 and dose into 'central' (a zero-order",
        "input of the fixed 0.0167 h = 1 min duration applies automatically); set",
        "ROUTE_IV = 0 and dose into 'depot' for the oral route."
      ),
      source_name        = "ROUTE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 20L,
    n_pregnant     = 20L,
    n_postpartum   = 15L,
    n_studies      = 1L,
    n_observations = "1571 plasma samples (920 during the pregnancy visit, 651 during the post-partum visit; Table 1). More than 45% of the artesunate observations were below the 1.2 ng/mL limit of quantification and were handled with the M3 likelihood method rather than being discarded.",
    age_range      = "not reported (adult women of childbearing age)",
    weight_range   = "40.0-64.0 kg pregnant (median 48.0); 37.0-52.0 kg post-partum (median 46.0) (Table 1)",
    weight_median  = "46 kg is the typical-value reference weight of Table 2",
    sex_female_pct = 100,
    race_ethnicity = "Karen and Burmese migrant and refugee women on the north-western border of Thailand and Myanmar (Shoklo Malaria Research Unit)",
    disease_state  = paste(
      "Uncomplicated Plasmodium falciparum malaria in the second or third",
      "trimester of pregnancy, with haematocrit not lower than 25%. Estimated",
      "gestational age 25.7 weeks (range 14.0-38.0) by dating ultrasound, 10",
      "women in the second and 10 in the third trimester; 7 classified mildly and",
      "13 moderately unwell at admission (Table 1). The same women were restudied",
      "3 months post-partum as healthy volunteers, the visit being postponed if",
      "malaria or any other illness was detected. 17 women delivered healthy",
      "singleton babies at a mean 39.2 weeks (range 35.6-41.5) and 3 were lost to",
      "follow-up."
    ),
    ga_range       = "14.0-38.0 weeks estimated gestational age at enrolment (median 25.7)",
    dose_range     = paste(
      "4 mg/kg artesunate daily for 7 days, with one of the seven doses given",
      "intravenously: group 1 received intravenous artesunate on admission",
      "followed by oral artesunate on the next 6 days; group 2 received oral",
      "artesunate on admission, intravenous artesunate on day 2, then oral",
      "artesunate for 5 days. Total artesunate dose 27.9 mg/kg (range 26.8-28.6)",
      "during the pregnancy visit and 27.4 mg/kg (4.08-29.0) post-partum",
      "(Table 1). Intravenous and oral artesunate were manufactured by the Guilin",
      "Pharmaceutical Factory and repackaged by Atlantic Pharmaceuticals; the same",
      "drug lots were reserved for the post-partum visit."
    ),
    sampling       = "Days 1 and 2: 0, 5, 15, 30, 60, 120, 180, 240 and 360 min after an intravenous dose, or 0, 15, 30, 60, 90, 120, 180, 240 and 360 min after an oral dose. Day 7: 0, 60, 120, 240 and 360 min after dose in all patients.",
    regions        = "Thailand-Myanmar border (Shoklo Malaria Research Unit clinics), April 2008 to March 2009",
    notes          = paste(
      "Re-analysis of a previously published pharmacokinetic study; the",
      "non-compartmental results of the same data appear in reference [8] of the",
      "source paper. Ethical approval TM-IR 029/2005 (Mahidol University) and",
      "OXTREC 007-05. Baseline demographics are Table 1 of Kloprogge 2015."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Absorption. Kloprogge 2015 Table 2. F and BIO are reported as
    # percentages and are converted to fractions here.
    # ---------------------------------------------------------------------
    lfdepot <- log(0.494)
    label("Absolute oral bioavailability of artesunate at the reference covariates (fraction)")  # Table 2 'F (%)' = 49.4 (%RSE 3.53); typical non-pregnant, non-acute patient of 46 kg
    logitfp <- logit(0.171)
    label("Logit of the first-pass fraction BIO, the share of the absorbed oral dose entering the circulation as artesunate (logit scale, -)")  # Table 2 'BIO (%)' = 17.1 (%RSE 6.98); Figure 1 defines bio as 'fraction of dose absorbed as AS'
    lmtt <- log(0.407)
    label("Mean transit time of the oral absorption chain (h)")  # Table 2 'MTT (h)' = 0.407 (%RSE 7.24)
    lka <- log(1.57)
    label("First-order absorption rate constant out of the final transit compartment (1/h)")  # Table 2 'ka (h-1)' = 1.57 (%RSE 8.59)
    ldur <- fixed(log(0.0167))
    label("Zero-order duration of the intravenous artesunate infusion (h)")  # Table 2 'DUR (h)' = 0.0167 (fixed); 1 min

    # ---------------------------------------------------------------------
    # Artesunate (parent) disposition, two compartments. Table 2.
    # ---------------------------------------------------------------------
    lvc <- log(8.80)
    label("Artesunate central volume of distribution (L)")  # Table 2 'VcART (l)' = 8.80 (%RSE 5.32)
    lq <- log(7.51)
    label("Artesunate inter-compartmental clearance (L/h)")  # Table 2 'QART (l h-1)' = 7.51 (%RSE 12.6)
    lvp <- log(2.43)
    label("Artesunate peripheral volume of distribution (L)")  # Table 2 'VpART (l)' = 2.43 (%RSE 13.8)
    lcl <- log(170)
    label("Artesunate elimination clearance, entirely metabolic conversion to dihydroartemisinin (L/h)")  # Table 2 'CLART (l h-1)' = 170 (%RSE 6.75)

    # ---------------------------------------------------------------------
    # Dihydroartemisinin (metabolite) disposition, two compartments. Table 2.
    # ---------------------------------------------------------------------
    lvc_dihydroart <- log(44.3)
    label("Dihydroartemisinin central volume of distribution (L)")  # Table 2 'VcDHA (l)' = 44.3 (%RSE 6.35)
    lq_dihydroart <- log(16.8)
    label("Dihydroartemisinin inter-compartmental clearance (L/h)")  # Table 2 'QDHA (l h-1)' = 16.8 (%RSE 9.7)
    lvp_dihydroart <- log(20.7)
    label("Dihydroartemisinin peripheral volume of distribution (L)")  # Table 2 'VpDHA (l)' = 20.7 (%RSE 9.35)
    lcl_dihydroart <- log(60.9)
    label("Dihydroartemisinin elimination clearance (L/h)")  # Table 2 'CLDHA (l h-1)' = 60.9 (%RSE 4.54)

    # ---------------------------------------------------------------------
    # Allometric exponents, both held fixed by the authors and applied to
    # every clearance and every volume parameter of both species.
    # ---------------------------------------------------------------------
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on all clearance parameters (unitless)")  # Results: allometry on all clearance, 'power fixed to 3/4'
    e_wt_vc <- fixed(1)
    label("Allometric exponent on all volume parameters (unitless)")  # Results: allometry on all volume parameters, 'power fixed to 1'

    # ---------------------------------------------------------------------
    # Covariate effects. Both act on oral bioavailability only, in the
    # fractional multiplicative form F = F_typ * (1 + effect * indicator),
    # which reproduces the paper's own post hoc F values (see the vignette
    # source-trace table).
    # ---------------------------------------------------------------------
    e_dis_malaria_acute_fdepot <- 0.866
    label("Fractional change in artesunate oral bioavailability during acute malaria (unitless)")  # Table 2 'Disease effect on F (%)' = 86.6 (%RSE 3.50)
    e_preg_fdepot <- -0.233
    label("Fractional change in artesunate oral bioavailability during pregnancy (unitless)")  # Table 2 'Pregnancy effect on F (%)' = -23.3 (%RSE 18.7); the minus sign is stated in Results and Conclusions ('pregnancy decreased oral bioavailability by 23%')

    # ---------------------------------------------------------------------
    # Inter-individual variability. Table 2 reports each random effect as a
    # %CV back-transformed from the NONMEM variance by the footnote formula
    # 100 * sqrt(exp(estimate) - 1), so the variance is recovered as
    # log(1 + (CV/100)^2). Table 2 leaves the IIV cell blank for VcART,
    # CLART, VcDHA, QDHA and VpDHA, so no eta is estimated on those five.
    # ---------------------------------------------------------------------
    etalfdepot ~ log(1 + 0.271^2)           # Table 2 'F (%)' IIV = 27.1 %CV (%RSE 20.9); variance 0.0709
    etalmtt ~ log(1 + 0.310^2)              # Table 2 'MTT (h)' IIV = 31.0 %CV (%RSE 20.2); variance 0.0918
    etalka ~ log(1 + 0.368^2)               # Table 2 'ka (h-1)' IIV = 36.8 %CV (%RSE 18.9); variance 0.1270
    etalq ~ log(1 + 0.0660^2)               # Table 2 'QART (l h-1)' IIV = 6.60 %CV (%RSE 23.2); variance 0.00435
    etalvp ~ log(1 + 0.633^2)               # Table 2 'VpART (l)' IIV = 63.3 %CV (%RSE 12.8); variance 0.3370
    etalcl_dihydroart ~ log(1 + 0.104^2)    # Table 2 'CLDHA (l h-1)' IIV = 10.4 %CV (%RSE 22.7); variance 0.0108

    # ---------------------------------------------------------------------
    # Residual error. Kloprogge 2015 Results: 'Residual variability was best
    # described using an additive error model on the log-transformed data',
    # which is a log-normal residual in nlmixr2's linear space, encoded with
    # lnorm(). The four Table 2 sigma values are NONMEM $SIGMA variance
    # estimates (they sit in the 'Population estimates' column with no
    # back-transformation applied, unlike the IIV column), so the SD is the
    # square root. The magnitudes are consistent with the sibling
    # Morris_2011_artesunate.R, whose Table 2 reports the same quantity as
    # 0.696 for artesunate and 0.174 for dihydroartemisinin.
    # ---------------------------------------------------------------------
    expSdIv <- sqrt(0.856)
    label("Log-scale additive residual SD for artesunate after intravenous dosing")  # Table 2 'sigma ARS i.v.' = 0.856 (variance); SD = 0.925
    expSdOral <- sqrt(1.16)
    label("Log-scale additive residual SD for artesunate after oral dosing")  # Table 2 'sigma ARS oral' = 1.16 (variance); SD = 1.077
    expSd_dihydroartIv <- sqrt(0.333)
    label("Log-scale additive residual SD for dihydroartemisinin after intravenous dosing")  # Table 2 'sigma DHA i.v.' = 0.333 (variance); SD = 0.577
    expSd_dihydroartOral <- sqrt(0.793)
    label("Log-scale additive residual SD for dihydroartemisinin after oral dosing")  # Table 2 'sigma DHA oral' = 0.793 (variance); SD = 0.891
  })

  model({
    # ---------------------------------------------------------------------
    # Absolute oral bioavailability of artesunate. The two covariates act
    # here and nowhere else, which is the paper's central finding: malaria
    # and pregnancy leave the disposition of intravenous artesunate
    # untouched and move only the oral bioavailability, in opposite
    # directions (Results, Discussion, Conclusions).
    # ---------------------------------------------------------------------
    f_oral <- exp(lfdepot + etalfdepot) *
      (1 + e_dis_malaria_acute_fdepot * DIS_MALARIA_ACUTE) *
      (1 + e_preg_fdepot * PREG)

    # First-pass fraction. fp is the share of the absorbed dose that reaches
    # the systemic circulation as artesunate; the complement (1 - fp) is
    # hydrolysed pre-systemically -- at gastric pH, by plasma esterases and
    # by hepatic CYP2A6 (Discussion) -- and appears directly as
    # dihydroartemisinin. Held on the logit scale so it cannot leave [0, 1],
    # following Bertrand_2011_S33138.R, which encodes the same
    # parent-presystemic-fraction construct.
    fp <- expit(logitfp)

    # ---------------------------------------------------------------------
    # Individual PK parameters with allometric weight scaling on 46 kg.
    # ---------------------------------------------------------------------
    mtt <- exp(lmtt + etalmtt)
    ka <- exp(lka + etalka)
    dur_iv <- exp(ldur)

    vc <- exp(lvc) * (WT / 46)^e_wt_vc
    q <- exp(lq + etalq) * (WT / 46)^e_wt_cl
    vp <- exp(lvp + etalvp) * (WT / 46)^e_wt_vc
    cl <- exp(lcl) * (WT / 46)^e_wt_cl

    vc_dihydroart <- exp(lvc_dihydroart) * (WT / 46)^e_wt_vc
    q_dihydroart <- exp(lq_dihydroart) * (WT / 46)^e_wt_cl
    vp_dihydroart <- exp(lvp_dihydroart) * (WT / 46)^e_wt_vc
    cl_dihydroart <- exp(lcl_dihydroart + etalcl_dihydroart) * (WT / 46)^e_wt_cl

    # ---------------------------------------------------------------------
    # Transit-chain rate constant. Figure 1 caption: 'ktr was calculated as
    # (transit compartments (n) + 1) / mean transit time', with n = 6 fixed
    # (Table 2). Seven states -- depot plus transit1..transit6 -- therefore
    # empty at ktr, and the eighth (transit7) is the absorption compartment
    # that empties at ka. This is the same chain idiom as
    # Bukkems_2021_raltegravir.R.
    # ---------------------------------------------------------------------
    ktr <- (6 + 1) / mtt

    # Two-compartment micro-constants for each species.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    kel_dihydroart <- cl_dihydroart / vc_dihydroart
    k12_dihydroart <- q_dihydroart / vc_dihydroart
    k21_dihydroart <- q_dihydroart / vp_dihydroart

    # ---------------------------------------------------------------------
    # ODE system.
    #
    # Oral:  depot -> transit1 -> ... -> transit6 (all at ktr) -> transit7,
    #        which empties at ka and splits: a fraction fp into artesunate
    #        central and (1 - fp) directly into dihydroartemisinin central.
    # I.v.:  dose records target `central` and are delivered as a zero-order
    #        input over dur_iv, bypassing the whole absorption chain.
    #
    # Artesunate is assumed to be completely metabolised to
    # dihydroartemisinin (Methods), so the whole of kel * central is
    # formation of dihydroartemisinin rather than true elimination. Amounts
    # are molar, and the conversion is 1:1 in moles, so no molecular-weight
    # factor appears at the formation step.
    # ---------------------------------------------------------------------
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot    - ktr * transit1
    d/dt(transit2) <-  ktr * transit1 - ktr * transit2
    d/dt(transit3) <-  ktr * transit2 - ktr * transit3
    d/dt(transit4) <-  ktr * transit3 - ktr * transit4
    d/dt(transit5) <-  ktr * transit4 - ktr * transit5
    d/dt(transit6) <-  ktr * transit5 - ktr * transit6
    d/dt(transit7) <-  ktr * transit6 - ka * transit7

    d/dt(central) <- ka * fp * transit7 -
      kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    d/dt(central_dihydroart) <- ka * (1 - fp) * transit7 +
      kel * central -
      kel_dihydroart * central_dihydroart -
      k12_dihydroart * central_dihydroart +
      k21_dihydroart * peripheral1_dihydroart
    d/dt(peripheral1_dihydroart) <- k12_dihydroart * central_dihydroart -
      k21_dihydroart * peripheral1_dihydroart

    # Oral bioavailability applies to the depot only; intravenous dose
    # records target `central` and are unaffected by it.
    f(depot) <- f_oral
    dur(central) <- dur_iv

    # ---------------------------------------------------------------------
    # Observations. Amounts are nmol and volumes are L, so both
    # concentrations are nmol/L. Multiply by the molecular weight and divide
    # by 1000 to compare against the ng/mL of Table 3: artesunate
    # 384.42 g/mol, dihydroartemisinin 284.35 g/mol.
    # ---------------------------------------------------------------------
    Cc <- central / vc
    Cc_dihydroart <- central_dihydroart / vc_dihydroart

    # Route-stratified residual magnitudes (Table 2). ROUTE_IV selects the
    # intravenous value; the oral value is the reference.
    expSd <- expSdIv * ROUTE_IV + expSdOral * (1 - ROUTE_IV)
    expSd_dihydroart <- expSd_dihydroartIv * ROUTE_IV +
      expSd_dihydroartOral * (1 - ROUTE_IV)

    Cc ~ lnorm(expSd)
    Cc_dihydroart ~ lnorm(expSd_dihydroart)
  })
}
