Comisar_2025_rimegepant <- function() {
  description <- paste(
    "Two-compartment population PK model for oral rimegepant with four",
    "transit absorption compartments and a terminal first-order (ka) step",
    "into the central compartment, jointly fitted to pooled adult and",
    "pediatric (6 to <12 years) phase 1 data with estimated allometric body",
    "weight exponents on clearances and volumes."
  )
  reference <- paste(
    "Comisar CM, Hughes JH, Mo G, Bhardwaj R, Jakate A, Lim CN, Liu J.",
    "Exposure Matching Using Population Pharmacokinetic Modeling and",
    "Simulation to Support Rimegepant Dose Selection for Pediatric Patients",
    "With Migraine. Clin Transl Sci. 2025;18(10):e70360.",
    "doi:10.1111/cts.70360.",
    "Covariate-model functional forms (power model for continuous",
    "covariates, theta_TV,REF * (1 + theta_x)^X for categorical covariates)",
    "and the F1 = (dose / 10 mg)^theta parameterisation are documented in",
    "the upstream adult-only model this analysis extends:",
    "Comisar CM, Hughes JH, Francis J, et al. Population Pharmacokinetic",
    "Modeling of the Oral Calcitonin Gene-Related Peptide Receptor",
    "Antagonist Rimegepant in Adults. CPT Pharmacometrics Syst Pharmacol.",
    "2025;14(8):1332-1345. doi:10.1002/psp4.70051."
  )
  vignette <- "Comisar_2025_rimegepant"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling on all four disposition parameters, centred at",
        "70 kg. The 70 kg reference is stated in the Comisar 2025 pediatric",
        "supplement (Supplemental methods for the sensitivity analysis:",
        "'allometric coefficients were fixed to literature values to",
        "maintain model stability in the absence of any information to",
        "center the covariate effect at 70 kg'). Unlike the upstream",
        "adult-only model, which fixed the exponents at the canonical 0.75",
        "and 1, this combined pediatric+adult model ESTIMATES both",
        "exponents: 0.575 (95% CI 0.413 to 0.737) on CL/F and Q/F and 1.18",
        "(0.988 to 1.372) on Vc/F and Vp/F. Analysis-population range 23.2",
        "to 134 kg (Table 1)."
      ),
      source_name        = "WT"
    ),
    FED = list(
      description        = "Fed state at the time of the dose",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = paste(
        "Per-dose-record indicator; the food-effect studies were crossover",
        "designs, so a single subject contributes both levels. Two",
        "multiplicative effects, both in Table 2: relative bioavailability",
        "F1 * (1 - 0.331 * FED) and transit rate ktr * (1 - 0.698 * FED).",
        "Meal definitions were study-specific (BHV3000-113 Part 2 food",
        "effect on the 75 mg ODT; BHV3000-120 a low-fat, low-calorie meal;",
        "supplement Table S1), so the general fed-vs-fasted FED canonical",
        "applies rather than FED_LOWFAT or FED_HIGHFAT. All 20 pediatric",
        "participants were dosed fasted."
      ),
      source_name        = "Fed"
    ),
    FORM_ODT = list(
      description        = "Orally disintegrating tablet formulation",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (immediate-release oral tablet)",
      notes              = paste(
        "Per-dose-record indicator. Multiplicative effect on the transit",
        "rate constant only: ktr * (1 + 0.470 * FORM_ODT) (Table 2, 'ODT on",
        "k tr'). The ODT is the marketed 75 mg presentation and the only",
        "formulation given to pediatric participants. Formulation does not",
        "affect bioavailability in this model. Reference category is the",
        "oral tablet (73.1% of upstream records); the capsule arm carries",
        "the separate FORM_CAPSULE indicator, so tablet is the state in",
        "which both indicators are 0."
      ),
      source_name        = "ODT"
    ),
    FORM_CAPSULE = list(
      description        = "Capsule formulation",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (immediate-release oral tablet)",
      notes              = paste(
        "Per-dose-record indicator. Multiplicative effect on the transit",
        "rate constant only: ktr * (1 + 2.19 * FORM_CAPSULE) (Table 2,",
        "'Capsule formulation on k tr'), i.e. a 3.19-fold faster transit",
        "than the tablet reference. Capsules were an early adult phase 1",
        "presentation (4.9% of upstream records) and were never given to",
        "pediatric participants. Paired with FORM_ODT: tablet is the",
        "reference when both indicators are 0. Bioavailability is",
        "unaffected, unlike the several capsule-vs-tablet models that put",
        "the effect on F."
      ),
      source_name        = "Capsule"
    ),
    HEPIMP_MOD = list(
      description        = "Moderate hepatic impairment (Child-Pugh B)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function, mild impairment, or not listed)",
      notes              = paste(
        "Multiplicative effect on apparent clearance: CL/F * (1 - 0.203 *",
        "HEPIMP_MOD) (Table 2), a 20.3% reduction. Mild hepatic impairment",
        "was tested and was NOT a significant covariate, so mild-impairment",
        "subjects sit in the reference category together with subjects",
        "having normal hepatic function. Six adults (1.4%) were moderately",
        "impaired; no pediatric participant had hepatic impairment",
        "(Table 1)."
      ),
      source_name        = "Moderate hepatic impairment"
    ),
    HEPIMP_SEV = list(
      description        = "Severe hepatic impairment (Child-Pugh C)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function, mild impairment, or not listed)",
      notes              = paste(
        "Multiplicative effect on apparent clearance: CL/F * (1 - 0.410 *",
        "HEPIMP_SEV) (Table 2), a 41.0% reduction. Six adults (1.4%) were",
        "severely impaired; no pediatric participant had hepatic",
        "impairment (Table 1). Severe hepatic impairment is one of only two",
        "covariates the upstream adult analysis judged clinically relevant."
      ),
      source_name        = "Severe hepatic impairment"
    ),
    CONMED_ITRACONAZOLE = list(
      description        = "Concomitant itraconazole (strong CYP3A4 / P-gp inhibitor)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant azole)",
      notes              = paste(
        "Two multiplicative effects, both in Table 2: CL/F * (1 - 0.744 *",
        "CONMED_ITRACONAZOLE), a 74.4% clearance reduction giving a 3.9-fold",
        "AUC increase, and ktr * (1 - 0.361 * CONMED_ITRACONAZOLE). Studied",
        "in a dedicated adult drug-drug-interaction arm (22 subjects",
        "upstream); no pediatric participant took itraconazole. This is the",
        "other covariate the upstream adult analysis judged clinically",
        "relevant."
      ),
      source_name        = "Itraconazole use"
    ),
    CONMED_FLUCONAZOLE = list(
      description        = "Concomitant fluconazole (moderate CYP3A4 inhibitor)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant azole)",
      notes              = paste(
        "Multiplicative effect on apparent clearance only: CL/F * (1 -",
        "0.429 * CONMED_FLUCONAZOLE) (Table 2), a 42.9% reduction giving a",
        "1.75-fold AUC increase. Studied in a dedicated adult",
        "drug-drug-interaction arm (23 subjects upstream); no pediatric",
        "participant took fluconazole. Unlike itraconazole, fluconazole has",
        "no effect on the transit rate constant in this model."
      ),
      source_name        = "Fluconazole use"
    ),
    DOSE_RIMEGEPANT_MG = list(
      description        = "Administered rimegepant dose level",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-dose-record dose amount in mg, driving the nonlinear relative",
        "bioavailability F1 = (DOSE_RIMEGEPANT_MG / 10)^0.192 (Table 2",
        "'Dose effect on F 1'; the reference dose of 10 mg and the power",
        "form are printed on the upstream structural schematic as",
        "'F1 = (Dose/10)^0.191'). Rimegepant shows greater-than-",
        "dose-proportional exposure over 10 to 150 mg, attributed to",
        "dose-dependent autoinhibition of CYP3A-mediated first-pass",
        "metabolism. Must be supplied as a data column separate from the",
        "event-table amt: a covariate column literally named DOSE (any",
        "casing) is consumed by rxode2's etTrans() and never reaches",
        "model()."
      ),
      source_name        = "Dose"
    ),
    DOSE_LOW = list(
      description        = "Low-dose cohort indicator (10 mg or 25 mg rimegepant)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (75 mg or 150 mg)",
      notes              = paste(
        "Per-dose-record indicator, 1 for the 10 mg and 25 mg dose levels",
        "only. Multiplicative effect on the transit rate constant: ktr *",
        "(1 + 0.596 * DOSE_LOW) (Table 2 '10/25mg dose effect on k tr',",
        "whose footnote reads 'Increase in the transit rate constant for",
        "10-25 mg doses', confirming the positive sign is an increase). A",
        "step indicator rather than a threshold on DOSE_RIMEGEPANT_MG,",
        "because the paper defines the covariate by the two studied dose",
        "levels, not by an inequality; the 6 pediatric participants who",
        "received 25 mg carry DOSE_LOW = 1."
      ),
      source_name        = "10/25 mg dose"
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "rimegepant", units = "mg", specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "rimegepant", units = "mg", specimen = "administration site", verified = TRUE),
    transit2    = list(analyte = "rimegepant", units = "mg", specimen = "administration site", verified = TRUE),
    transit3    = list(analyte = "rimegepant", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "rimegepant", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "rimegepant", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 443,
    n_studies      = 14,
    age_range      = "6-77 years",
    age_median     = "43 years",
    weight_range   = "23.2-134 kg",
    weight_median  = "73.0 kg",
    sex_female_pct = 31.4,
    race_ethnicity = c(White = 83.5, Black = 9.3, Asian = 7.2),
    disease_state  = paste(
      "Healthy adults plus adults with renal or hepatic impairment;",
      "children 6 to <12 years with a history of migraine"
    ),
    dose_range     = "10-150 mg single, daily, or every-other-day oral tablet, capsule, or ODT",
    hepatic_function = "405 normal or not listed, 6 mild, 6 moderate, 6 severe (all adults)",
    co_medication  = "45 participants took concomitant itraconazole or fluconazole in dedicated DDI studies (all adults)",
    notes          = paste(
      "Baseline characteristics from Comisar 2025 Table 1. 423 adults from",
      "13 phase 1 studies (11 in the upstream adult model plus BHV3000-113",
      "and BHV3000-120) contributed 14,063 observations, and 20 children",
      "aged 6 to <12 years with a history of migraine from study C4951008",
      "contributed 78 observations. Pediatric dosing was weight-banded:",
      "25 mg ODT for 15-30 kg (n = 6), 50 mg ODT for >30-50 kg (n = 9), and",
      "75 mg ODT for >50 kg (n = 5); pediatric median age 9 years (range",
      "6-11) and median weight 35.9 kg (range 23.2-61.8). No adolescent",
      "(12 to <18 years) PK data were collected; that gap was bridged by",
      "simulation."
    )
  )

  ini({
    # ---- Structural parameters: typical values at the 70 kg reference weight,
    #      10 mg fasted oral-tablet reference for the absorption covariates
    #      (Comisar 2025 Table 2, "Estimate (%RSE)" column).
    lcl  <- log(25.2); label("Apparent clearance CL/F at 70 kg (L/h)")                       # Table 2: CL/F 25.2 L/h (4.8% RSE); Table S3 95% CI 22.8-27.6
    lvc  <- log(113);  label("Apparent central volume Vc/F at 70 kg (L)")                    # Table 2: V1/F 113 L (5.2% RSE); Table S3 95% CI 101-125
    lvp  <- log(46.8); label("Apparent peripheral volume Vp/F at 70 kg (L)")                 # Table 2: V2/F 46.8 L (4.9% RSE); Table S3 95% CI 42.3-51.3
    lq   <- log(4.16); label("Apparent inter-compartmental clearance Q/F at 70 kg (L/h)")    # Table 2: Q/F 4.16 L/h (5.7% RSE); Table S3 95% CI 3.70-4.62
    lktr <- log(8.42); label("Transit absorption rate constant ktr, tablet fasted (1/h)")    # Table 2: k tr 8.42 1/h (5.1% RSE); Table S3 95% CI 7.58-9.26
    lka  <- log(3.05); label("First-order rate constant ka from transit3 to central (1/h)")  # Table 2: k a 3.05 1/h (15.5% RSE); Table S3 95% CI 2.12-3.98

    # F1 is a purely relative bioavailability anchored at 1 for the 10 mg
    # fasted reference dose; the entire dose-dependence is carried by
    # e_dose_rimegepant_mg_fdepot below. The upstream schematic prints the
    # whole model for F1 as "F1 = (Dose/10)^0.191" with no separate scale
    # term, so there is no estimated typical value to place here.
    lfdepot <- fixed(log(1)); label("Relative bioavailability F1 at the 10 mg fasted reference dose (fraction)")  # Comisar 2025 CPT:PSP Figure 2 schematic

    # ---- Allometric body weight exponents. Estimated in this combined
    #      pediatric+adult model (the upstream adult-only model fixed them at
    #      0.75 and 1); a single exponent is shared by the two clearances and
    #      a second single exponent by the two volumes.
    e_wt_cl_q  <- 0.575; label("Allometric exponent on WT/70 for CL/F and Q/F (unitless)")     # Table 2: Body weight effect on CL/F and Q/F 0.575 (14.4% RSE); 95% CI 0.413-0.737
    e_wt_vc_vp <- 1.18;  label("Allometric exponent on WT/70 for Vc/F and Vp/F (unitless)")    # Table 2: Body weight effect on V 1 /F and V 2 /F 1.18 (8.3% RSE); 95% CI 0.988-1.372

    # ---- Covariate effects on apparent clearance. Categorical covariates use
    #      the upstream supplement's form theta_TV,REF * (1 + theta_x)^X with
    #      X in {0, 1}, i.e. a fractional change: -0.744 is a 74.4% reduction.
    e_hepimp_mod_cl          <- -0.203; label("Fractional change in CL/F with moderate hepatic impairment (unitless)")  # Table 2: Moderate hepatic impairment on CL/F -0.203 (27.6% RSE)
    e_hepimp_sev_cl          <- -0.410; label("Fractional change in CL/F with severe hepatic impairment (unitless)")    # Table 2: Severe hepatic impairment on CL/F -0.410 (25.9% RSE)
    e_conmed_fluconazole_cl  <- -0.429; label("Fractional change in CL/F with concomitant fluconazole (unitless)")      # Table 2: Fluconazole use on CL/F -0.429 (2.9% RSE)
    e_conmed_itraconazole_cl <- -0.744; label("Fractional change in CL/F with concomitant itraconazole (unitless)")     # Table 2: Itraconazole use on CL/F -0.744 (1.3% RSE)

    # ---- Covariate effects on the transit absorption rate constant.
    e_fed_ktr                 <- -0.698; label("Fractional change in ktr when dosed fed (unitless)")                 # Table 2: Fed on k tr -0.698 (3.9% RSE)
    e_form_odt_ktr            <-  0.470; label("Fractional change in ktr for the ODT vs tablet (unitless)")          # Table 2: ODT on k tr 0.470 (18.4% RSE)
    e_form_capsule_ktr        <-  2.19;  label("Fractional change in ktr for the capsule vs tablet (unitless)")      # Table 2: Capsule formulation on k tr 2.19 (22.6% RSE)
    e_conmed_itraconazole_ktr <- -0.361; label("Fractional change in ktr with concomitant itraconazole (unitless)")  # Table 2: Itraconazole use on k tr -0.361 (28.3% RSE)
    e_dose_low_ktr            <-  0.596; label("Fractional change in ktr at the 10 mg and 25 mg dose levels (unitless)")  # Table 2: 10/25mg dose effect on k tr 0.596 (45.5% RSE); footnote a "Increase in the transit rate constant for 10-25 mg doses"

    # ---- Covariate effects on relative bioavailability. Dose is a continuous
    #      covariate and takes the supplement's power form (dose / 10 mg)^theta;
    #      fed status is categorical and takes the (1 + theta) form.
    e_dose_rimegepant_mg_fdepot <-  0.192; label("Power exponent on (dose/10 mg) for F1 (unitless)")  # Table 2: Dose effect on F 1 0.192 (12.4% RSE); form from CPT:PSP Figure 2 "F1 = (Dose/10)^0.191"
    e_fed_fdepot                <- -0.331; label("Fractional change in F1 when dosed fed (unitless)")  # Table 2: Fed on F 1 -0.331 (6.9% RSE)

    # ---- Inter-individual variability. Values are the variances and
    #      covariances tabulated in supplement Table S3 for the pooled (final)
    #      model, which reproduce Table 2's IIV% = 100 * omega and its
    #      percent-correlation entries (see the vignette Errata for the one
    #      3-part-per-1000 discrepancy on the CL:Vp correlation).
    etalcl + etalvc + etalvp ~ c(0.104,
                                 0.103, 0.168,
                                 0.0428, 0.0717, 0.0666)  # Table S3: omega2CL 0.104, omegaCL-Vc 0.103, omega2Vc 0.168, omegaCL-Vp 0.0428, omegaVc-Vp 0.0717, omega2Vp 0.0666
    etalktr ~ 0.300                                       # Table S3: omega2ktr 0.300 (Table 2 IIV 54.8%; sqrt(0.300) = 0.548)

    # ---- Residual unexplained variability, combined additive plus
    #      proportional on the linear concentration scale.
    addSd  <- 0.447; label("Additive residual error SD (ng/mL)")            # Table 2: Additive error (SD) 0.447 (3.2% RSE); Table S3 sigma2add 0.200, sqrt = 0.4472. Table 2 labels the unit "ng/L"; see vignette Errata -- the assay LLOQ is 0.5 ng/mL
    propSd <- 0.393; label("Proportional residual error SD (fraction)")     # Table 2: Proportional error 39.3 %CV (23.8% RSE); Table S3 sigma2prop 0.154, sqrt = 0.3924
  })

  model({
    # ---- Individual disposition parameters. Allometric scaling is centred at
    #      70 kg on all four parameters, with one shared exponent for the
    #      clearances and one for the volumes.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q *
      (1 + e_hepimp_mod_cl * HEPIMP_MOD) *
      (1 + e_hepimp_sev_cl * HEPIMP_SEV) *
      (1 + e_conmed_fluconazole_cl * CONMED_FLUCONAZOLE) *
      (1 + e_conmed_itraconazole_cl * CONMED_ITRACONAZOLE)
    q  <- exp(lq) * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp

    # ---- Individual absorption parameters. Formulation, food, itraconazole,
    #      and the low-dose step all act multiplicatively on ktr; ka has no
    #      covariates and no IIV.
    ka <- exp(lka)
    ktr <- exp(lktr + etalktr) *
      (1 + e_fed_ktr * FED) *
      (1 + e_form_odt_ktr * FORM_ODT) *
      (1 + e_form_capsule_ktr * FORM_CAPSULE) *
      (1 + e_conmed_itraconazole_ktr * CONMED_ITRACONAZOLE) *
      (1 + e_dose_low_ktr * DOSE_LOW)

    # ---- Relative bioavailability: nonlinear in dose (power form, 10 mg
    #      reference) and reduced when dosed fed.
    fdepot <- exp(lfdepot) *
      (DOSE_RIMEGEPANT_MG / 10)^e_dose_rimegepant_mg_fdepot *
      (1 + e_fed_fdepot * FED)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- ODE system. The paper's "four transit absorption compartments"
    #      are depot + transit1 + transit2 + transit3: the dose lands in the
    #      first of the four, three sequential transfers run at the shared
    #      rate ktr, and the fourth empties into the central compartment at
    #      the separate first-order rate ka (Comisar 2025 CPT:PSP Figure 2,
    #      which draws four stacked boxes joined by three ktr arrows with a
    #      single ka arrow leaving the stack).
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2 - ka  * transit3
    d/dt(central)     <-  ka  * transit3 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central  - k21 * peripheral1

    f(depot) <- fdepot

    # central is in mg and vc in L, so central/vc is mg/L = ug/mL; the paper
    # reports plasma rimegepant in ng/mL (assay LLOQ 0.5 ng/mL).
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
