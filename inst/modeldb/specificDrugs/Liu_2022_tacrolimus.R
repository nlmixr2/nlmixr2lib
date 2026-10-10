Liu_2022_tacrolimus <- function() {
  description <- paste0(
    "One-compartment population pharmacokinetic model for intravenous ",
    "(continuous 24 h infusion) and oral tacrolimus in pediatric hematopoietic ",
    "stem cell transplant (HSCT) recipients (Liu 2022, final model 1, n = 86): ",
    "first-order absorption with ka fixed at 4.48 1/h and an estimated oral ",
    "bioavailability; power effects of body weight and hematocrit on CL, ",
    "exponential effects of concomitant azole antifungals, concomitant ",
    "caspofungin and post-transplant day >= 28 on CL, and a power effect of ",
    "hematocrit on V."
  )
  reference <- paste0(
    "Liu XL, Guan YP, Wang Y, Huang K, Jiang FL, Wang J, Yu QH, Qiu KF, ",
    "Huang M, Wu JY, Zhou DH, Zhong GP, Yu XX. Population Pharmacokinetics ",
    "and Initial Dosage Optimization of Tacrolimus in Pediatric Hematopoietic ",
    "Stem Cell Transplant Patients. Front Pharmacol. 2022;13:891648. ",
    "doi:10.3389/fphar.2022.891648"
  )
  vignette <- "Liu_2022_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Power effect on CL, (WT / 17.8)^0.56 (Liu 2022 Eq. 4). The centring ",
        "value 17.8 kg is the one printed in Eq. 4; Table 1 reports a cohort ",
        "median of 17.4 kg (the paper's 'typical patient' in Table 4 also ",
        "uses 17.4 kg). Treated as time-varying in the source TDM dataset."
      ),
      source_name = "WT"
    ),
    HCT = list(
      description = "Hematocrit",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Time-varying. Liu 2022 reports hematocrit as a volume fraction ",
        "(Table 2 median 0.289, range 0.166-0.445); the canonical HCT column ",
        "is percent, so the paper's centring values are rescaled by 100 ",
        "inside model(): CL uses (HCT / 29.6)^-0.66 (Eq. 4 prints ",
        "Hct/0.296) and V uses (HCT / 28.9)^-0.66 (Eq. 5 prints Hct/0.289). ",
        "Pass HCT in percent (e.g. 28.9), not as a fraction. The paper's ",
        "Discussion gives HCT = 0.0029 * Hgb + 0.0038 (Hgb in g/L) for ",
        "deriving hematocrit from hemoglobin."
      ),
      source_name = "Hct"
    ),
    CONMED_AZOLE = list(
      description = "Concomitant azole antifungal therapy indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant azole antifungal)",
      notes = paste0(
        "Time-varying. Liu 2022 pooled voriconazole, itraconazole and ",
        "posaconazole into a single variable 'CZ' because few patients ",
        "received any one agent (Methods, Covariate Analysis). Effect ",
        "exp(-0.40 * CONMED_AZOLE) on CL (Eq. 4, Table 3 theta CZ,CL)."
      ),
      source_name = "CZ"
    ),
    CONMED_CASPOFUNGIN = list(
      description = "Concomitant caspofungin therapy indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant caspofungin)",
      notes = paste0(
        "Time-varying. Effect exp(0.16 * CONMED_CASPOFUNGIN) on CL (Eq. 4, ",
        "Table 3 theta CPFG,CL)."
      ),
      source_name = "CPFG"
    ),
    POD = list(
      description = "Post-transplant day",
      units = "days",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Time-varying days since the stem cell transplant. Liu 2022 enters it ",
        "only as a two-level categorical variable PTD (PTD = 1 before day 28, ",
        "PTD = 2 from day 28 on; Methods and the footnote to Eqs. 4-7 ",
        "'PTD = 2 is the post-transplant days >= 28 days'), so the model ",
        "derives the indicator (POD >= 28) internally and applies ",
        "exp(-0.32 * indicator) to CL (Eq. 4)."
      ),
      source_name = "PTD"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 86L,
    n_studies = 1L,
    n_observations = 578L,
    age_range = "1-16 years",
    age_median = "5 years",
    weight_range = "6.00-50.0 kg",
    weight_median = "17.4 kg",
    sex_female_pct = 38.4,
    disease_state = paste0(
      "Pediatric allogeneic hematopoietic stem cell transplant recipients ",
      "(bone marrow, peripheral blood or cord blood grafts). Diagnoses: ",
      "beta-thalassemia 27.9%, acute lymphatic leukemia 25.6%, acute ",
      "myelogenous leukemia 19.8%, aplastic anemia 11.6%, chronic myeloid ",
      "leukemia 3.5%, juvenile myelomonocytic leukemia 2.3%, other 9.3%. ",
      "Acute GVHD in 70.9%, chronic GVHD in 10.5%; 80.2% unrelated donors."
    ),
    dose_range = paste0(
      "Tacrolimus continuous 24 h IV infusion, median 0.024 mg/kg/day ",
      "(0.004-0.056), followed by oral tacrolimus, median 0.057 mg/kg/day ",
      "(0.008-0.205)."
    ),
    regions = "China (Sun Yat-sen Memorial Hospital, Guangzhou)",
    hematocrit = "median 0.289 (range 0.166-0.445) as a fraction (Table 2)",
    hemoglobin = "median 96.5 g/L (range 54-150) (Table 2)",
    notes = paste0(
      "Retrospective therapeutic-drug-monitoring data, January 2017 to ",
      "December 2020: 578 whole-blood trough concentrations (320 during IV ",
      "and 258 during oral dosing; median 5, range 1-21 per patient) ",
      "measured by EMIT (Viva-E). Fitted with FOCE-I in Phoenix NLME 7.0. ",
      "Demographics in Liu 2022 Table 1, laboratory values in Table 2. ",
      "Race/ethnicity is not reported (single centre in Guangzhou, ",
      "China). The CYP3A5-genotyped subpopulation (n = 24) model ",
      "is Liu_2022_tacrolimus_cyp3a5."
    )
  )

  ini({
    # Structural parameters (Liu 2022 Table 3, 'Final Model (n = 86)' column)
    lka <- fixed(log(4.48)); label("Absorption rate constant (1/h)") # Table 3 ka = 4.48 h-1, no RSE; Results: 'Ka was fixed at 4.48 according to the literature (Jusko et al., 1995; Wallin et al., 2009)'
    lcl <- log(2.42); label("Clearance at WT 17.8 kg and hematocrit 29.6% (L/h)") # Table 3 CL = 2.42 L/h (RSE 10.84%); Eq. 4
    lvc <- log(79.6); label("Volume of distribution at hematocrit 28.9% (L)") # Table 3 V = 79.6 L (RSE 16.51%); Eq. 5
    lfdepot <- log(0.19); label("Oral bioavailability (fraction)") # Table 3 F = 0.19 (RSE 13.01%)

    # Covariate effects (Liu 2022 Table 3 and Eqs. 4-5)
    e_wt_cl <- 0.56; label("Power exponent of body weight on CL (unitless)") # Table 3 theta WT,CL = 0.56 (RSE 25.65%)
    e_hct_cl <- -0.66; label("Power exponent of hematocrit on CL (unitless)") # Table 3 theta Hct,CL = -0.66 (RSE 24.78%)
    e_conmed_azole_cl <- -0.40; label("Exponential effect of concomitant azole antifungal on CL (unitless)") # Table 3 theta CZ,CL = -0.40 (RSE 12.54%)
    e_conmed_caspofungin_cl <- 0.16; label("Exponential effect of concomitant caspofungin on CL (unitless)") # Table 3 theta CPFG,CL = 0.16 (RSE 46.65%)
    e_pod_cl <- -0.32; label("Exponential effect of post-transplant day >= 28 on CL (unitless)") # Table 3 theta PTD,CL = -0.32 (RSE 38.32%)
    e_hct_vc <- -0.66; label("Power exponent of hematocrit on V (unitless)") # Table 3 theta Hct,V = -0.66 (RSE 24.61%)

    # IIV: exponential (Methods Eq. 1); Table 3 reports omega^2 (variances)
    etalcl ~ 0.06 # Table 3 omega2 CL = 0.06 (RSE 22.56%)
    etalvc ~ 0.66 # Table 3 omega2 V = 0.66 (RSE 46.27%)
    etalfdepot ~ 0.51 # Table 3 omega2 F = 0.51 (RSE 30.64%)

    # Residual error
    propSd <- 0.374; label("Proportional residual error (fraction)") # Table 3 sigma proportional = 37.4% (RSE 3.73%)
  })

  model({
    # Post-transplant day indicator: PTD = 2 (>= 28 days) vs PTD = 1 (Eqs. 4-7 footnote)
    ptd_late <- POD >= 28

    # Individual parameters (Eqs. 4 and 5); HCT in percent, so the printed
    # fraction centring values 0.296 and 0.289 become 29.6 and 28.9.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      (WT / 17.8)^e_wt_cl *
      (HCT / 29.6)^e_hct_cl *
      exp(e_conmed_azole_cl * CONMED_AZOLE) *
      exp(e_conmed_caspofungin_cl * CONMED_CASPOFUNGIN) *
      exp(e_pod_cl * ptd_late)
    vc <- exp(lvc + etalvc) * (HCT / 28.9)^e_hct_vc
    fdepot <- exp(lfdepot + etalfdepot)

    kel <- cl / vc

    # Oral doses go to depot; the continuous 24 h IV infusion goes to central.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    f(depot) <- fdepot

    # Dose in mg and V in L give mg/L; x1000 gives the whole-blood
    # concentration in ng/mL (ug/L), the unit of the EMIT assay and Figure 1.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
