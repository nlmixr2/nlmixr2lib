Ji_2022_pemigatinib <- function() {
  description <- paste(
    "Two-compartment population PK model for oral pemigatinib (a selective",
    "fibroblast growth factor receptor 1-3 inhibitor) in adults with",
    "advanced malignancies including cholangiocarcinoma (Ji 2022; N = 318",
    "from FIGHT-101, FIGHT-102 and FIGHT-202, 2968 concentration records).",
    "First-order absorption and linear elimination. Covariate effects:",
    "female sex lowers CL/F by 19% and raises ka 1.565-fold (male is the",
    "reference level); CL/F is 14.1% higher when no phosphate binder is",
    "coadministered (binder use is the reference level); Vc/F is 24.4%",
    "lower without a concomitant proton-pump inhibitor (PPI use is the",
    "reference level); body weight scales Vc/F and Vp/F as power functions",
    "of (WT / 73.3 kg). Correlated IIV on CL/F and Vc/F plus IIV on ka;",
    "residual error additive on the log scale (log-normal)."
  )
  reference <- paste(
    "Ji T, Chen X, Liu X, Yeleswaram S. Population Pharmacokinetics",
    "Analysis of Pemigatinib in Patients With Advanced Malignancies.",
    "Clin Pharmacol Drug Dev. 2022;11(4):454-466. doi:10.1002/cpdd.1038"
  )
  vignette <- "Ji_2022_pemigatinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline body weight (BWT in Ji 2022 Equations II and III).",
        "Power functions on Vc/F and Vp/F normalized to the cohort median",
        "73.3 kg (Table 2: median 73.3 kg, range 39.8-156.0 kg). No weight",
        "effect on CL/F, Q/F or ka."
      ),
      source_name = "BWT"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male) -- the paper's reference level",
      notes = paste(
        "1 = female, 0 = male. The source column SEXN is coded 1 = male,",
        "2 = female (Ji 2022 'Where:' block under Equations I-IV), so",
        "SEXN = 1 + SEXF and the paper's (1 - SEXN) term equals -SEXF.",
        "Equations I and IV therefore give female CL/F = 0.81 x male and",
        "female ka = 1.565 x male, which agrees with the Discussion",
        "('typical ka value is 56.5% higher for female patients'; 'typical",
        "CL/F value of female patients is 19% lower') and with Figure 2A",
        "(post hoc CL/F female vs male GMR 0.852). The Results/Abstract",
        "sentence attributing ka 1.49 1/h and CL/F 10.3 L/h to FEMALES",
        "contradicts the equations; those are the male, no-binder values."
      ),
      source_name = "SEXN"
    ),
    CONMED_PHOSBINDER = list(
      description = "Concomitant phosphate-binding agent use",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (phosphate binder used) -- the paper's reference level in Equation I",
      notes = paste(
        "BINDER in Ji 2022: 1 = used, 0 = not used, time-varying (defined",
        "for the jth subject on the ith visit). 117 of 318 patients (36.1%)",
        "took phosphate binders (Table 3), used to manage the on-target",
        "hyperphosphatemia of FGFR inhibition. Equation I applies",
        "(1 + 0.141 * (1 - BINDER)): CL/F is 9.00 L/h with a binder and",
        "1.141-fold higher without one."
      ),
      source_name = "BINDER"
    ),
    CONMED_PPI = list(
      description = "Concomitant proton-pump inhibitor use",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (PPI used) -- the paper's reference level in Equation II",
      notes = paste(
        "PPI in Ji 2022, time-varying (subscript ji in Equation II); coded",
        "1 = used, 0 = not used by analogy with BINDER (the 'Where:' block",
        "does not define PPI explicitly, but the Results text 'PPI",
        "coadministration increases typical Vc/F' fixes the direction).",
        "87 of 318 patients (26.9%) took a PPI (Table 3). Equation II",
        "applies (1 - 0.244 * (1 - PPI)): Vc/F is 161 L (at 73.3 kg) with a",
        "PPI and 122 L without one."
      ),
      source_name = "PPI"
    )
  )

  covariatesDataExcluded <- list(
    CONMED_CYP3A4_IND = list(
      description = "Concomitant weak (or moderate) CYP3A4 inducer use",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Retained on CL/F through forward and backward selection (+24.2% in",
        "the full multivariable model) but removed at model refinement",
        "because its 95% CI included zero (Ji 2022 Results 'Model",
        "Refinement' and Discussion). Not in the final model."
      )
    )
  )

  compartmentData <- list(
    depot = list(analyte = "pemigatinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pemigatinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "pemigatinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 318,
    n_studies = 3,
    n_observations = 2968,
    age_range = "21-79 years",
    age_median = "59 years",
    weight_range = "39.8-156.0 kg",
    weight_median = "73.3 kg",
    sex_female_pct = 55.0,
    race_ethnicity = c(White = 68.2, Japanese = 8.8, Asian = 7.2, Black = 6.3, Hispanic = 6.0, Other = 3.5),
    disease_state = paste(
      "Advanced malignancies (FIGHT-101, FIGHT-102 Japanese patients) and",
      "advanced/metastatic or surgically unresectable cholangiocarcinoma",
      "after at least one prior therapy (FIGHT-202)."
    ),
    dose_range = paste(
      "1-20 mg orally once daily (FIGHT-101 dose escalation), 9 and 13.5 mg",
      "(FIGHT-102), 13.5 mg (FIGHT-202); 2-weeks-on/1-week-off (84%) or",
      "continuous (16%) schedules."
    ),
    regions = "United States, Denmark, Japan, European Union, United Kingdom, Republic of Korea, Taiwan, Thailand",
    renal_function = "Normal 45.9%, mild impairment 42.1%, moderate impairment 11.9% (MDRD eGFR)",
    hepatic_function = "Normal 67.0%, mild impairment 29.6%, moderate impairment 3.5% (NCI criteria)",
    co_medication = paste(
      "Phosphate binders 36.1%, PPI 26.9%, H2-receptor antagonists 11.1%,",
      "diuretics 7.4%, weak/moderate CYP3A4 inhibitors 23.3%/4.1%, weak",
      "CYP3A4 inducers 9.1% (Table 3)."
    ),
    notes = "Baseline demographics from Ji 2022 Table 2; study designs Table 1; co-medications Table 3."
  )

  ini({
    # Structural parameters -- Ji 2022 Table 4 'Population Mean' and
    # Equations I-IV. Reference subject: male (SEXN = 1), on a phosphate
    # binder (BINDER = 1), on a PPI (PPI = 1), 73.3 kg.
    lka <- log(1.49)
    label("Absorption rate constant ka for a male (1/h)") # Table 4 ka 1.49 1/h; Equation IV
    lcl <- log(9.00)
    label("Apparent clearance CL/F for a male on a phosphate binder (L/h)") # Table 4 CL/F 9.00 L/h; Equation I
    lvc <- log(161)
    label("Apparent central volume Vc/F at 73.3 kg on a PPI (L)") # Table 4 Vc/F 161 L; Equation II
    lvp <- log(80.1)
    label("Apparent peripheral volume Vp/F at 73.3 kg (L)") # Table 4 Vp/F 80.1 L; Equation III
    lq <- log(16.0)
    label("Apparent intercompartmental clearance Q/F (L/h)") # Table 4 Q/F 16.0 L/h

    # Covariate effects
    e_phosbinder_cl <- 0.141
    label("Fractional CL/F increase when NOT on a phosphate binder (unitless)") # Table 4 'Phosphate binder on CL' 0.141; Equation I (1 + 0.141 x (1 - BINDER))
    e_sex_cl <- 0.190
    label("Sex coefficient on CL/F applied as (1 + e x (1 - SEXN)) (unitless)") # Table 4 'Sex (male vs female) on CL' 0.190; Equation I
    e_sex_ka <- 0.565
    label("Sex coefficient on ka applied as (1 - e x (1 - SEXN)) (unitless)") # Equation IV (1 - 0.565 x (1 - SEXN)); not tabulated in Table 4
    e_ppi_vc <- -0.244
    label("Fractional Vc/F change when NOT on a PPI (unitless)") # Table 4 'Proton pump inhibitor on Vc/F' -0.244; Equation II (1 - 0.244 x (1 - PPI))
    e_wt_vc <- 0.738
    label("Power exponent of (WT/73.3) on Vc/F (unitless)") # Table 4 'Body weight (median = 73.3 kg) on Vc/F' 0.738; Equation II
    e_wt_vp <- 1.22
    label("Power exponent of (WT/73.3) on Vp/F (unitless)") # Table 4 'Body weight on Vp/F' 1.22; Equation III

    # IIV -- Table 4 '%CV' is 100 x sqrt(omega^2) (the 95% CI of the ka
    # CV, 112-143, matches the %RSE of 12.3 on the variance scale only
    # under this reading). IIV on Vp/F and Q/F not estimated (NE).
    etalcl + etalvc ~ c(0.188356, 0.122, 0.123201) # Table 4: CL/F 43.4 %CV -> 0.434^2; 'Omega matrix for CL/F and Vc/F' covariance 0.122; Vc/F 35.1 %CV -> 0.351^2
    etalka ~ 1.6129 # Table 4: ka 127 %CV -> 1.27^2

    # Residual error: log-transformed data with an additive error (SD) on
    # the log scale, i.e. log-normal error on Cc.
    expSd <- 0.401
    label("Additive residual SD on the log scale (log ng/mL)") # Table 4 'RV,SD' 0.401
  })

  model({
    # 1. Covariate terms, written exactly as Ji 2022 Equations I-IV.
    #    SEXN is 1 = male, 2 = female; canonical SEXF is 1 = female.
    sexn <- 1 + SEXF

    # 2. Individual PK parameters
    ka <- exp(lka + etalka) * (1 - e_sex_ka * (1 - sexn))
    cl <- exp(lcl + etalcl) *
      (1 + e_phosbinder_cl * (1 - CONMED_PHOSBINDER)) *
      (1 + e_sex_cl * (1 - sexn))
    vc <- exp(lvc + etalvc) *
      (WT / 73.3)^e_wt_vc *
      (1 + e_ppi_vc * (1 - CONMED_PPI))
    vp <- exp(lvp) * (WT / 73.3)^e_wt_vp
    q <- exp(lq)

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Observation: dose mg / volume L = mg/L; x 1000 -> ng/mL (Figure 1
    #    axes are in ng/mL).
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd)
  })
}
