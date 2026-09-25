Zhu_2018_asunaprevir <- function() {
  description <- paste0(
    "Two-compartment population PK model for asunaprevir (ASV, BMS-650032; a ",
    "pangenotypic hepatitis C virus NS3 protease inhibitor) in 1239 adults with ",
    "chronic HCV genotype 1 or 4 infection, pooled from 3 Phase II and 2 Phase III ",
    "studies in which ASV was given with daclatasvir (DUAL) or with daclatasvir plus ",
    "peginterferon/ribavirin (QUAD). Absorption is sequential: a zero-order release ",
    "from the formulation over a duration D1 into the depot, followed by first-order ",
    "absorption (ka) into the central compartment; elimination is first order from ",
    "central. Apparent clearance (CL/F) carries a step-function auto-induction effect ",
    "that raises CL/F by 43% (50.8 to 72.5 L/h) from 48 h after the first dose onward, ",
    "attributed to CYP3A4 auto-induction; the induction dynamics themselves were not ",
    "estimable because the pooled dataset held few samples in the first 7 days. ",
    "CL/F additionally depends on age, sex, race (Black, Asian and Other relative to a ",
    "White reference), baseline AST, the on-treatment AST/baseline-AST ratio (a ",
    "time-varying hepatic-recovery term) and cirrhosis; Vc/F on sex and cirrhosis; ",
    "Vp/F on body weight; and ka, the zero-order duration D1 and relative ",
    "bioavailability on the formulation (Phase II tablet versus the Phase III / ",
    "commercial soft-gel capsule reference), with an additional bioavailability ",
    "increment at the 600 mg dose level that captures the observed more-than-dose-",
    "proportional exposure. Inter-individual variability is diagonal on CL/F, Vc/F, ",
    "Vp/F and ka, and residual variability is additive on log-transformed ",
    "concentrations (i.e. log-normal)."
  )

  reference <- paste(
    "Zhu L, Li H, Chan P, Eley T, Gandhi Y, Bifano M, Osawa M, Ueno T,",
    "Hughes E, AbuTarif M, Bertz R, Garimella T. (2018).",
    "Population Pharmacokinetic Analysis of Asunaprevir in Subjects with",
    "Hepatitis C Virus Infection. Infectious Diseases and Therapy 7(2):261-275.",
    "doi:10.1007/s40121-018-0197-y.",
    sep = " "
  )
  vignette <- "Zhu_2018_asunaprevir"

  # Doses are in mg and volumes in L, so `central / vc` is mg/L = ug/mL. Zhu 2018
  # reports exposures in ng.h/mL (Table 3), which is 1000x the model's native
  # ug.h/mL -- the vignette applies that factor at the comparison, not here.
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL/F normalised to the 55-year reference subject",
        "(Zhu 2018 Table 2 footnote a). Cohort median 57 years, range 18-79",
        "(Table 1)."
      ),
      source_name = "Age"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on Vp/F only, normalised to 70 kg. Weight was screened on",
        "CL/F in the univariate step but eliminated during backward elimination",
        "(Zhu 2018 Results, Covariate Models). Cohort median 70 kg, range 36-124",
        "(Table 1)."
      ),
      source_name = "Weight"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male; the reference subject of Zhu 2018 Table 2 footnote a)",
      notes = paste(
        "Zhu 2018 writes the effect as exp(theta * Female), so the source",
        "indicator is already female-coded and needs no transformation. Acts on",
        "both CL/F and Vc/F. Cohort 50.8% female (Table 1)."
      ),
      source_name = "Female"
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White / Caucasian is the reference race for the CL/F race effects)",
      notes = paste(
        "One of three mutually exclusive non-White indicators on CL/F",
        "(RACE_BLACK, RACE_ASIAN, RACE_OTHER); all three are 0 for the White",
        "reference. Cohort 6.0% Black (Table 1). The effect is small and its",
        "bootstrap 95% CI spans 0 (Table 2)."
      ),
      source_name = "Black"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White / Caucasian)",
      notes = paste(
        "Cohort 34.2% Asian, of which 21.5 percentage points were Japanese",
        "(Table 1). Zhu 2018 did not fit a separate Japanese indicator; the",
        "Japanese-versus-other-Asian exposure difference reported in Table 3",
        "arises from the age and sex distribution of the Japanese subset, not",
        "from a separate race parameter."
      ),
      source_name = "Asian"
    ),
    RACE_OTHER = list(
      description = "Race-category 'Other' indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White / Caucasian)",
      notes = paste(
        "Cohort 1.6% Other (Table 1). The effect's bootstrap 95% CI spans 0",
        "(Table 2), consistent with the small subgroup."
      ),
      source_name = "Otherrace"
    ),
    AST_BL = list(
      description = "Baseline (pre-treatment) serum aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed per subject; the anchor of the time-varying AST ratio term.",
        "Power effect on CL/F normalised to the 60 U/L reference subject (Zhu",
        "2018 Table 2 footnote a). Cohort median 51 U/L, range 13-595 (Table 1).",
        "Baseline ALT was also significant univariately but was dropped because",
        "it correlated with baseline AST at r^2 = 0.86 (Methods)."
      ),
      source_name = "Baseline AST"
    ),
    AST = list(
      description = "On-treatment serum aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "TIME-VARYING. Enters only as the ratio AST / AST_BL, so it equals 1 and",
        "contributes nothing at baseline. Captures the recovery of hepatic",
        "function during antiviral treatment: cohort mean AST fell by roughly 50%",
        "over the first 6 weeks (Results), and CL/F rises as AST falls. A dataset",
        "that carries no on-treatment AST measurements should set AST = AST_BL on",
        "every record, which reduces this term to 1."
      ),
      source_name = "AST"
    ),
    DIS_CIRRHOSIS = list(
      description = "Cirrhosis indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no cirrhosis; the reference subject of Zhu 2018 Table 2 footnote a)",
      notes = paste(
        "Acts on both CL/F and Vc/F. Only COMPENSATED cirrhosis with hepatic",
        "function no worse than Child-Pugh A was enrolled (Discussion), so the",
        "coefficient must not be extrapolated to Child-Pugh B or C. Cohort 19.5%",
        "cirrhotic, 80.4% not, 0.1% missing (Table 1)."
      ),
      source_name = "Cirrhosis"
    ),
    FORM_TABLET = list(
      description = "Phase II tablet formulation indicator",
      units = "(binary)",
      type = "binary",
      reference_category = paste(
        "0 (the Phase III / commercial soft-gel capsule, which anchors the",
        "relative bioavailability at 1)"
      ),
      notes = paste(
        "Per-dose-record. NOTE the comparator here is a soft-gel CAPSULE, not the",
        "non-tablet oral liquid named as this canonical's default reference",
        "category -- the same situation as Wada 2023 sparsentan. Acts on three",
        "absorption quantities at once: ka, the zero-order duration D1, and",
        "relative bioavailability. Cohort 30.3% tablet, 69.7% soft-gel (Table 1)."
      ),
      source_name = "Tablet"
    ),
    DOSE_600MG = list(
      description = "600 mg dose-level indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the 100 mg and 200 mg dose levels)",
      notes = paste(
        "Per-dose-record indicator that raises relative bioavailability at the",
        "600 mg level, encoding the more-than-dose-proportional exposure Zhu 2018",
        "observed on going from 200 mg BID to 600 mg BID of the Phase II tablet",
        "(Discussion). In the pooled dataset 600 mg was administered only as the",
        "tablet (7.1% of subjects, Table 1), so in practice DOSE_600MG = 1 implies",
        "FORM_TABLET = 1."
      ),
      source_name = "Dose(600 mg)"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "asunaprevir",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "asunaprevir",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "asunaprevir",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1239,
    n_studies = 5,
    age_range = "18-79 years",
    age_median = "57 years",
    weight_range = "36-124 kg",
    weight_median = "70 kg",
    sex_female_pct = 50.8,
    race_ethnicity = c(White = 58.3, Black = 6.0, Asian = 34.2, Other = 1.6),
    disease_state = paste(
      "Chronic hepatitis C virus genotype 1 (98.5%) or genotype 4 (1.5%)",
      "infection; 19.5% with compensated cirrhosis; mild hepatic impairment",
      "(Child-Pugh A) permitted; treatment-naive (33.5%), non-responder or",
      "null/partial responder (35.1%), or peginterferon/ribavirin",
      "ineligible/intolerant (31.4%)"
    ),
    dose_range = paste(
      "100, 200 or 600 mg once daily (4.3% of subjects) or twice daily (95.7%),",
      "as the Phase II tablet (30.3%) or the Phase III / commercial soft-gel",
      "capsule (69.7%); 69.7% received 100 mg BID soft-gel"
    ),
    regions = "Global, including a Japanese Phase II and a Japanese Phase III study",
    hepatic_function = paste(
      "Baseline AST median 51 U/L (range 13-595); baseline ALT median 60 U/L",
      "(range 7-475); both roughly twice the upper limit of normal on average and",
      "falling by about 50% over the first 6 weeks of treatment"
    ),
    renal_function = "Creatinine clearance median 101 mL/min (range 40-286)",
    co_medication = paste(
      "Daclatasvir in all subjects; peginterferon alfa and/or ribavirin in the",
      "QUAD and triple regimens"
    ),
    n_observations = 9496,
    notes = paste(
      "Demographics from Zhu 2018 Table 1 (n = 1239). The final PK analysis",
      "dataset held 9496 concentration records from 1236 subjects (Results);",
      "most subjects contributed 4-9 samples. The five contributing studies are",
      "listed in Supplemental Table 1 of the Electronic Supplementary Material."
    )
  )

  ini({
    # --- Structural parameters (Zhu 2018 Table 2, Fixed effects) ---------------
    # Reference subject: non-cirrhotic, male, 70 kg, 55 years old, White, baseline
    # AST 60 IU/L, receiving the soft-gel capsule, prior to induction
    # (Table 2 footnote a).
    lcl <- log(50.8); label("Apparent clearance CL/F, pre-induction (L/h)") # Zhu 2018 Table 2, 'CL/F (L/h)' = 50.8
    lvc <- log(47.6); label("Apparent central volume Vc/F (L)") # Zhu 2018 Table 2, 'Vc/F (L)' = 47.6
    lq <- log(21.6); label("Apparent intercompartmental clearance Q/F (L/h)") # Zhu 2018 Table 2, 'Q/F (L/h)' = 21.6
    lvp <- log(561); label("Apparent peripheral volume Vp/F (L)") # Zhu 2018 Table 2, 'Vp/F (L)' = 561
    lka <- log(0.484); label("First-order absorption rate constant ka (1/h)") # Zhu 2018 Table 2, 'Ka (1/h)' = 0.484
    ld1 <- log(1.12); label("Duration of the zero-order release into the depot, D1 (h)") # Zhu 2018 Table 2, 'D1 (h)' = 1.12

    # Relative bioavailability is referenced to the Phase III / commercial
    # soft-gel capsule at the 100 mg and 200 mg dose levels, for which Zhu 2018
    # estimated no F term -- every disposition parameter above is apparent
    # (/F) against that reference, so F is structurally anchored at 1 there.
    lfdepot <- fixed(log(1)); label("Relative bioavailability of the soft-gel capsule reference (unitless)") # Zhu 2018 Eq. for F: F = exp(-0.215*Tablet) * exp(0.65*Dose(600 mg)); the reference (Tablet = 0, 600 mg = 0) gives F = 1

    # --- Auto-induction of CL/F -----------------------------------------------
    # Modelled as a step function switching on 2 days (48 h) after the first
    # dose. exp(0.355) = 1.426, i.e. CL/F rises from 50.8 to 72.5 L/h, the +43%
    # the Abstract and Discussion report. A 6-day change point was also tried and
    # gave a worse OFV (Results, Structural Model).
    e_induction_cl <- 0.355; label("Log-scale effect of auto-induction on CL/F from 48 h after the first dose (unitless)") # Zhu 2018 Table 2, 'CL/F * induction' = 0.355

    # --- Covariate effects on CL/F (Zhu 2018 Eq. for CL/F) --------------------
    e_age_cl <- -0.341; label("Power exponent of age on CL/F, referenced to 55 years (unitless)") # Zhu 2018 Table 2, 'CL * age' = -0.341
    e_sexf_cl <- -0.117; label("Log-scale effect of female sex on CL/F (unitless)") # Zhu 2018 Table 2, 'CL * female' = -0.117
    e_race_black_cl <- 0.0386; label("Log-scale effect of Black race on CL/F versus White (unitless)") # Zhu 2018 Table 2, 'CL * Black Race' = 0.0386
    e_race_asian_cl <- -0.255; label("Log-scale effect of Asian race on CL/F versus White (unitless)") # Zhu 2018 Table 2, 'CL * Asian Race' = -0.255
    e_race_other_cl <- -0.0678; label("Log-scale effect of Other race on CL/F versus White (unitless)") # Zhu 2018 Table 2, 'CL * Other Race' = -0.0678
    e_ast_bl_cl <- -0.458; label("Power exponent of baseline AST on CL/F, referenced to 60 U/L (unitless)") # Zhu 2018 Eq. for CL/F, '(Baseline AST/60)^-0.458'; Table 2 prints the rounded 'CL * baseline AST' = -0.46
    e_ast_cl <- -0.291; label("Power exponent of the on-treatment AST / baseline AST ratio on CL/F (unitless)") # Zhu 2018 Eq. for CL/F, '(AST/Baseline AST)^-0.291'; Table 2 prints the rounded 'CL * AST' = -0.29
    e_dis_cirrhosis_cl <- -0.378; label("Log-scale effect of cirrhosis on CL/F (unitless)") # Zhu 2018 Table 2, 'CL * cirrhosis' = -0.378

    # --- Covariate effects on Vc/F (Zhu 2018 Eq. for Vc/F) -------------------
    e_sexf_vc <- -0.608; label("Log-scale effect of female sex on Vc/F (unitless)") # Zhu 2018 Table 2, 'Vc * female' = -0.608
    e_dis_cirrhosis_vc <- -0.835; label("Log-scale effect of cirrhosis on Vc/F (unitless)") # Zhu 2018 Table 2, 'Vc * cirrhosis' = -0.835

    # --- Covariate effect on Vp/F (Zhu 2018 Eq. for Vp/F) --------------------
    e_wt_vp <- 1.42; label("Power exponent of body weight on Vp/F, referenced to 70 kg (unitless)") # Zhu 2018 Table 2, 'Vp * weight' = 1.42

    # --- Covariate effects on absorption (Zhu 2018 Eqs. for Ka, D and F) -----
    e_form_tablet_ka <- -0.503; label("Log-scale effect of the tablet formulation on ka (unitless)") # Zhu 2018 Table 2, 'Ka * tablet' = -0.503
    e_form_tablet_d1 <- 0.864; label("Log-scale effect of the tablet formulation on D1 (unitless)") # Zhu 2018 Table 2, 'D1 * tablet' = 0.864
    e_form_tablet_fdepot <- -0.215; label("Log-scale effect of the tablet formulation on relative bioavailability (unitless)") # Zhu 2018 Table 2, 'Relative F * tablet' = -0.215
    e_dose_600mg_fdepot <- 0.65; label("Log-scale effect of the 600 mg dose level on relative bioavailability (unitless)") # Zhu 2018 Table 2, 'Relative F1 * 600 mg' = 0.65

    # --- Inter-individual variability -----------------------------------------
    # Zhu 2018 fitted a DIAGONAL omega block on CL/F, Vc/F, Vp/F and Ka
    # (Results, Structural Model), so no off-diagonal terms are carried. Table 2
    # prints the variance first and its square root in parentheses; the variance
    # is the value used here. Shrinkage was 13% on CL/F but 39%, 65% and 42% on
    # Vc/F, Vp/F and Ka respectively, so only the CL/F random effect is
    # well-informed by the data.
    etalcl ~ 0.168 # Zhu 2018 Table 2, random effect 'CL/F' variance = 0.168 (SD 0.41)
    etalvc ~ 2.19 # Zhu 2018 Table 2, random effect 'Vc/F' variance = 2.19 (SD 1.48); 39% shrinkage
    etalvp ~ 0.777 # Zhu 2018 Table 2, random effect 'Vp/F' variance = 0.777 (SD 0.881); 65% shrinkage
    etalka ~ 0.300 # Zhu 2018 Table 2, random effect 'Ka' variance = 0.300 (SD 0.548); 42% shrinkage

    # --- Residual variability --------------------------------------------------
    # "Residual variability was modeled using an additive error model with
    # log-transformed ASV concentrations" (Results, Structural Model), which is
    # nlmixr2's lnorm() error model with the SD on the natural-log scale.
    expSd <- 0.621; label("Residual SD, additive on log-transformed concentrations (log units)") # Zhu 2018 Table 2, residual error 'e' variance = 0.386, square root = 0.621
  })

  model({
    # 1. Derived terms ---------------------------------------------------------
    # Step-function auto-induction of CL/F, switching on 48 h (2 days) after the
    # FIRST dose. `max(0, tafd())` is 0 before any dose has been given and so
    # does not propagate NA into cl the way a bare tafd() would.
    induction <- (max(0, tafd()) > 48)

    # 2. Individual parameters -------------------------------------------------
    # CL/F, Zhu 2018 Eq. for CL/F (p. 266).
    cl <- exp(lcl + etalcl) *
      exp(e_induction_cl * induction) *
      (AGE / 55)^e_age_cl *
      exp(e_sexf_cl * SEXF) *
      exp(e_race_black_cl * RACE_BLACK) *
      exp(e_race_asian_cl * RACE_ASIAN) *
      exp(e_race_other_cl * RACE_OTHER) *
      (AST_BL / 60)^e_ast_bl_cl *
      (AST / AST_BL)^e_ast_cl *
      exp(e_dis_cirrhosis_cl * DIS_CIRRHOSIS)

    # Vc/F, Zhu 2018 Eq. for Vc/F (p. 266).
    vc <- exp(lvc + etalvc) *
      exp(e_sexf_vc * SEXF) *
      exp(e_dis_cirrhosis_vc * DIS_CIRRHOSIS)

    # Vp/F, Zhu 2018 Eq. for Vp/F (p. 267).
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vp

    # Q/F carried no covariates in the final model (Table 2).
    q <- exp(lq)

    # ka and the zero-order release duration D1, Zhu 2018 Eqs. for Ka and D
    # (p. 267).
    ka <- exp(lka + etalka) * exp(e_form_tablet_ka * FORM_TABLET)
    d1 <- exp(ld1) * exp(e_form_tablet_d1 * FORM_TABLET)

    # Relative bioavailability, Zhu 2018 Eq. for F (p. 267).
    fdepot <- exp(lfdepot) *
      exp(e_form_tablet_fdepot * FORM_TABLET) *
      exp(e_dose_600mg_fdepot * DOSE_600MG)

    # 3. Micro-constants -------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system ------------------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Absorption input ------------------------------------------------------
    # Sequential zero-order release then first-order absorption: the dose enters
    # `depot` at a constant rate over D1 hours and then leaves it first-order at
    # ka. Dose records MUST carry rate = -2 for rxode2 to honour dur(depot);
    # without it the dose is delivered as a bolus and D1 is silently ignored.
    dur(depot) <- d1
    f(depot) <- fdepot

    # 6. Observation and error -------------------------------------------------
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
