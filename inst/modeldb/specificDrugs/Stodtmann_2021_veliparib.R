Stodtmann_2021_veliparib <- function() {
  description <- "One-compartment population PK model for the oral PARP inhibitor veliparib (ABT-888) in adults with ovarian cancer, breast cancer or other solid tumours, pooled from 9 phase 1/2/3 studies (meta-analysis of individual-level data, n = 1470). First-order absorption without lag, with multiplicative fed and fasting effects on ka relative to an unknown-prandial-state reference; linear apparent clearance CL/F scaled by Cockcroft-Gault creatinine clearance (capped at 120 mL/min), serum albumin, male sex and strong CYP2D6 inhibitor co-medication; apparent volume Vc/F scaled by body weight, serum albumin and male sex. Separate combined additive-plus-proportional residual errors for the absorption phase (time after dose up to 2.5 h) and the elimination phase."
  reference <- paste(
    "Stodtmann S., Nuthalapati S., Eckert D., Kasichayanula S., Joshi R.,",
    "Bach B. A., Mensing S., Menon R., Xiong H. (2021). A Population",
    "Pharmacokinetic Meta-Analysis of Veliparib, a PARP Inhibitor, Across",
    "Phase 1/2/3 Trials in Cancer Patients.",
    "The Journal of Clinical Pharmacology 61(9):1195-1205.",
    "doi:10.1002/jcph.1875.",
    sep = " "
  )
  vignette <- "Stodtmann_2021_veliparib"

  # Stodtmann 2021 reports every parameter on a DAY time base (Table 4:
  # 'CL/F (L/day)', 'ka (1/day)') and the concentrations in ug/mL (Table 4
  # additive-error rows, Figures 1-2 axes 'mcg/mL'). Dose in mg over volume in
  # L is mg/L == ug/mL, so Cc = central / vc needs no scaling factor and the
  # additive residual SDs below are carried verbatim.
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(analyte = "veliparib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "veliparib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance estimated by the Cockcroft-Gault formula, NOT body-surface-area normalised",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Raw (non-BSA-normalised) Cockcroft-Gault creatinine clearance in mL/min",
        "(Stodtmann 2021 Table 2 footnote a: 'CrCL (based on Cockcroft-Gault",
        "formula) was tested both unrestricted and capped at 120 mL/min ... and",
        "the more significant improvement was taken forward'). The capped form was",
        "retained ('creatinine clearance (CrCL, capped at 120 mL/min)', Results,",
        "Significant Covariates) and the final equation prints min(CrCL, 120)/120,",
        "so the cap is applied INSIDE model() and a raw CRCL may be supplied.",
        "Reference 120 mL/min, which is also the cap, so the typical CL/F of 479",
        "L/day is that of a subject with CrCL >= 120 mL/min. Observed (Table 3,",
        "All Subjects): mean 102 mL/min, median 96.5, range 28.2-289. Baseline",
        "value (Table 3 is 'Patient Demographics and Baseline Factors'). Same raw",
        "Cockcroft-Gault usage of CRCL as the sibling Niu_2017_veliparib.R."
      ),
      source_name = "CrCL"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power term centred at the population median of 40 g/L on both CL/F",
        "(exponent 0.427) and Vc/F (exponent 0.260) (Stodtmann 2021 final",
        "equations and Table 4). Observed (Table 3, All Subjects): mean 39.6 g/L,",
        "median 40, range 20-52. Baseline value."
      ),
      source_name = "ALB"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power term centred at 70 kg on Vc/F only (exponent 0.505; Table 2",
        "reference value '70 kg' and the final Vc/F equation). Observed (Table 3,",
        "All Subjects): mean 69.6 kg, median 66.0, range 35.7-182. Baseline value."
      ),
      source_name = "WTKG"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (female)",
      notes = paste(
        "Stodtmann 2021 codes a MALE indicator (Table 2: 'Sex (male vs female)',",
        "reference Female; final equations '1.20^Male' on CL/F and '1.25^Male' on",
        "Vc/F). The canonical SEXF is the complement, so the model uses",
        "Male = 1 - SEXF; female (SEXF = 1) is the reference for which the",
        "typical values hold. The cohort was 97% female (1425 of 1470; Table 3)."
      ),
      source_name = "Male"
    ),
    CONMED_CYP2D6_INH = list(
      description = "Concomitant STRONG CYP2D6 inhibitor indicator, 1 = co-administered, 0 = not",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant strong CYP2D6 inhibitor)",
      notes = paste(
        "Only STRONG CYP2D6 inhibitors enter the 1 category (Table 2 comedication",
        "row 'strong inhibitors of CYP2D6'; final CL/F equation '0.885^CYP2D6',",
        "where CYP2D6 'is strong CYP2D6 inhibitor comedication'). Moderate and weak",
        "inhibitors are in the 0 category (Figure 3 legend: reference 'CYP2D6",
        "inhibitors other than strong'). The paper does not list which agents it",
        "classified as strong. Table 3 footnote a: 'Stated yes if at least 1",
        "observation occurred during comedication' -- a per-subject flag; 70 of 1470",
        "subjects (5%) were 'Yes'."
      ),
      source_name = "CYP2D6"
    ),
    FED = list(
      description = "Fed state at dosing, 1 = fed, 0 = fasting (read only when FED_MISSING = 0)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasting)",
      notes = paste(
        "Stodtmann 2021 enters the meal prior to the dose as a THREE-level factor",
        "'fasting vs fed vs unknown [reference]' on ka (Methods, Base Model). The",
        "printed typical ka of 59.4 1/day is therefore the UNKNOWN-prandial-state",
        "value; fasting multiplies it by 1.11 and fed by 0.356. Encoded with the",
        "FED_MISSING canonical for the unknown level: fed = FED * (1 - FED_MISSING),",
        "fasting = (1 - FED) * (1 - FED_MISSING). Only the food-effect subjects of",
        "phase 1 studies 3 and 4 had a recorded prandial state (fasting 27, fed 72;",
        "Table 3); the other 1396 were unknown."
      ),
      source_name = "Meal prior to dose"
    ),
    FED_MISSING = list(
      description = "Prandial state at dosing not recorded, 1 = unknown, 0 = recorded",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (prandial state recorded)",
      notes = paste(
        "The paper's REFERENCE level for the meal effect on ka ('unknown",
        "[reference]'), assigned to all subjects outside the phase 1 food-effect",
        "studies (1396 of 1470; Table 3). Set FED_MISSING = 1 to reproduce the",
        "published typical ka of 59.4 1/day and the fitted phase 2/3 population;",
        "set FED_MISSING = 0 with FED = 0 or 1 to simulate a defined fasting or fed",
        "dose."
      ),
      source_name = "Meal prior to dose = unknown"
    )
  )

  # Covariates Stodtmann 2021 SCREENED (Table 2) but did NOT retain.
  # Results / Discussion: "Race, age, region, cancer type, and concomitant use
  # of strong inhibitors of CYP3A4 and CYP2C19, strong inducers of CYP2C19 and
  # CYP3A4, and inhibitors of transporters (P-gp, [MATE]1/2, OCT2) were not
  # found to significantly impact veliparib pharmacokinetic parameters." AST,
  # ALT, total bilirubin and lean body weight were also screened (Table 2) and
  # are absent from the final equations. The CYP2C19 inhibitor / inducer and
  # the OCT2 and MATE1/2K inhibitor indicators have no register canonical and
  # are recorded in this comment only.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on CL/F and Vc/F (Table 2), not retained. Median 55 years, range 22-86 (Table 3)."
    ),
    RACE_BLACK = list(
      description = "Black race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened as 'Race (black vs other)' on CL/F and Vc/F (Table 2), not retained. 56 of 1470 subjects (4%) (Table 3)."
    ),
    REGION_JAPAN = list(
      description = "Japan region indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened as 'Region (Japan vs other)' on CL/F and Vc/F (Table 2), not retained. 77 of 1470 subjects (5%) (Table 3)."
    ),
    TUMTP_BREAST = list(
      description = "Breast cancer tumour-type indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened as part of 'Cancer type (breast vs ovarian vs other)' on CL/F and Vc/F (Table 2), not retained. 551 of 1470 subjects (38%) (Table 3)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened on CL/F (Table 2), not retained. Median 22 U/L, range 3-252 (Table 3)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened on CL/F (Table 2), not retained. Median 19 U/L, range 4-254 (Table 3)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened on CL/F (Table 2), not retained. Reported in mg/dL (not the register's umol/L): median 0.37, range 0.10-1.80 (Table 3)."
    ),
    LBM = list(
      description = "Lean body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened on Vc/F (Table 2), not retained (total body weight was). Median 41.3 kg, range 26.1-84.7 (Table 3)."
    ),
    CONMED_PGP_INH = list(
      description = "Concomitant P-glycoprotein inhibitor indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened on F, CL/F and Vc/F (Table 2), not retained. 25 of 1470 subjects (2%) (Table 3)."
    ),
    CONMED_CYP3A4_INH_STRONG = list(
      description = "Concomitant strong CYP3A4 inhibitor indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened on F, CL/F and Vc/F (Table 2), not retained. 11 of 1470 subjects (1%) (Table 3)."
    ),
    CONMED_CYP3A4_IND_STRONG = list(
      description = "Concomitant strong CYP3A4 inducer indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened on F, CL/F and Vc/F (Table 2), not retained. 6 of 1470 subjects (Table 3)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1470L,
    n_studies = 9L,
    age_range = "22-86 years",
    age_median = "55 years",
    weight_range = "35.7-182 kg",
    weight_median = "66.0 kg",
    sex_female_pct = 97,
    race_ethnicity = "White 75%, Asian 12%, Black 4%, Other 9% (Table 3)",
    disease_state = paste(
      "Adults with ovarian, fallopian tube or primary peritoneal cancer (58%), breast",
      "cancer (38%), or other advanced solid tumours (5%), receiving veliparib as",
      "monotherapy or with temozolomide, carboplatin/gemcitabine or",
      "carboplatin/paclitaxel."
    ),
    dose_range = paste(
      "Oral veliparib 10-400 mg twice daily, plus 40, 200 and 400 mg single doses in",
      "the phase 1 food-effect / crossover studies (Table 1). Phase 3 doses were",
      "120 mg BID (breast) and 150 or 300 mg BID (ovarian)."
    ),
    regions = "Multinational; Japan 5% (Table 3)",
    renal_function = paste(
      "Cockcroft-Gault CrCL median 96.5 mL/min (range 28.2-289). Normal 59%, mild",
      "impairment 32%, moderate 9%, severe 1 subject (Table 3)."
    ),
    hepatic_function = "Normal 76%, mild impairment 24% (Table 3)",
    notes = paste(
      "Individual-level meta-analysis of 9 AbbVie studies (6 phase 1, 1 phase 2, 2",
      "phase 3; Table 1). 9160 plasma concentrations after excluding 1.1% outliers",
      "(9262 with outliers); 2.9% below the LLOQ of about 1 ng/mL were imputed as",
      "LLOQ/2 (M5). Estimation by FOCE with interaction in NONMEM 7.4.3; 1000-sample",
      "bootstrap (all converged). Phase 2/3 sampling was sparse, giving eta shrinkage",
      "of 13% on CL/F and 28% on Vc/F (Discussion)."
    )
  )

  ini({
    # All values are FINAL estimates from Stodtmann 2021 Table 4 ('Population
    # Analysis, Estimate'), cross-checked against the printed final typical-
    # value equations (Results, Significant Covariates):
    #   CL/F = 479 * (min(CrCL,120)/120)^0.513 * (ALB/40)^0.427 * 1.20^Male * 0.885^CYP2D6  L/day
    #   Vc/F = 152 * (WTKG/70)^0.505 * (ALB/40)^0.260 * 1.25^Male  L
    #   ka   = 59.4 * 1.11^Fasting * 0.356^Fed  1/day
    lka <- log(59.4); label("First-order absorption rate constant ka, unknown prandial state (1/day)") # Table 4 'ka (1/day)' 59.4 (RSE 2.61%)
    lcl <- log(479); label("Apparent clearance CL/F, female, CrCL >= 120 mL/min, ALB 40 g/L, no strong CYP2D6 inhibitor (L/day)") # Table 4 'CL/F (L/day)' 479 (RSE 1.35%)
    lvc <- log(152); label("Apparent volume Vc/F, female, 70 kg, ALB 40 g/L (L)") # Table 4 'Vc/F (L)' 152 (RSE 1.10%)

    # Categorical effects are MULTIPLICATIVE FACTORS raised to the indicator
    # (theta^I), as printed in the final equations -- NOT the (1 + I*theta)
    # form of the Methods' generic categorical equation. The Results' effect
    # sizes confirm it: 1/1.20 - 1 = -16.7% AUC for males (paper -16.5%) and
    # 1/0.885 - 1 = +13.0% for strong CYP2D6 inhibitors (paper +13.0%); under
    # (1 + theta) the male effect would be -55%.
    e_fed_ka <- 0.356; label("Multiplicative factor on ka for a fed dose vs unknown prandial state (unitless)") # Table 4 'Fed on ka' 0.356 (RSE 3.93%)
    e_fast_ka <- 1.11; label("Multiplicative factor on ka for a fasting dose vs unknown prandial state (unitless)") # Table 4 'Fasting on ka' 1.11 (RSE 4.05%)
    e_crcl_cl <- 0.513; label("Power exponent on min(CRCL, 120)/120 for CL/F (unitless)") # Table 4 'Creatinine clearance on CL/F' 0.513 (RSE 5.98%)
    e_alb_cl <- 0.427; label("Power exponent on ALB/40 for CL/F (unitless)") # Table 4 'Albumin on CL/F' 0.427 (RSE 14.6%)
    e_male_cl <- 1.20; label("Multiplicative factor on CL/F for male sex (unitless)") # Table 4 'Male on CL/F' 1.20 (RSE 4.95%)
    e_cyp2d6inh_cl <- 0.885; label("Multiplicative factor on CL/F for strong CYP2D6 inhibitor co-medication (unitless)") # Table 4 'Strong inhibitors of CYP2D6 on CL/F' 0.885 (RSE 3.29%)
    e_wt_vc <- 0.505; label("Power exponent on WT/70 for Vc/F (unitless)") # Table 4 'Body weight on Vc/F' 0.505 (RSE 6.79%)
    e_alb_vc <- 0.260; label("Power exponent on ALB/40 for Vc/F (unitless)") # Table 4 'Albumin on Vc/F' 0.260 (RSE 23.9%)
    e_male_vc <- 1.25; label("Multiplicative factor on Vc/F for male sex (unitless)") # Table 4 'Male on Vc/F' 1.25 (RSE 5.44%)

    # Between-subject variability: variances (OMEGA). Table 4 footnote b:
    # '%CV is calculated as sqrt(exp[OMEGA(i,i)] - 1) x 100', and
    # sqrt(exp(0.085) - 1) = 29.8%, sqrt(exp(0.064) - 1) = 25.7% reproduce the
    # printed %CV, so the Estimate column is the variance. No BSV on ka (Results:
    # tested, not retained). No covariance reported; diagonal.
    etalcl ~ 0.085 # Table 4 'BSV on CL/F' 0.085 (29.8% CV)
    etalvc ~ 0.064 # Table 4 'BSV on Vc/F' 0.064 (25.7% CV)

    # Residual error: combined additive + proportional, with SEPARATE terms for
    # the absorption phase (before the 2.5 h Tmax) and the elimination phase
    # (Results: 'different error terms (proportional as well as additive) for
    # the absorption phase (before Tmax at 2.5 hours) and elimination phase
    # improved the OFV by 1128 points'). Read as standard deviations: the
    # additive rows carry concentration units 'ug/mL' (a variance would be
    # (ug/mL)^2), and the elimination-phase proportional RSE of 1.44% is below
    # the sqrt(2/9160) = 1.48% floor that a variance estimated from 9160
    # observations cannot beat.
    addSd_abs <- 0.004; label("Additive residual SD, absorption phase (ug/mL)") # Table 4 'Additive error in absorption phase (ug/mL)' 0.004 (RSE 4.78%)
    propSd_abs <- 0.208; label("Proportional residual SD, absorption phase (fraction)") # Table 4 'Proportional error in absorption phase' 0.208 (RSE 2.95%)
    addSd_elim <- 2.96e-7; label("Additive residual SD, elimination phase (ug/mL)") # Table 4 'Additive error in elimination phase (ug/mL)' 2.96 x 10^-7 (RSE 36.0%)
    propSd_elim <- 0.078; label("Proportional residual SD, elimination phase (fraction)") # Table 4 'Proportional error in elimination phase' 0.078 (RSE 1.44%)
  })

  model({
    # Covariate indicators in the paper's own orientation.
    male <- 1 - SEXF
    fed <- FED * (1 - FED_MISSING)
    fasting <- (1 - FED) * (1 - FED_MISSING)

    # CrCL is capped at 120 mL/min, which is also the reference value.
    crcl_capped <- min(CRCL, 120)

    ka <- exp(lka) * e_fast_ka^fasting * e_fed_ka^fed
    cl <- exp(lcl + etalcl) *
      (crcl_capped / 120)^e_crcl_cl *
      (ALB / 40)^e_alb_cl *
      e_male_cl^male *
      e_cyp2d6inh_cl^CONMED_CYP2D6_INH
    vc <- exp(lvc + etalvc) *
      (WT / 70)^e_wt_vc *
      (ALB / 40)^e_alb_vc *
      e_male_vc^male

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc

    # Phase-specific residual error. tad() is the time since the most recent
    # dose, in days here; the absorption phase is up to 2.5 h after the dose.
    tad_day <- tad()
    abs_phase <- tad_day <= 2.5 / 24
    addSd <- addSd_abs * abs_phase + addSd_elim * (1 - abs_phase)
    propSd <- propSd_abs * abs_phase + propSd_elim * (1 - abs_phase)
    Cc ~ add(addSd) + prop(propSd)
  })
}
