Schlender_2016_aging_physiology <- function() {
  description <- paste(
    "PBPK system-parameter (physiology) model of healthy human aging",
    "from 30 to 100 years (Schlender 2016 Clin Pharmacokinet). This is",
    "the SYSTEM layer of a whole-body PBPK approach, not a drug model:",
    "it has no drug, no dosing, no compartments and no ODEs, and every",
    "output is an algebraic function of the rxode2 time variable, which",
    "the model interprets as chronological age in YEARS. Outputs are",
    "body weight and height, thirteen age-varying organ / tissue masses",
    "plus a constant gastrointestinal mass, the cardiac index, ten organ",
    "blood flow rates, cardiac output, and the glomerular filtration",
    "rate. Anthropometry, masses and flows are the decade-knot values of",
    "Tables 2 and 3 with the linear interpolation between ten-year age",
    "bins that the paper specifies (Methods, Workflow for Elderly PBPK",
    "Model Development); outside 30 to 100 years the knot values are",
    "held flat. GFR is the paper's own contribution: a new reverse",
    "sigmoid hyperbolic AGING function (Equation 1) fitted with the",
    "Matlab lsqnonlin routine, giving specific GFR per 100 g of kidney,",
    "which this model multiplies by the age- and sex-specific kidney",
    "mass to return absolute GFR in mL/min. Sex is a required covariate",
    "because every trajectory in the paper is sex-specific.",
    "The drug-level PBPK models that Schlender 2016 uses to VERIFY this",
    "physiology (intravenous morphine and furosemide) are NOT included:",
    "they were built in PK-Sim v5.5, their tissue:plasma partition",
    "coefficients and cellular permeabilities are given only as the",
    "named Willmann et al. PK-Sim standard method rather than as values,",
    "and no project file was deposited, so the drug layer is not",
    "reproducible outside that platform. See the vignette for the full",
    "reasoning and for registry morphine / furosemide models that can be",
    "chained to this physiology.",
    sep = " "
  )
  reference <- paste(
    "Schlender JF, Meyer M, Thelen K, Krauss M, Willmann S, Eissing T,",
    "Jaehde U.",
    "Development of a Whole-Body Physiologically Based Pharmacokinetic",
    "Approach to Assess the Pharmacokinetics of Drugs in Elderly",
    "Individuals.",
    "Clin Pharmacokinet. 2016;55(12):1573-1589.",
    "doi:10.1007/s40262-016-0422-3.",
    sep = " "
  )
  vignette <- "Schlender_2016_aging_physiology"
  units <- list(
    time = "year (chronological age; valid 30 to 100)",
    dosing = "n/a (no exogenous dosing; PBPK system-physiology model)",
    concentration = paste(
      "n/a (no drug concentration). Outputs are body weight (kg),",
      "height (cm), organ masses (g), cardiac index (L/min/m^2),",
      "organ blood flows and cardiac output (mL/min), specific GFR",
      "(mL/min per 100 g kidney) and absolute GFR (mL/min)."
    )
  )

  covariateData <- list(
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Time-fixed and REQUIRED: Schlender 2016 reports every",
        "anthropometric, organ-mass and blood-flow trajectory as a",
        "separate female and male row (Tables 2 and 3), and the GFR",
        "aging half-time TA50 of Equation 1 is 59 years for females",
        "versus 54 years for males. The model selects the female branch",
        "when SEXF is 1 and the male branch when SEXF is 0; any other",
        "value linearly blends the two sexes and is not meaningful.",
        "The virtual populations used for the physiological-consistency",
        "checks were 5,000 individuals of each sex (Methods, Workflow",
        "for Elderly PBPK Model Development)."
      ),
      source_name = "Sex"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 10000L,
    n_studies = 116L,
    age_range = "30 to 100 years (the aging model; the underlying PK-Sim adult model extends below 30 y, and Fig. 2 shows organ weights from newborns, but the trajectories extracted here are defined on the 30-100 y decade knots of Tables 2 and 3)",
    weight_range = "Typical body weight 45.2 to 65.5 kg (female) and 55.6 to 74.4 kg (male) across the 30-100 y knots (Table 2)",
    height_range = "Typical height 147.9 to 163.0 cm (female) and 157.6 to 176.0 cm (male) across the 30-100 y knots (Table 2)",
    sex_female_pct = 50,
    race_ethnicity = "Healthy European Caucasians. Where European data were too sparse to describe an organ across the full age range, North American and Australian data were added (Methods, Workflow for Elderly PBPK Model Development; Discussion)",
    disease_state = paste(
      "Healthy aging. Inclusion criteria for the literature analysis",
      "required clear assignment of sex and race, a comprehensible",
      "statement of the method of analysis, and ruling out of effects",
      "of medication or disease on subject physiology (Methods,",
      "Collection of Data on Anthropometric, Anatomical and",
      "(Patho)physiological Changes with Age). Longitudinal studies",
      "were preferred. Organ volumes came mostly from autopsy series;",
      "in vivo measurements were used only for muscle, fat and blood",
      "flow."
    ),
    dose_range = "n/a (no exogenous drug PK is modelled by this system layer)",
    regions = "Europe, with North American and Australian data added for sparsely reported organs and for blood flow / cardiac output",
    n_studies_thompson = 19L,
    n_studies_additional = 97L,
    n_gfr_female = 1213L,
    n_gfr_male = 1081L,
    notes = paste(
      "n_studies is the 19 studies carried over from the Thompson et",
      "al. database plus the 97 studies added by this paper's own two",
      "literature searches (Methods). n_subjects is the size of the",
      "virtual cohort used for the physiological-consistency checks:",
      "a virtual 30-year-old was created in PK-Sim v5.5 and a virtual",
      "population of N = 5,000 was generated for EACH sex over the age",
      "range 30 to 100 years, i.e. 10,000 virtual subjects; it is not",
      "the number of real subjects behind the regressions. The GFR",
      "aging function of Equation 1 was fitted to body-surface-area-",
      "adjusted data from 1,213 females and 1,081 males, gathered",
      "mostly from prospective renal transplant donors assessed with",
      "exogenous markers (Results, Glomerular Filtration Rate).",
      "Secular trends were handled by generating a reference",
      "anthropometric measure for each 10-year time period over the",
      "past 60 years, used to normalise organ weights whenever a study",
      "reported no anthropometry (Methods).",
      "Blood-flow variability could not be robustly computed from the",
      "sparse perfusion literature, so the paper assumed a 5 %",
      "coefficient of variation for every blood flow rate; GFR",
      "variability was driven entirely by kidney-size variability.",
      "Neither is encoded here -- this model returns typical values."
    )
  )

  ini({
    # -----------------------------------------------------------------
    # Schlender 2016 Equation 1 (page 1581), the paper's own new
    # reverse sigmoid hyperbolic AGING function for kidney function,
    # optimised for both sexes with the Matlab lsqnonlin routine
    # starting at the age of 30 years:
    #
    #   Specific GFR_30-100y =
    #     F_PMA * ( 1 - Vmax * (age - 30 years)^Hill
    #                   / ( TA50^Hill + (age - 30 years)^Hill ) )
    #
    # Equation 1 is an embedded image in the published PDF and is
    # dropped by text converters; it was recovered from the publisher
    # equation graphic (Springer 40262_2016_422_Article_Equ1.gif) and
    # by rendering page 8 of the PDF at 400 dpi. All four constants are
    # stated in the prose immediately below the equation.
    #
    # Every value is a published point estimate of a fitted function
    # with no reported standard error, so all are fixed().
    # -----------------------------------------------------------------
    gfr_fpma <- fixed(26.6)
    label("Specific GFR after full maturation F_PMA (mL/min per 100 g kidney weight)")  # Schlender 2016 text below Eq. 1: the specific GFR from the Rhodin sigmoid hyperbolic MATURATION function after full maturation in the 40th week of gestation, 26.6 mL/min/100 g kidney weight
    gfr_vmax <- fixed(0.9)
    label("Maximal fractional decrease in specific GFR with aging Vmax (fraction of F_PMA)")  # Schlender 2016 text below Eq. 1: maximal decreasing rate factor Vmax = 0.9. The prose tags this with the unit 'mL/min/100 g kidney weight' carried over from F_PMA, but Eq. 1 multiplies Vmax by a unitless Hill ratio inside (1 - ...), so Vmax must be a dimensionless fraction; see vignette Errata
    gfr_ta50_female <- fixed(59)
    label("Aging half-time TA50 for females (years after the age of 30)")  # Schlender 2016 text below Eq. 1: aging half-time TA50 was 59 years for females
    gfr_ta50_male <- fixed(54)
    label("Aging half-time TA50 for males (years after the age of 30)")  # Schlender 2016 text below Eq. 1: aging half-time TA50 was 54 years for males
    gfr_hill <- fixed(1.5)
    label("Hill coefficient of the GFR aging function (unitless)")  # Schlender 2016 text below Eq. 1: the Hill coefficient was parameterized at 1.5
  })

  model({
    # The rxode2 time variable IS chronological age in years. The
    # trajectories are defined on the decade knots 30, 40, ..., 100 of
    # Schlender 2016 Tables 2 and 3.
    age <- time

    # ---- Linear interpolation between the ten-year age bins ----
    # The paper analysed the pooled data in 10-year age bins "as
    # recommended by the World Health Organization" and "interpolated
    # linearly" between them (Methods, Workflow for Elderly PBPK Model
    # Development), i.e. the physiology is a polygonal function of age.
    # Each rmp<N> is the clamped fraction of the decade ending at age N,
    # so a quantity is written as its value at 30 y plus the successive
    # decade increments. Clamping to [0, 1] reproduces the knot values
    # exactly and holds the trajectory flat outside 30-100 y rather than
    # extrapolating the last decade's slope, which the paper does not
    # support.
    rmp40 <- min(max((age - 30) / 10, 0), 1)
    rmp50 <- min(max((age - 40) / 10, 0), 1)
    rmp60 <- min(max((age - 50) / 10, 0), 1)
    rmp70 <- min(max((age - 60) / 10, 0), 1)
    rmp80 <- min(max((age - 70) / 10, 0), 1)
    rmp90 <- min(max((age - 80) / 10, 0), 1)
    rmp100 <- min(max((age - 90) / 10, 0), 1)

    # ---- Anthropometry (Schlender 2016 Table 2, rows 'Body weight' / 'Height') ----
    # Body weight, kg.
    bw_f <-
      60 + 3.3 * rmp40 + 1.9 * rmp50 + 0.3 * rmp60 - 3.3 * rmp70
        - 6.1 * rmp80 - 2.4 * rmp90 - 8.5 * rmp100
    bw_m <-
      73 + 1.4 * rmp40 - 1.4 * rmp50 - 1 * rmp60 - 0.9 * rmp70 - 3 * rmp80
        - 3.2 * rmp90 - 9.3 * rmp100
    bw <- SEXF * bw_f + (1 - SEXF) * bw_m

    # Body height, cm.
    ht_f <-
      163 - 2.1 * rmp40 - 0.9 * rmp50 - 1.7 * rmp60 - 2.9 * rmp70
        - 3.7 * rmp80 - 1.7 * rmp90 - 2.1 * rmp100
    ht_m <-
      176 - 0.9 * rmp40 - 1.2 * rmp50 - 2.7 * rmp60 - 3.8 * rmp70
        - 1.9 * rmp80 - 1.7 * rmp90 - 6.2 * rmp100
    ht <- SEXF * ht_f + (1 - SEXF) * ht_m

    # ---- Organ / tissue masses, g (Schlender 2016 Table 2, 'Organ masses') ----
    # Mass of blood pools (arterial + venous + portal-vein blood; Table 2 footnote b), g.
    mass_blood_f <-
      1899.4 + 21.7 * rmp40 + 40.1 * rmp50 - 14.3 * rmp60 - 99.2 * rmp70
        - 155.1 * rmp80 - 58.1 * rmp90 - 277.5 * rmp100
    mass_blood_m <-
      2264.6 + 4.9 * rmp40 - 34.6 * rmp50 - 54.7 * rmp60 - 68.9 * rmp70
        - 68.8 * rmp80 - 66 * rmp90 - 237 * rmp100
    mass_blood <- SEXF * mass_blood_f + (1 - SEXF) * mass_blood_m

    # Mass of bone, g.
    mass_bone_f <-
      9121.5 - 235.2 * rmp40 - 210.2 * rmp50 - 274.3 * rmp60
        - 590.2 * rmp70 - 830.2 * rmp80 - 362.4 * rmp90 - 355.5 * rmp100
    mass_bone_m <-
      11817.8 + 30.2 * rmp40 - 366.4 * rmp50 - 390 * rmp60 - 736.6 * rmp70
        - 717.9 * rmp80 - 668.1 * rmp90 - 668.1 * rmp100
    mass_bone <- SEXF * mass_bone_f + (1 - SEXF) * mass_bone_m

    # Mass of brain, g.
    mass_brain_f <-
      1357 - 4.7 * rmp40 - 4.7 * rmp50 - 32.5 * rmp60 - 27.8 * rmp70
        - 45.9 * rmp80 - 72 * rmp90 - 72 * rmp100
    mass_brain_m <-
      1508.8 - 1.9 * rmp40 - 4 * rmp50 - 39.6 * rmp60 - 29 * rmp70
        - 34.2 * rmp80 - 28 * rmp90 - 27.9 * rmp100
    mass_brain <- SEXF * mass_brain_f + (1 - SEXF) * mass_brain_m

    # Mass of fat, g.
    mass_fat_f <-
      19348 + 654.5 * rmp40 + 3812.9 * rmp50 + 3177.6 * rmp60
        + 891.7 * rmp70 - 2190.2 * rmp80 - 673.9 * rmp90 - 4849.2 * rmp100
    mass_fat_m <-
      14868 + 1441.7 * rmp40 + 366.8 * rmp50 + 1244.2 * rmp60
        + 1526.7 * rmp70 + 500.3 * rmp80 + 743.3 * rmp90 - 3608.7 * rmp100
    mass_fat <- SEXF * mass_fat_f + (1 - SEXF) * mass_fat_m

    # Mass of gonads, g.
    mass_gonads_f <-
      13.1 - 0.1 * rmp40 - 6.4 * rmp50 - 1.3 * rmp60 - 0.1 * rmp70
        - 0.1 * rmp80 + 0 * rmp90 - 0.2 * rmp100
    mass_gonads_m <-
      40.3 + 0.5 * rmp40 - 5.8 * rmp50 - 1.5 * rmp60 - 1.6 * rmp70
        - 1.3 * rmp80 - 0.3 * rmp90 - 0.8 * rmp100
    mass_gonads <- SEXF * mass_gonads_f + (1 - SEXF) * mass_gonads_m

    # Mass of heart, g.
    mass_heart_f <-
      328.4 + 11.6 * rmp40 + 15.6 * rmp50 + 21.5 * rmp60 + 22 * rmp70
        + 14.4 * rmp80 + 17.7 * rmp90 - 18.3 * rmp100
    mass_heart_m <-
      417.2 + 17.1 * rmp40 + 4.8 * rmp50 + 15.5 * rmp60 + 10.2 * rmp70
        - 23.4 * rmp80 - 3.5 * rmp90 - 12.8 * rmp100
    mass_heart <- SEXF * mass_heart_f + (1 - SEXF) * mass_heart_m

    # Mass of kidney, g.
    mass_kidney_f <-
      403.4 - 1.7 * rmp40 - 1.2 * rmp50 - 17.1 * rmp60 - 19.2 * rmp70
        - 39 * rmp80 - 16.4 * rmp90 - 13 * rmp100
    mass_kidney_m <-
      437.7 + 38.2 * rmp40 - 7.8 * rmp50 - 12.5 * rmp60 - 13.9 * rmp70
        - 46 * rmp80 - 29 * rmp90 - 10.8 * rmp100
    mass_kidney <- SEXF * mass_kidney_f + (1 - SEXF) * mass_kidney_m

    # Mass of liver, g.
    mass_liver_f <-
      1905.5 - 24.5 * rmp40 - 13.3 * rmp50 - 187.9 * rmp60 - 174.9 * rmp70
        - 81.5 * rmp80 - 93.9 * rmp90 - 140.1 * rmp100
    mass_liver_m <-
      2357.8 - 33.3 * rmp40 - 118.5 * rmp50 - 164.7 * rmp60 - 326.7 * rmp70
        - 297.3 * rmp80 - 93.3 * rmp90 - 122 * rmp100
    mass_liver <- SEXF * mass_liver_f + (1 - SEXF) * mass_liver_m

    # Mass of lung, g.
    mass_lung_f <-
      1009.5 + 11.9 * rmp40 + 17.4 * rmp50 - 15 * rmp60 - 18.4 * rmp70
        - 156.3 * rmp80 - 44.6 * rmp90 - 167.2 * rmp100
    mass_lung_m <-
      1294.3 + 40.4 * rmp40 + 2.7 * rmp50 + 11.6 * rmp60 - 77.2 * rmp70
        - 133.8 * rmp80 - 130.6 * rmp90 - 29.4 * rmp100
    mass_lung <- SEXF * mass_lung_f + (1 - SEXF) * mass_lung_m

    # Mass of muscle, g.
    mass_muscle_f <-
      20276.2 + 2782.1 * rmp40 - 1758.2 * rmp50 - 2328.9 * rmp60
        - 3224.8 * rmp70 - 2378.6 * rmp80 - 973.5 * rmp90 - 2333.7 * rmp100
    mass_muscle_m <-
      32338.6 - 20.4 * rmp40 - 1137.6 * rmp50 - 1580.6 * rmp60
        - 1037.9 * rmp70 - 2037.3 * rmp80 - 2867.8 * rmp90 - 4123 * rmp100
    mass_muscle <- SEXF * mass_muscle_f + (1 - SEXF) * mass_muscle_m

    # Mass of pancreas, g.
    mass_pancreas_f <-
      169.5 + 0.8 * rmp40 - 6.5 * rmp50 - 5.1 * rmp60 - 10.2 * rmp70
        - 16.6 * rmp80 - 3.9 * rmp90 - 5.4 * rmp100
    mass_pancreas_m <-
      190.3 + 0.5 * rmp40 - 7.2 * rmp50 - 5.4 * rmp60 - 12.8 * rmp70
        - 15.5 * rmp80 - 6.2 * rmp90 - 4.2 * rmp100
    mass_pancreas <- SEXF * mass_pancreas_f + (1 - SEXF) * mass_pancreas_m

    # Mass of skin, g.
    mass_skin_f <-
      2723.5 + 49.7 * rmp40 + 63 * rmp50 - 6.7 * rmp60 - 103.7 * rmp70
        - 166.7 * rmp80 - 65.3 * rmp90 - 270.3 * rmp100
    mass_skin_m <-
      3760.9 + 29.3 * rmp40 - 45.1 * rmp50 - 53.1 * rmp60 - 56.8 * rmp70
        - 98.2 * rmp80 - 97.3 * rmp90 - 317.1 * rmp100
    mass_skin <- SEXF * mass_skin_f + (1 - SEXF) * mass_skin_m

    # Mass of spleen, g.
    mass_spleen_f <-
      219.2 - 21.5 * rmp40 - 7 * rmp50 - 8.3 * rmp60 - 17.9 * rmp70
        - 14.8 * rmp80 - 50.3 * rmp90 - 20.7 * rmp100
    mass_spleen_m <-
      243.4 - 21.5 * rmp40 - 13.7 * rmp50 - 11.2 * rmp60 - 10.4 * rmp70
        - 26.2 * rmp80 - 27 * rmp90 - 25 * rmp100
    mass_spleen <- SEXF * mass_spleen_f + (1 - SEXF) * mass_spleen_m

    # Gastrointestinal organ mass is held CONSTANT across 30-100 y
    # (Schlender 2016 Table 2 footnote a: 1274.5 g female, 1304.8 g male).
    mass_gi <- SEXF * 1274.5 + (1 - SEXF) * 1304.8

    # ---- Cardiac index and organ blood flows (Schlender 2016 Table 3) ----
    # Cardiac index, L/min/m^2.
    cardiac_index_f <-
      3.34 - 0.01 * rmp40 - 0.25 * rmp50 - 0.19 * rmp60 - 0.14 * rmp70
        - 0.07 * rmp80 - 0.16 * rmp90 - 0.11 * rmp100
    cardiac_index_m <-
      3.23 - 0.01 * rmp40 - 0.13 * rmp50 - 0.23 * rmp60 - 0.28 * rmp70
        - 0.22 * rmp80 - 0.18 * rmp90 - 0.14 * rmp100
    cardiac_index <- SEXF * cardiac_index_f + (1 - SEXF) * cardiac_index_m

    # Blood flow to adipose, mL/min.
    flow_adipose_f <-
      501.1 + 27.1 * rmp40 + 49.8 * rmp50 + 77.3 * rmp60 + 29.4 * rmp70
        - 54.5 * rmp80 - 23.2 * rmp90 - 115.5 * rmp100
    flow_adipose_m <-
      324.7 + 27.4 * rmp40 + 5.2 * rmp50 + 24.1 * rmp60 + 26.6 * rmp70
        + 11.8 * rmp80 + 9.4 * rmp90 - 71.9 * rmp100
    flow_adipose <- SEXF * flow_adipose_f + (1 - SEXF) * flow_adipose_m

    # Blood flow to cerebral, mL/min.
    flow_cerebral_f <-
      707.8 - 20.9 * rmp40 - 20.8 * rmp50 - 34 * rmp60 - 30.9 * rmp70
        - 38.4 * rmp80 - 48.6 * rmp90 - 46.6 * rmp100
    flow_cerebral_m <-
      779.7 - 32.1 * rmp40 - 19.7 * rmp50 - 35 * rmp60 - 42.3 * rmp70
        - 26.8 * rmp80 - 29.4 * rmp90 - 58.3 * rmp100
    flow_cerebral <- SEXF * flow_cerebral_f + (1 - SEXF) * flow_cerebral_m

    # Blood flow to gonads, mL/min.
    flow_gonads_f <-
      1.2 + 0 * rmp40 - 0.6 * rmp50 - 0.1 * rmp60 + 0 * rmp70 + 0 * rmp80
        + 0 * rmp90 - 0.1 * rmp100
    flow_gonads_m <-
      3.2 + 0.1 * rmp40 - 0.5 * rmp50 - 0.1 * rmp60 - 0.1 * rmp70
        - 0.1 * rmp80 - 0.1 * rmp90 + 0 * rmp100
    flow_gonads <- SEXF * flow_gonads_f + (1 - SEXF) * flow_gonads_m

    # Blood flow to myocardial, mL/min.
    flow_myocardial_f <-
      295 + 16.2 * rmp40 + 2 * rmp50 + 49.9 * rmp60 + 13.2 * rmp70
        + 69.5 * rmp80 + 19 * rmp90 - 19.6 * rmp100
    flow_myocardial_m <-
      260.1 + 11 * rmp40 + 6.9 * rmp50 + 21.1 * rmp60 + 9.6 * rmp70
        + 26 * rmp80 - 2.6 * rmp90 - 9.8 * rmp100
    flow_myocardial <- SEXF * flow_myocardial_f + (1 - SEXF) * flow_myocardial_m

    # Blood flow to renal, mL/min.
    flow_renal_f <-
      1121 - 59.2 * rmp40 - 118 * rmp50 - 118 * rmp60 - 118 * rmp70
        - 118 * rmp80 - 118 * rmp90 - 118 * rmp100
    flow_renal_m <-
      1325 + 14.3 * rmp40 - 223 * rmp50 - 203 * rmp60 - 183 * rmp70
        - 163 * rmp80 - 143 * rmp90 - 123 * rmp100
    flow_renal <- SEXF * flow_renal_f + (1 - SEXF) * flow_renal_m

    # Blood flow to splanchnic (GI + pancreatic + splenic = portal; Table 3 footnote a), mL/min.
    flow_splanchnic_f <-
      1239 + 15.3 * rmp40 - 67.1 * rmp50 - 171.9 * rmp60 - 152.7 * rmp70
        - 91.1 * rmp80 - 92.4 * rmp90 - 108.6 * rmp100
    flow_splanchnic_m <-
      1235 + 6.9 * rmp40 - 4.2 * rmp50 - 181 * rmp60 - 243.6 * rmp70
        - 202.5 * rmp80 - 97.7 * rmp90 - 99.4 * rmp100
    flow_splanchnic <- SEXF * flow_splanchnic_f + (1 - SEXF) * flow_splanchnic_m

    # Blood flow to hepatic arterial, mL/min.
    flow_hepatic_f <-
      383.4 + 4.7 * rmp40 - 20.8 * rmp50 - 53.1 * rmp60 - 47.3 * rmp70
        - 28.2 * rmp80 - 28.6 * rmp90 - 33.6 * rmp100
    flow_hepatic_m <-
      423 + 2.4 * rmp40 - 1.5 * rmp50 - 62 * rmp60 - 83.4 * rmp70
        - 69.4 * rmp80 - 33.4 * rmp90 - 34.1 * rmp100
    flow_hepatic <- SEXF * flow_hepatic_f + (1 - SEXF) * flow_hepatic_m

    # Blood flow to muscle, mL/min.
    flow_muscle_f <-
      665.1 + 109.7 * rmp40 - 138.2 * rmp50 - 68.2 * rmp60 - 88.6 * rmp70
        - 74.2 * rmp80 - 36.6 * rmp90 - 70.2 * rmp100
    flow_muscle_m <-
      1105.7 - 18.9 * rmp40 - 50 * rmp50 - 61.8 * rmp60 - 52.7 * rmp70
        - 60.6 * rmp80 - 109.8 * rmp90 - 126 * rmp100
    flow_muscle <- SEXF * flow_muscle_f + (1 - SEXF) * flow_muscle_m

    # Blood flow to skeleton, mL/min.
    flow_skeleton_f <-
      294.9 - 7.6 * rmp40 - 6.8 * rmp50 - 8.9 * rmp60 - 19 * rmp70
        - 26.9 * rmp80 - 11.7 * rmp90 - 11.5 * rmp100
    flow_skeleton_m <-
      324.9 + 0.9 * rmp40 - 10.1 * rmp50 - 10.7 * rmp60 - 20.3 * rmp70
        - 19.7 * rmp80 - 18.4 * rmp90 - 18.4 * rmp100
    flow_skeleton <- SEXF * flow_skeleton_f + (1 - SEXF) * flow_skeleton_m

    # Blood flow to skin, mL/min.
    flow_skin_f <-
      295.7 + 5.4 * rmp40 + 6.8 * rmp50 - 0.7 * rmp60 - 11.3 * rmp70
        - 18.1 * rmp80 - 7.1 * rmp90 - 29.3 * rmp100
    flow_skin_m <-
      325.1 + 2.6 * rmp40 - 3.9 * rmp50 - 4.6 * rmp60 - 4.9 * rmp70
        - 8.5 * rmp80 - 8.4 * rmp90 - 27.4 * rmp100
    flow_skin <- SEXF * flow_skin_f + (1 - SEXF) * flow_skin_m
    # ---- Cardiac output ----
    # Schlender 2016 accounted for age-related changes in blood flow
    # distribution as changes in cardiac output, and states that "the
    # sum of all blood flows should be equal to CO" (Methods). Cardiac
    # output is recovered from the Table 3 cardiac index and the body
    # surface area implied by the Table 2 anthropometry. BSA uses the
    # Du Bois equation; the paper reports the cardiac index per m^2 of
    # BSA (Results, Cardiac Output Distribution) but does not state
    # which BSA formula it used, so the choice is an assumption -- see
    # the vignette Errata, which shows the Table 3 flows sum to within
    # 1-3 % of this cardiac output at every one of the 16 age-by-sex
    # knots.
    bsa <- 0.007184 * ht^0.725 * bw^0.425
    cardiac_output <- cardiac_index * bsa * 1000

    # ---- Glomerular filtration rate (Schlender 2016 Equation 1) ----
    # TA50 is sex-specific: 59 y (female) versus 54 y (male).
    gfr_ta50 <- SEXF * gfr_ta50_female + (1 - SEXF) * gfr_ta50_male
    # Equation 1 is defined from the age of 30 years onward. The age
    # offset is floored at 0 so the expression stays real below 30 y
    # (a negative base raised to the non-integer power 1.5 is NaN) and
    # returns the unaged F_PMA there.
    gfr_age <- max(age - 30, 0)
    gfr_specific <- gfr_fpma *
      (1 - gfr_vmax * gfr_age^gfr_hill /
        (gfr_ta50^gfr_hill + gfr_age^gfr_hill))
    # Absolute GFR. The paper parameterised the specific GFR function
    # "using the observed GFR values and assuming the mean kidney
    # weight at the respective age" (Results, Glomerular Filtration
    # Rate), so absolute GFR is recovered by multiplying by the Table 2
    # kidney mass. The /100 converts the per-100-g basis to per-gram.
    gfr <- gfr_specific * mass_kidney / 100
  })
}
