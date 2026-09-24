deKock_2018_sulfadoxinePyrimethamine <- function() {
  description <- paste(
    "Joint popPK model for the antimalarial fixed-dose combination of",
    "sulfadoxine and pyrimethamine, pooled from four published African",
    "trials at eight study sites in 801 patients (415 children, 386",
    "adults) with uncomplicated Plasmodium falciparum malaria",
    "(de Kock 2018). Each drug has one-compartment disposition with",
    "first-order absorption and elimination, fitted jointly so that the",
    "between-subject random effects on the two clearances are correlated",
    "(60%). Apparent clearance and volume are allometrically scaled with",
    "total body weight (exponents 0.75 and 1, reference WT = 18 kg, the",
    "population median), and clearance carries a sigmoidal maturation",
    "function of postgestational age (PGA50 = 8.12 months for sulfadoxine",
    "and 11.9 months for pyrimethamine). Underweight-for-age children",
    "have lower relative bioavailability through a 'hockey stick' effect",
    "of the weight-for-age Z-score that switches on below Z = -2:",
    "15.3% lower for sulfadoxine and 26.7% lower for pyrimethamine per",
    "Z-score unit below -2. Residual site-specific scaling is applied to",
    "the predicted concentrations (one group per drug, plus a separate",
    "Namaacha group for pyrimethamine), and pyrimethamine apparent",
    "clearance is 54.9% lower in the Bell et al. study, the only study",
    "that assayed liquid whole blood in a different laboratory. Relative",
    "bioavailability is fixed at 1 for both drugs and carries",
    "between-subject variability."
  )
  reference <- paste(
    "de Kock M, Tarning J, Workman L, Allen EN, Tekete MM, Djimde AA,",
    "Bell DJ, Ward SA, Barnes KI, Denti P.",
    "Population pharmacokinetic properties of sulfadoxine and",
    "pyrimethamine: a pooled analysis to inform optimal dosing in African",
    "children with uncomplicated malaria.",
    "Antimicrob Agents Chemother. 2018;62(5):e01370-17.",
    "doi:10.1128/AAC.01370-17.",
    sep = " "
  )
  vignette <- "deKock_2018_sulfadoxinePyrimethamine"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ug/mL for sulfadoxine (whole blood); ng/mL for pyrimethamine (whole blood)"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against de Kock 2018 Materials and Methods
  # ("All sulfadoxine-pyrimethamine concentrations were measured by use of
  # capillary whole-blood dried spots on filter paper ... except in the study
  # conducted by Bell et al., in which concentrations were measured in liquid
  # samples of either capillary or venous whole blood"). Unlike the sibling
  # de Kock 2017 IPTp model, this model is parameterised DIRECTLY on
  # whole-blood concentrations -- there is no hematocrit / plasma conversion
  # step -- so CL/F and V/F are apparent whole-blood quantities.
  compartmentData <- list(
    depot = list(
      analyte = "sulfadoxine",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(analyte = "sulfadoxine", units = "mg", specimen = "whole blood", verified = TRUE),
    depot_pyra = list(
      analyte = "pyrimethamine",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central_pyra = list(analyte = "pyrimethamine", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling on the population median 18 kg (de Kock 2018",
        "Table 2 footnote b), applied to the apparent volume of both drugs",
        "with exponent 1 and to the apparent clearance of both drugs with",
        "exponent 0.75 (Materials and Methods, 'The effect of body size was",
        "taken into account by using allometric scaling with total body",
        "weight to adjust all volumes with an exponent of 1 and flow rates",
        "... with an exponent of 0.75'). No height was recorded, so fat-free",
        "mass could not be tested as an alternative size descriptor. The",
        "cohort spans 5 to 80 kg (Table 1 median 18 kg, IQR 12 to 50)."
      ),
      source_name = "WT"
    ),
    PAGE = list(
      description = "Postgestational age (age measured from conception)",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying. Drives the sigmoidal maturation function on the",
        "apparent clearance of each drug (de Kock 2018 Materials and",
        "Methods: MAT = PGA^gamma / (PGA^gamma + PGA50^gamma)). The paper",
        "calls this covariate PGA, 'months after conception'; it is carried",
        "here on the canonical PAGE column because the two are the same",
        "quantity to within the roughly two-week (0.46 month) offset",
        "between conception and the last menstrual period, which is the",
        "convention the closely related antimalarial models",
        "Ali_2018_amodiaquine.R and Denti_2018_levofloxacin.R also use.",
        "Compute as postnatal age in months + 9 months of assumed term",
        "gestation; for adults, AGE_years * 12 + 9. The pooled data set",
        "contained no child under 12 months of age, so the maturation",
        "curve below 21 months PGA is an extrapolation -- reflected in the",
        "56% RSE on the sulfadoxine PGA50 (Discussion, final paragraph",
        "before Conclusions)."
      ),
      source_name = "PGA"
    ),
    WAZ = list(
      description = "Weight-for-age Z-score (WHO Multicentre Growth Reference Study standard)",
      units = "unitless (z-score; standard-deviation units)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline. The paper's only available measure of nutritional status;",
        "no height or mid-upper-arm circumference was recorded, so stunting",
        "and wasting could not be distinguished. Computed with the WHO",
        "igrowup.standard R macro against the WHO Multicentre Growth",
        "Reference Study curves (Materials and Methods, 'The nutritional",
        "status of children'). Enters a 'hockey stick' effect on relative",
        "bioavailability that is inert at WAZ >= -2 and falls linearly",
        "below it. Z-scores are defined only for children under 5 years of",
        "age; supply 0 for older children and adults so the effect is",
        "inert. No patient in the pooled analysis had a Z-score below",
        "-4.27 (Materials and Methods, simulation paragraph); the linear",
        "hockey-stick form turns the pyrimethamine bioavailability negative",
        "below WAZ = -5.74, so do not extrapolate past the observed range."
      ),
      source_name = "Z-score"
    ),
    REGION_MPUMALANGA = list(
      description = "Mpumalanga (South Africa) study-site indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the sulfadoxine reference sites Magude, Bancoumana, Bela Vista, Catuane and Chileka)",
      notes = paste(
        "1 = the Mpumalanga site of Barnes et al. (n = 122); 0 otherwise.",
        "One of the three sites sharing a single estimated -39.7% scaling",
        "of the OBSERVED sulfadoxine concentrations (de Kock 2018 Table 2,",
        "'Scaling on observations at site(s)', row 'Mpumalanga, Boane,",
        "Namaacha'). Carries no pyrimethamine effect -- Mpumalanga is in",
        "the pyrimethamine reference group."
      ),
      source_name = "SITE = 'Mpumalanga'"
    ),
    REGION_BOANE = list(
      description = "Boane (Mozambique) study-site indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the sulfadoxine reference sites Magude, Bancoumana, Bela Vista, Catuane and Chileka)",
      notes = paste(
        "1 = the Boane site of Allen et al. (n = 78); 0 otherwise. Shares",
        "the -39.7% sulfadoxine observation scaling with Mpumalanga and",
        "Namaacha (de Kock 2018 Table 2). Carries no pyrimethamine effect",
        "-- Boane is in the pyrimethamine reference group."
      ),
      source_name = "SITE = 'Boane'"
    ),
    REGION_NAMAACHA = list(
      description = "Namaacha (Mozambique) study-site indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Namaacha sites)",
      notes = paste(
        "1 = the Namaacha site; 0 otherwise. Namaacha contributed patients",
        "to both Barnes et al. (n = 91) and Allen et al. (n = 72) and the",
        "site is treated as one group in the covariate model. It is the",
        "only site that carries a scaling effect for BOTH drugs: it shares",
        "the -39.7% sulfadoxine observation scaling with Mpumalanga and",
        "Boane, and it is alone in the -22.0% pyrimethamine observation",
        "scaling group (de Kock 2018 Table 2)."
      ),
      source_name = "SITE = 'Namaacha'"
    ),
    REGION_BANCOUMANA = list(
      description = "Bancoumana (Mali) study-site indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the pyrimethamine reference sites Magude, Mpumalanga, Boane and Chileka)",
      notes = paste(
        "1 = the Bancoumana site of Tekete et al. (n = 114); 0 otherwise.",
        "One of the three sites sharing a single estimated +20.2% scaling",
        "of the OBSERVED pyrimethamine concentrations (de Kock 2018",
        "Table 2, row 'Bancoumana, Bela Vista, Catuane'). Carries no",
        "sulfadoxine effect -- Bancoumana is in the sulfadoxine reference",
        "group."
      ),
      source_name = "SITE = 'Bancoumana'"
    ),
    REGION_BELAVISTA = list(
      description = "Bela Vista (Mozambique) study-site indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the pyrimethamine reference sites Magude, Mpumalanga, Boane and Chileka)",
      notes = paste(
        "1 = the Bela Vista site of Barnes et al. (n = 65); 0 otherwise.",
        "Shares the +20.2% pyrimethamine observation scaling with",
        "Bancoumana and Catuane (de Kock 2018 Table 2). Carries no",
        "sulfadoxine effect -- Bela Vista is in the sulfadoxine reference",
        "group."
      ),
      source_name = "SITE = 'Bela Vista'"
    ),
    REGION_CATUANE = list(
      description = "Catuane (Mozambique) study-site indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the pyrimethamine reference sites Magude, Mpumalanga, Boane and Chileka)",
      notes = paste(
        "1 = the Catuane site of Allen et al. (n = 33); 0 otherwise. Shares",
        "the +20.2% pyrimethamine observation scaling with Bancoumana and",
        "Bela Vista (de Kock 2018 Table 2). Carries no sulfadoxine effect",
        "-- Catuane is in the sulfadoxine reference group. de Kock 2018",
        "spells the site both 'Catuane' (Table 2 footnotes, Results text)",
        "and 'Cutuane' (Table 1 column header); Catuane is the spelling",
        "used wherever the covariate effect itself is reported."
      ),
      source_name = "SITE = 'Catuane'"
    ),
    STUDY_BELL = list(
      description = "Bell et al. study-cohort indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the Barnes, Tekete and Allen studies)",
      notes = paste(
        "1 = patient from the Bell et al. study at the Chileka site in",
        "Malawi (n = 102); 0 otherwise. Multiplicatively reduces apparent",
        "pyrimethamine clearance by 54.9% (de Kock 2018 Table 2,",
        "'Difference from clearance in reference 5'). Encoded as a STUDY",
        "rather than a REGION indicator because the paper attributes the",
        "difference to the study's analytical method, not its geography:",
        "Bell et al. was 'the only one that assayed whole-blood liquid",
        "samples (capillary blood dried-spot samples were assayed in all",
        "the other studies), and its samples were assayed in a different",
        "lab' (Discussion). Study and site are perfectly confounded here",
        "-- Chileka is the only Bell et al. site -- so the paper could not",
        "separate matrix, assay method and population. Chileka is in the",
        "reference group for BOTH drugs' observation scaling, so a Bell",
        "et al. subject sets STUDY_BELL = 1 and every REGION_ indicator",
        "to 0."
      ),
      source_name = "Study = Bell et al. (reference 5)"
    )
  )

  # de Kock 2018 Materials and Methods, final covariate paragraph, and Results:
  # "No other predefined covariates (sex, baseline hemoglobin, dose [milligrams
  # per kilogram of body weight], concomitant medications, and baseline
  # parasitemia) were found to be significant, and these were therefore
  # excluded from the model." Sample blood matrix was also prespecified but is
  # perfectly confounded with STUDY_BELL and so is documented in that entry
  # rather than carried separately. Total dose in mg/kg was screened as a
  # continuous covariate and is likewise not carried as a column.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Prespecified and screened by the Wald-test full approach in de Kock 2018 but not significant; 398 / 801 (50%) of the cohort were male (Table 1)."
    ),
    HGB = list(
      description = "Baseline hemoglobin",
      units = "g/dL",
      type = "continuous",
      notes = "Prespecified and screened as a median-centred continuous covariate on the log-transformed PK parameters but not significant; cohort median 11 g/dL (IQR 10 to 12), Table 1."
    ),
    PARA = list(
      description = "Baseline Plasmodium falciparum parasitaemia",
      units = "parasites/uL",
      type = "continuous",
      notes = "Prespecified and screened as a median-centred continuous covariate but not significant; cohort geometric mean 21,700 parasites/uL (95% range 4,532 to 63,083), Table 1."
    ),
    CONMED_AMODIAQUINE = list(
      description = "Amodiaquine co-administration indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (sulfadoxine-pyrimethamine alone or with another partner drug)",
      notes = "Prespecified and screened but not significant; 60 / 801 (7%) received amodiaquine alongside sulfadoxine-pyrimethamine (Table 1)."
    ),
    CONMED_CHLOROQUINE = list(
      description = "Chloroquine co-administration indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (sulfadoxine-pyrimethamine alone or with another partner drug)",
      notes = "Prespecified and screened but not significant; 26 / 801 (3%) received chloroquine alongside sulfadoxine-pyrimethamine (Table 1)."
    ),
    CONMED_ARTESUNATE = list(
      description = "Artesunate co-administration indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (sulfadoxine-pyrimethamine alone or with another partner drug)",
      notes = "Prespecified and screened but not significant; 225 / 801 (28%) received artesunate alongside sulfadoxine-pyrimethamine (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 801L,
    n_children = 415L,
    n_adults = 386L,
    n_studies = 4L,
    n_sites = 8L,
    age_range = paste(
      "1 to over 20 years. Table 1 age strata: 32 (4%) under 2 years,",
      "383 (47%) over 2 to 5 years, 197 (25%) over 5 to 20 years, and",
      "189 (24%) 20 years and over. No patient under 12 months of age",
      "contributed pharmacokinetic data."
    ),
    weight_range = "Median 18 kg (IQR 12 to 50); per-site medians 11 to 55 kg (Table 1). The allometric reference is the pooled median, 18 kg.",
    sex_female_pct = 100 * (1 - 398 / 801),
    race_ethnicity = "Sub-Saharan African (Mozambique, South Africa, Mali, Malawi); not further reported.",
    disease_state = paste(
      "Nonpregnant patients with uncomplicated Plasmodium falciparum",
      "malaria. Median baseline hemoglobin 11 g/dL (IQR 10 to 12);",
      "geometric mean baseline parasitaemia 21,700 parasites/uL (95%",
      "range 4,532 to 63,083). Among the 383 children under 5 years of",
      "age with a nutrition score, 326 (85%) were of normal weight for",
      "age, 41 (11%) had -3 <= Z-score < -2, and 16 (4%) had a Z-score",
      "below -3 (Table 1)."
    ),
    dose_range = paste(
      "Single oral dose of the 500 mg sulfadoxine / 25 mg pyrimethamine",
      "tablet. Adults received 1,500 mg / 75 mg; children received at",
      "least 25 mg/kg sulfadoxine and 1.25 mg/kg pyrimethamine, given as",
      "half-tablet multiples by the per-study weight bands of Table 1.",
      "Sulfadoxine-pyrimethamine was given alone (250 children, 304",
      "adults) or with chloroquine (34 children), artesunate (85",
      "children, 113 adults), or amodiaquine (29 children)."
    ),
    sampling = paste(
      "Six to nine samples per patient, at least predose and on days 1,",
      "3, 7, 14, 21 and 28; several sites also sampled on days 2 and 42.",
      "8,981 samples were collected; 259 (2.88%) were excluded as",
      "biologically implausible outliers by a model-based normalised",
      "prediction distribution error screen, leaving 4,567 sulfadoxine",
      "and 4,155 pyrimethamine concentrations."
    ),
    regions = "Sub-Saharan Africa: Mozambique (Bela Vista, Namaacha, Boane, Catuane, Magude), South Africa (Mpumalanga), Mali (Bancoumana), Malawi (Chileka).",
    notes = paste(
      "Pooled individual patient data contributed to the Worldwide",
      "Antimalarial Resistance Network (WWARN) repository from four",
      "previously published trials: Barnes et al. (Bela Vista,",
      "Mpumalanga, Namaacha), Bell et al. (Chileka), Tekete et al.",
      "(Bancoumana) and Allen et al. (Boane, Catuane, Magude,",
      "Namaacha). Estimation used the SAEM algorithm in Monolix Suite",
      "2016R1. Concentrations below the limit of quantification were",
      "handled with Beal's M3 method at the two Mozambican sites that",
      "censored them (below 10 ng/mL pyrimethamine and 10 ug/mL",
      "sulfadoxine); elsewhere the laboratory's reported sub-LLOQ values",
      "were used as-is. 152 (18.9%) sulfadoxine and 125 (15.6%)",
      "pyrimethamine predose samples were detectable, and for those",
      "patients the fitted model initialised the disposition compartments",
      "to the observed predose concentration -- an estimation-time device",
      "for residual drug from a previous treatment that is not part of",
      "the forward-simulation model encoded here."
    )
  )

  ini({
    # =====================================================================
    # Sulfadoxine structural parameters. de Kock 2018 Table 2, "Sulfadoxine"
    # column. CL/F and V/F are the values for a FULLY MATURED patient at the
    # reference weight of 18 kg (Table 2 footnotes b and f), in the reference
    # site group, with a weight-for-age Z-score at or above -2.
    # Concentrations are whole blood; the apparent volume is therefore an
    # apparent whole-blood volume.
    # =====================================================================
    lcl <- log(0.0264)
    label("Sulfadoxine apparent clearance CL/F at 18 kg, fully matured (L/h)") # Table 2 sulfadoxine: CL/F = 0.0264 L/h, RSE 3%
    lvc <- log(5.29)
    label("Sulfadoxine apparent volume of distribution V/F at 18 kg (L)") # Table 2 sulfadoxine: V/F = 5.29 L, RSE 2%
    lka <- log(0.521)
    label("Sulfadoxine first-order absorption rate constant ka (1/h)") # Table 2 sulfadoxine: ka = 0.521 /h, RSE 16%
    lfdepot <- fixed(log(1))
    label("Sulfadoxine relative bioavailability F (fraction)") # Table 2 sulfadoxine: F = 1 fixed

    # =====================================================================
    # Pyrimethamine structural parameters. de Kock 2018 Table 2,
    # "Pyrimethamine" column, on the same reference conditions.
    # =====================================================================
    lcl_pyra <- log(0.829)
    label("Pyrimethamine apparent clearance CL/F at 18 kg, fully matured (L/h)") # Table 2 pyrimethamine: CL/F = 0.829 L/h, RSE 3%
    lvc_pyra <- log(91.4)
    label("Pyrimethamine apparent volume of distribution V/F at 18 kg (L)") # Table 2 pyrimethamine: V/F = 91.4 L, RSE 3%
    lka_pyra <- log(1.40)
    label("Pyrimethamine first-order absorption rate constant ka (1/h)") # Table 2 pyrimethamine: ka = 1.40 /h, RSE 80%
    lfdepot_pyra <- fixed(log(1))
    label("Pyrimethamine relative bioavailability F (fraction)") # Table 2 pyrimethamine: F = 1 fixed

    # =====================================================================
    # Allometric exponents. de Kock 2018 Materials and Methods: volumes are
    # adjusted with an exponent of 1 and flow rates with an exponent of 0.75,
    # citing Anderson and Holford. They are structural constants, not
    # estimates -- Table 2 reports no RSE for them -- and the same pair
    # applies to both drugs, so they are not drug-suffixed.
    # =====================================================================
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on the apparent clearance of both drugs (unitless)")
    e_wt_vc <- fixed(1)
    label("Allometric exponent on the apparent volume of both drugs (unitless)")

    # =====================================================================
    # Maturation of clearance. de Kock 2018 Materials and Methods:
    # MAT = PGA^gamma / (PGA^gamma + PGA50^gamma), a sigmoidal function of
    # postgestational age applied to clearance. Because Table 2's CL/F is
    # defined as the clearance 'for a fully matured child' (footnote f), the
    # factor is applied directly and is NOT renormalised to a reference age.
    # =====================================================================
    pma50_cl <- 8.12
    label("Sulfadoxine postgestational age at 50% of mature clearance (months)") # Table 2 sulfadoxine: PGA50 = 8.12 mo after conception, RSE 56%
    hill_cl <- 3.20
    label("Sulfadoxine maturation Hill coefficient (unitless)") # Table 2 sulfadoxine: gamma = 3.20, RSE 21%
    pma50_cl_pyra <- 11.9
    label("Pyrimethamine postgestational age at 50% of mature clearance (months)") # Table 2 pyrimethamine: PGA50 = 11.9 mo after conception, RSE 13%
    hill_cl_pyra <- 3.01
    label("Pyrimethamine maturation Hill coefficient (unitless)") # Table 2 pyrimethamine: gamma = 3.01, RSE 46%

    # =====================================================================
    # Nutritional-status effect on relative bioavailability. de Kock 2018
    # Materials and Methods encodes malnutrition as a 'hockey stick':
    # effect = (change in F per unit change in Z-score) x (Z-score + 2),
    # active only below the WHO underweight cut-off of -2. Table 2 reports
    # the coefficient as 'Change in F for each point in Z-score below -2 (%)'
    # = -15.3 for sulfadoxine and -26.7 for pyrimethamine.
    #
    # SIGN. Applied literally to the printed (Z-score + 2) the tabulated
    # negative coefficient would RAISE F in malnourished children. The
    # intended direction is unambiguous from the Abstract ('15.3% and 26.7%
    # lower bioavailabilities ... for each Z-score unit below -2') and the
    # Results ('children who had Z-scores of -3 having 15.3% and 26.7% lower
    # bioavailabilities ... than children with Z-scores of >= -2'), so the
    # covariate is carried as the POSITIVE depth below the knee,
    # max(0, -2 - WAZ), which reproduces the stated direction and magnitude
    # while keeping the tabulated coefficient verbatim.
    #
    # FORM. The LINEAR form (1 + e * depth) is used rather than the
    # compounding (1 + e)^depth that the covariate register prefers as a
    # default, because de Kock 2018 prints the equation and it is linear.
    # The register's preference is explicitly for papers that state a
    # per-unit percentage without printing the equation.
    # =====================================================================
    e_waz_fdepot <- -0.153
    label("Sulfadoxine fractional change in F per Z-score unit below -2 (unitless)") # Table 2 sulfadoxine: change in F per point of Z-score below -2 = -15.3%, RSE 31%
    e_waz_fdepot_pyra <- -0.267
    label("Pyrimethamine fractional change in F per Z-score unit below -2 (unitless)") # Table 2 pyrimethamine: change in F per point of Z-score below -2 = -26.7%, RSE 13%

    # =====================================================================
    # Study effect on pyrimethamine clearance. de Kock 2018 Table 2,
    # 'Difference from clearance in reference 5 (%)'; reference 5 is Bell
    # et al.
    # =====================================================================
    e_study_bell_cl_pyra <- -0.549
    label("Pyrimethamine fractional change in apparent clearance in the Bell et al. study (unitless)") # Table 2 pyrimethamine: difference from clearance in reference 5 = -54.9%, RSE 4%

    # =====================================================================
    # Site-specific multiplicative scaling of the PREDICTED concentrations
    # (de Kock 2018 Table 2, 'Scaling on observations at site(s) (%)'). Each
    # line is ONE estimated parameter shared by the sites its name lists, so
    # the group members are summed into a single 0/1 membership flag inside
    # model() rather than each carrying its own coefficient. The reference
    # group for sulfadoxine is Magude, Bancoumana, Bela Vista, Catuane and
    # Chileka (Table 2 footnote d); for pyrimethamine it is Magude,
    # Mpumalanga, Boane and Chileka (footnote e).
    # =====================================================================
    e_region_mpumalanga_boane_namaacha_cc <- -0.397
    label("Sulfadoxine observation-scaling effect for Mpumalanga, Boane and Namaacha (fraction)") # Table 2 sulfadoxine: scaling on observations, Mpumalanga / Boane / Namaacha = -39.7%, RSE 5%
    e_region_bancoumana_belavista_catuane_cc_pyra <- 0.202
    label("Pyrimethamine observation-scaling effect for Bancoumana, Bela Vista and Catuane (fraction)") # Table 2 pyrimethamine: scaling on observations, Bancoumana / Bela Vista / Catuane = 20.2%, RSE 19%
    e_region_namaacha_cc_pyra <- -0.220
    label("Pyrimethamine observation-scaling effect for Namaacha (fraction)") # Table 2 pyrimethamine: scaling on observations, Namaacha = -22.0%, RSE 23%

    # =====================================================================
    # Between-subject variability. de Kock 2018 Table 2 footnote c: 'BSV
    # values were assumed to be log-normally distributed and are reported
    # here as approximate percent coefficients of variation (CV%)', so the
    # internal log-scale variance is omega^2 = log(1 + CV^2). The two
    # clearances are correlated at 60.0% (Table 2, 'Correlation in CL of the
    # two drugs'); the correlation is taken to be between the random effects
    # on the log scale, following the sibling model
    # deKock_2017_sulfadoxinePyrimethamine.R.
    # =====================================================================

    # Table 2: BSV in CL = 33.9% CV (sulfadoxine, RSE 5%) and 29.0% CV
    # (pyrimethamine, RSE 5%); correlation in CL of the two drugs = 60.0%
    # (RSE 6%). Off-diagonal = rho * sqrt(var_sulfa * var_pyra).
    etalcl + etalcl_pyra ~ c(
      log(1 + 0.339^2),
      0.600 * sqrt(log(1 + 0.339^2) * log(1 + 0.290^2)),
      log(1 + 0.290^2)
    )

    # Table 2 sulfadoxine: BSV in F = 38.4% CV, RSE 12%
    etalfdepot ~ log(1 + 0.384^2)
    # Table 2 sulfadoxine: BSV in ka = 126% CV, RSE 21%
    etalka ~ log(1 + 1.26^2)
    # Table 2 sulfadoxine: BSV in V = 11.2% CV, RSE 23%
    etalvc ~ log(1 + 0.112^2)

    # Table 2 pyrimethamine: BSV in F = 36.1% CV, RSE 4%
    etalfdepot_pyra ~ log(1 + 0.361^2)
    # Table 2 pyrimethamine: BSV in ka = 171% CV, RSE 25%
    etalka_pyra ~ log(1 + 1.71^2)
    # Table 2 pyrimethamine: BSV in V = 15.5% CV, RSE 14%
    etalvc_pyra ~ log(1 + 0.155^2)

    # =====================================================================
    # Residual unexplained variability. de Kock 2018 Materials and Methods:
    # 'A combined error model with both additive and proportional components
    # was used'. Table 2 gives the additive term in the units of each drug's
    # observation.
    # =====================================================================
    addSd <- 3.79
    label("Sulfadoxine additive residual SD (ug/mL whole blood)") # Table 2 sulfadoxine: additive error = 3.79 ug/mL, RSE 5%
    propSd <- 0.171
    label("Sulfadoxine proportional residual SD (fraction)") # Table 2 sulfadoxine: proportional error = 17.1%, RSE 3%
    addSd_pyra <- 6.58
    label("Pyrimethamine additive residual SD (ng/mL whole blood)") # Table 2 pyrimethamine: additive error = 6.58 ng/mL, RSE 5%
    propSd_pyra <- 0.232
    label("Pyrimethamine proportional residual SD (fraction)") # Table 2 pyrimethamine: proportional error = 23.2%, RSE 2%
  })

  model({
    # ------------------------------------------------------------------
    # Allometric size factor. Reference weight 18 kg, the pooled median
    # (de Kock 2018 Table 2 footnote b).
    # ------------------------------------------------------------------
    wt_ratio <- WT / 18

    # ------------------------------------------------------------------
    # Sigmoidal maturation of clearance on postgestational age
    # (Materials and Methods: MAT = PGA^gamma / (PGA^gamma + PGA50^gamma)).
    # ------------------------------------------------------------------
    mat_cl <- PAGE^hill_cl / (PAGE^hill_cl + pma50_cl^hill_cl)
    mat_cl_pyra <- PAGE^hill_cl_pyra / (PAGE^hill_cl_pyra + pma50_cl_pyra^hill_cl_pyra)

    # ------------------------------------------------------------------
    # Nutritional-status 'hockey stick'. waz_below2 is the number of
    # Z-score units below the WHO underweight cut-off of -2 and is zero for
    # adequately nourished patients, so the factor is exactly 1 at
    # WAZ >= -2.
    # ------------------------------------------------------------------
    waz_below2 <- max(0, -2 - WAZ)
    waz_fdepot <- 1 + e_waz_fdepot * waz_below2
    waz_fdepot_pyra <- 1 + e_waz_fdepot_pyra * waz_below2

    # ------------------------------------------------------------------
    # Individual parameters. Sulfadoxine.
    # ------------------------------------------------------------------
    cl <- exp(lcl + etalcl) * wt_ratio^e_wt_cl * mat_cl
    vc <- exp(lvc + etalvc) * wt_ratio^e_wt_vc
    ka <- exp(lka + etalka)

    # ------------------------------------------------------------------
    # Individual parameters. Pyrimethamine, with the Bell et al. study
    # effect on apparent clearance.
    # ------------------------------------------------------------------
    cl_pyra <- exp(lcl_pyra + etalcl_pyra) * wt_ratio^e_wt_cl * mat_cl_pyra *
      (1 + e_study_bell_cl_pyra * STUDY_BELL)
    vc_pyra <- exp(lvc_pyra + etalvc_pyra) * wt_ratio^e_wt_vc
    ka_pyra <- exp(lka_pyra + etalka_pyra)

    # ------------------------------------------------------------------
    # One-compartment disposition with first-order absorption, one chain
    # per drug (Results: 'A one-compartment model with first-order
    # absorption and elimination provided the best fit for both
    # sulfadoxine and pyrimethamine').
    # ------------------------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (cl / vc) * central

    d/dt(depot_pyra) <- -ka_pyra * depot_pyra
    d/dt(central_pyra) <- ka_pyra * depot_pyra - (cl_pyra / vc_pyra) * central_pyra

    # Relative bioavailability is named rather than written inline so that it
    # is returned as an output column: AUC must equal fdepot * dose / cl for a
    # linear model, and that identity is only checkable per subject if the
    # drawn F is visible.
    fdepot <- exp(lfdepot + etalfdepot) * waz_fdepot
    fdepot_pyra <- exp(lfdepot_pyra + etalfdepot_pyra) * waz_fdepot_pyra

    f(depot) <- fdepot
    f(depot_pyra) <- fdepot_pyra

    # ------------------------------------------------------------------
    # Site-specific scaling of the predicted concentration. Each group
    # shares one estimated coefficient, and the member indicators are
    # mutually exclusive, so their sum is a 0/1 group-membership flag.
    # ------------------------------------------------------------------
    region_grp_sulfa <- REGION_MPUMALANGA + REGION_BOANE + REGION_NAMAACHA
    region_grp_pyra <- REGION_BANCOUMANA + REGION_BELAVISTA + REGION_CATUANE

    site_factor <- 1 + e_region_mpumalanga_boane_namaacha_cc * region_grp_sulfa
    site_factor_pyra <- 1 +
      e_region_bancoumana_belavista_catuane_cc_pyra * region_grp_pyra +
      e_region_namaacha_cc_pyra * REGION_NAMAACHA

    # ------------------------------------------------------------------
    # Whole-blood concentrations. Amounts are mg and volumes L, so
    # central / vc is mg/L = ug/mL, the unit de Kock 2018 reports
    # sulfadoxine in. Pyrimethamine is reported in ng/mL, hence the
    # factor of 1000.
    # ------------------------------------------------------------------
    Cc <- (central / vc) * site_factor
    Cc_pyra <- 1000 * (central_pyra / vc_pyra) * site_factor_pyra

    Cc ~ add(addSd) + prop(propSd)
    Cc_pyra ~ add(addSd_pyra) + prop(propSd_pyra)
  })
}
