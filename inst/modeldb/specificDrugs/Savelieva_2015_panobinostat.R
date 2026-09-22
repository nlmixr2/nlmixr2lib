Savelieva_2015_panobinostat <- function() {
  description <- "Three-compartment population PK model with first-order absorption and linear elimination for the pan-deacetylase inhibitor panobinostat in adults with advanced hematologic and solid tumors (Savelieva 2015, FIRST final model). Intravenous and oral data from 14 phase 1 and phase 2 studies were fitted simultaneously, so clearance and volumes are absolute rather than apparent and the absolute oral bioavailability is identified. Distribution is parameterized directly in first-order rate constants K23, K32, K24 and K42 as the authors did, not in intercompartmental clearances. Clearance and central volume carry power effects of body surface area and age plus multiplicative race factors; the absorption rate constant differs between the clinical service formulation and the final market image. A companion re-parameterization of the same data set, using weight-based allometry, intercompartmental clearances and an absorption lag, is Savelieva_2015_panobinostat_allometric."
  reference <- "Savelieva M, Woo MM, Schran H, Mu S, Nedelman J, Capdeville R. Population pharmacokinetics of intravenous and oral panobinostat in patients with hematologic and solid tumors. Eur J Clin Pharmacol. 2015;71(6):663-672. doi:10.1007/s00228-015-1846-7. Parameter estimates from Supplementary Table S2b; model code from Supplementary Table S2a."
  vignette <- "Savelieva_2015_panobinostat"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "panobinostat", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "panobinostat", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "panobinostat", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "panobinostat", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    BSA = list(
      description = "Baseline body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects on both CL and V2, centered on 1.9 m^2, the population median (Savelieva 2015 Supplementary Table S2a: '(BSA0/1.9)**THETA(10)' and '(BSA0/1.9)**THETA(11)'). Computed by the authors with the Gehan-George formula, printed in Methods, Construction of the population pharmacokinetics data set as BSA = 234.94 * (Weight^0.515 * Height^0.422) / 10000 with weight in kg and height in cm. That formula applied to the reported median weight 76.4 kg and median height 170 cm returns 1.915 m^2, which reproduces the stated median BSA of 1.9 m^2, so the centering value and the formula are mutually consistent. Observed quartiles 1.8 and 2.1 m^2 (Table 3). Body weight and BMI were also screened; BSA was the body-size covariate retained in this first final model, and the second final model replaces it with weight.",
      source_name = "BSA0"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects on both CL and V2, centered on 61 years, the population median (Savelieva 2015 Supplementary Table S2a: '(AGE0/61)**THETA(12)' and '(AGE0/61)**THETA(13)'). Observed range 16-88 years, quartiles 51 and 70 years (Table 3). The CL exponent is POSITIVE, so clearance rises and exposure falls with age; the paper quantifies this as AUC falling from 102 to 95 ng*h/mL across the age quartiles at median BSA, an effect it judges small relative to the 74 percent interindividual variability in clearance.",
      source_name = "AGE0"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Caucasian, when RACE_BLACK and RACE_OTHER are also 0)",
      notes = "Multiplicative factor applied as THETA^indicator on CL and on V2 (Savelieva 2015 Supplementary Table S2a: 'IF (RACE .EQ. 3) AS=1' then 'THETA(14)**(AS)'). The canonical 1 = Asian orientation matches the source coding, so no value flip is needed. The three race indicators are mutually exclusive; all three equal to 0 selects the Caucasian reference. 27 of the 581 patients were Asian (Table 1).",
      source_name = "AS (derived from RACE == 3)"
    ),
    RACE_BLACK = list(
      description = "Black race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Caucasian, when RACE_ASIAN and RACE_OTHER are also 0)",
      notes = "Multiplicative factor applied as THETA^indicator on CL and on V2 (Savelieva 2015 Supplementary Table S2a: 'IF (RACE .EQ. 2) BL=1' then 'THETA(16)**(BL)'). 34 of the 581 patients were Black (Table 1). The clearance factor is 1.010, i.e. essentially no effect on CL; the whole Black-vs-Caucasian exposure difference in Table 3 comes through the 1.241 factor on V2.",
      source_name = "BL (derived from RACE == 2)"
    ),
    RACE_OTHER = list(
      description = "Race category 'other' indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Caucasian, when RACE_ASIAN and RACE_BLACK are also 0)",
      notes = "Multiplicative factor applied as THETA^indicator on CL and on V2 (Savelieva 2015 Supplementary Table S2a: 'IF (RACE .EQ. 88) OT=1' then 'THETA(18)**(OT)'). 24 of the 581 patients fell in this category (Table 1). It carries the largest clearance effect of the three race indicators, a factor of 0.719 on CL, and Table 3 accordingly gives this group the highest typical exposure of the four race categories.",
      source_name = "OT (derived from RACE == 88)"
    ),
    FORM_PANO_CSF = list(
      description = "Panobinostat clinical service formulation versus final market image indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (final market image, FMI, the formulation intended for commercialization)",
      notes = "Selects the absorption rate constant between the two estimated values (Savelieva 2015 Supplementary Table S2a: 'KA=THETA(7)*(1-FORM)+THETA(8)*FORM', with THETA(7) labelled KA.FMI and THETA(8) labelled KA.CSF in Supplementary Table S2b). 1 = clinical service formulation (CSF), used in studies B2101, B2102 and B1101; 0 = final market image (FMI), used in every other oral study. 106 of the 494 orally dosed patients received the CSF. Formulation affected Ka but NOT bioavailability: the paper states explicitly that formulation 'had a significant impact on the Ka (but not on the bioavailability) of oral panobinostat', and the control stream applies a single THETA(9) to both formulations. This covariate is irrelevant to intravenous dosing, which bypasses the depot entirely.",
      source_name = "FORM"
    )
  )

  # Screened by the authors but NOT retained in either final model. Documented so
  # the provenance of the covariate search survives, without raising a
  # declared-but-unreferenced convention warning.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened as a body-size covariate on CL and V2 and not retained in THIS model, which uses BSA instead; median 76.4 kg, range 41-196.4 kg (Savelieva 2015 Results, Patients). Weight IS the body-size covariate of the second final model - see Savelieva_2015_panobinostat_allometric."
    ),
    HT = list(
      description = "Baseline height",
      units = "cm",
      type = "continuous",
      notes = "Collected at screening and used only as an input to the Gehan-George BSA formula, never entered as a covariate in its own right; median 170 cm, range 143-198 cm (Savelieva 2015 Results, Patients). Missing for 35 patients, who were assigned the population median."
    ),
    BMI = list(
      description = "Baseline body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened as a body-size covariate on CL and V2 and not retained (Savelieva 2015 Methods, Analysis of the effects of covariates). Missing for 35 patients, who were assigned the population median."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL and V2 and not retained; the paper reports that covariate analysis 'showed no impact on panobinostat clearance and volume by patients' sex'. 219 of 581 patients were female (Table 1)."
    ),
    CRCL = list(
      description = "Baseline creatinine clearance",
      units = "L/h",
      type = "continuous",
      notes = "Screened on CL and V2 and not retained; kidney function had no statistically significant effect. Computed by the authors with the Cockcroft-Gault formula from serum creatinine in micromoles per litre (Savelieva 2015 Methods, Construction of the population pharmacokinetics data set). Missing for 13 patients, who were assigned the population median."
    ),
    HEPIMP_MILD = list(
      description = "Mild hepatic impairment indicator",
      units = "(binary)",
      type = "binary",
      notes = "Part of the four-level liver-status covariate (normal 483, mild 91, moderate 6, severe 1; Table 1) graded on total bilirubin and AST against the upper limit of normal. Screened on CL and V2 and not retained. The paper cautions in its Discussion that this null result likely reflects the trial eligibility criteria, which generally required baseline bilirubin at or below 1.5x ULN and AST/ALT at or below 2x ULN, and that dedicated organ-impairment studies DID find a significant exposure increase with hepatic impairment."
    ),
    HEPIMP_MOD = list(
      description = "Moderate hepatic impairment indicator",
      units = "(binary)",
      type = "binary",
      notes = "Part of the same four-level liver-status covariate; only 6 of 581 patients, so the effect was not estimable with useful precision. Not retained. The US prescribing information nonetheless recommends a reduced starting dose in moderate hepatic impairment (Savelieva 2015 Discussion)."
    ),
    HEPIMP_SEV = list(
      description = "Severe hepatic impairment indicator",
      units = "(binary)",
      type = "binary",
      notes = "Part of the same four-level liver-status covariate; a single patient of 581. Not retained."
    ),
    TUMTP_OTHER = list(
      description = "Tumor-type indicator",
      units = "(binary)",
      type = "binary",
      notes = "Tumor type was screened on CL and V2 across the pooled hematologic and solid-tumor population and not retained (Savelieva 2015 Methods, Analysis of the effects of covariates, and Abstract). Placeholder entry standing for the paper's tumor-type covariate as a whole; the source does not publish the individual category codes it tested."
    ),
    CONMED_AZOLE = list(
      description = "Strong CYP3A4/5 inhibitor comedication indicator",
      units = "(binary)",
      type = "binary",
      notes = "One of five comedication groups screened on CL and not retained: drugs known to prolong QT, CYP2D6 substrates, strong CYP3A4/5 inhibitors, CYP3A4 inducers, and sensitive CYP3A4 substrates (Savelieva 2015 Methods, Analysis of the effects of covariates). Placeholder entry standing for that comedication screen. The Discussion attributes the null result to protocol guidance that recommended avoiding strong CYP3A4 inhibitors, and notes that a dedicated ketoconazole interaction study DID show a significant exposure increase."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 581,
    n_studies = 14,
    n_observations = 7834,
    age_range = "16-88 years",
    age_median = "61 years (quartiles 51 and 70 years)",
    weight_range = "41-196.4 kg",
    weight_median = "76.4 kg",
    height_range = "143-198 cm",
    height_median = "170 cm",
    bsa_median = "1.9 m^2 (quartiles 1.8 and 2.1 m^2)",
    sex_female_pct = 37.7,
    race_ethnicity = c(Caucasian = 85.4, Black = 5.9, Asian = 4.6, Other = 4.1),
    disease_state = "Advanced hematologic and solid tumors, including cutaneous T-cell lymphoma, chronic myeloid leukemia, multiple myeloma, Hodgkin lymphoma, non-Hodgkin lymphoma and advanced solid tumors",
    hepatic_function = "Liver status graded on total bilirubin and AST against the upper limit of normal: normal 483, mild 91, moderate 6, severe 1 (Table 1)",
    dose_range = "Intravenous 1.2-20 mg/m^2 daily under various intermittent regimens (studies A2101 and A2102, 87 patients); oral 10-80 mg/day in the phase 1 dose-escalation studies and 20-45 mg in the phase 2 and clinical pharmacology studies, most commonly 20 mg on days 1, 3 and 5 of each week (494 patients)",
    formulation = "Clinical service formulation (CSF) in oral studies B2101, B2102 and B1101 (106 patients); final market image (FMI) in every other oral study (388 patients)",
    regions = "International; study B1101 enrolled 13 Japanese patients and B1201 was conducted in Japan",
    notes = "Pooled from 14 open-label phase 1 and phase 2 studies listed in Supplementary Table S1 (A2101, A2102, B1101, B1201, B2101, B2102, B2109, B2110, B2111, B2201, B2202, B2203, B2211 and E2214). Baseline demographics: Table 1 and Fig. 1. The bioanalytical assay was linear from 0.5 to 500 ng/mL; the lower limit of quantification was 0.5 ng/mL in all studies except B2201 and B2203, where it was 0.1 ng/mL. Values below the limit of quantification, 6 percent of the total, were excluded. Missing baseline height (35 patients), BMI (35), creatinine clearance (13) and body weight (13) were imputed with population medians."
  )

  ini({
    # Structural parameters - typical values for the reference patient: a
    # Caucasian of BSA 1.9 m^2 and age 61 years, i.e. every covariate term equal
    # to 1. Values are from Savelieva 2015 Supplementary Table S2b; the
    # Interpretive Name column of that table supplies the units. Because the
    # intravenous and oral data were fitted jointly, CL and V2 are absolute, not
    # apparent.
    lcl  <- log(33.085); label("Clearance (L/h)")                                              # Suppl Table S2b, Theta 1 'CL (L/h)' = 33.085, PctSE 6.70, bootstrap Q10-Q90 30.467-36.286
    lvc  <- log(24.838); label("Central volume of distribution V2 (L)")                        # Suppl Table S2b, Theta 2 'V2 (L)' = 24.838, PctSE 9.84, bootstrap Q10-Q90 21.843-27.863

    # Distribution is parameterized directly in first-order rate constants, as
    # the authors did (Suppl Table S2a: 'K23=THETA(3)' ... 'K42=THETA(6)'), NOT
    # in intercompartmental clearances. The NONMEM compartment numbering is
    # 1 = depot, 2 = central, 3 = peripheral1, 4 = peripheral2, so K23 / K32 map
    # to the canonical k12 / k21 and K24 / K42 map to k13 / k31.
    lk12 <- log(1.810);  label("Rate constant, central to peripheral1, K23 (1/h)")             # Suppl Table S2b, Theta 3 'K23 (1/h)' = 1.810, PctSE 10.58, bootstrap Q10-Q90 1.584-2.070
    lk21 <- log(0.507);  label("Rate constant, peripheral1 to central, K32 (1/h)")             # Suppl Table S2b, Theta 4 'K32 (1/h)' = 0.507, PctSE 7.63, bootstrap Q10-Q90 0.464-0.559
    lk13 <- log(1.424);  label("Rate constant, central to peripheral2, K24 (1/h)")             # Suppl Table S2b, Theta 5 'K24 (1/h)' = 1.424, PctSE 9.37, bootstrap Q10-Q90 1.263-1.610
    lk31 <- log(0.040);  label("Rate constant, peripheral2 to central, K42 (1/h)")             # Suppl Table S2b, Theta 6 'K42 (1/h)' = 0.040, PctSE 4.75, bootstrap Q10-Q90 0.037-0.042

    # Absorption. Two separately estimated rate constants, selected by
    # formulation; the source switches between them linearly in the FORM
    # indicator rather than through a ratio.
    lka_fmi <- log(0.321); label("Absorption rate constant, final market image formulation (1/h)")     # Suppl Table S2b, Theta 7 'KA.FMI (1/h)' = 0.321, PctSE 8.49, bootstrap Q10-Q90 0.285-0.354
    lka_csf <- log(0.544); label("Absorption rate constant, clinical service formulation (1/h)")       # Suppl Table S2b, Theta 8 'KA.CSF (1/h)' = 0.544, PctSE 4.53, bootstrap Q10-Q90 0.515-0.578

    # Absolute oral bioavailability, shared by both oral formulations. It is
    # identified because intravenous and oral data were fitted together.
    lfdepot <- log(0.214); label("Absolute oral bioavailability (fraction)")                    # Suppl Table S2b, Theta 9 'F1' = 0.214, PctSE 6.90, bootstrap Q10-Q90 0.195-0.235; quoted in the Abstract as 21.4 percent

    # Covariate effects. Continuous covariates enter as a power of the
    # median-normalized value; race enters as a multiplicative factor raised to
    # its 0-or-1 indicator, so the factor applies when the indicator is 1 and
    # drops out otherwise (Suppl Table S2a).
    e_bsa_cl <- 1.002; label("Power exponent on BSA relative to 1.9 m^2 for clearance (unitless)")                          # Suppl Table S2b, Theta 10 'CL.BSA' = 1.002, PctSE 22.80, bootstrap Q10-Q90 0.724-1.319
    e_bsa_vc <- 1.359; label("Power exponent on BSA relative to 1.9 m^2 for central volume (unitless)")                     # Suppl Table S2b, Theta 11 'V2.BSA' = 1.359, PctSE 14.50, bootstrap Q10-Q90 1.068-1.576
    e_age_cl <- 0.176; label("Power exponent on age relative to 61 years for clearance (unitless)")                         # Suppl Table S2b, Theta 12 'CL.AGE' = 0.176, PctSE 56.05, bootstrap Q10-Q90 0.047-0.280
    e_age_vc <- 0.396; label("Power exponent on age relative to 61 years for central volume (unitless)")                    # Suppl Table S2b, Theta 13 'V2.AGE' = 0.396, PctSE 23.56, bootstrap Q10-Q90 0.270-0.499

    e_race_asian_cl <- 1.171; label("Multiplicative factor on clearance for Asian relative to Caucasian patients (unitless)")        # Suppl Table S2b, Theta 14 'CL.ASIAN' = 1.171, PctSE 9.16, bootstrap Q10-Q90 1.047-1.313
    e_race_asian_vc <- 1.373; label("Multiplicative factor on central volume for Asian relative to Caucasian patients (unitless)")   # Suppl Table S2b, Theta 15 'V2.ASIAN' = 1.373, PctSE 12.74, bootstrap Q10-Q90 1.169-1.609
    e_race_black_cl <- 1.010; label("Multiplicative factor on clearance for Black relative to Caucasian patients (unitless)")        # Suppl Table S2b, Theta 16 'CL.BLACK' = 1.010, PctSE 14.06, bootstrap Q10-Q90 0.854-1.195
    e_race_black_vc <- 1.241; label("Multiplicative factor on central volume for Black relative to Caucasian patients (unitless)")   # Suppl Table S2b, Theta 17 'V2.BLACK' = 1.241, PctSE 15.60, bootstrap Q10-Q90 1.017-1.490
    e_race_other_cl <- 0.719; label("Multiplicative factor on clearance for race category 'other' relative to Caucasian (unitless)")      # Suppl Table S2b, Theta 18 'CL.OTHER' = 0.719, PctSE 18.58, bootstrap Q10-Q90 0.559-0.908
    e_race_other_vc <- 1.127; label("Multiplicative factor on central volume for race category 'other' relative to Caucasian (unitless)") # Suppl Table S2b, Theta 19 'V2.OTHER' = 1.127, PctSE 12.69, bootstrap Q10-Q90 0.954-1.324

    # Interindividual variability. Savelieva 2015 Methods states that random
    # interindividual variability was estimated with an exponential
    # parameterization, and that eta_CL and eta_V were modelled as normally
    # distributed with mean zero and a full 2x2 covariance matrix. The tabulated
    # OMEGA entries are therefore VARIANCES on the log scale, which the paper's
    # own headline number confirms: sqrt(exp(0.439) - 1) = 0.742, the 74 percent
    # interindividual variability in clearance quoted in the Abstract.
    etalcl + etalvc ~ c(0.439,
                        0.178, 0.334)  # Suppl Table S2b, Omega 1,1 'OM.CL' = 0.439; 2,1 'OM.CLV2' = 0.178; 2,2 'OM.V2' = 0.334; CV 74 percent and 63 percent respectively

    # Residual error. Suppl Table S2a: 'Y=F*(1+EPS(1))+EPS(2)', a combined
    # proportional plus additive model. Suppl Table S2b reports the SIGMA
    # entries as variances ('VAR.PROP', 'VAR.ADD'), so each standard deviation
    # is the square root of the tabulated value.
    propSd <- 0.49193; label("Proportional residual error (fraction)")       # Suppl Table S2b, Sigma 1 'VAR.PROP' = 0.242, PctSE 4.45; sqrt(0.242)
    addSd  <- 0.11402; label("Additive residual error (ng/mL)")              # Suppl Table S2b, Sigma 2 'VAR.ADD' = 0.013, PctSE 80.33; sqrt(0.013)
  })

  model({
    # Covariate multipliers, transcribed from Savelieva 2015 Supplementary
    # Table S2a:
    #   CL1 = THETA(1)*EXP(ETA(1))*(BSA0/1.9)**THETA(10)
    #   CL  = CL1*(AGE0/61)**THETA(12)*THETA(14)**(AS)*THETA(16)**(BL)*THETA(18)**(OT)
    #   V21 = THETA(2)*EXP(ETA(2))*(BSA0/1.9)**THETA(11)
    #   V2  = V21*(AGE0/61)**THETA(13)*THETA(15)**(AS)*THETA(17)**(BL)*THETA(19)**(OT)
    # The race factors are raised to their 0-or-1 indicators, so each equals 1
    # for a Caucasian patient and the reference category needs no parameter.
    bsa_cl  <- (BSA / 1.9)^e_bsa_cl
    bsa_vc  <- (BSA / 1.9)^e_bsa_vc
    age_cl  <- (AGE / 61)^e_age_cl
    age_vc  <- (AGE / 61)^e_age_vc
    race_cl <- e_race_asian_cl^RACE_ASIAN * e_race_black_cl^RACE_BLACK * e_race_other_cl^RACE_OTHER
    race_vc <- e_race_asian_vc^RACE_ASIAN * e_race_black_vc^RACE_BLACK * e_race_other_vc^RACE_OTHER

    cl <- exp(lcl + etalcl) * bsa_cl * age_cl * race_cl
    vc <- exp(lvc + etalvc) * bsa_vc * age_vc * race_vc

    # Formulation-dependent absorption, written as the source's linear switch
    # 'KA = THETA(7)*(1-FORM) + THETA(8)*FORM'. The source additionally forces
    # KA = 0 and F1 = 1 for intravenous records; here that is achieved instead by
    # routing intravenous doses straight into central, which never touches the
    # depot, ka or f(depot).
    ka     <- exp(lka_fmi) * (1 - FORM_PANO_CSF) + exp(lka_csf) * FORM_PANO_CSF
    fdepot <- exp(lfdepot)

    # Rate constants. Elimination is first order from central; the four
    # distribution rate constants are estimated directly rather than derived
    # from intercompartmental clearances.
    kel <- cl / vc
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    k13 <- exp(lk13)
    k31 <- exp(lk31)

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central -
                          k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    f(depot) <- fdepot

    # Doses are in mg and volumes in L, so central/vc is mg/L; the factor 1000
    # converts to the ng/mL scale the paper reports and on which the additive
    # residual error is expressed. This reproduces the source's scaling
    # 'S2=V2/1000' (Suppl Table S2a), under which NONMEM computes A2/(V2/1000).
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
