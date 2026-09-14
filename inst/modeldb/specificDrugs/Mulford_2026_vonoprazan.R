Mulford_2026_vonoprazan <- function() {
  description <- paste0(
    "Two-compartment population PK model for the potassium-competitive ",
    "acid blocker vonoprazan, with delayed oral absorption through three ",
    "transit compartments parameterized by a mean transit time (all three ",
    "transit transfer rates fixed at 4/MTT) followed by first-order ",
    "absorption into the central compartment. Fitted to the pooled adult, ",
    "adolescent and child dataset (Mulford 2026 Model 3: 392 subjects, ",
    "8201 quantifiable plasma concentrations from eight adult studies in ",
    "healthy volunteers aged 18 to 54 years, one adolescent study ",
    "(VPED-102, NCT05343364, ages 12 to 17) and one child study ",
    "(VPED-103, NCT06106022, ages 6 to under 12) in patients with ",
    "gastroesophageal reflux disease). Relative bioavailability is ",
    "anchored at 1 (no intravenous data) and carries an over-proportional ",
    "dose-power effect and a female effect; absorption rate carries ",
    "female, age-power and day-after-first-dose effects; clearance and ",
    "central volume each carry a day-after-first-dose effect, and central ",
    "volume a body-weight power effect. Between-subject variability is ",
    "estimated on clearance and central volume (correlated), on mean ",
    "transit time and on absorption rate; residual variability is ",
    "combined proportional plus a fixed, negligible additive term. ",
    "Notably, Mulford 2026 found NO age or weight effect on clearance ",
    "across the 6-to-54-year range, so fixed 10 mg and 20 mg doses give ",
    "adults, adolescents and children comparable steady-state exposure. ",
    "For the organ-maturation extension used to project exposure in ",
    "infants and children aged 1 month to under 6 years, see ",
    "modellib('Mulford_2026_vonoprazan_maturation')."
  )
  reference <- paste(
    "Mulford DJ, Facius A, Witt G, Howden CW, Wagner T, Leifke E,",
    "Scarpignato C. Development and Use of a Population Pharmacokinetic",
    "Model for Characterizing the Pharmacokinetics of Vonoprazan in",
    "Pediatric Patients.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15:e70291.",
    "doi:10.1002/psp4.70291.",
    "Structural detail and covariate centering values taken from the",
    "Supporting Information Model Code (the final-model NONMEM control",
    "stream, PSP4-15-e70291-s001.docx).",
    sep = " "
  )
  vignette <- "Mulford_2026_vonoprazan"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against the Supporting Information Model
  # Code $MODEL block (GUT, TR1, TR2, TR3, CENTRAL, PERIPHERAL), which
  # doses AMT in mg into GUT and scales the central compartment to plasma
  # concentration with S5 = VC/1000.
  compartmentData <- list(
    depot       = list(analyte = "vonoprazan", units = "mg", specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "vonoprazan", units = "mg", specimen = "administration site", verified = TRUE),
    transit2    = list(analyte = "vonoprazan", units = "mg", specimen = "administration site", verified = TRUE),
    transit3    = list(analyte = "vonoprazan", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "vonoprazan", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vonoprazan", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters as a power model on central volume only, centered at 70 kg",
        "(Mulford 2026 Methods 2.2.2, and the Supporting Information Model",
        "Code line TVVC = TVVC * (WEIGHT/70)**THETA(11)). Mulford 2026",
        "deliberately did NOT carry a weight effect on clearance: Results",
        "3.3.3 reports that the dose-normalized steady-state AUC versus",
        "weight slopes had wide prediction intervals in both the adolescent",
        "and the child cohort and pointed in opposite directions, so no",
        "weight effect on clearance was added. Cohort weights span 20.7 to",
        "132 kg (Table S2); median 66 kg in adults, 64 kg in adolescents",
        "and 32 kg in children (Figure S2)."
      ),
      source_name        = "WEIGHT"
    ),
    AGE = list(
      description        = "Chronological age.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters as a power model on the absorption rate constant only,",
        "centered at 28 years per the Supporting Information Model Code",
        "line TVKA = TVKA * (AGE/28)**THETA(14). 28 years is essentially",
        "the adult cohort mean age of 28.2 years (Table S2). NOTE: Mulford",
        "2026 Methods 2.2.2 states the age centering value as 18 years,",
        "which contradicts the deposited control stream; the control-stream",
        "value is used here because it is the authors' executable code and",
        "because the same sentence's other two centering values (20 mg for",
        "dose, 70 kg for weight) match the control stream exactly. See the",
        "vignette Errata. As with weight, no age effect on clearance was",
        "retained (Results 3.3.3); the age effect is on absorption only,",
        "with a negative exponent, i.e. slower absorption in younger",
        "subjects. Cohort ages span 6 to 54 years."
      ),
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "The source column is SEX with SEX == 2 denoting female",
        "(Supporting Information Model Code: IF( SEX == 2 ) ... for both",
        "the relative-bioavailability and the absorption-rate effect), so",
        "SEXF = as.integer(SEX == 2). Carries two effects in opposite",
        "directions: females absorb vonoprazan more slowly but absorb more",
        "of it."
      ),
      source_name        = "SEX"
    ),
    DAY2 = list(
      description        = "Day-after-first-dose landmark indicator: 1 = the record falls on study day 2 or later, 0 = the record falls on study day 1 (the first dosing day).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (study day 1, the first dosing day)",
      notes              = paste(
        "Implements the Supporting Information Model Code branches",
        "IF(DAY>1) on absorption rate, clearance and central volume, i.e.",
        "DAY2 = as.integer(study_day >= 2). Mulford 2026 Results 3.1 is",
        "explicit that the day-1 / later-day step is a design artifact",
        "rather than a resolved time course: plasma concentrations were",
        "mainly observed on Day 1 and again on Day 7, so a continuous",
        "relationship could not be identified and a binary change for",
        "day > 1 was used instead. All three effects are negative, so",
        "day-2-and-later typical absorption rate, clearance and central",
        "volume are each roughly 11 to 21 percent below their day-1",
        "values. For a single-dose simulation set DAY2 = 0 throughout; for",
        "a repeated-dose simulation set DAY2 = 0 over the first 24 h and 1",
        "thereafter."
      ),
      source_name        = "DAY"
    ),
    DOSE_VONOPRAZAN_MG = list(
      description        = "Administered vonoprazan dose level carried on every record of the dosing interval, in mg.",
      units              = "mg",
      type               = "continuous",
      reference_category = "20 mg (the dose at which relative bioavailability equals its typical value)",
      notes              = paste(
        "Continuous power-model regressor on relative bioavailability,",
        "centered at 20 mg (Mulford 2026 Methods 2.2.2; Supporting",
        "Information Model Code line",
        "TVFREL = TVFREL * (DOSE/20)**THETA(10)). The positive exponent",
        "encodes the slightly over-proportional exposure the authors",
        "describe in Results 3.1. The source studies span 1 mg to 120 mg",
        "(Table S1), so the power model is supported over that range only.",
        "This column is the dose LEVEL and is distinct from the record's",
        "AMT: a subject on 20 mg once daily carries 20 on every record,",
        "including observation records."
      ),
      source_name        = "DOSE"
    )
  )

  covariatesDataExcluded <- list(
    RACE = list(
      description        = "Self-reported race category (Asian, Black, Other, White).",
      units              = "(categorical)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Carried in the analysis dataset (Supporting Information Model",
        "Code $INPUT includes RACE) and screened graphically in Figure S3,",
        "but not retained: the figure legend states race was formally",
        "tested during prior adult population pharmacokinetic model",
        "development and was neither statistically significant nor",
        "clinically meaningful. No point estimate is reported, so no",
        "effect can be reproduced."
      ),
      source_name        = "RACE"
    ),
    EGFR = list(
      description        = "Estimated glomerular filtration rate.",
      units              = "mL/min/1.73m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Present in the analysis dataset (Supporting Information Model",
        "Code $INPUT includes EGFR) but not tested or retained in the",
        "published covariate model, which screened only dose, time, sex,",
        "weight and age (Methods 2.2.2). Vonoprazan is extensively",
        "metabolized hepatically with only a small unchanged renal",
        "fraction, so a renal-function effect would not be expected. No",
        "point estimate is reported."
      ),
      source_name        = "EGFR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 392L,
    n_studies      = 10L,
    n_observations = 8201L,
    age_range      = "6-54 years",
    weight_range   = "20.7-132 kg",
    weight_median  = "66 kg (adults), 64 kg (adolescents), 32 kg (children)",
    disease_state  = "healthy volunteers (all eight adult studies) and patients with gastroesophageal reflux disease (the adolescent study VPED-102 and the child study VPED-103)",
    dose_range     = "1-120 mg oral; single dose, 10-40 mg once daily, and 20 mg twice daily",
    regions        = "Japan, Europe, China, United States",
    notes          = paste(
      "Baseline demographics from Mulford 2026 Table S2; per-study design,",
      "region, dose range and subject/observation counts from Table S1;",
      "pediatric study design from Table S3. Adults n = 354 (mean age 28.2",
      "years, SD 7.46; mean weight 68.3 kg, SD 10.7); adolescents n = 17",
      "(mean age 14.6 years, SD 1.69; mean weight 72.0 kg, SD 25.6,",
      "inflated by one 15-year-old weighing 132 kg); children n = 21 (mean",
      "age 8.76 years, SD 1.76; mean weight 36.2 kg, SD 12.4). The",
      "female fraction and the race distribution are not reported",
      "numerically; Figure S3 stratifies random effects by Asian, Black,",
      "Other and White. About 5 percent of concentrations were below the",
      "limit of quantification and were handled by the Beal M3 method",
      "(Table S4 reports the sensitivity analysis without M3)."
    )
  )

  ini({
    # Structural parameters. Mulford 2026 Table 1, rightmost (Children) column
    # = Model 3, the final model fitted to the pooled adult + adolescent +
    # child dataset; Results 3.1 directs the reader to this column for the
    # numerical estimates. Relative standard errors are quoted per parameter.

    # Relative bioavailability was anchored at 1 because no intravenous data
    # were available (Results 3.1 and $THETA row 1: '100 FIX').
    lfdepot <- fixed(log(1)); label("Typical relative oral bioavailability at a 20 mg dose in males (unitless)")

    lmtt <- log(0.762); label("Typical mean transit time of the absorption delay chain (h)")           # Table 1 Children, MTT TV, RSE 3.3
    lka  <- log(3.08);  label("Typical first-order absorption rate constant from the third transit compartment (1/h)") # Table 1 Children, ka TV, RSE 7.7
    lcl  <- log(118);   label("Typical apparent elimination clearance (L/h)")                          # Table 1 Children, CL TV, RSE 2.2
    lvc  <- log(751);   label("Typical apparent central volume of distribution (L)")                   # Table 1 Children, Vc TV, RSE 2.8
    lq   <- log(49.8);  label("Typical apparent distribution clearance (L/h)")                         # Table 1 Children, Q TV, RSE 4.8
    lvp  <- log(271);   label("Typical apparent peripheral volume of distribution (L)")                # Table 1 Children, Vp TV, RSE 2.2

    # Covariate effects. The two power exponents are applied as
    # (covariate/centre)^exponent; the three percent-change effects are applied
    # as (1 + effect * indicator) with the effect stored as a fraction.
    e_dose_fdepot <- 0.290;  label("Power exponent on DOSE_VONOPRAZAN_MG/20 for relative bioavailability (unitless)") # Table 1 Children, Frel dose-effect, RSE 4.8
    e_sexf_fdepot <- 0.351;  label("Fractional change in relative bioavailability for females (unitless)")            # Table 1 Children, Frel female-effect 35.1, RSE 10.8
    e_sexf_ka     <- -0.490; label("Fractional change in absorption rate constant for females (unitless)")            # Table 1 Children, ka female-effect -49.0, RSE 12.9
    e_age_ka      <- -0.693; label("Power exponent on AGE/28 for the absorption rate constant (unitless)")            # Table 1 Children, ka age-effect, RSE 10.6
    e_day2_ka     <- -0.208; label("Fractional change in absorption rate constant on study day 2 and later (unitless)") # Table 1 Children, ka Day >1-effect -20.8, RSE 9.9
    e_day2_cl     <- -0.115; label("Fractional change in elimination clearance on study day 2 and later (unitless)")    # Table 1 Children, CL Day >1-effect -11.5, RSE 3.5
    e_wt_vc       <- 0.668;  label("Power exponent on WT/70 for the central volume of distribution (unitless)")       # Table 1 Children, Vc weight-effect, RSE 6.6
    e_day2_vc     <- -0.116; label("Fractional change in central volume on study day 2 and later (unitless)")         # Table 1 Children, Vc Day >1-effect -11.6, RSE 11.4

    # Between-subject variability. Table 1 footnote a gives the BSV rows as
    # standard deviations and correlations 'as reported by NONMEM'; footnote b
    # gives the same quantity as a coefficient of variation derived by
    # 100 * sqrt(exp(omega^2) - 1). Both readings are printed, which pins the
    # scale: 0.526 as an SD gives 56.5 percent CV against the printed 56.4,
    # whereas reading 0.526 as a variance would give 84 percent. The
    # Supporting Information Model Code $OMEGA block confirms it directly for
    # the Model 2 fit, where 0.542^2 = 0.2939 matches the deposited 0.2935.
    # Variances below are the squares of the Table 1 Children BSV a entries;
    # the covariance is 0.913 * 0.366 * 0.389.
    etalcl + etalvc ~ c(0.133956,
                        0.129987, 0.151321)                                                            # Table 1 Children, CL BSV a 0.366, Cor 0.913, Vc BSV a 0.389
    etalmtt ~ 0.276676                                                                                 # Table 1 Children, MTT BSV a 0.526, RSE 4.1
    etalka  ~ 0.413449                                                                                 # Table 1 Children, ka BSV a 0.643, RSE 7.7

    # Residual variability. The Supporting Information Model Code $ERROR block
    # is W = sqrt((THETA(8)/100 * IPRED)^2 + THETA(9)^2) with $SIGMA 1 FIX,
    # i.e. a combined proportional-plus-additive standard deviation, which is
    # exactly what add() + prop() forms in nlmixr2. The additive term was
    # fixed at a negligible value ($THETA row 9: '0.001 FIX') so that the
    # residual is effectively proportional away from the limit of
    # quantification.
    propSd <- 0.247;        label("Proportional residual error (fraction)")   # Table 1 Children, Residual variability Prop 24.7, RSE 0.4
    addSd  <- fixed(0.001); label("Additive residual error (ng/mL)")          # Table 1 Children, Residual variability Add, reported as fixed
  })

  model({
    # 1. Derived covariate terms. Centering values: 20 mg for dose, 70 kg for
    #    weight (Methods 2.2.2, confirmed by the control stream) and 28 years
    #    for age (control stream line TVKA = TVKA * (AGE/28)**THETA(14); the
    #    Methods text's '18 years' is contradicted by the deposited code).
    dose_fdepot <- (DOSE_VONOPRAZAN_MG / 20)^e_dose_fdepot
    age_ka      <- (AGE / 28)^e_age_ka
    wt_vc       <- (WT / 70)^e_wt_vc

    # 2. Individual parameters. Relative bioavailability, distribution
    #    clearance and peripheral volume carry no between-subject variability
    #    (their $OMEGA diagonals are 0 FIX in the control stream).
    frel <- exp(lfdepot) * dose_fdepot * (1 + e_sexf_fdepot * SEXF)
    mtt  <- exp(lmtt + etalmtt)
    ka   <- exp(lka + etalka) * (1 + e_sexf_ka * SEXF) * age_ka * (1 + e_day2_ka * DAY2)
    cl   <- exp(lcl + etalcl) * (1 + e_day2_cl * DAY2)
    vc   <- exp(lvc + etalvc) * wt_vc * (1 + e_day2_vc * DAY2)
    q    <- exp(lq)
    vp   <- exp(lvp)

    # 3. Micro-constants. The chain carries three transit transfers at ktr and
    #    a fourth, terminal transfer into the central compartment at ka
    #    (Results 3.1: three delay compartments, all three transfer rates
    #    derived as 4/MTT, then a first-order absorption from the third delay
    #    into the central compartment). The control stream writes this as
    #    NT = 3; KTR = (NT+1)/MTT; K1T2 = K2T3 = K3T4 = KTR; K4T5 = KA.
    ktr <- 4 / mtt
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2 - ka  * transit3
    d/dt(central)     <-  ka  * transit3 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central  - k21 * peripheral1

    # 5. Bioavailability applies to the dosing (GUT) compartment: F1 = FREL.
    f(depot) <- frel

    # 6. Observation. The control stream scales the central compartment with
    #    S5 = VC/1000, converting an amount in mg and a volume in L to a
    #    concentration in ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
