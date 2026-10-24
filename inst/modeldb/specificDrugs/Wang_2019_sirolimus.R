Wang_2019_sirolimus <- function() {
  description <- "One-compartment population PK model with first-order absorption for oral sirolimus in Chinese children with kaposiform hemangioendothelioma, developed from routine therapeutic-drug-monitoring trough concentrations (Wang 2019). Apparent clearance CL/F rises exponentially with age and alanine aminotransferase and is 0.441-fold lower in girls than in boys; apparent volume V/F rises exponentially with the duration of sirolimus treatment. All covariates enter uncentered. The absorption rate constant is fixed at 0.485 per hour. Inter-individual variability is on CL/F and V/F and the residual error is exponential."
  reference <- paste(
    "Wang DD, Chen X, Li ZP. Population pharmacokinetics of sirolimus in",
    "pediatric patients with kaposiform hemangioendothelioma: A retrospective",
    "study. Oncol Lett. 2019;18(3):2412-2419. doi:10.3892/ol.2019.10562",
    sep = " "
  )
  vignette <- "Wang_2019_sirolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Wang 2019 Methods 'Analytical method': Emit 2000 Sirolimus Assay
  # (Siemens), a whole-blood immunoassay, on EDTA blood ('Drug
  # administration').
  compartmentData <- list(
    depot = list(analyte = "sirolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sirolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters CL/F uncentered as exp(0.215 * AGE) (Wang 2019 Results, final covariate model). Cohort 0.2-6 years, median 0.8, mean +/- SD 1.21 +/- 1.2 (Table I).",
      source_name = "age"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters CL/F uncentered as exp(0.0108 * ALT). Reported in IU/L, used interchangeably with U/L. Cohort 2-120, median 21, mean +/- SD 24.37 +/- 19.95 IU/L (Table I).",
      source_name = "ALT"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "The paper's covariate 'sex' enters CL/F as exp(-0.818 * sex) but its numeric coding is never printed. It is encoded here as sex = 1 (male) / 2 (female), i.e. sex = 1 + SEXF, the register's documented SEX alias (SEXF = as.integer(SEX == 2)). The direction is fixed by the Discussion ('males had a higher CL/F than females of the same age'); the 1/2 offset is chosen because only it reproduces the paper's own population-prediction level (Fig 2A, Fig 4) and agrees with the four other paediatric sirolimus popPK models' CL/F for a 7 kg child, including the same hospital's Chen 2020 KHE model. The literal 0/1 reading would put every male CL/F 2.27-fold higher. See the vignette's Assumptions and deviations section. Cohort 11 male / 6 female (Table I).",
      source_name = "sex"
    ),
    T_FIRSTDOSE = list(
      description = "Time elapsed since the first sirolimus dose (duration of treatment with sirolimus)",
      units = "h",
      type = "continuous",
      reference_category = NULL,
      notes = "The paper's DTT, in days, enters V/F uncentered as exp(0.0783 * DTT / 10). The column is carried in the canonical hours and divided by 24 inside model(). Time-varying within subject. Cohort DTT 7-749 days, median 117, mean +/- SD 154.2 +/- 149 (Table I).",
      source_name = "DTT"
    )
  )

  # Covariates screened in the stepwise search (Wang 2019 Methods
  # 'Covariate model') but not retained. Documentation only. Mean
  # corpuscular haemoglobin, mean corpuscular haemoglobin concentration,
  # and concomitant sulfamethoxazole, digoxin and ursodeoxycholic acid
  # were also screened; they have no registered canonical column and are
  # listed in population$notes.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened, not retained. The allometric model was assessed, but age plus sex gave the lower OFV (Discussion). Cohort 3.6-18 kg, median 7 (Table I)."
    ),
    DOSE_SIROLIMUS_MGD = list(
      description = "Daily sirolimus dose",
      units = "mg/day",
      type = "continuous",
      notes = "Screened ('DDS'), not retained. 0.16-1.5 mg/day, median 0.5 (Table I). Not a registered canonical; recorded only to document the screen."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened, not retained. 12-165 IU/L, median 35 (Table I)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Screened ('Cr'), not retained. 7-36 umol/L, median 18 (Table I)."
    ),
    BUN = list(
      description = "Serum urea",
      units = "mmol/L",
      type = "continuous",
      notes = "Screened ('UR'), not retained. 1-43 mmol/L, median 2.7 (Table I)."
    ),
    HCT = list(
      description = "Hematocrit",
      units = "%",
      type = "continuous",
      notes = "Screened, not retained. 22.7-43.9%, median 33.9 (Table I)."
    ),
    HGB = list(
      description = "Hemoglobin",
      units = "g/L",
      type = "continuous",
      notes = "Screened, not retained. 71.2-151 g/L, median 110 (Table I)."
    ),
    CONMED_METHYLPREDNISOLONE = list(
      description = "Concomitant methylprednisolone",
      units = "(binary)",
      type = "binary",
      notes = "Screened, not retained. 10 of 17 patients (Table II)."
    ),
    CONMED_CIMETIDINE = list(
      description = "Concomitant cimetidine",
      units = "(binary)",
      type = "binary",
      notes = "Screened, not retained. 1 of 17 patients (Table II)."
    ),
    CONMED_OMEPRAZOLE = list(
      description = "Concomitant omeprazole",
      units = "(binary)",
      type = "binary",
      notes = "Screened, not retained. 5 of 17 patients (Table II)."
    ),
    CONMED_FLUCONAZOLE = list(
      description = "Concomitant fluconazole",
      units = "(binary)",
      type = "binary",
      notes = "Screened, not retained. 1 of 17 patients (Table II)."
    ),
    CONMED_CEFOPERAZONE = list(
      description = "Concomitant cefoperazone/sulbactam",
      units = "(binary)",
      type = "binary",
      notes = "Screened ('cefoperazone/sulbactam'), not retained. 4 of 17 patients (Table II)."
    ),
    CONMED_VORICONAZOLE = list(
      description = "Concomitant voriconazole",
      units = "(binary)",
      type = "binary",
      notes = "Screened, not retained. 1 of 17 patients (Table II)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 17L,
    n_studies = 1L,
    n_observations = 133L,
    age_range = "0.2-6 years (median 0.8; mean +/- SD 1.21 +/- 1.2)",
    weight_range = "3.6-18 kg (median 7; mean +/- SD 7.99 +/- 3.04)",
    sex_female_pct = 35.3,
    race_ethnicity = "Chinese (single centre, Shanghai).",
    disease_state = "Paediatric kaposiform hemangioendothelioma treated with oral sirolimus.",
    dose_range = "Starting dose 0.8 mg/m^2 twice daily, then TDM-adjusted; observed 0.16-1.5 mg/day (median 0.5).",
    regions = "China (Children's Hospital of Fudan University, Shanghai), January 2016 - April 2018.",
    assay = "Emit 2000 Sirolimus Assay (Siemens), LLOQ 3.5 ng/mL, linear range 3.5-30 ng/mL.",
    sampling = "Retrospective real-world TDM; every concentration is a pre-dose trough, measured twice weekly or more often.",
    target_range = "Trough 10-15 ng/mL.",
    notes = "Covariates screened but not retained and without a registered canonical column: mean corpuscular haemoglobin (MCH), mean corpuscular haemoglobin concentration (MCHC), and concomitant sulfamethoxazole (16 of 17), digoxin (1 of 17) and ursodeoxycholic acid (1 of 17) (Methods, Table II)."
  )

  ini({
    # Final-model estimates, Wang 2019 Table IV. Every covariate enters
    # uncentered, so the CL/F intercept is the value at AGE = 0, ALT = 0 and
    # sex code 0, which no subject has; see the SEXF note.
    lcl <- log(3.19); label("CL/F equation intercept at AGE = 0, ALT = 0, sex code 0 (L/h)")  # Table IV: CL/F = 3.19 L/h (SE 23.7%)
    lvc <- log(165); label("V/F at zero duration of treatment (L)")                            # Table IV: V/F = 165 L (SE 26.9%)

    # Methods 'PPK modeling': Ka 'was fixed to 0.485 h-1 based on a previous
    # study (17)'; Table IV lists '0.485 (fixed)'.
    lka <- fixed(log(0.485)); label("Absorption rate constant Ka (1/h)")  # Table IV: Ka = 0.485 1/h, held constant

    # Results, final covariate model:
    #   CL/F = theta_CL/F * EXP(theta_AGE * age) * EXP(theta_ALT * ALT) * EXP(theta_SEX * sex)
    #   V/F  = theta_V/F * EXP[theta_DTT * (DTT/10)]
    e_age_cl <- 0.215; label("Exponential coefficient of AGE on CL/F (1/year)")                   # Table IV: theta AGE = 0.215 (SE 21.7%)
    e_alt_cl <- 0.0108; label("Exponential coefficient of ALT on CL/F (L/U)")                     # Table IV: theta ALT = 0.0108 (SE 18%)
    e_sexf_cl <- -0.818; label("Exponential coefficient of the sex code (1 + SEXF) on CL/F (unitless)")  # Table IV: theta SEX = -0.818 (SE 15.4%)
    e_t_firstdose_vc <- 0.0783; label("Exponential coefficient of treatment duration (days/10) on V/F (unitless)")  # Table IV: theta DTT = 0.0783 (SE 7.8%)

    # IIV, exponential (Methods equation 1). Table IV omega rows are read
    # as SDs of eta, the convention of this group's papers, and squared
    # here. The paper's own residual and prediction spreads (Fig 2A-C)
    # match the SD reading, not the variance reading; see the vignette.
    etalcl ~ 0.047089 # Table IV: omega CL/F = 0.217 (SE 20.7%), an SD; 0.217^2
    etalvc ~ 1.779556 # Table IV: omega V/F = 1.334 (SE 15.2%), an SD; 1.334^2

    # Residual error, Methods equation (3): Y = IPRED * exp(eps). Table IV
    # sigma 1 = 0.365 read as the SD of eps, by the same notation as the
    # omega rows.
    expSd <- 0.365; label("Exponential residual error SD (log scale)")  # Table IV: sigma 1 = 0.365 (SE 6.8%)
  })

  model({
    ka <- exp(lka)

    # The paper's sex code is 1 = male, 2 = female; see covariateData$SEXF.
    sex_code <- 1 + SEXF
    # DTT is in days in the paper; T_FIRSTDOSE is in hours.
    dtt_days <- T_FIRSTDOSE / 24

    cl <- exp(lcl + etalcl) *
      exp(e_age_cl * AGE) *
      exp(e_alt_cl * ALT) *
      exp(e_sexf_cl * sex_code)
    vc <- exp(lvc + etalvc) * exp(e_t_firstdose_vc * dtt_days / 10)

    kel <- cl / vc

    # One-compartment disposition with first-order absorption; F is not
    # identifiable from oral trough data and is absorbed into CL/F and V/F.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Doses in mg and volume in L give mg/L; x1000 for ng/mL.
    Cc <- central / vc * 1000
    Cc ~ lnorm(expSd)
  })
}
