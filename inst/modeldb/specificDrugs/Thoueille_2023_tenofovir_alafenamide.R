Thoueille_2023_tenofovir_alafenamide <- function() {
  description <- paste(
    "Joint two-analyte population PK model fitting tenofovir alafenamide and",
    "its active moiety tenofovir simultaneously in people living with HIV,",
    "with complete irreversible prodrug conversion and creatinine clearance,",
    "age, ethnicity and potent P-glycoprotein-inhibitor comedication on",
    "apparent tenofovir clearance."
  )
  reference <- paste(
    "Thoueille P, Alves Saldanha S, Desfontaine V, Kusejko K, Courlet P,",
    "Andre P, Cavassini M, Decosterd LA, Buclin T, Guidi M, and the Swiss HIV",
    "Cohort Study. Population pharmacokinetic modelling to characterize the",
    "effect of chronic kidney disease on tenofovir exposure after tenofovir",
    "alafenamide administration. J Antimicrob Chemother. 2023;78:1433-1443.",
    "doi:10.1093/jac/dkad103"
  )
  vignette <- "Thoueille_2023_tenofovir_ckd"
  # Doses were converted to nmol of tenofovir alafenamide and concentrations of
  # both analytes to nmol/mL for the population analysis (Methods, "PopPK
  # analysis"). Tenofovir alafenamide is the dosed parent and therefore keeps
  # the canonical unsuffixed compartment / parameter / output names; tenofovir
  # carries the registered `_tfv` metabolite suffix.
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/mL")

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Creatinine clearance estimated with the Cockcroft-Gault equation.",
        "NOTE: raw Cockcroft-Gault in mL/min, NOT normalized to 1.73 m^2 body",
        "surface area, so the values are not interchangeable with the",
        "BSA-normalized eGFR variants that also map to this canonical."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Reference value CL_CR-ref = 100 mL/min (Table S3 footnote). Enters as",
        "a linear centered fractional effect on tenofovir apparent clearance",
        "only; tenofovir alafenamide clearance carries no covariates. Study",
        "range 33-203 mL/min, median 93 mL/min (Table 1). No difference in",
        "tenofovir alafenamide concentrations was predicted between patients",
        "with CKD and those with normal renal function (Results,",
        "'Simulations'). Precedent for the raw, non-BSA-normalized",
        "Cockcroft-Gault form under this canonical: Delattre_2010_amikacin.R."
      ),
      source_name = "CLCR"
    ),
    AGE = list(
      description = "Age at baseline.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Centered on the study median Age_median = 51 years (Table S3",
        "footnote). Enters as a linear centered fractional effect on tenofovir",
        "apparent clearance. Study range 19-79 years (Table 1)."
      ),
      source_name = "Age"
    ),
    RACE_BLACK_HISPANIC = list(
      description = paste(
        "Composite ethnicity indicator: 1 = Black or Hispanic American,",
        "0 = White or Asian (reference)."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White or Asian)",
      notes = paste(
        "Same two-way pooled grouping as the tenofovir-alone models; Table S3",
        "labels the coefficient row 'theta_Black or Hispano American'.",
        "Model-building cohort composition was White 71%, Black 23%, Hispanic",
        "American 4%, Asian 2% with no 'Other' level (Table 1), so the",
        "indicator and its reference exactly partition the cohort."
      ),
      source_name = "Black or Hispano American"
    ),
    CONMED_PGP_INH = list(
      description = paste(
        "Concomitant potent P-glycoprotein-inhibitor indicator: 1 = receiving a",
        "potent P-gp inhibitor, 0 = not receiving one."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no potent P-glycoprotein inhibitor)",
      notes = paste(
        "ONLY the potent tier is pooled into the 1 category. Screened potent",
        "set (Table 1 footnote d): amiodarone, carvedilol, clarithromycin,",
        "fluoxetine, itraconazole, methadone, quetiapine, rilpivirine,",
        "risperidone, ritonavir, sertraline. Moderate P-gp inhibitors and P-gp",
        "inducers were screened separately and did NOT enter the final model.",
        "Cobicistat is itself a P-gp inhibitor but was deliberately excluded",
        "from this covariate and is carried separately as CONMED_COBICISTAT on",
        "relative bioavailability (Methods, 'Covariate model'); do not",
        "double-count it here."
      ),
      source_name = "Potent P-gp inhibitors"
    ),
    CONMED_COBICISTAT = list(
      description = paste(
        "Concomitant cobicistat coadministration indicator: 1 = tenofovir",
        "alafenamide 10 mg once daily given with cobicistat, 0 = tenofovir",
        "alafenamide 25 mg once daily without cobicistat."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (tenofovir alafenamide 25 mg once daily, no cobicistat; F fixed to 1)",
      notes = paste(
        "Cobicistat was forced into the model as a covariate on relative",
        "bioavailability F to absorb the 25 mg vs 10 mg dose difference",
        "(Methods, 'Structural model'), so in this dataset the indicator is",
        "perfectly confounded with dose level. The Discussion notes cobicistat",
        "may also act on tenofovir alafenamide cellular uptake by inhibiting",
        "P-gp-mediated efflux from PBMCs, a process this model cannot resolve",
        "separately from bioavailability. 35 of the 54 subjects with measured",
        "tenofovir alafenamide concentrations (65%) received cobicistat",
        "(Table S1)."
      ),
      source_name = "Cobicistat"
    )
  )

  # Covariates the authors screened but did NOT retain in the final model.
  # Documented here for provenance; none is referenced in model().
  covariatesDataExcluded <- list(
    CRCL_EPI = list(
      description = "CKD-EPI estimated glomerular filtration rate (BSA-normalized).",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = paste(
        "Screened and significant on tenofovir CL/F but fit less well than",
        "Cockcroft-Gault CLCR; only CLCR was retained (Results)."
      )
    ),
    WT = list(
      description = "Body weight at baseline.",
      units = "kg",
      type = "continuous",
      notes = "Screened and significant univariately, eliminated during backward deletion."
    ),
    HT = list(
      description = "Body height at baseline.",
      units = "cm",
      type = "continuous",
      notes = "Screened, not retained."
    ),
    BMI = list(
      description = "Body mass index at baseline.",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened, not retained."
    ),
    CREAT = list(
      description = "Serum creatinine at baseline.",
      units = "umol/L",
      type = "continuous",
      notes = "Screened, not retained standalone; enters indirectly via Cockcroft-Gault CLCR."
    ),
    SEXF = list(
      description = "Sex (1 = female).",
      units = "(binary)",
      type = "binary",
      notes = "Screened, not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 486,
    n_observations = "793 tenofovir and 100 tenofovir alafenamide concentrations (47 of the latter below the limit of quantification)",
    n_studies = 1,
    age_range = "19-79 years",
    age_median = "51 years",
    weight_range = "42-142 kg",
    weight_median = "74 kg",
    sex_female_pct = 30,
    race_ethnicity = c(White = 71, Black = 23, `Hispanic American` = 4, Asian = 2),
    disease_state = "HIV-1 infection (people living with HIV) on tenofovir alafenamide-based antiretroviral therapy",
    renal_function = "serum creatinine median 85 (range 42-256) umol/L; Cockcroft-Gault CLCR median 93 (range 33-203) mL/min; CKD-EPI eGFR median 87 (range 23-153) mL/min/1.73 m^2",
    dose_range = "tenofovir alafenamide 25 mg once daily, or 10 mg once daily with cobicistat",
    regions = "Switzerland",
    notes = paste(
      "Baseline demographics from Thoueille 2023 Table 1, model-building",
      "dataset column (n = 486 subjects). Only 54 of those subjects",
      "contributed tenofovir alafenamide measurements; their separate",
      "characteristics are in Table S1 (median age 62 years, range 22-76;",
      "median weight 75 kg; 83% male; 89% White; median CLCR 72 mL/min).",
      "Because tenofovir alafenamide has a plasma half-life of about 0.5 h,",
      "only samples drawn within 6 h of dosing were assayed for it, giving 100",
      "samples of which 47 were below the limit of quantification; those were",
      "handled with the M6 method. Steady state was assumed for all",
      "individuals."
    )
  )

  ini({
    # Structural parameters -- Thoueille 2023 Table S3 ("Final model" column).
    # k12, the depot-to-tenofovir-alafenamide absorption rate constant, was
    # fixed because too few samples were collected right after drug intake
    # (Methods, "Structural model"; Table S3 reports k12 = 2 with no RSE).
    lka <- fixed(log(2)); label("First-order absorption rate constant k12, depot to tenofovir alafenamide (1/h)")  # Table S3 (k12 = 2 h^-1, fixed)

    # Tenofovir alafenamide (dosed parent, canonical unsuffixed names).
    lcl <- log(211); label("Apparent clearance CL/F of tenofovir alafenamide (L/h)")                          # Table S3 (CL_TAF = 211 L/h)
    lvc <- log(49.8); label("Apparent central volume of distribution V/F of tenofovir alafenamide (L)")        # Table S3 (V_TAF = 49.8 L)

    # Tenofovir (active moiety, `_tfv` metabolite suffix).
    lcl_tfv <- log(42.2); label("Apparent clearance CL/F of tenofovir at reference covariates (L/h)")          # Table S3 (CL_TFV = 42.2 L/h)
    lvc_tfv <- log(2380); label("Apparent central volume of distribution V/F of tenofovir (L)")                # Table S3 (V_TFV = 2380 L)

    # F is a relative bioavailability anchored at 1 for the 25 mg
    # tenofovir alafenamide regimen without cobicistat; Table S3 reports
    # "1 FIX".
    lfdepot <- fixed(log(1)); label("Relative bioavailability F, tenofovir alafenamide 25 mg without cobicistat (unitless)")  # Table S3 (F = 1 FIX)

    # Covariate effects -- Table S3 "Final model:" equations. All act on
    # tenofovir clearance; tenofovir alafenamide clearance carries none.
    #   TVCL_TFV = CL_TFV * (1 + theta_CLCR * (CLCR - CLCR_ref) / CLCR_ref)
    #                     * (1 + theta_Age  * (Age - Age_median) / Age_median)
    #                     * (1 + theta_Black_or_Hispano_American)
    #                     * (1 + theta_Potent_P_gp_inhibitors)
    #     with CLCR_ref = 100 mL/min and Age_median = 51 years
    #   TVF      = 1 + theta_Cobicistat
    e_crcl_cl_tfv <- 0.707; label("Fractional CLCR effect on tenofovir CL/F per unit relative deviation from 100 mL/min")  # Table S3 (theta_CLCR = 0.707)
    e_age_cl_tfv <- -0.250; label("Fractional age effect on tenofovir CL/F per unit relative deviation from 51 years")     # Table S3 (theta_Age = -0.250)
    e_race_black_hispanic_cl_tfv <- 0.122; label("Fractional change in tenofovir CL/F for Black or Hispanic American vs White or Asian")  # Table S3 (theta_Black or Hispano American = 0.122)
    e_conmed_pgp_inh_cl_tfv <- -0.116; label("Fractional change in tenofovir CL/F with potent P-glycoprotein-inhibitor comedication")     # Table S3 (theta_Potent P-gp inhibitors = -0.116)
    e_cobi_fdepot <- 1.15; label("Fractional increase in F with cobicistat coadministration")                              # Table S3 (theta_Cobicistat = 1.15)

    # BSV was significant on F (dOFV = -420) and on tenofovir alafenamide
    # clearance (dOFV = -130); adding BSV on k12, CL_TFV, V_TFV or V_TAF did
    # not improve the fit (dOFV > -1) (Supplementary Methods, "Simultaneous
    # model development"). Parameters were assumed log-normally distributed,
    # so omega^2 = log(1 + CV^2):
    #   F       21.2% CV -> log(1 + 0.212^2) = 0.0439633
    #   CL_TAF  93.2% CV -> log(1 + 0.932^2) = 0.6252023
    etalfdepot ~ 0.0439633                                                                                    # Table S3 (BSV on F = 21.2% CV)
    etalcl ~ 0.6252023                                                                                        # Table S3 (BSV on CL_TAF = 93.2% CV)

    # Tenofovir alafenamide residual variability was best described by a mixed
    # (proportional plus additive) error model and tenofovir by an additive
    # model, with no correlation term between them (Supplementary Methods,
    # "Simultaneous model development").
    propSd <- 0.607; label("Proportional residual error for tenofovir alafenamide (fraction)")                 # Table S3 (sigma_propTAF = 60.7 CV%)
    addSd <- 2.35e-3; label("Additive residual error for tenofovir alafenamide (nmol/mL)")                     # Table S3 (sigma_addTAF = 2.35x10^-3 nmol/mL)
    addSd_tfv <- 9.92e-3; label("Additive residual error for tenofovir (nmol/mL)")                             # Table S3 (sigma_addTFV = 9.92x10^-3 nmol/mL)
  })

  model({
    ka <- exp(lka)

    # Tenofovir alafenamide: no covariates, BSV on clearance.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)

    # Tenofovir: multiplicative product of linear centered / fractional
    # covariate effects; references CLCR 100 mL/min and age 51 years.
    cl_tfv <- exp(lcl_tfv) *
      (1 + e_crcl_cl_tfv * (CRCL - 100) / 100) *
      (1 + e_age_cl_tfv * (AGE - 51) / 51) *
      (1 + e_race_black_hispanic_cl_tfv * RACE_BLACK_HISPANIC) *
      (1 + e_conmed_pgp_inh_cl_tfv * CONMED_PGP_INH)
    vc_tfv <- exp(lvc_tfv)

    # Tenofovir alafenamide is cleared ONLY by conversion into tenofovir
    # ("a complete and irreversible conversion of tenofovir alafenamide into
    # tenofovir was assumed"), so the whole k23 = CL_TAF / V_TAF flux is the
    # input to the tenofovir compartment. Both compartments hold nmol and the
    # conversion is 1:1 molar, so no stoichiometric factor is needed.
    k23 <- cl / vc
    kel_tfv <- cl_tfv / vc_tfv

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - k23 * central
    d/dt(central_tfv) <- k23 * central - kel_tfv * central_tfv

    # Relative bioavailability carries a between-subject random effect and the
    # cobicistat / dose-normalization effect.
    f(depot) <- exp(lfdepot + etalfdepot) * (1 + e_cobi_fdepot * CONMED_COBICISTAT)

    # Amounts are in nmol and volumes in L, so amount / volume is nmol/L; the
    # extra factor of 1000 converts to the nmol/mL scale the authors used for
    # the concentrations and for the residual-error estimates.
    Cc <- central / vc / 1000
    Cc_tfv <- central_tfv / vc_tfv / 1000

    Cc ~ add(addSd) + prop(propSd)
    Cc_tfv ~ add(addSd_tfv)
  })
}
