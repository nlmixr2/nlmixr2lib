Thoueille_2023_tenofovir_full <- function() {
  description <- paste(
    "One-compartment population PK model for tenofovir after oral tenofovir",
    "alafenamide in people living with HIV, with creatinine clearance, age,",
    "ethnicity and potent P-glycoprotein-inhibitor comedication on apparent",
    "clearance (full final covariate model)."
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
  # Doses were converted to nmol of tenofovir alafenamide and concentrations to
  # nmol/mL of tenofovir for the population analysis (Methods, "PopPK
  # analysis"). Complete, irreversible 1:1 molar conversion of tenofovir
  # alafenamide into tenofovir is assumed (Methods, "Structural model"), so a
  # dose expressed in nmol of tenofovir alafenamide enters the model as the
  # same number of nmol of tenofovir.
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot   = list(analyte = "tenofovir full", units = "nmol", specimen = "administration site", verified = FALSE),
    central = list(analyte = "tenofovir full", units = "nmol", specimen = "plasma", verified = FALSE)
  )

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
        "Reference value CL_CR-ref = 100 mL/min (Table 2 footnote). Enters as a",
        "linear centered fractional effect on CL/F. Study range 33-203 mL/min,",
        "median 93 mL/min (Table 1). CLCR was the single most influential",
        "covariate, explaining 53 of the 59 percentage points of between-subject",
        "variability on F that the full covariate model accounts for; the",
        "authors' sibling reduced model retains CLCR alone. Precedent for the",
        "raw, non-BSA-normalized Cockcroft-Gault form under this canonical:",
        "Delattre_2010_amikacin.R."
      ),
      source_name = "CLCR"
    ),
    AGE = list(
      description = "Age at baseline.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Centered on the study median Age_median = 51 years (Table 2 footnote).",
        "Enters as a linear centered fractional effect on CL/F. Study range",
        "19-79 years (Table 1). The effect is statistically significant but",
        "small: the paper reports CL/F falling to 36 L/h at 80 years versus 42",
        "L/h at the median age, and the simulated effect of age on tenofovir",
        "Cmin within each CKD stage was judged minor (Figure S5)."
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
        "The authors first fitted a separate clearance per ethnic group (rich",
        "model), then pooled to this two-way split; the reduced grouping cost",
        "only dOFV = +2 (P > 0.05) relative to the full ethnicity model",
        "(Results). Model-building cohort composition was White 71%, Black 23%,",
        "Hispanic American 4%, Asian 2% with no 'Other' level (Table 1), so the",
        "indicator and its reference exactly partition the cohort. Table 2",
        "labels the coefficient row 'theta_Black or Hispano American'. No",
        "clinically relevant difference in median tenofovir Cmin between ethnic",
        "groups was predicted (Results, 'Simulations')."
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
        "risperidone, ritonavir, sertraline; the actual positives in the",
        "model-building cohort were carvedilol and sertraline, 45 of 486",
        "subjects (9%). Moderate P-gp inhibitors (atazanavir, diltiazem,",
        "duloxetine, efavirenz) and P-gp inducers (rifabutine, nevirapine) were",
        "screened separately and did NOT enter the final model. Cobicistat is",
        "itself a P-gp inhibitor but was deliberately excluded from this",
        "covariate and is carried separately as CONMED_COBICISTAT on relative",
        "bioavailability (Methods, 'Covariate model'); do not double-count it",
        "here. The paper reports the effect as independent of the presence of",
        "cobicistat."
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
        "perfectly confounded with dose level: every cobicistat-treated",
        "subject received 10 mg and every other subject 25 mg. It therefore",
        "encodes the combined dose-normalization and boosting effect, not a",
        "pure drug-drug-interaction effect. 277 of 486 subjects (43%) received",
        "cobicistat (Table 1)."
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
        "Screened and statistically significant on CL/F (dOFV = -242,",
        "P < 0.001) but fit less well than Cockcroft-Gault CLCR (dOFV = -285);",
        "only CLCR was retained. Median 87, range 23-153 (Table 1)."
      )
    ),
    WT = list(
      description = "Body weight at baseline.",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Screened and significant in univariate analysis on CL/F (dOFV = -15,",
        "P < 0.001) but eliminated during backward deletion and application of",
        "the principle of parsimony. Median 74, range 42-142 (Table 1)."
      )
    ),
    HT = list(
      description = "Body height at baseline.",
      units = "cm",
      type = "continuous",
      notes = "Screened, not retained. Median 173, range 145-195 (Table 1)."
    ),
    BMI = list(
      description = "Body mass index at baseline.",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened, not retained. Median 24, range 15-51 (Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine at baseline.",
      units = "umol/L",
      type = "continuous",
      notes = paste(
        "Screened, not retained as a standalone covariate; it enters the model",
        "indirectly through the Cockcroft-Gault CLCR. Median 85, range 42-256",
        "(Table 1)."
      )
    ),
    SEXF = list(
      description = "Sex (1 = female).",
      units = "(binary)",
      type = "binary",
      notes = "Screened, not retained. 30% female in the model-building cohort (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 486,
    n_observations = 793,
    n_studies = 1,
    age_range = "19-79 years",
    age_median = "51 years",
    weight_range = "42-142 kg",
    weight_median = "74 kg",
    height_median = "173 cm",
    bmi_median = "24 kg/m^2",
    sex_female_pct = 30,
    race_ethnicity = c(White = 71, Black = 23, `Hispanic American` = 4, Asian = 2),
    disease_state = "HIV-1 infection (people living with HIV) on tenofovir alafenamide-based antiretroviral therapy",
    renal_function = "serum creatinine median 85 (range 42-256) umol/L; Cockcroft-Gault CLCR median 93 (range 33-203) mL/min; CKD-EPI eGFR median 87 (range 23-153) mL/min/1.73 m^2",
    co_medication = "potent P-glycoprotein inhibitors 9%, moderate P-gp inhibitors 2%, P-gp inducers 1%, OATP1B1/1B3 inhibitors 3%, nephrotoxic drugs 10%; cobicistat 43%",
    dose_range = "tenofovir alafenamide 25 mg once daily, or 10 mg once daily with cobicistat",
    regions = "Switzerland",
    notes = paste(
      "Baseline demographics from Thoueille 2023 Table 1, model-building",
      "dataset column (n = 486 subjects, 793 tenofovir concentrations).",
      "Subjects were enrolled in Swiss HIV Cohort Study project #815 or",
      "followed in the Lausanne routine therapeutic-drug-monitoring programme",
      "between January 2017 and January 2021. Mostly sparse sampling (1-2",
      "samples per subject) plus four subjects with detailed 0-24 h profiles.",
      "Steady state was assumed for all individuals. An independent external",
      "validation dataset of 83 subjects / 84 observations is described in the",
      "same table but was not used for estimation; external validation gave a",
      "non-significant mean prediction error of 3.6% (95% CI 0.2-7.1%) with a",
      "precision of 17%."
    )
  )

  ini({
    # Structural parameters -- Thoueille 2023 Table 2 ("Final" column).
    # ka was fixed, not estimated: the small number of samples collected soon
    # after dose intake could not support it, and tenofovir shows flip-flop
    # kinetics after tenofovir alafenamide dosing (Methods, "Structural
    # model"; Table 2 reports "2 FIX").
    lka <- fixed(log(2)); label("First-order absorption rate constant ka (1/h)")                          # Table 2 (ka = 2 h^-1 FIX)
    lcl <- log(42.2); label("Apparent clearance CL/F of tenofovir at reference covariates (L/h)")           # Table 2 (CL_TFV = 42.2 L/h)
    lvc <- log(2390); label("Apparent central volume of distribution V/F of tenofovir (L)")                 # Table 2 (V_TFV = 2390 L)
    # F is a relative bioavailability anchored at 1 for the 25 mg
    # tenofovir alafenamide regimen without cobicistat; Table 2 reports
    # "1 FIX".
    lfdepot <- fixed(log(1)); label("Relative bioavailability F, tenofovir alafenamide 25 mg without cobicistat (unitless)")  # Table 2 (F = 1 FIX)

    # Covariate effects -- Table 2 "Final model:" equations (identical algebraic
    # form to Table S3):
    #   TVCL_TFV = CL_TFV * (1 + theta_CLCR * (CLCR - CLCR_ref) / CLCR_ref)
    #                     * (1 + theta_Age  * (Age - Age_median) / Age_median)
    #                     * (1 + theta_Black_or_Hispano_American)
    #                     * (1 + theta_Potent_P_gp_inhibitors)
    #     with CLCR_ref = 100 mL/min and Age_median = 51 years
    #   TVF      = 1 + theta_Cobicistat
    e_crcl_cl <- 0.707; label("Fractional CLCR effect on CL/F per unit relative deviation from 100 mL/min")   # Table 2 (theta_CLCR = 0.707)
    e_age_cl <- -0.244; label("Fractional age effect on CL/F per unit relative deviation from 51 years")      # Table 2 (theta_Age = -0.244)
    e_race_black_hispanic_cl <- 0.119; label("Fractional change in CL/F for Black or Hispanic American vs White or Asian")  # Table 2 (theta_Black or Hispano American = 0.119)
    e_conmed_pgp_inh_cl <- -0.121; label("Fractional change in CL/F with potent P-glycoprotein-inhibitor comedication")     # Table 2 (theta_Potent P-gp inhibitors = -0.121)
    e_cobi_fdepot <- 1.15; label("Fractional increase in F with cobicistat coadministration")                 # Table 2 (theta_Cobicistat = 1.15)

    # BSV was supported only on F; adding BSV on ka, CL or V gave dOFV > -1
    # (Results, "Structural, statistical and covariate models"). Table 2
    # reports BSV on F as 21.1% CV. Parameters were assumed log-normally
    # distributed (Methods, "Structural model"), so the variance is
    # omega^2 = log(1 + CV^2) = log(1 + 0.211^2) = 0.0435584.
    etalfdepot ~ 0.0435584                                                                                   # Table 2 (BSV on F = 21.1% CV)

    # An additive error model best described tenofovir residual variability
    # (Results, "Structural, statistical and covariate models").
    addSd <- 0.0099; label("Additive residual error (nmol/mL)")                                               # Table 2 (sigma_add = 0.0099 nmol/mL)
  })

  model({
    ka <- exp(lka)
    # Multiplicative product of linear centered / fractional covariate effects
    # on apparent clearance; references CLCR 100 mL/min and age 51 years.
    cl <- exp(lcl) *
      (1 + e_crcl_cl * (CRCL - 100) / 100) *
      (1 + e_age_cl * (AGE - 51) / 51) *
      (1 + e_race_black_hispanic_cl * RACE_BLACK_HISPANIC) *
      (1 + e_conmed_pgp_inh_cl * CONMED_PGP_INH)
    vc <- exp(lvc)
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Relative bioavailability carries the only between-subject random effect
    # and the cobicistat / dose-normalization effect.
    f(depot) <- exp(lfdepot + etalfdepot) * (1 + e_cobi_fdepot * CONMED_COBICISTAT)

    # `central` is in nmol and `vc` in L, so `central / vc` is nmol/L; the
    # extra factor of 1000 converts to the nmol/mL scale the authors used for
    # the concentrations and for sigma_add.
    Cc <- central / vc / 1000
    Cc ~ add(addSd)
  })
}
