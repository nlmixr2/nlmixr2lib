Benkali_2010_tacrolimus <- function() {
  description <- "Two-compartment population PK model with Erlang-distributed transit absorption (3 transit compartments) for once-daily extended-release oral tacrolimus (Advagraf) in stable adult renal transplant recipients more than 6 months post-transplant who were switched from twice-daily ciclosporin (Benkali 2010), with a multiplicative CYP3A5*1-carrier (expresser) effect on apparent clearance and combined additive + proportional residual error."
  reference <- paste(
    "Benkali K, Rostaing L, Premaud A, Woillard JB, Saint-Marcoux F,",
    "Urien S, Kamar N, Marquet P, Rousseau A. Population pharmacokinetics",
    "and Bayesian estimation of tacrolimus exposure in renal transplant",
    "recipients on a new once-daily formulation.",
    "Clin Pharmacokinet. 2010;49(10):683-692.",
    sep = " "
  )
  vignette <- "Benkali_2010_tacrolimus"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CYP3A5_EXPR = list(
      description        = "CYP3A5 expresser indicator: 1 if the patient carries at least one functional CYP3A5*1 allele (genotype *1/*1 or *1/*3), 0 if homozygous *3/*3.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5 *3/*3 nonexpresser)",
      notes              = "Time-fixed germline genotype determined from the rs776746 (CYP3A5 6986A>G) polymorphism by TaqMan allelic-discrimination assay. The *1 (A) allele encodes functional CYP3A5; the *3 (G) allele creates a cryptic splice site and yields nonfunctional protein. In the Benkali 2010 cohort the genotype distribution was *1/*1 = 1, *1/*3 = 4, *3/*3 = 36 (Table I), so CYP3A5_EXPR = 1 for the 5 *1 carriers and 0 for the 36 nonexpressers. The covariate enters CL/F as `(1 + theta2)^CYP3A5_EXPR` with theta2 = 1.15 (Eq. 4 and Table II), giving expressers a CL/F approximately 2.15-fold higher than nonexpressers (close to the 2-fold value the paper reports in the Abstract / Discussion).",
      source_name        = "cyp"
    )
  )

  population <- list(
    n_subjects                = 41L,
    n_studies                 = 1L,
    n_observations            = 492L,
    age_range                 = "28-77 years",
    age_median                = "52 years",
    weight_range              = "45-116 kg",
    weight_median             = "68 kg",
    sex_female_pct            = 53.7,
    race_ethnicity            = "Not reported in source paper (single-country French cohort).",
    disease_state             = "Stable adult renal transplant recipients more than 6 months post-transplant who were switched from twice-daily ciclosporin to once-daily extended-release tacrolimus (Advagraf). All patients received concomitant mycophenolate mofetil; concomitant prednisolone median dose 2.5 mg (range 0-10 mg).",
    dose_range                = "Once-daily oral tacrolimus titrated to a target trough concentration of 4-8 ng/mL; the mean daily dose used in the bootstrap-VPC dose-normalisation was 4.63 mg.",
    regions                   = "France (Toulouse / Limoges).",
    cyp3a5_distribution       = "*1/*1 n = 1, *1/*3 n = 4, *3/*3 n = 36 (Table I); 5 expressers (CYP3A5_EXPR = 1) and 36 nonexpressers (CYP3A5_EXPR = 0).",
    sampling_window           = "Twelve blood samples per subject within a single 24-hour dosing interval at predose, 0.33, 0.66, 1, 1.5, 2, 3, 4, 6, 9, 12, and 24 hours post-dose (Methods section, Blood Collection).",
    baseline_demographics_other = "Median (range): haematocrit 38.5% (26.5-45.1); haemoglobin 13 g/dL (10.5-15.1); serum creatinine 113 umol/L (82-907); prednisolone dose 2.5 mg (0-10) (Table I).",
    notes                     = "Patients were randomly partitioned into a model-building dataset (n = 29) and a validation dataset (n = 12); the final model in Table II second-major-column was re-fit to all 41 patients (492 concentration-time observations). Tacrolimus quantified in whole blood by a validated LC-MS/MS assay with a 1 ug/L lower limit of quantification."
  )

  ini({
    # Final-model estimates from Table II "Final model obtained with the whole
    # dataset" column (n = 41), Benkali 2010, p. 688. Table II reports two
    # rate constants for the central <-> peripheral exchange (k56 = central
    # to peripheral, k65 = peripheral to central) rather than Q/F and V2/F,
    # so this file mirrors the source's rate-constant parameterisation
    # (kcp / kpc) instead of the canonical lq / lvp form. The implied
    # secondary parameters are Q/F = kcp * V1/F = 0.09 * 486 = 43.74 L/h and
    # V2/F = Q/F / kpc = 43.74 / 0.13 = 336.5 L (reported in vignette).

    # Erlang transit absorption: a chain of n = 3 transit compartments
    # connected by a common rate constant ktr (Methods + Table II). The
    # number of transits was selected by Benkali by increasing n until the
    # objective-function value stopped improving (p. 685, "Population
    # Pharmacokinetic Analysis").
    lktr <- log(3.3);  label("Erlang transit absorption rate constant ktr (1/h)")                              # Table II whole-dataset ktr = 3.3 1/h

    # Central volume of distribution (apparent, after oral dosing).
    lvc  <- log(486);  label("Apparent central volume V1/F (L)")                                                # Table II whole-dataset V1/F = 486 L

    # Apparent oral clearance for the CYP3A5 nonexpresser reference subject
    # (CYP3A5_EXPR = 0). The paper's Eq. 4 expresses CL/F as
    #   CL/F = theta1 * (1 + theta2)^CYP3A5_EXPR
    # with theta1 = 19 L/h (nonexpresser typical) and theta2 = 1.15
    # (expresser fold-effect parameter); see Eq. 4 and Table II.
    lcl  <- log(19);   label("Apparent oral clearance CL/F for CYP3A5 nonexpressers (L/h)")                     # Table II whole-dataset theta1 = 19 L/h

    # Multiplicative CYP3A5 expresser fold-effect on CL/F. Stored as the
    # full multiplier (1 + theta2) so the model() body can use the
    # `multiplier^CYP3A5_EXPR` form (matching `Bergmann_2014_tacrolimus`).
    e_cyp3a5_expr_cl <- 2.15; label("CYP3A5 expresser multiplicative factor on CL/F (1 + theta2; expressers have ~2.15-fold higher CL/F)")  # Table II whole-dataset theta2 = 1.15 -> multiplier 2.15

    # Central <-> peripheral rate constants (paper's k56 / k65). Carried as
    # log-rate-constant parameters because the source estimated them
    # directly rather than as Q/F and V2/F.
    lkcp <- log(0.09); label("Central -> peripheral rate constant k56 (1/h)")                                   # Table II whole-dataset k56 = 0.09 1/h
    lkpc <- log(0.13); label("Peripheral -> central rate constant k65 (1/h)")                                   # Table II whole-dataset k65 = 0.13 1/h

    # Inter-individual variability. Benkali 2010 used exponential IIV models
    # (Methods, Population Pharmacokinetic Analysis: "Interindividual
    # variability was described using exponential models"), so the reported
    # %CVs translate to log-scale variance via omega^2 = log(1 + CV^2):
    #   ktr   CV 52% -> log(1 + 0.52^2) = 0.2393
    #   V1/F  CV 53% -> log(1 + 0.53^2) = 0.2475
    #   CL/F  CV 35% -> log(1 + 0.35^2) = 0.1156
    #   k56   CV 54% -> log(1 + 0.54^2) = 0.2560
    # Table II does not report an IIV estimate for k65 (peripheral ->
    # central), so no eta is carried on lkpc.
    etalktr ~ 0.2393                                                                                            # Table II whole-dataset IIV ktr = 52% CV
    etalvc  ~ 0.2475                                                                                            # Table II whole-dataset IIV V1/F = 53% CV
    etalcl  ~ 0.1156                                                                                            # Table II whole-dataset IIV CL/F = 35% CV
    etalkcp ~ 0.2560                                                                                            # Table II whole-dataset IIV k56 = 54% CV

    # Combined additive + proportional residual error (Table II foot of
    # column "Final model obtained with the whole dataset"; corroborated in
    # Results: "low residual error (0.7 ng/mL and 9% for the additive and
    # the proportional parts, respectively)").
    propSd <- 0.09; label("Proportional residual error (fraction)")                                             # Table II whole-dataset proportional error = 9%
    addSd  <- 0.7;  label("Additive residual error (ng/mL)")                                                    # Table II whole-dataset additive error = 0.7 ng/mL
  })

  model({
    # Individual PK parameters. CYP3A5 expresser status enters CL/F as the
    # multiplicative factor (1 + theta2)^CYP3A5_EXPR per Eq. 4 of the
    # source paper.
    ktr <- exp(lktr + etalktr)
    vc  <- exp(lvc  + etalvc)
    cl  <- exp(lcl  + etalcl) * e_cyp3a5_expr_cl^CYP3A5_EXPR
    kcp <- exp(lkcp + etalkcp)
    kpc <- exp(lkpc)

    kel <- cl / vc

    # Two-compartment disposition with Erlang transit absorption. The chain
    # has three transit compartments (depot + transit1 + transit2) all
    # connected by the same ktr, matching the n = 3 Erlang structure that
    # Benkali identified as best (p. 688, Table II). Dose enters `depot`;
    # the final transit (`transit2`) feeds into `central` at rate ktr.
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot     - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1  - ktr * transit2
    d/dt(central)     <-  ktr * transit2  - kel * central - kcp * central + kpc * peripheral1
    d/dt(peripheral1) <-  kcp * central   - kpc * peripheral1

    # Tacrolimus is reported in whole blood in ng/mL. Dose units mg,
    # volumes L -> central / vc returns mg/L; multiply by 1000 to express
    # the prediction in ng/mL.
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
