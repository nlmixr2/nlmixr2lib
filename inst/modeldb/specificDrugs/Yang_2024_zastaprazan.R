Yang_2024_zastaprazan <- function() {
  description <- "Two-compartment population PK model with Erlang-type absorption through six sequential first-order transit steps and first-order elimination for zastaprazan (JP-1366), a potassium-competitive acid blocker, after oral dosing in patients with erosive gastroesophageal reflux disease (GERD) and in healthy volunteers (Yang 2024). Pooled analysis of 1590 plasma concentrations from 160 subjects across one phase 1 study in healthy volunteers (intensive sampling) and one phase 2 study in patients with erosive GERD (trough sampling), fit in NONMEM 7.4.4 by FOCEI. All six absorption steps share a single transit rate constant Ktr, so the absorption-time distribution is Erlang with shape 6 and rate Ktr. Disease status is the only retained covariate: apparent clearance is 41.4% lower in patients with GERD than in healthy volunteers, via the linear form CL/F = 29.4 * (1 - 0.414 * DIS_GERD). Typical values are for a healthy volunteer. Inter-individual variability is log-normal on Ktr, CL/F, Vc/F and Vp/F (none on Q/F); residual variability is proportional. CYP2C19 phenotype had no meaningful effect on zastaprazan PK, consistent with CYP3A-mediated rather than CYP2C19-mediated metabolism."
  reference   <- "Yang E, Hwang I, Ji SC, Kim J, Lee S. Population pharmacokinetic analysis of zastaprazan (JP-1366), a novel potassium-competitive acid blocker, in patients and healthy volunteers. CPT Pharmacometrics Syst Pharmacol. 2024;13(12):2150-2158. doi:10.1002/psp4.13228."
  vignette    <- "Yang_2024_zastaprazan"
  units       <- list(time = "h", dosing = "mg", concentration = "ug/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. The absorption chain is the six sequential
  # compartments drawn in Yang 2024 Figure 1 (boxes 1-6): `depot` is box 1
  # (the dosing compartment) and `transit1`-`transit5` are boxes 2-6.
  compartmentData <- list(
    depot       = list(analyte = "zastaprazan", units = "mg", specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "zastaprazan", units = "mg", specimen = "administration site", verified = TRUE),
    transit2    = list(analyte = "zastaprazan", units = "mg", specimen = "administration site", verified = TRUE),
    transit3    = list(analyte = "zastaprazan", units = "mg", specimen = "administration site", verified = TRUE),
    transit4    = list(analyte = "zastaprazan", units = "mg", specimen = "administration site", verified = TRUE),
    transit5    = list(analyte = "zastaprazan", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "zastaprazan", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "zastaprazan", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    DIS_GERD = list(
      description        = "Erosive gastroesophageal reflux disease indicator: 1 = patient with endoscopically confirmed erosive GERD enrolled in the phase 2 study (Study 2), 0 = healthy volunteer without any gastrointestinal disorder enrolled in the phase 1 study (Study 1).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteer; the 68-subject phase 1 cohort of Study 1)",
      notes              = "Time-fixed per subject. Maps 1:1 onto the paper's `disease status` covariate, which Yang 2024 Table 2 footnote a codes as PT = 1 for patient and PT = 0 for healthy volunteer, so the reported coefficient and the reported typical CL/F are carried through unchanged with no sign flip or re-baselining. The effect is LINEAR fractional rather than exponential: Table 2 footnote a prints CL/F (L/h) = 29.4 * (1 - 0.414 * PT), giving 17.2 L/h in patients versus 29.4 L/h in healthy volunteers (a 41.4% reduction). Yang 2024 attributes the lower apparent clearance to delayed gastric emptying in GERD raising bioavailability F rather than to a change in intrinsic clearance, and notes the same direction of effect for vonoprazan (another P-CAB) and lansoprazole (a PPI). Cohort composition 92 patients with erosive GERD / 68 healthy volunteers (Yang 2024 Table 1).",
      source_name        = "disease status (PT)"
    )
  )

  # Covariates Yang 2024 formally screened in its stepwise covariate search
  # (Methods, 'Population pharmacokinetic analysis': forward selection at
  # p < 0.05, backward elimination at p < 0.01) but did NOT retain in the
  # final model -- "Other tested candidates were not identified as
  # significant covariates in the model" (Results). No point estimate is
  # published for any of them, so they are documentation only and are not
  # referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a continuous candidate covariate and not retained. Total-cohort median 39 years (range 19-74); the patient and healthy cohorts differ markedly (patients median 58, range 24-74; healthy volunteers median 28, range 19-45; Yang 2024 Table 1), so age is strongly confounded with DIS_GERD in this dataset."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened as a continuous candidate covariate and not retained; no allometric scaling appears in the final model. Total-cohort median 71.0 kg (range 43.4-106.2; Yang 2024 Table 1)."
    ),
    HT = list(
      description = "Body height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened as a continuous candidate covariate and not retained. Yang 2024 Table 1 does not tabulate height."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a continuous candidate covariate and not retained. Yang 2024 Table 1 does not tabulate AST."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a continuous candidate covariate and not retained. Yang 2024 Table 1 does not tabulate ALT."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened as a continuous candidate covariate and not retained. Yang 2024 Table 1 does not tabulate albumin, and the paper does not state the units in which it was screened."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened as a continuous candidate covariate and not retained. Yang 2024 Table 1 does not tabulate total bilirubin, and the paper does not state the units in which it was screened."
    ),
    SEXF = list(
      description = "Female sex indicator: 1 = female, 0 = male",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a categorical candidate covariate and not retained. 31 of 160 subjects (19.4%) were female, all of them in the patient cohort (the phase 1 healthy cohort was 100% male; Yang 2024 Table 1), so sex is partly confounded with DIS_GERD."
    ),
    CYP2C19_IM = list(
      description = "CYP2C19 intermediate-metabolizer phenotype indicator: 1 = intermediate metabolizer (*1/*2, *1/*3), 0 = otherwise",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as part of the categorical CYP2C19-phenotype covariate and not retained: 'CYP2C19 phenotypes had no meaningful effect on zastaprazan PK' and 'No visual correlations between CYP2C19 phenotypes and PK parameters including CL/F were observed' (Yang 2024 Results and Figure S1). 69 of 160 subjects (43.1%) were intermediate metabolizers; the reference phenotype is the normal metabolizer (*1/*1, 70 subjects, 43.8%) and phenotype was unidentified for 4 subjects (Yang 2024 Table 1). Consistent with zastaprazan being metabolized mainly by CYP3A4/CYP3A5 rather than CYP2C19."
    ),
    CYP2C19_PM = list(
      description = "CYP2C19 poor-metabolizer phenotype indicator: 1 = poor metabolizer (*2/*2, *2/*3, *3/*3), 0 = otherwise",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as part of the categorical CYP2C19-phenotype covariate and not retained (see CYP2C19_IM). 17 of 160 subjects (10.6%) were poor metabolizers (Yang 2024 Table 1)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 160L,
    n_studies      = 2L,
    n_observations = 1590L,
    age_range      = "Total 19-74 years (median 39); patients with erosive GERD 24-74 years (median 58); healthy volunteers 19-45 years (median 28)",
    age_median     = "39 years",
    weight_range   = "Total 43.4-106.2 kg (median 71.0); patients with erosive GERD 43.4-106.2 kg (median 69.7); healthy volunteers 56.4-88.4 kg (median 73.6)",
    weight_median  = "71.0 kg",
    sex_female_pct = 19.4,
    race_ethnicity = "Not reported; Yang 2024 Table 1 tabulates no race or ethnicity",
    disease_state  = "92 patients with erosive gastroesophageal reflux disease (phase 2 Study 2) and 68 healthy volunteers without any gastrointestinal disorder (phase 1 Study 1)",
    dose_range     = "Study 1 (healthy volunteers): oral single doses of 5-60 mg or once-daily multiple doses of 5-40 mg for 7 consecutive days. Study 2 (patients with erosive GERD): oral once-daily multiple doses of 10-40 mg for 4 or 8 weeks. The approved regimen for erosive GERD is 20 mg once daily.",
    regions        = "Not explicitly stated; both studies were run by investigators at Seoul National University College of Medicine and Hospital and sponsored by Onconic Therapeutics Inc. (Seoul, Republic of Korea)",
    genotype       = "CYP2C19 phenotype identified for 156 of 160 subjects: 70 normal metabolizers (*1/*1, 43.8%), 69 intermediate metabolizers (*1/*2 or *1/*3, 43.1%), 17 poor metabolizers (*2/*2, *2/*3 or *3/*3, 10.6%); genotyped by Affymetrix DMET Plus microarray in Study 1 and by PCR in Study 2",
    notes          = "Demographics from Yang 2024 Table 1 (median and range for age and weight, n (%) for categorical variables); study designs and sampling schemes from the Methods 'Data' section and Table S1. Sampling was intensive in the healthy volunteers (1495 of the 1590 observations came from 68 subjects) and sparse trough-only in the patients (95 observations from 92 subjects). Eleven of the 92 patients (12.0%) were Helicobacter pylori positive and every healthy volunteer was negative; H. pylori status was examined graphically (Figure S2) rather than as a formal covariate candidate and did not influence zastaprazan exposure."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters - Yang 2024 Table 2, 'Final model
    # Estimates (RSE%)' column. Typical values are at the paper's own
    # reference state, a healthy volunteer (DIS_GERD = 0), so every value
    # below is the verbatim Table 2 estimate with no re-baselining.
    # ------------------------------------------------------------------
    lktr <- log(13.6)   ; label("Erlang transit rate constant Ktr, shared by all six absorption steps (1/h)")  # Yang 2024 Table 2: Ktr = 13.6 1/h (RSE 3.9%; bootstrap median 13.6, 95% CI 12.5-14.6)
    lcl  <- log(29.4)   ; label("Apparent clearance CL/F in a healthy volunteer (L/h)")                        # Yang 2024 Table 2: CL/F = 29.4 L/h (RSE 5.1%; bootstrap median 29.4, 95% CI 26.5-32.6)
    lvc  <- log(97.3)   ; label("Apparent central volume of distribution Vc/F (L)")                            # Yang 2024 Table 2: Vc/F = 97.3 L (RSE 5.0%; bootstrap median 97.1, 95% CI 87.9-107.3)
    lq   <- log(29.1)   ; label("Apparent inter-compartmental clearance Q/F (L/h)")                            # Yang 2024 Table 2: Q/F = 29.1 L/h (RSE 6.7%; bootstrap median 29.0, 95% CI 25.3-33.1)
    lvp  <- log(148.0)  ; label("Apparent peripheral volume of distribution Vp/F (L)")                         # Yang 2024 Table 2: Vp/F = 148.0 L (RSE 4.6%; bootstrap median 148.8, 95% CI 135.0-162.8)

    # ------------------------------------------------------------------
    # Covariate effect - the single retained covariate. Yang 2024 Table 2
    # footnote a prints the LINEAR fractional form
    #   CL/F (L/h) = 29.4 * (1 - 0.414 * PT);  PT = 1 (patient) or 0 (healthy)
    # so the coefficient is carried as -0.414 and applied as
    # (1 + e_dis_gerd_cl * DIS_GERD) in model().
    # ------------------------------------------------------------------
    e_dis_gerd_cl <- -0.414 ; label("Fractional effect of DIS_GERD (erosive GERD) on CL/F (unitless)")         # Yang 2024 Table 2 'Effect of disease status' = -0.414 (RSE 8.6%; bootstrap median -0.414, 95% CI -0.480 to -0.341); footnote a gives the linear form and the Results state CL/F is 41.4% lower in patients (17.2 L/h) than in healthy volunteers (29.4 L/h)

    # ------------------------------------------------------------------
    # Inter-individual variability - Yang 2024 Table 2 reports IIV as CV%
    # for an exponential (log-normal) random-effect model (Methods:
    # "IIV for each PK parameter was evaluated using the exponential error
    # model on the assumption of a log-normal distribution"), so each CV%
    # is back-transformed to a log-scale variance with the exact
    # log-normal relation omega^2 = log(1 + CV^2). No IIV was retained on
    # Q/F ("The final model included the IIV for most PK parameters except
    # for Q/F", Results), so there is no etalq term.
    # ------------------------------------------------------------------
    etalktr ~ log(1 + 0.261^2)   # Yang 2024 Table 2: IIV Ktr = 26.1 %CV (RSE 17.8%; bootstrap median 25.1, 95% CI 14.0-34.6) -> omega^2 = 0.06590
    etalcl  ~ log(1 + 0.382^2)   # Yang 2024 Table 2: IIV CL/F = 38.2 %CV (RSE 6.7%; bootstrap median 37.8, 95% CI 32.2-43.4) -> omega^2 = 0.13620
    etalvc  ~ log(1 + 0.274^2)   # Yang 2024 Table 2: IIV Vc/F = 27.4 %CV (RSE 13.8%; bootstrap median 27.0, 95% CI 18.3-35.0) -> omega^2 = 0.07239
    etalvp  ~ log(1 + 0.280^2)   # Yang 2024 Table 2: IIV Vp/F = 28.0 %CV (RSE 10.3%; bootstrap median 27.5, 95% CI 21.0-33.6) -> omega^2 = 0.07548

    # ------------------------------------------------------------------
    # Residual error - proportional only ("the residual variability was
    # described by proportional error model", Results).
    # ------------------------------------------------------------------
    propSd <- 0.349 ; label("Proportional residual error (fraction)")   # Yang 2024 Table 2: proportional residual error = 34.9% (RSE 3.6%; bootstrap median 34.8, 95% CI 32.3-37.3)
  })

  model({
    # 1. Individual PK parameters. Only CL/F carries a covariate; the
    #    disease effect is the linear fractional form of Yang 2024 Table 2
    #    footnote a, CL/F = 29.4 * (1 - 0.414 * PT), so a typical patient
    #    with erosive GERD has CL/F = 29.4 * 0.586 = 17.2 L/h.
    ktr <- exp(lktr + etalktr)
    cl  <- exp(lcl + etalcl) * (1 + e_dis_gerd_cl * DIS_GERD)
    vc  <- exp(lvc + etalvc)
    q   <- exp(lq)
    vp  <- exp(lvp + etalvp)

    # 2. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. ODE system. Yang 2024 Figure 1 draws the absorption as six
    #    sequential compartments (boxes 1-6) connected by five Ktr arrows,
    #    with a sixth Ktr arrow from box 6 into the central compartment.
    #    The oral dose enters box 1. Encoded canonically as
    #    depot (box 1) -> transit1..transit5 (boxes 2-6) -> central, which
    #    is six first-order transfers all governed by the single estimated
    #    Ktr, i.e. an Erlang(shape = 6, rate = Ktr) absorption-time
    #    distribution with mean absorption time 6 / Ktr = 0.441 h.
    #    Dose records target cmt = "depot".
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2 - ktr * transit3
    d/dt(transit4)    <-  ktr * transit3 - ktr * transit4
    d/dt(transit5)    <-  ktr * transit4 - ktr * transit5
    d/dt(central)     <-  ktr * transit5 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                   k12 * central - k21 * peripheral1

    # 4. Observation. Every volume and clearance is apparent (divided by
    #    the unknown F), so no separate bioavailability term is estimated
    #    and none is applied to the depot. Amounts are in mg and volumes
    #    in L, so central / vc is mg/L; multiply by 1000 to report ug/L as
    #    in the source paper (simulated Cmax,ss 110-445 ug/L and AUCtau
    #    560-2283 h*ug/L, Yang 2024 Results).
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
