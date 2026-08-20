Cojutti_2023_ceftobiprole <- function() {
  description <- "Three-compartment IV population PK model for ceftobiprole in adults with severe Gram-positive infections (real-life multicentre therapeutic drug monitoring cohort, Italy). Clearance rises exponentially with CKD-EPI estimated glomerular filtration rate; central volume V1 is larger in males. Supports probability-of-target-attainment analysis against free-trough or free-steady-state fCtrough/MIC and fCss/MIC targets."
  reference   <- "Cojutti PG, Giuliano S, Pascale R, Angelini J, Tascini C, Viale P, Pea F. Population Pharmacokinetic and Pharmacodynamic Analysis for Maximizing the Effectiveness of Ceftobiprole in the Treatment of Severe Methicillin-Resistant Staphylococcal Infections. Microorganisms. 2023;11(12):2964. doi:10.3390/microorganisms11122964"
  vignette    <- "Cojutti_2023_ceftobiprole"
  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central     = list(analyte = "ceftobiprole", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "ceftobiprole", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral2 = list(analyte = "ceftobiprole", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "CKD-EPI estimated glomerular filtration rate, BSA-normalized",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column eGFR. Cojutti 2023 Methods section 2.1: 'The CKD-EPI formula [18] was used to",
        "estimate patient eGFR.' Table 1 reports a cohort median (IQR) of 83.7 (50.5-101.7)",
        "mL/min/1.73 m^2; the Monte Carlo simulations in section 2.4 span five renal-function classes",
        "(<30, 30-50, 51-80, 81-130, >130 mL/min/1.73 m^2). Stored under the canonical CRCL column per",
        "inst/references/covariate-columns.md, which accepts a BSA-normalized CKD-EPI eGFR.",
        "IMPORTANT - the effect is applied UNCENTERED and EXPONENTIALLY,",
        "cl = exp(lcl) * exp(e_crcl_cl * CRCL), so exp(lcl) = 1.7 L/h is the CL intercept extrapolated",
        "to CRCL = 0 and not a CL at any physiological reference value. See the e_crcl_cl comment in",
        "ini() for the arithmetic that establishes this functional form against the paper's own",
        "reported median individual CL. eGFR is a repeated clinical-chemistry measurement over the TDM",
        "course, so it is time-varying within subject."
      ),
      source_name        = "eGFR"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female) is the source paper's reference category for V1",
      notes              = paste(
        "Source column Gender. Cojutti 2023 Table 1 reports male/female = 86/46 (34.8% female).",
        "Table 2 reports a single 'Gender on V1' coefficient of 0.39; the paper does not state which",
        "sex is the reference. Back-solving against the paper's own reported median individual V1 of",
        "19.7 L (Results section 3.2) identifies FEMALE as the reference: with V1_female = 14.77 L the",
        "male value is 14.77 * exp(0.39) = 21.8 L, and because the median subject in an 86-male /",
        "46-female cohort is male the predicted median is 21.8 L (+10.8% vs the reported 19.7 L). The",
        "opposite assignment (male as reference) predicts 14.77 L, which is 25% BELOW the reported",
        "value. The larger male V1 is also the physiologically expected direction. To preserve the",
        "paper's published V1 = 14.77 L verbatim while storing sex under the canonical SEXF (1 =",
        "female) convention, the effect is applied in model() as exp(e_sex_vc * (1 - SEXF)), following",
        "the Bajaj_2017_nivolumab.R precedent in this registry."
      ),
      source_name        = "Gender"
    )
  )

  # Covariates screened by Cojutti 2023 (Methods section 2.3: "the following
  # covariates were tested on the pharmacokinetic parameters: age, gender,
  # weight, height, and eGFR") but NOT retained by the forward/backward
  # selection. Only gender and eGFR survived into the final model, so age,
  # weight and height carry no coefficient and are documentation only.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in the forward/backward covariate search (Methods section 2.3) but not retained in the final model. Cohort median (IQR) 71.0 (61.8-79.0) years per Table 1."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened in the forward/backward covariate search (Methods section 2.3) but not retained in the final model. Cohort median (IQR) 73.5 (65.0-89.0) kg per Table 1. No allometric scaling is applied in the published model."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened in the forward/backward covariate search (Methods section 2.3) but not retained in the final model. Cojutti 2023 does not report a height summary in Table 1 (only the derived BMI, median 25.7 kg/m^2)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 132L,
    n_studies      = 1L,
    n_observations = 503L,
    age_range      = "Adults; median (IQR) 71.0 (61.8-79.0) years (Table 1)",
    weight_range   = "Median (IQR) 73.5 (65.0-89.0) kg (Table 1)",
    bmi_range      = "Median (IQR) 25.7 (22.5-30.1) kg/m^2 (Table 1)",
    sex_female_pct = 34.8,
    race_ethnicity = "Not reported (two Italian tertiary university hospital cohorts)",
    disease_state  = paste(
      "Adults with suspected or documented severe Gram-positive infection receiving ceftobiprole with",
      "therapeutic drug monitoring. Infection sites (Table 1): hospital-acquired pneumonia 38 (28.8%),",
      "endocarditis 27 (20.5%), bloodstream infection 22 (16.6%), community-acquired pneumonia 20",
      "(15.2%), bone and joint 9 (6.8%), skin and soft tissue 9 (6.8%), device-related 4 (3.0%), CNS 3",
      "(2.3%). A microbiological isolate was identified in 80/132 (60.6%)."
    ),
    renal_function = "CKD-EPI eGFR median (IQR) 83.7 (50.5-101.7) mL/min/1.73 m^2; serum creatinine median (IQR) 0.90 (0.68-1.36) mg/dL (Table 1)",
    dose_range     = paste(
      "eGFR-adjusted starting regimens administered as 3 h extended infusions (Methods section 2.1):",
      "500 mg q8h if eGFR >= 50, 500 mg q12h if eGFR 30-50, and 250 mg q12h if eGFR < 30",
      "mL/min/1.73 m^2, with subsequent TDM-guided adjustment. Median (IQR) daily dose 1500",
      "(1000-1500) mg; median (IQR) treatment duration 10.0 (2.0-81.0) days (Table 1)."
    ),
    co_medication  = "Monotherapy in 44/132 (33.3%); combination therapy in 88/132 (66.7%), most often with daptomycin (27/88), ampicillin (24/88), or fosfomycin (8/88)",
    regions        = "Italy (IRCCS Azienda Ospedaliero-Universitaria di Bologna and Azienda Sanitaria Universitaria Friuli Centrale, Udine)",
    notes          = paste(
      "Retrospective TDM cohort, January 2018 to December 2022 (Methods section 2.1). Baseline",
      "demographics per Table 1. Sampling was predominantly trough (Ctrough measured 5 days weekly)",
      "with additional samples at the end of and 1 h after the 3 h infusion where feasible; observed",
      "median (min-max) Ctrough was 7.4 (0.7-37.2) mg/L, end-of-infusion 18.4 (4.6-55.3) mg/L, and",
      "1 h post-infusion 14.8 (5.0-55.5) mg/L (Results section 3.2). Assay was LC-MS/MS with an LLOQ",
      "of 0.5 mg/L (Methods section 2.2). Table 1 labels serum albumin as 'g/L' with a median of 3.1,",
      "which is a unit-label error in the source (3.1 g/dL = 31 g/L); albumin was not among the",
      "screened covariates and is not used by this model."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural parameters: Cojutti 2023 Table 2, "Covariate Model Estimate
    # (% RSE)" column (the FINAL model). The paper's compartment numbering
    # maps onto the nlmixr2lib canonical names as
    #   V1 -> vc (central), Q2/V2 -> q/vp (peripheral1),
    #   Q3/V3 -> q2/vp2 (peripheral2).
    # Note that peripheral2 is the SMALL, FAST compartment in this model
    # (Q3 = 42.41 L/h into V3 = 4.76 L); the paper's numbering is preserved
    # rather than reordered by equilibration speed.
    # ---------------------------------------------------------------------

    # CL intercept. Because the eGFR effect is uncentered and exponential
    # (see e_crcl_cl below), exp(lcl) is the clearance extrapolated to
    # CRCL = 0 and is NOT a typical clearance at any physiological eGFR.
    # This is exactly why the reported CL falls from 3.43 L/h in the base
    # model to 1.7 L/h once the uncentered covariate enters.
    lcl  <- log(1.7);   label("Clearance intercept at CRCL = 0 (L/h)")                # Cojutti 2023 Table 2 covariate model: CL = 1.7 L/h (11.9% RSE)
    lvc  <- log(14.77); label("Central volume V1 in females (L)")                     # Cojutti 2023 Table 2 covariate model: V1 = 14.77 L (19.7% RSE); female is the reference category (see covariateData$SEXF)
    lq   <- log(6.11);  label("Intercompartmental clearance Q2 to peripheral1 (L/h)") # Cojutti 2023 Table 2 covariate model: Q2 = 6.11 L/h (13.6% RSE)
    lvp  <- log(46.16); label("Peripheral volume V2 (L)")                             # Cojutti 2023 Table 2 covariate model: V2 = 46.16 L (34.8% RSE)
    lq2  <- log(42.41); label("Intercompartmental clearance Q3 to peripheral2 (L/h)") # Cojutti 2023 Table 2 covariate model: Q3 = 42.41 L/h (61.6% RSE; the paper flags Q3 as the least precisely estimated structural parameter)
    lvp2 <- log(4.76);  label("Peripheral volume V3 (L)")                             # Cojutti 2023 Table 2 covariate model: V3 = 4.76 L (13.9% RSE)

    # ---------------------------------------------------------------------
    # Covariate effects: Cojutti 2023 Table 2, "Covariate effect" block.
    #
    # FUNCTIONAL FORM. Results section 3.2 says the covariates were added "as
    # a power function of ceftobiprole CL and V1, respectively", but the
    # reported coefficient of 0.011 CANNOT be a power exponent: any power
    # form gives CL = 1.7 * (83.7/ref)^0.011 ~= 1.70 L/h at the cohort median
    # eGFR, i.e. essentially no effect, which contradicts both the abstract
    # ("Estimated glomerular filtration rate significantly affected drug
    # clearance") and the paper's own reported median individual CL of
    # 4.04 L/h (Results section 3.2). The exponential (linear-on-log-scale)
    # form, which is Monolix's default for an untransformed continuous
    # covariate, reproduces it: 1.7 * exp(0.011 * 83.7) = 4.27 L/h, within
    # 5.7% of the reported 4.04 L/h. Candidate forms scored at the median
    # eGFR of 83.7 mL/min/1.73 m^2 against the reported 4.04 L/h were
    #   exponential 1.7*exp(0.011*eGFR)  = 4.27  (+5.7%)   <- adopted
    #   power       1.7*(eGFR/70)^0.011  = 1.70  (-57.9%)
    #   power       1.7*(eGFR/90)^0.011  = 1.70  (-57.9%)
    #   power raw   1.7*eGFR^0.011       = 1.79  (-55.8%)
    #   linear add  1.7 + 0.011*eGFR     = 2.62  (-35.1%)
    #   linear mult 1.7*(1+0.011*eGFR)   = 3.27  (-19.2%)
    # The control for this arithmetic is V2, which carries no covariate: its
    # typical value 46.16 L and its reported median individual value 46.6 L
    # agree to 1.0%, confirming that "median individual ~= typical value" is
    # a reliable discriminator in this paper. The residual 5.7% gap for CL is
    # consistent with eGFR being time-varying across the TDM course (the
    # exponential form back-solves to a median-over-records eGFR of 78.7 vs
    # the Table 1 baseline median of 83.7). Trusting the numbers over the
    # prose follows the standing policy for text-vs-value conflicts.
    e_crcl_cl <- 0.011; label("Exponential coefficient of CRCL on CL (per mL/min/1.73 m^2; uncentered)") # Cojutti 2023 Table 2 covariate model: eGFR on CL = 0.011 (12.6% RSE); bootstrap median 0.01 (95% CI 0.010-0.012)

    # Applied as exp(e_sex_vc * (1 - SEXF)) so that SEXF = 1 (female, the
    # source reference category) recovers V1 = 14.77 L verbatim and SEXF = 0
    # (male) gives 14.77 * exp(0.39) = 21.8 L. See covariateData$SEXF for the
    # back-solve that identifies female as the reference.
    e_sex_vc  <- 0.39;  label("Exponential coefficient of male sex on V1 (unitless; applied as (1 - SEXF))") # Cojutti 2023 Table 2 covariate model: Gender on V1 = 0.39 (54.8% RSE); bootstrap median 0.38 (95% CI 0.29-0.49)

    # ---------------------------------------------------------------------
    # Inter-individual variability. Cojutti 2023 Table 2 reports IIV as %CV
    # and the CL-V1 random-effect correlation as 0.6. Monolix parameters are
    # log-normal, so the log-scale variance is recovered with the standard
    # log-normal identity omega^2 = log(1 + CV^2):
    #   CL 47.4% -> 0.202676     V1 75.3% -> 0.449169
    #   Q2 35.0% -> 0.115558     V2 268.1% -> 2.102640
    #   Q3 235.2% -> 1.876698    V3 30.7% -> 0.090068
    # The alternative reading (that the paper printed 100 * omega rather than
    # a true CV) was tested and rejected: it is indistinguishable for CL, V1,
    # Q2 and V3, but for V2 it implies omega = 2.681, which over-predicts the
    # paper's own reported individual V2 range of 7.6-526.0 L by roughly 85x
    # (predicted 0.05-44520 L for n = 132), whereas omega = 1.450 predicts
    # 1.1-1894 L -- over-wide by only ~3.6x, which is what empirical-Bayes
    # shrinkage on a poorly identified peripheral volume produces.
    # Correlation: cov(etalcl, etalvc) = 0.6 * sqrt(0.202676 * 0.449169)
    #            = 0.6 * 0.450196 * 0.670200 = 0.181033.
    etalcl + etalvc ~ c(0.202676,
                        0.181033, 0.449169)  # Cojutti 2023 Table 2 covariate model: IIV on CL = 47.4 %CV, IIV on V1 = 75.3 %CV, correlation CL-V1 = 0.6 (24.8% RSE)
    etalq   ~ 0.115558                       # Cojutti 2023 Table 2 covariate model: IIV on Q2 = 35.0 %CV (34.5% RSE)
    etalvp  ~ 2.102640                       # Cojutti 2023 Table 2 covariate model: IIV on V2 = 268.1 %CV (19.4% RSE)
    etalq2  ~ 1.876698                       # Cojutti 2023 Table 2 covariate model: IIV on Q3 = 235.2 %CV (45.1% RSE); the paper flags the random effect on Q3 as imprecisely estimated
    etalvp2 ~ 0.090068                       # Cojutti 2023 Table 2 covariate model: IIV on V3 = 30.7 %CV (30.0% RSE)

    # ---------------------------------------------------------------------
    # Residual error. Cojutti 2023 Table 2 reports a combined error model
    # with a = 1.44 and b = 0.22; the table footnote states "a and b
    # represent additive and proportional residual error model,
    # respectively." The paper does not say whether Monolix's combined1
    # (sd = a + b*f) or combined2 (sd = sqrt(a^2 + (b*f)^2)) parameterization
    # was used; nlmixr2's add() + prop() is combined2. See the vignette
    # Assumptions and deviations section.
    addSd  <- 1.44; label("Additive residual error (mg/L)")            # Cojutti 2023 Table 2 covariate model: a = 1.44 (13.4% RSE)
    propSd <- 0.22; label("Proportional residual error (fraction)")    # Cojutti 2023 Table 2 covariate model: b = 0.22 (5.8% RSE)
  })

  model({
    # Derived sex term: Cojutti 2023 uses female as the V1 reference
    # category, so (1 - SEXF) reproduces the paper's male indicator while
    # keeping SEXF (1 = female) as the canonical storage convention.
    sex_male <- 1 - SEXF

    # Individual PK parameters. The eGFR effect on CL is exponential and
    # UNCENTERED (see the e_crcl_cl comment in ini()); the sex effect on V1
    # is exponential with female as the reference.
    cl  <- exp(lcl  + etalcl)  * exp(e_crcl_cl * CRCL)
    vc  <- exp(lvc  + etalvc)  * exp(e_sex_vc  * sex_male)
    q   <- exp(lq   + etalq)
    vp  <- exp(lvp  + etalvp)
    q2  <- exp(lq2  + etalq2)
    vp2 <- exp(lvp2 + etalvp2)

    # Three-compartment IV disposition micro-constants.
    kel <- cl  / vc
    k12 <- q   / vc
    k21 <- q   / vp
    k13 <- q2  / vc
    k31 <- q2  / vp2

    # ODE system. Ceftobiprole is given by intermittent extended infusion
    # (2-3 h) or by 24 h continuous infusion directly into the central
    # compartment; zero-order input is supplied through the event table
    # (dur() or rate on the dosing records), matching the paper's
    # "zero-order administration and first-order elimination from the
    # central compartment" (Methods section 2.3).
    d/dt(central)     <- -(kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <-    k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-    k13 * central - k31 * peripheral2

    # Total ceftobiprole plasma concentration in mg/L (dose in mg, volumes
    # in L). The PK/PD targets in the paper are expressed on the FREE
    # concentration, obtained as 0.84 * Cc (16% protein binding, Methods
    # section 2.4); that conversion is applied in the vignette rather than
    # in the model because the observed TDM data are total concentrations.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
