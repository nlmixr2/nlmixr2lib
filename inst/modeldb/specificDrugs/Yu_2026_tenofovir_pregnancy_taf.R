Yu_2026_tenofovir_pregnancy_taf <- function() {
  description <- paste(
    "Plasma tenofovir alafenamide (TAF) one-compartment first-order-absorption",
    "popPK model in non-pregnant, pregnant and postpartum women receiving oral TAF",
    "25 mg once daily. Non-pregnant (CONRAD 137 TAF arms) and pregnant / postpartum",
    "(IMPAACT P1026s and IMPAACT 2026 TAF arms, both without a pharmacokinetic",
    "booster) plasma TAF data were fitted jointly, with relative bioavailability F",
    "estimated separately for each of four pregnancy states (non-pregnant, fixed to",
    "1; second trimester; third trimester; 6-12 weeks postpartum) and clearance,",
    "volume and absorption shared across states. Doses and states are in umol and",
    "concentrations in umol/L (TAF 25 mg = 52.5 umol). This is the",
    "covariate-identification sub-model that supplied the pregnancy effect on TAF",
    "bioavailability used by the semi-mechanistic model Yu_2026_tenofovir.",
    sep = " "
  )
  reference <- paste(
    "Yu Y, Brooks KM, Doncel GF, Best BM, Marzinke MA, Mirochnick M, Anderson P,",
    "Myer L, Celum C, Heffron R, Coleman J, Joseph Davey D, Hendrix CW, Momper JD,",
    "Bies R, Scott RK. Development of a Semi-Mechanistic Population Pharmacokinetic",
    "Model for Predicting Tenofovir Disoproxil Fumarate and Tenofovir Alafenamide",
    "Exposure in Plasma and Cellular Matrices During Pregnancy and Postpartum.",
    "Clin Pharmacokinet. 2026;65(1):133-148. doi:10.1007/s40262-025-01589-y.",
    "Parameters from Electronic Supplementary Material Table S2.",
    sep = " "
  )
  vignette <- "Yu_2026_tenofovir"
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  covariateData <- list(
    EGA = list(
      description        = "Estimated gestational age of the mother",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used only to derive the paper's second- and third-trimester indicator",
        "variables TRI2 and TRI3; the model has no continuous EGA effect. Per the",
        "EGA register entry, a source paper that reports trimesters rather than",
        "weeks records the week values used per trimester here rather than",
        "introducing a trimester-indicator canonical. This model uses the standard",
        "obstetric boundaries TRI2 = 14 <= EGA < 28 weeks and TRI3 = EGA >= 28",
        "weeks. Yu 2026 does not print its own trimester cut-offs; the contributing",
        "IMPAACT P1026s protocol sampled the second trimester at 20-26 weeks and",
        "the third trimester at 30-38 weeks, both of which fall unambiguously",
        "inside these boundaries. Set EGA = 0 for non-pregnant and for postpartum",
        "subjects (the postpartum state is carried by TPP, not by EGA).",
        "No first-trimester data were available, so the model makes no prediction",
        "for 0 < EGA < 14 weeks; such a record is scored as non-pregnant.",
        sep = " "
      ),
      source_name        = "TRI2 / TRI3 (trimester indicator variables)"
    ),
    TPP = list(
      description        = "Time postpartum (time elapsed since delivery)",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used only to derive the paper's postpartum indicator PP; the model has no",
        "continuous TPP effect and applies a single step change for any TPP > 0.",
        "Yu 2026 sampled the postpartum state at 6-12 weeks after delivery, so the",
        "estimated postpartum bioavailability describes that window and should not",
        "be read as an immediately-post-delivery value. Set TPP = 0 during",
        "pregnancy and for non-pregnant subjects.",
        sep = " "
      ),
      source_name        = "PP (postpartum indicator variable)"
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "tenofovir alafenamide", units = "umol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tenofovir alafenamide", units = "umol", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species       = "human",
    n_subjects    = 53L,
    n_studies     = 3L,
    disease_state = paste(
      "Pregnant and postpartum women with HIV taking TAF 25 mg once daily without a",
      "pharmacokinetic booster (25 participants in the IMPAACT P1026s TAF arm and",
      "28 in the IMPAACT 2026 TAF arm, sampled in the second trimester, the third",
      "trimester and 6-12 weeks postpartum), jointly fitted with plasma TAF data",
      "from healthy HIV-negative non-pregnant women in CONRAD 137 (TAF 10 mg and",
      "25 mg arms).",
      sep = " "
    ),
    dose_range   = "TAF 25 mg orally once daily (52.5 umol); non-pregnant CONRAD 137 data also include a TAF 10 mg arm",
    regions      = "USA and international IMPAACT network sites; CONRAD 137 conducted in the USA",
    co_medication = paste(
      "Participants co-administered TAF with a pharmacokinetic booster (cobicistat",
      "or ritonavir) were excluded from the analysis, because P-glycoprotein",
      "inhibition raises TAF and TFV concentrations (Methods 2.1).",
      sep = " "
    ),
    notes = paste(
      "n_subjects counts the 25 + 28 = 53 pregnant IMPAACT participants (Table 1).",
      "The number of non-pregnant participants contributing plasma TAF to this",
      "particular fit is not stated separately; Table 1 reports 24 and 73",
      "non-pregnant women in the CONRAD 137 single- and multiple-dose phases (which",
      "may overlap), across TDF and TAF arms combined. Plasma TAF was heavily",
      "censored: 49.9-69.7% of samples were below the limit of quantification",
      "across the contributing arms, handled by the Beal M3 method. Age, weight,",
      "sex and race distributions are not reported in Yu 2026; the cohort is",
      "entirely female by design.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Relative bioavailability, estimated separately in each pregnancy
    # state with the non-pregnant state fixed to 1 as the structural
    # anchor. ESM Table S2 prints one absolute F per state, so each
    # state carries an explicit stratum suffix (parameter-names.md,
    # "Stratum-suffixed parameters") and none keeps the bare `lfdepot`.
    # Table S2 labels the F rows "(L/h)"; that unit is a copy-paste
    # error in the source table -- a bioavailability is unitless, and
    # the printed values (1, 0.827, 0.949, 1.18) are fractions that
    # reproduce the Results percentages exactly.
    # ------------------------------------------------------------------
    lfdepot_nonpreg <- fixed(log(1)); label("Relative bioavailability of oral TAF, non-pregnant (unitless)")          # ESM Table S2 (F Nonpregnant = 1, FIXED)
    lfdepot_tri2    <- log(0.827);    label("Relative bioavailability of oral TAF, second trimester (unitless)")      # ESM Table S2 (F 2nd Trimester = 0.827, RSE 16%); Results 3.2.2 "decrease of 17.3%"
    lfdepot_tri3    <- log(0.949);    label("Relative bioavailability of oral TAF, third trimester (unitless)")       # ESM Table S2 (F 3rd Trimester = 0.949, RSE 15%); Results 3.2.2 "decrease of ... 5.1%"
    lfdepot_pp      <- log(1.18);     label("Relative bioavailability of oral TAF, 6-12 weeks postpartum (unitless)") # ESM Table S2 (F Postpartum = 1.18, RSE 16%); Results 3.2.2 "increased by 18%"

    # Disposition parameters shared across all four pregnancy states.
    lka <- log(1.99); label("First-order absorption rate constant of oral TAF (1/h)")   # ESM Table S2 (Ka = 1.99 /h, RSE 4%)
    lcl <- log(145);  label("Apparent clearance CL/F of plasma TAF (L/h)")              # ESM Table S2 (CL = 145 L/h, RSE 13%)
    lvc <- log(46.2); label("Apparent volume of distribution V2/F of plasma TAF (L)")   # ESM Table S2 (V2 = 46.2 L, RSE 14%)

    # ------------------------------------------------------------------
    # Between-subject variability. ESM Table S2 reports BSV as %CV for
    # the exponential model P = TVP * exp(eta); omega^2 = log(CV^2 + 1).
    # ------------------------------------------------------------------
    etalka ~ 0.0220   # CV 14.9%: log(0.149^2 + 1)   (ESM Table S2, BSV on Ka, RSE 26%)
    etalcl ~ 0.3023   # CV 59.4%: log(0.594^2 + 1)   (ESM Table S2, BSV on CL, RSE 10%)
    etalvc ~ 0.6842   # CV 99.1%: log(0.991^2 + 1)   (ESM Table S2, BSV on V2, RSE 13%)

    # Residual error, plasma TAF. Table S2 reports a proportional term
    # only; unlike the semi-mechanistic model (Table 2) it carries no
    # additive TAF term.
    propSd <- 0.799; label("Proportional residual error, plasma TAF (fraction)")        # ESM Table S2 (PROP (TAF) = 79.9%, RSE 6%)
  })

  model({
    # 1. Pregnancy state. The paper's TRI2 / TRI3 / PP indicators are
    #    derived here from the canonical EGA and TPP covariate columns.
    #    Order matters: postpartum is tested first so that a record with
    #    TPP > 0 is scored postpartum regardless of any residual EGA.
    lfdepotState <-
      ifelse(TPP > 0, lfdepot_pp,
        ifelse(EGA >= 28, lfdepot_tri3,
          ifelse(EGA >= 14, lfdepot_tri2, lfdepot_nonpreg)))

    # 2. Individual parameters
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    # 3. Micro-constants
    kel <- cl / vc

    # 4. ODE system -- one-compartment with first-order oral absorption
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 5. Bioavailability -- the pregnancy effect enters here
    f(depot) <- exp(lfdepotState)

    # 6. Observation and error
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
