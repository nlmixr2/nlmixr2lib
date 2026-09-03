Salem_2013_cyp_ontogeny <- function() {
  description <- paste(
    "Drug-elimination-pathway ontogeny (system-parameter) model for",
    "birth to 20 years (Salem 2013 J Clin Pharmacol). Ten algebraic",
    "functions give the activity or protein expression of an",
    "elimination pathway relative to its own adult value: hepatic",
    "CYP1A2, CYP2B6, CYP2C8, CYP2C9, CYP2C18/19, CYP2D6, CYP2E1 and",
    "CYP3A4/5, enterocytic (gut) CYP3A4/5, and renal function. The",
    "nine CYP pathways follow the sigmoidal (hyperbolic) maturation",
    "form of Equation 1, fraction of adult = Fbirth + (AdultMax -",
    "Fbirth) * Age^n / (Agemid^n + Age^n), with Age in years; renal",
    "function follows a quadratic in body surface area normalised to",
    "an adult glomerular filtration rate of 120 mL/min. This is a",
    "SYSTEM layer, not a drug model: it has no drug, no dosing, no",
    "compartments and no ODEs, and every output is an algebraic",
    "function of the rxode2 time variable, which the model interprets",
    "as postnatal age in YEARS. Ratios of two pathway functions",
    "(Equation 2) quantify differential ontogeny and hence the",
    "age-dependence of the fraction metabolised by each route, which",
    "is what drives the age-dependence of metabolic drug-drug",
    "interactions. Valid from birth to 20 years. No between-subject",
    "variability is encoded: the paper propagates uncertainty by",
    "bootstrapping the residual variance about each regression line,",
    "and the per-pathway residual variances are not reported.",
    sep = " "
  )
  reference <- paste(
    "Salem F, Johnson TN, Barter ZE, Leeder JS, Rostami-Hodjegan A.",
    "Age related changes in fractional elimination pathways for drugs:",
    "assessing the impact of variable ontogeny on metabolic",
    "drug-drug interactions.",
    "J Clin Pharmacol. 2013;53(8):857-865.",
    "doi:10.1002/jcph.100.",
    sep = " "
  )
  vignette <- "Salem_2013_cyp_ontogeny"
  units <- list(
    time          = "year (postnatal age; valid 0 to 20)",
    dosing        = "n/a (no exogenous dosing; enzyme / renal ontogeny system model)",
    concentration = paste(
      "n/a (no drug concentration). Every output is a unitless",
      "fraction of the corresponding adult value."
    )
  )

  covariateData <- list(
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Drives the renal-function output only. Salem 2013",
        "Supplementary Table S1 gives renal function as a quadratic in",
        "BSA divided by 120 mL/min, the adult glomerular filtration",
        "rate: renal = ((-0.61604 * BSA^2) + (99.054 * BSA) - 17.74) /",
        "120. The paper states (Methods, Comparison of the Ontogeny of",
        "Drug Elimination Pathways) that 'for the renal model we used a",
        "linear model based on the relationship between glomerular",
        "filtration rate (GFR) and body surface area'. No age-to-BSA",
        "growth function is supplied by the paper, so BSA must be",
        "provided by the user; the nine CYP outputs do not use it.",
        "The quadratic is negative below BSA = 0.18 m^2 and reaches 1",
        "(adult GFR) near BSA = 1.42 m^2; keep BSA inside the paediatric",
        "-to-adult range 0.2 to 2.0 m^2."
      ),
      source_name        = "BSA"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Salem 2013 Supplementary Table S1 reports separate renal",
        "'fractional expression at birth relative to adult' summary",
        "values for males (0.15) and females (0.14), but the renal",
        "equation printed in the same row is a function of body surface",
        "area only and carries no sex term, and no sex term appears in",
        "any CYP row. Sex is therefore documented here rather than",
        "implemented. The Simcyp DDI simulations set the fraction of",
        "females to 0.5 in every age band (Methods, Hypothetical Age",
        "Related Changes in DDI)."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 10000L,
    n_studies      = 16L,
    age_range      = "Birth to 20 years (ontogeny functions evaluated over 0 to 20 years at 0.1-year increments; DDI simulations at age bands 1 day, 7 days, 1 month, 1 year, 2 years and 20 years)",
    weight_range   = "Not reported (Simcyp Paediatric v12 virtual population defaults)",
    sex_female_pct = 50,
    race_ethnicity = "Not reported",
    disease_state  = paste(
      "Healthy paediatric and adult virtual subjects. The ontogeny",
      "functions themselves are regressions on pooled in vitro human",
      "liver-bank expression and catalytic-activity data (Salem 2013",
      "Supplementary Table S1 references 1 to 16, e.g. Sonnier 1998",
      "CYP1A2, Croom 2009 CYP2B6, Treluyer 1997 and Koukouritaki 2004",
      "CYP2C, Stevens 2008 CYP2D6, Johnsrud 2003 CYP2E1, Lacroix 1997",
      "and Stevens 2003 CYP3A, Johnson 2001 enterocytic CYP3A4,",
      "Johnson 2006 renal), updated from the Johnson 2006 models. The",
      "per-study subject counts are not restated in this paper."
    ),
    dose_range     = "n/a (no exogenous drug PK is modelled by this system layer)",
    regions        = "Not reported",
    n_bootstrap    = 10000L,
    notes          = paste(
      "n_subjects is the size of the virtual cohort used for the",
      "hypothetical DDI simulations, not of the ontogeny regressions:",
      "'One hundred simulations with 10 trials and 10 subjects were",
      "performed with fraction of females set as 0.5 for each of three",
      "age bands, 1 day (0.00274-0.00276 years), 1 year (0.99-1.01",
      "years), and 20 years (19.99-20.01 years)' (Methods), i.e. 100 x",
      "10 x 10 = 10,000 virtual subjects per age band. Those",
      "simulations were run in Simcyp Paediatric v12 with proprietary",
      "system parameters (liver size, liver blood flow, enzyme",
      "abundance, fraction unbound) and are NOT reproducible outside",
      "that platform; only the ontogeny functions of Supplementary",
      "Table S1 are extracted here. Confidence intervals around each",
      "ontogeny line were obtained by 10,000 bootstrap iterations per",
      "pathway per age point (Methods, Equations 3 and 4) using the",
      "residual variance S0^2 about the regression; the per-pathway",
      "S0 values are not reported, so no variability is encoded."
    )
  )

  ini({
    # -----------------------------------------------------------------
    # Salem 2013 Equation 1 (page 858), per pathway:
    #
    #   fraction of adult = (AdultMax - Fbirth) * Age^n / (Agemid^n + Age^n)
    #                       + Fbirth
    #
    # where AdultMax is the maximum adult relative expression (1 at
    # younger adults, "but with relevant corrections if the reference
    # group were older"), Fbirth is the relative expression at birth,
    # Agemid is the age at which expression is midway between the birth
    # and adult values, and n is the sigmoidicity exponent. Age is in
    # years. Naming follows the registered system-model family used by
    # Deferm_2025_lactation_physiology.R: <pathway>_max / _birth / _t50
    # / _hill.
    #
    # The per-pathway constants come from Supplementary Table S1, whose
    # "Hyperbolic function" column is a set of ten embedded Microsoft
    # Equation 3.0 objects inside the supplementary .doc (OLE streams
    # ObjectPool/_1425384659 .. _1425384668 / "Equation Native"). Each
    # object was decoded and is quoted VERBATIM in the comment on its own
    # parameter block; the model reproduces the printed equation exactly.
    #
    # POLICY NOTE. Supplementary Table S1 also carries two descriptive
    # summary columns ("Fractional expression at birth relative to adult"
    # and "Time to half adult expression (Y)"). Those columns are NOT
    # reproducible from the printed equations for most rows and are
    # treated as approximate annotation, not as parameter sources. Taking
    # "time to half adult expression" as the age at which the function
    # reaches 0.5, the printed equations give 0.008 (column 0.008) for
    # CYP2C8, 0.19 (0.2) for CYP2C18/19, 0.57 (0.6) for hepatic CYP3A4/5,
    # 0.09 (0.08) for CYP2D6, 0.12 (0.1) for CYP2E1, 0.38 (0.3) for gut
    # CYP3A4/5, 0.007 (0.005) for CYP2C9, 1.32 (1.8) for CYP1A2 and 0.80
    # (1.2) for CYP2B6 - i.e. the column disagrees with the equation on
    # every row, by up to 40 %. Where a summary column and the printed
    # equation conflict, the EQUATION is used. See the vignette "Errata"
    # section, which tabulates every disagreement.
    #
    # Every value is a published point estimate of a fitted regression
    # with no reported standard error, so all are fixed().
    # -----------------------------------------------------------------

    # ---- Hepatic CYP1A2 ----
    # Supplementary Table S1 prints, verbatim:
    #   CYP1A2 = (1.05 - 0.08) * Age^1.1 / (1.69^1.1 + Age^1.1) + 0.08
    # Used as printed. Three other places in the paper are in tension
    # with the Fbirth = 0.08 term and are recorded in the vignette
    # "Errata" section rather than being used to overwrite it:
    #   (a) the Table S1 birth column for this row reads "Negligible";
    #   (b) Figure 1 annotates "CYP1A2 Birth level= 0";
    #   (c) the Discussion (and Supplementary Table S4) give a 20.9-fold
    #       CYP2C9-to-CYP1A2 disparity at day 1; this equation gives 5.0.
    # No alternative CYP1A2 parameter set is printed anywhere in the
    # paper, so the printed equation is the only sourced model. Fitting
    # a replacement curve to the Figure 1 trace would substitute three
    # unpublished values for three published ones and would require
    # AdultMax near 1.57, which contradicts the paper's own definition
    # of AdultMax ("1 at younger adults but with relevant corrections if
    # the reference group were older") and is far outside the 1.05-1.074
    # range of every other AdultMax in Table S1.
    cyp1a2_max <- fixed(1.05)
    label("CYP1A2 maximum adult relative expression AdultMax (fraction of adult)")  # Supplementary Table S1 CYP1A2 equation, numerator leading constant
    cyp1a2_birth <- fixed(0.08)
    label("CYP1A2 relative expression at birth Fbirth (fraction of adult)")  # Supplementary Table S1 CYP1A2 equation, subtracted and added constant. The row's summary column reads "Negligible" and Figure 1 annotates "Birth level= 0"; see vignette Errata
    cyp1a2_t50 <- fixed(1.69)
    label("CYP1A2 age at half-maximal maturation Agemid (year)")  # Supplementary Table S1 CYP1A2 equation denominator base
    cyp1a2_hill <- fixed(1.1)
    label("CYP1A2 sigmoidicity exponent n (unitless)")  # Supplementary Table S1 CYP1A2 equation exponent

    # ---- Hepatic CYP2B6 ----
    # Supplementary Table S1 prints, verbatim:
    #   CYP2B6 = (1 - 0.1) * Age / (1 + Age) + 0.1
    # i.e. AdultMax = 1, Fbirth = 0.1, Agemid = 1 y, n = 1, used as
    # printed. The decoded equation object for this row does contain
    # three superscript template boxes (one on the numerator "Age", one
    # on the denominator "1" and one on the denominator "Age") but ALL
    # THREE ARE EMPTY, so the equation typesets - and a reader of the
    # published supplement sees - a bare "Age" with no exponent, which
    # is n = 1. The distinction is visible in the OLE structure: this
    # row carries 4 templates (1 fraction + 3 empty superscripts), the
    # same count as CYP1A2 / CYP2C18-19 / CYP2E1 / CYP3A4-5 which have
    # filled exponents, whereas CYP2C8, CYP2D6 and gut CYP3A4/5 - the
    # rows that genuinely have no exponent - carry only 1 template.
    # The row's summary columns (Fbirth 0.15, time to half adult
    # expression 1.2 y) are not reachable from the printed constants
    # under ANY exponent: with Fbirth = 0.1 and Agemid = 1 the function
    # reaches 0.5 at Age = 0.8^(1/n), which is < 1 y for every n > 0 and
    # so can never equal 1.2 y. That confirms the columns are
    # approximate annotation. See vignette "Errata".
    cyp2b6_max <- fixed(1)
    label("CYP2B6 maximum adult relative expression AdultMax (fraction of adult)")  # Supplementary Table S1 CYP2B6 equation numerator leading constant
    cyp2b6_birth <- fixed(0.1)
    label("CYP2B6 relative expression at birth Fbirth (fraction of adult)")  # Supplementary Table S1 CYP2B6 equation, subtracted and added constant. The row's summary column reads 0.15; see vignette Errata
    cyp2b6_t50 <- fixed(1)
    label("CYP2B6 age at half-maximal maturation Agemid (year)")  # Supplementary Table S1 CYP2B6 equation denominator base. The row's summary column reads 1.2 y; see vignette Errata
    cyp2b6_hill <- fixed(1)
    label("CYP2B6 sigmoidicity exponent n (unitless)")  # Supplementary Table S1 CYP2B6 equation: bare Age (the three superscript boxes in the equation object are empty)

    # ---- Hepatic CYP2C8 ----
    # Supplementary Table S1: CYP2C8 = (1 - 0.3) * Age / (0.02 + Age) + 0.3
    cyp2c8_max <- fixed(1)
    label("CYP2C8 maximum adult relative expression AdultMax (fraction of adult)")  # Supplementary Table S1 CYP2C8 equation
    cyp2c8_birth <- fixed(0.3)
    label("CYP2C8 relative expression at birth Fbirth (fraction of adult)")  # Supplementary Table S1 CYP2C8 equation and column (0.3)
    cyp2c8_t50 <- fixed(0.02)
    label("CYP2C8 age at half-maximal maturation Agemid (year)")  # Supplementary Table S1 CYP2C8 equation denominator
    cyp2c8_hill <- fixed(1)
    label("CYP2C8 sigmoidicity exponent n (unitless)")  # Supplementary Table S1 CYP2C8 equation: bare Age

    # ---- Hepatic CYP2C9 ----
    # Supplementary Table S1:
    #   CYP2C9 = (1 - 0.17) * Age^0.53 / (0.016^0.53 + Age^0.53) + 0.17
    cyp2c9_max <- fixed(1)
    label("CYP2C9 maximum adult relative expression AdultMax (fraction of adult)")  # Supplementary Table S1 CYP2C9 equation
    cyp2c9_birth <- fixed(0.17)
    label("CYP2C9 relative expression at birth Fbirth (fraction of adult)")  # Supplementary Table S1 CYP2C9 equation and column (0.17)
    cyp2c9_t50 <- fixed(0.016)
    label("CYP2C9 age at half-maximal maturation Agemid (year)")  # Supplementary Table S1 CYP2C9 equation denominator
    cyp2c9_hill <- fixed(0.53)
    label("CYP2C9 sigmoidicity exponent n (unitless)")  # Supplementary Table S1 CYP2C9 equation exponent

    # ---- Hepatic CYP2C18/19 ----
    # Supplementary Table S1:
    #   CYP2C19 = (1 - 0.3) * Age^2.44 / (0.28^2.44 + Age^2.44) + 0.3
    # The Table S1 row is labelled "CYP2C18/19"; the equation object is
    # labelled "CYP2C19". Parameters are named cyp2c19 here.
    cyp2c19_max <- fixed(1)
    label("CYP2C18/19 maximum adult relative expression AdultMax (fraction of adult)")  # Supplementary Table S1 CYP2C18/19 equation
    cyp2c19_birth <- fixed(0.3)
    label("CYP2C18/19 relative expression at birth Fbirth (fraction of adult)")  # Supplementary Table S1 CYP2C18/19 equation and column (0.3)
    cyp2c19_t50 <- fixed(0.28)
    label("CYP2C18/19 age at half-maximal maturation Agemid (year)")  # Supplementary Table S1 CYP2C18/19 equation denominator
    cyp2c19_hill <- fixed(2.44)
    label("CYP2C18/19 sigmoidicity exponent n (unitless)")  # Supplementary Table S1 CYP2C18/19 equation exponent

    # ---- Hepatic CYP2D6 ----
    # Supplementary Table S1: CYP2D6 = (1.0 - 0.036) * Age / (0.1 + Age) + 0.036
    cyp2d6_max <- fixed(1)
    label("CYP2D6 maximum adult relative expression AdultMax (fraction of adult)")  # Supplementary Table S1 CYP2D6 equation (printed 1.0)
    cyp2d6_birth <- fixed(0.036)
    label("CYP2D6 relative expression at birth Fbirth (fraction of adult)")  # Supplementary Table S1 CYP2D6 equation and column (0.036); Figure 1 annotates "CYP2D6 Birth level= 0.04"
    cyp2d6_t50 <- fixed(0.1)
    label("CYP2D6 age at half-maximal maturation Agemid (year)")  # Supplementary Table S1 CYP2D6 equation denominator
    cyp2d6_hill <- fixed(1)
    label("CYP2D6 sigmoidicity exponent n (unitless)")  # Supplementary Table S1 CYP2D6 equation: bare Age

    # ---- Hepatic CYP2E1 ----
    # Supplementary Table S1:
    #   CYP2E1 = (1.074 - 0.086) * Age^0.496 / (0.226^0.496 + Age^0.496) + 0.086
    cyp2e1_max <- fixed(1.074)
    label("CYP2E1 maximum adult relative expression AdultMax (fraction of adult)")  # Supplementary Table S1 CYP2E1 equation
    cyp2e1_birth <- fixed(0.086)
    label("CYP2E1 relative expression at birth Fbirth (fraction of adult)")  # Supplementary Table S1 CYP2E1 equation and column (0.086)
    cyp2e1_t50 <- fixed(0.226)
    label("CYP2E1 age at half-maximal maturation Agemid (year)")  # Supplementary Table S1 CYP2E1 equation denominator
    cyp2e1_hill <- fixed(0.496)
    label("CYP2E1 sigmoidicity exponent n (unitless)")  # Supplementary Table S1 CYP2E1 equation exponent

    # ---- Hepatic CYP3A4/5 ----
    # Supplementary Table S1:
    #   CYP3A4/5 = 1.061 * Age^0.78 / (0.66^0.78 + Age^0.78)
    # No additive birth term is printed, i.e. Fbirth = 0, matching the
    # Table S1 column ("Negligible") and the Figure 1 annotation
    # "CYP3A4 Birth level= 0".
    cyp3a4_max <- fixed(1.061)
    label("Hepatic CYP3A4/5 maximum adult relative expression AdultMax (fraction of adult)")  # Supplementary Table S1 hepatic CYP3A4/5 equation
    cyp3a4_birth <- fixed(0)
    label("Hepatic CYP3A4/5 relative expression at birth Fbirth (fraction of adult)")  # Supplementary Table S1 hepatic CYP3A4/5 equation has no additive term; column reads "Negligible"; Figure 1 annotates "CYP3A4 Birth level= 0"
    cyp3a4_t50 <- fixed(0.66)
    label("Hepatic CYP3A4/5 age at half-maximal maturation Agemid (year)")  # Supplementary Table S1 hepatic CYP3A4/5 equation denominator
    cyp3a4_hill <- fixed(0.78)
    label("Hepatic CYP3A4/5 sigmoidicity exponent n (unitless)")  # Supplementary Table S1 hepatic CYP3A4/5 equation exponent

    # ---- Enterocytic (gut) CYP3A4/5 ----
    # Supplementary Table S1:
    #   CYP3A4/5 (gut) = (1 - 0.42) * Age / (2.357 + Age) + 0.42
    cyp3a4gut_max <- fixed(1)
    label("Gut CYP3A4/5 maximum adult relative expression AdultMax (fraction of adult)")  # Supplementary Table S1 gut CYP3A4/5 equation
    cyp3a4gut_birth <- fixed(0.42)
    label("Gut CYP3A4/5 relative expression at birth Fbirth (fraction of adult)")  # Supplementary Table S1 gut CYP3A4/5 equation and column (0.42)
    cyp3a4gut_t50 <- fixed(2.357)
    label("Gut CYP3A4/5 age at half-maximal maturation Agemid (year)")  # Supplementary Table S1 gut CYP3A4/5 equation denominator
    cyp3a4gut_hill <- fixed(1)
    label("Gut CYP3A4/5 sigmoidicity exponent n (unitless)")  # Supplementary Table S1 gut CYP3A4/5 equation: bare Age

    # ---- Renal function ----
    # Supplementary Table S1 (bottom row):
    #   Renal function = ((-0.61604 * BSA^2) + (99.054 * BSA) - 17.74) / 120
    # A quadratic in body surface area (m^2) giving GFR in mL/min,
    # divided by the adult GFR of 120 mL/min to put it on the same
    # fraction-of-adult scale as the CYP pathways.
    renal_bsa2 <- fixed(-0.61604)
    label("Quadratic coefficient of body surface area in the renal GFR polynomial (mL/min per m^4)")  # Supplementary Table S1 renal equation
    renal_bsa1 <- fixed(99.054)
    label("Linear coefficient of body surface area in the renal GFR polynomial (mL/min per m^2)")  # Supplementary Table S1 renal equation
    renal_int <- fixed(-17.74)
    label("Intercept of the renal GFR polynomial (mL/min)")  # Supplementary Table S1 renal equation
    renal_adult <- fixed(120)
    label("Adult glomerular filtration rate used to normalise the renal polynomial (mL/min)")  # Supplementary Table S1 renal equation denominator
  })

  model({
    # The rxode2 time variable IS postnatal age in years (paper symbol
    # "Age"). Solve over 0 to 20 to trace the published range; Salem
    # 2013 evaluates the functions "with intervals of 0.1 year up to 20
    # years of age" (Methods, Equation 4 description).
    age_y <- time

    # ---- Hepatic CYP maturation, Equation 1 ----
    ont_cyp1a2 <- cyp1a2_birth +
      (cyp1a2_max - cyp1a2_birth) * age_y^cyp1a2_hill /
        (cyp1a2_t50^cyp1a2_hill + age_y^cyp1a2_hill)
    ont_cyp2b6 <- cyp2b6_birth +
      (cyp2b6_max - cyp2b6_birth) * age_y^cyp2b6_hill /
        (cyp2b6_t50^cyp2b6_hill + age_y^cyp2b6_hill)
    ont_cyp2c8 <- cyp2c8_birth +
      (cyp2c8_max - cyp2c8_birth) * age_y^cyp2c8_hill /
        (cyp2c8_t50^cyp2c8_hill + age_y^cyp2c8_hill)
    ont_cyp2c9 <- cyp2c9_birth +
      (cyp2c9_max - cyp2c9_birth) * age_y^cyp2c9_hill /
        (cyp2c9_t50^cyp2c9_hill + age_y^cyp2c9_hill)
    ont_cyp2c19 <- cyp2c19_birth +
      (cyp2c19_max - cyp2c19_birth) * age_y^cyp2c19_hill /
        (cyp2c19_t50^cyp2c19_hill + age_y^cyp2c19_hill)
    ont_cyp2d6 <- cyp2d6_birth +
      (cyp2d6_max - cyp2d6_birth) * age_y^cyp2d6_hill /
        (cyp2d6_t50^cyp2d6_hill + age_y^cyp2d6_hill)
    ont_cyp2e1 <- cyp2e1_birth +
      (cyp2e1_max - cyp2e1_birth) * age_y^cyp2e1_hill /
        (cyp2e1_t50^cyp2e1_hill + age_y^cyp2e1_hill)
    ont_cyp3a4 <- cyp3a4_birth +
      (cyp3a4_max - cyp3a4_birth) * age_y^cyp3a4_hill /
        (cyp3a4_t50^cyp3a4_hill + age_y^cyp3a4_hill)

    # ---- Enterocytic (gut) CYP3A4/5 maturation, Equation 1 ----
    ont_cyp3a4gut <- cyp3a4gut_birth +
      (cyp3a4gut_max - cyp3a4gut_birth) * age_y^cyp3a4gut_hill /
        (cyp3a4gut_t50^cyp3a4gut_hill + age_y^cyp3a4gut_hill)

    # ---- Renal function ----
    # Body-surface-area quadratic normalised to the adult GFR. Unlike
    # the CYP outputs this is a function of BSA, not of age; the paper
    # supplies no age-to-BSA growth function, so the user must provide
    # BSA (see covariateData$BSA).
    ont_renal <- (renal_bsa2 * BSA^2 + renal_bsa1 * BSA + renal_int) /
      renal_adult
  })
}
