Fasanmade_2009_infliximab <- function() {
  description <- "Two-compartment population PK model of infliximab (anti-TNF-alpha) in patients with ulcerative colitis (Fasanmade 2009)"
  reference <- "Fasanmade AA, Adedokun OJ, Blank M, Zhou H, Davis HM. Population pharmacokinetic analysis of infliximab in patients with ulcerative colitis. Eur J Clin Pharmacol. 2009;65(12):1211-1228. doi:10.1007/s00228-009-0718-4"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on Vc; normalized as WT/77 per Fasanmade 2009 (reference: 77 kg).",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL; normalized as ALB/4.1 per Fasanmade 2009 (reference: 4.1 g/dL).",
      source_name        = "ALB"
    ),
    ADA_POS = list(
      description        = "Anti-drug antibody positivity (antibodies to infliximab)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative)",
      notes              = "Time-invariant in the source analysis: a subject is coded 1 if antibodies were detected at any visit through week 42/54, otherwise 0. Source paper uses the column label 'ATI' (antibodies to infliximab); renamed to the canonical ADA_POS per covariate-columns.md.",
      source_name        = "ATI"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Source paper uses 'SEX' with 1 = female / 0 = male, matching the canonical SEXF encoding; column renamed to SEXF per covariate-columns.md.",
      source_name        = "SEX"
    )
  )

  population <- list(
    n_subjects     = 482L,
    n_studies      = 2L,
    age_range      = "18-75 years (adults)",
    age_median     = "41 years",
    weight_range   = "36-159 kg",
    weight_median  = "77 kg",
    sex_female_pct = 39.2,
    race_ethnicity = "Predominantly White; smaller proportions of Black, Asian, and Other. Race was not retained as a covariate in the final model.",
    disease_state  = "Moderate-to-severe ulcerative colitis (pooled from the ACT 1 and ACT 2 phase III trials).",
    dose_range     = "5 mg/kg and 10 mg/kg IV infusion over 2 hours at weeks 0, 2, 6, then every 8 weeks through week 46.",
    regions        = "Multi-regional (North America, Europe, and others).",
    notes          = "Immunogenicity incidence (ATI-positive at any time through week 42/54) was ~6.8% across pooled ACT 1/ACT 2 subjects. Concomitant immunomodulators (azathioprine, 6-mercaptopurine) did not have a statistically significant effect on infliximab CL in the final model. Reference covariate values: WT = 77 kg, ALB = 4.1 g/dL, male, ADA-negative."
  )

  ini({
    # Structural parameters (typical values for 77-kg male, albumin 4.1 g/dL,
    # ADA-negative) from Fasanmade 2009 Table 3 (final population PK model).
    lcl <- log(0.407); label("Clearance (CL, L/day)")               # Fasanmade 2009 Table 3
    lvc <- log(3.29);  label("Central volume of distribution (Vc, L)")   # Fasanmade 2009 Table 3
    lvp <- log(4.13);  label("Peripheral volume of distribution (Vp, L)")# Fasanmade 2009 Table 3
    lq  <- log(7.14);  label("Intercompartmental clearance (Q, L/day)")  # Fasanmade 2009 Table 3

    # Covariate effects from Fasanmade 2009 Table 3.
    e_alb_cl <- -1.54;  label("Power exponent of albumin on CL (unitless)")              # Fasanmade 2009 Table 3
    e_ada_cl <-  0.471; label("Fractional change in CL for ADA-positive (unitless)")     # Fasanmade 2009 Table 3
    e_sex_cl <- -0.236; label("Fractional change in CL for females (unitless)")          # Fasanmade 2009 Table 3
    e_wt_vc  <-  0.538; label("Allometric exponent of body weight on Vc (unitless)")     # Fasanmade 2009 Table 3
    e_sex_vc <- -0.137; label("Fractional change in Vc for females (unitless)")          # Fasanmade 2009 Table 3

    # Inter-individual variability (diagonal). omega^2 on the log scale;
    # sqrt(exp(omega^2) - 1) gives the implied %CV.
    etalcl ~ 0.131   # implies ~37.7% CV; Fasanmade 2009 Table 3
    etalvc ~ 0.048   # implies ~22.1% CV; Fasanmade 2009 Table 3

    # Residual error (combined additive + proportional); Fasanmade 2009 Table 3.
    propSd <- 0.403;  label("Proportional residual error (fraction)")
    addSd  <- 0.0413; label("Additive residual error (ug/mL)")
  })
  model({
    # Individual PK parameters. Reference subject: 77-kg male, albumin 4.1 g/dL,
    # ADA-negative. Covariate forms per Fasanmade 2009 Table 3:
    #   - ALB on CL: power model (ALB/4.1)^e_alb_cl
    #   - ADA on CL, SEX on CL/Vc: multiplicative (1 + e * COV)
    #   - WT on Vc: power model (WT/77)^e_wt_vc
    cl <- exp(lcl + etalcl) *
      (ALB / 4.1)^e_alb_cl *
      (1 + e_ada_cl * ADA_POS) *
      (1 + e_sex_cl * SEXF)

    vc <- exp(lvc + etalvc) *
      (WT / 77)^e_wt_vc *
      (1 + e_sex_vc * SEXF)

    vp <- exp(lvp)
    q  <- exp(lq)

    # Two-compartment model (IV infusion administered to the central compartment).
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL
    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
