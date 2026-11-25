Dudkowski_2017_alogliptin <- function() {
  description <- "Two-compartment population PK model with first-order absorption and allometric body-weight scaling for oral alogliptin, with separate pediatric and adult absorption rate, apparent clearance and EC50, linked to a direct sigmoid-Emax model of plasma DPP-4 inhibition, in children, adolescents and adults with type 2 diabetes mellitus (Dudkowski 2017)"
  reference <- paste(
    "Dudkowski C, Tsai M, Liu J, Zhao Z, Schmidt E, Xie J.",
    "The pharmacokinetics and pharmacodynamics of alogliptin in children,",
    "adolescents, and adults with type 2 diabetes mellitus.",
    "Eur J Clin Pharmacol. 2017;73(3):279-288.",
    "doi:10.1007/s00228-016-2175-1"
  )
  vignette <- "Dudkowski_2017_alogliptin"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL (alogliptin plasma concentration); % (DPP-4 inhibition)"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline weight; allometric scaling (WT/70)^0.75 on CL/F and Q/F and (WT/70)^1 on V2/F and V3/F (Methods, PK/PD modeling and simulations).",
      source_name = "WT"
    ),
    CHILD = list(
      description = "Pediatric indicator: 1 = pediatric subject (10 to <18 years, children and adolescents pooled), 0 = adult (18-65 years)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Switches the typical KA, CL/F and EC50 between the pediatric and adult values of Supplemental Table 1. The paper pools its children (10 to <14 y) and adolescents (14 to <18 y) into one pediatric group for the model; no separate adolescent effect was estimated.",
      source_name = "pediatric vs adult subjects (Supplemental Table 1 row labels)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "alogliptin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "alogliptin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "alogliptin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 46L,
    n_studies = 1L,
    studies = "Phase 1 randomized, open-label, single-dose, multicenter study NCT00957268.",
    age_range = "10-65 years (children 10 to <14 y, adolescents 14 to <18 y, adults 18-65 y)",
    age_mean = "12.0-12.4 y (children), 15.1-15.4 y (adolescents), 51.3 y (adults) (Table 1)",
    weight_mean = "86.6-116.3 kg across the five dose/age arms; 92.3 kg in adults (Table 1)",
    sex_female_pct = 73.9,
    race_ethnicity = c(Black = 69.6, White = 30.4),
    disease_state = "Type 2 diabetes mellitus without renal impairment (CrCl >= 50 mL/min) and, in pediatric subjects, without clinically significant comorbidity; 39 of 46 on concomitant metformin.",
    dose_range = "Single oral dose of 12.5 or 25 mg (pediatric, randomized 1:1) or 25 mg (adults), after an overnight fast",
    regions = "United States",
    n_pediatric = 24L,
    n_adult = 22L,
    notes = "Table 1: 9 children (5 at 12.5 mg, 4 at 25 mg), 15 adolescents (8 at 12.5 mg, 7 at 25 mg) and 22 gender- and race-matched adults (25 mg). Female 34/46 and Black 32/46 computed from the Table 1 counts. Regions are not stated in the paper; the US is inferred from the US-CFR compliance statement and Takeda Development Center Americas sponsorship."
  )

  ini({
    # Structural PK parameters (Supplemental Table 1), typical values at WT = 70 kg.
    # The adult values are the reference; the pediatric values are switched in with CHILD.
    lka <- log(1.66)
    label("Absorption rate constant KA, adults (1/h)") # Suppl. Table 1 'KA in adult subjects' = 1.66 (%SEM 41.2)
    lka_ped <- log(0.590)
    label("Absorption rate constant KA, pediatric subjects (1/h)") # Suppl. Table 1 'KA in pediatric subjects' = 0.590 (%SEM 19.0)
    lcl <- log(13.0)
    label("Apparent clearance CL/F, adults, 70 kg (L/h)") # Suppl. Table 1 'CL/F in adult subjects' = 13.0 (%SEM 3.53)
    lcl_ped <- log(14.4)
    label("Apparent clearance CL/F, pediatric subjects, 70 kg (L/h)") # Suppl. Table 1 'CL/F in pediatric subjects' = 14.4 (%SEM 4.86)
    lvc <- log(125)
    label("Apparent central volume V2/F, 70 kg (L)") # Suppl. Table 1 'V2/F' = 125 (%SEM 8.83)
    lvp <- log(117)
    label("Apparent peripheral volume V3/F, 70 kg (L)") # Suppl. Table 1 'V3/F' = 117 (%SEM 7.16)
    lq <- log(11.7)
    label("Apparent intercompartmental clearance Q/F, 70 kg (L/h)") # Suppl. Table 1 'Q/F' = 11.7 (%SEM 21.3)

    # Allometric exponents, fixed (Methods equation: b = 0.75 for CL/F, b = 1 for volume)
    e_wt_cl_q <- fixed(0.75)
    label("Allometric exponent on CL/F and Q/F (unitless)") # Methods, PK/PD modeling: 'b = 0.75 for CL/F'
    e_wt_vc_vp <- fixed(1)
    label("Allometric exponent on V2/F and V3/F (unitless)") # Methods, PK/PD modeling: 'b = 1 for Vz/F'

    # PD: sigmoid Emax on DPP-4 inhibition (%), driven by plasma alogliptin (ng/mL)
    e0 <- fixed(0)
    label("Baseline DPP-4 inhibition E0 (%)") # Suppl. Table 1 'E0' = 0 FIXED; Methods 'fixed to 0'
    emax <- 99.3
    label("Maximal DPP-4 inhibition Emax (%)") # Suppl. Table 1 'Emax' = 99.3 (%SEM 1.78)
    lec50 <- log(5.45)
    label("Alogliptin concentration at half-maximal inhibition EC50, adults (ng/mL)") # Suppl. Table 1 'EC50 in adult subjects' = 5.45 (%SEM 8.13)
    lec50_ped <- log(6.38)
    label("Alogliptin concentration at half-maximal inhibition EC50, pediatric subjects (ng/mL)") # Suppl. Table 1 'EC50 in pediatric subjects' = 6.38 (%SEM 8.64)
    lhill <- log(0.857)
    label("Sigmoid Emax shape factor gamma (unitless)") # Suppl. Table 1 'Shape factor' = 0.857 (%SEM 10.1)

    # IIV (Suppl. Table 1). %CV is read as 100*sqrt(omega^2), the convention the same
    # table uses for its residual row (variance 0.0820 printed as 28.6 %CV).
    # KA: 'NE'; V3/F, Q/F, Emax and gamma: '0 %CV FIXED' -> no eta.
    # One eta on CL/F and one on EC50 are shared by both age groups (same %CV and %SEM printed for both rows).
    etalcl ~ 0.034225 # Suppl. Table 1 CL/F IIV 18.5 %CV -> 0.185^2
    etalvc ~ 0.090601 # Suppl. Table 1 V2/F IIV 30.1 %CV -> 0.301^2
    etalec50 ~ 0.060516 # Suppl. Table 1 EC50 IIV 24.6 %CV -> 0.246^2

    # Residual error (Suppl. Table 1; values printed as variances)
    propSd <- 0.28636
    label("Proportional residual error on alogliptin concentration (fraction)") # Suppl. Table 1 'PK Proportional Residual Variability' = 0.0820 (28.6 %CV); sqrt(0.0820)
    addSd <- fixed(0.1)
    label("Additive residual error on alogliptin concentration (ng/mL)") # Suppl. Table 1 'PK Additive Residual Variability' = 0.0100 FIXED (0.100 SD)
    addSd_Dpp4Inh <- 1.9157
    label("Additive residual error on DPP-4 inhibition (%)") # Suppl. Table 1 'PD Additive Residual Variability' = 3.67 (1.92 SD); sqrt(3.67)
  })

  model({
    # Typical values switch between adult (CHILD = 0) and pediatric (CHILD = 1) estimates
    ka <- exp(lka * (1 - CHILD) + lka_ped * CHILD)
    cl <- exp(lcl * (1 - CHILD) + lcl_ped * CHILD + etalcl) * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    vp <- exp(lvp) * (WT / 70)^e_wt_vc_vp
    q <- exp(lq) * (WT / 70)^e_wt_cl_q
    ec50 <- exp(lec50 * (1 - CHILD) + lec50_ped * CHILD + etalec50)
    hill <- exp(lhill)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # mg / L -> ng/mL
    Cc <- 1000 * central / vc

    # Methods: E = E0 + Emax * Conc^gamma / (EC50^gamma + Conc^gamma)
    Dpp4Inh <- e0 + emax * Cc^hill / (ec50^hill + Cc^hill)

    Cc ~ add(addSd) + prop(propSd)
    Dpp4Inh ~ add(addSd_Dpp4Inh)
  })
}
