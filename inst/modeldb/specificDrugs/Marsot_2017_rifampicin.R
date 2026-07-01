Marsot_2017_rifampicin <- function() {
  description <- paste(
    "One-compartment population pharmacokinetic model for oral",
    "rifampicin in adult patients with staphylococcal osteoarticular",
    "infections (Marsot 2017; 62 patients, 103 steady-state plasma",
    "concentrations from routine therapeutic drug monitoring at",
    "300 mg three times daily). Absorption uses a single-transit",
    "compartment chain (depot -> transit1 -> central) with a fixed",
    "first-order absorption rate constant ka = 1.15 1/h; first-order",
    "elimination from the central compartment. Coadministration of",
    "oral fusidic acid (500 mg three times daily in 16 % of the",
    "cohort) was retained as the sole significant covariate and",
    "reduces both apparent oral clearance and apparent central",
    "volume of distribution: typical CL/F is 5.1 L/h (with fusidic",
    "acid) vs 13.7 L/h (without) and typical V/F is 23.8 L (with)",
    "vs 61.1 L (without), an interaction attributed to CYP3A4",
    "inhibition and altered plasma protein binding by fusidic acid.",
    "Inter-individual variability is reported on CL/F (72.9 % CV)",
    "and V/F (59.1 % CV); residual variability is additive with",
    "standard deviation 2.256 mg/L."
  )
  reference <- paste(
    "Marsot A, Menard A, Dupouey J, Muziotti C, Guilhaumou R, Blin O.",
    "Population pharmacokinetics of rifampicin in adult patients with",
    "osteoarticular infections: interaction with fusidic acid.",
    "Br J Clin Pharmacol. 2017 May;83(5):1039-1047.",
    "doi:10.1111/bcp.13178."
  )
  vignette <- "Marsot_2017_rifampicin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CONMED_FUSIDIC = list(
      description        = paste(
        "1 = patient is coadministered oral fusidic acid (500 mg",
        "three times daily) alongside rifampicin during the",
        "observation interval; 0 = no concomitant fusidic acid.",
        "Time-fixed per subject in the Marsot 2017 cohort (each",
        "patient's antibiotic regimen was stable across the",
        "steady-state sampling window)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant fusidic acid)",
      notes              = paste(
        "Marsot 2017 Table 2 encodes fusidic acid coadministration",
        "as separate typical values for CL/F and V/F rather than as",
        "a multiplicative theta: CL/F = 5.1 L/h with fusidic acid,",
        "13.7 L/h without; V/F = 23.8 L with, 61.1 L without.",
        "The model encodes this as multiplicative effects",
        "e_conmed_fusidic_cl = log(5.1 / 13.7) and",
        "e_conmed_fusidic_vc = log(23.8 / 61.1) applied on the log",
        "scale so the CONMED_FUSIDIC = 0 reference recovers the",
        "13.7 L/h and 61.1 L typical values entered as lcl and",
        "lvc. Fusidic acid coadministration frequency in the",
        "studied cohort is 10 / 62 = 16.1 % (Marsot 2017 Table 1).",
        "Marsot 2017 Discussion attributes the interaction to",
        "fusidic acid inhibition of CYP3A4 and to displacement of",
        "rifampicin from plasma protein binding sites, but the",
        "exact mechanism remains uncertain."
      ),
      source_name        = "fusidic acid (Table 1 / Table 2)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at study entry (years). Tested in the covariate screen but not retained in the final model.",
      units       = "year",
      type        = "continuous",
      notes       = "Marsot 2017 Results paragraph 3: 'The other covariables tested on the model, including gender, age, height and body weight, did not improve the pharmacokinetic model.' Cohort age range 20-89 years, mean 57.4 (Table 1)."
    ),
    WT = list(
      description = "Total body weight at study entry (kg). Tested in the covariate screen but not retained in the final model.",
      units       = "kg",
      type        = "continuous",
      notes       = "Marsot 2017 Results paragraph 3 (see AGE notes). Cohort weight range 46-119 kg, mean 72.3 (Table 1). No allometric scaling in the final model."
    ),
    HT = list(
      description = "Height at study entry (cm). Tested in the covariate screen but not retained in the final model.",
      units       = "cm",
      type        = "continuous",
      notes       = "Marsot 2017 Results paragraph 3 (see AGE notes)."
    ),
    SEXF = list(
      description = "1 = female, 0 = male. Tested in the covariate screen but not retained in the final model.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Marsot 2017 Results paragraph 3 (see AGE notes). Female proportion 16 / 62 = 25.8 % (Table 1)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 62L,
    n_studies      = 1L,
    n_observations = 103L,
    age_range      = "20-89 years; mean 57.4 (Table 1)",
    weight_range   = "46-119 kg; mean 72.3 (Table 1)",
    sex_female_pct = 25.8,
    race_ethnicity = "not reported",
    disease_state  = paste(
      "Adult outpatients with staphylococcal osteoarticular",
      "infections (OAI) followed at the Infectious and Tropical",
      "Diseases Unit of Conception Hospital (Marseille, France).",
      "75.8 % of patients had infections associated with prosthetic",
      "materials (prosthesis 51.6 %, osteosynthesis 25.8 %; 21.0 %",
      "no material). Main suspected pathogens: Staphylococcus",
      "aureus (64.5 %), coagulase-negative staphylococci (35.5 %),",
      "Streptococcus spp. (6.5 %). Only steady-state samples",
      "(less than 4 weeks after treatment start) were analysed."
    ),
    dose_range     = paste(
      "Oral rifampicin 300 mg three times daily (median 12.4",
      "mg/kg). Treatment duration: 3 months (no prosthetic",
      "material) or 6 months (prosthetic material) per Marsot 2017",
      "Materials and method 'Data sources'."
    ),
    regions        = "France (Marseille)",
    notes          = paste(
      "Retrospective TDM cohort (August 2012-August 2015). Plasma",
      "rifampicin quantified by HPLC-UV (LLOQ 0.5 mg/L, linear",
      "range 0.5-20 mg/L). Concomitant antibiotics beyond",
      "rifampicin (n, % of cohort): ofloxacin 32 (51.6),",
      "fusidic acid 10 (16.1), clindamycin 6 (9.7), teicoplanin 6",
      "(9.7), ciprofloxacin 5 (8.1), vancomycin 5 (8.1),",
      "amoxicillin 5 (8.1), cotrimoxazole 2 (3.2), ceftazidime 2",
      "(3.2), cloxacillin 1 (1.6). Only fusidic acid was retained",
      "as a significant covariate in the final model (Marsot 2017",
      "Table 1 and Results paragraph 2)."
    )
  )

  ini({
    # ============================================================
    # Structural PK -- one-compartment with a single-transit
    # absorption chain (depot -> transit1 -> central) and
    # first-order elimination. Point estimates from Marsot 2017
    # Table 2 ('Estimate Mean' column, final model). The CL/F and
    # V/F typical values entered here are the reference category
    # 'without fusidic acid' values; the effect of fusidic acid
    # coadministration is applied as the multiplicative
    # e_conmed_fusidic_cl / e_conmed_fusidic_vc log-scale offsets
    # below so that CONMED_FUSIDIC = 0 recovers Table 2's
    # 'Without fusidic acid' values (13.7 L/h and 61.1 L).
    # ============================================================
    lcl <- log(13.7);  label("Apparent oral clearance without fusidic acid (CL/F, L/h)")               # Marsot 2017 Table 2 row "CL/F (L h-1) Without fusidic acid" = 13.7 (RSE 26.3 %)
    lvc <- log(61.1);  label("Apparent central volume of distribution without fusidic acid (V/F, L)")  # Marsot 2017 Table 2 row "V/F (L)   Without fusidic acid" = 61.1 (RSE 56.5 %)

    # Absorption rate constant fixed at 1.15 1/h. Marsot 2017
    # Results paragraph 2: "The absorption rate was fixed at 1.15
    # h^-1 in our study" and Discussion: "In the presented model,
    # the absorption rate constant was fixed at 1.15 h^-1 to
    # improve the fit, as presented in previous studies [14, 28]."
    # This is the transit-chain rate constant applied at both
    # depot -> transit1 and transit1 -> central transitions.
    lka <- fixed(log(1.15));  label("Absorption / transit-chain rate constant (ka, 1/h; fixed)")

    # Bioavailability anchor. Oral-only data; absolute F is not
    # identifiable, so the typical value is fixed at 1.0. Marsot
    # 2017 parameterises exclusively in the CL/F and V/F apparent
    # forms and does not estimate F.
    lfdepot <- fixed(log(1));  label("Relative bioavailability anchor (F, fixed at 1)")

    # ============================================================
    # Covariate effect of fusidic acid coadministration. Marsot
    # 2017 Table 2 reports separate typical values for the two
    # strata rather than a multiplicative theta:
    #   CL/F with fusidic acid    = 5.1  L/h
    #   CL/F without fusidic acid = 13.7 L/h  (reference)
    #   V/F  with fusidic acid    = 23.8 L
    #   V/F  without fusidic acid = 61.1 L    (reference)
    # Encoded here as log-scale offsets so that
    #   exp(lcl + e_conmed_fusidic_cl) = 13.7 * exp(log(5.1/13.7)) = 5.1
    #   exp(lvc + e_conmed_fusidic_vc) = 61.1 * exp(log(23.8/61.1)) = 23.8
    # ============================================================
    e_conmed_fusidic_cl <- log(5.1 / 13.7);   label("Power-form effect of CONMED_FUSIDIC on CL/F (log with/without ratio, unitless)")   # Marsot 2017 Table 2 CL/F 'With' 5.1 vs 'Without' 13.7
    e_conmed_fusidic_vc <- log(23.8 / 61.1);  label("Power-form effect of CONMED_FUSIDIC on V/F  (log with/without ratio, unitless)")   # Marsot 2017 Table 2 V/F  'With' 23.8 vs 'Without' 61.1

    # ============================================================
    # Inter-individual variability (natural-log-scale variance,
    # per Marsot 2017 Table 2 footnote "omega, variance of eta_ki"
    # and the exponential IIV model P_i = theta_k * exp(eta_ki)
    # described in Materials and method 'Population
    # pharmacokinetic analysis'). Consistency check against the
    # abstract's summary "72.9 % (49.5, 86.0 %) and 59.1 %
    # (5.5, 105.4 %)" for CL/F and V/F respectively:
    #   sqrt(0.531) = 0.729 -> 72.9 % CV  (matches abstract)
    #   sqrt(0.349) = 0.591 -> 59.1 % CV  (matches abstract)
    # So the reported values 0.531 and 0.349 are variances on
    # the log scale and translate directly to CV via sqrt().
    # No off-diagonal correlations are reported.
    # ============================================================
    etalcl ~ 0.531
    # Marsot 2017 Table 2 row "omega CL/F" = 0.531 (RSE 38.3 %); sqrt(0.531) = 0.729 = 72.9 % CV.

    etalvc ~ 0.349
    # Marsot 2017 Table 2 row "omega V/F"  = 0.349 (RSE 238.1 %; poorly estimated but retained in the final model per Marsot 2017 Discussion paragraph 5). sqrt(0.349) = 0.591 = 59.1 % CV.

    # ============================================================
    # Residual unexplained variability. Marsot 2017 Materials and
    # method 'Population pharmacokinetic analysis': "Ultimately,
    # an additive error model was used to model random residual
    # variability according to: C_ij = CP_ij + epsilon_ij". Table 2
    # reports sigma_add = 2.256 mg/L (RSE 12.4 %), consistent with
    # the abstract's "Residual variability was 2.3 mg l^-1 (1.6,
    # 2.6 mg l^-1)". The 2.256 mg/L value is the additive-error
    # standard deviation in linear concentration units (mg/L);
    # the Table 2 footnote's "sigma, variance of epsilon_ij"
    # description is inconsistent with the reported units and
    # with the abstract's linear-space variability -- the value
    # is used here as an SD in mg/L (see vignette Assumptions
    # and deviations for the reconciliation). The Results
    # paragraph 2 description "additive model on natural
    # log-transformed data" appears to describe the goodness-of-
    # fit diagnostic used to select the residual, not the
    # linear-space equation implemented in NONMEM.
    # ============================================================
    addSd <- 2.256;  label("Additive residual error standard deviation (mg/L)")  # Marsot 2017 Table 2 row "sigma_add (mg l-1)" = 2.256 (RSE 12.4 %)
  })

  model({
    # Individual PK parameters. Fusidic acid coadministration
    # (CONMED_FUSIDIC = 1) applies the log-scale offset that
    # multiplies the typical CL/F by 5.1/13.7 and the typical
    # V/F by 23.8/61.1, exactly reproducing Table 2's stratified
    # typical values.
    cl <- exp(lcl + e_conmed_fusidic_cl * CONMED_FUSIDIC + etalcl)
    vc <- exp(lvc + e_conmed_fusidic_vc * CONMED_FUSIDIC + etalvc)
    ka <- exp(lka)

    # Micro-constant
    kel <- cl / vc

    # ODE system: oral single-transit absorption chain
    # (depot -> transit1 -> central), both transitions at rate ka
    # per Marsot 2017 Methods (transit compartment feeding an
    # 'absorption compartment' at rate ka). See also Svensson 2016
    # rifampicin for the same Savic-style N = 1 transit-chain
    # implementation.
    d/dt(depot)    <- -ka * depot
    d/dt(transit1) <-  ka * depot    - ka * transit1
    d/dt(central)  <-  ka * transit1 - kel * central

    # Bioavailability anchor (F fixed at 1 for oral-only apparent
    # CL/V parameterisation).
    f(depot) <- exp(lfdepot)

    # Plasma concentration of rifampicin (mg/L).
    Cc <- central / vc

    # Additive residual error (mg/L) per Marsot 2017 Methods:
    # C_ij = CP_ij + epsilon_ij.
    Cc ~ add(addSd)
  })
}
