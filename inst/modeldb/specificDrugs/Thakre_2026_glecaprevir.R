Thakre_2026_glecaprevir <- function() {
  description <- "One-compartment population PK model with first-order absorption and first-order elimination for oral glecaprevir (given with pibrentasvir) in pediatric patients aged 3 to <18 years with chronic hepatitis C virus genotype 1-6 infection (DORA study, Thakre 2026). Body weight enters CL/F and V2/F allometrically at a 70 kg reference with estimated exponents; the absorption rate is parameterised as a multiple of each subject's own elimination rate constant (ka = 3.65 x CL/V2) so that absorption is always faster than elimination. The pediatric coated-granule formulation (Part 2) carries a relative bioavailability of 0.8 versus the adult tablet, fixed from an adult phase 1 bioavailability study. No other covariate was retained. Time is in days."
  reference <- "Thakre N, Beck D, Saeed AM, Eckert D, Chen MJ, Liu W, Lon HK, Mostafa NM, Mensing S. Population Pharmacokinetic Analysis of Glecaprevir and Pibrentasvir in Pediatric Patients >= 3 to < 18 Years of Age with Genotypes 1-6 Chronic Hepatitis C Virus (HCV) Infection. Clin Pharmacokinet. 2026;65(7):1061-1072. doi:10.1007/s40262-026-01637-1"
  vignette <- "Thakre_2026_glecaprevir_pibrentasvir"
  units <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "glecaprevir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "glecaprevir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power (allometric) scaling of CL/F and V2/F centered on the standard 70 kg body weight, with both exponents estimated; part of the structural base model (Thakre 2026 Methods 'Model Development'; ESM Fig. S3 legend: CL/F = exp(TVCL + log(WTKG/70)*WTKGonCL + ETA1)). Analysis population median 36.45 kg, range 12.70-108.90 kg (ESM Table S1).",
      source_name = "WTKG"
    ),
    FORM_GRANULE = list(
      description = "Pediatric coated-granule formulation indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (adult co-formulated GLE/PIB 100 mg/40 mg film-coated tablet)",
      notes = "1 = DORA Part 2 children (3 to <12 years) dosed with the pediatric coated-granule formulation (separate GLE and PIB coated granules in the intensive-PK portion, or one packet of co-packaged GLE/PIB coated granules mixed with soft food in the non-intensive portion); 0 = DORA Part 1 adolescents dosed with the adult tablet. Multiplies F1 by the fixed relative bioavailability 0.800 (ratio of pediatric to adult formulation; Thakre 2026 Table 2 'Formulation (Part 2) on F1' and Methods 'Model Development'). The paper applies one F1 ratio to every Part 2 dose and does not distinguish the two pediatric granule presentations.",
      source_name = "Formulation (Part 2)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested on CL/F and V2/F as a power function normalised to the study median (10 years) after the allometric weight scaling; not significant (Thakre 2026 Results 'Model Results'; ESM Supplementary Results 1.2.2 and Fig. S3)."
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      notes = "Tested as an allometric alternative to body weight (reference 1.7 m^2; analogous results, weight kept as more convenient) and as a power covariate normalised to the study median (1.18 m^2) on top of the weight scaling; not significant (Thakre 2026 Methods and Results; ESM 1.2.3 and Fig. S5)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Biological sex tested as a multiplicative categorical effect on CL/F and V2/F; not significant (Thakre 2026 Methods 'Model Development' and Results 'Model Results')."
    ),
    RACE_WHITE = list(
      description = "White race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Race (white vs others) tested as a multiplicative categorical effect. Race on V2/F was the only covariate significant in forward inclusion but was removed in backward elimination (P > 0.001), so the base and final glecaprevir models are identical (Thakre 2026 Results 'Model Results')."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 126L,
    n_studies = 1L,
    age_range = "3-17 years",
    age_median = "10 years",
    weight_range = "12.70-108.90 kg",
    weight_median = "36.45 kg",
    sex_female_pct = 55.56,
    race_ethnicity = c(White = 70.63, Black = 5.56, Asian = 15.87, Other = 7.94),
    disease_state = "Chronic hepatitis C virus infection, genotypes 1-4 in the analysed patients (GT1 75.4%, GT2 4.0%, GT3 16.7%, GT4 4.0%); no cirrhosis; 3 patients (2.4%) HCV-HIV coinfected; 89.7% treatment naive.",
    dose_range = "Once-daily oral glecaprevir with pibrentasvir. Part 1 adolescents (12 to <18 years): adult tablets GLE/PIB 300 mg/120 mg. Part 2 children: pediatric coated granules dosed by weight band, initial 40/15 mg dose ratio (120/45, 160/60 or 200/75 mg) and final 50/20 mg ratio (150/60 mg for 12 to <20 kg, 200/80 mg for 20 to <30 kg, 250/100 mg for 30 to <45 kg). Treatment 8, 12 or 16 weeks.",
    regions = "Multinational, including Japan (10.3% of patients).",
    notes = "DORA phase 2/3 study (NCT03067129): Part 1 cohort 1 adolescents (n = 47) and Part 2 cohorts 2-4 children aged 9 to <12, 6 to <9 and 3 to <6 years (n = 79 analysed). 693 glecaprevir concentrations from intensive (day 1 at 4 h; week 2 at 0, 2, 4, 6, 12 h) and sparse sampling. Demographics from Thakre 2026 ESM Table S1."
  )

  ini({
    # Structural parameters, typical values for a 70 kg patient. Table 2
    # footnote a: CL/F and V2/F were estimated as log-transformed THETAs and
    # back-transformed for the table, so they enter on the log scale here.
    lcl <- log(1830); label("Apparent clearance CL/F at 70 kg (L/day)") # Thakre 2026 Table 2: CL/F = 1830 L/day (SEE 388; 95% CI 1220-2740)
    lvc <- log(358); label("Apparent central volume V2/F at 70 kg (L)") # Thakre 2026 Table 2: V2/F = 358 L (SEE 97.0; 95% CI 216-593)
    lkaratio <- log(3.65); label("Log ratio of absorption rate ka to elimination rate CL/V2 (unitless)") # Thakre 2026 Table 2: Absorption factor = 3.65 (SEE 0.438), footnote b ka estimated via a factor on ke = CL/V2; derived typical ka 18.6 1/day

    # Allometric body-weight exponents (estimated), reference 70 kg
    e_wt_cl <- 1.04; label("Allometric exponent of body weight on CL/F (unitless)") # Thakre 2026 Table 2: Body weight on CL/F = 1.04 (SEE 0.211)
    e_wt_vc <- 1.19; label("Allometric exponent of body weight on V2/F (unitless)") # Thakre 2026 Table 2: Body weight on V2/F = 1.19 (SEE 0.246)

    # Relative bioavailability of the pediatric granules vs the adult tablet
    e_form_granule_fdepot <- fixed(0.8); label("Relative bioavailability F1 of the pediatric granule formulation vs the adult tablet (ratio)") # Thakre 2026 Table 2: Formulation (Part 2) on F1 = 0.800 FIX; Methods: from an adult phase 1 bioavailability study

    # IIV (variances). The paper states CL/F and V2/F carry correlated random
    # effects but prints no covariance, so the block is encoded as diagonal.
    etalcl ~ 1.39 # Thakre 2026 Table 2: IIV on CL/F variance 1.39 (174 %CV), shrinkage 3.8%
    etalvc ~ 1.52 # Thakre 2026 Table 2: IIV on V2/F variance 1.52 (189 %CV), shrinkage 5.1%

    # Residual error: Table 2 reports sigma^2 = 0.568, so SD = sqrt(0.568).
    propSd <- 0.7537; label("Proportional residual error (fraction)") # Thakre 2026 Table 2: sigma1^2 (proportional) = 0.568 (SEE 0.0320); sqrt(0.568) = 0.7537
  })

  model({
    # Individual parameters; body weight power-scaled at 70 kg
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    kel <- cl / vc
    # Absorption rate constant is a fixed multiple of the subject's own
    # elimination rate constant (Thakre 2026 Table 2 footnote b), which keeps
    # ka > ke and excludes flip-flop kinetics.
    karatio <- exp(lkaratio)
    ka <- karatio * kel

    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - kel * central

    # F1 = 1 for the adult tablet; the pediatric granules carry the fixed ratio
    f(depot) <- 1 + (e_form_granule_fdepot - 1) * FORM_GRANULE

    # mg / L = ug/mL; x 1000 for ng/mL
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
