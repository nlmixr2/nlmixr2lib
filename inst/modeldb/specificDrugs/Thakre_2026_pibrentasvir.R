Thakre_2026_pibrentasvir <- function() {
  description <- "Two-compartment population PK model with first-order absorption, an absorption lag time and first-order elimination for oral pibrentasvir (given with glecaprevir) in pediatric patients aged 3 to <18 years with chronic hepatitis C virus genotype 1-6 infection (DORA study, Thakre 2026). Body weight enters CL/F and V2/F allometrically at a 70 kg reference with estimated exponents; Q/F and V3/F are not weight-scaled. The absorption rate is parameterised as a multiple of each subject's own elimination rate constant (ka = 2.47 x CL/V2) so that absorption is always faster than elimination. The pediatric coated-granule formulation (Part 2) carries a relative bioavailability of 1.2 versus the adult tablet, fixed from an adult phase 1 bioavailability study. No other covariate was retained. Time is in days."
  reference <- "Thakre N, Beck D, Saeed AM, Eckert D, Chen MJ, Liu W, Lon HK, Mostafa NM, Mensing S. Population Pharmacokinetic Analysis of Glecaprevir and Pibrentasvir in Pediatric Patients >= 3 to < 18 Years of Age with Genotypes 1-6 Chronic Hepatitis C Virus (HCV) Infection. Clin Pharmacokinet. 2026;65(7):1061-1072. doi:10.1007/s40262-026-01637-1"
  vignette <- "Thakre_2026_glecaprevir_pibrentasvir"
  units <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "pibrentasvir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pibrentasvir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "pibrentasvir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power (allometric) scaling of CL/F and V2/F centered on the standard 70 kg body weight, with both exponents estimated; part of the structural base model (Thakre 2026 Methods 'Model Development'; ESM Fig. S4 legend: CL/F = exp(TVCL + log(WTKG/70)*WTKGonCL + ETA1)). Q/F and V3/F carry no weight term (Table 3). Analysis population median 36.45 kg, range 12.70-108.90 kg (ESM Table S1).",
      source_name = "WTKG"
    ),
    FORM_GRANULE = list(
      description = "Pediatric coated-granule formulation indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (adult co-formulated GLE/PIB 100 mg/40 mg film-coated tablet)",
      notes = "1 = DORA Part 2 children (3 to <12 years) dosed with the pediatric coated-granule formulation (separate GLE and PIB coated granules in the intensive-PK portion, or one packet of co-packaged GLE/PIB coated granules mixed with soft food in the non-intensive portion); 0 = DORA Part 1 adolescents dosed with the adult tablet. Multiplies F1 by the fixed relative bioavailability 1.20 (ratio of pediatric to adult formulation; Thakre 2026 Table 3 'Formulation (Part 2) on F1' and Methods 'Model Development'). The paper applies one F1 ratio to every Part 2 dose and does not distinguish the two pediatric granule presentations.",
      source_name = "Formulation (Part 2)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested on CL/F and V2/F as a power function normalised to the study median (10 years) after the allometric weight scaling; not significant (Thakre 2026 Results 'Model Results'; ESM Supplementary Results 1.2.2 and Fig. S4)."
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      notes = "Tested as an allometric alternative to body weight (reference 1.7 m^2; analogous results, weight kept as more convenient) and as a power covariate normalised to the study median (1.18 m^2) on top of the weight scaling; not significant (Thakre 2026 Methods and Results; ESM 1.2.3 and Fig. S6)."
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
      notes = "Race (white vs others) tested as a multiplicative categorical effect; not significant for pibrentasvir (Thakre 2026 Results 'Model Results')."
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
    dose_range = "Once-daily oral pibrentasvir with glecaprevir. Part 1 adolescents (12 to <18 years): adult tablets GLE/PIB 300 mg/120 mg. Part 2 children: pediatric coated granules dosed by weight band, initial 40/15 mg dose ratio (120/45, 160/60 or 200/75 mg) and final 50/20 mg ratio (150/60 mg for 12 to <20 kg, 200/80 mg for 20 to <30 kg, 250/100 mg for 30 to <45 kg). Treatment 8, 12 or 16 weeks.",
    regions = "Multinational, including Japan (10.3% of patients).",
    notes = "DORA phase 2/3 study (NCT03067129): Part 1 cohort 1 adolescents (n = 47) and Part 2 cohorts 2-4 children aged 9 to <12, 6 to <9 and 3 to <6 years (n = 79 analysed). 693 pibrentasvir concentrations from intensive (day 1 at 4 h; week 2 at 0, 2, 4, 6, 12 h) and sparse sampling. Demographics from Thakre 2026 ESM Table S1."
  )

  ini({
    # Structural parameters, typical values for a 70 kg patient. Table 3
    # footnote a: CL/F, V2/F, Q/F and V3/F were estimated as log-transformed
    # THETAs and back-transformed for the table.
    lcl <- log(2170); label("Apparent clearance CL/F at 70 kg (L/day)") # Thakre 2026 Table 3: CL/F = 2170 L/day (SEE 231; 95% CI 1760-2660)
    lvc <- log(528); label("Apparent central volume V2/F at 70 kg (L)") # Thakre 2026 Table 3: V2/F = 528 L (SEE 86.2; 95% CI 386-723)
    lq <- log(765); label("Apparent intercompartmental clearance Q/F (L/day)") # Thakre 2026 Table 3: Q/F = 765 L/day (SEE 92.0; 95% CI 606-966)
    lvp <- log(3750); label("Apparent peripheral volume V3/F (L)") # Thakre 2026 Table 3: V3/F = 3750 L (SEE 1421; 95% CI 1900-7390)
    ltlag <- log(0.0556); label("Absorption lag time ALAG1 (day)") # Thakre 2026 Table 3: ALAG1 = 0.0556 days (SEE 0.00242; 95% CI 0.0510-0.0605)
    lkaratio <- log(2.47); label("Log ratio of absorption rate ka to elimination rate CL/V2 (unitless)") # Thakre 2026 Table 3: Absorption factor = 2.47 (SEE 0.422), footnote b ka estimated via a factor on ke = CL/V2; derived typical ka 10.1 1/day

    # Allometric body-weight exponents (estimated), reference 70 kg
    e_wt_cl <- 0.44; label("Allometric exponent of body weight on CL/F (unitless)") # Thakre 2026 Table 3: Body weight on CL/F = 0.440 (SEE 0.111)
    e_wt_vc <- 0.757; label("Allometric exponent of body weight on V2/F (unitless)") # Thakre 2026 Table 3: Body weight on V2/F = 0.757 (SEE 0.116)

    # Relative bioavailability of the pediatric granules vs the adult tablet
    e_form_granule_fdepot <- fixed(1.2); label("Relative bioavailability F1 of the pediatric granule formulation vs the adult tablet (ratio)") # Thakre 2026 Table 3: Formulation (Part 2) on F1 = 1.20 FIX; Methods: from an adult phase 1 bioavailability study

    # IIV (variances). The paper states CL/F and V2/F carry correlated random
    # effects but prints no covariance, so the block is encoded as diagonal.
    etalcl ~ 0.359 # Thakre 2026 Table 3: IIV on CL/F variance 0.359 (65.7 %CV), shrinkage 4.9%
    etalvc ~ 0.17 # Thakre 2026 Table 3: IIV on V2/F variance 0.170 (43.0 %CV), shrinkage 8.9%

    # Residual error: Table 3 reports sigma^2 = 0.229, so SD = sqrt(0.229).
    propSd <- 0.4785; label("Proportional residual error (fraction)") # Thakre 2026 Table 3: sigma1^2 (proportional) = 0.229 (SEE 0.0164); sqrt(0.229) = 0.4785
  })

  model({
    # Individual parameters; body weight power-scaled at 70 kg on CL/F and V2/F
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    # Absorption rate constant is a fixed multiple of the subject's own
    # elimination rate constant ke = CL/V2 (Thakre 2026 Table 3 footnote b),
    # which keeps ka > ke and excludes flip-flop kinetics.
    karatio <- exp(lkaratio)
    ka <- karatio * kel

    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1

    alag(depot) <- exp(ltlag)
    # F1 = 1 for the adult tablet; the pediatric granules carry the fixed ratio
    f(depot) <- 1 + (e_form_granule_fdepot - 1) * FORM_GRANULE

    # mg / L = ug/mL; x 1000 for ng/mL
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
