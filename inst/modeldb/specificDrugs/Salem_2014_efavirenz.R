Salem_2014_efavirenz <- function() {
  description <- "One-compartment population PK model for oral efavirenz in HIV-1-infected children (Salem 2014). Allometric body-weight scaling on apparent clearance (fixed exponent 0.75) and apparent volume of distribution (fixed exponent 1.0) referenced to 70 kg; sigmoid Emax maturation of CL/F with postnatal age (TM50 = 4.6 months, Hill = 3.4); 51% reduction in CL/F for CYP2B6-516 T/T homozygotes; Emax maturation of relative bioavailability for the oral liquid (suspension or solution) formulations vs the capsule reference (mature F = 0.79; TM50 = 10.6 months; Hill fixed at 1)."
  reference <- paste(
    "Salem AH, Fletcher CV, Brundage RC.",
    "Pharmacometric characterization of efavirenz developmental pharmacokinetics",
    "and pharmacogenetics in HIV-infected children.",
    "Antimicrob Agents Chemother. 2014;58(1):136-143.",
    "doi:10.1128/AAC.01738-13."
  )
  vignette <- "Salem_2014_efavirenz"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying (children were followed up to 4 years, with PK assessments at weeks 2, 6, 56, and 112; Salem 2014 Methods 'Patient population and study design'). Drives fixed-exponent allometric scaling on CL/F (exponent 0.75) and V/F (exponent 1.0) with standard reference weight 70 kg per Salem 2014 Methods 'Development of the covariate model' paragraph 3.",
      source_name        = "WT"
    ),
    PNA = list(
      description        = "Postnatal age",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying (per Salem 2014 Methods, AGE is the postnatal age of the subject in months). Drives the sigmoid Emax maturation of CL/F (TM50 = 4.6 months, Hill = 3.4) and the Emax maturation of the relative bioavailability of the oral liquid formulations vs the capsule reference (TM50 = 10.6 months, Hill fixed at 1).",
      source_name        = "AGE"
    ),
    SNP_CYP2B6_RS3745274_T_COUNT = list(
      description        = "Count of CYP2B6 c.516G>T (rs3745274, p.Q172H) T-alleles per subject (0/1/2). 0 = GG homozygous wild-type, 1 = GT heterozygous, 2 = TT homozygous variant.",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (germline genotype). Only the homozygous variant (T/T, count = 2) carries a CL/F effect in Salem 2014 -- no difference was observed between G/T and G/G (Salem 2014 Results paragraph 3). The packaged model derives the binary T/T indicator inline as (SNP_CYP2B6_RS3745274_T_COUNT >= 2). Cohort allele-genotype frequencies (n = 96; Salem 2014 Table 1): GG 36%, GT 30%, TT 13%, Missing 23%. Mixture modelling on the 22 missing-genotype subjects did not suggest a hidden T/T fraction (Salem 2014 Results paragraph 3); for downstream simulation those subjects can be assigned the non-T/T phenotype (count = 0 or 1).",
      source_name        = "CYP2B6-G516T"
    ),
    FORM_CAPSULE = list(
      description        = "Formulation indicator at the dose record: 1 = capsule (structural F = 1 reference), 0 = oral liquid (suspension or solution; F decreases below 1 with decreasing age via an Emax function).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (capsule; F fixed structurally to 1 per Salem 2014 Methods 'Development of the covariate model' paragraph 5)",
      notes              = "Per-dose-record covariate (children could switch formulations during follow-up as they aged out of the liquid formulation). Salem 2014 explicitly states that no difference in relative bioavailability was observed between the suspension and the solution (Results paragraph 4), so the same liquid F profile applies to both. The packaged model gates the age-dependent Emax bioavailability and its associated IIV (etaltvf_liq) so they apply only when FORM_CAPSULE = 0; FORM_CAPSULE = 1 yields F = 1 exactly (no IIV on the capsule arm).",
      source_name        = "FORM"
    )
  )

  population <- list(
    species            = "human",
    n_subjects         = 96,
    n_studies          = 1,
    age_range          = "2-202 months postnatal (~ 0.17-16.8 years)",
    age_median         = "66 months postnatal (~ 5.5 years)",
    weight_range       = "4.8-96.4 kg",
    weight_median      = "18.7 kg",
    bsa_range          = "0.27-2.07 m^2 (median 0.75; Mosteller formula)",
    sex_female_pct     = 59,
    race_ethnicity     = c(`Non-Hispanic white` = 13, `Non-Hispanic black` = 56, Hispanic = 29, `Native American` = 1, Other = 1),
    cyp2b6_freq        = "CYP2B6-516G>T (rs3745274) cohort frequencies (n = 96; Salem 2014 Table 1): G/G 36%, G/T 30%, T/T 13%, Missing 23%. MDR1-C3435T frequencies were collected (C/C 33%, C/T 35%, T/T 7%, Missing 26%) but the MDR1 polymorphism had no effect on CL/F or V/F and is not encoded in the model.",
    disease_state      = "HIV-1-infected children enrolled in the Pediatric AIDS Clinical Trials Group 382 (PACTG382) study; combination antiretroviral therapy of efavirenz plus nelfinavir plus at least one nucleoside reverse transcriptase inhibitor. Inclusion criteria: < 16 years old at enrolment with plasma HIV-1 RNA > 400 copies/mL by reverse transcription-PCR (Amplicor Monitor assay).",
    dose_range         = "AUC-controlled design targeting 24 h steady-state AUC0-24 between 190 and 380 uM*h. Initial allometric dosing: capsule mg/day = (weight[kg]/70)^0.7 * 600 mg; oral suspension or solution mg/day = (weight[kg]/70)^0.7 * 720 mg (the 20% dose top-up anticipates the lower liquid bioavailability). Doses were rounded to the nearest 25 mg, and adjusted proportionately up to a 200 mg maximum increase if the measured AUC fell outside target.",
    regions            = "United States and Puerto Rico (18 participating PACTG sites)",
    n_observations     = "3172 plasma efavirenz concentrations across 96 subjects, sampled before-dose and 2, 5, 8, 12 h post-dose at weeks 2, 6, and (if dose was adjusted at week 6) week 10. Additional 6 h and 24 h samples were collected from cohort I (ages 3-16 years) only. Follow-up AUC samples were collected at weeks 56 and 112 to monitor for changes during growth and development. Concentrations were quantified by validated HPLC at a PACTG Pharmacology Laboratory, LLOQ 0.020 mg/L, assay variability 1-4.5%.",
    study              = "PACTG382 (open-label phase I/II two-cohort study). Cohort I: ages 3-16 years at enrolment. Cohort II: ages 2 months to 8 years at enrolment.",
    notes              = "Baseline demographics from Salem 2014 Table 1. The study was approved by the institutional review boards of all 18 participating sites; informed written consent was obtained from parents or guardians."
  )

  ini({
    # ---- Structural PK parameters (Salem 2014 Table 2 final model) ----
    lka <- log(0.84); label("First-order absorption rate constant ka (1/h)")                          # Salem 2014 Table 2 final ka = 0.84 h^-1 (RSE 12.6%; 90% CI 0.68-1.05)
    lcl <- log(11.2); label("Apparent mature clearance CL/F (L/h) at WT = 70 kg, non-T/T genotype")   # Salem 2014 Table 2 final TVCL = 11.2 L/h (RSE 6.8%; 90% CI 9.9-12.5); reference 70 kg and non-T/T CYP2B6-516 genotype (the mature CL value standardised to a 70 kg adult)
    lvc <- log(468);  label("Apparent volume of distribution V/F (L) at WT = 70 kg")                   # Salem 2014 Table 2 final TVV = 468.0 L (RSE 8.7%; 90% CI 400.8-535.2); reference 70 kg

    # ---- Maturation parameters for CL/F (sigmoid Emax in postnatal age, months) ----
    # Maturation form: cl_mat = PNA^hill_cl / (PNA^hill_cl + tm50_cl^hill_cl)
    # Salem 2014 Results paragraph 2 + equation following Results paragraph 3:
    # "CL/F = 11.2 * (WT/70)^0.75 * [AGE^3.4 / (AGE^3.4 + 4.6^3.4)] * 0.49^GTflag"
    tm50_cl  <- 4.6;        label("Maturation half-life for CL/F (months); age at 50% mature CL/F")   # Salem 2014 Table 2 final TM50,CL = 4.6 months (RSE 8.6%; 90% CI 3.9-5.3)
    lhill_cl <- log(3.4);   label("Hill coefficient for the sigmoid Emax maturation of CL/F (unitless)") # Salem 2014 Table 2 final H = 3.4 (RSE 8.1%; 90% CI 2.9-3.9)

    # ---- CYP2B6-516 T/T genotype effect on CL/F ----
    # CL/F is multiplied by 0.49 when the subject is CYP2B6-516 T/T (51% reduction).
    # Encoded as a log-shift on lcl that is gated by the binary T/T indicator
    # derived inline from SNP_CYP2B6_RS3745274_T_COUNT >= 2.
    e_cyp2b6_tt_cl <- log(0.49); label("Log-ratio of T/T CL/F vs non-T/T (unitless); 51% reduction") # Salem 2014 Table 2 final CYP2B6 T/T GT = 0.49 (RSE 11.5%; 90% CI 0.40-0.58); log(0.49) = -0.7133

    # ---- Fixed allometric exponents on body weight (Salem 2014 Methods 'Development of the covariate model' paragraph 3) ----
    # "The exponents in the allometric model were fixed to 0.75 and 1 for CL/F and V/F, respectively."
    e_wt_cl <- fixed(0.75); label("Allometric exponent of (WT/70) on CL/F (unitless; fixed)")        # Salem 2014 Methods 'Development of the covariate model' paragraph 3
    e_wt_vc <- fixed(1.0);  label("Allometric exponent of (WT/70) on V/F (unitless; fixed)")         # Salem 2014 Methods 'Development of the covariate model' paragraph 3

    # ---- Bioavailability: oral liquid (suspension or solution) vs capsule (reference) ----
    # Salem 2014 Results paragraph 4 + equation following: "F_solution_and_suspension = 0.79 * [AGE/(AGE+10.6)]"
    # The capsule is the structural reference (F = 1). The liquid F is an Emax function of postnatal age
    # with mature asymptote tvf_liq = 0.79 and TM50,F = 10.6 months. The Hill coefficient for the F maturation
    # was fixed at 1 because "the Hill factor was not significantly different from 1, and hence the maturation
    # model was reduced to an E max model and the Hill factor was fixed at 1" (Salem 2014 Results paragraph 4).
    ltvf_liq <- log(0.79); label("Log mature relative bioavailability for liquid (suspension/solution) vs capsule (fraction)") # Salem 2014 Table 2 final TVF = 0.79 (RSE 12.5%; 90% CI 0.63-0.95); log(0.79) = -0.2357
    tm50_f   <- 10.6;      label("Maturation half-life for liquid bioavailability (months); age at 50% mature F") # Salem 2014 Table 2 final TM50,F = 10.6 months (RSE 38.7%; 90% CI 3.8-17.4)

    # ---- IIV (diagonal omega; exponential errors, log-normal per Methods paragraph 2 of Development) ----
    # Salem 2014 Methods 'Development of the population pharmacokinetic base model' paragraph 2:
    # CL/F = TVCL * EXP(ETA); individual PK parameters log-normally distributed.
    # CV-to-variance conversion: omega^2 = log(CV^2 + 1).
    etalcl      ~ 0.189669 # 45.7% CV CL/F (Salem 2014 Table 2 final, RSE 28.4%, 90% CI 35.0-56.4); log(1 + 0.457^2) = 0.189669
    etalvc      ~ 0.174034 # 43.6% CV V/F (Salem 2014 Table 2 final, RSE 30.5%, 90% CI 32.5-54.7); log(1 + 0.436^2) = 0.174034
    etaltvf_liq ~ 0.147731 # 39.9% CV liquid TVF (Salem 2014 Table 2 final, RSE 32.8%, 90% CI 29.1-50.7); log(1 + 0.399^2) = 0.147731
    # Note: Salem 2014 also reports inter-occasion variability on CL/F of 30.0% CV (Salem 2014 Table 2;
    # IOV CL/F RSE 13.7%, 90% CI 26.6-33.4) in addition to the diagonal IIV. The packaged model does NOT
    # encode IOV structurally -- the source paper does not define an operational occasion column for the
    # model-library use case, and the nlmixr2lib convention (Andrews 2017 / Brooks 2021 tacrolimus
    # precedent) is to omit IOV when no occasion mapping is defined. Downstream users who want to simulate
    # IOV can add an OCC indicator and a per-occasion eta on lcl in rxode2. See vignette Assumptions and
    # deviations.

    # ---- Residual error (combined proportional + additive) ----
    # Salem 2014 Table 2 footer: "The residual unexplained variability in efavirenz observed concentrations
    # was described by a proportional error of 25% and an additive standard deviation of 0.25 ug/mL."
    propSd <- 0.25;  label("Proportional residual error (fraction)")                                  # Salem 2014 Table 2 footer: proportional residual error 25%
    addSd  <- 0.25;  label("Additive residual error (mg/L; equivalent to ug/mL)")                     # Salem 2014 Table 2 footer: additive residual SD 0.25 ug/mL = 0.25 mg/L (numerically equal in concentration units)
  })

  model({
    # 1. Derive the binary CYP2B6-516 T/T indicator from the SNP allele count.
    #    Salem 2014 Results paragraph 3: "Children with the CYP2B6-516-T/T genotype were found
    #    to have a CL/F 51% lower than that of the other children. No difference in CL/F was
    #    observed between CYP2B6-516-G/T and CYP2B6-516-G/G genotypes." The T/T indicator is
    #    therefore the homozygous-variant level only (T-allele count = 2).
    is_tt <- (SNP_CYP2B6_RS3745274_T_COUNT >= 2)

    # 2. Maturation of CL/F: sigmoid Emax in postnatal age (months).
    #    Maturation = PNA^hill_cl / (PNA^hill_cl + tm50_cl^hill_cl)
    hill_cl       <- exp(lhill_cl)
    maturation_cl <- PNA^hill_cl / (PNA^hill_cl + tm50_cl^hill_cl)

    # 3. Individual PK parameters with fixed-exponent allometric weight scaling
    #    (reference WT = 70 kg), maturation, and CYP2B6 T/T genotype effect.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl + e_cyp2b6_tt_cl * is_tt) * (WT / 70)^e_wt_cl * maturation_cl
    vc <- exp(lvc + etalvc)                          * (WT / 70)^e_wt_vc

    # 4. Micro-constants
    kel <- cl / vc

    # 5. ODE system: one-compartment model with first-order absorption (oral)
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 6. Bioavailability: capsule is the structural reference (F = 1); liquid F is an
    #    Emax function of postnatal age. The IIV etaltvf_liq applies on the log-scale
    #    of tvf_liq, and the FORM_CAPSULE indicator in fdepot then selects the F = 1
    #    capsule reference (eta value computed but unused) or the liquid Emax F.
    tvf_liq  <- exp(ltvf_liq + etaltvf_liq)
    f_liquid <- tvf_liq * PNA / (PNA + tm50_f)
    fdepot   <- FORM_CAPSULE * 1 + (1 - FORM_CAPSULE) * f_liquid
    f(depot) <- fdepot

    # 7. Observation and error
    #    Concentration: dose in mg, volume in L -> mg/L (numerically equal to ug/mL).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
