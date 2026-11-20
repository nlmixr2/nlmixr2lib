Dickinson_2014_nevirapine <- function() {
  description <- "One-compartment population PK model for oral nevirapine in HIV-infected Malawian adults (Dickinson 2014), with first-order absorption, combined CYP2B6 983T>C (rs28399499) / 516G>T (rs3745274) genotype effects entering additively on linear-scale CL/F, and inter-individual plus inter-occasion variability on CL/F."
  reference <- paste(
    "Dickinson L, Chaponda M, Carr DF, van Oosterhout JJ, Kumwenda J,",
    "Lalloo DG, Pirmohamed M, Heyderman RS, Khoo SH. Population",
    "pharmacokinetic and pharmacogenetic analysis of nevirapine in",
    "hypersensitive and tolerant HIV-infected patients from Malawi.",
    "Antimicrob Agents Chemother. 2014;58(2):706-712.",
    "doi:10.1128/AAC.02069-13.",
    sep = " "
  )
  vignette <- "Dickinson_2014_nevirapine"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "nevirapine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "nevirapine", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    SNP_CYP2B6_RS28399499_C_COUNT = list(
      description = "Count of CYP2B6 c.983T>C (rs28399499) C-alleles per subject (0/1/2). 0 = TT homozygous wild-type, 1 = TC heterozygous, 2 = CC homozygous variant.",
      units = "(count, 0/1/2)",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed (germline genotype). Combined with the 516G>T count into the three genotype groups of the Dickinson 2014 final (reduced) model: reference = 983TT with 516GG or GT; 983TT/516TT; 983TC with 516GG or GT. The 983TC/516GG-or-GT shift (-1.08 L/h) applies when this count is 1. No 983CC homozygotes were observed (Dickinson 2014 Results 'Pharmacogenetic analysis', TC 17.2%), so count 2 is outside the model's support and receives no 983 shift. The heterozygous 983TC with 516TT combination is not described in the source; the model adds both shifts for it (additive extrapolation).",
      source_name = "GEN indicators for CYP2B6 c.983T>C (SNP3)"
    ),
    SNP_CYP2B6_RS3745274_T_COUNT = list(
      description = "Count of CYP2B6 c.516G>T (rs3745274) T-alleles per subject (0/1/2). 0 = GG homozygous wild-type, 1 = GT heterozygous, 2 = TT homozygous variant.",
      units = "(count, 0/1/2)",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed (germline genotype). In the Dickinson 2014 reduced final model 516GT was pooled with 516GG (5.3% difference, Results 'Multiple SNP analysis'), so only the homozygous 516TT group carries a shift (-0.696 L/h), applied when this count is 2. Cohort frequencies GT 48.7%, TT 16.2%.",
      source_name = "GEN indicators for CYP2B6 c.516G>T (SNP4)"
    ),
    OCC = list(
      description = "Integer-valued sampling-occasion index for inter-occasion variability on CL/F",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Three occasion slots, the largest number any subject contributed in the source: rich-sampling patients had one intensive sampling day plus 1 to 2 additional random samples (Dickinson 2014 Results 'Pharmacokinetic sampling'), and 219 occasions were counted across the 180 patients (Results 'Multiple SNP analysis'). All slots share the one estimated IOV magnitude. Set OCC = 1 throughout for a single-occasion simulation, or OCC = 0 to switch inter-occasion variability off.",
      source_name = "occasion (kappa_ij subscript j)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at start of ART",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL/F and not significant (dOFV -0.02 to -3.1; Dickinson 2014 Results 'Pharmacokinetic and covariate model'). Hypersensitivity reaction and viral hepatitis were screened the same way with the same result."
    ),
    WT = list(
      description = "Body weight at start of ART",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL/F and not significant (dOFV -0.02 to -3.1; Dickinson 2014 Results 'Pharmacokinetic and covariate model'). Hypersensitivity reaction and viral hepatitis were screened the same way with the same result."
    ),
    BMI = list(
      description = "Body mass index at start of ART",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened linearly and as a categorical cut-off at 18.5 kg/m^2; not significant (Dickinson 2014 Results 'Pharmacokinetic and covariate model')."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened on CL/F and not significant (dOFV -0.02 to -3.1; Dickinson 2014 Results 'Pharmacokinetic and covariate model'). Hypersensitivity reaction and viral hepatitis were screened the same way with the same result."
    ),
    CD4_ABS = list(
      description = "CD4 cell count at start of ART",
      units = "cells/mm^3",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened linearly and as a categorical cut-off at 250 cells/mm^3; not significant (Dickinson 2014 Results 'Pharmacokinetic and covariate model')."
    ),
    SNP_CYP2B6_RS2279343_G_COUNT = list(
      description = "Count of CYP2B6 c.785A>G (rs2279343) G-alleles (0/1/2)",
      units = "(count, 0/1/2)",
      type = "continuous",
      reference_category = NULL,
      notes = "Significant univariately (785GG -22.9% CL/F) and in a pairwise model with 983T>C, but dropped at backwards elimination (Dickinson 2014 Results 'Multiple SNP analysis')."
    ),
    SNP_CYP3A5_RS10264272_A_COUNT = list(
      description = "Count of CYP3A5*6 (rs10264272) variant alleles (0/1/2)",
      units = "(count, 0/1/2)",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested univariately; no significant association with CL/F (Dickinson 2014 Results 'Pharmacogenetic analysis')."
    ),
    SNP_CYP3A5_RS776746_G_COUNT = list(
      description = "Count of CYP3A5*3 (rs776746) variant alleles (0/1/2)",
      units = "(count, 0/1/2)",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested univariately; no significant association with CL/F (Dickinson 2014 Results 'Pharmacogenetic analysis')."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 180,
    n_studies = 1,
    n_observations = 383,
    age_range = "21-62 years",
    age_median = "34 years",
    weight_range = "35-94 kg",
    weight_median = "54 kg",
    bmi_median = "20 kg/m^2 (range 15-38)",
    cd4_median = "156 cells/mm^3 (range 1-812)",
    sex_female_pct = 56.1,
    race_ethnicity = "Black African (Malawian)",
    disease_state = "HIV-infected, antiretroviral-naive adults starting nevirapine-based ART (WHO stage III/IV and/or CD4 < 250 cells/mm^3); 25 developed nevirapine hypersensitivity (18 rash, 3 jaundice, 4 Stevens-Johnson syndrome); 23 coinfected with hepatitis B and/or C.",
    dose_range = "200 mg nevirapine orally twice daily (fixed-dose combination with lamivudine and stavudine); 3 patients were in the 2-week 200 mg once-daily lead-in.",
    regions = "Malawi (Queen Elizabeth Central Hospital, Blantyre)",
    cyp2b6_freq = "983T>C: TC 17.2%, CC 0%. 516G>T: GT 48.7%, TT 16.2%. Final-model groups (Figure 2): 983TT/516GG or GT n = 102, 983TT/516TT n = 25, 983TC/516GG or GT n = 27, genotype missing n = 26.",
    notes = "Pharmacokinetic substudy of a 1,117-patient prospective cohort, selected to represent sex, hypersensitivity and viral-hepatitis phenotypes. 40 tolerant patients had rich sampling (0, 1, 2, 5, 9 h post-dose, fed, directly observed); 140 had random sampling. Samples 0.42-28 h post-dose; median duration of therapy at sampling 6 weeks (1-26)."
  )

  ini({
    lcl <- log(3.02); label("Apparent oral clearance CL/F, reference genotype (L/h)") # Table 1 final model CL/F = 3.02 L/h (RSE 4.44%); reference = 983TT with 516GG or GT
    lvc <- log(114); label("Apparent volume of distribution V/F (L)") # Table 1 final model V/F = 114 L (RSE 22.8%)
    lka <- log(0.578); label("First-order absorption rate constant ka (1/h)") # Table 1 final model ka = 0.578 1/h (RSE 26.1%)

    e_983tt_516tt_cl <- -0.696; label("Additive shift on CL/F for CYP2B6 983TT/516TT (L/h)") # Table 1 final model theta1 = -0.696 (RSE 36.6%); -23.0% vs reference (Results 'Multiple SNP analysis')
    e_983tc_cl <- -1.08; label("Additive shift on CL/F for CYP2B6 983TC with 516GG or GT (L/h)") # Table 1 final model theta2 = -1.08 (RSE 22.3%); -35.8% vs reference (Results 'Multiple SNP analysis')

    etalcl ~ 0.104657 # Table 1 final model IIV CL/F 33.2% CV (RSE 27.4%); log(1 + 0.332^2)

    etaiov_lcl_1 ~ 0.096312 # Table 1 final model IOV CL/F 31.8% CV (RSE 24.2%); log(1 + 0.318^2); occasion 1
    etaiov_lcl_2 ~ fixed(0.096312) # Table 1 final model IOV CL/F, same magnitude for occasion 2
    etaiov_lcl_3 ~ fixed(0.096312) # Table 1 final model IOV CL/F, same magnitude for occasion 3

    propSd <- 0.131; label("Proportional residual error (fraction)") # Table 1 final model proportional residual error 13.1% (RSE 17.0%)
  })

  model({
    # Genotype groups of the reduced final model: 516GT is pooled with the
    # 516GG reference, so only 516TT carries a shift; 983TC carries its own.
    hom_516 <- (SNP_CYP2B6_RS3745274_T_COUNT == 2)
    het_983 <- (SNP_CYP2B6_RS28399499_C_COUNT == 1)

    # Inter-occasion variability, multiplexed by the occasion index; OCC = 0
    # zeroes every indicator.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    iov_cl <- oc1 * etaiov_lcl_1 + oc2 * etaiov_lcl_2 + oc3 * etaiov_lcl_3

    # CL/F_ij = (CL0 + theta1*GEN1 + theta2*GEN2) * exp(eta_i + kappa_ij)
    tvcl <- exp(lcl) + e_983tt_516tt_cl * hom_516 + e_983tc_cl * het_983
    cl <- tvcl * exp(etalcl + iov_cl)
    vc <- exp(lvc)
    ka <- exp(lka)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
