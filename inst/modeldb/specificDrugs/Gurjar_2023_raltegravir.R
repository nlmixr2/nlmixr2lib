Gurjar_2023_raltegravir <- function() {
  description <- "Two-compartment first-order-absorption population PK model for oral raltegravir 400 mg twice daily in treatment-naive HIV-1-infected adults of the NEAT001/ANRS143 trial (darunavir/ritonavir background), estimated from sparse single samples at weeks 4 and 24 with NONMEM $PRIOR (NWPRI) informative priors on Q/F, Vp/F, ka and Vc/F from Arab-Alameddine 2012; interindividual variability on CL/F only, proportional residual error, and no retained covariates (weight, age, sex, ethnicity, UGT1A1*28 and SLC22A6 genotypes all screened and rejected) (Gurjar 2023)."
  reference <- "Gurjar R, Dickinson L, Carr D, Stohr W, Bonora S, Owen A, D'Avolio A, Cursley A, De Castro N, Fatkenheuer G, Vandekerckhove L, Di Perri G, Pozniak A, Schwimmer C, Raffi F, Boffito M, and the NEAT001/ANRS143 Study Group. Influence of UGT1A1 and SLC22A6 polymorphisms on the population pharmacokinetics and pharmacodynamics of raltegravir in HIV-infected adults: a NEAT001/ANRS143 substudy. Pharmacogenomics J. 2023;23:14-20. doi:10.1038/s41397-022-00293-5"
  vignette <- "Gurjar_2023_raltegravir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened univariably as allometric scaling on CL/F, Q/F (exponent 0.75) and Vc/F, Vp/F (exponent 1) with reference 70 kg; dOFV -5.3 for 4 d.f., not significant (Supplementary Table S1). Not retained. Cohort median 72 (41-135) kg (Table 1).",
      source_name = "WEIGHT0"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened univariably as a linear effect on CL/F centred on the median 36 years; dOFV -1.4, not significant (Supplementary Table S1). Not retained.",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened univariably as a power-of-indicator effect on CL/F (reference male); dOFV -2.1, not significant (Supplementary Table S1). Not retained.",
      source_name = "SEX"
    ),
    RACE_BLACK = list(
      description = "Black race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Caucasian reference)",
      notes = "Screened with RACE_ASIAN and RACE_OTHER as indicator effects on CL/F (reference Caucasian); dOFV -1.0 for 3 d.f., not significant (Supplementary Table S1). Not retained.",
      source_name = "ETH1/ETH2/ETH3"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Caucasian, Black or Other)",
      notes = "Screened on CL/F jointly with RACE_BLACK / RACE_OTHER, and separately against a pooled Caucasian/Black/Other reference (dOFV -0.9); not significant (Supplementary Table S1). Not retained.",
      source_name = "ETH1/ETH2/ETH3"
    ),
    RACE_OTHER = list(
      description = "Other race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Caucasian reference)",
      notes = "Screened on CL/F jointly with RACE_BLACK / RACE_ASIAN; not significant (Supplementary Table S1). Not retained.",
      source_name = "ETH1/ETH2/ETH3"
    ),
    UGT1A1_STAR28_HOM = list(
      description = "UGT1A1 low-activity genotype indicator (in this cohort, UGT1A1*28/*28 only)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal or reduced UGT1A1 activity: *1/*1, *1/*36, *1/*28, *28/*36, *36/*37)",
      notes = "The paper groups UGT1A1 (rs8175347) diplotypes by CPIC activity as normal / reduced / low (low = *28/*28, *28/*37, *37/*37; only *28/*28 occurred, n = 40). Screened on CL/F with normal as reference (REDUCED, LOW and MISSING indicators; dOFV -4.8 for 3 d.f.) and with normal+reduced pooled as reference (LOW and MISSING; dOFV -3.6 for 2 d.f.); neither significant (Supplementary Table S1). Results text reports a non-significant 21% lower CL/F in the low-activity group. Not retained.",
      source_name = "UGT1A1"
    ),
    SNP_SLC22A6_RS4149170 = list(
      description = "SLC22A6 (OAT1) 453G>A (rs4149170) variant indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (GG)",
      notes = "Screened on CL/F as AG, AA and MISSING indicators (reference GG; dOFV -3.8 for 3 d.f.) and as AA vs pooled GG/AG (dOFV -3.8 for 2 d.f.); not significant (Supplementary Table S1). Not retained. Name follows the SNP_<GENE>_RS<rsid> family; not registered because the model does not use it.",
      source_name = "OAT1AG/OAT1AA"
    ),
    SNP_SLC22A6_RS11568626 = list(
      description = "SLC22A6 (OAT1) 728C>T (rs11568626) variant carrier indicator (CT or TT)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CC)",
      notes = "Screened on CL/F as a pooled CT/TT indicator plus a MISSING indicator (reference CC); dOFV -3.6 for 2 d.f., not significant (Supplementary Table S1). Not retained. Name follows the SNP_<GENE>_RS<rsid> family; not registered because the model does not use it.",
      source_name = "rs11568626"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "raltegravir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "raltegravir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "raltegravir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 349L,
    n_studies = 1L,
    n_observations = 602L,
    age_range = "20-71 years",
    age_median = "37 years",
    weight_range = "41-135 kg",
    weight_median = "72 kg",
    sex_female_pct = 12.3,
    race_ethnicity = c(Caucasian = 82.5, Black = 12.6, Asian = 2.3, Other = 2.6),
    disease_state = "Treatment-naive HIV-1-infected adults (plasma HIV RNA > 1000 copies/mL, CD4 < 500 cells/mm^3 unless symptomatic) randomised to the raltegravir + darunavir/ritonavir NRTI-sparing arm of NEAT001/ANRS143. Median baseline CD4 340 (5-780) cells/mm^3; median HIV RNA 4.82 (3.11-6.31) log10 copies/mL.",
    dose_range = "Oral raltegravir 400 mg twice daily with darunavir/ritonavir 800/100 mg once daily.",
    regions = "Europe (78 sites in 15 countries).",
    co_medication = "Darunavir/ritonavir (all subjects).",
    genotypes = "UGT1A1 (rs8175347): normal 109, reduced 115, low (*28/*28) 40, *36/*36 1, missing 84. SLC22A6 453G>A: GG 216, GA 68, AA 9, missing 56. SLC22A6 728C>T: CC 285, CT 6, TT 2, missing 56.",
    notes = "Table 1 of Gurjar 2023. Single plasma samples at week 4 (n = 313) and week 24 (n = 289), 0.17-16.0 h post-dose; LLQ 0.0117 mg/L. NONMEM 7.3 FOCE-I with $PRIOR NWPRI on Q/F, Vp/F, ka and Vc/F (prior means 8.5 L/h, 113 L, 0.21 1/h, 223 L from Arab-Alameddine 2012); CL/F estimated without a prior."
  )

  ini({
    # Table 2 'Fixed effects'; NONMEM control stream in the Supplementary Information
    # (ADVAN4 TRANS4, THETA(1..5) = Q, V3, KA, V2, CL). Q/F, Vp/F, ka and Vc/F were
    # estimated under informative NWPRI priors; CL/F had no prior.
    lcl <- log(55.8); label("Apparent oral clearance CL/F (L/h)") # Table 2: CL/F = 55.8 L/h (RSE 4.1%)
    lvc <- log(194); label("Apparent central volume of distribution Vc/F (L)") # Table 2: Vc/F = 194 L (RSE 6.5%)
    lq <- log(13.0); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2: Q/F = 13.0 L/h (RSE 4.0%)
    lvp <- log(117); label("Apparent peripheral volume of distribution Vp/F (L)") # Table 2: Vp/F = 117 L (RSE 0.6%)
    lka <- log(1.12); label("First-order absorption rate constant ka (1/h)") # Table 2: ka = 1.12 1/h (RSE 13.0%)

    # Table 2 'Random effects': IIV CL/F 62.7% (RSE 12.1%), read as a CV; log(1 + 0.627^2) = 0.33156.
    # The control stream fixes the IIV on Vc, Q, Vp and ka to 0, so they are omitted here.
    etalcl ~ 0.33156

    # Table 2 'Residual error': proportional 69.9% (RSE 7.0%); Y = F*(1 + ERR(1)).
    propSd <- 0.699; label("Proportional residual error (fraction)") # Table 2: proportional = 69.9%
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg and volumes in L give mg/L, the paper's concentration unit.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
