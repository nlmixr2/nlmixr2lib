Stillemans_2022_atorvastatin <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption for",
    "atorvastatin acid in adult ambulatory patients at cardiovascular risk",
    "(Stillemans 2022), fitted jointly to a sparsely sampled real-life",
    "investigation cohort and two richly sampled support datasets. Apparent",
    "clearance CL/F carries a separate typical value and a separate",
    "log-normal random effect for the investigation cohort and for the",
    "support datasets (STUDY_ATORVA_SUPPORT); in the investigation cohort",
    "CL/F is reduced by the SLCO1B1 c.521T>C (rs4149056) genotype,",
    "CL/F = theta_CL * (1 + theta_SLCO1B1), with separate heterozygous and",
    "homozygous-variant effects. Q/F, Vc/F and Vp/F carry log-normal IIV;",
    "ka is fixed to 2.5 1/h with no IIV. Exponential residual error."
  )
  reference <- paste(
    "Stillemans G, Paquot A, Muccioli GG, Hoste E, Panin N, Asberg A,",
    "Balligand JL, Haufroid V, Elens L. Atorvastatin population",
    "pharmacokinetics in a real-life setting: Influence of genetic",
    "polymorphisms and association with clinical response.",
    "Clin Transl Sci. 2022;15(3):667-679.",
    "doi:10.1111/cts.13185.",
    sep = " "
  )
  vignette <- "Stillemans_2022_atorvastatin"

  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL"
  )

  covariateData <- list(
    STUDY_ATORVA_SUPPORT = list(
      description = paste(
        "Dataset indicator: 1 = subject from one of the two richly sampled",
        "support datasets (Lemahieu 2005, Hermann 2006); 0 = subject from",
        "the sparsely sampled real-life investigation cohort."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (investigation cohort)",
      notes = paste(
        "Stillemans 2022 Results, Covariate analysis: 'clearance was defined",
        "with two fixed and two random effects: one pair for the",
        "investigation cohort and one pair for the rest of the dataset' so",
        "that covariate relationships could be estimated in the",
        "investigation cohort (the only one with covariate data) without",
        "affecting the support-dataset estimates. Selects between",
        "lcl_invest / etalcl_invest (STUDY_ATORVA_SUPPORT = 0) and",
        "lcl_support / etalcl_support (STUDY_ATORVA_SUPPORT = 1). The",
        "SLCO1B1 effect applies only when STUDY_ATORVA_SUPPORT = 0 (genotype",
        "was coded as missing in the support datasets). Set to 0 to",
        "simulate the real-life ambulatory population."
      ),
      source_name = "cohort (investigation vs support)"
    ),
    SNP_SLCO1B1_RS4149056_HET = list(
      description = "SLCO1B1 c.521T>C heterozygous (521T/C) indicator.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (521T/T wild type, with SNP_SLCO1B1_RS4149056_HOM = 0)",
      notes = paste(
        "CL/F in the investigation cohort = theta_CL * (1 + theta_SLCO1B1),",
        "theta_SLCO1B1 = -0.402 for 521TC (Stillemans 2022 Table 3 and its",
        "footnote). 16 of 70 investigation-cohort patients (22.9%) were TC",
        "(Table 2). Only effect retained after backward elimination",
        "(Results, Covariate analysis)."
      ),
      source_name = "SLCO1B1 521TC"
    ),
    SNP_SLCO1B1_RS4149056_HOM = list(
      description = "SLCO1B1 c.521T>C homozygous-variant (521C/C) indicator.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (521T/T wild type, with SNP_SLCO1B1_RS4149056_HET = 0)",
      notes = paste(
        "theta_SLCO1B1 = -0.041 for 521CC (Stillemans 2022 Table 3; RSE",
        "469.7%, 95% CI -0.422 to 0.339). The paper notes the effect was",
        "'more pronounced in heterozygotes (and the effect in CC",
        "homozygotes was estimated with very poor precision)' (Results,",
        "Covariate analysis); 5 of 70 patients (7.1%) were CC (Table 2)."
      ),
      source_name = "SLCO1B1 521CC"
    )
  )

  # Screened in the stepwise covariate search on CL/F but not retained in the
  # final model (Stillemans 2022 Methods, Covariate analysis; Results,
  # Covariate analysis). OATP2B1-inhibitor co-medication (dOFV -4.3) and sex
  # (dOFV -4) entered the forward step but were removed at backward
  # elimination (alpha = 0.01); no final-model coefficient is reported for any
  # of them. Race, dose, the other genotypes (CYP3A4*22, CYP3A5*3, ABCB1
  # c.1199G>A and c.3435C>T, ABCC1 c.2012G>T, SLCO1B1 c.388A>G, SLCO2B1
  # c.935G>A, SLCO1B3 c.334T>G) and the other interaction classes were also
  # screened without being retained.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Entered the forward search on CL/F (dOFV -4; women had lower",
        "CL/F, Figure 1b) but removed at backward elimination; no",
        "final-model estimate reported."
      )
    ),
    CONMED_OATP1B_INH = list(
      description = "Co-medication with an OATP inhibitor (here: OATP2B1 inhibitor).",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "The paper's screened class is OATP2B1 inhibitors (usually",
        "L-thyroxine; 7 of 70 patients, Table 1). Entered the forward",
        "search on CL/F (dOFV -4.3; lower CL/F, Figure 1a) but removed at",
        "backward elimination; no final-model estimate reported."
      )
    ),
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      notes = "Screened on CL/F; not retained."
    ),
    BMI = list(
      description = "Body mass index.",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened on CL/F; not retained."
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "atorvastatin",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "atorvastatin",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "atorvastatin",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 111L,
    n_studies = 3L,
    n_observations_investigation = 132L,
    age_median = "53.8 years (IQR 21.6), investigation cohort",
    bmi_median = "26.0 kg/m^2 (IQR 5.4), investigation cohort",
    sex_female_pct = 50,
    race_ethnicity = c(White = 94.3, Hispanic = 2.9, Other = 2.9),
    disease_state = paste(
      "Investigation cohort: ambulatory adults with hypercholesterolemia",
      "at risk of cardiovascular disease (91.4% primary prevention),",
      "normal hepatic function, de novo (27.1%), recently switched (17.1%)",
      "or long-term (55.7%) atorvastatin therapy. Support datasets: healthy",
      "volunteers (Lemahieu 2005) and patients with statin-induced myopathy",
      "plus healthy controls (Hermann 2006)."
    ),
    dose_range = paste(
      "Oral atorvastatin 5-80 mg q24h in the investigation cohort",
      "(20 mg 31.4%, 40 mg 27.1%, 10 mg 20.0%, 80 mg 14.3%, 5 mg 5.7%,",
      "30 mg 1.4%); 40 mg q24h in Lemahieu 2005 and 10 mg q24h in Hermann",
      "2006 (Supplementary Table S1)."
    ),
    regions = "Belgium (investigation cohort, Cliniques Universitaires Saint-Luc, Brussels); support-dataset sites not stated in the paper",
    notes = paste(
      "Three datasets (Stillemans 2022 Supplementary Table S1):",
      "investigation cohort n = 70 patients, sparse opportunistic sampling",
      "(one sample per visit, 1-4 visits, 132 samples, post-intake delay",
      "2.2-40 h, median 14.6 h; NCT03604471); Lemahieu 2005 n = 13 healthy",
      "volunteers, 6 h profile from the atorvastatin-alone control phase of",
      "a calcineurin-inhibitor interaction study; Hermann 2006 n = 28 (13",
      "patients with statin-induced myopathy and 15 healthy controls), 24 h",
      "profile. Demographics (Tables 1 and 2) are reported for the",
      "investigation cohort only. SLCO1B1 c.521T>C in the investigation",
      "cohort: TT 44 (62.9%), TC 16 (22.9%), CC 5 (7.1%), missing 5 (7.1%)."
    )
  )

  ini({
    # Structural parameters (Stillemans 2022 Table 3, final model). All
    # clearances and volumes are apparent (X/F) values (Table 3 note).
    # Two CL/F typical values -- one per cohort -- estimated in the same
    # NONMEM run (Results, Covariate analysis), so both carry a stratum suffix.
    lcl_invest  <- log(535);  label("Apparent clearance CL/F, investigation cohort, SLCO1B1 521TT (L/h)") # Table 3: theta CL investigation = 535 L/h (RSE 6.2%, 95% CI 470-600)
    lcl_support <- log(400);  label("Apparent clearance CL/F, support datasets (L/h)")                     # Table 3: theta CL support = 400 L/h (RSE 7.4%, 95% CI 342-458)
    lq          <- log(1690); label("Apparent inter-compartmental clearance Q/F (L/h)")                    # Table 3: theta Q = 1690 L/h (RSE 20.3%, 95% CI 1018-2362)
    lvc         <- log(1960); label("Apparent central volume Vc/F (L)")                                    # Table 3: theta Vc = 1960 L (RSE 21.7%, 95% CI 1125-2795)
    lvp         <- log(3900); label("Apparent peripheral volume Vp/F (L)")                                 # Table 3: theta Vp = 3900 L (RSE 15.7%, 95% CI 2700-5100)
    lka         <- fixed(log(2.5)); label("First-order absorption rate constant ka (1/h)")                 # Table 3: theta ka = 2.5 1/h (fixed); Results, Structural model: fixed to a literature value, IIV set to zero

    # SLCO1B1 c.521T>C effect on CL/F in the investigation cohort:
    # CL/F = theta_CL * (1 + theta_SLCO1B1), theta_SLCO1B1 = 0 for 521TT
    # (Table 3 footnote).
    e_snp_slco1b1_rs4149056_het_cl <- -0.402; label("Fractional change in CL/F for SLCO1B1 521T/C vs 521T/T (unitless)") # Table 3: theta SLCO1B1 521TC = -0.402 (RSE 15.2%, 95% CI -0.522 to -0.282)
    e_snp_slco1b1_rs4149056_hom_cl <- -0.041; label("Fractional change in CL/F for SLCO1B1 521C/C vs 521T/T (unitless)") # Table 3: theta SLCO1B1 521CC = -0.041 (RSE 469.7%, 95% CI -0.422 to 0.339)

    # IIV. Table 3 labels every omega row '(SD)', but each printed estimate is
    # the NONMEM OMEGA variance: the midpoint of each row's 95% CI equals the
    # square root of the printed estimate (CL invest sqrt(0.0677) = 0.260 vs
    # (0.176 + 0.344) / 2 = 0.260; CL support 0.442 vs 0.442; Q 0.846 vs 0.845;
    # Vc 1.086 vs 1.090; Vp 0.713 vs 0.713), i.e. the CI is on the SD scale
    # and the estimate is not. The Results text also gives the CL/F variance
    # as 0.127 before and 0.0677 after the SLCO1B1 covariate. Diagonal omega
    # (a full block improved OFV but was not used; Results, Structural model).
    etalcl_invest  ~ 0.0677 # Results, Covariate analysis: IIV of CL/F lowered from 0.127 to 0.0677; Table 3 omega CL investigation = 0.068 (RSE 16.5%, SD-scale CI 0.176-0.344, shrinkage 44.3%)
    etalcl_support ~ 0.195  # Table 3: omega CL support = 0.195 (RSE 9.8%, SD-scale CI 0.357-0.527, shrinkage 40.1%)
    etalq          ~ 0.715  # Table 3: omega Q = 0.715 (RSE 19.3%, SD-scale CI 0.526-1.164, shrinkage 52.4%)
    etalvc         ~ 1.18   # Table 3: omega Vc = 1.18 (RSE 16.1%, SD-scale CI 0.745-1.435, shrinkage 40.3%)
    etalvp         ~ 0.508  # Table 3: omega Vp = 0.508 (RSE 19.5%, SD-scale CI 0.441-0.985, shrinkage 38.8%)

    # Residual error: exponential (Results, Structural model). Table 3 prints
    # sigma = 0.085 labelled '(SD)' with 95% CI 0.067-0.103 centred on 0.085
    # and RSE 5.5%. Read as the SIGMA variance: an SD-scale RSE of 5.5% is a
    # variance-scale RSE of 11%, giving 0.085 +/- 1.96 * 0.0094 = 0.067-0.103,
    # which is the printed CI; reading 0.085 as the SD gives 0.076-0.094.
    # ini() takes the SD, so expSd = sqrt(0.085).
    expSd <- 0.291548; label("Exponential (log-normal) residual error SD (log scale)") # sqrt(0.085); Table 3: sigma = 0.085 (RSE 5.5%, 95% CI 0.067-0.103, shrinkage 13.9%)
  })

  model({
    # Individual parameters. The cohort indicator selects which CL/F typical
    # value and which random effect apply (Results, Covariate analysis).
    snp_cl <- 1 +
      e_snp_slco1b1_rs4149056_het_cl * SNP_SLCO1B1_RS4149056_HET +
      e_snp_slco1b1_rs4149056_hom_cl * SNP_SLCO1B1_RS4149056_HOM
    cl_invest <- exp(lcl_invest + etalcl_invest) * snp_cl
    cl_support <- exp(lcl_support + etalcl_support)
    cl <- cl_invest * (1 - STUDY_ATORVA_SUPPORT) +
      cl_support * STUDY_ATORVA_SUPPORT

    ka <- exp(lka)
    q <- exp(lq + etalq)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg and volumes in L give mg/L; x 1000 gives ng/mL, the unit of
    # the paper's observations (Results, Final model evaluation: 'above 10
    # ng ml-1').
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd)
  })
}
