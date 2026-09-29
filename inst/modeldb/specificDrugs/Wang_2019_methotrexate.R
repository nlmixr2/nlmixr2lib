Wang_2019_methotrexate <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and",
    "linear elimination for low-dose oral methotrexate in Chinese adults with",
    "rheumatoid arthritis (Wang 2019; n = 71 patients, 85 mostly trough",
    "concentrations). Clearance depends on the SLCO1B1 (OATP1B1) c.388A>G",
    "(rs2306283) genotype through multiplicative factors for the 388A/G",
    "heterozygote and 388G/G homozygote relative to 388A/A. Volume,",
    "absorption rate and bioavailability are fixed. Exponential",
    "between-subject variability on clearance only, and combined",
    "proportional plus additive residual error.",
    sep = " "
  )
  reference <- paste(
    "Wang Z, Zhang N, Chen C, Chen S, Xu J, Zhou Y, Zhao X, Cui Y (2019).",
    "Influence of the OATP Polymorphism on the Population Pharmacokinetics of",
    "Methotrexate in Chinese Patients.",
    "Curr Drug Metab 20(7):592-600. doi:10.2174/1389200220666190701094756.",
    sep = " "
  )
  vignette <- "Wang_2019_methotrexate"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "methotrexate", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "methotrexate", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    SNP_SLCO1B1_RS2306283_HET = list(
      description = "SLCO1B1 c.388A>G (rs2306283) heterozygous 388A/G indicator (1 = A/G, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Source column RS2306283 (a.k.a. RS230) coded 1 = 388G/G, 2 = 388A/G,",
        "3 = 388A/A (Abstract). Derive SNP_SLCO1B1_RS2306283_HET =",
        "(RS2306283 == 2). Multiplies clearance by 0.647 (Table 3 theta6);",
        "the Discussion describes this as CL/F 'diminished by 32.3% in",
        "patients carrying the OATP1B1-388AG' genotype. Cohort frequency",
        "22-23 of 71 (Tables 1 and 2)."
      ),
      source_name = "RS2306283"
    ),
    SNP_SLCO1B1_RS2306283_HOM = list(
      description = "SLCO1B1 c.388A>G (rs2306283) homozygous-variant 388G/G indicator (1 = G/G, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Derive SNP_SLCO1B1_RS2306283_HOM = (RS2306283 == 1) from the source",
        "coding above. Multiplies clearance by 0.805 (Table 3 theta5); the",
        "Discussion quotes a 17.8% CL/F reduction for 388GG versus 388AA.",
        "G/G is the MAJORITY genotype in this Chinese cohort (42-43 of 71,",
        "Tables 1 and 2), so the reference 388A/A stratum is only 6 subjects."
      ),
      source_name = "RS2306283"
    )
  )

  covariatesDataExcluded <- list(
    RBC = list(
      description = "Red blood cell count",
      units = "10^12/L",
      type = "continuous",
      notes = "Entered CL/F in the SCM forward step but was removed in backward deletion (Section 3.3)."
    ),
    SNP_SLCO1B1_RS4149056 = list(
      description = "SLCO1B1 c.521T>C (rs4149056) genotype",
      units = "(categorical)",
      type = "categorical",
      notes = "Entered CL/F in the SCM forward step but was removed in backward deletion (Section 3.3)."
    ),
    SNP_SLC19A1_RS1051266 = list(
      description = "SLC19A1 rs1051266 G>A genotype",
      units = "(categorical)",
      type = "categorical",
      notes = "Entered CL/F in the SCM forward step but was removed in backward deletion (Section 3.3)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 71L,
    n_studies = 1L,
    n_observations = 85L,
    age_range = "mean 48.0 years (SD 15.2)",
    weight_range = "mean 59.4 kg (SD 10.7)",
    height = "mean 1.61 m (SD 0.06)",
    bsa = "mean 1.6 m^2 (SD 0.16)",
    sex_female_pct = 100 * 11 / 71,
    race_ethnicity = "Chinese (100%)",
    disease_state = "Rheumatoid arthritis on low-dose methotrexate therapy",
    dose_range = "Low-dose oral methotrexate; the dose and dosing schedule are not reported.",
    regions = "China (Peking University First Hospital, Beijing).",
    genotypes = "SLCO1B1 rs2306283: GG 43, AG 22, AA 6 (Table 2).",
    co_medication = "Leflunomide 17/71, folic acid 46/71, hydroxychloroquine 23/71, calcium 23/71, prednisolone 11/71 (Table 1).",
    notes = paste(
      "Table 1 prints 'GEND(Male/Female) 60/11'; sex_female_pct follows that",
      "print literally although a predominantly male RA cohort is atypical",
      "and the columns may be transposed. Most samples were trough",
      "concentrations (Discussion), which is why V/F, Ka and F were fixed."
    )
  )

  ini({
    lcl <- log(7.75); label("Clearance for SLCO1B1 388A/A (L/h)") # Table 3 theta1 'CL/F (L/h)' = 7.75
    lvc <- fixed(log(32.8)); label("Volume of distribution (L)") # Table 3 theta2 'V/F (L)' = 32.8 FIX
    lka <- fixed(log(1.69)); label("Absorption rate constant (1/h)") # Table 3 theta3 'Ka' = 1.69 FIX (unit printed 'h/L', a typo for 1/h)
    lfdepot <- fixed(log(0.704)); label("Oral bioavailability (fraction)") # Table 3 theta4 'F' = 0.704 FIX; final-model equation 'F = 0.704'

    e_snp_slco1b1_rs2306283_hom_cl <- 0.805; label("Clearance multiplier for SLCO1B1 388G/G (unitless)") # Table 3 theta5 'RS2306283 on CL' = 0.805; equation 'RS230 = 1'
    e_snp_slco1b1_rs2306283_het_cl <- 0.647; label("Clearance multiplier for SLCO1B1 388A/G (unitless)") # Table 3 theta6 'RS2306283 on CL' = 0.647; equation 'RS230 = 2'

    etalcl ~ fixed(0.167) # Table 3 'omega CL' = 0.167 FIX, NONMEM OMEGA variance; omega V, Ka, F are 0 FIX and omitted

    propSd <- sqrt(0.713); label("Proportional residual error (fraction)") # Table 3 'sigma1 (pro)' = 0.713, a SIGMA variance (bootstrap-CI width test)
    addSd <- sqrt(2.83); label("Additive residual error (ng/mL)") # Table 3 'sigma2 (add)' = 2.83, a SIGMA variance (bootstrap-CI width test)
  })

  model({
    # Final-model equation (Section 3.3) printed as
    # 'RS230 = 1 : CL = 7.75 x e^(0.167 x 0.805)', '= 2 : ... 0.647', '= 3 : ... 1'.
    # Encoded as CL = theta1 x theta_genotype x exp(eta), the reading
    # reproduced by the base-model CL/F (5.98 L/h) and the Discussion's
    # quoted genotype reductions; see the vignette Errata.
    cl <- exp(lcl + etalcl) *
      e_snp_slco1b1_rs2306283_hom_cl^SNP_SLCO1B1_RS2306283_HOM *
      e_snp_slco1b1_rs2306283_het_cl^SNP_SLCO1B1_RS2306283_HET
    vc <- exp(lvc)
    ka <- exp(lka)
    fdepot <- exp(lfdepot)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    f(depot) <- fdepot

    # Amounts in mg, volume in L -> mg/L; x 1000 gives ng/mL (assay units, LLOQ 0.5 ng/mL).
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
