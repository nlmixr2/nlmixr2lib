Stillemans_2021_darunavir <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order absorption and elimination for oral cobicistat- or ritonavir-boosted darunavir in adults with HIV-1 infection (Stillemans 2021). Apparent clearance decreases exponentially with alpha-1 acid glycoprotein (AAG), is lower in women and in CYP3A5 nonexpressers (*3/*3); apparent volume decreases exponentially with AAG and is higher in SLCO3A1 rs8027174 G>T carriers."
  reference <- paste(
    "Stillemans G, Belkhir L, Vandercam B, Vincent A, Haufroid V, Elens L.",
    "Exploration of Reduced Doses and Short-Cycle Therapy for",
    "Darunavir/Cobicistat in Patients with HIV Using Population",
    "Pharmacokinetic Modeling and Simulations. Clin Pharmacokinet.",
    "2021;60:177-189. doi:10.1007/s40262-020-00920-z. PMCID: PMC7862523.",
    sep = " "
  )
  vignette <- "Stillemans_2021_darunavir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    AAG = list(
      description = "Plasma alpha-1 acid glycoprotein concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters CL/F and V/F through the exponential form",
        "P = theta * exp(theta_cov * (AAG - median)) (Methods section 2.6",
        "equation; ESM Online Resource 2 lists both AAG effects as",
        "'Exponential'), with theta_cov = -0.61 on CL/F and -0.68 on V/F",
        "per g/L (Table 3). The source centres on the cohort median but",
        "does NOT print that median anywhere (Table 1, the ESM and the",
        "companion validation paper doi:10.1007/s00228-020-03036-2 omit",
        "AAG). The model file centres on 1 g/L, a round value inside the",
        "adult reference interval; see the vignette Assumptions section",
        "for the size of the resulting typical-value shift. The only AAG",
        "values the source prints are the two high-AAG outliers, 2.66 and",
        "2.2 g/L (Results 3.4; ESM Online Resource 3). Missing values were",
        "substituted by the most recent measurement or the population",
        "median (Methods 2.6)."
      ),
      source_name = "AAG"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Enters CL/F as exp(-0.21 * SEXF): ESM Online Resource 2 lists the",
        "sex effect on CL with an 'Exponential' equation and Table 3",
        "reports 'Female sex on CL' = -0.21, giving a female/male ratio of",
        "0.81. The Results prose rounds this to 'CL was 21% lower in",
        "female individuals'. 33.1% female (42 of 127; Table 1)."
      ),
      source_name = "sex"
    ),
    CYP3A5_EXPR = list(
      description = "CYP3A5 expresser status: 1 = at least one CYP3A5*1 allele (*1/*1 or *1/*3), 0 = *3/*3 nonexpresser",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (nonexpresser, register convention); the published typical value is for expressers (CYP3A5_EXPR = 1)",
      notes = paste(
        "The source parameter is 'CYP3A5*3 on CL' = -0.16 (Table 3), a",
        "categorical effect P = theta * (1 + theta_cov) (Methods 2.6)",
        "applied to nonexpressers; subjects were stratified 'simply as",
        "expressors vs non-expressors in the case of CYP3A5*3' (Methods",
        "2.6) and the Results state 'CL being 19% higher in expressors",
        "compared with non-expressors' (1 / 0.84 = 1.19). The indicator",
        "is therefore inverted relative to the canonical column: CL/F is",
        "multiplied by (1 + e_cyp3a5_cl * (1 - CYP3A5_EXPR)). Genotypes:",
        "*1/*1 33, *1/*3 33, *3/*3 58, missing 3 (Table 2); missing",
        "genotypes were set to the most frequent category."
      ),
      source_name = "CYP3A5*3 (rs776746)"
    ),
    SNP_SLCO3A1_RS8027174 = list(
      description = "SLCO3A1 rs8027174 (g.91941607G>T) variant indicator: 1 = at least one T allele, 0 = G/G",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (G/G)",
      notes = paste(
        "Enters V/F as (1 + 0.81 * SNP_SLCO3A1_RS8027174) (Table 3",
        "'SLCO3A1 G>T on V' = 0.81, ESM Online Resource 2 'Categorical';",
        "Results 'V was increased by 81% in individuals harboring SLCO3A1",
        "rs8027174'). Genotypes G/G 105, G/T 18, T/T 0, missing 4",
        "(Table 2), so the indicator is heterozygote-vs-G/G in this cohort."
      ),
      source_name = "SLCO3A1 G>T (rs8027174)"
    )
  )

  covariatesDataExcluded <- list(
    RACE_BLACK = list(
      description = "African race indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Retained on V/F only in the rich-sampling (dataset B) submodel",
        "(V 70% lower in Caucasian vs African subjects, Results 3.3); not",
        "retained in the final model on the combined dataset."
      )
    )
  )

  compartmentData <- list(
    depot = list(analyte = "darunavir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "darunavir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 127,
    n_studies = 1,
    n_observations = 405,
    age_median = "55 years (IQR 13)",
    weight_median = "73 kg (IQR 17)",
    sex_female_pct = 33.1,
    race_ethnicity = "Caucasian 52.8%, African 43.3%, other 3.9% (Table 1)",
    disease_state = "HIV-1 infection, adult outpatients on boosted darunavir (median treatment duration 4.2 years; 78.7% with viral load < 40 copies/mL)",
    dose_range = "darunavir 800 mg q24h (91.3%), 600 mg q12h (7.9%) or 1200 mg q24h (0.8%); boosted with cobicistat (85.8%) or ritonavir (14.2%)",
    regions = "Belgium (Cliniques Universitaires Saint-Luc, Brussels)",
    genotypes = "CYP3A5 *1/*1 26.0%, *1/*3 26.0%, *3/*3 45.7%; SLCO3A1 rs8027174 G/T 14.2% (Table 2)",
    notes = paste(
      "Prospective observational TDM study (NCT03101644). Dataset A: 309",
      "sparse samples, one per visit at random post-intake times (1-29 h,",
      "median 11.9 h). Dataset B: 96 rich samples from 12 adherent patients",
      "(pre-dose, 0.5, 1, 2, 3, 4, 5, 6 h after witnessed intake). The",
      "final model was fitted to the pooled data with a NONMEM PRIOR on all",
      "fixed and random effects (except residual error) from the dataset-B",
      "fit (Methods 2.5). Baseline demographics in Table 1; genotypes in",
      "Table 2."
    )
  )

  ini({
    # Table 3 'Final population pharmacokinetic model', Estimate column
    # (NONMEM 7.4.3, combined dataset, after removal of 7 outliers).
    lka <- log(0.68)
    label("Absorption rate constant ka (1/h)") # Table 3, 'k a (h -1)' = 0.68 (RSE 73.5%)
    lcl <- log(12.9)
    label("Apparent clearance CL/F for a male CYP3A5 expresser at the AAG centring value (L/h)") # Table 3, 'CL/F (L.h -1)' = 12.9 (RSE 5.3%)
    lvc <- log(152)
    label("Apparent volume V/F for an SLCO3A1 G/G subject at the AAG centring value (L)") # Table 3, 'V/F (L)' = 152 (RSE 21.1%)

    # Covariate effects (Table 3 'Covariate model'); equation forms from
    # Methods 2.6 and ESM Online Resource 2.
    e_aag_cl <- -0.61
    label("Exponential AAG coefficient on CL/F (per g/L)") # Table 3, 'AAG on CL' = -0.61; ESM Online Resource 2 'Exponential'
    e_aag_vc <- -0.68
    label("Exponential AAG coefficient on V/F (per g/L)") # Table 3, 'AAG on V' = -0.68; ESM Online Resource 2 'Exponential'
    e_sexf_cl <- -0.21
    label("Exponential female-sex coefficient on CL/F (unitless)") # Table 3, 'Female sex on CL' = -0.21; ESM Online Resource 2 'Exponential'
    e_cyp3a5_cl <- -0.16
    label("Fractional change in CL/F for CYP3A5 nonexpressers (*3/*3) vs expressers (unitless)") # Table 3, 'CYP3A5*3 on CL' = -0.16; ESM Online Resource 2 'Categorical'
    e_snp_slco3a1_rs8027174_vc <- 0.81
    label("Fractional change in V/F for SLCO3A1 rs8027174 T carriers vs G/G (unitless)") # Table 3, 'SLCO3A1 G>T on V' = 0.81; ESM Online Resource 2 'Categorical'

    # IIV: exponential (log-normal) random effects, Table 3 reports the SD;
    # variance = SD^2. Diagonal (the block matrix was dropped for the
    # combined dataset, Results 3.2).
    etalka ~ 0.36 # Table 3, 'omega ka (SD)' = 0.60 -> 0.60^2
    etalcl ~ 0.0484 # Table 3, 'omega CL (SD)' = 0.22 -> 0.22^2
    etalvc ~ 0.1089 # Table 3, 'omega V (SD)' = 0.33 -> 0.33^2

    # Residual error: 'mixed (additive plus exponential)' (Results 3.2).
    propSd <- 0.281
    label("Exponential (proportional) residual SD (fraction)") # Table 3, 'sigma exponential (SD)' = 0.281
    addSd <- 0.641
    label("Additive residual SD (mg/L)") # Table 3, 'sigma additive (SD)' = 0.641
  })

  model({
    # AAG centring value (g/L). The source centres on the unreported cohort
    # median; 1 g/L is a maintainers' assumption (see covariateData$AAG).
    aag_ref <- 1

    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) *
      exp(e_aag_cl * (AAG - aag_ref)) *
      exp(e_sexf_cl * SEXF) *
      (1 + e_cyp3a5_cl * (1 - CYP3A5_EXPR))
    vc <- exp(lvc + etalvc) *
      exp(e_aag_vc * (AAG - aag_ref)) *
      (1 + e_snp_slco3a1_rs8027174_vc * SNP_SLCO3A1_RS8027174)

    kel <- cl / vc

    # One-compartment model with first-order absorption and elimination
    # (Results 3.2); no lag time in the final combined-dataset model.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
