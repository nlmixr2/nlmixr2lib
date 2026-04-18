Soehoel_2022_tralokinumab <- function() {
  description <- "Two-compartment population PK model for tralokinumab (Soehoel 2022) in adults with moderate-to-severe atopic dermatitis, with SC first-order absorption and allometric body-weight effects."
  reference <- "Soehoel A, Larsen MS, Timmermann S. Population Pharmacokinetics of Tralokinumab in Adult Subjects With Moderate to Severe Atopic Dermatitis. Clinical Pharmacology in Drug Development. 2022;11(8):910-921. doi:10.1002/cpdd.1113"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # From Table 2 footnotes
  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric-style effect on CL/Q and Vc/Vp with reference weight 75 kg.",
      source_name        = "WT"
    ),
    nonECZTRA = list(
      description        = "Indicator for non-ECZTRA trial enrollment",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ECZTRA trial)",
      notes              = "1 = any study other than ECZTRA; 0 = ECZTRA study. Mixed-case preserved from source per covariate-columns.md; future models should rename to NON_ECZTRA.",
      source_name        = "nonECZTRA"
    ),
    dilution = list(
      description        = "Indicator for diluted drug product (study D2213C00001)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not diluted)",
      notes              = "1 = drug diluted as in study D2213C00001; 0 = not diluted (typical). Lower-case preserved from source per covariate-columns.md; future models should rename to DILUTION.",
      source_name        = "dilution"
    )
  )

  population <- list(
    n_subjects     = 2561,
    n_studies      = 10,
    age_range      = "18-92 years",
    age_median     = "38 years (131 subjects >=65 years)",
    weight_range   = "36-165 kg",
    weight_median  = "74.5 kg",
    sex_female_pct = 44.9,
    race_ethnicity = "White 1721 (67%), Asian 560 (22%), Black/African American 183 (7%), Other 4%; Hispanic/Latino 222 (8.7%), Not Hispanic/Latino 2339 (91%)",
    disease_state  = "Pooled across indications: atopic dermatitis 2066 (81%), asthma 441 (17%), healthy 54 (2%); baseline EASI median 27.5 (range 12-72) in AD subjects",
    dose_range     = "Single-dose phase 1 through multi-dose phase 3; labelled AD regimen is a 600 mg SC loading dose followed by 300 mg SC every 2 weeks (Q4W extensions evaluated)",
    regions        = "Not explicitly stated in Soehoel 2022 Table 1; pooled from 10 multinational trials (3 phase 3 ECZTRA, 4 phase 2, 3 phase 1)",
    notes          = "Demographics from Soehoel 2022 Table 1 (pooled analysis population, n=2561). 2204 subjects (86%) enrolled in ECZTRA phase 3 atopic-dermatitis trials; 49 subjects (2%) received diluted drug product in study D2213C00001."
  )

  ini({
    lka <- log(0.184); label("Absorption rate (1/day)")
    lvc <- log(2.71); label("Central volume of distribution (L)")
    lcl <- log(0.149); label("Clearance (L/day)")
    lvp <- log(1.44); label("Peripheral volume of distribution (L)")
    lq <- log(0.159); label("Intercompartmental clearance (L/day)")
    lfdepot <- log(0.761); label("Subcutaneous bioavailability (fraction)")
    CcaddSd <- 0.238; label("Additive residual error (ug/mL)")
    CcpropSd <- 0.216; label("Proportional residual error (fraction)")

    e_wt_vcvp <- 0.783; label("Effect of body weight on central and peripheral volumes (unitless)")
    e_wt_clq <- 0.873; label("Effect of body weight on clearance and intercompartmental clearance (unitless)")
    e_nonECZTRA_cl <- 0.344; label("Effect of non-ECZTRA trials on clearance (unitless)")
    e_nonECZTRA_vc <- 0.258; label("Effect of non-ECZTRA trials on central volume (unitless)")
    e_f_dilution <- 0.354; label("Effect of dilution on bioavailability (unitless)")
    e_ka_dilution <- -0.519; label("Effect of dilution trials on absorption rate (unitless)")

    # Variance-covariance matrix for (etalvc, etalcl); Table 2 reports
    # IIV via CV% = sqrt(exp(omega^2) - 1) * 100%, so omega^2 = log(1 + CV^2).
    # CV_V2 = 40.1% -> omega^2_V2 = log(1 + 0.401^2) = 0.148971
    # CV_CL = 31.3% -> omega^2_CL = log(1 + 0.313^2) = 0.093459
    # correlation 0.61 -> cov = 0.61 * sqrt(0.148971 * 0.093459) = 0.071977
    etalvc + etalcl ~ c(0.148971,
                        0.071977, 0.093459)
  })
  model({
    fdepot <- exp(lfdepot)*(1 + e_f_dilution*dilution)
    ka <- exp(lka)*(1 + e_ka_dilution*dilution)
    cl <- exp(lcl + etalcl)*(WT/75)^e_wt_clq * (1 + e_nonECZTRA_cl*nonECZTRA)
    vc <- exp(lvc + etalvc)*(WT/75)^e_wt_vcvp * (1 + e_nonECZTRA_vc*nonECZTRA)
    q <- exp(lq)*(WT/75)^e_wt_clq
    vp <- exp(lvp)*(WT/75)^e_wt_vcvp

    # No unit conversion is required to change mg/L (dosing amount/central
    # volume unit) to ug/mL (measurement unit)
    Cc <- linCmt()
    f(depot) <- fdepot
    Cc ~ add(CcaddSd) + prop(CcpropSd)
  })
}

# IIV variance/covariance derivation (Soehoel 2022 Table 2 footnote: IIV was
# computed as sqrt(exp(omega^2) - 1), i.e. CV on the log-normal scale). In
# ini(), the "etalvc + etalcl ~ c(...)" form stores the variance-covariance
# matrix directly (variances on the diagonal, covariance off-diagonal):
#   omega^2_V2  = log(1 + 0.401^2) = 0.148971
#   omega^2_CL  = log(1 + 0.313^2) = 0.093459
#   cov(V2, CL) = 0.61 * sqrt(0.148971 * 0.093459) = 0.071977
#
# A prior revision of this file stored sqrt(.) values (SDs) in the
# variance-covariance triple and used an off-diagonal of
# sqrt(0.61 * omega_V2 * omega_CL), both of which were incorrect.
