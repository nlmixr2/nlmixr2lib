Schipani_2012_lopinavir <- function() {
  description <- "Population PK model for boosted lopinavir (lopinavir/ritonavir 400/100 mg) in HIV-infected adults from the Liverpool Therapeutic Drug Monitoring Registry. One-compartment with first-order absorption; apparent clearance is modified additively by body weight (deviation from median 72 kg) and by SLCO1B1 521T>C (rs4149056) genotype, encoded via the paired SLCO1B1_HAP15_HET / SLCO1B1_HAP15_HOM indicators (the source paper genotyped only 521T>C so *5- and *15-haplotype carriers are pooled, per the canonical's documented pooling rule)."
  reference <- paste(
    "Schipani A, Egan D, Dickinson L, Davies G, Boffito M, Youle M,",
    "Khoo SH, Back DJ, Owen A.",
    "Estimation of the effect of SLCO1B1 polymorphisms on lopinavir",
    "plasma concentration in HIV-infected adults.",
    "Antivir Ther. 2012;17(5):861-868.",
    "doi:10.3851/IMP2095.",
    sep = " "
  )
  vignette <- "Schipani_2012_lopinavir"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centered at the cohort median 72 kg. Paper covariate equation (p. 1311 Eq. 2): TVCL = theta0 + theta1 * (WT - WTmedian), where theta0 = 5.67 L/h is the typical CL/F at median weight (Table 1 Final Model). Range 45-117 kg.",
      source_name        = "BW"
    ),
    SLCO1B1_HAP15_HET = list(
      description        = "SLCO1B1 *15 haplotype heterozygote indicator: 1 if subject carries the 521C variant heterozygously (521TC), 0 otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (521TT homozygous wild-type)",
      notes              = "Time-fixed germline genotype. Schipani 2012 genotyped only the SLCO1B1 521T>C (rs4149056) SNP and did not phase the 388A>G (rs2306283) variant, so *5 (521T>C alone) and *15 (388A>G + 521T>C) carriers are pooled in this indicator; recorded under SLCO1B1_HAP15_HET per the canonical's documented pooling rule for Ide-style extractions (see inst/references/covariate-columns.md). Paper effect equation (p. 1311): CL = CL0 + theta_HET * HET + theta_HOM * HOM, with HET = 1 for 521TC. Distribution in the Liverpool cohort (Methods p. 863): 73 of 375 (19.5%) heterozygous 521TC.",
      source_name        = "HET (= 1 for 521TC)"
    ),
    SLCO1B1_HAP15_HOM = list(
      description        = "SLCO1B1 *15 haplotype homozygote indicator: 1 if subject carries the 521C variant homozygously (521CC), 0 otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (521TT homozygous wild-type)",
      notes              = "Time-fixed germline genotype. As for SLCO1B1_HAP15_HET, the Schipani 2012 cohort was genotyped at 521T>C only and the indicator pools *5/*5 and *15/*15 (and mixed *5/*15) homozygotes for the 521C allele. Paper effect equation (p. 1311): CL = CL0 + theta_HET * HET + theta_HOM * HOM, with HOM = 1 for 521CC. Distribution in the Liverpool cohort (Methods p. 863): 7 of 375 (1.9%) homozygous 521CC.",
      source_name        = "HOM (= 1 for 521CC)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 375L,
    n_observations = 594L,
    n_studies      = 1L,
    age_range      = "19-66 years",
    age_median     = "40 years",
    weight_range   = "45-117 kg",
    weight_median  = "72 kg",
    sex_female_pct = 18,
    race_ethnicity = "not collected (Discussion limitation: 'lack of ethnicity data')",
    disease_state  = "HIV-positive adults on lopinavir/ritonavir 400/100 mg tablets twice daily, all with HIV viral load <50 copies/mL at the time of sampling.",
    dose_range     = "Lopinavir/ritonavir 400/100 mg tablets twice daily (boosted lopinavir).",
    regions        = "United Kingdom (Liverpool Therapeutic Drug Monitoring Registry, Liverpool; external validation cohort from the Royal Free NHS Trust, London).",
    sampling_window = "Sparse TDM sampling at random time points post-dose; only a limited number of samples in the absorption phase. External validation set: 42 observations from 6 patients.",
    genotype_distribution = "SLCO1B1 521T>C (rs4149056), Methods p. 863: 295 of 375 (78%) 521TT homozygous wild-type, 73 (20%) 521TC heterozygous, 7 (2%) 521CC homozygous variant. Minor allele frequency 11%; SNP in Hardy-Weinberg equilibrium.",
    notes          = "TDM registry data; exclusion criteria were pregnancy, undetectable plasma lopinavir concentrations (suggesting non-adherence), and concomitant use of known enzyme inducers. LLOQ for lopinavir = 95 ng/mL (0.095 mg/L); observed plasma concentration range 0.114-22.432 mg/L. The paper's simulations of dose-reduction scenarios use a ritonavir-effect sequential model carried over from Schipani 2011 (PMID 22128223) and is not part of this packaged model."
  )

  ini({
    # One-compartment first-order absorption model (Schipani 2012 Table 1
    # Final Model). The paper's parameterisation is additive on CL:
    #   TVCL = CL0 + theta_WT * (WT - WTmedian) + theta_HET * HET + theta_HOM * HOM
    #   CLi  = TVCL * exp(eta_CL)                    (exponential IIV)
    # where CL0 = 5.67 L/h is the typical CL/F at the cohort median weight
    # of 72 kg in 521TT (wild-type) subjects. Body weight effect was
    # significant (dOFV = 20.5, p < 0.001); SLCO1B1 521T>C effect was
    # significant (dOFV = 15, p < 0.001). Age and sex were screened and
    # rejected during backward elimination (Results p. 864).

    # Absorption rate constant (apparent first-order, no lag-time, no
    # IIV retained: "lag time did not significantly improve the fit";
    # "inter-individual variability was supported only for apparent
    # clearance").
    lka <- log(0.20); label("Absorption rate constant ka (1/h)")  # Table 1 Final Model ka = 0.20 1/h (RSE 6%)

    # Typical CL/F at median weight (72 kg) for 521TT homozygous
    # wild-type (CL0 in the paper's covariate equation).
    lcl <- log(5.67); label("Typical apparent clearance CL/F at WT = 72 kg, 521TT (L/h)")  # Table 1 Final Model CL/F = 5.67 L/h (RSE 4%)

    # Apparent central volume of distribution. No IIV retained.
    lvc <- log(45.5); label("Apparent central volume V/F (L)")  # Table 1 Final Model V/F = 45.5 L (RSE 14%)

    # Additive body-weight effect on CL/F: change in L/h per kg
    # deviation from median 72 kg. Text on p. 864: "CL/F increased
    # by 0.5 L/h with body weight increases of 10 kg" (i.e.,
    # ~0.0457 L/h per kg). NB this is an additive effect on the
    # linear CL scale, not the canonical multiplicative allometric form.
    e_wt_cl <- 0.0457; label("Additive WT effect on CL (L/h per kg from 72 kg)")  # Table 1 Final Model factor BW on CL = 0.0457 (RSE 22%)

    # Additive SLCO1B1 521T>C genotype effects on CL/F (paired
    # heterozygote / homozygote indicators against the 521TT reference).
    # Both effects are negative: variant 521C carriage reduces hepatic
    # OATP1B1 uptake -> lower apparent CL/F.
    e_slco1b1_hap15_het_cl <- -0.791; label("Additive 521TC genotype effect on CL (L/h vs 521TT)")  # Table 1 Final Model factor T/C on CL = -0.791 (RSE 36%)
    e_slco1b1_hap15_hom_cl <- -2.09;  label("Additive 521CC genotype effect on CL (L/h vs 521TT)")  # Table 1 Final Model factor C/C on CL = -2.09 (RSE 29%)

    # IIV on CL/F only (exponential, theta_i = theta * exp(eta)).
    # Paper reports IIV CL/F = 37% (Table 1 Final Model, RSE 44%).
    # For an exponential / log-normal random effect with CV = 37%:
    #   omega^2 = log(1 + 0.37^2) = log(1.1369) = 0.1284
    etalcl ~ 0.1284  # Table 1 Final Model IIV CL/F = 37% CV (= log(1 + 0.37^2) = 0.1284)

    # Residual variability: proportional only ("residual variability
    # was best described by a purely proportional structure").
    propSd <- 0.331; label("Proportional residual error (fraction)")  # Table 1 Final Model proportional residual = 33.1% (RSE 13%)
  })
  model({
    # Individual PK parameters. No IIV on ka or V/F per the paper
    # (only CL/F carries an eta).
    ka <- exp(lka)
    vc <- exp(lvc)

    # Typical CL is additive in the covariates; individual CL is the
    # typical CL multiplied by exp(eta) (exponential IIV).
    tvcl <- exp(lcl) +
            e_wt_cl                  * (WT - 72) +
            e_slco1b1_hap15_het_cl   * SLCO1B1_HAP15_HET +
            e_slco1b1_hap15_hom_cl   * SLCO1B1_HAP15_HOM
    cl <- tvcl * exp(etalcl)

    kel <- cl / vc

    # ODE system: depot -> central -> elimination.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Plasma lopinavir concentration. Dose is in mg and Vc is in L,
    # so central/vc has units mg/L (= ug/mL), matching the paper's
    # plasma concentrations.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
