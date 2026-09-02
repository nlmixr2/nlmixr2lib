Vujkovic_2018_efavirenz <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption for",
    "oral efavirenz in HIV-infected black African adults initiating",
    "efavirenz-based antiretroviral therapy in Botswana (Vujkovic 2018).",
    "Only apparent oral clearance CL/F was estimated: the apparent volume",
    "of distribution (150 L) and the absorption rate constant (0.18 /h)",
    "were held fixed at the values obtained from the same cohort in the",
    "earlier Gross 2017 analysis. CL/F carries a three-level CYP2B6 516G>T",
    "(rs3745274) genotype effect (GG 9.439, GT 7.233, TT 4.033 L/h per",
    "70 kg) encoded as log-ratio multiplicative shifts on the 516GG",
    "wild-type reference, together with 0.75-exponent allometric scaling of",
    "CL/F to a 70 kg reference weight. Residual variability is combined",
    "additive plus proportional."
  )
  reference <- paste(
    "Vujkovic M, Bellamy SL, Zuppa AF, Gastonguay MR, Moorthy GS,",
    "Ratshaa B, Han X, Steenhoff AP, Mosepele M, Strom BL, Bisson GP,",
    "Aplenc R, Gross R (2018).",
    "Polymorphisms in cytochrome P450 are associated with extensive",
    "efavirenz pharmacokinetics and CNS toxicities in an HIV cohort in",
    "Botswana.",
    "The Pharmacogenomics Journal 18(5):678-688.",
    "doi:10.1038/s41397-018-0028-2.",
    "The fixed Vd and Ka values are inherited from the same cohort's",
    "earlier analysis, Gross R et al. (2017) AIDS 31(15):2107-2113,",
    "doi:10.1097/QAD.0000000000001593 (paper reference 13).",
    sep = " "
  )
  vignette <- "Vujkovic_2018_efavirenz"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(
      analyte = "efavirenz", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "efavirenz", units = "mg",
      specimen = "plasma", verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Drives 0.75-exponent allometric scaling of apparent oral clearance CL/F, normalised to a 70 kg reference weight. Vujkovic 2018 Methods 'Population PK Model Development' paragraph 2: 'Weight was included in an allometric model and fixed at 0.75 CL/F and normalized to a reference weight of 70 kg.' The apparent volume of distribution was held fixed at 150 L and is NOT weight-scaled in this model. The paper does not report the cohort body-weight distribution (Vujkovic 2018 Table 1 lists body mass index, median 22.0 kg/m2, IQR 19.8-25.1, but neither weight nor height), so downstream users must supply WT themselves; at WT = 70 kg the allometric factor is exactly 1 and the published typical CL/F values apply directly.",
      source_name        = "subject weight"
    ),
    SNP_CYP2B6_RS3745274_T_COUNT = list(
      description        = "Count of CYP2B6 c.516G>T (rs3745274, p.Q172H) T-alleles per subject (0/1/2). 0 = GG homozygous wild-type (normal metaboliser), 1 = GT heterozygous (slow), 2 = TT homozygous variant (very slow).",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (germline genotype). Vujkovic 2018 Results 'Population Pharmacokinetics' paragraph 1 reports the three genotype-specific typical apparent clearances directly: 'CYP2B6 516G>T genotypes demonstrated an impact on EFV clearance for normal (GG, 9.439 L/hr/70kg), slow (GT, 7.233 L/hr/70kg), and very slow (TT, 4.033 L/hr/70kg) genotypes respectively, with a MOFV of 2085.' The count column is decomposed in model() into two mutually exclusive indicators, ('count == 1') and ('count == 2'), following the Schipani 2011 / Siccardi 2012 precedent; the paper's estimates are not a per-allele dose-response (the TT shift is not twice the GT shift), so a single per-allele slope would not reproduce the published values. Cohort genotype counts, derived from Vujkovic 2018 Table 3 by summing over the 983T>C and extensive-SNP strata: GG n = 309 (39%), GT n = 375 (47%), TT n = 115 (14%); 814 of 941 enrolled participants were successfully genotyped for rs3745274.",
      source_name        = "CYP2B6 516G>T (rs3745274)"
    )
  )

  covariatesDataExcluded <- list(
    SNP_CYP2B6_RS28399499_C_COUNT = list(
      description = "Count of CYP2B6 c.983T>C (rs28399499) C-alleles per subject (0/1/2)",
      units       = "(count, 0/1/2)",
      type        = "continuous",
      notes       = "Screened univariately on CL/F and strongly significant (Vujkovic 2018 Table 2: wild-type 7.984, heterozygote 4.691, homozygote 1.440 L/h per 70 kg; MOFV 1862 vs 2075 for the pre-SNP model), but the paper did NOT carry it into a joint multivariate PK model. Its role in the paper is to define the composite CYP2B6 metaboliser group (Table 3) used in the downstream pharmacodynamic analyses, not to extend the PK model. Carrying the univariate estimate alongside the retained rs3745274 effect would double-count the shared CYP2B6 signal; see the vignette Assumptions and deviations section."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested as a covariate on CL/F and not retained. Vujkovic 2018 Results 'Population Pharmacokinetics' paragraph 1: 'History of tuberculosis, gender, and drug adherence did not significantly improve MOFV'. No point estimate is reported."
    ),
    TB_HX = list(
      description = "Suspected or confirmed tuberculosis co-occurrence (1 = yes, 0 = no)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested as a covariate on CL/F and not retained (Vujkovic 2018 Results 'Population Pharmacokinetics' paragraph 1). Prevalence in the cohort was 3.8% (31 of 941; Vujkovic 2018 Table 1). No point estimate is reported. Name is descriptive only -- no canonical register entry was created because the covariate is not used in model()."
    ),
    ADHERENCE = list(
      description = "Efavirenz medication possession ratio over the first month of therapy",
      units       = "(fraction)",
      type        = "continuous",
      notes       = "Tested as a covariate on CL/F and not retained (Vujkovic 2018 Results 'Population Pharmacokinetics' paragraph 1). Derived from pharmacy refill data as (doses dispensed - doses returned) / (days between the initial fill and the fill closest to the target date); Vujkovic 2018 Methods 'Medication Adherence'. No point estimate is reported. Name is descriptive only -- no canonical register entry was created because the covariate is not used in model()."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 742L,
    n_studies      = 1L,
    age_range      = "21 years and older (inclusion criterion)",
    age_median     = "37 years (IQR 33-44; Vujkovic 2018 Table 1)",
    weight_range   = "not reported",
    weight_median  = "not reported; body mass index median 22.0 kg/m2 (IQR 19.8-25.1), Vujkovic 2018 Table 1. The model's allometric reference weight is 70 kg.",
    sex_female_pct = 51.0,
    race_ethnicity = "black African (study inclusion criterion (ii): 'black African origin'); all participants recruited in and around Gaborone, Botswana",
    disease_state  = "HIV-1 infection, antiretroviral-naive at enrolment, initiating a first three-drug regimen containing efavirenz plus two nucleoside reverse-transcriptase inhibitors. Median baseline CD4 count 196 cells/mm3 (IQR 112-256) and median plasma viral load 4.9 log10 copies/mL (IQR 4.2-5.4); 3.8% had a history of tuberculosis (Vujkovic 2018 Table 1).",
    dose_range     = "efavirenz 600 mg once daily (fixed dose; study inclusion criterion (v))",
    regions        = "Botswana (clinics in and around Gaborone), enrolment June 2009 to November 2013",
    notes          = "941 participants were enrolled, 814 were successfully genotyped for CYP2B6 516G>T (rs3745274), and 742 of those contributed at least one steady-state plasma efavirenz concentration at month 1; 562 also contributed a month-6 sample (Vujkovic 2018 Results 'Population Pharmacokinetics' paragraph 1 and Table 1). Samples were single midpoint draws per visit collected once steady state had been reached (efavirenz reaches steady state within 6-10 days); the date and time of both the blood draw and the last dose were recorded, so the analysis dataset carries an observed time after dose rather than a nominal one. Assay was HPLC-MS/MS in negative ionisation mode with multiple reaction monitoring, lower limit of quantification 1 ng/mL (Vujkovic 2018 Methods 'Drug assay'). Median observed plasma efavirenz was 2.17 ug/mL (IQR 1.62-3.88) at month 1 and 2.05 ug/mL at month 6. Fitted in NONMEM VII with ADVAN2 and FOCE-I with interaction."
  )

  ini({
    # =========================================================================
    # Structural disposition. Vujkovic 2018 Methods 'Population PK Model
    # Development' paragraph 1: 'The baseline EFV PK model parameters were
    # adopted from a previous study in the same cohort, where clearance
    # (CL/F) was the only parameter to be estimated while holding the
    # parameters of volume of distribution (Vd) and absorption rate
    # constant (Ka) constant at 150 L and 0.18 h-1 respectively [13].'
    # Reference 13 is Gross 2017, AIDS 31(15):2107-2113.
    #
    # Because Vd and Ka are inherited constants rather than estimates from
    # this dataset, the model reproduces steady-state EXPOSURE (CL/F drives
    # AUC and hence the midpoint concentrations that were fitted) far more
    # reliably than the within-interval concentration-time SHAPE. See the
    # vignette Assumptions and deviations section.
    # =========================================================================
    lcl <- log(9.439)
    label("Apparent oral clearance CL/F for the CYP2B6 516GG wild-type reference (L/h per 70 kg)")
    # Vujkovic 2018 Results 'Population Pharmacokinetics' paragraph 1: 'normal (GG, 9.439 L/hr/70kg)'

    lvc <- fixed(log(150))
    label("Apparent volume of distribution Vd/F (L), inherited from Gross 2017 (paper reference 13)")
    # Vujkovic 2018 Methods 'Population PK Model Development' paragraph 1: 'holding the parameters of volume of distribution (Vd) and absorption rate constant (Ka) constant at 150 L and 0.18 h-1'

    lka <- fixed(log(0.18))
    label("First-order absorption rate constant ka (1/h), inherited from Gross 2017 (paper reference 13)")
    # Vujkovic 2018 Methods 'Population PK Model Development' paragraph 1 (same sentence as lvc)

    # =========================================================================
    # CYP2B6 516G>T (rs3745274) effect on CL/F. The paper prints the three
    # genotype-specific typical clearances rather than the two shifts, so
    # the shifts below are exact log-ratios of published values against the
    # 516GG reference and reproduce all three numbers exactly:
    #   GG: exp(lcl)                    = 9.439 L/h per 70 kg
    #   GT: exp(lcl + e_516gt_cl)       = 7.233 L/h per 70 kg
    #   TT: exp(lcl + e_516tt_cl)       = 4.033 L/h per 70 kg
    # A single per-allele slope cannot reproduce these (the TT log-shift is
    # -0.8503, not twice the GT log-shift of -0.2662), so the count column
    # is decomposed into two mutually exclusive indicators in model(),
    # following the Schipani 2011 / Siccardi 2012 precedent.
    # =========================================================================
    e_516gt_cl <- log(7.233 / 9.439)
    label("Log-ratio of CYP2B6 516GT heterozygote CL/F vs the 516GG wild-type reference (unitless)")
    # Vujkovic 2018 Results 'Population Pharmacokinetics' paragraph 1: 'slow (GT, 7.233 L/hr/70kg)'; log(7.233/9.439) = -0.2662

    e_516tt_cl <- log(4.033 / 9.439)
    label("Log-ratio of CYP2B6 516TT homozygous-variant CL/F vs the 516GG wild-type reference (unitless)")
    # Vujkovic 2018 Results 'Population Pharmacokinetics' paragraph 1: 'very slow (TT, 4.033 L/hr/70kg)'; log(4.033/9.439) = -0.8503

    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on CL/F with body weight (unitless)")
    # Vujkovic 2018 Methods 'Population PK Model Development' paragraph 2: 'Weight was included in an allometric model and fixed at 0.75 CL/F and normalized to a reference weight of 70 kg.'

    # =========================================================================
    # Inter-individual variability on CL/F. Vujkovic 2018 Results
    # 'Population Pharmacokinetics' paragraph 1 reports the final baseline
    # covariate model 'with an inter-individual variability (omega2 CL) of
    # 40.5% coefficient variation'. The stated UNIT (% coefficient of
    # variation) is taken over the omega-squared SYMBOL, so 40.5% is read as
    # a CV and converted to the log-scale variance nlmixr2 expects:
    # omega^2 = log(1 + CV^2) = log(1 + 0.405^2) = 0.15189.
    # Vujkovic 2018 Methods 'Population PK Model Development' paragraph 1:
    # 'We used an exponential variance model to describe the unexplained
    # random variability of parameters across individuals.'
    #
    # The paper also retained inter-occasion variability on CL/F between the
    # month-1 and month-6 visits (Table 2 footnote: 'Clearance estimates are
    # corrected for the covariates: CYP2B6 G516T genotype (rs3745274),
    # allometrically scaled weight, and inter-occasional variability between
    # visits'). Its magnitude is never reported for the final model, so it
    # is NOT encoded here -- consistent with the between-occasion-variability
    # convention already used in Bienczak_2016_efavirenz.R and
    # Svensson_2018_bedaquiline.R when a between-subject term is reported on
    # the same parameter. See the vignette Assumptions and deviations
    # section.
    # =========================================================================
    etalcl ~ 0.15189
    # Vujkovic 2018 Results 'Population Pharmacokinetics' paragraph 1, IIV on CL 40.5% CV; log(1 + 0.405^2) = 0.15189

    # =========================================================================
    # Residual unexplained variability. Vujkovic 2018 Results 'Population
    # Pharmacokinetics' paragraph 1: '... proportional residual variability
    # of 0.0831% coefficient variation, and residual additive variability of
    # 0.0699 SD.' Both numbers are taken on the SD scale that nlmixr2's
    # prop() / add() expect: 0.0831 as a proportional CV (8.31%) and 0.0699
    # as an additive SD in mg/L. Vujkovic 2018 Methods 'Population PK Model
    # Development' paragraph 2: 'A composite additive and proportional error
    # model was used to describe random residual variability.'
    # =========================================================================
    propSd <- 0.0831
    label("Proportional residual error (fraction)")
    # Vujkovic 2018 Results 'Population Pharmacokinetics' paragraph 1: 'proportional residual variability of 0.0831'

    addSd <- 0.0699
    label("Additive residual error (mg/L)")
    # Vujkovic 2018 Results 'Population Pharmacokinetics' paragraph 1: 'residual additive variability of 0.0699 SD'
  })

  model({
    # --- 1. CYP2B6 516G>T genotype indicators -------------------------------
    # Mutually exclusive; the 516GG wild-type reference is implicit (both
    # indicators zero).
    is_516gt <- (SNP_CYP2B6_RS3745274_T_COUNT == 1)
    is_516tt <- (SNP_CYP2B6_RS3745274_T_COUNT == 2)

    # --- 2. Individual PK parameters ----------------------------------------
    # Apparent oral clearance: 516GG typical CL/F shifted multiplicatively by
    # genotype, scaled allometrically to a 70 kg reference weight, with
    # exponential between-subject variability on the log scale.
    cl <- exp(lcl + e_516gt_cl * is_516gt + e_516tt_cl * is_516tt + etalcl) *
      (WT / 70)^e_wt_cl

    # Apparent volume of distribution and absorption rate constant: inherited
    # constants, not weight-scaled and carrying no reported variability.
    vc <- exp(lvc)
    ka <- exp(lka)

    # --- 3. Micro-constants --------------------------------------------------
    kel <- cl / vc

    # --- 4. ODE system -------------------------------------------------------
    # NONMEM ADVAN2: one-compartment linear disposition with first-order
    # absorption from an oral depot. Doses are given to `depot`
    # (cmt = "depot"); bioavailability is subsumed into the apparent
    # parameters CL/F and Vd/F and is not separately identifiable.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # --- 5. Observation and residual error -----------------------------------
    # Dose in mg / volume in L -> concentration in mg/L, which is numerically
    # identical to the ug/mL used in Vujkovic 2018 Table 1.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
