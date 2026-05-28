Zhao_2018_omeprazole <- function() {
  description <- "Population PK-pharmacogenetic model for oral omeprazole and its two metabolites 5-hydroxy-omeprazole and omeprazole sulfone in Caucasian neonates and young infants (Zhao 2018). One-compartment parent disposition with first-order absorption (Ka modulated by ABCB1 C3435T genotype) is followed by parallel formation into two one-compartment metabolites with apparent volume V_M/F fixed to 1 L; the omeprazole-to-5-hydroxy-omeprazole formation clearance (CLOMZ-M1) is modulated by CYP2C19 metabolizer phenotype (poor / intermediate / extensive-or-ultrarapid) and a postnatal-age power function, while the omeprazole-to-omeprazole-sulfone formation clearance (CLOMZ-M2) and the metabolite apparent eliminations carry no covariates. Linear omeprazole elimination was estimated as negligible (< 0.0001 L/h) and is therefore not included in the final structural model."
  reference   <- "Zhao W, Leroux S, Biran V, Jacqz-Aigrain E. Developmental pharmacogenetics of CYP2C19 in neonates and young infants: omeprazole as a probe drug. Br J Clin Pharmacol. 2018;84(5):997-1005. doi:10.1111/bcp.13526"
  vignette    <- "Zhao_2018_omeprazole"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    PNA = list(
      description        = "Postnatal age (chronological since birth).",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Zhao 2018 reports postnatal age in DAYS (cohort median 38 days, range 7-87 days). The canonical PNA column is in MONTHS, so the model `model()` block converts to months internally via PNA / 30.4375 only at the reference-ratio step: the paper's `F_PNA = (PNA_days / 38)^0.472` is reparameterised as `F_PNA = (PNA_months / 1.249)^0.472` because both numerator and denominator carry the same units factor and cancel. Users should supply PNA in months in the dataset.",
      source_name        = "PNA"
    ),
    CYP2C19_IM = list(
      description        = "CYP2C19 intermediate-metabolizer phenotype indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (extensive or ultrarapid metabolizer; *1/*1, *1/*17, *17/*17 -- both CYP2C19_IM = 0 and CYP2C19_PM = 0)",
      notes              = "1 = subject has CYP2C19 IM phenotype (*1/*2 or *2/*17 in Zhao 2018; one functional and one loss-of-function allele); 0 = otherwise. Cohort distribution: IM 21.6% (11/51), PM 3.9% (2/51), EM/UM reference 74.5% (38/51). Source paper pooled EM and UM into a single reference because typical-value CLOMZ-M1 was indistinguishable between the two strata.",
      source_name        = "CYP2C19_IM"
    ),
    CYP2C19_PM = list(
      description        = "CYP2C19 poor-metabolizer phenotype indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (extensive, ultrarapid, or intermediate metabolizer)",
      notes              = "1 = subject has CYP2C19 PM phenotype (*2/*2 in Zhao 2018; two loss-of-function alleles); 0 = otherwise. Paired with `CYP2C19_IM` to encode the three-level EM/UM (reference) / IM / PM phenotype with two binary indicators.",
      source_name        = "CYP2C19_PM"
    ),
    ABCB1_C3435T_HET = list(
      description        = "ABCB1 (rs1045642) C3435T heterozygote indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (C/C homozygous wild-type, when paired with `ABCB1_C3435T_MUT = 0`)",
      notes              = "1 = subject has ABCB1 C/T genotype at rs1045642; 0 = otherwise (C/C wild-type or T/T homozygous variant). Cohort distribution: C/C 49.0% (25/51), C/T 43.1% (22/51), T/T 7.8% (4/51).",
      source_name        = "ABCB1_C3435T_HET"
    ),
    ABCB1_C3435T_MUT = list(
      description        = "ABCB1 (rs1045642) C3435T homozygous-variant indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (C/C homozygous wild-type, when paired with `ABCB1_C3435T_HET = 0`)",
      notes              = "1 = subject has ABCB1 T/T genotype at rs1045642; 0 = otherwise. Paired with `ABCB1_C3435T_HET` to encode the three-level C/C (reference) / C/T / T/T genotype with two binary indicators.",
      source_name        = "ABCB1_C3435T_MUT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 51L,
    n_studies      = 1L,
    n_observations = "73 omeprazole + paired 5-hydroxy-omeprazole and omeprazole-sulfone plasma concentrations (Zhao 2018 Results paragraph 2). LLOQ 10 ng/mL omeprazole, 25 ng/mL each metabolite; sub-LLOQ values imputed as LLOQ/2 in the source modelling.",
    age_range      = "Neonates and young infants, PNA 7-87 days (gestational age at birth 24-41 weeks; postmenstrual age covers preterm-newborn through young-infant range)",
    age_median     = "PNA median 38 days (Table 1)",
    weight_range   = "0.78-3.80 kg (current weight; birth gestational age range 24-41 weeks)",
    weight_median  = "2.13 kg (Table 1)",
    sex_female_pct = 52.9,
    race_ethnicity = c(White = 100),
    disease_state  = "Caucasian neonates and young infants with gastroesophageal reflux enrolled in an omeprazole dose-finding study (ClinicalTrials.gov NCT01657578).",
    dose_range     = "Oral omeprazole 1-3 mg/kg/day once daily in the morning (cohort median 3.6 mg/day, range 0.8-7.3 mg/day; per-kg median 2.0 mg/kg/day, range 1.0-2.5 mg/kg/day).",
    regions        = "France (Robert Debre University Hospital, Paris; Comite de Protection des Personnes Saint Louis ethics committee).",
    notes          = "Demographics from Zhao 2018 Table 1. Pharmacokinetic samples obtained 0.5-4 h and 4-12 h after the first dose. Genotyping by TaqMan allelic discrimination for CYP2C19*2 (rs4244285), CYP2C19*17 (rs12248560), and ABCB1 rs1045642 / rs2032582 / rs1128503; the final model uses only the rs1045642 C3435T SNP on Ka (the other ABCB1 SNPs were tested but not retained). The model is fit to a single-dose pharmacokinetic dataset; multiple-dose extrapolation should consider that ABCB1 ontogeny and CYP2C19 ontogeny were observed to evolve over the first months of life so steady-state predictions outside the studied PNA range carry extrapolation risk."
  )

  ini({
    # Final-model estimates from Zhao 2018 Table 2 ("Final model" column).
    # Bootstrap medians (5th-95th percentile) shown in Table 2 are
    # comparable to the final estimates and indicate model stability.
    # The source paper reports IIV as %CV on an exponential model
    # eta_i ~ N(0, omega^2); the conversion to log-normal variance for
    # nlmixr2 is omega^2 = log(1 + CV^2). Residual error is additive on
    # the linear scale with a single SD applied to all three outputs
    # (Zhao 2018 Methods "Population pharmacokinetic-pharmacogenetic
    # modelling" paragraph 3 and Table 2 single-row "Residual additive").

    # Structural parameters -- parent omeprazole.
    lka     <- log(0.0497);  label("Absorption rate constant Ka at ABCB1 C/C wild-type reference (1/h)")  # Zhao 2018 Table 2: Ka = theta1 = 0.0497 1/h (RSE 29.4%)
    lvc     <- log(0.513);   label("Apparent central volume of distribution V1/F of omeprazole (L)")     # Zhao 2018 Table 2: V1/F = 0.513 L (RSE 25.2%)

    # Structural parameters -- 5-hydroxy-omeprazole (5-OH-OMZ) metabolite.
    # V2/F fixed to 1 L; the formation parameter CLOMZ-M1 = theta2 * F_PNA
    # * F_CYP2C19 was estimated. The paper interprets CLOMZ-M1 as the
    # ratio CL_form_5oh / V_M (with V_M = 1 L it is numerically equal to
    # the formation clearance in L/h); see vignette Assumptions and
    # deviations for the dimensional analysis. The eta on the formation
    # parameter pairs with `lkmet_5oh` as the structural name; the
    # `kmet` paper-named-parameter token is used because it canonically
    # tags a metabolite-formation rate / formation-clearance scalar in
    # parent-plus-metabolite popPK models.
    lvc_5oh    <- fixed(log(1));   label("Apparent volume of distribution V2/F of 5-hydroxy-omeprazole (L) -- fixed")   # Zhao 2018 Table 2: V2/F = 1 L FIX
    lkmet_5oh  <- log(0.658);      label("Formation parameter CLOMZ-M1 reference (theta2) for 5-hydroxy-omeprazole (L/h at PNA = 38 days / 1.249 months and EM/UM reference)")  # Zhao 2018 Table 2: theta2 = 0.658 L/h (RSE 21.4%)
    lcl_5oh    <- log(0.846);      label("Apparent elimination clearance CLM1/F of 5-hydroxy-omeprazole (L/h)")         # Zhao 2018 Table 2: CLM1/F = 0.846 L/h (RSE 18.4%)

    # Structural parameters -- omeprazole sulfone (M2) metabolite.
    lvc_sfn    <- fixed(log(1));   label("Apparent volume of distribution V3/F of omeprazole sulfone (L) -- fixed")     # Zhao 2018 Table 2: V3/F = 1 L FIX
    lkmet_sfn  <- log(0.140);      label("Formation parameter CLOMZ-M2 for omeprazole sulfone (L/h)")                   # Zhao 2018 Table 2: CLOMZ-M2 = 0.140 L/h (RSE 8.4%)
    lcl_sfn    <- log(0.130);      label("Apparent elimination clearance CLM2/F of omeprazole sulfone (L/h)")           # Zhao 2018 Table 2: CLM2/F = 0.130 L/h (RSE 27.3%)

    # Covariate effects.
    # Ka shifters for ABCB1 C3435T (rs1045642) genotype; multiplicative
    # via power-of-binary-indicator (the indicator equals 0 or 1, so
    # `e^IND = 1` if IND = 0 and `e` if IND = 1).
    e_abcb1_c3435t_het_ka <- 1.86;   label("ABCB1 C/T heterozygote multiplicative factor on Ka (unitless)")     # Zhao 2018 Table 2: F_ABCB1 HET = 1.86 (RSE 44.4%)
    e_abcb1_c3435t_mut_ka <- 6.93;   label("ABCB1 T/T homozygous-variant multiplicative factor on Ka (unitless)") # Zhao 2018 Table 2: F_ABCB1 MUT = 6.93 (RSE 44.6%)

    # CLOMZ-M1 covariate effects.
    # CYP2C19 phenotype multiplicative factors on the formation
    # parameter for 5-hydroxy-omeprazole (EM/UM reference).
    e_cyp2c19_im_kmet_5oh <- 0.449;  label("CYP2C19 IM multiplicative factor on CLOMZ-M1 (unitless)")              # Zhao 2018 Table 2: F_CYP2C19 IM = 0.449 (RSE 25.6%)
    e_cyp2c19_pm_kmet_5oh <- 0.125;  label("CYP2C19 PM multiplicative factor on CLOMZ-M1 (unitless)")              # Zhao 2018 Table 2: F_CYP2C19 PM = 0.125 (RSE 44.5%)
    e_pna_kmet_5oh        <- 0.472;  label("Postnatal-age power exponent on CLOMZ-M1 (unitless; PNA reference 38 days = 1.249 months)")  # Zhao 2018 Table 2: F_PNA = (PNA / 38)^theta3, theta3 = 0.472 (RSE 29.2%)

    # IIV (exponential model). Source paper reports %CV in Table 2;
    # log-normal variance via omega^2 = log(1 + CV^2). No correlation
    # block is reported in the paper (independent etas).
    # Ka: log(1 + 1.30^2) = 0.98954
    etalka         ~ 0.98954
    # V1/F: log(1 + 0.99^2) = 0.68318
    etalvc         ~ 0.68318
    # CLOMZ-M1: log(1 + 0.526^2) = 0.24432
    etalkmet_5oh   ~ 0.24432
    # CLOMZ-M2: log(1 + 0.359^2) = 0.12127
    etalkmet_sfn   ~ 0.12127
    # CLM1/F: log(1 + 0.559^2) = 0.27182
    etalcl_5oh     ~ 0.27182
    # CLM2/F: log(1 + 0.688^2) = 0.38735
    etalcl_sfn     ~ 0.38735

    # Residual error. Zhao 2018 Table 2 reports a single additive
    # residual SD of 55.6 ng/mL applied to all three outputs (combined
    # row). Conversion to model units (mg/L): 55.6 ng/mL = 0.0556 mg/L.
    # The same value seeds the three per-output residual-error
    # parameters; in the original NONMEM fit a single SIGMA(1,1)
    # multiplexed across all observations.
    addSd     <- 0.0556;   label("Additive residual SD for omeprazole (mg/L)")               # Zhao 2018 Table 2: 55.6 ng/mL = 0.0556 mg/L (RSE 50.2%)
    addSd_5oh <- 0.0556;   label("Additive residual SD for 5-hydroxy-omeprazole (mg/L)")     # Zhao 2018 Table 2: same shared 55.6 ng/mL value
    addSd_sfn <- 0.0556;   label("Additive residual SD for omeprazole sulfone (mg/L)")       # Zhao 2018 Table 2: same shared 55.6 ng/mL value
  })

  model({
    # ----- Derived covariate terms -----
    # PNA enters as a power function relative to the cohort reference.
    # Zhao 2018 expresses PNA in days (reference 38 d); the canonical
    # nlmixr2lib PNA column is in months, so the reference is rescaled:
    # 38 days / 30.4375 days/month = 1.249 months. The exponent and
    # the dynamic relationship are unchanged.
    pna_ref_months <- 38 / 30.4375
    f_pna <- (PNA / pna_ref_months) ^ e_pna_kmet_5oh

    # ABCB1 C3435T multiplicative factor on Ka.
    f_abcb1 <- (e_abcb1_c3435t_het_ka ^ ABCB1_C3435T_HET) *
               (e_abcb1_c3435t_mut_ka ^ ABCB1_C3435T_MUT)

    # CYP2C19 phenotype multiplicative factor on CLOMZ-M1.
    f_cyp2c19 <- (e_cyp2c19_im_kmet_5oh ^ CYP2C19_IM) *
                 (e_cyp2c19_pm_kmet_5oh ^ CYP2C19_PM)

    # ----- Individual parameters -----
    ka       <- exp(lka + etalka) * f_abcb1
    vc       <- exp(lvc + etalvc)
    vc_5oh   <- exp(lvc_5oh)
    vc_sfn   <- exp(lvc_sfn)
    kmet_5oh <- exp(lkmet_5oh + etalkmet_5oh) * f_pna * f_cyp2c19
    kmet_sfn <- exp(lkmet_sfn + etalkmet_sfn)
    cl_5oh   <- exp(lcl_5oh + etalcl_5oh)
    cl_sfn   <- exp(lcl_sfn + etalcl_sfn)

    # ----- ODE system -----
    # Parent omeprazole loss to the two metabolites; no direct linear
    # omeprazole elimination (CLOMZ/F was estimated as negligible,
    # < 0.0001 L/h, and is excluded from the final structural model
    # per Zhao 2018 Results paragraph 2). The formation rate of the
    # metabolite per unit time is `kmet_<metab> * vc_<metab> * (central
    # / vc)` (units: L/h * L * mg/L / L = mg/h); with V_M = 1 L the
    # `vc_<metab>` factor is numerically 1 and the expression collapses
    # to `kmet_<metab> * (central / vc)`. The explicit `vc_<metab>`
    # factor is retained so dimensional analysis is transparent and the
    # equation remains correct should V_M ever be relaxed.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          (kmet_5oh * vc_5oh + kmet_sfn * vc_sfn) * (central / vc)
    d/dt(central_5oh) <-  kmet_5oh * vc_5oh * (central / vc) -
                          cl_5oh * (central_5oh / vc_5oh)
    d/dt(central_sfn) <-  kmet_sfn * vc_sfn * (central / vc) -
                          cl_sfn * (central_sfn / vc_sfn)

    # ----- Observations -----
    Cc     <- central / vc
    Cc_5oh <- central_5oh / vc_5oh
    Cc_sfn <- central_sfn / vc_sfn

    Cc     ~ add(addSd)
    Cc_5oh ~ add(addSd_5oh)
    Cc_sfn ~ add(addSd_sfn)
  })
}
