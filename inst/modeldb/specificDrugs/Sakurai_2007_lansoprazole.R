Sakurai_2007_lansoprazole <- function() {
  description <- "Two-compartment population PK model for intravenously administered lansoprazole in 56 healthy Japanese adult males (Sakurai 2007). Volumes (V1, V2) and clearances (CL, Q) scale linearly with body weight via per-kg reference values; systemic clearance is stratified by CYP2C19 metabolizer phenotype using two binary indicators (homoEM reference; heteroEM and PM groups carry multiplicative factors of 0.612 and 0.212 respectively). Inter-individual variability is log-normal on V1, CL, V2 (no IIV on Q); residual error is combined proportional plus additive."
  reference   <- "Sakurai Y, Hirayama M, Hashimoto M, Tanaka T, Hasegawa S, Irie S, Ashida K, Kayano Y, Taguchi M, Hashimoto Y. Population pharmacokinetics and proton pump inhibitory effects of intravenous lansoprazole in healthy Japanese males. Biol Pharm Bull. 2007;30(12):2238-2243. doi:10.1248/bpb.30.2238"
  vignette    <- "Sakurai_2007_lansoprazole"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight (current).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-kg linear scaling on V1, V2, CL, Q (Sakurai 2007 Equations 1-4). Cohort mean (S.D.) 61.8 (6.0) kg. The per-kg parameter values reported in Table 1 are multiplied by the subject's body weight to obtain individual totals; no reference / normalisation weight is used.",
      source_name        = "WT"
    ),
    CYP2C19_IM = list(
      description        = "CYP2C19 intermediate-metabolizer phenotype indicator (Sakurai 2007 heteroEM group).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (extensive metabolizer; CYP2C19*1/*1, when both CYP2C19_IM = 0 and CYP2C19_PM = 0)",
      notes              = "1 = subject is a CYP2C19 intermediate metabolizer (Sakurai 2007 heteroEM group, CYP2C19*1/*2 or *1/*3, one functional and one loss-of-function allele); 0 = otherwise. Paired with CYP2C19_PM to encode the three-level EM (homoEM, reference) / IM (heteroEM) / PM phenotype with two binary indicators. Cohort distribution (Methods 'Subjects' and Results paragraph 2): homoEM 28.6% (16/56), heteroEM 57.1% (32/56), PM 14.3% (8/56). Time-fixed per subject (germline CYP2C19 genotype determined by allele-specific PCR or PCR-RFLP).",
      source_name        = "CYP2C19_IM"
    ),
    CYP2C19_PM = list(
      description        = "CYP2C19 poor-metabolizer phenotype indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (extensive or intermediate metabolizer)",
      notes              = "1 = subject is a CYP2C19 poor metabolizer (Sakurai 2007 PM group, CYP2C19*2/*2, *2/*3, or *3/*3; two loss-of-function alleles); 0 = otherwise. Paired with CYP2C19_IM to encode the three-level EM (reference) / IM / PM phenotype with two binary indicators. Time-fixed per subject (germline genotype).",
      source_name        = "CYP2C19_PM"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 56L,
    n_studies       = 1L,
    n_observations  = "1069 serum lansoprazole concentrations across single-dose and multiple-dose intravenous trials (Sakurai 2007 Results paragraph 2). HPLC-UV assay with lower limit of quantification 10 ng/mL in 500 uL of human serum; linearity 5-2000 ng/mL.",
    age_range       = "20-34 years",
    age_median      = "24.0 years (mean +/- S.D.: 24.0 +/- 4.0)",
    weight_range    = "approximately 55.8-67.8 kg (1 SD around the mean; individual extremes not reported)",
    weight_median   = "61.8 kg (mean +/- S.D.: 61.8 +/- 6.0)",
    sex_female_pct  = 0,
    race_ethnicity  = c(Japanese = 100),
    disease_state   = "Healthy Japanese adult males enrolled in single-dose and multiple-dose intravenous lansoprazole pharmacokinetic trials.",
    dose_range      = "Lansoprazole 15 mg or 30 mg by 30-min intravenous drip infusion; 30 mg as a single intravenous bolus; 30 mg twice-daily 30-min infusions for 5 days in a 16-subject multiple-dose subset.",
    regions         = "Japan (Ohsaki Clinic [Tokyo]; Sekino Clinical Pharmacology Clinic [Tokyo]; Kyushu Clinical Pharmacology Research Clinic [Fukuoka]).",
    notes           = "All subjects were genotyped for CYP2C19*1, *2, *3 by allele-specific PCR (SNP Typing Kit, Toyobo) or PCR-RFLP. Cohort genotype counts: 16 homoEM (CYP2C19*1/*1), 32 heteroEM (*1/*2 or *1/*3), 8 PM (*2/*2, *2/*3, or *3/*3). Multiple study arms were pooled in the popPK fit: single 30-min IV infusion of 30 mg (all 56 subjects); twice-daily 30-min IV infusion of 30 mg for 5 d (16 subjects: 8 PMs, 7 heteroEMs, 1 homoEM); single 30 mg IV bolus (8 subjects: 3 heteroEMs, 5 homoEMs); single 30-min IV infusion of 15 mg (8 subjects: 5 heteroEMs, 3 homoEMs). NONMEM V (level 1.1) on ThinkCentre A50p; first-order conditional estimation with ADVAN3 and TRANS4 PREDPP subroutines."
  )

  ini({
    # All parameter estimates from Sakurai 2007 Table 1, "2-Compartment" column.
    # Per-kg reference values are multiplied by individual WT inside model() to
    # obtain total V1, V2, CL, Q (Sakurai 2007 Equations 1-4). The paper reports
    # the residual additive SD in ng/mL (Table 1: 6.40 ng/mL); it is converted to
    # mg/L here for unit consistency with the rest of the file (6.40 ng/mL =
    # 0.00640 mg/L). Sakurai 2007 reports omegas as the standard deviation of
    # the log-normal eta (Equations 1-3: "variance of omega^2"); variances for
    # nlmixr2 ini() are therefore the squared values.

    # Structural parameters (per-kg reference values)
    lvc <- log(0.110)  ; label("Central volume per kg of body weight V1/WT (L/kg)")                             # Sakurai 2007 Table 1 theta1 = 0.110 L/kg (95% CI 0.104-0.116)
    lcl <- log(0.179)  ; label("Systemic clearance per kg in homoEM reference CL_homoEM/WT (L/h/kg)")           # Sakurai 2007 Table 1 theta2 = 0.179 L/h/kg (95% CI 0.163-0.195)
    lvp <- log(0.201)  ; label("Peripheral volume per kg of body weight V2/WT (L/kg)")                          # Sakurai 2007 Table 1 theta5 = 0.201 L/kg (95% CI 0.179-0.223)
    lq  <- log(0.0882) ; label("Inter-compartmental clearance per kg Q/WT (L/h/kg)")                            # Sakurai 2007 Table 1 theta6 = 0.0882 L/h/kg (95% CI 0.0838-0.0926)

    # CYP2C19 covariate effects on CL (multiplicative factors; homoEM = 1.0)
    e_cyp2c19_im_cl <- 0.612 ; label("CYP2C19 IM (heteroEM) multiplicative factor on CL (unitless)")            # Sakurai 2007 Table 1 theta3 = 0.612 (heteroEM/homoEM CL ratio; 95% CI 0.548-0.676)
    e_cyp2c19_pm_cl <- 0.212 ; label("CYP2C19 PM multiplicative factor on CL (unitless)")                       # Sakurai 2007 Table 1 theta4 = 0.212 (PM/homoEM CL ratio; 95% CI 0.178-0.246)

    # IIV (log-normal); variance = (Table 1 omega)^2
    # omega_V1 = 0.149 (95% CI 0.126-0.169) -> variance 0.022201
    etalvc ~ 0.022201
    # omega_CL = 0.221 (95% CI 0.191-0.247) -> variance 0.048841
    etalcl ~ 0.048841
    # omega_V2 = 0.0880 (95% CI 0.0525-0.113) -> variance 0.007744
    etalvp ~ 0.007744
    # Q has no IIV (Sakurai 2007 Equation 4: Qi = theta6 * WT, no eta_Q term)

    # Residual error -- combined proportional + additive (Sakurai 2007 Equation 5:
    # C_ij = C_pred,ij * (1 + eps_CV,ij) + eps_ADD,ij, with variances sigma^2_CV
    # and sigma^2_ADD reported as SDs in Table 1).
    propSd <- 0.0604  ; label("Proportional residual SD (unitless fraction)")  # Sakurai 2007 Table 1 sigma_CV = 0.0604 (95% CI 0.0521-0.0677)
    addSd  <- 0.00640 ; label("Additive residual SD (mg/L)")                   # Sakurai 2007 Table 1 sigma_ADD = 6.40 ng/mL = 0.00640 mg/L (95% CI 0.00-10.2 ng/mL)
  })

  model({
    # CYP2C19 phenotype multiplicative factor on CL.
    # Sakurai 2007 Equation 2: CL_i = theta2 * theta3 * theta4 * WT * exp(eta_CL),
    # with theta3 fixed to 1 for non-heteroEM subjects and theta4 fixed to 1 for
    # non-PM subjects. With the canonical binary indicators (CYP2C19_IM = 1 for
    # heteroEM; CYP2C19_PM = 1 for PM; both = 0 for homoEM reference), the
    # power-of-binary-indicator form theta^IND is 1 when IND = 0 and theta when
    # IND = 1, reproducing the paper's per-subject parameter selection.
    f_cyp2c19 <- (e_cyp2c19_im_cl ^ CYP2C19_IM) *
                 (e_cyp2c19_pm_cl ^ CYP2C19_PM)

    # Individual parameters (per-kg scaled by WT; CYP2C19 effect on CL only).
    vc <- exp(lvc + etalvc) * WT
    cl <- exp(lcl + etalcl) * f_cyp2c19 * WT
    vp <- exp(lvp + etalvp) * WT
    q  <- exp(lq) * WT

    # Two-compartment intravenous PK (no depot; dosing is into central).
    Cc <- linCmt()
    Cc ~ add(addSd) + prop(propSd)
  })
}
