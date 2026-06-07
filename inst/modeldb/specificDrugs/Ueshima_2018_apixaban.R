Ueshima_2018_apixaban <- function() {
  description <- paste0(
    "One-compartment population pharmacokinetic and pharmacogenomic model ",
    "for oral apixaban in Japanese adult patients with atrial fibrillation ",
    "(Ueshima 2018). Apparent oral clearance CL/F is the sum of an apparent ",
    "renal arm (power on creatinine clearance, CCR/70) and an apparent ",
    "non-renal arm carrying two recessive-/dominant-style pharmacogenomic ",
    "factors: CYP3A5 *3 carrier (genotype *1/*3 or *3/*3) reduces non-renal ",
    "CL/F by a factor of 0.312, and ABCG2 421A/A (rs2231142 homozygous ",
    "variant) reduces non-renal CL/F by a factor of 0.341. Apparent volume ",
    "of distribution Vd/F = 24.7 L (no significant covariates). Absorption ",
    "rate constant ka was fixed at 0.42 1/h from a prior publication ",
    "(Frost 2013 Br J Clin Pharmacol, reference 13 in the paper) because ",
    "the sparse trough-and-2-point-postdose sampling design lacked enough ",
    "absorption-phase data to identify ka."
  )
  reference <- paste0(
    "Ueshima S, Hira D, Kimura Y, Fujii R, Tomitsuka C, Yamane T, ",
    "Tabuchi Y, Ozawa T, Itoh H, Ohno S, Horie M, Terada T, Katsura T. ",
    "Population pharmacokinetics and pharmacogenomics of apixaban in ",
    "Japanese adult patients with atrial fibrillation. ",
    "Br J Clin Pharmacol. 2018;84(6):1301-1312. doi:10.1111/bcp.13561."
  )
  vignette <- "Ueshima_2018_apixaban"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance (Cockcroft-Gault)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Time-fixed (baseline) per subject. Ueshima 2018 Methods: 'Ccr was ",
        "calculated using the Cockcroft-Gault equation' (raw mL/min; not ",
        "BSA-normalised). Cohort range 30.6-145.5 mL/min, median 69.8 ",
        "mL/min (Table 1). Reference value 70 mL/min in the CL/F equation ",
        "(CCR/70)^theta4; the paper rounds the cohort median 69.8 to 70 ",
        "for the centring constant. The raw Cockcroft-Gault scale follows ",
        "the existing `CLCR` raw-CrCl precedent in the register; the ",
        "canonical CRCL description notes that BSA-normalised vs raw is ",
        "paper-dependent and recorded in per-model notes."
      ),
      source_name        = "Ccr"
    ),
    CYP3A5_STAR1_HOM = list(
      description        = "CYP3A5*1/*1 homozygote indicator (two functional CYP3A5*1 alleles at rs776746)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5*1/*3 heterozygote or CYP3A5*3/*3 nonexpresser; i.e., NOT CYP3A5*1/*1)",
      notes              = paste0(
        "Time-fixed per subject (germline genotype). 1 = subject is ",
        "CYP3A5*1/*1 homozygote; 0 = otherwise (the union of *3/*3 ",
        "nonexpressers and *1/*3 heterozygotes). Ueshima 2018 cohort: ",
        "4/81 (4.9%) *1/*1, 30/81 (37.0%) *1/*3, 47/81 (58.1%) *3/*3 ",
        "(Table 2). The source paper encodes its CYP3A5 indicator with ",
        "the OPPOSITE orientation: paper CYP3A5_ind = 1 if subject ",
        "carries any *3 allele (*1/*3 or *3/*3) and 0 if *1/*1, i.e., ",
        "paper CYP3A5_ind = 1 - CYP3A5_STAR1_HOM. The model() block ",
        "applies the effect to (1 - CYP3A5_STAR1_HOM) so the canonical ",
        "0/1 encoding preserves the paper's reported multiplicative ",
        "factor 0.312 unchanged. Paired indicator CYP3A5_STAR1_HET is ",
        "NOT used because Ueshima 2018 pools *1/*3 with *3/*3 (recessive-",
        "expresser model) rather than estimating a distinct ",
        "*1/*3-vs-*3/*3 effect (the cohort had only n = 4 *1/*1 subjects ",
        "and likely could not have identified a three-level effect)."
      ),
      source_name        = "CYP3A5 *3 carrier indicator (paper Methods Eq.: CYP3A5 = 1 if *1/*3 or *3/*3, 0 if *1/*1)"
    ),
    SNP_ABCG2_RS2231142_HOM = list(
      description        = "ABCG2 rs2231142 (Q141K / c.421C>A) homozygous-variant indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (homozygous wild-type 421C/C or heterozygous 421C/A; i.e., NOT 421A/A)",
      notes              = paste0(
        "Time-fixed per subject (germline genotype). 1 = subject is ",
        "homozygous variant 421A/A; 0 = otherwise (the union of 421C/C ",
        "homozygous wild-type and 421C/A heterozygous variant carriers). ",
        "Ueshima 2018 cohort: 9/81 (11.1%) 421A/A, 33/81 (40.7%) 421C/A, ",
        "39/81 (48.2%) 421C/C (Table 2). Recessive-model encoding: ",
        "heterozygotes are pooled with wild-type homozygotes because ",
        "Ueshima 2018 reported that only the 421A/A stratum had a ",
        "distinct typical-value effect on apixaban non-renal clearance ",
        "(paper Methods Eq.: ABCG2 = 1 if 421A/A, 0 if 421C/C or 421C/A; ",
        "Table 4 theta6 = 0.341 multiplier on the non-renal arm)."
      ),
      source_name        = "ABCG2 421A/A genotype indicator (paper Methods Eq.: ABCG2 = 1 if 421A/A, 0 if 421C/C or 421C/A)"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 81L,
    n_studies        = 1L,
    n_observations   = 276L,
    age_range        = "40.5-84.9 years (median 68.1)",
    age_median       = "68.1 years",
    weight_range     = "41.0-92.2 kg (median 65.0)",
    weight_median    = "65.0 kg",
    sex_female_pct   = 24.7,
    race_ethnicity   = c(Asian = 100),
    disease_state    = paste0(
      "Japanese adult inpatients and outpatients with non-valvular ",
      "atrial fibrillation receiving oral apixaban at Shiga University ",
      "of Medical Science Hospital from February 2015 to May 2016. ",
      "Co-morbidities: hypertension (44%), diabetes mellitus (42%), ",
      "heart failure (25%), dyslipidemia (22%). Concomitant CYP3A4 / ",
      "P-glycoprotein inhibitors (amiodarone, n = 17; verapamil, n = 6) ",
      "and inducers (rifampicin, n = 1 -- excluded from analysis) were ",
      "permitted but did not survive covariate selection. Patients with ",
      "poor drug compliance were excluded."
    ),
    dose_range       = paste0(
      "Oral apixaban (Eliquis, Bristol-Myers Squibb / Pfizer) 5, 10, or ",
      "20 mg/day total, given as twice-daily 2.5, 5, or 10 mg tablets. ",
      "Inpatients sampled trough + 0.5-2 h + 9-12 h post-dose; ",
      "outpatients sampled a single timepoint 0.3-16 h post-dose per ",
      "hospital visit. Apixaban plasma concentrations measured by ",
      "LC/MS/MS (lower limit of quantification 2.5 ng/mL; linear range ",
      "2.5-500 ng/mL)."
    ),
    regions          = "Japan (single-centre observational study at Shiga University of Medical Science Hospital)",
    creatinine_clearance = "30.6-145.5 mL/min, median 69.8 (Cockcroft-Gault; raw mL/min, not BSA-normalised)",
    serum_creatinine = "0.41-1.34 mg/dL, median 0.88 (Table 1)",
    ast_range        = "13-97 IU/L, median 23",
    alt_range        = "5-115 IU/L, median 19",
    cyp3a5_genotype  = c(`*1/*1` = 4.9, `*1/*3` = 37.0, `*3/*3` = 58.1),
    abcg2_421_genotype = c(`421C/C` = 48.2, `421C/A` = 40.7, `421A/A` = 11.1),
    notes            = paste0(
      "Software: NONMEM 7.3.0 (FOCE-I); Perl-speaks-NONMEM v4.7.0 for ",
      "bootstrap. Baseline demographics in Table 1; final-model parameter ",
      "estimates in Table 4. Covariates tested and NOT retained: AST, ALT, ",
      "body weight, age, ABCB1 1236C>T / 2677G>T/A / 3435C>T genotypes, ",
      "co-morbidities (hypertension, diabetes, heart failure, ",
      "dyslipidemia), and concomitant CYP3A4 / P-glycoprotein ",
      "inhibitors (amiodarone, verapamil). No covariates were retained ",
      "on Vd/F."
    )
  )

  ini({
    # Structural parameters - Ueshima 2018 Table 4 final estimates.

    # Apparent oral absorption rate constant ka was FIXED at the literature
    # value of 0.42 1/h from a prior healthy-subject popPK analysis (Frost
    # 2013, reference 13 in the paper) because Ueshima 2018 had sparse
    # trough + 2-point post-dose sampling and could not identify ka from
    # the absorption phase (Methods: 'ka was fixed at the reported value
    # of 0.42 h^-1 because of a lack of data on the absorption phase').
    lka <- fixed(log(0.42)); label("Absorption rate constant ka (1/h); fixed at literature value")  # Ueshima 2018 Table 4 theta3 = 0.42 h^-1 (fixed)

    # Base apparent oral clearance theta1 = 1.53 L/h. This is the typical
    # value of CL/F for a reference patient with CCR = 70 mL/min,
    # CYP3A5*1/*1, and ABCG2 421C/C or 421C/A. In the paper's final
    # equation the same theta1 multiplies both the renal arm (CCR/70)^theta4
    # and the non-renal arm theta5^CYP3A5 * theta6^ABCG2, so at the
    # reference covariate values both arms equal theta1 = 1.53 L/h and the
    # total reference CL/F is 2 * theta1 = 3.06 L/h (paper Results: 'The
    # population mean of CL/F for a typical patient (Ccr value of 70 ml
    # min-1) with the CYP3A5 *1/*1 and ABCG2 421C/C or C/A genotypes was
    # estimated to be 3.06 l h-1').
    lcl <- log(1.53); label("Per-arm apparent oral clearance theta1 (L/h) at the reference patient (CCR=70, CYP3A5*1/*1, ABCG2 421C/C or 421C/A)")  # Ueshima 2018 Table 4 theta1 = 1.53 L/h

    # Apparent volume of distribution theta2 = 24.7 L. No covariates
    # retained on Vd/F (Ueshima 2018 Results: 'no covariates affected the
    # population mean of the apparent volume of distribution').
    lvc <- log(24.7); label("Apparent volume of distribution Vd/F (L)")  # Ueshima 2018 Table 4 theta2 = 24.7 L

    # Power exponent on (CCR / 70) for the renal arm CLR/F = theta1 *
    # (CCR/70)^theta4. theta4 = 0.700 (95% CI 0.471-0.929; bootstrap
    # median 0.714).
    e_crcl_cl_renal <- 0.700; label("Power exponent on (CRCL / 70) for the renal arm of CL/F (unitless)")  # Ueshima 2018 Table 4 theta4 = 0.700

    # Multiplicative factor on the non-renal arm for CYP3A5 *3 carriers
    # (genotype *1/*3 or *3/*3 -- i.e., NOT CYP3A5*1/*1). Paper Table 4
    # theta5 = 0.312 (95% CI 0.273-0.351; bootstrap median 0.342). Encoded
    # as the log of the multiplicative factor so the model() block can
    # apply it via exp(e_cyp3a5_star1_hom_cl_nonren * (1 -
    # CYP3A5_STAR1_HOM)): when the subject is *1/*1 (CYP3A5_STAR1_HOM = 1)
    # the exponent (1 - 1) = 0 collapses the factor to 1, and when the
    # subject is a *3 carrier (CYP3A5_STAR1_HOM = 0) the exponent is 1 and
    # the factor reduces to exp(log(0.312)) = 0.312.
    e_cyp3a5_star1_hom_cl_nonren <- log(0.312); label("Log of CL/F non-renal arm multiplicative factor for CYP3A5 *3 carriers (genotype *1/*3 or *3/*3) relative to *1/*1 reference")  # Ueshima 2018 Table 4 theta5 = 0.312

    # Multiplicative factor on the non-renal arm for ABCG2 421A/A
    # homozygotes (rs2231142 homozygous variant). Paper Table 4 theta6 =
    # 0.341 (95% CI 0.160-0.522; bootstrap median 0.478). Encoded as the
    # log of the multiplicative factor; in model() applied via
    # exp(e_snp_abcg2_rs2231142_hom_cl_nonren * SNP_ABCG2_RS2231142_HOM):
    # when SNP_ABCG2_RS2231142_HOM = 1 (421A/A) the factor reduces to
    # exp(log(0.341)) = 0.341; otherwise the factor is 1.
    e_snp_abcg2_rs2231142_hom_cl_nonren <- log(0.341); label("Log of CL/F non-renal arm multiplicative factor for ABCG2 421A/A homozygous variant relative to 421C/C or 421C/A reference")  # Ueshima 2018 Table 4 theta6 = 0.341

    # Inter-individual variability - Ueshima 2018 Table 4.
    # IIV is expressed in the paper as %CV (square root of omega^2 on the
    # log scale, conventional NONMEM reporting): eta_1 = 26.6% on CL/F
    # (variance = 0.266^2 = 0.0708), eta_2 = 56.6% on Vd/F (variance =
    # 0.566^2 = 0.3204). The off-diagonal covariance (Methods: 'The
    # covariance between inter-individual variability for CL/F and that
    # for Vd/F was 11.6%') corresponds to 0.116 on the omega-matrix scale;
    # the implied correlation 0.116 / (0.266 * 0.566) = 0.770 matches the
    # paper's reported correlation coefficient between individual CL/F and
    # Vd/F (Results: 'the correlation coefficient between individual CL/F
    # and Vd/F was 0.770').
    etalcl + etalvc ~ c(0.0708,
                        0.116, 0.3204)  # Ueshima 2018 Table 4: Var(CL/F) = 0.266^2 = 0.0708; Cov(CL/F, Vd/F) = 0.116; Var(Vd/F) = 0.566^2 = 0.3204

    # Residual unexplained variability - Ueshima 2018 Methods Eq.: an
    # exponential error model C_ij_obs = C_ij_pred * exp(epsilon), with
    # epsilon ~ N(0, sigma^2). Table 4 reports the residual as 34.0%
    # (RSE 12.0%; 95% CI 28.0-40.0). The 34.0% is the standard deviation
    # of epsilon on the log scale, encoded directly as the log-normal
    # error magnitude expSd = 0.34.
    expSd <- 0.340; label("Log-normal residual error SD (unitless; SD of log(Cobs / Cpred))")  # Ueshima 2018 Table 4 epsilon = 34.0% CV
  })

  model({
    # 1. Derived covariate indicators.
    # The paper's CYP3A5 indicator is 1 if the subject carries any *3
    # allele (*1/*3 or *3/*3) and 0 if *1/*1. The canonical column
    # CYP3A5_STAR1_HOM uses the OPPOSITE orientation (1 if and only if
    # *1/*1); apply the inversion here so the encoded effect parameter
    # (log of the source-paper multiplicative factor 0.312) is used
    # unchanged.
    cyp3a5_star3_carrier <- 1 - CYP3A5_STAR1_HOM

    # 2. Typical-value CL/F arms.
    # Renal arm: theta1 * (CCR / 70)^theta4. Approaches zero for CRCL ->
    # 0 (anuric patient) which underestimates true CL/F because the paper
    # does not include the (small) renal clearance of unbound apixaban via
    # alternative pathways; the simulation should be restricted to the
    # CCR range of the fitting cohort (30.6-145.5 mL/min).
    cl_renal_typ <- exp(lcl) * (CRCL / 70)^e_crcl_cl_renal

    # Non-renal arm: theta1 * theta5^CYP3A5 * theta6^ABCG2. The factor is
    # the product of two power-of-binary-indicator multipliers,
    # implemented as exp(log(factor) * indicator) so the same effect
    # parameter can be applied both in the reference and in the perturbed
    # state.
    cl_nonren_typ <- exp(lcl) *
      exp(e_cyp3a5_star1_hom_cl_nonren * cyp3a5_star3_carrier) *
      exp(e_snp_abcg2_rs2231142_hom_cl_nonren * SNP_ABCG2_RS2231142_HOM)

    # 3. Individual CL/F and Vd/F with correlated exponential IIV.
    # IIV is applied multiplicatively to the TOTAL CL/F (not to each arm
    # separately), per the standard NONMEM exponential-IIV convention.
    cl <- (cl_renal_typ + cl_nonren_typ) * exp(etalcl)
    vc <- exp(lvc + etalvc)
    ka <- exp(lka)

    # 4. One-compartment ODE system with first-order oral absorption.
    kel <- cl / vc
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 5. Observation and error.
    # Convert mg/L (central / vc with dose in mg, vc in L) to ng/mL by
    # multiplying by 1000 (1 mg/L = 1 ug/mL = 1000 ng/mL). Exponential /
    # log-normal residual error per the paper's Methods.
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd)
  })
}
