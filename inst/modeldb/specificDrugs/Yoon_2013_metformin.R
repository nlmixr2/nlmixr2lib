Yoon_2013_metformin <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order oral absorption and an absorption-lag time for a single 500 mg oral dose of metformin in healthy Korean male adults (Yoon 2013). Body weight enters V/F as a linear-in-deviation covariate; two organic-cation-transporter (OCT) genetic polymorphisms enter CL/F as multiplicative fractional shifts, OCT2 c.808G>T (SLC22A2, A270S, widely known as rs316019) and OCTN1 c.917C>T (SLC22A4, T306I). Variant carriers are pooled (heterozygotes and homozygous variants together vs the homozygous wild-type reference) per the source paper's dominant-model encoding. Slow first-order absorption (ka = 0.248 1/h) relative to elimination (kel = cl/vc = 1.21 1/h) produces the apparent flip-flop kinetics reported for metformin."
  reference   <- "Yoon H, Cho H-Y, Yoo H-D, Kim S-M, Lee Y-B. Influences of organic cation transporter polymorphisms on the population pharmacokinetics of metformin in healthy subjects. AAPS J. 2013;15(2):571-580. doi:10.1208/s12248-013-9460-z"
  vignette    <- "Yoon_2013_metformin"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight at enrolment as reported by Yoon 2013 Methods 'Subjects and Study Design' (mean +/- SD 67.74 +/- 8.24 kg; range 53.1-95.6 kg across 96 healthy Korean males).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-in-deviation covariate on V/F: V/F = 112 * (1 + theta_WT * (WT - 70)). The paper does not explicitly print the reference weight; 70 kg is used here as the rounded-standard convention (see the vignette Assumptions and deviations section for the sensitivity of predicted V/F to alternative choices of 67.74 kg cohort mean vs 70 kg round). Time-fixed per subject.",
      source_name        = "WT"
    ),
    SNP_SLC22A2_808GT = list(
      description        = "Binary genotype indicator for the SLC22A2 (OCT2) c.808G>T single-nucleotide polymorphism (protein change A270S; widely reported in the pharmacogenomic literature as rs316019). 1 = at least one variant (T) allele present (heterozygous 808GT or homozygous 808TT carrier); 0 = homozygous wild-type 808GG.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (homozygous wild-type 808GG)",
      notes              = "Time-fixed per subject (germline genotype). Yoon 2013 pooled heterozygous and homozygous variant carriers into a single variant group (dominant-model encoding) because only heterozygous 808GT carriers were observed in the Korean cohort (n = 15 GT, 0 TT out of 96; Table I). Enters CL/F multiplicatively as a fractional shift: CL/F = 136 * (1 - 0.248 * SNP_SLC22A2_808GT) * (1 - 0.234 * SNP_SLC22A4_917CT). Variant carriers have 24.8% lower apparent oral clearance than wild-type homozygotes at the reference OCTN1 status (Yoon 2013 Table III).",
      source_name        = "OCT2_808GT"
    ),
    SNP_SLC22A4_917CT = list(
      description        = "Binary genotype indicator for the SLC22A4 (OCTN1) c.917C>T single-nucleotide polymorphism (protein change T306I). 1 = at least one variant (T) allele present (heterozygous 917CT or homozygous 917TT carrier); 0 = homozygous wild-type 917CC.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (homozygous wild-type 917CC)",
      notes              = "Time-fixed per subject (germline genotype). Yoon 2013 pooled heterozygous and homozygous variant carriers into a single variant group (dominant-model encoding); Table I reports 52 CT plus 34 TT out of 96 in the Korean cohort (carrier rate 86 of 96 = 90%). Enters CL/F multiplicatively as a fractional shift alongside SNP_SLC22A2_808GT: CL/F = 136 * (1 - 0.248 * SNP_SLC22A2_808GT) * (1 - 0.234 * SNP_SLC22A4_917CT). Variant carriers have 23.4% lower apparent oral clearance than wild-type homozygotes at the reference OCT2 status (Yoon 2013 Table III).",
      source_name        = "OCTN1_917CT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 96L,
    n_studies      = 5L,
    age_range      = "19-31 years",
    age_median     = "22.41 years (mean; SD 2.43)",
    weight_range   = "53.1-95.6 kg",
    weight_median  = "67.74 kg (mean; SD 8.24)",
    sex_female_pct = 0,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Healthy",
    dose_range     = "500 mg single oral dose (Glucophage tablet 500 mg; Boehringer Ingelheim Korea)",
    regions        = "Korea (Gwangju; Chonnam National University)",
    pharmacogenomics = "Germline genotyping for OCT1 (289C>A, 350C>T, 616C>T, 1022C>T), OCT2 (596C>T, 602C>T, 808G>T), OCT3 (1199C>T, 1267G>T), OCTN1 (917C>T, 1507C>T), MATE1 (191G>A, 373C>T, 929C>T, 983A>C, 1421A>G), and MATE2-K (130G>A, 192G>T, 632_633GC>TT) by direct sequencing plus SNaPshot Multiplex assay (Yoon 2013 Methods 'Genotype Analysis'; allele frequencies in Table I). Only OCT2 808G>T and OCTN1 917C>T were retained as significant covariates on CL/F in the final population model (Table III).",
    notes          = "Ninety-six healthy Korean male volunteers pooled from five independent single-dose bioequivalence studies conducted at the Institute of Bioequivalence and Bridging Study, Chonnam National University; only the reference-formulation arm of each BE study was used. Serum sampling schedule: pre-dose and 0.5, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 10, 12, and 24 h post-dose (13 samples per subject). Bioanalysis by validated HPLC-UV (LLOQ 10 ng/mL; linear range 10-2000 ng/mL). No dropouts across the five studies. Baseline demographics from Yoon 2013 Methods 'Subjects and Study Design' (age 19-31 y, mean 22.41 +/- 2.43 y; weight 53.1-95.6 kg, mean 67.74 +/- 8.24 kg)."
  )

  ini({
    # Structural PK parameters (Yoon 2013 Table III final-model point
    # estimates; %RSE shown in the source; typical-value reference is a
    # 70-kg wild-type-both subject at OCT2-808GG and OCTN1-917CC).
    lka   <- log(0.248); label("First-order oral absorption rate constant ka (1/h)")               # Yoon 2013 Table III TV(ka) = 0.248 1/h (RSE 3.00%)
    lcl   <- log(136);   label("Apparent oral clearance CL/F at the OCT2-808GG + OCTN1-917CC reference (L/h)") # Yoon 2013 Table III TV(CL/F) = 136 L/h (RSE 18.4%)
    lvc   <- log(112);   label("Apparent central volume of distribution V/F at the 70 kg reference (L)")       # Yoon 2013 Table III TV(V/F) = 112 L (RSE 6.88%)
    ltlag <- log(0.182); label("Absorption lag time (h)")                                                       # Yoon 2013 Table III TV(Tlag) = 0.182 h (RSE 9.78%)

    # Covariate effects. All three are estimated (not fixed).
    #   CL/F = TV(CL/F) * (1 + theta_OCT2 * SNP_SLC22A2_808GT)
    #                  * (1 + theta_OCTN1 * SNP_SLC22A4_917CT) * exp(etalcl)
    # with theta_OCT2 = -0.248 and theta_OCTN1 = -0.234 (negative signs
    # because variant carriers have LOWER CL/F than wild-type).
    #   V/F = TV(V/F) * (1 + theta_WT * (WT - 70)) * exp(etalvc)
    # with theta_WT = 0.0183.
    e_snp_oct2_cl  <- -0.248;  label("Fractional shift in CL/F for OCT2-808G>T variant carriers (unitless)")  # Yoon 2013 Table III theta_OCT2 = -0.248 (RSE 15.8%)
    e_snp_octn1_cl <- -0.234;  label("Fractional shift in CL/F for OCTN1-917C>T variant carriers (unitless)") # Yoon 2013 Table III theta_OCTN1 = -0.234 (RSE 61.1%)
    e_wt_vc        <-  0.0183; label("Linear slope of V/F on (WT - 70 kg) (1/kg)")                             # Yoon 2013 Table III theta_BW = 0.0183 (RSE 42.8%)

    # IIV - exponential log-normal between-subject variability on all four
    # structural parameters (Yoon 2013 Methods 'Population Pharmacokinetic
    # Analysis'; variances w^2 reported directly on the log-scale in
    # Table III, so ini() uses the reported w^2 without transformation).
    # Yoon 2013 Table III does not report off-diagonal covariances; all
    # etas are treated as independent.
    etalcl   ~ 0.0800   # Yoon 2013 Table III w^2 CL/F = 0.0800 (RSE 13.4%; shrinkage 1.20%)
    etalvc   ~ 0.295    # Yoon 2013 Table III w^2 V/F  = 0.295  (RSE 15.4%; shrinkage 11.7%)
    etalka   ~ 0.0572   # Yoon 2013 Table III w^2 ka   = 0.0572 (RSE 9.63%; shrinkage 12.0%)
    etaltlag ~ 0.281    # Yoon 2013 Table III w^2 lag  = 0.281  (RSE 24.6%; shrinkage 24.8%)

    # Residual error - proportional error model (Yoon 2013 Methods; final
    # residual variance sigma^2_pro = 0.0291 corresponding to a proportional
    # SD of sqrt(0.0291) = 0.1706 ~= 17%).
    propSd <- 0.1706; label("Proportional residual SD on Cc (fraction)")  # Yoon 2013 Table III sigma^2_pro = 0.0291 (RSE 8.25%)
  })

  model({
    # 1. Reference values for the linear-in-deviation covariate on V/F.
    ref_wt <- 70

    # 2. Individual PK parameters. The two SNP effects on CL/F combine
    # multiplicatively as fractional shifts (Yoon 2013 Results 'Population
    # Pharmacokinetic Model'; final equation reported symbolically in the
    # paragraph following 'Model 7 (Supplementary Table S2)'). The BW
    # effect on V/F is linear-in-deviation centred on 70 kg (see the
    # vignette Assumptions and deviations for the reference-weight choice).
    ka   <- exp(lka   + etalka)
    cl   <- exp(lcl   + etalcl)   * (1 + e_snp_oct2_cl  * SNP_SLC22A2_808GT) *
                                    (1 + e_snp_octn1_cl * SNP_SLC22A4_917CT)
    vc   <- exp(lvc   + etalvc)   * (1 + e_wt_vc * (WT - ref_wt))
    tlag <- exp(ltlag + etaltlag)

    # 3. Micro-constant
    kel <- cl / vc

    # 4. ODE system (one-compartment, first-order absorption with lag)
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 5. Absorption lag applied to the depot compartment
    alag(depot) <- tlag

    # 6. Observation and residual error. Dose in mg, V/F in L give
    # central/vc in mg/L, matching Yoon 2013's reported concentration
    # units (Table II C_max in mg/L).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
