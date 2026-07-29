Mukonzo_2009_efavirenz <- function() {
  description <- "Two-compartment population PK model for single-dose oral efavirenz in 121 healthy Ugandan adults, with sequential zero-order followed by first-order absorption to the central compartment. Apparent oral clearance CL/F is reduced by 21% in homozygous CYP2B6*6 (rs3745274 T/T) and by 20% in homozygous CYP2B6*11 (rs35303484 G/G) carriers (multiplicative fractional effects). Relative bioavailability Frel is increased by 26% in ABCB1 rs3842 mutant carriers (heterozygote or homozygote). Apparent peripheral volume Vp/F is 2.08-fold higher in women than in men. Concentrations are reported in mg/L (1 mg/L efavirenz = 3.168 micromol/L)."
  reference <- paste(
    "Mukonzo JK, Roshammar D, Waako P, Andersson M, Fukasawa T, Milani L,",
    "Svensson JO, Ogwal-Okeng J, Gustafsson LL, Aklillu E.",
    "A novel polymorphism in ABCB1 gene, CYP2B6*6 and sex predict",
    "single-dose efavirenz population pharmacokinetics in Ugandans.",
    "Br J Clin Pharmacol. 2009;68(5):690-699.",
    "doi:10.1111/j.1365-2125.2009.03516.x."
  )
  vignette <- "Mukonzo_2009_efavirenz"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Time-fixed per subject. Mukonzo 2009 Table 3 reports a multiplicative factor of 2.08 (95% CI 1.64, 2.52) on the apparent peripheral volume Vp/F for females relative to males. Cohort distribution (Results paragraph 1): 121 healthy adults of whom 57% (n = 69) female.",
      source_name        = "sex"
    ),
    SNP_CYP2B6_RS3745274_T_COUNT = list(
      description        = "Count of CYP2B6 c.516G>T (rs3745274, p.Q172H, CYP2B6*6) T-alleles per subject (0 / 1 / 2). 0 = GG homozygous wild-type, 1 = GT heterozygous, 2 = TT homozygous variant.",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (germline genotype). Mukonzo 2009 Table 1 reports CYP2B6 c.516G>T (rs3745274) and c.785A>G (rs2279343) as in complete linkage disequilibrium in the Ugandan cohort, jointly defining the CYP2B6*6 haplotype. Either SNP alone deterministically identifies *6 status; the canonical column uses rs3745274 (T-allele count) to match existing nlmixr2lib precedent (Dhoro 2015, Schipani 2011, Olagunju 2018). Mukonzo 2009 Table 3 reports a -20.9% multiplicative shift on CL/F for homozygous mutant (T/T, count = 2) only; heterozygotes are pooled with wild-type and receive no shift (Results: 'Homozygous CYP2B6*6 (G516T, A785G) ... displayed 21 ... percent lower apparent oral clearance'). Cohort allele frequencies (Table 1, n = 121): 516G>T 35.6%, 785A>G 36.4%.",
      source_name        = "CYP2B6 516G>T (rs3745274) and 785A>G (rs2279343); CYP2B6*6 haplotype"
    ),
    SNP_CYP2B6_RS35303484_G_COUNT = list(
      description        = "Count of CYP2B6 c.136A>G (rs35303484, p.M46V, CYP2B6*11) G-alleles per subject (0 / 1 / 2). 0 = AA homozygous wild-type, 1 = AG heterozygous, 2 = GG homozygous variant.",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (germline genotype). Mukonzo 2009 Table 1 reports CYP2B6 c.136A>G (rs35303484) as defining the CYP2B6*11 phenotypic-null allele. Table 3 reports a -19.9% multiplicative shift on CL/F for homozygous mutant (G/G, count = 2) only; heterozygotes are pooled with wild-type and receive no shift (Results: 'CYP2B6*11 polymorphism resulted in a 20% lower efavirenz clearance following single-dose administration'). Cohort allele frequency (Table 1, n = 121): 136A>G 13.6%.",
      source_name        = "CYP2B6 136A>G (rs35303484); CYP2B6*11"
    ),
    SNP_ABCB1_RS3842 = list(
      description        = "ABCB1 c.4036A>G (rs3842) mutant-allele carrier indicator (1 = G allele present, heterozygous or homozygous; 0 = AA wild-type)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (AA homozygous wild-type)",
      notes              = "Time-fixed (germline genotype). Mukonzo 2009 Table 1 reports rs3842 (c.4036A>G) as a previously uncharacterised polymorphism in the ABCB1 3' UTR. Table 3 reports a +25.7% multiplicative shift on relative bioavailability Frel for mutant carriers regardless of zygosity (Results: 'Mutant homozygote and heterozygote individuals for ABCB1 rs3842 exhibited 26% greater efavirenz bioavailability than wild-type carriers'). Encoded as binary (any G allele = 1) to match the carrier-grouping in the source paper. Cohort allele frequency (Table 1, n = 121): 4036A>G 16.8%.",
      source_name        = "ABCB1 4036A>G (rs3842)"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 121L,
    n_studies        = 1L,
    n_observations   = 402L,
    age_range        = "mean 26.5 years (SD 8.2)",
    weight_range     = "mean 57.5 kg (SD 5.9)",
    sex_female_pct   = 57,
    race_ethnicity   = c(Ugandan = 100),
    disease_state    = "healthy adult volunteers",
    dose_range       = "600 mg single oral dose efavirenz (Stocrin)",
    regions          = "Uganda (Kampala)",
    sampling_window  = "intensive (n = 32): 0, 1, 2, 4, 8, 24, 48, and 72 h post-dose; sparse (n = 89): 4 and 24 h post-dose",
    biochem_baseline = "mean serum albumin 41.0 g/L (SD 8.9); ALT 10.8 U/L (SD 9.7); urea 4.14 mmol/L (SD 9.0); creatinine 108.4 umol/L (SD 37.4) -- Results paragraph 1",
    cyp2b6_freq      = "CYP2B6*6 (rs3745274/rs2279343, complete LD): 516G>T 35.6%, 785A>G 36.4%. CYP2B6*11 (rs35303484): 136A>G 13.6% (Mukonzo 2009 Table 1).",
    abcb1_freq       = "ABCB1 rs3842 (4036A>G): 16.8% (Mukonzo 2009 Table 1).",
    notes            = "Adult Ugandan healthy volunteers receiving a single 600 mg oral efavirenz dose; plasma concentrations measured by HPLC-UV (LLOQ 0.35 micromol/L). NONMEM VI FOCE-INTER estimation with stepwise univariate covariate inclusion (P < 0.05) and backward elimination (P < 0.01). HIV / hepatitis B serology and liver / renal function tests confirmed health; subjects were drug-free for one week pre-dose."
  )

  ini({
    # ---- Structural fixed-effect estimates (Mukonzo 2009 Table 3 final model) ----
    # Two-compartment apparent-oral disposition; zero-order input over duration D
    # followed by first-order absorption rate ka from the depot to the central
    # compartment ("two-compartment pharmacokinetic model with zero-order input
    # to the dose compartment followed by sequential first-order absorption to
    # the central compartment" -- Mukonzo 2009 Results paragraph "Pharmacokinetic
    # modelling"). Typical values are reported for a wild-type extensive
    # metaboliser male (the reference category for all four covariates).
    lka <- log(0.146); label("Absorption rate constant ka (1/h)")                                       # Mukonzo 2009 Table 3: ka = 0.146 1/h (95% CI 0.0558, 0.236)
    lcl <- log(4.00) ; label("Apparent oral clearance CL/F at wild-type extensive metaboliser (L/h)")   # Mukonzo 2009 Table 3: CL/F = 4.00 L/h (95% CI 3.47, 4.53)
    lvc <- log(19.1) ; label("Apparent central volume of distribution Vc/F (L)")                        # Mukonzo 2009 Table 3: Vc/F = 19.1 L (95% CI 7.46, 30.7)
    lvp <- log(155)  ; label("Apparent peripheral volume of distribution Vp/F at male reference (L)")   # Mukonzo 2009 Table 3: Vp/F = 155 L (95% CI 131, 179); typical for males
    lq  <- log(13.7) ; label("Apparent intercompartmental clearance Q/F (L/h)")                         # Mukonzo 2009 Table 3: Q/F = 13.7 L/h (95% CI 6.1, 21.3)
    ld  <- log(1.07) ; label("Duration of zero-order input D (h)")                                      # Mukonzo 2009 Table 3: D = 1.07 h (95% CI 0.758, 1.38)

    # Relative bioavailability anchor: Mukonzo 2009 Table 3 reports Frel = "1 FIX"
    # for the wild-type ABCB1 genotype; covariate effects shift mutant carriers
    # multiplicatively from this anchor.
    lfdepot <- fixed(log(1)); label("Relative bioavailability Frel at ABCB1 rs3842 wild-type reference (fixed)")  # Mukonzo 2009 Table 3: Frel = 1 FIX

    # ---- Covariate effects ----
    # All four effects are fractional / multiplicative shifts from the reference
    # category typical value. Mukonzo 2009 simulation cross-check: a homozygous-
    # mutant subject for CYP2B6*6 + CYP2B6*11 + ABCB1 rs3842 has predicted
    # AUC = 600 * (1 + 0.257) / (4 * (1 - 0.209) * (1 - 0.199))
    #     = 600 * 1.257 / 2.535 = 297.5 mg.h/L = 942 micromol.h/L (vs reported
    # 943 micromol.h/L) confirming the multiplicative-on-CL parameterisation.
    # Wild-type cross-check: AUC = 600 / 4 = 150 mg.h/L = 475 micromol.h/L
    # (vs reported 475 micromol.h/L).
    e_2b6_6_cl       <- -0.209; label("Fractional change in CL/F for homozygous CYP2B6*6 (rs3745274 T/T) vs reference (unitless)")  # Mukonzo 2009 Table 3: Effect of CYP2B6*6 = -0.209 (95% CI -0.386, -0.032)
    e_2b6_11_cl      <- -0.199; label("Fractional change in CL/F for homozygous CYP2B6*11 (rs35303484 G/G) vs reference (unitless)")# Mukonzo 2009 Table 3: Effect of CYP2B6*11 = -0.199 (95% CI -0.329, -0.0691)
    e_rs3842_fdepot  <-  0.257; label("Fractional change in Frel for ABCB1 rs3842 mutant carrier vs AA wild-type (unitless)")       # Mukonzo 2009 Table 3: Effect of ABCB1 (rs3842) = 0.257 (95% CI 0.0873, 0.427)
    e_sexf_vp        <-  2.08 ; label("Ratio of Vp/F in females (SEXF = 1) to males (SEXF = 0); applied as Vp = Vp_male * (e_sexf_vp^SEXF)")  # Mukonzo 2009 Table 3: Effect of sex = 2.08 (95% CI 1.64, 2.52)

    # ---- IIV (exponential / log-normal on each structural parameter) ----
    # Mukonzo 2009 Table 3 reports between-subject variability as CV%. The
    # log-normal variance is omega^2 = log(CV^2 + 1).
    etalcl     ~ 0.019408   # Mukonzo 2009 Table 3: omega(CL) = 14.0% CV; omega^2 = log(0.140^2 + 1) = 0.019408
    etalvc     ~ 0.688200   # Mukonzo 2009 Table 3: omega(Vc) = 99.5% CV; omega^2 = log(0.995^2 + 1) = 0.688200
    etalvp     ~ 0.074979   # Mukonzo 2009 Table 3: omega(Vp) = 27.9% CV; omega^2 = log(0.279^2 + 1) = 0.074979
    etalq      ~ 0.098055   # Mukonzo 2009 Table 3: omega(Q)  = 32.1% CV; omega^2 = log(0.321^2 + 1) = 0.098055
    etalka     ~ 0.038077   # Mukonzo 2009 Table 3: omega(ka) = 19.7% CV; omega^2 = log(0.197^2 + 1) = 0.038077
    etald      ~ 0.395940   # Mukonzo 2009 Table 3: omega(D1) = 69.7% CV; omega^2 = log(0.697^2 + 1) = 0.395940
    etalfdepot ~ 0.034729   # Mukonzo 2009 Table 3: omega(Frel)= 18.8% CV; omega^2 = log(0.188^2 + 1) = 0.034729

    # ---- Residual error ----
    # The Methods paragraph specifies an "intercept-slope" (combined additive +
    # proportional) residual error model; the Results paragraph after Table 3
    # states "In the final model the additive part of the combined residual
    # error model was insignificantly small", so only the proportional slope
    # is encoded here. Table 3 reports sigma_prop = 13.9% CV which is read as
    # the proportional SD on the linear concentration scale.
    propSd <- 0.139; label("Proportional residual error (fraction; the additive term was insignificantly small in the final fit)")  # Mukonzo 2009 Table 3: sigma_prop = 13.9% CV (95% CI 9.62, 17.1)
  })

  model({
    # 1. Genotype indicators derived from the per-allele count columns. Only
    #    homozygous mutants receive the CYP2B6 effects (Mukonzo 2009 Results).
    is_2b6_6_hom  <- (SNP_CYP2B6_RS3745274_T_COUNT  == 2)
    is_2b6_11_hom <- (SNP_CYP2B6_RS35303484_G_COUNT == 2)

    # 2. Individual PK parameters with multiplicative covariate effects.
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) *
          (1 + e_2b6_6_cl  * is_2b6_6_hom) *
          (1 + e_2b6_11_cl * is_2b6_11_hom)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp) * (e_sexf_vp ^ SEXF)
    q  <- exp(lq  + etalq)
    d  <- exp(ld  + etald)

    # 3. Micro-constants for explicit two-compartment ODEs.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 4. ODE system: zero-order input into depot over duration D, sequential
    #    first-order absorption to central at rate ka, two-compartment
    #    disposition with linear elimination from central.
    d / dt(depot)       <- -ka * depot
    d / dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d / dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # 5. Zero-order release into the depot (duration D), then first-order
    #    transfer to central at rate ka. Bioavailability is anchored at 1
    #    for wild-type and scaled by (1 + 0.257) for ABCB1 rs3842 carriers;
    #    log-normal IIV on Frel multiplies the typical value.
    dur(depot) <- d
    f(depot)   <- (1 + e_rs3842_fdepot * SNP_ABCB1_RS3842) * exp(lfdepot + etalfdepot)

    # 6. Observation: efavirenz concentration in mg/L (dose mg / Vc L).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
