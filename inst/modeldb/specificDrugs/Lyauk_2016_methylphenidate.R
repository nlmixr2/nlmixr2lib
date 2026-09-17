Lyauk_2016_methylphenidate <- function() {
  description <- "Two-compartment population PK model for d-methylphenidate in healthy adults after a single oral dose of racemic methylphenidate, with Savic transit-compartment absorption and CES1 pharmacogenetic (rs71647871, rs115629050, CES1A2 diplotype) plus body-weight and sex covariate effects"
  reference <- paste(
    "Lyauk YK, Stage C, Bergmann TK, Ferrero-Milliani L, Bjerre D,",
    "Thomsen R, Dalhoff KP, Rasmussen HB, Jurgens G.",
    "Population Pharmacokinetics of Methylphenidate in Healthy Adults",
    "Emphasizing Novel and Known Effects of Several Carboxylesterase 1",
    "(CES1) Variants. Clin Transl Sci. 2016;9(6):337-345.",
    "doi:10.1111/cts.12423. PMCID: PMC5351003.",
    sep = " "
  )
  vignette <- "Lyauk_2016_methylphenidate"

  # The analysed observations are d-MPH plasma concentrations in ug/L
  # (numerically identical to ng/mL); the administered dose is the RACEMIC
  # dl-MPH amount in mg. The NONMEM control stream in Supplementary Material 1
  # carries `;DV == ug/l`, `;TIME == hours` and
  # `;AMT 10 mg Ritalin (dl-Methylphenidate)`, and scales the central
  # compartment with `S2 = V2/1000` (mg / L -> ug/L).
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometrically scaled with exponents FIXED at 0.75 on CL/F and Q/F",
        "and at 1 on Vdcentral/F and Vdperipheral/F, normalised to a 70 kg",
        "reference subject (Table 1 'Weight exponent' rows; Results paragraph",
        "3 states estimation of the exponents produced no significant change",
        "in OFV or parameter estimates). Cohort mean 72.5 kg (SD 12.9),",
        "Supplementary Table S1; the paper's own large-population simulation",
        "(Supplementary Material 2) truncates body weight to 46.1-113.8 kg.",
        "Time-invariant (single-dose studies)."
      ),
      source_name = "WT"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = male",
      notes = paste(
        "Acts on the mean transit time MTT as 1 + 0.925 * SEXF, i.e. females",
        "have an approximately 1.93-fold longer MTT than males. The source",
        "NONMEM column is `SEX` with `IF(SEX.EQ.0) MTTSEX = 1` (male, the",
        "most common category) and `IF(SEX.EQ.1) MTTSEX = (1 + THETA(16))`",
        "(Supplementary Material 1), so SEX maps onto the canonical SEXF",
        "directly with NO value inversion. Cohort 60 of 122 female (49.2%),",
        "Supplementary Table S1."
      ),
      source_name = "SEX"
    ),
    SNP_CES1_RS71647871 = list(
      description = "CES1 rs71647871 (c.428G>A, p.Gly143Glu) variant-allele carrier indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = rs71647871 GG (wild-type homozygote)",
      notes = paste(
        "1 = GA heterozygote (no GG variant homozygotes were observed; the",
        "cohort contained 110 wild-type, 6 heterozygous and 0 homozygous",
        "subjects, Supplementary Table S2). Reduces CL/F by 58.7%, the",
        "largest single covariate effect in the model. NOTE the value",
        "INVERSION relative to the source NONMEM data set, whose column",
        "`X4A` codes 1 = wild-type and 0 = the GA variant",
        "(`IF(X4A.EQ.1) CLX4A = 1` / `IF(X4A.EQ.0) CLX4A = (1 + THETA(13))`,",
        "Supplementary Material 1): SNP_CES1_RS71647871 = 1 - X4A, which is",
        "also the orientation of the paper's own printed typical-value",
        "equation (Results, where the variable `rs71647871` takes the value",
        "one for individuals possessing the reported genotype). Pair with",
        "SNP_CES1_RS71647871_MISSING, which flags ungenotyped subjects.",
        "Time-invariant (germline genotype)."
      ),
      source_name = "X4A"
    ),
    SNP_CES1_RS71647871_MISSING = list(
      description = "CES1 rs71647871 genotype-missing indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = rs71647871 genotype measured",
      notes = paste(
        "1 = genotype could not be determined. Missing in 5% of the studied",
        "population (Results paragraph 3). Carries its own estimated CL/F",
        "effect under the EXTRA (EST) full-maximum-likelihood method for data",
        "missing not at random (Keizer 2012; Johansson & Karlsson 2013), which",
        "adds one parameter per covariate containing missing data rather than",
        "imputing a genotype. Source coding is the sentinel value 999 in the",
        "`X4A` column (`IF(X4A.EQ.999) CLX4A = (1 + THETA(12))`,",
        "Supplementary Material 1). Member of the established <COV>_MISSING",
        "family alongside SNP_CYP2C19_RS3814637_MISSING and",
        "SNP_CYP3A4_RS2242480_MISSING. Mutually exclusive with",
        "SNP_CES1_RS71647871 = 1."
      ),
      source_name = "X4A (sentinel 999)"
    ),
    SNP_CES1_RS115629050 = list(
      description = "CES1 rs115629050 (p.Ala270Ser) variant-allele carrier indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = rs115629050 GG (wild-type homozygote)",
      notes = paste(
        "1 = TG heterozygote (64 wild-type, 7 heterozygous, 0 homozygous of",
        "the 71 successfully genotyped subjects, Supplementary Table S2).",
        "Reduces CL/F by 40.3%. This is the paper's novel finding - the first",
        "report of an rs115629050 effect on methylphenidate PK - and the",
        "Discussion notes that in-vitro work on CES1-metabolised ACE",
        "inhibitors did not reproduce it, so the effect should be treated as",
        "requiring confirmation. Source column `X7A` codes 0 = wild-type and",
        "1 = the TG variant, so this canonical matches the source orientation",
        "directly (NO inversion; contrast SNP_CES1_RS71647871 above).",
        "Time-invariant (germline genotype)."
      ),
      source_name = "X7A"
    ),
    SNP_CES1_RS115629050_MISSING = list(
      description = "CES1 rs115629050 genotype-missing indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = rs115629050 genotype measured",
      notes = paste(
        "1 = genotype could not be determined. Missing in 42% of the studied",
        "population (Results paragraph 3) - by far the largest missing",
        "fraction in the model, because the genetic assay cannot resolve",
        "sequence downstream of CES1A1 exon 5 when CES1A2 is present and",
        "rs115629050 lies in exon 7 (Discussion). Handled by the EXTRA (EST)",
        "method; source coding is the sentinel 999 in `X7A`",
        "(`IF(X7A.EQ.999) CLX7A = (1 + THETA(14))`, Supplementary Material 1).",
        "Its estimated effect (+9.0%) is the only covariate coefficient in the",
        "model whose bootstrap 95% CI spans zero (-0.0623 to 0.262, Table 1),",
        "i.e. ungenotyped subjects behave much like wild-type. Mutually",
        "exclusive with SNP_CES1_RS115629050 = 1."
      ),
      source_name = "X7A (sentinel 999)"
    ),
    CES1_HAPA2_HET = list(
      description = "CES1A2-bearing haplotype heterozygote indicator (exactly one CES1A2 copy)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no CES1A2 copy (CES1A1-CES1P1 / CES1A1-CES1P1, the most common diplotype)",
      notes = paste(
        "A CES1 haplotype carries either CES1P1 or the hybrid gene CES1A2 in",
        "addition to CES1A1 (Supplementary Methods 1.1), so the diplotype is a",
        "three-level germline dosage: no CES1A2 copy, one copy, or two copies.",
        "This indicator flags the one-copy (heterozygous) group and reduces",
        "CL/F by 18.2%; pair with CES1_HAPA2_HOM for the two-copy group, both",
        "zero being the no-CES1A2 reference. 82 / 34 / 5 subjects carried",
        "zero / one / two copies (Supplementary Table S2). Source column",
        "`NCOP` carries the total CES1A1 + CES1A2 copy number, so",
        "CES1_HAPA2_HET = as.integer(NCOP == 3) (`IF(NCOP.EQ.3) CLNCOP =",
        "(1 + THETA(9))`, Supplementary Material 1). Member of the",
        "<GENE>_HAP<haplotype>_HET / _HOM family established by",
        "SLCO1B1_HAP15_HET / SLCO1B1_HAP15_HOM. Time-invariant."
      ),
      source_name = "NCOP (== 3)"
    ),
    CES1_HAPA2_HOM = list(
      description = "CES1A2-bearing haplotype homozygote indicator (two CES1A2 copies)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no CES1A2 copy (CES1A1-CES1P1 / CES1A1-CES1P1, the most common diplotype)",
      notes = paste(
        "Flags the two-copy (homozygous CES1A1-CES1A2 / CES1A1-CES1A2) group",
        "and reduces CL/F by 41.0%; paired with CES1_HAPA2_HET, both zero",
        "being the no-CES1A2 reference. Only 5 of 121 genotyped subjects",
        "carried two copies (Supplementary Table S2), and the paper's",
        "Discussion flags that the CES1A2 direction conflicts with an",
        "irinotecan study and with an oseltamivir study that found no",
        "diplotype effect, so the effect requires confirmation. Source coding",
        "CES1_HAPA2_HOM = as.integer(NCOP == 4) (`IF(NCOP.EQ.4) CLNCOP =",
        "(1 + THETA(10))`, Supplementary Material 1). Member of the",
        "<GENE>_HAP<haplotype>_HET / _HOM family established by",
        "SLCO1B1_HAP15_HET / SLCO1B1_HAP15_HOM. Time-invariant."
      ),
      source_name = "NCOP (== 4)"
    ),
    CES1_HAPA2_MISSING = list(
      description = "CES1A2 diplotype (copy-number) missing indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = CES1A2 copy number measured",
      notes = paste(
        "1 = CES1A2 copy number could not be determined; missing in 0.8%",
        "(1 of 122) of the studied population (Results paragraph 3,",
        "Supplementary Table S2). Handled by the EXTRA (EST) method with its",
        "own estimated CL/F coefficient; source coding is the sentinel 999 in",
        "`NCOP` (`IF(NCOP.EQ.999) CLNCOP = (1 + THETA(11))`, Supplementary",
        "Material 1). The estimate (-53.5%) is large and precise despite",
        "resting on a single subject, so it should be read as a nuisance",
        "parameter absorbing that subject rather than as a transferable",
        "effect. Mutually exclusive with CES1_HAPA2_HET and CES1_HAPA2_HOM."
      ),
      source_name = "NCOP (sentinel 999)"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "methylphenidate",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "dexmethylphenidate",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "dexmethylphenidate",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 122,
    n_studies = 2,
    age_range = "mean 23.5 years (SD 2.6)",
    age_median = "23.5 years (mean; median not reported)",
    weight_range = "mean 72.5 kg (SD 12.9)",
    weight_median = "72.5 kg (mean; median not reported)",
    sex_female_pct = 49.2,
    race_ethnicity = c(White = 100),
    disease_state = "healthy volunteers",
    dose_range = "single oral 10 mg racemic (dl) methylphenidate (Ritalin) immediate-release tablet",
    regions = "Denmark",
    n_observations = 503,
    height_mean_cm = 176.4,
    bmi_mean = 23.2,
    notes = paste(
      "Pooled from two single-dose studies at Bispebjerg Hospital,",
      "Copenhagen (NCT02135263 and NCT02147535). Study I: 44 subjects,",
      "rich sampling predose and at 0.5, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 10, 24",
      "and 33 h. Study II: 78 subjects, a SINGLE sample at 3 h postdose.",
      "503 d-MPH plasma concentrations in total. All participants were",
      "Caucasian and were genotyped for CES1 from saliva-derived DNA.",
      "Baseline demographics: Supplementary Table S1; CES1 variant",
      "distribution: Supplementary Table S2. Participants fasted overnight",
      "and received methylphenidate 1 h after a standardised breakfast",
      "(two participants did not receive the standard meal)."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural parameters. All values are the final estimates from the
    # 'NONMEM estimate (RSE%)' column of Table 1. The $THETA block of the
    # control stream in Supplementary Material 1 holds the INITIAL estimates
    # (e.g. CL/F 233, V2 98.4, Q 69.1, V3 250, MTT 0.508) and is deliberately
    # NOT used here.
    #
    # Every disposition parameter is APPARENT (/F): only oral data were
    # collected, and the dose is the RACEMIC dl-MPH amount while the
    # observation is the d-enantiomer only, so F additionally absorbs the
    # d-fraction of the administered racemate.
    # ---------------------------------------------------------------------
    lmtt <- log(0.505)
    label("Mean transit time through the absorption transit chain (h)")
    # Table 1 theta 1 'MTT (/h)' = 0.505, RSE 13.7%; bootstrap 0.521, nonparametric 95% CI 0.371-0.707

    lntr <- fixed(log(3))
    label("Number of absorption transit compartments (unitless)")
    # Table 1 theta 2 'No. of transit compartments (fixed)' = 3. Results paragraph 2: originally estimated as 3.05 but fixed to 3.0 to aid model stability and convergence, with no significant change in fit or parameter estimates. Control stream $THETA '(3) FIX'

    lka <- log(0.418)
    label("First-order absorption rate constant out of the absorption compartment (1/h)")
    # Table 1 theta 3 'ka (/h)' = 0.418, RSE 18.8%; note the bootstrap RSE is 134% and the 95% CI 0.370-1.304 is markedly asymmetric, so ka is the least well determined structural parameter

    lcl <- log(233.0)
    label("Apparent oral clearance CL/F for a 70 kg wild-type-CES1 subject (L/h)")
    # Table 1 theta 4 'CL/F (L/h)' = 233.0, RSE 3.6%; bootstrap 231.7, 95% CI 214.0-251.8

    lvc <- log(97.6)
    label("Apparent central volume of distribution Vdcentral/F for a 70 kg subject (L)")
    # Table 1 theta 5 'Vdcentral/F (L)' = 97.6, RSE 28.6%; bootstrap 101.8, 95% CI 47.6-288.7

    lq <- log(70.1)
    label("Apparent intercompartmental clearance Q/F for a 70 kg subject (L/h)")
    # Table 1 theta 6 'Q/F (L/h)' = 70.1, RSE 47.8%; bootstrap 74.3, 95% CI 49.9-175.3

    lvp <- log(252)
    label("Apparent peripheral volume of distribution Vdperipheral/F for a 70 kg subject (L)")
    # Table 1 theta 7 'Vdperipheral/F (L)' = 252, RSE 30.9%; bootstrap 260.7, 95% CI 198.2-453.2

    # ---------------------------------------------------------------------
    # Allometric body-weight exponents. Both were FIXED, not estimated:
    # Table 1 lists them in dedicated '(fixed)' rows with no RSE and no
    # bootstrap CI, and Results paragraph 3 states that estimating them
    # produced no significant change in OFV or parameter estimates. One
    # exponent is shared by CL/F and Q/F, the other by the two volumes, so
    # each is encoded once under the registered shared-exponent name.
    # ---------------------------------------------------------------------
    e_wt_cl_q <- fixed(0.75)
    label("Allometric exponent of body weight on CL/F and Q/F (unitless)")
    # Table 1 rows 'Weight exponent on CL/F (fixed)' = 0.75 and 'Weight exponent on Q/F (fixed)' = 0.75; control stream `TVCL = THETA(1)*CLCOV*(WT/70)**0.75` and `Q=THETA(5)*(WT/70)**0.75`

    e_wt_vc_vp <- fixed(1)
    label("Allometric exponent of body weight on Vdcentral/F and Vdperipheral/F (unitless)")
    # Table 1 rows 'Weight exponent on Vdcentral/F (fixed)' = 1 and 'Weight exponent on Vdperipheral/F (fixed)' = 1; control stream `TVV2 = THETA(2)*(WT/70)` and `V3=THETA(6)*(WT/70)`

    # ---------------------------------------------------------------------
    # Covariate effects. Every effect is a LINEAR PROPORTIONAL shift of the
    # form (1 + theta * indicator), per the typical-value equations printed
    # in the Results section and the `( 1 + THETA(n))` blocks of the control
    # stream. The genetic effects multiply each other on CL/F; sex acts on
    # MTT alone.
    # ---------------------------------------------------------------------
    e_sexf_mtt <- 0.925
    label("Proportional change in MTT for female sex (fraction)")
    # Table 1 theta 8 'Female gender on MTT' = 0.925, RSE 44.3%; bootstrap 0.852, 95% CI 0.297-1.831. Results equation MTT = theta1 * (1 + Gender * theta8)

    e_snp_ces1_rs71647871_cl <- -0.587
    label("Proportional change in CL/F for CES1 rs71647871 GA carriers (fraction)")
    # Table 1 theta 9 'rs71647871 (GA) on CL/F' = -0.587, RSE 6.9%; bootstrap -0.585, 95% CI -0.677 to -0.473

    e_snp_ces1_rs71647871_missing_cl <- -0.157
    label("Proportional change in CL/F for subjects with an undetermined CES1 rs71647871 genotype (fraction)")
    # Table 1 theta 10 'Missing data rs71647871 on CL/F' = -0.157, RSE 46.9%; bootstrap -0.161, 95% CI -0.330 to -0.00164

    e_ces1_hapa2_het_cl <- -0.182
    label("Proportional change in CL/F for carriers of one CES1A2 copy (fraction)")
    # Table 1 theta 11 'One CES1A2 on CL/F' = -0.182, RSE 31.9%; bootstrap -0.184, 95% CI -0.307 to -0.0475

    e_ces1_hapa2_hom_cl <- -0.410
    label("Proportional change in CL/F for carriers of two CES1A2 copies (fraction)")
    # Table 1 theta 12 'Two CES1A2 on CL/F' = -0.410, RSE 22.5%; bootstrap -0.408, 95% CI -0.583 to -0.192

    e_ces1_hapa2_missing_cl <- -0.535
    label("Proportional change in CL/F for subjects with an undetermined CES1A2 diplotype (fraction)")
    # Table 1 theta 13 'Missing data CES1A diplotype on CL/F' = -0.535, RSE 9.0%; bootstrap -0.535, 95% CI -0.663 to -0.440

    e_snp_ces1_rs115629050_cl <- -0.403
    label("Proportional change in CL/F for CES1 rs115629050 TG carriers (fraction)")
    # Table 1 theta 14 'rs115629050 (TG) on CL/F' = -0.403, RSE 26.1%; bootstrap -0.399, 95% CI -0.714 to -0.180

    e_snp_ces1_rs115629050_missing_cl <- 0.090
    label("Proportional change in CL/F for subjects with an undetermined CES1 rs115629050 genotype (fraction)")
    # Table 1 theta 15 'Missing data rs115629050 on CL/F' = 0.090, RSE 79.6%; bootstrap 0.0913, 95% CI -0.0623 to 0.262 (the only covariate coefficient whose CI spans zero)

    # ---------------------------------------------------------------------
    # Inter-individual variability. Results paragraph 2: IIV was estimated
    # for CL/F, Vdcentral/F and MTT 'in a full variance-covariance matrix,
    # containing both diagonal and nondiagonal elements', matching the
    # $OMEGA BLOCK(3) of the control stream whose row order is CL, V2, MTT.
    #
    # SCALE OF THE PRINTED 'IIV (%CV)' COLUMN. Table 1 reports CL/F 21.6,
    # Vdcentral/F 90.1 and MTT 62.1 %CV. Those are omega on the log scale
    # expressed as a percentage, i.e. %CV = 100 * sqrt(omega^2) - NOT the
    # log-normal 100 * sqrt(exp(omega^2) - 1). The control stream's INITIAL
    # $OMEGA diagonal settles it: 0.0478 / 0.78 / 0.383 give 21.9 / 88.3 /
    # 61.9 under sqrt(omega^2) against the printed 21.6 / 90.1 / 62.1, but
    # 22.1 / 108.7 / 68.3 under the exponential form. Every initial-to-final
    # THETA shift in this run is under 1.5%, and only the sqrt reading keeps
    # the omega shifts in that range (the exponential reading would demand a
    # 17% move on Vdcentral/F alone). Variances below are therefore
    # (%CV / 100)^2: 0.216^2, 0.901^2, 0.621^2.
    #
    # OFF-DIAGONAL PROVENANCE. The final off-diagonal covariances are NOT
    # published - Table 1 prints only the three diagonal %CV values. The
    # covariances used here are the INITIAL estimates from the $OMEGA
    # BLOCK(3) of Supplementary Material 1 (cov(CL,V2) = 0.0669,
    # cov(CL,MTT) = 0.01, cov(V2,MTT) = 0.01), combined with the final
    # variances above. The resulting correlations (0.34, 0.075, 0.018) are
    # all admissible and the block is positive definite, but the two smaller
    # values are the conventional NONMEM 0.01 starting value and should be
    # treated as weakly determined. See the vignette's Assumptions and
    # deviations section.
    # ---------------------------------------------------------------------
    etalcl + etalvc + etalmtt ~ c(
      0.046656,
      0.066900, 0.811801,
      0.010000, 0.010000, 0.385641
    )
    # diagonal from Table 1 'IIV (%CV)' rows CL/F 21.6, Vdcentral/F 90.1, MTT 62.1 as (%CV/100)^2; off-diagonals are the initial estimates of $OMEGA BLOCK(3) in Supplementary Material 1

    propSd <- 0.184
    label("Proportional residual error on d-MPH plasma concentration (fraction)")
    # Table 1 'Residual variability / Proportional error' = 0.184, RSE 9.6%; bootstrap 0.181, 95% CI 0.144-0.218. Control stream $ERROR uses `W = THETA(4)*IPRED` with `Y = IPRED + W*EPS(1)` and $SIGMA 1 FIX, so THETA(4) is the proportional residual SD itself
  })

  model({
    # -------------------------------------------------------------------
    # 1. Covariate multipliers
    #
    # Results typical-value equations:
    #   MTT  = theta1 * (1 + Gender * theta8)
    #   CL/F = theta4 * (1 + rs71647871 * theta9)
    #                 * (1 + One CES1A2 * theta11)
    #                 * (1 + Two CES1A2 * theta12)
    #                 * (1 + rs115629050 * theta14)
    #                 * (WT/70)^0.75
    #
    # The printed CL/F equation lists only the observed-genotype terms. The
    # three EXTRA-method missing-genotype terms (theta 10, 13, 15) are
    # tabulated in Table 1 and appear in the control stream's `CLCOV =
    # CLNCOP*CLX4A*CLX7A` product, where each of the three factors selects
    # the wild-type (1), the variant, or the missing-data multiplier. They
    # are included here so the encoding reproduces the paper's own fitted
    # model rather than only its reference-genotype special case; for a
    # fully genotyped subject all three indicators are 0 and the product
    # collapses to the printed equation exactly.
    # -------------------------------------------------------------------
    mttcov <- 1 + e_sexf_mtt * SEXF

    clcov <-
      (1 + e_snp_ces1_rs71647871_cl * SNP_CES1_RS71647871) *
      (1 + e_snp_ces1_rs71647871_missing_cl * SNP_CES1_RS71647871_MISSING) *
      (1 + e_ces1_hapa2_het_cl * CES1_HAPA2_HET) *
      (1 + e_ces1_hapa2_hom_cl * CES1_HAPA2_HOM) *
      (1 + e_ces1_hapa2_missing_cl * CES1_HAPA2_MISSING) *
      (1 + e_snp_ces1_rs115629050_cl * SNP_CES1_RS115629050) *
      (1 + e_snp_ces1_rs115629050_missing_cl * SNP_CES1_RS115629050_MISSING)

    # -------------------------------------------------------------------
    # 2. Individual parameters. IIV is carried by CL/F, Vdcentral/F and MTT
    # only; ka, Q/F and Vdperipheral/F are typical values with no eta, per
    # Table 1 and the control stream ($PK assigns ETA(1) to CL, ETA(2) to
    # V2 and ETA(3) to MTT, with `KA = THETA(3)`, `Q=THETA(5)*...` and
    # `V3=THETA(6)*...` carrying none).
    # -------------------------------------------------------------------
    cl <- exp(lcl + etalcl) * clcov * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q <- exp(lq) * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp) * (WT / 70)^e_wt_vc_vp
    ka <- exp(lka)
    mtt <- exp(lmtt + etalmtt) * mttcov
    ntr <- exp(lntr)

    # -------------------------------------------------------------------
    # 3. Micro-constants. Figure 1 defines ktr as (number of transit
    # compartments + 1) divided by the mean transit time, which the control
    # stream implements as `KTR = (NN+1)/MTT`. With NN fixed at 3 this is
    # 4/MTT.
    # -------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    ktr <- (ntr + 1) / mtt

    # -------------------------------------------------------------------
    # 4. Savic (2007) analytical transit-absorption input, transcribed
    # verbatim from $DES / $PK of Supplementary Material 1:
    #
    #   LNFAC   = LOG(2.5066)+(NN+0.5)*LOG(NN)-NN
    #   LNDK    = LOG(DOSE+0.00001)+LOG(KTR)
    #   X       = KTR*T
    #   DADT(1) = EXP(LNDK+NN*LOG(X+0.00001)-X-LNFAC) - KA*A(1)
    #
    # LNFAC is STIRLING'S approximation of log(NN!), not the exact
    # lgamma(NN+1), and it is reproduced here rather than replaced by
    # rxode2's builtin transit() (which uses the exact log-gamma). At
    # NN = 3 Stirling gives 3! ~ 5.8355 against the true 6, so this input
    # function delivers 6/5.8355 = 1.0281 times the administered dose over
    # the whole profile. That 2.8% excess is baked into the published CL/F,
    # so substituting the exact factorial would be locally more correct and
    # globally inconsistent with Table 1. The vignette's mass-balance gate
    # therefore targets Dose * 1.028065 rather than Dose.
    #
    # The two 0.00001 offsets are the paper's own guards against log(0) at
    # the moment of dosing and are kept. podo(depot) returns the raw dose
    # amount before bioavailability, matching NONMEM's PODO, and tad(depot)
    # is time after the most recent dose, matching T in these single-dose
    # studies.
    # -------------------------------------------------------------------
    tad_dose <- tad(depot)
    x_tr <- ktr * tad_dose
    lfact_ntr <- log(2.5066) + (ntr + 0.5) * log(ntr) - ntr
    transit_in <- exp(
      log(podo(depot) + 0.00001) + log(ktr) +
        ntr * log(x_tr + 0.00001) - x_tr - lfact_ntr
    )

    # -------------------------------------------------------------------
    # 5. ODE system. $MODEL declares COMP=(ABS,DEFDOSE), COMP=(CENT,DEFOBS)
    # and COMP=(PERI); ABS is the absorption compartment fed by the transit
    # input and emptied at ka, mapped here onto the canonical `depot`.
    # -------------------------------------------------------------------
    d/dt(depot) <- transit_in - ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # -------------------------------------------------------------------
    # 6. `F1 = 0` in $PK suppresses the dose bolus into the absorption
    # compartment so that the analytical transit term above is the sole
    # input pathway. podo(depot) is unaffected by f(depot), so the transit
    # term still sees the full administered amount.
    # -------------------------------------------------------------------
    f(depot) <- 0

    # -------------------------------------------------------------------
    # 7. Observation. The control stream sets `S2 = V2/1000`, i.e. the
    # central amount in mg divided by the volume in L is multiplied by 1000
    # to give the observed d-MPH concentration in ug/L (= ng/mL).
    # -------------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
