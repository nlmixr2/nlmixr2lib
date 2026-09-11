Li_2025_lacosamide_cyp2c19 <- function() {
  description <- "One-compartment population PK model with first-order absorption for oral lacosamide in Chinese children with epilepsy, with a body-weight power function and CYP2C19*2 (rs12769205) genotype on apparent clearance; Model II of Li 2025, for the clinical scenario in which CYP2C19 genotype is available"
  reference <- paste(
    "Li Y, Guo HL, Fan L, Wang J, Hu YH, Zhang YY, Qiu JC, Chen J, Wu CF,",
    "Zhang G, Lu XP, Chen F (2025).",
    "PopPK modeling supports BW band dosing of lacosamide for pediatric epilepsy.",
    "npj Genomic Medicine 10:80.",
    "doi:10.1038/s41525-025-00519-y. PMCID: PMC12394691.",
    sep = " "
  )
  vignette <- "Li_2025_lacosamide"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters apparent clearance as a power function normalized to the cohort",
        "median of 30 kg: CL/F = 1.7 * (WT/30)^0.319 * theta_rs12769205",
        "(Results Eq. 3; Table 2 'Covariate model structure'). The 30 kg",
        "reference is the median body weight of the model-development group",
        "(Table 1: BW 30 kg, range 10-80 kg), which is what Supplementary",
        "Eq. 13 defines as COVmedian. Supplementary Table S5 confirms body",
        "weight remained 'the covariate exerting the greatest influence on",
        "CL/F' after the genotype term was added. Body weight does NOT enter",
        "V/F: the V-BW allometric term was rejected for parameter instability",
        "(RSE 327.4%; Results, paragraph following Eq. 4), so V/F is a single",
        "population value of 26.4 L for every subject.",
        sep = " "
      ),
      source_name        = "BW"
    ),
    SNP_CYP2C19_RS12769205_GA = list(
      description        = "CYP2C19*2 rs12769205 heterozygous GA genotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (AA genotype), jointly with SNP_CYP2C19_RS12769205_GG = 0",
      notes              = paste(
        "1 = subject's reported rs12769205 genotype is GA; 0 = otherwise. Used",
        "together with SNP_CYP2C19_RS12769205_GG; both are 0 for the AA",
        "reference group, which Li 2025 calls the wild type. Time-fixed per",
        "subject (germline genotype). Genotype counts in the model-development",
        "group were AA/GA/GG = 56/62/15 (Table 1). Carriers have 0.879-fold the",
        "apparent clearance of AA subjects, i.e. about 12 percent lower CL/F and",
        "correspondingly higher trough concentrations. rs12769205 is the",
        "CYP2C19*2 tag single-nucleotide polymorphism genotyped on the Agena",
        "MassARRAY platform (Methods, 'Genotyping'); it was the only one of nine",
        "screened polymorphisms across ABCB1, ABCC2, CYP2C9 and CYP2C19 to",
        "affect CL/F significantly (dOFV 16.465, p < 0.001; Supplementary",
        "Table S5).",
        sep = " "
      ),
      source_name        = "rs12769205"
    ),
    SNP_CYP2C19_RS12769205_GG = list(
      description        = "CYP2C19*2 rs12769205 homozygous GG genotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (AA genotype), jointly with SNP_CYP2C19_RS12769205_GA = 0",
      notes              = paste(
        "1 = subject's reported rs12769205 genotype is GG; 0 = otherwise. Used",
        "together with SNP_CYP2C19_RS12769205_GA; both are 0 for the AA",
        "reference group. Time-fixed per subject (germline genotype). Only 15 of",
        "133 subjects in the model-development group were GG (Table 1), so this",
        "coefficient rests on a small subgroup. GG subjects have 0.736-fold the",
        "apparent clearance of AA subjects; the Discussion states this as",
        "'homozygous carriers exhibited a 26.4% reduction in CL/F (i.e.,",
        "0.736 x wild-type; p < 0.001)', which is consistent with the loss of",
        "function annotated for CYP2C19*2 in PharmGKB.",
        sep = " "
      ),
      source_name        = "rs12769205"
    )
  )

  # Covariates that Li 2025 screened but did NOT retain in Model II.
  # Documented here for provenance only; none is referenced in model().
  covariatesDataExcluded <- list(
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Screened on CL/F with a proportional model (Supplementary Eq. 16) and",
        "rejected at forward inclusion: dOFV -0.173, p > 0.05 (Supplementary",
        "Table S5). Cohort was 52 of 133 female (Table 1).",
        sep = " "
      ),
      source_name        = "SEX"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Offered to covariate screening (Methods, 'Covariate model') and not",
        "retained on any parameter. Median 7.5 years, range 1-18 years",
        "(Table 1). Age enters only indirectly through body weight; two explicit",
        "maturation parameterizations (Supplementary Eqs. 8-11) were fitted and",
        "both abandoned for excessive RSE and shrinkage (> 80%; Supplementary",
        "S.4.4), so the final model carries no maturation term.",
        sep = " "
      ),
      source_name        = "Age"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Three eGFR formulations were screened on CL/F -- Shull's creatinine",
        "equation and two cystatin-C equations (Supplementary Eqs. 18-20;",
        "Table 1 rows eGFR a/b/c) -- and all three were rejected at forward",
        "inclusion (Supplementary Table S5). Supplementary S.4.4 states that",
        "'renal function did not affect CL/F of LCM'.",
        sep = " "
      ),
      source_name        = "eGFR"
    ),
    RBC = list(
      description        = "Red blood cell (erythrocyte) count",
      units              = "10^12/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "One of five covariates that reached the forward-inclusion threshold on",
        "CL/F (dOFV > 3.84, p < 0.05) -- with sodium channel blockers, serum",
        "creatinine, potassium and mean corpuscular hemoglobin concentration --",
        "but none survived backward elimination (dOFV < 10.83, p > 0.001;",
        "Results 'Model development'; Supplementary Table S5). Representative of",
        "that group; the full screened panel is tabulated in Supplementary",
        "Table S5 and listed in the vignette rather than enumerated here. The",
        "eight non-CYP2C19*2 polymorphisms screened alongside rs12769205",
        "(ABCB1 rs1045642, rs2032582, rs3789243; ABCC2 rs3740066, rs717620;",
        "CYP2C9 rs1057910; CYP2C19 rs4986893 and rs3758581) were likewise all",
        "rejected.",
        sep = " "
      ),
      source_name        = "RBC"
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "lacosamide", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "lacosamide", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species          = "human",
    n_subjects       = 133,
    n_observations   = 347,
    n_studies        = 1,
    age_range        = "1-18 years",
    age_median       = "7.5 years",
    weight_range     = "10-80 kg",
    weight_median    = "30 kg",
    sex_female_pct   = 39.1,
    race_ethnicity   = c(Asian = 100),
    disease_state    = "epilepsy diagnosed by ILAE criteria; focal or generalized seizures",
    dose_range       = "oral tablet twice daily; per-dose 2.0-8 mg/kg in children under 50 kg and 75-200 mg in children 50 kg and over",
    concentration_range = "0.70-11.90 mg/L (steady-state trough)",
    target_range     = "2-7 mg/L steady-state trough",
    genotype_distribution = c(AA = 56, GA = 62, GG = 15),
    regions          = "China (single centre, Children's Hospital of Nanjing Medical University)",
    notes            = paste(
      "Same cohort as the body-weight-only sibling Li_2025_lacosamide.R:",
      "a retrospective real-world therapeutic-drug-monitoring cohort collected",
      "between June 2021 and March 2023, in which 190 children contributing 493",
      "concentrations were randomly split 70/30 into a model-development group",
      "(133 children, 347 concentrations -- the numbers recorded above) and an",
      "external validation group (57 children, 146 concentrations). Baseline",
      "demographics are Table 1. genotype_distribution gives the rs12769205",
      "AA/GA/GG counts in the development group. All samples are pre-dose",
      "troughs drawn 30 min before the next maintenance dose after at least 3",
      "days on an unchanged regimen, which is why ka is fixed and V/F carries no",
      "inter-individual variability. Race is recorded as 100 percent Asian",
      "because the cohort is a Chinese single-centre pediatric population; the",
      "paper reports no formal race or ethnicity breakdown. The Discussion notes",
      "that the CYP2C19*2 allele frequency exceeds 30 percent in Asians and",
      "cautions that these genotype coefficients should be 'contextually",
      "adapted' for populations with different pharmacogenetic profiles, where",
      "CYP2C19*17 rather than *2 may predominate.",
      sep = " "
    )
  )

  ini({
    # ================================================================
    # Structural parameters -- Li 2025 Table 2, "Model II" block,
    # Estimate column. Same one-compartment first-order-absorption
    # structure as Model I (NONMEM ADVAN2 TRANS2, Methods 'Base
    # model'), refitted with the CYP2C19*2 genotype term on CL/F, so
    # the typical CL/F and V/F both shift upward relative to Model I.
    #
    # CL/F and V/F are APPARENT values (Methods, 'Base model': "as
    # bioavailability (F) could not be determined, CL and Vd were
    # expressed as apparent values").
    #
    # UNITS OF CL/F. Table 2's row label is "CL/F (L/h)"; the Results
    # sentence following Eq. 4 misprints it as "1.7 L/h/kg", the same
    # slip Model I's text carries. Table 2 governs and the per-kg
    # reading is arithmetically impossible against the observed
    # concentration range; see the vignette Errata.
    # ================================================================
    lka <- fixed(log(2.45)); label("Absorption rate constant ka (1/h)")                     # Table 2 Model II: "Ka 2.45 (fixed)", the same value fixed in Model I. Methods 'Base model': ka "could not be reliably estimated. Therefore, ka was fixed at 2.45 h-1, based on values reported in a previous study" (the paper's reference 15). No RSE reported.
    lcl <- log(1.7);         label("Apparent clearance CL/F for the AA reference genotype at the reference body weight of 30 kg (L/h)") # Table 2 Model II: CL/F = 1.7 L/h, RSE 4.7% (bootstrap median 1.69, 95% CI 1.47-1.91). Also the leading coefficient of Results Eq. 3. This is the AA (theta = 1) value.
    lvc <- log(26.4);        label("Apparent volume of distribution V/F (L)")                # Table 2 Model II: V/F = 26.4 L, RSE 13.4% (bootstrap median 26.1, 95% CI 16.8-40.0). Also Results Eq. 4, which carries no covariate term.

    # ================================================================
    # Covariate effect on CL/F -- body weight, power function
    # normalized to the cohort median of 30 kg (Supplementary Eq. 13).
    # Table 2's "Covariate model structure" footnote prints the
    # instantiated relationship as
    #   CL/F = 1.7 * (BW/30)^0.319 * theta_rs12769205,
    # which is Results Eq. 3.
    # ================================================================
    e_wt_cl <- 0.319; label("Power exponent on (WT/30 kg) for CL/F (unitless)")              # Table 2 Model II: "CL-BW" = 0.319, RSE 13.3% (bootstrap median 0.317, 95% CI 0.231-0.408). Also the exponent printed in Results Eq. 3.

    # ================================================================
    # Covariate effects on CL/F -- CYP2C19*2 rs12769205 genotype.
    #
    # These two parameters are MULTIPLICATIVE FACTORS, not log-scale or
    # exponential coefficients. Table 2 lists them as "CL-rs12769205*GA
    # = 0.879" and "CL-rs12769205*GG = 0.736", and Results Eq. 3 applies
    # them as a bare product, "* theta_rs12769205", with the gloss
    # "theta = 1 for CYP2C19*2 AA". So AA is the reference with an
    # implicit factor of 1 and the tabulated values are the fold-change
    # in CL/F relative to AA. Supplementary Table S5 records the term as
    # entering through Eq. 16, Pi = TV(P)*(1 + theta*COV), which is the
    # same multiplicative family (theta_GA = -0.121, theta_GG = -0.264
    # in that parameterization); the fold-change form is used here
    # because it is what Table 2 prints.
    #
    # WHICH VALUES. The paper states the genotype factors twice and
    # disagrees with itself: Table 2 gives 0.879 / 0.736 while the
    # Results gloss under Eq. 4 gives 0.887 / 0.744. Table 2 is used
    # here, on three grounds. (1) The bootstrap columns of Table 2
    # independently report medians of 0.880 and 0.736. (2) The
    # Discussion derives "a 26.4% reduction in CL/F (i.e., 0.736 x
    # wild-type)", and 1 - 0.736 = 0.264 exactly, so 0.736 is
    # arithmetically self-confirming while 0.744 is not. (3) The two
    # gloss values are each exactly 0.008 above the tabulated ones -- a
    # single systematic slip in one sentence rather than two independent
    # estimates. Table 2 therefore carries the vote 3-1 for GG and 2-1
    # for GA. See the vignette Errata.
    #
    # Sign and monotonicity check: both factors are below 1 and decrease
    # with G-allele count (AA 1 > GA 0.879 > GG 0.736), i.e. each copy
    # of the *2 loss-of-function allele lowers apparent clearance. This
    # is the direction the Discussion asserts and the direction
    # annotated for CYP2C19*2 in PharmGKB. Consistent.
    # ================================================================
    e_snp_cyp2c19_rs12769205_ga_cl <- 0.879; label("Multiplicative factor on CL/F for CYP2C19*2 rs12769205 GA vs AA (unitless)") # Table 2 Model II: "CL-rs12769205*GA" = 0.879, RSE 4.0% (bootstrap median 0.880, 95% CI 0.818-0.950). The Results gloss under Eq. 4 says 0.887 instead; see the block comment above and the vignette Errata.
    e_snp_cyp2c19_rs12769205_gg_cl <- 0.736; label("Multiplicative factor on CL/F for CYP2C19*2 rs12769205 GG vs AA (unitless)") # Table 2 Model II: "CL-rs12769205*GG" = 0.736, RSE 5.5% (bootstrap median 0.736, 95% CI 0.661-0.825). Confirmed by the Discussion's "26.4% reduction in CL/F (i.e., 0.736 x wild-type)". The Results gloss under Eq. 4 says 0.744 instead; see the block comment above.

    # ================================================================
    # Inter-individual variability -- Li 2025 Table 2, "IIV in CL/F,
    # %CV" row of the Model II block. Exponential (log-normal) IIV on
    # CL/F only, on the same reasoning documented in
    # Li_2025_lacosamide.R: Supplementary Tables S4/S5 label the base
    # model "Eq.1 & Eq.7" (Eq. 1 being the exponential IIV form), and
    # Supplementary S.4.4 states that the IIV on V/F "was not
    # informative enough to be estimated and then was excluded from the
    # model" because only trough samples were available.
    #
    # Adding the genotype term reduced unexplained IIV in CL/F from
    # 18.0 %CV (Model I) to 15.6 %CV, which is the paper's central
    # pharmacogenetic claim.
    #
    # %CV -> variance uses the plain square, following the same table's
    # residual rows (0.554^2 = 0.307 and 0.174^2 = 0.0301, the two
    # variances quoted in Results 'Model evaluation'). The exact
    # log-normal alternative log(1 + 0.156^2) = 0.0241 differs by 0.6%
    # in omega; see the vignette Errata.
    # ================================================================
    etalcl ~ 0.024336  # Table 2 Model II: IIV in CL/F = 15.6 %CV, RSE 12% (bootstrap median 15.1, 95% CI 10.6-18.9); eta-shrinkage 22.6%. Variance = 0.156^2 = 0.024336.

    # ================================================================
    # Residual error -- Li 2025 Table 2, Model II block. The mixed
    # error model is Supplementary Eq. 7, Y = IPRED*(1 + eps1) + eps2,
    # giving Var(Y) = (sigma1*IPRED)^2 + sigma2^2, which is exactly
    # nlmixr2's prop(propSd) + add(addSd) combination. Both components
    # are unchanged from Model I.
    #
    # Table 2's Model II block prints "Proportional error, %CV" twice:
    # once as 17.4 with RSE 21.2% in the parameter position, and once as
    # 12.7 in the row immediately below the additive-error row, where
    # Model I's block instead reads "sigma-shrinkage (%) 14.1". The 12.7
    # entry is the sigma-shrinkage, mislabelled: Results 'Model
    # evaluation' says "both residual error components showed reduced
    # epsilon-shrinkage to 12.7% (vs. 14.1% previously)". So the Model
    # II proportional error is 17.4 %CV, identical to Model I, and 12.7
    # is not a second error estimate. See the vignette Errata.
    #
    # Concentrations are in ug/mL throughout the paper, numerically
    # identical to the declared mg/L, so the additive SD transfers
    # unchanged.
    # ================================================================
    propSd <- 0.174;        label("Proportional residual error (fraction)")                  # Table 2 Model II: "Proportional error, %CV" = 17.4, RSE 21.2% (bootstrap median 17.3, 95% CI 13.6-21.3). Results 'Model evaluation': the proportional error variance "remained stable at 0.0301" (= 0.174^2 to rounding), i.e. unchanged from Model I.
    addSd  <- fixed(0.554); label("Additive residual error SD on Cc (mg/L)")                 # Table 2 Model II: "Additive error, SD 0.554 (fixed)". No RSE reported, consistent with a fixed value. Results 'Model evaluation': the additive error variance "remained stable at ... 0.307" (= 0.554^2).
  })

  model({
    # ---- Derived covariate term -----------------------------------
    # theta_rs12769205 of Results Eq. 3: 1 for the AA reference
    # genotype, 0.879 for GA, 0.736 for GG. The two indicators are
    # mutually exclusive and both 0 for AA, so writing the factor as
    # 1 + (f_GA - 1)*I(GA) + (f_GG - 1)*I(GG) reproduces all three
    # levels while keeping Table 2's published fold-changes literal in
    # ini().
    genoCl <- 1 +
      (e_snp_cyp2c19_rs12769205_ga_cl - 1) * SNP_CYP2C19_RS12769205_GA +
      (e_snp_cyp2c19_rs12769205_gg_cl - 1) * SNP_CYP2C19_RS12769205_GG

    # ---- Individual PK parameters ---------------------------------
    # Body weight scales CL/F as a power function normalized to the
    # cohort median of 30 kg and the genotype factor multiplies it
    # (Results Eq. 3). V/F carries no covariate and no eta (Results
    # Eq. 4); ka is fixed.
    ka  <- exp(lka)
    cl  <- exp(lcl + etalcl) * (WT / 30)^e_wt_cl * genoCl
    vc  <- exp(lvc)

    # ---- Micro-constant -------------------------------------------
    kel <- cl / vc

    # ---- ODE system ------------------------------------------------
    # One compartment with first-order absorption and elimination,
    # matching the NONMEM ADVAN2 TRANS2 subroutine named in Methods,
    # 'Base model'. Oral tablets are dosed into depot; because F is not
    # identifiable, no bioavailability term is applied and CL/F and V/F
    # absorb it.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ---- Observation and error -------------------------------------
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
