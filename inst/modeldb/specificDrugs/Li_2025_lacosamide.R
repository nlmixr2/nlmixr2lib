Li_2025_lacosamide <- function() {
  description <- "One-compartment population PK model with first-order absorption for oral lacosamide in Chinese children with epilepsy, with a body-weight power function on apparent clearance; Model I of Li 2025, for the clinical scenario in which CYP2C19 genotype is unavailable"
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
        "median of 30 kg: CL/F = 1.51 * (WT/30)^0.294 (Results Eq. 1; Table 2",
        "'Covariate model structure'). The 30 kg reference is the median body",
        "weight of the model-development group (Table 1: BW 30 kg, range",
        "10-80 kg), which is what Supplementary Eq. 13 defines as COVmedian.",
        "Body weight does NOT enter V/F in the final model: the V-BW allometric",
        "term was screened and rejected for parameter instability (RSE 327.4%",
        "versus < 15% for the base model; Results, paragraph following Eq. 4),",
        "so V/F is a single population value of 23.7 L for every subject.",
        "Body weight was recorded at each therapeutic-drug-monitoring visit, so",
        "it is time-varying in principle, but the model was fit to steady-state",
        "trough samples and the paper reports a single baseline weight per",
        "subject; supply the weight in effect at the time of the dose.",
        sep = " "
      ),
      source_name        = "BW"
    )
  )

  # Covariates that Li 2025 screened but did NOT retain in Model I. Documented
  # here for provenance only; none is referenced in model().
  covariatesDataExcluded <- list(
    SNP_CYP2C19_RS12769205_GA = list(
      description        = "CYP2C19*2 rs12769205 heterozygous GA genotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (AA genotype), jointly with SNP_CYP2C19_RS12769205_GG = 0",
      notes              = paste(
        "DELIBERATELY excluded from Model I. Li 2025 built two final models for",
        "the two real-world clinical scenarios defined by genotype availability",
        "(Methods, 'Simulation of dosing regimen'): Model I uses body weight",
        "only, and Model II adds CYP2C19*2 rs12769205. This file is Model I; the",
        "genotype-informed sibling is Li_2025_lacosamide_cyp2c19.R. rs12769205",
        "was a highly significant predictor of CL/F when tested (dOFV 16.465,",
        "p < 0.001; Supplementary Table S5), so its absence here reflects the",
        "intended use case, not a failed covariate screen.",
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
        "DELIBERATELY excluded from Model I; see the",
        "SNP_CYP2C19_RS12769205_GA notes and Li_2025_lacosamide_cyp2c19.R.",
        sep = " "
      ),
      source_name        = "rs12769205"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Screened on CL/F with a proportional model (Supplementary Eq. 16) and",
        "rejected at forward inclusion: dOFV -0.173, p > 0.05 (Supplementary",
        "Table S4). Cohort was 52 of 133 female (Table 1).",
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
        "Listed among the demographics offered to covariate screening (Methods,",
        "'Covariate model') and not retained on any parameter. Median 7.5 years,",
        "range 1-18 years (Table 1). Age enters the model only indirectly,",
        "through body weight; the Discussion attributes the positive BW-CL/F",
        "relationship partly to age-related maturation of renal function. Two",
        "explicit maturation parameterizations (a BW simple-exponent model and a",
        "fixed 0.75 allometric model, Supplementary Eqs. 8-11) were fitted and",
        "both were abandoned for excessive RSE and shrinkage (> 80%;",
        "Supplementary S.4.4), so the final model carries no maturation term.",
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
        "Three separate eGFR formulations were screened on CL/F -- Shull's",
        "creatinine equation and two cystatin-C equations (Supplementary",
        "Eqs. 18-20; Table 1 rows eGFR a/b/c) -- and all three were rejected at",
        "forward inclusion (dOFV -0.134, -1.087 and -0.123; Supplementary",
        "Table S4). Supplementary S.4.4 states outright that 'renal function did",
        "not affect CL/F of LCM'. Serum creatinine itself did reach the forward",
        "inclusion threshold (dOFV -16.341) but was dropped at backward",
        "elimination.",
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
        "Results 'Model development'; Supplementary Table S4). Representative of",
        "that group; the full screened panel (about 30 demographic, comedication,",
        "hematology and biochemistry covariates) is tabulated in Supplementary",
        "Table S4 and listed in the vignette rather than enumerated here.",
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
    regions          = "China (single centre, Children's Hospital of Nanjing Medical University)",
    notes            = paste(
      "Retrospective real-world therapeutic-drug-monitoring cohort collected",
      "between June 2021 and March 2023. 190 children contributing 493",
      "concentrations were randomly split 70/30 into a model-development group",
      "(133 children, 347 concentrations -- the numbers recorded above and the",
      "data this model was fitted to) and an external validation group (57",
      "children, 146 concentrations). Baseline demographics are Table 1;",
      "additional screened characteristics are Supplementary Table S1. All",
      "samples are pre-dose troughs drawn 30 min before the next maintenance",
      "dose after at least 3 days on an unchanged regimen (Methods, 'Sample",
      "collection and concentration measurement'), which is why absorption and",
      "distribution are only weakly identified: ka is fixed and V/F carries no",
      "inter-individual variability. Race is recorded as 100 percent Asian",
      "because the cohort is a Chinese single-centre pediatric population; the",
      "paper reports no formal race or ethnicity breakdown. Concomitant",
      "antiseizure medications were present (valproate 19.5 percent,",
      "levetiracetam 21.1 percent, sodium channel blockers 9.7 percent) but none",
      "was retained as a covariate, and the Discussion cautions against applying",
      "the model to children on enzyme-inducing comedication or sodium channel",
      "blockers.",
      sep = " "
    )
  )

  ini({
    # ================================================================
    # Structural parameters -- Li 2025 Table 2, "Model I" block,
    # Estimate column. A one-compartment model with first-order
    # absorption and elimination, implemented in NONMEM as ADVAN2
    # TRANS2 (Methods, 'Base model').
    #
    # CL/F and V/F are APPARENT values: Methods 'Base model' states
    # that "as bioavailability (F) could not be determined, CL and Vd
    # were expressed as apparent values".
    #
    # UNITS OF CL/F. Table 2's row label is "CL/F (L/h)" and the
    # Discussion says "Model I estimated the CL/F of 1.51 L/h", but the
    # Results sentence introducing Eq. 1 misprints it as "1.51 L/h/kg".
    # Table 2 and the Discussion agree on L/h, and the per-kg reading is
    # arithmetically impossible: at 1.51 L/h/kg a 30 kg child would
    # clear 45.3 L/h and the steady-state trough on this cohort's median
    # regimen would be about 0.12 mg/L, 30-fold below every observed
    # concentration (Table 1 range 0.70-11.90 mg/L). L/h is used here;
    # see the vignette Errata.
    # ================================================================
    lka <- fixed(log(2.45)); label("Absorption rate constant ka (1/h)")                     # Table 2 Model I: "Ka 2.45 (fixed)". Methods 'Base model': "due to limited sampling data during the absorption phase, the ka could not be reliably estimated. Therefore, ka was fixed at 2.45 h-1, based on values reported in a previous study" (the paper's reference 15). No RSE is reported, consistent with a fixed value.
    lcl <- log(1.51);        label("Apparent clearance CL/F at the reference body weight of 30 kg (L/h)") # Table 2 Model I: CL/F = 1.51 L/h, RSE 3.4% (bootstrap median 1.51, 95% CI 1.36-1.67). Also the leading coefficient of Results Eq. 1.
    lvc <- log(23.7);        label("Apparent volume of distribution V/F (L)")                # Table 2 Model I: V/F = 23.7 L, RSE 12.1% (bootstrap median 23.6, 95% CI 16.7-38.1). Also Results Eq. 2, which carries no covariate term.

    # ================================================================
    # Covariate effect on CL/F -- body weight, power function
    # normalized to the cohort median of 30 kg.
    # Supplementary Eq. 13 defines the power form as
    #   Pi = TV(P) * (COV / COVmedian)^theta
    # and Table 2's "Covariate model structure" footnote prints the
    # instantiated relationship as CL/F = 1.51 * (BW/30)^0.294, which
    # is Results Eq. 1.
    # Sign check: the Discussion states "a significant positive
    # correlation between BW and CL/F was observed, indicating that
    # higher BW necessitates increased LCM doses". Positive exponent.
    # Consistent.
    # ================================================================
    e_wt_cl <- 0.294; label("Power exponent on (WT/30 kg) for CL/F (unitless)")              # Table 2 Model I: "CL-BW" = 0.294, RSE 14.1% (bootstrap median 0.296, 95% CI 0.207-0.380). Also the exponent printed in Results Eq. 1.

    # ================================================================
    # Inter-individual variability -- Li 2025 Table 2, "IIV in CL/F,
    # %CV" row. Exponential (log-normal) IIV on CL/F only.
    #
    # WHICH IIV MODEL. Supplementary S.4.4 prose says the IIV "on CL and
    # V were estimated using exponential models (Eq. 3)", but its own
    # Eq. 3 is the PROPORTIONAL form Pi = TV(P)*(1 + eta) while Eq. 1 is
    # the exponential form Pi = TV(P)*exp(eta). Supplementary Tables S4
    # and S5 both label the base model "Eq.1 & Eq.7", i.e. exponential
    # IIV with a mixed residual error, and the main text (Results,
    # 'Model development') independently says "exponential models". The
    # exponential form is used here; the "(Eq. 3)" cross-reference is a
    # typo.
    #
    # WHY NO IIV ON V/F. The main text says IIV was estimated on both
    # CL/F and V/F, but Table 2 reports an IIV row for CL/F only, and
    # Supplementary S.4.4 resolves it: "Due to only the C0 samples of
    # LCM being available, the IIV for Vd/F was not informative enough
    # to be estimated and then was excluded from the model." Trough-only
    # sampling. No eta on V/F here, and none on the fixed ka.
    #
    # %CV -> VARIANCE CONVENTION. Table 2 reports IIV as a percent CV
    # and does not print omega^2. The same table's residual rows fix the
    # convention: the reported "Proportional error, %CV = 17.4" and
    # "Additive error, SD = 0.554" correspond to the variances of 0.0301
    # and 0.307 quoted in Results 'Model evaluation', and 0.554^2 =
    # 0.3069 -> 0.307 while 0.174^2 = 0.0303 -> 0.0301. Both are the
    # plain square, so the table's %CV entries are sqrt(variance)*100
    # and omega^2 = 0.180^2. The exact log-normal alternative,
    # omega^2 = log(1 + 0.180^2) = 0.0319, changes omega by 0.8% and no
    # validation gate in the vignette; see the vignette Errata.
    # ================================================================
    etalcl ~ 0.0324  # Table 2 Model I: IIV in CL/F = 18.0 %CV, RSE 19.1% (bootstrap median 17.8, 95% CI 14.1-21.4); eta-shrinkage 18.2%. Variance = 0.180^2 = 0.0324.

    # ================================================================
    # Residual error -- Li 2025 Table 2, "Proportional error" and
    # "Additive error" rows. The mixed error model is Supplementary
    # Eq. 7, Y = IPRED*(1 + eps1) + eps2, selected over the additive,
    # exponential and proportional alternatives (Eqs. 4-6) for the
    # lowest OFV and AIC (Results 'Model development'; Supplementary
    # S.4.4). Eq. 7 gives Var(Y) = (sigma1*IPRED)^2 + sigma2^2, which is
    # exactly nlmixr2's prop(propSd) + add(addSd) combination.
    #
    # Concentrations are reported in ug/mL throughout the paper, which
    # is numerically identical to the mg/L declared in units above, so
    # the additive SD transfers unchanged.
    # ================================================================
    propSd <- 0.174;        label("Proportional residual error (fraction)")                  # Table 2 Model I: "Proportional error, %CV" = 17.4, RSE 20.7% (bootstrap median 17.3, 95% CI 13.8-20.7); sigma-shrinkage 14.1%. Results 'Model evaluation' quotes the corresponding variance as 0.0301 (= 0.174^2 to rounding).
    addSd  <- fixed(0.554); label("Additive residual error SD on Cc (mg/L)")                 # Table 2 Model I: "Additive error, SD 0.554 (fixed)". No RSE reported, consistent with a fixed value. Results 'Model evaluation' quotes the corresponding variance as 0.307 (= 0.554^2).
  })

  model({
    # ---- Individual PK parameters ---------------------------------
    # Body weight scales CL/F as a power function normalized to the
    # cohort median of 30 kg (Results Eq. 1). V/F carries no covariate
    # and no eta (Results Eq. 2); ka is fixed.
    ka  <- exp(lka)
    cl  <- exp(lcl + etalcl) * (WT / 30)^e_wt_cl
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
