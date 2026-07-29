# One-compartment parent-plus-metabolite population PK model for oral
# clozapine and its primary active metabolite norclozapine
# (N-desmethylclozapine) in Chinese adult inpatients with refractory
# schizophrenia (Li 2012, Acta Pharmacol Sin 33:1409-1416;
# doi:10.1038/aps.2012.71).

Li_2012_clozapine <- function() {
  description <- paste(
    "One-compartment parent-plus-metabolite population PK model for oral",
    "clozapine and its primary active metabolite norclozapine",
    "(N-desmethylclozapine) in 162 Chinese adult inpatients (74 male, 88",
    "female; 35.5 +/- 10.6 years) with refractory schizophrenia on",
    "maintenance oral clozapine therapy (Li 2012). First-order absorption",
    "(Ka fixed at 1.3 1/h from prior rich-data clozapine PK studies) into a",
    "single central compartment with first-order elimination; a fixed",
    "fraction (KF = 0.66) of the absorbed clozapine dose is converted in",
    "the parent central compartment to norclozapine and feeds a separate",
    "one-compartment metabolite compartment with its own apparent",
    "clearance and apparent volume. Two binary covariates were retained in",
    "the final forward-and-backward-selected model: current-smoker status",
    "increases apparent clearance of both species (clozapine by 45%,",
    "norclozapine by 54.3%), and male sex increases apparent clearance of",
    "both species (clozapine by 20.8%, norclozapine by 24.2%); the typical",
    "values reported in Table 2 are for the female-nonsmoker reference",
    "stratum. A combined additive-plus-proportional residual error model",
    "is reported separately for clozapine and norclozapine. The model was",
    "internally validated using normalized prediction distribution errors",
    "(NPDE).",
    sep = " "
  )
  reference <- paste(
    "Li LJ, Shang DW, Li WB, Guo W, Wang XP, Ren YP, Li AN, Fu PX,",
    "Ji SM, Lu W, Wang CY (2012). Population pharmacokinetics of",
    "clozapine and its primary metabolite norclozapine in Chinese",
    "patients with schizophrenia. Acta Pharmacol Sin 33(11):1409-1416.",
    "doi:10.1038/aps.2012.71.",
    sep = " "
  )
  vignette <- "Li_2012_clozapine"
  units    <- list(time = "hour", dosing = "umol", concentration = "umol/L")

  covariateData <- list(
    SEXF = list(
      description        = "Biological sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male) in the canonical column. The Li 2012 paper's reference category is female (Table 2 typical values are for female nonsmokers; Methods example formula uses SEX = 0 for female, 1 for male).",
      notes              = paste(
        "Li 2012 encodes sex as a male-indicator (1 = male, 0 = female) with",
        "female as the reference category (Methods 'For male subjects, the",
        "theta_cov term is added to the population estimate'). To store under",
        "the canonical SEXF (1 = female, 0 = male) while preserving Li's",
        "female-reference CL/F and CLM typical values, the effect is applied",
        "in model() as (1 + e_sex_cl * (1 - SEXF)) and (1 + e_sex_cl_norcloz *",
        "(1 - SEXF)), so SEXF = 1 (female) yields factor 1 and SEXF = 0 (male)",
        "yields the paper's +20.8% (parent) / +24.2% (metabolite) increment."
      ),
      source_name        = "SEX (1 = male, 0 = female)"
    ),
    SMOKE = list(
      description        = "Current-smoker status indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-smoker).",
      notes              = paste(
        "Li 2012 Methods: 'Persons who had smoked 5 or more cigarettes per",
        "day within the last week were defined as smokers.' Self-reported and",
        "nurse-confirmed; no biochemical verification (Discussion lists this",
        "as a limitation). Distribution in the 162-subject cohort (Table 1):",
        "72 female nonsmokers, 2 female smokers, 40 male nonsmokers, 48 male",
        "smokers. Applied multiplicatively on apparent clearance of both",
        "clozapine and norclozapine."
      ),
      source_name        = "smoking status"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age in years",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened as a continuous covariate on apparent clearance and volume",
        "of distribution of both clozapine and norclozapine via the",
        "power-law form TVCL = theta_CL * (AGE/AGEAVE)^theta_AGE (Methods).",
        "Not retained in the final model: 'Other covariates such as weight",
        "and age did not significantly influence the PK parameters of",
        "clozapine and norclozapine' (Results); the chi-square OFV threshold",
        "for forward inclusion was P < 0.01 (df = 1). Documented here to",
        "preserve the FREM-style covariate screen without carrying a",
        "convention warning for a declared-but-not-used covariate."
      )
    ),
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened alongside AGE as a continuous covariate on apparent",
        "clearance and volume of distribution; not retained in the final",
        "model (Results). Cohort range and central tendency are not",
        "tabulated separately by Li 2012; Table 1 reports only age,",
        "sex / smoking strata, and concentration summaries."
      )
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 162L,
    n_observations  = 1617L,
    n_studies       = 1L,
    age_range       = "18-59 years (mean 35.5, SD 10.6; Table 1)",
    age_median      = "mean 35.5 years (median not reported)",
    weight_range    = "not tabulated in Table 1; WT and BSA were collected but only Methods narrates the screen",
    weight_median   = "not reported",
    sex_female_pct  = 54.3,
    race_ethnicity  = c(Asian = 100),
    disease_state   = paste(
      "DSM-IV-diagnosed refractory schizophrenia, hospitalized at multiple",
      "mental-health sites in China. All patients on maintenance oral",
      "clozapine therapy (Jiangsu Nhwa pharmaceutical, Xuzhou); most",
      "dosed twice or three times per day. Compliance was confirmed by",
      "repeated serum-clozapine measurements during the study."
    ),
    dose_range      = paste(
      "Patient-clinician-prescribed oral clozapine; specific per-subject",
      "daily-dose distribution is not tabulated in the publication.",
      "Sampling: 20% of total samples taken 0.5-4.5 h after the last dose,",
      "remaining samples taken 8.5-15.5 h after the last dose (sparse",
      "design typical of therapeutic-drug-monitoring data)."
    ),
    regions         = "China (multiple mental-health sites; lead site Beijing An Ding Hospital)",
    smoke_strata    = "female nonsmoker n=72, female smoker n=2, male nonsmoker n=40, male smoker n=48 (Table 1)",
    notes           = paste(
      "Refractory schizophrenia, Chinese-ethnicity inpatient cohort.",
      "Clozapine and norclozapine concentrations were quantified by HPLC",
      "with UV detection at 254 nm using desipramine as internal standard;",
      "lower limit of detection 0.08 umol/L for both species; CV < 5%",
      "(Methods 'Determination of clozapine and norclozapine concentrations').",
      "Population PK fit with NONMEM 7.1 ADVAN5; FOCE without interaction",
      "was selected over FOCE-I after both were tested. Internal validation",
      "used the NPDE method (R package 'npde'). Concomitant medications",
      "were recorded but not formally analyzed for drug-drug interactions",
      "(Discussion limitations)."
    )
  )

  ini({
    # ---- Structural typical values (Li 2012 Table 2) ----
    # All clozapine and norclozapine typical values in Table 2 are the
    # estimates for the female-nonsmoker reference stratum
    # (SEXF = 1 / SMOKE = 0); covariate adjustments are applied
    # multiplicatively in model() so a male and / or a smoker is driven
    # away from these reference values.

    # Absorption rate constant. Fixed at 1.3 1/h based on prior
    # rich-data clozapine PK studies (Methods Results paragraph 1:
    # "The first-order absorption rate constant (Ka) for clozapine was
    # fixed at 1.3 h-1 based on several pharmacokinetic studies that
    # obtained rich data describing the pharmacokinetics of clozapine
    # in patients [ref 12]").
    lka <- fixed(log(1.3)); label("Absorption rate constant Ka (1/h); fixed per Methods Results paragraph 1")  # Li 2012 Table 2: Ka = 1.3 (Fixed)

    # Apparent oral clearance of clozapine (CL/F) and apparent volume of
    # distribution of clozapine (V/F) at the female-nonsmoker reference.
    lcl <- log(21.9);  label("Apparent oral clearance of clozapine CL/F at the female-nonsmoker reference (L/h)")  # Li 2012 Table 2: CL/F = 21.9 L/h (RSE 6%)
    lvc <- log(526);   label("Apparent volume of distribution of clozapine V/F (L)")                              # Li 2012 Table 2: V/F = 526 L (RSE 10%)

    # Fixed fraction of absorbed clozapine converted to norclozapine in
    # the parent central compartment (Methods Results paragraph 1: "The
    # fraction of the absorbed dose of clozapine converted into
    # norclozapine (KF) was fixed at 0.66 in published papers and was
    # validated by the ratio of the mean amount of norclozapine to the
    # mean amount of clozapine at steady-state in these articles
    # [refs 22, 23]"). Linear-scale (in [0, 1]), wrapped in fixed().
    kf  <- fixed(0.66); label("Fraction of absorbed clozapine converted to norclozapine (unitless, in [0, 1]); fixed")  # Li 2012 Table 2: KF = 0.66 (Fixed)

    # Apparent oral clearance of norclozapine (CLM) and apparent volume
    # of distribution of norclozapine (VM) at the female-nonsmoker
    # reference. CLM and VM are "apparent" in the same parent-bioavailability
    # sense because norclozapine is not directly dosed; the values absorb
    # any norclozapine-formation bioavailability that is not separately
    # identifiable from the data.
    lcl_norcloz <- log(32.7); label("Apparent clearance of norclozapine CLM at the female-nonsmoker reference (L/h)")  # Li 2012 Table 2: CLM = 32.7 L/h (RSE 5.6%)
    lvc_norcloz <- log(624);  label("Apparent volume of distribution of norclozapine VM (L)")                          # Li 2012 Table 2: VM = 624 L (RSE 5.5%)

    # ---- Covariate effects ----
    # Linear multiplicative form per Methods Eq. for covariate model:
    #   TVCL = theta_CL * (1 + theta_smoke * SMOKE + theta_sex * MALE)
    # where MALE = (1 - SEXF). The paper's reported coefficients (0.45,
    # 0.208, 0.543, 0.242) match the percent increases stated in
    # Results ("Smoking was associated with increases in the clearance
    # of clozapine and norclozapine of 45% and 54.3%, respectively. The
    # clearance of clozapine and norclozapine were 20.8% and 24.2%
    # greater, respectively, in males than in females"), confirming the
    # multiplicative-fractional parameterization.

    e_smoke_cl         <- 0.45;  label("Smoking effect on clozapine CL/F (fraction; applied as (1 + e_smoke_cl * SMOKE))")          # Li 2012 Table 2: theta_Smoking (cloz) = 0.45 (RSE 34.9%)
    e_sex_cl           <- 0.208; label("Male-sex effect on clozapine CL/F (fraction; applied as (1 + e_sex_cl * (1 - SEXF)))")      # Li 2012 Table 2: theta_Gender (cloz) = 0.208 (RSE 44.6%)
    e_smoke_cl_norcloz <- 0.543; label("Smoking effect on norclozapine CLM (fraction; applied as (1 + e_smoke_cl_norcloz * SMOKE))")          # Li 2012 Table 2: theta_Smoking (norcloz) = 0.543 (RSE 35.7%)
    e_sex_cl_norcloz   <- 0.242; label("Male-sex effect on norclozapine CLM (fraction; applied as (1 + e_sex_cl_norcloz * (1 - SEXF)))")      # Li 2012 Table 2: theta_Gender (norcloz) = 0.242 (RSE 49.2%)

    # ---- Inter-individual variability (Li 2012 Table 2) ----
    # The paper reports IIV in the "Interindividual variability%" column
    # as a CV%. For log-normal IIV the relationship between CV and the
    # NONMEM OMEGA variance is omega^2 = log(CV^2 + 1); the values below
    # use that conversion.
    #   CL/F:  CV = 42.9% -> omega^2 = log(0.429^2 + 1) = 0.16898
    #   V/F:   CV = 65.7% -> omega^2 = log(0.657^2 + 1) = 0.35892
    #   CLM:   CV = 42.1% -> omega^2 = log(0.421^2 + 1) = 0.16325
    #   VM:    CV = 75.6% -> omega^2 = log(0.756^2 + 1) = 0.45235
    # IIV is reported only on CL/F, V/F, CLM, and VM (no IIV on Ka,
    # which is fixed, and no IIV on KF, which is also fixed).
    etalcl         ~ 0.16898  # Li 2012 Table 2: CL/F  IIV CV% = 42.9 -> omega^2 = log(0.429^2 + 1)
    etalvc         ~ 0.35892  # Li 2012 Table 2: V/F   IIV CV% = 65.7 -> omega^2 = log(0.657^2 + 1)
    etalcl_norcloz ~ 0.16325  # Li 2012 Table 2: CLM   IIV CV% = 42.1 -> omega^2 = log(0.421^2 + 1)
    etalvc_norcloz ~ 0.45235  # Li 2012 Table 2: VM    IIV CV% = 75.6 -> omega^2 = log(0.756^2 + 1)

    # ---- Residual error (Li 2012 Table 2 Residual Variability) ----
    # Combined additive + proportional error model fitted separately for
    # clozapine and norclozapine (Methods "Combined additive and
    # proportional error model" and Table 2). Reported as:
    #   sigma_1 (clozapine additive)        SD = 0.162 umol/L
    #   sigma_2 (clozapine proportional)    CV = 26.6%
    #   sigma_3 (norclozapine additive)     SD = 0.117 umol/L
    #   sigma_4 (norclozapine proportional) CV = 16.9%
    # In nlmixr2 propSd / addSd are on the SD scale; the values below
    # are taken directly from Table 2.
    addSd         <- 0.162; label("Additive residual error for clozapine (umol/L)")          # Li 2012 Table 2: sigma_1 SD = 0.162 umol/L
    propSd        <- 0.266; label("Proportional residual error for clozapine (fraction)")    # Li 2012 Table 2: sigma_2 CV = 26.6%
    addSd_norcloz <- 0.117; label("Additive residual error for norclozapine (umol/L)")        # Li 2012 Table 2: sigma_3 SD = 0.117 umol/L
    propSd_norcloz <- 0.169; label("Proportional residual error for norclozapine (fraction)") # Li 2012 Table 2: sigma_4 CV = 16.9%
  })

  model({
    # Derived sex term. Li 2012 encodes sex as a male-indicator (1 =
    # male, 0 = female) with female as the reference category; (1 -
    # SEXF) reproduces the paper's male = 1 column while keeping SEXF
    # (1 = female) as the canonical storage convention. See
    # covariateData[[SEXF]]$notes.
    sex_male <- 1 - SEXF

    # Apparent oral clearance of clozapine with smoking and male-sex
    # multiplicative effects on the female-nonsmoker reference TVCL,
    # plus log-normal IIV on the linear-scale typical value.
    cl <- exp(lcl + etalcl) *
      (1 + e_smoke_cl * SMOKE + e_sex_cl * sex_male)
    vc <- exp(lvc + etalvc)
    ka <- exp(lka)

    # Apparent clearance and apparent volume of distribution of
    # norclozapine with smoking and male-sex effects on CLM.
    cl_norcloz <- exp(lcl_norcloz + etalcl_norcloz) *
      (1 + e_smoke_cl_norcloz * SMOKE + e_sex_cl_norcloz * sex_male)
    vc_norcloz <- exp(lvc_norcloz + etalvc_norcloz)

    # Disposition. d/dt(depot): first-order oral absorption of
    # clozapine. d/dt(central): clozapine in the central compartment,
    # eliminated at total rate (cl / vc); a fraction KF of that
    # elimination flux feeds the norclozapine central compartment.
    # d/dt(central_norcloz): norclozapine accrued at rate kf * (cl /
    # vc) * central and eliminated at rate (cl_norcloz / vc_norcloz).
    # Mass conversion assumes 1:1 molar stoichiometry (consistent with
    # the molar dosing / observation unit pair umol / umol_per_L);
    # clozapine MW 326.83 g/mol and norclozapine MW 312.80 g/mol
    # differ by ~4%, so a mass-unit fit would carry the same KF with
    # only a 4% bias absorbed into the apparent CLM / VM estimates.
    d/dt(depot)           <- -ka * depot
    d/dt(central)         <-  ka * depot - (cl / vc) * central
    d/dt(central_norcloz) <-  kf * (cl / vc) * central - (cl_norcloz / vc_norcloz) * central_norcloz

    # Plasma concentrations. Both observation variables are in umol/L
    # to match the source paper's reporting (Methods 'Determination of
    # clozapine and norclozapine concentrations' and Table 1 mean
    # concentrations).
    Cc         <- central         / vc
    Cc_norcloz <- central_norcloz / vc_norcloz

    # Combined additive + proportional residual error model, fitted
    # separately for parent and metabolite (Methods 'Combined additive
    # and proportional error model' and Table 2 Residual Variability).
    Cc         ~ add(addSd)         + prop(propSd)
    Cc_norcloz ~ add(addSd_norcloz) + prop(propSd_norcloz)
  })
}
