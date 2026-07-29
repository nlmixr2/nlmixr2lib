# Adult one-compartment mixture-model population PK of oral risperidone and
# its active metabolite 9-OH-risperidone in the Clinical Antipsychotic
# Trials of Intervention Effectiveness (CATIE) cohort, with three CYP2D6
# metabolizer subpopulations (PM / IM / EM) and a power age effect on
# metabolite clearance (Feng 2008, Br J Clin Pharmacol 66(5):629-639;
# doi:10.1111/j.1365-2125.2008.03276.x).

Feng_2008_risperidone <- function() {
  description <- paste(
    "Adult one-compartment parent-plus-metabolite population PK model for",
    "oral risperidone and its active metabolite 9-OH-risperidone in 490",
    "subjects pooled across the CATIE-AD (n = 110, behavioural symptoms of",
    "Alzheimer disease, mean age 78.3 years) and CATIE-SZ (n = 380,",
    "schizophrenia, mean age 40.6 years) trials (Feng 2008). First-order",
    "absorption with Ka fixed at 1.7 1/h into a single central compartment",
    "with first-order elimination; the fraction of risperidone metabolized",
    "to 9-OH-risperidone (KF) feeds a single metabolite compartment whose",
    "apparent volume of distribution is set equal to the parent apparent",
    "volume per the paper's identifiability constraint. A mixture model",
    "with three CYP2D6 metabolizer subpopulations (poor PM, intermediate",
    "IM, extensive EM) yields subpopulation-specific apparent oral",
    "clearances (CL/F) and metabolite formation fractions (KF); CL/F in IM",
    "(36 L/h) and KF in IM (1) are fixed per Table 3 to stabilize the",
    "mixture estimation. Age is the only retained subject-level covariate,",
    "acting on 9-OH-risperidone apparent clearance (CLM/F) via a power",
    "model with exponent -0.378 referenced at a nominal median age of 45",
    "years. Inter-individual variability is reported separately for CL/F",
    "in PM and EM (no IIV is reported for CL/F in IM or for CLM/F), for Ka",
    "(despite a fixed typical value), and for the shared Vd/F; combined",
    "additive-plus-proportional residual error is reported separately for",
    "risperidone and 9-OH-risperidone plasma concentrations.",
    sep = " "
  )
  reference <- paste(
    "Feng Y, Pollock BG, Coley K, Marder S, Miller D, Kirshner M,",
    "Aravagiri M, Schneider L, Bies RR (2008). Population pharmacokinetic",
    "analysis for risperidone using highly sparse sampling measurements",
    "from the CATIE study. British Journal of Clinical Pharmacology",
    "66(5):629-639. doi:10.1111/j.1365-2125.2008.03276.x.",
    sep = " "
  )
  vignette <- "Feng_2008_risperidone"
  units    <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Combined cohort age 18-93 years (mean 49.1, SD 18.8; Table 1).",
        "CATIE-AD subjects 57-93 years (mean 78.3, SD 6.7); CATIE-SZ",
        "subjects 18-65 years (mean 40.6, SD 11.2). Reference age 45 years",
        "is a nominal round value approximating the median of the combined",
        "cohort per Methods 'Final model development' ('the reference value",
        "for a covariate was specified to be a nominal value that",
        "approximates the median for the covariate'). Age enters CLM/F via",
        "a power model with exponent -0.378 (Table 3). The Discussion",
        "reports POSTHOC CLM averages of 6.1 L/h at age 45 and 4.9 L/h at",
        "age 70 that are not exactly reproducible from the Table 3 typical",
        "value (8.83 L/h) and exponent (-0.378) with reference 45; see the",
        "vignette Assumptions and deviations section."
      ),
      source_name        = "AGE"
    ),
    CYP2D6_PM = list(
      description        = "CYP2D6 poor-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (intermediate or extensive metabolizer; both CYP2D6_PM and CYP2D6_EM = 0 indicates IM)",
      notes              = paste(
        "1 = subject is a CYP2D6 poor metabolizer, 0 otherwise. Paired with",
        "CYP2D6_EM to encode the three-level PM / IM / EM phenotype; IM is",
        "the implicit reference (both CYP2D6_PM = 0 and CYP2D6_EM = 0).",
        "Feng 2008 inferred CYP2D6 phenotype from the NONMEM mixture-model",
        "posterior rather than from external genotyping; the final-model",
        "mixture-fraction estimates are 41.2% PM (Table 3 P1), 52.4% EM",
        "(Table 3 P2), and 6.4% IM (1 - P1 - P2). To simulate the cohort,",
        "draw the phenotype jointly from the (P1, P2, 1 - P1 - P2)",
        "categorical distribution. Same encoding convention as",
        "Sherwin_2012_risperidone.R."
      ),
      source_name        = "P1 (mixture-model PM subpopulation fraction); inferred per subject from NONMEM mixture posterior"
    ),
    CYP2D6_EM = list(
      description        = "CYP2D6 extensive-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (intermediate or poor metabolizer; both CYP2D6_PM and CYP2D6_EM = 0 indicates IM)",
      notes              = paste(
        "1 = subject is a CYP2D6 extensive metabolizer, 0 otherwise. Paired",
        "with CYP2D6_PM; IM is the implicit reference (both indicators = 0).",
        "In Feng 2008 the mixture-model assignment estimated 41.2% PM",
        "(Table 3 P1), 52.4% EM (Table 3 P2), and 6.4% IM (1 - P1 - P2)."
      ),
      source_name        = "P2 (mixture-model EM subpopulation fraction); inferred per subject from NONMEM mixture posterior"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 490L,
    n_observations = 2472L,
    n_studies      = 2L,
    age_range      = "18-93 years (mean 49.1, SD 18.8; combined median approximately 45)",
    age_median     = "approximately 45 years (mean 49.1 reported; CATIE-SZ dominates the lower half of the combined sample)",
    weight_range   = "42.7-187.7 kg (mean 84.1, SD 22.5)",
    weight_median  = "84.1 kg (mean reported, not median)",
    sex_female_pct = 32.4,
    race_ethnicity = c(
      White                     = 66.9,
      `Black or African-American` = 28.6,
      Asian                     = 2.4,
      `American Indian`         = 1.0,
      `Two or more races`       = 0.8,
      `Native Hawaiian`         = 0.2
    ),
    disease_state  = paste(
      "Pooled cohort of two CATIE substudies: 110 outpatients with",
      "Alzheimer disease and behavioural disturbance (CATIE-AD, mean age",
      "78.3 years) and 380 patients with schizophrenia aged 18-65 years",
      "(CATIE-SZ). All received oral risperidone as part of a multi-arm",
      "trial of atypical antipsychotics."
    ),
    dose_range     = paste(
      "Oral risperidone tablet, total daily dose 0.5-6.0 mg (CATIE-AD",
      "0.5-3.5 mg/day; CATIE-SZ 0.75-6.0 mg/day). 313 subjects dosed once",
      "daily; 177 subjects dosed twice daily (Results)."
    ),
    regions        = "USA (multicentre CATIE network)",
    cyp2d6_distribution = paste(
      "Mixture-model assignment in the final model (Table 3): 41.2% PM",
      "(P1), 52.4% EM (P2), 6.4% IM (1 - P1 - P2). Phenotype was inferred",
      "from the NONMEM mixture posterior, not from external CYP2D6",
      "genotyping. The Discussion notes the unusually high PM fraction",
      "(~41%) versus the ~5-10% PM fraction expected from CYP2D6 allele",
      "frequencies and attributes the difference partly to concomitant",
      "CYP2D6 inhibitors (paroxetine / fluoxetine) re-classifying subjects",
      "into the PM stratum and partly to variable medication adherence."
    ),
    covariates_screened = paste(
      "Age, weight, sex, smoking status, race, and concomitant",
      "fluoxetine / paroxetine were tested as covariates on CL/F and",
      "Vd/F (Table 2). After the three-component mixture model was",
      "introduced, none of these were retained on parent risperidone CL/F",
      "or on Vd/F. Only age on 9-OH-risperidone CLM/F was retained in",
      "the final model (Methods 'Final model development'; Table 2",
      "section 3-3 model M26 vs M25, dOFV = -72.2, p < 0.005)."
    ),
    notes          = paste(
      "Random plasma samples (1-6 per subject) collected over the CATIE",
      "follow-up window: CATIE-AD samples at weeks 2, 4, 12 or at",
      "medication-switch points; CATIE-SZ random samples every 3 months",
      "for up to 6 samples per subject (Methods). Risperidone and",
      "9-OH-risperidone quantified by LC-MS/MS with a 0.1 ng/mL limit of",
      "detection. NONMEM v5 level 1.1 was used with the FO estimation",
      "method; FOCEI was unstable on this sparse data and produced biased",
      "predictions (Methods). A two-component mixture (PM + EM only) was",
      "also tested and rejected (dOFV +186.6 versus three-component, AIC",
      "p < 0.005)."
    )
  )

  ini({
    # Structural absorption. Ka is fixed at 1.7 1/h per Table 3 (the
    # Ka SE column shows 'NA'); inter-individual variability is
    # estimated even though the typical value is fixed
    # (Table 3 w_ka% = 53.7, SE 89.3%).
    lka <- fixed(log(1.7)); label("Absorption rate constant Ka (1/h); typical value fixed per Table 3")  # Feng 2008 Table 3: K_a = 1.7 (Fixed)

    # Shared apparent volume of distribution. The paper sets VM = V
    # for the metabolite due to identifiability (Results 'Base model':
    # "Due to the identifiability problem associated with KF and V_M,
    # V_M was set to the same value as V").
    lvc <- log(444); label("Apparent central volume of distribution Vd/F (L); shared with metabolite VdM/F")  # Feng 2008 Table 3: V, V_M = 444 L (SE 17.8%)

    # Subpopulation-specific apparent oral clearances of risperidone.
    # Estimated for PM and EM; CL/F in IM is fixed at 36 L/h per
    # Table 3 (SE column 'NA') to stabilize the three-component
    # mixture estimation.
    lcl_pm <- log(12.9); label("Apparent oral clearance in CYP2D6 PMs, CL/F (L/h)")  # Feng 2008 Table 3: CL in PM = 12.9 L/h (SE 6.5%)
    lcl_em <- log(65.4); label("Apparent oral clearance in CYP2D6 EMs, CL/F (L/h)")  # Feng 2008 Table 3: CL in EM = 65.4 L/h (SE 9.9%)
    lcl_im <- fixed(log(36)); label("Apparent oral clearance in CYP2D6 IMs, CL/F (L/h); fixed")  # Feng 2008 Table 3: CL in IM = 36 L/h (Fixed)

    # Metabolite apparent clearance referenced at age 45 years
    # (nominal round value approximating the combined-cohort median;
    # Methods 'Final model development' specifies "a nominal value
    # that approximates the median for the covariate"). Power-model
    # age effect with exponent -0.378 is the only retained
    # subject-level covariate effect.
    lclm <- log(8.83); label("Apparent clearance of 9-OH-risperidone CLM/F at reference age 45 (L/h)")  # Feng 2008 Table 3: CLM = 8.83 L/h (SE 42.6%)
    e_age_clm <- -0.378; label("Power exponent for age on CLM/F (Methods Eq. P_TV = P1*(R/R_ref)^P2)")  # Feng 2008 Table 3: Age on CLM = -0.378 (SE 34.7%)

    # Subpopulation-specific fraction of risperidone metabolized to
    # 9-OH-risperidone. KF in IM is fixed at 1 per Table 3 (SE 'NA')
    # following the Results note: "The model with estimation of all
    # three parameters was unstable, including utilizing the
    # literature published values (0.05, 0.2, 0.3 for PM, IM and EM,
    # respectively) for the three groups, thus KF for subjects in the
    # IM group was fixed to stabilize the model estimation of KF in
    # the PM and EM populations."
    kf_pm <- 0.96;     label("Fraction of risperidone metabolized to 9-OH in CYP2D6 PMs (unitless, 0-1)")  # Feng 2008 Table 3: KF_PM = 0.96 (SE 42.8%)
    kf_em <- 0.595;    label("Fraction of risperidone metabolized to 9-OH in CYP2D6 EMs (unitless, 0-1)")  # Feng 2008 Table 3: KF_EM = 0.595 (SE 40.0%)
    kf_im <- fixed(1); label("Fraction of risperidone metabolized to 9-OH in CYP2D6 IMs; fixed")  # Feng 2008 Table 3: KF_IM = 1 (Fixed)

    # Inter-individual variability (NONMEM OMEGA, variance scale on
    # log eta). Table 3 reports omega as a CV-style percent (column
    # header w%) following the small-CV linearisation CV ~ omega so
    # omega^2 = (w% / 100)^2; this matches the Sherwin 2012
    # convention. No omega is reported in Table 3 for CL/F in IM,
    # for CLM/F, or for the KF parameters; those are deliberately
    # omitted in the final mixture model (the text mentions IIV on
    # CL, V, and CLM in the base-model section but Table 3 lists
    # only the four omegas below for the final mixture model).
    etalka    ~ 0.288  # Feng 2008 Table 3: w_ka% = 53.7 (SE 89.3%); omega^2 = 0.537^2 = 0.288
    etalvc    ~ 0.130  # Feng 2008 Table 3: w_V,VM% = 36.1 (SE 24.4%); omega^2 = 0.361^2 = 0.130
    etalcl_pm ~ 0.920  # Feng 2008 Table 3: w_CL_PM% = 95.9 (SE 39.5%); omega^2 = 0.959^2 = 0.920
    etalcl_em ~ 0.320  # Feng 2008 Table 3: w_CL_EM% = 56.6 (SE 16.8%); omega^2 = 0.566^2 = 0.320

    # Residual error -- combined additive + constant CV (proportional)
    # error model for risperidone and 9-OH-risperidone separately
    # (Methods Equation 2 and surrounding text). Table 3 reports each
    # sigma directly on its native scale (% for proportional, ng/mL
    # for additive; the table calls additive units micrograms per
    # litre which equals nanograms per millilitre and matches the
    # Table 1 observed concentration units).
    propSd     <- 0.639; label("Proportional residual error for risperidone (fraction; SD)")  # Feng 2008 Table 3: sigma_1 % = 63.9 (SE 12.5%)
    addSd      <- 4.29;  label("Additive residual error for risperidone (ng/mL; SD)")  # Feng 2008 Table 3: sigma_3 (ug/L) = 4.29 (SE 104.9%)
    propSd_9oh <- 0.379; label("Proportional residual error for 9-OH-risperidone (fraction; SD)")  # Feng 2008 Table 3: sigma_2 % = 37.9 (SE 35.4%)
    addSd_9oh  <- 0.88;  label("Additive residual error for 9-OH-risperidone (ng/mL; SD)")  # Feng 2008 Table 3: sigma_4 (ug/L) = 0.88 (SE 38.7%)
  })

  model({
    # CYP2D6 phenotype indicator-gating; IM is the implicit reference
    # (both CYP2D6_PM = 0 and CYP2D6_EM = 0 -> IM). Same encoding as
    # Sherwin_2012_risperidone.R.
    ind_im <- (1 - CYP2D6_PM - CYP2D6_EM)

    # Subpopulation-specific apparent oral clearances of risperidone.
    # IIV is applied only to PM and EM (no omega for IM in Table 3);
    # the IM clearance is the deterministic fixed typical value.
    cl_pm <- exp(lcl_pm + etalcl_pm)
    cl_em <- exp(lcl_em + etalcl_em)
    cl_im <- exp(lcl_im)

    # Indicator-gated active CL and KF: only one of the three
    # indicator products is 1 for any subject, so cl equals the
    # corresponding phenotype's cl_* (and similarly for kf).
    cl <- cl_pm * CYP2D6_PM + cl_em * CYP2D6_EM + cl_im * ind_im
    kf <- kf_pm * CYP2D6_PM + kf_em * CYP2D6_EM + kf_im * ind_im

    # Metabolite apparent clearance with age power covariate. No IIV
    # reported in Table 3 for CLM/F so this is deterministic for a
    # given AGE. The Methods covariate-model form P_TV = P1 * (R /
    # R_ref)^P2 with R_ref = 45 years (nominal median).
    clm <- exp(lclm) * (AGE / 45)^e_age_clm

    # Shared apparent volume of distribution. The metabolite Vd is
    # constrained to equal the parent Vd per the Results 'Base model'
    # identifiability note.
    vc     <- exp(lvc + etalvc)
    vc_9oh <- vc

    # Absorption: typical-value Ka fixed, individual Ka has IIV.
    ka <- exp(lka + etalka)

    # Disposition. d/dt(depot): first-order absorption out of the gut.
    # d/dt(central): risperidone eliminated at rate (cl / vc).
    # d/dt(central_9oh): metabolite formed at rate kf * (cl / vc) *
    # central (mass-fraction conversion as in the source NONMEM
    # mixture model; no molar correction is applied -- MW of
    # risperidone 410.5 and 9-OH-risperidone 426.5 differ by ~4%).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - (cl / vc) * central
    d/dt(central_9oh) <-  kf * (cl / vc) * central - (clm / vc_9oh) * central_9oh

    # Plasma concentrations in ng/mL. Dose in mg, volumes in L:
    # central / vc has units mg/L = ug/mL; multiply by 1000 to obtain
    # ng/mL, matching Table 1 observed concentration units.
    Cc     <- 1000 * central     / vc
    Cc_9oh <- 1000 * central_9oh / vc_9oh

    # Combined additive + proportional residual error, separate for
    # parent and metabolite (Table 3 sigma_1..sigma_4).
    Cc     ~ prop(propSd)     + add(addSd)
    Cc_9oh ~ prop(propSd_9oh) + add(addSd_9oh)
  })
}
