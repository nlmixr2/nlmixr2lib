Lioger_2017_rituximab <- function() {
  description <- "Two-compartment population PK model of rituximab in rheumatoid arthritis patients, parameterised with first-order distribution (k12, k21) and elimination (k10) rate constants (per Lioger 2017 Methods 'Structural pharmacokinetics model design': rate-constant parameterisation gave lower AIC, shrinkages, and Vc/elimination correlation than CL/Vc). Final covariate model: BSA, sex, and rituximab treatment course on V1; time-varying CD19+ B-cell count and IgG serum concentration on k10 (Lioger 2017 Table 2)."
  reference <- "Lioger B, Edupuganti SR, Mulleman D, Passot C, Desvignes C, Bejan-Angoulvant T, Thibault G, Gouilleux-Gruart V, Melet J, Paintaud G, Ternant D. Antigenic burden and serum IgG concentrations influence rituximab pharmacokinetics in rheumatoid arthritis patients. Br J Clin Pharmacol. 2017 Sep;83(9):1773-1781. doi:10.1111/bcp.13270. PMID:28230269."
  vignette <- "Lioger_2017_rituximab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central     = list(analyte = "rituximab", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "rituximab", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Median-centred power effect on V1 with exponent 1.0; reference 1.8 m^2 (Lioger 2017 Table 1 median). BSA computation formula not specified in the paper.",
      source_name        = "BSA"
    ),
    SEXF = list(
      description        = "Sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "categorical",
      reference_category = "0 (male)",
      notes              = "Lioger 2017 uses female as the reference category (women are 82.8% of the cohort). The paper's categorical-covariate parameterisation is Ln(V1_TV) = Ln(V1_female) + beta_sex_on_male * SEXM, with beta_sex_on_male = 0.37 (V1 higher in men). Because the canonical column is SEXF (1 = female), the model derives SEXM = 1 - SEXF inline and applies exp(e_sexmale_vc * SEXM). Typical V1 = 4.1 L is the female-reference value from Table 2.",
      source_name        = "SEX"
    ),
    OCC = list(
      description        = "Rituximab treatment course (RTC) number, integer 1..5 (course 4 and 5 were pooled in the paper's categorical analysis but the final model uses RTC as a continuous power covariate)",
      units              = "(count)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Lioger 2017 Table 2 retains RTC as a continuous power covariate on V1 (RTC parameterisation as continuous gave lower AIC than as categorical, 5236.1 vs. 5252.5). Median-centred at RTC = 1 (first course). Same column is also the occasion index the paper used for IOV (see population$notes and 'Assumptions and deviations' in the vignette).",
      source_name        = "RTC"
    ),
    IGG = list(
      description        = "Pre-therapeutic (pre-course) serum immunoglobulin G concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying: measured by nephelometry prior to each rituximab infusion (Methods, Biological data). Median-centred power effect on k10 with exponent 0.50; reference 10 g/L (Lioger 2017 Table 1 medians across the five courses are 10.6, 10.4, 9.2, 8.9, 10.2 g/L; the pooled median is not published as a single number, so a rounded 10 g/L is used per the standing 'undefined reference/centering value -> rounded standard' extraction policy).",
      source_name        = "IgG"
    ),
    CD19_ABS = list(
      description        = "Pre-course CD19+ B-cell count (surrogate of target-antigen burden)",
      units              = "cells/mm^3",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying: measured by flow cytometry prior to each rituximab infusion (Methods, Biological data). Median-centred power effect on k10 with exponent 0.035; reference 100 cells/mm^3 (Lioger 2017 Table 1 per-course medians are 214, 44, 33, 126, 142 cells/mm^3 -- course 1 is the pre-treatment baseline and courses 2-5 fall after B-cell depletion, so the pooled median is much lower than the course-1 value; a rounded 100 cells/mm^3 is used per the standing 'undefined reference/centering value -> rounded standard' extraction policy).",
      source_name        = "CD19"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 64L,
    n_studies      = 1L,
    n_observations = 674L,
    n_courses      = 125L,
    age_range      = "38-84 years",
    age_median     = "59 years",
    sex_female_pct = 82.8,
    bsa_range      = "1.35-2.29 m^2",
    bsa_median     = "1.8 m^2",
    disease_state  = "Rheumatoid arthritis fulfilling American College of Rheumatology criteria. Median disease duration 1.4 years (range 0.27-4.6 years), median initial DAS28 5.24 (range 2.1-8.35). 68.8% rheumatoid factor positive, 82.9% anti-citrullinated protein antibody positive, 79.6% past anti-TNF use.",
    dose_range     = "1000 mg IV infusion on days 1 and 15 per course, 1-5 courses per patient (median 2). Retreatment on an as-needed basis after week 24 based on symptoms of relapse.",
    concomitant_medications = "Corticosteroids 75%, methotrexate 48.4%.",
    regions        = "Single centre, University Hospital of Tours, France, July 2007 - October 2010; follow-up completed November 2012.",
    notes          = "Reported medians from Lioger 2017 Table 1. Rituximab concentrations were quantified by a validated ELISA derived from Blasco et al. 2007 (LOD 0.061 mg/L, LLOQ 0.20 mg/L, ULOQ 9.0 mg/L). 27.9% of samples were below the quantitation limit and interval-censored between 0 and 0.061 mg/L in the source estimation. Source NLME software: MONOLIX 4.3.2 with SAEM (K1 = 700, K2 = 300). Retrospective single-centre therapeutic-drug-monitoring cohort; ethical approval waived as data were routine clinical practice. IOV on V1 (gamma_V1 = 0.33) and k10 (gamma_k10 = 0.08) improved the fit (dOFV = 56.3) but is not encoded structurally here -- see the vignette 'Assumptions and deviations' section."
  )

  ini({
    # Structural parameters -- Lioger 2017 chose rate-constant parameterisation
    # (k10/k12/k21) over CL/Vc parameterisation (Methods 'Structural
    # pharmacokinetics model design': rate-constant parameterisation gave lower
    # AIC 4524.4 vs. 4645.5, lower shrinkages Shk10 = 31% vs. ShCL = 49%, and
    # decreased Vc/elimination-parameter correlation rk10 = 0.15 vs. rCL = 0.32).
    lvc <- log(4.1);  label("Typical central volume of distribution (V1, L) for a median-BSA female at first course") # Lioger 2017 Table 2 final V1 = 4.1 (RSE 7%)
    lk10 <- log(0.11); label("Typical elimination rate constant (k10, 1/day)")                                       # Lioger 2017 Table 2 final k10 = 0.11 (RSE 8%)
    lk12 <- log(0.42); label("Typical central -> peripheral rate constant (k12, 1/day)")                             # Lioger 2017 Table 2 final k12 = 0.42 (RSE 6%)
    lk21 <- log(0.25); label("Typical peripheral -> central rate constant (k21, 1/day)")                             # Lioger 2017 Table 2 final k21 = 0.25 (RSE 4%)

    # Covariate effects on V1 (Lioger 2017 Table 2 "final estimate" column).
    # Continuous covariates use median-centred power form V1 = V1_typ *
    # (COV / med(COV))^beta (Methods 'Covariates'). Categorical sex enters
    # as Ln(V1_TV) = Ln(V1_female) + beta_sex_on_male * SEXM (women are the
    # reference category; men have larger V1).
    e_bsa_vc     <- 1.0;  label("Power exponent of BSA on V1 (unitless; V1 scales with (BSA/1.8)^e_bsa_vc)")             # Lioger 2017 Table 2 BSA on V1 = 1.0 (RSE 38%)
    e_sexmale_vc <- 0.37; label("Additive shift on log-V1 for males vs. female reference (unitless; V1_male = V1_female * exp(e_sexmale_vc))") # Lioger 2017 Table 2 Sex on V1 = 0.37 (RSE 34%)
    e_rtc_vc     <- 0.41; label("Power exponent of rituximab treatment course on V1 (unitless; V1 scales with (RTC/1)^e_rtc_vc)") # Lioger 2017 Table 2 RTC on V1 = 0.41 (RSE 19%)

    # Covariate effects on k10 (Lioger 2017 Table 2 "final estimate" column).
    # Both IgG and CD19+ counts are pre-course (time-varying) measurements
    # entering k10 as median-centred power effects. The positive IgG exponent
    # is hypothesised (Discussion) to reflect FcRn saturation at high
    # endogenous-IgG levels; the positive CD19 exponent reflects
    # target-mediated drug disposition via the CD20+ B-cell pool.
    e_igg_k10    <- 0.50;  label("Power exponent of serum IgG on k10 (unitless; k10 scales with (IGG/10)^e_igg_k10)")     # Lioger 2017 Table 2 IgG on k10 = 0.50 (RSE 19%)
    e_cd19_k10   <- 0.035; label("Power exponent of CD19+ count on k10 (unitless; k10 scales with (CD19_ABS/100)^e_cd19_k10)") # Lioger 2017 Table 2 CD19 on k10 = 0.035 (RSE 27%)

    # Inter-individual variability (Lioger 2017 Table 2 "omega" rows; the
    # legend states omega = interindividual standard deviation on the
    # exponential-random-effect model theta_i = theta_TV * exp(eta_i)).
    # Nlmixr2 stores variance, so eta variance = omega^2.
    # IIV on k12 and k21 were poorly identifiable in the source and set to 0
    # (Results paragraph 2); not encoded here.
    etalvc  ~ 0.0484  # Lioger 2017 Table 2 omega_V1 = 0.22 (RSE 25%); variance = 0.22^2
    etalk10 ~ 0.0576  # Lioger 2017 Table 2 omega_k10 = 0.24 (RSE 13%); variance = 0.24^2

    # Residual error -- combined additive + proportional on the linear scale
    # (Methods 'Interindividual, intercourse and residual models':
    # Yo,ij = Yp,ij * (1 + prop,ij) + add,ij). Table 2 reports the SDs
    # sigma_add and sigma_prop directly.
    addSd  <- 0.31; label("Additive residual error (SD, mg/L)")           # Lioger 2017 Table 2 sigma_add = 0.31 (RSE 4%)
    propSd <- 0.27; label("Proportional residual error (SD, fraction)")   # Lioger 2017 Table 2 sigma_prop = 0.27 (RSE 2%)
  })

  model({
    # Sex indicator: paper uses female as the reference category. Derive
    # SEXM = 1 - SEXF inline so the source-paper coefficient (positive for
    # males) can be applied without changing sign.
    sex_male <- 1 - SEXF

    # Individual parameters. Reference subject: female, BSA = 1.8 m^2, first
    # course (RTC = 1), IgG = 10 g/L, CD19+ = 100 cells/mm^3.
    vc  <- exp(lvc + etalvc) * (BSA / 1.8)^e_bsa_vc * exp(e_sexmale_vc * sex_male) * (OCC / 1)^e_rtc_vc
    k10 <- exp(lk10 + etalk10) * (IGG / 10)^e_igg_k10 * (CD19_ABS / 100)^e_cd19_k10
    k12 <- exp(lk12)
    k21 <- exp(lk21)

    # Two-compartment ODE with first-order distribution and elimination from
    # the central compartment (Lioger 2017 Methods 'Structural
    # pharmacokinetics model design').
    d/dt(central)     <- -k10 * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                    k12 * central - k21 * peripheral1

    # Observed concentration in mg/L (dose in mg, V1 in L).
    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
