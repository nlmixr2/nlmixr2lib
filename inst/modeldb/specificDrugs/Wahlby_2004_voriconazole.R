Wahlby_2004_voriconazole <- function() {
  description <- "Pediatric (2-11 years) two-compartment population PK model for intravenous voriconazole in 35 children, demonstrating Wahlby 2004's extended covariate-model formulation. All disposition parameters scale linearly with body weight. Final-model clearance depends on the time-varying log-ratio (log(ALP/ALP_BASE), 'log(DALKP)' in the source) and on log(ALT) with individual variability in both covariate-effect coefficients (Wahlby 2004 Eq 3 demonstrated). A binary CYP2C19 non-extensive-metabolizer indicator (PM + heterozygous-EM versus homozygous-EM) multiplicatively modifies CL. Underlying structural PK comes from Walsh TJ et al. (Antimicrob Agents Chemother 2004;48(6):2166-2172) and the Karlsson 1995 (J Pharmacokin Biopharm 1998;26(2):207-246) sigma-IIV residual-error pattern is approximated in this entry."
  reference <- paste(
    "Wahlby U, Thomson AH, Milligan PA, Karlsson MO.",
    "Models for time-varying covariates in population pharmacokinetic-pharmacodynamic analysis.",
    "Br J Clin Pharmacol 2004;58(4):367-377.",
    "doi:10.1111/j.1365-2125.2004.02170.x.",
    "Structural PK model carried from Walsh TJ, Karlsson MO, Driscoll T, Arguedas AG et al.",
    "Pharmacokinetics and safety of intravenous voriconazole in children after single- or multiple-dose administration.",
    "Antimicrob Agents Chemother 2004;48(6):2166-2172;",
    "Karlsson 1995 / Karlsson 1998 residual-error pattern: Karlsson MO, Jonsson EN, Wiltse CG, Wade JR.",
    "Assumption testing in population pharmacokinetic models: illustrated with an analysis of moxonidine data from congestive heart failure patients.",
    "J Pharmacokin Biopharm 1998;26(2):207-246. Encodes Wahlby 2004 Table 7 Final-Model column.",
    sep = " "
  )
  vignette <- "Wahlby_2004_time_varying_covariates"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight, time-varying (Walsh 2004 / Wahlby 2004 voriconazole pediatric cohort uses linear WT scaling on all disposition parameters).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear (not allometric) WT scaling on CL, V1, Q, V2 per Methods. The pediatric cohort spans ages 2-11 years.",
      source_name        = "WT"
    ),
    ALT = list(
      description        = "Serum alanine aminotransferase, time-varying.",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Table 3 reports median 25 IU/L overall (range 7-535). Reference value in the log(ALT) effect: 25 IU/L (Eq 7).",
      source_name        = "ALT"
    ),
    ALP = list(
      description        = "Serum alkaline phosphatase, time-varying.",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Table 3 reports median 130 IU/L (range 47-761). Used together with ALP_BASE in the final model to form log(ALP/ALP_BASE) = the paper's log(DALKP) term.",
      source_name        = "ALKP"
    ),
    ALP_BASE = list(
      description        = "Per-subject baseline alkaline phosphatase (BALKP in the source paper), time-fixed.",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-subject constant. Median across 34 subjects with baseline value = 136 IU/L (Table 3, BALKP row). One subject had missing ALKP records imputed to the median baseline value per the source paper. The final model uses log(ALP / ALP_BASE) for the within-subject delta (paper's log(DALKP) term).",
      source_name        = "BALKP"
    ),
    CYP2C19_NON_EM = list(
      description        = "Composite CYP2C19 non-homozygous-extensive-metabolizer indicator (1 = poor or heterozygous-extensive metabolizer, 0 = homozygous extensive metabolizer).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0",
      notes              = "Wahlby 2004 codes this column as 'PM' but groups both true PMs (homozygous loss-of-function) and heterozygous-EMs (one functional, one loss-of-function allele) under value 1. This is a composite carrier indicator rather than the strict CYP2C19_PM category. Distinct from CYP2C19_PM (which is the strict homozygous-PM canonical) and from CYP2C19_S2_CARRIER (which is *2-specific).",
      source_name        = "PM"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 35L,
    n_studies      = 2L,
    age_range      = "2-11 years",
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    disease_state  = "Pediatric subjects (aged 2-11 years) receiving intravenous voriconazole on 1-5 occasions; the cohort pooled two pediatric studies described in Walsh TJ et al. 2004 (Antimicrob Agents Chemother).",
    dose_range     = "Intravenous voriconazole infusions over 0.2-7.6 days follow-up; specific dose ranges are reported by Walsh 2004 and are not reproduced in Wahlby 2004.",
    n_observations = 355L,
    follow_up      = "0.2-7.6 days per subject; on average 10 voriconazole plasma concentrations per individual.",
    regions        = NA_character_,
    notes          = "Pediatric population. Per Methods, disposition parameters are normalised to body weight (linear, not allometric). The Karlsson 1995-style sigma-IIV residual-error pattern reported by Wahlby 2004 (omega_sigma = 0.77) is documented but is approximated in this entry by a homoscedastic proportional residual error; see the vignette Errata."
  )

  ini({
    # Structural typical-value parameters (Wahlby 2004 Table 7 Final-Model column).
    # Disposition parameters are reported per kg body weight; CL, Q in L/h/kg, V1, V2 in L/kg.
    lcl <- log(0.36); label("Population typical CL per kg WT (L/h/kg)")  # Table 7: theta_CL = 0.36 L/h/kg (RSE 11%)
    lvc <- log(0.78); label("Population typical V1 per kg WT (L/kg)")     # Table 7: theta_V1 = 0.78 L/kg (RSE 23%)
    lq  <- log(0.62); label("Population typical Q per kg WT (L/h/kg)")    # Table 7: theta_Q  = 0.62 L/h/kg (RSE 25%)
    lvp <- log(1.62); label("Population typical V2 per kg WT (L/kg)")     # Table 7: theta_V2 = 1.62 L/kg (RSE 11%)

    # Covariate effects on CL (Wahlby 2004 Eq 7 final-model form: log(ALT) effect, log(DALKP)
    # effect, and a CYP2C19 non-EM multiplicative discount)
    e_logalt_cl   <- 0.75; label("log(ALT) effect on CL: coefficient in 1 - e_logalt_cl_i * (log(ALT) - log(25))")        # Table 7: theta_log(ALT)-CL = 0.75 (RSE 26%); reported constrained per footnote a (0.66*0.75-0.33 = 0.23 lower bound)
    e_logdalkp_cl <- 0.59; label("log(ALP/ALP_BASE) effect on CL: coefficient in 1 - e_logdalkp_cl_i * log(ALP/ALP_BASE)") # Table 7: theta_log(DALKP)-CL = 0.59 (RSE 64%)
    e_pm_cl       <- 0.46; label("CYP2C19_NON_EM effect on CL: multiplicative factor (1 - e_pm_cl * CYP2C19_NON_EM)")     # Table 7: theta_PM-CL = 0.46 (RSE 31%); encoded as 1 - theta_PM * PM per the biological constraint that PM = 0 must return factor 1

    # Inter-individual variability on CL (the standard parameter-IIV term)
    etalcl ~ 0.61^2  # Table 7: omega_CL = 0.61 (RSE 11% relative to variance term)

    # IIV in the covariate-effect coefficients (Wahlby 2004 Eq 3 demonstrated)
    etae_logalt_cl   ~ 0.30^2  # Table 7: omega_log(ALT) = 0.30 (RSE 61% rel. to variance); footnote b: 'Approximate SD from logit transformation' - the source paper used a logit constraint to keep CL positive; this entry uses the reported SD directly
    etae_logdalkp_cl ~ 1.94^2  # Table 7: omega_log(DALKP) = 1.94 (RSE 21% relative to variance term)

    # Residual error (Wahlby 2004 Table 7 Final-Model column)
    # The source paper used log-transformed data with additive residual error whose magnitude
    # varies between individuals (Karlsson 1995 / Karlsson 1998 pattern, omega_sigma = 0.77).
    # This entry uses a homoscedastic proportional residual error as an approximation. See
    # the vignette Errata for the deviation rationale.
    propSd <- 0.42; label("Proportional residual error (fraction; approximates the log-transformed additive form)")  # Table 7: sigma = 0.42 (RSE 20%, log-additive); for moderate sigma, log-additive ~= proportional in linear space
  })

  model({
    # Individual covariate-effect coefficients (Wahlby 2004 Eq 3)
    e_logalt_cl_i   <- e_logalt_cl   * exp(etae_logalt_cl)
    e_logdalkp_cl_i <- e_logdalkp_cl * exp(etae_logdalkp_cl)

    # log(DALKP) is interpreted as log(ALP / ALP_BASE) - the within-subject log-ratio,
    # which is positive when ALP increases above its baseline and negative when ALP falls
    # below baseline. This avoids the log-of-negative problem that a literal ALP - ALP_BASE
    # delta would have when ALP decreases.
    log_dalkp <- log(ALP / ALP_BASE)

    # CL with the Wahlby 2004 Eq 7 final-model substitutions. The CYP2C19 non-EM term is
    # encoded as (1 - e_pm_cl * CYP2C19_NON_EM): with theta_PM = 0.46, a homozygous-EM
    # (CYP2C19_NON_EM = 0) gets factor 1, and a non-EM (CYP2C19_NON_EM = 1) gets factor
    # 0.54 - i.e. 46 percent lower CL. Reading the published equation literally
    # ([PM * (1 - theta_PM)]) would zero out CL when PM = 0, which is biologically
    # nonsensical; this multiplicative-discount form is the standard NONMEM equivalent.
    cl <- exp(lcl + etalcl) * WT *
      (1 - e_logalt_cl_i   * (log(ALT) - log(25)) -
           e_logdalkp_cl_i * log_dalkp) *
      (1 - e_pm_cl * CYP2C19_NON_EM)

    # Linear WT scaling on the rest of the disposition parameters per Methods
    vc <- exp(lvc) * WT
    q  <- exp(lq)  * WT
    vp <- exp(lvp) * WT

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment IV voriconazole
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
