Honda_2005_carvedilol <- function() {
  description <- paste(
    "One-compartment population PK model for orally administered racemic",
    "carvedilol in 23 healthy Japanese volunteers, with R- and S-enantiomer",
    "whole-blood concentrations measured by chiral HPLC at 2 h and 6 h after",
    "a single 5- or 10-mg oral dose (Honda 2005). NONMEM ADVAN1/TRANS2 with",
    "very rapid absorption: the racemic dose is split equally between two",
    "parallel central compartments (central_r, central_s) with no separate",
    "absorption depot. CL/F and V/F scale linearly with body weight; an S/R",
    "ratio theta_3 (CL/F) and theta_4 (V/F) parameterise the stereoselective",
    "difference. One subject-level eta on CL/F and one on V/F are shared",
    "between enantiomers (correlated block IIV, rho ~ 0.90). Power-variance",
    "residual error with fixed exponent 1/2 (Honda Eq. 3), shared between",
    "R- and S-enantiomer observations. CYP2D6*10 genotype is not in the",
    "structural model; Honda 2005 reports the *10-carrier effect only as a",
    "post-hoc stratification of the individual Bayes estimates (Figs. 3-4)."
  )
  reference <- paste(
    "Honda M, Nozawa T, Igarashi N, Inoue H, Arakawa R, Ogura Y,",
    "Okabe H, Taguchi M, Hashimoto Y.",
    "Effect of CYP2D6*10 on the pharmacokinetics of R- and S-carvedilol",
    "in healthy Japanese volunteers.",
    "Biol Pharm Bull. 2005;28(8):1476-1479.",
    "doi:10.1248/bpb.28.1476"
  )
  vignette <- "Honda_2005_carvedilol"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Honda 2005 Methods 'Subjects and Study Protocols': cohort body",
        "weight 47-86 kg (mean 64.7 kg). The 5-mg subgroup mean (SD) was",
        "57.9 (6.7) kg (n=9); the 10-mg subgroup was 69.1 (9.9) kg (n=14).",
        "Body weight enters Eq. 1 and Eq. 2 linearly (allometric exponent",
        "= 1): CL/F_i = theta_1 * theta_3^S * WT_i * (1 + eta_CL/F_i) and",
        "V/F_i = theta_2 * theta_4^S * WT_i * (1 + eta_V/F_i). The",
        "underlying theta_1 and theta_2 are reported in per-kg units",
        "(L/h/kg and L/kg). Time-fixed at baseline in the published model."
      ),
      source_name        = "WT"
    )
  )

  covariatesDataExcluded <- list(
    CYP2D6_PM = list(
      description        = paste(
        "Honda 2005 post-hoc stratification indicator: 1 = subject carries",
        "at least one CYP2D6*10 allele; 0 = subject is CYP2D6*1/*1 or",
        "*1/*2 (the *10-non-carrier reference). NOT a structural-model",
        "covariate. The paper reports significantly lower (CL/F)/WT and",
        "(V/F)/WT in *10-allele carriers vs the reference complement",
        "(Honda 2005 Figs. 3-4; p<0.001 for R-carvedilol CL/F), but the",
        "effect is shown only as a t-test on the individual Bayes",
        "estimates -- the NONMEM structural model itself does not include",
        "CYP2D6 as a covariate."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP2D6*1/*1 or *1/*2)",
      notes              = paste(
        "CYP2D6 alleles determined by PCR-RFLP (*1, *10, *14) and",
        "allele-specific PCR (*2); long-PCR for *5 (Honda 2005 Methods",
        "'Genotyping of CYP2D6'). Cohort distribution: 5 *1/*1, 1 *1/*2,",
        "12 *1/*10, 3 *2/*10, 2 *10/*10. No null alleles (*5 or *14)",
        "detected. The broader canonical CYP2D6_PM (poor-metabolizer",
        "indicator) in inst/references/covariate-columns.md is used here",
        "as documentation only; downstream users who want the *10-carrier",
        "stratification effect can post-process the Bayes-estimate output",
        "rather than re-fit with CYP2D6 as a covariate."
      ),
      source_name        = "(Honda 2005 reports per-subject CYP2D6 genotypes; no NONMEM column name)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 23L,
    n_studies      = 1L,
    age_range      = "22-44 years",
    age_median     = "29.1 years (cohort mean)",
    weight_range   = "47-86 kg",
    weight_median  = "64.7 kg (cohort mean)",
    sex_female_pct = round(4 / 23 * 100, 1),
    race_ethnicity = "Japanese (single-centre Toyama, Japan)",
    disease_state  = "Healthy adult volunteers (all participants were physicians or pharmacists; no concomitant medications reported).",
    dose_range     = "Single oral dose of racemic carvedilol (Artist tablet): 5 mg as two 2.5-mg tablets (n=9; 7 men, 2 women) or 10 mg as one 10-mg tablet (n=14; 12 men, 2 women). Doses were self-selected, taken with water at least 2 h before a meal, after an overnight fast.",
    regions        = "Japan (Toyama Medical and Pharmaceutical University)",
    n_observations = 46L,
    notes          = paste(
      "Sampling: 5 mL of whole blood at 2 h and 6 h post-dose (two-",
      "time-point sparse design). Both R- and S-carvedilol were assayed",
      "in every sample, so each subject contributes two R observations",
      "+ two S observations (~92 enantiomer-resolved concentrations in",
      "total). Chiral HPLC method per Mihara et al. with assay CVs <=",
      "4.1% at 0.5-3 ng/mL; LOQ 0.2 ng/mL. Parameter estimates obtained",
      "with the Bayesian NONMEM post-hoc option (POSTHOC) starting from",
      "the population estimates. CYP2D6 genotype distribution in the",
      "cohort: 5 *1/*1, 1 *1/*2, 12 *1/*10, 3 *2/*10, 2 *10/*10 (no *5",
      "or *14 null alleles detected); see covariatesDataExcluded."
    )
  )

  ini({
    # Honda 2005 Table 1 final NONMEM estimates (95% CI = +/-1.96 * S.E.).
    # The S indicator in Honda Eq. 1 / 2 is fixed to one for S-carvedilol
    # and to zero for R-carvedilol; the encoding below maps Honda's
    # theta_3^S and theta_4^S onto log-additive S/R ratios lratiocl /
    # lratiov so that R-enantiomer parameters use lclf / lvf alone and
    # S-enantiomer parameters use lclf + lratiocl and lvf + lratiov.
    lclf     <- log(1.01) ; label("R-carvedilol CL/F per kg body weight (theta_1, L/h/kg)") # Honda 2005 Table 1: theta_1 = 1.01 L/h/kg (95% CI 0.84-1.18)
    lvf      <- log(2.53) ; label("R-carvedilol V/F per kg body weight (theta_2, L/kg)")    # Honda 2005 Table 1: theta_2 = 2.53 L/kg (95% CI 2.04-3.02)
    lratiocl <- log(2.13) ; label("Log S/R ratio for CL/F (theta_3; unitless)")              # Honda 2005 Table 1: theta_3 = 2.13 (95% CI 1.64-2.62)
    lratiov  <- log(2.94) ; label("Log S/R ratio for V/F (theta_4; unitless)")               # Honda 2005 Table 1: theta_4 = 2.94 (95% CI 1.98-3.90)

    # Correlated block IIV between eta_CL/F and eta_V/F. Honda 2005
    # Eq. 1 and Eq. 2 use a single eta per subject for CL/F (shared
    # across R and S enantiomers) and a single eta per subject for V/F
    # (also shared across R and S). Lower-triangle order
    # c(var(etalclf), cov, var(etalvf)) per nlmixr2 block syntax.
    # Honda reports omega_CL/F,V/F = 0.130, giving correlation
    # rho = 0.130 / sqrt(0.130 * 0.161) = 0.898 (paper says 0.899).
    etalclf + etalvf ~ c(0.130, 0.130, 0.161)
    # Honda 2005 Table 1: omega^2(CL/F) = 0.130 (95% CI 0.073-0.187);
    #                     omega^2(V/F)  = 0.161 (95% CI 0.054-0.268);
    #                     omega(CL/F,V/F) = 0.130 (95% CI 0.058-0.202).

    # Power-variance residual error. Honda 2005 Eq. 3 writes
    #   Cb_ij = Cb*_ij + Cb*_ij^(1/2) * eps_ij,  var(eps) = sigma^2
    # encoded in nlmixr2 as 'Cc ~ pow(propSd, powExp)' with SD(Cc - Cb*) =
    # propSd * Cb*^powExp. propSd = sqrt(sigma^2); powExp = 1/2 (fixed).
    # Concentrations are interpreted on the ng/mL scale (Honda Fig. 1).
    # Honda 2005 fits ONE sigma shared between R- and S-enantiomer
    # observations (one $SIGMA entry, single sigma^2 = 0.0584). nlmixr2
    # requires a distinct residual-SD parameter per `~` endpoint, so the
    # shared sigma is encoded as two separate but numerically identical
    # parameters (propSd_r = propSd_s = sqrt(0.0584)); a re-fit that
    # estimates the two would break Honda's structural shared-sigma
    # assumption -- see vignette 'Assumptions and deviations'.
    propSd_r <- sqrt(0.0584) ; label("Power-error SD coefficient for R-carvedilol ((ng/mL)^(1/2), = sqrt(Honda sigma^2))")  # Honda 2005 Table 1: sigma^2 = 0.0584 (95% CI 0.0074-0.1094); shared between R and S
    propSd_s <- sqrt(0.0584) ; label("Power-error SD coefficient for S-carvedilol ((ng/mL)^(1/2), = sqrt(Honda sigma^2))")  # Honda 2005 Table 1: sigma^2 = 0.0584 (95% CI 0.0074-0.1094); shared between R and S
    powExp_r <- fixed(0.5)   ; label("Power-error exponent for R-carvedilol (Honda Eq. 3: Cb*^(1/2); unitless, FIXED)")      # Honda 2005 Eq. 3 (text after Eq. 2)
    powExp_s <- fixed(0.5)   ; label("Power-error exponent for S-carvedilol (Honda Eq. 3: Cb*^(1/2); unitless, FIXED)")      # Honda 2005 Eq. 3 (text after Eq. 2)
  })

  model({
    # NOTE on the IIV encoding. Honda 2005 Eq. 1 and Eq. 2 write the
    # individual CL/F and V/F as CL/F_typ * (1 + eta) and V/F_typ *
    # (1 + eta) -- a linear (1 + eta) form on untransformed CL/F and
    # V/F. The model below uses the nlmixr2 standard exp(lparam +
    # etalparam) log-normal IIV form as an approximation. For the modest
    # reported variances (omega^2 = 0.130 and 0.161, CV ~ 36% and ~40%),
    # the linear (1 + eta) form and the log-normal exp(eta) form agree
    # to <~5% on per-subject CL/F across the bulk of the distribution
    # but the log-normal form avoids unphysical negative CL/F in the
    # extreme tail of eta. See the validation vignette
    # ('Assumptions and deviations' section) for details.

    # Typical CL/F and V/F (units L/h and L) for R- and S-enantiomers.
    # WT enters linearly per Eq. 1-2 (allometric exponent = 1). The
    # same etalclf drives both R and S CL/F; the same etalvf drives
    # both R and S V/F (one subject-level pair of etas, not two).
    # The S/R ratios exp(lratiocl) and exp(lratiov) are kept outside
    # the mu-referenced exp(... + eta) block so nlmixr2's mu-reference
    # parser sees a single fixed-effect parameter per eta-bearing
    # exponential (otherwise it errors "2+ population parameters in a
    # single mu-referenced expression").
    cl_r <- exp(lclf + etalclf) * WT
    cl_s <- exp(lclf + etalclf) * exp(lratiocl) * WT
    v_r  <- exp(lvf  + etalvf)  * WT
    v_s  <- exp(lvf  + etalvf)  * exp(lratiov)  * WT

    kel_r <- cl_r / v_r
    kel_s <- cl_s / v_s

    # One-compartment with very rapid absorption: the racemic oral dose
    # is delivered as a bolus to each central compartment (no depot is
    # carried; NONMEM ADVAN1/TRANS2 with effectively instantaneous
    # absorption). The user data set provides two dose records per
    # administration: cmt = central_r with amt = 0.5 * racemic_dose_mg,
    # and cmt = central_s with amt = 0.5 * racemic_dose_mg.
    d/dt(central_r) <- -kel_r * central_r
    d/dt(central_s) <- -kel_s * central_s

    # Honda 2005 observed blood concentrations are in ng/mL (Fig. 1 and
    # Table 1's sigma^2 = 0.0584 (ng/mL)), while doses are in mg and
    # V/F is in L; convert mg/L -> ng/mL via the factor 1000.
    Cc_r <- 1000 * central_r / v_r
    Cc_s <- 1000 * central_s / v_s

    # Power-variance residual error (Honda's single sigma encoded as
    # paired propSd_r = propSd_s = sqrt(0.0584); see ini() note).
    Cc_r ~ pow(propSd_r, powExp_r)
    Cc_s ~ pow(propSd_s, powExp_s)
  })
}
