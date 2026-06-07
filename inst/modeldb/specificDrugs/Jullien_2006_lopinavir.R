Jullien_2006_lopinavir <- function() {
  description <- "One-compartment population PK model for oral lopinavir (boosted with ritonavir) in HIV-infected children from birth to 18 years, with the absorption and elimination rate constants constrained to a single shared rate constant k = CL/F divided by V/F (Jullien 2006, simplified parameterisation per Wahlby 2002). Body weight is allometrically scaled on CL/F and V/F (reference 27 kg), nevirapine coadministration increases CL/F by 34%, and male sex increases CL/F by 39% in children older than 12 years."
  reference <- "Jullien V, Urien S, Hirt D, Delaugerre C, Rey E, Teglas JP, Vaz P, Rouzioux C, Chaix ML, Macassa E, Firtion G, Pons G, Blanche S, Treluyer JM. Population Analysis of Weight-, Age-, and Sex-Related Differences in the Pharmacokinetics of Lopinavir in Children from Birth to 18 Years. Antimicrob Agents Chemother. 2006 Nov;50(11):3548-55. doi:10.1128/AAC.00943-05"
  vignette <- "Jullien_2006_lopinavir"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight (time-varying); reference 27 kg (cohort median per Jullien 2006 Table 1).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on CL/F and V/F with reference 27 kg (the cohort median BW). Source column name was BW in the NONMEM dataset; canonical name WT is used here.",
      source_name        = "BW"
    ),
    AGE = list(
      description        = "Subject age in years (time-varying).",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used as a 12-year threshold to gate the sex effect on CL/F (the sex effect applies only when AGE > 12). The age-stratified effect was identified by Jullien 2006 (Results 'Population pharmacokinetics' paragraph 2 and final covariate submodel).",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Biological sex indicator: 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Jullien 2006 used the inverse coding S = 1 for boys / 0 for girls. The canonical SEXF column inverts the source values (SEXF = 1 - S). The published 39% increase in CL/F for boys older than 12 years is preserved by applying the effect as exp(e_sexf_cl * (1 - SEXF) * (AGE > 12)) inside model().",
      source_name        = "SEX"
    ),
    CONMED_NVP = list(
      description        = "Concomitant nevirapine indicator: 1 = nevirapine coadministered with lopinavir, 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no nevirapine)",
      notes              = "Nevirapine is an NNRTI and a CYP3A inducer; coadministration increased lopinavir CL/F by 34% in Jullien 2006 (Results 'Population pharmacokinetics' paragraph 2, final covariate submodel). Source column was the in-equation indicator N (1 if combined with nevirapine, 0 otherwise). This canonical entry is the nevirapine analog of the registered CONMED_EFV indicator.",
      source_name        = "N"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 157L,
    n_studies      = 1L,
    n_observations = 541L,
    age_range      = "3 days to 18 years",
    age_median     = "10.2 years",
    weight_range   = "2 to 73 kg",
    weight_median  = "27.6 kg",
    sex_female_pct = 42.7,
    race_ethnicity = "Mixed (paediatric HIV cohorts from Cochin-Saint-Vincent de Paul and Necker-Enfants Malades, Paris); paper does not stratify by race.",
    disease_state  = "HIV-1 infection (or maternal-fetal transmission prophylaxis) on lopinavir/ritonavir-containing combination antiretroviral therapy.",
    dose_range     = "Twice-daily oral lopinavir/ritonavir (mean dose 279 mg lopinavir, range 30-532 mg; 109 mg/kg/day mean, range 4.4-29.4 mg/kg per dose). Liquid formulation used in younger children; solid oral formulation used in older children.",
    regions        = "France (Paris).",
    notes          = "Retrospective therapeutic-drug-monitoring cohort. Concomitant ART: at least one nucleoside reverse-transcriptase inhibitor in 90% of samples, one protease inhibitor in 10%, one non-nucleoside reverse-transcriptase inhibitor in 23%. Nevirapine combined with lopinavir in 16% of samples; efavirenz in 8%; amprenavir in 4%. Baseline demographics from Jullien 2006 Table 1. Median 3 samples per patient (range 1-14)."
  )

  ini({
    # Structural parameters (Jullien 2006 Table 3, final model original data set).
    # The paper uses a one-compartment model in which a single first-order rate
    # constant k governs both absorption and elimination (k = ka = kel) per the
    # Wahlby 2002 simplified parameterisation. The estimated parameters are CL/F
    # and V/F; the rate constant is derived inside model() as kel = cl/vc.
    lcl <- log(2.58); label("Apparent clearance at WT=27 kg (CL/F, L/h)")             # Table 3 row 'Final model, Mean': TV(CL/F) = 2.58 L/h
    lvc <- log(24.6); label("Apparent volume of distribution at WT=27 kg (V/F, L)")   # Table 3 row 'Final model, Mean': TV(V/F) = 24.6 L

    # Allometric exponents on body weight, reference WT = 27 kg (cohort median).
    # Estimated -- reported with SEs in Table 3; NOT held fixed.
    e_wt_cl <- 0.46; label("Allometric exponent on CL/F with body weight (unitless)")  # Table 3: theta_BW on CL/F
    e_wt_vc <- 0.72; label("Allometric exponent on V/F with body weight (unitless)")   # Table 3: theta_BW on V/F

    # Nevirapine coadministration effect on CL/F. Paper reports the multiplicative
    # factor as 1.34 (Table 3 theta_nevirapine and the final covariate submodel
    # equation). Encoded as a log-additive coefficient: exp(e_nvp_cl) = 1.34.
    e_nvp_cl <- log(1.34); label("Effect of nevirapine coadministration on log CL/F (multiplier 1.34)")  # Table 3: theta_nevirapine

    # Sex effect on CL/F among adolescents (AGE > 12). Paper reports the
    # multiplicative factor as 1.39 (Table 3 theta_SEX (if age > 12 yr) and the
    # final covariate submodel equation). Paper convention: S = 1 for boys, so
    # the +39% applies to boys older than 12. Encoded with the SEXF canonical
    # (SEXF = 1 - S), so the effect is applied to (1 - SEXF) * (AGE > 12) inside
    # model() and exp(e_sexf_cl) = 1.39 multiplies CL/F for adolescent boys only.
    e_sexf_cl <- log(1.39); label("Effect of male sex on log CL/F when AGE > 12 yr (multiplier 1.39)")  # Table 3: theta_SEX (if age > 12 yr)

    # Inter-individual variability block on log CL/F and log V/F.
    # Jullien 2006 Table 3 reports omega^2(CL/F) = 0.0957, omega^2(V/F) = 0.180,
    # and cov(CL,V) = 0.131 (final model, original data set means). Encoded as a
    # 2x2 block with variances on the log scale (exponential / log-normal IIV).
    etalcl + etalvc ~ c(0.0957,
                        0.131, 0.180)  # Table 3: omega^2_CL/F, cov_CL,V, omega^2_V/F

    # Residual error. Paper reports a combined exponential (sigma_1^2 = 0.138 on
    # log scale) and additive (sigma_2^2 = 1.83 mg^2/L^2) error model. nlmixr2's
    # native simulation-time residual does not support add(...) + lnorm(...)
    # combinations; following the Ruhs 2012 methotrexate precedent, the
    # exponential component is encoded as the small-sigma-equivalent proportional
    # SD (propSd = sqrt(0.138)) and the additive component is the SD on the
    # linear concentration scale (addSd = sqrt(1.83)). The combined-error model
    # is approximately equivalent at the small-CV limit; see the vignette
    # Assumptions and deviations section for the exact equivalence and
    # numerical tolerance.
    propSd <- sqrt(0.138); label("Proportional residual SD (small-sigma equivalent of exponential, log-scale CV ~ 37%)")  # Table 3: sigma_1^2 = 0.138
    addSd  <- sqrt(1.83);  label("Additive residual SD (mg/L)")                                                            # Table 3: sigma_2^2 = 1.83 mg^2/L^2
  })

  model({
    # Individual PK parameters. CL/F is allometrically scaled with body weight
    # (reference 27 kg) and modified by the nevirapine coadministration indicator
    # and the adolescent-male sex indicator (effective only when AGE > 12).
    cl_size <- (WT / 27)^e_wt_cl
    cl_nvp  <- exp(e_nvp_cl * CONMED_NVP)
    cl_sex  <- exp(e_sexf_cl * (1 - SEXF) * (AGE > 12))
    cl <- exp(lcl + etalcl) * cl_size * cl_nvp * cl_sex
    vc <- exp(lvc + etalvc) * (WT / 27)^e_wt_vc

    # Simplified one-compartment kinetics with a single shared rate constant
    # k = CL/V used for both absorption and elimination (Wahlby 2002, cited as
    # reference 31 by Jullien 2006). The classical 1-cmt model with first-order
    # absorption was overparameterised for this dataset (paper Results
    # 'Population pharmacokinetics' paragraph 1); the constrained ka = kel form
    # produced a 20-point further reduction in the objective function.
    kel <- cl / vc
    ka  <- kel

    d/dt(depot)   <- -ka  * depot
    d/dt(central) <-  ka  * depot - kel * central

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
