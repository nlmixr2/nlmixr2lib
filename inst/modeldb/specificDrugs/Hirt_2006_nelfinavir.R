Hirt_2006_nelfinavir <- function() {
  description <- paste(
    "Population PK model for oral nelfinavir and its active metabolite",
    "M8 (hydroxy-tert-butylamide) in 182 pediatric HIV-1 infected children",
    "aged 3 days to 17 years (Hirt 2006). One-compartment model for",
    "nelfinavir (depot + central) with first-order absorption (Ka) and",
    "linear elimination via apparent total clearance CL_T/F; the active",
    "metabolite M8 is described by a single compartment (central_m8) with",
    "apparent volume FIXED to 1 L (not identifiable). Only the fraction",
    "F_MT (~2.5%) of nelfinavir's total clearance enters M8; the remainder",
    "is lost to non-M8 elimination pathways. Body-weight scaling is linear",
    "(per-kg parameterisation of V/F and CL/F); both V/F and CL/F decrease",
    "with age via a shared power exponent of -0.29 relative to the median",
    "age 8.2 years. M8 elimination rate KM0 is increased ~1.9-fold by",
    "concomitant administration of an enzyme-inducing NNRTI (efavirenz or",
    "nevirapine, pooled under the indicator CONMED_NNRTI_IND consistent",
    "with the paper's finding that the two drugs' inducer effects on KM0",
    "were not significantly different and the two were never administered",
    "simultaneously). Inter-individual variability on V/F, CL/F, and KM0;",
    "correlation 0.45 between IIVs on CL/F and KM0; additive residual",
    "error in mg/L for each output."
  )
  reference <- "Hirt D, Urien S, Jullien V, Firtion G, Rey E, Pons G, Blanche S, Treluyer JM. Age-related effects on nelfinavir and M8 pharmacokinetics: a population study with 182 children. Antimicrob Agents Chemother. 2006;50(3):910-916. doi:10.1128/aac.50.3.910-916.2006"
  vignette  <- "Hirt_2006_nelfinavir"
  units     <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/L"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight; multiplier on the per-kg V/F and CL/F.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-kg parameterisation: V/F (L) = (V/F per kg) * WT and CL/F (L/h) = (CL/F per kg) * WT. Linear (exponent 1) scaling was retained over the theoretical allometric exponents 0.75 on CL and 1 on V because no significant OFV or goodness-of-fit improvement was observed (Results, 'Nelfinavir pharmacokinetic model building'). Cohort range 1.7-70 kg; median 21 kg.",
      source_name        = "BW"
    ),
    AGE = list(
      description        = "Subject age in years; drives a shared power-form age effect on V/F and CL/F.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form (AGE/8.2)^-0.29 applied to BOTH V/F and CL/F. Reference age 8.2 years is the median of the 182-child cohort. Range in the source: 3 days (i.e., 0.0082 years) to 17 years. For neonates whose age is recorded in days, convert via AGE = days/365.25 before applying the model.",
      source_name        = "AGE"
    ),
    CONMED_NNRTI_IND = list(
      description        = "Pooled enzyme-inducing NNRTI coadministration indicator (1 = subject coadministered efavirenz or nevirapine, 0 = otherwise).",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Pooled because the source paper found the inducer effects of efavirenz and nevirapine on KM0 were statistically indistinguishable, and the two drugs were never administered simultaneously in the cohort (efavirenz n=10/53 obs; nevirapine n=33/133 obs). Drives the multiplicative form KM0 = 1.88 * (1 + 0.91 * CONMED_NNRTI_IND), consistent with CYP3A4 induction of M8 elimination. Saquinavir (n=10/48 obs) was tested but not significant; ritonavir (n=3/11 obs) was excluded from the final analysis (unstable model). For a future paper that needs to distinguish efavirenz vs. nevirapine effects separately, use CONMED_EFV and a (future) CONMED_NVP rather than this pooled indicator.",
      source_name        = "NNI"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 182L,
    n_studies       = 1L,
    age_range       = "3 days to 17 years",
    age_median      = "8.2 years",
    weight_range    = "1.7-70 kg",
    weight_median   = "21 kg",
    sex_female_pct  = 47.8,
    disease_state   = "Pediatric HIV-1 infection on oral nelfinavir as part of combination antiretroviral therapy. Nelfinavir given only as 250-mg tablets (crumbled in water and added to milk or food for children unable to swallow); powder formulation not used due to large administration volume, unpleasant consistency, and dissolution difficulties.",
    dose_range      = "Per-administration mean (SD) by age cohort (Table 1): <2 mo BID 147 (61) mg / TID 130 (70) mg; 2 mo-2 yr BID 504 (137) mg / TID 233 (89) mg; 2-7 yr TID 416 (108) mg; >=8 yr BID 849 (278) mg / TID 655 (131) mg. Typical mg/kg doses 25-60 mg/kg per administration.",
    regions         = "France (multicentre; Paris-area pediatric HIV clinics including Hopital Cochin-Saint-Vincent-de-Paul and Hopital Necker-Enfants Malades).",
    n_observations  = "742 nelfinavir + 557 M8 plasma concentrations (median 3 nelfinavir and 2 M8 samples per patient). Most samples at steady state (>= 10 days of treatment); 18 plasma samples in neonates younger than 10 days were not at steady state.",
    age_cohorts     = "Three pre-specified age cohorts for the FDA-recommendation evaluation: <2 mo (n=25); 2 mo-2 yr (n=36); 2-13 yr (n=121).",
    co_medication   = "Coadministered antiretrovirals known to influence nelfinavir or M8: efavirenz (n=10 subjects, 53 samples), nevirapine (n=33, 133 samples), ritonavir (n=3, 11 samples; excluded from final analysis), saquinavir (n=10, 48 samples; tested, not significant).",
    notes           = "Demographics from Hirt 2006 Table 1 and the Results 'Demographic data' paragraph. Therapeutic drug monitoring data; ethics committee approval not required under French TDM regulations at the time of the study. Concentrations < LOQ were set to half the LOQ (0.1 ug/mL). NONMEM v V level 1.1 with FOCE+I. Bootstrap validation: 1000 resamples (Table 2)."
  )

  ini({
    # Structural PK parameters. All clearance and volume terms are per
    # kg body weight; the per-kg values multiply by WT inside model() to
    # give whole-subject values. Reference age 8.2 years is the cohort
    # median.

    # Nelfinavir absorption rate (no IIV: ISV on Ka was not identifiable
    # and excluding it had no influence on OFV; Hirt 2006 Results,
    # 'Nelfinavir pharmacokinetic model building' paragraph).
    lka <- log(0.48)
    label("Nelfinavir first-order absorption rate constant (1/h)")  # Hirt 2006 Table 2: Ka mean 0.48 1/h (RSE 32%)

    # Nelfinavir apparent total clearance per kg (CL_T/F) at AGE = 8.2 yr.
    lcl <- log(0.93)
    label("Nelfinavir apparent total clearance per kg at AGE = 8.2 yr (L/h/kg)")  # Hirt 2006 Table 2: CL_T mean 0.93 L/h/kg (RSE 4.6%)

    # Nelfinavir apparent central volume per kg (V/F) at AGE = 8.2 yr.
    lvc <- log(6.86)
    label("Nelfinavir apparent central volume of distribution per kg at AGE = 8.2 yr (L/kg)")  # Hirt 2006 Table 2: V mean 6.86 L/kg (RSE 28%)

    # Nelfinavir-to-M8 apparent formation clearance fraction. No IIV:
    # ISV on F_MT was not significant and exclusion had no influence on
    # OFV; Hirt 2006 Results, '(ii) M8 pharmacokinetic model building'.
    lfmt_m8 <- log(0.025)
    label("Nelfinavir-to-M8 apparent formation clearance fraction (unitless)")  # Hirt 2006 Table 2: F_MT mean 0.025 (RSE 15%)

    # M8 apparent terminal elimination rate constant (K_M0 in the paper).
    lkel_m8 <- log(1.88)
    label("M8 apparent elimination rate constant at CONMED_NNRTI_IND = 0 (1/h)")  # Hirt 2006 Table 2: K_M0 mean 1.88 1/h (RSE 16%)

    # M8 apparent distribution volume FIXED to 1 L. K_M0 and Vm are not
    # jointly identifiable (only their ratio CLm0/Vm is); the paper
    # fixes Vm = 1 so K_M0 absorbs the physical M8 volume (Appendix:
    # "K_M0 = CLm0/Vm, with Vm = 1").
    lvc_m8 <- fixed(log(1))
    label("M8 apparent distribution volume (L) -- FIXED to 1 per Hirt 2006 (not identifiable)")  # Hirt 2006 Appendix: 'K_M0 = CLm0/Vm with Vm = 1'

    # Shared age exponent on V/F and CL/F. The paper applies the SAME
    # age effect to both V and CL (a 101-U combined OFV decrease vs.
    # the age-on-CL-alone submodel); Results, 'Nelfinavir pharmacokinetic
    # model building'.
    e_age_cl_vc <- -0.29
    label("Shared power exponent of (AGE/8.2 yr) on CL/F and V/F (unitless)")  # Hirt 2006 Table 2: theta_AGE mean -0.29 (RSE 12%)

    # Pooled NNRTI-inducer effect on M8 elimination (efavirenz or
    # nevirapine; never co-administered together in the cohort). Applied
    # multiplicatively as K_M0 * (1 + theta_NNI * CONMED_NNRTI_IND);
    # CONMED_NNRTI_IND = 1 inflates K_M0 by 1 + 0.91 = 1.91, consistent
    # with the 1.9-fold lower M8 concentrations reported in Discussion.
    e_conmed_nnrti_ind_kel_m8 <- 0.91
    label("Additive-form multiplicative effect of CONMED_NNRTI_IND on M8 elimination rate (unitless)")  # Hirt 2006 Table 2: theta_NNI mean 0.91 (RSE 25%)

    # Inter-individual variability. CV% from Hirt 2006 Table 2 converted
    # to log-normal omega^2 via omega^2 = log(1 + CV^2).
    #   ISV V/F   = 109%   -> omega^2 = log(1 + 1.09^2)   = 0.7831
    #   ISV CL_T  = 39.1%  -> omega^2 = log(1 + 0.391^2)  = 0.1423
    #   ISV K_M0  = 49.2%  -> omega^2 = log(1 + 0.492^2)  = 0.2169
    #   r(CL, KM0) = 0.45  -> cov = 0.45 * sqrt(0.1423) * sqrt(0.2169)
    #                            = 0.45 * 0.3772 * 0.4657
    #                            = 0.0790
    # IIV on V/F is independent; IIVs on CL/F and KM0 are correlated
    # (Hirt 2006 Results, '(iii) Nelfinavir-M8 pharmacokinetic model
    # building': covariance term between total clearance and M8
    # elimination rate produced a 11-U OFV decrease).
    etalvc ~ 0.7831  # Hirt 2006 Table 2: ISV V 109% (RSE 53%)
    etalcl + etalkel_m8 ~ c(0.1423,
                            0.0790, 0.2169)  # Hirt 2006 Table 2: ISV CL_T 39.1% (RSE 22%), K_M0 49.2% (RSE 26%), r(CL,KM0) = 0.45 (RSE 36%)

    # Residual error: additive on each analyte concentration. Hirt 2006
    # Results, 'Nelfinavir pharmacokinetic model building': "intersubject
    # and residual variabilities were best described by exponential and
    # additive error models, respectively."
    addSd    <- 1.65
    label("Nelfinavir additive residual SD (mg/L)")  # Hirt 2006 Table 2: sigma_NELFI mean 1.65 ug/mL = mg/L (RSE 8.9%)
    addSd_m8 <- 0.63
    label("M8 additive residual SD (mg/L)")  # Hirt 2006 Table 2: sigma_M8 mean 0.63 ug/mL = mg/L (RSE 21%)
  })

  model({
    # 1. Individual PK parameters. Body-weight scaling is linear
    # (per-kg parameters multiplied by WT); age scaling is power-form
    # with the SAME exponent on CL/F and V/F.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * WT * (AGE / 8.2)^e_age_cl_vc
    vc <- exp(lvc + etalvc) * WT * (AGE / 8.2)^e_age_cl_vc

    # 2. M8 parameters. F_MT is a fraction (unitless), no IIV. vc_m8 is
    # FIXED to 1 L (not identifiable). kel_m8 has an additive-form
    # multiplicative NNRTI-inducer effect.
    fmt_m8 <- exp(lfmt_m8)
    vc_m8  <- exp(lvc_m8)
    kel_m8 <- exp(lkel_m8 + etalkel_m8) *
              (1 + e_conmed_nnrti_ind_kel_m8 * CONMED_NNRTI_IND)

    # 3. Nelfinavir elimination rate constant from CL_T/F and V/F.
    kel <- cl / vc

    # 4. ODE system (Hirt 2006 Appendix). Nelfinavir absorption is
    # first-order from the depot; total clearance removes nelfinavir
    # at rate kel*central, of which only the F_MT fraction is observed
    # as M8 formation. M8 is eliminated terminally at rate kel_m8.
    d/dt(depot)      <- -ka * depot
    d/dt(central)    <-  ka * depot - kel * central
    d/dt(central_m8) <-  fmt_m8 * kel * central - kel_m8 * central_m8

    # 5. Observation variables. Concentrations are compartment amounts
    # divided by apparent volumes. Because vc_m8 is FIXED to 1 L per
    # the source paper's identifiability normalisation, the M8
    # "concentration" Cc_m8 is numerically equal to the central_m8
    # compartment amount in mg.
    Cc    <- central    / vc
    Cc_m8 <- central_m8 / vc_m8

    Cc    ~ add(addSd)
    Cc_m8 ~ add(addSd_m8)
  })
}
