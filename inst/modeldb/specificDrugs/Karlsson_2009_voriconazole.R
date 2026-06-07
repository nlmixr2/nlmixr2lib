Karlsson_2009_voriconazole <- function() {
  description <- paste(
    "Two-compartment population pharmacokinetic model with Michaelis-Menten",
    "elimination for voriconazole in pediatric patients aged 2 to <12 years",
    "(Karlsson 2009), pooled from three open-label intravenous and oral",
    "studies; first-order oral absorption with bioavailability, no lag time;",
    "all disposition parameters proportional to body weight; CYP2C19",
    "metabolizer status (heterozygous extensive metabolizers pooled with poor",
    "metabolizers) and alanine aminotransferase as covariates on clearance;",
    "residual error stratified by CYP2C19 metabolizer group"
  )
  reference <- paste(
    "Karlsson MO, Lutsar I, Milligan PA. Population pharmacokinetic analysis",
    "of voriconazole plasma concentration data from pediatric studies.",
    "Antimicrob Agents Chemother. 2009 Mar;53(3):935-944.",
    "doi:10.1128/AAC.00751-08"
  )
  vignette <- "Karlsson_2009_voriconazole"
  units    <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Karlsson 2009 Methods: 'The parameters CL, Vc, Q, and Vp were",
        "modeled as directly proportional to weight.' All four disposition",
        "parameters scale linearly with WT (no reference-weight",
        "normalization, no allometric exponent). Cohort weight range was",
        "10.8 to 54.9 kg with an average of 22.8 kg (paper Methods 'Study",
        "data')."
      ),
      source_name        = "WT"
    ),
    ALT = list(
      description        = "Alanine aminotransferase",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power covariate on CL: CL = TVCL * (ALT / 26.5)^(-0.0931). The",
        "reference ALT 26.5 IU/L is the population median used by Karlsson",
        "2009 for the typical subject (paper Table 2 footnote a, page 941:",
        "'Typical subject characteristics were as follows: EM, 26.5 IU/liter",
        "ALT, not taking CYP2C9 inhibitors'). Coefficient -0.0931 is",
        "reported as 'fractional change per unit change of log(ALT)' in",
        "Table 1 (final-model column). Karlsson 2009 used time-varying ALT",
        "rather than baseline ALT, because the time-varying form was more",
        "predictive (Discussion paragraph on ALT covariate); per-subject",
        "baseline ALT was 7-242 IU/L with average 40.7 IU/L (paper",
        "Methods)."
      ),
      source_name        = "ALT"
    ),
    CYP2C19_IM = list(
      description        = "CYP2C19 intermediate-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (EM, UM, PM, or RM); when paired with CYP2C19_PM = 0, reference is the extensive metabolizer (EM) phenotype (Karlsson 2009: 58 of 82 patients)",
      notes              = paste(
        "Karlsson 2009 only had three CYP2C19 poor metabolizers in the",
        "dataset, so the 21 heterozygous extensive metabolizers (HEMs) and",
        "the 3 PMs were pooled into a single 'HEM/PM' group estimated as",
        "-35.5% CL relative to the 58 EMs (paper Methods 'Pharmacokinetic",
        "analysis' and Table 1 final-model row 'Decrease in CL for",
        "HEMs/PMs'). The model file reproduces the pooling by attaching the",
        "identical coefficient (-0.355) to both CYP2C19_IM and CYP2C19_PM.",
        "The same pooled HEM/PM indicator also stratifies the residual",
        "error (propSd_em = 0.573 vs propSd_pmim = 0.299, Table 1",
        "final-model rows 'Residual error (EMs)' and 'Residual error",
        "(HEMs/PMs)')."
      ),
      source_name        = "CYP2C19 == 'HEM'"
    ),
    CYP2C19_PM = list(
      description        = "CYP2C19 poor-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (EM, UM, IM, or RM); when paired with CYP2C19_IM = 0, reference is the extensive metabolizer (EM) phenotype",
      notes              = paste(
        "Karlsson 2009 had only 3 PMs (of 82 total patients) and pooled",
        "them with the 21 HEMs into a single 'HEM/PM' covariate group; the",
        "PM-specific effect was therefore not separately estimated. The",
        "model file applies the same coefficient (-0.355) to CYP2C19_PM as",
        "to CYP2C19_IM to faithfully reproduce the paper's pooled HEM/PM",
        "estimate. The pooled group also stratifies the residual error",
        "(see CYP2C19_IM notes and Table 1 final-model column)."
      ),
      source_name        = "CYP2C19 == 'PM'"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 82L,
    n_studies      = 3L,
    n_observations = 1274L,
    age_range      = "2 to <12 years",
    weight_range   = "10.8-54.9 kg (mean 22.8)",
    weight_mean    = "22.8 kg",
    sex_female_pct = 42.7,
    cyp2c19_dist   = c(EM_pct = 70.7, HEM_pct = 25.6, PM_pct = 3.7),
    race_ethnicity = c(White = 69.5, Black = 7.3, Asian = 4.9, Other = 18.3),
    disease_state  = paste(
      "Pediatric patients aged 2 to <12 years with serious fungal",
      "infections, predominantly immunocompromised due to underlying",
      "hematologic malignancy or bone marrow transplantation. Disease",
      "distribution: 31 leukemia (37.8%), 39 bone marrow transplantation",
      "(47.6%), 2 lymphoma (2.4%), 1 aplastic anemia (1.2%), 9 other",
      "(11.0%). Mucositis was present in 57 patients (69.5%)."
    ),
    dose_range     = paste(
      "Pooled from three studies. Study A (n=11): single 3 or 4 mg/kg i.v.",
      "doses. Study B (n=28): loading 6 mg/kg BID i.v. on day 1, then 3",
      "mg/kg BID days 2-4 and 4 mg/kg BID days 4-8. Study C (n=43, two",
      "cohorts): cohort 1 received 6 mg/kg BID i.v. day 1, 4 mg/kg BID",
      "i.v. days 2-4, 6 mg/kg BID i.v. days 5-8, then 4 mg/kg BID p.o.",
      "days 9-12; cohort 2 received 6 mg/kg BID i.v. days 1-4, 8 mg/kg",
      "BID i.v. days 5-8, then 6 mg/kg BID p.o. days 9-12. i.v. infusions",
      "0.5-3.5 h. Average 15.5 plasma samples per subject."
    ),
    alt_range      = "7-242 IU/L (mean 40.7; population median 26.5 per Table 2 footnote)",
    alp_range      = "46-309 IU/L (mean 135)",
    notes          = paste(
      "Karlsson 2009 Methods: 5 of the original 87 patients were excluded",
      "due to missing CYP2C19 genotype (4 from study B and 1 from study C),",
      "yielding the n=82 analyzed cohort. NONMEM V level 1.1 with FOCE",
      "INTER. Concomitant CYP3A4, CYP2C19, and CYP2C9 inhibitors and CYP450",
      "inducers were tested as covariates; only CYP2C9 inhibitor effect was",
      "borderline-significant and was retained in the alternative final",
      "model but excluded from the final model used for simulations (paper",
      "Results 'Model building' and Discussion paragraph on the alternative",
      "final model)."
    )
  )

  ini({
    # Structural parameters per kg body weight (Karlsson 2009 Table 1,
    # final-model column, page 938-939). All disposition parameters are
    # directly proportional to weight; F and Ka are not weight-scaled.
    lka     <- log(0.849); label("Absorption rate constant (1/h)")                                     # Karlsson 2009 Table 1 final-model: Ka = 0.849 /h
    lcl     <- log(0.582); label("CL per kg in CYP2C19 EMs at reference ALT 26.5 IU/L (L/h/kg)")        # Karlsson 2009 Table 1 final-model: CL in EMs = 0.582 L/h/kg
    lvc     <- log(0.807); label("Central volume per kg (L/kg)")                                       # Karlsson 2009 Table 1 final-model: Vc = 0.807 L/kg
    lq      <- log(0.609); label("Inter-compartmental clearance per kg (L/h/kg)")                      # Karlsson 2009 Table 1 final-model: Q = 0.609 L/h/kg
    lvp     <- log(2.17);  label("Peripheral volume per kg (L/kg)")                                    # Karlsson 2009 Table 1 final-model: Vp = 2.17 L/kg
    lkm     <- log(3.03);  label("Michaelis-Menten constant (ug/mL = mg/L)")                           # Karlsson 2009 Table 1 final-model: Km = 3030 ng/mL = 3.03 ug/mL
    lfdepot <- log(0.446); label("Oral bioavailability (fraction)")                                    # Karlsson 2009 Table 1 final-model: F = 44.6 (%)

    # Covariate effects on CL.
    # CYP2C19 HEM/PM pooled: -35.5% on CL relative to EMs. Karlsson 2009 had
    # 3 PMs in the dataset and merged them with the 21 HEMs into a single
    # 'HEM/PM' group; the model file reproduces the pooling by setting the
    # IM and PM coefficients equal.
    e_cyp2c19_im_cl <- fixed(-0.355); label("CL fractional change for CYP2C19 IM, pooled with PM by Karlsson (unitless)")   # Karlsson 2009 Table 1 final-model: Decrease in CL for HEMs/PMs = 35.5%
    e_cyp2c19_pm_cl <- fixed(-0.355); label("CL fractional change for CYP2C19 PM, pooled with IM by Karlsson (unitless)")   # Karlsson 2009 Table 1 final-model: same coefficient as HEM (pooled)
    # ALT power effect: CL = TVCL * (ALT / 26.5)^(-0.0931). Reference 26.5
    # IU/L is the population median from Table 2 footnote 'a': "Typical
    # subject characteristics were as follows: EM, 26.5 IU/liter ALT, not
    # taking CYP2C9 inhibitors." Coefficient is reported in Table 1 as
    # 'fractional change per unit change of log(ALT)'.
    e_alt_cl <- fixed(-0.0931); label("Power exponent of (ALT / 26.5) on CL (unitless)")               # Karlsson 2009 Table 1 final-model: Log(ALT) on CL = -0.0931

    # IIV. Karlsson 2009 reports correlated IIV on CL, Km, and F (Table 1
    # final-model column 'CV(...)' and 'cor(...)' rows). CVs are converted
    # to log-normal variance using omega^2 = log(CV^2 + 1):
    #   CV(CL) = 52.8% -> omega^2 = log(0.528^2 + 1) = 0.2459
    #   CV(Km) = 131%  -> omega^2 = log(1.31^2  + 1) = 0.9992
    #   CV(F)  = 69.7% -> omega^2 = log(0.697^2 + 1) = 0.3960
    # Covariances derived from reported correlations:
    #   cor(CL,Km) = -0.685 -> cov = -0.685 * sqrt(0.2459 * 0.9992) = -0.3395
    #   cor(CL,F)  =  0.660 -> cov =  0.660 * sqrt(0.2459 * 0.3960) =  0.2059
    #   cor(Km,F)  = -0.646 -> cov = -0.646 * sqrt(0.9992 * 0.3960) = -0.4063
    # Karlsson 2009 also reports interoccasion variability on CL (CV 43%)
    # and an IIV-on-residual-error term (CV 43%) -- both are omitted here
    # because nlmixr2 lacks a clean syntactic analog for either without an
    # explicit OCC column or the eta-on-sigma construct, and the final
    # model used for the paper's simulations remains identifiable with the
    # subject-level CL/Km/F IIV preserved. See vignette Assumptions and
    # deviations.
    etalcl + etalkm + etalfdepot ~ c(
       0.2459,
      -0.3395, 0.9992,
       0.2059, -0.4063, 0.3960
    )                                                                                                  # Karlsson 2009 Table 1 final-model: CV(CL)=52.8, CV(Km)=131, CV(F)=69.7; cor(CL,Km)=-0.685, cor(CL,F)=0.66, cor(Km,F)=-0.646

    # Residual error stratified by CYP2C19 metabolizer status. Karlsson
    # 2009 reports as 'additive on log-concentration', equivalent to
    # proportional on linear scale.
    propSd_em   <- 0.573; label("Proportional residual SD for CYP2C19 EMs (fraction)")                  # Karlsson 2009 Table 1 final-model: Residual error (EMs) = 57.3
    propSd_pmim <- 0.299; label("Proportional residual SD for CYP2C19 HEM/PM (fraction)")               # Karlsson 2009 Table 1 final-model: Residual error (HEMs/PMs) = 29.9
  })

  model({
    # Pooled HEM/PM indicator (Karlsson 2009 pooled IM and PM into a single
    # group). Mutually exclusive in canonical encoding, so the sum is 0 for
    # EMs and 1 for HEM or PM subjects.
    cyp2c19_pmim <- CYP2C19_IM + CYP2C19_PM

    # Individual PK parameters. CL has a power ALT effect (reference 26.5
    # IU/L) and a pooled HEM/PM linear-deviation effect (-35.5%). All
    # disposition parameters scale linearly with body weight (no reference
    # weight normalization, no allometric exponent).
    ka  <- exp(lka)
    cl_em_per_kg <- exp(lcl + etalcl) * (ALT / 26.5)^e_alt_cl
    cl_per_kg    <- cl_em_per_kg *
      (1 + e_cyp2c19_im_cl * CYP2C19_IM + e_cyp2c19_pm_cl * CYP2C19_PM)
    cl  <- cl_per_kg * WT
    vc  <- exp(lvc) * WT
    q   <- exp(lq)  * WT
    vp  <- exp(lvp) * WT
    km  <- exp(lkm + etalkm)
    f_oral <- exp(lfdepot + etalfdepot)

    # Michaelis-Menten elimination. Karlsson 2009 parameterised in CL and
    # Km; Vmax is then derived as Vmax = CL * Km (paper Methods, 'Vmax was
    # obtained as CL multiplied by Km'). Rewriting the MM term in amounts:
    #   Vmax * C / (Km + C) = Vmax * A_central / (Km * Vc + A_central)
    # so the ODE term below is mass-per-time with consistent units.
    vmax <- cl * km

    # 2-compartment with first-order oral absorption and MM elimination.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          vmax * central / (km * vc + central) -
                          q * (central / vc - peripheral1 / vp)
    d/dt(peripheral1) <-  q * (central / vc - peripheral1 / vp)

    # Oral bioavailability (44.6%) applied to depot. I.v. doses bypass
    # depot and enter central directly via the data's `cmt` column.
    f(depot) <- f_oral

    # Plasma voriconazole concentration. dose mg / vc L = mg/L = ug/mL,
    # matching the Karlsson 2009 Km reported in ng/mL.
    Cc <- central / vc

    # CYP2C19-stratified proportional residual error: EMs use propSd_em,
    # HEM/PM cohort uses propSd_pmim.
    prop_strat <- propSd_em * (1 - cyp2c19_pmim) + propSd_pmim * cyp2c19_pmim
    Cc ~ prop(prop_strat)
  })
}
