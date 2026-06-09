# Population pharmacokinetic model for oral lumefantrine in pregnant and
# non-pregnant Tanzanian women with uncomplicated Plasmodium falciparum
# malaria (Mosha 2014, Antimicrobial Agents and Chemotherapy 58(8):
# 4583-4592; doi:10.1128/AAC.02595-14).

Mosha_2014_lumefantrine <- function() {
  description <- paste(
    "Population PK model for oral lumefantrine in 33 pregnant (2nd or 3rd",
    "trimester) and 22 non-pregnant women with uncomplicated Plasmodium",
    "falciparum malaria in Rufiji, Tanzania after standard fixed-dose",
    "artemether-lumefantrine (Mosha 2014). One-compartment disposition",
    "with first-order absorption and ka fixed at 0.54 1/h. Relative",
    "bioavailability F1 is fixed at 1 (structural anchor) with a",
    "categorical pregnancy effect of -33% on F1 (linear-deviation form)",
    "and log-normal IIV around the typical F1 (65% CV). The published",
    "model used a logit transformation on individual F1 to constrain",
    "individuals to (0, 1); this encoding uses log-normal IIV on F",
    "(matching the established Kloprogge 2013 / 2018 lumefantrine",
    "precedents in nlmixr2lib). Structural CL and Vc do not carry IIV in",
    "the final model; the F1 IIV absorbs the joint CL/Vc variability via",
    "the AUC = D x F / CL relationship. The desbutyl-lumefantrine (DLF)",
    "metabolite arm of the published joint model is not encoded; see the",
    "validation vignette for the rationale.",
    sep = " "
  )
  reference <- paste(
    "Mosha D, Guidi M, Mwingira F, Abdulla S, Mercier T, Decosterd LA,",
    "Csajka C, Genton B (2014).",
    "Population pharmacokinetics and clinical response for",
    "artemether-lumefantrine in pregnant and nonpregnant women with",
    "uncomplicated Plasmodium falciparum malaria in Tanzania.",
    "Antimicrobial Agents and Chemotherapy 58(8):4583-4592.",
    "doi:10.1128/AAC.02595-14.",
    sep = " "
  )
  vignette <- "Mosha_2014_artemether_lumefantrine"
  units    <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    PREG = list(
      description        = "Pregnancy status indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = pregnant (second or third trimester), 0 = non-pregnant.",
        "Time-fixed per subject. Mosha 2014 enrolled 33 pregnant and 22",
        "non-pregnant women (Table 1). The paper applies pregnancy as a",
        "linear-deviation effect on relative bioavailability F1:",
        "F1_typical = TVF1 * (1 + e_preg_f * PREG) with TVF1 = 1 and",
        "e_preg_f = -0.33, so pregnant women have a 33% lower typical",
        "F1 than non-pregnant women (the published model also retains a",
        "diarrhoea effect on F1 of -84%, but only 2 women in the cohort",
        "had diarrhoea and the diarrhoea covariate is not encoded here",
        "-- see vignette Assumptions and deviations). The paper",
        "additionally retains a pregnancy effect on the LF-to-DLF",
        "formation rate K23 (+80%); since the DLF compartment is not",
        "encoded in this model, that K23 covariate is also not carried",
        "(vignette Errata).",
        sep = " "
      ),
      source_name        = "PREG"
    )
  )

  covariatesDataExcluded <- list(
    DIARR = list(
      description = "Acute diarrhoea on treatment indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Mosha 2014 reports a categorical diarrhoea effect on F1 of",
        "-84% (linear deviation, theta_diarrF1 = -0.84 with 95% CI",
        "-0.95 to -0.44). Only 2 women in the cohort had diarrhoea on",
        "treatment, all pregnant; the effect is statistically",
        "significant (Delta-OFV = -15, p = 1.1e-4) but rests on a very",
        "small subsample. The canonical covariate-column register",
        "(inst/references/covariate-columns.md) does not currently",
        "list DIARR, and the operator stop-and-ask gate prevents this",
        "skill from introducing a new canonical without prior",
        "registration. The diarrhoea effect is therefore not encoded",
        "in the model file; it is documented here in",
        "covariatesDataExcluded so the paper's covariate-screen",
        "provenance is preserved.",
        sep = " "
      )
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 55L,
    n_studies       = 1L,
    n_pregnant      = 33L,
    n_nonpregnant   = 22L,
    n_observations_lf = 265L,
    age_range       = "18-41 years (pregnant median 25; non-pregnant median 21.5; Table 1)",
    weight_range    = "40-80 kg (pregnant median 52; non-pregnant median 48.5; Table 1)",
    sex_female_pct  = 100,
    disease_state   = paste(
      "Uncomplicated Plasmodium falciparum malaria; pregnant women in",
      "the second or third trimester (median gestational age 27 weeks,",
      "range 14-37) compared concurrently with non-pregnant women",
      "(controls) recruited from the same clinic during the same period."
    ),
    dose_range      = paste(
      "Coartem (Novartis): 20 mg artemether + 120 mg lumefantrine per",
      "tablet; four tablets per dose (480 mg lumefantrine) given twice",
      "daily for 3 days (oral; dose times 0, 8, 24, 36, 48, 60 hours)",
      "with 200 mL milk containing 4.5 g fat to optimise lumefantrine",
      "bioavailability."
    ),
    regions         = "Tanzania (Rufiji district, Kibiti health center; April-September 2012)",
    notes           = paste(
      "Demographics from Mosha 2014 Table 1. 89.5% (51 / 57) of recruited",
      "patients had detectable residual lumefantrine at baseline (mean",
      "37.3 ng/mL); the published joint model estimates a residual",
      "previous-treatment dose F0 = 2.7 mg with 87% CV IIV (Table 3,",
      "text states 3.2 mg). F0 is a study-specific artefact representing",
      "lumefantrine carry-over from prior antimalarial use and is not",
      "encoded in this forward-simulation model; users who wish to",
      "reproduce the published baseline-positive cohort should pre-load",
      "the depot. Companion artemether extraction from the same cohort:",
      "modellib('Mosha_2014_artemether')."
    )
  )

  ini({
    # Structural population mean parameters from Mosha 2014 Table 3
    # "Lumefantrine" section, "Population pharmacokinetics analysis"
    # column. The paper reports clearance and volume estimates on the
    # linear scale; log() is applied here for the nlmixr2 internal scale.
    lcl     <- log(2.8)      ; label("Apparent elimination clearance, CL/F (L/h)")        # Mosha 2014 Table 3: CL = 2.8 (RSE 12%; bootstrap 2.8, 95% CI 2.2-3.6)
    lvc     <- log(134)      ; label("Apparent central volume of distribution, Vc/F (L)") # Mosha 2014 Table 3: Vc = 134 (RSE 14%; bootstrap 134, 95% CI 101-174)
    lka     <- fixed(log(0.54)) ; label("Absorption rate constant, ka (1/h); fixed")       # Mosha 2014 Table 3: Ka fixed at 0.54 1/h (Methods: "the absorption rate constants ... were thus fixed ... to achieve peak plasma AM and LF concentrations, respectively, 2 h and 6 to 8 h after drug intake")

    # Relative bioavailability anchored at 1 (structural, fixed by the
    # source paper; Table 3 reports F1 = "Fixed at 1"). All structural
    # variability around F1 is captured via the etalfdepot IIV (65% CV),
    # which the paper attributes to the joint correlation between CL and
    # Vc absorbed into F1 (Results: "The assignment of interpatient
    # variability on LF bioavailability F1 (fixed at 1) accounting for
    # the correlation between CL and Vc and their variability resulted
    # in additional improvement of the model fit").
    lfdepot <- fixed(log(1)) ; label("Relative bioavailability F1 (unitless); fixed at 1") # Mosha 2014 Table 3: F1 = "Fixed at 1"

    # Pregnancy covariate effect on F1, applied as a linear-deviation
    # form F1_typical = TVF1 * (1 + e_preg_f * PREG) per the paper's
    # generic covariate-model equation theta = theta_a * (1 + theta_b * X)
    # (Methods, Covariate model). Theta_PregF1 = -0.33 corresponds to a
    # 33% reduction in typical F1 in pregnant women relative to the
    # non-pregnant reference.
    e_preg_f <- -0.33 ; label("Pregnancy effect on F1: F1_pregnant / F1_nonpregnant - 1 = -0.33") # Mosha 2014 Table 3: theta_PregF1 = -0.33 (RSE 37%; bootstrap -0.31, 95% CI -0.52 to -0.05)

    # IIV. Mosha 2014 Table 3 column "IIV (%)" reports the F1 IIV as
    # 65% CV. The internal variance for log-normal IIV is recovered by
    # omega^2 = log((CV/100)^2 + 1) = log(0.65^2 + 1) = 0.35272 -- the
    # same convention used by Kloprogge 2013 / 2018. The published model
    # parameterised individual F1 on a logit scale (logit(F1) = logit(
    # F1_typical) + eta) to constrain individuals to (0, 1); this
    # encoding uses log-normal IIV on F directly, matching the
    # established lumefantrine precedents in nlmixr2lib. See vignette
    # Assumptions and deviations.
    etalfdepot ~ 0.35272  # Mosha 2014 Table 3: BSV on F1 = 65% CV (RSE 50%; bootstrap 61%, 95% CI 43-77)

    # Residual error. Mosha 2014 retained a proportional residual error
    # model for LF on the linear-concentration scale (Results,
    # "Residual intrapatient variability was best described using a
    # proportional and mixed-error model for LF and DLF, respectively").
    # Table 3 reports sigma_prop,LF = 51% CV. The proportional error
    # parameter propSd is therefore the fractional SD on the linear
    # concentration scale.
    propSd  <- 0.51 ; label("Proportional residual SD for LF plasma concentration (fraction)") # Mosha 2014 Table 3: sigma_prop,LF = 51% CV (RSE 32%; bootstrap 51%, 95% CI 45-56)
  })

  model({
    # Individual PK parameters. Mosha 2014 retained no body-weight,
    # body-mass-index, age, or gestational-age covariates in the final
    # LF model (Results, "None of the remaining covariates influenced
    # LF and DLF pharmacokinetics") so no allometric scaling is applied.
    # CL and Vc do not carry IIV; the F1 IIV (lognormal, on the depot
    # compartment) captures the joint CL/Vc variability.
    ka       <- exp(lka)
    cl       <- exp(lcl)
    vc       <- exp(lvc)

    # Typical relative bioavailability with the linear-deviation
    # pregnancy effect on F1.
    fdepot_typ <- exp(lfdepot) * (1 + e_preg_f * PREG)

    # Elimination rate constant.
    kel <- cl / vc

    # ODE system: first-order absorption from depot into a
    # one-compartment disposition.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Relative bioavailability F applied to the depot compartment. The
    # log-normal IIV (etalfdepot) is layered on top of the typical
    # covariate-adjusted F. See ini() comment regarding the logit ->
    # lognormal encoding deviation from the source paper.
    f(depot) <- fdepot_typ * exp(etalfdepot)

    # Lumefantrine plasma concentration in ng/mL. With doses in mg and
    # Vc in L, central / vc has units mg/L = ug/mL; multiply by 1000 to
    # obtain ng/mL, which matches the units the paper reports (Figures
    # 2B and 4, day-7 threshold values of 175, 280, 600 ng/mL).
    Cc <- 1000 * central / vc

    # Proportional residual error on the linear-concentration scale.
    Cc ~ prop(propSd)
  })
}
