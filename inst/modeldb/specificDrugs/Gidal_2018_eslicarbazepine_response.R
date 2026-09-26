Gidal_2018_eslicarbazepine_response <- function() {
  description <- paste0(
    "Logistic-regression exposure-efficacy model for the PROBABILITY OF ",
    "RESPONSE (at least a 50% reduction in seizure frequency from ",
    "baseline during the maintenance phase) in adults with focal-onset ",
    "seizures taking adjunctive eslicarbazepine acetate (ESL) (Gidal ",
    "2018, phase 3 trials 2093-301, 2093-302 and 2093-304). The linear ",
    "predictor is -1.30 + 0.735 * (Cav-ss / 10205)^0.609 ",
    "- 0.668*WesternEurope (Gidal 2018 Appendix S1 Eq. E-10, ",
    "Table S-8). The exposure term is a POWER function of the average ",
    "steady-state eslicarbazepine concentration, centred at the ",
    "population median of 10,205 ng/mL; a constant effect, a linear ",
    "effect, a log-linear effect and an Emax effect were all screened ",
    "and the power form was selected. Because the exponent 0.609 is ",
    "below 1 the exposure-response curve is concave: most of the benefit ",
    "is realised well below the median concentration and the curve is ",
    "shallow above it, which is the basis for the paper's conclusion ",
    "that routine plasma-concentration monitoring is not useful for ",
    "guiding ESL dose selection. Western European patients carry a ",
    "0.668 lower logit at every exposure. All eight published predicted ",
    "probabilities -- placebo, 400, 800 and 1,200 mg in each of the two ",
    "region groups -- reproduce from this equation to two decimal ",
    "places; see the vignette source trace. There is no PK layer and no ",
    "ODE: exposure enters as the static per-patient column CAV, an ",
    "empirical-Bayes prediction from ",
    "modellib('Gidal_2018_eslicarbazepine'). No between-subject random ",
    "effect and no residual error are estimated (Bernoulli likelihood). ",
    "Hosmer-Lemeshow chi-squared 6.03 on 8 df (p = 0.6442), area under ",
    "the ROC curve 0.62, minimum objective function 1,345.498."
  )
  reference <- paste(
    "Gidal BE, Jacobson MP, Ben-Menachem E, Carreno M, Blum D,",
    "Soares-da-Silva P, Falcao A, Rocha F, Moreira J, Grinnell T,",
    "Ludwig E, Fiedler-Kelly J, Passarell J, Sunkaraneni S.",
    "Exposure-safety and efficacy response relationships and population",
    "pharmacokinetics of eslicarbazepine acetate.",
    "Acta Neurol Scand. 2018;138(3):203-211. doi:10.1111/ane.12950.",
    "Parameter table and equation are in Appendix S1 (supporting",
    "information), Table S-8 and Equation E-10.",
    "Exposure metric produced by modellib('Gidal_2018_eslicarbazepine').",
    sep = " "
  )
  vignette <- "Gidal_2018_eslicarbazepine_exposure_response"
  units <- list(
    time = "n/a (static landmark responder analysis over the maintenance phase; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the covariate CAV)",
    concentration = "prob_response (probability of at least a 50% reduction in seizure frequency, 0-1; also logit_response)"
  )

  covariateData <- list(
    CAV = list(
      description = paste(
        "Individual predicted average steady-state eslicarbazepine plasma",
        "concentration over the once-daily dosing interval."
      ),
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Empirical-Bayes prediction from the population PK model of the",
        "same paper; compute it as dose / (24 * CL/F) with",
        "modellib('Gidal_2018_eslicarbazepine'), equivalently AUC0-24 / 24.",
        "CENTRED at the population median of 10,205 ng/mL, which is the",
        "median for the pooled active arms and -- per the main text --",
        "close to the median on ESL 800 mg once daily. This is the single",
        "most useful anchor in the paper for checking exposure units, and",
        "it is consistent with the PK model: 800 mg / (24 h * 2.43 L/h)",
        "= 13,717 ng/mL for a patient on no comedication, falling to about",
        "9,500 ng/mL for the roughly half of the population taking",
        "carbamazepine 800 mg/day. Set to 0 for placebo patients; the",
        "power term is then exactly 0 and the intercept is the placebo",
        "logit, which reproduces the published placebo probabilities of",
        "0.21 (non-Western-Europe) and 0.12 (Western Europe)."
      ),
      source_name = "C_av-ss_i (Eq. E-10)"
    ),
    REGION_WESTERNEUROPE = list(
      description = "Western European study site; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all non-Western-European sites pooled: North America, Latin America and Rest of World)",
      notes = paste(
        "The only covariate retained in this model, and the reference",
        "group is the pooled remainder rather than a single named region",
        "-- a two-way split, unlike the four-way split the companion",
        "standardized-seizure-frequency and weekly-seizure-count models",
        "use and unlike the Europe-referenced split the TEAE models use.",
        "Do not transfer a region encoding between models of this family.",
        "The negative coefficient predicts a lower response probability in",
        "Western Europe at every exposure; Gidal 2018 Discussion",
        "attributes this to demographic and clinical differences between",
        "the groups and notes that the study was not designed to evaluate",
        "it. The same directional finding appears independently in the",
        "standardized-seizure-frequency model, where Western Europe",
        "shrinks the Emax."
      ),
      source_name = "WEU_i (Eq. E-10)"
    )
  )

  population <- list(
    species = "human",
    n_studies = 3L,
    n_subjects = 1152L,
    age_median = "37 years",
    disease_state = "adults with focal-onset seizures, at least four in the 4 weeks before screening despite 1-3 concomitant antiepileptic drugs",
    dose_range = "eslicarbazepine acetate 400, 800 or 1,200 mg orally once daily, or placebo",
    regions = "Western Europe versus the pooled non-Western-European remainder",
    median_cavss = "10,205 ng/mL (the centring value used by Eq. E-10)",
    notes = paste(
      "The responder endpoint is dichotomous: a patient is a responder",
      "when the number of seizures fell by at least 50% from baseline",
      "during the maintenance phase. Published predicted probabilities --",
      "Western Europe 0.12 / 0.18 / 0.22 / 0.26 and non-Western-Europe",
      "0.21 / 0.30 / 0.35 / 0.40 for placebo / 400 / 800 / 1,200 mg --",
      "all reproduce from Eq. E-10. The paper notes that this model agreed",
      "closely with the separately fitted standardized-seizure-frequency",
      "model. The number of subjects is the safety analysis set size,",
      "which Appendix S1 does not restate separately for the responder",
      "analysis."
    )
  )

  ini({
    # ==================================================================
    # Gidal 2018 Appendix S1 Table S-8 and Equation E-10:
    #
    #   logit(p) = -1.30 + [0.735 * (Cav-ss / 10205)^0.609]
    #              - 0.668*WEU
    #
    # Four rows, reported once in Table S-8 (footnote a: 'Parameter
    # estimates on the logit scale') and once in Eq. E-10; both agree
    # exactly. Cav-ss is scaled by its median 10,205 ng/mL but NOT
    # subtracted, so the intercept is the logit at zero exposure, i.e.
    # the placebo logit.
    # ==================================================================

    logit_ref <- -1.30 ; label("Logit of the probability of response for a non-Western-European patient at zero eslicarbazepine exposure (unitless logit)")  # Table S-8, 'Placebo effect' -1.30, 9.2% SEM; Eq. E-10. Check: expit(-1.30) = 0.214, the 0.21 non-Western-Europe placebo probability the paper reports
    e_cav_logit <- 0.735 ; label("Log-odds gain in the probability of response at the median average steady-state eslicarbazepine concentration of 10,205 ng/mL (unitless logit)")  # Table S-8, 'Intercept for the eslicarbazepine effect' 0.735, 19.5% SEM; Eq. E-10
    e_cav_logit_pow <- 0.609 ; label("Power exponent on (average steady-state eslicarbazepine concentration / 10,205 ng/mL) in the exposure term (unitless)")  # Table S-8, 'Power for eslicarbazepine effect' 0.609, 31.5% SEM; Eq. E-10. Below 1, so the exposure-response curve is concave
    e_region_westerneurope_logit <- -0.668 ; label("Log-odds shift for a Western European versus a non-Western-European study site (unitless logit)")  # Table S-8, 'Effect of Western European region' -0.668, 32.9% SEM; Eq. E-10. Check: expit(-1.30 - 0.668) = 0.123, the 0.12 Western-Europe placebo probability the paper reports

    # ----- No between-subject variability, no residual error -----
    # Bernoulli likelihood: the source estimates no sigma and no random
    # effects. The tiny fixed additive residual exists only so rxode2 has
    # an error model to attach to the typical-value probability.
    addSd_prob_response <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value response probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ----- Linear predictor (Gidal 2018 Eq. E-10) -----
    # CAV is scaled by its median but not centred, so the exposure term
    # is exactly 0 at zero exposure and the intercept is the placebo
    # logit.
    logit_response <- logit_ref +
      e_cav_logit * (CAV / 10205)^e_cav_logit_pow +
      e_region_westerneurope_logit * REGION_WESTERNEUROPE

    prob_response <- expit(logit_response)

    # ----- Observation -----
    prob_response ~ add(addSd_prob_response)
  })
}
