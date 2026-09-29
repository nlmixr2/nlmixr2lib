Gidal_2018_eslicarbazepine_weekly_seizures <- function() {
  description <- paste0(
    "Zero-inflated Poisson exposure-efficacy model for the mean WEEKLY ",
    "SEIZURE COUNT over up to 18 weeks of treatment in adults with ",
    "focal-onset seizures taking adjunctive eslicarbazepine acetate ",
    "(ESL) (Gidal 2018, phase 3 trials 2093-301, 2093-302 and 2093-304). ",
    "Unlike the paper's other efficacy models this one carries a TIME ",
    "COURSE: the mean count lambda is a baseline count multiplied by ",
    "[1 - 0.560 * ((1 - 0.390) * DRUG + 0.390 * week / 23)], where ",
    "DRUG = Cav-ss / (Cav-ss + 9450) and the baseline is ",
    "2.17 + 0.751*WesternEurope + 0.958*LatinAmerica ",
    "+ 1.12*NorthAmerica - 0.0243*(age - 37), scaled by exp(eta) ",
    "(Gidal 2018 Appendix S1 Eqs. E-11, E-12 and E-13, Table S-9). The ",
    "maximum achievable reduction from baseline is 56%, of which a ",
    "time-driven placebo component accounts for 39% and eslicarbazepine ",
    "exposure for the remaining 61%. The EC50 of 9,450 ng/mL is the ",
    "value quoted in the main text as 9.5 ug/mL and described there as ",
    "similar to the median Cav-ss on ESL 800 mg once daily, so about ",
    "half the maximal drug effect is reached at that dose. The time term ",
    "rises LINEARLY in week number scaled by 23, the maximum week in the ",
    "dataset, and is not saturating; it is only valid inside that ",
    "window. Rest of World is the region reference. A zero-inflation ",
    "factor of 0.0809 was estimated and is carried as a parameter so the ",
    "marginal expected count can be derived, since rxode2 has no ",
    "zero-inflated Poisson likelihood. There is no PK layer and no ODE: ",
    "exposure enters as the column CAV, an empirical-Bayes prediction ",
    "from modellib('Gidal_2018_eslicarbazepine'). Minimum objective ",
    "function 101,991.667."
  )
  reference <- paste(
    "Gidal BE, Jacobson MP, Ben-Menachem E, Carreno M, Blum D,",
    "Soares-da-Silva P, Falcao A, Rocha F, Moreira J, Grinnell T,",
    "Ludwig E, Fiedler-Kelly J, Passarell J, Sunkaraneni S.",
    "Exposure-safety and efficacy response relationships and population",
    "pharmacokinetics of eslicarbazepine acetate.",
    "Acta Neurol Scand. 2018;138(3):203-211. doi:10.1111/ane.12950.",
    "Parameter table and equations are in Appendix S1 (supporting",
    "information), Table S-9 and Equations E-11, E-12 and E-13.",
    "Exposure metric produced by modellib('Gidal_2018_eslicarbazepine').",
    sep = " "
  )
  vignette <- "Gidal_2018_eslicarbazepine_exposure_response"
  units <- list(
    time = "week",
    dosing = "n/a (no dose events; exposure enters as the covariate CAV)",
    concentration = "seizweek (mean number of seizures per week, the Poisson rate lambda; the zero-inflation-adjusted marginal mean seizweek_marg is derived)"
  )

  covariateData <- list(
    CAV = list(
      description = paste(
        "Individual predicted average steady-state eslicarbazepine plasma",
        "concentration over the once-daily dosing interval, for the week",
        "in question."
      ),
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Empirical-Bayes prediction from the population PK model of the",
        "same paper; compute it as dose / (24 * CL/F) with",
        "modellib('Gidal_2018_eslicarbazepine'), equivalently AUC0-24 / 24.",
        "Gidal 2018 describes it as the WEEKLY Cav-ss in the ith patient",
        "at the jth time, so it can vary across the titration period",
        "before settling at the maintenance value. Set to 0 for placebo",
        "patients and during the pre-treatment period: the drug term",
        "Cav-ss / (Cav-ss + 9450) is then exactly 0 and Eq. E-12 collapses",
        "to the placebo Eq. E-11, which is why the two printed equations",
        "can be encoded as one. NOT centred and NOT scaled -- it enters",
        "raw into the Emax denominator."
      ),
      source_name = "C_av-ss_ij (Eqs. E-12 and E-13)"
    ),
    PLACEBO = list(
      description = "Randomised placebo-arm membership; 1 = placebo, 0 = active eslicarbazepine acetate.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (active treatment arm)",
      notes = paste(
        "Gidal 2018 prints Eq. E-11 for placebo patients and Eq. E-12 for",
        "active patients; the two differ ONLY by the presence of the DRUG",
        "term. This model gates that term with (1 - PLACEBO) so both",
        "printed equations are reproduced by a single expression. Note",
        "that this is a WEAKER switch than the one in the companion",
        "standardized-seizure-frequency model: here the time-driven",
        "placebo component applies to BOTH arms (it is the shared",
        "39% share of the maximum effect), whereas in the SSF model the",
        "placebo shift and the drug effect are mutually exclusive.",
        "Setting CAV = 0 for a placebo patient makes PLACEBO redundant,",
        "but both are kept so the encoding matches the printed equations",
        "literally."
      ),
      source_name = "implicit in the split between Eqs. E-11 and E-12"
    ),
    REGION_WESTERNEUROPE = list(
      description = "Western European study site; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      reference_category = "0; with REGION_LATINAMERICA and REGION_NORTHAMERICA also 0 this selects the REST OF WORLD reference group",
      notes = paste(
        "Acts on the baseline weekly seizure count only (+0.751 seizures",
        "per week). The three region indicators are mutually exclusive;",
        "all three 0 selects Rest of World. The region split here matches",
        "the companion standardized-seizure-frequency model and differs",
        "from the splits used by the probability-of-response model",
        "(Western Europe versus pooled remainder) and by the TEAE models",
        "(Europe as reference)."
      ),
      source_name = "WEU_i (Eqs. E-11 and E-12)"
    ),
    REGION_LATINAMERICA = list(
      description = "Latin American study site; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      reference_category = "0; see REGION_WESTERNEUROPE for the shared Rest-of-World reference",
      notes = "Acts on the baseline weekly seizure count only (+0.958 seizures per week). Mutually exclusive with the other two region indicators.",
      source_name = "LA_i (Eqs. E-11 and E-12)"
    ),
    REGION_NORTHAMERICA = list(
      description = "North American study site; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      reference_category = "0; see REGION_WESTERNEUROPE for the shared Rest-of-World reference",
      notes = paste(
        "Acts on the baseline weekly seizure count only (+1.12 seizures",
        "per week), the largest of the three regional shifts. Study",
        "2093-304 was the North American trial. Mutually exclusive with",
        "the other two region indicators."
      ),
      source_name = "NA_i (Eqs. E-11 and E-12)"
    ),
    AGE = list(
      description = "Patient age.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters the baseline count linearly, centred at the population",
        "median of 37 years; the negative slope means older patients are",
        "predicted to have fewer seizures per week at baseline. The",
        "figures of Gidal 2018 are all drawn at this median age. The",
        "slope (-0.0243 seizures per week per year) is on the natural",
        "count scale here, whereas the companion standardized-seizure-",
        "frequency model's age slope (-0.00922) is on the log scale; the",
        "two are not interchangeable."
      ),
      source_name = "AGE_i (Eqs. E-11 and E-12)"
    )
  )

  population <- list(
    species = "human",
    n_studies = 3L,
    n_subjects = 1152L,
    age_median = "37 years (the centring value used by Eqs. E-11 and E-12)",
    disease_state = "adults with focal-onset seizures, at least four in the 4 weeks before screening despite 1-3 concomitant antiepileptic drugs",
    dose_range = "eslicarbazepine acetate 400, 800 or 1,200 mg orally once daily, or placebo",
    regions = "Rest of World (reference), Western Europe, Latin America and North America",
    observation_window = "weekly epochs for up to 18 weeks of treatment; the week index NMWK runs to a maximum of 23 in the dataset and is the value the time term is scaled by",
    notes = paste(
      "Weekly seizure counts were modelled by Poisson regression within",
      "weekly epochs, following Ette & Williams (2007) and the pregabalin",
      "exposure-response analysis of Miller 2003. A posterior predictive",
      "check of the percentage of responders at week 14 is shown in",
      "Figure S-3 for each of the four dose groups. The number of subjects",
      "is the safety analysis set size, which Appendix S1 does not restate",
      "separately for the weekly-seizure-count analysis."
    )
  )

  ini({
    # ==================================================================
    # Gidal 2018 Appendix S1 Table S-9 and Equations E-11, E-12, E-13:
    #
    #   lambda_ij = (2.17 + 0.751*WEU + 0.958*LA + 1.12*NA
    #                     - 0.0243*(AGE - 37)) * exp(eta_i)
    #               * [1 - 0.56 * ((1 - 0.39)*DRUG_ij
    #                              + 0.39 * NMWK_ij / 23)]
    #   DRUG_ij = Cav-ss_ij / (Cav-ss_ij + 9450)
    #
    # Eq. E-11 (placebo) is the same expression with the DRUG term
    # absent, which is what (1 - PLACEBO) reproduces.
    # ==================================================================

    # ----- Baseline weekly seizure count (additive covariate shifts) ---
    lrbase <- log(2.17) ; label("Mean baseline seizure frequency for a 37-year-old Rest-of-World patient (seizures/week)")  # Table S-9, 'Mean baseline seizure frequency (per week)' 2.17, 3.6% SEM; Eqs. E-11 and E-12
    e_region_westerneurope_rbase <- 0.751 ; label("Additive shift in the baseline weekly seizure count for a Western European versus a Rest-of-World study site (seizures/week)")  # Table S-9, 'Additive shift for Western Europe on baseline seizure frequency' 0.751, 26.5% SEM
    e_region_latinamerica_rbase <- 0.958 ; label("Additive shift in the baseline weekly seizure count for a Latin American versus a Rest-of-World study site (seizures/week)")  # Table S-9, 'Additive shift for Latin America on baseline seizure frequency' 0.958, 18.8% SEM
    e_region_northamerica_rbase <- 1.12 ; label("Additive shift in the baseline weekly seizure count for a North American versus a Rest-of-World study site (seizures/week)")  # Table S-9, 'Additive shift for North America on baseline seizure frequency' 1.12, 18.1% SEM
    e_age_rbase <- -0.0243 ; label("Change in the baseline weekly seizure count per year of age above 37 years (seizures/week/year)")  # Table S-9, 'Slope for age on baseline seizure frequency (seizures per week per year of age)' -0.0243, 20.6% SEM

    # ----- Maximum effect and its time / exposure split -----
    emax <- 0.560 ; label("Maximum fractional reduction in the weekly seizure count attributable to time and eslicarbazepine exposure combined (unitless fraction)")  # Table S-9, 'Emax due to time and Cav-ss (maximum fractional reduction in weekly seizure frequency)' 0.560, 8.9% SEM. Matches the main text's 'maximum reduction from baseline of 56% during treatment with ESL'
    f_time <- 0.390 ; label("Share of the maximum effect attributable to the time-driven placebo component; the remainder is attributable to eslicarbazepine exposure (unitless fraction)")  # Table S-9, 'Fraction of maximum effect due to time' 0.390, 10.7% SEM. Matches the main text's '39% ... and eslicarbazepine Cav-ss accounted for the remaining 61%'
    lec50 <- log(9450) ; label("Average steady-state eslicarbazepine concentration giving half the maximum drug component of the effect (ng/mL)")  # Table S-9, 'EC50 (ng/mL)' 9,450, 31.1% SEM; Eq. E-13. This is the value the main text quotes as 9.5 ug/mL -- NOT the 3,530 ng/mL EC50 of the companion standardized-seizure-frequency model

    # ----- Zero inflation -----
    f_zeroinfl <- 0.0809 ; label("Probability that an observation is a structural zero, outside the Poisson count process (unitless fraction)")  # Table S-9, 'Zero-inflation factor' 0.0809, 3.5% SEM. Carried as a parameter because rxode2 has no zero-inflated Poisson likelihood; see model() and the vignette Errata item 5

    # ----- Interindividual variability -----
    # Table S-9 reports the baseline IIV as a VARIANCE (footnote a);
    # sqrt(0.654) = 0.809, the printed 80.9 %CV.
    etalrbase ~ 0.654  # Table S-9, baseline seizure frequency IIV 0.654, 4.4% SEM; footnote a: 'The estimate provided in the table (0.654) is a variance term. The corresponding %CV = 80.9%'

    # ----- No residual error in the source -----
    # The source likelihood is zero-inflated Poisson, which has no
    # residual-error parameter. The tiny fixed additive residual below
    # exists only so rxode2 has an error model to attach to the
    # typical-value count; it is NOT a published quantity. See the
    # vignette's Assumptions and deviations.
    addSd_seizweek <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value weekly seizure count; the source likelihood is zero-inflated Poisson (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ----- Constants -----
    # 23 is the maximum week index NMWK observed in the dataset; Gidal
    # 2018 scales the time term by it 'to prevent negative predictions'.
    # The term is linear, not saturating, so the model is only valid for
    # week <= 23.
    nmwk_max <- 23

    # ----- Baseline weekly seizure count (Eqs. E-11 and E-12) -----
    # The covariate shifts are absolute seizures per week added to the
    # 2.17 typical value; the exponential IIV then scales the sum.
    seizweek_base <- (exp(lrbase) +
      e_region_westerneurope_rbase * REGION_WESTERNEUROPE +
      e_region_latinamerica_rbase * REGION_LATINAMERICA +
      e_region_northamerica_rbase * REGION_NORTHAMERICA +
      e_age_rbase * (AGE - 37)) *
      exp(etalrbase)

    # ----- Drug component (Eq. E-13) -----
    ec50 <- exp(lec50)
    drug_eff <- CAV / (CAV + ec50)

    # ----- Time-driven placebo component, common to both arms ---------
    time_eff <- time / nmwk_max

    # ----- Mean weekly seizure count -----
    # (1 - PLACEBO) removes the drug term for placebo patients, turning
    # Eq. E-12 into Eq. E-11.
    seizweek <- seizweek_base *
      (1 - emax * ((1 - f_time) * drug_eff * (1 - PLACEBO) +
        f_time * time_eff))

    # ----- Marginal expected count under the zero-inflation mixture ----
    # seizweek above is the Poisson rate of the count component; a
    # fraction f_zeroinfl of observations are structural zeros, so the
    # marginal expectation is (1 - f_zeroinfl) * lambda. Which of the two
    # the paper's figures plot is not stated; see the vignette Errata.
    seizweek_marg <- (1 - f_zeroinfl) * seizweek

    # ----- Observation -----
    seizweek ~ add(addSd_seizweek)
  })
}
