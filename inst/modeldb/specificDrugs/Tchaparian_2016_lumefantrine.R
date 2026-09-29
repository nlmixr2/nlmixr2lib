# Population PK model of lumefantrine in capillary whole blood of Ugandan
# children aged 6 months to 2 years treated with artemether-lumefantrine for
# uncomplicated Plasmodium falciparum malaria. Source: Tchaparian E, Sambol NC,
# Arinaitwe E, et al. J Infect Dis. 2016;214(8):1243-1251.
# doi:10.1093/infdis/jiw338 (PMC5034953).

Tchaparian_2016_lumefantrine <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption and",
    "linear elimination for oral lumefantrine in capillary whole blood of 101",
    "Ugandan children aged 6 months to 2 years (207 malaria episodes, 806",
    "concentrations) treated with artemether-lumefantrine for uncomplicated",
    "malaria (Tchaparian 2016). Absorption feeds a depot compartment at a rate",
    "constant held at the literature value 0.45 1/h because the study had too",
    "few absorption-phase samples to estimate it. Body-weight allometric",
    "scaling is applied unconditionally to every disposition parameter",
    "(exponent 0.75 shared by CL/F and Q/F, exponent 1 shared by V1/F and",
    "V2/F), referenced to 8.43 kg, the typical weight of a 1-year-old child.",
    "Age is the single retained covariate and acts on relative bioavailability",
    "as the power function F = (AGE_months / 12)^0.596, so a 6-month-old has",
    "F = 0.66 and a 2-year-old F = 1.51 relative to the 1-year-old reference",
    "- an increase with maturation, which the authors note is the opposite of",
    "what CYP3A4 and P-glycoprotein ontogeny would predict and attribute to",
    "immature biliary function limiting the dissolution of this highly",
    "lipophilic drug in the youngest children. Because every disposition term",
    "is apparent (divided by F), correlation among CL/F, V1/F, Q/F and V2/F is",
    "induced through variability on F itself rather than through an estimated",
    "covariance block. Both interindividual and interoccasion variability are",
    "carried on all five parameters, with the occasion being a malaria episode",
    "(up to 5 per child); the Q/F and V2/F variances were constrained to 50",
    "percent CV by the authors because the data were too sparse to estimate",
    "them. The residual is combined additive plus proportional. The fit was a",
    "Monolix 4.3 SAEM run in which the 23 percent of samples below the limit",
    "of detection were retained as left-censored.",
    sep = " "
  )
  reference <- paste(
    "Tchaparian E, Sambol NC, Arinaitwe E, McCormack SA, Bigira V, Wanzira H,",
    "Muhindo M, Creek DJ, Sukumar N, Blessborn D, Tappero JW, Kakuru A,",
    "Bergqvist Y, Aweeka FT, Parikh S.",
    "Population pharmacokinetics and pharmacodynamics of lumefantrine in young",
    "Ugandan children treated with artemether-lumefantrine for uncomplicated",
    "malaria. J Infect Dis. 2016;214(8):1243-1251. doi:10.1093/infdis/jiw338.",
    "Parameter values are from Table 2; the covariate model for F is the",
    "display equation in Results, 'Population Pharmacokinetic Model'; the",
    "statement that the Q/F and V2/F variances were constrained to 50 percent",
    "CV is in the Supplement, 'Population Pharmacokinetic analysis'.",
    sep = " "
  )
  vignette <- "Tchaparian_2016_lumefantrine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE: checked against Tchaparian 2016
  # Methods ('Sample Collection and Analysis': 100 uL of capillary whole blood
  # dried on tartaric-acid-pretreated filter paper) and Table 2, whose
  # parameter names are explicitly capillary whole-blood terms.
  compartmentData <- list(
    depot = list(analyte = "lumefantrine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "lumefantrine", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "lumefantrine", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight recorded at the time of malaria diagnosis; time-fixed within a treatment episode but re-recorded at each new episode.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling was applied unconditionally in every model the",
        "authors fitted, including the base model with no other covariates",
        "(Methods, 'Pharmacokinetic Analysis'). The exponent is 0.75 on the two",
        "clearance terms and 1 on the two volume terms, both held at the",
        "canonical theoretical values rather than estimated. The reference",
        "weight is 8.43 kg, the typical weight of a 1-year-old child; the",
        "authors used the cohort median of 9.0 kg during model development and",
        "switched to 8.43 kg for the final model 'for ease of interpretation'",
        "(Supplement, 'Population Pharmacokinetic analysis'), so the Table 2",
        "point estimates are anchored at 8.43 kg. Cohort median 9.1 kg (range",
        "6.1-13.0 kg; Table 1). All children weighed under 14 kg and therefore",
        "received the same 1-tablet artemether-lumefantrine dose.",
        sep = " "
      ),
      source_name = "WT"
    ),
    AGE = list(
      description = "Subject age at the time of malaria diagnosis, the only covariate retained in the final model; it acts on relative bioavailability.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "AGE is supplied in the canonical unit of years; the model converts to",
        "months internally (AGE * 12) because Tchaparian 2016 writes the",
        "covariate equation on the month scale with a reference age of 12",
        "months: F = (AGE_months / 12)^0.596 (Results, 'Population",
        "Pharmacokinetic Model'). Cohort median 14.4 months, range 6.6-22.2",
        "months (Table 1), so the supported range is roughly 0.55-2.0 years.",
        "The effect is a genuine ontogeny signal that survived allometric",
        "scaling of both clearance and volume, and its direction is positive:",
        "older children absorb more. Extrapolating this power function far",
        "beyond 2 years is not supported by the data.",
        sep = " "
      ),
      source_name = "Age"
    ),
    OCC = list(
      description = "Integer occasion index identifying which malaria episode a record belongs to, used to multiplex the interoccasion-variability random effects.",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "An occasion is one treated malaria episode. Children were eligible for",
        "repeat PK sampling at every episode during longitudinal follow-up, and",
        "Table 1 reports 1-5 episodes per child in the population PK data set",
        "(101 children with 1 episode, 56 with 2, 34 with 3, 12 with 4 and 4",
        "with 5; median 1.9). Five occasions are therefore encoded, which spans",
        "the observed maximum. Monolix reports a single interoccasion magnitude",
        "per parameter rather than one per occasion, so occasion 1 carries the",
        "estimated variance and occasions 2-5 repeat it via fixed() - the",
        "Monolix analogue of the NONMEM $OMEGA BLOCK(1) SAME idiom used by",
        "Stoschus_2025_phenobarbital.R and Ding_2026_vancomycin.R. Supply",
        "OCC = 1 for a single-episode simulation.",
        sep = " "
      ),
      source_name = "occasion"
    )
  )

  # Covariates that Tchaparian 2016 screened but did not retain in the final
  # model. Documented here rather than in covariateData because the model()
  # block never references them; see Methods, 'Pharmacokinetic Analysis'.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not retained. 47.3 percent of population PK episodes were in male children (Table 1)."
    ),
    WAZ = list(
      description = "Weight-for-age z-score computed against the WHO Child Growth Standards.",
      units = "unitless (z-score; standard-deviation units)",
      type = "continuous",
      notes = paste(
        "Screened; not retained on any PK parameter. Median -0.74 (IQR -3.63 to",
        "2.08) and 26 of 207 episodes (12.6 percent) in underweight children",
        "(Table 1 footnote a). The companion recurrence model",
        "Tchaparian_2016_lumefantrine_recurrence does use WAZ, dichotomised at",
        "-2, as an adjustment covariate.",
        sep = " "
      )
    ),
    HGB = list(
      description = "Blood hemoglobin concentration at diagnosis.",
      units = "g/dL",
      type = "continuous",
      notes = "Screened; not retained. Median 10.0 g/dL (range 5.6-15.9; Table 1)."
    ),
    PARA = list(
      description = "Asexual Plasmodium parasite density at diagnosis, log-transformed for covariate screening.",
      units = "parasites/uL",
      type = "continuous",
      notes = paste(
        "Screened; not retained in the final model. The Supplement notes that",
        "'in earlier stages of the model building, there was a suggestion that",
        "parasite density may have an influence upon either or both volumes of",
        "distribution; this covariate was not significant, however, in the last",
        "step.' Geometric mean 15,568 parasites/uL (Table 1). Note that the",
        "sibling lumefantrine models Kloprogge_2018_lumefantrine.R and",
        "Ding_2026_lumefantrine.R DO retain a parasitaemia effect on F.",
        sep = " "
      )
    ),
    HIV_POS = list(
      description = "HIV-infected indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not retained. Only 9 HIV-infected children contributing 16 episodes (Table 1), which the authors state limited power to detect an effect."
    ),
    CONMED_TMPSMX = list(
      description = "Daily trimethoprim-sulfamethoxazole prophylaxis indicator.",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened; not retained as a population PK covariate, which the authors",
        "flag explicitly in the Discussion ('it should be noted, however, that",
        "TMP-SMZ use was not identified as a significant covariate in our",
        "population pharmacokinetic analysis') even though day 7 concentrations",
        "were higher in children receiving it (median 243 vs 206 ng/mL,",
        "P = .018). 47 children over 80 episodes received prophylaxis (Table 1).",
        "It IS retained in the companion recurrence model",
        "Tchaparian_2016_lumefantrine_recurrence, where it interacts with",
        "lumefantrine exposure.",
        sep = " "
      )
    ),
    BREASTFEED_EXCLUSIVE = list(
      description = "Breast-feeding indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not retained. 95 of 207 population PK episodes involved breast-feeding children (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 101L,
    n_studies = 1L,
    n_episodes = 207L,
    n_observations = 806L,
    age_range = "6.6-22.2 months (median 14.4); enrolment window 6 months to 2 years",
    weight_range = "6.1-13.0 kg (median 9.1)",
    sex_female_pct = 52.7,
    race_ethnicity = "Ugandan children resident in Tororo district; race / ethnicity not otherwise reported",
    disease_state = paste(
      "Uncomplicated Plasmodium falciparum malaria diagnosed on a positive",
      "thick blood smear plus documented or recent fever. A small number of",
      "Plasmodium ovale (5 episodes) and Plasmodium malariae (9 episodes)",
      "infections were retained after species was ruled out as a significant",
      "influence on PK (Supplement). 9 children were HIV-infected and received",
      "nevirapine-, lamivudine- and stavudine- or zidovudine-based",
      "antiretroviral therapy; 47 children received daily",
      "trimethoprim-sulfamethoxazole prophylaxis because they were HIV-infected",
      "or HIV-exposed.",
      sep = " "
    ),
    dose_range = paste(
      "Artemether-lumefantrine (Coartem, Novartis) 20 mg artemether plus 120 mg",
      "lumefantrine (1 tablet) twice daily for 3 days - 6 doses, 720 mg total",
      "lumefantrine - for every child, because all body weights were under",
      "14 kg. Median total body-weight-adjusted lumefantrine dose 81.1 mg/kg",
      "(range 55.4-118.0). Morning doses were given in clinic as crushed",
      "tablets dispersed in water followed by 150 mL of reconstituted cow's",
      "milk containing about 5 g of fat, or by breast-feeding, to standardise",
      "the food effect on lumefantrine absorption; evening doses were given at",
      "home with milk powder. The full 6-dose course was taken in 99.5 percent",
      "of analysed treatments.",
      sep = " "
    ),
    regions = "Tororo district, eastern Uganda (perennial high-transmission area, entomological inoculation rate up to 562 infective bites per person-year)",
    notes = paste(
      "Sampling was sparse and by design: one sample per episode on each of day",
      "0 (pre-dose), day 2 (before the fifth dose), day 3 (12 h after the last",
      "dose, i.e. 72 h), day 7 (108 h after the last dose, i.e. 168 h) and",
      "day 14 (276 h after the last dose, i.e. 336 h). The assay was HPLC on",
      "dried capillary whole blood spotted onto tartaric-acid-pretreated filter",
      "paper, with LLOQ 132 ng/mL and LOD 52 ng/mL. 188 of 806 samples (23.3",
      "percent) fell below the LOD and were retained as left-censored rather",
      "than discarded, which is what lets the model describe the 14-day tail at",
      "all; only 25 of 197 day-14 samples were above the LOD. 56 of 862",
      "candidate concentrations (6.4 percent) were excluded, mostly 44 samples",
      "from 9 episodes with detectable pre-first-dose concentrations.",
      "Because concentrations are capillary WHOLE BLOOD on filter paper rather",
      "than venous plasma, the apparent volumes and the exposure metrics are",
      "not directly comparable with the plasma-based lumefantrine models in",
      "this library (Hoglund_2015_lumefantrine, Kloprogge_2018_lumefantrine,",
      "Ding_2026_lumefantrine); the authors devote a Discussion paragraph to",
      "this. Companion pharmacodynamic model from the same paper:",
      "modellib('Tchaparian_2016_lumefantrine_recurrence').",
      sep = " "
    )
  )

  ini({
    # ================================================================
    # Structural parameters, Tchaparian 2016 Table 2, 'Point Estimate'
    # column. Every disposition term is apparent (divided by F) and is
    # anchored at the final model's reference child: 8.43 kg and 12
    # months, i.e. the typical 1-year-old. Values are on the linear
    # scale in the source and log() is applied here for the nlmixr2
    # internal log scale.
    # ================================================================

    lcl <- log(2.19)
    label("Apparent lumefantrine clearance CL/F at 8.43 kg (L/h)")
    # Table 2: CL/F = 2.19 L/h, RSE 8 percent. Results text confirms
    # 'the prediction of CL/F in a typical 12-month-old child is 2.19 L/h'.

    lvc <- log(83.2)
    label("Apparent lumefantrine central volume of distribution V1/F at 8.43 kg (L)")
    # Table 2: V1/F = 83.2 L/8.43 kg, RSE 7 percent

    lq <- log(0.23)
    label("Apparent lumefantrine intercompartmental clearance Q/F at 8.43 kg (L/h)")
    # Table 2: Q/F = 0.23 L/h, RSE 35 percent

    lvp <- log(441)
    label("Apparent lumefantrine peripheral volume of distribution V2/F at 8.43 kg (L)")
    # Table 2: V2/F = 441 L/8.43 kg, RSE 74 percent

    lka <- fixed(log(0.45))
    label("First-order lumefantrine absorption rate constant ka (1/h)")
    # Table 2 footnote c marks ka as not estimated. Methods: 'owing to
    # limited data in the absorption phase, the rate constant of
    # absorption (ka) was fixed to the previously reported value of
    # 0.45 h-1' (the paper's reference 24).

    lfdepot <- fixed(log(1))
    label("Relative bioavailability F at the reference age of 12 months (unitless)")
    # Table 2 footnote c: F = 1, not estimated. This is the standard
    # anchor of an apparent-parameter (CL/F, V/F) parameterisation; the
    # age covariate below is a deviation from it, and the variability
    # on F is what induces correlation among the four apparent
    # disposition terms (Supplement, 'Population Pharmacokinetic
    # analysis').

    # ---------------- Covariate effects ----------------

    e_wt_cl_q <- fixed(0.75)
    label("Allometric exponent on body weight shared by CL/F and Q/F (unitless)")
    # Table 2 footnote a: 'Model: Pk * [WT/8.43]^0.75, in which WT is
    # total body weight, and 8.43 is the typical weight for a
    # 12-month-old child'. Methods states the exponent was imposed
    # unconditionally, not estimated, and no RSE is reported for it.

    e_wt_vc_vp <- fixed(1)
    label("Allometric exponent on body weight shared by V1/F and V2/F (unitless)")
    # Methods, 'Pharmacokinetic Analysis': allometric scaling was
    # applied 'with respect to the 2 clearance parameters by
    # multiplying each parameter by [weight/reference weight]^0.75 and,
    # for the 2 volume terms, by multiplying each parameter by
    # [weight/reference weight]'. Imposed, not estimated.

    e_age_fdepot <- 0.596
    label("Power exponent on (age in months / 12) for relative bioavailability F (unitless)")
    # Table 2: 'Age effect' = 0.596, RSE 38 percent, with footnote b
    # giving the covariate model as F * [age/12]^theta_x. Results
    # writes the equation out: F_i = [1 * Age_i/12]^0.596.

    # ================================================================
    # Interindividual variability, Table 2 'IIV, percent (RSE, percent)'
    # column. The Supplement fixes the scale of this column beyond
    # doubt: 'IIV and IOV of Q/F and V2/F were each fixed to be 50
    # percent (CV)', and the Results text repeats it for clearance
    # ('2.19 L/h, with 8 percent (CV) interindividual variability').
    # The entries are therefore coefficients of variation of a
    # log-normal distribution, not variances, so each is converted to
    # the internal log-scale variance by omega^2 = log(CV^2 + 1).
    # ================================================================

    etalcl ~ 0.0063796
    # Table 2, CL/F row, IIV = 8 CV (RSE 174); log(0.08^2 + 1) = 0.0063796

    etalvc ~ 0.0228411
    # Table 2, V1/F row, IIV = 15.2 CV (RSE 46); log(0.152^2 + 1) = 0.0228411

    etalq ~ fixed(0.2231436)
    # Table 2, Q/F row, IIV = 50 CV, footnote c; constrained by the
    # authors because the data were too sparse to estimate it
    # (Supplement); log(0.5^2 + 1) = 0.2231436

    etalvp ~ fixed(0.2231436)
    # Table 2, V2/F row, IIV = 50 CV, footnote c; constrained by the
    # authors on the same grounds as Q/F; log(0.5^2 + 1) = 0.2231436

    etalfdepot ~ 0.1149354
    # Table 2, F row, IIV = 34.9 CV (RSE 23); log(0.349^2 + 1) = 0.1149354

    # ================================================================
    # Interoccasion variability, Table 2 'IOV (percent)' column, same
    # CV-to-variance conversion. The occasion is a malaria episode.
    # Monolix estimates ONE magnitude per parameter that is shared by
    # every occasion, so occasion 1 carries the value and occasions 2-5
    # repeat it through fixed() - the same encoding the registered
    # Monolix precedent Stoschus_2025_phenobarbital.R uses for the
    # NONMEM $OMEGA BLOCK(1) SAME idiom.
    # ================================================================

    etaiov_cl_1 ~ 0.0086118
    # Table 2, CL/F row, IOV = 9.3 CV (RSE 105); log(0.093^2 + 1) = 0.0086118
    etaiov_cl_2 ~ fixed(0.0086118) # shared magnitude, occasion 2
    etaiov_cl_3 ~ fixed(0.0086118) # shared magnitude, occasion 3
    etaiov_cl_4 ~ fixed(0.0086118) # shared magnitude, occasion 4
    etaiov_cl_5 ~ fixed(0.0086118) # shared magnitude, occasion 5

    etaiov_vc_1 ~ 0.0073688
    # Table 2, V1/F row, IOV = 8.6 CV (RSE 133); log(0.086^2 + 1) = 0.0073688
    etaiov_vc_2 ~ fixed(0.0073688) # shared magnitude, occasion 2
    etaiov_vc_3 ~ fixed(0.0073688) # shared magnitude, occasion 3
    etaiov_vc_4 ~ fixed(0.0073688) # shared magnitude, occasion 4
    etaiov_vc_5 ~ fixed(0.0073688) # shared magnitude, occasion 5

    etaiov_q_1 ~ fixed(0.2231436)
    # Table 2, Q/F row, IOV = 50 CV, footnote c; constrained by the
    # authors (Supplement); log(0.5^2 + 1) = 0.2231436
    etaiov_q_2 ~ fixed(0.2231436) # shared magnitude, occasion 2
    etaiov_q_3 ~ fixed(0.2231436) # shared magnitude, occasion 3
    etaiov_q_4 ~ fixed(0.2231436) # shared magnitude, occasion 4
    etaiov_q_5 ~ fixed(0.2231436) # shared magnitude, occasion 5

    etaiov_vp_1 ~ fixed(0.2231436)
    # Table 2, V2/F row, IOV = 50 CV, footnote c; constrained by the
    # authors (Supplement); log(0.5^2 + 1) = 0.2231436
    etaiov_vp_2 ~ fixed(0.2231436) # shared magnitude, occasion 2
    etaiov_vp_3 ~ fixed(0.2231436) # shared magnitude, occasion 3
    etaiov_vp_4 ~ fixed(0.2231436) # shared magnitude, occasion 4
    etaiov_vp_5 ~ fixed(0.2231436) # shared magnitude, occasion 5

    etaiov_fdepot_1 ~ 0.3022031
    # Table 2, F row, IOV = 59.4 CV (RSE 8); log(0.594^2 + 1) = 0.3022031
    etaiov_fdepot_2 ~ fixed(0.3022031) # shared magnitude, occasion 2
    etaiov_fdepot_3 ~ fixed(0.3022031) # shared magnitude, occasion 3
    etaiov_fdepot_4 ~ fixed(0.3022031) # shared magnitude, occasion 4
    etaiov_fdepot_5 ~ fixed(0.3022031) # shared magnitude, occasion 5

    # ================================================================
    # Residual error. Combined additive plus proportional (Supplement:
    # 'a combined additive plus proportional error model was used for
    # residual variability'). Results states the point estimates are
    # standard deviations, so they map directly onto propSd / addSd
    # without a variance-to-SD conversion.
    # ================================================================

    propSd <- 0.378
    label("Proportional residual SD for capillary whole-blood lumefantrine concentration (SD on the linear scale)")
    # Table 2: Proportional = 37.8 percent, RSE 5 percent. Results:
    # 'the point estimate of residual error (standard deviation) was
    # 37.8 percent (proportional) and 19.2 ng/mL (additive)'.

    addSd <- 19.2
    label("Additive residual SD for capillary whole-blood lumefantrine concentration (ng/mL)")
    # Table 2: Additive = 19.2 ng/mL, RSE 18 percent. For orientation
    # the assay LOD is 52 ng/mL and the LLOQ 132 ng/mL.
  })

  model({
    # ---- Reference values of the final model (Table 2 footnote a and
    # the Results covariate equation). The authors moved from the
    # cohort median 9.0 kg to the 1-year-old's 8.43 kg when they fitted
    # the final model, so these two constants and the Table 2 point
    # estimates belong together.
    WT_REF <- 8.43 # kg, typical weight of a 12-month-old child
    AGE_REF_MO <- 12 # months, reference age

    # ---- Occasion indicators. OCC is the malaria-episode index; the
    # five indicators are mutually exclusive and any OCC outside 1-5
    # zeroes all of them, leaving that record with interindividual
    # variability only.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)

    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5
    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2 + oc3 * etaiov_vc_3 +
      oc4 * etaiov_vc_4 + oc5 * etaiov_vc_5
    iov_q <- oc1 * etaiov_q_1 + oc2 * etaiov_q_2 + oc3 * etaiov_q_3 +
      oc4 * etaiov_q_4 + oc5 * etaiov_q_5
    iov_vp <- oc1 * etaiov_vp_1 + oc2 * etaiov_vp_2 + oc3 * etaiov_vp_3 +
      oc4 * etaiov_vp_4 + oc5 * etaiov_vp_5
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 +
      oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4 + oc5 * etaiov_fdepot_5

    # ---- Individual apparent disposition parameters. Allometric
    # scaling is on total body weight referenced to 8.43 kg, with the
    # 0.75 exponent shared by the two clearances and the exponent of 1
    # shared by the two volumes (Methods, 'Pharmacokinetic Analysis').
    cl <- exp(lcl + etalcl + iov_cl) * (WT / WT_REF)^e_wt_cl_q
    vc <- exp(lvc + etalvc + iov_vc) * (WT / WT_REF)^e_wt_vc_vp
    q <- exp(lq + etalq + iov_q) * (WT / WT_REF)^e_wt_cl_q
    vp <- exp(lvp + etalvp + iov_vp) * (WT / WT_REF)^e_wt_vc_vp
    ka <- exp(lka)

    # ---- Relative bioavailability. Results gives the covariate model
    # verbatim as F_i = [1 * Age_i/12]^0.596 with age in MONTHS, so the
    # canonical AGE column (years) is converted here. A 6-month-old has
    # F = 0.66 and a 2-year-old F = 1.51, reproducing the 0.7 and 1.5
    # the Discussion quotes.
    AGE_MO <- AGE * 12
    fdepot <- exp(lfdepot + etalfdepot + iov_fdepot) *
      (AGE_MO / AGE_REF_MO)^e_age_fdepot

    # ---- Disposition micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- Two-compartment disposition with first-order absorption.
    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    # ---- Capillary whole-blood concentration in ng/mL. The dose is in
    # mg and vc in L, so central / vc is mg/L = ug/mL; the factor of
    # 1000 converts to the ng/mL scale that Table 2 and Figures 1-2
    # report (same convention as Kay_2020_lumefantrine.R).
    Cc <- 1000 * central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
