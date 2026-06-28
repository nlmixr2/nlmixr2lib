vanIersel_2018_posaconazole <- function() {
  description <- "Population PK model for the delayed-release solid oral tablet formulation of posaconazole in adult healthy volunteers and patients at high risk for invasive fungal disease (van Iersel 2018). One-compartment disposition with sequential zero-order then first-order absorption: each oral dose loads into the depot compartment as a zero-order infusion of duration D1, after which depot drains to central with first-order rate constant ka and central eliminates with first-order rate constant CL/V. The random effect on D1 is the same as the random effect on ka multiplied by a correlation factor (cor_kad1 = -0.586). Covariates retained in the final model are body weight on relative bioavailability (allometric power exponent), tablet formulation A/B versus C/D on bioavailability, AML/MDS disease state on bioavailability, fed status on absorption rate, and single-dose-versus-multiple-dose record indicator on clearance. Residual variability is log-additive with separate magnitudes for phase 1 versus phase 3 studies."
  reference   <- "van Iersel MLPS, Rossenu S, de Greef R, Waskin H. A population pharmacokinetic model for a solid oral tablet formulation of posaconazole. Antimicrob Agents Chemother. 2018;62(7):e02465-17. doi:10.1128/AAC.02465-17."
  vignette    <- "vanIersel_2018_posaconazole"
  units       <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject in the source analysis. Allometric-style power effect on relative bioavailability F1 with reference 74.9 kg (the simulation-cohort median used in van Iersel 2018 Fig S4 and matching the overall ~75 kg mean across the six contributing studies, Table 1).",
      source_name        = "WT"
    ),
    FED = list(
      description        = "Fed-state indicator at the dosing event",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted at dosing; the typical-value ka reference in van Iersel 2018 Table 2)",
      notes              = "Per-dose-record indicator switching the absorption rate constant ka via a multiplicative (1 + 0.530 * FED) effect; 53% higher ka after a meal (van Iersel 2018 Table 2 row 'Food intake on ka').",
      source_name        = "Food status"
    ),
    MULTI_DOSE_PT = list(
      description        = "Multiple-dose-record indicator (1 = dose record from the multiple-dose phase of a study; 0 = single-dose record)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (single-dose record)",
      notes              = "Per-dose-record indicator switching apparent clearance via a multiplicative (1 + 0.750 * MULTI_DOSE_PT) effect; 75% higher CL in multiple-dose records (van Iersel 2018 Table 2 row 'Dosing regimen on CL'). The dose-record-level encoding follows the MULTI_DOSE_PT canonical (Goel 2016 Sonidegib precedent).",
      source_name        = "Dosing regimen (single dose vs multiple doses)"
    ),
    FORM_POSA_AB = list(
      description        = "Posaconazole prototype-tablet formulation indicator (1 = tablet formulation A or B; 0 = tablet formulation C or D, the later / marketed formulations including the commercial formulation D)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (tablet C or D; van Iersel 2018 Results 'Posaconazole exposure in the clinical setting and impact of covariates on exposure' states tablet formulation D is the marketed image and the predominant formulation used across the data set and was taken as the reference formulation for the relative bioavailability estimate)",
      notes              = "Per-dose-record indicator switching relative bioavailability F1 via a multiplicative (1 + 0.247 * FORM_POSA_AB) effect; 24.7% higher F1 for the prototype A/B formulations (van Iersel 2018 Table 2 row 'Tablet formulation A/B on F1'). Drug-specific member of the FORM_<drug>_<contrast> family (mirrors FORM_ASV_LIQUID for asunaprevir, FORM_SAR_DP2 for sarilumab).",
      source_name        = "Tablet formulation (A/B vs C/D)"
    ),
    DIS_MDS_AML = list(
      description        = "Combined AML/MDS disease-state indicator (1 = patient with acute myeloid leukemia or myelodysplastic syndrome; 0 = healthy volunteer or HSCT recipient)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-AML/MDS subjects: healthy volunteers and HSCT recipients in the van Iersel 2018 pooled analysis)",
      notes              = "Time-fixed per subject. Switches relative bioavailability F1 via a multiplicative (1 + (-0.165) * DIS_MDS_AML) effect; 16.5% lower F1 in AML/MDS patients versus the non-AML/MDS reference (van Iersel 2018 Table 2 row 'AML/MDS on F1'). The combined MDS+AML encoding follows the DIS_MDS_AML canonical (Ogasawara 2020 durvalumab precedent).",
      source_name        = "Disease state (AML/MDS)"
    ),
    STUDY_POSA_PHASE3 = list(
      description        = "Posaconazole phase 3 study indicator (1 = subject enrolled in the phase 3 patient study P05615 of the van Iersel 2018 pooled analysis; 0 = phase 1 studies P04975, P05637, P07764, P07783, P07691 of healthy volunteers)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (phase 1 studies in healthy volunteers)",
      notes              = "Per-subject (study-fixed) indicator switching the log-additive residual-error magnitude between expSd_p1 (0.42, phase 1 studies) and expSd_p3 (0.322, phase 3 study); van Iersel 2018 Results 'Model development' and Discussion 'The residual error model with two random-effect parameters'. Drug-specific paper-anchored canonical following the STUDY_<DRUG>_PHASE<N> family (mirrors STUDY_NIPOCALIMAB_PHASE1 for nipocalimab, STUDY_FARLETUZUMAB_PHASE2 for farletuzumab).",
      source_name        = "Study phase (phase 1 vs phase 3)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 335L,
    n_studies      = 6L,
    age_range      = "31.4-51.0 years (study means; van Iersel 2018 Table 1)",
    weight_range   = "73.9-79.6 kg (study means; van Iersel 2018 Table 1)",
    weight_median  = "74.9 kg (simulation-cohort median used as the body-weight reference in Fig S4 and in the model file's allometric F1 normalisation)",
    sex_female_pct = "14-50% across studies (van Iersel 2018 Table 1; 38% in the phase 3 study P05615)",
    race_ethnicity = "Predominantly white (26-100% across studies; van Iersel 2018 Table 1; van Iersel 2018 Discussion notes the limited race distribution as a potential limitation)",
    disease_state  = "Pooled cohort of 104 healthy volunteers (studies P04975, P05637, P07764, P07783, P07691) and 231 patients (study P05615) at high risk for invasive fungal disease: 125 AML, 82 MDS, 6 GVHD (HSCT recipients), 18 missing.",
    dose_range     = "100-400 mg oral posaconazole as the delayed-release solid tablet (formulations A, B, C, or D; D = commercial formulation), single or multiple doses across the six studies.",
    regions        = "Not reported in detail",
    n_observations = "5,756 plasma posaconazole concentration measurements (van Iersel 2018 Results 'Model development').",
    notes          = "Phase 1 studies P04975/P05637/P07764/P07783/P07691 contributed 848 + 1,253 + 335 + 505 + 675 = 3,616 samples in healthy volunteers; the phase 3 patient study P05615 contributed 2,140 samples. Five samples with CWRES > 6 were treated as outliers during base-model development but reintroduced in the final model used here (Table 2 final-model column with outliers). PK sampling schemes per study are summarised in van Iersel 2018 Table 3."
  )

  ini({
    # ====================================================================
    # Structural PK parameters (van Iersel 2018 Table 2 final-model
    # column 'Estimate (% RSE) for final model with outliers') --
    # reference covariates: WT = 74.9 kg, fasted (FED = 0), single-dose
    # record (MULTI_DOSE_PT = 0), tablet C or D (FORM_POSA_AB = 0),
    # non-AML/MDS subject (DIS_MDS_AML = 0). One-compartment disposition
    # with sequential zero-order then first-order absorption: each oral
    # dose loads into depot as a zero-order infusion of duration D1, then
    # depot drains to central with first-order rate constant ka, then
    # central eliminates with first-order rate constant CL/V.
    # ====================================================================
    lcl <- log(9.70);  label("Apparent clearance CL/F at reference covariates (L/h)")  # van Iersel 2018 Table 2: CL/F = 9.70 L/h, RSE 5.00%
    lvc <- log(393);   label("Apparent central volume of distribution V/F (L)")         # van Iersel 2018 Table 2: V    = 393 L,  RSE 2.77%
    lka <- log(0.853); label("First-order absorption rate constant ka at reference covariates (1/h)")  # van Iersel 2018 Table 2: ka  = 0.853 1/h, RSE 7.75%
    ld1 <- log(2.54);  label("Duration of zero-order absorption into the depot compartment D1 (h)")     # van Iersel 2018 Table 2: D1  = 2.54 h,   RSE 3.45%

    # Bioavailability anchor. F1 = 1 at the reference covariate set
    # (tablet D in healthy / HSCT recipients at WT = 74.9 kg), per
    # van Iersel 2018 Results 'Posaconazole exposure ... and impact of
    # covariates': 'tablet formulation D is the marketed image and the
    # predominant formulation used across the data set and was taken as
    # the reference formulation for the relative bioavailability estimate'.
    lfdepot <- fixed(log(1)); label("Reference relative bioavailability F1 (unitless; tablet D in non-AML/MDS subjects at WT = 74.9 kg, fixed anchor)")

    # Correlation factor coupling the random effect on D1 to etalka.
    # van Iersel 2018 Methods 'Covariate model': 'To account for the high
    # degree of correlation between ka and D1, the model was modified so
    # that the random effect for D1 was the same as that for the ka
    # parameter multiplied by a correlation factor.' Encoded as
    #     d1_i = exp(ld1 + cor_kad1 * etalka)
    # so D1 borrows its individual-level deviation directly from etalka
    # and there is no separate etald1 estimated in the final model.
    cor_kad1 <- -0.586; label("Correlation factor linking etalka to the individual-level deviation of D1 (unitless)")  # van Iersel 2018 Table 2: COR = -0.586, RSE 15.2%

    # Covariate effects (van Iersel 2018 Table 2 final-model column).
    # Body weight on F1 is encoded as an allometric-style power
    # exponent on (WT / 74.9); the remaining four covariate effects are
    # encoded as multiplicative (1 + theta * cov) shifts on their
    # respective typical-value parameters. The NONMEM .lst file was not
    # available; the (1 + theta * cov) form is the convention used here
    # and is documented in the vignette Assumptions and deviations.
    e_wt_fdepot          <- -1.03;   label("Allometric power exponent of body weight on relative bioavailability F1 (unitless; F1 ~ (WT/74.9)^e_wt_fdepot)")  # van Iersel 2018 Table 2: 'Body wt on F1' = -1.03, RSE 8.91%
    e_md_cl              <-  0.750;  label("Multiplicative effect of MULTI_DOSE_PT on apparent clearance CL/F (unitless; CL ~ (1 + e_md_cl * MULTI_DOSE_PT))")  # van Iersel 2018 Table 2: 'Dosing regimen on CL' = 0.750, RSE 5.87%
    e_formposaab_fdepot  <-  0.247;  label("Multiplicative effect of FORM_POSA_AB on relative bioavailability F1 (unitless; F1 ~ (1 + e_formposaab_fdepot * FORM_POSA_AB))")  # van Iersel 2018 Table 2: 'Tablet formulation A/B on F1' = 0.247, RSE 21.5%
    e_fed_ka             <-  0.530;  label("Multiplicative effect of FED on first-order absorption rate ka (unitless; ka ~ (1 + e_fed_ka * FED))")  # van Iersel 2018 Table 2: 'Food intake on ka' = 0.530, RSE 17.3%
    e_dis_mds_aml_fdepot <- -0.165;  label("Multiplicative effect of DIS_MDS_AML on relative bioavailability F1 (unitless; F1 ~ (1 + e_dis_mds_aml_fdepot * DIS_MDS_AML))")  # van Iersel 2018 Table 2: 'AML/MDS on F1' = -0.165, RSE 25.8%

    # Inter-individual variability. van Iersel 2018 Table 2 reports IIV
    # as CV%; convert to log-scale variance via omega^2 = log(1 + CV^2)
    # per the lognormal IIV convention. The IIV magnitudes here are the
    # final-model (with outliers) values; IOV on ka, D1, and F1 (also
    # reported in Table 2) is NOT encoded because the source analysis
    # does not define an operational occasion column for the simulation
    # use case (Brooks 2021 / Andrews 2017 nlmixr2lib precedent); the
    # IOV magnitudes are recorded in the vignette Assumptions and
    # deviations section.
    etalcl     ~ log(1 + 0.379^2)  # van Iersel 2018 Table 2: IIV for CL = 37.9% CV, RSE 13.1%; omega^2 = log(1 + 0.379^2) = 0.1343
    etalka     ~ log(1 + 0.575^2)  # van Iersel 2018 Table 2: IIV for ka = 57.5% CV, RSE 29.3%; omega^2 = log(1 + 0.575^2) = 0.2858
    etalfdepot ~ log(1 + 0.242^2)  # van Iersel 2018 Table 2: IIV for F1 = 24.2% CV, RSE 26.7%; omega^2 = log(1 + 0.242^2) = 0.0569

    # Residual error -- log-additive on Cc (NONMEM 'additive on log-
    # transformed observation' per van Iersel 2018 Methods 'Pharmaco-
    # kinetic model development': 'log-transformed, both-sides
    # approach'). Encoded as ~ lnorm(expSd) so the SD applies directly
    # on log(Cc). Two study-phase-specific SDs are switched at
    # observation time via the STUDY_POSA_PHASE3 indicator.
    expSd_p1 <- 0.42;  label("Log-scale residual SD on Cc for phase 1 studies (van Iersel 2018 NONMEM log-additive)")  # van Iersel 2018 Table 2: 'SD (phase 1 studies)' = 0.42, RSE 8.69%
    expSd_p3 <- 0.322; label("Log-scale residual SD on Cc for phase 3 study (van Iersel 2018 NONMEM log-additive)")    # van Iersel 2018 Table 2: 'SD (phase 3 study)'   = 0.322, RSE 10.3%
  })

  model({
    # ====================================================================
    # 1. Derived covariate terms
    # ====================================================================
    # Body weight is normalised to the simulation-cohort median 74.9 kg
    # (van Iersel 2018 Fig S4 'median 74.9 kg'; matches the ~75 kg mean
    # of the six studies' Table 1 weights).
    wt_ratio <- WT / 74.9

    # ====================================================================
    # 2. Individual PK parameters
    # ====================================================================
    cl <- exp(lcl + etalcl) * (1 + e_md_cl * MULTI_DOSE_PT)
    vc <- exp(lvc)
    ka <- exp(lka + etalka) * (1 + e_fed_ka * FED)

    # D1 borrows its individual-level deviation from etalka via
    # cor_kad1 (van Iersel 2018 Methods 'Covariate model').
    d1 <- exp(ld1 + cor_kad1 * etalka)

    # Relative bioavailability F1 with the four retained covariate
    # effects (allometric on WT, additive shifts on FORM_POSA_AB and
    # DIS_MDS_AML).
    fdepot <- exp(lfdepot + etalfdepot) *
              wt_ratio^e_wt_fdepot *
              (1 + e_formposaab_fdepot * FORM_POSA_AB) *
              (1 + e_dis_mds_aml_fdepot * DIS_MDS_AML)

    # ====================================================================
    # 3. Micro-constants
    # ====================================================================
    kel <- cl / vc

    # ====================================================================
    # 4. ODE system -- 1-compartment with sequential zero-then-first-
    # order absorption. Each oral dose record in the event table targets
    # cmt = depot with amt = nominal dose; dur(depot) <- d1 imposes the
    # zero-order infusion duration into depot; depot then drains to
    # central with first-order rate ka; central eliminates with first-
    # order rate kel.
    # ====================================================================
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ====================================================================
    # 5. Bioavailability and zero-order duration
    # ====================================================================
    f(depot)   <- fdepot
    dur(depot) <- d1

    # ====================================================================
    # 6. Observation and residual error
    # ====================================================================
    # Cc (ng/mL) = central (mg) / vc (L) * 1e6 (mg to ng) / 1e3 (L to mL)
    #            = central / vc * 1000.
    Cc <- (central / vc) * 1000

    # Study-phase-specific log-additive residual SD: phase 1 SD applies
    # when STUDY_POSA_PHASE3 = 0; phase 3 SD applies when = 1.
    expSd <- expSd_p3 * STUDY_POSA_PHASE3 + expSd_p1 * (1 - STUDY_POSA_PHASE3)
    Cc ~ lnorm(expSd)
  })
}
