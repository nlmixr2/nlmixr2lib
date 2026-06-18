Svensson_2012_nevirapine <- function() {
  description <- paste(
    "One-compartment population PK model for oral nevirapine in HIV-infected",
    "South African adults (Svensson 2012, the multi-source 'mega-model'",
    "integration paper) with first-order absorption through two transit",
    "compartments (Savic 2007 parameterisation, ktr = (NTRANS+1)/MTT shared",
    "across all transit-rate steps), a two-population mixture on apparent oral",
    "clearance CL/F (fast eliminators 3.12 L/h at 82.7% probability vs slow",
    "eliminators 1.45 L/h at 17.3% probability, the slow class associated by",
    "the paper Discussion with CYP2B6 516TT homozygotes), Anderson-Holford",
    "allometric scaling (fat-free mass at exponent 0.75 for CL/F with reference",
    "FFM 42 kg corresponding to a 70 kg / 1.6 m woman, body weight at exponent",
    "1 for V/F with reference WT 70 kg, both exponents fixed), a fed/fasted",
    "binary covariate on absorption mean transit time MTT (2.46 h fed vs 0.596",
    "h fasted, a 4.1-fold slowing of absorption with food) and a concomitant",
    "tuberculosis-treatment (rifampicin + isoniazid +/- ethambutol) effect on",
    "bioavailability F (39% decrease; F = 0.613 when on TB treatment vs",
    "F = 1.0 reference, with an additional 34.1% between-subject variability",
    "specific to the TB-treatment effect). Residual error is proportional",
    "(8.41% CV)."
  )
  reference <- paste(
    "Svensson E, van der Walt JS, Barnes KI, Cohen K, Kredo T, Huitema A,",
    "Nachega JB, Karlsson MO, Denti P (2012).",
    "Integration of data from multiple sources for simultaneous modelling",
    "analysis: experience from nevirapine population pharmacokinetics.",
    "Br J Clin Pharmacol 74(3):465-476.",
    "doi:10.1111/j.1365-2125.2012.04205.x",
    sep = " "
  )
  vignette <- "Svensson_2012_nevirapine"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline; cohort median 67-72 kg across the three",
        "South African studies (range 43-128 kg per Table 1). Used for the",
        "Anderson-Holford allometric scaling of apparent volume of",
        "distribution V/F with exponent 1.0 (fixed) and reference 70 kg per",
        "Svensson 2012 Methods 'Modelling' (Equation 2 cites the customary",
        "0.75/1.0 fixed exponents) and Results 'Population pharmacokinetics'",
        "paragraph 2 ('For BW a reference value 70 kg was used, while for FFM",
        "the reference was 42 kg, corresponding to a woman of 70 kg and 1.6",
        "m.'). Note that FFM (not WT) is used for CL/F allometric scaling per",
        "the same paragraph."
      ),
      source_name        = "BW"
    ),
    FFM = list(
      description        = "Fat-free mass derived from body weight, height, and sex via the Janmahasatian formula",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used for the Anderson-Holford allometric scaling of apparent oral",
        "clearance CL/F with exponent 0.75 (fixed) and reference 42 kg per",
        "Svensson 2012 Results 'Population pharmacokinetics' paragraph 2",
        "('for FFM the reference was 42 kg, corresponding to a woman of 70",
        "kg and 1.6 m'). FFM was the best size descriptor for CL/F over body",
        "weight or normal fat weight (NFW; paragraph 2 same section). For",
        "subjects with missing height (studies 2 and 4 in the source paper),",
        "the authors developed an exponential imputation model (Equation 5):",
        "FFM = FFMmax * (1 - exp(-K * WT)) with parameters FFMmax = 104 kg",
        "/ K = 0.0107 for South African men, FFMmax = 69.6 kg / K = 0.0131",
        "for South African women, FFMmax = 106 kg / K = 0.0106 for P3M-",
        "matched men, FFMmax = 73.2 kg / K = 0.0126 for P3M-matched women.",
        "Users with directly measured WT, height, and sex should compute FFM",
        "via the Janmahasatian (2005) formula and supply it as a covariate;",
        "users without height can use the paper's imputation model."
      ),
      source_name        = "FFM"
    ),
    FED = list(
      description        = "Fed-vs-fasted dose-record indicator (per-occasion)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted dose; e.g. morning dose in the rich-sampling protocol of study 1)",
      notes              = paste(
        "Per-occasion indicator: 1 = nevirapine dose taken under fed",
        "conditions (e.g. evening dose with a meal in study 1's rich-",
        "sampling protocol), 0 = nevirapine dose taken under fasted",
        "conditions (e.g. morning dose). Drives a 4.1-fold slowing of the",
        "absorption mean transit time MTT: 2.46 h fed vs 0.596 h fasted",
        "per Svensson 2012 Table 2 'MTT (h) fed = 2.46' and 'MTT (h) fasted",
        "= 0.596' and Results 'Population pharmacokinetics' paragraph 4",
        "(food on MTT, dropOFV = -90.4 points). The paper found no effect",
        "of food on extent of absorption F (only rate), and the combined",
        "mega-model design separated the food effect from a diurnal effect",
        "by pooling studies with different day/night sampling protocols.",
        "Users simulating a typical fed dose should set FED = 1; a typical",
        "fasted dose, FED = 0."
      ),
      source_name        = "FED (per the paper's narrative; the underlying NONMEM dataset is not on disk)"
    ),
    TB_POS = list(
      description        = "Active tuberculosis with concomitant TB-treatment regimen indicator (rifampicin + isoniazid +/- ethambutol)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no active TB; no concomitant TB treatment)",
      notes              = paste(
        "Time-fixed per subject during the TB-treatment window; effectively",
        "indicates patients who are on the continuation phase of TB therapy",
        "(rifampicin + isoniazid, with ethambutol in a small minority).",
        "TB-coinfected patients in study 1 contributed paired samples: one",
        "while on TB treatment (TB_POS = 1) and one 3 months after",
        "concluding TB treatment (TB_POS = 0); studies 2 and 3 did not",
        "include patients on TB treatment. Drives a 39% decrease in",
        "nevirapine bioavailability F: F = 0.613 when TB_POS = 1 vs F = 1.0",
        "reference per Svensson 2012 Table 2 'F (%) when TB treatment =",
        "61.30' and Results 'Population pharmacokinetics' paragraph 4 (TB",
        "effect on F, dropOFV = -84.6 points). The TB-treatment effect on F",
        "also carries a 34.1% between-subject variability per Table 2",
        "'BSV F when TB treatment (%) = 34.10', encoded here as a separate",
        "log-normal eta gated on TB_POS. The paper attributes the F effect",
        "biologically to rifampicin-mediated CYP3A4 induction in the gut",
        "wall and to first-pass metabolism (Discussion paragraph 5)."
      ),
      source_name        = "TB (paper narrative description of TB-treatment status)"
    ),
    MIX_SLOW_ELIM_NVP = list(
      description        = "Per-subject latent mixture-model class indicator for the slow-vs-fast CL/F sub-population in the Svensson 2012 nevirapine model",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fast-eliminator subpopulation; majority class at 82.7% of the source cohort, typical CL/F = 3.12 L/h)",
      notes              = paste(
        "1 = subject classified to the slow-eliminator subpopulation",
        "(minority class, 17.3% of the source cohort, typical CL/F = 1.45",
        "L/h); 0 = subject classified to the fast-eliminator subpopulation",
        "(majority class, 82.7%, typical CL/F = 3.12 L/h). Not a measured",
        "clinical covariate -- the mixture assignment is a posterior class",
        "assignment from the NONMEM mixture block (Svensson 2012 Results",
        "'Population pharmacokinetics' paragraph 1: 'A mixture model",
        "estimating two typical values of CL/F and the probability of",
        "belonging to the low-clearance population improved the predictive",
        "performance of the model, decreasing BSV in CL from 38 to 25%,",
        "and prevented these outlying data from unduly biasing the",
        "estimation of the typical value of CL/F of the majority of the",
        "population.'). The estimated population probability of belonging",
        "to the slow class is 17.30% (95% CI from the 45.30% RSE: very",
        "imprecise; Table 2 'Probability (%) of belonging to pop. 2 =",
        "17.30 (45.30)'). The paper Discussion paragraph 4 attributes the",
        "slow class biologically to CYP2B6 516TT homozygotes ('The 516G > T",
        "mutation in CYP2B6 is associated with substantial loss of",
        "metabolic function, and the frequency of 516TT homozygotes in a",
        "South African population has been reported to be 13-23%. The",
        "estimate of the probability of belonging to the low-clearance",
        "population (17.3%) agrees well with the reported proportion of",
        "CYP2B6-516-TT homozygotes in the South African population,",
        "implying biological plausibility of two populations with",
        "different clearance rates.'). For typical-value simulation set",
        "MIX_SLOW_ELIM_NVP = 0 to reproduce the majority fast-eliminator",
        "phenotype; set MIX_SLOW_ELIM_NVP = 1 to reproduce the rarer slow",
        "phenotype. For population simulation, draw MIX_SLOW_ELIM_NVP ~",
        "Bernoulli(0.173) per subject. Drug-specific name (_NVP suffix)",
        "follows the precedent of CYP2B6-driven slow-eliminator mixtures",
        "in nevirapine PK; future fast/slow eliminator mixture extractions",
        "for other drugs should register their own canonical or extend a",
        "more general MIX_SLOW_ELIM."
      ),
      source_name        = "NONMEM $MIXTURE class assignment (component 1 = fast, component 2 = slow); MIX_SLOW_ELIM_NVP = as.integer(MIXTURE == 2)"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 115L,
    n_studies        = 3L,
    age_range        = "21-60 years (across the three South African studies; Study 1 median 34, Study 2 median 32, Study 3 median 32)",
    weight_range     = "43-128 kg (across the three South African studies; Study 1 median 67, Study 2 median 72, Study 3 median 60)",
    sex_female_pct   = 83.5,
    race_ethnicity   = "African (Cape Town, South Africa)",
    disease_state    = paste(
      "HIV-1 infection on antiretroviral therapy regimens including 200 mg",
      "nevirapine twice daily at steady state. Patients in Study 1 were",
      "stratified by tuberculosis co-infection: TB-coinfected patients were",
      "sampled both during the continuation phase of anti-TB treatment",
      "(rifampicin + isoniazid, with ethambutol in a small subset) and 3",
      "months after concluding TB treatment. Studies 2 and 3 had no TB-",
      "coinfected patients; Study 3 was an interaction study with the",
      "antimalarial artemether-lumefantrine (no malaria co-infection)."
    ),
    dose_range       = paste(
      "200 mg oral nevirapine twice daily at steady state across all three",
      "model-development studies. Concomitant antiretroviral therapy",
      "(other NRTI/NNRTI backbone agents) varied between and within",
      "studies and was not retained as a PK covariate in the final model."
    ),
    regions          = "South Africa (Cape Town: Study 1, Study 2, Study 3)",
    notes            = paste(
      "Pooled mega-model cohort from three South African studies (49 + 50",
      "+ 16 = 115 patients after exclusions; 1270 plasma nevirapine",
      "concentration samples). Study 1 (TB-interaction): rich sampling on",
      "n_ID = 25 (predose, 0.5, 1, 1.5, 2, 4, 6, 10 h after morning fasted",
      "dose; 0, 0.25, 0.75, 2, 10, 12 h after evening fed dose) plus sparse",
      "sampling on n_ID = 24 (two samples 0-12 h). Study 2 (DOT-HAART",
      "TDM): one sample per visit 0-12 h, one to four visits per subject 6",
      "months apart. Study 3 (artemether-lumefantrine interaction): rich",
      "sampling per visit (predose, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8 h morning",
      "fasted; predose, 1.5, 2, 3, 4, 5, 6, 8, 10, 12 h evening fed) plus",
      "5 sparse samples 14, 24, 96, 120, 144 h after first morning dose.",
      "Baseline demographics from Svensson 2012 Table 1 columns Study 1-3.",
      "An external validation cohort of 173 Dutch HIV-infected adults",
      "(Study 4) was used for VPC-based generalisation testing but not for",
      "parameter estimation; the cohort represented here is the 115",
      "subjects of the model-development set only."
    )
  )

  ini({
    # =========================================================================
    # Structural disposition (Svensson 2012 Table 2 'Parameter estimates').
    # CL/F parameters are reported at the cohort-anchor combination of FFM =
    # 42 kg and WT = 70 kg per Results 'Population pharmacokinetics'
    # paragraph 2. V/F has no IIV reported and is the same across the two
    # mixture sub-populations.
    # =========================================================================
    lcl_fast <- log(3.12);  label("Typical apparent oral CL/F in the fast-eliminator subpopulation (L/h, at FFM 42 kg reference)")  # Svensson 2012 Table 2 'CL/F pop. 1 (l h-1) = 3.12 (RSE 5.10%)'
    lcl_slow <- log(1.45);  label("Typical apparent oral CL/F in the slow-eliminator subpopulation (L/h, at FFM 42 kg reference)")  # Svensson 2012 Table 2 'CL/F pop. 2 (l h-1) = 1.45 (RSE 14.70%)'
    lvc      <- log(105);   label("Typical apparent volume of distribution V/F (L, at WT 70 kg reference)")                         # Svensson 2012 Table 2 'V/F (l) = 105 (RSE 4.90%)'

    # =========================================================================
    # Transit-compartment absorption (Svensson 2012 Methods 'Modelling',
    # citing Savic 2007 reference [35]). The paper fixes NTRANS = 2 transit
    # compartments and reports two MTT typical values selected per dose
    # occasion by the fed/fasted indicator (2.46 h fed vs 0.596 h fasted).
    # The Savic 2007 standard formulation places the (NTRANS + 1) = 3
    # first-order rate constants between depot, transit_1, transit_2, and
    # central all equal to ktr = (NTRANS + 1) / MTT, with MTT = mean transit
    # time from depot to central. The paper does not report a separate ka,
    # consistent with this shared-rate Savic implementation.
    # =========================================================================
    lmtt_fasted <- log(0.596);  label("Typical absorption mean transit time MTT in the fasted state (h, FED = 0)")  # Svensson 2012 Table 2 'MTT (h) fasted = 0.596 (RSE 8.70%)'
    e_fed_mtt   <- log(2.46 / 0.596);  label("Log-additive multiplicative shift on MTT for the fed state (FED = 1)")  # Svensson 2012 Table 2 'MTT (h) fed = 2.46 (RSE 7.50%)' / 'MTT (h) fasted = 0.596'; e_fed_mtt = log(2.46/0.596) = 1.418 reproduces the 4.13-fold slowing of MTT with food
    nn_fix      <- fixed(2);  label("Number of Savic-style transit compartments NTRANS (integer, unitless)")  # Svensson 2012 Results 'Population pharmacokinetics' paragraph 1: 'In the final model, the number of transit compartments was fixed to two without significant loss of goodness of fit.'

    # =========================================================================
    # Bioavailability (Svensson 2012 Table 2 'F (%) when TB treatment').
    # Reference F = 1.0 fixed (oral-only data; absolute F not identifiable);
    # concomitant TB treatment decreases F by 39%, encoded as a log-additive
    # multiplicative shift on the depot bioavailability.
    # =========================================================================
    lfdepot     <- fixed(log(1.0));  label("Reference bioavailability F (fixed at 1.0; oral-only data, absolute F not identifiable)")  # Svensson 2012 Methods 'Modelling' paragraph 6 and Discussion paragraph 5: oral-only data implies that an absolute F is not estimable; the 'F (%) when TB treatment = 61.30' value is a RELATIVE bioavailability ratio against the reference fasted-no-TB state
    e_tb_pos_fdepot <- log(0.613);   label("Log-additive multiplicative shift on F for concomitant TB treatment (TB_POS = 1)")  # Svensson 2012 Table 2 'F (%) when TB treatment = 61.30 (RSE 8.70%)'; e_tb_pos_fdepot = log(0.613) = -0.490 reproduces the 39% decrease

    # =========================================================================
    # Allometric exponents (Svensson 2012 Results 'Population pharmacokinetics'
    # paragraph 2: 'Estimating the allometric coefficients instead of using
    # 3/4 for CL and 1 for V fixed did not improve the OFV convincingly (-6.9
    # points, 2 degrees of freedom) and did not decrease the variability in
    # the parameters, hence scaling with the fixed coefficient was found to
    # be sufficient.'). Both exponents wrapped in fixed() per nlmixr2lib
    # convention.
    # =========================================================================
    e_ffm_cl <- fixed(0.75);  label("Allometric exponent of FFM on CL/F (unitless; fixed at 0.75 per Anderson-Holford)")  # Svensson 2012 Results 'Population pharmacokinetics' paragraph 2 (allometric coefficients fixed to customary 3/4 and 1)
    e_wt_vc  <- fixed(1.00);  label("Allometric exponent of WT on V/F (unitless; fixed at 1 per Anderson-Holford)")  # Svensson 2012 Results 'Population pharmacokinetics' paragraph 2

    # =========================================================================
    # Inter-individual / inter-occasion variability (Svensson 2012 Table 2).
    # The source paper reports CV% on the SD scale; nlmixr2lib uses
    # variances on the log scale via omega^2 = log(1 + CV^2). Following the
    # convention used in Bienczak_2016_nevirapine.R, Svensson_2018_
    # bedaquiline.R, and Svensson_2016_rifampicin.R, BOV is dropped when a
    # BSV term is reported on the same parameter and folded in as a BSV-
    # equivalent when only BOV is reported:
    #
    #   - CL/F: BSV 24.9% kept as etalcl_fast (shared across both mixture
    #     sub-populations; the source paper reports a single BSV after the
    #     mixture, applied to whichever class the subject was assigned to);
    #     no BOV reported.
    #   - F: BOV 26.9% folded in as BSV-equivalent etalfdepot. An additional
    #     34.1% BSV is reported specifically on the TB-treatment effect on
    #     F and is encoded as a separate eta gated on TB_POS (etalfdepot_tb).
    #   - MTT: BOV 64.0% folded in as BSV-equivalent etalmtt (no separate
    #     BSV on MTT reported).
    #
    # See vignette Assumptions and deviations for the BSV-vs-BOV folding
    # rationale and the cross-mixture-class CL eta convention.
    # =========================================================================
    etalcl_fast    ~ 0.06015  # Svensson 2012 Table 2 'BSV CL/F (%) = 24.90'; omega^2 = log(1 + 0.249^2) = 0.06015 (shared across mixture sub-populations)
    etalmtt        ~ 0.34331  # Svensson 2012 Table 2 'BOV MTT (%) = 64.00' folded in as BSV-equivalent; omega^2 = log(1 + 0.640^2) = 0.34331 (no separate BSV on MTT reported)
    etalfdepot     ~ 0.06986  # Svensson 2012 Table 2 'BOV F (%) = 26.90' folded in as BSV-equivalent; omega^2 = log(1 + 0.269^2) = 0.06986
    etalfdepot_tb  ~ 0.11000  # Svensson 2012 Table 2 'BSV F when TB treatment (%) = 34.10'; omega^2 = log(1 + 0.341^2) = 0.11000 (gated by TB_POS in model())

    # =========================================================================
    # Residual error (Svensson 2012 Table 2 'Proportional error (%) = 8.41').
    # Proportional-only on linear-scale concentration; no additive component
    # supported by the data.
    # =========================================================================
    propSd <- 0.0841;  label("Proportional residual error (fraction)")  # Svensson 2012 Table 2 'Proportional error (%) = 8.41 (RSE 5.40%)'
  })

  model({
    # --- 1. Derived covariate terms ----------------------------------------
    # Mixture-class gating (Sherwin 2012 risperidone idiom adapted to a
    # two-class binary). Only one of (1 - MIX_SLOW_ELIM_NVP) and
    # MIX_SLOW_ELIM_NVP is 1 per subject; the typical CL is selected
    # accordingly. The BSV eta etalcl_fast is shared across both classes
    # (single 24.9% BSV after the mixture, per Table 2).
    cl_typ <- exp(lcl_fast) * (1 - MIX_SLOW_ELIM_NVP) +
              exp(lcl_slow) * MIX_SLOW_ELIM_NVP

    # --- 2. Individual PK parameters ---------------------------------------
    # Apparent oral clearance: mixture-class typical value scaled
    # allometrically to fat-free mass (FFM reference 42 kg), with a single
    # log-normal BSV shared across the two mixture sub-populations.
    cl <- cl_typ * exp(etalcl_fast) * (FFM / 42)^e_ffm_cl

    # Apparent volume of distribution: allometric on body weight (WT
    # reference 70 kg). No IIV reported in Svensson 2012 Table 2.
    vc <- exp(lvc) * (WT / 70)^e_wt_vc

    # Mean transit time: fasted reference shifted multiplicatively for fed
    # occasions (log-additive form so that the typical fed MTT = 0.596 *
    # exp(log(2.46 / 0.596)) = 2.46 h matches Table 2). BSV-equivalent
    # eta from the folded-in BOV term.
    mtt <- exp(lmtt_fasted + e_fed_mtt * FED + etalmtt)
    nn  <- nn_fix
    ktr <- (nn + 1) / mtt  # Savic 2007 standard formulation: (NTRANS + 1) first-order rate constants at rate ktr, all shared; ktr = (NTRANS + 1) / MTT where MTT is the depot-to-central mean transit time

    # Bioavailability: reference F = 1.0 modulated by the TB-treatment
    # covariate, with a TB-conditional eta on the bioavailability path. The
    # BOV-folded-in etalfdepot acts on all dose records; the
    # TB-conditional etalfdepot_tb acts only when TB_POS = 1.
    fdepot <- exp(lfdepot + e_tb_pos_fdepot * TB_POS +
                  etalfdepot + etalfdepot_tb * TB_POS)

    # Micro-constant for systemic elimination.
    kel <- cl / vc

    # --- 3. ODE system -----------------------------------------------------
    # depot -> transit1 -> transit2 -> central -> elimination, all transit
    # steps at the shared Savic rate ktr = (NTRANS + 1) / MTT.
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot    - ktr * transit1
    d/dt(transit2) <-  ktr * transit1 - ktr * transit2
    d/dt(central)  <-  ktr * transit2 - kel * central

    # --- 4. Bioavailability ------------------------------------------------
    f(depot) <- fdepot

    # --- 5. Observation and residual error ---------------------------------
    # Dose in mg, V/F in L => Cc in mg/L (matches the source paper unit
    # convention and Table 1 / Figure 2 plot axes).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
