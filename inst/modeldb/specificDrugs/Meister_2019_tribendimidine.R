Meister_2019_tribendimidine <- function() {
  description <- "Pooled population PK model for the two tribendimidine metabolites dADT (deacetylated amidantel, the anthelminthically active species) and adADT (acetylated dADT) in Opisthorchis viverrini-infected Lao adolescents and adults given single oral doses of 25-600 mg (Meister 2019). Pools the two phase 2a ascending-dose trials of Vanobberghen 2016 (68 patients) with a phase 2b trial (125 patients). A Savic transit-compartment absorption model with a non-integer, log-normally distributed number of transit compartments (5.27 typical) feeds a one-compartment dADT disposition model, from which a fixed 65% of elimination is routed to a one-compartment adADT model (the remaining 35% is assumed renal). Allometric body-weight scaling (fixed 0.75 on clearances, 1 on volumes, reference 52 kg), a linear age effect on dADT clearance, and two absorption covariates: the 200-mg-versus-50-mg tablet formulation and the breaking of the enteric coating of a split 50-mg tablet. Systematic whole-blood and dried-blood-spot matrix conversion factors, each metabolite carrying its own matrix-specific residual error. Fitted on natural-log-transformed molar concentrations, so amounts are nmol and concentrations nmol/L."
  reference <- "Meister I, Assawasuwannakit P, Vanobberghen F, Penny MA, Odermatt P, Sayasone S, Huwyler J, Tarning J, Keiser J. Pooled population pharmacokinetic analysis of tribendimidine for the treatment of Opisthorchis viverrini infections. Antimicrob Agents Chemother. 2019;63(4):e01391-18. doi:10.1128/AAC.01391-18"
  vignette <- "Meister_2019_tribendimidine"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. The absorption chain is analytic (rxode2 `transit()`),
  # so the six explicit transit states of the predecessor model
  # (Vanobberghen_2016_tribendimidine) collapse into `depot`.
  compartmentData <- list(
    depot = list(
      analyte = "tribendimidine / dADT in transit (dosed prodrug; never measured)",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "dADT (deacetylated amidantel)",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    central_adadt = list(
      analyte = "adADT (acetylated dADT)",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters dADT clearance only, as the linear form (1 + slope * (AGE - 45)). Meister 2019 Table 2 footnote d defines the effect as a 'linear covariate relationship between age and CL/F dADT centered on the median age of 45 years', and Table 2 footnote a states the printed estimates are for 'a typical patient at 45 years of age'. Unlike the predecessor model (Vanobberghen_2016_tribendimidine), Meister 2019 retained NO age effect on adADT clearance. Because the form is linear rather than exponential it is only valid over the ages actually studied (15 to 79 years); extrapolating dADT clearance beyond roughly 129 years would make it non-positive.",
      source_name = "AGE"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling on both clearances (exponent fixed at 0.75) and both volumes (exponent fixed at 1), normalized to 52 kg. NOT a row of Meister 2019 Table 2 because the exponents were fixed a priori rather than estimated, but stated twice in the source: Materials and Methods 'Covariate analysis' ('Total body weight was implemented a priori as an allometric function, centered on the median body weight, on clearance and volume parameters simultaneously using a fixed exponent of 0.75 for clearance and 1 for volume'), and Table 2 footnote a, which fixes the reference patient at 52 kg. The pooled phase 2a median weight is 52 kg and the phase 2b median is 54 kg (Table 1); 52 kg is the value footnote a attaches to the printed estimates and is therefore the centering value used here.",
      source_name = "WEIGHT"
    ),
    FORM_TRI_TAB200 = list(
      description = "Tribendimidine tablet-strength formulation indicator (1 = 200-mg enteric-coated tablets; 0 = 50-mg enteric-coated tablets, the reference formulation)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (whole 50-mg enteric-coated tablets)",
      notes = "Multiplicative effect of +42.9% on mean absorption transit time relative to whole 50-mg tablets (Meister 2019 Table 2, 'Formulation on MTT (%)' = 42.9). The Results give the direction explicitly: 'a 42.9% slower mean absorption transit time for the 200-mg formulation than for the 50-mg formulation', which makes the 50-mg tablet the reference level. The authors attribute the delay to the 200-mg tablets floating in the stomach (Discussion, citing the in vitro physicochemical characterisation in reference 16). Unlike Vanobberghen 2016, Meister 2019 retained NO formulation effect on either central volume.",
      source_name = "FORM"
    ),
    FORM_TRI_SPLIT50 = list(
      description = "Split-tablet indicator for tribendimidine 50-mg enteric-coated tablets (1 = the administered 50-mg tablet was broken in half, destroying the enteric coating; 0 = whole tablet)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (whole tablet)",
      notes = "Multiplicative effect of -79.4% on mean absorption transit time relative to whole 50-mg tablets (Meister 2019 Table 2, 'Split 50-mg tablets on MTT (%)' = -79.4). The Results give the direction explicitly: 'a 79.4% faster mean absorption transit time for broken tablets than for whole 50-mg tablets'. In the source data only the 25-mg dose level of the second phase 2a trial used split tablets, so this indicator is 1 only when FORM_TRI_TAB200 is 0; the two effects were never observed in combination and the multiplicative composition used here is an extrapolation outside that cell. Registered as the sibling indicator that the FORM_TRI_TAB200 register entry anticipated: Vanobberghen 2016's split-tablet interaction model did not converge, so that earlier extraction pooled split tablets into its reference level.",
      source_name = "SPLIT"
    ),
    SAMPLE_WHOLEBLOOD = list(
      description = "Per-observation sampling-matrix indicator (1 = the concentration was measured in venous whole blood; 0 = venous plasma, the reference matrix)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (venous plasma)",
      notes = "Selects both the systematic matrix conversion factor applied to the predicted concentration and the matrix-specific residual error. Meister 2019 Table 2 footnote a fixes the reference matrix: the printed estimates are 'with drug concentrations measured in plasma'. Whole-blood samples were collected in the two phase 2a trials only; the phase 2b trial used dried blood spots exclusively (Table 1). Must be 0 on any record for which SAMPLE_DBS is 1.",
      source_name = "MATRIX"
    ),
    SAMPLE_DBS = list(
      description = "Per-observation sampling-matrix indicator (1 = the concentration was measured in a dried blood spot; 0 = venous plasma, the reference matrix)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (venous plasma)",
      notes = "Selects both the systematic matrix conversion factor applied to the predicted concentration and the matrix-specific residual error. Dried blood spots were collected from a fingertip capillary draw onto DMPK-C cards in all three trials and were the only matrix in the phase 2b trial (Meister 2019 Materials and Methods, 'PK sampling and analysis'). Meister 2019 reports that 'the difference between drug concentrations measured in whole blood and dried blood spots was not statistically significant', so this level differs from SAMPLE_WHOLEBLOOD only in its point estimates and its residual magnitude. Must be 0 on any record for which SAMPLE_WHOLEBLOOD is 1.",
      source_name = "MATRIX"
    )
  )

  # Screened during covariate model building but not retained in the final
  # model, so they are documented rather than implemented.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Listed among the covariates tested in Meister 2019 Materials and Methods 'Covariate analysis' ('Other covariates tested included age, sex, formulation, and the effect of breaking the administered tablets') but not retained by the forward-selection / backward-elimination procedure. 54% of the phase 2b participants and 51% of the pooled phase 2a participants were female (Table 1)."
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units = "mL/min",
      type = "continuous",
      notes = "Explicitly excluded from the analysis rather than screened and rejected: Meister 2019 Materials and Methods 'Covariate analysis' states 'Creatinine clearance data were not available for all patients and could not be considered in this analysis'. The rural setting of the phase 2b trial did not permit biochemical assessment (Materials and Methods, 'Patients, treatment, and study procedures'). The authors nonetheless attribute the retained age effect on dADT clearance to declining renal function with age (Discussion)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 193L,
    n_studies = 3L,
    age_range = "median 42 years (range 15-65) in the pooled phase 2a trials; median 48 years (range 15-79) in the phase 2b trial",
    weight_range = "median 52 kg (range 38-67) in the pooled phase 2a trials; median 54 kg (range 32-85) in the phase 2b trial",
    sex_female_pct = 53.4,
    renal_function = "Not assessed. Creatinine clearance was unavailable for part of the pooled data set and was not tested as a covariate.",
    disease_state = "Adolescents and adults aged 15 years and older with confirmed Opisthorchis viverrini infection, diagnosed by duplicate Kato-Katz thick smears on two stool samples. In the phase 2b trial 90% (n = 113) had a light infection intensity at baseline and the mean egg burden was 145.3 eggs per gram. Across all three trials 151/191 (79%) were cured at 21 days; 93% of the phase 2b participants were cured and the average egg reduction rate exceeded 99%.",
    dose_range = "Single oral doses of 25, 50, 100, 200, 400 and 600 mg tribendimidine. Phase 2a trial 1 (n = 31) gave 200, 400 and 600 mg as 200-mg enteric-coated tablets (n = 13, 9, 9); phase 2a trial 2 (n = 37) gave 25, 50, 100 and 200 mg as 50-mg enteric-coated tablets (n = 9, 9, 9, 10), the 25-mg dose being a split 50-mg tablet; the phase 2b trial (n = 125) gave 400 mg as 200-mg tablets (n = 123), with 2 patients wrongly dosed at 200 mg and modelled at the dose actually received.",
    regions = "Champasack district, Lao People's Democratic Republic",
    studies = "Two phase 2a single-ascending-dose trials (previously reported by Vanobberghen 2016 and Duthaler 2016) pooled with the PK substudy of a noninferiority randomized controlled phase 2b trial of tribendimidine 400 mg versus praziquantel 75 mg/kg conducted February-April 2014; ISRCTN Registry no. ISRCTN96948551.",
    notes = "Sampling differed by trial (Table 1). The phase 2a trials used dense venous sampling at 0, 1, 2, 3, 4, 4.5, 5, 6, 8, 10 and 24 h for whole blood and plasma, plus dried blood spots at 0, 1, 3, 4.5, 6, 10 h or 0, 2, 4, 5, 8, 24 h. The phase 2b trial used a WinPOPT-optimised sparse scheme of five dried blood spots per patient at 0.32, 2.00, 7.75, 8.00 and 30.0 h (supplemental Tables S1 and S2). Both metabolites were quantified by LC-MS/MS over an analytical range of 1 to 2,000 ng/mL. Concentrations were natural-log transformed and the two metabolites fitted sequentially by the PPP&D method; values below the LLOQ (13% for dADT and 21% for adADT across all matrices) were handled by Beal's M3 method. Parameter precision comes from a 1,000-replicate nonparametric bootstrap stratified by study."
  )

  ini({
    # ---- Absorption -------------------------------------------------------
    # Savic 2007 transit-compartment model implemented analytically. Meister
    # 2019 Results, 'Pooled population PK modeling': 'a full Stirling
    # approximation implementation of the transit compartment absorption model
    # (5.27 theoretical compartments) where the absorption rate constant was
    # assumed to equal the transit rate constant'. The number of compartments
    # is therefore an estimated, non-integer, log-normally distributed
    # parameter rather than a fixed chain length, which is why this model uses
    # rxode2's analytic transit() rather than the explicit six-state chain of
    # the predecessor Vanobberghen_2016_tribendimidine.
    #
    # lmtt is the reference-level MTT, i.e. a WHOLE 50-mg tablet. Both
    # covariate effects below are stated in the Results relative to the 50-mg
    # whole tablet, which makes that the reference level of the formulation
    # covariate. Table 2 footnote a instead describes the typical patient as
    # receiving 'the 200-mg formulation as whole tablets'; see the vignette
    # Errata for why the reference-level reading is the one implemented.
    lmtt <- log(3.18); label("Mean absorption transit time MTT for a whole 50-mg tablet (h)") # Table 2, 'MTT (h)' = 3.18 (95% CI 3.05-3.73)
    lntr <- log(5.27); label("Number of theoretical transit compartments (unitless)") # Table 2, 'No. of transit compartments' = 5.27 (95% CI 5.16-8.01)

    # Relative bioavailability, fixed to unity with between-subject
    # variability allowed (Materials and Methods, 'Structural and stochastic
    # model development': 'Relative bioavailability was fixed to unity, but
    # between-subject variability (BSV) was allowed').
    lfdepot <- fixed(log(1)); label("Relative bioavailability F (fraction)") # Table 2, 'F (%)' = 100 (fixed)

    # ---- dADT (deacetylated amidantel) disposition ------------------------
    # Apparent values, at the reference patient of Table 2 footnote a.
    lcl <- log(15.8); label("Apparent dADT clearance CL/F at 52 kg and 45 years (L/h)") # Table 2, dADT 'CL/FdADT (liters/h)' = 15.8 (95% CI 15.1-17.9)
    lvc <- log(88.8); label("Apparent dADT central volume Vc/F at 52 kg (L)") # Table 2, dADT 'V/FdADT (liters)' = 88.8 (95% CI 84.4-103)

    # ---- adADT (acetylated dADT) disposition ------------------------------
    # Apparent values that also absorb the assumed 65% metabolic fraction, so
    # they are CL/(F * fm) and V/(F * fm).
    lcl_adadt <- log(65.8); label("Apparent adADT clearance CL/F at 52 kg (L/h)") # Table 2, adADT 'CL/FadADT (liters/h)' = 65.8 (95% CI 59.6-87.0)
    lvc_adadt <- log(15.7); label("Apparent adADT central volume Vc/F at 52 kg (L)") # Table 2, adADT 'V/FadADT (liters)' = 15.7 (95% CI 11.5-18.1)

    # Fraction of dADT elimination routed to adADT. Assumed, not estimated.
    fm <- fixed(0.65); label("Fraction of dADT elimination converted to adADT (unitless)") # Materials and Methods, 'Structural and stochastic model development': 'The percentage of dADT cleared renally was fixed at 35%, based on a previously reported value (9), and the remaining 65% was assumed to be metabolized to adADT'

    # ---- Allometric exponents, fixed a priori -----------------------------
    e_wt_cl <- fixed(0.75); label("Allometric exponent on both clearances (unitless)") # Materials and Methods, 'Covariate analysis': 'a fixed exponent of 0.75 for clearance and 1 for volume'
    e_wt_vc <- fixed(1); label("Allometric exponent on both central volumes (unitless)") # Materials and Methods, 'Covariate analysis', as above

    # ---- Covariate effects ------------------------------------------------
    e_age_cl <- -0.0119; label("Linear age effect on dADT CL/F, per year older (fraction)") # Table 2, 'Age on CL/FdADT (%)' = -1.19 (95% CI -1.35 to -0.88), i.e. -1.19% per year
    e_tab200_mtt <- 0.429; label("200-mg-tablet effect on MTT, relative to a whole 50-mg tablet (fraction longer)") # Table 2, 'Formulation on MTT (%)' = 42.9 (95% CI 16.0-61.4)
    e_split50_mtt <- -0.794; label("Split-50-mg-tablet effect on MTT, relative to a whole 50-mg tablet (fraction shorter)") # Table 2, 'Split 50-mg tablets on MTT (%)' = -79.4 (95% CI -100 to -60.2)

    # Systematic sampling-matrix conversion factors, relative to the venous
    # plasma reference matrix of Table 2 footnote a.
    e_wb_cc <- -0.145; label("Whole-blood matrix conversion factor on dADT concentration, relative to plasma (fraction)") # Table 2, dADT 'Whole blood-to-plasma matrix conversion factor (%)' = -14.5 (95% CI -19.4 to -7.41)
    e_dbs_cc <- -0.137; label("DBS matrix conversion factor on dADT concentration, relative to plasma (fraction)") # Table 2, dADT 'DBS-to-plasma matrix conversion factor (%)' = -13.7 (95% CI -19.0 to -4.55)
    e_wb_cc_adadt <- 0.05; label("Whole-blood matrix conversion factor on adADT concentration, relative to plasma (fraction)") # Table 2, adADT 'Whole blood-to-plasma matrix conversion factor (%)' = 5.00 (95% CI 2.16-9.51)
    e_dbs_cc_adadt <- 0.07; label("DBS matrix conversion factor on adADT concentration, relative to plasma (fraction)") # Table 2, adADT 'DBS-to-plasma matrix conversion factor (%)' = 7.00 (95% CI -1.45 to 17.7)

    # ---- Between-subject variability --------------------------------------
    # Table 2 footnote b prints BSV as sqrt(exp(variance)) - 1, which inverts
    # to variance = log(1 + CV^2). Table 2 reports no correlations, so the
    # omega matrix is diagonal.
    etalfdepot ~ 0.132235 # Table 2, F '% CV for BSV' = 37.6 (95% CI 32.1-45.8); log(1 + 0.376^2) = 0.132235
    etalmtt ~ 0.256718 # Table 2, MTT '% CV for BSV' = 54.1 (95% CI 49.9-66.2); log(1 + 0.541^2) = 0.256718
    etalntr ~ 2.302585 # Table 2, transit-compartment '% CV for BSV' = 300 (95% CI 204-397); log(1 + 3.00^2) = 2.302585
    etalcl ~ 0.038075 # Table 2, dADT CL/F '% CV for BSV' = 19.7 (95% CI 16.6-22.5); log(1 + 0.197^2) = 0.038075
    etalvc ~ 0.064927 # Table 2, dADT V/F '% CV for BSV' = 25.9 (95% CI 20.1-28.5); log(1 + 0.259^2) = 0.064927
    etalcl_adadt ~ 0.852541 # Table 2, adADT CL/F '% CV for BSV' = 116 (95% CI 105-172); log(1 + 1.16^2) = 0.852541
    etalvc_adadt ~ 0.087836 # Table 2, adADT V/F '% CV for BSV' = 30.3 (95% CI 25.9-54.3); log(1 + 0.303^2) = 0.087836

    # ---- Residual unexplained variability ---------------------------------
    # Meister 2019 fitted natural-log-transformed concentrations with an
    # additive error on that scale, which is an exponential residual error on
    # the arithmetic scale, i.e. lnorm() in nlmixr2. Table 2 prints each RUV
    # as a % CV under the same footnote b transform as the BSV rows, so the
    # log-scale SD is sqrt(log(1 + CV^2)). One RUV per sampling matrix per
    # metabolite.
    expSdPlasma <- 0.460477; label("Log-scale residual SD for dADT in plasma (unitless)") # Table 2, dADT 'RUVplasma (% CV)' = 48.6 (95% CI 44.1-59.4); sqrt(log(1 + 0.486^2)) = 0.460477
    expSdWholeBlood <- 0.564023; label("Log-scale residual SD for dADT in whole blood (unitless)") # Table 2, dADT 'RUVwhole blood (% CV)' = 61.2 (95% CI 52.7-72.8); sqrt(log(1 + 0.612^2)) = 0.564023
    expSdDbs <- 0.622523; label("Log-scale residual SD for dADT in dried blood spots (unitless)") # Table 2, dADT 'RUVDBS (% CV)' = 68.8 (95% CI 54.0-81.7); sqrt(log(1 + 0.688^2)) = 0.622523
    expSdPlasma_adadt <- 0.379875; label("Log-scale residual SD for adADT in plasma (unitless)") # Table 2, adADT 'RUVplasma (CV%)' = 39.4 (95% CI 31.4-44.3); sqrt(log(1 + 0.394^2)) = 0.379875
    expSdWholeBlood_adadt <- 0.445891; label("Log-scale residual SD for adADT in whole blood (unitless)") # Table 2, adADT 'RUVwhole blood (% CV)' = 46.9 (95% CI 38.5-53.6); sqrt(log(1 + 0.469^2)) = 0.445891
    expSdDbs_adadt <- 0.444166; label("Log-scale residual SD for adADT in dried blood spots (unitless)") # Table 2, adADT 'RUVDBS (% CV)' = 46.7 (95% CI 37.1-63.7); sqrt(log(1 + 0.467^2)) = 0.444166
  })

  model({
    # ---- Covariate factors ------------------------------------------------
    # Age is linear and centered on the median 45 years. Both absorption
    # covariates are referenced to a WHOLE 50-mg tablet; only the 25-mg dose
    # level of phase 2a trial 2 used split tablets, so the two indicators were
    # never 1 simultaneously in the source data.
    cl_age <- 1 + e_age_cl * (AGE - 45)
    mtt_form <- (1 + e_tab200_mtt * FORM_TRI_TAB200) * (1 + e_split50_mtt * FORM_TRI_SPLIT50)

    # Sampling-matrix conversion factors. Venous plasma is the reference
    # matrix and carries a factor of exactly 1.
    matrix_cc <- 1 + e_wb_cc * SAMPLE_WHOLEBLOOD + e_dbs_cc * SAMPLE_DBS
    matrix_cc_adadt <- 1 + e_wb_cc_adadt * SAMPLE_WHOLEBLOOD + e_dbs_cc_adadt * SAMPLE_DBS

    # ---- Individual parameters --------------------------------------------
    # Allometric scaling applies to disposition only; neither absorption
    # parameter carries a weight term.
    mtt <- exp(lmtt + etalmtt) * mtt_form
    ntr <- exp(lntr + etalntr)
    fdepot <- exp(lfdepot + etalfdepot)

    cl <- exp(lcl + etalcl) * (WT / 52)^e_wt_cl * cl_age
    vc <- exp(lvc + etalvc) * (WT / 52)^e_wt_vc
    cl_adadt <- exp(lcl_adadt + etalcl_adadt) * (WT / 52)^e_wt_cl
    vc_adadt <- exp(lvc_adadt + etalvc_adadt) * (WT / 52)^e_wt_vc

    # ---- Micro-constants ---------------------------------------------------
    # The absorption rate constant equals the transit rate constant (Results,
    # 'Pooled population PK modeling'), and rxode2's transit() defines
    # ktr = (n + 1) / MTT, so ka takes that same value.
    ka <- (ntr + 1) / mtt
    kel <- cl / vc
    kel_adadt <- cl_adadt / vc_adadt

    # ---- ODE system --------------------------------------------------------
    # Amounts are in nmol and volumes in L, so concentrations come out in
    # nmol/L -- the molar scale the model was fitted on. Tribendimidine
    # degrades to dADT without enzymatic involvement and the dosed prodrug is
    # never measured, so the absorption chain is written in dADT-equivalent
    # moles at 1:1 stoichiometry and no molecular-weight factor appears
    # anywhere in the system. Event tables must therefore dose in nmol; the
    # vignette converts a milligram tribendimidine dose for the user.
    #
    # rxode2's transit(ntr, mtt, fdepot) returns the analytic Savic gamma-PDF
    # input rate into depot from the most recent dose, applying fdepot as the
    # bioavailability factor. f(depot) <- 0 suppresses the ordinary dose bolus
    # so transit() is the only input pathway; depot then empties at ka.
    d/dt(depot) <- transit(ntr, mtt, fdepot) - ka * depot

    # dADT: one-compartment disposition. Total elimination is kel; the split
    # into a renal share (1 - fm) and a metabolic share (fm) does not change
    # the dADT profile, only how much mass reaches adADT.
    d/dt(central) <- ka * depot - kel * central

    # adADT: one-compartment disposition formed from dADT at 1:1 molar
    # stoichiometry. Formation-rate-limited, so the terminal slope of adADT
    # tracks the dADT half-life.
    d/dt(central_adadt) <- fm * kel * central - kel_adadt * central_adadt

    f(depot) <- 0

    # ---- Observation and residual error ------------------------------------
    Cc <- central / vc * matrix_cc
    Cc_adadt <- central_adadt / vc_adadt * matrix_cc_adadt

    # One residual magnitude per sampling matrix, selected per observation
    # record. Plasma is the reference and is selected when neither indicator
    # is set.
    expSdCc <- expSdPlasma * (1 - SAMPLE_WHOLEBLOOD - SAMPLE_DBS) +
      expSdWholeBlood * SAMPLE_WHOLEBLOOD +
      expSdDbs * SAMPLE_DBS
    expSdCcAdadt <- expSdPlasma_adadt * (1 - SAMPLE_WHOLEBLOOD - SAMPLE_DBS) +
      expSdWholeBlood_adadt * SAMPLE_WHOLEBLOOD +
      expSdDbs_adadt * SAMPLE_DBS

    Cc ~ lnorm(expSdCc)
    Cc_adadt ~ lnorm(expSdCcAdadt)
  })
}
