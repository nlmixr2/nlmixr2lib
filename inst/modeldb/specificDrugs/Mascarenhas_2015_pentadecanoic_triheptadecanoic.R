Mascarenhas_2015_pentadecanoic_triheptadecanoic <- function() {
  description <- "Simultaneous two-analyte population PK model for the malabsorption blood test (MBT) in healthy subjects and subjects with cystic fibrosis and pancreatic insufficiency (Mascarenhas 2015). The MBT co-administers pentadecanoic acid (PA, C15:0), a free fatty acid absorbed without pancreatic lipase, and triheptadecanoic acid (THA), a triglyceride that must be hydrolysed by lipase to release heptadecanoic acid (HA, C17:0) before absorption; the postdose difference between the two analytes measures pancreatic-based fat absorption. Each analyte has its own 1-compartment disposition, its own Savic 2007 analytical transit-compartment absorption chain, and its own estimated fasting baseline concentration added to the model prediction. Allometric body-weight scaling on CL/F and V/F at a 70 kg reference (exponents fixed at 0.75 and 1). Healthy subjects are the bioavailability reference (F = 1); subjects with cystic fibrosis take a relative bioavailability that depends on whether pancreatic enzymes were administered and, for heptadecanoic acid only, on the timing of the enzyme dose relative to the meal. The heptadecanoic-acid random effects are modelled as scaled multiples of the pentadecanoic-acid random effects (correlations 0.974-0.999), and between-occasion variability is carried on mean transit time and on bioavailability for both analytes."
  reference <- "Mascarenhas MR, Mondick J, Barrett JS, Wilson M, Stallings VA, Schall JI. Malabsorption Blood Test: Assessing Fat Absorption in Patients With Cystic Fibrosis and Pancreatic Insufficiency. J Clin Pharmacol. 2015;55(8):854-865. doi:10.1002/jcph.484"
  vignette <- "Mascarenhas_2015_pentadecanoic_triheptadecanoic"

  # Concentrations are micromol/L and V/F is in L, so dose amounts must be
  # supplied in micromoles. Converting the paper's gram doses:
  #   depot       (PA) : umol = g_PA  / 242.40 * 1e6
  #   depot_hepta (HA) : umol = g_THA / 849.42 * 1e6 * 3
  # The factor of 3 is because each triheptadecanoic acid molecule carries
  # three heptadecanoic acids ("a triglyceride with 3 heptadecanoic acids",
  # Mascarenhas 2015 MBT Preparation). See the vignette for the check that
  # the final selected MBT (5.0 g PA + 5.5 g THA) is near-equimolar under
  # this conversion (20 627 umol PA vs 19 425 umol HA).
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  # Both analytes are dosed; declared explicitly because buildModelDb()'s
  # heuristic only recognises a compartment literally named "depot" and
  # would otherwise omit depot_hepta from the registry's dosing column.
  dosing <- c("depot", "depot_hepta")

  compartmentData <- list(
    depot         = list(analyte = "pentadecanoic acid", units = "umol", specimen = "administration site", verified = TRUE),
    central       = list(analyte = "pentadecanoic acid", units = "umol", specimen = "plasma", verified = TRUE),
    depot_hepta   = list(analyte = "heptadecanoic acid", units = "umol", specimen = "administration site", verified = TRUE),
    central_hepta = list(analyte = "heptadecanoic acid", units = "umol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling of CL/F (exponent fixed 0.75) and V/F (exponent fixed 1.0) for both analytes, normalized to a 70 kg reference weight. Mascarenhas 2015 Results, Population Pharmacokinetic Analysis: 'described by an allometric model, where the typical value of a model parameter was described as a function of individual body weight (WTi), normalized by a reference weight, which was 70 kg.' Cohort weights ranged 30.7-101.6 kg (Table 1).",
      source_name        = "WT"
    ),
    DIS_CF = list(
      description        = "Cystic fibrosis with pancreatic insufficiency indicator (1 = CF subject, 0 = healthy comparison subject)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy comparison subject; relative bioavailability fixed at 1 for both analytes)",
      notes              = "Selects which relative-bioavailability branch applies. Healthy subjects are the F = 1 reference by construction (Mascarenhas 2015 Results: 'Separate bioavailability fractions (relative to healthy subjects) were estimated for subjects with CF both with and without enzyme administration'). All CF subjects in this analysis had pancreatic insufficiency confirmed by fecal elastase 1 < 200 ug/g stool (Methods, Subjects).",
      source_name        = "CF"
    ),
    CONMED_PANCRELIPASE = list(
      description        = "Pancreatic enzyme replacement (pancrelipase) administered with the MBT on this occasion (1 = enzymes taken, 0 = no enzymes)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no pancreatic enzymes administered with the MBT)",
      notes              = "Occasion-level, not subject-level: the No Enzymes Protocol gave the same 6 CF subjects one MBT with and one without enzymes. Only meaningful when DIS_CF = 1; healthy subjects took no enzymes and are handled by the DIS_CF = 0 reference branch. The product was Creon 20 (delayed-release, enteric-coated) at 4-7 capsules = 80 000-140 000 lipase units, taken with the MBT and again with lunch (Methods, MBT Preparation). Mascarenhas 2015 Results reports that doses above the 80 000-unit standard 'did not yield an increase in bioavailability', so a binary indicator rather than a continuous lipase-unit column is the correct encoding for this dataset.",
      source_name        = "ENZ"
    ),
    CONMED_PANCRELIPASE_PRE30 = list(
      description        = "Pancreatic enzymes taken 30 minutes BEFORE the MBT (1 = yes, 0 = no)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (with the all-zero reference across the three timing indicators meaning enzymes taken WITH the MBT)",
      notes              = "Indicator variable t1 of the Mascarenhas 2015 final-model equation block. Acts on heptadecanoic-acid bioavailability only; the paper found no corresponding effect on pentadecanoic acid ('PA was absorbed similarly to healthy subjects regardless of timing of enzymes', Discussion), and the PA column of the equation block carries no timing term. Only meaningful when DIS_CF = 1 and CONMED_PANCRELIPASE = 1.",
      source_name        = "T1"
    ),
    CONMED_PANCRELIPASE_POST30 = list(
      description        = "Pancreatic enzymes taken 30 minutes AFTER the MBT (1 = yes, 0 = no)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (all-zero across the three timing indicators = enzymes taken WITH the MBT)",
      notes              = "Indicator variable t2 of the Mascarenhas 2015 final-model equation block. Heptadecanoic acid only; see CONMED_PANCRELIPASE_PRE30.",
      source_name        = "T2"
    ),
    CONMED_PANCRELIPASE_POST60 = list(
      description        = "Pancreatic enzymes taken 60 minutes AFTER the MBT (1 = yes, 0 = no)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (all-zero across the three timing indicators = enzymes taken WITH the MBT)",
      notes              = "Indicator variable t3 of the Mascarenhas 2015 final-model equation block. Heptadecanoic acid only; see CONMED_PANCRELIPASE_PRE30. This arm was discontinued after 9 of the planned 16 subjects once an interim review showed markedly reduced concentrations, which is why its effect (0.78) is the least precisely estimated of the three timing factors (95%CI 0.491-1.13).",
      source_name        = "T3"
    ),
    OCC = list(
      description        = "Occasion index for between-occasion variability (1-4)",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "The MBT was repeated on up to 4 separate occasions at least 5 days apart (Timing of Enzymes Protocol: 4 occasions; Reproducibility Protocol: 3; No Enzymes Protocol: 2). Decomposed inside model() into binary indicators oc1..oc4 that multiplex the between-occasion etas on mean transit time and on relative bioavailability, carried independently for each analyte. For a single-occasion simulation pass OCC = 1 so the first IOV eta applies.",
      source_name        = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a post hoc power covariate on CL/F for both analytes but NOT retained in the final model. Mascarenhas 2015 Results: 'Adding these effect models resulted in no drop in the objective function value when compared with the reported model. Estimated exponents (95%CI) for the age effect were -0.421 (-1.11, 0.267) for PA CL/F and ... (-0.995, 0.467) for HA CL/F. Although the 95%CIs include the null value of zero, the low precision of the estimates indicates insufficient information in the data set to accurately characterize any age effects on CL/F.' Documented here to preserve the provenance of the covariate screen; both confidence intervals span zero, so no age term is implemented."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 60L,
    n_studies      = 5L,
    age_range      = "9.9-49.9 years (CF 9.9-24.2; healthy 19.0-49.9)",
    age_median     = "CF 15.6 years; healthy 28.2 years",
    weight_range   = "30.7-101.6 kg (CF 30.7-74.5; healthy 60.0-101.6)",
    weight_median  = "CF 50.5 kg; healthy 71.6 kg",
    sex_female_pct = 50.0,
    disease_state  = "33 subjects with cystic fibrosis and pancreatic insufficiency (fecal elastase 1 < 200 ug/g stool, FEV1 >= 40% predicted, no fibrosing colonopathy or significant bowel resection) and 27 healthy comparison subjects (BMI 21-30 kg/m2, no chronic illness affecting nutrient absorption).",
    dose_range     = "Single oral MBT test meal (8 oz, 550 kcal, 32 g fat) containing 2.5 or 5.0 g pentadecanoic acid and 5.0, 5.5 or 8.0 g triheptadecanoic acid. Healthy Orlistat Protocol 2.5 g PA + 8.0 g THA; healthy Timing Protocol 5.0 g PA + 5.5 g THA; CF No Enzymes Protocol 2.5 g PA + 5.0 g (n=3) or 8.0 g (n=3) THA; CF Timing and Reproducibility Protocols 5.0 g PA + 5.5 g THA. 5.0 g PA + 5.5 g THA was selected as the final MBT.",
    regions        = "United States (Children's Hospital of Philadelphia and Pennsylvania Presbyterian Medical Center)",
    co_medication  = "CF subjects received Creon 20 pancrelipase, 4-7 capsules (80 000-140 000 lipase units): 52% took 80 000 units, 15% 100 000, 21% 120 000, 12% 140 000. Enzymes were taken with the MBT and again with the 6-hour lunch meal, except on the no-enzyme occasion and on the timing-arm occasions where the dose was displaced by -30, +30 or +60 minutes.",
    studies        = "Pooled analysis across 5 protocols: CF No Enzymes (n=6, 2 occasions), CF Timing of Enzymes (n=16, 3-4 occasions), CF Reproducibility (n=11, 3 occasions), healthy Orlistat (n=15, pre-Orlistat MBT only), healthy Timing of Enzymes comparison group (n=12).",
    notes          = "Baseline demographics from Mascarenhas 2015 Table 1. Plasma sampled hourly over 8 hours starting with a premeal baseline sample, after a 12-hour overnight fast; subjects abstained from alcohol and dairy for 24 hours beforehand. Estimation in NONMEM VII level 2.0 with FOCE-INT. Parameter precision from a 500-replicate nonparametric bootstrap stratified by subject status and enzyme administration method. sex_female_pct is the pooled figure across both cohorts (CF 45%, healthy 56%)."
  )

  ini({
    # -----------------------------------------------------------------
    # Structural parameters -- pentadecanoic acid (PA, C15:0).
    # PA is the unsuffixed parent analyte. All values are the "PA
    # Estimate" column of Mascarenhas 2015 Table 2. Concentrations are
    # micromol/L. The extracted text layer of Table 2 reads "mmol/L"
    # because the Symbol-font micro sign is mis-decoded as a Latin "m";
    # rendering the Table 2 units column as an image shows "umol/L"
    # directly. The same page proves the text layer is lossy in a second
    # way -- it renders the Ka unit "h^-1" as "h1", dropping the minus
    # sign. Two independent checks agree with the rendered glyph: the
    # assay QC concentrations in Methods (Sample Analysis) are in mg/dL,
    # and the highest PA control of 6.70 mg/dL is
    # 6.70 * 10 / 242.40 g/mol = 276 umol/L, matching the scale of the
    # Figure 2 profiles; and a fasting plasma fatty-acid baseline of
    # 24.9 mmol/L is not physiological.
    # -----------------------------------------------------------------
    lcl    <- log(9.66);  label("Apparent oral clearance CL/F of pentadecanoic acid at 70 kg (L/h)")            # Table 2 PA CL/F = 9.66 L/h (RSE 18.5%, 95%CI 8.14-12.4)
    lvc    <- log(24.5);  label("Apparent central volume of distribution V/F of pentadecanoic acid at 70 kg (L)") # Table 2 PA V/F = 24.5 L (RSE 27.6%, 95%CI 4.89-32.5)
    lrbase <- log(24.9);  label("Fasting baseline plasma pentadecanoic acid concentration (umol/L)")            # Table 2 PA C0 = 24.9 umol/L (RSE 12.1%, 95%CI 22.7-27.4)
    lmtt   <- log(0.817); label("Mean absorption transit time MTT of pentadecanoic acid (h)")                   # Table 2 PA MTT = 0.817 h (RSE 7.3%, 95%CI 0.751-2.84)
    lka    <- log(0.266); label("First-order absorption rate constant Ka of pentadecanoic acid (1/h)")          # Table 2 PA Ka = 0.266 1/h (RSE 17.6%, 95%CI 0.206-0.440)
    lnn    <- log(6.96);  label("Number of absorption transit compartments N for pentadecanoic acid (continuous, dimensionless)") # Table 2 PA N = 6.96 (RSE 39.1%, 95%CI 0.195-18.2)

    # -----------------------------------------------------------------
    # Structural parameters -- heptadecanoic acid (HA, C17:0), the
    # lipase-dependent co-analyte released from triheptadecanoic acid.
    # "HA Estimate" column of Table 2.
    # -----------------------------------------------------------------
    lcl_hepta    <- log(16.3);  label("Apparent oral clearance CL/F of heptadecanoic acid at 70 kg (L/h)")            # Table 2 HA CL/F = 16.3 L/h (RSE 15.8%, 95%CI 13.9-20.7)
    lvc_hepta    <- log(14.3);  label("Apparent central volume of distribution V/F of heptadecanoic acid at 70 kg (L)") # Table 2 HA V/F = 14.3 L (RSE 64.1%, 95%CI 1.35-37.2)
    lrbase_hepta <- log(27.4);  label("Fasting baseline plasma heptadecanoic acid concentration (umol/L)")            # Table 2 HA C0 = 27.4 umol/L (RSE 9.5%, 95%CI 25.3-29.7)
    lmtt_hepta   <- log(3.52);  label("Mean absorption transit time MTT of heptadecanoic acid (h)")                   # Table 2 HA MTT = 3.52 h (RSE 12.2%, 95%CI 3.21-4.33)
    lka_hepta    <- log(0.307); label("First-order absorption rate constant Ka of heptadecanoic acid (1/h)")          # Table 2 HA Ka = 0.307 1/h (RSE 17.2%, 95%CI 0.253-0.563)
    lnn_hepta    <- log(7.08);  label("Number of absorption transit compartments N for heptadecanoic acid (continuous, dimensionless)") # Table 2 HA N = 7.08 (RSE 20.6%, 95%CI 5.79-8.49)

    # -----------------------------------------------------------------
    # Allometric exponents. Paper-fixed, no uncertainty reported.
    # -----------------------------------------------------------------
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL/F for both analytes (unitless)") # Results, Population Pharmacokinetic Analysis: allometric model, exponent 0.75 on CL/F, reference weight 70 kg
    e_wt_vc <- fixed(1.0);  label("Allometric exponent on V/F for both analytes (unitless)")  # Results, Population Pharmacokinetic Analysis: allometric model, exponent 1.0 on V/F, reference weight 70 kg

    # -----------------------------------------------------------------
    # Relative bioavailability. Healthy subjects are the reference
    # (F = 1 for both analytes, fixed by construction, not estimated).
    # The CF branches below are relative to that reference.
    # -----------------------------------------------------------------
    e_cf_fdepot           <- 1.07;   label("Relative bioavailability of pentadecanoic acid in CF without enzymes (fraction of healthy)")   # Table 2 PA F_CF = 1.07 (RSE 26.0%, 95%CI 0.827-1.42) = theta7
    e_cf_enz_fdepot       <- 0.877;  label("Relative bioavailability of pentadecanoic acid in CF with enzymes taken WITH the MBT (fraction of healthy)") # Table 2 PA F_CF,ENZ = 0.877 (RSE 30.0%, 95%CI 0.720-1.09) = theta8
    e_cf_fdepot_hepta     <- 0.0292; label("Relative bioavailability of heptadecanoic acid in CF without enzymes (fraction of healthy)")    # Table 2 HA F_CF = 0.0292 (RSE 67.8%, 95%CI 0.0192-0.0459) = theta15
    e_cf_enz_fdepot_hepta <- 0.606;  label("Relative bioavailability of heptadecanoic acid in CF with enzymes taken WITH the MBT (fraction of healthy)") # Table 2 HA F_CF,ENZ = 0.606 (RSE 26.9%, 95%CI 0.482-0.823) = theta16

    # Enzyme-timing effects on heptadecanoic acid bioavailability, as
    # multiplicative factors relative to enzymes taken WITH the MBT.
    # No corresponding PA terms: the PA column of the equation block
    # carries no timing effect (Discussion: "PA was absorbed similarly
    # to healthy subjects regardless of timing of enzymes").
    e_pre30_fdepot_hepta  <- 0.911; label("Heptadecanoic acid bioavailability factor, enzymes 30 min BEFORE the MBT (relative to enzymes with the MBT)") # Table 2 HA F_T1 = 0.911 (RSE 21.0%, 95%CI 0.710-1.12) = theta17
    e_post30_fdepot_hepta <- 0.829; label("Heptadecanoic acid bioavailability factor, enzymes 30 min AFTER the MBT (relative to enzymes with the MBT)")  # Table 2 HA F_T2 = 0.829 (RSE 25.0%, 95%CI 0.664-0.979) = theta18
    e_post60_fdepot_hepta <- 0.78;  label("Heptadecanoic acid bioavailability factor, enzymes 60 min AFTER the MBT (relative to enzymes with the MBT)")  # Table 2 HA F_T3 = 0.78 (RSE 26.7%, 95%CI 0.491-1.13) = theta19

    # -----------------------------------------------------------------
    # Random-effect sharing factors. Mascarenhas 2015 Results: "PA and
    # HA random-effects parameters were highly correlated for CL/F
    # (0.974), V/F (0.987), and baseline concentrations (0.999).
    # Therefore, random-effects parameters for HA were modeled as a
    # fraction of the corresponding PA random effects." The equation
    # block writes these as exp(eta_1 * theta20) etc. Same idiom as
    # Prytula_2016_tacrolimus.R (q_eta_scale) and
    # Hirt_2009_efavirenz.R (vc_eta_scale).
    # -----------------------------------------------------------------
    cl_hepta_eta_scale    <- 0.994; label("Scaling factor relating the heptadecanoic-acid CL/F random effect to the pentadecanoic-acid CL/F random effect (unitless)")       # Table 2 theta20 = 0.994 (RSE 9.6%, 95%CI 0.888-1.17)
    vc_hepta_eta_scale    <- 1.3;   label("Scaling factor relating the heptadecanoic-acid V/F random effect to the pentadecanoic-acid V/F random effect (unitless)")         # Table 2 theta21 = 1.3 (RSE 65.5%, 95%CI 0.870-2.43)
    rbase_hepta_eta_scale <- 0.998; label("Scaling factor relating the heptadecanoic-acid baseline random effect to the pentadecanoic-acid baseline random effect (unitless)") # Table 2 theta22 = 0.998 (RSE 19.9%, 95%CI 0.890-1.10)

    # -----------------------------------------------------------------
    # Between-subject variability: a full 3x3 block on the PA CL/F,
    # V/F and baseline etas ("A full covariance matrix was estimated to
    # capture random-effects correlation", Methods). Values are the
    # "Interindividual variance" rows of Table 2, in lower-triangle
    # order (CL, CL-V, V, CL-S0, V-S0, S0) -- exactly the order the
    # table prints them.
    #
    # The two covariances involving the baseline are NEGATIVE. Text
    # extraction of this PDF drops the leading minus sign (the same
    # defect that turns the Ka unit "h^-1" into "h1" and the F_T1
    # description "enzymes at -0.5 h" into "enzymes at 0.5 h"), so the
    # extracted values read as +0.0894 and +0.151. Rendering Table 2 as
    # an image shows the printed cells directly: estimates -0.0894 and
    # -0.151, RSEs -51.7% and -53.8%, and 95%CIs (-0.162, -0.0259) and
    # (-0.231, -0.0899). The text layer alone also settles it without
    # any rendering, because those two CIs extract as "(0.162, 0.0259)"
    # and "(0.231, 0.0899)" -- DESCENDING, and therefore impossible for
    # a 2.5th-97.5th bootstrap quantile pair. The signs are
    # mechanistically sensible (higher clearance -> lower baseline
    # concentration) and the resulting matrix is positive-definite
    # (leading minors 0.315, 0.0812, determinant 0.00607).
    # -----------------------------------------------------------------
    etalcl + etalvc + etalrbase ~ c(0.315,
                                    0.162, 0.341,
                                    -0.0894, -0.151, 0.143)
    # Table 2: V2CL = 0.315 (RSE 47.6%, 95%CI 0.188-0.416); sqrt = 56.1% = printed BSV
    # Table 2: V2CL-V = 0.162 (RSE 51.7%, 95%CI -0.0666-0.241)
    # Table 2: V2V = 0.341 (RSE 55.4%, 95%CI 0.230-0.991); sqrt = 58.4% = printed BSV
    # Table 2: V2CL-S0 = -0.0894 (RSE -51.7%, 95%CI -0.162 to -0.0259)
    # Table 2: V2V-S0 = -0.151 (RSE -53.8%, 95%CI -0.231 to -0.0899)
    # Table 2: V2S0 = 0.143 (RSE 44.9%, 95%CI 0.0892-0.192); sqrt = 37.8% = printed BSV

    # -----------------------------------------------------------------
    # Between-occasion variability on MTT and on relative
    # bioavailability, for each analyte. Encoded as the
    # occasion-indicator expansion (one eta per occasion, occasions
    # 2-4 fixed to the occasion-1 magnitude to represent NONMEM
    # $OMEGA BLOCK(1) SAME), because rxode2 parses but cannot SIMULATE
    # the `eta ~ var | occ` multi-level syntax. Up to 4 occasions
    # (Timing of Enzymes Protocol).
    #
    # Each variance below is confirmed by the printed BOV percentage:
    # sqrt(0.205) = 45.3%, sqrt(0.102) = 31.9%, sqrt(0.0925) = 30.4%,
    # sqrt(0.136) = 36.9% -- matching Results: "For mean absorption
    # transit time (MTT), between-occasion variability was estimated to
    # be 45.3% and 31.9%, and for bioavailability it was 30.4% and
    # 36.9% for PA and HA, respectively."
    # -----------------------------------------------------------------
    etaiov_mtt_1 ~ 0.205        ; # Table 2 PA V2MTT = 0.205 (95%CI 0.0471-0.307); sqrt = 45.3% = printed PA BOV
    etaiov_mtt_2 ~ fixed(0.205) ; # NONMEM $OMEGA BLOCK(1) SAME
    etaiov_mtt_3 ~ fixed(0.205) ; # NONMEM $OMEGA BLOCK(1) SAME
    etaiov_mtt_4 ~ fixed(0.205) ; # NONMEM $OMEGA BLOCK(1) SAME

    etaiov_fdepot_1 ~ 0.0925        ; # Table 2 PA V2F = 0.0925 (95%CI 0.0565-0.127); sqrt = 30.4% = printed PA BOV
    etaiov_fdepot_2 ~ fixed(0.0925) ; # NONMEM $OMEGA BLOCK(1) SAME
    etaiov_fdepot_3 ~ fixed(0.0925) ; # NONMEM $OMEGA BLOCK(1) SAME
    etaiov_fdepot_4 ~ fixed(0.0925) ; # NONMEM $OMEGA BLOCK(1) SAME

    etaiov_mtt_hepta_1 ~ 0.102        ; # Table 2 HA V2MTT = 0.102 (95%CI 0.0725-0.127); sqrt = 31.9% = printed HA BOV
    etaiov_mtt_hepta_2 ~ fixed(0.102) ; # NONMEM $OMEGA BLOCK(1) SAME
    etaiov_mtt_hepta_3 ~ fixed(0.102) ; # NONMEM $OMEGA BLOCK(1) SAME
    etaiov_mtt_hepta_4 ~ fixed(0.102) ; # NONMEM $OMEGA BLOCK(1) SAME

    etaiov_fdepot_hepta_1 ~ 0.136        ; # Table 2 HA V2F = 0.136 (95%CI 0.0411-0.189); sqrt = 36.9% = printed HA BOV
    etaiov_fdepot_hepta_2 ~ fixed(0.136) ; # NONMEM $OMEGA BLOCK(1) SAME
    etaiov_fdepot_hepta_3 ~ fixed(0.136) ; # NONMEM $OMEGA BLOCK(1) SAME
    etaiov_fdepot_hepta_4 ~ fixed(0.136) ; # NONMEM $OMEGA BLOCK(1) SAME

    # -----------------------------------------------------------------
    # Residual error. Methods specify a proportional error model for
    # both analytes; Table 2 reports the residual VARIANCE, so the SD
    # is its square root (and matches the printed percentage).
    # -----------------------------------------------------------------
    propSd       <- 0.3137; label("Proportional residual error for pentadecanoic acid (fraction)") # Table 2 PA delta2_prop = 0.0984 (RSE 12.0%, 95%CI 0.0877-0.118); sqrt(0.0984) = 0.3137 = printed 31.4%
    propSd_hepta <- 0.2782; label("Proportional residual error for heptadecanoic acid (fraction)")  # Table 2 HA delta2_prop = 0.0774 (RSE 10.1%, 95%CI 0.0674-0.0887); sqrt(0.0774) = 0.2782 = 27.8%
  })

  model({
    # ---------------------------------------------------------------
    # 1. Occasion indicators and the between-occasion random effects
    #    they multiplex. OCC = 1..4; any other value zeroes all four
    #    indicators and leaves the subject at the typical MTT and F.
    # ---------------------------------------------------------------
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)

    iov_mtt    <- oc1 * etaiov_mtt_1    + oc2 * etaiov_mtt_2 +
                  oc3 * etaiov_mtt_3    + oc4 * etaiov_mtt_4
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 +
                  oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4

    iov_mtt_hepta    <- oc1 * etaiov_mtt_hepta_1    + oc2 * etaiov_mtt_hepta_2 +
                        oc3 * etaiov_mtt_hepta_3    + oc4 * etaiov_mtt_hepta_4
    iov_fdepot_hepta <- oc1 * etaiov_fdepot_hepta_1 + oc2 * etaiov_fdepot_hepta_2 +
                        oc3 * etaiov_fdepot_hepta_3 + oc4 * etaiov_fdepot_hepta_4

    # ---------------------------------------------------------------
    # 2. Individual parameters. Allometric size on CL/F and V/F for
    #    both analytes at the 70 kg reference. The heptadecanoic-acid
    #    etas are the pentadecanoic-acid etas scaled by theta20-22
    #    (correlations 0.974-0.999 in the source fit), so the two
    #    analytes share three between-subject random effects rather
    #    than carrying six.
    # ---------------------------------------------------------------
    wt70 <- WT / 70

    cl    <- exp(lcl + etalcl) * wt70^e_wt_cl
    vc    <- exp(lvc + etalvc) * wt70^e_wt_vc
    rbase <- exp(lrbase + etalrbase)
    mtt   <- exp(lmtt + iov_mtt)
    ka    <- exp(lka)
    nn    <- exp(lnn)

    cl_hepta    <- exp(lcl_hepta + etalcl * cl_hepta_eta_scale) * wt70^e_wt_cl
    vc_hepta    <- exp(lvc_hepta + etalvc * vc_hepta_eta_scale) * wt70^e_wt_vc
    rbase_hepta <- exp(lrbase_hepta + etalrbase * rbase_hepta_eta_scale)
    mtt_hepta   <- exp(lmtt_hepta + iov_mtt_hepta)
    ka_hepta    <- exp(lka_hepta)
    nn_hepta    <- exp(lnn_hepta)

    # ---------------------------------------------------------------
    # 3. Relative bioavailability.
    #
    #    Healthy subjects (DIS_CF = 0) are the reference, F = 1.
    #    CF subjects take one of two branches depending on whether
    #    enzymes were administered on this occasion. The
    #    between-occasion eta on F sits on the enzyme branch only,
    #    following both the equation block (which places
    #    eta_IOV,F on F_CF+Enz, not on F_CF) and Methods ("Random
    #    effects for between-occasion variability in absorption mean
    #    transit time and bioavailability FOR CF SUBJECTS WITH
    #    ADMINISTRATION OF ENZYMES were also included").
    #
    #    For heptadecanoic acid the enzyme branch is additionally
    #    multiplied by the three timing factors in power form, so the
    #    reference timing (enzymes with the MBT; all three indicators
    #    zero) contributes a factor of exactly 1.
    # ---------------------------------------------------------------
    fdepot <- (1 - DIS_CF) +
      DIS_CF * ((1 - CONMED_PANCRELIPASE) * e_cf_fdepot +
                CONMED_PANCRELIPASE * e_cf_enz_fdepot * exp(iov_fdepot))

    fdepot_hepta <- (1 - DIS_CF) +
      DIS_CF * ((1 - CONMED_PANCRELIPASE) * e_cf_fdepot_hepta +
                CONMED_PANCRELIPASE * e_cf_enz_fdepot_hepta *
                  e_pre30_fdepot_hepta^CONMED_PANCRELIPASE_PRE30 *
                  e_post30_fdepot_hepta^CONMED_PANCRELIPASE_POST30 *
                  e_post60_fdepot_hepta^CONMED_PANCRELIPASE_POST60 *
                  exp(iov_fdepot_hepta))

    # ---------------------------------------------------------------
    # 4. ODE system: two independent 1-compartment dispositions, each
    #    fed by its own Savic 2007 analytical transit chain.
    #    rxode2's transit(n, mtt, bio) is compartment-aware -- it reads
    #    podo() for the compartment whose d/dt() contains the call and
    #    fires only when a dose targets that compartment -- so the two
    #    chains stay independent. f(<depot>) <- 0 suppresses the bolus
    #    so the transit chain alone delivers each dose.
    #
    #    Dose amounts are in micromoles; see the `units` note above for
    #    the gram-to-micromole conversion, including the factor of 3
    #    for triheptadecanoic acid -> heptadecanoic acid.
    # ---------------------------------------------------------------
    d/dt(depot)         <- transit(nn, mtt, fdepot) - ka * depot
    d/dt(central)       <- ka * depot - cl / vc * central

    d/dt(depot_hepta)   <- transit(nn_hepta, mtt_hepta, fdepot_hepta) - ka_hepta * depot_hepta
    d/dt(central_hepta) <- ka_hepta * depot_hepta - cl_hepta / vc_hepta * central_hepta

    f(depot)       <- 0
    f(depot_hepta) <- 0

    # ---------------------------------------------------------------
    # 5. Observations. The estimated fasting baseline is an ADDITIVE
    #    constant offset, not a decaying initial condition: these are
    #    endogenous odd-chain fatty acids continuously present in
    #    plasma, and Methods describes the baseline as accounting "for
    #    variability in concentration of fats at baseline prior to MBT
    #    administration". The alternative reading (central(0) = C0*V,
    #    decaying at kel) is falsified by Figure 2, where the
    #    heptadecanoic-acid median is flat at ~28 umol/L from t = 0 to
    #    t = 1 h; with kel_HA = 16.3/14.3 = 1.14 /h a decaying initial
    #    condition would have fallen to ~8.8 umol/L by 1 h.
    # ---------------------------------------------------------------
    Cc       <- central / vc + rbase
    Cc_hepta <- central_hepta / vc_hepta + rbase_hepta

    Cc       ~ prop(propSd)
    Cc_hepta ~ prop(propSd_hepta)
  })
}
