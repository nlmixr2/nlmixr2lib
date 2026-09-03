Keunecke_2020_regorafenib_phase1 <- function() {
  description <- "Joint parent + metabolite population PK model for oral regorafenib and its two pharmacologically active metabolites M-2 (N-oxide, BAY 75-7495) and M-5 (N-oxide N-desmethyl, BAY 81-8752) in patients with advanced solid tumours, fitted to densely sampled phase 1 data (Keunecke 2020 Table 2). Each of the three analytes carries its own two-compartment disposition plus a gallbladder reservoir, giving enterohepatic circulation: drug moves from the central compartment into the gallbladder at a constant first-order rate kbm, and is released back essentially instantaneously (kehc = 100 1/h) during three 30-minute post-prandial windows per day anchored to breakfast, lunch and dinner. M-2 is formed pre-systemically, so the dose splits at absorption into a two-transit parent chain and a three-transit M-2 chain; M-5 is formed systemically from the parent central compartment as a fraction fm_m5 of parent clearance. Metabolite volumes, inter-compartmental flow and gallbladder rate constants are fixed to the parent values because no intravenous or direct-metabolite data were available. All parent volumes and clearances are apparent, relative to oral bioavailability AND to the parent fraction (1 - fm_m2), exactly as tabulated by the authors; the model therefore delivers the full nominal dose to the parent chain via f(depot) = 1 / (1 - fm_m2)."
  reference <- paste(
    "Keunecke A, Hoefman S, Drenth H-J, Zisowsky J, Cleton A, Ploeger BA.",
    "Population pharmacokinetics of regorafenib in solid tumours:",
    "Exposure in clinical practice considering enterohepatic circulation",
    "and food intake.",
    "Br J Clin Pharmacol. 2020;86(12):2362-2376.",
    "doi:10.1111/bcp.14334.",
    sep = " "
  )
  vignette <- "Keunecke_2020_regorafenib"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The phase 1 model carries no covariate effects; sex and BMI enter only in
  # the phase 3 covariate model (Keunecke_2020_regorafenib_phase3).
  covariateData <- list()

  compartmentData <- list(
    depot            = list(analyte = "regorafenib", units = "mg", specimen = "administration site", verified = TRUE),
    transit1         = list(analyte = "regorafenib", units = "mg", specimen = "administration site", verified = TRUE),
    transit2         = list(analyte = "regorafenib", units = "mg", specimen = "administration site", verified = TRUE),
    central          = list(analyte = "regorafenib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1      = list(analyte = "regorafenib", units = "mg", specimen = "plasma", verified = TRUE),
    gallbladder      = list(analyte = "regorafenib", units = "mg", specimen = "bile", verified = TRUE),
    transit1_m2      = list(analyte = "regorafenib M-2 (N-oxide, BAY 75-7495)", units = "mg", specimen = "administration site", verified = TRUE),
    transit2_m2      = list(analyte = "regorafenib M-2 (N-oxide, BAY 75-7495)", units = "mg", specimen = "administration site", verified = TRUE),
    transit3_m2      = list(analyte = "regorafenib M-2 (N-oxide, BAY 75-7495)", units = "mg", specimen = "administration site", verified = TRUE),
    central_m2       = list(analyte = "regorafenib M-2 (N-oxide, BAY 75-7495)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_m2   = list(analyte = "regorafenib M-2 (N-oxide, BAY 75-7495)", units = "mg", specimen = "plasma", verified = TRUE),
    gallbladder_m2   = list(analyte = "regorafenib M-2 (N-oxide, BAY 75-7495)", units = "mg", specimen = "bile", verified = TRUE),
    central_m5       = list(analyte = "regorafenib M-5 (N-oxide N-desmethyl, BAY 81-8752)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_m5   = list(analyte = "regorafenib M-5 (N-oxide N-desmethyl, BAY 81-8752)", units = "mg", specimen = "plasma", verified = TRUE),
    gallbladder_m5   = list(analyte = "regorafenib M-5 (N-oxide N-desmethyl, BAY 81-8752)", units = "mg", specimen = "bile", verified = TRUE)
  )

  population <- list(
    species     = "human",
    n_subjects  = 62L,
    n_studies   = 2L,
    age_range   = "not reported for the phase 1 subset",
    disease_state = "Adults with advanced solid tumours.",
    dose_range  = "Regorafenib 160 mg once daily (4 x 40 mg tablets), 3 weeks on / 1 week off; study 15823 additionally gave a single 160 mg dose on Day 1.",
    regions     = "USA (study 14814) and Mainland China (study 15823).",
    notes       = "Keunecke 2020 section 2.1.1 and Appendix Table A1. Study 14814 (NCT01339104, USA, n = 44, 774 PK observations) sampled cycle 1 day 21 predose, up to 10 h post-dose and 24 h post-dose. Study 15823 (NCT02398513, Mainland China, n = 18, 566 PK observations) sampled cycle 0 day 1 and cycle 1 day 21, predose and up to 96 h post-dose. 1340 parent regorafenib observations in total. Absorption (ka, its IIV, and the number of transit compartments) was determined from single-dose study 15823 data alone and then fixed; the parent structural parameters were fixed to the parent-only phase 1 model (Table 1) when the metabolite data were added. Baseline demographics for the phase 1 subset are not tabulated separately in the paper; Appendix Table A2 reports covariate distributions for the phase 3 studies only."
  )

  ini({
    # ---------------------------------------------------------------------
    # Absorption. Two transit compartments for the parent chain and three for
    # the pre-systemically formed M-2 chain (Figure 1), all sharing the single
    # rate constant ka. ka and its IIV were estimated from single-dose study
    # 15823 data and then FIXED for every subsequent step (section 2.3.1).
    # ---------------------------------------------------------------------
    lka <- fixed(log(0.482)); label("First-order absorption / transit rate constant (1/h)")  # Keunecke 2020 Table 2: ka = 0.482 1/h, Fixed

    # ---------------------------------------------------------------------
    # Pre-systemic fraction of dose forming M-2, carried on the LOGIT scale
    # exactly as the authors report it. Table 2 footnote: "FRM2 and FRM5 on
    # logit scale corresponds with 41% and 25%, respectively"; inverse-logit
    # of -0.355 is 0.4121 and of -1.09 is 0.2516, which reproduces that note.
    # ---------------------------------------------------------------------
    logitfm_m2 <- -0.355; label("Logit of the pre-systemic fraction of dose forming M-2 (unitless)")  # Keunecke 2020 Table 2: FRM2 = -0.355, SE 0.0946, RSE 26.7%, 95% CI -0.540 to -0.169
    logitfm_m5 <- -1.09;  label("Logit of the fraction of parent clearance forming M-5 (unitless)")   # Keunecke 2020 Table 2: FRM5 = -1.09, SE 0.212, RSE 19.4%, 95% CI -1.51 to -0.679

    # ---------------------------------------------------------------------
    # Parent disposition. Every value is APPARENT: the paper tabulates
    # CL_P/(1-FRM2), VC_P/(1-FRM2), Q_P/(1-FRM2) and VP_P/(1-FRM2), each
    # additionally relative to oral bioavailability F_oral (Table 2 footnote:
    # "Volumes, clearances and intercompartmental flow are relative to oral
    # bioavailability (F_oral), which is not listed for clarity"). They were
    # fixed here to the values obtained from the parent-only phase 1 fit
    # (Table 1) when the metabolite data were added (section 2.3.2).
    # ---------------------------------------------------------------------
    lcl <- fixed(log(4.02)); label("Apparent parent clearance CL_P/(1-FRM2)/F (L/h)")                       # Keunecke 2020 Table 2: CL_P/(1-FRM2) = 4.02 L/h, Fixed (estimated in Table 1: SE 0.266, RSE 6.60%)
    lvc <- fixed(log(10.7)); label("Apparent parent central volume VC_P/(1-FRM2)/F (L)")                    # Keunecke 2020 Table 2: VC_P/(1-FRM2) = 10.7 L, Fixed (estimated in Table 1: SE 1.92, RSE 18.1%)
    lq  <- fixed(log(11.0)); label("Apparent parent inter-compartmental clearance Q_P/(1-FRM2)/F (L/h)")    # Keunecke 2020 Table 2: Q_P/(1-FRM2) = 11.0 L/h, Fixed (estimated in Table 1: SE 1.37, RSE 12.5%)
    lvp <- fixed(log(162));  label("Apparent parent peripheral volume VP_P/(1-FRM2)/F (L)")                 # Keunecke 2020 Table 2: VP_P/(1-FRM2) = 162 L, Fixed (estimated in Table 1: SE 26.2, RSE 16.1%)

    # ---------------------------------------------------------------------
    # Enterohepatic circulation (Figure 1 and section 2.3.1). Transfer from
    # the central compartment into the gallbladder is first-order and
    # independent of food; the reverse release is gated by food intake.
    # kbm is the paper's k_CG and kehc is the paper's k_GE, fixed to 100 1/h
    # to make emptying "complete and (almost) instantaneous".
    # ---------------------------------------------------------------------
    lkbm  <- fixed(log(0.141)); label("Central-to-gallbladder transfer rate constant k_CG (1/h)")   # Keunecke 2020 Table 2: k_CG = 0.141 1/h, Fixed (estimated in Table 1: SE 0.0190, RSE 13.5%)
    lkehc <- fixed(log(100));   label("Gallbladder emptying rate constant k_GE (1/h)")              # Keunecke 2020 Table 2: k_GE = 100 1/h, Fixed; section 2.3.1 "complete and (almost) instantaneous emptying (assumed first-order rate constant 100 h-1)"

    # Meal clock times (hours since midnight) and the duration of each
    # gallbladder-emptying window. Three releases per day (section 2.3.1).
    tmeal1 <- fixed(8);   label("Clock time of breakfast, the first gallbladder-emptying trigger (h since midnight)")  # Keunecke 2020 Table 2: Breakfast = 08:00, Fixed
    tmeal2 <- fixed(12);  label("Clock time of lunch, the second gallbladder-emptying trigger (h since midnight)")     # Keunecke 2020 Table 2: Lunch = 12:00, Fixed
    tmeal3 <- fixed(18);  label("Clock time of dinner, the third gallbladder-emptying trigger (h since midnight)")     # Keunecke 2020 Table 2: Dinner = 18:00, Fixed
    dge    <- fixed(0.5); label("Duration of gallbladder emptying after each meal (h)")                                # Keunecke 2020 Table 2: DGE = 0.500 h, Fixed; section 2.3.1 "the duration of transfer ... was approximately 30 min"

    # ---------------------------------------------------------------------
    # Metabolite clearances. Apparent, and the only metabolite disposition
    # parameters that could be estimated: "Since no PK data are available
    # after direct administration of M-2 and M-5, the volume of distribution
    # of M-2 and M-5 cannot be estimated, and it is assumed that their volume
    # is the same as that of parent regorafenib" (section 2.3.2). Volumes,
    # inter-compartmental flow and k_CG for both metabolites are therefore
    # NOT separate parameters -- model() reuses the parent's vc, vp, q, kbm.
    # ---------------------------------------------------------------------
    lcl_m2 <- log(2.45);  label("Apparent M-2 clearance CL_M-2/F (L/h)")  # Keunecke 2020 Table 2: CL_M-2 = 2.45 L/h, SE 0.122, RSE 4.98%, 95% CI 2.21 to 2.69
    lcl_m5 <- log(0.746); label("Apparent M-5 clearance CL_M-5/F (L/h)")  # Keunecke 2020 Table 2: CL_M-5 = 0.746 L/h, SE 0.215, RSE 28.8%, 95% CI 0.325 to 1.17

    # ---------------------------------------------------------------------
    # Inter-individual variability. The appendix states the IIV model is
    # exponential (theta_i = TV_theta * exp(eta_i)), so the tabulated
    # omega^2 values are variances on the log scale and go into ini()
    # unchanged. FRM2 / FRM5 IIV is on the logit scale.
    # ---------------------------------------------------------------------
    # Keunecke 2020 Table 2, row "omega^2 k a" = 0.127, held constant (carried
    # from the single-dose absorption step, section 2.3.1). Eta traces are kept
    # on their own line above the declaration: rxode2 rewrites a trailing
    # comment on an eta line into a label(), which then duplicates fixed().
    etalka ~ fixed(0.127)

    # Full 3x3 clearance block. The published matrix is very nearly singular
    # (determinant 9.2e-05; corr(CL_M-2, CL_M-5) = 0.91), which is the
    # numerical difficulty the authors describe: the CL_P / CL_M-5 covariance
    # is the one parameter with RSE > 50% and a 95% CI spanning zero, and
    # section 3.3.2 records that the covariance step failed for the
    # corresponding phase 3 base model. It is still positive definite, so it
    # is carried here exactly as published; the phase 3 model resolves the
    # problem by sharing one eta between the two metabolites.
    #
    # Source trace for the block below, kept on its own lines above the
    # declaration: a trailing comment INSIDE a multi-line omega c(...) is
    # rewritten into a bare `;` and makes the model unparseable.
    #   Keunecke 2020 Table 2, row "omega^2 CL_P"        = 0.117 (SE 0.0231, RSE 19.7%)
    #   Keunecke 2020 Table 2, row "omega CL_P,CL_M-2"   = 0.122 (RSE 21.4%)
    #   Keunecke 2020 Table 2, row "omega^2 CL_M-2"      = 0.267 (RSE 18.9%)
    #   Keunecke 2020 Table 2, row "omega CL_P,CL_M-5"   = 0.157 (RSE 52.4%, 95% CI -0.00412 to 0.319)
    #   Keunecke 2020 Table 2, row "omega CL_M-2,CL_M-5" = 0.656 (RSE 18.4%)
    #   Keunecke 2020 Table 2, row "omega^2 CL_M-5"      = 1.95 (RSE 24.0%)
    etalcl + etalcl_m2 + etalcl_m5 ~ c(
      0.117,
      0.122, 0.267,
      0.157, 0.656, 1.95
    )

    # Keunecke 2020 Table 2, row "omega^2 FRM2" = 0.156 (SE 0.0759, RSE 48.7%)
    etalogitfm_m2 ~ 0.156
    # Keunecke 2020 Table 2, row "omega^2 FRM5" = 0.841 (SE 0.258, RSE 30.7%)
    etalogitfm_m5 ~ 0.841

    # ---------------------------------------------------------------------
    # Residual error. Proportional for all three analytes; the metabolites
    # additionally carry an additive term fixed at 0.001 mg/L (1 ug/L, just
    # below the 2.0 ug/L LLOQ) that acts as a numerical stabiliser.
    # ---------------------------------------------------------------------
    propSd    <- 0.406;         label("Parent proportional residual error (fraction)")  # Keunecke 2020 Table 2: Parent Prop. error SD = 0.406, SE 0.0208, RSE 5.11%
    addSd_m2  <- fixed(0.001);  label("M-2 additive residual error (mg/L)")             # Keunecke 2020 Table 2: M-2 Add. error SD = 0.001, Fixed
    propSd_m2 <- 0.380;         label("M-2 proportional residual error (fraction)")     # Keunecke 2020 Table 2: M-2 Prop. error SD = 0.380, SE 0.0231, RSE 6.07%
    addSd_m5  <- fixed(0.001);  label("M-5 additive residual error (mg/L)")             # Keunecke 2020 Table 2: M-5 Add. error SD = 0.001, Fixed
    propSd_m5 <- 0.455;         label("M-5 proportional residual error (fraction)")     # Keunecke 2020 Table 2: M-5 Prop. error SD = 0.455, SE 0.0256, RSE 5.63%
  })

  model({
    # 1. Individual parameters.
    ka    <- exp(lka + etalka)
    cl    <- exp(lcl + etalcl)
    vc    <- exp(lvc)
    q     <- exp(lq)
    vp    <- exp(lvp)
    kbm   <- exp(lkbm)
    kehc  <- exp(lkehc)

    cl_m2 <- exp(lcl_m2 + etalcl_m2)
    cl_m5 <- exp(lcl_m5 + etalcl_m5)

    # Formation fractions, inverse-logit back onto (0, 1). The individual
    # logit-scale value is collected on its own line first so that the eta
    # stays mu-referenced.
    logitfm_m2_ind <- logitfm_m2 + etalogitfm_m2
    logitfm_m5_ind <- logitfm_m5 + etalogitfm_m5
    fm_m2 <- 1 / (1 + exp(-logitfm_m2_ind))
    fm_m5 <- 1 / (1 + exp(-logitfm_m5_ind))

    # 2. Food-intake switch. Gallbladder emptying is triggered by each of the
    #    three daily meals and lasts dge hours (Figure 1: "Not during meal" /
    #    "During meal"). `t` is absolute time, so the model expects an event
    #    table whose origin is midnight; tod is the resulting time of day.
    tod   <- t - 24 * floor(t / 24)
    meal  <- (tod >= tmeal1) * (tod < tmeal1 + dge) +
             (tod >= tmeal2) * (tod < tmeal2 + dge) +
             (tod >= tmeal3) * (tod < tmeal3 + dge)

    # 3. Parent regorafenib: depot -> 2 transit compartments -> central,
    #    two-compartment disposition, plus the gallbladder cycle. This is the
    #    Figure 1 caption ODE system for A3 (central), A4 (gallbladder) and
    #    A5 (peripheral), with the meal switch selecting between the
    #    "Not during meal" and "During meal" forms.
    d/dt(depot)       <- -ka * depot
    d/dt(transit1)    <-  ka * (1 - fm_m2) * depot - ka * transit1
    d/dt(transit2)    <-  ka * transit1 - ka * transit2
    d/dt(central)     <-  ka * transit2 -
                            kbm * central -
                            (cl / vc) * central -
                            (q / vc) * central + (q / vp) * peripheral1 +
                            meal * kehc * gallbladder
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1
    d/dt(gallbladder) <-  kbm * central - meal * kehc * gallbladder

    # 4. M-2, formed pre-systemically (section 2.3.2), through its own
    #    three-transit absorption chain (Figure 1, left branch). Its volumes,
    #    inter-compartmental flow and gallbladder rate constant are the
    #    parent's, per the paper's identifiability assumption.
    d/dt(transit1_m2)    <-  ka * fm_m2 * depot - ka * transit1_m2
    d/dt(transit2_m2)    <-  ka * transit1_m2 - ka * transit2_m2
    d/dt(transit3_m2)    <-  ka * transit2_m2 - ka * transit3_m2
    d/dt(central_m2)     <-  ka * transit3_m2 -
                               kbm * central_m2 -
                               (cl_m2 / vc) * central_m2 -
                               (q / vc) * central_m2 + (q / vp) * peripheral1_m2 +
                               meal * kehc * gallbladder_m2
    d/dt(peripheral1_m2) <-  (q / vc) * central_m2 - (q / vp) * peripheral1_m2
    d/dt(gallbladder_m2) <-  kbm * central_m2 - meal * kehc * gallbladder_m2

    # 5. M-5, formed systemically from the parent central compartment as the
    #    fraction fm_m5 of parent clearance (Figure 1: the parent's CL splits
    #    into CL_Parent * FRM5 towards M-5 and CL_Parent * (1 - FRM5) out).
    d/dt(central_m5)     <-  fm_m5 * (cl / vc) * central -
                               kbm * central_m5 -
                               (cl_m5 / vc) * central_m5 -
                               (q / vc) * central_m5 + (q / vp) * peripheral1_m5 +
                               meal * kehc * gallbladder_m5
    d/dt(peripheral1_m5) <-  (q / vc) * central_m5 - (q / vp) * peripheral1_m5
    d/dt(gallbladder_m5) <-  kbm * central_m5 - meal * kehc * gallbladder_m5

    # 6. Bioavailability. The tabulated parent parameters are already divided
    #    by (1 - FRM2), i.e. they describe a parent chain that receives the
    #    FULL nominal dose. Scaling the depot by 1 / (1 - fm_m2) before the
    #    Figure 1 split restores that convention, and simultaneously puts the
    #    M-2 branch on the same apparent scale so that the metabolite volumes
    #    really are the parent's tabulated 10.7 L and 162 L.
    f(depot) <- 1 / (1 - fm_m2)

    # 7. Observations. All three concentrations use the parent's apparent
    #    central volume (metabolite volumes fixed to the parent values).
    Cc    <- central / vc
    Cc_m2 <- central_m2 / vc
    Cc_m5 <- central_m5 / vc

    Cc    ~ prop(propSd)
    Cc_m2 ~ add(addSd_m2) + prop(propSd_m2)
    Cc_m5 ~ add(addSd_m5) + prop(propSd_m5)
  })
}
