Palmer_2025_moxifloxacin <- function() {
  description <- "Two-compartment population PK model with a four-step first-order transit-chain absorption for oral moxifloxacin in children with rifampicin-resistant tuberculosis (Palmer 2025, CATALYST trial). Apparent disposition parameters are allometrically scaled on body weight to a 16 kg reference child with fixed exponents 0.75 (CL/F, Q/F) and 1 (Vc/F, Vp/F); bioavailability is fixed to 1 because only oral data were available. Interindividual variability is carried on CL/F alone, while interoccasion variability sits on bioavailability and on the mean absorption time MAT, with the IOV magnitude inflated 1.70-fold on dosing occasions that were not directly observed. A dispersible 100 mg paediatric tablet and the standard crushed 400 mg adult tablet were found bioequivalent, so the final model carries no formulation term."
  reference <- paste(
    "Palmer M, Zou Y, Hesseling AC, van der Laan L, Courtney I, Kinikar AA,",
    "Sonkawade N, Paradkar M, Kulkarni V, Casalme DJO, Frias MVG, Draper H,",
    "Wiesner L, Karlsson MO, Denti P, Svensson EM, Garcia-Prats AJ.",
    "Population pharmacokinetics and dosing of dispersible moxifloxacin",
    "formulation in children with rifampicin-resistant tuberculosis.",
    "Br J Clin Pharmacol. 2025;91(6):1853-1864. doi:10.1002/bcp.70005.",
    "Structural model adapted from Radtke KK et al.,",
    "Clin Infect Dis. 2022;74(8):1372-1381.",
    sep = " "
  )
  vignette <- "Palmer_2025_moxifloxacin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "The only covariate retained in the final model. Allometric scaling on all four apparent disposition parameters relative to a 16 kg reference child, with fixed exponents 0.75 on CL/F and Q/F and 1 on Vc/F and Vp/F (Table 2 footnote a). Body weight was selected in preference to fat-free mass as the body-size descriptor (Results paragraph 2). Observed range 6.9-42.1 kg; the paper's own dosing simulations extend to 3-31 kg and it warns that weights below 7 kg are extrapolated beyond the CATALYST data.",
      source_name = "WT"
    ),
    OCC = list(
      description = "Integer-valued occasion index used to multiplex the interoccasion-variability etas on bioavailability and on mean absorption time",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "The CATALYST design gives four occasions per child: the predose sample at the PK1 visit arises from the previous day's at-home dose (OCC = 1), the directly observed standard-formulation dose given at the PK1 visit (OCC = 2), the predose sample at the PK2 visit arising from the previous day's at-home dose (OCC = 3), and the directly observed dispersible-formulation dose given at the PK2 visit (OCC = 4). Methods 2.2: 'Predose samples resulting from prior doses were handled as a separate occasion from the actual postdose samples.' The source NONMEM control stream (Supplementary Material 9) writes a single shared IOV variance per parameter across occasions via an IF (OCC.EQ.N) block, so occasions 2-4 are encoded here as fixed etas equal to the estimated occasion-1 variance, mirroring the NONMEM $OMEGA BLOCK(1) SAME idiom that nlmixr2 has no shortcut for. Decomposed inside model() into the binary indicators oc1..oc4.",
      source_name = "OCC"
    ),
    SELFADMIN = list(
      description = "Binary indicator that the dosing occasion was not directly observed, i.e. the dose was given at home by the caregiver rather than under supervision at the study visit",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (directly observed dose administered at the PK visit)",
      notes = "Source column UNOBS, defined in the Supplementary Material 9 NONMEM control stream as the '; flag for PK samples from unobserved dosing events' driving UNOBS_EFF = 1 + THETA(10) * UNOBS. Polarity matches the SELFADMIN canonical exactly (1 = unobserved / self-administered, 0 = directly observed), so no inversion is required. Unlike the founding Wallender 2021 usage, the effect here is not on the typical bioavailability but on the MAGNITUDE of the interoccasion variability: the control stream writes F1 = 1 * EXP(IOVF * UNOBS_EFF) and MAT = TVMAT * EXP(IOVMAT * UNOBS_EFF), so an unobserved occasion multiplies the IOV standard deviation on both parameters by the estimated 1.70-fold factor. The paper attributes this to caregivers reporting imprecise dosing times for at-home doses (Discussion paragraph 5). In the CATALYST design SELFADMIN = 1 exactly on OCC 1 and 3.",
      source_name = "UNOBS"
    )
  )

  covariatesDataExcluded <- list(
    HAZ = list(
      description = "Height-for-age z-score, UK-WHO growth chart",
      units = "(z-score)",
      type = "continuous",
      notes = "Carried over from the Radtke model as a clearance covariate and re-tested in CATALYST, but not retained: dOFV < 0.01 and the estimated effect was -0.96% (95% CI -5.8 to 3.9)% per z-score unit, against Radtke's -9.8 (95% CI -16 to -3.8)%. Discussion paragraph 6 attributes the discrepancy to a different growth reference and to Indian and Philippine children having lower average height-for-age than the international standard."
    ),
    WAZ = list(
      description = "Weight-for-age z-score, UK-WHO growth chart",
      units = "(z-score)",
      type = "continuous",
      notes = "Explored as a malnutrition-status covariate on clearance (Methods 2.2) and not retained in the final model; no point estimate is reported."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested on clearance as both a sigmoidal maturation function (dOFV < 0.01) and a proportional effect (dOFV = -1.6, estimated at +8% in children under 2 years); neither was significant and neither was retained. A published renal-maturation function added as a fixed component worsened the OFV by 13 points (Discussion paragraph 5)."
    ),
    FFM = list(
      description = "Fat-free mass, calculated with a paediatric formula",
      units = "kg",
      type = "continuous",
      notes = "Tested against body weight as the allometric body-size descriptor (Methods 2.2). Body weight was selected; no FFM-based parameter estimates are reported."
    ),
    HIV_POS = list(
      description = "HIV-positive status",
      units = "(binary)",
      type = "binary",
      notes = "A clearance covariate in the Radtke model. Only 1 of the 36 CATALYST children was living with HIV, so the effect could not be meaningfully estimated, and fixing it to the previously reported value worsened the model fit; it was therefore not included (Results paragraph 2)."
    ),
    FORM_DISPERSIBLE = list(
      description = "Dispersible 100 mg scored paediatric tablet versus the crushed standard 400 mg adult tablet",
      units = "(binary)",
      type = "binary",
      notes = "The whole point of the trial, and formally estimated stepwise on bioavailability, on MAT, and on their variabilities -- but not retained, because none of the differences was significant. The bioavailability ratio (dispersible vs standard) was 105% with a log-likelihood-profiling 90% CI of 95-115%, inside the 80-125% bioequivalence limits, and the MAT ratio was 106% (90% CI 90-126%). The final model in Table 2 therefore has no formulation term, and the paper's conclusion is that dosing recommendations can be identical for the two formulations. FORM_DISPERSIBLE is not a ratified canonical -- the name is recorded here for documentation only and no covariate column is required to use this model."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "moxifloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "moxifloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "moxifloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "moxifloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "moxifloxacin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "moxifloxacin", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 36L,
    n_studies = 1L,
    n_observations = 384L,
    age_range = "Median 4.8 years (range 0.4-15); 4 children were under 1 year and 1 was under 6 months at the pharmacokinetic visit (Table 1; Discussion paragraph 5).",
    weight_range = "Median 15.6 kg (range 6.9-42.1); enrolment was stratified into two parallel weight cohorts, 16 participants under 15 kg of whom 6 were under 10 kg (Table 1; Results paragraph 1).",
    sex_female_pct = 61,
    disease_state = "Children being treated for rifampicin-resistant tuberculosis in routine TB programmes, on a regimen containing both clofazimine and a fluoroquinolone, within 16 weeks of treatment start. TB was microbiologically confirmed in 17 (47.2%) and unconfirmed in 19 (52.8%). Median height-for-age z-score -1.1 (range -4.1 to 0.46) and weight-for-age z-score -1.6 (range -4.7 to 1.2). One participant (3%) was living with HIV.",
    dose_range = "Daily oral moxifloxacin at the WHO weight-band doses (WHO 2022 guideline Table 4), 80-400 mg once daily. Each child was sampled twice: at visit PK1 on the standard non-dispersible 400 mg tablet crushed and suspended in water, and at visit PK2, 1-14 days later, on the first dose of the dispersible 100 mg scored tablet dispersed in water.",
    regions = "South Africa (n = 20), the Philippines (n = 10), India (n = 6).",
    notes = "CATALYST, an open-label multisite trial (Pan African Clinical Trials Registry 202012756409365). Intensive sampling at 0, 1, 2, 4, 8 and 24 h postdose on each of the two visit days, with children fasted at least 4 h before the predose draw and then fed breakfast before dosing. Moxifloxacin was assayed by validated LC-MS/MS with an LLOQ of 0.0628 ug/mL. The paper is internally inconsistent about the observation count: Results paragraph 1 states that 438 observations were obtained of which 4 predose samples were below the quantification limit and excluded, while the Table 1 header gives obs = 384 for the analysis dataset. The Table 1 figure is recorded here; see the validation vignette Errata. Fitted in NONMEM 7.5 with FOCE-I; parameter precision by sampling importance resampling, and the formulation-ratio confidence intervals by log-likelihood profiling. sex_female_pct is the complement of the 14/36 (39%) male count in Table 1."
  )

  ini({
    # ---- Structural fixed effects (Palmer 2025 Table 2, 'Estimate' column) ----
    # All disposition parameters are APPARENT (X/F): Table 2 footnote a records
    # that F is inestimable from oral-only data and was fixed at 1, and that the
    # typical values are stated for a 16 kg reference child.
    lcl  <- log(6.90);  label("Apparent oral clearance CL/F at 16 kg (L/h)")                     # Table 2: CL/F = 6.90 L/h (RSE 3.5%, 95% CI 6.44-7.38)
    lvc  <- log(61.1);  label("Apparent central volume of distribution Vc/F at 16 kg (L)")       # Table 2: Vc/F = 61.1 L (RSE 3.7%, 95% CI 57.1-65.6)
    lq   <- log(0.860); label("Apparent intercompartmental clearance Q/F at 16 kg (L/h)")        # Table 2: Q/F = 0.860 L/h (RSE 19%, 95% CI 0.542-1.20)
    lvp  <- log(44.5);  label("Apparent peripheral volume of distribution Vp/F at 16 kg (L)")    # Table 2: Vp/F = 44.5 L (RSE 36%, 95% CI 18.7-81.8)
    lmat <- log(1.01);  label("Mean absorption time MAT across the four transit steps (h)")      # Table 2: MAT = 1.01 h (RSE 5.5%, 95% CI 0.909-1.13)
    lfdepot <- fixed(log(1)); label("Typical oral bioavailability F (unitless)")                 # Table 2: F = 1 (fixed); footnote a - not estimable with only oral data

    # ---- Allometric exponents (Palmer 2025 Table 2 footnote a, both fixed) ----
    # "Allometric scaling with body weight (WT) is applied following: (WT/16)^theta,
    #  with a fixed exponent (theta) of 0.75 for CL and Q, and 1 for Vc and Vp."
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL/F and Q/F (unitless)")              # Table 2 footnote a: exponent fixed at 0.75 for the clearances
    e_wt_vc <- fixed(1.0);  label("Allometric exponent on Vc/F and Vp/F (unitless)")             # Table 2 footnote a: exponent fixed at 1 for the volumes

    # ---- Unobserved-dosing inflation of the IOV magnitude ----
    # Supplementary Material 9 $PK: UNOBS_EFF = 1 + THETA(10) * UNOBS, and both
    # F1 = 1 * EXP(IOVF * UNOBS_EFF) and MAT = TVMAT * EXP(IOVMAT * UNOBS_EFF).
    # The multiplier therefore scales the IOV STANDARD DEVIATION, not a typical
    # value. Table 2 prints the multiplier itself (1 + THETA(10)), so the value
    # transcribed here is the printed fold change and it is applied as
    # e_selfadmin_iovscale^SELFADMIN, which is 1 at a directly observed occasion.
    e_selfadmin_iovscale <- 1.70; label("Fold change in the IOV standard deviation on F and on MAT at an unobserved dosing occasion (unitless)") # Table 2: fold change of unobserved dosing event on IOV of F and MAT = 1.70 (RSE 26%, 95% CI 0.970-2.68)

    # ---- Interindividual variability (Palmer 2025 Table 2 'IIV/IOV' block) ----
    # Table 2 footnote c defines the reported %CV as "the square root of estimated
    # variances * 100%", i.e. the printed number is 100 * omega where omega is the
    # standard deviation on the log scale. The nlmixr2 eta slot holds the VARIANCE,
    # so omega^2 = (CV / 100)^2 directly -- this is NOT the log(1 + CV^2) transform
    # that applies when a paper reports sqrt(exp(omega^2) - 1).
    #
    # IIV sits on clearance only; Results paragraph 2: "The model included IIV for
    # clearance and IOV for bioavailability and MAT."
    etalcl ~ 0.018496   # Table 2: IIV on CL = 13.6 %CV -> omega = 0.136, omega^2 = 0.136^2 = 0.018496

    # ---- Interoccasion variability (Palmer 2025 Table 2 'IIV/IOV' block) ----
    # One shared variance per parameter across the four CATALYST occasions
    # (Supplementary Material 9 IF (OCC.EQ.N) block). Occasion 1 carries the
    # estimated variance and occasions 2-4 are fixed to the same value, which is
    # how nlmixr2 expresses NONMEM's $OMEGA BLOCK(1) SAME.
    etaiov_fdepot_1 ~ 0.094864          # Table 2: IOV on F = 30.8 %CV -> omega = 0.308, omega^2 = 0.094864
    etaiov_fdepot_2 ~ fixed(0.094864)   # shared IOV variance across occasions
    etaiov_fdepot_3 ~ fixed(0.094864)
    etaiov_fdepot_4 ~ fixed(0.094864)
    etaiov_mat_1    ~ 0.178084          # Table 2: IOV on MAT = 42.2 %CV -> omega = 0.422, omega^2 = 0.178084
    etaiov_mat_2    ~ fixed(0.178084)   # shared IOV variance across occasions
    etaiov_mat_3    ~ fixed(0.178084)
    etaiov_mat_4    ~ fixed(0.178084)

    # ---- Residual unexplained variability (Palmer 2025 Table 2) ----
    # Supplementary Material 9 $ERROR: Y = IPRED + IPRED * EPS(1) + EPS(2),
    # i.e. combined proportional plus additive on the linear concentration scale.
    propSd <- 0.136;  label("Proportional residual error (fraction)")   # Table 2: proportional error = 13.6% (RSE 6.2%); footnote c -> SD scale
    addSd  <- 0.0435; label("Additive residual error (mg/L)")           # Table 2: additive error = 0.0435 mg/L (RSE 13%); footnote d - "reported in the standard deviation scale"
  })

  model({
    # ---- Occasion decomposition for the IOV etas ----
    # OCC = 1 PK1 predose (previous day's at-home dose), 2 PK1 observed dose,
    #       3 PK2 predose (previous day's at-home dose), 4 PK2 observed dose.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)

    # An unobserved (home-administered) dosing occasion inflates the IOV standard
    # deviation on BOTH F and MAT by e_selfadmin_iovscale. At SELFADMIN = 0 the
    # multiplier is exactly 1, so a directly observed occasion is unmodified.
    iov_scale <- e_selfadmin_iovscale^SELFADMIN

    iov_fdepot <- (oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 +
                   oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4) * iov_scale
    iov_mat    <- (oc1 * etaiov_mat_1    + oc2 * etaiov_mat_2 +
                   oc3 * etaiov_mat_3    + oc4 * etaiov_mat_4)    * iov_scale

    # ---- Individual PK parameters, allometrically scaled to a 16 kg child ----
    cl  <- exp(lcl + etalcl) * (WT / 16)^e_wt_cl
    vc  <- exp(lvc)          * (WT / 16)^e_wt_vc
    q   <- exp(lq)           * (WT / 16)^e_wt_cl
    vp  <- exp(lvp)          * (WT / 16)^e_wt_vc
    mat <- exp(lmat + iov_mat)

    # Transit-chain rate constant. Supplementary Material 9 $PK sets KTR = 4/MAT
    # and then K14 = K45 = K56 = K62 = KTR, i.e. the dose compartment and the
    # three transit compartments each empty at rate ktr into the next, with the
    # last discharging directly into central. Four sequential first-order steps
    # at rate ktr give a total mean absorption time of 4/ktr = MAT.
    ktr <- 4 / mat

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---- ODE system ----
    # depot -> transit1 -> transit2 -> transit3 -> central <-> peripheral1
    # $MODEL NCOMP=6 COMP=(DOS,DEFDOSE) COMP=(MOXC,DEFOBS) COMP=(MOXP)
    #                COMP=(TRANSI) COMP=(TRANSI2) COMP=(TRANSI3)
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2 - ktr * transit3
    d/dt(central)     <-  ktr * transit3 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                   k12 * central - k21 * peripheral1

    # Bioavailability applies to the dose compartment (NONMEM F1 on COMP 1 = DOS).
    # The typical value is fixed at 1 and carries interoccasion variability only.
    f(depot) <- exp(lfdepot + iov_fdepot)

    # Dose in mg, Vc/F in L -> mg/L, the units of the additive residual error and
    # of the paper's reported concentrations.
    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
