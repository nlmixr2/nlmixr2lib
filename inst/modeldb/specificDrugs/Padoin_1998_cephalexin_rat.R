Padoin_1998_cephalexin_rat <- function() {
  description <- "Preclinical (rat, male Wistar). Two-compartment population PK model for cephalexin after intra-arterial (IA) or oral (gastric-tube) administration in rats, with first-order absorption and a competitive drug-drug interaction from coadministered oral quinapril that lowers cephalexin Ka (paper Table 4: 0.249 to 0.177 1/h; ~29% lower) and CL (paper Table 4: 0.810 to 0.640 L/h/kg; ~21% lower) when both drugs are given by the oral route. The paper parameterizes the disposition as {CL, Vc, CL_D, Vss = Vc + Vp}; this implementation uses the canonical {CL, Vc, Q, Vp} parameterization with the typical-value Vp derived as Vss - Vc = 1.23 - 0.416 = 0.814 L/kg. Intra-arterial quinapril or intra-arterial cephalexin produced no detectable interaction on cephalexin elimination, attributed by the authors to the much higher cephalexin renal concentrations outcompeting quinapril at the renal anionic transport carrier (and to Ka being irrelevant for IA dosing)."
  reference <- paste(
    "Padoin C, Tod M, Perret G, Petitjean O.",
    "Analysis of the pharmacokinetic interaction between cephalexin",
    "and quinapril by a nonlinear mixed-effect model.",
    "Antimicrob Agents Chemother. 1998;42(6):1463-1469.",
    "doi:10.1128/aac.42.6.1463.",
    sep = " "
  )
  vignette <- "Padoin_1998_cephalexin_rat"
  units <- list(
    time          = "hour",
    dosing        = "mg/kg",
    concentration = "mg/L"
  )

  covariateData <- list(
    CONMED_QPRL_ORAL = list(
      description        = "Indicator for concomitant oral (gastric-tube) quinapril coadministration: 1 = subject received oral quinapril (Padoin 1998 group 5: cephalexin GT + quinapril GT), 0 = otherwise (Padoin 1998 groups 1 (cephalexin IA only), 2 (cephalexin IA + quinapril IA), 3 (cephalexin IA + quinapril GT), 4 (cephalexin GT only)).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no oral quinapril coadministration; reference pools groups 1-4 because the paper found no DDI on cephalexin CL or Ka in any of these conditions)",
      notes              = "Padoin 1998 dose / regimen: single 0.8 mg/kg oral dose of quinapril via gastric tube, administered 15 min before the 50 mg/kg oral cephalexin dose (also via gastric tube), in male Wistar rats (n = 8 per group). Time-fixed per subject (parallel-group design). The paper's final model specifies CL = 0.640 L/h/kg and Ka = 0.177 1/h when CONMED_QPRL_ORAL = 1 vs the reference CL = 0.810 L/h/kg and Ka = 0.249 1/h when CONMED_QPRL_ORAL = 0. The reference category pools group 3 (cephalexin IA + quinapril GT) into the no-DDI condition because the paper specifies CL_j = CL1 = 0.810 L/h/kg for groups 1-4 and only group 5 receives the lower CL; Ka is irrelevant in groups 1-3 (IA dosing skips the depot compartment).",
      source_name        = "(derived from the paper's group indicator; the paper's final model writes 'Ka_j = Ka1 if quinapril is not given; Ka_j = Ka2 if quinapril is given via a GT' and 'CL_j = CL1 for groups 1 to 4; CL_j = CL2 for group 5', collapsing both effects into the single binary CONMED_QPRL_ORAL canonical column for downstream simulation)"
    )
  )

  population <- list(
    species        = "rat (male Wistar, IFACREDO, L'Arbresle, France)",
    n_subjects     = 40L,
    n_studies      = 1L,
    age_range      = "not reported (animals fasted 18 h before the experiment, water freely available)",
    weight_range   = "250-280 g body weight at experiment start",
    sex_female_pct = 0,
    race_ethnicity = NA,
    disease_state  = "healthy rats prepared with an indwelling carotid artery catheter (installed 24 h pre-experiment under thiopental anesthesia, 50 mg/kg IP) for intra-arterial drug administration and serial arterial blood sampling; oral doses delivered via gastric tube as a 2% methylcellulose suspension",
    dose_range     = "single 50 mg/kg dose of cephalexin (IA or via gastric tube) +/- single 0.8 mg/kg dose of quinapril (IA or via gastric tube, administered 15 min before cephalexin when given orally)",
    regions        = "preclinical (in-vivo rat); Hopital Avicenne, Bobigny, France",
    notes          = "Five parallel groups of n = 8 rats each (total n = 40): group 1 cephalexin IA only; group 2 cephalexin IA + quinapril IA; group 3 cephalexin IA + quinapril via GT; group 4 cephalexin via GT only; group 5 cephalexin via GT + quinapril via GT (the only group with the demonstrated DDI). Arterial blood samples (0.15 mL) collected at 0, 5, 15, 30, 45, 60, 90 min and 2, 3, 4, 5, 6 h after cephalexin administration; blood volume replaced with twice-volume isotonic saline after 30 min. Cephalexin assayed by HPLC with UV detection at 262 nm (Spherisorb C18 column); calibration linear 2-100 mg/L, LLOQ 2.0 mg/L, interassay precision 7-10% CV (paper Methods, Analytical methods). Cephalexin plasma protein binding determined ex vivo (n = 5 separate rats) was constant at fu = 0.82 +/- 0.08 at 5, 30, and 120 min post-dose. Population fit by NONMEM IV.2.0 with the first-order conditional estimation method (METHOD=COND); both groups 1-4 and group 5 were fit jointly to give the final model in paper Table 4."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural PK parameters. Paper Table 4 reports the final-model
    # population means for the reduced (no-DDI) condition (Ka1, CL1) and
    # the oral-quinapril DDI condition (Ka2, CL2). The reference values
    # below correspond to CONMED_QPRL_ORAL = 0 (groups 1-4); the
    # multiplicative covariate effects (below) recover the CL2 and Ka2
    # values when CONMED_QPRL_ORAL = 1.
    lka     <- log(0.249);   label("Cephalexin absorption rate constant Ka1 (1/h) - reference (no oral quinapril)")                   # Padoin 1998 Table 4: Ka_1 = 0.249 1/h, SE 0.056
    lcl     <- log(0.810);   label("Cephalexin elimination clearance CL1 (L/h/kg) - reference (groups 1-4)")                         # Padoin 1998 Table 4: CL_1 = 0.810 L/h/kg, SE 0.110
    lvc     <- log(0.416);   label("Central volume of distribution Vc (L/kg)")                                                       # Padoin 1998 Table 4: Vc = 0.416 L/kg, SE 0.057
    lq      <- log(0.363);   label("Distributional clearance Q (= paper's CL_D) between central and peripheral1 (L/h/kg)")           # Padoin 1998 Table 4: CL_D = 0.363 L/h/kg, SE 0.071
    lvp     <- log(0.814);   label("Peripheral1 volume Vp (L/kg); derived as Vss - Vc = 1.23 - 0.416 = 0.814")                       # Derived from Padoin 1998 Table 4: Vss = 1.23 L/kg (SE 0.238) and Vc = 0.416 L/kg (SE 0.057)
    lfdepot <- log(0.89);    label("Bioavailability F of cephalexin after oral (gastric-tube) dosing (fraction)")                     # Padoin 1998 Table 4: F = 0.89, SE 0.25; the same F applies to both oral groups (paper Results: "without modification of the extent of absorption")

    # ------------------------------------------------------------------
    # Covariate effects of oral quinapril coadministration (CONMED_QPRL_ORAL).
    # Encoded as multiplicative exponential on log-scale so that
    # CONMED_QPRL_ORAL = 0 recovers the reference Ka1 / CL1 and
    # CONMED_QPRL_ORAL = 1 recovers Ka2 / CL2 from paper Table 4 exactly.
    e_conmed_qprl_oral_ka <- log(0.177 / 0.249); label("Multiplicative log-effect of CONMED_QPRL_ORAL on Ka (log(0.177/0.249) = -0.3413; ~29% lower Ka under oral quinapril coadministration)") # Padoin 1998 Table 4: Ka_2 = 0.177 1/h vs Ka_1 = 0.249 1/h
    e_conmed_qprl_oral_cl <- log(0.640 / 0.810); label("Multiplicative log-effect of CONMED_QPRL_ORAL on CL (log(0.640/0.810) = -0.2356; ~21% lower CL under oral quinapril coadministration)") # Padoin 1998 Table 4: CL_2 = 0.640 L/h/kg vs CL_1 = 0.810 L/h/kg

    # ------------------------------------------------------------------
    # Interindividual variability. Paper Results states the IIV form is
    # log-normal: Pj = P# * exp(eta_j), eta_j ~ N(0, var). The variances
    # reported in Table 4 are entered directly as the log-scale variance
    # of the corresponding eta. Var(eta_Ka) and Var(eta_F) were not
    # statistically different from zero and were fixed to zero in the
    # final model (paper Results); they are therefore omitted from this
    # block. Var(eta_CL_D) and Var(eta_Vss) carry large SEs (CIs include
    # zero) but the paper retained them because fixing to zero worsened
    # the fit per the likelihood ratio test.
    etalcl ~ 0.382  # Padoin 1998 Table 4: Var(eta_CL) = 0.382, SE 0.109
    etalvc ~ 0.783  # Padoin 1998 Table 4: Var(eta_Vc) = 0.783, SE 0.173
    etalq  ~ 2.34   # Padoin 1998 Table 4: Var(eta_CL_D) = 2.34, SE 0.997
    etalvp ~ 2.38   # Approximation: paper Table 4 reports Var(eta_Vss) = 2.38 (SE 1.12); since Vp dominates Vss (Vp/Vss = 0.814/1.23 = 66%), Var(eta_Vp) is approximated by Var(eta_Vss). See vignette Assumptions and deviations.

    # ------------------------------------------------------------------
    # Residual error. Paper Methods specifies the additive-proportional
    # form Cij = C_hat_ij + epsilon_i * C_hat_ij = C_hat_ij * (1 + epsilon),
    # epsilon ~ N(0, sigma2). Table 4 reports sigma2_e = 0.033 (SE 0.015);
    # the SD is sqrt(0.033) = 0.1817 (~18.2% CV).
    propSd <- sqrt(0.033); label("Proportional residual error SD on cephalexin plasma concentrations (fraction; ~0.1817 = ~18.2% CV)") # Padoin 1998 Table 4: sigma2_e = 0.033, SE 0.015; SD computed as sqrt(0.033)
  })

  model({
    # ---- Individual structural parameters ----
    # Note: paper Results fixed Var(eta_Ka) and Var(eta_F) to zero (not
    # significantly different from zero); the typical Ka and F values are
    # therefore applied without an eta term.
    ka      <- exp(lka + e_conmed_qprl_oral_ka * CONMED_QPRL_ORAL)
    cl      <- exp(lcl + e_conmed_qprl_oral_cl * CONMED_QPRL_ORAL + etalcl)
    vc      <- exp(lvc + etalvc)
    vp      <- exp(lvp + etalvp)
    q       <- exp(lq  + etalq)
    fdepot  <- exp(lfdepot)

    # ---- Micro-constants ----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---- ODE system: standard 2-compartment with first-order absorption ----
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ---- Bioavailability on the depot compartment (oral route only) ----
    f(depot) <- fdepot

    # ---- Observation: cephalexin plasma concentration (mg/L) ----
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
