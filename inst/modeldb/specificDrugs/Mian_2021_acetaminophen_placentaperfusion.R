Mian_2021_acetaminophen_placentaperfusion <- function() {
  description <- "Ex vivo (human term placenta, dual-side recirculating single-cotyledon perfusion). Seven-compartment mechanistic in silico cotyledon perfusion model of acetaminophen placental transfer: maternal reservoir, maternal perfusate in the intervillous space, intervillous interstitial space, trophoblasts, intravillous interstitial space, fetal perfusate in the villous capillaries, and fetal reservoir. Perfusate recirculates between each reservoir and its half of the cotyledon; drug crosses the endothelium into the interstitial spaces (Schmitt interstitial partition coefficients), the apical trophoblast membrane with asymmetric influx/efflux permeability factors f_in and f_out, and the basolateral trophoblast membrane, with a fitted trophoblast:perfusate partition coefficient. Fitted in MoBi to ex vivo perfusion data (Mian 2021 Equations 1-13). The paper's whole-body maternal-fetal PBPK model (PK-Sim / MoBi), into which the fitted placental parameters were copied, is NOT reproduced here; only the self-contained ex vivo placental-transfer model is."
  reference <- paste(
    "Mian P, Nolan B, van den Anker JN, van Calsteren K, Allegaert K, Lakhi N,",
    "Dallmann A. Mechanistic Coupling of a Novel in silico Cotyledon Perfusion",
    "Model and a Physiologically Based Pharmacokinetic Model to Predict Fetal",
    "Acetaminophen Pharmacokinetics at Delivery. Front Pediatr. 2021;9:733520.",
    "doi:10.3389/fped.2021.733520. PMCID: PMC8496351.",
    "Structure: Figure 3, Table 1 and Equations 1-13; fixed inputs: In silico",
    "Cotyledon Perfusion Model section; fitted values: Results. Compartment",
    "volume fractions, cotyledon perfusate volumes and surface-area scaling are",
    "not printed in the paper and are taken from the authors' MoBi 9.1 project",
    "that the paper states is shared on the Open Systems Pharmacology GitHub",
    "(github.com/Open-Systems-Pharmacology/Pregnancy-Models,",
    "CotyledonPerfusionModel/CotyledonPerfusionModel.mbp3, commit 73b0acd4f0,",
    "2021-09-23). Observed ex vivo data: Conings S et al. (reference 20 of the",
    "source paper), as deposited in the same MoBi project.",
    sep = " "
  )
  vignette <- "Mian_2021_acetaminophen_placentaperfusion"
  units <- list(time = "min", dosing = "umol", concentration = "umol/L")

  # Acetaminophen is added either to the maternal reservoir (maternal-to-fetal
  # experiments) or to the fetal reservoir (fetal-to-maternal experiments).
  dosing <- c("maternal_reservoir", "fetal_reservoir")

  # The seven states are the compartments of the ex vivo perfusion circuit and
  # the perfused cotyledon (Figure 3, Table 1); none of them is a whole-body
  # compartment, so they are declared paper-specific rather than registered.
  paper_specific_compartments <- c(
    "maternal_reservoir",
    "maternal_perfusate",
    "maternal_interstitial",
    "trophoblast",
    "fetal_interstitial",
    "fetal_perfusate",
    "fetal_reservoir"
  )

  # The reservoirs and the cotyledon perfusate compartments hold the perfusion
  # medium (a buffer with bovine serum albumin), which is not a biological
  # matrix in the specimen vocabulary, hence "not applicable"; the interstitial
  # spaces and the trophoblasts are placental tissue.
  compartmentData <- list(
    maternal_reservoir = list(analyte = "acetaminophen", units = "umol", specimen = "not applicable", verified = TRUE),
    maternal_perfusate = list(analyte = "acetaminophen", units = "umol", specimen = "not applicable", verified = TRUE),
    maternal_interstitial = list(analyte = "acetaminophen", units = "umol", specimen = "tissue", verified = TRUE),
    trophoblast = list(analyte = "acetaminophen", units = "umol", specimen = "tissue", verified = TRUE),
    fetal_interstitial = list(analyte = "acetaminophen", units = "umol", specimen = "tissue", verified = TRUE),
    fetal_perfusate = list(analyte = "acetaminophen", units = "umol", specimen = "not applicable", verified = TRUE),
    fetal_reservoir = list(analyte = "acetaminophen", units = "umol", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "ex vivo (human term placenta, isolated single-cotyledon dual perfusion)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    n_experiments = 14L,
    n_observations = 455L,
    system = paste(
      "Recirculating dual perfusion of one cotyledon per placenta (data of Conings et al.,",
      "reference 20 of the source paper). Perfusate albumin (bovine serum albumin)",
      "40 mg/mL maternal and 30 mg/mL fetal."
    ),
    dose_range = paste(
      "acetaminophen at an initial reservoir concentration of about 10 mg/L (Figure 4 legend),",
      "added to the fetal reservoir in 4 experiments (Figure 4A-D) and to the maternal",
      "reservoir in 10 experiments (Figure 4E-N); 17.0-18.9 umol per experiment"
    ),
    perfusion = "maternal flow 14 mL/min, fetal flow 6 mL/min; maternal and fetal reservoir volumes 280 and 284 mL",
    disease_state = "not applicable (ex vivo placental tissue from term deliveries)",
    regions = "Belgium (Leuven)",
    notes = paste(
      "The fitted parameters were copied into a whole-body maternal-fetal PBPK model that",
      "was evaluated against umbilical-vein concentrations after oral acetaminophen 1000 mg",
      "at delivery (Nitsche et al., n = 34; Mehraban et al., n = 43; Table 2). That",
      "platform PBPK model is out of scope for this file."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Perfusion circuit and cotyledon volumes. Experimental settings, not
    # estimated. Volumes in L, flows in L/min.
    # ------------------------------------------------------------------
    q_mat <- fixed(0.014); label("Maternal perfusate flow rate, Q_M (L/min)")                      # Text below Eq 2: '14 mL/min and 6 mL/min for the flow rate in the maternal and fetal system'
    q_fet <- fixed(0.006); label("Fetal perfusate flow rate, Q_F (L/min)")                         # Text below Eq 2, same sentence
    v_mres <- fixed(0.28); label("Maternal reservoir volume (L)")                                  # Text below Eq 2: '280 and 284 mL for the maternal and fetal reservoir volume'
    v_fres <- fixed(0.284); label("Fetal reservoir volume (L)")                                    # Text below Eq 2, same sentence
    v_mcot <- fixed(0.023); label("Intervillous (maternal) cotyledon fraction volume (L)")         # Text below Eq 4: 'The volumes of the intervillous and intravillous cotyledon fraction were assumed to be 23 and 35 mL'
    v_fcot <- fixed(0.035); label("Intravillous (fetal) cotyledon fraction volume (L)")            # Text below Eq 4, same sentence
    v_mperf <- fixed(0.02); label("Maternal perfusate volume in the cotyledon (L)")                # Not printed; authors' MoBi project, 'Cotyledon (maternal fraction)|Perfusate|Volume' = 0.02 L
    v_fperf <- fixed(0.02); label("Fetal perfusate volume in the cotyledon (L)")                   # Not printed; authors' MoBi project, 'Cotyledon (fetal fraction)|Perfusate|Volume' = 0.02 L
    fvas_mcot <- fixed(0.674); label("Vascular fraction of the intervillous cotyledon fraction (fraction)")   # Not printed; authors' MoBi project, 'Cotyledon (maternal fraction)|Fraction vascular' = 0.674 (used only to scale the endothelial surface area)
    fvas_fcot <- fixed(0.168); label("Vascular fraction of the intravillous cotyledon fraction (fraction)")   # Not printed; authors' MoBi project, 'Cotyledon (fetal fraction)|Fraction vascular' = 0.168
    fint_mcot <- fixed(0.1018125); label("Interstitial fraction of the intervillous cotyledon fraction (fraction)") # Not printed; authors' MoBi project, 'Cotyledon (maternal fraction)|Fraction interstitial' = 0.1018125
    fint_fcot <- fixed(0.4065634); label("Interstitial fraction of the intravillous cotyledon fraction (fraction)") # Not printed; authors' MoBi project, 'Cotyledon (fetal fraction)|Fraction interstitial' = 0.4065634; trophoblast fraction = 1 - fvas_fcot - fint_fcot

    # ------------------------------------------------------------------
    # Protein binding and interstitial composition (Equation 5, Schmitt).
    # ------------------------------------------------------------------
    fu <- fixed(0.84); label("Fraction unbound in maternal perfusate (fraction)")                  # Text below Eq 4: 'This resulted in an unbound fraction of 0.84 and 0.88 in maternal and fetal perfusate' (MoBi project: 0.8408769)
    fu_fetus <- fixed(0.88); label("Fraction unbound in fetal perfusate (fraction)")               # Text below Eq 4, same sentence (MoBi project: 0.8757135)
    fwater_int <- fixed(0.935); label("Fractional water content of the interstitial space (fraction)") # Text below Eq 5: f_water_int 'assumed to be the same than for the intervillous interstitial space [0.935 (28)]'
    fwater_perf <- fixed(0.926); label("Fractional water content of the perfusate (fraction)")     # Text below Eq 5: f_water_perf 'similar to the fractional volume content reported for plasma [0.926 (28)]'
    fprot_int_perf <- fixed(0.37); label("Interstitial-to-perfusate protein content ratio (ratio)") # Text below Eq 5: ratio f_protein_int / f_protein_perf 'assumed to be the same as in adult tissue [0.37 (26)]'

    # ------------------------------------------------------------------
    # Permeabilities (dm/min) and surface areas (dm^2); P x SA is in L/min.
    # ------------------------------------------------------------------
    p_endo <- fixed(10); label("Endothelial permeability, P_endo (dm/min)")                        # Text below Eq 4: P_endo 'set to a value of 100 cm/min' = 10 dm/min
    p_tro <- fixed(0.00429); label("Trophoblast membrane permeability, P (dm/min)")                # Text below Eq 11: 'resulting in a value of 4.29 x 10^-2 cm/min for acetaminophen' = 4.29e-3 dm/min
    sa_villi_placenta <- fixed(1178); label("Surface area of all fetal villi in the term placenta (dm^2)") # Text below Eq 11: 'the absolute surface area of all fetal villi in the term placenta, ~1178 dm2 (27)'
    n_cotyledon <- fixed(35); label("Number of cotyledons per term placenta (count)")              # Text below Eq 11: 'the average number of cotyledons in the placenta which varies around 35 at term (34)'
    k_sa_endo <- fixed(9500); label("Endothelial surface area per unit vascular volume (1/dm)")    # Not printed; authors' MoBi project, 'SA proportionality factor' = 9500 1/dm (SA_perf:int = k x fraction vascular x volume)

    # ------------------------------------------------------------------
    # Fitted placental transfer parameters (Results, scenario 4: asymmetric
    # transfer plus a common placental partition coefficient; fitted value
    # +/- 95% confidence interval, Monte-Carlo optimisation in MoBi).
    # ------------------------------------------------------------------
    lfin <- log(0.060); label("Apical influx permeability factor, f_in, maternal perfusate to trophoblast (unitless)")   # Results: 'The fitted values +/- 95% confidence intervals for f in and fout were 0.060 +/- 0.0058' (MoBi project: 0.0595778)
    lfout <- log(0.051); label("Apical efflux permeability factor, f_out, trophoblast to maternal perfusate (unitless)") # Results: 'and 0.051 +/- 0.0061, respectively' (MoBi project: 0.0507334)
    lkp_trophoblast <- log(4.31); label("Trophoblast-to-perfusate partition coefficient, K_FM_cell:perf = K_F_cell:perf (unitless)") # Results: 'The fitted value +/- 95% confidence interval for the placental partition coefficients ... was 4.31 +/- 0.57' (MoBi project: 4.313041)

    # The optimisation was a deterministic least-squares fit to pooled ex
    # vivo data; no between-placenta variability or residual-error model was
    # estimated or reported, so this file carries typical values only.
  })

  model({
    fin <- exp(lfin)
    fout <- exp(lfout)
    kp_trophoblast <- exp(lkp_trophoblast)

    # Compartment volumes (L). The two cotyledon perfusate volumes are set
    # directly; the interstitial and trophoblast volumes are fractions of the
    # intervillous (maternal) and intravillous (fetal) cotyledon volumes.
    v_mint <- fint_mcot * v_mcot
    v_fint <- fint_fcot * v_fcot
    v_tro <- (1 - fvas_fcot - fint_fcot) * v_fcot

    # Surface areas (dm^2). Endothelial surface areas scale with the vascular
    # volume; the basolateral (interstitial/trophoblast) surface area uses the
    # volume-based geometric-similarity scaling of the authors' MoBi project
    # ((V [mL] / 1.2)^0.75 x 7.54 x 100); the apical surface area is the
    # whole-placenta villous surface divided by the number of cotyledons.
    sa_mperf_int <- k_sa_endo * fvas_mcot * v_mcot
    sa_fperf_int <- k_sa_endo * fvas_fcot * v_fcot
    sa_fint_cell <- (v_fcot * 1000 / 1.2)^0.75 * 7.54 * 100
    sa_villi <- sa_villi_placenta / n_cotyledon

    # Interstitial-to-perfusate partition coefficients (Equation 5, Schmitt),
    # with the maternal or fetal fraction unbound.
    kp_mint <- (fwater_int + fprot_int_perf * (1 / fu - fwater_perf)) * fu
    kp_fint <- (fwater_int + fprot_int_perf * (1 / fu_fetus - fwater_perf)) * fu_fetus

    # Concentrations (umol/L).
    c_mres <- maternal_reservoir / v_mres
    c_mperf <- maternal_perfusate / v_mperf
    c_mint <- maternal_interstitial / v_mint
    c_tro <- trophoblast / v_tro
    c_fint <- fetal_interstitial / v_fint
    c_fperf <- fetal_perfusate / v_fperf
    c_fres <- fetal_reservoir / v_fres

    # Intercompartmental fluxes (umol/min), Equations 1-11.
    j_mres_mperf <- q_mat * (c_mres - c_mperf)                                                     # Eq 1
    j_fres_fperf <- q_fet * (c_fres - c_fperf)                                                     # Eq 2
    j_mperf_mint <- fu * p_endo * sa_mperf_int * (c_mperf - c_mint / kp_mint)                      # Eq 3
    j_fperf_fint <- fu_fetus * p_endo * sa_fperf_int * (c_fperf - c_fint / kp_fint)                # Eq 4
    j_fint_tro <- p_tro * sa_fint_cell *
      (c_fint * fu_fetus / kp_fint - c_tro * fu / kp_trophoblast)                                  # Eq 10
    j_mperf_tro <- p_tro * sa_villi * fu *
      (fin * c_mperf - fout * c_tro / kp_trophoblast)                                              # Eq 11

    # Equations 12-13 (the E matrix): amount balances for the seven states.
    d/dt(maternal_reservoir) <- -j_mres_mperf
    d/dt(maternal_perfusate) <- j_mres_mperf - j_mperf_mint - j_mperf_tro
    d/dt(maternal_interstitial) <- j_mperf_mint
    d/dt(trophoblast) <- j_mperf_tro + j_fint_tro
    d/dt(fetal_interstitial) <- j_fperf_fint - j_fint_tro
    d/dt(fetal_perfusate) <- j_fres_fperf - j_fperf_fint
    d/dt(fetal_reservoir) <- -j_fres_fperf

    # Observed quantities: acetaminophen in the maternal and fetal reservoirs
    # (Figure 4), plus the trophoblast concentration discussed in the text.
    Cmaternal <- c_mres
    Cfetal <- c_fres
    Ctrophoblast <- c_tro
  })
}
