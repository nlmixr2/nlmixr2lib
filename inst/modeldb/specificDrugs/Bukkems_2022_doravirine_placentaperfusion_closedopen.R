Bukkems_2022_doravirine_placentaperfusion_closedopen <- function() {
  description <- "Ex vivo (human term placenta, dual-side closed-open single-cotyledon perfusion). Five-compartment mechanistic placenta model of doravirine placental transfer fitted to closed-open perfusions, in which the dosed circulation recirculates and the other circulation is single-pass: maternal reservoir, maternal part of the placenta, placental barrier, fetal part of the placenta and fetal reservoir, all well stirred. On the closed side perfusate recirculates between a 200 mL reservoir and its half of the cotyledon; on the open side drug-free perfusate enters the cotyledon and the effluent passes through a 3 mL collecting reservoir to waste. Drug crosses the maternal-facing and the fetal-facing placental barrier by bidirectional passive diffusion (CLpdm and CLpdf) acting on unbound concentrations (FU, FUp). Placental sub-volumes and both clearances scale linearly with cotyledon weight (reference 42 g = 44.02 mL). Log-normal between-placenta variability on CLpdm and CLpdf fixed at 100% CV; combined additive-plus-proportional residual error, both components specific to the perfusion direction. Fitted in NONMEM 7.4 (Bukkems 2022 Section 2.5; Electronic Supplementary Material Online Resources 16, 17 and 19). This is the paper's sensitivity analysis of the perfusion configuration; the closed-closed fit used for the primary PBPK predictions is Bukkems_2022_doravirine_placentaperfusion."
  reference <- paste(
    "Bukkems VE, van Hove H, Roelofsen D, Freriksen JJM, van Ewijk-Beneken Kolmer EWJ,",
    "Burger DM, van Drongelen J, Svensson EM, Greupink R, Colbers A. Prediction of",
    "Maternal and Fetal Doravirine Exposure by Integrating Physiologically Based",
    "Pharmacokinetic Modeling and Human Placenta Perfusion Experiments.",
    "Clin Pharmacokinet. 2022;61:1129-1141. doi:10.1007/s40262-022-01127-0.",
    "PMCID: PMC9349081. Structure: Section 2.5, Equations 1-7 of the main paper",
    "adapted to the open circulation as drawn in Online Resource 16A; fixed inputs:",
    "Online Resource 17; final estimates: Online Resource 19 and Section 3.5;",
    "cotyledon weights: Online Resource 14.",
    sep = " "
  )
  vignette <- "Bukkems_2022_doravirine_placentaperfusion"
  units <- list(time = "min", dosing = "ug", concentration = "ug/mL")

  # Doravirine is added to the closed (recirculating) reservoir: the maternal
  # reservoir in maternal-to-fetal experiments, the fetal reservoir in
  # fetal-to-maternal experiments.
  dosing <- c("maternal_reservoir", "fetal_reservoir")

  # The five states are the compartments of the ex vivo perfusion circuit and
  # of the perfused cotyledon (Online Resource 16A).
  paper_specific_compartments <- c(
    "maternal_reservoir",
    "maternal_placenta",
    "placental_barrier",
    "fetal_placenta",
    "fetal_reservoir"
  )

  compartmentData <- list(
    maternal_reservoir = list(analyte = "doravirine", units = "ug", specimen = "not applicable", verified = TRUE),
    maternal_placenta = list(analyte = "doravirine", units = "ug", specimen = "tissue", verified = TRUE),
    placental_barrier = list(analyte = "doravirine", units = "ug", specimen = "tissue", verified = TRUE),
    fetal_placenta = list(analyte = "doravirine", units = "ug", specimen = "tissue", verified = TRUE),
    fetal_reservoir = list(analyte = "doravirine", units = "ug", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT_COTYLEDON = list(
      description = "Wet weight of the perfused placental cotyledon",
      units = "g",
      type = "continuous",
      reference_category = NULL,
      notes = "Scales the volumes of the maternal part of the placenta, the placental barrier and the fetal part of the placenta (11.55%, 11.05% and 8.25% of the cotyledon volume, standardised to 42 g = 44.02 mL; Online Resource 17) and, per Online Resource 19 footnote a ('A typical cotyledon weighs 42g, assuming to be equal to 44.02 mL'), the transfer clearances. Closed-open cotyledon weights were 22.5-41.4 g (median 33.15 g; Online Resource 14).",
      source_name = "cotyledon weight"
    ),
    PERF_DIR_FTM = list(
      description = "Perfusion direction: 1 = doravirine added to the fetal circulation (fetal-to-maternal experiment); 0 = added to the maternal circulation (maternal-to-fetal experiment)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "In the closed-open configuration the dosed circulation is the closed one, so the indicator also selects the open side: PERF_DIR_FTM = 0 closes the maternal circulation (200 mL maternal reservoir) and opens the fetal one (3 mL fetal collecting reservoir); PERF_DIR_FTM = 1 the reverse (Online Resource 16: 'these models were run simultaneously while the relevant compartments were turned on based on the dosing compartment covariate'; Online Resource 17: 'VMR, mL 200 or 3 ... based on dosing compartment'). Also selects the additive and proportional residual errors (Online Resource 19).",
      source_name = "dosing compartment"
    )
  )

  population <- list(
    species = "ex vivo (human term placenta, isolated single-cotyledon dual perfusion)",
    n_subjects = 6L,
    n_studies = 1L,
    n_experiments = 6L,
    system = paste(
      "Closed-open dual-side perfusion of one cotyledon per placenta: the dosed circulation",
      "recirculates (200 mL reservoir) and the other is single-pass (Online Resource 13);",
      "maternal flow 12 mL/min, fetal flow 6 mL/min (Online Resource 17)."
    ),
    dose_range = paste(
      "doravirine added to the closed circulation, maternal-to-fetal in 3 experiments and",
      "fetal-to-maternal in 3 (Section 2.5); the dosed concentration is not restated for",
      "this configuration (0.96 mg/L in the closed-closed experiments)"
    ),
    gestational_age_range = "38-40 weeks at delivery (Online Resource 14)",
    cotyledon_weight_range = "22.5-41.4 g (Online Resource 14)",
    disease_state = "not applicable (ex vivo placental tissue from term deliveries; 5 of 6 by caesarean section)",
    regions = "the Netherlands (Nijmegen)",
    notes = paste(
      "Sensitivity analysis of the perfusion configuration: these estimates were also imputed",
      "into the Simcyp pregnancy PBPK model (Figure 5G-H), which is out of scope for this file."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Experimental settings and placental physiology (Online Resource 17).
    # Not estimated. Volumes in mL, flows in mL/min.
    # ------------------------------------------------------------------
    q_mat <- fixed(12)
    label("Maternal perfusate flow rate, Q_M (mL/min)") # Online Resource 17: 'QM, mL/min 12, Experimental condition'
    q_fet <- fixed(6)
    label("Fetal perfusate flow rate, Q_F (mL/min)") # Online Resource 17: 'QF, mL/min 6, Experimental condition'
    v_res_closed <- fixed(200)
    label("Volume of the closed (recirculating) reservoir (mL)") # Online Resource 17: 'VMR, mL 200 or 3' and 'VFR, mL 200 or 3', 'based on dosing compartment'
    v_res_open <- fixed(3)
    label("Volume of the open-side collecting reservoir (mL)") # Online Resource 17, same rows
    v_cot <- fixed(44.02)
    label("Volume of the reference 42 g cotyledon (mL)") # Online Resource 17: 'Standardized for placenta of 42g = 44.02 mL'
    f_vmp <- fixed(0.1155)
    label("Maternal part of the placenta as a fraction of cotyledon volume (fraction)") # Online Resource 17: 'VMP, mL 5.08, 11.55% of total placental volume'
    f_vpb <- fixed(0.1105)
    label("Placental barrier as a fraction of cotyledon volume (fraction)") # Online Resource 17: 'VPB, mL 4.86, 11.05% of total placental volume'
    f_vfp <- fixed(0.0825)
    label("Fetal part of the placenta as a fraction of cotyledon volume (fraction)") # Online Resource 17: 'VFP, mL 3.63, 8.25% of total placental volume'
    fu <- fixed(0.529)
    label("Ex vivo fraction unbound in the perfusion buffer, FU (fraction)") # Online Resource 17: 'FU 0.529, Measured'
    fu_pb <- fixed(0.01)
    label("Fraction unbound in the placental barrier, FUp (fraction)") # Online Resource 17: 'FUp 0.01, Estimated with Simcyp PBPK simulator version 20'

    # ------------------------------------------------------------------
    # Estimated intrinsic placental transfer clearances (Online Resource
    # 19), for the typical 42 g cotyledon.
    # ------------------------------------------------------------------
    lcl_pdm <- log(11.0)
    label("Passive diffusion clearance over the maternal-facing barrier, CLpdm (mL/min)") # Online Resource 19: 'CLpdm, mL/min 11.0 (95%CI from SIR 5.6 - 22.7)'; Section 3.5
    lcl_pdf <- log(4.3)
    label("Passive diffusion clearance over the fetal-facing barrier, CLpdf (mL/min)") # Online Resource 19: 'CLpdf, mL/min 4.3 (2.3 - 8.3)'; Section 3.5

    # Between-placenta variability fixed at 100% CV; variance = log(1 + 1^2).
    etalcl_pdm ~ fixed(0.6931472) # Online Resource 19: 'IIV CLpdm, % Fixed op 100', footnote b: CV = sqrt(exp(variance) - 1)
    etalcl_pdf ~ fixed(0.6931472) # Online Resource 19: 'IIV CLpdf, % Fixed op 100', footnote b: CV = sqrt(exp(variance) - 1)

    # ------------------------------------------------------------------
    # Residual error (Online Resource 19). Proportional SDs back-transformed
    # from the reported CV with footnote b: SD = sqrt(log(1 + CV^2)).
    # ------------------------------------------------------------------
    addSd_mtf <- 0.00188
    label("Additive residual error after dosing in the maternal compartment (ug/mL)") # Online Resource 19: 'Additive residual error after dosing in maternal compartment, ug/mL 0.00188'
    addSd_ftm <- 0.00007
    label("Additive residual error after dosing in the fetal compartment (ug/mL)") # Online Resource 19: 'Additive residual error after dosing in fetal compartment, ug/mL 0.00007'
    propSd_mtf <- 0.04498
    label("Proportional residual error after dosing in the maternal compartment (fraction)") # Online Resource 19: 'Proportional residual error after dosing in maternal compartment,% 4.5'; sqrt(log(1 + 0.045^2)) = 0.04498
    propSd_ftm <- 0.08684
    label("Proportional residual error after dosing in the fetal compartment (fraction)") # Online Resource 19: 'Proportional residual error after dosing in fetal compartment, % 8.7'; sqrt(log(1 + 0.087^2)) = 0.08684
  })

  model({
    # Which circulation is single-pass: the one that was not dosed.
    open_mat <- PERF_DIR_FTM
    open_fet <- 1 - PERF_DIR_FTM
    v_mres <- v_res_closed * (1 - open_mat) + v_res_open * open_mat
    v_fres <- v_res_closed * (1 - open_fet) + v_res_open * open_fet

    # Cotyledon size: placental sub-volumes and transfer clearances are
    # proportional to the perfused cotyledon weight (42 g = 44.02 mL).
    size_cot <- WT_COTYLEDON / 42
    v_mp <- f_vmp * v_cot * size_cot
    v_pb <- f_vpb * v_cot * size_cot
    v_fp <- f_vfp * v_cot * size_cot

    cl_pdm <- exp(lcl_pdm + etalcl_pdm) * size_cot
    cl_pdf <- exp(lcl_pdf + etalcl_pdf) * size_cot

    # Total concentrations (ug/mL).
    c_mr <- maternal_reservoir / v_mres
    c_mp <- maternal_placenta / v_mp
    c_pb <- placental_barrier / v_pb
    c_fp <- fetal_placenta / v_fp
    c_fr <- fetal_reservoir / v_fres

    # Main-paper Equations 1, 2, 4, 6 and 7 with the open circulation of
    # Online Resource 16A: on the open side the perfusate entering the
    # cotyledon is drug-free and the collecting reservoir drains to waste at
    # the same flow; on the closed side the reservoir recirculates as in the
    # closed-closed model. Fluxes are flow or clearance times concentration;
    # the states are amounts (ug).
    d/dt(maternal_reservoir) <- q_mat * (c_mp - c_mr)
    d/dt(maternal_placenta) <- q_mat * ((1 - open_mat) * c_mr - c_mp) +
      cl_pdm * (c_pb * fu_pb - c_mp * fu)
    d/dt(placental_barrier) <- cl_pdm * (c_mp * fu - c_pb * fu_pb) +
      cl_pdf * (c_fp * fu - c_pb * fu_pb)
    d/dt(fetal_placenta) <- q_fet * ((1 - open_fet) * c_fr - c_fp) +
      cl_pdf * (c_pb * fu_pb - c_fp * fu)
    d/dt(fetal_reservoir) <- q_fet * (c_fp - c_fr)

    # Observed: total doravirine in the closed reservoir and in the open-side
    # effluent (Online Resources 15 and 18). One residual-error model for both
    # reservoirs (Online Resource 19 has no reservoir dimension) with
    # direction-specific additive and proportional parts; each endpoint needs
    # its own error variables, so the shared values are copied to both.
    Cmaternal <- c_mr
    Cfetal <- c_fr
    addSd_dir <- addSd_mtf * (1 - PERF_DIR_FTM) + addSd_ftm * PERF_DIR_FTM
    propSd_dir <- propSd_mtf * (1 - PERF_DIR_FTM) + propSd_ftm * PERF_DIR_FTM
    addSd_maternal <- addSd_dir
    addSd_fetal <- addSd_dir
    propSd_maternal <- propSd_dir
    propSd_fetal <- propSd_dir
    Cmaternal ~ add(addSd_maternal) + prop(propSd_maternal)
    Cfetal ~ add(addSd_fetal) + prop(propSd_fetal)
  })
}
