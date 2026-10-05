Bukkems_2022_doravirine_placentaperfusion <- function() {
  description <- "Ex vivo (human term placenta, dual-side recirculating closed-closed single-cotyledon perfusion). Five-compartment mechanistic placenta model of doravirine placental transfer: maternal reservoir, maternal part of the placenta, placental barrier, fetal part of the placenta and fetal reservoir, all well stirred. Perfusate recirculates between each reservoir and its half of the cotyledon at the experimental flow; drug crosses the maternal-facing and the fetal-facing placental barrier by bidirectional passive diffusion (CLpdm and CLpdf) acting on the unbound concentrations, with the ex vivo perfusate fraction unbound FU and the placental-barrier fraction unbound FUp. Placental sub-volumes and both transfer clearances scale linearly with the perfused cotyledon weight (reference 42 g = 44.02 mL). Log-normal between-placenta variability on CLpdm and CLpdf fixed at 100% CV; combined additive-plus-proportional residual error whose proportional part differs between maternal-to-fetal and fetal-to-maternal experiments. Fitted in NONMEM 7.4 (Bukkems 2022 Equations 1, 2, 4, 6, 7; Tables 1 and 2). The paper's Simcyp whole-body maternal-fetal pregnancy PBPK model, into which the fitted transfer clearances were imputed, is NOT reproduced here; only the self-contained ex vivo placental-transfer model is."
  reference <- paste(
    "Bukkems VE, van Hove H, Roelofsen D, Freriksen JJM, van Ewijk-Beneken Kolmer EWJ,",
    "Burger DM, van Drongelen J, Svensson EM, Greupink R, Colbers A. Prediction of",
    "Maternal and Fetal Doravirine Exposure by Integrating Physiologically Based",
    "Pharmacokinetic Modeling and Human Placenta Perfusion Experiments.",
    "Clin Pharmacokinet. 2022;61:1129-1141. doi:10.1007/s40262-022-01127-0.",
    "PMCID: PMC9349081. Structure: Section 2.2, Figure 1A and Equations 1, 2, 4, 6",
    "and 7; fixed inputs: Table 1; final estimates: Table 2; cotyledon weights:",
    "Electronic Supplementary Material Online Resource 4.",
    sep = " "
  )
  vignette <- "Bukkems_2022_doravirine_placentaperfusion"
  units <- list(time = "min", dosing = "ug", concentration = "ug/mL")

  # Doravirine is added either to the maternal reservoir (maternal-to-fetal
  # experiments) or to the fetal reservoir (fetal-to-maternal experiments).
  dosing <- c("maternal_reservoir", "fetal_reservoir")

  # The five states are the compartments of the ex vivo perfusion circuit and
  # of the perfused cotyledon (Figure 1A); none is a whole-body compartment,
  # so they are declared paper-specific rather than registered.
  paper_specific_compartments <- c(
    "maternal_reservoir",
    "maternal_placenta",
    "placental_barrier",
    "fetal_placenta",
    "fetal_reservoir"
  )

  # The reservoirs hold the perfusion buffer (human albumin in buffer), which
  # is not a biological matrix in the specimen vocabulary; the three placental
  # compartments are placental tissue.
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
      notes = "Scales the volumes of the maternal part of the placenta, the placental barrier and the fetal part of the placenta (11.55%, 11.05% and 8.25% of the cotyledon volume; Table 1), with 42 g taken as 44.02 mL (Section 2.2: 'Absolute volumes of MP, PB, and FP were scaled with the individual cotyledon weight of each perfusion experiment and were standarized to a typical cotyledon volume of 44.02 mL'). CLpdm and CLpdf are reported 'For the typical cotyledon weighing 42 g' (Table 2 footnote a) and were imputed into the PBPK model per mL of placenta (37.2 mL/min / 44.02 mL = 0.0507 L/h/mL; Section 3.3), so they are scaled linearly with cotyledon weight as well. Closed-closed cotyledon weights were 18.9-64.4 g (median 32.55 g; Online Resource 4).",
      source_name = "cotyledon weight"
    ),
    PERF_DIR_FTM = list(
      description = "Perfusion direction: 1 = doravirine added to the fetal circulation (fetal-to-maternal experiment); 0 = added to the maternal circulation (maternal-to-fetal experiment)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Selects the proportional residual error: 7.5% CV after dosing in the maternal compartment (PERF_DIR_FTM = 0) and 12.6% CV after dosing in the fetal compartment (PERF_DIR_FTM = 1) (Table 2). The structural model is the same in both directions; the dose record goes to maternal_reservoir or fetal_reservoir accordingly.",
      source_name = "dosing compartment"
    )
  )

  population <- list(
    species = "ex vivo (human term placenta, isolated single-cotyledon dual perfusion)",
    n_subjects = 8L,
    n_studies = 1L,
    n_experiments = 8L,
    system = paste(
      "Closed-closed (both circulations recirculating) dual-side perfusion of one cotyledon",
      "per placenta, 180 min, maternal buffer human albumin 29 g/L and fetal buffer 32 g/L;",
      "maternal flow 12 mL/min, fetal flow 6 mL/min, 200 mL per reservoir; antipyrine",
      "100 mg/L as the overlap control (Section 2.1)."
    ),
    dose_range = paste(
      "doravirine 0.96 mg/L in the dosed reservoir (192 ug in 200 mL), added to the",
      "maternal circulation in 4 experiments and to the fetal circulation in 4"
    ),
    gestational_age_range = "37-40 weeks at delivery (Online Resource 4)",
    cotyledon_weight_range = "18.9-64.4 g (Online Resource 4)",
    disease_state = "not applicable (ex vivo placental tissue from term deliveries; 7 of 8 by caesarean section)",
    regions = "the Netherlands (Nijmegen)",
    notes = paste(
      "The fitted transfer clearances were converted to 0.0507 and 0.0075 L/h/mL placenta",
      "and imputed into a Simcyp V20 whole-body maternal-fetal pregnancy PBPK model for",
      "doravirine 100 mg QD and BID at 26, 32 and 40 weeks of gestation; that platform",
      "PBPK model is out of scope for this file. A separate closed-open perfusion fit is",
      "packaged as Bukkems_2022_doravirine_placentaperfusion_closedopen."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Experimental settings and placental physiology for the typical 42 g
    # (44.02 mL) cotyledon (Table 1). Not estimated. Volumes mL, flows mL/min;
    # the placental sub-volumes are fractions of the cotyledon volume.
    # ------------------------------------------------------------------
    q_mat <- fixed(12)
    label("Maternal perfusate flow rate, Q_M (mL/min)") # Table 1: 'Q M, mL/min 12, Experimental condition'
    q_fet <- fixed(6)
    label("Fetal perfusate flow rate, Q_F (mL/min)") # Table 1: 'Q F, mL/min 6, Experimental condition'
    v_mres <- fixed(200)
    label("Maternal reservoir volume, V_MR (mL)") # Table 1: 'V MR, mL 200, Experimental condition'
    v_fres <- fixed(200)
    label("Fetal reservoir volume, V_FR (mL)") # Table 1: 'V FR, mL 200, Experimental condition'
    v_cot <- fixed(44.02)
    label("Volume of the reference 42 g cotyledon (mL)") # Table 1 and Section 2.2: 'Standardized for placenta of 42 g = 44.02 mL'
    f_vmp <- fixed(0.1155)
    label("Maternal part of the placenta as a fraction of cotyledon volume (fraction)") # Table 1: 'V MP, mL 5.08, 11.55% of total placental volume'
    f_vpb <- fixed(0.1105)
    label("Placental barrier as a fraction of cotyledon volume (fraction)") # Table 1: 'V PB, mL 4.86, 11.05% of total placental volume'
    f_vfp <- fixed(0.0825)
    label("Fetal part of the placenta as a fraction of cotyledon volume (fraction)") # Table 1: 'V FP, mL 3.63, 8.25% of total placental volume'
    fu <- fixed(0.529)
    label("Ex vivo fraction unbound in the perfusion buffer, FU (fraction)") # Table 1: 'FU 0.529, Measured'
    fu_pb <- fixed(0.01)
    label("Fraction unbound in the placental barrier, FUp (fraction)") # Table 1: 'FU p 0.01, Estimated with Simcyp PBPK simulator version 20'

    # ------------------------------------------------------------------
    # Estimated intrinsic placental transfer clearances (Table 2), typical
    # 42 g cotyledon.
    # ------------------------------------------------------------------
    lcl_pdm <- log(37.2)
    label("Passive diffusion clearance over the maternal-facing barrier, CLpdm (mL/min)") # Table 2: 'CL pdm, mL/min 37.2 (95% CI from SIR 19.8-73.7)'
    lcl_pdf <- log(5.5)
    label("Passive diffusion clearance over the fetal-facing barrier, CLpdf (mL/min)") # Table 2: 'CL pdf, mL/min 5.5 (3.1-9.8)'

    # Between-placenta variability fixed at 100% CV; variance = log(1 + 1^2).
    etalcl_pdm ~ fixed(0.6931472) # Table 2: 'IIV CL pdm, % Fixed to 100', footnote b: CV = sqrt(exp(variance) - 1)
    etalcl_pdf ~ fixed(0.6931472) # Table 2: 'IIV CL pdf, % Fixed to 100', footnote b: CV = sqrt(exp(variance) - 1)

    # ------------------------------------------------------------------
    # Residual error (Table 2). Proportional SDs back-transformed from the
    # reported CV with footnote b: SD = sqrt(log(1 + CV^2)).
    # ------------------------------------------------------------------
    addSd <- 0.000007
    label("Additive residual error (ug/mL)") # Table 2: 'Additive residual error, ug/mL 0.000007'
    propSd_mtf <- 0.07490
    label("Proportional residual error after dosing in the maternal compartment (fraction)") # Table 2: 'Proportional residual error after dosing in maternal compartment, % 7.5'; sqrt(log(1 + 0.075^2)) = 0.07490
    propSd_ftm <- 0.12551
    label("Proportional residual error after dosing in the fetal compartment (fraction)") # Table 2: 'Proportional residual error after dosing in fetal compartment, % 12.6'; sqrt(log(1 + 0.126^2)) = 0.12551
  })

  model({
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

    # Equations 1, 2, 4, 6 and 7 (diffusion-only transfer model). The paper
    # writes each balance as dN/dt = (flux terms in N) / V; the N on the
    # right-hand side are concentrations, so each flux is a flow or clearance
    # times a concentration and the states are amounts (ug).
    d/dt(maternal_reservoir) <- q_mat * (c_mp - c_mr) # Eq 1
    d/dt(maternal_placenta) <- q_mat * (c_mr - c_mp) +
      cl_pdm * (c_pb * fu_pb - c_mp * fu) # Eq 2
    d/dt(placental_barrier) <- cl_pdm * (c_mp * fu - c_pb * fu_pb) +
      cl_pdf * (c_fp * fu - c_pb * fu_pb) # Eq 4
    d/dt(fetal_placenta) <- q_fet * (c_fr - c_fp) +
      cl_pdf * (c_pb * fu_pb - c_fp * fu) # Eq 6
    d/dt(fetal_reservoir) <- q_fet * (c_fp - c_fr) # Eq 7

    # Observed: total doravirine in the maternal and fetal reservoirs
    # (Figure 3, Online Resource 7). One residual-error model for both
    # reservoirs (Table 2 has no reservoir dimension) with a direction-specific
    # proportional part; each endpoint needs its own error variables, so the
    # shared values are copied to both.
    Cmaternal <- c_mr
    Cfetal <- c_fr
    addSd_maternal <- addSd
    addSd_fetal <- addSd
    propSd_dir <- propSd_mtf * (1 - PERF_DIR_FTM) + propSd_ftm * PERF_DIR_FTM
    propSd_maternal <- propSd_dir
    propSd_fetal <- propSd_dir
    Cmaternal ~ add(addSd_maternal) + prop(propSd_maternal)
    Cfetal ~ add(addSd_fetal) + prop(propSd_fetal)
  })
}
