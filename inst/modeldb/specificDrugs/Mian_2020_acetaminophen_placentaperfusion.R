Mian_2020_acetaminophen_placentaperfusion <- function() {
  description <- "Ex vivo (human term placenta, dual-side recirculating single-cotyledon perfusion). Four-compartment mechanistic model of acetaminophen maternal-to-fetal transfer across an isolated perfused cotyledon: a maternal reservoir and a fetal reservoir, each recirculating through its own half of the cotyledon (maternal and fetal parts), which exchange drug by symmetric passive diffusion (transcotyledon diffusion clearance D_cot) and partition against the perfusate with coefficient K_f,m; a first-order placental elimination K_pe acts on the fetal part. Fitted in MONOLIX to ex vivo perfusion data (Mian 2020 Equations 5-8, Table 1). The paper's whole-body fetal-maternal PBPK model (MoBi / PK-Sim), into which D_cot and K_f,m were up-scaled, is NOT reproduced here; only the transferable ex vivo placental-transfer model is."
  reference <- paste(
    "Mian P, Allegaert K, Conings S, Annaert P, Tibboel D, Pfister M,",
    "van Calsteren K, van den Anker JN, Dallmann A.",
    "Integration of Placental Transfer in a Fetal-Maternal Physiologically Based",
    "Pharmacokinetic Model to Characterize Acetaminophen Exposure and Metabolic",
    "Clearance in the Fetus. Clin Pharmacokinet. 2020;59(7):911-925.",
    "doi:10.1007/s40262-020-00861-7. PMCID: PMC7329787.",
    "Perfusion set-up: Section 2.5.1.1; model structure and fixed volumes:",
    "Section 2.5.1.2, Figure 3 and Equations 5-8; estimates: Table 1;",
    "model fit: Figure 4.",
    sep = " "
  )
  vignette <- "Mian_2020_acetaminophen_placentaperfusion"
  units <- list(time = "min", dosing = "ug", concentration = "ug/mL")

  # In the ex vivo experiment acetaminophen is added to the maternal reservoir,
  # so doses land in `maternal_reservoir`, not in `depot` / `central`.
  dosing <- c("maternal_reservoir")

  # The four states are the perfusion circuit's reservoirs and the two halves
  # of the perfused cotyledon (Figure 3); none of them is a whole-body
  # compartment, so they are declared paper-specific rather than registered.
  paper_specific_compartments <- c(
    "maternal_reservoir",
    "maternal_cotyledon",
    "fetal_cotyledon",
    "fetal_reservoir"
  )

  # The reservoirs hold the perfusion medium, which is not a biological matrix
  # in the specimen vocabulary, hence "not applicable"; the two cotyledon
  # states are the perfused placental tissue.
  compartmentData <- list(
    maternal_reservoir = list(analyte = "acetaminophen", units = "ug", specimen = "not applicable", verified = TRUE),
    maternal_cotyledon = list(analyte = "acetaminophen", units = "ug", specimen = "tissue", verified = TRUE),
    fetal_cotyledon = list(analyte = "acetaminophen", units = "ug", specimen = "tissue", verified = TRUE),
    fetal_reservoir = list(analyte = "acetaminophen", units = "ug", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "ex vivo (human term placenta, isolated single-cotyledon dual perfusion)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    system = paste(
      "Recirculating (closed-closed) dual perfusion of one intact cotyledon per placenta,",
      "started within 30 min of delivery; the chorionic artery and vein were cannulated",
      "(Section 2.5.1.1; set-up published as reference 18 of the source paper)."
    ),
    dose_range = "acetaminophen 10 ug/mL added to the maternal reservoir (maternal-to-fetal transfer experiment)",
    perfusion = "maternal flow 14 mL/min, fetal flow 6 mL/min; mean maternal and fetal reservoir volumes 280 and 284 mL",
    sampling = "maternal and fetal reservoirs at 0, 3, 6, 10, 15, 20 and 30 min, every 15 min to 150 min, then every 30 min to 210 min",
    ethics = "University Hospitals Leuven ethics board s54819; EudraCT 2012-004580-51; NCT02622802",
    disease_state = "not applicable (ex vivo placental tissue from term deliveries)",
    regions = "Belgium (Leuven)",
    notes = paste(
      "The number of perfused placentas is not stated in this paper (it is in the",
      "earlier perfusion report, reference 18). The paper's in vivo arm -- 34 women and",
      "their newborns at a median 39 weeks gestation after oral acetaminophen 1000 mg at",
      "delivery (reference 19) -- evaluates the whole-body fetal-maternal PBPK model, not",
      "this ex vivo model, and is out of scope for this file."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Perfusion circuit settings (Section 2.5.1.1) and cotyledon volumes
    # (Section 2.5.1.2). All are experimental constants that were not
    # estimated; change them to describe a different perfusion set-up.
    # ------------------------------------------------------------------
    q_mat <- fixed(14); label("Maternal perfusion flow rate, Q_m (mL/min)")                 # Section 2.5.1.1: 'The fetal (Qf) and maternal (Qm) circulations were established at a flow of 6 and 14 mL/min, respectively'
    q_fet <- fixed(6); label("Fetal perfusion flow rate, Q_f (mL/min)")                     # Section 2.5.1.1, same sentence
    v_mres <- fixed(280); label("Maternal reservoir volume, V_m (mL)")                      # Section 2.5.1.1: 'The mean maternal and fetal reservoir volumes (Vm and Vf) were 280 and 284 mL, respectively'
    v_fres <- fixed(284); label("Fetal reservoir volume, V_f (mL)")                         # Section 2.5.1.1, same sentence
    v_mcot <- fixed(23); label("Maternal part of the cotyledon volume, V_mp (mL)")          # Section 2.5.1.2: 'The maternal cotyledon volume (Vmp) was assumed to be 23 mL' (reference 44); 'Vmp and Vfp ... were fixed'
    v_fcot <- fixed(35); label("Fetal part of the cotyledon volume, V_fp (mL)")             # Section 2.5.1.2: 'and the fetal cotyledon volume (Vfp) 35 mL' (total cotyledon 58 mL)

    # ------------------------------------------------------------------
    # Estimated placental-transfer parameters (Table 1, 'mean (residual
    # standard error)'; the bracketed numbers are %, and Section 3.2 calls the
    # 229 for K_pe a 'high imprecision').
    # ------------------------------------------------------------------
    ldcot <- log(36); label("Transcotyledon passive diffusion clearance, D_cot (mL/min)")   # Table 1: D_cot = 36 (81) mL/min; same value in both directions (Section 2.5.1.2: 'no polarity was assumed')
    lkpe <- log(0.0126); label("Placental elimination rate constant from the fetal part, K_pe (1/min)") # Table 1: K_pe = 0.0126 (229) /min
    lkfm <- log(0.737); label("Partition coefficient between the fetal and maternal compartment, K_f,m (unitless)") # Table 1: K_f,m = 0.737 (35)

    # Between-placenta variability and the residual-error model were
    # estimated (Section 2.5.1.2: IIV 'was tested for significance on all
    # parameters except Vmp and Vfp'; proportional, constant and mixed error
    # models were investigated) but the final variances and error magnitudes
    # are not reported, so this file carries the typical values only.
  })

  model({
    dcot <- exp(ldcot)
    kpe <- exp(lkpe)
    kfm <- exp(lkfm)

    # Concentrations (ug/mL) from amounts (ug) and volumes (mL).
    c_mres <- maternal_reservoir / v_mres
    c_mcot <- maternal_cotyledon / v_mcot
    c_fcot <- fetal_cotyledon / v_fcot
    c_fres <- fetal_reservoir / v_fres

    # Equations 5-8. The paper writes each equation as
    # dN/dt = (flux terms in N) / V, i.e. the N on its right-hand sides are
    # concentrations; multiplying through by V gives the amount balances
    # below. K_f,m divides the cotyledon concentration on BOTH the maternal
    # and the fetal outflow, as printed.
    d/dt(maternal_reservoir) <- q_mat * c_mcot / kfm - q_mat * c_mres                          # Eq 5
    d/dt(maternal_cotyledon) <- q_mat * c_mres - q_mat * c_mcot / kfm -
      dcot * c_mcot + dcot * c_fcot                                                            # Eq 6
    d/dt(fetal_cotyledon) <- q_fet * c_fres - q_fet * c_fcot / kfm +
      dcot * c_mcot - dcot * c_fcot - kpe * v_fcot * c_fcot                                    # Eq 7
    d/dt(fetal_reservoir) <- q_fet * c_fcot / kfm - q_fet * c_fres                             # Eq 8

    # Observed quantities: acetaminophen in the maternal and fetal reservoirs.
    Cmaternal <- c_mres
    Cfetal <- c_fres
  })
}
