Fan_2026_cyclophosphamide_mouse_qsp <- function() {
  description <- "Preclinical (mouse). QSP model of donor T-cell kinetics after MHC-haploidentical hematopoietic cell transplant (B6C3F1 -> B6D2F1) with posttransplantation cyclophosphamide (PTCy). Three layers: (1) the Campagne 2019 mouse PK model of cyclophosphamide, 4-hydroxycyclophosphamide (4HCY) and carboxyethylphosphoramide mustard (CEPM), each with central, peripheral and brain extracellular-fluid compartments; (2) a sigmoid-Emax effect of central 4HCY delayed through four signal-transduction transit compartments that stimulates removal of T cells; (3) a Ganusov-Tomura recirculation model for cytotoxic (CD8+), helper (CD4+Foxp3-) and regulatory (CD4+CD25+Foxp3+) T cells between blood, spleen, lymph nodes, liver, intestine and lung, with density-dependent proliferation, blood-only death, a delayed sigmoid thymic output and saturable homing capacity in spleen, lymph node and liver (Fan 2026). Naive-pooled fit (no IIV); the paper's neural-network GVHD-score layer is not included. Caution: Emax, EC50, gamma and tau were calibrated with the authors' code, whose time stepping gave each dose as 4 x the nominal amount over the 0.025-day step ending at the nominal dose time; to reproduce the published T-cell curves give that input (see the vignette). A nominal 9-min input of the same dose produces far less T-cell depletion."
  reference <- paste(
    "Fan HH, McCune JS, Kanakry CG, Mager DE. Hybrid Systems Modeling of Donor",
    "T-Cell Responses and Graft-Versus-Host-Disease After Posttransplantation",
    "Cyclophosphamide Administration. CPT Pharmacometrics Syst Pharmacol.",
    "2026;15(4):e70252. doi:10.1002/psp4.70252.",
    "PK layer: Campagne O, et al. CNS Penetration of Cyclophosphamide and",
    "Metabolites in Mice Bearing Group 3 Medulloblastoma and Non-Tumor Bearing",
    "Mice. J Pharm Pharm Sci. 2019;22(1):612-629. doi:10.18433/jpps30608.",
    sep = " "
  )
  vignette <- "Fan_2026_cyclophosphamide_mouse_qsp"
  # Time is in days since donor-graft infusion (day 0): the thymic-output
  # function lambda(t) reads the solver time directly, so a simulation must
  # start at the day of transplant. The PK layer is body-weight normalised:
  # dose `amt` is umol of cyclophosphamide per kg body weight (mg/kg divided by
  # the molar mass 0.261086 mg/umol), the PK states hold umol/kg and the PK
  # concentrations are umol/L. The T-cell states hold cells.
  units <- list(time = "day", dosing = "umol", concentration = "umol/L")

  paper_specific_compartments <- c(
    "tc_blood",
    "th_blood",
    "treg_blood",
    "tc_spleen1",
    "th_spleen1",
    "treg_spleen1",
    "tc_spleen2",
    "th_spleen2",
    "treg_spleen2",
    "tc_lnode1",
    "th_lnode1",
    "treg_lnode1",
    "tc_lnode2",
    "th_lnode2",
    "treg_lnode2",
    "tc_liver",
    "th_liver",
    "treg_liver",
    "tc_gut1",
    "th_gut1",
    "treg_gut1",
    "tc_gut2",
    "th_gut2",
    "treg_gut2",
    "tc_lung",
    "th_lung",
    "treg_lung"
  )

  compartmentData <- list(
    central = list(analyte = "cyclophosphamide", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cyclophosphamide", units = "umol", specimen = "plasma", verified = TRUE),
    ecf = list(analyte = "cyclophosphamide", units = "umol", specimen = "brain ISF", verified = TRUE),
    central_4ohctx = list(
      analyte = "4-hydroxycyclophosphamide",
      units = "umol",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_4ohctx = list(
      analyte = "4-hydroxycyclophosphamide",
      units = "umol",
      specimen = "plasma",
      verified = TRUE
    ),
    ecf_4ohctx = list(
      analyte = "4-hydroxycyclophosphamide",
      units = "umol",
      specimen = "brain ISF",
      verified = TRUE
    ),
    central_cepm = list(
      analyte = "carboxyethylphosphoramide mustard",
      units = "umol",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_cepm = list(
      analyte = "carboxyethylphosphoramide mustard",
      units = "umol",
      specimen = "plasma",
      verified = TRUE
    ),
    ecf_cepm = list(
      analyte = "carboxyethylphosphoramide mustard",
      units = "umol",
      specimen = "brain ISF",
      verified = TRUE
    ),
    transit1 = list(analyte = "4HCY removal signal", units = "1/day", specimen = "not applicable", verified = TRUE),
    transit2 = list(analyte = "4HCY removal signal", units = "1/day", specimen = "not applicable", verified = TRUE),
    transit3 = list(analyte = "4HCY removal signal", units = "1/day", specimen = "not applicable", verified = TRUE),
    transit4 = list(analyte = "4HCY removal signal", units = "1/day", specimen = "not applicable", verified = TRUE),
    tc_blood = list(analyte = "donor CD8+ T cells", units = "cells", specimen = "whole blood", verified = TRUE),
    th_blood = list(analyte = "donor CD4+Foxp3- T cells", units = "cells", specimen = "whole blood", verified = TRUE),
    treg_blood = list(
      analyte = "donor CD4+CD25+Foxp3+ T cells",
      units = "cells",
      specimen = "whole blood",
      verified = TRUE
    ),
    tc_spleen1 = list(analyte = "donor CD8+ T cells", units = "cells", specimen = "tissue", verified = TRUE),
    th_spleen1 = list(analyte = "donor CD4+Foxp3- T cells", units = "cells", specimen = "tissue", verified = TRUE),
    treg_spleen1 = list(
      analyte = "donor CD4+CD25+Foxp3+ T cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    tc_spleen2 = list(analyte = "donor CD8+ T cells", units = "cells", specimen = "tissue", verified = TRUE),
    th_spleen2 = list(analyte = "donor CD4+Foxp3- T cells", units = "cells", specimen = "tissue", verified = TRUE),
    treg_spleen2 = list(
      analyte = "donor CD4+CD25+Foxp3+ T cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    tc_lnode1 = list(analyte = "donor CD8+ T cells", units = "cells", specimen = "tissue", verified = TRUE),
    th_lnode1 = list(analyte = "donor CD4+Foxp3- T cells", units = "cells", specimen = "tissue", verified = TRUE),
    treg_lnode1 = list(
      analyte = "donor CD4+CD25+Foxp3+ T cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    tc_lnode2 = list(analyte = "donor CD8+ T cells", units = "cells", specimen = "tissue", verified = TRUE),
    th_lnode2 = list(analyte = "donor CD4+Foxp3- T cells", units = "cells", specimen = "tissue", verified = TRUE),
    treg_lnode2 = list(
      analyte = "donor CD4+CD25+Foxp3+ T cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    tc_liver = list(analyte = "donor CD8+ T cells", units = "cells", specimen = "tissue", verified = TRUE),
    th_liver = list(analyte = "donor CD4+Foxp3- T cells", units = "cells", specimen = "tissue", verified = TRUE),
    treg_liver = list(analyte = "donor CD4+CD25+Foxp3+ T cells", units = "cells", specimen = "tissue", verified = TRUE),
    tc_gut1 = list(analyte = "donor CD8+ T cells", units = "cells", specimen = "tissue", verified = TRUE),
    th_gut1 = list(analyte = "donor CD4+Foxp3- T cells", units = "cells", specimen = "tissue", verified = TRUE),
    treg_gut1 = list(analyte = "donor CD4+CD25+Foxp3+ T cells", units = "cells", specimen = "tissue", verified = TRUE),
    tc_gut2 = list(analyte = "donor CD8+ T cells", units = "cells", specimen = "tissue", verified = TRUE),
    th_gut2 = list(analyte = "donor CD4+Foxp3- T cells", units = "cells", specimen = "tissue", verified = TRUE),
    treg_gut2 = list(analyte = "donor CD4+CD25+Foxp3+ T cells", units = "cells", specimen = "tissue", verified = TRUE),
    tc_lung = list(analyte = "donor CD8+ T cells", units = "cells", specimen = "tissue", verified = TRUE),
    th_lung = list(analyte = "donor CD4+Foxp3- T cells", units = "cells", specimen = "tissue", verified = TRUE),
    treg_lung = list(analyte = "donor CD4+CD25+Foxp3+ T cells", units = "cells", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight of the transplanted mouse",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. Used only to convert the blood T-cell amount to a concentration through the blood volume Vb = 77 uL/g x body weight (Fan 2026 Eq. 12). The paper describes each dose cohort's measured body-weight course with a piecewise-polynomial spline (Methods S2; coefficients in Data S2 and the authors' code); supply that course, converted from g to kg, as a time-varying covariate.",
      source_name = "W(t) (g)"
    )
  )

  population <- list(
    species = "mouse (B6C3F1 donor -> B6D2F1 recipient, MHC-haploidentical)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    disease_state = "T-cell-replete MHC-haploidentical murine hematopoietic cell transplant after total body irradiation; recipients received 10e6 T-cell-depleted bone marrow cells and 40e6 splenocytes on day 0 (Wachsmuth 2019 preclinical data).",
    dose_range = "Vehicle (0 mg/kg) or cyclophosphamide 5, 10, 25, 50 or 100 mg/kg/day intraperitoneally on days 3 and 4 after transplant, simulated as a 9-min zero-order input into the central compartment.",
    regions = "USA (preclinical data from the NCI laboratory of Kanakry; modeling at the University at Buffalo).",
    notes = "Each mouse contributed one destructive time point: donor T-cell subsets counted by flow cytometry in blood (cells/uL), spleen, peripheral lymph nodes and liver (cells) on days 7, 21 and 200. Cyclophosphamide and metabolite concentrations were not measured in this study; the PK layer is the published Campagne 2019 model run with its typical values. Migration, death and proliferation parameters come from the literature (Table 1 sources); m21, m31, m41 and hreg were estimated from the vehicle group, the drug-effect parameters were calibrated to day-7 observations, and the thymic-output parameters were estimated from days 21 and 200 in all groups. Naive-pooled estimation (Nelder-Mead in scipy) with a proportional residual error."
  )

  ini({
    # ---- PK layer: Campagne 2019 mouse model as listed in Fan 2026 Table S2.
    # Clearances are printed in L/h/kg and converted to L/day/kg in model().
    lcl <- fixed(log(4.4)); label("Cyclophosphamide central clearance CL (L/h/kg)") # Table S2: CY CL = 4.4
    lvc <- fixed(log(0.65)); label("Cyclophosphamide central volume Vc (L/kg)") # Table S2: CY Vc = 0.65
    lq <- fixed(log(0.18)); label("Cyclophosphamide peripheral clearance Q (L/h/kg)") # Table S2: CY Q = 0.18
    lvp <- fixed(log(0.062)); label("Cyclophosphamide peripheral volume Vp (L/kg)") # Table S2: CY Vp = 0.062
    lclin <- fixed(log(4.3e-4)); label("Cyclophosphamide ECF influx clearance CLin (L/h/kg)") # Table S2: CY CLin = 4.3e-4
    lclef <- fixed(log(2.4e-3)); label("Cyclophosphamide ECF efflux clearance CLef (L/h/kg)") # Table S2: CY CLef = 2.4e-3
    fu <- fixed(0.26); label("Cyclophosphamide plasma fraction unbound (fraction)") # Table S2: CY Fu = 0.26

    lcl_4ohctx <- fixed(log(11)); label("4HCY central clearance CL (L/h/kg)") # Table S2: 4HCY CL = 11
    lvc_4ohctx <- fixed(log(2.4)); label("4HCY central volume Vc (L/kg)") # Table S2: 4HCY Vc = 2.4
    lq_4ohctx <- fixed(log(0.075)); label("4HCY peripheral clearance Q (L/h/kg)") # Table S2: 4HCY Q = 0.075
    lvp_4ohctx <- fixed(log(0.22)); label("4HCY peripheral volume Vp (L/kg)") # Table S2: 4HCY Vp = 0.22
    lclin_4ohctx <- fixed(log(2.3e-4)); label("4HCY ECF influx clearance CLin (L/h/kg)") # Table S2: 4HCY CLin = 2.3e-4
    lclef_4ohctx <- fixed(log(3.3e-3)); label("4HCY ECF efflux clearance CLef (L/h/kg)") # Table S2: 4HCY CLef = 3.3e-3
    fu_4ohctx <- fixed(0.39); label("4HCY plasma fraction unbound (fraction)") # Table S2: 4HCY Fu = 0.39

    lcl_cepm <- fixed(log(6.4)); label("CEPM central clearance CL (L/h/kg)") # Table S2: CEPM CL = 6.4
    lvc_cepm <- fixed(log(1.1)); label("CEPM central volume Vc (L/kg)") # Table S2: CEPM Vc = 1.1
    lq_cepm <- fixed(log(1.6)); label("CEPM peripheral clearance Q (L/h/kg)") # Table S2: CEPM Q = 1.6
    lvp_cepm <- fixed(log(1.7)); label("CEPM peripheral volume Vp (L/kg)") # Table S2: CEPM Vp = 1.7
    lclin_cepm <- fixed(log(8.8e-5)); label("CEPM ECF influx clearance CLin (L/h/kg)") # Table S2: CEPM CLin = 8.8e-5
    lclef_cepm <- fixed(log(1.5e-3)); label("CEPM ECF efflux clearance CLef (L/h/kg)") # Table S2: CEPM CLef = 1.5e-3
    fu_cepm <- fixed(0.31); label("CEPM plasma fraction unbound (fraction)") # Table S2: CEPM Fu = 0.31

    lvecf <- fixed(log(0.001)); label("Brain extracellular-fluid volume Vecf, all three compounds (L/kg)") # Table S2: Vecf = 0.001

    # ---- Drug effect: sigmoid Emax of central 4HCY through four transit
    # compartments (Eqs. 2-4). Calibrated by hand to day-7 T-cell counts.
    lemax <- fixed(log(350)); label("Maximum 4HCY-stimulated T-cell removal rate Emax (1/day)") # Table 1: Emax = 350 (calibrated)
    lec50 <- fixed(log(85)); label("Central 4HCY concentration at half-maximal effect EC50 (umol/L)") # Table 1: EC50 = 85 (calibrated)
    lhill <- fixed(log(4)); label("Hill coefficient of the 4HCY effect (unitless)") # Table 1: gamma = 4 (calibrated)
    ltau <- fixed(log(0.4)); label("Mean transit time of each signal-transduction compartment tau (day)") # Table 1: tau = 0.4 (calibrated)

    # ---- T-cell migration rate constants between blood (1) and spleen (2),
    # lymph node (3), liver (4), intestine (5), lung (6); 1/h, converted to
    # 1/day in model(). Shared by the three T-cell subsets.
    lk_blood_spleen <- fixed(log(5.4)); label("Blood-to-spleen migration rate constant m12 (1/h)") # Table 1: m12 = 5.4 (Ganusov and Tomura)
    lk_spleen_blood <- log(0.530681); label("Spleen-to-blood migration rate constant m21 (1/h)") # Table 1: m21 = 0.531 (RSE 14.4%); full precision from the authors' notebook
    lk_blood_lnode <- fixed(log(0.7)); label("Blood-to-lymph-node migration rate constant m13 (1/h)") # Table 1: m13 = 0.7 (Ganusov and Tomura)
    lk_lnode_blood <- log(0.882723); label("Lymph-node-to-blood migration rate constant m31 (1/h)") # Table 1: m31 = 0.883 (RSE 14.3%); full precision from the authors' notebook
    lk_blood_liver <- fixed(log(8.61)); label("Blood-to-liver migration rate constant m14 (1/h)") # Table 1: m14 = 8.61 (Ganusov and Tomura)
    lk_liver_blood <- log(0.236853); label("Liver-to-blood migration rate constant m41 (1/h)") # Table 1: m41 = 0.237 (RSE 15.6%); full precision from the authors' notebook
    lk_blood_gut <- fixed(log(5.5)); label("Blood-to-intestine migration rate constant m15 (1/h)") # Table 1: m15 = 5.5 (Ganusov and Tomura)
    lk_gut_blood <- fixed(log(0.1)); label("Intestine-to-blood migration rate constant m51 (1/h)") # Table 1: m51 = 0.1 (Ganusov and Tomura)
    lk_blood_lung <- fixed(log(23.64)); label("Blood-to-lung migration rate constant m16 (1/h)") # Table 1: m16 = 23.64 (Ganusov and Tomura)
    lk_lung_blood <- fixed(log(1.7)); label("Lung-to-blood migration rate constant m61 (1/h)") # Table 1: m61 = 1.7 (Ganusov and Tomura)

    # ---- Death (blood only) and density-dependent proliferation (Eq. 9)
    lkdeath_tc <- fixed(log(0.018)); label("Death rate constant of cytotoxic T cells dc (1/day)") # Table 1: dc = 0.018 (Westera 2013, Borghans 2018)
    lkdeath_th <- fixed(log(0.023)); label("Death rate constant of helper T cells dh (1/day)") # Table 1: dh = 0.023 (Westera 2013, Borghans 2018)
    lkdeath_treg <- fixed(log(0.014)); label("Death rate constant of regulatory T cells dreg (1/day)") # Table 1: dreg = 0.014 (Milanez-Almeida 2015)
    lkprol_tc <- fixed(log(1)); label("Maximum proliferation rate constant of cytotoxic T cells rc (1/day)") # Table 1: rc = 1 (den Braber)
    lkprol_th <- fixed(log(1)); label("Maximum proliferation rate constant of helper T cells rh (1/day)") # Table 1: rh = 1 (den Braber)
    lkprol_treg <- fixed(log(1)); label("Maximum proliferation rate constant of regulatory T cells rreg (1/day)") # Table 1: rreg = 1 (den Braber)
    lh50_tc <- fixed(log(3.0e4)); label("Cytotoxic T-cell amount at half-maximal proliferation hc (cells)") # Table 1: hc = 3.0e4 (den Braber)
    lh50_th <- fixed(log(2.4e4)); label("Helper T-cell amount at half-maximal proliferation hh (cells)") # Table 1: hh = 2.4e4 (den Braber)
    lh50_treg <- log(4289.14); label("Regulatory T-cell amount at half-maximal proliferation hreg (cells)") # Table 1: hreg = 4.3e3 (RSE 24.6%); full precision 0.428914e4 from the authors' notebook

    # ---- Thymic output lambda_i(t) (Eq. 10)
    llambda_b <- log(28.7); label("Time after transplant at which thymic output reaches half of lambda0 (day)") # Table 1: lambda_b = 28.7 (RSE 4.64%)
    llambda_r <- fixed(log(20)); label("Steepness of the thymic-output rise lambda_r (unitless)") # Table 1: lambda_r = 20 (calibrated)
    llambda0_tc <- log(4.65186605e9); label("Zero-order thymic output of cytotoxic and helper T cells lambda0 (cells/day)") # Table 1: lambda0,c = 4.65e9 (RSE 6.84%), lambda0,h assumed equal; full precision from the authors' notebook
    llambda0_treg <- log(6.865977e7); label("Zero-order thymic output of regulatory T cells lambda0,reg (cells/day)") # Table 1: lambda0,reg = 6.87e7 (RSE 9.68%); full precision from the authors' notebook

    # ---- Homing capacities (Eq. 11): amount giving 50% inhibition of entry.
    # 'Maximum of observations' in Table 1; full precision from the notebook.
    cap_spleen_tc <- fixed(3.312715e6); label("Spleen capacity for cytotoxic T cells x2c0 (cells)") # Table 1: x2c0 = 3.31e6
    cap_spleen_th <- fixed(7.367487e6); label("Spleen capacity for helper T cells x2h0 (cells)") # Table 1: x2h0 = 7.37e6
    cap_spleen_treg <- fixed(67607.33262); label("Spleen capacity for regulatory T cells x2reg0 (cells)") # Table 1: x2reg0 = 6.76e4
    cap_lnode_tc <- fixed(8.110854e5); label("Lymph-node capacity for cytotoxic T cells x3c0 (cells)") # Table 1: x3c0 = 8.11e5
    cap_lnode_th <- fixed(1.411981e6); label("Lymph-node capacity for helper T cells x3h0 (cells)") # Table 1: x3h0 = 1.41e6
    cap_lnode_treg <- fixed(1.050772e5); label("Lymph-node capacity for regulatory T cells x3reg0 (cells)") # Table 1: x3reg0 = 1.05e5
    cap_liver_tc <- fixed(3.973779e6); label("Liver capacity for cytotoxic T cells x4c0 (cells)") # Table 1: x4c0 = 3.97e6
    cap_liver_th <- fixed(2.726987e6); label("Liver capacity for helper T cells x4h0 (cells)") # Table 1: x4h0 = 2.73e6
    cap_liver_treg <- fixed(72987.76315); label("Liver capacity for regulatory T cells x4reg0 (cells)") # Table 1: x4reg0 = 7.30e4
    hill_cap <- fixed(6); label("Exponent of the homing-capacity inhibition function (unitless)") # Eq. 11 text: exponent set to 6

    # ---- Graft and blood volume
    graft_splenocytes <- fixed(40e6); label("Donor splenocytes infused on day 0 (cells)") # Study Data: 40 x 10^6 splenocytes
    frac_tcell_graft <- fixed(0.343); label("Fraction of infused splenocytes that are T cells (fraction)") # Authors' notebook: total_Tcells = 40e6 x 34.3/100
    ratio_cd4_cd8 <- fixed(0.598); label("CD4:CD8 ratio of the infused T cells (ratio)") # Authors' notebook: R = 0.598
    frac_treg_cd4 <- fixed(0.004); label("Treg fraction of infused CD4+ T cells (fraction)") # Authors' notebook: Rreg = 0.4/100
    vblood_per_wt <- fixed(77); label("Blood volume per gram of body weight (uL/g)") # Eq. 12: Vb = W x 77 uL/g

    # ---- Residual error (Eq. 14). The paper estimated ONE proportional
    # parameter xi shared by every T-cell output; rxode2 needs a separate
    # parameter per endpoint, so each output carries the same value.
    propSd_tc_blood_conc <- 0.821; label("Proportional residual error of tc_blood_conc, the shared xi (fraction)") # Table 1: xi = 0.821 (RSE 10.4%)
    propSd_th_blood_conc <- 0.821; label("Proportional residual error of th_blood_conc, the shared xi (fraction)") # Table 1: xi = 0.821 (RSE 10.4%)
    propSd_treg_blood_conc <- 0.821; label("Proportional residual error of treg_blood_conc, the shared xi (fraction)") # Table 1: xi = 0.821 (RSE 10.4%)
    propSd_spleen_tc <- 0.821; label("Proportional residual error of spleen_tc, the shared xi (fraction)") # Table 1: xi = 0.821 (RSE 10.4%)
    propSd_spleen_th <- 0.821; label("Proportional residual error of spleen_th, the shared xi (fraction)") # Table 1: xi = 0.821 (RSE 10.4%)
    propSd_spleen_treg <- 0.821; label("Proportional residual error of spleen_treg, the shared xi (fraction)") # Table 1: xi = 0.821 (RSE 10.4%)
    propSd_lnode_tc <- 0.821; label("Proportional residual error of lnode_tc, the shared xi (fraction)") # Table 1: xi = 0.821 (RSE 10.4%)
    propSd_lnode_th <- 0.821; label("Proportional residual error of lnode_th, the shared xi (fraction)") # Table 1: xi = 0.821 (RSE 10.4%)
    propSd_lnode_treg <- 0.821; label("Proportional residual error of lnode_treg, the shared xi (fraction)") # Table 1: xi = 0.821 (RSE 10.4%)
    propSd_tc_liver <- 0.821; label("Proportional residual error of tc_liver, the shared xi (fraction)") # Table 1: xi = 0.821 (RSE 10.4%)
    propSd_th_liver <- 0.821; label("Proportional residual error of th_liver, the shared xi (fraction)") # Table 1: xi = 0.821 (RSE 10.4%)
    propSd_treg_liver <- 0.821; label("Proportional residual error of treg_liver, the shared xi (fraction)") # Table 1: xi = 0.821 (RSE 10.4%)
  })

  model({
    # ---- PK layer (Fan 2026 Methods S1 = Campagne 2019; per-kg amounts) ----
    cl <- exp(lcl) * 24
    vc <- exp(lvc)
    q <- exp(lq) * 24
    vp <- exp(lvp)
    clin <- exp(lclin) * 24
    clef <- exp(lclef) * 24
    cl_4ohctx <- exp(lcl_4ohctx) * 24
    vc_4ohctx <- exp(lvc_4ohctx)
    q_4ohctx <- exp(lq_4ohctx) * 24
    vp_4ohctx <- exp(lvp_4ohctx)
    clin_4ohctx <- exp(lclin_4ohctx) * 24
    clef_4ohctx <- exp(lclef_4ohctx) * 24
    cl_cepm <- exp(lcl_cepm) * 24
    vc_cepm <- exp(lvc_cepm)
    q_cepm <- exp(lq_cepm) * 24
    vp_cepm <- exp(lvp_cepm)
    clin_cepm <- exp(lclin_cepm) * 24
    clef_cepm <- exp(lclef_cepm) * 24
    vecf <- exp(lvecf)

    Cc <- central / vc
    Cc_4ohctx <- central_4ohctx / vc_4ohctx
    Cc_cepm <- central_cepm / vc_cepm
    Cecf <- ecf / vecf
    Cecf_4ohctx <- ecf_4ohctx / vecf
    Cecf_cepm <- ecf_cepm / vecf

    # All of CY clearance forms 4HCY and all of 4HCY clearance forms CEPM
    # (molar, Fm = 1), as in the authors' code.
    d/dt(central) <- -cl * Cc - q * (Cc - peripheral1 / vp) - fu * clin * Cc + clef * Cecf
    d/dt(peripheral1) <- q * (Cc - peripheral1 / vp)
    d/dt(ecf) <- fu * clin * Cc - clef * Cecf
    d/dt(central_4ohctx) <- cl * Cc - cl_4ohctx * Cc_4ohctx -
      q_4ohctx * (Cc_4ohctx - peripheral1_4ohctx / vp_4ohctx) -
      fu_4ohctx * clin_4ohctx * Cc_4ohctx + clef_4ohctx * Cecf_4ohctx
    d/dt(peripheral1_4ohctx) <- q_4ohctx * (Cc_4ohctx - peripheral1_4ohctx / vp_4ohctx)
    d/dt(ecf_4ohctx) <- fu_4ohctx * clin_4ohctx * Cc_4ohctx - clef_4ohctx * Cecf_4ohctx
    d/dt(central_cepm) <- cl_4ohctx * Cc_4ohctx - cl_cepm * Cc_cepm -
      q_cepm * (Cc_cepm - peripheral1_cepm / vp_cepm) -
      fu_cepm * clin_cepm * Cc_cepm + clef_cepm * Cecf_cepm
    d/dt(peripheral1_cepm) <- q_cepm * (Cc_cepm - peripheral1_cepm / vp_cepm)
    d/dt(ecf_cepm) <- fu_cepm * clin_cepm * Cc_cepm - clef_cepm * Cecf_cepm

    # ---- Drug effect (Eqs. 2-4): transit4 is the removal rate M4 (1/day) ----
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)
    tau <- exp(ltau)
    Cpos_4ohctx <- max(Cc_4ohctx, 0)
    edrug <- emax * Cpos_4ohctx^hill / (ec50^hill + Cpos_4ohctx^hill)
    d/dt(transit1) <- (edrug - transit1) / tau
    d/dt(transit2) <- (transit1 - transit2) / tau
    d/dt(transit3) <- (transit2 - transit3) / tau
    d/dt(transit4) <- (transit3 - transit4) / tau

    # ---- T-cell kinetics (Eqs. 5-11) ----
    k_blood_spleen <- exp(lk_blood_spleen) * 24
    k_spleen_blood <- exp(lk_spleen_blood) * 24
    k_blood_lnode <- exp(lk_blood_lnode) * 24
    k_lnode_blood <- exp(lk_lnode_blood) * 24
    k_blood_liver <- exp(lk_blood_liver) * 24
    k_liver_blood <- exp(lk_liver_blood) * 24
    k_blood_gut <- exp(lk_blood_gut) * 24
    k_gut_blood <- exp(lk_gut_blood) * 24
    k_blood_lung <- exp(lk_blood_lung) * 24
    k_lung_blood <- exp(lk_lung_blood) * 24
    kdeath_tc <- exp(lkdeath_tc)
    kdeath_th <- exp(lkdeath_th)
    kdeath_treg <- exp(lkdeath_treg)
    kprol_tc <- exp(lkprol_tc)
    kprol_th <- exp(lkprol_th)
    kprol_treg <- exp(lkprol_treg)
    h50_tc <- exp(lh50_tc)
    h50_th <- exp(lh50_th)
    h50_treg <- exp(lh50_treg)
    lambda_b <- exp(llambda_b)
    lambda_r <- exp(llambda_r)
    lambda0_tc <- exp(llambda0_tc)
    lambda0_th <- lambda0_tc
    lambda0_treg <- exp(llambda0_treg)

    # Thymic output (Eq. 10); t is days since transplant
    thymus <- (1 - exp(-2 * t / lambda_b)) / (1 + exp(lambda_r * (1 - t / lambda_b)))
    lambda_tc <- lambda0_tc * thymus
    lambda_th <- lambda0_th * thymus
    lambda_treg <- lambda0_treg * thymus

    # Saturable homing into spleen, lymph node and liver (Eq. 11); the spleen
    # and lymph-node amounts include both residence sub-compartments.
    spleen_tc <- tc_spleen1 + tc_spleen2
    spleen_th <- th_spleen1 + th_spleen2
    spleen_treg <- treg_spleen1 + treg_spleen2
    lnode_tc <- tc_lnode1 + tc_lnode2
    lnode_th <- th_lnode1 + th_lnode2
    lnode_treg <- treg_lnode1 + treg_lnode2
    i_spleen_tc <- 1 - spleen_tc^hill_cap / (cap_spleen_tc^hill_cap + spleen_tc^hill_cap)
    i_spleen_th <- 1 - spleen_th^hill_cap / (cap_spleen_th^hill_cap + spleen_th^hill_cap)
    i_spleen_treg <- 1 - spleen_treg^hill_cap / (cap_spleen_treg^hill_cap + spleen_treg^hill_cap)
    i_lnode_tc <- 1 - lnode_tc^hill_cap / (cap_lnode_tc^hill_cap + lnode_tc^hill_cap)
    i_lnode_th <- 1 - lnode_th^hill_cap / (cap_lnode_th^hill_cap + lnode_th^hill_cap)
    i_lnode_treg <- 1 - lnode_treg^hill_cap / (cap_lnode_treg^hill_cap + lnode_treg^hill_cap)
    i_liver_tc <- 1 - tc_liver^hill_cap / (cap_liver_tc^hill_cap + tc_liver^hill_cap)
    i_liver_th <- 1 - th_liver^hill_cap / (cap_liver_th^hill_cap + th_liver^hill_cap)
    i_liver_treg <- 1 - treg_liver^hill_cap / (cap_liver_treg^hill_cap + treg_liver^hill_cap)

    # Cytotoxic (CD8+) T cells. Proliferation r / (1 + x / h) uses each
    # compartment's own amount; the second residence sub-compartments of
    # spleen, lymph node and intestine have neither proliferation nor drug
    # effect (Eq. 8); death acts in blood only (Eq. 5).
    d/dt(tc_blood) <- lambda_tc +
      k_spleen_blood * tc_spleen2 - k_blood_spleen * i_spleen_tc * tc_blood +
      k_lnode_blood * tc_lnode2 - k_blood_lnode * i_lnode_tc * tc_blood +
      k_liver_blood * tc_liver - k_blood_liver * i_liver_tc * tc_blood +
      k_gut_blood * tc_gut2 - k_blood_gut * tc_blood +
      k_lung_blood * tc_lung - k_blood_lung * tc_blood -
      kdeath_tc * tc_blood + kprol_tc / (1 + tc_blood / h50_tc) * tc_blood - transit4 * tc_blood
    d/dt(tc_spleen1) <- k_blood_spleen * i_spleen_tc * tc_blood - k_spleen_blood * tc_spleen1 +
      kprol_tc / (1 + tc_spleen1 / h50_tc) * tc_spleen1 - transit4 * tc_spleen1
    d/dt(tc_spleen2) <- k_spleen_blood * (tc_spleen1 - tc_spleen2)
    d/dt(tc_lnode1) <- k_blood_lnode * i_lnode_tc * tc_blood - k_lnode_blood * tc_lnode1 +
      kprol_tc / (1 + tc_lnode1 / h50_tc) * tc_lnode1 - transit4 * tc_lnode1
    d/dt(tc_lnode2) <- k_lnode_blood * (tc_lnode1 - tc_lnode2)
    d/dt(tc_liver) <- k_blood_liver * i_liver_tc * tc_blood - k_liver_blood * tc_liver +
      kprol_tc / (1 + tc_liver / h50_tc) * tc_liver - transit4 * tc_liver
    d/dt(tc_gut1) <- k_blood_gut * tc_blood - k_gut_blood * tc_gut1 +
      kprol_tc / (1 + tc_gut1 / h50_tc) * tc_gut1 - transit4 * tc_gut1
    d/dt(tc_gut2) <- k_gut_blood * (tc_gut1 - tc_gut2)
    d/dt(tc_lung) <- k_blood_lung * tc_blood - k_lung_blood * tc_lung +
      kprol_tc / (1 + tc_lung / h50_tc) * tc_lung - transit4 * tc_lung

    # Helper (CD4+Foxp3-) T cells
    d/dt(th_blood) <- lambda_th +
      k_spleen_blood * th_spleen2 - k_blood_spleen * i_spleen_th * th_blood +
      k_lnode_blood * th_lnode2 - k_blood_lnode * i_lnode_th * th_blood +
      k_liver_blood * th_liver - k_blood_liver * i_liver_th * th_blood +
      k_gut_blood * th_gut2 - k_blood_gut * th_blood +
      k_lung_blood * th_lung - k_blood_lung * th_blood -
      kdeath_th * th_blood + kprol_th / (1 + th_blood / h50_th) * th_blood - transit4 * th_blood
    d/dt(th_spleen1) <- k_blood_spleen * i_spleen_th * th_blood - k_spleen_blood * th_spleen1 +
      kprol_th / (1 + th_spleen1 / h50_th) * th_spleen1 - transit4 * th_spleen1
    d/dt(th_spleen2) <- k_spleen_blood * (th_spleen1 - th_spleen2)
    d/dt(th_lnode1) <- k_blood_lnode * i_lnode_th * th_blood - k_lnode_blood * th_lnode1 +
      kprol_th / (1 + th_lnode1 / h50_th) * th_lnode1 - transit4 * th_lnode1
    d/dt(th_lnode2) <- k_lnode_blood * (th_lnode1 - th_lnode2)
    d/dt(th_liver) <- k_blood_liver * i_liver_th * th_blood - k_liver_blood * th_liver +
      kprol_th / (1 + th_liver / h50_th) * th_liver - transit4 * th_liver
    d/dt(th_gut1) <- k_blood_gut * th_blood - k_gut_blood * th_gut1 +
      kprol_th / (1 + th_gut1 / h50_th) * th_gut1 - transit4 * th_gut1
    d/dt(th_gut2) <- k_gut_blood * (th_gut1 - th_gut2)
    d/dt(th_lung) <- k_blood_lung * th_blood - k_lung_blood * th_lung +
      kprol_th / (1 + th_lung / h50_th) * th_lung - transit4 * th_lung

    # Regulatory (CD4+CD25+Foxp3+) T cells
    d/dt(treg_blood) <- lambda_treg +
      k_spleen_blood * treg_spleen2 - k_blood_spleen * i_spleen_treg * treg_blood +
      k_lnode_blood * treg_lnode2 - k_blood_lnode * i_lnode_treg * treg_blood +
      k_liver_blood * treg_liver - k_blood_liver * i_liver_treg * treg_blood +
      k_gut_blood * treg_gut2 - k_blood_gut * treg_blood +
      k_lung_blood * treg_lung - k_blood_lung * treg_blood -
      kdeath_treg * treg_blood + kprol_treg / (1 + treg_blood / h50_treg) * treg_blood - transit4 * treg_blood
    d/dt(treg_spleen1) <- k_blood_spleen * i_spleen_treg * treg_blood - k_spleen_blood * treg_spleen1 +
      kprol_treg / (1 + treg_spleen1 / h50_treg) * treg_spleen1 - transit4 * treg_spleen1
    d/dt(treg_spleen2) <- k_spleen_blood * (treg_spleen1 - treg_spleen2)
    d/dt(treg_lnode1) <- k_blood_lnode * i_lnode_treg * treg_blood - k_lnode_blood * treg_lnode1 +
      kprol_treg / (1 + treg_lnode1 / h50_treg) * treg_lnode1 - transit4 * treg_lnode1
    d/dt(treg_lnode2) <- k_lnode_blood * (treg_lnode1 - treg_lnode2)
    d/dt(treg_liver) <- k_blood_liver * i_liver_treg * treg_blood - k_liver_blood * treg_liver +
      kprol_treg / (1 + treg_liver / h50_treg) * treg_liver - transit4 * treg_liver
    d/dt(treg_gut1) <- k_blood_gut * treg_blood - k_gut_blood * treg_gut1 +
      kprol_treg / (1 + treg_gut1 / h50_treg) * treg_gut1 - transit4 * treg_gut1
    d/dt(treg_gut2) <- k_gut_blood * (treg_gut1 - treg_gut2)
    d/dt(treg_lung) <- k_blood_lung * treg_blood - k_lung_blood * treg_lung +
      kprol_treg / (1 + treg_lung / h50_treg) * treg_lung - transit4 * treg_lung

    # Infused donor T cells start in blood; every other state starts at 0
    graft_tcells <- graft_splenocytes * frac_tcell_graft
    tc_blood(0) <- graft_tcells / (1 + ratio_cd4_cd8)
    th_blood(0) <- (1 - frac_treg_cd4) * ratio_cd4_cd8 / (1 + ratio_cd4_cd8) * graft_tcells
    treg_blood(0) <- frac_treg_cd4 * ratio_cd4_cd8 / (1 + ratio_cd4_cd8) * graft_tcells

    # ---- Observations: blood in cells/uL (Eq. 12), tissues in cells ----
    vblood <- WT * 1000 * vblood_per_wt
    tc_blood_conc <- tc_blood / vblood
    th_blood_conc <- th_blood / vblood
    treg_blood_conc <- treg_blood / vblood

    tc_blood_conc ~ prop(propSd_tc_blood_conc)
    th_blood_conc ~ prop(propSd_th_blood_conc)
    treg_blood_conc ~ prop(propSd_treg_blood_conc)
    spleen_tc ~ prop(propSd_spleen_tc)
    spleen_th ~ prop(propSd_spleen_th)
    spleen_treg ~ prop(propSd_spleen_treg)
    lnode_tc ~ prop(propSd_lnode_tc)
    lnode_th ~ prop(propSd_lnode_th)
    lnode_treg ~ prop(propSd_lnode_treg)
    tc_liver ~ prop(propSd_tc_liver)
    th_liver ~ prop(propSd_th_liver)
    treg_liver ~ prop(propSd_treg_liver)
  })
}
