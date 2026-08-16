vandenMaagdenberg_2025_a2ar_tme_qsp <- function() {
  description <- paste(
    "QSP. Preclinical (mouse). Tumour-microenvironment quantitative systems",
    "pharmacology model coupling adenosine / A2A-receptor (A2AR) signalling to",
    "CD8 T-cell dynamics, PD-L1 checkpoint expression and logistic tumour",
    "growth in syngeneic mice, driven by a de novo generated A2AR inhibitor",
    "given intravenously and an anti-PD-L1 monoclonal antibody given",
    "intraperitoneally. Adapted by van den Maagdenberg 2025 from the AZD4635",
    "systems model of Voronova 2021: the original two-compartment oral AZD4635",
    "PK is replaced by a one-compartment IV model for a generated small",
    "molecule, whose central volume (vc_a2ari), elimination rate (kel_a2ari),",
    "A2AR binding affinity (kd_a2ari) and molecular weight (mw_a2ari) are the",
    "four compound-specific inputs the generative-design workflow supplies.",
    "Seven ODE states plus tumour volume: intratumoral adenosine sets A2AR",
    "occupancy, which suppresses the immune activation rate; the A2AR inhibitor",
    "competes with adenosine for the receptor; the anti-PD-L1 antibody frees",
    "PD-L1-bound signalling; the resulting immune activation rate drives naive",
    "CD8 T-cell influx, proliferation and differentiation into cytotoxic",
    "T-lymphocytes, which kill tumour cells. Study-level covariates are fixed",
    "to the Voronova 2021 MCA205-2 syngeneic-mouse study, the arm the authors",
    "selected because its low intratumoral adenosine makes it the most",
    "A2AR-inhibitor-sensitive. Deterministic mechanism model: the authors ran",
    "every simulation with the two structural random effects (on sl and sr)",
    "set to zero and fitted no residual error, so no etas and no error model",
    "are encoded. Default compound-specific values are the",
    "low-elimination-rate / high-affinity example molecule of Figure 8B.",
    sep = " "
  )
  reference <- paste(
    "van den Maagdenberg HW, de Mol van Otterloo J, van Hasselt JGC,",
    "van der Graaf PH, van Westen GJP. Integrating Pharmacokinetics and",
    "Quantitative Systems Pharmacology Approaches in Generative Drug Design.",
    "J Chem Inf Model. 2025;65(10):4783-4796. doi:10.1021/acs.jcim.5c00107.",
    "Model equations Supporting Information Figure S9; variable definitions",
    "Table S4; parameter values Table S5. Systems structure and all system",
    "parameter values adapted from Voronova V, Peskov K, Kosinsky Y, et al.",
    "Evaluation of Combination Strategies for the A2AR Inhibitor AZD4635",
    "Across Tumor Microenvironment Conditions via a Systems Pharmacology",
    "Model. Front Immunol. 2021;12:617316. doi:10.3389/fimmu.2021.617316.",
    "Full-precision parameter values and the four Figure 8B example-compound",
    "property sets are taken from the authors' own deposited implementation",
    "(https://github.com/CDDLeiden/PK-in-generative-drug-design,",
    "05_QSP_modelling/models/TME_model_1cmp.R) and data archive",
    "(doi:10.5281/zenodo.15082627); see the vignette source-trace table.",
    sep = " "
  )
  vignette <- "vandenMaagdenberg_2025_a2ar_tme_qsp"

  # Mechanistic tumour-microenvironment states. `tumor` is canonical; the
  # remaining seven are paper-mechanistic and do not generalise, so they are
  # declared here rather than added to the compartment register.
  paper_specific_compartments <- c(
    "depot_mab", "central_mab", "central_a2ari",
    "tcell_naive", "tcell_ctl", "pdl1", "adenosine"
  )

  units <- list(
    time   = "day",
    dosing = paste(
      "mg for both agents. The anti-PD-L1 antibody is dosed into depot_mab",
      "(intraperitoneal); the A2AR inhibitor is dosed into central_a2ari",
      "(intravenous bolus). The paper's regimen for a 0.025 kg mouse is",
      "0.125 mg antibody (5 mg/kg) every 3.5 days and 1.25 mg A2AR inhibitor",
      "(50 mg/kg) every 0.5 days, both starting on day 7.",
      sep = " "
    ),
    concentration = paste(
      "nM for both drug observables (Cc_mab, Cc_a2ari) and for intratumoral",
      "adenosine. Tumour volume (tumor) is in uL, which equals mm^3 -- the",
      "unit the paper's figures use. T-cell states are cell counts and pdl1",
      "is a normalised PD-L1 level.",
      sep = " "
    )
  )

  compartmentData <- list(
    depot_mab     = list(analyte = "anti-PD-L1 monoclonal antibody", units = "mg", specimen = "administration site", verified = TRUE),
    central_mab   = list(analyte = "anti-PD-L1 monoclonal antibody", units = "nmol", specimen = "plasma", verified = TRUE),
    central_a2ari = list(analyte = "A2A receptor inhibitor", units = "mg", specimen = "plasma", verified = TRUE),
    tumor         = list(analyte = "tumor", units = "uL", specimen = "not applicable", verified = TRUE),
    tcell_naive   = list(analyte = "naive CD8 T cell", units = "cells", specimen = "not applicable", verified = TRUE),
    tcell_ctl     = list(analyte = "cytotoxic CD8 T-lymphocyte", units = "cells", specimen = "not applicable", verified = TRUE),
    pdl1          = list(analyte = "PD-L1", units = "nM", specimen = "not applicable", verified = TRUE),
    adenosine     = list(analyte = "adenosine", units = "nM", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    SYNGENEIC_MODEL = list(
      description = paste(
        "Syngeneic mouse tumour model identity. Voronova 2021 estimated a",
        "separate covariate set (on sl, sr, vado and tvin) for each of its four",
        "studies (CT26, MC38, MCA205 study 1, MCA205 study 2). van den",
        "Maagdenberg 2025 selected the MCA205-2 set and used it throughout, so",
        "this model hard-codes those four coefficients as fixed parameters",
        "rather than carrying a switch covariate.",
        sep = " "
      ),
      units              = NA_character_,
      type               = "categorical",
      reference_category = "MCA205-2",
      notes              = paste(
        "Not retained as a covariate column: the source paper reports only the",
        "MCA205-2 coefficients (Methods, 'Quantitative Systems Pharmacology",
        "Model'). The other three studies' coefficients are in Voronova 2021",
        "and would have to come from that paper, not this one.",
        sep = " "
      )
    )
  )

  population <- list(
    species    = "mouse (MCA205 syngeneic tumour model)",
    n_subjects = NA_integer_,
    n_studies  = 1L,
    disease_state = "subcutaneous MCA205 syngeneic tumour",
    dose_range = paste(
      "Anti-PD-L1 monoclonal antibody 5 mg/kg intraperitoneally twice weekly",
      "(five doses) plus A2AR inhibitor 50 mg/kg intravenously twice daily",
      "(31 doses); both start on day 7 after tumour inoculation and run to",
      "day 22, with the system simulated for 30 days.",
      sep = " "
    ),
    notes = paste(
      "The systems parameters were estimated by Voronova 2021 against pooled",
      "individual tumour-size dynamics from four syngeneic-mouse studies (116",
      "mice in three efficacy studies plus 40 mice in an MCA205 dose-finding",
      "study). van den Maagdenberg 2025 did not re-estimate any system",
      "parameter; it re-used the Voronova MCA205 study-2 covariate set and",
      "substituted QSPR-predicted, allometrically scaled PK and A2AR affinity",
      "for de novo generated molecules in place of AZD4635. Human QSPR",
      "predictions were scaled to the 0.025 kg mouse with allometric exponents",
      "of 0.65 for clearance and 0.95 for volume (main text eq 1), and the",
      "mouse elimination rate is then clearance divided by volume.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Compound-specific inputs. Table S5 lists these four as "-": they are
    # supplied per generated molecule, not estimated. Defaults below are the
    # "low elimination rate, high affinity" example compound of Figure 8B,
    # whose printed legend values (kel 5.9 1/day, Ki 0.9 nM) reproduce from
    # the authors' deposited generated-compound table. Override all four via
    # rxSolve(params=) to simulate a different molecule.
    # ------------------------------------------------------------------
    vc_a2ari <- fixed(0.2242577481014); label("Central volume of distribution of the A2AR inhibitor in the mouse (L)")  # Table S5 VcARinh (reported as "-"); Zenodo 15082627 A2AR_VDSSmax_0 generated compound, scaled by eq 1 with exponent 0.95
    kel_a2ari <- fixed(5.90479524688); label("Elimination rate constant of the A2AR inhibitor (1/day)")  # Table S5 kelARinh (reported as "-"); Figure 8B legend "kel: 5.9 1/day"
    kd_a2ari <- fixed(0.947208676066); label("A2AR binding affinity (Ki) of the A2AR inhibitor (nM)")  # Table S5 KdARinh (reported as "-"); Figure 8B legend "Ki: 0.9 nM"
    mw_a2ari <- fixed(454.484985352); label("Molecular weight of the A2AR inhibitor (g/mol)")  # Table S5 MW (reported as "-"); RDKit MolWt of the same Figure 8B compound

    # ------------------------------------------------------------------
    # Anti-PD-L1 monoclonal antibody PK and binding (Table S5)
    # ------------------------------------------------------------------
    ka_mab <- fixed(8); label("Intraperitoneal absorption rate constant of the anti-PD-L1 antibody (1/day)")  # Table S5 kainput1 = 8 1/d
    convf_mab <- fixed(6.66); label("Unit conversion for the antibody dose, mg to nmol (nmol/mg)")  # Table S5 convF1 = 6.66 (1/0.15; MW 150 kDa)
    vc_mab <- fixed(0.003); label("Central volume of distribution of the anti-PD-L1 antibody (L)")  # Table S5 Vc = 0.003 L
    kel_mab <- fixed(0.1); label("Elimination rate constant of the anti-PD-L1 antibody (1/day)")  # Table S5 kelmAb = 0.1 1/d
    kd_mab <- fixed(30); label("Anti-PD-L1 antibody / PD-L1 binding affinity (nM)")  # Table S5 Kd1 = 30 nM

    # ------------------------------------------------------------------
    # Adenosine / A2AR receptor module (Table S5)
    # ------------------------------------------------------------------
    kd_ado <- fixed(1182); label("Binding affinity of adenosine for the A2A receptor (nM)")  # Table S5 Kdado = 1182 nM
    scf_ado <- fixed(1); label("Scaling factor from total intratumoral to extracellular adenosine (unitless)")  # Table S5 sCf = 1
    tvado <- fixed(100000); label("Typical adenosine level in the tumour before the study covariate (nM)")  # Table S5 Vado = 100000 nM (before exp(betaMCA205))
    kado_tumor <- fixed(80); label("Tumour volume at half-maximal adenosine accumulation (uL)")  # Table S5 Kado = 80
    ic50_a2ar <- fixed(1.82315256991075); label("A2AR occupancy giving a 50% decrease in immune cell activity (unitless)")  # Table S5 IC50 = 1.823
    vmax_supr <- fixed(0.7); label("Maximal adenosine-mediated suppression of effector T cells and APCs (unitless)")  # Table S5 Vmaxsupr = 0.7

    # ------------------------------------------------------------------
    # T-cell and PD-L1 dynamics (Table S5)
    # ------------------------------------------------------------------
    kln <- fixed(209.674023531827); label("Maximal influx rate of naive effector T cells into the tumour (cells/day)")  # Table S5 kLn = 209.674 cells/d
    tsl <- fixed(4.55575391256157); label("Typical T-cell tumour-infiltration constant under antigen exposure before the study covariate (uL/day)")  # Table S5 sL = 4.556 uL/d (before exp(betaCIV151))
    tsr <- fixed(57.0519720073586); label("Typical sensitivity of cellular immunosuppression to systemic antigen before the study covariate (uL/day)")  # Table S5 sR = 57.052 uL/d (before exp(betaMCA205))
    kel_tn <- fixed(0.2); label("Loss rate constant of naive effector T cells (1/day)")  # Table S5 kel = 0.2 1/d (renamed from the paper's bare `kel` to avoid the canonical PK elimination-rate name)
    kapo_ctl <- fixed(2.0); label("Apoptosis rate constant of cytotoxic T-lymphocytes (1/day)")  # Table S5 kapo = 2.0 1/d
    kpro_tn <- fixed(3); label("Maximal T-cell proliferation rate (1/day)")  # Table S5 kpro = 3 1/d
    kdif_tn <- fixed(3.2); label("Maximal T-cell differentiation rate (1/day)")  # Table S5 kdif = 3.2 1/d
    kp_pdl1 <- fixed(1279.91714283438); label("Effector T-cell count giving half-maximal PD-L1 up-regulation (cells)")  # Table S5 Kp = 1279.917 cells

    # ------------------------------------------------------------------
    # Tumour growth and kill (Table S5)
    # ------------------------------------------------------------------
    tvmax_tumor <- fixed(3500); label("Maximal (carrying-capacity) tumour volume (uL)")  # Table S5 TVmax = 3500 uL
    r_tumor <- fixed(0.521743095907443); label("Logistic tumour growth rate constant (1/day)")  # Table S5 r = 0.522 1/d
    beff <- fixed(0.001); label("Rate of tumour cell kill per cytotoxic T-lymphocyte (1/(day*cell))")  # Table S5 beff = 0.001 1/(d*cell)
    kdeath_spont <- fixed(0.01); label("Slow spontaneous tumour cell death rate constant (1/day)")  # Table S5 d = 0.01 1/d (renamed from the paper's bare `d` to avoid colliding with the d/dt() operator)
    ttvin <- fixed(2.05452388904701); label("Typical initial tumour volume before the study covariate (uL)")  # Table S5 TVin = 2.055 uL (before exp(betaMCA205))

    # ------------------------------------------------------------------
    # Study covariate coefficients, fixed to the Voronova 2021 MCA205-2
    # syngeneic-mouse study. Values from the paper's Methods, which are more
    # precise than the rounded exponents printed in Table S5.
    # ------------------------------------------------------------------
    sl_cov <- fixed(0); label("MCA205-2 study coefficient on sl (log scale)")  # Methods: "sLcov = 0"; Table S5 sL * exp(0)
    sr_cov <- fixed(0.5308); label("MCA205-2 study coefficient on sr (log scale)")  # Methods: "sRcov = 0.5308"; Table S5 prints the rounded exp(0.531)
    vado_cov <- fixed(-3); label("MCA205-2 study coefficient on vado (log scale)")  # Methods: "Vadocov = -3"; Table S5 Vado * exp(-3)
    tvin_cov <- fixed(0.69); label("MCA205-2 study coefficient on the initial tumour volume (log scale)")  # Methods: "TVincov = 0.69"; Table S5 TVin * exp(0.69)
  })

  model({
    # ------------------------------------------------------------------
    # 1. Study-covariate-adjusted system parameters.
    #    The paper writes these as sL = TsL * exp(eta.sL + sL_cov) etc.;
    #    every simulation in the paper sets eta.sL = eta.sR = 0, so the etas
    #    are not encoded here (see description and vignette Errata).
    # ------------------------------------------------------------------
    sl <- tsl * exp(sl_cov)
    sr <- tsr * exp(sr_cov)
    vado <- tvado * exp(vado_cov)
    tvin <- ttvin * exp(tvin_cov)

    # ------------------------------------------------------------------
    # 2. Drug disposition.
    #    Anti-PD-L1 antibody: first-order intraperitoneal absorption into a
    #    one-compartment central pool, with the mg -> nmol conversion applied
    #    on transfer (SI Figure S9 eqs 1a-1b).
    #    A2AR inhibitor: one-compartment IV bolus (SI Figure S9 eq 1d).
    # ------------------------------------------------------------------
    d/dt(depot_mab) <- -ka_mab * depot_mab
    d/dt(central_mab) <- ka_mab * depot_mab * convf_mab - kel_mab * central_mab
    d/dt(central_a2ari) <- -kel_a2ari * central_a2ari

    # nmol / L = nM
    Cc_mab <- central_mab / vc_mab                                # SI Figure S9 eq 1c
    # (mg / L) / (g/mol) * 1e6 = nM
    Cc_a2ari <- central_a2ari / vc_a2ari / mw_a2ari * 1e6         # SI Figure S9 eq 1e

    # ------------------------------------------------------------------
    # 3. Receptor occupancy and immune suppression.
    # ------------------------------------------------------------------
    pdl1free <- pdl1 / (1 + Cc_mab / kd_mab)                      # SI Figure S9 eq 1f

    # Competitive occupancy of the A2A receptor by adenosine, with the
    # inhibitor competing for the same site.
    a2ar_occup <-
      (adenosine * scf_ado / kd_ado) /
      (1 + adenosine * scf_ado / kd_ado + Cc_a2ari / kd_a2ari)    # SI Figure S9 eq 1g

    ado_suppr <- vmax_supr * a2ar_occup / (a2ar_occup + ic50_a2ar)  # SI Figure S9 eq 1h

    # ------------------------------------------------------------------
    # 4. Antigen release and immune activation.
    # ------------------------------------------------------------------
    tkr <- beff * tcell_ctl + kdeath_spont                        # SI Figure S9 eq 1i
    ag <- tkr * tumor * (1 - ado_suppr)                           # SI Figure S9 eq 1j
    prfunc <- (1 - pdl1free) * (1 - ag / (ag + sr)) * (1 - ado_suppr)  # SI Figure S9 eq 1k
    tninf <- kln * ag / (ag + sl)                                 # SI Figure S9 eq 1l
    cd8tot <- tcell_ctl + tcell_naive                             # SI Figure S9 eq 1m
    isc <- ag / (ag + sr)                                         # SI Figure S9 eq 1n

    # ------------------------------------------------------------------
    # 5. Growth-arrest switch. Reproduces the source implementation exactly:
    #    tumour growth is switched off once treatment has begun (day 7) and
    #    the tumour has fallen below 10 uL. Table S5 lists xf as "-" and
    #    describes it in words; the conditional is from the authors'
    #    TME_model_1cmp.R.
    # ------------------------------------------------------------------
    if (time > 7 & tumor < 10) {
      xf <- 0
    } else {
      xf <- 1
    }

    # ------------------------------------------------------------------
    # 6. Initial condition. Only the tumour is seeded; the T-cell, PD-L1 and
    #    adenosine states start at zero and are driven up by the system, as
    #    in the source implementation.
    # ------------------------------------------------------------------
    tumor(0) <- tvin

    # ------------------------------------------------------------------
    # 7. Model reactions and ODE system (SI Figure S9 eqs 1o-1x).
    # ------------------------------------------------------------------
    tum_gr <- tumor * r_tumor * (1 - tumor / tvmax_tumor) * xf    # eq 1o
    tum_kill <- (beff * tcell_ctl + kdeath_spont) * tumor         # eq 1p
    ctl_dynamic1 <- tninf + kpro_tn * prfunc * tcell_naive - kel_tn * tcell_naive  # eq 1q
    ctl_dynamic2 <- kdif_tn * prfunc * tcell_naive                # eq 1r
    ctl_dynamic5 <- kapo_ctl * tcell_ctl                          # eq 1s

    d/dt(tumor) <- tum_gr - tum_kill                              # eq 1t
    d/dt(tcell_naive) <- ctl_dynamic1 - ctl_dynamic2              # eq 1u
    d/dt(tcell_ctl) <- ctl_dynamic2 - ctl_dynamic5                # eq 1v
    d/dt(pdl1) <- tcell_ctl / (tcell_ctl + kp_pdl1) - pdl1        # eq 1w
    d/dt(adenosine) <- vado * tumor / (tumor + kado_tumor) - adenosine  # eq 1x
  })
}
