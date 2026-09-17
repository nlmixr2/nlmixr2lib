Kletting_2015_antiCD66_pbpk_model1 <- function() {
  description <- paste(
    "PBPK (whole-body, 86 ODE states, SAAM II). Physiologically based",
    "pharmacokinetic model for the murine anti-CD66 monoclonal antibody",
    "BW 250/183, DTPA-conjugated and radiolabelled with 111In for",
    "pre-therapeutic imaging or with 90Y for radioimmunotherapy of acute",
    "leukaemia before stem-cell transplantation. Six parallel circulations are",
    "carried: the antibody is split into fully, half and non-immunoreactive",
    "species by the probability r_im that one antibody arm is immunoreactive,",
    "and each of those three is carried as a radiolabelled and an unlabelled",
    "copy, coupled only by physical decay, which converts labelled species into",
    "their unlabelled counterparts in every compartment. Fully immunoreactive",
    "antibody binds CD66 both monovalently and bivalently, half immunoreactive",
    "antibody only monovalently, and all species compete for the same finite",
    "antigen pool in red marrow, liver, spleen and blood. Antigen numbers in red",
    "marrow and blood are not estimated but constrained to the fitted liver and",
    "spleen numbers through granulocyte-pool ratios. Model 1 constrains the blood antigen number to the unweighted liver plus spleen sum (S1 Text Eq 8).",
    "Liver, spleen, red marrow, GI tract, an interstitial space and the main",
    "vascular compartment are perfused by plasma flow; degraded antibody follows",
    "an Eger/Houston four-state submodel. Nine parameters were fitted per patient",
    "to gamma-camera and serum data from 27 patients with acute leukaemia; the",
    "shipped values are the cohort means of Table 1 for Model 1. Amounts",
    "are nmol and time is minutes.",
    sep = " "
  )
  reference <- paste(
    "Kletting P, Maass C, Reske S, Beer AJ, Glatting G. Physiologically Based",
    "Pharmacokinetic Modeling Is Essential in 90Y-Labeled Anti-CD66",
    "Radioimmunotherapy. PLoS One. 2015;10(5):e0127934.",
    "doi:10.1371/journal.pone.0127934. Model equations, fixed parameters and",
    "data-assignment equations are in supplement S1 Text (Eqs 1-43, Tables A",
    "and B); per-patient administered amounts and fitting results are in",
    "supplement S1 Table.",
    sep = " "
  )
  vignette <- "Kletting_2015_antiCD66_pbpk"
  units <- list(time = "min", dosing = "nmol", concentration = "nmol/L")

  paper_specific_compartment_pattern <- "^(ab|agbi|agmono|ex|metap|metaex1|metaex2|cleared)_"

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Sets the number of CD66-expressing cells in the red marrow,",
        "N_cells,RM = 188e8 per kg * BW (S1 Text Table B), from which the liver,",
        "spleen and blood cell numbers and hence every bivalent enhancement factor",
        "alpha_i follow. Measured individually but not published per patient; the",
        "demonstration value is the 70 kg reference adult."
      ),
      source_name = "BW"
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters only the red-marrow readout. The UlmDos height-corrected scaling",
        "factor from the L2-L4 lumbar-spine region of interest to the entire red",
        "marrow is K = (170 / HT) / 0.06665 (S1 Text Table A, used in Eq 20). It",
        "affects no state, only the abMarrowRoi observable. Not published per",
        "patient; the demonstration value is 170 cm, which makes K = 1 / 0.06665."
      ),
      source_name = "height"
    )
  )

  compartmentData <- list(
    ab_fa_plasma_unlab = list(analyte = "anti-CD66 antibody", units = "nmol", specimen = "serum", verified = TRUE),
    ab_fa_plasma_lab = list(
      analyte = "radiolabelled anti-CD66 antibody",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    ab_fa_liver_unlab = list(analyte = "anti-CD66 antibody", units = "nmol", specimen = "tissue", verified = TRUE),
    ab_fa_liver_lab = list(
      analyte = "radiolabelled anti-CD66 antibody",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    ab_fa_spleen_unlab = list(analyte = "anti-CD66 antibody", units = "nmol", specimen = "tissue", verified = TRUE),
    ab_fa_spleen_lab = list(
      analyte = "radiolabelled anti-CD66 antibody",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    ab_fa_marrow_unlab = list(analyte = "anti-CD66 antibody", units = "nmol", specimen = "tissue", verified = TRUE),
    ab_fa_marrow_lab = list(
      analyte = "radiolabelled anti-CD66 antibody",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    ab_fa_gi_unlab = list(analyte = "anti-CD66 antibody", units = "nmol", specimen = "tissue", verified = TRUE),
    ab_fa_gi_lab = list(
      analyte = "radiolabelled anti-CD66 antibody",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    ab_fa_int_unlab = list(analyte = "anti-CD66 antibody", units = "nmol", specimen = "tissue", verified = TRUE),
    ab_fa_int_lab = list(
      analyte = "radiolabelled anti-CD66 antibody",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    ab_ha_plasma_unlab = list(analyte = "anti-CD66 antibody", units = "nmol", specimen = "serum", verified = TRUE),
    ab_ha_plasma_lab = list(
      analyte = "radiolabelled anti-CD66 antibody",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    ab_ha_liver_unlab = list(analyte = "anti-CD66 antibody", units = "nmol", specimen = "tissue", verified = TRUE),
    ab_ha_liver_lab = list(
      analyte = "radiolabelled anti-CD66 antibody",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    ab_ha_spleen_unlab = list(analyte = "anti-CD66 antibody", units = "nmol", specimen = "tissue", verified = TRUE),
    ab_ha_spleen_lab = list(
      analyte = "radiolabelled anti-CD66 antibody",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    ab_ha_marrow_unlab = list(analyte = "anti-CD66 antibody", units = "nmol", specimen = "tissue", verified = TRUE),
    ab_ha_marrow_lab = list(
      analyte = "radiolabelled anti-CD66 antibody",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    ab_ha_gi_unlab = list(analyte = "anti-CD66 antibody", units = "nmol", specimen = "tissue", verified = TRUE),
    ab_ha_gi_lab = list(
      analyte = "radiolabelled anti-CD66 antibody",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    ab_ha_int_unlab = list(analyte = "anti-CD66 antibody", units = "nmol", specimen = "tissue", verified = TRUE),
    ab_ha_int_lab = list(
      analyte = "radiolabelled anti-CD66 antibody",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    ab_na_plasma_unlab = list(analyte = "anti-CD66 antibody", units = "nmol", specimen = "serum", verified = TRUE),
    ab_na_plasma_lab = list(
      analyte = "radiolabelled anti-CD66 antibody",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    ab_na_liver_unlab = list(analyte = "anti-CD66 antibody", units = "nmol", specimen = "tissue", verified = TRUE),
    ab_na_liver_lab = list(
      analyte = "radiolabelled anti-CD66 antibody",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    ab_na_spleen_unlab = list(analyte = "anti-CD66 antibody", units = "nmol", specimen = "tissue", verified = TRUE),
    ab_na_spleen_lab = list(
      analyte = "radiolabelled anti-CD66 antibody",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    ab_na_marrow_unlab = list(analyte = "anti-CD66 antibody", units = "nmol", specimen = "tissue", verified = TRUE),
    ab_na_marrow_lab = list(
      analyte = "radiolabelled anti-CD66 antibody",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    ab_na_gi_unlab = list(analyte = "anti-CD66 antibody", units = "nmol", specimen = "tissue", verified = TRUE),
    ab_na_gi_lab = list(
      analyte = "radiolabelled anti-CD66 antibody",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    ab_na_int_unlab = list(analyte = "anti-CD66 antibody", units = "nmol", specimen = "tissue", verified = TRUE),
    ab_na_int_lab = list(
      analyte = "radiolabelled anti-CD66 antibody",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    agbi_fa_plasma_unlab = list(
      analyte = "anti-CD66 antibody bound to two CD66 antigens",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    agbi_fa_plasma_lab = list(
      analyte = "radiolabelled anti-CD66 antibody bound to two CD66 antigens",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    agbi_fa_liver_unlab = list(
      analyte = "anti-CD66 antibody bound to two CD66 antigens",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    agbi_fa_liver_lab = list(
      analyte = "radiolabelled anti-CD66 antibody bound to two CD66 antigens",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    agbi_fa_spleen_unlab = list(
      analyte = "anti-CD66 antibody bound to two CD66 antigens",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    agbi_fa_spleen_lab = list(
      analyte = "radiolabelled anti-CD66 antibody bound to two CD66 antigens",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    agbi_fa_marrow_unlab = list(
      analyte = "anti-CD66 antibody bound to two CD66 antigens",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    agbi_fa_marrow_lab = list(
      analyte = "radiolabelled anti-CD66 antibody bound to two CD66 antigens",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    agmono_fa_plasma_unlab = list(
      analyte = "anti-CD66 antibody bound to one CD66 antigen",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    agmono_fa_plasma_lab = list(
      analyte = "radiolabelled anti-CD66 antibody bound to one CD66 antigen",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    agmono_fa_liver_unlab = list(
      analyte = "anti-CD66 antibody bound to one CD66 antigen",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    agmono_fa_liver_lab = list(
      analyte = "radiolabelled anti-CD66 antibody bound to one CD66 antigen",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    agmono_fa_spleen_unlab = list(
      analyte = "anti-CD66 antibody bound to one CD66 antigen",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    agmono_fa_spleen_lab = list(
      analyte = "radiolabelled anti-CD66 antibody bound to one CD66 antigen",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    agmono_fa_marrow_unlab = list(
      analyte = "anti-CD66 antibody bound to one CD66 antigen",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    agmono_fa_marrow_lab = list(
      analyte = "radiolabelled anti-CD66 antibody bound to one CD66 antigen",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    agmono_ha_plasma_unlab = list(
      analyte = "anti-CD66 antibody bound to one CD66 antigen",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    agmono_ha_plasma_lab = list(
      analyte = "radiolabelled anti-CD66 antibody bound to one CD66 antigen",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    agmono_ha_liver_unlab = list(
      analyte = "anti-CD66 antibody bound to one CD66 antigen",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    agmono_ha_liver_lab = list(
      analyte = "radiolabelled anti-CD66 antibody bound to one CD66 antigen",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    agmono_ha_spleen_unlab = list(
      analyte = "anti-CD66 antibody bound to one CD66 antigen",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    agmono_ha_spleen_lab = list(
      analyte = "radiolabelled anti-CD66 antibody bound to one CD66 antigen",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    agmono_ha_marrow_unlab = list(
      analyte = "anti-CD66 antibody bound to one CD66 antigen",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    agmono_ha_marrow_lab = list(
      analyte = "radiolabelled anti-CD66 antibody bound to one CD66 antigen",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    ex_fa_unlab = list(
      analyte = "degraded anti-CD66 antibody, extravascular delay space",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    metap_fa_unlab = list(
      analyte = "anti-CD66 antibody degradation product",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    metaex1_fa_unlab = list(
      analyte = "anti-CD66 antibody degradation product",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    metaex2_fa_unlab = list(
      analyte = "anti-CD66 antibody degradation product",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    ex_fa_lab = list(
      analyte = "radiolabelled degraded anti-CD66 antibody, extravascular delay space",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    metap_fa_lab = list(
      analyte = "radiolabelled anti-CD66 antibody degradation product",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    metaex1_fa_lab = list(
      analyte = "radiolabelled anti-CD66 antibody degradation product",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    metaex2_fa_lab = list(
      analyte = "radiolabelled anti-CD66 antibody degradation product",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    ex_ha_unlab = list(
      analyte = "degraded anti-CD66 antibody, extravascular delay space",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    metap_ha_unlab = list(
      analyte = "anti-CD66 antibody degradation product",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    metaex1_ha_unlab = list(
      analyte = "anti-CD66 antibody degradation product",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    metaex2_ha_unlab = list(
      analyte = "anti-CD66 antibody degradation product",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    ex_ha_lab = list(
      analyte = "radiolabelled degraded anti-CD66 antibody, extravascular delay space",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    metap_ha_lab = list(
      analyte = "radiolabelled anti-CD66 antibody degradation product",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    metaex1_ha_lab = list(
      analyte = "radiolabelled anti-CD66 antibody degradation product",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    metaex2_ha_lab = list(
      analyte = "radiolabelled anti-CD66 antibody degradation product",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    ex_na_unlab = list(
      analyte = "degraded anti-CD66 antibody, extravascular delay space",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    metap_na_unlab = list(
      analyte = "anti-CD66 antibody degradation product",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    metaex1_na_unlab = list(
      analyte = "anti-CD66 antibody degradation product",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    metaex2_na_unlab = list(
      analyte = "anti-CD66 antibody degradation product",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    ex_na_lab = list(
      analyte = "radiolabelled degraded anti-CD66 antibody, extravascular delay space",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    metap_na_lab = list(
      analyte = "radiolabelled anti-CD66 antibody degradation product",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    metaex1_na_lab = list(
      analyte = "radiolabelled anti-CD66 antibody degradation product",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    metaex2_na_lab = list(
      analyte = "radiolabelled anti-CD66 antibody degradation product",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    cleared_unlab = list(
      analyte = "antibody-derived material excreted from the body",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    ),
    cleared_lab = list(
      analyte = "radiolabelled antibody-derived material excreted from the body",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 27,
    n_studies = 2,
    disease_state = "acute leukaemia (21 acute myeloid, 6 acute lymphoblastic); radioimmunotherapy to intensify conditioning before stem-cell transplantation",
    dose_range = paste(
      "Pre-therapeutic imaging: 0.5 +/- 0.1 mg anti-CD66 antibody",
      "(1 mg = 6.7 nmol; 3.3 +/- 0.6 nmol in S1 Table) carrying a mean 111In",
      "activity of 130 +/- 16 MBq, intravenous bolus. Therapy: 1.3 +/- 0.5 mg",
      "(8.7 +/- 3.1 nmol) carrying 3.2 +/- 0.9 GBq of 90Y, intravenous bolus about",
      "eight days later."
    ),
    regions = "Germany (Ulm University)",
    notes = paste(
      "Both study protocols were approved by the Ethics Committee of Ulm",
      "University. Age, sex, weight and height distributions are not reported.",
      "Each patient was fitted individually (27 separate fits, SAAM II 2.2)",
      "rather than by a population model, so there is no between-subject random",
      "effect to extract: Table 1 reports the mean and SD of the 27 individual",
      "estimates and the means are shipped here as the typical values. Data were",
      "weighted with a fractional standard deviation of 0.1 (S1 Table), and the",
      "total serum volume additionally carried a Bayesian prior with mean V_m,1",
      "from the 5 min serum sample and SD (V_m,1 - V_BSA). The between-patient",
      "SDs of the individual estimates for Model 1 were AgL 0.26, AgS 0.19, exL 0.089, exS 0.048, fRM 0.21 %, cRM 0.33, lambda_db 1.7e-5 /min, r_im 0.090 and V 0.62 L;",
      "they are estimation-error-inflated individual spreads, not a fitted omega,",
      "and are recorded here rather than encoded as inter-individual variability."
    )
  )

  ini({
    # ---- Fitted parameters, Table 1 cohort means (Model 1) ---------------
    ag0Liver <- fixed(0.31) ; label("Total CD66 antigen number in the liver, Ag0,L (nmol)")      # Table 1
    ag0Spleen <- fixed(0.22) ; label("Total CD66 antigen number in the spleen, Ag0,S (nmol)")    # Table 1
    exLiver <- fixed(0.235) ; label("Liver fraction of the extravascular delay compartment, exL (unitless)")    # Table 1
    exSpleen <- fixed(0.107) ; label("Spleen fraction of the extravascular delay compartment, exS (unitless)") # Table 1
    fRm <- fixed(0.0067) ; label("Fraction of total plasma flow reaching the red marrow, fRM (unitless)")      # Table 1 reports fRM in per cent
    cRm <- fixed(1.22) ; label("Individual correction of the L2-L4 to whole red marrow scaling factor, cRM (unitless)") # Table 1
    lambdaDb <- fixed(6.8e-5) ; label("Degradation rate of bound antibody, lambda_db (1/min)")       # Table 1 reports lambda_db in units of 1e-5 /min
    rIm <- fixed(0.801) ; label("Probability that one antibody arm is immunoreactive, r_im (unitless)")  # Table 1
    vSerum <- fixed(2.99) ; label("Total serum volume, V (L)")                                        # Table 1

    # ---- Binding (S1 Text Table B) ---------------------------------------
    konMono <- fixed(0.006) ; label("Monovalent association rate onto CD66, k_on,mono (L/nmol/min)")      # Table B, ref [3]; fixed to a typical value because sensitivity analysis showed no material effect
    koff <- fixed(0.06) ; label("Dissociation rate, k_off (1/min)")                                       # Table B, ref [12]; fixed to a typical value
    enhancement <- fixed(1672000) ; label("Bivalent enhancement factor, E = k_on,bi / k_on,mono (1/cm)")  # Table B, ref [6] Kaufman and Jain 1992
    rCell <- fixed(6) ; label("Radius of a CD66-expressing cell, r_cell (um)")                            # Table B, ref [15]
    nCellRmPerKg <- fixed(1.88e10) ; label("CD66-expressing cells in the red marrow per kg body weight (1/kg)") # Table B, refs [7,13]: 188e8/kg * BW
    agRatioRmBlood <- fixed(38) ; label("Red marrow to circulating CD66-positive cell ratio (unitless)")  # S1 Text Eq 4, ref [7]
    agBloodOrganFrac <- fixed(0.9) ; label("Liver plus spleen resident CD66-positive cells as a fraction of the circulating pool (unitless)") # S1 Text Eqs 5-7, ref [8]

    # ---- Physical decay (S1 Text Table B) --------------------------------
    # Table B gives 1.72e-4 /min for 111In and 1.80e-4 /min for 90Y. The shipped
    # value is 111In, the pre-therapeutic label these parameters were fitted to;
    # switch it with rxode2::ini(mod, lambdaPhy = 1.80e-4) for the 90Y therapy.
    lambdaPhy <- fixed(1.72e-4) ; label("Physical decay constant of the radiolabel, lambda_phy (1/min)")

    # ---- Distribution, degradation and clearance (S1 Text Table B) -------
    kIn <- fixed(0.0017) ; label("Transport rate from plasma to the interstitial space, k_in (1/min)")    # Table B, ref [1] Eger 1987
    kOut <- fixed(0.005) ; label("Transport rate from the interstitial space to plasma, k_out (1/min)")   # Table B, ref [1]
    lambdaDu <- fixed(3.9e-4) ; label("Degradation rate of unbound antibody, lambda_du (1/min)")          # Table B, ref [1]
    lambdaClex <- fixed(3.9e-5) ; label("Clearance from the extravascular delay space, lambda_clex (1/min)")  # Table B, ref [1]
    lambdaMetaex1 <- fixed(0.39) ; label("Degradation product, vascular to extravascular pool 1 (1/min)")     # Table B, ref [17] Houston 1979
    lambdaMetaex2 <- fixed(0.17) ; label("Degradation product, extravascular pool 1 to vascular (1/min)")     # Table B, ref [17]
    lambdaMetaex3 <- fixed(0.018) ; label("Degradation product, extravascular pool 1 to pool 2 (1/min)")      # Table B, ref [17]
    lambdaMetaex4 <- fixed(0.013) ; label("Degradation product, extravascular pool 2 to pool 1 (1/min)")      # Table B, ref [17]
    lambdaCl <- fixed(0.047) ; label("Clearance of degradation product from the body, lambda_cl (1/min)")     # Table B, ref [17]

    # ---- Volume and flow fractions (S1 Text Eq 1 and Table B) ------------
    fVolLiver <- fixed(0.1) ; label("Antibody distribution volume of the liver as a fraction of V (unitless)")      # S1 Text Eq 1, ref [5] Leggett 1995
    fVolSpleen <- fixed(0.014) ; label("Antibody distribution volume of the spleen as a fraction of V (unitless)")  # S1 Text Eq 1, ref [5]
    fVolGi <- fixed(0.076) ; label("Antibody distribution volume of the GI tract as a fraction of V (unitless)")    # S1 Text Eq 1, ref [5]
    fVolMarrow <- fixed(0.04) ; label("Antibody distribution volume of the red marrow as a fraction of V (unitless)") # S1 Text Eq 1, ref [5]
    fFlowTotal <- fixed(1.23) ; label("Total plasma flow per unit serum volume, F / V (1/min)")  # Table B footnote: F = 6500 mL/min at V = 5300 mL for the average normal adult
    fFlowLiver <- fixed(0.065) ; label("Hepatic arterial plasma flow as a fraction of F (unitless)")   # Table B, ref [5]
    fFlowSpleen <- fixed(0.03) ; label("Splenic plasma flow as a fraction of F (unitless)")            # Table B, ref [5]
    fFlowGi <- fixed(0.16) ; label("GI tract plasma flow as a fraction of F (unitless)")               # Table B, ref [5]

    # ---- Readout parameters (S1 Text Table A and Eqs 20, 21) -------------
    intLiver <- fixed(0.04) ; label("Liver fraction of the nonsaturable interstitial compartment, intL (unitless)")    # Table A, ref [1]
    intSpleen <- fixed(0.04) ; label("Spleen fraction of the nonsaturable interstitial compartment, intS (unitless)")  # Table A, ref [1]
    vRoiL2L4 <- fixed(0.03) ; label("Blood volume of the arteries and veins overlapping the L2-L4 region of interest (L)") # Table B: 30 mL, ref [16]
    kRefHeight <- fixed(170) ; label("Reference-man height in the UlmDos L2-L4 scaling factor (cm)")        # Table A, ref [11] Glatting 2005
    kRefFraction <- fixed(0.06665) ; label("Reference-man L2-L4 fraction of the total red marrow (unitless)") # Table A, ref [11]

    # ---- Dose split (S1 Text Eqs 10-16) ----------------------------------
    # The radiolabelled fraction of the injected antibody: a * f_l,PT for the
    # pre-therapeutic injection and b * f_l,T for therapy. Table B gives
    # radiochemical purities a = 0.94 and b = 0.96 and labelled antibody
    # fractions f_l,PT = 2.3 per cent and f_l,T = 21 per cent, so the two values
    # are 0.0216 and 0.2016. The shipped value is the pre-therapeutic one, the
    # data the model was fitted to.
    fracLabeled <- fixed(0.0216) ; label("Radiolabelled fraction of the injected antibody (unitless)")
  })

  model({
    # ================= Volumes, flows and antigen pools ==================
    # Antibody distribution volumes, S1 Text Eq 1 (L).
    vLiver <- fVolLiver * vSerum
    vSpleen <- fVolSpleen * vSerum
    vGi <- fVolGi * vSerum
    vMarrow <- fVolMarrow * vSerum
    # Table B: V_P = V - V_L - V_S - V_GI - V_RM.
    vPlasma <- vSerum - vLiver - vSpleen - vGi - vMarrow

    # Plasma flows, S1 Text Table B (L/min).
    fTotal <- fFlowTotal * vSerum
    fLiver <- fFlowLiver * fTotal
    fSpleen <- fFlowSpleen * fTotal
    fGi <- fFlowGi * fTotal
    fMarrow <- fRm * fTotal

    # Total antigen numbers (nmol). Only the liver and spleen numbers are
    # fitted; blood and red marrow follow from the granulocyte-pool ratios.
    # Eq 8 (Model 1). Eq 8 prints '* 0.9', but the sentence it follows states
    # that liver plus spleen equal 90 % of the circulating pool, and Table 1's
    # means are a linear map of the fitted pair: AgL 0.31 + AgS 0.22 = 0.53
    # against AgB 0.58, which only reproduces under division (0.589) and not
    # under multiplication (0.477). Encoded as '/ 0.9'.
    ag0Plasma <- (ag0Liver + ag0Spleen) / agBloodOrganFrac
    ag0Marrow <- agRatioRmBlood * ag0Plasma  # S1 Text Eq 4

    # CD66-expressing cell numbers, S1 Text Eq 4 and Table B. The liver and
    # spleen marginating-granulocyte numbers are proportional to their antigen
    # numbers, N_cells,i = Ag0,i / Ag0,RM * N_cells,RM.
    nCellMarrow <- nCellRmPerKg * WT
    nCellPlasma <- nCellMarrow / agRatioRmBlood
    nCellLiver <- nCellMarrow * ag0Liver / ag0Marrow
    nCellSpleen <- nCellMarrow * ag0Spleen / ag0Marrow

    # Bivalent enhancement, S1 Text Eq 3: alpha_i = E / (4 pi r_cell^2 N_cells,i),
    # which is k_on,bi * [Ag]_s / (k_on,mono * Ag_i). r_cell is tabulated in um
    # and E in 1/cm, so the cell surface is taken in cm2 and alpha_i comes out in
    # 1/cm3; the factor 1000 converts it to 1/L, the volume unit k_on,mono is
    # expressed in.
    surfCell <- 4 * pi * (rCell * 1e-4)^2
    alphaMarrow <- 1000 * enhancement / (surfCell * nCellMarrow)
    alphaPlasma <- 1000 * enhancement / (surfCell * nCellPlasma)
    alphaLiver <- 1000 * enhancement / (surfCell * nCellLiver)
    alphaSpleen <- 1000 * enhancement / (surfCell * nCellSpleen)

    # Free antigen, S1 Text Eq 2. A bivalently bound antibody occupies two.
    agPlasma <- ag0Plasma - agmono_fa_plasma_unlab - agmono_fa_plasma_lab -
      2 * agbi_fa_plasma_unlab - 2 * agbi_fa_plasma_lab -
      agmono_ha_plasma_unlab - agmono_ha_plasma_lab
    agLiver <- ag0Liver - agmono_fa_liver_unlab - agmono_fa_liver_lab -
      2 * agbi_fa_liver_unlab - 2 * agbi_fa_liver_lab -
      agmono_ha_liver_unlab - agmono_ha_liver_lab
    agSpleen <- ag0Spleen - agmono_fa_spleen_unlab - agmono_fa_spleen_lab -
      2 * agbi_fa_spleen_unlab - 2 * agbi_fa_spleen_lab -
      agmono_ha_spleen_unlab - agmono_ha_spleen_lab
    agMarrow <- ag0Marrow - agmono_fa_marrow_unlab - agmono_fa_marrow_lab -
      2 * agbi_fa_marrow_unlab - 2 * agbi_fa_marrow_lab -
      agmono_ha_marrow_unlab - agmono_ha_marrow_lab

    # ========================== Dose splitting ==========================
    # S1 Text Eqs 10-16. Dose the same total nmol of antibody into each of the
    # six main vascular compartments; bioavailability applies the radiolabelled
    # fraction and the immunoreactivity split r^2 : 2 r (1 - r) : (1 - r)^2.
    f(ab_fa_plasma_lab) <- fracLabeled * rIm^2
    f(ab_ha_plasma_lab) <- fracLabeled * 2 * rIm * (1 - rIm)
    f(ab_na_plasma_lab) <- fracLabeled * (1 - rIm)^2
    f(ab_fa_plasma_unlab) <- (1 - fracLabeled) * rIm^2
    f(ab_ha_plasma_unlab) <- (1 - fracLabeled) * 2 * rIm * (1 - rIm)
    f(ab_na_plasma_unlab) <- (1 - fracLabeled) * (1 - rIm)^2

    # ============================== ODEs ================================
    # ---- Free antibody, main vascular compartment (Eqs 28, 34, 39) -----
    # Eq 39 prints k_out * Ab_ha,int and lambda_phy * Ab*_ha,P in the
    # non-immunoreactive unlabelled equation; its own labelled twin uses the
    # na states, so the ha subscripts are a transcription slip and are
    # encoded as na.
    d/dt(ab_fa_plasma_unlab) <- -2 * konMono * agPlasma * ab_fa_plasma_unlab / vPlasma + koff * agmono_fa_plasma_unlab +
      (fLiver + fSpleen + fGi) / vLiver * ab_fa_liver_unlab +
      fMarrow / vMarrow * ab_fa_marrow_unlab -
      (fLiver + fSpleen + fGi + fMarrow) / vPlasma * ab_fa_plasma_unlab -
      (lambdaDu + kIn) * ab_fa_plasma_unlab + kOut * ab_fa_int_unlab + lambdaPhy * ab_fa_plasma_lab
    d/dt(ab_fa_plasma_lab) <- -2 * konMono * agPlasma * ab_fa_plasma_lab / vPlasma + koff * agmono_fa_plasma_lab +
      (fLiver + fSpleen + fGi) / vLiver * ab_fa_liver_lab +
      fMarrow / vMarrow * ab_fa_marrow_lab -
      (fLiver + fSpleen + fGi + fMarrow) / vPlasma * ab_fa_plasma_lab -
      (lambdaDu + kIn) * ab_fa_plasma_lab + kOut * ab_fa_int_lab - lambdaPhy * ab_fa_plasma_lab
    d/dt(ab_ha_plasma_unlab) <- -2 * konMono * agPlasma * ab_ha_plasma_unlab / vPlasma + koff * agmono_ha_plasma_unlab +
      (fLiver + fSpleen + fGi) / vLiver * ab_ha_liver_unlab +
      fMarrow / vMarrow * ab_ha_marrow_unlab -
      (fLiver + fSpleen + fGi + fMarrow) / vPlasma * ab_ha_plasma_unlab -
      (lambdaDu + kIn) * ab_ha_plasma_unlab + kOut * ab_ha_int_unlab + lambdaPhy * ab_ha_plasma_lab
    d/dt(ab_ha_plasma_lab) <- -2 * konMono * agPlasma * ab_ha_plasma_lab / vPlasma + koff * agmono_ha_plasma_lab +
      (fLiver + fSpleen + fGi) / vLiver * ab_ha_liver_lab +
      fMarrow / vMarrow * ab_ha_marrow_lab -
      (fLiver + fSpleen + fGi + fMarrow) / vPlasma * ab_ha_plasma_lab -
      (lambdaDu + kIn) * ab_ha_plasma_lab + kOut * ab_ha_int_lab - lambdaPhy * ab_ha_plasma_lab
    d/dt(ab_na_plasma_unlab) <- (fLiver + fSpleen + fGi) / vLiver * ab_na_liver_unlab +
      fMarrow / vMarrow * ab_na_marrow_unlab -
      (fLiver + fSpleen + fGi + fMarrow) / vPlasma * ab_na_plasma_unlab -
      (lambdaDu + kIn) * ab_na_plasma_unlab + kOut * ab_na_int_unlab + lambdaPhy * ab_na_plasma_lab
    d/dt(ab_na_plasma_lab) <- (fLiver + fSpleen + fGi) / vLiver * ab_na_liver_lab +
      fMarrow / vMarrow * ab_na_marrow_lab -
      (fLiver + fSpleen + fGi + fMarrow) / vPlasma * ab_na_plasma_lab -
      (lambdaDu + kIn) * ab_na_plasma_lab + kOut * ab_na_int_lab - lambdaPhy * ab_na_plasma_lab

    # ---- Free antibody, liver (Eqs 25, 31, 36) -------------------------
    d/dt(ab_fa_liver_unlab) <- -2 * konMono * agLiver * ab_fa_liver_unlab / vLiver + koff * agmono_fa_liver_unlab +
      fLiver / vPlasma * ab_fa_plasma_unlab +
      fSpleen / vSpleen * ab_fa_spleen_unlab + fGi / vGi * ab_fa_gi_unlab -
      (fLiver + fSpleen + fGi) / vLiver * ab_fa_liver_unlab + lambdaPhy * ab_fa_liver_lab
    d/dt(ab_fa_liver_lab) <- -2 * konMono * agLiver * ab_fa_liver_lab / vLiver + koff * agmono_fa_liver_lab +
      fLiver / vPlasma * ab_fa_plasma_lab +
      fSpleen / vSpleen * ab_fa_spleen_lab + fGi / vGi * ab_fa_gi_lab -
      (fLiver + fSpleen + fGi) / vLiver * ab_fa_liver_lab - lambdaPhy * ab_fa_liver_lab
    d/dt(ab_ha_liver_unlab) <- -2 * konMono * agLiver * ab_ha_liver_unlab / vLiver + koff * agmono_ha_liver_unlab +
      fLiver / vPlasma * ab_ha_plasma_unlab +
      fSpleen / vSpleen * ab_ha_spleen_unlab + fGi / vGi * ab_ha_gi_unlab -
      (fLiver + fSpleen + fGi) / vLiver * ab_ha_liver_unlab + lambdaPhy * ab_ha_liver_lab
    d/dt(ab_ha_liver_lab) <- -2 * konMono * agLiver * ab_ha_liver_lab / vLiver + koff * agmono_ha_liver_lab +
      fLiver / vPlasma * ab_ha_plasma_lab +
      fSpleen / vSpleen * ab_ha_spleen_lab + fGi / vGi * ab_ha_gi_lab -
      (fLiver + fSpleen + fGi) / vLiver * ab_ha_liver_lab - lambdaPhy * ab_ha_liver_lab
    d/dt(ab_na_liver_unlab) <- fLiver / vPlasma * ab_na_plasma_unlab +
      fSpleen / vSpleen * ab_na_spleen_unlab + fGi / vGi * ab_na_gi_unlab -
      (fLiver + fSpleen + fGi) / vLiver * ab_na_liver_unlab + lambdaPhy * ab_na_liver_lab
    d/dt(ab_na_liver_lab) <- fLiver / vPlasma * ab_na_plasma_lab +
      fSpleen / vSpleen * ab_na_spleen_lab + fGi / vGi * ab_na_gi_lab -
      (fLiver + fSpleen + fGi) / vLiver * ab_na_liver_lab - lambdaPhy * ab_na_liver_lab

    # ---- Free antibody, spleen and red marrow (Eqs 24, 30, 35) ---------
    # The unlabelled half of Eq 35 prints a minus on the inflow from the main
    # vascular compartment; its labelled twin and Eq 36 both print a plus, so
    # the minus is a transcription slip and the inflow is encoded positive.
    d/dt(ab_fa_spleen_unlab) <- -2 * konMono * agSpleen * ab_fa_spleen_unlab / vSpleen + koff * agmono_fa_spleen_unlab +
      fSpleen / vPlasma * ab_fa_plasma_unlab - fSpleen / vSpleen * ab_fa_spleen_unlab + lambdaPhy * ab_fa_spleen_lab
    d/dt(ab_fa_spleen_lab) <- -2 * konMono * agSpleen * ab_fa_spleen_lab / vSpleen + koff * agmono_fa_spleen_lab +
      fSpleen / vPlasma * ab_fa_plasma_lab - fSpleen / vSpleen * ab_fa_spleen_lab - lambdaPhy * ab_fa_spleen_lab
    d/dt(ab_ha_spleen_unlab) <- -2 * konMono * agSpleen * ab_ha_spleen_unlab / vSpleen + koff * agmono_ha_spleen_unlab +
      fSpleen / vPlasma * ab_ha_plasma_unlab - fSpleen / vSpleen * ab_ha_spleen_unlab + lambdaPhy * ab_ha_spleen_lab
    d/dt(ab_ha_spleen_lab) <- -2 * konMono * agSpleen * ab_ha_spleen_lab / vSpleen + koff * agmono_ha_spleen_lab +
      fSpleen / vPlasma * ab_ha_plasma_lab - fSpleen / vSpleen * ab_ha_spleen_lab - lambdaPhy * ab_ha_spleen_lab
    d/dt(ab_na_spleen_unlab) <- fSpleen / vPlasma * ab_na_plasma_unlab - fSpleen / vSpleen * ab_na_spleen_unlab + lambdaPhy * ab_na_spleen_lab
    d/dt(ab_na_spleen_lab) <- fSpleen / vPlasma * ab_na_plasma_lab - fSpleen / vSpleen * ab_na_spleen_lab - lambdaPhy * ab_na_spleen_lab
    d/dt(ab_fa_marrow_unlab) <- -2 * konMono * agMarrow * ab_fa_marrow_unlab / vMarrow + koff * agmono_fa_marrow_unlab +
      fMarrow / vPlasma * ab_fa_plasma_unlab - fMarrow / vMarrow * ab_fa_marrow_unlab + lambdaPhy * ab_fa_marrow_lab
    d/dt(ab_fa_marrow_lab) <- -2 * konMono * agMarrow * ab_fa_marrow_lab / vMarrow + koff * agmono_fa_marrow_lab +
      fMarrow / vPlasma * ab_fa_plasma_lab - fMarrow / vMarrow * ab_fa_marrow_lab - lambdaPhy * ab_fa_marrow_lab
    d/dt(ab_ha_marrow_unlab) <- -2 * konMono * agMarrow * ab_ha_marrow_unlab / vMarrow + koff * agmono_ha_marrow_unlab +
      fMarrow / vPlasma * ab_ha_plasma_unlab - fMarrow / vMarrow * ab_ha_marrow_unlab + lambdaPhy * ab_ha_marrow_lab
    d/dt(ab_ha_marrow_lab) <- -2 * konMono * agMarrow * ab_ha_marrow_lab / vMarrow + koff * agmono_ha_marrow_lab +
      fMarrow / vPlasma * ab_ha_plasma_lab - fMarrow / vMarrow * ab_ha_marrow_lab - lambdaPhy * ab_ha_marrow_lab
    d/dt(ab_na_marrow_unlab) <- fMarrow / vPlasma * ab_na_plasma_unlab - fMarrow / vMarrow * ab_na_marrow_unlab + lambdaPhy * ab_na_marrow_lab
    d/dt(ab_na_marrow_lab) <- fMarrow / vPlasma * ab_na_plasma_lab - fMarrow / vMarrow * ab_na_marrow_lab - lambdaPhy * ab_na_marrow_lab

    # ---- Free antibody, gastrointestinal tract (Eqs 26, 32, 37) --------
    d/dt(ab_fa_gi_unlab) <- fGi / vPlasma * ab_fa_plasma_unlab - fGi / vGi * ab_fa_gi_unlab + lambdaPhy * ab_fa_gi_lab
    d/dt(ab_fa_gi_lab) <- fGi / vPlasma * ab_fa_plasma_lab - fGi / vGi * ab_fa_gi_lab - lambdaPhy * ab_fa_gi_lab
    d/dt(ab_ha_gi_unlab) <- fGi / vPlasma * ab_ha_plasma_unlab - fGi / vGi * ab_ha_gi_unlab + lambdaPhy * ab_ha_gi_lab
    d/dt(ab_ha_gi_lab) <- fGi / vPlasma * ab_ha_plasma_lab - fGi / vGi * ab_ha_gi_lab - lambdaPhy * ab_ha_gi_lab
    d/dt(ab_na_gi_unlab) <- fGi / vPlasma * ab_na_plasma_unlab - fGi / vGi * ab_na_gi_unlab + lambdaPhy * ab_na_gi_lab
    d/dt(ab_na_gi_lab) <- fGi / vPlasma * ab_na_plasma_lab - fGi / vGi * ab_na_gi_lab - lambdaPhy * ab_na_gi_lab

    # ---- Free antibody, interstitial space (Eqs 27, 33, 38) ------------
    d/dt(ab_fa_int_unlab) <- kIn * ab_fa_plasma_unlab - kOut * ab_fa_int_unlab + lambdaPhy * ab_fa_int_lab
    d/dt(ab_fa_int_lab) <- kIn * ab_fa_plasma_lab - kOut * ab_fa_int_lab - lambdaPhy * ab_fa_int_lab
    d/dt(ab_ha_int_unlab) <- kIn * ab_ha_plasma_unlab - kOut * ab_ha_int_unlab + lambdaPhy * ab_ha_int_lab
    d/dt(ab_ha_int_lab) <- kIn * ab_ha_plasma_lab - kOut * ab_ha_int_lab - lambdaPhy * ab_ha_int_lab
    d/dt(ab_na_int_unlab) <- kIn * ab_na_plasma_unlab - kOut * ab_na_int_unlab + lambdaPhy * ab_na_int_lab
    d/dt(ab_na_int_lab) <- kIn * ab_na_plasma_lab - kOut * ab_na_int_lab - lambdaPhy * ab_na_int_lab

    # ---- Bivalently bound fully immunoreactive antibody (Eq 22) --------
    d/dt(agbi_fa_plasma_unlab) <- konMono * alphaPlasma * agPlasma * agmono_fa_plasma_unlab -
      2 * koff * agbi_fa_plasma_unlab - lambdaDb * agbi_fa_plasma_unlab + lambdaPhy * agbi_fa_plasma_lab
    d/dt(agbi_fa_plasma_lab) <- konMono * alphaPlasma * agPlasma * agmono_fa_plasma_lab -
      2 * koff * agbi_fa_plasma_lab - lambdaDb * agbi_fa_plasma_lab - lambdaPhy * agbi_fa_plasma_lab
    d/dt(agbi_fa_liver_unlab) <- konMono * alphaLiver * agLiver * agmono_fa_liver_unlab -
      2 * koff * agbi_fa_liver_unlab - lambdaDb * agbi_fa_liver_unlab + lambdaPhy * agbi_fa_liver_lab
    d/dt(agbi_fa_liver_lab) <- konMono * alphaLiver * agLiver * agmono_fa_liver_lab -
      2 * koff * agbi_fa_liver_lab - lambdaDb * agbi_fa_liver_lab - lambdaPhy * agbi_fa_liver_lab
    d/dt(agbi_fa_spleen_unlab) <- konMono * alphaSpleen * agSpleen * agmono_fa_spleen_unlab -
      2 * koff * agbi_fa_spleen_unlab - lambdaDb * agbi_fa_spleen_unlab + lambdaPhy * agbi_fa_spleen_lab
    d/dt(agbi_fa_spleen_lab) <- konMono * alphaSpleen * agSpleen * agmono_fa_spleen_lab -
      2 * koff * agbi_fa_spleen_lab - lambdaDb * agbi_fa_spleen_lab - lambdaPhy * agbi_fa_spleen_lab
    d/dt(agbi_fa_marrow_unlab) <- konMono * alphaMarrow * agMarrow * agmono_fa_marrow_unlab -
      2 * koff * agbi_fa_marrow_unlab - lambdaDb * agbi_fa_marrow_unlab + lambdaPhy * agbi_fa_marrow_lab
    d/dt(agbi_fa_marrow_lab) <- konMono * alphaMarrow * agMarrow * agmono_fa_marrow_lab -
      2 * koff * agbi_fa_marrow_lab - lambdaDb * agbi_fa_marrow_lab - lambdaPhy * agbi_fa_marrow_lab

    # ---- Monovalently bound fully immunoreactive antibody (Eq 23) ------
    d/dt(agmono_fa_plasma_unlab) <- 2 * konMono * agPlasma * ab_fa_plasma_unlab / vPlasma -
      konMono * alphaPlasma * agPlasma * agmono_fa_plasma_unlab - koff * agmono_fa_plasma_unlab +
      2 * koff * agbi_fa_plasma_unlab - lambdaDb * agmono_fa_plasma_unlab + lambdaPhy * agmono_fa_plasma_lab
    d/dt(agmono_fa_plasma_lab) <- 2 * konMono * agPlasma * ab_fa_plasma_lab / vPlasma -
      konMono * alphaPlasma * agPlasma * agmono_fa_plasma_lab - koff * agmono_fa_plasma_lab +
      2 * koff * agbi_fa_plasma_lab - lambdaDb * agmono_fa_plasma_lab - lambdaPhy * agmono_fa_plasma_lab
    d/dt(agmono_fa_liver_unlab) <- 2 * konMono * agLiver * ab_fa_liver_unlab / vLiver -
      konMono * alphaLiver * agLiver * agmono_fa_liver_unlab - koff * agmono_fa_liver_unlab +
      2 * koff * agbi_fa_liver_unlab - lambdaDb * agmono_fa_liver_unlab + lambdaPhy * agmono_fa_liver_lab
    d/dt(agmono_fa_liver_lab) <- 2 * konMono * agLiver * ab_fa_liver_lab / vLiver -
      konMono * alphaLiver * agLiver * agmono_fa_liver_lab - koff * agmono_fa_liver_lab +
      2 * koff * agbi_fa_liver_lab - lambdaDb * agmono_fa_liver_lab - lambdaPhy * agmono_fa_liver_lab
    d/dt(agmono_fa_spleen_unlab) <- 2 * konMono * agSpleen * ab_fa_spleen_unlab / vSpleen -
      konMono * alphaSpleen * agSpleen * agmono_fa_spleen_unlab - koff * agmono_fa_spleen_unlab +
      2 * koff * agbi_fa_spleen_unlab - lambdaDb * agmono_fa_spleen_unlab + lambdaPhy * agmono_fa_spleen_lab
    d/dt(agmono_fa_spleen_lab) <- 2 * konMono * agSpleen * ab_fa_spleen_lab / vSpleen -
      konMono * alphaSpleen * agSpleen * agmono_fa_spleen_lab - koff * agmono_fa_spleen_lab +
      2 * koff * agbi_fa_spleen_lab - lambdaDb * agmono_fa_spleen_lab - lambdaPhy * agmono_fa_spleen_lab
    d/dt(agmono_fa_marrow_unlab) <- 2 * konMono * agMarrow * ab_fa_marrow_unlab / vMarrow -
      konMono * alphaMarrow * agMarrow * agmono_fa_marrow_unlab - koff * agmono_fa_marrow_unlab +
      2 * koff * agbi_fa_marrow_unlab - lambdaDb * agmono_fa_marrow_unlab + lambdaPhy * agmono_fa_marrow_lab
    d/dt(agmono_fa_marrow_lab) <- 2 * konMono * agMarrow * ab_fa_marrow_lab / vMarrow -
      konMono * alphaMarrow * agMarrow * agmono_fa_marrow_lab - koff * agmono_fa_marrow_lab +
      2 * koff * agbi_fa_marrow_lab - lambdaDb * agmono_fa_marrow_lab - lambdaPhy * agmono_fa_marrow_lab

    # ---- Monovalently bound half immunoreactive antibody (Eq 29) -------
    # As printed, Eq 29 also carries -k_on,mono * alpha_i * Ag_i * AgAb_ha,mono
    # and +2 k_off * AgAb_ha,bi. A half immunoreactive antibody has one active
    # arm and cannot crosslink, no AgAb_ha,bi state exists, and Eq 2 carries no
    # ha bivalent term; keeping the alpha loss would drain bound antibody into
    # nothing and break mass balance. Both terms are dropped as copies of the
    # fully immunoreactive Eq 23. The printed factor 2 on the association is
    # kept as published.
    d/dt(agmono_ha_plasma_unlab) <- 2 * konMono * agPlasma * ab_ha_plasma_unlab / vPlasma -
      koff * agmono_ha_plasma_unlab - lambdaDb * agmono_ha_plasma_unlab + lambdaPhy * agmono_ha_plasma_lab
    d/dt(agmono_ha_plasma_lab) <- 2 * konMono * agPlasma * ab_ha_plasma_lab / vPlasma -
      koff * agmono_ha_plasma_lab - lambdaDb * agmono_ha_plasma_lab - lambdaPhy * agmono_ha_plasma_lab
    d/dt(agmono_ha_liver_unlab) <- 2 * konMono * agLiver * ab_ha_liver_unlab / vLiver -
      koff * agmono_ha_liver_unlab - lambdaDb * agmono_ha_liver_unlab + lambdaPhy * agmono_ha_liver_lab
    d/dt(agmono_ha_liver_lab) <- 2 * konMono * agLiver * ab_ha_liver_lab / vLiver -
      koff * agmono_ha_liver_lab - lambdaDb * agmono_ha_liver_lab - lambdaPhy * agmono_ha_liver_lab
    d/dt(agmono_ha_spleen_unlab) <- 2 * konMono * agSpleen * ab_ha_spleen_unlab / vSpleen -
      koff * agmono_ha_spleen_unlab - lambdaDb * agmono_ha_spleen_unlab + lambdaPhy * agmono_ha_spleen_lab
    d/dt(agmono_ha_spleen_lab) <- 2 * konMono * agSpleen * ab_ha_spleen_lab / vSpleen -
      koff * agmono_ha_spleen_lab - lambdaDb * agmono_ha_spleen_lab - lambdaPhy * agmono_ha_spleen_lab
    d/dt(agmono_ha_marrow_unlab) <- 2 * konMono * agMarrow * ab_ha_marrow_unlab / vMarrow -
      koff * agmono_ha_marrow_unlab - lambdaDb * agmono_ha_marrow_unlab + lambdaPhy * agmono_ha_marrow_lab
    d/dt(agmono_ha_marrow_lab) <- 2 * konMono * agMarrow * ab_ha_marrow_lab / vMarrow -
      koff * agmono_ha_marrow_lab - lambdaDb * agmono_ha_marrow_lab - lambdaPhy * agmono_ha_marrow_lab

    # ---- Degraded-antibody submodel, Eger/Houston (Eqs 40-43) ---------
    d/dt(ex_fa_unlab) <- lambdaDu * ab_fa_plasma_unlab - lambdaClex * ex_fa_unlab + lambdaPhy * ex_fa_lab
    d/dt(ex_fa_lab) <- lambdaDu * ab_fa_plasma_lab - lambdaClex * ex_fa_lab - lambdaPhy * ex_fa_lab
    d/dt(ex_ha_unlab) <- lambdaDu * ab_ha_plasma_unlab - lambdaClex * ex_ha_unlab + lambdaPhy * ex_ha_lab
    d/dt(ex_ha_lab) <- lambdaDu * ab_ha_plasma_lab - lambdaClex * ex_ha_lab - lambdaPhy * ex_ha_lab
    d/dt(ex_na_unlab) <- lambdaDu * ab_na_plasma_unlab - lambdaClex * ex_na_unlab + lambdaPhy * ex_na_lab
    d/dt(ex_na_lab) <- lambdaDu * ab_na_plasma_lab - lambdaClex * ex_na_lab - lambdaPhy * ex_na_lab

    d/dt(metap_fa_unlab) <- lambdaClex * ex_fa_unlab - lambdaCl * metap_fa_unlab -
      lambdaMetaex1 * metap_fa_unlab + lambdaMetaex2 * metaex1_fa_unlab +
      lambdaDb * (agmono_fa_plasma_unlab + agbi_fa_plasma_unlab +
        agmono_fa_liver_unlab + agbi_fa_liver_unlab +
        agmono_fa_spleen_unlab + agbi_fa_spleen_unlab +
        agmono_fa_marrow_unlab + agbi_fa_marrow_unlab) + lambdaPhy * metap_fa_lab
    d/dt(metap_fa_lab) <- lambdaClex * ex_fa_lab - lambdaCl * metap_fa_lab -
      lambdaMetaex1 * metap_fa_lab + lambdaMetaex2 * metaex1_fa_lab +
      lambdaDb * (agmono_fa_plasma_lab + agbi_fa_plasma_lab +
        agmono_fa_liver_lab + agbi_fa_liver_lab +
        agmono_fa_spleen_lab + agbi_fa_spleen_lab +
        agmono_fa_marrow_lab + agbi_fa_marrow_lab) - lambdaPhy * metap_fa_lab
    d/dt(metap_ha_unlab) <- lambdaClex * ex_ha_unlab - lambdaCl * metap_ha_unlab -
      lambdaMetaex1 * metap_ha_unlab + lambdaMetaex2 * metaex1_ha_unlab +
      lambdaDb * (agmono_ha_plasma_unlab + agmono_ha_liver_unlab + agmono_ha_spleen_unlab + agmono_ha_marrow_unlab) + lambdaPhy * metap_ha_lab
    d/dt(metap_ha_lab) <- lambdaClex * ex_ha_lab - lambdaCl * metap_ha_lab -
      lambdaMetaex1 * metap_ha_lab + lambdaMetaex2 * metaex1_ha_lab +
      lambdaDb * (agmono_ha_plasma_lab + agmono_ha_liver_lab + agmono_ha_spleen_lab + agmono_ha_marrow_lab) - lambdaPhy * metap_ha_lab
    d/dt(metap_na_unlab) <- lambdaClex * ex_na_unlab - lambdaCl * metap_na_unlab -
      lambdaMetaex1 * metap_na_unlab + lambdaMetaex2 * metaex1_na_unlab + lambdaPhy * metap_na_lab
    d/dt(metap_na_lab) <- lambdaClex * ex_na_lab - lambdaCl * metap_na_lab -
      lambdaMetaex1 * metap_na_lab + lambdaMetaex2 * metaex1_na_lab - lambdaPhy * metap_na_lab

    d/dt(metaex1_fa_unlab) <- lambdaMetaex1 * metap_fa_unlab -
      (lambdaMetaex2 + lambdaMetaex3) * metaex1_fa_unlab +
      lambdaMetaex4 * metaex2_fa_unlab + lambdaPhy * metaex1_fa_lab
    d/dt(metaex1_fa_lab) <- lambdaMetaex1 * metap_fa_lab -
      (lambdaMetaex2 + lambdaMetaex3) * metaex1_fa_lab +
      lambdaMetaex4 * metaex2_fa_lab - lambdaPhy * metaex1_fa_lab
    d/dt(metaex1_ha_unlab) <- lambdaMetaex1 * metap_ha_unlab -
      (lambdaMetaex2 + lambdaMetaex3) * metaex1_ha_unlab +
      lambdaMetaex4 * metaex2_ha_unlab + lambdaPhy * metaex1_ha_lab
    d/dt(metaex1_ha_lab) <- lambdaMetaex1 * metap_ha_lab -
      (lambdaMetaex2 + lambdaMetaex3) * metaex1_ha_lab +
      lambdaMetaex4 * metaex2_ha_lab - lambdaPhy * metaex1_ha_lab
    d/dt(metaex1_na_unlab) <- lambdaMetaex1 * metap_na_unlab -
      (lambdaMetaex2 + lambdaMetaex3) * metaex1_na_unlab +
      lambdaMetaex4 * metaex2_na_unlab + lambdaPhy * metaex1_na_lab
    d/dt(metaex1_na_lab) <- lambdaMetaex1 * metap_na_lab -
      (lambdaMetaex2 + lambdaMetaex3) * metaex1_na_lab +
      lambdaMetaex4 * metaex2_na_lab - lambdaPhy * metaex1_na_lab

    d/dt(metaex2_fa_unlab) <- lambdaMetaex3 * metaex1_fa_unlab - lambdaMetaex4 * metaex2_fa_unlab + lambdaPhy * metaex2_fa_lab
    d/dt(metaex2_fa_lab) <- lambdaMetaex3 * metaex1_fa_lab - lambdaMetaex4 * metaex2_fa_lab - lambdaPhy * metaex2_fa_lab
    d/dt(metaex2_ha_unlab) <- lambdaMetaex3 * metaex1_ha_unlab - lambdaMetaex4 * metaex2_ha_unlab + lambdaPhy * metaex2_ha_lab
    d/dt(metaex2_ha_lab) <- lambdaMetaex3 * metaex1_ha_lab - lambdaMetaex4 * metaex2_ha_lab - lambdaPhy * metaex2_ha_lab
    d/dt(metaex2_na_unlab) <- lambdaMetaex3 * metaex1_na_unlab - lambdaMetaex4 * metaex2_na_unlab + lambdaPhy * metaex2_na_lab
    d/dt(metaex2_na_lab) <- lambdaMetaex3 * metaex1_na_lab - lambdaMetaex4 * metaex2_na_lab - lambdaPhy * metaex2_na_lab

    # ---- Cleared from the body (bookkeeping, not a published state) ----
    # lambda_cl is the only route out of the system (Table B, 'clearance from
    # body'); accumulating it makes total antibody mass exactly conserved and
    # gives the vignette a mass-balance gate.
    d/dt(cleared_unlab) <- lambdaCl * (metap_fa_unlab + metap_ha_unlab + metap_na_unlab) + lambdaPhy * cleared_lab
    d/dt(cleared_lab) <- lambdaCl * (metap_fa_lab + metap_ha_lab + metap_na_lab) - lambdaPhy * cleared_lab

    # ========================= Observed signals =========================
    # S1 Text Eqs 17-21. Every readout counts radiolabelled material only.
    metaLab <- metap_fa_lab + metap_ha_lab + metap_na_lab
    exLab <- ex_fa_lab + ex_ha_lab + ex_na_lab
    intLab <- ab_fa_int_lab + ab_ha_int_lab + ab_na_int_lab
    boundPlasmaLab <- agbi_fa_plasma_lab + agmono_fa_plasma_lab + agmono_ha_plasma_lab
    freePlasmaLab <- ab_fa_plasma_lab + ab_ha_plasma_lab + ab_na_plasma_lab

    # Eq 17: liver region of interest (nmol).
    abLiverRoi <- agbi_fa_liver_lab + agmono_fa_liver_lab + ab_fa_liver_lab +
      agmono_ha_liver_lab + ab_ha_liver_lab + ab_na_liver_lab +
      exLiver * exLab + (vLiver / vSerum) * (metaLab + boundPlasmaLab) +
      intLiver * intLab

    # Eq 18: spleen region of interest (nmol).
    abSpleenRoi <- agbi_fa_spleen_lab + agmono_fa_spleen_lab + ab_fa_spleen_lab +
      agmono_ha_spleen_lab + ab_ha_spleen_lab + ab_na_spleen_lab +
      exSpleen * exLab + (vSpleen / vSerum) * (metaLab + boundPlasmaLab) +
      intSpleen * intLab

    # Eq 19: serum concentration of radiolabelled material (nmol/L).
    Cc <- (freePlasmaLab + (vPlasma / vSerum) * metaLab) / vPlasma

    # Total red marrow content: the bracket of Eq 20 without the blood
    # background x. This is the quantity the red-marrow time-integrated activity
    # coefficient of Table 2 is computed from.
    abMarrow <- agbi_fa_marrow_lab + agmono_fa_marrow_lab + ab_fa_marrow_lab +
      agmono_ha_marrow_lab + ab_ha_marrow_lab + ab_na_marrow_lab +
      (vMarrow / vSerum) * (metaLab + boundPlasmaLab)

    # Eq 21: activity of the arteries and veins overlapping the L2-L4 region of
    # interest.
    xBlood <- (vRoiL2L4 / vPlasma) *
      (freePlasmaLab + (vPlasma / vSerum) * (metaLab + boundPlasmaLab))

    # Eq 20 with Table A: the L2-L4 gamma-camera reading the model predicts.
    kScale <- (kRefHeight / HT) / kRefFraction
    abMarrowRoi <- cRm * (abMarrow + xBlood) / kScale

    # Whole-body retained radiolabelled material (nmol). Not a published
    # equation: the whole-body gamma-camera signal is every radiolabelled state
    # that has not yet been excreted.
    abWholeBody <- ab_fa_plasma_lab + ab_fa_liver_lab + ab_fa_spleen_lab + ab_fa_marrow_lab +
      ab_fa_gi_lab + ab_fa_int_lab + ab_ha_plasma_lab + ab_ha_liver_lab +
      ab_ha_spleen_lab + ab_ha_marrow_lab + ab_ha_gi_lab + ab_ha_int_lab +
      ab_na_plasma_lab + ab_na_liver_lab + ab_na_spleen_lab + ab_na_marrow_lab +
      ab_na_gi_lab + ab_na_int_lab + agbi_fa_plasma_lab + agbi_fa_liver_lab +
      agbi_fa_spleen_lab + agbi_fa_marrow_lab + agmono_fa_plasma_lab + agmono_fa_liver_lab +
      agmono_fa_spleen_lab + agmono_fa_marrow_lab + agmono_ha_plasma_lab + agmono_ha_liver_lab +
      agmono_ha_spleen_lab + agmono_ha_marrow_lab + ex_fa_lab + ex_ha_lab +
      ex_na_lab + metap_fa_lab + metap_ha_lab + metap_na_lab +
      metaex1_fa_lab + metaex1_ha_lab + metaex1_na_lab + metaex2_fa_lab +
      metaex2_ha_lab + metaex2_na_lab

    # Total antibody-derived mass in the system, for the mass-balance gate. The
    # only route out is lambda_cl, so this is conserved and equals the injected
    # amount.
    abMassBalance <- abWholeBody + cleared_lab +
      ab_fa_plasma_unlab + ab_fa_liver_unlab + ab_fa_spleen_unlab + ab_fa_marrow_unlab +
      ab_fa_gi_unlab + ab_fa_int_unlab + ab_ha_plasma_unlab + ab_ha_liver_unlab +
      ab_ha_spleen_unlab + ab_ha_marrow_unlab + ab_ha_gi_unlab + ab_ha_int_unlab +
      ab_na_plasma_unlab + ab_na_liver_unlab + ab_na_spleen_unlab + ab_na_marrow_unlab +
      ab_na_gi_unlab + ab_na_int_unlab + agbi_fa_plasma_unlab + agbi_fa_liver_unlab +
      agbi_fa_spleen_unlab + agbi_fa_marrow_unlab + agmono_fa_plasma_unlab + agmono_fa_liver_unlab +
      agmono_fa_spleen_unlab + agmono_fa_marrow_unlab + agmono_ha_plasma_unlab + agmono_ha_liver_unlab +
      agmono_ha_spleen_unlab + agmono_ha_marrow_unlab + ex_fa_unlab + ex_ha_unlab +
      ex_na_unlab + metap_fa_unlab + metap_ha_unlab + metap_na_unlab +
      metaex1_fa_unlab + metaex1_ha_unlab + metaex1_na_unlab + metaex2_fa_unlab +
      metaex2_ha_unlab + metaex2_na_unlab + cleared_unlab
  })
}
