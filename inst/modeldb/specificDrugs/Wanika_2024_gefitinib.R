Wanika_2024_gefitinib <- function() {
  description <- "In vitro (NSCLC cell lines; H292, HCC827, PC9, A549, H1650). Mechanistic in vitro PK/PD model of gefitinib uptake into non-small cell lung cancer cells and inhibition of tumour cell viability over a 72-h exposure (Wanika 2024 Eqs 1-5, Table 2 gefitinib rows). The drug is applied to the culture medium at t = 0 (dose = medium concentration in uM) and exchanges first-order with the pooled intracellular space; the medium-to-cell and cell-to-medium terms are scaled by the medium:cell volume ratio, and the cell volume shrinks in proportion to viability (Eq 3), which concentrates the drug in the remaining cells. Viability follows an indirect-response model with zero-order production kin, inhibited by an Imax function of the intracellular concentration, and first-order loss kout. The per-curve Imax and IC50 are KNOWN inputs (digitised from each source study) supplied as covariates; the initial cell volume cancels out of the equations. Fitted in Monolix to 82 digitised concentration-viability curves; lognormal IIV on all five structural parameters, additive residual error. Sibling models: Wanika_2024_alectinib, Wanika_2024_crizotinib, Wanika_2024_erlotinib (same structure, fitted separately)."
  reference <- paste(
    "Wanika L, Evans ND, Johnson M, Tomkinson H, Chappell MJ (2024).",
    "In vitro PK/PD modeling of tyrosine kinase inhibitors in non-small cell",
    "lung cancer cell lines.",
    "Clinical and Translational Science 17(3):e13714.",
    "doi:10.1111/cts.13714.",
    sep = " "
  )
  vignette <- "Wanika_2024_tki_nsclc_invitro"

  units <- list(
    time = "h",
    dosing = "umol/L (drug concentration applied to the culture medium at t = 0)",
    concentration = "umol/L (medium and intracellular drug concentration); % (tumour cell viability, the PD readout)"
  )

  compartmentData <- list(
    medium = list(analyte = "gefitinib", units = "umol/L", specimen = "administration site", verified = TRUE),
    intracellular = list(analyte = "gefitinib", units = "umol/L", specimen = "tumor", verified = TRUE),
    viability = list(analyte = "tumour cell viability", units = "%", specimen = "not applicable", verified = TRUE)
  )

  # medium / intracellular are the culture-medium and pooled-cell drug
  # concentration states of an in vitro uptake model (paper TKI_M / TKI_C);
  # neither is a registered canonical compartment, so both are declared
  # paper-specific here, as Mody_2023_doxorubicin_dexrazoxane_jimt1.R does
  # for its medium states.
  paper_specific_compartments <- c("medium", "intracellular")

  covariateData <- list(
    IC50_VIABILITY = list(
      description = "Alectinib concentration producing 50% inhibition of tumour cell viability in the source study's 72-h assay for this cell line (a known per-curve input, not estimated)",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Digitised by the paper's authors from each source study and supplied as a Monolix regressor (Wanika 2024 Table S4, column 'IC50 (uM)'; range 0.0052-8.42 uM). Enters Eq 4 as the half-maximal intracellular concentration.",
      source_name = "IC50"
    ),
    IMAX_VIABILITY = list(
      description = "Maximum fractional inhibition of tumour cell viability by gefitinib in the source study's 72-h assay for this cell line (0-1; a known per-curve input, not estimated)",
      units = "(fraction)",
      type = "continuous",
      reference_category = NULL,
      notes = "Digitised by the paper's authors from each source study and supplied as a Monolix regressor (Wanika 2024 Table S4, column 'IMAX'; range 0.753-0.983). Table 1 of the paper lists the unit as %, but the supplied values and Eq 4 use a fraction.",
      source_name = "IMAX"
    )
  )

  population <- list(
    species = "in vitro (NSCLC cell lines)",
    n_subjects = 82L,
    n_studies = 7L,
    disease_state = "Non-small cell lung cancer cell lines (EGFR-mutant HCC827, PC9 and H1650; KRAS-mutant A549; H292) exposed to gefitinib for 72 h.",
    dose_range = "0.001-50 umol/L gefitinib in the culture medium.",
    cell_line = "H292, HCC827, PC9, A549, H1650",
    n_observations = "82 cell-viability observations at 72 h (each curve also carries the 100% viability record at t = 0).",
    notes = "Each 'subject' is one concentration level of one digitised concentration-viability curve (Wanika 2024 Table S4, IDs 1-82), taken with WebPlotDigitizer from 7 published in vitro studies. IIV therefore spans cell line, study and experimental conditions."
  )

  ini({
    # Wanika 2024 Table 2, gefitinib rows (Monolix SAEM). Means are the
    # population values; 'Variance' is the lognormal IIV variance (omega^2).
    lk_medium_intracellular <- log(0.856)
    label("First-order uptake rate constant from medium into cells KDRUGIN (1/h)") # Table 2: 0.856 (RSE 3.927%)
    lk_intracellular_medium <- log(0.313)
    label("First-order efflux rate constant from cells back into medium KDRUGOUT (1/h)") # Table 2: 0.313 (RSE 17.533%)
    lvr <- log(3.192)
    label("Medium-to-initial-cell volume ratio VR = VM / VC(0) (unitless)") # Table 2: 3.192 (RSE 15.527%)
    lkin <- log(0.941)
    label("Zero-order tumour cell viability production rate KIN (%/h)") # Table 2: 0.941 (RSE 4.420%)
    lkout <- log(0.034)
    label("First-order tumour cell viability loss rate constant KOUT (1/h)") # Table 2: 0.034 (RSE 10.693%)
    lrbase <- fixed(log(100))
    label("Initial tumour cell viability TumorCellViability(0) (%)") # Table 1: initial value 100 (%)

    etalk_medium_intracellular ~ 0.113 # Table 2: variance 0.113 (shrinkage 100%)
    etalk_intracellular_medium ~ 0.355 # Table 2: variance 0.355 (shrinkage 100%)
    etalvr ~ 0.198 # Table 2: variance 0.198 (shrinkage 100%)
    etalkin ~ 0.068 # Table 2: variance 0.068 (shrinkage 99.3%)
    etalkout ~ 0.760 # Table 2: variance 0.760 (shrinkage -0.90%)

    addSd_viability <- 0.737
    label("Additive residual SD on tumour cell viability a (%)") # Table 2: a = 0.737 (RSE 7.807%); Eq 9 constant error model
  })

  model({
    k_medium_intracellular <- exp(lk_medium_intracellular + etalk_medium_intracellular)
    k_intracellular_medium <- exp(lk_intracellular_medium + etalk_intracellular_medium)
    vr <- exp(lvr + etalvr)
    kin <- exp(lkin + etalkin)
    kout <- exp(lkout + etalkout)
    rbase <- exp(lrbase)

    # Eq 3: VC(t) = VC(0) * viability(t) / viability(0). With Eq 5,
    # VM = VC(0) * VR, the medium:cell volume ratio is VM / VC(t) =
    # vr / vcRel, so the absolute initial cell volume VC(0) cancels.
    vcRel <- viability / rbase

    # Eqs 1-2: concentration states; each transfer term is scaled by the
    # volume ratio so that drug amount is conserved between the two spaces.
    d/dt(medium) <- -k_medium_intracellular * medium +
      k_intracellular_medium * (vcRel / vr) * intracellular
    d/dt(intracellular) <- k_medium_intracellular * (vr / vcRel) * medium -
      k_intracellular_medium * intracellular

    # Eq 4: indirect response, zero-order production inhibited by the
    # intracellular concentration through the known per-curve Imax / IC50.
    d/dt(viability) <- kin * (1 - IMAX_VIABILITY * intracellular / (IC50_VIABILITY + intracellular)) -
      kout * viability
    viability(0) <- rbase

    viability ~ add(addSd_viability)
  })
}
