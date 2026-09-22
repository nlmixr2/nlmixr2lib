Lu_2016_pinatuzumab_polatuzumab <- function() {
  description <- "Integrated two-analyte population PK model for the MMAE antibody-drug conjugates pinatuzumab vedotin and polatuzumab vedotin in patients with relapsed/refractory B-cell non-Hodgkin lymphoma (Lu 2016). Antibody-conjugated MMAE (acMMAE; output Cc) and total antibody (Tab; output Cc_tab) are each described by a linear two-compartment model that shares CL, Q, Vc, Vp and all random effects, the only structural difference being an additional first-order deconjugation loss kdec from the acMMAE central compartment. The two ADCs are fitted jointly in one run: CL, Vc, kdec and the acMMAE assay cross-calibration slope take molecule-specific values selected by the TRT_PINATUZUMAB_VEDOTIN / TRT_POLATUZUMAB_VEDOTIN indicators, while Q, Vp, the IIV and the residual errors are shared. States hold molar amounts (nmol); the two observables are returned in ng/mL using the molecular weights carried in the published control stream. No covariates were assessed. Simulation requires dosing central and central_tab simultaneously with the same molar antibody amount; f(central) applies the mean drug-to-antibody ratio of 3.585."
  reference <- "Lu D, Gibiansky L, Agarwal P, Dere RC, Li C, Chu Y-W, Hirata J, Joshi A, Jin JY, Girish S. Integrated Two-Analyte Population Pharmacokinetic Model for Antibody-Drug Conjugates in Patients: Implications for Reducing Pharmacokinetic Sampling. CPT Pharmacometrics Syst Pharmacol. 2016;5(12):665-673. doi:10.1002/psp4.12137"
  vignette <- "Lu_2016_pinatuzumab_polatuzumab"
  units <- list(time = "h", dosing = "nmol", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Lu 2016 Figure 2 defines A1/A2 as the molar amounts of
  # Tab and A3/A4 as the molar amounts of acMMAE; the supplement S1 control
  # stream declares the dose unit as nM (nmol) and the concentration unit as
  # nmol/L before the molecular-weight conversion in $ERROR. Tab was assayed by
  # ELISA in serum and acMMAE by LC-MS/MS in plasma (Lu 2016 Bioanalytical
  # methods), so the two analytes carry different specimens.
  compartmentData <- list(
    central = list(analyte = "acMMAE", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "acMMAE", units = "nmol", specimen = "plasma", verified = TRUE),
    central_tab = list(analyte = "total antibody", units = "nmol", specimen = "serum", verified = TRUE),
    peripheral1_tab = list(analyte = "total antibody", units = "nmol", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    TRT_PINATUZUMAB_VEDOTIN = list(
      description = "Study-treatment indicator: 1 = the subject received pinatuzumab vedotin (anti-CD22 MMAE ADC), 0 = did not",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the subject received polatuzumab vedotin)",
      notes = "Mutually exclusive with TRT_POLATUZUMAB_VEDOTIN; exactly one of the two must be 1 on every record. Selects the pinatuzumab vedotin values of CL, Vc, kdec and the acMMAE assay cross-calibration slope. Lu 2016 supplement S1 carries this as the integer column MOL, with CL = (THETA(1)*(2-MOL) + THETA(7)*(MOL-1)) and similarly for V1, KDEC and CORR, so MOL = 1 selects the first (pinatuzumab vedotin) THETA block.",
      source_name = "MOL (MOL = 1)"
    ),
    TRT_POLATUZUMAB_VEDOTIN = list(
      description = "Study-treatment indicator: 1 = the subject received polatuzumab vedotin (anti-CD79b MMAE ADC), 0 = did not",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the subject received pinatuzumab vedotin)",
      notes = "Mutually exclusive with TRT_PINATUZUMAB_VEDOTIN; exactly one of the two must be 1 on every record. Selects the polatuzumab vedotin values of CL, Vc, kdec and the acMMAE assay cross-calibration slope, i.e. MOL = 2 in the Lu 2016 supplement S1 control stream.",
      source_name = "MOL (MOL = 2)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 154,
    n_studies = 2,
    disease_state = "Relapsed or refractory B-cell non-Hodgkin's lymphoma; patients with chronic lymphocytic leukaemia were excluded from the analysis datasets.",
    dose_range = "Pinatuzumab vedotin 0.1 to 3.2 mg/kg Q3W (dose expansion at 2.4 mg/kg) and polatuzumab vedotin 0.1 to 2.4 mg/kg Q3W (dose expansion at 2.4 mg/kg), each as a single agent; plus 1.8 and 2.4 mg/kg Q3W in combination with rituximab.",
    studies = "Dataset 1 (the phase I development dataset behind Model 1) pooled DCT4862g (NCT01209130, pinatuzumab vedotin: 65 single-agent subjects plus 15 rituximab-combination subjects) and DCS4968g (NCT01290549, polatuzumab vedotin: 66 single-agent subjects plus 8 rituximab-combination subjects). Lu 2016 Supplemental Table 1.",
    covariates = "None. Lu 2016 states that a relatively simple model without covariates was preferred given the goals of the analysis, so no covariate assessment was performed; the only subject-level classifier in the model is which ADC the patient received.",
    notes = "Both ADCs link monomethyl auristatin E to the antibody through the same protease-labile MC-VC-PABC linker and share a mean drug-to-antibody ratio (mDAR) of 3.585 measured by hydrophobic interaction chromatography. Sampling in the phase I studies was intensive (more than 35 serum and plasma samples per patient for the single-agent cohorts, more than 30 for the rituximab-combination cohorts). Lu 2016 reports four fits of this one structural model: Model 1 (phase I data, dataset 1), and Models 2, 3 and 4, which refit it to progressively reduced phase II Tab sampling schemes in order to demonstrate that Tab sampling can be cut or eliminated. This file carries the Model 1 estimates, which are the ones Lu 2016 presents as the final integrated model developed from the intensive phase I data.",
    dosing_note = "Each ADC administration generates input to BOTH analytes. To simulate, provide TWO dose events per administration sharing the same amt, time and ii/addl: one with cmt = 'central_tab' and one with cmt = 'central', where amt is the molar antibody dose in nmol. The model applies f(central) <- mdar = 3.585 so the acMMAE compartment receives mDAR * DTab, reproducing the Lu 2016 Figure 2 initial conditions A1(0) = DTab and A3(0) = mDAR * DTab. Convert a mass dose to nmol with amt = dose_mg * 1e6 / 146455 using the antibody molecular weight carried in the supplement S1 control stream. Lu 2016 computed its published steady-state exposures with bolus dosing; use rate = 0 to reproduce them, or supply rate / dur for the clinical infusions."
  )

  ini({
    # ============================================================
    # Structural parameters, Lu 2016 Table 1, Model 1 column
    # (fit to the phase I dataset 1). Values are the final
    # estimates from the paper's results table, NOT the initial
    # estimates in the supplement S1 $THETA block.
    # Molecule-specific parameters carry a stratum suffix because
    # Lu 2016 estimates each of them twice in one joint fit.
    # ============================================================
    lcl_pina <- log(0.0292)
    label("Pinatuzumab vedotin proteolytic clearance (L/h)")
    lcl_pola <- log(0.0355)
    label("Polatuzumab vedotin proteolytic clearance (L/h)")
    lq <- log(0.0267)
    label("Intercompartmental clearance, shared by both ADCs (L/h)")
    lvc_pina <- log(5.35)
    label("Pinatuzumab vedotin central volume (L)")
    lvc_pola <- log(5.00)
    label("Polatuzumab vedotin central volume (L)")
    lvp <- log(8.21)
    label("Peripheral volume, shared by both ADCs (L)")
    lkdec_pina <- log(0.00855)
    label("Pinatuzumab vedotin deconjugation rate of a single MMAE molecule from the conjugate (1/h)")
    lkdec_pola <- log(0.00647)
    label("Polatuzumab vedotin deconjugation rate of a single MMAE molecule from the conjugate (1/h)")

    # Lu 2016 Eq. 5, Cpredicted = CORR * Cmodel. CORR reconciles the
    # mDAR measured in the dosing solution by HIC, which sets the
    # acMMAE initial condition, against the molar Tab:acMMAE ratio
    # implied by the ELISA and LC-MS/MS PK assays; without it the
    # model underestimated acMMAE and overestimated Tab at every
    # timepoint. Unitless multiplicative assay gain with no
    # intercept, hence cal_slope_<assay> and no cal_int_<assay>.
    cal_slope_acmmae_pina <- 1.31
    label("Pinatuzumab vedotin acMMAE assay cross-calibration slope (unitless)")
    cal_slope_acmmae_pola <- 1.45
    label("Polatuzumab vedotin acMMAE assay cross-calibration slope (unitless)")

    # ============================================================
    # Inter-individual variability, Lu 2016 Table 1, Model 1
    # column. The table reports variances (the column labels are
    # omega^2), which is also what the supplement S1 $OMEGA block
    # specifies, and each eta enters as EXP(ETA(n)) on the linear
    # parameter. Diagonal; Lu 2016 fits no OMEGA BLOCK. The IIV is
    # shared across both ADCs per Lu 2016 Methods.
    # ============================================================
    etalcl ~ 0.488 # Lu 2016 Table 1, row 'Random effect on CL: omega^2 CL' = 0.488
    etalq ~ 0.269 # Lu 2016 Table 1, row 'Random effect on Q: omega^2 Q' = 0.269
    etalvc ~ 0.0519 # Lu 2016 Table 1, row 'Random effect on VC: omega^2 VC' = 0.0519
    etalvp ~ 0.707 # Lu 2016 Table 1, row 'Random effect on VP: omega^2 VP' = 0.707
    etalkdec ~ 0.053 # Lu 2016 Table 1, row 'Random effect on kdec: omega^2 kdec' = 0.053

    # ============================================================
    # Residual error. Supplement S1 $ERROR applies Y = TY*(1+EPS)
    # separately to the two analytes, i.e. proportional in linear
    # space, and Lu 2016 Table 1 reports the EPS variances. The
    # SDs below are the square roots of those variances. Shared by
    # both ADCs per Lu 2016 Methods.
    # ============================================================
    propSd <- 0.177200
    label("acMMAE proportional residual SD (fraction)")
    propSd_tab <- 0.241868
    label("Total antibody proportional residual SD (fraction)")
  })
  model({
    # ------------------------------------------------------------
    # Constants carried in the Lu 2016 supplement S1 control stream.
    # mdar is the mean drug-to-antibody ratio of the dosing solution
    # measured by hydrophobic interaction chromatography, identical
    # for both ADCs (Lu 2016 Figure 2 legend and Methods).
    # mw_tab and mw_mmae convert a molar concentration in nmol/L to
    # a mass concentration in ng/mL; they are the $ERROR block's
    # coDose = 146455/1000 (antibody, g/mol) and the 0.718 factor
    # applied to A(3)/V1 (MMAE, g/mol).
    # ------------------------------------------------------------
    mdar <- 3.585 # unitless MMAE molecules per antibody
    mw_tab <- 146.455 # (ng/mL) per (nmol/L) for the antibody
    mw_mmae <- 0.718 # (ng/mL) per (nmol/L) for MMAE

    # ------------------------------------------------------------
    # Individual parameters. The molecule indicators select the
    # ADC-specific estimates; Q, Vp and every eta are shared.
    # ------------------------------------------------------------
    cl <- exp(lcl_pina * TRT_PINATUZUMAB_VEDOTIN + lcl_pola * TRT_POLATUZUMAB_VEDOTIN + etalcl)
    vc <- exp(lvc_pina * TRT_PINATUZUMAB_VEDOTIN + lvc_pola * TRT_POLATUZUMAB_VEDOTIN + etalvc)
    kdec <- exp(lkdec_pina * TRT_PINATUZUMAB_VEDOTIN + lkdec_pola * TRT_POLATUZUMAB_VEDOTIN + etalkdec)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)
    corr <- cal_slope_acmmae_pina * TRT_PINATUZUMAB_VEDOTIN +
      cal_slope_acmmae_pola * TRT_POLATUZUMAB_VEDOTIN

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ------------------------------------------------------------
    # ODE system, Lu 2016 Eqs. 1-4 (supplement S1 $DES). Both
    # analytes use the same k10, k12 and k21; deconjugation is the
    # only additional loss and acts on the acMMAE central
    # compartment alone. States hold molar amounts in nmol.
    # ------------------------------------------------------------
    d/dt(central) <- -(kel + k12) * central - kdec * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    d/dt(central_tab) <- -(kel + k12) * central_tab + k21 * peripheral1_tab
    d/dt(peripheral1_tab) <- k12 * central_tab - k21 * peripheral1_tab

    # Lu 2016 Figure 2 initial conditions: A1(0) = DTab and
    # A3(0) = mDAR * DTab. Dosing both analyte compartments with the
    # same molar antibody amount and scaling the acMMAE input by
    # mdar reproduces them.
    f(central) <- mdar

    # ------------------------------------------------------------
    # Observations, in ng/mL. Lu 2016 Eq. 5 applies the assay
    # cross-calibration slope to the acMMAE prediction. The paper's
    # text for Eq. 5 writes Cmodel = A1/VC, which is a typo for
    # A3/VC: A1 is the Tab central amount, and the supplement S1
    # $ERROR block computes ACMMAE = A(3)/V1*0.718 with
    # TY = CORR*ACMMAE for the acMMAE observation type. See the
    # vignette Errata section.
    # ------------------------------------------------------------
    Cc <- corr * mw_mmae * central / vc
    Cc_tab <- mw_tab * central_tab / vc

    Cc ~ prop(propSd)
    Cc_tab ~ prop(propSd_tab)
  })
}
