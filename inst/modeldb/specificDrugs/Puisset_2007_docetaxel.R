Puisset_2007_docetaxel <- function() {
  description <- "Friberg-style semi-mechanistic myelosuppression PD model for docetaxel-induced neutropenia in adult cancer patients (Puisset 2007). PD-only: docetaxel plasma concentration is supplied as the time-varying CP_MGL covariate (mg/L) and drives a linear drug effect E_drug = Slope * CP_MGL on the proliferating compartment. The five-compartment Friberg PD chain (one proliferating pool, three transit compartments, one circulating ANC compartment) and feedback (Circ0 / circ)^gamma reproduce the structure of Friberg LE et al. (2002) J Clin Oncol 20(24):4713-4721. Three covariates retained in the published final covariate model act multiplicatively on Slope: alpha-1 acid glycoprotein (AAG) as a power form (AAG / 1.29)^(-0.72), prior chemotherapy >= 2 lines (PRIOR_CHEMO_LINES_GE2) as a 1.69-fold multiplier, and treatment centre Toulouse vs Paris (STUDY_TOULOUSE) as a 1.82-fold multiplier (the centre effect is acknowledged by the authors to most likely reflect a between-centre HPLC-assay bias on docetaxel concentration rather than a clinical PD covariate). The upstream docetaxel PK was held fixed at Baille 1997 / Bruno 1996 individual posthoc profiles during the published PD fit; users couple this model with their preferred docetaxel popPK (e.g. modellib('Ozawa_2007_docetaxel') or modellib('Netterberg_2017_docetaxel')) to drive CP_MGL."
  reference <- paste(
    "Puisset F, Alexandre J, Treluyer J-M, Raoul V, Roche H, Goldwasser F,",
    "Chatelut E. (2007).",
    "Clinical pharmacodynamic factors in docetaxel toxicity.",
    "British Journal of Cancer 97(3):290-296.",
    "doi:10.1038/sj.bjc.6603872.",
    "PD structure follows Friberg LE et al. (2002) J Clin Oncol 20(24):4713-4721",
    "(see modellib('Friberg_2002_paclitaxel') for the leukocyte sibling).",
    "Upstream docetaxel PK was supplied via Bayesian individual posthoc profiles",
    "from Baille P et al. (1997) Clin Cancer Res 3:1535-1538 (built on the Bruno R et al.",
    "(1996) J Pharmacokinet Biopharm 24:153-172 docetaxel popPK structural model);",
    "this PD-only model therefore consumes plasma docetaxel concentration via the",
    "time-varying CP_MGL covariate rather than encoding the upstream PK explicitly.",
    sep = " "
  )
  vignette <- "Puisset_2007_docetaxel"

  units <- list(
    time          = "hour",
    dosing        = "mg",
    concentration = "mg/L",
    anc           = "10^9 cells/L",
    aag           = "g/L"
  )

  covariateData <- list(
    AAG = list(
      description        = "Serum alpha-1 acid glycoprotein concentration, time-fixed at baseline.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Multiplicative power-form effect on the drug-effect slope: slope = TVSlope * (AAG / 1.29)^e_aag_slope (Puisset 2007 Table 3 final model). Reference value 1.29 g/L is the cohort mean (Table 1: All patients n=92, mean 1.29, range 0.46-2.98 g/L). Higher AAG decreases Slope (lower neutropenia sensitivity per unit total plasma docetaxel) because docetaxel is highly bound to AAG and only the unbound fraction drives the pharmacodynamic effect. Note that increased AAG simultaneously decreases docetaxel clearance, so the net effect on ANC nadir is modest (Discussion, Figure 6A). Time-fixed per subject.",
      source_name        = "AAG"
    ),
    PRIOR_CHEMO_LINES_GE2 = list(
      description        = "Binary indicator: 1 = patient received at least two prior chemotherapy lines before docetaxel (3rd-line-or-later), 0 = 0 or 1 prior chemotherapy lines (1st- or 2nd-line).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (0 or 1 prior chemotherapy lines)",
      notes              = "Multiplicative effect on the drug-effect slope: slope = ... * e_prior_chemo_lines_ge2_slope^PRIOR_CHEMO_LINES_GE2 with e_prior_chemo_lines_ge2_slope = 1.69 (Puisset 2007 Table 3 final model). Patients with >= 2 prior chemotherapy lines have a 69% higher Slope than less-pretreated patients, consistent with cumulative bone-marrow depletion. Source NM-TRAN column `PTT2` (Methods, page 291: 'PTT2 = 0 or = 1 if patients had less than two lines, or at least two lines of chemotherapy before docetaxel, respectively'). Cohort distribution (Table 1, All patients n=92, 0 / 1 / >=2 prior lines = 27 / 44 / 21).",
      source_name        = "PTT2"
    ),
    STUDY_TOULOUSE = list(
      description        = "Binary indicator: 1 = patient enrolled at the Institut Claudius-Regaud (Toulouse, n=37), 0 = Hopital Cochin (Paris, n=55). Captures a between-centre analytical bias acknowledged by the authors rather than a clinical PD covariate.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Paris cohort)",
      notes              = "Multiplicative effect on the drug-effect slope: slope = ... * e_study_toulouse_slope^STUDY_TOULOUSE with e_study_toulouse_slope = 1.82 (Puisset 2007 Table 3 final model). Patients enrolled in Toulouse have an 82% higher Slope than the Paris cohort, attributed by the authors (Discussion, page 293) to an HPLC-assay underestimation of plasma docetaxel in Toulouse and a coarser PK sampling schedule (4 sampling times vs 3) that biased the Bayesian individual CL estimates upward. A simultaneous +70% CEN effect appears on docetaxel CL in the same paper (Results: CL = 31.3 * (AAG/1.29)^(-0.412) * 1.7^CEN, page 293). Source NM-TRAN column `CEN` (Methods, page 291: 'CEN = 0 or = 1 if data corresponded to Paris or Toulouse, respectively'). For typical-value simulation, leave STUDY_TOULOUSE = 0 (Paris reference) since the effect is a study-level artefact, not a generalisable PD covariate.",
      source_name        = "CEN"
    ),
    CP_MGL = list(
      description        = "Time-varying instantaneous docetaxel plasma concentration supplied per event row as the PD driver.",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Drug effect input: edrug = slope * CP_MGL with slope in 1/(mg/L) units. CP_MGL is supplied as a time-varying covariate column rather than computed from a coupled PK model, reflecting the Puisset 2007 sequential-PK-then-PD design (Methods, page 291: 'Individual-specific PK parameters were obtained by analysis of docetaxel plasma concentrations vs time using the Bayesian estimation method previously proposed by Baille et al (1997). Complete plasma docetaxel concentration vs time profiles used for the PK/PD analysis were generated from this analysis'). Population mean (range) docetaxel clearance is 40.0 (15.9-74.4) L/h and AUC is 4.1 (1.9-8.7) mg.h/L for a typical dose of 70-100 mg/m^2 IV over 1 h (Results, page 293). Set CP_MGL = 0 outside the drug-exposure window. The vignette validates the model by coupling it with the in-package modellib('Ozawa_2007_docetaxel') docetaxel PK to generate CP_MGL trajectories; users may substitute any other docetaxel popPK source.",
      source_name        = NA_character_
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 92L,
    n_studies      = 2L,
    age_range      = "46-77 years",
    age_median     = "60.7 years (mean across the pooled cohort)",
    weight_range   = "39-106 kg",
    weight_median  = "70.8 kg (mean across the pooled cohort)",
    sex_female_pct = 48.9,
    race_ethnicity = NULL,
    disease_state  = "Adult cancer patients receiving single-agent docetaxel monochemotherapy. Primary diseases (Table 1, All patients n=92): breast 31, prostate 27, lung 15, other 19. Performance status (ECOG) distribution 0 / 1 / 2 / 3 = 13 / 61 / 16 / 2.",
    dose_range     = "Docetaxel (Taxotere, Aventis Pharma) 70-100 mg/m^2 IV over 1 h, single-agent. Per-dose-level counts (Table 1, All patients): 70 mg/m^2 n=1; 75 mg/m^2 n=29; 85 mg/m^2 n=39; 100 mg/m^2 n=23. All patients received antiallergic and antiemetic premedication; no prophylactic G-CSF.",
    regions        = "France (Institut Claudius-Regaud, Toulouse, n=37; Hopital Cochin, Paris, n=55)",
    prior_chemo    = "Prior chemotherapy lines distribution (Table 1, All patients): 0 / 1 / >=2 = 27 / 44 / 21.",
    baseline_labs  = "Baseline characteristics (Table 1, All patients): age 60.7 (46-77) y; body weight 70.8 (39-106) kg; BSA 1.77 (1.35-2.20) m^2; alpha-1 acid glycoprotein 1.29 (0.46-2.98) g/L; serum albumin 36 (16-45) g/L. Liver-function-test elevation (>1.5xULN / <=1.5xULN): AST 12 / 80; ALT 8 / 84. Toulouse vs Paris differed significantly in AAG (1.51 vs 1.13 g/L, p<0.01), liver-function elevations (p<0.01), performance-status distribution (p<0.01), and docetaxel dose distribution (p<0.001) but not in age, body weight, BSA, serum albumin, sex ratio, prior-chemotherapy distribution, or primary-disease mix.",
    notes          = "Pooled cohort from two sequential single-centre clinical trials. Toulouse trial (n=37) correlated docetaxel CL with dexamethasone CL as a CYP3A probe (Puisset et al 2007, Cancer Chemother Pharmacol 54:265-272). Paris trial (n=55) related docetaxel toxicity to CYP3A / MDR1 / GST genetic polymorphisms (Tran et al 2006, Clin Pharmacol Ther 79:570-580). The two protocols differed in blood-sampling schedules for the PK (Toulouse: pre-dose + 0.5, 2, 6 h after end of infusion; Paris: pre-dose, end of infusion, 6 h after end of infusion) and in the docetaxel HPLC assay (interday CV: Toulouse 15.7%, Paris 6.0%). All 92 patients received weekly complete blood counts after the first cycle; no patient received prophylactic G-CSF. The PD analysis used a HYBRID NONMEM estimation method (FOCE for Circ0; FO for MTT and Slope)."
  )

  ini({
    # ------------------------------------------------------------------
    # Friberg myelosuppression structural PD parameters (Puisset 2007 Table 2 -- 'Present study' column).
    # Values are population means estimated by NONMEM (HYBRID method: FOCE for Circ0, FO for MTT and Slope).
    # ------------------------------------------------------------------
    lcirc0 <- log(4.41);  label("Baseline circulating ANC Circ0 (10^9 cells/L)")                                  # Table 2 'Present study': estimate 4.41, RSE 6.3%
    lmtt   <- log(96.2);  label("Mean transit time MTT through the proliferation -> circulation chain (h)")      # Table 2 'Present study': estimate 96.2, RSE 8.8%
    lgamma <- log(0.146); label("Feedback exponent gamma on (Circ0 / circ) (unitless)")                          # Table 2 'Present study': estimate 0.146, RSE 42.5%

    # ------------------------------------------------------------------
    # Drug-effect slope and its three retained covariate effects (Puisset 2007 Table 3 'Final model' row).
    # lslope encodes the typical-value Slope for a reference subject at AAG = 1.29 g/L (cohort mean),
    # PRIOR_CHEMO_LINES_GE2 = 0 (< 2 prior chemo lines), and STUDY_TOULOUSE = 0 (Paris reference).
    # The covariate effects multiply the typical Slope inside model().
    #
    # Important typo note (paper Errata): Table 3's final-model row prints
    # 'Slope = theta1 * (AGE / 1.29)^theta2 * theta3^PTT2 * theta4^CEN' but the reference value 1.29
    # is the cohort-mean AAG (Table 1), not AGE (cohort mean 60.7 years). The intermediate-model row
    # immediately above carries both (AGE / 60.7) and (AAG / 1.29) separately; the Analysis-of-
    # covariates text on page 293 states 'The absolute theta value corresponding to AGE was not
    # significantly different from 0. The final covariate model was based on PTT2, AAG, and CEN'.
    # The Discussion confirms 'the figure for patients treated at the Toulouse Centre was 82%, and
    # was nearly inversely proportional to the AAG level'. The 'AGE' in the Table 3 final-model row
    # is therefore a typographical error for AAG; this model encodes the intended (AAG / 1.29)^theta2
    # form. The vignette Errata reproduces this audit.
    # ------------------------------------------------------------------
    lslope <- log(7.40);  label("Typical-value drug-effect slope at AAG=1.29 g/L, PRIOR_CHEMO_LINES_GE2=0, STUDY_TOULOUSE=0 (1/(mg/L))")  # Table 3 final model: theta1 = 7.40 (+/- 1.22)

    e_aag_slope                   <- -0.72; label("Power-form exponent of AAG/1.29 on Slope (unitless, negative)")  # Table 3 final model: theta2 = -0.72 (+/- 0.18)
    e_prior_chemo_lines_ge2_slope <-  1.69; label("Multiplicative effect of PRIOR_CHEMO_LINES_GE2 on Slope (unitless)")  # Table 3 final model: theta3 = 1.69 (+/- 0.32)
    e_study_toulouse_slope        <-  1.82; label("Multiplicative effect of STUDY_TOULOUSE on Slope (unitless; HPLC-assay between-centre bias)")  # Table 3 final model: theta4 = 1.82 (+/- 0.46)

    # ------------------------------------------------------------------
    # Inter-individual variability (Puisset 2007 Methods page 291: 'The interindividual variability
    # was estimated for Circ0, MTT, and Slope according to an exponential model. A log-normal
    # distribution of the parameters was assumed.'). No IIV on gamma. Reported as CV%; converted
    # to log-eta variance via omega^2 = log(1 + CV^2):
    #   Circ0 IIV  49.2%  -> log(1 + 0.492^2) = 0.2169
    #   MTT   IIV  25.5%  -> log(1 + 0.255^2) = 0.0630
    #   Slope IIV  44.0%  -> log(1 + 0.440^2) = 0.1770  (final-model IIV, Table 3 -- reduced from
    #                                                    the structural-only 96.5% in Table 2 by
    #                                                    the AAG + PRIOR_CHEMO_LINES_GE2 +
    #                                                    STUDY_TOULOUSE covariate effects)
    # ------------------------------------------------------------------
    etalcirc0 ~ 0.2169   # Table 2 'Present study': IIV Circ0 = 49.2% CV (RSE 18.3%)
    etalmtt   ~ 0.0630   # Table 2 'Present study': IIV MTT   = 25.5% CV (RSE 40.2%)
    etalslope ~ 0.1770   # Table 3 final model:    IIV Slope = 44%   CV

    # ------------------------------------------------------------------
    # Residual error (Puisset 2007 Methods page 291: 'A proportional model for residual variability
    # was used.'). Reported as CV%; mapped to the propSd fraction directly.
    # Table 2's 'Additive residual error (ANC 10^9 / L)' row reads 'Not applicable' for the present
    # study (additive component used by Friberg 2002 was not retained).
    # ------------------------------------------------------------------
    propSd <- 0.356;     label("Proportional residual SD on ANC (fraction)")                                     # Table 2 'Present study': proportional residual error = 35.6% (RSE 21.3%)
  })

  model({
    # ------------------------------------------------------------------
    # Individual PD parameters
    # ------------------------------------------------------------------
    circ0 <- exp(lcirc0 + etalcirc0)
    mtt   <- exp(lmtt   + etalmtt)
    gamma <- exp(lgamma)

    # Slope carries IIV plus three covariate effects (Table 3 final model).
    slope <-
      exp(lslope + etalslope) *
      (AAG / 1.29)^e_aag_slope *
      e_prior_chemo_lines_ge2_slope^PRIOR_CHEMO_LINES_GE2 *
      e_study_toulouse_slope^STUDY_TOULOUSE

    # Transit-rate constant: 3 transit compartments plus 1 proliferation pool, so
    # ktr = (n_transit + 1) / MTT = 4 / MTT (Friberg 2002 convention; same as
    # Friberg_2002_paclitaxel and Netterberg_2017_docetaxel).
    ktr <- 4 / mtt

    # ------------------------------------------------------------------
    # Drug effect (linear): E_drug = Slope * Conc (Puisset 2007 Methods, Figure 1 caption).
    # CP_MGL is the time-varying docetaxel plasma concentration (mg/L), supplied per event row.
    # ------------------------------------------------------------------
    edrug <- slope * CP_MGL

    # ------------------------------------------------------------------
    # Friberg myelosuppression chain
    # Compartment ordering follows the convention used by Friberg_2002_paclitaxel and
    # Netterberg_2017_docetaxel:
    #   circ        = Circ (observed circulating ANC)        (CMT 1 by ODE order)
    #   precursor1  = Prol (self-renewing proliferating pool) (CMT 2)
    #   precursor2  = Transit 1                              (CMT 3)
    #   precursor3  = Transit 2                              (CMT 4)
    #   precursor4  = Transit 3                              (CMT 5)
    # The proliferation equation
    #   dProl/dt = ktr * Prol * (1 - E_drug) * (Circ0 / Circ)^gamma - ktr * Prol
    # combines the self-renewal term, the linear cell-loss term (1 - E_drug), and the feedback
    # amplification (Circ0 / Circ)^gamma. Transit compartments are first-order with rate ktr
    # (Methods page 291: 'the rate constants were assumed to be equal: kprol, ktr, and kcirc').
    # ------------------------------------------------------------------
    d/dt(circ)       <-  ktr * precursor4 - ktr * circ
    d/dt(precursor1) <-  ktr * precursor1 * (1 - edrug) * (circ0 / circ)^gamma - ktr * precursor1
    d/dt(precursor2) <-  ktr * precursor1 - ktr * precursor2
    d/dt(precursor3) <-  ktr * precursor2 - ktr * precursor3
    d/dt(precursor4) <-  ktr * precursor3 - ktr * precursor4

    # Initial conditions: all five chain compartments start at the per-subject baseline ANC (circ0),
    # the steady-state assumption in the absence of drug (Methods page 291).
    circ(0)       <- circ0
    precursor1(0) <- circ0
    precursor2(0) <- circ0
    precursor3(0) <- circ0
    precursor4(0) <- circ0

    # ------------------------------------------------------------------
    # Observation: absolute neutrophil count ANC (10^9 cells/L), proportional residual error
    # ------------------------------------------------------------------
    ANC <- circ
    ANC ~ prop(propSd)
  })
}
