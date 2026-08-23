vanderGaag_2025_osimertinib_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body target-site, 26 states, hand-coded deSolve). Osimertinib",
    "disposition and EGFR target engagement in advanced-stage NSCLC, covering",
    "both a 2.7 ug intravenous [11C]C-osimertinib microdose and an 80 mg oral",
    "therapeutic dose. Thirteen perfusion-limited tissues (adipose, brain, gut,",
    "heart, kidney, liver, lung, lung tumor, pancreas, muscle, skin, spleen and",
    "a lumped extracellular-water rest-of-body) exchange with venous and",
    "arterial blood. The lung and the lung tumor are perfused in parallel from",
    "venous blood and drain into arterial blood; gut and spleen drain into the",
    "liver through the portal vein. Tissue-to-plasma partition coefficients are",
    "not fitted: they are computed inside the model from the Rodgers and Rowland",
    "moderate-to-strong-base equations, so each Kp emerges from the drug's pKa,",
    "logD, fraction unbound and blood-to-plasma ratio together with the tissue's",
    "fractional water, neutral-lipid, neutral-phospholipid, acidic-phospholipid",
    "and lysosomal composition. Lysosomal sequestration is resolved per cell",
    "type in liver and lung, and the tumor carries an acidified extracellular",
    "water (pH 6.7) plus a 100 percent residual-cell (immune-deprived)",
    "composition. Superimposed on this backbone is saturable 1:1 EGFR binding in",
    "the ten EGFR-expressing tissues: bound drug is carried as its own state per",
    "tissue, association is proportional to the free-EGFR fraction, and",
    "wild-type on/off rates are used everywhere except the tumor, which uses the",
    "50/50 wild-type / L858R-T790M mixture. Elimination is hepatic (well-stirred",
    "intrinsic clearance on unbound drug) plus renal. This is a bottom-up model:",
    "nothing was estimated from the clinical data, so there is no between-subject",
    "variability and the paper reports no residual-error model -- the propSd term",
    "is a placeholder. See the vignette Errata."
  )
  reference <- paste(
    "van der Gaag S, Jordens T, Yaqub M, Grijseels R, van Valkengoed DW,",
    "de Langen EN, van den Broek R, Thijssen VLJL, de Langen AJ,",
    "Kouwenhoven MCM, Bahce I, Westerman BA, Hendrikse NH, Bartelink IH.",
    "Physiologically Based Pharmacokinetic Model of Tyrosine Kinase Inhibitors",
    "to Predict Target Site Penetration, with PET-Guided Verification.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(5):918-928.",
    "doi:10.1002/psp4.70006.",
    "Organ volumes, intracellular-water volumes and blood flows from",
    "Appendix S2 Table S1; body/pH/haematocrit constants from Table S2;",
    "tissue fractional volumes from Table S3; lysosomal fractions, acidic",
    "phospholipid and EGFR concentrations from Table S4; liver and lung cell-type",
    "composition from Table S5; osimertinib physicochemical and clearance",
    "parameters from Appendix S3 Table S6; tumor and EGFR binding parameters from",
    "main-text Table 1; Kpu equations from Appendix S4 Suppl Eq 1-6; whole-body",
    "ODEs from Appendix S5 Suppl Eq 7-20; EGFR-binding and free-receptor-fraction",
    "equations from main-text Equations 1-2.",
    "The tumor extracellular-water pH term and the tumor 100% residual-cell",
    "composition are inherited from the explicitly cited predecessor model",
    "Bartelink IH, et al. Physiologically Based Pharmacokinetic (PBPK) Modeling",
    "to Predict PET Image Quality of Three Generations EGFR TKI in",
    "Advanced-Stage NSCLC Patients. Pharmaceuticals (Basel). 2022;15(7):796.",
    "doi:10.3390/ph15070796 (Supplement section VC and supplement note 9).",
    sep = " "
  )
  vignette <- "vanderGaag_2025_osimertinib_pbpk"
  units    <- list(time = "h", dosing = "nmol", concentration = "nM")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Every state holds an amount of osimertinib in nmol.
  # The `<organ>_complex` states hold osimertinib covalently bound to EGFR in
  # the intracellular water of that organ (1:1 stoichiometry, main text
  # section 2.2.2). verified = TRUE: the state list is transcribed from
  # Appendix S5 Suppl Eq 7-20, which writes out every differential equation.
  compartmentData <- list(
    gut_lumen        = list(analyte = "osimertinib", units = "nmol", specimen = "administration site", verified = TRUE),
    venous           = list(analyte = "osimertinib", units = "nmol", specimen = "whole blood", verified = TRUE),
    arterial         = list(analyte = "osimertinib", units = "nmol", specimen = "whole blood", verified = TRUE),
    brain            = list(analyte = "osimertinib", units = "nmol", specimen = "tissue", verified = TRUE),
    spleen           = list(analyte = "osimertinib", units = "nmol", specimen = "tissue", verified = TRUE),
    other            = list(analyte = "osimertinib", units = "nmol", specimen = "tissue", verified = TRUE),
    adipose          = list(analyte = "osimertinib", units = "nmol", specimen = "tissue", verified = TRUE),
    adipose_complex  = list(analyte = "osimertinib-EGFR complex", units = "nmol", specimen = "tissue", verified = TRUE),
    gut              = list(analyte = "osimertinib", units = "nmol", specimen = "tissue", verified = TRUE),
    gut_complex      = list(analyte = "osimertinib-EGFR complex", units = "nmol", specimen = "tissue", verified = TRUE),
    heart            = list(analyte = "osimertinib", units = "nmol", specimen = "tissue", verified = TRUE),
    heart_complex    = list(analyte = "osimertinib-EGFR complex", units = "nmol", specimen = "tissue", verified = TRUE),
    kidney           = list(analyte = "osimertinib", units = "nmol", specimen = "tissue", verified = TRUE),
    kidney_complex   = list(analyte = "osimertinib-EGFR complex", units = "nmol", specimen = "tissue", verified = TRUE),
    liver            = list(analyte = "osimertinib", units = "nmol", specimen = "tissue", verified = TRUE),
    liver_complex    = list(analyte = "osimertinib-EGFR complex", units = "nmol", specimen = "tissue", verified = TRUE),
    lung             = list(analyte = "osimertinib", units = "nmol", specimen = "tissue", verified = TRUE),
    lung_complex     = list(analyte = "osimertinib-EGFR complex", units = "nmol", specimen = "tissue", verified = TRUE),
    tumor            = list(analyte = "osimertinib", units = "nmol", specimen = "tumor", verified = TRUE),
    tumor_complex    = list(analyte = "osimertinib-EGFR complex", units = "nmol", specimen = "tumor", verified = TRUE),
    pancreas         = list(analyte = "osimertinib", units = "nmol", specimen = "tissue", verified = TRUE),
    pancreas_complex = list(analyte = "osimertinib-EGFR complex", units = "nmol", specimen = "tissue", verified = TRUE),
    muscle           = list(analyte = "osimertinib", units = "nmol", specimen = "tissue", verified = TRUE),
    muscle_complex   = list(analyte = "osimertinib-EGFR complex", units = "nmol", specimen = "tissue", verified = TRUE),
    skin             = list(analyte = "osimertinib", units = "nmol", specimen = "tissue", verified = TRUE),
    skin_complex     = list(analyte = "osimertinib-EGFR complex", units = "nmol", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 4L,
    n_studies      = 2L,
    age_range      = NA_character_,
    weight_range   = "single 87.5 kg reference individual (Table S2); the model is deterministic and is not stratified by body size",
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = paste(
      "advanced-stage (stage 3-4) non-small cell lung cancer with an activating",
      "EGFR mutation (exon 19 deletion or exon 21 L858R) and progression on a",
      "first-generation EGFR-TKI; the tumor compartment assumes a 50/50 mix of",
      "wild-type and L858R/T790M-mutated EGFR alleles"
    ),
    dose_range     = paste(
      "2.7 ug (280 MBq) intravenous [11C]C-osimertinib microdose, and 80 mg",
      "oral osimertinib once daily"
    ),
    regions        = "Netherlands (Amsterdam UMC)",
    notes          = paste(
      "This is a bottom-up whole-body PBPK model: every parameter is a",
      "literature or in-vitro value and none was estimated from the clinical",
      "data, so the model has no between-subject variability. Verification used",
      "two independent data sets. (1) Microdose: four NSCLC patients each",
      "underwent two dynamic [11C]C-osimertinib PET scans (60 or 120 min) within",
      "a two-week window, one on and one off concomitant first- or",
      "second-generation TKI; because concomitant TKI did not change uptake the",
      "two scans were pooled (n = 8 dynamic scans). Three of the four also had a",
      "whole-body static PET at 60 min, giving brain, kidney, liver, lung, spleen",
      "and tumor concentrations. Study CCMO NL64722.031.18. (2) Therapeutic:",
      "digitised plasma concentration-time profiles from 10 healthy volunteers",
      "given a single 80 mg oral dose and at steady state, extracted from phase 1",
      "study D5160C00001 with WebPlotDigitizer 4.6. Spherical tumor volumes came",
      "from the median baseline diameters of 12 NSCLC lesions in four patients.",
      "The authors' acceptance criterion was a 3-fold prediction-error window",
      "rather than the usual 2-fold, justified by the small sample size."
    )
  )

  ini({
    # =================================================================
    # All parameters are fixed(). This is a bottom-up PBPK model: the
    # paper estimated nothing from the clinical data (main text section
    # 2.1, "developed using a bottom-up approach"). Values are the
    # paper's tabulated model inputs.
    #
    # The paper works in minutes; rate constants and flows are converted
    # to hours here by an explicit "* 60" so the source-trace comment can
    # quote the published per-minute number unchanged.
    # =================================================================

    # ---- Osimertinib physicochemical properties (Appendix S3 Table S6) ----
    pka <- fixed(9.0)
    label("Acid dissociation constant (unitless)")                                       # Table S6: pKa 9.0
    logd <- fixed(3.2)
    label("Octanol:water distribution coefficient at pH 7.4, log10 (unitless)")          # Table S6: LogD 3.2
    logpv <- fixed(2.22)
    label("Vegetable oil:water partition coefficient, log10, adipose only (unitless)")   # Table S6: LogPV 2.22 = 1.115*LogD - 1.35
    fu <- fixed(0.0535)
    label("Fraction unbound in plasma (fraction)")                                       # Table S6: Fu 5.35%
    bp <- fixed(0.77)
    label("Blood-to-plasma concentration ratio (unitless)")                              # Table S6: bp 0.77 = 1 / plasma-to-blood 1.30
    kabc <- fixed(38.5)
    label("Association constant of basic compounds to acidic phospholipids (unitless)")  # Table S6: Ka_bc 38.5

    # ---- Absorption, distribution and elimination (Table S6) ----
    lka <- fixed(log(4.0e-3 * 60))
    label("First-order absorption rate constant from gut lumen (1/h)")                   # Table S6: Ka 4.0e-3 /min
    lfdepot <- fixed(log(0.70))
    label("Oral bioavailability (fraction)")                                             # Table S6: F 70%
    lclint <- fixed(log(3.63 * 60))
    label("Hepatic intrinsic clearance, unbound (L/h)")                                  # Table S6: CL_int 3.63 L/min (Suppl Eq 16, from CL_H = 0.17 L/min = 86% of total)
    lclr <- fixed(log(392e-5 * 60))
    label("Renal clearance (L/h)")                                                       # Table S6: CL_R 392e-5 L/min

    # ---- EGFR binding kinetics (main-text Table 1) ----
    # Wild-type rates are used in every healthy EGFR-expressing tissue;
    # the tumor uses the 50/50 wild-type / L858R-T790M mixture, which is
    # the arithmetic mean of the two pure sets (0.174/0.186 -> 0.18 and
    # 162/10.2 -> 86.1), consistent with a 50% mutated allele frequency.
    lkon <- fixed(log(0.174 * 60))
    label("EGFR association rate constant, wild-type (1/(nM*h))")                        # Table 1: Kon WT EGFR 0.174 /nM/min
    lkoff <- fixed(log(162 * 60))
    label("EGFR dissociation rate constant, wild-type (1/h)")                            # Table 1: Koff WT EGFR 162 /min
    lkon_tumor <- fixed(log(0.18 * 60))
    label("EGFR association rate constant, tumor 50/50 WT/L858R-T790M (1/(nM*h))")       # Table 1: Kon,tumor 0.18 /nM/min
    lkoff_tumor <- fixed(log(86.1 * 60))
    label("EGFR dissociation rate constant, tumor 50/50 WT/L858R-T790M (1/h)")           # Table 1: Koff,tumor 86.1 /min

    # ---- pH of the body-water spaces (Appendix S2 Table S2) ----
    ph_plasma <- fixed(7.4)
    label("pH of plasma (unitless)")                                                     # Table S2: pHp 7.4
    ph_iw <- fixed(7.0)
    label("pH of intracellular water (unitless)")                                        # Table S2: pHiw 7.0
    ph_ew_tumor <- fixed(6.7)
    label("pH of tumor extracellular water (unitless)")                                  # Table 1: PH EW 6.7 (healthy tissue uses plasma pH, Table S2 pHew 7.4)

    # ---- Organ total volumes, L (Appendix S2 Table S1) ----
    lv_venous <- fixed(log(5.6 * 2 / 3))
    label("Venous blood volume (L)")                                                     # Table S1: Blood 5.6 L; 2/3 venous per the Jones 2006 framework (split not reported - see vignette Errata)
    lv_arterial <- fixed(log(5.6 * 1 / 3))
    label("Arterial blood volume (L)")                                                   # Table S1: Blood 5.6 L; 1/3 arterial per the Jones 2006 framework (split not reported - see vignette Errata)
    lv_adipose <- fixed(log(10.4))
    label("Adipose volume (L)")                                                          # Table S1: Adipose 10.4
    lv_brain <- fixed(log(1.49))
    label("Brain volume (L)")                                                            # Table S1: Brain 1.49
    lv_gut <- fixed(log(1.5))
    label("Gut volume (L)")                                                              # Table S1: Gut 1.5
    lv_heart <- fixed(log(0.31))
    label("Heart volume (L)")                                                            # Table S1: Heart 0.31
    lv_kidney <- fixed(log(0.31))
    label("Kidney volume (L)")                                                           # Table S1: Kidney 0.31
    lv_liver <- fixed(log(2.52))
    label("Liver volume (L)")                                                            # Table S1: Liver 2.52
    lv_lung <- fixed(log(0.87))
    label("Lung volume (L)")                                                             # Table S1: Lung 0.87
    lv_tumor <- fixed(log(0.05))
    label("Lung tumor virtual volume (L)")                                               # Table S1: Lung-tumor (virtual volume) 0.05; Table 1 virtual tumor volume 5% of lung, = 3% tumor x 1.7-fold cell density
    lv_pancreas <- fixed(log(0.21))
    label("Pancreas volume (L)")                                                         # Table S1: Pancreas 0.21
    lv_muscle <- fixed(log(33.9))
    label("Muscle volume (L)")                                                           # Table S1: Muscle 33.9
    lv_skin <- fixed(log(5.63))
    label("Skin volume (L)")                                                             # Table S1: Skin 5.63
    lv_spleen <- fixed(log(0.20))
    label("Spleen volume (L)")                                                           # Table S1: Spleen 0.20
    lv_other <- fixed(log(15.4))
    label("Lumped rest-of-body extracellular-water volume (L)")                          # Table S1: EW 15.4 (weight minus the sum of all total tissue volumes)

    # ---- Organ intracellular-water volumes, L (Table S1) ----
    # Only the EGFR-expressing tissues need an IW volume: it converts the
    # tissue EGFR concentration to the concentration in the water space
    # where binding happens (main-text Equation 1).
    lviw_adipose <- fixed(log(0.18))
    label("Adipose intracellular-water volume (L)")                                      # Table S1: Adipose IW 0.18
    lviw_gut <- fixed(log(0.71))
    label("Gut intracellular-water volume (L)")                                          # Table S1: Gut IW 0.71
    lviw_heart <- fixed(log(0.14))
    label("Heart intracellular-water volume (L)")                                        # Table S1: Heart IW 0.14
    lviw_kidney <- fixed(log(0.14))
    label("Kidney intracellular-water volume (L)")                                       # Table S1: Kidney IW 0.14
    lviw_liver <- fixed(log(1.41))
    label("Liver intracellular-water volume (L)")                                        # Table S1: Liver IW 1.41
    lviw_lung <- fixed(log(0.38))
    label("Lung intracellular-water volume (L)")                                         # Table S1: Lung IW 0.38
    lviw_tumor <- fixed(log(0.02))
    label("Lung tumor intracellular-water volume (L)")                                   # Table S1: Lung-tumor IW 0.02
    lviw_pancreas <- fixed(log(0.14))
    label("Pancreas intracellular-water volume (L)")                                     # Table S1: Pancreas IW 0.14
    lviw_muscle <- fixed(log(21.4))
    label("Muscle intracellular-water volume (L)")                                       # Table S1: Muscle IW 21.4
    lviw_skin <- fixed(log(1.64))
    label("Skin intracellular-water volume (L)")                                         # Table S1: Skin IW 1.64

    # ---- Organ blood flows, L/min -> L/h (Table S1) ----
    lq_adipose <- fixed(log(0.33 * 60))
    label("Adipose blood flow (L/h)")                                                    # Table S1: Adipose Q 0.33 L/min
    lq_brain <- fixed(log(0.78 * 60))
    label("Brain blood flow (L/h)")                                                      # Table S1: Brain Q 0.78 L/min
    lq_gut <- fixed(log(0.98 * 60))
    label("Gut blood flow (L/h)")                                                        # Table S1: Gut Q 0.98 L/min
    lq_heart <- fixed(log(0.26 * 60))
    label("Heart blood flow (L/h)")                                                      # Table S1: Heart Q 0.26 L/min
    lq_hepatic_artery <- fixed(log(0.42 * 60))
    label("Hepatic artery blood flow (L/h)")                                             # Table S1: Hepatic artery Q 0.42 L/min
    lq_kidney <- fixed(log(1.24 * 60))
    label("Kidney blood flow (L/h)")                                                     # Table S1: Kidney Q 1.24 L/min
    lq_liver <- fixed(log(1.59 * 60))
    label("Total hepatic blood flow (L/h)")                                              # Table S1: Liver Q 1.59 L/min (= gut 0.98 + spleen 0.20 + hepatic artery 0.42)
    lq_lung <- fixed(log(6.31 * 60))
    label("Lung blood flow (L/h)")                                                       # Table S1: Lung Q 6.31 L/min
    lq_tumor <- fixed(log(0.20 * 60))
    label("Lung tumor blood flow (L/h)")                                                 # Table S1: Lung-tumor Q 0.20 L/min; Table 1 tumor blood flow 3% of lung
    lq_pancreas <- fixed(log(0.32 * 60))
    label("Pancreas blood flow (L/h)")                                                   # Table S1: Pancreas Q 0.32 L/min
    lq_muscle <- fixed(log(1.11 * 60))
    label("Muscle blood flow (L/h)")                                                     # Table S1: Muscle Q 1.11 L/min
    lq_skin <- fixed(log(0.33 * 60))
    label("Skin blood flow (L/h)")                                                       # Table S1: Skin Q 0.33 L/min
    lq_spleen <- fixed(log(0.20 * 60))
    label("Spleen blood flow (L/h)")                                                     # Table S1: Spleen Q 0.20 L/min
    lq_other <- fixed(log(0.24 * 60))
    label("Lumped rest-of-body blood flow (L/h)")                                        # Table S1: EW Q 0.24 L/min

    # ---- Tissue EGFR concentrations, nM (Appendix S2 Table S4) ----
    # Expressed per unit of TOTAL tissue volume (main-text Equation 3);
    # brain, spleen and blood cells scored "Not Detected" in the Human
    # Protein Atlas and carry no EGFR, so they have no binding states.
    legfr_adipose <- fixed(log(1.82))
    label("Adipose EGFR concentration (nM)")                                             # Table S4: Adipose EGFR 1.82
    legfr_gut <- fixed(log(40.6))
    label("Gut EGFR concentration (nM)")                                                 # Table S4: Gut EGFR 40.6
    legfr_heart <- fixed(log(48.4))
    label("Heart EGFR concentration (nM)")                                               # Table S4: Heart EGFR 48.4
    legfr_kidney <- fixed(log(205))
    label("Kidney EGFR concentration (nM)")                                              # Table S4: Kidney EGFR 205
    legfr_liver <- fixed(log(235))
    label("Liver EGFR concentration (nM)")                                               # Table S4: Liver EGFR 235
    legfr_lung <- fixed(log(1616))
    label("Lung EGFR concentration (nM)")                                                # Table S4: Lung EGFR 1616
    legfr_tumor <- fixed(log(382))
    label("Lung tumor EGFR concentration (nM)")                                          # Table S4: Lung-Tumor EGFR 382 (main-text Table 1 rounds to 383)
    legfr_pancreas <- fixed(log(281))
    label("Pancreas EGFR concentration (nM)")                                            # Table S4: Pancreas EGFR 281
    legfr_muscle <- fixed(log(67.7))
    label("Muscle EGFR concentration (nM)")                                              # Table S4: Muscle EGFR 67.7
    legfr_skin <- fixed(log(120))
    label("Skin EGFR concentration (nM)")                                                # Table S4: Skin EGFR 120

    # ---- Residual error ----
    # The paper fits nothing and reports no residual-error model, but an
    # nlmixr2 model definition requires a residual-error term. propSd is a
    # fixed placeholder for syntactic completeness only and must NOT be
    # read as an estimate. Same convention as Aoki_2024_bosentan_pbpk,
    # Mi_2023_cefquinome_pbpk and An_2012_mitoxantrone_*_pbpk.
    propSd <- fixed(0.10)
    label("Proportional residual error placeholder (fraction)")                          # not reported in van der Gaag 2025; placeholder only
  })

  model({
    # =================================================================
    # 1. Back-transform the log-scale parameters
    # =================================================================
    ka        <- exp(lka)
    clint     <- exp(lclint)
    clr       <- exp(lclr)
    kon       <- exp(lkon)
    koff      <- exp(lkoff)
    kon_t     <- exp(lkon_tumor)
    koff_t    <- exp(lkoff_tumor)

    v_venous   <- exp(lv_venous)
    v_arterial <- exp(lv_arterial)
    v_adipose  <- exp(lv_adipose)
    v_brain    <- exp(lv_brain)
    v_gut      <- exp(lv_gut)
    v_heart    <- exp(lv_heart)
    v_kidney   <- exp(lv_kidney)
    v_liver    <- exp(lv_liver)
    v_lung     <- exp(lv_lung)
    v_tumor    <- exp(lv_tumor)
    v_pancreas <- exp(lv_pancreas)
    v_muscle   <- exp(lv_muscle)
    v_skin     <- exp(lv_skin)
    v_spleen   <- exp(lv_spleen)
    v_other    <- exp(lv_other)

    viw_adipose  <- exp(lviw_adipose)
    viw_gut      <- exp(lviw_gut)
    viw_heart    <- exp(lviw_heart)
    viw_kidney   <- exp(lviw_kidney)
    viw_liver    <- exp(lviw_liver)
    viw_lung     <- exp(lviw_lung)
    viw_tumor    <- exp(lviw_tumor)
    viw_pancreas <- exp(lviw_pancreas)
    viw_muscle   <- exp(lviw_muscle)
    viw_skin     <- exp(lviw_skin)

    q_adipose  <- exp(lq_adipose)
    q_brain    <- exp(lq_brain)
    q_gut      <- exp(lq_gut)
    q_heart    <- exp(lq_heart)
    q_ha       <- exp(lq_hepatic_artery)
    q_kidney   <- exp(lq_kidney)
    q_liver    <- exp(lq_liver)
    q_lung     <- exp(lq_lung)
    q_tumor    <- exp(lq_tumor)
    q_pancreas <- exp(lq_pancreas)
    q_muscle   <- exp(lq_muscle)
    q_skin     <- exp(lq_skin)
    q_spleen   <- exp(lq_spleen)
    q_other    <- exp(lq_other)

    egfr_adipose  <- exp(legfr_adipose)
    egfr_gut      <- exp(legfr_gut)
    egfr_heart    <- exp(legfr_heart)
    egfr_kidney   <- exp(legfr_kidney)
    egfr_liver    <- exp(legfr_liver)
    egfr_lung     <- exp(legfr_lung)
    egfr_tumor    <- exp(legfr_tumor)
    egfr_pancreas <- exp(legfr_pancreas)
    egfr_muscle   <- exp(legfr_muscle)
    egfr_skin     <- exp(legfr_skin)

    # =================================================================
    # 2. Tissue composition -- Rodgers and Rowland fractional volumes
    #
    # Appendix S2 Table S3 (fnl, fnp, few, fiw), Table S4 (flys, acidic
    # phospholipid concentration AP in mg/g) and Table S5 (liver and lung
    # cell-type composition). These are drug-independent physiological
    # constants of the Rodgers and Rowland tissue-composition scheme, so
    # they live here rather than in ini().
    #
    # Column order below is: fnl, fnp, few, fiw.
    # =================================================================
    # Neutral lipids
    fnl_adipose  <- 0.85     # Table S3
    fnl_brain    <- 0.04     # Table S3
    fnl_gut      <- 0.04     # Table S3
    fnl_heart    <- 0.01     # Table S3
    fnl_kidney   <- 0.04     # Table S3
    fnl_liver    <- 0.03     # Table S3
    fnl_lung     <- 0.009    # Table S3
    fnl_pancreas <- 0.04     # Table S3
    fnl_muscle   <- 0.02     # Table S3
    fnl_skin     <- 0.06     # Table S3
    fnl_spleen   <- 0.02     # Table S3
    # Neutral phospholipids
    fnp_adipose  <- 0.002    # Table S3
    fnp_brain    <- 0.002    # Table S3
    fnp_gut      <- 0.01     # Table S3
    fnp_heart    <- 0.01     # Table S3
    fnp_kidney   <- 0.01     # Table S3
    fnp_liver    <- 0.02     # Table S3
    fnp_lung     <- 0.003    # Table S3
    fnp_pancreas <- 0.009    # Table S3
    fnp_muscle   <- 0.007    # Table S3
    fnp_skin     <- 0.004    # Table S3
    fnp_spleen   <- 0.02     # Table S3
    # Extracellular water
    few_adipose  <- 0.14     # Table S3
    few_brain    <- 0.16     # Table S3
    few_gut      <- 0.28     # Table S3
    few_heart    <- 0.32     # Table S3
    few_kidney   <- 0.27     # Table S3
    few_liver    <- 0.16     # Table S3
    few_lung     <- 0.34     # Table S3
    few_pancreas <- 0.12     # Table S3
    few_muscle   <- 0.12     # Table S3
    few_skin     <- 0.38     # Table S3
    few_spleen   <- 0.21     # Table S3
    # Intracellular water
    fiw_adipose  <- 0.02     # Table S3
    fiw_brain    <- 0.61     # Table S3
    fiw_gut      <- 0.48     # Table S3
    fiw_heart    <- 0.46     # Table S3
    fiw_kidney   <- 0.47     # Table S3
    fiw_liver    <- 0.56     # Table S3
    fiw_lung     <- 0.43     # Table S3
    fiw_pancreas <- 0.66     # Table S3
    fiw_muscle   <- 0.63     # Table S3
    fiw_skin     <- 0.29     # Table S3
    fiw_spleen   <- 0.53     # Table S3
    # Lysosomal fractional volume
    flys_adipose  <- 2.7e-4  # Table S4
    flys_brain    <- 0.01    # Table S4
    flys_gut      <- 0.01    # Table S4
    flys_heart    <- 0.003   # Table S4
    flys_kidney   <- 0.02    # Table S4
    flys_pancreas <- 0.013   # Table S4 (not measured; mean of all measured organs)
    flys_muscle   <- 0.001   # Table S4
    flys_skin     <- 0.002   # Table S4
    flys_spleen   <- 0.05    # Table S4
    # Acidic phospholipids, mg/g
    ap_adipose  <- 0.4       # Table S4
    ap_brain    <- 0.4       # Table S4
    ap_gut      <- 2.41      # Table S4
    ap_heart    <- 3.15      # Table S4
    ap_kidney   <- 2.45      # Table S4
    ap_liver    <- 4.56      # Table S4
    ap_lung     <- 0.57      # Table S4
    ap_pancreas <- 1.67      # Table S4
    ap_muscle   <- 3.39      # Table S4
    ap_skin     <- 1.32      # Table S4
    ap_spleen   <- 3.18      # Table S4

    # =================================================================
    # 3. Unbound tissue-to-plasma partition coefficients (Kpu)
    #
    # Appendix S4 Suppl Eq 1-5, the Rodgers and Rowland equations for
    # moderate-to-strong bases. Kpu is the sum of four terms -- ionised
    # drug in intracellular water (Suppl Eq 1), drug in extracellular
    # water, binding to neutral lipids and neutral phospholipids (Suppl
    # Eq 2), and electrostatic binding to acidic phospholipids (Suppl Eq
    # 3) -- plus lysosomal sequestration (Suppl Eq 4).
    #
    # Suppl Eq 2/5 print "LogP" where the octanol:water partition
    # coefficient P itself is meant; the Bartelink 2022 supplement
    # (section III, note 5) states explicitly that "the octanol/water
    # partition coefficient (P) is included for binding affinity to
    # neutral lipids and phospholipids". P = 10^LogD accordingly, and
    # adipose uses the vegetable-oil:water coefficient 10^LogPV.
    # =================================================================
    pp   <- 10^logd                                        # octanol:water partition coefficient
    ppv  <- 10^logpv                                       # vegetable oil:water, adipose only
    dp   <- 1 + 10^(pka - ph_plasma)                       # denominator, unprotonated fraction in plasma
    diw  <- 1 + 10^(pka - ph_iw)
    xiw  <- diw / dp                                       # Suppl Eq 1: IW-to-plasma ionisation ratio
    aiw  <- 10^(pka - ph_iw)

    # Lysosomal pH terms. Liver and lung are resolved per cell type
    # (Table S5); every other tissue uses the whole-tissue lysosomal pH
    # of 5.3 (Table S2).
    l530 <- (1 + 10^(pka - 5.30)) / diw                    # Table S2: pHlys 5.3
    a530 <- 10^(pka - 5.30)
    l510 <- (1 + 10^(pka - 5.10)) / diw                    # Table S5: type II, residual, fat-storing, endothelial cells
    a510 <- 10^(pka - 5.10)
    l475 <- (1 + 10^(pka - 4.75)) / diw                    # Table S5: alveolar macrophages
    a475 <- 10^(pka - 4.75)
    l472 <- (1 + 10^(pka - 4.72)) / diw                    # Table S5: hepatocytes
    a472 <- 10^(pka - 4.72)
    l494 <- (1 + 10^(pka - 4.94)) / diw                    # Table S5: Kupffer cells
    a494 <- 10^(pka - 4.94)

    # Suppl Eq 2: binding to neutral lipids and neutral phospholipids
    lip_adipose  <- (ppv * fnl_adipose  + (0.3 * ppv + 0.7) * fnp_adipose)  / dp
    lip_brain    <- (pp  * fnl_brain    + (0.3 * pp  + 0.7) * fnp_brain)    / dp
    lip_gut      <- (pp  * fnl_gut      + (0.3 * pp  + 0.7) * fnp_gut)      / dp
    lip_heart    <- (pp  * fnl_heart    + (0.3 * pp  + 0.7) * fnp_heart)    / dp
    lip_kidney   <- (pp  * fnl_kidney   + (0.3 * pp  + 0.7) * fnp_kidney)   / dp
    lip_liver    <- (pp  * fnl_liver    + (0.3 * pp  + 0.7) * fnp_liver)    / dp
    lip_lung     <- (pp  * fnl_lung     + (0.3 * pp  + 0.7) * fnp_lung)     / dp
    lip_pancreas <- (pp  * fnl_pancreas + (0.3 * pp  + 0.7) * fnp_pancreas) / dp
    lip_muscle   <- (pp  * fnl_muscle   + (0.3 * pp  + 0.7) * fnp_muscle)   / dp
    lip_skin     <- (pp  * fnl_skin     + (0.3 * pp  + 0.7) * fnp_skin)     / dp
    lip_spleen   <- (pp  * fnl_spleen   + (0.3 * pp  + 0.7) * fnp_spleen)   / dp

    # Suppl Eq 1 + few + Eq 2 + Eq 3, then + Suppl Eq 4 scaled by
    # flys * fcell_type * xiw. Tissues with a single (whole-tissue) cell
    # type take fcell_type = 1.
    kpu_adipose <- xiw * fiw_adipose + few_adipose + lip_adipose +
      kabc * ap_adipose * aiw / dp +
      xiw * (l530 * fiw_adipose + kabc * ap_adipose * a530 / diw + lip_adipose) *
      flys_adipose * 1
    kpu_brain <- xiw * fiw_brain + few_brain + lip_brain +
      kabc * ap_brain * aiw / dp +
      xiw * (l530 * fiw_brain + kabc * ap_brain * a530 / diw + lip_brain) *
      flys_brain * 1
    kpu_gut <- xiw * fiw_gut + few_gut + lip_gut +
      kabc * ap_gut * aiw / dp +
      xiw * (l530 * fiw_gut + kabc * ap_gut * a530 / diw + lip_gut) *
      flys_gut * 1
    kpu_heart <- xiw * fiw_heart + few_heart + lip_heart +
      kabc * ap_heart * aiw / dp +
      xiw * (l530 * fiw_heart + kabc * ap_heart * a530 / diw + lip_heart) *
      flys_heart * 1
    kpu_kidney <- xiw * fiw_kidney + few_kidney + lip_kidney +
      kabc * ap_kidney * aiw / dp +
      xiw * (l530 * fiw_kidney + kabc * ap_kidney * a530 / diw + lip_kidney) *
      flys_kidney * 1
    kpu_pancreas <- xiw * fiw_pancreas + few_pancreas + lip_pancreas +
      kabc * ap_pancreas * aiw / dp +
      xiw * (l530 * fiw_pancreas + kabc * ap_pancreas * a530 / diw + lip_pancreas) *
      flys_pancreas * 1
    kpu_muscle <- xiw * fiw_muscle + few_muscle + lip_muscle +
      kabc * ap_muscle * aiw / dp +
      xiw * (l530 * fiw_muscle + kabc * ap_muscle * a530 / diw + lip_muscle) *
      flys_muscle * 1
    kpu_skin <- xiw * fiw_skin + few_skin + lip_skin +
      kabc * ap_skin * aiw / dp +
      xiw * (l530 * fiw_skin + kabc * ap_skin * a530 / diw + lip_skin) *
      flys_skin * 1
    kpu_spleen <- xiw * fiw_spleen + few_spleen + lip_spleen +
      kabc * ap_spleen * aiw / dp +
      xiw * (l530 * fiw_spleen + kabc * ap_spleen * a530 / diw + lip_spleen) *
      flys_spleen * 1

    # Liver: four cell types (Table S5) -- hepatocytes (flys 0.008,
    # fcell 0.78, pH 4.72), Kupffer (0.14, 0.02, 4.94), fat-storing
    # (0.002, 0.01, 5.1) and endothelial (0.07, 0.03, 5.1).
    kpu_liver <- xiw * fiw_liver + few_liver + lip_liver +
      kabc * ap_liver * aiw / dp +
      xiw * (l472 * fiw_liver + kabc * ap_liver * a472 / diw + lip_liver) * 0.008 * 0.78 +
      xiw * (l494 * fiw_liver + kabc * ap_liver * a494 / diw + lip_liver) * 0.14  * 0.02 +
      xiw * (l510 * fiw_liver + kabc * ap_liver * a510 / diw + lip_liver) * 0.002 * 0.01 +
      xiw * (l510 * fiw_liver + kabc * ap_liver * a510 / diw + lip_liver) * 0.07  * 0.03

    # Lung: three cell types (Table S5) -- alveolar macrophages (flys
    # 0.08, fcell 0.04, pH 4.75), type II cells (0.03, 0.08, 5.1) and
    # residual cells (0.01, 0.88, 5.1).
    kpu_lung <- xiw * fiw_lung + few_lung + lip_lung +
      kabc * ap_lung * aiw / dp +
      xiw * (l475 * fiw_lung + kabc * ap_lung * a475 / diw + lip_lung) * 0.08 * 0.04 +
      xiw * (l510 * fiw_lung + kabc * ap_lung * a510 / diw + lip_lung) * 0.03 * 0.08 +
      xiw * (l510 * fiw_lung + kabc * ap_lung * a510 / diw + lip_lung) * 0.01 * 0.88

    # Tumor: same fractional composition as lung (Table S3) but with two
    # NSCLC hallmarks. (i) The extracellular water is acidified to pH 6.7
    # (Table 1), which raises the EW term by (1 + 10^(pKa - pHewt)) / dp;
    # the Bartelink 2022 supplement section VC derives this and shows the
    # intracellular terms are unchanged because the two pH factors
    # cancel. (ii) Immune deprivation: the tumor is simulated as 100%
    # residual cells, i.e. the alveolar-macrophage and type II cell
    # fractions are removed (Bartelink 2022 supplement note 9).
    kpu_tumor <- xiw * fiw_lung +
      (1 + 10^(pka - ph_ew_tumor)) / dp * few_lung + lip_lung +
      kabc * ap_lung * aiw / dp +
      xiw * (l510 * fiw_lung + kabc * ap_lung * a510 / diw + lip_lung) * 0.01 * 1.0

    # The lumped rest-of-body compartment is pure extracellular water
    # (Table S1 footnote d), so few = 1 and every other term is zero.
    kpu_other <- 1

    # Suppl Eq 6: Kp = Kpu * Fu. The ODEs divide by Kp/BP throughout, so
    # the composite tissue-to-BLOOD partition coefficient is formed once.
    kb_adipose  <- kpu_adipose  * fu / bp
    kb_brain    <- kpu_brain    * fu / bp
    kb_gut      <- kpu_gut      * fu / bp
    kb_heart    <- kpu_heart    * fu / bp
    kb_kidney   <- kpu_kidney   * fu / bp
    kb_liver    <- kpu_liver    * fu / bp
    kb_lung     <- kpu_lung     * fu / bp
    kb_tumor    <- kpu_tumor    * fu / bp
    kb_pancreas <- kpu_pancreas * fu / bp
    kb_muscle   <- kpu_muscle   * fu / bp
    kb_skin     <- kpu_skin     * fu / bp
    kb_spleen   <- kpu_spleen   * fu / bp
    kb_other    <- kpu_other    * fu / bp

    # =================================================================
    # 4. Emergent blood concentrations
    #
    # cven and cart are blood concentrations in the two blood pools;
    # out_<tissue> is the blood concentration leaving that tissue,
    # C_tissue / (Kp_tissue / BP), as written throughout Suppl Eq 7-20.
    # =================================================================
    cven <- venous   / v_venous
    cart <- arterial / v_arterial
    out_adipose  <- adipose  / v_adipose  / kb_adipose
    out_brain    <- brain    / v_brain    / kb_brain
    out_gut      <- gut      / v_gut      / kb_gut
    out_heart    <- heart    / v_heart    / kb_heart
    out_kidney   <- kidney   / v_kidney   / kb_kidney
    out_liver    <- liver    / v_liver    / kb_liver
    out_lung     <- lung     / v_lung     / kb_lung
    out_tumor    <- tumor    / v_tumor    / kb_tumor
    out_pancreas <- pancreas / v_pancreas / kb_pancreas
    out_muscle   <- muscle   / v_muscle   / kb_muscle
    out_skin     <- skin     / v_skin     / kb_skin
    out_spleen   <- spleen   / v_spleen   / kb_spleen
    out_other    <- other    / v_other    / kb_other

    # =================================================================
    # 5. Saturable EGFR binding (main-text Equations 1-2, Suppl Eq 11-12)
    #
    # tot_<tissue> is the total EGFR amount in the tissue, [EGFR] * V,
    # in nmol. Because binding is 1:1 and covalent, the amount of EGFR
    # already occupied equals the amount of drug in the complex state, so
    # the free-EGFR fraction is 1 - complex / tot (Equation 2). The
    # association term is mass action between the free-EGFR amount and
    # the drug concentration in intracellular water, which the paper
    # writes as A_tissue / V_iw:
    #
    #   kon * A_tissue * ([EGFR] * V / V_iw) * fEGFRunbound
    #
    # As occupancy approaches 1 the free fraction goes to 0 and binding
    # stops, so the expression is self-limiting.
    # =================================================================
    tot_adipose  <- egfr_adipose  * v_adipose
    tot_gut      <- egfr_gut      * v_gut
    tot_heart    <- egfr_heart    * v_heart
    tot_kidney   <- egfr_kidney   * v_kidney
    tot_liver    <- egfr_liver    * v_liver
    tot_lung     <- egfr_lung     * v_lung
    tot_tumor    <- egfr_tumor    * v_tumor
    tot_pancreas <- egfr_pancreas * v_pancreas
    tot_muscle   <- egfr_muscle   * v_muscle
    tot_skin     <- egfr_skin     * v_skin

    bind_adipose  <- kon   * adipose  * (tot_adipose  / viw_adipose)  * (1 - adipose_complex  / tot_adipose)
    bind_gut      <- kon   * gut      * (tot_gut      / viw_gut)      * (1 - gut_complex      / tot_gut)
    bind_heart    <- kon   * heart    * (tot_heart    / viw_heart)    * (1 - heart_complex    / tot_heart)
    bind_kidney   <- kon   * kidney   * (tot_kidney   / viw_kidney)   * (1 - kidney_complex   / tot_kidney)
    bind_liver    <- kon   * liver    * (tot_liver    / viw_liver)    * (1 - liver_complex    / tot_liver)
    bind_lung     <- kon   * lung     * (tot_lung     / viw_lung)     * (1 - lung_complex     / tot_lung)
    bind_tumor    <- kon_t * tumor    * (tot_tumor    / viw_tumor)    * (1 - tumor_complex    / tot_tumor)
    bind_pancreas <- kon   * pancreas * (tot_pancreas / viw_pancreas) * (1 - pancreas_complex / tot_pancreas)
    bind_muscle   <- kon   * muscle   * (tot_muscle   / viw_muscle)   * (1 - muscle_complex   / tot_muscle)
    bind_skin     <- kon   * skin     * (tot_skin     / viw_skin)     * (1 - skin_complex     / tot_skin)

    # =================================================================
    # 6. Differential equations (Appendix S5 Suppl Eq 7-20)
    # =================================================================
    # Absorption (Suppl Eq 13-14). The published Suppl Eq 13 is written
    # as -Ka * Dose * F, i.e. against the administered dose rather than
    # the amount remaining, which would empty the lumen linearly and then
    # drive it negative. The intended and standard reading is first-order
    # loss from the lumen with the absorbed fraction F reaching the gut;
    # F is applied as a bioavailability on the dosing compartment so that
    # dosing gut_lumen with a plain amount behaves correctly. See the
    # vignette Errata.
    d/dt(gut_lumen) <- -ka * gut_lumen
    f(gut_lumen)    <- exp(lfdepot)

    # Venous blood (Suppl Eq 9). The lung and the lung tumor are both
    # perfused from venous blood, in parallel; every other tissue that
    # does not drain through the portal vein returns here.
    d/dt(venous) <- -(q_lung + q_tumor) * cven +
      q_adipose * out_adipose + q_brain * out_brain + q_heart * out_heart +
      q_kidney * out_kidney + q_liver * out_liver + q_muscle * out_muscle +
      q_other * out_other + q_pancreas * out_pancreas + q_skin * out_skin

    # Arterial blood (Suppl Eq 10). Fed by the lung and tumor outflow and
    # drained by every arterially perfused tissue, including the gut, the
    # spleen and the hepatic artery.
    d/dt(arterial) <- q_lung * out_lung + q_tumor * out_tumor -
      (q_adipose + q_brain + q_gut + q_heart + q_ha + q_kidney + q_muscle +
         q_other + q_pancreas + q_skin + q_spleen) * cart

    # Non-eliminating tissues without EGFR (Suppl Eq 7)
    d/dt(brain)  <- q_brain  * (cart - out_brain)
    d/dt(spleen) <- q_spleen * (cart - out_spleen)
    d/dt(other)  <- q_other  * (cart - out_other)

    # Non-eliminating tissues with EGFR (Suppl Eq 11-12)
    d/dt(adipose)          <- q_adipose * (cart - out_adipose) - bind_adipose + koff * adipose_complex
    d/dt(adipose_complex)  <- bind_adipose - koff * adipose_complex
    d/dt(heart)            <- q_heart * (cart - out_heart) - bind_heart + koff * heart_complex
    d/dt(heart_complex)    <- bind_heart - koff * heart_complex
    d/dt(muscle)           <- q_muscle * (cart - out_muscle) - bind_muscle + koff * muscle_complex
    d/dt(muscle_complex)   <- bind_muscle - koff * muscle_complex
    d/dt(skin)             <- q_skin * (cart - out_skin) - bind_skin + koff * skin_complex
    d/dt(skin_complex)     <- bind_skin - koff * skin_complex
    d/dt(pancreas)         <- q_pancreas * (cart - out_pancreas) - bind_pancreas + koff * pancreas_complex
    d/dt(pancreas_complex) <- bind_pancreas - koff * pancreas_complex

    # Gut (Suppl Eq 14): arterial inflow, absorption from the lumen and
    # EGFR binding; the outflow is collected by the liver.
    d/dt(gut)         <- ka * gut_lumen + q_gut * (cart - out_gut) - bind_gut + koff * gut_complex
    d/dt(gut_complex) <- bind_gut - koff * gut_complex

    # Liver (Suppl Eq 17-18): portal inflow from gut and spleen plus the
    # hepatic artery, minus hepatic venous outflow, minus well-stirred
    # intrinsic clearance acting on the unbound fraction.
    d/dt(liver) <- q_gut * out_gut + q_spleen * out_spleen + q_ha * cart -
      q_liver * out_liver - clint * out_liver * fu -
      bind_liver + koff * liver_complex
    d/dt(liver_complex) <- bind_liver - koff * liver_complex

    # Kidney (Suppl Eq 19-20): renal clearance of unbound drug.
    d/dt(kidney) <- q_kidney * (cart - out_kidney) - clr * fu * out_kidney -
      bind_kidney + koff * kidney_complex
    d/dt(kidney_complex) <- bind_kidney - koff * kidney_complex

    # Lung (Suppl Eq 8) and lung tumor: both perfused from venous blood
    # and draining into arterial blood (Suppl Eq 9-10). Suppl Eq 8 is the
    # base-model lung equation; lung and tumor are both EGFR-expressing
    # (main text section 2.2.2), so the Suppl Eq 11-12 binding terms
    # apply to them as well. The tumor uses the 50/50 wild-type /
    # L858R-T790M on and off rates.
    d/dt(lung)          <- q_lung * (cven - out_lung) - bind_lung + koff * lung_complex
    d/dt(lung_complex)  <- bind_lung - koff * lung_complex
    d/dt(tumor)         <- q_tumor * (cven - out_tumor) - bind_tumor + koff_t * tumor_complex
    d/dt(tumor_complex) <- bind_tumor - koff_t * tumor_complex

    # =================================================================
    # 7. Observations
    #
    # Cc is the venous PLASMA concentration, which is what the clinical
    # concentration-time profiles of section 2.3.2 are measured in; the
    # venous state holds BLOOD, so it is divided by the blood-to-plasma
    # ratio. Tissue concentrations are total (free + EGFR-bound) drug per
    # unit total tissue volume, matching what PET measures. occ_<tissue>
    # is fractional EGFR occupancy, the quantity plotted in Figure 4C.
    # =================================================================
    Cc          <- cven / bp
    Cc_arterial <- cart / bp
    Cblood      <- cven

    Ct_adipose  <- (adipose  + adipose_complex)  / v_adipose
    Ct_brain    <- brain / v_brain
    Ct_gut      <- (gut      + gut_complex)      / v_gut
    Ct_heart    <- (heart    + heart_complex)    / v_heart
    Ct_kidney   <- (kidney   + kidney_complex)   / v_kidney
    Ct_liver    <- (liver    + liver_complex)    / v_liver
    Ct_lung     <- (lung     + lung_complex)     / v_lung
    Ct_tumor    <- (tumor    + tumor_complex)    / v_tumor
    Ct_pancreas <- (pancreas + pancreas_complex) / v_pancreas
    Ct_muscle   <- (muscle   + muscle_complex)   / v_muscle
    Ct_skin     <- (skin     + skin_complex)     / v_skin
    Ct_spleen   <- spleen / v_spleen

    occ_adipose  <- adipose_complex  / tot_adipose
    occ_gut      <- gut_complex      / tot_gut
    occ_heart    <- heart_complex    / tot_heart
    occ_kidney   <- kidney_complex   / tot_kidney
    occ_liver    <- liver_complex    / tot_liver
    occ_lung     <- lung_complex     / tot_lung
    occ_tumor    <- tumor_complex    / tot_tumor
    occ_pancreas <- pancreas_complex / tot_pancreas
    occ_muscle   <- muscle_complex   / tot_muscle
    occ_skin     <- skin_complex     / tot_skin

    Cc ~ prop(propSd)
  })
}
