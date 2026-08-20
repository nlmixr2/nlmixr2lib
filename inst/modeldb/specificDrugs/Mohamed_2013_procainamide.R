Mohamed_2013_procainamide <- function() {
  description <- paste(
    "Joint parent-metabolite two-compartment popPK model for procainamide",
    "and its major active N-acetyl metabolite N-acetylprocainamide (NAPA)",
    "in a single 40-year-old male with chronic kidney disease stage 5",
    "receiving continuous renal replacement therapy (CRRT) after aortic",
    "valve replacement (Mohamed 2013). Structural model (paper equations",
    "1-3): procainamide central and peripheral compartments with linear",
    "distribution clearance (Cld) and two parallel first-order elimination",
    "arms from the central compartment - Cl_other combining all non-",
    "metabolism procainamide clearance pathways (CRRT, residual renal, and",
    "any additional non-NAPA-forming metabolic routes) and Cl_f,napa the",
    "formation clearance of NAPA from procainamide; the metabolite NAPA is",
    "described as a single central compartment fed by the Cl_f,napa flux",
    "out of the procainamide central compartment and eliminated linearly",
    "with total clearance Cl_napa (CRRT plus residual renal). The NAPA",
    "apparent volume of distribution (Vc)_N was FIXED at 100 L, derived as",
    "the patient's body weight (70 kg) times a literature Vd of 1.5 L/kg",
    "for NAPA in functionally anephric patients on hemodialysis (paper ref",
    "9), because the fraction of procainamide elimination going through",
    "NAPA formation could not be identified independently from the",
    "concentration data. All other structural parameters were successfully",
    "estimated by first-order conditional estimation in NONMEM VII. RUV",
    "was described by a proportional error model (paper equation 4) on",
    "both procainamide and NAPA concentrations; the estimated residual-",
    "error magnitude is not reported in Mohamed 2013 Table 1, so this",
    "library implementation encodes a nominal 15 percent CV proportional",
    "SD on each output and documents the omission in the vignette Errata.",
    "Because the source is a single-subject case report, no between-",
    "subject variability is identifiable and no eta parameters are",
    "encoded; the model is a typical-value description tied to this",
    "specific patient's disposition on CRRT."
  )
  reference <- paste(
    "Mohamed AN, Abdelhady AM, Spencer D, Sowinski KM, Tisdale JE,",
    "Overholser BR. Pharmacokinetic modeling and simulation of",
    "procainamide and N-acetylprocainamide during continuous renal",
    "replacement therapy. Am J Kidney Dis. 2013 Jun;61(6):1046-1048.",
    "doi:10.1053/j.ajkd.2013.02.358."
  )
  vignette <- "Mohamed_2013_procainamide"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/L"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    central      = list(analyte = "procainamide", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1  = list(analyte = "procainamide", units = "mg", specimen = "plasma", verified = FALSE),
    central_napa = list(analyte = "N-acetylprocainamide (NAPA)", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list()

  population <- list(
    species        = "human",
    n_subjects     = 1L,
    n_studies      = 1L,
    n_observations = "Serial plasma samples for procainamide and NAPA collected over 120 hours as part of routine clinical care (Case Report section, paper Figure 1 caption).",
    age_range      = "40 years (single patient)",
    weight_range   = "70 kg (single patient)",
    sex_female_pct = 0,
    race_ethnicity = "White (single 40-year-old white male; Case Report first paragraph).",
    disease_state  = "CKD stage 5 (serum creatinine 10 mg/dL on ICU admission) with residual renal function (~1400 mL urine per 24 h) on continuous renal replacement therapy (CRRT) at 2000 mL/h via a PRISMA system with an AN69/M100 hemofilter set and 150 mL/min blood flow, following aortic valve replacement complicated by monomorphic ventricular tachycardia treated with procainamide loading and maintenance infusions.",
    dose_range     = "IV procainamide: 1000 mg loading dose infused at 50 mg/min for 20 min; maintenance infusions of 2 mg/min for 6 h, 4 mg/min for 8 h, 2 mg/min for 7 h, and 1 mg/min for ~13 h (approximately 36 h total therapy). Simulated regimens in Figure 2 use a 20 mg/min loading for 30 min followed by 48 h of 1, 2, or 4 mg/min maintenance.",
    regions        = "United States (Indiana University Health Methodist Hospital, Indianapolis, IN).",
    notes          = "Single-subject case report; the popPK model is a typical-value description of this specific patient's disposition. Between-subject variability is not identifiable from a single subject and is therefore not encoded (paper Table 1 does not report omega values). CRRT was maintained continuously at 2000 mL/h during the entire procainamide administration period; the estimated Cl_other therefore combines residual renal (~58 mL/h from the ~1400 mL/24 h urine output) with CRRT hemofiltration removal of procainamide and cannot be decomposed further without procainamide dialysate concentrations."
  )

  ini({
    # All structural parameter values are from Mohamed 2013 Table 1
    # (final compartmental PK model estimated from the patient's serial
    # procainamide and NAPA plasma concentrations by first-order
    # conditional estimation in NONMEM VII). RSE values shown in the
    # comments are the paper's reported relative standard errors on the
    # estimate; none of the estimated parameters carry between-subject
    # variability because the dataset is a single subject.

    # Procainamide non-formation clearance (paper Cl_other) - combines
    # CRRT hemofiltration, residual renal excretion, and any procainamide
    # elimination routes other than N-acetylation to NAPA (paper "Cl_other
    # combines all other procainamide clearance pathways (CRRT, renal,
    # metabolic)"). Encoded as `lcl_nonmet` per the parent-metabolite
    # canonical (parent's non-formation elimination arm; sibling of
    # `lcl_met`).
    lcl_nonmet <- log(3.54)          ; label("Procainamide non-formation clearance Cl_other (L/h): CRRT + residual renal + non-NAPA-forming metabolic")           # Table 1: Cl_other = 3.54 L/h, RSE 44.3%

    # Procainamide-to-NAPA formation clearance (paper Cl_f,napa) - the
    # metabolic-formation flux from parent to metabolite via hepatic
    # N-acetylation. Encoded as `lcl_met` per the parent-metabolite
    # canonical (metabolic-formation arm; sibling of `lcl_nonmet`).
    lcl_met    <- log(3.70)          ; label("Procainamide-to-NAPA formation clearance Cl_f,napa (L/h)")                                                          # Table 1: Cl_f,napa = 3.70 L/h, RSE 26.3%

    # NAPA total elimination clearance (paper Cl_napa) - combines CRRT
    # hemofiltration and residual renal excretion of NAPA. Encoded as
    # `lcl_napa` per the metabolite canonical (`lcl_<metabsuffix>`).
    lcl_napa   <- log(2.96)          ; label("NAPA total elimination clearance Cl_napa (L/h): CRRT + residual renal")                                             # Table 1: Cl_napa = 2.96 L/h, RSE 32.1%

    # Procainamide inter-compartmental (distribution) clearance between
    # its central and peripheral compartments (paper Cld).
    lq         <- log(19.1)          ; label("Procainamide distribution clearance Cld (L/h)")                                                                     # Table 1: Cld = 19.1 L/h, RSE 24.6%

    # Procainamide central volume of distribution (paper (Vc)_P). The
    # apparent volume is relative to the total dose because parent is
    # administered IV (F = 1); no separate F term is estimated. The RSE
    # is high (82.5%) because the central volume is largely informed by
    # the earliest post-loading samples of a single subject.
    lvc        <- log(30.7)          ; label("Procainamide central volume (Vc)_P (L)")                                                                            # Table 1: (Vc)_P = 30.7 L, RSE 82.5%

    # Procainamide peripheral volume of distribution (paper (Vp)_P).
    lvp        <- log(169)           ; label("Procainamide peripheral volume (Vp)_P (L)")                                                                         # Table 1: (Vp)_P = 169 L, RSE 26.9%

    # NAPA central volume of distribution (paper (Vc)_N). FIXED at 100 L,
    # computed as the patient's body weight (70 kg) times a literature
    # NAPA Vd of ~1.5 L/kg in functionally anephric patients on
    # hemodialysis (paper reference 9). Fixing was required because the
    # fraction of procainamide elimination going through NAPA formation
    # could not be identified independently from the concentration data.
    lvc_napa   <- fixed(log(100))    ; label("NAPA central volume (Vc)_N (L), 100 L (= 70 kg x 1.5 L/kg from paper reference 9)")                        # Table 1: (Vc)_N = 100 L FIXED

    # Residual error - Mohamed 2013 equation 4 defines a proportional
    # error model (yi = yhat_i * (1 + eps_i)) on both procainamide and
    # NAPA plasma concentrations. Table 1 does not report the estimated
    # value of the proportional-error variance; a nominal 15% CV is used
    # here as a plausible value for downstream simulations, encoded as
    # a typical-value fallback per the operator's standing "unreported
    # RUV -> typical-value" policy. See the vignette Assumptions and
    # deviations section.
    propSd       <- 0.15             ; label("Procainamide proportional residual SD (fraction); nominal 15% CV - value not reported in Table 1")                  # not reported in Table 1; nominal 15% CV
    propSd_napa  <- 0.15             ; label("NAPA proportional residual SD (fraction); nominal 15% CV - value not reported in Table 1")                          # not reported in Table 1; nominal 15% CV
  })

  model({
    # -----------------------------------------------------------------
    # Individual structural parameters. No IIV is encoded because the
    # source dataset is a single subject and between-subject variability
    # is not identifiable; typical-value replication uses `rxode2::zeroRe`
    # in the vignette in any case, but keeping the individual-parameter
    # form here would just multiply typical values by exp(0) throughout.
    # -----------------------------------------------------------------
    cl_nonmet <- exp(lcl_nonmet)
    cl_met    <- exp(lcl_met)
    cl_napa   <- exp(lcl_napa)
    q         <- exp(lq)
    vc        <- exp(lvc)
    vp        <- exp(lvp)
    vc_napa   <- exp(lvc_napa)

    # -----------------------------------------------------------------
    # Two-compartment procainamide + one-compartment NAPA metabolite
    # ODE system (Mohamed 2013 equations 1-3, transcribed to
    # nlmixr2 rate-form). The parent central compartment loses mass to
    # non-formation elimination (Cl_other), to metabolite formation
    # (Cl_f,napa), and to the peripheral compartment (Cld). The parent
    # peripheral compartment exchanges only with the central via Cld.
    # The metabolite central compartment is fed by the parent formation
    # flux (Cl_f,napa / (Vc)_P * X1) and eliminated with total
    # clearance Cl_napa.
    #
    # State units: `central` and `peripheral1` in mg (procainamide);
    # `central_napa` in mg (NAPA). Concentrations `Cc` and `Cc_napa`
    # in mg/L (matching the paper's Figure 1 y-axis and the
    # therapeutic-range thresholds of 4-12 mg/L procainamide and
    # 10-30 mg/L combined procainamide + NAPA).
    # -----------------------------------------------------------------
    d/dt(central)      <- -(cl_nonmet + cl_met + q) / vc * central + q / vp * peripheral1
    d/dt(peripheral1)  <-  q / vc * central - q / vp * peripheral1
    d/dt(central_napa) <-  cl_met / vc * central - cl_napa / vc_napa * central_napa

    # -----------------------------------------------------------------
    # Plasma concentration outputs (mg/L). The parent uses the bare
    # canonical `Cc`; the metabolite uses the `_napa` suffix per the
    # metabolite-canonical naming rule.
    # -----------------------------------------------------------------
    Cc      <- central      / vc
    Cc_napa <- central_napa / vc_napa

    Cc      ~ prop(propSd)
    Cc_napa ~ prop(propSd_napa)
  })
}
