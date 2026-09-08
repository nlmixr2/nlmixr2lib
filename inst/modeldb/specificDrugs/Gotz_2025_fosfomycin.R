Gotz_2025_fosfomycin <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous fosfomycin in critically",
    "ill adults with and without kidney replacement therapy (KRT), pooled from four",
    "prospective observational studies (45 patients, 727 concentrations). Total",
    "clearance is the sum of two parallel arms. The body clearance arm is a power",
    "function of MDRD-estimated glomerular filtration rate normalised to the",
    "48.4 mL/min/1.73 m^2 cohort median, switched off entirely in anuric patients",
    "(24-hour urine output < 100 mL) and in the study B cohort, which showed no",
    "fosfomycin elimination between KRT sessions. The dialysis clearance arm is a",
    "power function of the dialysate flow rate normalised to 42 mL/min and is gated",
    "to periods when KRT is actually running; the same term covers continuous KRT",
    "(Q_D = 42 mL/min) and prolonged-intermittent KRT (Q_D = 250 mL/min), which",
    "differ only in Q_D and in the duration of the KRT period. The peripheral volume",
    "of distribution expands linearly with time since the first dose, by 0.07% per",
    "minute, and this expansion is restricted to anuric patients. Interindividual",
    "variability is diagonal on body clearance and on both volumes; residual error is",
    "combined proportional plus additive (Gotz 2025).",
    sep = " "
  )
  reference <- paste(
    "Gotz KM, Kreuer S, Volz AK, Parker SL, Roberts JA, Dimopoulos G, Dimski T,",
    "Kindgen-Milles D, Beuche LKV, Kielstein JT, Lehr T. Population pharmacokinetics",
    "of intravenous fosfomycin: dose optimization for critically ill patients with and",
    "without kidney replacement therapy. Antimicrob Agents Chemother.",
    "2025;69(6):e01779-24. doi:10.1128/aac.01779-24",
    sep = " "
  )
  vignette <- "Gotz_2025_fosfomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE: checked against Gotz 2025, which
  # administers fosfomycin as 4-8 g (mg-scale) i.v. doses and quantifies
  # "fosfomycin serum concentrations" (study A) or "fosfomycin plasma
  # concentrations" (studies B-D) by HPLC/MS or LC-MS/MS (Methods, Patients and
  # study design). Both states hold amounts of unchanged fosfomycin; the drug is
  # not metabolized (Introduction).
  compartmentData <- list(
    central     = list(analyte = "fosfomycin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "fosfomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Estimated glomerular filtration rate calculated with the MDRD equation",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "BSA-NORMALIZED, creatinine-based eGFR from the four-variable MDRD equation, which",
        "Table S1 of the supplement gives as 175 * SCr^-1.154 * age^-0.203 (* 0.742 if female)",
        "(* 1.212 if black), with SCr in mg/dL and age in years. Time-varying: recalculated",
        "from repeated serum creatinine measurements, a median (range) of 4 (0-10) per patient",
        "(Methods, Data analysis). The 48.4 mL/min/1.73 m^2 normalisation constant is the",
        "entire-cohort median in Table 1, not a rounded standard.",
        "Gotz 2025 deliberately chose eGFR_MDRD over the alternatives it also screened:",
        "Cockcroft-Gault eCrCL was rejected because the KRT patients were considerably heavier",
        "than the non-KRT patients so a weight-containing equation would overestimate their",
        "kidney function, CKD-EPI was rejected because its accuracy is best above",
        "60 mL/min/1.73 m^2, and measured 24-hour urinary CrCL was available in only two of",
        "the four pooled studies (Discussion). Absolute (non-BSA-normalized) eGFR_MDRD in",
        "mL/min was also tested and did NOT improve the fit (Results), so the model is",
        "calibrated against the BSA-NORMALIZED value and the canonical mL/min/1.73 m^2 units",
        "are the correct ones here. This is a DIFFERENT kidney-function marker from the one",
        "in the sibling model Huppe_2023_fosfomycin.R, which uses a raw, non-normalized",
        "MEASURED urinary creatinine clearance; the two are not interchangeable."
      ),
      source_name        = "eGFR_MDRD"
    ),
    URINE_VOL_24H = list(
      description        = "24-hour urine output, used as the anuria gate",
      units              = "mL/24h",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Gates TWO separate parts of the model, in OPPOSITE directions, and is the only",
        "covariate in this model that does so.",
        "(1) Body clearance is switched off in anuric patients: Methods, Model development,",
        "'CL_body was fixed to 0 in patients presenting 24-hour urine output < 100 mL',",
        "printed again in Table 2 footnote b as '(x 0 if 24 - h urine output < 100 mL)'.",
        "(2) The time-dependent expansion of the peripheral volume applies ONLY to anuric",
        "patients: Results, Population pharmacokinetic model, 'As this variation in the",
        "distribution of fosfomycin may primarily occur in patients experiencing fluid",
        "retention ... the increasing V_P was restricted to patients experiencing anuria.'",
        "Table 2 footnote c omits this second gate; the restriction is confirmed against the",
        "paper's own Figure S6 simulations - see the vignette source-trace section and Errata.",
        "Ten of 45 patients (22.2%) were anuric (Results, Patients). Entire-cohort median",
        "(IQR) 24-hour urine output 410 (0-1400) mL, ranging from 25 mL in study C to",
        "3400 mL in study D (Table 1). The 100 mL/24h cutoff is the paper's own definition of",
        "anuria (Table 1 footnote a) and coincides with the preserved-diuresis cutoff already",
        "used by the sibling model Huppe_2023_fosfomycin.R.",
        "The bare ANURIA canonical is deliberately NOT used: the source ascertains anuria from",
        "a urine VOLUME against a stated cutoff rather than tabulating a yes/no flag, and the",
        "sibling fosfomycin model encodes the identical cutoff from the identical column."
      ),
      source_name        = "24-hour urine output"
    ),
    RRT_CRRT_ACTIVE = list(
      description        = "Kidney-replacement-therapy-active indicator (1 while KRT is running, 0 otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no KRT running)",
      notes              = paste(
        "Time-varying WITHIN subject. Methods, Model development: 'CL_KRT was fixed to 0 for",
        "patients not receiving KRT or time periods between KRT sessions.' The pooled data",
        "set supports this because blood samples were available between KRT sessions",
        "(study B) and on non-dialysis days (study C) - the Discussion identifies this as",
        "what makes CL_body and CL_KRT separately identifiable.",
        "Covers BOTH modalities in this model. Gotz 2025 fits ONE dialysis-clearance",
        "equation to continuous KRT (CKRT, study C, CVVHD) and to prolonged-intermittent KRT",
        "(PIKRT, studies A-B, Genius batch dialysis); the modalities are distinguished only by",
        "their dialysate flow rate DFR and by how long the KRT period lasts, not by a separate",
        "parameter. The canonical's own definition spans 'a continuous or extended",
        "extracorporeal renal-replacement-therapy modality', explicitly including sustained",
        "low-efficiency dialysis and extended daily diafiltration, which is the class PIKRT",
        "belongs to; so the single RRT_CRRT_ACTIVE column is the right gate for both and",
        "RRT_HEMODIAL_ACTIVE is not additionally needed.",
        "Simulation settings (Methods, Simulations): CKRT for 48 h from the start of",
        "treatment, PIKRT for 8 h on the second day of treatment. 33 of 45 patients (73.3%)",
        "underwent KRT - 18 PIKRT and 15 CKRT."
      ),
      source_name        = "KRT"
    ),
    DFR = list(
      description        = "Dialysate flow rate through the extracorporeal circuit",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying within subject; the single driver of the dialysis clearance arm, which",
        "enters as a power function normalised to 42 mL/min (the entire-cohort median in",
        "Table 1). Entire-cohort median (IQR) 42 (33-50) mL/min, with a very wide between-",
        "modality spread: 190 (28-190) mL/min in study A and 250 (240-250) mL/min in study B",
        "(both PIKRT) against 33 (33-42) mL/min in study C (CKRT), P < 0.001 (Table 1).",
        "The paper's Q_D-vs-CL_KRT figure (Fig. 3C) is labelled 'Q_D (L/h)' but is in fact",
        "plotted in mL/min: CL_KRT = 2.0 * (Q_D/42)^0.587 evaluates to 1.64 L/h at Q_D = 30",
        "and 5.96 L/h at Q_D = 270, which reproduces that panel's endpoints exactly, whereas",
        "the same values read as L/h do not. The canonical mL/min is therefore correct and",
        "the axis label is a publication erratum - see the vignette Errata.",
        "Gotz 2025 deliberately did NOT use the Michaels equation that the sibling model",
        "Huppe_2023_fosfomycin.R uses, even though it evaluated it: the mass transfer-area",
        "coefficient K0A is dialyzer-specific and varies with Q_D for small molecules, and the",
        "pooled studies used different dialyzers, 'which would have limited the generalizability",
        "of our model findings' (Discussion). Blood flow rate Q_B was screened alongside Q_D",
        "and not retained, so BFR is not a covariate of this model.",
        "Meaningful only while RRT_CRRT_ACTIVE = 1; the arm is gated off otherwise."
      ),
      source_name        = "Q_D"
    ),
    STUDY_GERECKE = list(
      description        = "Study B (Gerecke 2021 PIKRT cohort) membership indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any of studies A, C, D, or a new patient being simulated)",
      notes              = paste(
        "A hard, non-estimated 0/1 gate that zeroes the body-clearance arm, exactly as",
        "printed in Table 2 footnote b: '(x 0 if study B)'. Results, Population",
        "pharmacokinetic model: 'CL_body was fixed to 0 for study B since these patients",
        "presented no fosfomycin elimination without KRT (Fig. 2), which led to individual",
        "CL_body estimates close to zero.'",
        "Study B is Gerecke LKV et al., J Antimicrob Chemother 2021;77:169-173,",
        "doi:10.1093/jac/dkab357 (reference 24 of Gotz 2025; reference 2 of the supplement):",
        "eight patients on prolonged-intermittent KRT with a Genius 90 batch dialysis system.",
        "It contributed 8 of the 45 patients (17.8%).",
        "This gate is NOT the same as, and is not implied by, the anuria gate: study B's",
        "median (IQR) 24-hour urine output was 700 (500-1000) mL, so none of its patients",
        "meets the < 100 mL anuria criterion, and both gates must be carried separately.",
        "Set to 0 when simulating a new patient - the paper's own Monte Carlo simulations",
        "(Methods, Simulations; Fig. S6) do not invoke it. Its purpose here is to reproduce",
        "the published fit faithfully, and it is the reason the paper also reports a reduced",
        "data set (n = 37) without study B, whose estimates Table 2 shows to be consistent."
      ),
      source_name        = "study B"
    )
  )

  covariatesDataExcluded <- list(
    BFR = list(
      description = "Blood flow rate through the extracorporeal circuit",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened on CL_KRT together with Q_D (Methods, Model development: 'we examined the effects of the dialyzate flow rate (Q_D) and Q_B on CL_KRT') and not retained; Q_D was 'the key variable for CL_KRT' (Discussion). Entire-cohort median (IQR) 100 (100-150) mL/min (Table 1). Retained by the sibling model Huppe_2023_fosfomycin.R, which uses the Michaels equation and therefore needs both flow rates."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened on all model parameters (Methods, Model development) and not retained. Entire-cohort median (IQR) 80 (70-90) kg (Table 1). Enters the model only indirectly, through the BSA used to compute eGFR_MDRD."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened on all model parameters and not retained. Entire-cohort median (IQR) 63 (57-75) years (Table 1). Enters the model only indirectly, as a term of the MDRD equation."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on all model parameters and not retained. 11 of 45 patients (24%) were female (Table 1). Enters the model only indirectly, as the 0.742 female factor of the MDRD equation."
    ),
    BSA = list(
      description = "Body surface area, Mosteller equation",
      units       = "m^2",
      type        = "continuous",
      notes       = "Screened on all model parameters and not retained. Entire-cohort median (IQR) 1.96 (1.83-2.05) m^2 (Table 1). Used only to convert relative eGFR_MDRD to absolute eGFR_MDRD, a substitution that did not significantly improve the model (Results)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened on all model parameters and not retained as a direct covariate. Entire-cohort median (IQR) 1.3 (0.94-2.1) mg/dL (Table 1). Enters the model only indirectly, as the principal term of the MDRD equation."
    ),
    POTASSIUM = list(
      description = "Serum potassium",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened on all model parameters and not retained. Entire-cohort median (IQR) 4.3 (4-4.5) mmol/L, missing for 44.4% of patients (studies B and D), which also prevented any evaluation of hypokalemia as an adverse effect (Methods, Data analysis)."
    ),
    SODIUM = list(
      description = "Serum sodium",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened on all model parameters and not retained. Entire-cohort median (IQR) 140 (140-150) mmol/L, missing for 44.4% of patients. Of clinical interest because intravenous fosfomycin is given as the disodium salt (Discussion)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 45L,
    n_studies      = 4L,
    n_observations = 727L,
    age_range      = "median (IQR) 63 (57-75) years",
    weight_range   = "median (IQR) 80 (70-90) kg",
    sex_female_pct = 24,
    race_ethnicity = "Not reported.",
    disease_state  = paste(
      "Critically ill adults with acute kidney injury or chronic kidney disease. 33 of 45",
      "(73.3%) underwent kidney replacement therapy: 18 prolonged-intermittent KRT",
      "(studies A and B) and 15 continuous KRT (study C). Twelve patients (study D) required",
      "no KRT. Ten patients (22.2%) were anuric, defined as 24-hour urine output < 100 mL."
    ),
    renal_function = paste(
      "Entire cohort median (IQR): eGFR_MDRD 48.4 (33.3-77.2) mL/min/1.73 m^2; eGFR_CKD-EPI",
      "51.3 (33.9-80.4) mL/min/1.73 m^2; Cockcroft-Gault eCrCL 56.3 (45.3-89.5) mL/min;",
      "measured 24-hour urinary CrCL 0 (0-31.1) mL/min; serum creatinine 1.3 (0.94-2.1) mg/dL;",
      "24-hour urine output 410 (0-1400) mL. eGFR_MDRD category counts at first dose",
      "(number anuric in parentheses): >= 90, 7 (1); 60-89, 8 (2); 45-59, 6 (3); 30-44, 12 (1);",
      "15-29, 4 (3); < 15, 0 (0)."
    ),
    rrt_settings   = paste(
      "Entire cohort median (IQR) blood flow rate 100 (100-150) mL/min and dialysate flow rate",
      "42 (33-50) mL/min, differing sharply by modality (P < 0.001): study A 190 (28-190)",
      "mL/min for both, study B 250 (240-250) mL/min for both, study C 100 (100-100) mL/min",
      "blood flow and 33 (33-42) mL/min dialysate flow. Study A used the Genius system with a",
      "Polyflux140H polyamix hemofilter and 8-hour PIKRT episodes; study B used the GENIUS 90",
      "batch dialysis system with Polyflux 17L, FX 60 or Polyflux 170H dialysers and 6-hour",
      "PIKRT episodes; study C used the multiFiltrate Ci-Ca with Ultraflux AV 1000S",
      "polysulfone hemofilters for CVVHD of physician-determined duration (Supplement,",
      "Kidney replacement therapy modalities)."
    ),
    dose_range     = paste(
      "5 g three times daily intravenously in patients with KRT; 4 g four times daily or 6 g",
      "three times daily in patients without KRT. Infusion duration 30-60 min in studies A, B",
      "and D and 120 min in study C. The Monte Carlo simulations additionally explored 4, 5",
      "and 8 g three times daily, 8 g twice daily and 4 g four times daily."
    ),
    regions        = "Germany (studies A, B, C), Greece (study D).",
    notes          = paste(
      "Pooled analysis of four prospective observational studies: A = Dimski 2021",
      "(doi:10.1038/s41598-021-91423-9, n = 10, PIKRT), B = Gerecke 2021",
      "(doi:10.1093/jac/dkab357, n = 8, PIKRT), C = Huppe 2023",
      "(doi:10.1038/s41598-023-45084-5, n = 15, CVVHD), D = Parker 2015",
      "(doi:10.1128/AAC.01321-15, n = 12, no KRT). Study C is the source of the sibling model",
      "Huppe_2023_fosfomycin.R already in this library, so that model is a proper subset of",
      "this one's data.",
      "Missing data were handled by median imputation for entirely missing patient",
      "characteristics and by last-observation-carried-forward for continuous laboratory data",
      "(Methods, Data analysis). Study B contributed no height, BMI, BSA or laboratory",
      "markers; serum potassium and sodium were missing for 44.4% of patients and measured",
      "24-hour urinary CrCL for 40.0%.",
      "Estimation was by FOCE with interaction in NONMEM 7.4. The final model was evaluated",
      "by a prediction-corrected VPC on 1000 replicates (Fig. S4, S5). All final estimates",
      "had RSE <= 32% (Table 2), and a sensitivity re-fit excluding study B (n = 37) gave",
      "consistent estimates."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural parameters. All values are the "Final model / Full data set
    # (n = 45)" column of Gotz 2025 Table 2; the parenthesised number in each
    # trace is that table's RSE. Nothing here is fixed - every parameter was
    # estimated, and the paper reports RSEs <= 32% throughout (Results).
    # ---------------------------------------------------------------------

    # Body clearance intercept, i.e. the clearance of a non-anuric patient
    # outside study B whose eGFR_MDRD sits at the 48.4 mL/min/1.73 m^2 cohort
    # median. "The typical CL_body was estimated at 1.6 L/h" (Results).
    lcl_body <- log(1.6)
    label("Typical body clearance CL_body at the cohort-median eGFR_MDRD (L/h)")           # Table 2, Body clearance (theta_CL) = 1.6 L/h (RSE 20%)

    # Power exponent of the eGFR_MDRD effect on body clearance. Table 2
    # footnote b: CL_body = theta_CL * (eGFR_MDRD/48.4)^theta_KF * exp(eta_CL).
    # Verified against Fig. 3B, whose curve runs from 1.06 L/h at
    # eGFR_MDRD = 30 to 6.43 L/h at 240 - exactly what this exponent gives.
    e_crcl_cl_body <- 0.869
    label("Power exponent on (CRCL/48.4) for body clearance (unitless)")                   # Table 2, Kidney function (theta_KF) = 0.869 (RSE 25%)

    # Central volume of distribution. No covariate effect was retained on it.
    lvc <- log(23.1)
    label("Typical central volume of distribution Vc (L)")                                 # Table 2, Volume of distribution / Central (theta_VC) = 23.1 L (RSE 10%)

    # Peripheral volume of distribution at the time of the first dose.
    lvp <- log(15.4)
    label("Typical peripheral volume of distribution Vp at time of first dose (L)")        # Table 2, Volume of distribution / Peripheral (theta_VP) = 15.4 L (RSE 18%)

    # Linear slope of the peripheral volume against time since the first dose,
    # PER MINUTE. Table 2 footnote c: Vp = theta_VP * (1 + TSFD * theta_T) *
    # exp(eta_VP), with the same footnote defining TSFD as "time since the first
    # dose in minutes". Results states the effect as "increased linearly by
    # 0.07% per minute", i.e. 0.0007/min, which agrees. Verified against
    # Fig. 3D, which runs from 15.4 L at time 0 to 93 L at 120 h
    # (15.4 * (1 + 7200 * 0.0007) = 93.0).
    e_tsfd_vp <- 0.0007
    label("Linear slope on peripheral volume per minute since the first dose (1/min)")     # Table 2, Time after the first dose (theta_T) = 0.0007 min^-1 (RSE 28%)

    # Intercompartmental clearance. No covariate effect and no IIV were
    # retained on it.
    lq <- log(12.0)
    label("Typical intercompartmental clearance Q (L/h)")                                  # Table 2, Intercompartmental clearance (theta_Q) = 12.0 L/h (RSE 17%)

    # Dialysis clearance intercept, i.e. the added clearance while KRT runs at
    # the 42 mL/min cohort-median dialysate flow rate. Table 2 names the
    # parameter theta_KRT in the row and theta_DIAL in footnote d; they are the
    # same parameter.
    lcl_hemodialysis <- log(2.0)
    label("Typical dialysis clearance CL_KRT at a dialysate flow rate of 42 mL/min (L/h)") # Table 2, Dialysis clearance (theta_KRT) = 2.0 L/h (RSE 15%)

    # Power exponent of the dialysate-flow-rate effect on dialysis clearance.
    # Table 2 footnote d: CL_KRT = theta_DIAL * (Q_D/42)^theta_QD. Verified
    # against Fig. 3C, whose curve runs from 1.64 L/h to 5.96 L/h over an
    # x-range of 30-270 - which this exponent reproduces exactly when that axis
    # is read in mL/min rather than the L/h it is labelled with (see the DFR
    # covariate notes and the vignette Errata).
    e_dfr_cl_hemodialysis <- 0.587
    label("Power exponent on (DFR/42) for dialysis clearance (unitless)")                  # Table 2, Q_D (theta_QD) = 0.587 (RSE 16%)

    # ---------------------------------------------------------------------
    # Interindividual variability. Table 2 reports these three as %CV under a
    # "Random effects: IIV" heading, and its footnote e gives the conversion
    # from the NONMEM variance as sqrt(exp(omega^2) - 1). The variances below
    # are therefore log(CV^2 + 1). The footnote's trailing "omega = variance"
    # is a typo for "omega^2 = variance": as printed it contradicts the very
    # formula it annotates, which squares omega inside the exponential. The
    # formula itself is legible only in the publisher's equation artwork
    # (the inline image behind Table 2 footnote e) - both pdftotext and the
    # markdown conversion silently drop the radical sign and the superscript.
    # nlmixr2 takes the VARIANCE on the eta line, so these are the NONMEM
    # OMEGA diagonal elements.
    #
    # Exponential IIV throughout: "Interindividual variability (IIV) was
    # explored using exponential random-effects models" (Methods, Data
    # analysis), consistent with the exp(eta_CL) and exp(eta_VP) terms printed
    # in Table 2 footnotes b and c.
    etalcl_body ~ 0.53983
    label("IIV on body clearance (variance; 84.6 %CV)")                                    # Table 2, Body clearance (eta_CL) = 84.6 %CV (RSE 14%); log(0.846^2 + 1)

    etalvc ~ 0.44437
    label("IIV on central volume of distribution (variance; 74.8 %CV)")                    # Table 2, Volume of distribution / Central (eta_VC) = 74.8 %CV (RSE 12%); log(0.748^2 + 1)

    # The base-model IIV on Vp was 247.0 %CV; the time-since-first-dose effect
    # "explained substantial parts of the associated IIV, reducing it from
    # 247.0 %CV to 71.1 %CV" (Results). Only the final 71.1 %CV is encoded.
    etalvp ~ 0.40914
    label("IIV on peripheral volume of distribution (variance; 71.1 %CV)")                 # Table 2, Volume of distribution / Peripheral (eta_VP) = 71.1 %CV (RSE 32%); log(0.711^2 + 1)

    # ---------------------------------------------------------------------
    # Residual unexplained variability: combined proportional plus additive.
    # These two rows sit under a separate "Random effects: residual variability"
    # heading and carry NO footnote-e superscript, so the proportional term is a
    # plain CV on the linear scale and needs no back-transformation, and the
    # additive term is already an SD in the concentration units of the assay.
    # ---------------------------------------------------------------------
    propSd <- 0.147
    label("Proportional residual error (fraction, i.e. 14.7 %CV)")                         # Table 2, Random effects: residual variability / Proportional = 14.7 %CV (RSE 20%)

    addSd <- 21.9
    label("Additive residual error (mg/L)")                                                # Table 2, Random effects: residual variability / Additive = 21.9 mg/L (RSE 19%)
  })

  model({
    # ---- 1. Derived covariate terms ------------------------------------------
    # Anuria gate, from the paper's own < 100 mL/24 h definition (Table 1
    # footnote a). It acts in BOTH directions: anuric patients lose the body
    # clearance arm and gain the expanding peripheral volume.
    anuric   <- (URINE_VOL_24H < 100)
    diuresis <- 1 - anuric

    # Time since the first dose, converted from the model's hours to the minutes
    # in which theta_T was estimated (Table 2 footnote c defines TSFD as "time
    # since the first dose in minutes"). TSFD is taken to be the model's own time
    # variable, so an event table MUST place the first dose at time = 0 or the
    # peripheral-volume trajectory will be wrong.
    tsfd_min <- t * 60

    # ---- 2. Individual parameters --------------------------------------------
    # Body clearance arm. Table 2 footnote b, in full:
    #   CL_body = theta_CL * (eGFR_MDRD/48.4)^theta_KF * exp(eta_CL)
    #             * 0 if study B * 0 if 24-h urine output < 100 mL
    # Note (CRCL/48.4)^e_crcl_cl_body is 0 at CRCL = 0, so a patient simulated at
    # eGFR_MDRD = 0 - the paper's own anuria scenario in Fig. S6 - already has no
    # body clearance before the two gates are applied.
    cl_body <-
      diuresis * (1 - STUDY_GERECKE) *
      exp(lcl_body + etalcl_body) * (CRCL / 48.4)^e_crcl_cl_body

    # Dialysis clearance arm. Table 2 footnote d:
    #   CL_KRT = theta_DIAL * (Q_D/42)^theta_QD
    # gated by whether KRT is running at all (Methods, Model development).
    # DFR is in mL/min, matching the 42 mL/min normalisation constant. Gotz 2025
    # deliberately avoided the Michaels equation used by the sibling
    # Huppe_2023_fosfomycin.R, so no dialyzer-specific K0A appears here and there
    # is no singular rational expression to guard - the arm is a plain power
    # function and is finite for every non-negative DFR.
    cl_hemodialysis <-
      RRT_CRRT_ACTIVE * exp(lcl_hemodialysis) * (DFR / 42)^e_dfr_cl_hemodialysis

    # "A two-compartment model incorporating concomitant body and dialysis
    # clearance" (Results); Fig. 3A draws the two arms summing to total clearance.
    cl_total <- cl_body + cl_hemodialysis

    vc <- exp(lvc + etalvc)

    # Peripheral volume, expanding linearly with time since the first dose in
    # anuric patients only (Results; see the URINE_VOL_24H covariate notes for
    # why the gate is present even though Table 2 footnote c omits it).
    vp <- exp(lvp + etalvp) * (1 + anuric * e_tsfd_vp * tsfd_min)

    q <- exp(lq)

    # ---- 3. Micro-constants --------------------------------------------------
    kel <- cl_total / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- 4. ODE system -------------------------------------------------------
    # Standard two-compartment disposition with intravenous input into the
    # central compartment (Fig. 3A). Written on amounts.
    #
    # NOTE on the degenerate zero-clearance limit, which this model genuinely
    # predicts and which the paper genuinely simulates: an anuric patient
    # (eGFR_MDRD = 0) who is not currently on KRT has cl_total EXACTLY zero.
    # Concentrations still plateau rather than diverge, because the peripheral
    # volume grows with time since the first dose - this is the dark-orange
    # "eGFR_MDRD = 0" trace of Fig. S6.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ---- 5. Observation and error --------------------------------------------
    # Dose in mg and vc in L give mg/L, the units the paper reports throughout
    # (MICs of 32-256 mg/L; additive residual error of 21.9 mg/L).
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
