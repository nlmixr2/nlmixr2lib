Park_2025_efineptakin_alfa <- function() {
  description <- "Joint population PK/PD model of rhIL-7-hyFc (efineptakin alfa), a hybrid-Fc-fused long-acting recombinant human interleukin-7, given by intramuscular injection to adults with locally advanced or metastatic solid tumours (Park 2025, final sequential PK-PD model). PK is two-compartment with linear disposition and a DOUBLE-PEAK absorption phase built from TWO parallel depot compartments: the dose is split between them by a logistic determinant (BIOF) and depot 2 is delayed by an absorption lag. Clearance is PIECEWISE-CONSTANT IN TIME -- a hard step from a sex-dependent early value to a single sex-independent late value at an estimated breakpoint TCLchange = 350 h after treatment initiation (NONMEM MTIME). Because the ELISA cannot separate endogenous interleukin-7 from drug, the observed serum concentration is the model-predicted drug concentration PLUS a per-subject endogenous IL-7 baseline. PD is a Friberg-type lymphopoiesis chain -- one progenitor pool, three maturation transit compartments and a circulating pool observed as the absolute lymphocyte count (ALC) -- with power-law homeostatic feedback (Circ0 / Circ)^gamma and Emax STIMULATION (not suppression) of progenitor proliferation by drug. Estimated on 402 serum concentrations and 256 ALC observations from 75 treatment cycles in 35 patients (NCT03478995 phase 1b). Two outputs: total serum IL-7 (ng/mL) and ALC (cells/uL)."
  reference <- paste(
    "Park S, Lee SM, Pak KC, Byun M-S, Choi D, Lim H-S. (2025).",
    "Population Pharmacokinetic and Pharmacodynamic Modeling Analysis of rhIL-7-hyFc,",
    "a Hybrid Fc-Fused Long-Acting Interleukin-7, to Support Optimal Dosing Regimens",
    "in Patients with Solid Cancer.",
    "Drug Des Devel Ther 19:9303-9320. doi:10.2147/DDDT.S564085. PMC12539412.",
    "The lymphopoiesis structure is Friberg LE, Henningsson A, Maas H, Nguyen L, Karlsson MO. (2002).",
    "Model of chemotherapy-induced myelosuppression with parameter consistency across drugs.",
    "J Clin Oncol 20(24):4713-4721. doi:10.1200/JCO.2002.02.140.",
    sep = " "
  )
  vignette <- "Park_2025_efineptakin_alfa"
  # Declared explicitly: buildModelDb()'s fallback heuristic matches only the
  # literal names `depot` and `central`, so without this it would record
  # `dosing = depot,central` -- true of `depot`, but silently 1 short (it misses
  # `depot2`) and wrong about `central`, which is never dosed here. Both depots
  # receive the full amount; the f() values do the splitting.
  dosing <- c("depot", "depot2")
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL",
    alc           = "cells/uL"
  )
  # UNITS, established arithmetically (the paper never states the dose base).
  # Doses are reported in mg/kg and as absolute amounts ("0.06 to 1.7 mg/kg;
  # equivalent to 2.8-133.5 mg"), volumes in L (Vc 1290 L, Vp 15000 L) and
  # concentrations in ng/mL (EC50 0.066 ng/mL, LLOQ 0.031 ng/mL). mg / L is
  # ug/mL, so the central concentration must be multiplied by 1000 to reach
  # ng/mL. Corroboration against the paper's own steady-state simulation, which
  # pins the concentration base exactly. At steady state every dosing interval
  # begins well past TCLchange, so the late arm alone governs and the TYPICAL
  # Cavg is Dose / (CL_late * tau). At the 61.7 kg pooled cohort mean weight
  # that is 0.034 / 0.408 / 0.816 ng/mL for 0.05 / 0.6 / 1.2 mg/kg every 12
  # weeks, against the reported 0.05 / 0.54 / 1.09. The gap is a CONSTANT
  # 1.34-1.36x at every dose, and it is not slack: the paper reports the MEAN
  # across a simulated cohort, and Cavg is proportional to 1/CL, so with a
  # log-normal CL the cohort mean is inflated over the typical value by exactly
  # exp(omega^2 / 2) = exp(0.613 / 2) = 1.3587. Applying it gives
  # 0.046 / 0.554 / 1.108 against 0.05 / 0.54 / 1.09 -- agreement to 2.6%, 2.6%
  # and 1.6%. A CONSTANT ratio explained by an independently reported variance
  # confirms the concentration base; a dose-dependent one would have indicated a
  # clearance error. The vignette reproduces this check by simulation.

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot       = list(analyte = "rhIL-7-hyFc (efineptakin alfa)", units = "mg", specimen = "administration site", verified = TRUE),
    depot2      = list(analyte = "rhIL-7-hyFc (efineptakin alfa)", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "rhIL-7-hyFc (efineptakin alfa)", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "rhIL-7-hyFc (efineptakin alfa)", units = "mg", specimen = "serum", verified = TRUE),
    precursor1  = list(analyte = "lymphocytes", units = "cells/uL", specimen = "not applicable", verified = TRUE),
    precursor2  = list(analyte = "lymphocytes", units = "cells/uL", specimen = "not applicable", verified = TRUE),
    precursor3  = list(analyte = "lymphocytes", units = "cells/uL", specimen = "not applicable", verified = TRUE),
    precursor4  = list(analyte = "lymphocytes", units = "cells/uL", specimen = "not applicable", verified = TRUE),
    circ        = list(analyte = "lymphocytes", units = "cells/uL", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    SEXF = list(
      description        = "Female sex indicator (1 = female, 0 = male).",
      units              = "(binary)",
      type               = "categorical",
      reference_category = "male (SEXF = 0)",
      notes              = "The ONLY covariate retained in the final PK or PD model, and it acts on the EARLY clearance arm only. Park 2025 Results: 'Sex was the only covariate influencing CL< TCLchange. The estimated clearance values were 4.82 L/h in females, 12.40 L/h in males (both CL<TCLchange), and 45.80 L/h CL> TCLchange in both sexes.' Male is taken as the reference because it carries the plain lcl; the female effect is therefore e_sexf_cl = log(4.82 / 12.40). No covariate acts on the LATE arm (lcl_late) -- the paper is explicit that the post-breakpoint clearance is the same in both sexes -- so the sex effect switches off at TCLchange along with the early arm. Cohort split: 19 male / 16 female of 35 (Table 1). The paper reports the sex effect on clearance had limited clinical consequence for the PD endpoint: predicted ALC at 1.2 mg/kg every 12 weeks was 2723 cells/uL in males versus 2709 cells/uL in females.",
      source_name        = "SEX"
    )
  )

  # Covariates SCREENED but not retained. Park 2025 evaluated 26 demographic and
  # clinical covariates by stepwise covariate modelling (forward p < 0.005,
  # backward p > 0.001) on both the PK and the PK-PD model, and every one except
  # sex was rejected. Only those with an existing entry in
  # inst/references/covariate-columns.md are listed here as structured entries;
  # the remainder of the screen is recorded verbatim in population$notes rather
  # than asserted against canonical names that do not exist in the register.
  covariatesDataExcluded <- list(
    AGE           = list(description = "Age", units = "years", type = "continuous", reference_category = NULL, notes = "Screened in the 26-covariate stepwise analysis; not retained. Cohort mean 57.9 years, CV 16.2% (Table 1); range 40-75 years; 31% (11/35) aged over 60."),
    WT            = list(description = "Body weight", units = "kg", type = "continuous", reference_category = NULL, notes = "Screened; not retained. NOTE this model carries NO allometric scaling -- clearance and volume are absolute, and mg/kg doses must be converted to an absolute mg amount before simulation. Cohort means 67.1 kg in males (range 49-110) and 55.3 kg in females (range 42-68) (Table 1 and Results text)."),
    HT            = list(description = "Height", units = "cm", type = "continuous", reference_category = NULL, notes = "Screened; not retained. Cohort means 169.7 cm male, 157.9 cm female (Table 1)."),
    SMOKE_CURRENT = list(description = "Current-smoker indicator", units = "(binary)", type = "categorical", reference_category = "non-smoker", notes = "Screened; not retained. Cohort: 1 smoker, 9 ex-smokers, 25 non-smokers (Table 1)."),
    WBC           = list(description = "White blood cell count", units = "10^9/L", type = "continuous", reference_category = NULL, notes = "Screened; not retained. Cohort mean 7.77, CV 40.3% (Table 1)."),
    RBC           = list(description = "Red blood cell count", units = "10^12/L", type = "continuous", reference_category = NULL, notes = "Screened; not retained. Cohort mean 3.80, CV 14.2% (Table 1)."),
    HGB           = list(description = "Hemoglobin", units = "g/dL", type = "continuous", reference_category = NULL, notes = "Screened; not retained. Cohort mean 11.86, CV 13.2% (Table 1)."),
    HCT           = list(description = "Hematocrit", units = "%", type = "continuous", reference_category = NULL, notes = "Screened; not retained. Cohort mean 35.61, CV 12.9% (Table 1).")
  )

  population <- list(
    species        = "human",
    n_subjects     = 35L,
    n_studies      = 1L,
    age_range      = "40-75 years (mean 57.9, CV 16.2%)",
    weight_range   = "42-110 kg (male mean 67.1, range 49-110; female mean 55.3, range 42-68)",
    height_range   = "male mean 169.7 cm, female mean 157.9 cm",
    sex_female_pct = 45.7,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Histologically confirmed, locally advanced, recurrent or metastatic solid tumours that were incurable and had failed or were unsuitable for standard therapy. ECOG performance status 0-1, age at least 19 years, adequate haematologic and end-organ function. Anti-tumour therapy within 3 weeks, or prior immune checkpoint inhibitor or immunomodulatory antibody within 12 weeks, were exclusions.",
    dose_range     = "Intramuscular rhIL-7-hyFc 0.06, 0.12, 0.24, 0.48, 0.72, 0.96, 1.2 and 1.7 mg/kg (equivalent to 2.8-133.5 mg absolute) every 3 weeks in 29 patients, plus 1.2 mg/kg every 6 weeks in 6 patients. 1.2 mg/kg was the maximum tolerated dose; one dose-limiting grade 3 or higher hypersensitivity occurred at 1.7 mg/kg.",
    regions        = "Republic of Korea (Yonsei Cancer Center, Asan Medical Center, Seoul St. Mary's Hospital)",
    observations   = "402 serum rhIL-7-hyFc concentration measurements and 256 ALC observations across 75 treatment cycles in 35 patients. All 35 contributed cycle-1 data; 26 (74%) contributed cycle-2 data; 6 contributed cycle-3 data; only 1 patient extended beyond three cycles. PK sampling at 0, 0.5, 6, 12, 24, 48, 72, 168 and 336 h in cycle 1, then 504, 672, 840 and 1008 h; ALC weekly to 1008 h then at 1512, 2016, 2520, 3024, 3528, 4032 and 4536 h. Assay LLOQ 0.031 ng/mL (Quantikine HS ELISA).",
    baseline_alc   = "Median 1268 cells/uL (range 341-2453); Table 1 mean 1314.4 cells/uL (CV 40.3%). 12 of 35 patients (34%) had ALC above 1500 cells/uL and 2 had severe lymphopenia below 500 cells/uL.",
    baseline_il7   = "Median endogenous serum IL-7 0.05 ng/mL (range 0.02-0.23). The ELISA measures total IL-7 and cannot separate endogenous from drug-derived, so this baseline is carried in the observation model rather than subtracted from the data.",
    notes          = "Estimation used NONMEM 7.4.3 with ADVAN13 and FOCE-I; PK and PD were fitted SEQUENTIALLY, the PD model taking individual empirical Bayes PK estimates as input. Precision from a 1000-replicate nonparametric bootstrap (PK 988/1000 and PD 997/1000 successful). Twenty-six covariates were screened by stepwise covariate modelling and only sex was retained: sex, age, body weight, height, smoking status, alcohol consumption, red blood cell count, hemoglobin, hematocrit, white blood cell count, neutrophil %, lymphocyte %, monocyte %, eosinophil %, basophil %, platelet count, mean corpuscular hemoglobin, mean corpuscular hemoglobin concentration, absolute neutrophil count, absolute lymphocyte count, prothrombin time, activated partial thromboplastin time, PT in international normalization ratio, thyroid-stimulating hormone, T3 and free T4. Anti-drug antibodies developed across all doses and dosing intervals but did not affect safety and were NOT modelled as a covariate on clearance -- the time-dependent clearance step is structural, and the paper reports the clearance increase was not dose-dependent (no correlation between Bayesian individual clearance and dose). Long-term data are thin: 27 of 35 patients (77%) received only one or two doses."
  )

  ini({
    # =========================================================================
    # PHARMACOKINETICS -- Park 2025 Table 2 (final PK model, point-estimate
    # column; the bootstrap-median column is NOT used). All IIV terms in Table 2
    # are reported as "omega^2 (CV %)" pairs and are on the VARIANCE scale: the
    # table footnote gives CV% = 100 * sqrt(exp(omega) - 1), and every printed
    # pair reproduces exactly under that reading (0.609 -> 91.6; 0.377 -> 67.7;
    # 0.198 -> 46.8; 1.81 -> 226.1; 0.613 -> 91.9). Read as standard deviations
    # they would not. So the tabulated numbers go into the eta blocks unchanged.
    # =========================================================================

    # -- Absorption. Two parallel depots reproduce the observed double peak:
    # "Introduction of two depot compartments in the absorption model, where
    # each of the two doses is administered divided by two fractions with a time
    # interval, well described the double peaks during the absorption phase."
    # ONE shared Ka serves both depots (Table 2 abbreviations: "Ka, absorption
    # rate constant from depot compartment 1 and 2 to central compartment").
    lka         <- log(0.341); label("ka: first-order absorption rate constant, shared by both depots (1/h)") # Park 2025 Table 2: Ka = 0.341 1/h
    logitfdepot <- 0.346;      label("BIOF: logit-scale determinant of the dose split between depot 1 and depot 2 (unitless)") # Park 2025 Table 2: BIOF = 0.346
    ltlag       <- log(22.5);  label("ALAG2: absorption lag time on depot 2 only (h)") # Park 2025 Table 2: ALAG2 = 22.5 h

    # BIOA is FIXED at 1 (Table 2 point estimate 1, with "NA" in both bootstrap
    # columns, and the abbreviations line states "BIOA, bioavailability of
    # rhIL-7-hyFc which is fixed at 1"). It is an identifiability anchor for an
    # intramuscular-only data set: no intravenous arm was studied, so Vc, Vp, Q
    # and both clearance arms are all apparent (per-F) quantities. The IIV on it
    # is estimated and is the paper's "IIV in bioavailability".
    lfdepot     <- fixed(log(1)); label("BIOA: total bioavailability across both depots (fraction)") # Park 2025 Table 2: BIOA = 1, fixed

    # -- Disposition. Two compartments, linear.
    lvc <- log(1290.0);  label("Vc/F: apparent central volume of distribution (L)")                # Park 2025 Table 2: V_C = 1290.0 L
    lvp <- log(15000.0); label("Vp/F: apparent peripheral volume of distribution (L)")             # Park 2025 Table 2: V_P = 15000.0 L
    lq  <- log(10.80);   label("Q/F: apparent inter-compartmental clearance (L/h)")                # Park 2025 Table 2: Q = 10.80 L/h

    # -- Clearance: a HARD STEP IN TIME, not a smooth time-course.
    # "The PK model indicated that clearance of rhIL-7-hyFc increases
    # approximately 350 hours (14 days) after treatment initiation ... captured
    # by implementing a stepwise increase in CL at an estimated time point
    # (TCLchange) within each patient in the model, using the MTIME feature in
    # NONMEM."
    #
    # NAMING. `tclchange` and `cl_late` are new canonical parameter names
    # registered with this model (inst/references/parameter-names.md). The
    # naming is deliberately ASYMMETRIC: the EARLY arm is an ordinary clearance
    # and so is carried by the plain canonical `lcl`, which lets the sex
    # covariate take the standard two-token `e_sexf_cl` shape. There is no
    # `lcl_early`; that would be a second canonical for a concept the register
    # already covers. The existing `lcl_ss` / `lcl_time` pair does NOT apply --
    # that is an ADDITIVE decomposition CL(t) = CL_ss + CL_time(t), whereas this
    # is a piecewise-constant switch between two mutually exclusive arms.
    ltclchange <- log(350.0); label("TCLchange: time after treatment initiation at which clearance steps up (h)") # Park 2025 Table 2: TCLchange = 350.0 h
    lcl        <- log(12.40); label("CL/F before TCLchange, male reference (L/h)")                                # Park 2025 Table 2: CL_M,<TCLchange = 12.40 L/h
    lcl_late   <- log(45.80); label("CL/F after TCLchange, both sexes (L/h)")                                     # Park 2025 Table 2: CL_>TCLchange = 45.80 L/h

    # Sex effect on the EARLY clearance arm only. Male is the reference (it
    # carries lcl), so the female multiplier is 4.82 / 12.40 = 0.3887, i.e. a
    # ~2.6-fold LOWER early clearance in females. The ratio is written out
    # inline so the arithmetic is auditable.
    e_sexf_cl <- log(4.82 / 12.40); label("Effect of female sex on log CL/F before TCLchange (unitless)") # Park 2025 Table 2: CL_F,<TCLchange = 4.82 L/h vs CL_M,<TCLchange = 12.40 L/h

    # -- Endogenous interleukin-7 baseline. The ELISA "does not distinguish
    # between endogenous IL-7 and exogenous rhIL-7-hyFc", so "baseline
    # (pre-dose) endogenous IL-7 levels were retained in the analysis without
    # substitution or correction, under the assumption that baseline levels
    # remained constant for each individual." It is therefore an additive
    # offset on the OBSERVATION, not a compartment and not a drug property.
    # Park 2025 reports no estimated typical value; the Results text gives the
    # cohort median directly.
    lbl_il7 <- log(0.05); label("Baseline endogenous serum interleukin-7 concentration (ng/mL)") # Park 2025 Results: "The median baseline serum IL-7 concentration was 0.05 ng/mL (range: 0.02-0.23)"

    # -- PK inter-individual variability (Table 2, variance scale; see header).
    etalka      ~ 0.609 # Park 2025 Table 2: IIV_Ka   = 0.609 (CV 91.6%)
    etalfdepot  ~ 0.377 # Park 2025 Table 2: IIV_BIOA = 0.377 (CV 67.7%)
    etalvc      ~ 0.198 # Park 2025 Table 2: IIV_Vc   = 0.198 (CV 46.8%)
    etalq       ~ 1.81  # Park 2025 Table 2: IIV_Q    = 1.81  (CV 226.1%)
    # A SINGLE combined random effect on clearance. "IOV was tested for each
    # fixed effect parameter; inclusion of IOV in clearance (CL) significantly
    # improved model fit. However, due to data limitations, IIV and IOV in CL
    # could not be estimated separately, and a combined random effect parameter
    # was estimated instead." It is one eta, so it multiplies BOTH clearance
    # arms (see model()); it cannot be split into an IIV and an IOV part here,
    # and encoding it as occasion-varying would overstate what was estimated.
    etalcl      ~ 0.613 # Park 2025 Table 2: IIV+IOV_CL = 0.613 (CV 91.9%)
    # IIV on the endogenous baseline. Park 2025 observation equation (page 9309)
    # carries a random effect omega2 on the baseline concentration and states
    # "omega_1 and epsilon_1 share the same variance, as do omega_2 and
    # epsilon_2, and the random effect omega_1 and omega_2 persist across all
    # observations for each individual." The variance of omega_2 is therefore
    # the variance of epsilon_2, which Table 2 reports as a 27.5% CV. Converted
    # to a log-normal omega^2 by the standard identity so the baseline stays
    # strictly positive -- see the vignette Errata for this deviation from the
    # paper's additive-on-one form.
    etalbl_il7  ~ log(0.275^2 + 1) # Park 2025 page 9309 + Table 2: var(omega_2) = var(epsilon_2), epsilon proportional = 0.275 as CV -> omega^2 = 0.07290

    # -- PK residual error. Table 2 gives ONE residual term, proportional,
    # 0.275, footnoted "a proportional residual error represented as CV". The
    # paper's observation equation also carries an ADDITIVE term epsilon_1, but
    # its value is reported nowhere in the paper, so it is omitted rather than
    # invented (recorded in the vignette Errata). The proportional term applies
    # to the TOTAL measured IL-7 (baseline + predicted drug), which is exactly
    # the epsilon_2 * (C_BSL,obs + C_pred) term of the published equation.
    propSd_Cil7 <- 0.275; label("Proportional residual error on total serum IL-7 (fraction)") # Park 2025 Table 2: epsilon (proportional) = 0.275, reported as CV

    # =========================================================================
    # PHARMACODYNAMICS -- Park 2025 Table 3 (final PD model, point-estimate
    # column). Same variance-scale reading as Table 2, and it reproduces the
    # printed CV% here too (0.127 -> 36.8; 3.900 -> 695.7; and 0.11, carrying
    # the rounding of a value near 0.1137, -> 34.7).
    # =========================================================================
    lmtt   <- log(146.0); label("MTT: mean transit time for lymphocyte maturation (h)")                              # Park 2025 Table 3: MTT = 146.0 h
    lgamma <- log(0.130); label("gamma: exponent of the (Circ0 / Circ) homeostatic-feedback term (unitless)")        # Park 2025 Table 3: gamma = 0.130
    lkout  <- log(0.051); label("Kcirc: first-order degradation rate constant of circulating lymphocytes (1/h)")     # Park 2025 Table 3: Kcirc = 0.051 1/h
    lemax  <- log(0.155); label("Emax: maximum fractional stimulation of progenitor proliferation (unitless)")       # Park 2025 Table 3: Emax = 0.155
    lec50  <- log(0.066); label("EC50: drug concentration giving half-maximal stimulation (ng/mL)")                  # Park 2025 Table 3: EC50 = 0.066 ng/mL

    # Baseline circulating ALC. This is NOT an estimated parameter: the PD model
    # was fitted sequentially and used each subject's own observed pre-dose ALC,
    # so Table 3 has no Circ0 row. The cohort mean from Table 1 is used as the
    # typical value. It is preferred over the Results-text median (1268
    # cells/uL) because inverting the paper's own simulation pins the simulated
    # cohort baseline nearer the mean: "2439-2898 cells/uL at 0.6 mg/kg
    # (increases of +1088 and +1535 from baseline, respectively)" implies
    # 2439 - 1088 = 1351 and 2898 - 1535 = 1363 cells/uL. NO IIV is attached --
    # the 40.3% in Table 1 is the observed CV of the data, not an estimated
    # omega, and inventing an omega from it would fabricate a variance the
    # paper never estimated.
    lrbase_alc <- log(1314.4); label("Circ0: baseline absolute lymphocyte count (cells/uL)") # Park 2025 Table 1: ALC mean = 1314.4 cells/uL (CV 40.3%)

    # -- PD inter-individual variability (Table 3, variance scale).
    etalmtt  ~ 0.127 # Park 2025 Table 3: IIV_MTT  = 0.127 (CV 36.8%)
    # WARNING -- THIS VARIANCE MAKES STOCHASTIC SIMULATION EXPLOSIVE. It is
    # transcribed exactly as published and is NOT to be quietly shrunk, but any
    # Monte Carlo use needs care. The saturating ceiling on the ALC response is
    # (1 + Emax)^(1 / gamma), and 1 / gamma = 1 / 0.130 = 7.69, so a log-normal
    # Emax with omega^2 = 3.900 (SD 1.975 on the log scale) is raised to the
    # 7.69th power: a +2 SD subject has Emax = 0.155 * exp(3.95) = 8.05 and a
    # ceiling of 9.05^7.69 = 2e7 times baseline. A 100-subject draw reaches ALC
    # around 1e8 cells/uL at the 95th percentile. The paper itself reports this
    # term as barely identified (bootstrap 95% CI 0.005-7.774, spanning three
    # orders of magnitude) and its own Monte Carlo means (2439-3165 cells/uL)
    # are reproduced by the TYPICAL-VALUE trajectory, not by a draw from this
    # variance -- see the vignette, which validates against typical values and
    # suppresses this eta for any cohort plot.
    etalemax ~ 3.900 # Park 2025 Table 3: IIV_Emax = 3.900 (CV 695.7%)
    etalec50 ~ 0.11  # Park 2025 Table 3: IIV_EC50 = 0.11  (CV 34.7%)
    # There is deliberately NO eta on gamma, Kcirc or Circ0: Table 3 lists no
    # IIV row for them, and "Including the IOV on any fixed-effect parameter did
    # not improve the model fit" for the PD model.

    # -- PD residual error, combined additive plus proportional. Table 3
    # footnotes are explicit about the scales: "epsilon_1, additive residual
    # error represented as standard deviation" and "epsilon_2, proportional
    # residual error represented as CV", so both go in unchanged. The additive
    # term is in cells/uL (Table 1 labels the ALC row "ALC, 1/uL"), which makes
    # 2.450 cells/uL negligible beside a ~1300 cells/uL baseline -- the
    # proportional term carries essentially all of the residual variability.
    addSd_ALC  <- 2.450; label("Additive residual SD on ALC (cells/uL)")      # Park 2025 Table 3: epsilon_1 (additive) = 2.450, as SD
    propSd_ALC <- 0.195; label("Proportional residual SD on ALC (fraction)")  # Park 2025 Table 3: epsilon_2 (proportional) = 0.195, as CV
  })

  model({
    # ---- 1. Individual PK parameters ----
    ka      <- exp(lka + etalka)
    vc      <- exp(lvc + etalvc)
    vp      <- exp(lvp)                 # population-only: Table 2 lists no IIV_Vp row
    q       <- exp(lq  + etalq)
    bioa    <- exp(lfdepot + etalfdepot) # iBIOA: individual total bioavailability, typical value fixed at 1
    tlag    <- exp(ltlag)
    bl_il7  <- exp(lbl_il7 + etalbl_il7)

    # ---- 2. Time-varying clearance: a piecewise-constant step ----
    # Reproduces the NONMEM MTIME construction. `t` is time since the start of
    # the record, which for these simulations is time since treatment
    # initiation -- the quantity the paper's TCLchange is measured from ("350
    # hours (14 days) after treatment initiation"), NOT time since the most
    # recent dose. The boolean-multiplier idiom follows Yoneyama_2017_emicizumab.
    #
    # The sex effect and the single combined CL random effect are applied
    # differently on purpose: e_sexf_cl acts on the EARLY arm only (the paper
    # states the late arm is identical in both sexes), while etalcl multiplies
    # BOTH arms because only one CL random effect was estimated.
    tclchange <- exp(ltclchange)
    cl_early  <- exp(lcl + etalcl + e_sexf_cl * SEXF)
    cl_late   <- exp(lcl_late + etalcl)
    cl        <- cl_early * (t < tclchange) + cl_late * (t >= tclchange)

    # ---- 3. Micro-constants ----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---- 4. PK ODEs: two parallel depots feeding a two-compartment system ----
    d/dt(depot)       <- -ka * depot
    d/dt(depot2)      <- -ka * depot2
    d/dt(central)     <-  ka * depot + ka * depot2 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ---- 5. Dose split and lag ----
    # Park 2025 Table 2 footnote, verbatim: "BIOF, determinant of relative
    # bioavailabilities between depot compartment 1 (F1) and 2 (F2), which is
    # expressed as following logistic model:
    #   F1 = iBIOA * e^BIOF / (1 + e^BIOF),  F2 = iBIOA * (1 - e^BIOF / (1 + e^BIOF))
    # where iBIOA is individual BIOA".
    # At the point estimate BIOF = 0.346 this splits the dose 0.5857 / 0.4143.
    # Both depots must be dosed with the FULL amount; the f() values do the
    # splitting. The lag applies to depot 2 only, which is what separates the
    # two absorption peaks.
    fr1        <- exp(logitfdepot) / (1 + exp(logitfdepot))
    f(depot)   <- bioa * fr1
    f(depot2)  <- bioa * (1 - fr1)
    alag(depot2) <- tlag

    # ---- 6. Drug concentration ----
    # Drug-derived serum concentration, ng/mL (see the units note above). This,
    # not the total, drives the PD: Park 2025 defines the stimulatory effect as
    # "E = Emax * CONC / (EC50 + CONC), where CONC is the serum concentration of
    # rhIL-7-hyFc". The drug-only reading is also structurally forced -- if CONC
    # included the endogenous baseline then E would be non-zero at t = 0
    # (0.155 * 0.05 / (0.066 + 0.05) = 0.067) and the lymphopoiesis chain would
    # not be at steady state at its own stated baseline, so the model would
    # drift away from Circ0 with no drug present.
    Cc <- central / vc * 1000

    # ---- 7. Individual PD parameters ----
    mtt       <- exp(lmtt + etalmtt)
    gamma     <- exp(lgamma)
    kout      <- exp(lkout)
    emax      <- exp(lemax + etalemax)
    ec50      <- exp(lec50 + etalec50)
    rbase_alc <- exp(lrbase_alc)

    # ktr from MTT. Park 2025 Table 3 footnote: "MTT = (n+1) / Ktr, where Ktr is
    # transit compartments and rate constant, and n is 3, the number of transit
    # compartments." So ktr = 4 / 146.0 = 0.027397 1/h -- numerically distinct
    # from Kcirc = 0.051 1/h, which is what proves the two are different
    # parameters (see the dCirc/dt note below).
    ktr <- 4 / mtt

    # ---- 8. Stimulatory drug effect ----
    edrug <- emax * Cc / (ec50 + Cc)

    # ---- 9. Friberg-type lymphopoiesis chain ----
    # Compartment naming follows the house Friberg precedent
    # (Gebhard_2023_mercaptopurine_anc, Kawamura_2018_eribulin,
    # Guo_2022_PF_06939999): Prol -> precursor1, Transit1..3 ->
    # precursor2..precursor4, Circ -> circ.
    #
    # THREE DEPARTURES FROM THE EQUATIONS AS TYPESET ON PAGE 9310. Each is
    # forced by the paper's own numbers; none is a guess. All three are
    # restated in the vignette Errata.
    #
    # (i) THE FEEDBACK TERM IS PRINTED UPSIDE DOWN. The typeset equation reads
    #     dProl/dt = KProl * Prol * (1 + E) * (Circ / Circ0)^gamma - Ktr * Prol,
    #     with Circ in the NUMERATOR. That direction is arithmetically
    #     incompatible with the paper's own results. Setting dProl/dt = 0 with
    #     KProl = Ktr gives Circ_ss / Circ0 = (1 + Emax)^(-1/gamma) =
    #     1.155^(-1/0.130) = 0.33 -- i.e. the drug would cut ALC to a THIRD of
    #     baseline. The paper reports the opposite everywhere: "a
    #     dose-dependent increase in ALC", "approximately 1.8- to 2.3-fold
    #     increase from baseline", "+1088 to +1802 cells/uL". Inverting the
    #     ratio gives Circ_ss / Circ0 = 1.155^(1/0.130) = 3.03, a maximum
    #     3.03-fold rise which brackets the reported 1.8-2.3-fold time-averages
    #     exactly as a sub-Emax average should. The inverted form is also what
    #     the cited source model uses (Friberg 2002), what the Discussion
    #     describes ("homeostatic feedback ... constrain lymphocyte
    #     expansion"), and what every Friberg model already in this library
    #     uses. Implemented as (Circ0 / Circ)^gamma.
    #
    # (ii) THE LAST ODE PRINTS Ktr WHERE Kcirc BELONGS. The typeset dCirc/dt
    #      reads "= Ktr * Transit3 - Ktr * Circ". Taken literally, Kcirc --
    #      an ESTIMATED parameter with a point estimate (0.051 1/h), a
    #      bootstrap median (0.039) and a 95% CI (0.026-0.158) -- would appear
    #      nowhere in the model at all. Both the Table 3 footnote and the
    #      Figure 2 legend define Kcirc as the "degradation rate constant of
    #      circulating lymphocytes", and 0.051 differs from ktr = 0.0274, so
    #      they are distinct quantities. (The prose beneath the equations
    #      compounds the error by calling Kcirc "the transit rate constant",
    #      contradicting both the table footnote and the figure legend.)
    #      Implemented with kout = Kcirc on the circulating pool.
    #
    # (iii) KProl IS NOT REPORTED, BUT IS STRUCTURALLY FORCED. Table 3 has no
    #       KProl row. Requiring the undosed system to sit at its own baseline,
    #       dProl/dt = 0 at Circ = Circ0 and E = 0, gives Prol * (KProl - Ktr)
    #       = 0 and hence KProl = Ktr identically. This is a derivation, not an
    #       assumption, and it is the same constraint the Friberg family
    #       imposes throughout this library.
    fb <- (rbase_alc / circ)^gamma
    d/dt(precursor1) <- ktr * precursor1 * (1 + edrug) * fb - ktr * precursor1
    d/dt(precursor2) <- ktr * (precursor1 - precursor2)
    d/dt(precursor3) <- ktr * (precursor2 - precursor3)
    d/dt(precursor4) <- ktr * (precursor3 - precursor4)
    d/dt(circ)       <- ktr * precursor4 - kout * circ

    # ---- 10. Initial conditions ----
    # The lymphopoiesis chain starts at the drug-free steady state implied by
    # the parameters themselves, so an undosed simulation holds flat at Circ0.
    # Solving the chain at rest: every precursor equals Prol0, and
    # dCirc/dt = ktr * Prol0 - kout * Circ0 = 0 forces
    # Prol0 = Circ0 * kout / ktr. At the typical values that is
    # 1314.4 * 0.051 / 0.027397 = 2447 cells/uL.
    precursor1(0) <- rbase_alc * kout / ktr
    precursor2(0) <- rbase_alc * kout / ktr
    precursor3(0) <- rbase_alc * kout / ktr
    precursor4(0) <- rbase_alc * kout / ktr
    circ(0)       <- rbase_alc

    # ---- 11. Observations ----
    # Two endpoints. The measured serum quantity is TOTAL IL-7, because the
    # Quantikine HS ELISA "does not distinguish between endogenous IL-7 and
    # exogenous rhIL-7-hyFc". Park 2025 page 9309 gives the observation model as
    #   C_obs = C_pred + (C_BSL,obs + omega_1 + omega_2 * C_BSL,obs)
    #                  + epsilon_1 + epsilon_2 * (C_BSL,obs + C_pred)
    # so the proportional residual applies to the SUM, which is what
    # `Cil7 ~ prop(propSd_Cil7)` encodes. Cc is exposed without residual error
    # as the drug-only concentration, since that is the quantity the PD model
    # consumes and the quantity the paper's exposure-response summaries report.
    Cil7 <- Cc + bl_il7
    ALC  <- circ
    Cil7 ~ prop(propSd_Cil7)
    ALC  ~ add(addSd_ALC) + prop(propSd_ALC)
  })
}
