Panetta_2024_palbociclib_plt <- function() {
  description <- "Semi-mechanistic PK/PD model of palbociclib-induced thrombocytopenia in children and young adults (4.9-21.6 years) with recurrent, progressive or refractory brain tumors treated on the Pediatric Brain Tumor Consortium phase 1 study PBTC-042 (Panetta 2024, Table 3B). Structurally identical to the neutrophil model modellib('Panetta_2024_palbociclib_anc') -- section 2.4.1 states 'we used a similar structural model for platelet dynamics' -- but separately estimated on the platelet counts, so none of the parameters are shared. The palbociclib PK layer is the one-compartment first-order-absorption model with an absorption lag time and an AST power effect on CL/F reported in Table 2 (also available on its own as modellib('Panetta_2024_palbociclib')), converted from hours to days. Plasma palbociclib inhibits the proliferation rate of a bone-marrow precursor pool through a fractional Imax term IC50^n / (IC50^n + C^n); the pool feeds three maturation transit compartments and then the circulating platelet compartment, which exerts negative feedback on proliferation through a KM / (KM + PLT) term. KM and the drug-free initial conditions are fixed by the paper's steady-state requirement (Equations 6 and 7). The transit rate constant is parameterised as the estimated fraction fk_bp of the proliferation rate kin, on a logit scale, which guarantees kin > k_bp for every simulated subject. The circulating-platelet elimination rate kout and the Hill coefficient n are fixed. Volumes and clearances are per m^2 of body-surface area and doses are supplied in mg/m^2, so BSA cancels and no BSA covariate is needed."
  reference <- paste(
    "Panetta JC, Selvo NS, Van Mater D, Stewart CF. (2024).",
    "Population Pharmacokinetic and Pharmacodynamic Study of Palbociclib in Children and Young Adults",
    "with Recurrent, Progressive, or Refractory Brain Tumors.",
    "Pharmaceutics 16(12):1528. doi:10.3390/pharmaceutics16121528.",
    "The myelosuppression structure was previously described in",
    "Panetta JC, Schaiquevich P, Santana VM, Stewart CF. (2008) Using pharmacokinetic and",
    "pharmacodynamic modeling and simulation to evaluate importance of schedule in topotecan therapy",
    "for pediatric neuroblastoma. Clin Cancer Res 14(1):318-325. doi:10.1158/1078-0432.CCR-07-1243,",
    "itself adapted from a temozolomide myelosuppression model.",
    "The underlying phase 1 trial is PBTC-042 (NCT02255461).",
    sep = " "
  )
  vignette <- "Panetta_2024_palbociclib"
  units <- list(
    time          = "day",
    dosing        = "mg/m^2",
    concentration = "ng/mL",
    plt           = "10^9/L"
  )
  # Unit notes are identical to the neutrophil model: the PK parameters are
  # converted from hours to days inline, and doses are BSA-normalized so
  # central / vc is mg/L and Cc multiplies by 1000 to reach ng/mL. Equation 1
  # is then driven by Cc in ng/mL: Panetta 2024 labels IC50 and C_palbociclib
  # "nM", but that reading contradicts the paper's own Figure 4 and Table 4 --
  # see the IC50 entry in ini().

  compartmentData <- list(
    depot      = list(analyte = "palbociclib", units = "mg/m^2", specimen = "administration site", verified = TRUE),
    central    = list(analyte = "palbociclib", units = "mg/m^2", specimen = "plasma", verified = TRUE),
    precursor1 = list(analyte = "platelets", units = "10^9/L", specimen = "tissue", verified = TRUE),
    precursor2 = list(analyte = "platelets", units = "10^9/L", specimen = "tissue", verified = TRUE),
    precursor3 = list(analyte = "platelets", units = "10^9/L", specimen = "tissue", verified = TRUE),
    precursor4 = list(analyte = "platelets", units = "10^9/L", specimen = "tissue", verified = TRUE),
    circ       = list(analyte = "platelets", units = "10^9/L", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    AST = list(
      description        = "Baseline serum aspartate aminotransferase activity.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL/F centred on the cohort median AST of 25 U/L (Table 1), the only covariate retained by the pharmacokinetic covariate analysis (Table 2, AST column; Table S1). See modellib('Panetta_2024_palbociclib') for the full covariate-screen description.",
      source_name        = "AST"
    ),
    OCC = list(
      description        = "Treatment-course index used by the inter-occasion random effects.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Panetta 2024 estimated the pharmacodynamic parameters on courses 1 and 2 (section 3.5), so OCC = 1 is course 1 and OCC = 2 is course 2. The paper reports inter-occasion variability on kin, k_bp and IC50 (Table 3B) but never prints an occasion column name; OCC is the nlmixr2lib canonical. The CL/F inter-occasion term of the PK layer is keyed to the same OCC here: in the PK analysis its occasions were the course 1 day 1 and course 1 day 21 sampling visits, roughly three weeks apart and therefore comparable in span to a treatment course. For a single-occasion simulation set OCC = 1 throughout.",
      source_name        = NA_character_
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 31L,
    n_studies      = 1L,
    age_range      = "4.9-21.6 years (median 12.8)",
    weight_range   = "23.8-110.4 kg (median 52.5)",
    bsa_range      = "0.8-2.4 m^2 (median 1.5)",
    sex_female_pct = 36,
    race_ethnicity = c(Caucasian = 58, Black = 12, Other = 30),
    disease_state  = "Recurrent, progressive or refractory pediatric brain tumors. Stratum I (not heavily pretreated) versus stratum II (heavily pretreated); simulations at 75 mg/m^2 predicted grade 3 or greater thrombocytopenia in 15% of stratum II versus 0% of stratum I patients (Figure S5).",
    dose_range     = "Oral palbociclib 50, 75 or 95 mg/m^2 once daily for 21 days of each 28-day course.",
    regions        = "United States (Pediatric Brain Tumor Consortium)",
    observations   = "189 platelet counts from 31 patients during course 1 and 99 during course 2; 55 further counts from six patients during courses 3 and 4 were held out to test prediction. Of the 33 patients in the PK analysis, 2 were excluded from the PD analysis because they had only a single pre-treatment platelet observation.",
    notes          = "Estimated in Monolix 2023R1 by SAEM with the palbociclib pharmacokinetics fixed to each individual's post hoc estimates (section 2.4.2). Observed thrombocytopenia was mild: only three of 31 patients reached grade 1 or worse, and the median (5th-95th percentile) nadir from the individual model estimates was 169 (59, 240) x 10^9/L. Unlike the neutrophil model, the courses 3 and 4 platelet predictions were biased, over-predicting depletion by a median of 8.0% (p = 4 x 10^-4)."
  )

  ini({
    # =========================================================================
    # PALBOCICLIB PK LAYER. Panetta 2024 Table 2, AST column (the final model),
    # converted from hours to days. Section 2.4.2: "The palbociclib
    # pharmacokinetics were fixed to each individual's post hoc estimated
    # parameters".
    # =========================================================================
    ltlag <- log(0.8 / 24); label("Tlag: absorption lag time (day)")                    # Table 2, AST column: Tlag = 0.8 h  = 0.033333 day
    lka   <- log(0.46 * 24); label("ka: first-order absorption rate constant (1/day)")  # Table 2, AST column: ka = 0.46 /h  = 11.04 /day
    lvc   <- log(653); label("V/F: apparent oral volume per BSA (L/m^2)")               # Table 2, AST column: V/F = 653 L/m^2
    lcl   <- log(36.5 * 24); label("CL/F: apparent oral clearance per BSA (L/day/m^2)") # Table 2, AST column: CL/F = 36.5 L/h/m^2 = 876 L/day/m^2

    e_ast_cl <- -0.3; label("Exponent of the (AST / 25 U/L) power effect on CL/F (unitless)")  # Table 2, AST column: AST on CL/F = -0.3 (RSE 39.3%); p = 0.0066

    # =========================================================================
    # PLATELET PD LAYER. Panetta 2024 Table 3B.
    # =========================================================================
    lrbase <- log(267.40); label("Ncirc0: drug-free steady-state circulating platelet count (10^9/L)")  # Table 3B: Ncirc0 = 267.40 x 10^9/L (RSE 4.2%)
    lkin   <- log(0.93); label("kin: proliferation rate of the bone-marrow platelet precursor pool (1/day)")  # Table 3B: kin = 0.93 /day (RSE 3.3%)
    # IC50 UNITS. Panetta 2024 labels IC50 and C_palbociclib "nM" (Table 3B and
    # the Equation 1-5 glossary), but the value is only self-consistent with
    # the paper's own PK model and simulations when read as ng/mL, the unit the
    # PK model actually produces. This model is the check that settles it: read
    # as ng/mL it reproduces Figure 4D, 4E and 4F at all three dosages (typical
    # nadirs 191, 163 and 143 x 10^9/L against roughly 205, 175 and 160 read
    # off the published median curves) and it keeps the typical subject well
    # above the 100 x 10^9/L threshold that Table 4 says only 6, 8 and 12% of
    # individuals cross. Applying the nM conversion (palbociclib is C24H29N7O2,
    # 447.54 g/mol, so 1 ng/mL = 2.2344 nM) means the tabulated 559.14 nM is
    # only 559.14 / 2.2344 = 250.2 ng/mL on the scale the PK model produces,
    # i.e. a MORE potent drug: typical nadirs become 128, 89 and
    # 67 x 10^9/L, putting the typical subject below the threshold at every
    # dosage and contradicting both the figure and the table. See
    # modellib('Panetta_2024_palbociclib_anc') for the third check, against
    # Supplementary Figure S2_5. This is a reporting error in the source,
    # recorded in the vignette's assumptions section.
    lic50  <- log(559.14); label("IC50: palbociclib concentration halving marrow proliferation (ng/mL)")  # Table 3B: IC50 = 559.14 (RSE 11.2%), printed as nM; see the units note above

    # As in the neutrophil model, k_bp is not estimated directly: section 2.4.1
    # defines k_bp = fk_bp * kin with fk_bp in [0, 1] and a logit distribution
    # (section 2.4.2), so that kin > k_bp always holds. Table 3B tabulates the
    # derived typical k_bp = 0.52 /day, so
    #   fk_bp = 0.52 / 0.93 = 0.559140,  logit(0.559140) = 0.237618.
    # The implied mean transit time from marrow to circulation is
    # 4 / k_bp = 4 / 0.52 = 7.69 days.
    logitfktr <- qlogis(0.52 / 0.93)
    label("logit(fk_bp): logit of the transit rate constant expressed as a fraction of kin (unitless)")  # Table 3B: k_bp = 0.52 /day (RSE 5.9%) with kin = 0.93 /day

    lkout <- fixed(log(4.80)); label("kout: elimination rate of circulating platelets (1/day)")  # Table 3B: kout = 4.80 /day, reported as Fixed
    lhill <- fixed(log(1));    label("n: Hill coefficient of the palbociclib marrow effect (unitless)")  # Table 3B: n = 1, reported as Fixed

    # =========================================================================
    # Inter-individual and inter-occasion variability. Table 3B reports these
    # as CV%, and section 2.4.2 states the random effects are log-normal (all
    # parameters except fk_bp), so omega^2 = log(CV^2 + 1). The inversion is
    # written out inline so it is auditable.
    # =========================================================================
    etalrbase ~ log(0.214^2 + 1)  # Table 3B: IIV Ncirc0 = 21.4% (RSE 14.4%) -> omega^2 = 0.04478
    etalkin   ~ log(0.076^2 + 1)  # Table 3B: IIV kin    = 7.6%  (RSE 39.4%) -> omega^2 = 0.00576
    etalic50  ~ log(0.330^2 + 1)  # Table 3B: IIV IC50   = 33.0% (RSE 33.6%) -> omega^2 = 0.10306

    # fk_bp is logit-distributed. Table 3B's CV is converted with the delta
    # method for the logit transform, sd_logit = CV / (1 - p) with
    # p = 0.559140:
    #   IIV: 0.076 / 0.440860 = 0.17239,  IOV: 0.059 / 0.440860 = 0.13383.
    # Round-trip check: expit(0.237618 +/- 0.17239) gives k_bp between 0.484
    # and 0.556 /day, a 7.6% CV.
    etalogitfktr ~ 0.17239^2  # Table 3B: IIV k_bp = 7.6% (RSE 99.1%)

    etaiov_kin_1  ~ log(0.041^2 + 1)
    label("Inter-occasion variability in kin, occasion 1 (log-scale variance)")   # Table 3B: IOV kin  = 4.1%  (RSE 43.8%) -> 0.00168
    etaiov_kin_2  ~ fix(log(0.041^2 + 1))                                         # same IOV magnitude on occasion 2
    etaiov_ic50_1 ~ log(0.196^2 + 1)
    label("Inter-occasion variability in IC50, occasion 1 (log-scale variance)")  # Table 3B: IOV IC50 = 19.6% (RSE 72.8%) -> 0.03769
    etaiov_ic50_2 ~ fix(log(0.196^2 + 1))                                         # same IOV magnitude on occasion 2
    etaiov_fktr_1 ~ 0.13383^2
    label("Inter-occasion variability in logit(fk_bp), occasion 1 (logit-scale variance)")  # Table 3B: IOV k_bp = 5.9% (RSE 57.2%)
    etaiov_fktr_2 ~ fix(0.13383^2)                                                 # same IOV magnitude on occasion 2

    # PK-layer variability carried over from Table 2.
    etaltlag ~ 0.67^2  # Table 2, AST column: IIV Tlag = 0.67 (reported as a standard deviation)
    etalka   ~ 0.87^2  # Table 2, AST column: IIV ka   = 0.87
    etalvc   ~ 0.04^2  # Table 2, AST column: IIV V/F  = 0.04
    etalcl   ~ 0.26^2  # Table 2, AST column: IIV CL/F = 0.26
    etaiov_cl_1 ~ 0.12^2
    label("Inter-occasion variability in CL/F, occasion 1 (log-scale variance)")  # Table 2, AST column: IOV CL/F = 0.12
    etaiov_cl_2 ~ fix(0.12^2)                                                     # same IOV magnitude on occasion 2

    # Residual error. Section 2.4.2: "A proportional residual error model was
    # used". Table 3B reports it as b = 0.20. The palbociclib concentrations
    # are not an endpoint of this fit -- the PK was fixed to individual post hoc
    # estimates -- so only the platelet count carries an error model here.
    propSd_PLT <- 0.20; label("Proportional residual error on platelet count (fraction)")  # Table 3B: residual error b = 0.20 (RSE 6.0%)
  })

  model({
    # ---- 1. Occasion indicators ----
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_cl   <- oc1 * etaiov_cl_1   + oc2 * etaiov_cl_2
    iov_kin  <- oc1 * etaiov_kin_1  + oc2 * etaiov_kin_2
    iov_ic50 <- oc1 * etaiov_ic50_1 + oc2 * etaiov_ic50_2
    iov_fktr <- oc1 * etaiov_fktr_1 + oc2 * etaiov_fktr_2

    # ---- 2. Individual parameters ----
    tlag <- exp(ltlag + etaltlag)
    ka   <- exp(lka + etalka)
    vc   <- exp(lvc + etalvc)
    cl   <- exp(lcl + etalcl + iov_cl) * (AST / 25)^e_ast_cl

    rbase <- exp(lrbase + etalrbase)
    kin   <- exp(lkin + etalkin + iov_kin)
    ic50  <- exp(lic50 + etalic50 + iov_ic50)
    kout  <- exp(lkout)
    hill  <- exp(lhill)
    fktr  <- expit(logitfktr + etalogitfktr + iov_fktr)
    ktr   <- fktr * kin

    # ---- 3. Palbociclib PK: one compartment, first-order absorption ----
    kel <- cl / vc
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central
    alag(depot) <- tlag
    Cc <- 1000 * central / vc

    # ---- 4. PK-PD link ----
    # Equation 1 is driven directly by Cc, the plasma palbociclib concentration
    # in ng/mL. See the IC50 entry in ini() for why no nM conversion is applied
    # despite Panetta 2024 labelling IC50 and C_palbociclib "nM".

    # ---- 5. Steady-state-derived constants (Equation 6) ----
    # KM = k_bp / (kin - k_bp) * Ncirc0 = fk_bp / (1 - fk_bp) * Ncirc0 after
    # kin cancels. At the typical values that is
    # 0.559140 / 0.440860 * 267.40 = 339.14 x 10^9/L.
    km_fb <- fktr / (1 - fktr) * rbase

    # ---- 6. Platelet dynamics: Equations 1-5 applied to platelets ----
    # Equation 1 as printed omits the trailing N_PBM factor on the bracketed
    # proliferation term; that it belongs there is forced by the paper's own
    # Equations 6 and 7 (see modellib('Panetta_2024_palbociclib_anc') for the
    # derivation).
    feedback <- km_fb / (km_fb + circ)
    inhib    <- ic50^hill / (ic50^hill + Cc^hill)
    d/dt(precursor1) <- kin * feedback * inhib * precursor1 - ktr * precursor1  # Equation 1, N_PBM
    d/dt(precursor2) <- ktr * (precursor1 - precursor2)                        # Equation 2, N_T1
    d/dt(precursor3) <- ktr * (precursor2 - precursor3)                        # Equation 3, N_T2
    d/dt(precursor4) <- ktr * (precursor3 - precursor4)                        # Equation 4, N_T3
    d/dt(circ)       <- ktr * precursor4 - kout * circ                         # Equation 5, N_circ

    # ---- 7. Drug-free steady-state initial conditions (Equation 7) ----
    precursor1(0) <- kout / ktr * rbase
    precursor2(0) <- kout / ktr * rbase
    precursor3(0) <- kout / ktr * rbase
    precursor4(0) <- kout / ktr * rbase
    circ(0)       <- rbase

    # ---- 8. Observation and error ----
    PLT <- circ
    PLT ~ prop(propSd_PLT)
  })
}
