Panetta_2024_palbociclib <- function() {
  description <- "One-compartment population PK model with first-order absorption, an absorption lag time and first-order elimination for oral palbociclib in children and young adults (4.9-21.6 years) with recurrent, progressive or refractory brain tumors treated on the Pediatric Brain Tumor Consortium phase 1 study PBTC-042 (NCT02255461). This is the paper's FINAL model, i.e. the AST column of Table 2: apparent oral clearance falls as baseline aspartate aminotransferase rises, through a power model centred on the cohort median AST of 25 U/L with exponent -0.3 (CL/F 48.2 L/h/m^2 at AST = 10 U/L versus 24.1 L/h/m^2 at AST = 100 U/L). Inter-individual variability is diagonal on Tlag, ka, V/F and CL/F, with additional inter-occasion variability on CL/F across the two pharmacokinetic sampling occasions (course 1 day 1 and course 1 day 21). Residual error is proportional. Volumes and clearances are per m^2 of body-surface area and doses are supplied in mg/m^2, so BSA cancels and no BSA covariate is needed. The companion myelosuppression models driven by this PK layer are modellib('Panetta_2024_palbociclib_anc') and modellib('Panetta_2024_palbociclib_plt')."
  reference <- paste(
    "Panetta JC, Selvo NS, Van Mater D, Stewart CF. (2024).",
    "Population Pharmacokinetic and Pharmacodynamic Study of Palbociclib in Children and Young Adults",
    "with Recurrent, Progressive, or Refractory Brain Tumors.",
    "Pharmaceutics 16(12):1528. doi:10.3390/pharmaceutics16121528.",
    "The underlying phase 1 trial is PBTC-042 (NCT02255461).",
    sep = " "
  )
  vignette <- "Panetta_2024_palbociclib"
  units <- list(time = "h", dosing = "mg/m^2", concentration = "ng/mL")
  # Unit note. Panetta 2024 reports V/F in L/m^2 and CL/F in L/h/m^2 (Table 2)
  # and states in section 2.3 that "the palbociclib dosage (i.e., the
  # palbociclib dosage normalized to the patient body surface area (BSA)) was
  # used as the model input". Doses are therefore supplied in mg/m^2, BSA
  # cancels out of the model entirely, and central / vc carries units of
  # mg/L. Plasma concentrations are reported in ng/mL throughout (LLOQ
  # 1 ng/mL, section 2.2; Figure 2 axes), so Cc multiplies by the 1000 ng/mL
  # per mg/L conversion. A clinical dose in mg is converted with
  #   amt [mg/m^2] = dose [mg] / BSA [m^2].

  compartmentData <- list(
    depot   = list(analyte = "palbociclib", units = "mg/m^2", specimen = "administration site", verified = TRUE),
    central = list(analyte = "palbociclib", units = "mg/m^2", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    AST = list(
      description        = "Baseline serum aspartate aminotransferase activity.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Panetta 2024 section 2.3.1: continuous covariates were modelled with a power model scaled to the population median covariate value. The cohort median AST is 25 U/L (range 11-93 U/L, Table 1), which is the centring value. AST was the only covariate retained by the stepwise forward inclusion (Table S1); age, body weight, BSA, gender, serum creatinine, serum albumin, eGFR, total bilirubin and eleven concomitant medications were screened and rejected. The exponent is negative, so clearance falls as AST rises.",
      source_name        = "AST"
    ),
    OCC = list(
      description        = "Pharmacokinetic sampling-occasion index used by the inter-occasion random effect on CL/F.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Two occasions in PBTC-042 (section 2.2): OCC = 1 is the course 1 day 1 single-dose study (pre-dose and 0.5, 1, 2, 4, 8, 24 and 48 h) and OCC = 2 is the course 1 day 21 steady-state study (pre-dose and 1, 2, 4, 8 and 24 h). Panetta 2024 does not print an occasion column name; OCC is the nlmixr2lib canonical. For a single-occasion simulation set OCC = 1 throughout.",
      source_name        = NA_character_
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 33L,
    n_studies      = 1L,
    age_range      = "4.9-21.6 years (median 12.8)",
    weight_range   = "23.8-110.4 kg (median 52.5)",
    height_range   = "109-185.1 cm (median 147.8)",
    bsa_range      = "0.8-2.4 m^2 (median 1.5)",
    sex_female_pct = 36,
    race_ethnicity = c(Caucasian = 58, Black = 12, Other = 30),
    disease_state  = "Recurrent, progressive or refractory pediatric brain tumors. Patients were split into two strata: stratum I (21 patients, 64%) were not heavily pretreated and stratum II (12 patients, 36%) were heavily pretreated.",
    dose_range     = "Oral palbociclib 50 mg/m^2 (6 patients), 75 mg/m^2 (21 patients) or 95 mg/m^2 (6 patients) once daily for 21 days of each 28-day course. The course 1 day 2 dose was held so that a 48 h single-dose sample could be drawn.",
    regions        = "United States (Pediatric Brain Tumor Consortium)",
    observations   = "387 palbociclib plasma concentrations from 33 patients. 34 patients enrolled; one was excluded for taking phenytoin and omeprazole against protocol. 33 patients had course 1 day 1 sampling and 26 of them also had course 1 day 21 sampling. Ten of 397 samples were excluded (four late samples in two patients whose 24 h concentration exceeded their 8 h concentration, and six spurious day-21 samples in one patient).",
    notes          = "Baseline laboratory medians (range), Table 1: albumin 4.2 g/dL (3.2-4.9), AST 25 U/L (11-93), bilirubin 0.4 mg/dL (0.1-1.2), creatinine 0.5 mg/dL (0.2-1.2). Other race comprises Asian (3), Hispanic (4), American Indian / Alaskan (1) and unknown (2). Estimated in Monolix 2023R1 by SAEM with M3 censoring of the below-quantitation-limit samples; precision was assessed by 500-replicate non-parametric bootstrap (Rsmlx 2.0.4)."
  )

  ini({
    # =========================================================================
    # Structural parameters. Values are the AST column of Panetta 2024 Table 2
    # ("Final pharmacokinetic population parameter estimates"), which is the
    # final model -- Figure 2's diagnostics and the bootstrap in section 3.2
    # both refer to the model that includes AST. The Base Model column
    # (Tlag 0.8 h, ka 0.48 /h, V/F 664.5 L/m^2, CL/F 36.8 L/h/m^2) is the
    # covariate-free model quoted in the Abstract and is NOT used here.
    # =========================================================================
    ltlag <- log(0.8);   label("Tlag: absorption lag time (h)")                      # Table 2, AST column: Tlag = 0.8 h (RSE 12.1%); bootstrap mean 0.81 [0.65, 0.97]
    lka   <- log(0.46);  label("ka: first-order absorption rate constant (1/h)")     # Table 2, AST column: ka = 0.46 /h (RSE 17.6%); bootstrap mean 0.51 [0.39, 0.70]
    lvc   <- log(653);   label("V/F: apparent oral volume per BSA (L/m^2)")          # Table 2, AST column: V/F = 653 L/m^2 (RSE 3.5%); bootstrap mean 668.5 [617.1, 757.1]
    lcl   <- log(36.5);  label("CL/F: apparent oral clearance per BSA (L/h/m^2)")    # Table 2, AST column: CL/F = 36.5 L/h/m^2 (RSE 5.2%); bootstrap mean 33.5 [30.5, 37.9]

    # Covariate effect. Power model centred on the cohort median AST of
    # 25 U/L (Table 1); section 2.3.1 states continuous covariates were
    # "modeled using a power model scaled to the population median covariate
    # value". The centring value is confirmed arithmetically by section 3.2,
    # which reports CL/F falling from 48.2 to 24.1 L/h/m^2 between AST = 10
    # and AST = 100 U/L:
    #   36.5 * (10 / 25)^-0.3 = 48.05  and  36.5 * (100 / 25)^-0.3 = 24.08.
    e_ast_cl <- -0.3; label("Exponent of the (AST / 25 U/L) power effect on CL/F (unitless)")  # Table 2, AST column: AST on CL/F = -0.3 (RSE 39.3%); bootstrap mean -0.24 [-0.12, -0.36]; p = 0.0066

    # =========================================================================
    # Inter-individual variability. The Table 2 footnote states "IIV and IOV
    # are reported as a standard deviation", and section 2.3 states the random
    # effects are log-normal, so the tabulated numbers are omega on the log
    # scale and the ini() entries are their squares. The squaring is written
    # out inline so it is auditable.
    # =========================================================================
    etaltlag ~ 0.67^2  # Table 2, AST column: IIV Tlag = 0.67 (RSE 14.1%); bootstrap mean 0.69 [0.54, 0.88]
    etalka   ~ 0.87^2  # Table 2, AST column: IIV ka   = 0.87 (RSE 14.8%); bootstrap mean 0.86 [0.34, 1.19]
    etalvc   ~ 0.04^2  # Table 2, AST column: IIV V/F  = 0.04 (RSE 86%);   bootstrap mean 0.06 [0.02, 0.26]
    etalcl   ~ 0.26^2  # Table 2, AST column: IIV CL/F = 0.26 (RSE 16.3%); bootstrap mean 0.24 [0.15, 0.30]

    # Inter-occasion variability on CL/F across the two sampling occasions.
    # Panetta 2024 reports a single IOV standard deviation, so both occasions
    # share it; the second is fixed to the first (the Monolix analogue of
    # NONMEM's $OMEGA BLOCK(1) SAME), following the two-occasion precedent in
    # Chen_2023_nemonoxacin.R.
    etaiov_cl_1 ~ 0.12^2
    label("Inter-occasion variability in CL/F, occasion 1 (log-scale variance)")  # Table 2, AST column: IOV CL/F = 0.12 (RSE 25.4%); bootstrap mean 0.13 [0.09, 0.18]
    etaiov_cl_2 ~ fix(0.12^2)                                                    # same IOV standard deviation on occasion 2

    # Residual error. Section 2.3: "A proportional residual error model was
    # used"; the tabulated value is the proportional standard deviation.
    propSd <- 0.32; label("Proportional residual error (fraction)")  # Table 2, AST column: residual proportional error = 0.32 (RSE 4.7%); bootstrap mean 0.32 [0.26, 0.39]
  })

  model({
    # ---- 1. Occasion indicators ----
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2

    # ---- 2. Individual parameters ----
    tlag <- exp(ltlag + etaltlag)
    ka   <- exp(lka + etalka)
    vc   <- exp(lvc + etalvc)
    cl   <- exp(lcl + etalcl + iov_cl) * (AST / 25)^e_ast_cl

    # ---- 3. Micro-constants ----
    kel <- cl / vc

    # ---- 4. ODE system: one compartment, first-order absorption ----
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ---- 5. Absorption lag time ----
    alag(depot) <- tlag

    # ---- 6. Observation and error ----
    # central / vc is in mg/L because doses are mg/m^2 and vc is L/m^2;
    # 1 mg/L = 1000 ng/mL, which is the unit Panetta 2024 reports.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
