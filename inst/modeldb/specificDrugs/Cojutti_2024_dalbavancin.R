Cojutti_2024_dalbavancin <- function() {
  description <- paste(
    "Simultaneously fitted two-compartment intravenous population PK model for dalbavancin and",
    "indirect-response (turnover) PD model for C-reactive protein (C-RP) in adults receiving",
    "long-term dalbavancin for documented or suspected staphylococcal osteoarticular infections",
    "(prosthetic joint infection, spondylodiscitis, osteomyelitis, infected pseudoarthrosis,",
    "septic arthritis). Dalbavancin clearance rises exponentially with CKD-EPI estimated",
    "glomerular filtration rate; the effect is UNCENTERED, so exp(lcl) is the non-renal clearance",
    "intercept extrapolated to eGFR = 0 rather than a clearance at any physiological renal",
    "function. Total plasma dalbavancin inhibits C-RP production with FULL inhibition (the",
    "printed equation carries no Imax term, i.e. Imax is structurally 1) and an IC50 of 0.70",
    "mg/L. The C-RP baseline R0 was not estimated: it is fixed to each individual's own",
    "pre-treatment C-RP value, supplied through the CRP covariate column, and the production rate",
    "is derived as kin = R0 * kout.",
    sep = " "
  )
  reference <- paste(
    "Cojutti PG, Tedeschi S, Zamparini E, Viale P, Pea F. Population Pharmacokinetics and",
    "Pharmacodynamics of Dalbavancin and C-Reactive Protein in Patients with Staphylococcal",
    "Osteoarticular Infections. Clin Pharmacokinet. 2024;63(9):1271-1282.",
    "doi:10.1007/s40262-024-01410-2",
    sep = " "
  )
  vignette <- "Cojutti_2024_dalbavancin"

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Checked against Cojutti 2024 ESM Figure S1 (schematic of
  # the PK/PD model) and Equation 1.
  compartmentData <- list(
    central     = list(analyte = "dalbavancin", units = "mg",    specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "dalbavancin", units = "mg",    specimen = "plasma", verified = TRUE),
    # The C-RP turnover state holds a CONCENTRATION, not an amount: Equation 1 is
    # written directly on R = "C-RP concentration in plasma" with no volume term,
    # so kin carries units of mg/dL per h.
    crp         = list(analyte = "C-reactive protein", units = "mg/dL", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "CKD-EPI estimated glomerular filtration rate, BSA-normalized",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column eGFR. Cojutti 2024 Methods 2.1 computed eGFR by three formulas (CKD-EPI,",
        "Cockcroft-Gault, MDRD) and kept 'the formula having the best performance in estimating",
        "dalbavancin clearance'; Results 3.2 reports that the winner was CKD-EPI: 'The only",
        "covariate significantly associated with dalbavancin CL in the base two-compartment model",
        "was eGFR, estimated by means of the CKD-EPI formula.' Methods 2.1 also states eGFR 'was",
        "normalized to 1.73 m^2, in agreement with what we have just done in our previous",
        "dalbavancin population PK model', and Table 1 reports the cohort value in",
        "mL/min/1.73 m^2 (median 93, range 33-144), so the column is stored under the canonical",
        "BSA-normalized CRCL. IMPORTANT - the effect is applied UNCENTERED and EXPONENTIALLY:",
        "Results 3.2 prints the relationship verbatim as CL = 0.030 x e^(0.0042 x eGFR), so",
        "exp(lcl) is the clearance intercept at eGFR = 0 (the paper's non-renal clearance) and",
        "NOT a clearance at any physiological reference value. The arithmetic is closed by the",
        "paper's own numbers: 0.031 x exp(0.0042 x 93) = 0.0458 L/h against the reported total",
        "(non-renal plus renal) clearance of 0.045 L/h from the individual posterior estimates",
        "(Results 3.2), a 1.8% agreement that no centered or power form reproduces. Only CL",
        "carries a covariate; V1, Q and V2 do not. The Monte Carlo simulations in Methods 2.3",
        "span four renal-function classes (eGFR 0-29, 30-59, 60-89 and 90-120 mL/min). The same",
        "uncentered exponential CRCL-on-CL form appears in this group's ceftobiprole model, see",
        "modellib('Cojutti_2023_ceftobiprole')."
      ),
      source_name        = "eGFR"
    ),
    CRP = list(
      description        = "Pre-treatment (time-zero) plasma C-reactive protein concentration, used as the individual C-RP baseline R0",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TIME-FIXED per subject, and NOT a covariate on a PK parameter. Cojutti 2024 Methods 2.2:",
        "'C-RP values at time zero (R0), i.e. when starting dalbavancin treatment, were considered",
        "as the baseline values for modeling[.] As C-RP concentration at baseline is independent",
        "from dalbavancin concentration, R0 was fixed to the value of C-RP for the i-th individual",
        "at time zero, therefore only the inhibition production of C-RP is estimated. R0 is",
        "equivalent to the kin/kout ratio (R0 = kin/kout).' The Table 2 footnote repeats it: 'Only",
        "kout is estimated. kin derives from a parameter transformation (kin = R0 x kout). R0 was",
        "fixed to the C-RP value at time zero of each individual.' Accordingly this column enters",
        "model() twice - as the initial condition crp(0) and inside kin - and R0 appears nowhere",
        "in Table 2 because it was never estimated. UNITS ARE mg/dL, not the register's default",
        "mg/L: Table 1 reports 'C-RP baseline value, mg/dL 2.67 (1.1-30.6)' and the paper's",
        "efficacy target throughout the abstract and Tables 3-4 is C-RP < 1 mg/dL. Supplying this",
        "column in mg/L would rescale every C-RP prediction by 10. SEPARATELY, Methods 2.2 notes",
        "that 'C-RP was also tested as a continuous covariate on dalbavancin CL' - that screening",
        "test was NOT retained (Results 3.2 keeps only eGFR on CL), so no e_crp_cl term exists.",
        "Assay: routine clinical laboratory C-RP measured at baseline and at each TDM instance",
        "(Methods 2.1); the paper does not name the assay platform. The final model was fit to 211",
        "C-RP observations, a median of 4 (range 1-8) subsequent assessments per patient (Table 1)."
      ),
      source_name        = "C-RP"
    )
  )

  # Covariates screened by Cojutti 2024 (Methods 2.2: "The following clinical
  # covariates were then tested on the PK parameters: sex, weight, height, serum
  # creatinine, and eGFR") but NOT retained in the final model. Results 3.2:
  # "The only covariate significantly associated with dalbavancin CL in the base
  # two-compartment model was eGFR". Documentation only - none of these is
  # referenced in model().
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on the PK parameters (Methods 2.2) but not retained. Table 1 reports M/F = 31/14 (68.9%/31.1%)."
    ),
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened on the PK parameters (Methods 2.2) but not retained; the published model applies NO allometric scaling. Table 1 median 78 kg, range 50-110 kg (BMI median 27.4, range 18.8-42.9 kg/m^2)."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened on the PK parameters (Methods 2.2) but not retained. Height was collected (Methods 2.1) but Table 1 reports only the derived BMI, not height itself."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened on the PK parameters (Methods 2.2) but not retained; the renal-function signal entered through the derived CKD-EPI eGFR instead. Collected as laboratory data (Methods 2.1) but not summarised in Table 1."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 45L,
    n_studies      = 1L,
    n_observations = "175 dalbavancin plasma concentrations and 211 C-RP concentrations, fitted simultaneously (Results 3.2)",
    age_range      = "18-80 years",
    age_median     = "61 years",
    weight_range   = "50-110 kg",
    weight_median  = "78 kg",
    bmi_range      = "18.8-42.9 kg/m^2 (median 27.4)",
    sex_female_pct = 31.1,
    race_ethnicity = "Not reported; single-centre Italian cohort.",
    disease_state  = paste(
      "Adults with documented or suspected staphylococcal osteoarticular infection treated with",
      "dalbavancin monotherapy after completing an initial 2-week in-hospital daptomycin-based",
      "combination regimen. Infection types (Table 1): prosthetic joint infection 23 (51.1%; 11",
      "hip, 8 knee, 4 implant-related), spondylodiscitis 8 (17.8%), infected pseudoarthrosis",
      "non-unions 6 (13.3%), osteomyelitis 4 (8.9%), septic arthritis 4 (8.9%). A microbiological",
      "isolate was identified in 40/45 (88.9%; 35 monomicrobial, 5 polymicrobial): 19",
      "methicillin-resistant coagulase-negative staphylococci, 10 MRSA, 10 methicillin-susceptible",
      "CoNS and 6 MSSA. Test of cure was positive in 41/45 (91.1%)."
    ),
    renal_function = "CKD-EPI eGFR median 93 mL/min/1.73 m^2, range 33-144 (Table 1). The Monte Carlo target-attainment analysis extrapolated to four classes: eGFR 0-29, 30-59, 60-89 and 90-120 mL/min.",
    baseline_crp   = "Median 2.67 mg/dL, range 1.1-30.6 mg/dL (Table 1). The Conclusion limits the model's applicability to patients with a baseline C-RP below 30.6 mg/dL.",
    dose_range     = paste(
      "All patients started on two 1500 mg intravenous doses one week apart (days 1 and 8);",
      "further TDM-guided 1500 mg doses were added case by case. Median total dose per treatment",
      "course 3000 mg, range 3000-7500 mg (Table 1). 25/45 (55.6%) received exactly two doses,",
      "16/45 (35.6%) three doses and 4/45 (8.9%) four or more. ESM Figure S6 lists the exact",
      "per-patient schedule; the later doses fall between day 29 and day 111. The paper does not",
      "state the infusion duration."
    ),
    regions        = "Italy (IRCCS Azienda Ospedaliero-Universitaria di Bologna)",
    notes          = paste(
      "Retrospective single-centre study, January 2021 to August 2023 (Ethics Committee",
      "897/2021/Oss/AOUBo). Inclusion required a bone-and-joint infection diagnosis, prior",
      "daptomycin-based combination therapy, and at least two C-RP concentrations of which one at",
      "the start of dalbavancin treatment. Total plasma dalbavancin was measured by a validated",
      "LC-MS/MS method with an LLOQ of 0.5 mg/L; intra- and interday quality-control coefficients",
      "of variation were 9.0-14.0% and 4.8-14.2%. TDM began 21-35 days after treatment start and",
      "was repeated whenever feasible (median 4 assessments per patient, range 1-8). Estimation",
      "was SAEM in Monolix 2023R1, with 1000 non-parametric bootstrap iterations run through the",
      "Rsmlx package. The structural PK model and its initial values were taken from this group's",
      "earlier 69-patient dalbavancin popPK analysis, but every population PK parameter was",
      "RE-ESTIMATED here, so this file has no upstream-model dependency. Model fit: R^2 of observed",
      "versus population- and individual-predicted dalbavancin concentrations 0.85 and 0.94; for",
      "C-RP the population-predicted fit was poor (R^2 = 0.22) and the individual-predicted fit",
      "high (R^2 = 0.94). ESM Figure S1 is the model schematic; Figures S2-S7 are goodness-of-fit,",
      "residual and individual-fit plots and contain no additional parameter values."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural PK parameters -- Cojutti 2024 Table 2, "PK/PD model
    # Value" column (the FINAL simultaneously fitted PK/PD model; there is
    # no separate base-model column). Monolix reports these on the linear
    # scale; they are log-transformed here per the nlmixr2lib convention,
    # which is also the scale Monolix estimated them on ("All individual
    # parameters were considered to be lognormally distributed", Methods
    # 2.2). The paper's V1 / V2 map onto the canonical vc / vp.
    #
    # CL is the INTERCEPT of the uncentered exponential eGFR relationship,
    # i.e. the non-renal clearance extrapolated to eGFR = 0 -- Results 3.2:
    # "the non-renal CL of dalbavancin was 0.031 L/h, while the total CL
    # (non-renal and renal) obtained from the individual posterior
    # estimates was 0.045 L/h". Note that Results 3.2 and the Table 2
    # footnote both print the relationship with a rounded intercept,
    # CL = 0.030 x e^(0.0042 x eGFR); the value used here is Table 2's
    # 0.031, which is the estimate carrying an RSE and a bootstrap
    # interval, and which the same Results paragraph states in words as
    # the non-renal clearance. The two differ by 3% and both reproduce the
    # reported 0.045 L/h at the median eGFR of 93 (0.0458 vs 0.0443).
    # ------------------------------------------------------------------
    lcl  <- log(0.031); label("Non-renal clearance intercept at CRCL = 0 (CL, L/h)")   # Cojutti 2024 Table 2: CL = 0.031 L/h (RSE 15.4%; bootstrap median 0.028, 5th-95th 0.02-0.03)
    lvc  <- log(5.93);  label("Central volume of distribution V1 (L)")                 # Cojutti 2024 Table 2: V1 = 5.93 L (RSE 4.8%; bootstrap median 5.98, 5th-95th 5.43-6.32)
    lq   <- log(0.038); label("Intercompartmental clearance Q (L/h)")                  # Cojutti 2024 Table 2: Q = 0.038 L/h (RSE 34.9%; bootstrap median 0.04, 5th-95th 0.02-0.19)
    lvp  <- log(9.55);  label("Peripheral volume of distribution V2 (L)")              # Cojutti 2024 Table 2: V2 = 9.55 L (RSE 14.2%; bootstrap median 9.08, 5th-95th 6.81-25.59)

    # ------------------------------------------------------------------
    # Covariate effect. Applied as cl = exp(lcl) * exp(e_crcl_cl * CRCL),
    # UNCENTERED, per the verbatim relationship printed in Results 3.2 and
    # repeated in the Table 2 footnote: "CL = 0.030 x e^(0.0042 x eGFR)".
    # Adding eGFR on CL dropped the objective function value by 23.18
    # points and the AIC and BIC by 21.17 and 19.71 (Results 3.2).
    # ------------------------------------------------------------------
    e_crcl_cl <- 0.0042; label("Exponential coefficient of CRCL on CL (per mL/min/1.73 m^2; uncentered)")  # Cojutti 2024 Table 2: beta_eGFR = 0.0042 (RSE 32.8%; bootstrap median 0.005, 5th-95th 0.003-0.007)

    # ------------------------------------------------------------------
    # Structural PD parameters (indirect turnover model with full
    # inhibition of C-RP production, Cojutti 2024 Equation 1). Only kout is
    # estimated; kin is the derived transformation kin = R0 * kout and R0
    # is the individual's own baseline C-RP (the CRP covariate column), so
    # neither appears in Table 2.
    # ------------------------------------------------------------------
    lkout  <- log(0.0037); label("First-order elimination rate constant of plasma C-RP (kout, 1/h)")           # Cojutti 2024 Table 2: kout = 0.0037 1/h (RSE 11.6%; bootstrap median 0.0038, 5th-95th 0.003-0.005)
    lic50  <- log(0.70);   label("Total dalbavancin concentration causing half-maximal inhibition of C-RP production (IC50, mg/L)")  # Cojutti 2024 Table 2: IC50 = 0.70 mg/L (RSE 35.0%; bootstrap median 0.58, 5th-95th 0.35-1.16)

    # ------------------------------------------------------------------
    # Inter-individual variability. Cojutti 2024 Table 2 heads this block
    # "SD of the random effects" and Methods 2.2 states that "All
    # individual parameters were considered to be lognormally distributed,
    # random effects were normally distributed, and an exponential model
    # was used for describing the individual parameter estimates". These
    # are therefore Monolix omegas -- standard deviations of the normal
    # random effect on the log scale -- so the nlmixr2 variance is
    # omega^2 directly. They are NOT %CV values, so the
    # omega^2 = log(CV^2 + 1) conversion used for CV-reporting papers does
    # NOT apply here. (The sibling model modellib('Cojutti_2023_ceftobiprole')
    # DOES report %CV and does carry that conversion; the two papers use
    # different Table conventions and the header wording is what
    # distinguishes them.)
    # ------------------------------------------------------------------
    etalcl   ~ 0.0144  # 0.12^2;  Cojutti 2024 Table 2: omega CL   = 0.12 (bootstrap median 0.11, 5th-95th 0.05-0.18). The printed RSE of 334% is out of scale with every other entry and is very likely a misprint for 3.34% or 33.4%; it does not affect the point estimate encoded here.
    etalvc   ~ 0.0196  # 0.14^2;  Cojutti 2024 Table 2: omega V1   = 0.14 (RSE 42.7%; bootstrap median 0.13, 5th-95th 0.09-0.27)
    etalq    ~ 0.49    # 0.70^2;  Cojutti 2024 Table 2: omega Q    = 0.70 (RSE 39.9%; bootstrap median 0.57, 5th-95th 0.29-1.33)
    etalvp   ~ 0.2209  # 0.47^2;  Cojutti 2024 Table 2: omega V2   = 0.47 (RSE 32.5%; bootstrap median 0.60, 5th-95th 0.39-1.40)
    etalkout ~ 0.3969  # 0.63^2;  Cojutti 2024 Table 2: omega kout = 0.63 (RSE 13.4%; bootstrap median 0.62, 5th-95th 0.46-0.83)
    etalic50 ~ 2.25    # 1.50^2;  Cojutti 2024 Table 2: omega IC50 = 1.5  (RSE 17.9%; bootstrap median 1.38, 5th-95th 1.09-1.83). This is a very large random effect (95% of individual IC50 values span roughly 0.06-8.3 mg/L); Discussion attributes it to the heterogeneity of osteoarticular infection type and disease history.

    # ------------------------------------------------------------------
    # Residual variability. Table 2 lists exactly two entries, b1 and b2,
    # and the table footnote defines them as "b1 and b2 proportional
    # residual errors of the PK and PD models, respectively" -- so both
    # outputs carry a proportional error and neither carries an additive
    # term, even though Methods 2.2 says constant, proportional and
    # combined error models were all tested.
    # ------------------------------------------------------------------
    propSd     <- 0.27; label("Proportional residual error on dalbavancin plasma concentration (fraction)")  # Cojutti 2024 Table 2: b1 = 0.27 (RSE 8.9%; bootstrap median 0.25, 5th-95th 0.17-0.29)
    propSd_crp <- 0.32; label("Proportional residual error on plasma C-RP (fraction)")                        # Cojutti 2024 Table 2: b2 = 0.32 (RSE 6.3%; bootstrap median 0.32, 5th-95th 0.24-0.35)
  })

  model({
    # 1. Individual PK parameters. The eGFR effect on CL is exponential and
    #    UNCENTERED (Results 3.2: CL = 0.030 x e^(0.0042 x eGFR)); V1, Q and
    #    V2 carry no covariate.
    cl <- exp(lcl + etalcl) * exp(e_crcl_cl * CRCL)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    # 2. Individual PD parameters.
    kout <- exp(lkout + etalkout)
    ic50 <- exp(lic50 + etalic50)

    # 3. C-RP baseline. Methods 2.2 and the Table 2 footnote: R0 was NOT
    #    estimated but fixed to the individual's own time-zero C-RP value,
    #    and kin is the derived transformation kin = R0 * kout. R0 is
    #    therefore supplied as data through the CRP covariate column (mg/dL).
    rbase <- CRP
    kin   <- rbase * kout

    # 4. Two-compartment IV disposition micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 5. Total dalbavancin plasma concentration in mg/L (dose in mg, volumes
    #    in L). Defined before the ODE block because it drives the C-RP
    #    production-inhibition term; ESM Figure S1 shows the inhibition arrow
    #    originating from the CENTRAL compartment concentration C(t), and
    #    Equation 1 names the driver "Cp[,] the dalbavancin total plasma
    #    concentration". There is no effect compartment.
    Cc <- central / vc

    # 6. ODE system. Dalbavancin is given intravenously, so doses enter
    #    `central` directly and there is no absorption compartment. The C-RP
    #    equation is Cojutti 2024 Equation 1 verbatim,
    #      dR/dt = kin * (1 - Cp / (IC50 + Cp)) - kout * R,
    #    which has no Imax symbol: inhibition is structurally FULL (Imax = 1),
    #    described in Results 3.2 and the Discussion as an "indirect turnover
    #    Imax model with full inhibition of the C-RP production". The state
    #    starts at the individual baseline, which is the steady state of the
    #    undrugged system because kin = rbase * kout.
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-    k12       * central - k21 * peripheral1
    d/dt(crp)         <-  kin * (1 - Cc / (ic50 + Cc)) - kout * crp
    crp(0)            <-  rbase

    # 7. Observations. Cc is total (not free) plasma dalbavancin -- the paper
    #    fitted the LC-MS/MS total-concentration assay and states all of its
    #    thresholds (8.04 mg/L for fAUC24/MIC target attainment, 14.5 mg/L for
    #    C-RP production inhibition) on the total concentration. crp is plasma
    #    C-RP in mg/dL.
    Cc  ~ prop(propSd)
    crp ~ prop(propSd_crp)
  })
}
