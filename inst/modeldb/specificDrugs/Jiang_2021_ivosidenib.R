Jiang_2021_ivosidenib <- function() {
  description <- paste(
    "Two-compartment population PK model for oral ivosidenib (AG-120) in",
    "adults with IDH1-mutant advanced hematologic malignancies (mostly",
    "relapsed or refractory AML) from the phase 1 AG120-C-001 study",
    "(Jiang 2021). Sequential zero-order release into the depot followed",
    "by first-order absorption; first-order elimination. The model is",
    "parameterised on steady-state apparent parameters at 500 mg once",
    "daily, with a step change between the first dose and repeated",
    "dosing (a 0.50-fold change in relative bioavailability and a",
    "1.66-fold change in clearance, MULTI_DOSE_PT), a less-than-dose-",
    "proportional power effect of dose on relative bioavailability,",
    "baseline-albumin and albumin-ratio power effects on CL/F and Vc/F,",
    "a baseline body-weight power effect on Vc/F, and multiplicative",
    "CYP3A4-inhibitor effects on CL/F (voriconazole, fluconazole,",
    "posaconazole, other moderate/strong and mild inhibitors).",
    "The concentration-QTcF model of the same paper is packaged",
    "separately as Jiang_2021_ivosidenib_QTcF."
  )
  reference <- paste(
    "Jiang X, Wada R, Poland B, Kleijn HJ, Fan B, Liu G, Liu H, Kapsalis S,",
    "Yang H, Le K. Population pharmacokinetic and exposure-response analyses",
    "of ivosidenib in patients with IDH1-mutant advanced hematologic",
    "malignancies. Clin Transl Sci. 2021;14(3):942-953.",
    "doi:10.1111/cts.12959"
  )
  vignette <- "Jiang_2021_ivosidenib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline (time-fixed) weight; power effect on Vc/F only,",
        "(WT / 74.2)^0.92. The 74.2 kg reference is the typical patient",
        "of the Jiang 2021 Figure 2b forest plot ('weight = 74.2 kg'),",
        "which equals the Table S1 overall median weight."
      ),
      source_name = "Wt (baseline body weight)"
    ),
    ALB_BASE = list(
      description = "Per-subject baseline serum albumin (time-fixed)",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters CL/F as (ALB_BASE / 37)^0.82 and Vc/F as",
        "(ALB_BASE / 37)^0.73. Reference 37 g/L is the typical patient of",
        "the Jiang 2021 Figure 2 forest plot ('albumin = 37 g/L') and the",
        "Table S1 overall median. Paired with the time-varying ALB column",
        "through the albumin ratio ALB / ALB_BASE (Wahlby-style",
        "baseline-plus-within-subject decomposition)."
      ),
      source_name = "Baseline albumin"
    ),
    ALB = list(
      description = "Time-varying serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Used only through the albumin ratio ALB / ALB_BASE ('albumin at a",
        "given time divided by baseline albumin (reflects within-patient",
        "variability)', Jiang 2021 Table 1 footnote), which enters CL/F with",
        "exponent 0.99 and Vc/F with exponent 1.1. A user who supplies",
        "ALB = ALB_BASE on every record gets ratio 1 and no within-subject",
        "effect."
      ),
      source_name = "Albumin (time-varying); 'albumin ratio' = albumin / baseline albumin"
    ),
    DOSE = list(
      description = "Ivosidenib dose per administration",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters relative bioavailability as (DOSE / 500)^-0.49 (Jiang 2021",
        "Table 1 'Dose-Frel exponent'); the 500 mg reference is the",
        "approved once-daily dose, confirmed by the Figure 2c forest-plot",
        "arithmetic (300 mg -> 71 and 800 mg -> 118 ug*h/mL from 93).",
        "The paper does not state whether the 100 mg twice-daily cohort",
        "(n = 4) was coded as 100 mg or 200 mg; this model uses the amount",
        "per administration. Supply DOSE on every record (dose and",
        "observation rows) and place the column after amt in the event",
        "table."
      ),
      source_name = "Dose"
    ),
    MULTI_DOSE_PT = list(
      description = "Repeated-dosing indicator: 0 = from the first ivosidenib dose until the second dose; 1 = from the second dose onward",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (repeated dosing; the steady-state apparent parameters are the model's reference)",
      notes = paste(
        "Step change between first-dose and steady-state PK (Jiang 2021",
        "Figure 1a: 'Frel(day 1), Frel(ss)' and 'CL(day 1), CL(ss)').",
        "At MULTI_DOSE_PT = 0 the dose's relative bioavailability is",
        "1/0.50 = 2-fold and clearance is 1/1.66-fold the steady-state",
        "value, so first-dose CL/F = 5.39 / (2 * 1.66) = 1.63 L/h as",
        "printed. Supplementary Figure S1c states that the factors 'begin",
        "to influence the pharmacokinetic curve' at the start of",
        "continuous once-daily dosing, and the Figure 1b VPC groups the",
        "escalation day -3 single dose with the expansion cycle 1 day 1",
        "dose as the 'first dose'. Code the first dose record and every",
        "record before the second dose as 0, and the second dose and",
        "every later record as 1. The paper estimates no onset time",
        "course (auto-induction is described as complete within about 2",
        "weeks, but the model is a step)."
      ),
      source_name = "single-dose to steady-state factor (no column name printed)"
    ),
    CONMED_VORICONAZOLE = list(
      description = "Concomitant voriconazole indicator (1 = coadministered)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no voriconazole)",
      notes = "Multiplicative 0.64-fold change in CL/F (Jiang 2021 Table 1). Time-varying per concentration record (Table S2 counts records, not patients).",
      source_name = "Voriconazole"
    ),
    CONMED_FLUCONAZOLE = list(
      description = "Concomitant fluconazole indicator (1 = coadministered)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no fluconazole)",
      notes = "Multiplicative 0.59-fold change in CL/F (Jiang 2021 Table 1). Time-varying per concentration record.",
      source_name = "Fluconazole"
    ),
    CONMED_POSACONAZOLE = list(
      description = "Concomitant posaconazole indicator (1 = coadministered)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no posaconazole)",
      notes = "Multiplicative 0.65-fold change in CL/F (Jiang 2021 Table 1). Time-varying per concentration record.",
      source_name = "Posaconazole"
    ),
    CONMED_CYP3A4_INH_STRONG = list(
      description = "Concomitant strong CYP3A4 inhibitor indicator (1 = coadministered)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no strong CYP3A4 inhibitor)",
      notes = paste(
        "Jiang 2021 pooled every moderate or strong CYP3A4 inhibitor other",
        "than voriconazole, fluconazole and posaconazole into one 'other",
        "strong or moderate CYP3A4 inhibitors' group with a 0.92-fold CL/F",
        "effect (Table 1). This model forms that group as",
        "max(CONMED_CYP3A4_INH_MOD, CONMED_CYP3A4_INH_STRONG) and switches",
        "it off when any of the three named azoles is flagged, so a record",
        "coded to the register definition (voriconazole and posaconazole",
        "are strong inhibitors) is not double-counted. The paper does not",
        "name the agents in the 'other' group."
      ),
      source_name = "Other moderate/strong CYP3A4 inhibitors"
    ),
    CONMED_CYP3A4_INH_MOD = list(
      description = "Concomitant moderate CYP3A4 inhibitor indicator (1 = coadministered)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no moderate CYP3A4 inhibitor)",
      notes = paste(
        "Pooled with CONMED_CYP3A4_INH_STRONG into the paper's 'other",
        "strong or moderate CYP3A4 inhibitors' group (0.92-fold CL/F);",
        "switched off when fluconazole, voriconazole or posaconazole is",
        "flagged, since those carry their own coefficients. See the",
        "CONMED_CYP3A4_INH_STRONG notes."
      ),
      source_name = "Other moderate/strong CYP3A4 inhibitors"
    ),
    CONMED_CYP3A4_INH_WEAK = list(
      description = "Concomitant mild (weak) CYP3A4 inhibitor indicator (1 = coadministered)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no mild CYP3A4 inhibitor)",
      notes = "Jiang 2021 'mild CYP3A inhibitors' group (Table S2 'Other mild CYP3A4 inhibitors'); multiplicative 1.04-fold change in CL/F (Table 1). The paper does not name the agents.",
      source_name = "Mild CYP3A4 inhibitors"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened; 'uncorrelated to ivosidenib CL/F' (Jiang 2021 Results 'Final population PK model')."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not retained (Jiang 2021 Results)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened; not retained (Jiang 2021 Results)."
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units = "mL/min",
      type = "continuous",
      notes = "Screened (continuous and renal-impairment category); no effect on CL/F. Only two patients had severe renal impairment."
    ),
    HEPIMP = list(
      description = "Hepatic impairment category (NCI ODWG)",
      units = "(categorical)",
      type = "categorical",
      notes = "Screened together with ALT, AST and bilirubin; no correlation with CL/F in patients with mild hepatic impairment."
    ),
    CONMED_PPI = list(
      description = "Concomitant proton-pump inhibitor indicator",
      units = "(binary)",
      type = "binary",
      notes = "Pantoprazole, the most frequently used proton-pump inhibitor, was screened on its own and 'did not affect ivosidenib CL/F' (Jiang 2021 Results)."
    ),
    CONMED_H2RA = list(
      description = "Concomitant H2-receptor antagonist indicator",
      units = "(binary)",
      type = "binary",
      notes = "Famotidine, the most frequently used H2-receptor antagonist, was screened on its own and 'did not affect ivosidenib CL/F' (Jiang 2021 Results)."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "ivosidenib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ivosidenib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ivosidenib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 253L,
    n_studies = 1L,
    age_range = "18-89 years",
    age_median = "68 years",
    weight_range = "37.7-150.4 kg",
    weight_median = "74.2 kg",
    sex_female_pct = 46,
    race_ethnicity = c(White = 69, Black = 6, Asian = 2, Other = 4, Missing = 20),
    disease_state = "IDH1-mutant advanced hematologic malignancies: relapsed or refractory AML (80%), untreated AML (14%), other (6%)",
    dose_range = "100 mg twice daily, or 300, 500, 800 or 1200 mg once daily orally in continuous 28-day cycles (225 of 255 patients at 500 mg once daily); escalation patients also received a single dose on day -3",
    regions = "Multinational phase 1 (AG120-C-001, NCT02074839)",
    notes = paste(
      "4656 plasma concentrations from 253 patients (Jiang 2021 Results",
      "'Baseline patient characteristics'). Demographics from Table S1",
      "(N = 255 before exclusion of two patients with missing",
      "administration times): 54% male, ECOG PS 0-1 in 78%, baseline",
      "albumin median 37 g/L (15-48), CrCl median 83 mL/min/1.73 m^2."
    )
  )

  ini({
    # Steady-state apparent parameters at the reference patient (500 mg
    # once daily, repeated dosing, ALB_BASE = ALB = 37 g/L, WT = 74.2 kg,
    # no CYP3A4 inhibitor). Jiang 2021 Table 1.
    lcl <- log(5.39); label("Steady-state apparent clearance CL/F (L/h)") # Table 1 'Steady-state CL/F, L/h' = 5.39 (RSE 4%)
    lvc <- log(234); label("Steady-state apparent central volume Vc/F (L)") # Table 1 'Steady-state Vc/F, L' = 234 (RSE 7%)
    lq <- log(15.8); label("Steady-state apparent intercompartmental clearance Q/F (L/h)") # Table 1 'Steady-state Q/F, L/h' = 15.8 (RSE 19%)
    lvp <- log(151); label("Steady-state apparent peripheral volume Vp/F (L)") # Table 1 'Steady-state Vp/F, L' = 151 (RSE 22%)
    lka <- log(1.38); label("First-order absorption rate constant ka (1/h)") # Table 1 'ka, 1/h' = 1.38 (RSE 10%)
    ld1 <- log(0.27); label("Duration of zero-order release into the depot (h)") # Table 1 'Tlag, h' = 0.27 (RSE 11%), defined in the footnote as 'zero-order release duration (lag-time)'

    # First-dose versus repeated-dosing step (MULTI_DOSE_PT)
    e_multi_dose_pt_f <- 0.50; label("Fold change in relative bioavailability at steady state versus the first dose (unitless)") # Table 1 'Steady-state fold change in Frel' = 0.50 (RSE 7%)
    e_md_cl <- 1.66; label("Fold change in clearance at steady state versus the first dose (unitless)") # Table 1 'Steady-state fold change in CL' = 1.66 (RSE 11%)

    # Dose nonlinearity on relative bioavailability
    e_dose_fdepot <- -0.49; label("Power exponent of DOSE/500 mg on relative bioavailability (unitless)") # Table 1 'Dose-Frel exponent' = -0.49 (RSE 19%)

    # Body weight and albumin
    e_wt_vc <- 0.92; label("Power exponent of WT/74.2 kg on Vc/F (unitless)") # Table 1 'Wt-Vc/F exponent' = 0.92 (RSE 13%)
    e_alb_base_cl <- 0.82; label("Power exponent of ALB_BASE/37 g/L on CL/F (unitless)") # Table 1 'Baseline albumin-CL/F exponent' = 0.82 (RSE 20%)
    e_alb_ratio_cl <- 0.99; label("Power exponent of the albumin ratio ALB/ALB_BASE on CL/F (unitless)") # Table 1 'Albumin ratio-CL/F exponent' = 0.99 (RSE 19%)
    e_alb_base_vc <- 0.73; label("Power exponent of ALB_BASE/37 g/L on Vc/F (unitless)") # Table 1 'Baseline albumin-Vc/F exponent' = 0.73 (RSE 28%)
    e_alb_ratio_vc <- 1.1; label("Power exponent of the albumin ratio ALB/ALB_BASE on Vc/F (unitless)") # Table 1 'Albumin ratio-Vc/F exponent' = 1.1 (RSE 38%)

    # CYP3A4 inhibitors: multiplicative fold changes on CL/F
    e_conmed_voriconazole_cl <- 0.64; label("Fold change in CL/F with voriconazole (unitless)") # Table 1 'Fold change in CL with voriconazole' = 0.64 (RSE 6%)
    e_conmed_fluconazole_cl <- 0.59; label("Fold change in CL/F with fluconazole (unitless)") # Table 1 'Fold change in CL with fluconazole' = 0.59 (RSE 6%)
    e_conmed_posaconazole_cl <- 0.65; label("Fold change in CL/F with posaconazole (unitless)") # Table 1 'Fold change in CL with posaconazole' = 0.65 (RSE 12%)
    e_conmed_cyp3a4_inh_modstrong_cl <- 0.92; label("Fold change in CL/F with other moderate or strong CYP3A4 inhibitors (unitless)") # Table 1 'Fold change in CL with other moderate/strong CYP3A inhibitors' = 0.92 (RSE 17%)
    e_conmed_cyp3a4_inh_weak_cl <- 1.04; label("Fold change in CL/F with mild CYP3A4 inhibitors (unitless)") # Table 1 'Fold change in CL with mild CYP3A inhibitors' = 1.04 (RSE 6%)

    # IIV: Table 1 'CV%' taken as omega * 100 (the footnote reports RSE
    # on these standard-deviation terms as RSE of variance / 2), so
    # omega^2 = (CV%/100)^2. No off-diagonal elements are reported.
    etalcl ~ 0.1225 # Table 1 CL/F BSV CV% = 35 -> 0.35^2
    etalvc ~ 0.2209 # Table 1 Vc/F BSV CV% = 47 -> 0.47^2
    etalka ~ 1.1664 # Table 1 ka BSV CV% = 108 -> 1.08^2

    # Residual error: log-additive (additive on log-transformed data)
    expSd <- 0.26; label("Log-additive residual error SD (log scale)") # Table 1 'Log-additive CV%' = 26 (RSE 3%)
  })

  model({
    # 1. Derived covariate terms
    # The paper's 'other strong or moderate CYP3A4 inhibitors' group
    # excludes the three azoles, which carry their own coefficients.
    azole <- max(CONMED_VORICONAZOLE, max(CONMED_FLUCONAZOLE, CONMED_POSACONAZOLE))
    inh_other <- max(CONMED_CYP3A4_INH_MOD, CONMED_CYP3A4_INH_STRONG) * (1 - azole)
    cyp3a4_cl <- e_conmed_voriconazole_cl^CONMED_VORICONAZOLE *
      e_conmed_fluconazole_cl^CONMED_FLUCONAZOLE *
      e_conmed_posaconazole_cl^CONMED_POSACONAZOLE *
      e_conmed_cyp3a4_inh_modstrong_cl^inh_other *
      e_conmed_cyp3a4_inh_weak_cl^CONMED_CYP3A4_INH_WEAK

    alb_ratio <- ALB / ALB_BASE

    # First-dose step: the reference is repeated dosing (MULTI_DOSE_PT = 1);
    # before the second dose CL is 1/1.66-fold and the dose's relative
    # bioavailability is 1/0.50-fold the steady-state value.
    md_cl <- e_md_cl^(MULTI_DOSE_PT - 1)
    md_f <- e_multi_dose_pt_f^(MULTI_DOSE_PT - 1)

    # 2. Individual parameters
    cl <- exp(lcl + etalcl) * md_cl * (ALB_BASE / 37)^e_alb_base_cl *
      alb_ratio^e_alb_ratio_cl * cyp3a4_cl
    vc <- exp(lvc + etalvc) * (WT / 74.2)^e_wt_vc *
      (ALB_BASE / 37)^e_alb_base_vc * alb_ratio^e_alb_ratio_vc
    q <- exp(lq)
    vp <- exp(lvp)
    ka <- exp(lka + etalka)
    d1 <- exp(ld1)

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODEs
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Bioavailability and zero-order release. Dose records must carry
    # rate = -2 for rxode2 to use the modelled duration d1.
    f(depot) <- (DOSE / 500)^e_dose_fdepot * md_f
    dur(depot) <- d1

    # 6. Observation (dose mg / volume L = mg/L; x 1000 -> ng/mL)
    Cc <- central / vc * 1000
    Cc ~ lnorm(expSd)
  })
}
