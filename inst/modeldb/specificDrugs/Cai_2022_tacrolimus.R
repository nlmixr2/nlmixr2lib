Cai_2022_tacrolimus <- function() {
  description <- "Theory-based one-compartment population PK model for oral immediate-release tacrolimus in adult Chinese recipients early after first liver transplantation (Cai 2022 theory-based final model). Parameters are expressed on the PLASMA scale: the measured whole-blood trough is linked to the plasma concentration through saturable binding to erythrocytes (Bmax 418 ug/L, KD 3.8 ug/L, fixed from the literature) scaled by hematocrit. Apparent plasma clearance CL/F depends on AST (power), the patient's own total daily tacrolimus dose (Emax form), recipient CYP3A5 expresser status and triazole antifungal co-therapy; apparent plasma volume V/F rises exponentially with postoperative day. Absorption rate fixed at 4.48 1/h."
  reference <- "Cai XJ, Li RD, Li JH, Tao YF, Zhang QB, Shen CH, Zhang XF, Wang ZX, Jiao Z. Prospective population pharmacokinetic study of tacrolimus in adult recipients early after liver transplantation: A comparison of Michaelis-Menten and theory-based pharmacokinetic models. Front Pharmacol. 2022;13:1031969. doi:10.3389/fphar.2022.1031969"
  vignette <- "Cai_2022_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    AST = list(
      description = "Serum aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. Power effect (AST / 45)^-0.216 on CL/F (Cai 2022 Section 3.2.1 final-model equation; 45 U/L = Table 1 cohort median).",
      source_name = "AST"
    ),
    DOSE_TAC_MGD = list(
      description = "Patient's own total daily tacrolimus dose",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying; the sum of the two q12h doses (Cai 2022 'DD', 'daily tacrolimus dose'; Table 1 mean 2.73 +/- 1.56, median 2.5, range 0.25-8 mg/day). Enters CL/F through the Emax factor 3.51 * DD / (2.44 + DD) (Section 3.2.1 final-model equation), so CL/F rises with dose and approaches 0 as DD approaches 0. Keep this column consistent with the dose records: for q12h dosing DOSE_TAC_MGD is twice the per-administration amt.",
      source_name = "DD"
    ),
    CYP3A5_EXPR = list(
      description = "Recipient CYP3A5 expresser status (1 = *1/*1 or *1/*3; 0 = *3/*3)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP3A5*3/*3 nonexpresser)",
      notes = "Recipient (not donor) genotype at rs776746. CL/F multiplied by 1.168 when CYP3A5_EXPR = 1 (Table 4 theta CYP3A5*1 = 0.168; Eq. 11 scale model). Cohort: 12 *1/*1, 71 *1/*3, 93 *3/*3 (Table 2), i.e. 47.2 percent expressers.",
      source_name = "CYP3A5*1"
    ),
    CONMED_AZOLE = list(
      description = "Concomitant triazole antifungal (voriconazole or fluconazole) co-therapy",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no triazole antifungal)",
      notes = "Time-varying. Cai 2022 'TAF' (triazole antifungal agents; Table 4 row labelled theta VRCZ). CL/F multiplied by 0.425 when CONMED_AZOLE = 1 (Table 4 theta = -0.575; Eq. 11 scale model; Section 3.2.1 'x 0.425 (if combined with TAF)'). Table 1 counts 35 voriconazole and 21 fluconazole samples.",
      source_name = "TAF"
    ),
    POD = list(
      description = "Postoperative day after liver transplantation",
      units = "days",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. V/F = 11000 * exp(0.887 * POD / 15) (Section 3.2.1 final-model equation). The printed normaliser is 15 days, not the Table 1 median of 10 days; the printed equation is followed. Cohort range 2-72 days (Table 1).",
      source_name = "POD"
    ),
    HCT = list(
      description = "Hematocrit",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. Not a covariate on any PK parameter; used only to convert the model's plasma concentration to whole blood through Cai 2022 Eq. 6, Cwb = Cp + Cp * HCT * Bmax / (Cp + KD). The paper states HCT 'as %' but the equation is only coherent with HCT as a volume fraction (with HCT = 31 the whole-blood/plasma ratio would be ~3400 rather than the ~34 that maps the observed 6 ng/mL whole-blood troughs onto the ~0.2 ng/mL plasma scale of Figures 1A and 2A), so the model divides HCT by 100, as Storset_2014_tacrolimus.R does for the same equation. Table 1 median 31.3 percent (range 15.4-53.6).",
      source_name = "HCT"
    )
  )

  covariatesDataExcluded <- list(
    TBILI = list(
      description = "Total serum bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Section 2.5.2) but retained only in the companion Michaelis-Menten model Cai_2022_tacrolimus_mm.R, not in this theory-based model.",
      source_name = "TBIL"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 176L,
    n_observations = 906L,
    n_studies = 1L,
    age_range = "18-74 years",
    age_median = "51 years (mean 50.68, SD 9.71)",
    weight_range = "40-104 kg",
    weight_median = "64 kg (mean 64.54, SD 11.85)",
    sex_female_pct = 14.8,
    race_ethnicity = "Chinese (single centre, Huashan Hospital, Shanghai)",
    disease_state = "Adult recipients of a first liver transplant from donation after cardiac death, early postoperative period (POD 2-72, median 10) until discharge; most recipients hepatitis B, C or E positive.",
    dose_range = "Immediate-release oral tacrolimus (Prograf) q12h, starting 0.5-1 mg q12h and titrated to whole-blood C0 8-12 ng/mL in the first 3 months; daily dose 0.25-8 mg/day (median 2.5).",
    regions = "China",
    notes = "Cai 2022 Table 1 (demographics) and Table 2 (genotypes). 150 male / 26 female. Sparse TDM troughs only (median five per patient) measured by EMIT in whole blood, collected after at least three doses at a constant dose rate. Concomitant steroids in all patients; mycophenolate mofetil when eGFR < 60 mL/min/1.73 m^2."
  )

  ini({
    lka <- fixed(log(4.48)); label("Absorption rate constant ka (1/h)") # Section 2.5.1 and Table 4: Ka fixed at 4.48 1/h from Zhu 2015
    lcl <- log(234); label("Apparent plasma clearance CL/F at AST 45 U/L, CYP3A5 nonexpresser, no triazole antifungal, before the dose factor (L/h)") # Table 4 CLpl/F = 234 L/h (RSE 45.7%)
    lvc <- log(11000); label("Apparent plasma volume V/F at POD = 0 (L)") # Table 4 Vpl/F = 11000 L (RSE 10.8%)

    e_ast_cl <- -0.216; label("Power exponent of (AST / 45 U/L) on CL/F (unitless)") # Table 4 theta AST = -0.216 (RSE 19.8%)
    e_dose_tac_cl_emax <- 3.51; label("Maximum multiplicative factor of the daily-dose Emax term on CL/F (unitless)") # Table 4 theta DDmax = 3.51 (RSE 55.8%)
    e_dose_tac_cl_ed50 <- 2.44; label("Daily dose giving half of the Emax factor on CL/F (mg/day)") # Table 4 theta DD50 = 2.44 (RSE 25%)
    e_cyp3a5_expr_cl <- 0.168; label("Fractional change in CL/F for CYP3A5 expressers (unitless)") # Table 4 theta CYP3A5*1 = 0.168 (RSE 37.1%)
    e_conmed_azole_cl <- -0.575; label("Fractional change in CL/F with triazole antifungal co-therapy (unitless)") # Table 4 theta VRCZ = -0.575 (RSE 17%)
    e_pod_vc <- 0.887; label("Exponential coefficient of POD / 15 on V/F (unitless)") # Table 4 theta POD = 0.887 (RSE 12.3%)

    # Table 4 BSV 'omega CLpl/F (%)' = 29.7 and 'omega Vpl/F (%)' = 59.2, read as
    # CV percent of a log-normal eta: omega^2 = log(CV^2 + 1).
    etalcl ~ 0.0845 # Table 4 omega CLpl/F = 29.7 percent; log(0.297^2 + 1) = 0.0845
    etalvc ~ 0.3004 # Table 4 omega Vpl/F = 59.2 percent; log(0.592^2 + 1) = 0.3004

    propSd <- 0.308; label("Proportional residual error on the plasma concentration (fraction)") # Table 4 delta1 = 30.8 percent (RSE 7.7%)
  })

  model({
    # Erythrocyte-binding constants fixed from the literature (Cai 2022 Section
    # 2.5.2.1, citing Jusko 1995 and Storset 2014a): Bmax in ug/L erythrocytes,
    # KD in ug/L plasma. 1 ug/L = 1 ng/mL.
    bmax_rbc <- 418
    kd_rbc <- 3.8

    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      (AST / 45)^e_ast_cl *
      e_dose_tac_cl_emax * DOSE_TAC_MGD / (e_dose_tac_cl_ed50 + DOSE_TAC_MGD) *
      (1 + e_cyp3a5_expr_cl * CYP3A5_EXPR) *
      (1 + e_conmed_azole_cl * CONMED_AZOLE)
    vc <- exp(lvc + etalvc) * exp(e_pod_vc * POD / 15)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Plasma concentration: mg / L * 1000 = ng/mL. The model was fitted to
    # whole-blood troughs converted to this plasma scale (Figures 1A and 2A).
    Cc <- 1000 * central / vc

    # Whole-blood concentration via saturable erythrocyte binding (Eq. 6);
    # HCT is supplied in percent and converted to a volume fraction here.
    Cwb <- Cc + Cc * (HCT / 100) * bmax_rbc / (Cc + kd_rbc)

    Cc ~ prop(propSd)
  })
}
