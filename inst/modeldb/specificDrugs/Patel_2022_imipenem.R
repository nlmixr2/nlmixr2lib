Patel_2022_imipenem <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for imipenem in adult",
    "healthy participants and patients with complicated intra-abdominal",
    "infection, complicated urinary tract infection, or hospital-acquired /",
    "ventilator-associated bacterial pneumonia (Patel 2022; 1,197",
    "participants pooled across 12 phase I-III studies of the",
    "imipenem/cilastatin/relebactam fixed-dose combination, 6,100",
    "quantifiable imipenem plasma concentrations). An update of the Bhagunde",
    "2019 model with the phase III RESTORE-IMI 2 (HABP/VABP) and PN017 data.",
    "Zero-order infusion into a central compartment with first-order linear",
    "elimination and linear distribution to one peripheral compartment.",
    "Cockcroft-Gault creatinine clearance and body weight enter clearance as",
    "power functions centred on 105.5 mL/min and 75 kg, and pneumonia (HABP",
    "or VABP) lowers clearance by 38%. Body weight, pneumonia and mechanical",
    "ventilation within pneumonia enter the central volume. Inter-individual",
    "variability is log-normal on CL, V1 and V2 with an estimated CL-V1",
    "correlation, and residual error is proportional. The companion",
    "relebactam model fitted in the same NONMEM run is Patel_2022_relebactam."
  )
  reference <- paste(
    "Patel M, Bellanti F, Daryani NM, Noormohamed N, Hilbert DW, Young K,",
    "Kulkarni P, Copalu W, Gheyas F, Rizk ML (2022).",
    "Population pharmacokinetic/pharmacodynamic assessment of",
    "imipenem/cilastatin/relebactam in patients with",
    "hospital-acquired/ventilator-associated bacterial pneumonia.",
    "Clin Transl Sci 15(2):396-408.",
    "doi:10.1111/cts.13158. PMCID PMC8841461.",
    sep = " "
  )
  vignette <- "Patel_2022_imipenem_relebactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Imipenem is given as imipenem/cilastatin with relebactam; cilastatin was
  # not modelled ("Cilastatin concentrations were not considered in the model
  # because cilastatin is primarily excreted renally unchanged and any impact
  # of altered cilastatin levels would be captured by the effect on imipenem
  # PK", Supplementary Methods, 'Data sources').
  #
  # STATE UNITS. The final-model NONMEM control stream embedded in the
  # supplement ('finalmodel.txt', Supplementary Methods 'Final Model') carries
  # AMT in nmol and DV in nmol/L ("AMT = nmol; DV = nmol/L; TIME = h / CL and
  # Q = L/hr; V1 & V2 = L"). The model is linear and CL / Q / V1 / V2 are in
  # L/h and L, so the parameters do not depend on the amount unit: dosing in
  # mg gives Cc in mg/L (= ug/mL). Divide by the imipenem molar mass 299.35
  # g/mol and multiply by 1000 to reproduce the paper's uM exposure metrics.
  compartmentData <- list(
    central = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Creatinine clearance calculated with the Cockcroft-Gault equation,",
        "raw mL/min and NOT normalised to 1.73 m^2 body surface area",
        "(Discussion: 'the Cockcroft-Gault equation was used to derive",
        "CrCl')."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed baseline value per subject; an exploratory analysis found",
        "no need for a time-varying CrCl effect (Results, 'Covariate",
        "analysis'). Power effect on CL centred on 105.5 mL/min: Table 3",
        "footnote c 'CL (L/h) = 12.68 x (CrCl/105.5)^0.48 x ...' and the",
        "control stream CL_IPCRCL = ((CRCL/105.5)**THETA(9)). 105.5 mL/min is",
        "the pooled-cohort median (Table 2 prints 106, range 8-452). Missing",
        "CrCl was imputed with the population median (Supplementary Methods);",
        "the control stream codes the missing value -99 as a covariate factor",
        "of 1. Simulated AUC0-24 fold changes versus normal renal function",
        "(CrCl 90-150 mL/min) were 1.23, 1.59 and 2.18 for mild, moderate and",
        "severe renal impairment (Results, 'Final model')."
      ),
      source_name = "CRCL"
    ),
    WT = list(
      description = "Total body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed baseline value per subject. Estimated power effects on",
        "CL (0.29) and V1 (1.03), both centred on 75 kg (Table 3 footnote c",
        "'(WT/75)'; control stream TVCL_IP = THETA(1)*((WT/75)**THETA(17)),",
        "TVV1_IP = THETA(2)*((WT/75)**THETA(18))). 75 kg is the cohort median",
        "(Table 2, range 27-180 kg). Missing weight was imputed with the",
        "population median (Supplementary Methods). Standard allometric",
        "exponents (0.75 / 1) were tested only in an alternative model",
        "(Table S4) and were not adopted."
      ),
      source_name = "WT"
    ),
    DIS_HABP = list(
      description = paste(
        "Index infection is hospital-acquired bacterial pneumonia (1 = yes,",
        "0 = no), ventilated or nonventilated."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = paste(
        "0 -- the pooled reference group of healthy participants,",
        "participants with cIAI and participants with cUTI (Table 3 footnote",
        "f)"
      ),
      notes = paste(
        "The paper's 'Pneumonia' flag is the control stream's INFC2 = 3",
        "(INFC2 0/1/2 = healthy/cIAI/cUTI share a covariate factor of 1).",
        "It covers HABP and VABP together, so one coefficient is applied to",
        "DIS_HABP + DIS_VABP in model(); the two columns are kept separate",
        "per the register's DIS_VABP entry. The pneumonia effect is",
        "(1 - 0.38) on CL and (1 - 0.39) on V1. cIAI and cUTI were first",
        "estimated as separate infection types and then pooled with healthy",
        "participants because the separate estimates were imprecise (Results,",
        "'Final model'). Table 2: pneumonia 278 of 1,197 (23.2%); of the 261",
        "PN014 participants with pneumonia, 139 nonventilated HABP, 30",
        "ventilated HABP and 92 VABP."
      ),
      source_name = "INFC2"
    ),
    DIS_VABP = list(
      description = paste(
        "Index infection is ventilator-associated bacterial pneumonia (1 =",
        "yes, 0 = no)."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = paste(
        "0 -- the pooled healthy / cIAI / cUTI reference group (Table 3",
        "footnote f)"
      ),
      notes = paste(
        "Shares the single 'Pneumonia' coefficient with DIS_HABP (applied to",
        "DIS_HABP + DIS_VABP). A VABP patient is by definition ventilated, so",
        "model() applies the ventilation effect on V1 to every DIS_VABP = 1",
        "subject regardless of MECH_VENT; see MECH_VENT."
      ),
      source_name = "INFC2"
    ),
    MECH_VENT = list(
      description = paste(
        "Receiving mechanical ventilation (1 = yes, 0 = no), as recorded for",
        "participants with pneumonia in RESTORE-IMI 2 (PN014)."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (nonventilated pneumonia; Table 3 footnote g)",
      notes = paste(
        "Time-fixed at baseline (Discussion: 'Covariates were assumed to",
        "remain at baseline values, which may not be the case for all",
        "covariates (e.g., CrCl and ventilation status)'). Table 3 footnote",
        "g: 'The ventilation effect on participants with pneumonia receiving",
        "ventilation. The reference group was participants with pneumonia.'",
        "The control stream column is VENT2 (V1_IPVENT2 = 1 + THETA(16) when",
        "VENT2 = 1). Encoded in model() as ventilated pneumonia =",
        "DIS_VABP + DIS_HABP * MECH_VENT, so the effect reaches ventilated",
        "HABP and all VABP patients and is never applied outside pneumonia",
        "(ventilation was collected only in the pneumonia study). Ventilated",
        "pneumonia comprised ventilated HABP (30) and VABP (92) of the 261",
        "PN014 pneumonia participants (Table 2). Effect on V1 only; ventilation",
        "did not affect CL (Results, 'Covariate analysis')."
      ),
      source_name = "VENT2"
    )
  )

  # Screened but not retained in the final model. No usable point estimate is
  # published for any of them.
  covariatesDataExcluded <- list(
    DIS_HEALTHY = list(
      description = "Healthy participant (1) versus participant with an infection (0).",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "The Bhagunde 2019 base model carried a health-status effect on",
        "imipenem V1; it was removed before the stepwise covariate search",
        "('the model was modified by removing the HLTH effect (participant",
        "with infection vs. healthy participant) covariate on IPM V1',",
        "Results, 'Base model'). In the final model healthy participants are",
        "part of the pooled reference group with cIAI and cUTI. Table 2:",
        "231 healthy participants (19.3%)."
      ),
      source_name = "HLTH"
    ),
    DIS_CIAI = list(
      description = "Complicated intra-abdominal infection (1 = yes, 0 = no).",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Estimated as a separate infection type during finalisation, then",
        "pooled into the reference group because of imprecision (Results,",
        "'Final model'). Table 2: 308 (25.7%)."
      ),
      source_name = "INFC2"
    ),
    DIS_CUTI = list(
      description = "Complicated urinary tract infection (1 = yes, 0 = no).",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Pooled into the reference group together with cIAI and healthy",
        "participants (Results, 'Final model'). Table 2: 380 (31.7%)."
      ),
      source_name = "INFC2"
    ),
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      notes = paste(
        "Assessed; 'no trend was observed in the covariate-EBE plots'",
        "(Results, 'Covariate analysis'). Age enters the Cockcroft-Gault",
        "CrCl that was retained. Table 2: median 55, range 18-96 years."
      ),
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Sex, female indicator.",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Assessed and not retained (Results, 'Covariate analysis'). The",
        "dataset column is MALE (control stream $INPUT), so SEXF = 1 - MALE.",
        "Table 2: 464 female (38.8%)."
      ),
      source_name = "MALE"
    ),
    RACE_BLACK = list(
      description = "Race, Black / African American indicator.",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Race was significant at forward inclusion but did not meet the",
        "backward-elimination criterion (Results, 'Covariate analysis').",
        "Table 2: 36 (3.0%)."
      ),
      source_name = "RACE"
    ),
    RACE_ASIAN = list(
      description = "Race, Asian indicator.",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "One level of the race covariate that was dropped at backward",
        "elimination; see RACE_BLACK. Table 2: Asian non-Japanese 23 (1.9%),",
        "Japanese 123 (10.3%), Japanese status unknown 18 (1.5%)."
      ),
      source_name = "RACE"
    ),
    RACE_JAPANESE = list(
      description = "Japanese heritage indicator.",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Recorded as a separate dataset column JAP (control stream $INPUT);",
        "part of the race screen dropped at backward elimination. Table 2:",
        "123 Japanese participants (10.3%), largely from PN012, PN019 and",
        "PN017 (Table S1)."
      ),
      source_name = "JAP"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1197,
    n_observations = 6100,
    n_studies = 12,
    age_range = "18-96 years",
    age_median = "55 years",
    weight_range = "27-180 kg",
    weight_median = "75 kg",
    sex_female_pct = 38.8,
    race_ethnicity = c(
      White = 78.8,
      Black = 3.0,
      `Asian (non-Japanese)` = 1.9,
      `Asian (Japanese)` = 10.3,
      `Asian (Japanese status unknown)` = 1.5,
      Other = 4.5
    ),
    disease_state = paste(
      "Pooled healthy participants (231, 19.3%) and adults with complicated",
      "intra-abdominal infection (308, 25.7%), complicated urinary tract",
      "infection (380, 31.7%) or hospital-acquired / ventilator-associated",
      "bacterial pneumonia (278, 23.2%)"
    ),
    renal_function = paste(
      "Cockcroft-Gault CrCl 8-452 mL/min, median 106 mL/min, over 1,194",
      "participants with a recorded value: <15 mL/min 0.4%, 15 to <30 1.9%,",
      "30 to <60 12.0%, 60 to <90 23.2%, 90 to <150 46.6%, 150 to <180",
      "11.0%, 180 to <210 3.0%, 210 to <250 1.0%, >=250 0.9%"
    ),
    dose_range = paste(
      "Imipenem 250-500 mg as a 30-minute intravenous infusion, single dose",
      "or every 6 hours; the clinical regimen is",
      "imipenem/cilastatin/relebactam 500/500/250 mg every 6 hours in normal",
      "renal function, reduced to 400, 300 and 200 mg imipenem for CrCl",
      "60-89, 30-59 and 15-29 mL/min (Table S2)"
    ),
    regions = "Global, including Japanese phase I and phase III studies",
    notes = paste(
      "Table 2 ('Clinical and demographic data for study participants') and",
      "Table S1. The Bhagunde 2019 dataset (855 participants, 10 studies)",
      "plus the phase III studies PN014 (RESTORE-IMI 2, HABP/VABP, n = 261)",
      "and PN017 (Japanese cIAI / cUTI, n = 81). The 1,197 participants",
      "provided 6,100 imipenem and 6,531 relebactam quantifiable",
      "concentrations; BLQ samples (10.2% of imipenem overall) were",
      "excluded, and observations with |CWRES| > 6 were removed as outliers."
    )
  )

  ini({
    # =========================================================================
    # Structural parameters: Table 3 ('Final imipenem and relebactam model
    # parameter estimates'), imipenem NONMEM column, taken at the unrounded
    # precision printed in Table 3 footnote c (the body of the table prints
    # 12.7 / 11.4 / 7.79 / 23.1). Typical values are for the reference
    # subject: CrCl 105.5 mL/min, 75 kg, healthy / cIAI / cUTI, not
    # ventilated. Structure: 'two-compartment, zero-order intravenous infusion
    # models with first-order linear elimination' (Results, 'Base model');
    # control stream $SUBROUTINE ADVAN3 TRANS4 with S1 = V1_IP.
    # =========================================================================
    lcl <- log(12.68)
    label("Clearance at the covariate reference (L/h)") # Table 3 footnote c 12.68; Table 3 CL 12.7 (RSE 1.7), 95% CI 12.3-13.1; bootstrap 12.7
    lvc <- log(11.39)
    label("Central volume of distribution at the covariate reference (L)") # Table 3 footnote c 11.39; Table 3 V1 11.4 (RSE 3.8), 95% CI 10.5-12.3; bootstrap 11.5
    lvp <- log(7.79)
    label("Peripheral volume of distribution (L)") # Table 3 footnote c 7.79; Table 3 V2 7.79 (RSE 5.7), 95% CI 6.90-8.68; bootstrap 7.76
    lq <- log(23.07)
    label("Intercompartmental clearance (L/h)") # Table 3 footnote c 23.07; Table 3 Q 23.1 (RSE 10.9), 95% CI 18.0-28.1; bootstrap 22.9

    # =========================================================================
    # Covariate effects, Table 3 imipenem NONMEM column. Continuous covariates
    # are power functions of (CRCL / 105.5) and (WT / 75). Pneumonia and
    # ventilation are proportional shifts (1 + theta * flag), as in the
    # control stream (CL_IPINFC2 = 1 + THETA(10) for INFC2 = 3; V1_IPINFC2 =
    # 1 + THETA(13); V1_IPVENT2 = 1 + THETA(16) for VENT2 = 1).
    # =========================================================================
    e_crcl_cl <- 0.48
    label("Power exponent of (CRCL / 105.5) on CL (unitless)") # Table 3 'Covariates on CL / CrCl (power)' 0.48 (RSE 4.0), 95% CI 0.44-0.52
    e_wt_cl <- 0.29
    label("Power exponent of (WT / 75) on CL (unitless)") # Table 3 'Covariates on CL / WT (power)' 0.29 (RSE 17.7), 95% CI 0.19-0.39
    e_habp_vabp_cl <- -0.38
    label("Proportional shift in CL for pneumonia, HABP or VABP (fraction)") # Table 3 'Covariates on CL / Pneumonia' -0.38 (RSE 7.7), 95% CI -0.44 to -0.32
    e_wt_vc <- 1.03
    label("Power exponent of (WT / 75) on V1 (unitless)") # Table 3 'Covariates on V1 / WT (power)' 1.03 (RSE 9.7), 95% CI 0.83-1.23
    e_habp_vabp_vc <- -0.39
    label("Proportional shift in V1 for pneumonia, HABP or VABP (fraction)") # Table 3 'Covariates on V1 / Pneumonia' -0.39 (RSE 15.5), 95% CI -0.52 to -0.27
    e_mech_vent_vc <- 0.23
    label("Proportional shift in V1 for ventilated versus nonventilated pneumonia (fraction)") # Table 3 'Covariates on V1 / Ventilation' 0.23 (RSE 46.2), 95% CI 0.02-0.45; bootstrap 0.24; sign per table and control stream, see vignette

    # =========================================================================
    # Inter-individual variability. Table 3 footnote h: '%CV = sqrt(omega^2)
    # x 100', so omega^2 = (CV/100)^2 directly (NOT log(CV^2 + 1)). The
    # control stream $OMEGA values confirm the scale: BSV_CL_IP 0.281 ->
    # sqrt 0.530 = 53.0% (Table 3 53.0%). Footnote i: correlation =
    # omega_ij / sqrt(omega_ii * omega_jj), so the covariance is
    # 0.97 x 0.530 x 0.863 = 0.4436683. Q carries no IIV ($OMEGA 0 FIX ;
    # [BSV_Q_IP]).
    # =========================================================================
    etalcl + etalvc ~ c(
      0.2809000, # Table 3 imipenem 'BSV in CL' CV 53.0 (RSE 9.4, shrinkage 7.2): 0.530^2
      0.4436683, # Table 3 imipenem 'Corr CL ~ V1' 0.97 (RSE 11.5): 0.97 x 0.530 x 0.863
      0.7447690 # Table 3 imipenem 'BSV in V1' CV 86.3 (RSE 14.0, shrinkage 8.7): 0.863^2
    )
    etalvp ~ 0.4032250 # Table 3 imipenem 'BSV in V2' CV 63.5 (RSE 14.7, shrinkage 35.6): 0.635^2

    # =========================================================================
    # Residual error: proportional ('a proportional error model described the
    # residual error', Results, 'Base model'). The control stream writes
    # F*(1+ERR(1)) + ERR(2) with the additive ERR(2) variance at 0 FIX, and
    # $SIGMA 0.086 -> sqrt 0.293, in line with the Table 3 29.5%.
    # =========================================================================
    propSd <- 0.295
    label("Proportional residual error (fraction)") # Table 3 imipenem 'Residual error, proportional' 29.5 (RSE 7.4, shrinkage 14.1), 95% CI 27.2-31.6
  })

  model({
    # -----------------------------------------------------------------------
    # 1. Covariate indicators. The paper's 'Pneumonia' flag (control stream
    #    INFC2 = 3) covers HABP and VABP. The ventilation flag applies only
    #    within pneumonia (Table 3 footnote g): ventilated HABP (MECH_VENT =
    #    1) and every VABP patient, who is ventilated by definition.
    # -----------------------------------------------------------------------
    pneumonia <- DIS_HABP + DIS_VABP
    vent_pneumonia <- DIS_VABP + DIS_HABP * MECH_VENT

    # -----------------------------------------------------------------------
    # 2. Individual disposition parameters. Table 3 footnote c and the
    #    control stream $PK block:
    #      CL = 12.68 * (CRCL/105.5)^0.48 * (WT/75)^0.29 * (1 - 0.38*PNEU) * exp(eta1)
    #      V1 = 11.39 * (WT/75)^1.03 * (1 - 0.39*PNEU) * (1 + 0.23*VENT) * exp(eta2)
    #      V2 = 7.79 * exp(eta3);  Q = 23.07
    #    The typeset footnote c V1 equation carries three typesetting slips:
    #    a stray 'Pneumonia' exponent, '+exp(eta2)' for 'x exp(eta2)', and
    #    'Flag -0.23 Ventilation' with a minus sign. The Table 3 estimate
    #    (+0.23, 95% CI 0.02-0.45), the bootstrap (+0.24, 0.02-0.48) and the
    #    control stream (V1_IPVENT2 = 1 + THETA(16), initial estimate +0.215)
    #    all give a positive ventilation effect, which is what is encoded.
    #    Supply CRCL in raw Cockcroft-Gault mL/min and WT in kg.
    # -----------------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (CRCL / 105.5)^e_crcl_cl * (WT / 75)^e_wt_cl *
      (1 + e_habp_vabp_cl * pneumonia)
    vc <- exp(lvc + etalvc) * (WT / 75)^e_wt_vc *
      (1 + e_habp_vabp_vc * pneumonia) * (1 + e_mech_vent_vc * vent_pneumonia)
    vp <- exp(lvp + etalvp)
    q <- exp(lq)

    # -----------------------------------------------------------------------
    # 3. Micro-rate constants (ADVAN3 TRANS4: CL / V1 / Q / V2).
    # -----------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # -----------------------------------------------------------------------
    # 4. Two-compartment intravenous disposition. The zero-order infusion
    #    rate comes from the event table (30-minute infusion clinically).
    # -----------------------------------------------------------------------
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # -----------------------------------------------------------------------
    # 5. Observation: total plasma imipenem (S1 = V1_IP). Unbound-drug PK/PD
    #    metrics (fT>MIC) need the plasma unbound fraction applied downstream;
    #    Patel 2022 does not print it.
    # -----------------------------------------------------------------------
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
