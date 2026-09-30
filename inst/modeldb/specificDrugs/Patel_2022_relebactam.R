Patel_2022_relebactam <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for relebactam, a",
    "beta-lactamase inhibitor given with imipenem/cilastatin, in adult healthy",
    "participants and patients with complicated intra-abdominal infection,",
    "complicated urinary tract infection, or hospital-acquired /",
    "ventilator-associated bacterial pneumonia (Patel 2022; 1,197",
    "participants pooled across 12 phase I-III studies, 6,531 quantifiable",
    "relebactam plasma concentrations). An update of the Bhagunde 2019 model",
    "with the phase III RESTORE-IMI 2 (HABP/VABP) and PN017 data. Zero-order",
    "infusion into a central compartment with first-order linear elimination",
    "and linear distribution to one peripheral compartment. Cockcroft-Gault",
    "creatinine clearance enters clearance as a power function centred on",
    "105.5 mL/min, and pneumonia (HABP or VABP) lowers clearance by 43%. Body",
    "weight (centred on 75 kg), pneumonia and mechanical ventilation within",
    "pneumonia enter the central volume. Inter-individual variability is",
    "log-normal on CL, V1 and V2 with an estimated CL-V1 correlation, and",
    "residual error is proportional. The companion imipenem model fitted in",
    "the same NONMEM run is Patel_2022_imipenem."
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

  # Relebactam (MK-7655) is co-administered with imipenem/cilastatin; the
  # assayed analyte is relebactam in plasma.
  #
  # STATE UNITS. The final-model NONMEM control stream embedded in the
  # supplement ('finalmodel.txt', Supplementary Methods 'Final Model') carries
  # AMT in nmol and DV in nmol/L ("AMT = nmol; DV = nmol/L; TIME = h / CL and
  # Q = L/hr; V1 & V2 = L"). The model is linear and CL / Q / V1 / V2 are in
  # L/h and L, so the parameters do not depend on the amount unit: dosing in
  # mg gives Cc in mg/L (= ug/mL). Divide by the relebactam molar mass 348.37
  # g/mol and multiply by 1000 to reproduce the paper's uM exposure metrics.
  compartmentData <- list(
    central = list(analyte = "relebactam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "relebactam", units = "mg", specimen = "plasma", verified = TRUE)
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
        "footnote c 'CL (L/h) = 7.23 x (CrCl/105.5)^0.75 x ...' and the",
        "control stream CL_RLCRCL = ((CRCL/105.5)**THETA(11)). CrCl on",
        "relebactam CL gave the largest OFV change of the analysis (-596.6).",
        "Missing CrCl was imputed with the population median (Supplementary",
        "Methods); the control stream codes the missing value -99 as a",
        "covariate factor of 1. Simulated AUC0-24 fold changes versus normal",
        "renal function (CrCl 90-150 mL/min) were 1.39, 2.05 and 3.35 for",
        "mild, moderate and severe renal impairment (Results, 'Final model')."
      ),
      source_name = "CRCL"
    ),
    WT = list(
      description = "Total body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed baseline value per subject. An estimated power effect on",
        "V1 only (0.65), centred on 75 kg (Table 3 footnote c '(WT / 75)';",
        "control stream TVV1_RL = THETA(6)*((WT/75)**THETA(19))). Unlike",
        "imipenem, relebactam CL carries no weight term: Table 3 prints NA",
        "for 'Covariates on CL / WT (power)' and the control stream has",
        "TVCL_RL = THETA(5). Missing weight was imputed with the population",
        "median (Supplementary Methods)."
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
        "(1 - 0.43) on CL and (1 - 0.29) on V1. Table 2: pneumonia 278 of",
        "1,197 (23.2%); of the 261 PN014 participants with pneumonia, 139",
        "nonventilated HABP, 30 ventilated HABP and 92 VABP."
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
        "remain at baseline values'). Table 3 footnote g: 'The ventilation",
        "effect on participants with pneumonia receiving ventilation. The",
        "reference group was participants with pneumonia.' The control",
        "stream column is VENT2 (V1_RLVENT2 = 1 + THETA(15) when VENT2 = 1).",
        "Encoded in model() as ventilated pneumonia = DIS_VABP + DIS_HABP *",
        "MECH_VENT, so the effect reaches ventilated HABP and all VABP",
        "patients and is never applied outside pneumonia. Effect on V1 only;",
        "ventilation did not affect CL (Results, 'Covariate analysis')."
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
        "Part of the pooled reference group with cIAI and cUTI in the final",
        "model: after reparametrising with healthy participants as reference,",
        "the separate infection-type estimates were imprecise and the model",
        "was simplified 'by grouping participants with cIAI/cUTI with healthy",
        "participants to form the reference group' (Results, 'Final model').",
        "Table 2: 231 healthy participants (19.3%)."
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
    n_observations = 6531,
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
      "Relebactam 25-1,150 mg as an intravenous infusion, single dose or",
      "every 6 hours (Table S1); the clinical regimen is",
      "imipenem/cilastatin/relebactam 500/500/250 mg every 6 hours in normal",
      "renal function, reduced to 200, 150 and 100 mg relebactam for CrCl",
      "60-89, 30-59 and 15-29 mL/min (Table S2)"
    ),
    regions = "Global, including Japanese phase I and phase III studies",
    notes = paste(
      "Table 2 ('Clinical and demographic data for study participants') and",
      "Table S1. The Bhagunde 2019 dataset (855 participants, 10 studies)",
      "plus the phase III studies PN014 (RESTORE-IMI 2, HABP/VABP, n = 261)",
      "and PN017 (Japanese cIAI / cUTI, n = 81). The 1,197 participants",
      "provided 6,100 imipenem and 6,531 relebactam quantifiable",
      "concentrations; BLQ samples (9.4% of relebactam overall) were",
      "excluded, and observations with |CWRES| > 6 were removed as outliers."
    )
  )

  ini({
    # =========================================================================
    # Structural parameters: Table 3 ('Final imipenem and relebactam model
    # parameter estimates'), relebactam NONMEM column, taken at the precision
    # printed in Table 3 footnote c (the body of the table prints 7.23 / 11.2
    # / 6.15 / 10.9). Typical values are for the reference subject: CrCl
    # 105.5 mL/min, 75 kg, healthy / cIAI / cUTI, not ventilated. Structure:
    # 'two-compartment, zero-order intravenous infusion models with
    # first-order linear elimination' (Results, 'Base model'); control stream
    # $SUBROUTINE ADVAN3 TRANS4 with S1 = V1_RL.
    # =========================================================================
    lcl <- log(7.23)
    label("Clearance at the covariate reference (L/h)") # Table 3 footnote c 7.23; Table 3 CL 7.23 (RSE 1.6), 95% CI 7.00-7.47; bootstrap 7.23
    lvc <- log(11.21)
    label("Central volume of distribution at the covariate reference (L)") # Table 3 footnote c 11.21; Table 3 V1 11.2 (RSE 2.7), 95% CI 10.6-11.8; bootstrap 11.2
    lvp <- log(6.15)
    label("Peripheral volume of distribution (L)") # Table 3 footnote c 6.15; Table 3 V2 6.15 (RSE 3.8), 95% CI 5.68-6.62; bootstrap 6.16
    lq <- log(10.93)
    label("Intercompartmental clearance (L/h)") # Table 3 footnote c 10.93; Table 3 Q 10.9 (RSE 7.8), 95% CI 9.22-12.6; bootstrap 10.9

    # =========================================================================
    # Covariate effects, Table 3 relebactam NONMEM column. (CRCL / 105.5) and
    # (WT / 75) enter as power functions; pneumonia and ventilation as
    # proportional shifts (1 + theta * flag), as in the control stream
    # (CL_RLINFC2 = 1 + THETA(12) for INFC2 = 3; V1_RLINFC2 = 1 + THETA(14);
    # V1_RLVENT2 = 1 + THETA(15) for VENT2 = 1). There is no weight effect on
    # relebactam CL (Table 3 'NA'; TVCL_RL = THETA(5)).
    # =========================================================================
    e_crcl_cl <- 0.75
    label("Power exponent of (CRCL / 105.5) on CL (unitless)") # Table 3 'Covariates on CL / CrCl (power)' 0.75 (RSE 4.2), 95% CI 0.68-0.81
    e_habp_vabp_cl <- -0.43
    label("Proportional shift in CL for pneumonia, HABP or VABP (fraction)") # Table 3 'Covariates on CL / Pneumonia' -0.43 (RSE 5.0), 95% CI -0.48 to -0.39
    e_wt_vc <- 0.65
    label("Power exponent of (WT / 75) on V1 (unitless)") # Table 3 'Covariates on V1 / WT (power)' 0.65 (RSE 10.7), 95% CI 0.51-0.79
    e_habp_vabp_vc <- -0.29
    label("Proportional shift in V1 for pneumonia, HABP or VABP (fraction)") # Table 3 'Covariates on V1 / Pneumonia' -0.29 (RSE 16.2), 95% CI -0.38 to -0.19
    e_mech_vent_vc <- 0.36
    label("Proportional shift in V1 for ventilated versus nonventilated pneumonia (fraction)") # Table 3 'Covariates on V1 / Ventilation' 0.36 (RSE 32.0), 95% CI 0.13-0.58; bootstrap 0.36

    # =========================================================================
    # Inter-individual variability. Table 3 footnote h: '%CV = sqrt(omega^2)
    # x 100', so omega^2 = (CV/100)^2 directly (NOT log(CV^2 + 1)). The
    # control stream $OMEGA values confirm the scale: BSV_CL_RL 0.188 ->
    # sqrt 0.434 = 43.4% (Table 3 43.6%). Footnote i: correlation =
    # omega_ij / sqrt(omega_ii * omega_jj), so the covariance is
    # 0.62 x 0.436 x 0.561 = 0.1516495. Q carries no IIV ($OMEGA 0 FIX ;
    # [BSV_Q_RL]).
    # =========================================================================
    etalcl + etalvc ~ c(
      0.1900960, # Table 3 relebactam 'BSV in CL' CV 43.6 (RSE 8.1, shrinkage 13.9): 0.436^2
      0.1516495, # Table 3 relebactam 'Corr CL ~ V1' 0.62 (RSE 12.2): 0.62 x 0.436 x 0.561
      0.3147210 # Table 3 relebactam 'BSV in V1' CV 56.1 (RSE 10.1, shrinkage 18.8): 0.561^2
    )
    etalvp ~ 0.3445690 # Table 3 relebactam 'BSV in V2' CV 58.7 (RSE 26.2, shrinkage 49.5): 0.587^2

    # =========================================================================
    # Residual error: proportional (Results, 'Base model'). The control stream
    # writes F*(1+ERR(3)) + ERR(4) with the additive ERR(4) variance at 0 FIX,
    # and $SIGMA 0.050 -> sqrt 0.224, in line with the Table 3 22.6%.
    # =========================================================================
    propSd <- 0.226
    label("Proportional residual error (fraction)") # Table 3 relebactam 'Residual error, proportional' 22.6 (RSE 6.7, shrinkage 14.0), 95% CI 21.0-24.0
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
    #      CL = 7.23 * (CRCL/105.5)^0.75 * (1 - 0.43*PNEU) * exp(eta5)
    #      V1 = 11.21 * (WT/75)^0.65 * (1 - 0.29*PNEU) * (1 + 0.36*VENT) * exp(eta6)
    #      V2 = 6.15 * exp(eta7);  Q = 10.93
    #    Footnote c prints a stray 'x 0.75' after the CrCl term of the
    #    relebactam CL equation; there is no such factor in the control
    #    stream (CL_RL = THETA(5) * CL_RLCRCL * CL_RLINFC2 * EXP(ETA(5))), and
    #    it would contradict the tabulated typical CL of 7.23 L/h.
    #    Supply CRCL in raw Cockcroft-Gault mL/min and WT in kg.
    # -----------------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (CRCL / 105.5)^e_crcl_cl *
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
    # 5. Observation: total plasma relebactam (S1 = V1_RL). Unbound-drug
    #    PK/PD metrics (fAUC0-24/MIC) need the plasma unbound fraction applied
    #    downstream; Patel 2022 does not print it.
    # -----------------------------------------------------------------------
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
