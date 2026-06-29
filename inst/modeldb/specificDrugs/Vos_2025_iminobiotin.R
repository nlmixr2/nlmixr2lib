Vos_2025_iminobiotin <- function() {
  description <- "Two-compartment intravenous population PK model for the selective neuronal and inducible nitric oxide synthase inhibitor 2-iminobiotin (2-IB) in adults with acute large-vessel-occlusion (LVO) ischemic stroke treated with endovascular thrombectomy (Vos 2025). Central, peripheral, and inter-compartmental clearance are fixed at the upstream TIBOHCA-trial values (Vc = 10.2 L, Q = 15.0 L/h, Vp = 10.4 L); clearance is the only estimated structural parameter. Two typical clearance values are reported, one for patients who did not receive concomitant intravenous thrombolysis (alteplase) (9.29 L/h) and one for patients who did (15.3 L/h, a +65% increase). Baseline estimated glomerular filtration rate enters as an allometric power-form covariate on CL with exponent 0.817 and reference 90 mL/min/1.73 m^2. IIV is estimated on CL only (omega^2 = 0.046, ~22% CV); Vc had no IIV because the volume was fixed. Residual error is a single proportional component with SD 12.3%. The 24-hour continuous infusion uses an eGFR-stratified pump-speed table (Supplemental Table S2) so that all subjects regardless of renal function target the same average exposure (AUC_avg_4h ~ 365 ng*h/mL)."
  reference <- "Vos EM, Peeters-Scholte CMPCD, Boiten J, Hund HM, Jellema K, Kloppenborg RP, et al. Safety, Tolerability, and Pharmacokinetics of the Neuroprotectant 2-Iminobiotin in Patients With Large-Vessel Occlusion Ischemic Stroke Treated With Endovascular Thrombectomy. Stroke. 2025;56(8):1991-1999. doi:10.1161/STROKEAHA.125.050560. PMID:40270284. Structural Vc / Q / Vp inherited from the upstream TIBOHCA out-of-hospital cardiac arrest analysis (van den Heuvel 2024 / earlier Peeters-Scholte work referenced as Vos 2025 reference 11)."
  vignette <- "Vos_2025_iminobiotin"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CRCL = list(
      description        = "Baseline estimated glomerular filtration rate (eGFR) on admission",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject (admission eGFR; clinical-laboratory creatinine-based MDRD estimate per standard Dutch hospital practice). Power-form effect on CL with reference 90 mL/min/1.73 m^2: cl = exp(lcl) * (1 + e_alteplase_cl * CONMED_ALTEPLASE) * (CRCL / 90)^e_crcl_cl. Vos 2025 Supplemental Table S8 reports the covariate as 'COVeGFR0 CL' (baseline eGFR effect on CL) with estimate 0.817 (17.2% RSE). The eGFR reference value of 90 mL/min/1.73 m^2 is the canonical adult reference (sidecar Q4 operator decision; the paper does not state the reference explicitly).",
      source_name        = "eGFR0"
    ),
    CONMED_ALTEPLASE = list(
      description        = "Concomitant intravenous alteplase (recombinant tissue plasminogen activator, r-tPA) coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant alteplase / IVT)",
      notes              = "1 = subject received concomitant intravenous thrombolysis with alteplase before / concurrent with 2-IB infusion; 0 = no IVT. In Vos 2025 the cohort prevalence was 65% (13 of 20 in the 2-IB arm; 11 of 18 in the PK-evaluable subset). The +65% multiplicative increase in clearance (TVCL 9.29 L/h no-IVT vs 15.3 L/h with IVT; Vos 2025 Supplemental Table S8 'CL (L/hr) + IVT') has no firmly established mechanism: the Vos 2025 Discussion ruled out plasmin-mediated cleavage because 2-IB is a small molecule rather than a peptide (unlike nerinetide in the ESCAPE-NA1 trial). Source column 'IVT' in the Vos 2025 NONMEM dataset.",
      source_name        = "IVT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 18L,
    n_studies      = 1L,
    age_range      = "47-89 years (mITT 2-IB arm: median 78, IQR 70-82; overall median 76, IQR 66-82)",
    age_median     = "78 years (2-IB arm; 76 across the mITT pool)",
    weight_range   = "not reported in Vos 2025",
    weight_median  = "not reported in Vos 2025",
    sex_female_pct = 50,
    race_ethnicity = "not reported (single-center Dutch trial; cohort predominantly Northern European)",
    disease_state  = "Acute ischemic stroke due to large-vessel occlusion of the anterior circulation (MCA-M1, proximal MCA-M2, or internal carotid artery), treated with endovascular thrombectomy with or without intravenous alteplase (IVT) within 6 hours of stroke onset",
    dose_range     = "Intravenous loading dose of 3 mL (2.25 mg) of 0.75 mg/mL 2-IB solution over 1 minute (groups A, B intravenously; group C intra-arterially), followed by 1.3 mL/h (0.975 mg/h) intravenous continuous infusion for 4 h (phase 2), followed by an eGFR-adjusted intravenous continuous infusion for 20 h (phase 3) per Vos 2025 Supplemental Table S2 dosing scheme. Pump speeds in phase 3 range from 0.45 mL/h (eGFR 20-29) to 3.5 mL/h (eGFR 200-220), targeting AUC_avg_4h ~ 365 ng*h/mL across the eGFR strata.",
    regions        = "the Netherlands (single-center: Haaglanden Medical Center, The Hague)",
    trial_design   = "Single-center, randomized, double-blind, placebo-controlled phase 2a trial (EudraCT 2021-002162-40 / Dutch trial register 51194). n=20 patients per arm in the mITT population; 18 evaluable for PK in the 2-IB arm (2 missing). Randomization 1:1 across three timing-of-start treatment groups (A: study drug at admission intravenously; B: study drug after reperfusion intravenously; C: study drug after reperfusion intra-arterially).",
    co_medication  = "Concomitant intravenous thrombolysis (alteplase / r-tPA) administered in 65% (13/20) of 2-IB-arm patients per the standard-of-care window; encoded as the CONMED_ALTEPLASE binary covariate.",
    notes          = "Severe renal impairment excluded (eGFR <= 20 mL/min/1.73 m^2 or requiring dialysis). Median NIHSS at admission 19 (IQR 14-21) in the 2-IB arm. Baseline demographics from Vos 2025 Table 1 (mITT) and Supplemental Table S4 (ITT). The cohort age skew (median 78 years, IQR 70-82) reflects the LVO-eligible-for-EVT population which is enriched for elderly patients with atrial fibrillation; ~10-30% prevalence of hypertension / cardiovascular co-morbidities."
  )

  ini({
    # Typical clearance for the IVT-negative reference (no concomitant alteplase),
    # at the eGFR reference of 90 mL/min/1.73 m^2. The IVT effect is applied
    # multiplicatively in model(). Estimated.
    lcl <- log(9.29);    label("Typical clearance CL at CRCL=90 and CONMED_ALTEPLASE=0 (L/h)")  # Vos 2025 Supplemental Table S8: "CL (L/hr) - IVT" = 9.29 (RSE 5.8%)

    # Structural distribution parameters fixed from the upstream TIBOHCA model
    # (Vos 2025 reference 11). Footnote on Supplemental Table S8 reads "10.2 fixed",
    # "15.0 fixed", "10.4 fixed".
    lvc <- fixed(log(10.2)); label("Central volume of distribution Vc (L; fixed from upstream TIBOHCA)")           # Vos 2025 Supplemental Table S8: "Vcentral (L) 10.2 fixed"
    lq  <- fixed(log(15.0)); label("Inter-compartmental clearance Q (L/h; fixed from upstream TIBOHCA)")           # Vos 2025 Supplemental Table S8: "Q (L/hr) 15.0 fixed"
    lvp <- fixed(log(10.4)); label("Peripheral volume of distribution Vp (L; fixed from upstream TIBOHCA)")        # Vos 2025 Supplemental Table S8: "Vperipheral (L) 10.4 fixed"

    # Concomitant intravenous alteplase (IVT) effect on CL: multiplicative
    # increase from TVCL = 9.29 (no-IVT) to TVCL = 15.3 (with IVT), i.e.
    # +(15.3/9.29 - 1) = +0.6469 fractional change in CL.
    e_alteplase_cl <- 0.6469;  label("Fractional change in CL with concomitant alteplase (unitless)")  # Vos 2025 Supplemental Table S8: derived from "CL (L/hr) + IVT" 15.3 vs no-IVT 9.29, i.e. 15.3/9.29 - 1

    # Baseline eGFR allometric power exponent on CL (Vos 2025 Supplemental
    # Table S8: "COVeGFR0 CL"). Reference eGFR = 90 mL/min/1.73 m^2
    # (sidecar Q4 = A; the paper does not state the reference).
    e_crcl_cl <- 0.817;        label("Power-law exponent of (CRCL/90) on CL (unitless)")              # Vos 2025 Supplemental Table S8: "COVeGFR0 CL" = 0.817 (RSE 17.2%)

    # Inter-individual variability on CL only (Vc, Q, Vp are fixed therefore
    # have no IIV; ETA Vcentral reported as a dash in Vos 2025 Supplemental
    # Table S8). omega^2 = 0.046 -> ~21.7% CV. Sidecar Q1 = A; the bolded
    # Random-effect-parameters | 0.225 (0.09) row above ETA1 CL in the
    # same table is a redundant CV / shrinkage summary, not a separate
    # estimated omega.
    etalcl ~ 0.046  # Vos 2025 Supplemental Table S8 row ETA1 CL = 0.046 (RSE 37.5 percent)

    # Proportional residual error: SD = 0.123 directly (12.3% CV).
    # Sidecar Q2 = B confirms 0.123 is the residual SD, not the variance.
    propSd <- 0.123;     label("Proportional residual SD (fraction)")                                # Vos 2025 Supplemental Table S8: "Residual error" = 0.123 (RSE 24.5%); sidecar Q2 = B (SD scale)
  })

  model({
    # Concomitant alteplase fractional multiplier on CL (binary; 1 + 0.6469
    # when CONMED_ALTEPLASE = 1).
    alteplase_eff_cl <- 1 + e_alteplase_cl * CONMED_ALTEPLASE

    # Power-law eGFR effect on CL at reference 90 mL/min/1.73 m^2.
    crcl_eff_cl <- (CRCL / 90)^e_crcl_cl

    # Individual PK parameters.
    cl <- exp(lcl + etalcl) * alteplase_eff_cl * crcl_eff_cl
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # Micro-constants for explicit two-compartment ODEs.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment intravenous disposition. No depot: 2-IB is administered
    # intravenously (loading + continuous infusion) or intra-arterially
    # (group C loading); both inputs enter `central` directly via the event
    # table.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Observation in ng/mL. With dose in mg and Vc in L, central / vc yields
    # mg/L = ug/mL; the 1000 factor converts to ng/mL to match the units the
    # paper reports.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
