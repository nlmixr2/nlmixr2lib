Jian_2026_peginterferon_alfa_2b_wbc_et <- function() {
  description <- "Sequential PK/PD model for the white-blood-cell-count response to subcutaneous peginterferon alfa-2b (Pegbing) extrapolated from chronic hepatitis B to adults with essential thrombocythemia (Jian 2026 WBC model, ET extrapolation). The Jian 2025 quasi-equilibrium TMDD population PK model (all parameters fixed) drives a Friberg semi-mechanistic myelosuppression model (proliferating pool, three transit compartments, circulating white blood cells, feedback exponent gamma fixed at 0.16) in which total plasma drug inhibits stem-cell proliferation through an Emax function. All drug-effect and system parameters are the CHB estimates; only the typical baseline count and its IIV are replaced by real-world ET values, and the PK layer is evaluated for healthy subjects."
  reference <- "Jian W, Yin Y, Chen R, Xue J, Gu J, Du Z, He R, Zhou T. Drug Repositioning of Pegbing for Essential Thrombocythemia: Insights from PK/PD Modeling and Extrapolation. Clin Pharmacol Ther. 2026;119(1):219-227. doi:10.1002/cpt.70079 (PK layer: Jian W, et al. Population Pharmacokinetic Model of Pegbing in Healthy Subjects and Chronic Hepatitis B Patients. CPT Pharmacometrics Syst Pharmacol. 2025;14:2014-2025. doi:10.1002/psp4.70104)"
  vignette <- "Jian_2026_peginterferon_alfa_2b_myelosuppression"
  units <- list(
    time = "h",
    dosing = "ug",
    concentration = "ug/L total plasma peginterferon alfa-2b (Cc); white-blood-cell count (circ) in 10^9/L"
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault equation, raw (NOT BSA-normalized), serum creatinine in umol/L: male CrCL = (140 - Age) * Weight / (0.818 * Creatinine); female 0.85 times that (Jian 2026 Table 1 footnotes e and f; identical to Jian 2025).",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Acts only in the inherited Jian 2025 PK layer: linear effect on CL/F centred at 110.38 mL/min, CL = CL_TV * (1 + 0.00504 * (CRCL - 110.38)). No covariate was retained on any WBC parameter (Jian 2026 Table 2 has no covariate row). Jian 2026 reports no creatinine clearance for the ET cohort (Table 1 RWD column is blank) and does not say what value the extrapolation used; the 110.38 mL/min reference makes the CrCL term neutral.",
      source_name = "Creatinine clearance rate"
    ),
    DIS_HEALTHY = list(
      description = "Health-status indicator inherited from the Jian 2025 PK layer: 1 = healthy volunteer, 0 = chronic hepatitis B patient.",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (healthy volunteer)",
      notes = "SET TO 1 FOR EVERY SUBJECT IN THIS MODEL. Jian 2026 Methods 'PK/PD relationship extrapolation': 'it was assumed that ET patients exhibited similar PK characteristics to healthy subjects ... with the covariate of healthy status fixed to healthy'. The column is kept so the inherited PK layer matches Jian_2025_peginterferon_alfa_2b exactly; the CHB fit (Jian_2026_peginterferon_alfa_2b_wbc_chb) sets it to 0.",
      source_name = "Health status (healthy v.s. CHB)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "peginterferon alfa-2b", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "peginterferon alfa-2b", units = "ug", specimen = "plasma", verified = TRUE),
    total_target = list(
      analyte = "type I interferon receptor (IFNAR)",
      units = "ug/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    prol = list(
      analyte = "proliferating marrow stem cells (leukocyte lineage)",
      units = "10^9/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit1 = list(
      analyte = "maturing leukocyte precursors",
      units = "10^9/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit2 = list(
      analyte = "maturing leukocyte precursors",
      units = "10^9/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit3 = list(
      analyte = "maturing leukocyte precursors",
      units = "10^9/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    circ = list(analyte = "white blood cells", units = "10^9/L", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 100L,
    n_studies = 1L,
    sex_female_pct = 56.3,
    race_ethnicity = "Not reported; Chinese real-world cohort",
    disease_state = "Essential thrombocythemia (Philadelphia-negative myeloproliferative neoplasm) treated off-label with Pegbing. Baseline platelet count median 764 x 10^9/L (range 193-2295); WBC median 8.71 x 10^9/L (range 2.95-30.45) (Jian 2026 Table 1, RWD; baseline data available for 71 of the 100 patients).",
    dose_range = "Real-world: mostly Pegbing 180 ug SC once weekly, some 90 or 135 ug depending on disease severity. In silico trial: 90, 180 or 270 ug weekly for 52 weeks.",
    regions = "China",
    parameter_source = "Drug-effect, transit and feedback parameters are the CHB estimates (Phase II NCT01143662, n = 154 CHB patients); only the typical baseline and its IIV come from the ET real-world data. The ET data were used for external VPC validation (Jian 2026 Figure 3b-c) and did not re-estimate any drug-effect parameter.",
    notes = "Sex percentage from the 71 ET patients with baseline data (Jian 2026 Table 1 footnote a). The Discussion states the RWD came from 100 ET patients who received Pegbing and consented; follow-up visits were irregular, with at least one per patient."
  )

  ini({
    # ==================================================================
    # PHARMACOKINETIC LAYER -- carried unchanged from Jian 2025 and FIXED.
    # Jian 2026 Methods 'Population PK/PD modeling': "The population
    # pharmacokinetic (PopPK) model of Pegbing ... has been previously
    # published [ref 22 = Jian 2025]. Exposure of Pegbing was simulated
    # using the PopPK model with individual parameter estimates." The fit
    # was therefore sequential and every PK value below is fixed. Full
    # source-trace in Jian_2025_peginterferon_alfa_2b.R.
    # ==================================================================
    lcl <- fixed(log(0.136)); label("Apparent linear clearance of free drug CL/F at CrCL 110.38 mL/min (L/h)") # Jian 2025 Table 2 final model CL/F = 0.136 L/h
    lvc <- fixed(log(2.95)); label("Apparent central volume of distribution V/F (L)") # Jian 2025 Table 2 final model V/F = 2.95 L
    lka <- fixed(log(0.0101)); label("First-order SC absorption rate constant, healthy-subject reference (1/h)") # Jian 2025 Table 2 final model ka = 0.0101 /h
    lrbase_target <- fixed(log(0.328)); label("Baseline total IFN-receptor level R0 (ug/L, drug-equivalent)") # Jian 2025 Table 2 final model R0 = 0.328 ug/L
    lkint <- fixed(log(0.0827)); label("Internalization rate constant of the drug-receptor complex (1/h)") # Jian 2025 Table 2 final model kint = 0.0827 /h
    lkdeg <- fixed(log(0.544)); label("First-order IFN-receptor degradation rate constant (1/h)") # Jian 2025 Table 2 final model kdeg = 0.544 /h (fixed there too)
    lkd <- fixed(log(0.0493)); label("Quasi-equilibrium dissociation constant KD (ug/L)") # Jian 2025 Table 2 final model KD = 0.0493 ug/L (fixed there too)
    lkdes <- fixed(log(0.0068)); label("Exponential decay rate of IFN-receptor synthesis kR (1/day)") # Jian 2025 Table 2 final model kR = 0.0068 /day; /24 applied in model()
    e_crcl_cl <- fixed(0.00504); label("Linear coefficient of CrCL on CL/F, centred at 110.38 mL/min (per mL/min)") # Jian 2025 Table 2 = 0.00504; linear form per the Jian 2025 Data S2 control stream
    e_chb_ka <- fixed(-0.405); label("Linear fractional effect of CHB disease status on ka (unitless)") # Jian 2025 Table 2 theta CHB on ka = -0.405

    etalcl ~ fixed(0.118) # Jian 2025 Data S2 OMEGA 1 = 0.118 (Table 2 IIV CL/F 34.4%)
    etalvc ~ fixed(0.232) # Jian 2025 Data S2 OMEGA 2 = 0.232 (Table 2 IIV V/F 48.2%)
    etalka ~ fixed(0.197) # Jian 2025 Data S2 OMEGA 3 = 0.197 (Table 2 IIV ka 44.4%)
    etalkint ~ fixed(0.272) # Jian 2025 Data S2 OMEGA 5 = 0.272 (Table 2 IIV kint 52.2%)

    # ==================================================================
    # WHITE-BLOOD-CELL LAYER -- Jian 2026 Table 2, 'WBC model' block.
    # ==================================================================
    # ET EXTRAPOLATION: Jian 2026 Methods 'The typical value and IIV of
    # baseline PLT and WBC were replaced with medians and geometric variances
    # derived from RWD'; Results '... as well as a higher WBC median (8.74 vs.
    # 5.25) and larger geometric variance (0.235 vs. 0.025)'. Every other WBC
    # parameter is the CHB estimate (Table 2).
    lrbase <- log(8.74); label("Baseline white-blood-cell count WBC0, ET real-world median (10^9/L)") # Jian 2026 Results 'Extrapolation of the PK/PD relationship': ET median WBC 8.74 x10^9/L (Table 1 prints the RWD median as 8.71; the Results value used in the extrapolation is taken); replaces the CHB Table 2 WBC0 of 5.25
    lmtt <- log(212); label("Mean transit time MTT = 4 / ktr (h)") # Jian 2026 Table 2 WBC model: MTT = 212 h (RSE 8.1%; bootstrap median 212, 95% CI 188-241)
    lgamma <- fixed(log(0.16)); label("Feedback exponent gamma on (WBC0 / WBC) (unitless)") # Jian 2026 Table 2 WBC model: gamma = 0.16 fixed (footnote a: fixed from literature, refs 25-26)
    lemax <- log(0.0952); label("Maximum fractional inhibition of stem-cell proliferation Emax (unitless)") # Jian 2026 Table 2 WBC model: Emax = 0.0952 (RSE 9.1%; bootstrap median 0.0962, 95% CI 0.0847-0.113)
    lec50 <- log(1940); label("Total Pegbing concentration at half-maximal inhibition EC50 (pg/mL)") # Jian 2026 Table 2 WBC model: EC50 = 1,940 pg/mL (RSE 35.3%; bootstrap median 2,048, 95% CI 1,142-3,521)

    # IIV on EC50: Table 2 footnote prints 'calculated as SQRT(omega2-1)', which is
    # impossible as written. The percentages are 100 * sqrt(omega^2): the
    # Results text gives the CHB 'geometric variance' of baseline WBC as
    # 0.025 = 0.159^2, the PLT block's 0.057 = 0.239^2 discriminates against
    # the log-normal CV reading, and the same authors' Jian 2025 control
    # stream confirms 100 * sqrt(omega^2) for the identical footnote.
    # The WBC block's baseline IIV row is misprinted 'IIV PLT0' in Table 2.
    etalrbase ~ 0.235 # Jian 2026 Results 'Extrapolation of the PK/PD relationship': ET geometric variance of baseline WBC = 0.235 (CHB 0.025 = Table 2 IIV WBC0 15.9% squared, so the geometric variance is omega^2)
    etalec50 ~ 1.292^2 # Jian 2026 Table 2 WBC model: IIV EC50 = 129.2% (RSE 16.5%, shrinkage 27%) -> 1.292^2

    # Residual: combined proportional + additive (Supplementary Method 1,
    # Eq. 10, 'mixed model was finally selected').
    propSd <- 0.193; label("Proportional residual SD on WBC count (fraction)") # Jian 2026 Table 2 WBC model: RUVprop = 19.3% (RSE 9.3%, shrinkage 5%; bootstrap median 19.3, 95% CI 18.2-20.4)
    addSd <- 0.004; label("Additive residual SD on WBC count (10^9/L)") # Jian 2026 Table 2 WBC model: RUVadd = 0.004 x10^9/L (RSE 3.4%, shrinkage 5%; bootstrap median 0.004, 95% CI 0.003-0.004)
  })

  model({
    # ================= Inherited Jian 2025 PK layer ===================
    cl <- exp(lcl + etalcl) * (1 + e_crcl_cl * (CRCL - 110.38))
    vc <- exp(lvc + etalvc)
    chb <- 1 - DIS_HEALTHY
    ka <- exp(lka + etalka) * (1 + e_chb_ka * chb)
    rbase_target <- exp(lrbase_target)
    kint <- exp(lkint + etalkint)
    kdeg <- exp(lkdeg)
    kd <- exp(lkd)
    kdes <- exp(lkdes)

    # Receptor synthesis decays from the start of treatment (Jian 2025 Data
    # S2: KSYN = REC0*KDEG*EXP(-THETA(9)*TIME)); kdes is per day, t in hours.
    ksyn <- kdeg * rbase_target * exp(-kdes * t / 24)

    ctot <- central / vc
    qecore <- ctot - total_target - kd
    cfree <- 0.5 * (qecore + sqrt(qecore * qecore + 4 * kd * ctot))
    bound <- ctot - cfree

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl * cfree - kint * bound * vc
    d/dt(total_target) <- ksyn - kdeg * total_target - (kint - kdeg) * bound
    total_target(0) <- rbase_target

    # Total plasma Pegbing (the ELISA-measured quantity of Jian 2025), ug/L.
    Cc <- ctot

    # ================= WBC layer (Jian 2026 Supplementary Method 1) ======
    rbase <- exp(lrbase + etalrbase)
    mtt <- exp(lmtt)
    gamma <- exp(lgamma)
    emax <- exp(lemax)
    ec50 <- exp(lec50 + etalec50)

    # Supplementary Method 1: kprol = ktr = ke, MTT = 4 / ktr.
    ktr <- 4 / mtt

    # Eq. 6: E = 1 - Emax * C / (EC50 + C). EC50 is in pg/mL and Cc in ug/L
    # (1 ug/L = 1000 pg/mL). The paper names the driver only as 'the
    # concentration of Pegbing'; the total concentration the PopPK model
    # predicts and the ELISA measures is used.
    cpg <- Cc * 1000
    edrug <- 1 - emax * cpg / (ec50 + cpg)

    # Eqs. 1-5 (Eq. 5 prints dPLT/dt on the left for the Circ state).
    d/dt(prol) <- ktr * prol * (rbase / circ)^gamma * edrug - ktr * prol
    d/dt(transit1) <- ktr * (prol - transit1)
    d/dt(transit2) <- ktr * (transit1 - transit2)
    d/dt(transit3) <- ktr * (transit2 - transit3)
    d/dt(circ) <- ktr * transit3 - ktr * circ

    prol(0) <- rbase
    transit1(0) <- rbase
    transit2(0) <- rbase
    transit3(0) <- rbase
    circ(0) <- rbase

    circ ~ add(addSd) + prop(propSd)
  })
}
