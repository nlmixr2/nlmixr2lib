Sikma_2020_tacrolimus_unbound_plasma <- function() {
  description <- "Two-compartment population PK model for whole-blood (Cc), unbound plasma (Cupc), and total plasma (Ctpc) tacrolimus in 30 adult thoracic-organ (10 heart + 20 lung) transplant recipients during the first 6 postoperative days (Sikma 2020). First-order oral absorption with ka, F, and the within-PK fixed-parameter variabilities inherited from a previously estimated tacrolimus model; non-linear saturable binding of tacrolimus to erythrocytes (UPC = WBC * Kd / (Bmax * HCT - WBC)) with the maximum erythrocyte binding capacity Bmax scaled by hematocrit, and a linear non-specific plasma binding constant Nplasma linking unbound to total plasma (TPC = Nplasma * UPC)."
  reference <- "Sikma MA, Van Maarseveen EM, Hunault CC, Moreno JM, Van de Graaf EA, Kirkels JH, Verhaar MC, Grutters JC, Kesecioglu J, De Lange DW, Huitema ADR. Unbound Plasma, Total Plasma, and Whole-Blood Tacrolimus Pharmacokinetics Early After Thoracic Organ Transplantation. Clin Pharmacokinet. 2020;59(6):771-780. doi:10.1007/s40262-019-00854-1"
  vignette <- "Sikma_2020_tacrolimus_unbound_plasma"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL (whole blood) and pg/mL (unbound and total plasma)")

  covariateData <- list(
    HCT = list(
      description        = "Hematocrit, expressed as a fraction of total blood volume",
      units              = "fraction",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying covariate. Last-observation-carried-forward per Sikma 2020 Methods (Mixed-Effects Modeling): hematocrit was introduced into the model by multiplying Bmax (and Nplasma in a sensitivity branch) with the observed hematocrit. Study-population medians by day were 0.31, 0.28, 0.27, 0.27, 0.27, 0.28 across days 1-6 post-transplant (Sikma 2020 Table 1). Source paper reports HCT as a fraction (0-1) directly in the binding equation, not as percent; the canonical-register HCT entry's units (%) are explicitly overridden here so that the published values for Bmax (2700) and Kd (0.142) reproduce the paper's UPC equation directly. To use a dataset that records HCT in percent, multiply the column by 0.01 before passing it to this model.",
      source_name        = "Ht"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 30L,
    n_studies       = 1L,
    n_observations  = "1180 tacrolimus concentrations across whole-blood, total plasma, and unbound plasma matrices; 119 whole-blood 0-12 h profiles (median 5 per patient, range 1-6) and 90 total-and-unbound plasma profiles (median 3 per patient, range 0-6); 46 of 1180 (3.9%) below the lower limit of quantification were discarded.",
    age_range       = "34-60 years (Q1-Q3 of 30 patients)",
    age_median      = "43 years",
    weight_range    = "61-86 kg (Q1-Q3)",
    weight_median   = "73.5 kg",
    height_median   = "173.5 cm",
    sex_female_pct  = 50.0,
    race_ethnicity  = "Not reported in source paper (single-centre Utrecht, Netherlands cohort).",
    disease_state   = "Adult thoracic organ transplant recipients (10 heart, 20 lung; 18 of 20 lung were double-lung transplantations) admitted to the intensive care unit during the first 6 postoperative days. Heart indications: 5 ischemic and 5 non-ischemic dilated cardiomyopathy. Lung indications: cystic fibrosis (10), chronic obstructive pulmonary disease (3), idiopathic pulmonary arterial hypertension (2), idiopathic pulmonary fibrosis (2), bronchiectasis (1), Langerhans cell histiocytosis (1), sarcoidosis (1).",
    dose_range      = "Oral tacrolimus (Prograft) twice daily, starting at 0.1 mg/kg/dose for lung recipients and 2 mg/dose for heart recipients on the day of transplantation; thereafter titrated to a target whole-blood trough of 9-15 ng/mL (12-hour post-dose, 6 a.m.).",
    regions         = "Netherlands (University Medical Center Utrecht).",
    ecmo_frequency  = "8 of 30 patients (27%) received postoperative ECMO with a median duration of 4 days (Q1-Q3 2-6).",
    baseline_labs   = "Median (Q1-Q3) day-1 hematocrit 0.31 (0.28-0.35); day-1 albumin 26.2 g/L (22.5-29.3); day-1 HDL 0.84 mmol/L (0.70-1.06); day-1 alpha-1-acid glycoprotein 0.89 g/L (0.76-1.18); day-1 pH 7.39 (7.33-7.43); day-1 packed-red-blood-cell transfusion volume 275 mL (275-550). Hematocrit fell from a day-1 median of 0.31 to 0.27-0.28 by days 3-6 (Sikma 2020 Table 1).",
    immunosuppression = "Triple therapy with tacrolimus + a cell-cycle blocker + an interleukin-2 inhibitor + corticosteroids per institutional protocol.",
    notes           = "All thoracic organ recipients admitted to the Utrecht ICU between June 2013 and March 2015 (NTR 3912 / EudraCT 2012-001909-24). Blood samples taken at 0, 2 (or 3 in cystic fibrosis), 6, and 12 hours after administration; collected between 6 p.m. and 6 a.m. Unbound tacrolimus quantified by LC-MS/MS (Stienstra method) with assay range 1.00-200 pg/mL (LLOQ 1 pg/mL in ultrafiltrate, 100 pg/mL in plasma); whole-blood tacrolimus by HPLC-MS/MS with LLOQ 0.5 ng/mL (assay range 1-50 ng/mL)."
  )

  ini({
    # Final-model parameter estimates from Sikma 2020 Table 3. NONMEM 7.3.0
    # fit by SAEM with importance-sampling likelihood and sampling-importance-
    # resampling (SIR) precision; reported 95% CIs are SIR-based.
    #
    # Absorption parameters (Ka, IIVs on Ka / V1 / Q / V2 / F, IOV on Ka, F)
    # were fixed at values inherited from a prior tacrolimus model in this
    # programme; the lead paper does not list the upstream citation. Ka and F
    # are encoded with `fixed()`; the "10% Fixed" IIVs on V1, Ka, Q, V2, F
    # and the "3% Fixed" IIV on Kd are kept as `~ fixed(...)` so the
    # parameter file preserves the source's structural choices. Inter-occasion
    # variability on Ka (98.3% Fixed) and F (65%) is documented in vignette
    # Assumptions and deviations but is not encoded in the model file because
    # it requires per-occasion data structure that nlmixr2lib does not
    # standardise for popPK extractions.

    # ---- Structural disposition (Sikma 2020 Table 3) ------------------------
    lcl  <- log(20.9)         ; label("Whole-blood apparent clearance CL/F (L/h)")                              # Sikma 2020 Table 3 CL = 20.9 (95% CI 16.8-24.7)
    lvc  <- log(220)          ; label("Whole-blood apparent central volume V1/F (L)")                           # Sikma 2020 Table 3 V1 = 220 (95% CI 187-246)
    lq   <- log(72.0)         ; label("Whole-blood apparent inter-compartmental clearance Q/F (L/h)")           # Sikma 2020 Table 3 Q = 72.0 (CI as printed 529-767 is an apparent typo; point estimate 72.0 is used)
    lvp  <- log(469)          ; label("Whole-blood apparent peripheral volume V2/F (L)")                        # Sikma 2020 Table 3 V2 = 469 (95% CI 399-579)

    # ---- Absorption / bioavailability (fixed from upstream model) -----------
    lka     <- fixed(log(0.579)) ; label("First-order absorption rate constant Ka (1/h; FIXED from upstream tacrolimus model)") # Sikma 2020 Table 3 Ka = 0.579 Fixed
    lfdepot <- fixed(log(1))     ; label("Oral bioavailability F (FIXED at 1)")                                  # Sikma 2020 Table 3 F = 1 Fixed

    # ---- Erythrocyte and plasma binding (Sikma 2020 Table 3) ----------------
    # Bmax is the maximum binding capacity to erythrocytes. The model
    # multiplies Bmax by HCT (fraction) to give the effective whole-blood
    # capacity. Source paper labels Bmax as pg/mL and Kd as pg/mL in Table 3,
    # but the published equation UPC = WBC * Kd / (Bmax * Ht - WBC) is only
    # dimensionally consistent if Bmax and Kd are in the same concentration
    # unit as WBC. With WBC in ng/mL, Bmax = 2700 (ng/mL) * Ht 0.30 = 810
    # ng/mL gives a positive denominator at typical WBC = 9.5 ng/mL and
    # reproduces the simulated Figure 4 (UPC 1.06-2.14 pg/mL across HCT
    # 0.25-0.50 at WBC 9 ng/mL). The numeric values 2700 and 0.142 are
    # therefore carried as-is; the units interpretation is documented in
    # vignette Assumptions and deviations.
    lbmax    <- log(2700)        ; label("Maximum erythrocyte binding capacity Bmax (ng/mL whole blood equivalent; Sikma 2020 Table 3 label pg/mL)")    # Sikma 2020 Table 3 Bmax = 2700 (95% CI 1750-3835)
    lkd      <- log(0.142)       ; label("Erythrocyte binding dissociation constant Kd (ng/mL whole blood equivalent; Sikma 2020 Table 3 label pg/mL)") # Sikma 2020 Table 3 Kd = 0.142 (95% CI 0.087-0.195)
    lnplasma <- log(137)         ; label("Linear non-specific binding constant Nplasma = TPC / UPC (unitless ratio)")                                    # Sikma 2020 Table 3 Nplasma = 137 (95% CI 120-152)

    # ---- Inter-individual variability ---------------------------------------
    # Source paper reports IIV as %CV; converted to log-scale variance via
    # omega^2 = log(1 + CV^2). "Fixed" entries from Table 3 are wrapped in
    # `fixed(...)` so the source's structural choice is preserved.
    #   CL      CV 42.1%       -> log(1 + 0.421^2)  = 0.16320
    #   V1      CV 10% Fixed   -> log(1 + 0.10^2)   = 0.00995
    #   Q       CV 10% Fixed   -> log(1 + 0.10^2)   = 0.00995
    #   V2      CV 10% Fixed   -> log(1 + 0.10^2)   = 0.00995
    #   Ka      CV 10% Fixed   -> log(1 + 0.10^2)   = 0.00995
    #   F       CV 10% Fixed   -> log(1 + 0.10^2)   = 0.00995
    #   Bmax    CV 27%         -> log(1 + 0.27^2)   = 0.07034
    #   Kd      CV  3% Fixed   -> log(1 + 0.03^2)   = 0.000899
    #   Nplasma CV 29%         -> log(1 + 0.29^2)   = 0.08076
    etalcl      ~ 0.16320                                       # Sikma 2020 Table 3 IIV CL = 42.1% (95% CI 30-60)
    etalvc      ~ fixed(0.00995)                                # Sikma 2020 Table 3 IIV V1 = 10% Fixed
    etalq       ~ fixed(0.00995)                                # Sikma 2020 Table 3 IIV Q = 10% Fixed
    etalvp      ~ fixed(0.00995)                                # Sikma 2020 Table 3 IIV V2 = 10% Fixed
    etalka      ~ fixed(0.00995)                                # Sikma 2020 Table 3 IIV Ka = 10% Fixed
    etalfdepot  ~ fixed(0.00995)                                # Sikma 2020 Table 3 IIV F = 10% Fixed
    etalbmax    ~ 0.07034                                       # Sikma 2020 Table 3 IIV Bmax = 27% (95% CI 19-36)
    etalkd      ~ fixed(0.000899)                               # Sikma 2020 Table 3 IIV Kd = 3% Fixed
    etalnplasma ~ 0.08076                                       # Sikma 2020 Table 3 IIV Nplasma = 29% (95% CI 22-41)

    # ---- Residual error (Sikma 2020 Table 3) --------------------------------
    # Proportional residual on the linear concentration scale, separate per
    # matrix. Source paper additionally reports correlated residual error
    # (R(WBC,UPC) = 0.26, R(WBC,TPC) = 0.51, R(UPC,TPC) = 0.51) via the L2
    # NONMEM data option; residual-error correlations across distinct
    # observation outputs are not part of the standard nlmixr2lib residual
    # model and are documented in vignette Assumptions and deviations rather
    # than encoded here.
    propSd      <- 0.167   ; label("Proportional residual SD for whole-blood concentration Cc (fraction)")  # Sikma 2020 Table 3 RUV WBC = 16.7% (95% CI 15.8-17.6)
    propSd_Cupc <- 0.363   ; label("Proportional residual SD for unbound plasma concentration Cupc (fraction)") # Sikma 2020 Table 3 RUV UPC = 36.3% (95% CI 33.9-40.4)
    propSd_Ctpc <- 0.316   ; label("Proportional residual SD for total plasma concentration Ctpc (fraction)")   # Sikma 2020 Table 3 RUV TPC = 31.6% (95% CI 28.6-34.2)
  })

  model({
    # ---- Individual structural parameters ----------------------------------
    cl      <- exp(lcl      + etalcl)
    vc      <- exp(lvc      + etalvc)
    q       <- exp(lq       + etalq)
    vp      <- exp(lvp      + etalvp)
    ka      <- exp(lka      + etalka)
    fdepot  <- exp(lfdepot  + etalfdepot)
    bmax    <- exp(lbmax    + etalbmax)
    kd      <- exp(lkd      + etalkd)
    nplasma <- exp(lnplasma + etalnplasma)

    # Micro-constants for the two-compartment disposition
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---- ODE system: two-compartment first-order oral absorption -----------
    # Whole-blood concentration drives the kinetics (CL/F, V1/F, Q/F, V2/F);
    # unbound and total plasma are computed algebraically from WBC under the
    # binding-equilibrium assumption (Sikma 2020 Methods, Mixed-Effects
    # Modeling: "whole-blood, total plasma, and unbound plasma concentrations
    # were assumed to be in equilibrium all the time").
    d/dt(depot)       <- -ka  * depot
    d/dt(central)     <-  ka  * depot     - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central   - k21 * peripheral1

    # Bioavailability fixed at 1 (Sikma 2020 Table 3); IIV carried via etalfdepot.
    f(depot) <- fdepot

    # ---- Observation variables --------------------------------------------
    # Cc:   whole-blood concentration in ng/mL. Dose mg, volume L -> central
    #       / vc is mg/L; multiply by 1000 to express in ng/mL.
    # Cupc: unbound plasma concentration in pg/mL. The saturable equilibrium
    #       UPC = WBC * Kd / (Bmax * HCT - WBC) returns UPC in the same
    #       concentration unit as WBC (ng/mL); multiply by 1000 for pg/mL.
    # Ctpc: total plasma concentration in pg/mL. TPC = Nplasma * UPC; both in
    #       pg/mL.
    Cc   <- central / vc * 1000
    Cupc <- Cc * kd / (bmax * HCT - Cc) * 1000
    Ctpc <- nplasma * Cupc

    Cc   ~ prop(propSd)
    Cupc ~ prop(propSd_Cupc)
    Ctpc ~ prop(propSd_Ctpc)
  })
}
