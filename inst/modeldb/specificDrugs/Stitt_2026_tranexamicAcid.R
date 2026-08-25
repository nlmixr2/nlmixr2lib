Stitt_2026_tranexamicAcid <- function() {
  description <- "Two-compartment intravenous population PK model for tranexamic acid (TXA) with first-order elimination and allometric body-weight scaling (exponents fixed at 0.75 on the clearances and 1 on the volumes, 70 kg reference), estimated in adults with severe traumatic injury (TAMPITI trial) and extrapolated by Stitt 2026 to children with trauma-related bleeding; platelet count, near-infrared-spectroscopy skeletal-muscle tissue oxygen saturation and interleukin-8 act on clearance (Stitt 2026)."
  reference   <- paste(
    "Stitt G, Downes K, Zuppa A, Leeper C, Watt K, Spinella P.",
    "Tranexamic acid dosing in pediatric trauma: dose simulation based on",
    "population pharmacokinetic modeling in adult trauma patients.",
    "Transfusion. 2026;66(Suppl. 1):S257-S265. doi:10.1111/trf.70047.",
    "The adult structural, covariate, IIV and residual-error estimates encoded",
    "here were first published in Stitt G, Spinella PC, Bocchicchio GV,",
    "Roberts I, Downes KJ, Zuppa AF. Population pharmacokinetic modelling and",
    "simulation of tranexamic acid in adult trauma patients.",
    "Br J Clin Pharmacol. 2024;90:1932-1941. doi:10.1111/bcp.16075",
    "(Stitt 2026 reference 22); Stitt 2026 reproduces the complete final",
    "parameter set in its Table 1 and Equations 1-4 and adds the allometric",
    "scaling layer.",
    sep = " "
  )
  vignette    <- "Stitt_2026_tranexamicAcid"
  units       <- list(time = "min", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Actual body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric power scaling against a 70 kg reference weight, applied to all four disposition parameters: (WT/70)^0.75 on CL and Q, (WT/70)^1 on V1 and V2 (Stitt 2026 Equations 1-4, and Supporting Information 'Adult Population PK Model' / 'Scaling Adult Population PK Model to Children', which state that theta_allometric was FIXED at 0.75 for clearance parameters and 1 for volume parameters and that the same 70 kg reference weight was used for the paediatric extrapolation). Time-fixed at the admission value. The TAMPITI adults on whom the model was estimated had a median weight of 80.1 kg (Stitt 2026 Table 2); the virtual paediatric cohort spans 17.7-58.2 kg.",
      source_name        = "WT"
    ),
    PLT = list(
      description        = "Platelet count at time 0 (admission)",
      units              = "10^9 cells/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on clearance, (PLT/196)^0.468 (Stitt 2026 Equation 1). Stitt 2026 reports platelet count in K/uL, which is numerically identical to the canonical 10^9 cells/L, so no value transformation is needed. Baseline only (the admission value at time 0; Stitt 2026 Table 2 footnote a). The centring constant 196 is the value printed inside Equation 1; Table 2 gives the TAMPITI median admission PLT as 197 K/uL. The 0.5% difference is numerically immaterial and the equation value is used per the standing trust-the-printed-equation policy. Higher platelet count predicts higher TXA clearance.",
      source_name        = "PLT"
    ),
    STO2 = list(
      description        = "Skeletal-muscle tissue oxygen saturation measured by near-infrared spectroscopy (NIRS)",
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on clearance, (STO2/88)^-0.29 (Stitt 2026 Equation 1). Baseline only. IMPORTANT non-equivalence that the paper itself flags (Discussion, S262): the adult reference 88% is the median ADMISSION NIRS value in TAMPITI participants, whereas the paediatric range 51-80% used for the virtual cohort is the median LOWEST NIRS value recorded in the paediatric literature sources -- the two are not the same quantity, and the lower paediatric values are part of why the extrapolated paediatric clearance is fast. Lower tissue oxygen saturation predicts higher TXA clearance.",
      source_name        = "NIRS"
    ),
    IL8 = list(
      description        = "Interleukin-8 (CXCL8) concentration at time 0 (admission)",
      units              = "pg/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on clearance, entered UN-NORMALIZED as IL8^-0.0873 (Stitt 2026 Equation 1) -- i.e. with an implicit reference of 1 pg/mL, unlike the PLT and STO2 terms in the same equation, which are median-normalized. This is not a transcription slip: encoding the raw form reproduces the paper's own Table 3 adult AUC0-4h and AUC0-8h to +1.9%, whereas a median-normalized (IL8/20.3) reading raises typical adult clearance from 162 to 211 mL/min and undershoots the published AUCs by roughly 15%. Because the term is un-normalized, the column units must be pg/mL exactly for the published exponent to apply. Baseline only (the admission value at time 0). Higher IL-8 predicts LOWER TXA clearance; note that the Stitt 2026 Discussion prose (S262) states the opposite direction, contradicting both its own Table 1 (-0.0887) and Equation 1 (-0.0873) -- the negative exponent is used here.",
      source_name        = "IL-8"
    )
  )

  compartmentData <- list(
    central     = list(analyte = "tranexamic acid", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tranexamic acid", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 94L,
    n_studies      = 1L,
    n_observations = 597L,
    weight_median  = "80.1 kg (TAMPITI median at time 0)",
    disease_state  = "Adults with severe traumatic injury enrolled in the Tranexamic Acid Mechanisms and Pharmacokinetics in Traumatic Injury (TAMPITI) trial",
    dose_range     = "2 g or 4 g intravenous TXA bolus (Stitt 2026 Methods, Adult population PK model)",
    notes          = "Estimation population: 597 plasma TXA samples from 94 TAMPITI participants (placebo participants excluded), per Stitt 2026 Methods 'Adult population PK model'. Median covariate values at time 0 (Stitt 2026 Table 2): weight 80.1 kg, platelet count 197 K/uL, NIRS 88%, IL-8 20.3 pg/mL. The final estimates are reproduced in Stitt 2026 Table 1 and Equations 1-4; they were first reported in Stitt 2024 (Br J Clin Pharmacol 90:1932-1941, doi:10.1111/bcp.16075). Stitt 2026 itself contains NO new observed pharmacokinetic data: it is a simulation study that allometrically scales this adult model to children. The virtual paediatric population it simulates (Stitt 2026 Table 2 and Supporting Information) covers age 4.7-15.4 years, weight 17.7-58.2 kg, platelet count 205-433 K/uL, NIRS 51-80% and IL-8 6.6-50.2 pg/mL, with each covariate drawn independently and uniformly from its range; paediatric weights were derived from the MATIC-1 trauma-cohort age distribution via CDC 50th-percentile weight-for-age. There are no observed paediatric trauma TXA PK data anywhere in the source, so every paediatric prediction from this model is an extrapolation outside the estimation range and, as the paper states, requires clinical validation."
  )

  ini({
    # Structural parameters. Stitt 2026 Equation 1-4 (Results, p.S260), which
    # print the final typical-value equations in full. Table 1 lists the same
    # estimates with RSE and 95% CI. The equations are reported in mL/min (CL,
    # Q) and mL (V1, V2) at the 70 kg reference weight; they are divided by
    # 1000 here so that the model works in L and L/min and therefore returns
    # mg/L from a mg dose, which is the concentration unit used throughout the
    # paper's figures and NCA tables.
    #
    # Where Table 1 and Equation 1 disagree the equation value is used, per the
    # standing trust-the-printed-equation policy; both are recorded below. Both
    # differences are far too small to be numerically discriminating.
    lcl <- log(190 / 1000)
    label("Clearance at WT 70 kg, PLT 196 x 10^9/L, StO2 88%, IL-8 1 pg/mL (L/min)")            # Stitt 2026 Equation 1: 190 mL/min. Table 1 prints 192 mL/min/70 kg (RSE 10.2%, 95% CI 154-230).
    lvc <- log(17300 / 1000)
    label("Central volume of distribution V1 at WT 70 kg (L)")                                  # Stitt 2026 Equation 2: 17,300 mL. Table 1: 17,300 mL/70 kg (RSE 4.94%, 95% CI 15,600-19,000).
    lq  <- log(80.1 / 1000)
    label("Inter-compartmental clearance Q at WT 70 kg (L/min)")                                # Stitt 2026 Equation 3: 80.1 mL/min. Table 1: 80.1 mL/min/70 kg (RSE 7.12%, 95% CI 68.9-91.3).
    lvp <- log(11400 / 1000)
    label("Peripheral volume of distribution V2 at WT 70 kg (L)")                               # Stitt 2026 Equation 4: 11,400 mL. Table 1: 11,400 mL/70 kg (RSE 4.31%, 95% CI 10,400-12,400).

    # Allometric exponents. Both are FIXED, not estimated: they carry no RSE or
    # 95% CI in Table 1, and the Supporting Information states outright that
    # "theta_allometric was fixed at 0.75 for clearance parameters and 1 for
    # volume parameters". A single shared exponent is applied to each pair, so
    # the shared-exponent naming convention (e_<cov>_<param1>_<param2>) is used.
    e_wt_cl_q  <- fixed(0.75)
    label("Allometric exponent on (WT/70) shared by CL and Q (unitless)")                       # Stitt 2026 Equations 1 and 3; Supporting Information "Adult Population PK Model" and "Scaling Adult Population PK Model to Children"
    e_wt_vc_vp <- fixed(1)
    label("Allometric exponent on (WT/70) shared by V1 and V2 (unitless)")                      # Stitt 2026 Equations 2 and 4; Supporting Information (same sections)

    # Covariate effects on clearance. All three were estimated (Table 1 reports
    # an RSE and a 95% CI for each). Note the asymmetry in Equation 1: PLT and
    # NIRS are median-normalized but IL-8 is NOT -- see covariateData$IL8$notes
    # for the arithmetic that confirms the un-normalized form against the
    # paper's own Table 3 output.
    e_plt_cl  <-  0.468
    label("Power exponent on (PLT/196) for CL (unitless)")                                      # Stitt 2026 Equation 1; Table 1 "PLT on CL" 0.468 (RSE 14.6%, 95% CI 0.335-0.601)
    e_sto2_cl <- -0.29
    label("Power exponent on (STO2/88) for CL (unitless)")                                      # Stitt 2026 Equation 1; Table 1 "NIRS on CL" -0.29 (RSE 15.7%, 95% CI -0.379 to -0.201)
    e_il8_cl  <- -0.0873
    label("Power exponent on un-normalized IL8 (pg/mL) for CL (unitless)")                      # Stitt 2026 Equation 1: -0.0873. Table 1 prints -0.0887 (RSE 25.9%, 95% CI -0.134 to -0.0436).

    # Inter-individual variability, on the VARIANCE scale. Stitt 2026 Table 1
    # reports one "Inter-individual variability" number per parameter, and
    # Equations 1-4 typeset the exponential IIV term by substituting that same
    # number into the exponent (e^0.106, e^0.0688, e^0.879, e^0.0589), which
    # fixes the column as the quantity in the variance slot but not its scale.
    # Two independent checks settle the scale as omega^2 rather than omega:
    #   (a) The RSE on the CL IIV is 14.7%. The asymptotic RSE floor for a
    #       VARIANCE estimated from N = 94 subjects is sqrt(2/94) = 14.6%; for
    #       an SD it would be sqrt(1/(2*94)) = 7.3%. The reported value sits on
    #       the variance floor.
    #   (b) The published adult Cmax interquartile range (Table 3: 98.9-139.7
    #       mg/L) has Q3/Q1 = 1.413. Cmax after a 1-minute bolus is essentially
    #       Dose/V1, so for log-normal V1 the ratio is exp(2 * 0.6745 * omega).
    #       With omega^2 = 0.0688 that predicts 1.424 (0.8% from observed);
    #       reading 0.0688 as an SD predicts 1.097, which is impossible.
    # The paper reports no off-diagonal covariances, so the block is diagonal.
    etalcl ~ 0.106   # Stitt 2026 Table 1, IIV on CL (RSE 14.7%, 95% CI 0.0754-0.137; shrinkage 2.7%); also Equation 1 e^0.106
    etalvc ~ 0.0688  # Stitt 2026 Table 1, IIV on V1 (RSE 31.3%, 95% CI 0.0267-0.111; shrinkage 13.4%); also Equation 2 e^0.0688
    etalq  ~ 0.879   # Stitt 2026 Table 1, IIV on Q  (RSE 29.5%, 95% CI 0.371-1.39; shrinkage 14.7%); also Equation 3 e^0.879
    etalvp ~ 0.0589  # Stitt 2026 Table 1, IIV on V2 (RSE 36.5%, 95% CI 0.0168-0.101; shrinkage 15.1%); also Equation 4 e^0.0589

    # Residual error. Table 1 lists the proportional residual error in the same
    # column as the IIV terms, so it is likewise a variance and the nlmixr2 SD
    # is its square root: sqrt(0.0238) = 0.154, i.e. 15.4% proportional error.
    # This is independently corroborated by the published Cmax values, which
    # sit ~16% above the typical-value Dose/V1 in BOTH the adult (117.1 vs
    # 101.0 mg/L) and the paediatric 25 mg/kg (114.6 vs 101.2 mg/L) arms: an
    # NCA Cmax taken as the maximum over a noisy sampling grid is biased upward
    # by roughly one to one-and-a-half residual SDs, which 15.4% explains and
    # a 2.4% SD reading (which would give ~+3%) cannot.
    propSd <- sqrt(0.0238)
    label("Proportional residual error (fraction)")                                             # Stitt 2026 Table 1, Residual error / Proportional 0.0238 (RSE 4.5%, 95% CI 0.0217-0.0259; shrinkage 21.5%)
  })

  model({
    # Individual parameters, exactly as printed in Stitt 2026 Equations 1-4.
    # All four covariate terms act on clearance; the volumes and the
    # inter-compartmental clearance carry allometric weight scaling only.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q *
      (PLT / 196)^e_plt_cl * (STO2 / 88)^e_sto2_cl * IL8^e_il8_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq  + etalq)  * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp

    # Micro-constants for the explicit two-compartment system.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Intravenous administration only: doses enter the central compartment
    # directly (Stitt 2026 simulates IV boluses given over 1-10 minutes).
    # No depot and no bioavailability term.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Dose in mg over vc in L gives mg/L, the concentration unit used in every
    # Stitt 2026 figure and NCA table.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
