Alli_2023_cefazolin <- function() {
  description <- "Three-compartment joint total/unbound population PK model for intravenous cefazolin surgical prophylaxis in adults undergoing elective cardiac valve replacement with cardiopulmonary bypass (Alli 2023). The dose enters a central compartment whose concentration is the TOTAL plasma cefazolin; it exchanges first-order with a peripheral compartment and with a separate unbound-drug compartment (volume Vt) whose concentration is the UNBOUND plasma cefazolin. The return rate constant from the unbound compartment to the central (total) compartment is the protein-binding term PB = (ktc * Bmax + NS) * (ALB / 40), so higher serum albumin pulls drug back into the total pool. Clearance acts on the central compartment and scales inversely with serum creatinine, CL = TVCL * (90 / SeCr). Estimated with the Pmetrics non-parametric adaptive grid (NPAG); the non-parametric distribution is approximated here by log-normal IIV from the reported CV%."
  reference <- "Alli A, Paruk F, Roger C, Lipman J, Calleemalay D, Wallis SC, Scribante J, Richards GA, Roberts JA, Ballotta A. Peri-operative pharmacokinetics of cefazolin prophylaxis during valve replacement surgery. PLoS One. 2023;18(9):e0291425. doi:10.1371/journal.pone.0291425"
  vignette <- "Alli_2023_cefazolin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(
      analyte = "cefazolin (total)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE,
      notes = "The 'dosing compartment (volume: Vc)' of Alli 2023 Figure 2. Dividing by Vc gives the TOTAL plasma concentration, which is what the paper's total-concentration observations were fitted to. Elimination (Total clearance, Cl) acts here."
    ),
    peripheral1 = list(
      analyte = "cefazolin",
      units = "mg",
      specimen = "tissue",
      verified = TRUE,
      notes = "The 'Peripheral compartment' of Alli 2023 Figure 2, exchanging with the central compartment through kcp / kpc. Unobserved; its volume is not a model parameter (the paper derives Vp = Vc * kcp / kpc after the fit)."
    ),
    peripheral2 = list(
      analyte = "cefazolin (unbound)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE,
      notes = "The 'Unbound drug compartment (volume: Vt)' of Alli 2023 Figure 2. Structurally the third compartment of the paper's 'three-compartment linear model'; dividing by Vt gives the UNBOUND plasma concentration fitted to the ultracentrifugation (Centrifree) measurements. Filled from central by kct and returned to central by the albumin-scaled protein-binding rate constant PB."
    )
  )

  covariateData <- list(
    CREAT = list(
      description = "Serum creatinine concentration",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Pre-operative value, time-fixed per subject. Enters clearance as the inverse-proportional term CL = TVCL * (90 / SeCr), printed in the Alli 2023 Results ('The final clearance model was described as follows'). The text calls 90 umol/L 'the population mean value', but Table 1 reports a cohort mean of 76.7 umol/L (SD 18.6, range 42-102); the equation's 90 is used as printed. The paper's dosing simulations use SeCr 50, 75 and 100 umol/L.",
      source_name = "SeCr"
    ),
    ALB = list(
      description = "Serum albumin concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Pre-operative value, time-fixed per subject. Scales the protein-binding return rate constant from the unbound to the central (total) compartment, PB = (ktc * Bmax + NS) * (ALB / 40), printed in the Alli 2023 Results. The reference 40 g/L is called 'the population mean value' in the text; Table 1 reports 41.6 g/L (SD 4.57, range 33-47). Higher albumin raises PB, shifting drug from the unbound into the total pool. The paper's dosing simulations use 30, 35 and 40 g/L.",
      source_name = "ALB"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 16L,
    n_studies = 1L,
    n_samples = "195 total and 64 unbound plasma cefazolin concentrations",
    age_range = "18-59 years (mean 44.0, SD 12.0)",
    weight_range = "51-97 kg (mean 75.1, SD 15.6)",
    sex_female_pct = 50,
    disease_state = "Adults (18-60 years eligible) undergoing cardiopulmonary bypass (CPB) for elective cardiac valve replacement, receiving cefazolin surgical prophylaxis. CPB duration mean 165 min (range 105-278). Excluded: pregnancy, Jehovah's Witnesses, dialysis, cefazolin in the previous 72 h.",
    renal_function = "Serum creatinine mean 76.7 umol/L (SD 18.6, range 42-102); Cockcroft-Gault creatinine clearance mean 81.5 mL/min (range 42.5-124.63).",
    albumin = "Serum albumin mean 41.6 g/L (SD 4.57, range 33-47).",
    dose_range = "2 g cefazolin IV (Figure 2: 1-min bolus infusion) 30 min before skin incision, re-dosed 2 g after 4 h if surgery was ongoing.",
    regions = "South Africa (Charlotte Maxeke Johannesburg Academic Hospital)",
    notes = "Demographics from Alli 2023 Table 1. Prospective convenience-sampled PK study; plasma sampled pre-dose and 2-300 min after the dose from the central venous catheter (plus CPB arterial-port samples at 90-180 min). Total and unbound (ultracentrifugation, Centrifree 30 kDa) cefazolin by UHPLC-MS/MS, total 1-500 mg/L, unbound 0.1-500 mg/L. CPB flow 60-70 mL/kg/min at 30-32 C. Model fitted in Pmetrics (NPAG)."
  )

  ini({
    # Structural parameters: Alli 2023 Table 2 'Mean' column of the Pmetrics
    # NPAG population distribution (final covariate model). The mean is the
    # value the Results text quotes (CL 3.23 L/h, Vc 3.38 L, Vt 6.08 L,
    # NS 15.14); the Table 2 median is noted per line. The Table 2 values
    # are the typical values at the covariate references SeCr = 90 umol/L
    # and ALB = 40 g/L.
    lcl <- log(3.23)
    label("Clearance from the central (total) compartment at SeCr 90 umol/L (L/h)") # Table 2 Cl mean 3.23 (SD 1.16, CV 36.06%, median 3.00)
    lvc <- log(3.38)
    label("Central (total-drug) volume Vc (L)") # Table 2 Vc mean 3.38 (SD 0.73, CV 21.59%, median 3.41)
    lvp2 <- log(6.08)
    label("Volume of the unbound-drug compartment Vt (L)") # Table 2 Vt mean 6.08 (SD 2.03, CV 33.42%, median 6.28)
    lk12 <- log(5.92)
    label("Central-to-peripheral rate constant kcp (1/h)") # Table 2 kcp mean 5.92 (SD 3.88, CV 65.58%, median 4.59)
    lk21 <- log(2.24)
    label("Peripheral-to-central rate constant kpc (1/h)") # Table 2 kpc mean 2.24 (SD 1.16, CV 51.81%, median 2.35)
    lk13 <- log(144.27)
    label("Central-to-unbound-compartment rate constant kct (1/h)") # Table 2 kct mean 144.27 (SD 14.67, CV 10.17%, median 143.84)
    lk31 <- log(195.91)
    label("Unbound-to-central rate constant ktc, scaled by Bmax inside PB (1/h)") # Table 2 ktc mean 195.91 (SD 25.34, CV 12.94%, median 188.90)
    lbmax <- log(1.01)
    label("Maximal binding capacity Bmax multiplying ktc in PB (g/L, as tabulated)") # Table 2 Bmax mean 1.01 (SD 0.46, CV 45.23%, median 0.87)
    lkns <- log(15.14)
    label("Non-saturable protein-binding term NS added in PB (1/h)") # Table 2 NS mean 15.14 (SD 3.75, CV 24.79%, median 15.20); no unit printed

    # IIV. Pmetrics NPAG estimates a discrete non-parametric distribution;
    # Table 2 summarises it by mean, SD and CV%. Each CV% is carried into a
    # log-normal random effect through omega^2 = log(CV^2 + 1). This is a
    # parametric approximation of the non-parametric distribution; no
    # covariances are reported, so the etas are independent.
    etalcl ~ 0.122246 # Table 2 Cl CV 36.06%
    etalvc ~ 0.045559 # Table 2 Vc CV 21.59%
    etalvp2 ~ 0.105881 # Table 2 Vt CV 33.42%
    etalk12 ~ 0.357726 # Table 2 kcp CV 65.58%
    etalk21 ~ 0.237778 # Table 2 kpc CV 51.81%
    etalk13 ~ 0.010290 # Table 2 kct CV 10.17%
    etalk31 ~ 0.016606 # Table 2 ktc CV 12.94%
    etalbmax ~ 0.186127 # Table 2 Bmax CV 45.23%
    etalkns ~ 0.059640 # Table 2 NS CV 24.79%

    # Residual error: Pmetrics assay-error polynomial SD = C0 + C1 * obs,
    # printed in the Results with C0 = 1.0 mg/L and C1 = 0.1 for total
    # plasma data and C0 = 0.1 mg/L and C1 = 0.1 for unbound data. The
    # polynomial is multiplied by an estimated noise factor gamma whose
    # value the paper does not report; it is carried here unscaled
    # (gamma = 1). The linear sum of the two terms is nlmixr2's combined1
    # form.
    addSd <- 1.0
    label("Additive residual error on total Cc (mg/L)") # Results, assay error C0 = 1.0 mg/L (total)
    propSd <- 0.1
    label("Proportional residual error on total Cc (fraction)") # Results, assay error C1 = 0.1 (total)
    addSd_Cu <- 0.1
    label("Additive residual error on unbound Cu (mg/L)") # Results, assay error C0 = 0.1 mg/L (unbound)
    propSd_Cu <- 0.1
    label("Proportional residual error on unbound Cu (fraction)") # Results, assay error C1 = 0.1 (unbound)
  })

  model({
    # Covariate references printed in the Results equations.
    creat_ref <- 90
    alb_ref <- 40

    # Individual parameters.
    cl <- exp(lcl + etalcl) * (creat_ref / CREAT)
    vc <- exp(lvc + etalvc)
    vp2 <- exp(lvp2 + etalvp2)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    k13 <- exp(lk13 + etalk13)
    k31 <- exp(lk31 + etalk31)
    bmax <- exp(lbmax + etalbmax)
    kns <- exp(lkns + etalkns)

    # Protein-binding rate constant from the unbound to the total (central)
    # compartment, PB = (ktc * Bmax + NS) * (ALB / 40) (Results). Figure 2
    # draws Bmax and NS acting on the unbound -> central arrow, so PB
    # replaces ktc as that arrow's rate constant.
    kpb <- (k31 * bmax + kns) * (ALB / alb_ref)

    kel <- cl / vc

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 -
      k13 * central + kpb * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - kpb * peripheral2

    # Total plasma concentration from the central compartment; unbound
    # plasma concentration from the unbound-drug compartment.
    Cc <- central / vc
    Cu <- peripheral2 / vp2

    Cc ~ add(addSd) + prop(propSd) + combined1()
    Cu ~ add(addSd_Cu) + prop(propSd_Cu) + combined1()
  })
}
