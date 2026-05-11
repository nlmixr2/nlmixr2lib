Quartino_2016_trastuzumab <- function() {
  description <- "Two-compartment population PK model with parallel linear and Michaelis-Menten nonlinear elimination from the central compartment and first-order subcutaneous absorption (with bioavailability) for trastuzumab (Herceptin) administered IV or as a fixed 600 mg manual-syringe SC dose in women with HER2-positive early breast cancer; covariates body weight (on CL, Vc, Vp) and ALT (on CL) (Quartino 2016, HannaH study)"
  reference <- "Quartino AL, Hillenbach C, Li J, Li H, Wada DR, Visich J, Li C, Heinzmann D, Jin JY, Lum BL. Population pharmacokinetic and exposure-response analysis for trastuzumab administered using a subcutaneous 'manual syringe' injection or intravenously in women with HER2-positive early breast cancer. Cancer Chemother Pharmacol. 2016;77(1):77-88. doi:10.1007/s00280-015-2922-5"
  vignette <- "Quartino_2016_trastuzumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effects on linear CL (exponent 1.04), Vc (exponent 0.443), and Vp (exponent 0.500); reference 68 kg per Quartino 2016 Table 1 and the in-text covariate equations CLi = 0.111 * (WTi/68)^1.04 * (ALTi/19)^0.144, Vci = 2.91 * (WTi/68)^0.443, Vpi = 3.06 * (WTi/68)^0.500.",
      source_name        = "WT"
    ),
    ALT = list(
      description        = "Baseline serum alanine aminotransferase activity",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effect on linear CL only (exponent 0.144); reference 19 IU/L per Quartino 2016 in-text covariate equation CLi = 0.111 * (WTi/68)^1.04 * (ALTi/19)^0.144. No ALT effect on Vc or Vp in the final model.",
      source_name        = "ALT"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 592L,
    n_studies        = 1L,
    n_observations   = 15761L,
    study            = "HannaH (NCT00950300); phase III, randomized, international, open-label, neoadjuvant-adjuvant trial.",
    age_range        = "Not reported in the main publication; baseline demographics in Online Resource 3 (not on disk).",
    weight_range     = "Spans at least <58 kg to >85 kg per the body-weight quartile analysis in Results; reference (median-equivalent) is 68 kg per Table 1.",
    weight_median    = "68 kg (reference subject of the Quartino 2016 covariate model)",
    sex_female_pct   = 100,
    disease_state    = "HER2-positive operable, locally advanced, or inflammatory early breast cancer (EBC).",
    dose_range       = "Subcutaneous (SC): fixed 600 mg trastuzumab q3w by manual handheld syringe (5 min). Intravenous (IV): 8 mg/kg loading dose, 6 mg/kg maintenance q3w (90 min then 30 min infusions). Both arms received 8 cycles neoadjuvant (with concomitant chemotherapy) + 10 cycles adjuvant trastuzumab monotherapy.",
    regions          = "International, multicentre HannaH study (sites in multiple countries).",
    reference_subject = "68 kg, ALT 19 IU/L per Quartino 2016 covariate equations and the 'typical patient' used in Figure 3 typical-subject simulations.",
    notes            = "Baseline demographics in Quartino 2016 Online Resource 3 (not on disk). HannaH enrolled 596 patients; 595 received at least one dose (297 SC, 298 IV). After outlier and BLQ handling the PK analysis dataset contained 592 patients and 15,761 trastuzumab serum concentrations (Quartino 2016 Results 'Patient population and PK samples'). All patients are women with HER2-positive EBC. Co-medication with anthracycline-, taxane-, or other neoadjuvant chemotherapy during cycles 1-8 (Patients and methods). ATA- and AHA-positivity rates (5.0% and 5.7%) were too low to detect immunogenicity effects on PK."
  )

  ini({
    # Structural parameters (Quartino 2016 Table 1, FOCEI final estimates).
    # Reference subject: 68 kg female with EBC, ALT 19 IU/L. Concentration in
    # the central compartment is Cc = central / vc, with dose in mg and
    # volumes in L -> Cc in mg/L (= ug/mL). Vmax is mg/day and Km is mg/L so
    # both elimination fluxes Vmax * Cc / (Km + Cc) and kel * central are in
    # mg/day. SC bioavailability F applies only to depot-entered doses.
    lcl     <- log(0.111); label("Linear CL for the reference subject (L/day)")                   # Quartino 2016 Table 1
    lvc     <- log(2.91);  label("Central volume of distribution Vc for the reference subject (L)") # Quartino 2016 Table 1
    lq      <- log(0.445); label("Intercompartmental clearance Q (L/day)")                        # Quartino 2016 Table 1
    lvp     <- log(3.06);  label("Peripheral volume of distribution Vp for the reference subject (L)") # Quartino 2016 Table 1
    lvmax   <- log(11.9);  label("Maximum nonlinear (Michaelis-Menten) elimination rate Vmax (mg/day)") # Quartino 2016 Table 1
    lkm     <- log(33.9);  label("Michaelis-Menten constant Km (mg/L = ug/mL)")                   # Quartino 2016 Table 1
    lka     <- log(0.404); label("First-order SC absorption rate Ka (1/day)")                     # Quartino 2016 Table 1
    lfdepot <- log(0.771); label("Subcutaneous bioavailability F (fraction)")                     # Quartino 2016 Table 1

    # Covariate exponents on the structural parameters (Quartino 2016 Table 1
    # and the in-text covariate equations
    #   CLi = 0.111 * (WTi/68)^1.04 * (ALTi/19)^0.144
    #   Vci = 2.91  * (WTi/68)^0.443
    #   Vpi = 3.06  * (WTi/68)^0.500
    # ). Reference WT = 68 kg, reference ALT = 19 IU/L.
    e_wt_cl  <- 1.04;  label("Power exponent of WT on linear CL (unitless; reference 68 kg)")    # Quartino 2016 Table 1
    e_wt_vc  <- 0.443; label("Power exponent of WT on Vc (unitless; reference 68 kg)")           # Quartino 2016 Table 1
    e_wt_vp  <- 0.500; label("Power exponent of WT on Vp (unitless; reference 68 kg)")           # Quartino 2016 Table 1
    e_alt_cl <- 0.144; label("Power exponent of ALT on linear CL (unitless; reference 19 IU/L)") # Quartino 2016 Table 1

    # Inter-individual variability. Quartino 2016 Table 1 reports omega as
    # %CV on log-normal parameters; convert via omega^2 = log(CV^2 + 1).
    # IIV is reported on F (13.0% CV), linear CL (30.0% CV), Vc (19.1% CV),
    # and Vp (50.4% CV) only; no IIV reported on Ka, Vmax, Km, or Q
    # (dashes in the Table 1 BSV column). The paper states
    # "Between-subject variability was modeled using a log-normal variance
    # model" (Results), so all four etas are log-normal; the paper does not
    # report a covariance / correlation block, so the diagonal is used.
    etalfdepot ~ 0.016759  # F  13.0% CV -- Quartino 2016 Table 1; omega^2 = log(0.130^2 + 1)
    etalcl     ~ 0.086178  # CL 30.0% CV -- Quartino 2016 Table 1; omega^2 = log(0.300^2 + 1)
    etalvc     ~ 0.035831  # Vc 19.1% CV -- Quartino 2016 Table 1; omega^2 = log(0.191^2 + 1)
    etalvp     ~ 0.226351  # Vp 50.4% CV -- Quartino 2016 Table 1; omega^2 = log(0.504^2 + 1)

    # Residual error: combined proportional + additive (Quartino 2016 Table 1
    # and Methods 'Outliers were identified using ... a proportional plus
    # additive residual error model'; the combined form was retained in the
    # final model per Results). Table 1 reports proportional variability
    # 23.9% and additive variability 4.48 ug/mL as SDs (standard NONMEM
    # additive-plus-proportional convention).
    propSd <- 0.239; label("Proportional residual error (fraction)")                              # Quartino 2016 Table 1
    addSd  <- 4.48;  label("Additive residual error (ug/mL)")                                    # Quartino 2016 Table 1
  })

  model({
    # Per-subject covariate multipliers (Quartino 2016 in-text covariate
    # equations).
    cl_cov <- (WT / 68)^e_wt_cl * (ALT / 19)^e_alt_cl
    vc_cov <- (WT / 68)^e_wt_vc
    vp_cov <- (WT / 68)^e_wt_vp

    # Individual parameters
    cl   <- exp(lcl + etalcl) * cl_cov
    vc   <- exp(lvc + etalvc) * vc_cov
    vp   <- exp(lvp + etalvp) * vp_cov
    q    <- exp(lq)
    vmax <- exp(lvmax)
    km   <- exp(lkm)
    ka   <- exp(lka)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    Cc <- central / vc

    # Two-compartment model with first-order SC absorption from depot to
    # central, plus parallel linear and Michaelis-Menten elimination from
    # central (Quartino 2016 Results: "a two-compartment model with parallel
    # linear and nonlinear (Michaelis-Menten) elimination from the central
    # compartment ... SC absorption was modeled as a first-order process").
    # IV doses are entered directly into central and bypass the depot.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - vmax * Cc / (km + Cc) - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Subcutaneous bioavailability applies to depot-entered doses (SC route);
    # IV doses entered into central bypass this factor.
    f(depot) <- exp(lfdepot + etalfdepot)

    Cc ~ add(addSd) + prop(propSd)
  })
}
