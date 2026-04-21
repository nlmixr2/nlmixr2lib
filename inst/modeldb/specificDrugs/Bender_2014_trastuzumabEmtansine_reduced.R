Bender_2014_trastuzumabEmtansine_reduced <- function() {
  description <- "Reduced three-compartment population PK model for trastuzumab emtansine (T-DM1) and naked trastuzumab (DAR0) in cynomolgus monkeys (default) and rats (Bender 2014): single lumped T-DM1 conjugate species deconjugates into DAR0 via a single deconjugation clearance; both species share V1/V2/V3 and distributional clearances."
  reference <- "Bender B, Leipold DD, Xu K, Shen BQ, Tibbitts J, Friberg LE. A mechanistic pharmacokinetic model elucidating the disposition of trastuzumab emtansine (T-DM1), an antibody-drug conjugate (ADC) for treatment of metastatic breast cancer. AAPS J. 2014;16(5):994-1008. doi:10.1208/s12248-014-9618-3"
  vignette <- "Bender_2014_trastuzumabEmtansine_reduced"
  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "ug/mL"
  )

  covariateData <- list()

  population <- list(
    n_subjects     = 18,
    n_studies      = 1,
    age_range      = "not reported",
    weight_range   = "not reported",
    sex_female_pct = NA_real_,
    species        = "Cynomolgus monkey (Macaca fascicularis); rat parameter set also provided",
    disease_state  = "Healthy preclinical animals (no tumour xenograft)",
    dose_range     = "Cynomolgus: 30 mg/kg IV single (n=4) and 10 mg/kg IV q3w x4 infusion (n=14); Rat: 0.3, 3, 20 mg/kg IV single (n=24) plus 10 mg/kg IV stability cohorts (n=10)",
    regions        = "Preclinical (Genentech/Roche in-house studies)",
    scope_note     = "Preclinical-only (no human PK). All existing specificDrugs/ entries are human except Grimm_2023_trontinemab and this model; confirmed by the operator as an intentional preclinical addition.",
    model_variant  = "Reduced T-DM1 model (Bender 2014 Table III). A separate Bender_2014_trastuzumabEmtansine_mechanistic model file encodes the DAR0-DAR7 catenary mechanistic variant (Bender 2014 Table II).",
    species_parameters = list(
      cynomolgus = list(
        source   = "Bender 2014 Table III, cynomolgus columns",
        CL_TT    = "19.9 mL/day (IIV CV 19.8%)",
        V1       = "154 mL (IIV CV 7.65%)",
        CLd2     = "56.8 mL/day (no IIV)",
        V2       = "50.0 mL (no IIV)",
        CLd3     = "60.4 mL/day (no IIV)",
        V3       = "84.7 mL (IIV CV 27.3%)",
        CL_DEC   = "22.0 mL/day (IIV CV 11.9%)",
        res_err  = "9.64% proportional"
      ),
      rat = list(
        source   = "Bender 2014 Table III, rat columns",
        CL_TT    = "2.37 mL/day (IIV CV 24.6%)",
        V1       = "10.7 mL (IIV CV 19.8%)",
        CLd2     = "59.7 mL/day (no IIV)",
        V2       = "2.52 mL (IIV CV 68.0%)",
        CLd3     = "13.9 mL/day (no IIV)",
        V3       = "15.5 mL (IIV CV 16.5%)",
        CL_DEC   = "2.24 mL/day (IIV CV 15.3%)",
        res_err  = "10.9% proportional",
        notes    = "Rat fit adds IIV on V2; to switch from the cynomolgus default to the rat parameter set, call ini() on the resulting model with the rat values above (converted from mL/day to L/day by /1000)."
      )
    ),
    notes          = "NONMEM 7.2, FOCE-INTER; residual error reported as additive on log-scale (proportional in linear space). Derived parameters (t1/2 of TT, t1/2 of T-DM1, CL_T-DM1 = CL_TT + CL_DEC) are quoted in Table III but are computed, not fit."
  )

  ini({
    # Structural parameters - cynomolgus monkey (default ini values)
    # Paper uses absolute volumes in mL and clearances in mL/day; here V is
    # stored in L (mL/1000) and CL in L/day so that dose in mg and V in L give
    # Cc = central / vc directly in mg/L = ug/mL.
    lcl    <- log(0.0199);  label("Total trastuzumab clearance (CL_TT, L/day)")             # Bender 2014 Table III, cyno: 19.9 mL/day
    lvc    <- log(0.154);   label("Central volume shared by T-DM1 and DAR0 (V1, L)")         # Bender 2014 Table III, cyno: 154 mL
    lqd2   <- log(0.0568);  label("Distributional clearance to peripheral1 (CLd2, L/day)")   # Bender 2014 Table III, cyno: 56.8 mL/day
    lvp    <- log(0.0500);  label("Peripheral volume 1 shared by T-DM1 and DAR0 (V2, L)")    # Bender 2014 Table III, cyno: 50.0 mL
    lqd3   <- log(0.0604);  label("Distributional clearance to peripheral2 (CLd3, L/day)")   # Bender 2014 Table III, cyno: 60.4 mL/day
    lvp2   <- log(0.0847);  label("Peripheral volume 2 shared by T-DM1 and DAR0 (V3, L)")    # Bender 2014 Table III, cyno: 84.7 mL
    lcldec <- log(0.0220);  label("Deconjugation clearance from T-DM1 to DAR0 (CL_DEC, L/day)") # Bender 2014 Table III, cyno: 22.0 mL/day

    # IIV - cynomolgus Table III; log-normal, omega^2 = log(CV^2 + 1)
    etalcl    ~ 0.038454   # cyno CL_TT CV 19.8%  (rat CV 24.6% -> 0.05878)
    etalvc    ~ 0.005835   # cyno V1    CV 7.65%  (rat CV 19.8% -> 0.038454)
    etalvp2   ~ 0.071878   # cyno V3    CV 27.3%  (rat CV 16.5% -> 0.026858)
    etalcldec ~ 0.014062   # cyno CL_DEC CV 11.9% (rat CV 15.3% -> 0.023151)

    # Residual error - cynomolgus; paper reports a single residual magnitude
    # applied to both T-DM1 and total trastuzumab observations
    CcpropSd  <- 0.0964; label("Proportional residual error on T-DM1 concentration (fraction)")            # Bender 2014 Table III, cyno: 9.64%
    CttpropSd <- 0.0964; label("Proportional residual error on total trastuzumab concentration (fraction)") # Bender 2014 Table III, cyno: 9.64%
  })

  model({
    cl     <- exp(lcl + etalcl)
    vc     <- exp(lvc + etalvc)
    qd2    <- exp(lqd2)
    vp     <- exp(lvp)
    qd3    <- exp(lqd3)
    vp2    <- exp(lvp2 + etalvp2)
    cldec  <- exp(lcldec + etalcldec)

    # T-DM1 three-compartment distribution with linear antibody elimination
    # (CL_TT) and deconjugation loss (CL_DEC) from the central compartment;
    # DAR0 inherits the same V1/V2/V3 and distributional clearances and is
    # fed by the deconjugation flux out of the T-DM1 central compartment.
    d/dt(central)           <- -(cl + cldec) / vc * central -
                                qd2 / vc * central + qd2 / vp  * peripheral1 -
                                qd3 / vc * central + qd3 / vp2 * peripheral2
    d/dt(peripheral1)       <-  qd2 / vc * central - qd2 / vp  * peripheral1
    d/dt(peripheral2)       <-  qd3 / vc * central - qd3 / vp2 * peripheral2

    d/dt(dar0_central)      <-  cldec / vc * central - cl / vc * dar0_central -
                                qd2 / vc * dar0_central + qd2 / vp  * dar0_peripheral1 -
                                qd3 / vc * dar0_central + qd3 / vp2 * dar0_peripheral2
    d/dt(dar0_peripheral1)  <-  qd2 / vc * dar0_central - qd2 / vp  * dar0_peripheral1
    d/dt(dar0_peripheral2)  <-  qd3 / vc * dar0_central - qd3 / vp2 * dar0_peripheral2

    # Observations: T-DM1 = conjugated ADC central concentration;
    # total trastuzumab (TT) = T-DM1 + DAR0 central concentration
    Cc  <- central / vc
    Ctt <- (central + dar0_central) / vc

    Cc  ~ prop(CcpropSd)
    Ctt ~ prop(CttpropSd)
  })
}
