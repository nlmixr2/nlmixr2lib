Bender_2014_trastuzumabEmtansine_mechanistic <- function() {
  description <- "Mechanistic DAR0-DAR7 catenary deconjugation PK model for trastuzumab emtansine (T-DM1) in cynomolgus monkeys (default) and rats (Bender 2014): each DAR moiety distributes into a shared three-compartment backbone and deconjugates sequentially toward naked trastuzumab (DAR0); uses five shared upper-chain rate constants (k7->3) plus separate k_2->1 and k_1->0."
  reference <- "Bender B, Leipold DD, Xu K, Shen BQ, Tibbitts J, Friberg LE. A mechanistic pharmacokinetic model elucidating the disposition of trastuzumab emtansine (T-DM1), an antibody-drug conjugate (ADC) for treatment of metastatic breast cancer. AAPS J. 2014;16(5):994-1008. doi:10.1208/s12248-014-9618-3"
  vignette <- "Bender_2014_trastuzumabEmtansine_mechanistic"
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
    scope_note     = "Preclinical-only (no human PK). Companion to the simpler Bender_2014_trastuzumabEmtansine_reduced model file (Bender 2014 Table III).",
    model_variant  = "Mechanistic DAR0-DAR7 model (Bender 2014 Table II). Each DAR species has its own three-compartment distribution sharing V1/V2/V3/CLd2/CLd3 with all other DAR species; a catenary chain of first-order deconjugation rates drives DAR_n -> DAR_(n-1).",
    dose_products  = list(
      note         = "Paper reports two preclinical dose products. The DAR-specific fractions below are seeded into the DAR0_central..DAR7_central compartments at t=0 when simulating. Default vignette uses DAR 3.1 (cynomolgus studies).",
      T_DM1_DAR3_1 = list(
        DAR7 = 0.02, DAR6 = 0.05, DAR5 = 0.10, DAR4 = 0.19,
        DAR3 = 0.26, DAR2 = 0.23, DAR1 = 0.13, DAR0 = 0.02
      ),
      T_DM1_DAR1_5 = list(
        DAR7 = 0.00, DAR6 = 0.00, DAR5 = 0.01, DAR4 = 0.04,
        DAR3 = 0.13, DAR2 = 0.26, DAR1 = 0.35, DAR0 = 0.21
      )
    ),
    species_parameters = list(
      cynomolgus = list(
        source   = "Bender 2014 Table II, cynomolgus columns",
        CL_TT    = "17.4 mL/day (IIV CV 24.8%)",
        k_plasma = "0.0939 /day",
        V1       = "148 mL (IIV CV 11.7%)",
        CLd2     = "25.5 mL/day (no IIV)",
        V2       = "57.2 mL (IIV CV 46.8%)",
        CLd3     = "81.2 mL/day (no IIV)",
        V3       = "127 mL (no IIV)",
        k_7_to_3 = "0.341 /day (no IIV)  [shared k_7->6 = k_6->5 = k_5->4 = k_4->3 = k_3->2]",
        k_2_to_1 = "0.255 /day (no IIV)",
        k_1_to_0 = "0.0939 /day (no IIV)",
        res_err  = "15.3% proportional"
      ),
      rat = list(
        source   = "Bender 2014 Table II, rat columns",
        CL_TT    = "2.42 mL/day (IIV CV 24.0%)",
        k_plasma = "0.156 /day",
        V1       = "11.0 mL (IIV CV 18.5%)",
        CLd2     = "49.0 mL/day (no IIV)",
        V2       = "3.44 mL (IIV CV 49.4%)",
        CLd3     = "12.0 mL/day (no IIV)",
        V3       = "16.7 mL (IIV CV 16.8%)",
        k_7_to_3 = "0.543 /day (IIV CV 21.8%)",
        k_2_to_1 = "0.388 /day (no IIV)",
        k_1_to_0 = "0.114 /day (IIV CV 15.1%)",
        res_err  = "11.1% proportional",
        notes    = "Rat fit has additional IIV terms on V3, k_7->3, and k_1->0 not active in the cynomolgus default. To switch to rat, update ini() fixed-effect values (mL/day -> L/day and mL -> L by /1000)."
      )
    ),
    notes          = "NONMEM 7.2, FOCE-INTER; residual error reported as additive on log-scale (proportional in linear space). CL_TT is a fit parameter composed of CL_in_vivo (in vivo antibody clearance) and k_plasma * V1 (plasma degradation, constrained by in vitro plasma stability data). In the in vivo ODE the sum CL_in_vivo/V1 + k_plasma reduces to CL_TT/V1, so this model uses CL_TT directly; k_plasma is retained as a reported parameter because it is needed to simulate the in vitro plasma-stability experiments (set CL_in_vivo = CL_TT - k_plasma * V1 and zero the distributional clearances)."
  )

  ini({
    # Structural parameters - cynomolgus monkey (default ini values)
    # Volumes stored in L (mL / 1000) and clearances in L/day so that dose in
    # mg and V in L give concentrations in mg/L = ug/mL directly.
    lcl      <- log(0.0174);  label("Total trastuzumab clearance (CL_TT, L/day)")           # Bender 2014 Table II, cyno: 17.4 mL/day
    lkplasma <- log(0.0939);  label("Plasma antibody degradation rate constant (k_plasma, 1/day)") # Bender 2014 Table II, cyno: 0.0939 /day
    lvc      <- log(0.148);   label("Central volume shared by all DAR species (V1, L)")      # Bender 2014 Table II, cyno: 148 mL
    lqd2     <- log(0.0255);  label("Distributional clearance to peripheral1 (CLd2, L/day)") # Bender 2014 Table II, cyno: 25.5 mL/day
    lvp      <- log(0.0572);  label("Peripheral volume 1 shared by all DAR species (V2, L)") # Bender 2014 Table II, cyno: 57.2 mL
    lqd3     <- log(0.0812);  label("Distributional clearance to peripheral2 (CLd3, L/day)") # Bender 2014 Table II, cyno: 81.2 mL/day
    lvp2     <- log(0.127);   label("Peripheral volume 2 shared by all DAR species (V3, L)") # Bender 2014 Table II, cyno: 127 mL
    lkhi     <- log(0.341);   label("Shared upper-chain deconjugation rate k_7->6 = k_6->5 = k_5->4 = k_4->3 = k_3->2 (1/day)") # Bender 2014 Table II, cyno: 0.341 /day
    lk21     <- log(0.255);   label("DAR2 deconjugation rate constant k_2->1 (1/day)")       # Bender 2014 Table II, cyno: 0.255 /day
    lk10     <- log(0.0939);  label("DAR1 deconjugation rate constant k_1->0 (1/day)")       # Bender 2014 Table II, cyno: 0.0939 /day

    # IIV - cynomolgus Table II; log-normal, omega^2 = log(CV^2 + 1)
    etalcl ~ 0.05969  # cyno CL_TT CV 24.8% (rat CV 24.0% -> 0.05612)
    etalvc ~ 0.01359  # cyno V1    CV 11.7% (rat CV 18.5% -> 0.03384)
    etalvp ~ 0.19800  # cyno V2    CV 46.8% (rat CV 49.4% -> 0.21770)

    # Residual error - cynomolgus; paper reports a single residual magnitude
    # applied to T-DM1 and total trastuzumab observations
    CcpropSd  <- 0.153; label("Proportional residual error on T-DM1 concentration (fraction)")            # Bender 2014 Table II, cyno: 15.3%
    CttpropSd <- 0.153; label("Proportional residual error on total trastuzumab concentration (fraction)") # Bender 2014 Table II, cyno: 15.3%
  })

  model({
    cl      <- exp(lcl + etalcl)
    kplasma <- exp(lkplasma)
    vc      <- exp(lvc + etalvc)
    qd2     <- exp(lqd2)
    vp      <- exp(lvp + etalvp)
    qd3     <- exp(lqd3)
    vp2     <- exp(lvp2)
    khi     <- exp(lkhi)
    k21     <- exp(lk21)
    k10     <- exp(lk10)

    # In vivo antibody clearance: CL_TT - k_plasma * V1 (Bender 2014 Table II
    # footnote b). In vivo the sum (cl_in_vivo/vc + kplasma) equals cl/vc, so
    # the ODE elimination term simplifies to cl/vc. k_plasma is retained for
    # in vitro plasma-stability simulations.
    cl_in_vivo <- cl - kplasma * vc
    kel        <- cl / vc

    # Three-compartment distribution parameters, shared across all DAR species
    k12 <- qd2 / vc
    k21d <- qd2 / vp
    k13 <- qd3 / vc
    k31d <- qd3 / vp2

    # DAR7: no incoming deconjugation (no DAR8)
    d/dt(dar7_central)     <- -kel * dar7_central - khi * dar7_central -
                               k12 * dar7_central + k21d * dar7_peripheral1 -
                               k13 * dar7_central + k31d * dar7_peripheral2
    d/dt(dar7_peripheral1) <-  k12 * dar7_central - k21d * dar7_peripheral1
    d/dt(dar7_peripheral2) <-  k13 * dar7_central - k31d * dar7_peripheral2

    # DAR6..DAR3: upper chain, shared rate khi both outgoing and incoming
    d/dt(dar6_central)     <-  khi * dar7_central - kel * dar6_central - khi * dar6_central -
                               k12 * dar6_central + k21d * dar6_peripheral1 -
                               k13 * dar6_central + k31d * dar6_peripheral2
    d/dt(dar6_peripheral1) <-  k12 * dar6_central - k21d * dar6_peripheral1
    d/dt(dar6_peripheral2) <-  k13 * dar6_central - k31d * dar6_peripheral2

    d/dt(dar5_central)     <-  khi * dar6_central - kel * dar5_central - khi * dar5_central -
                               k12 * dar5_central + k21d * dar5_peripheral1 -
                               k13 * dar5_central + k31d * dar5_peripheral2
    d/dt(dar5_peripheral1) <-  k12 * dar5_central - k21d * dar5_peripheral1
    d/dt(dar5_peripheral2) <-  k13 * dar5_central - k31d * dar5_peripheral2

    d/dt(dar4_central)     <-  khi * dar5_central - kel * dar4_central - khi * dar4_central -
                               k12 * dar4_central + k21d * dar4_peripheral1 -
                               k13 * dar4_central + k31d * dar4_peripheral2
    d/dt(dar4_peripheral1) <-  k12 * dar4_central - k21d * dar4_peripheral1
    d/dt(dar4_peripheral2) <-  k13 * dar4_central - k31d * dar4_peripheral2

    d/dt(dar3_central)     <-  khi * dar4_central - kel * dar3_central - khi * dar3_central -
                               k12 * dar3_central + k21d * dar3_peripheral1 -
                               k13 * dar3_central + k31d * dar3_peripheral2
    d/dt(dar3_peripheral1) <-  k12 * dar3_central - k21d * dar3_peripheral1
    d/dt(dar3_peripheral2) <-  k13 * dar3_central - k31d * dar3_peripheral2

    # DAR2: receives from DAR3 at khi, loses to DAR1 at k21
    d/dt(dar2_central)     <-  khi * dar3_central - kel * dar2_central - k21 * dar2_central -
                               k12 * dar2_central + k21d * dar2_peripheral1 -
                               k13 * dar2_central + k31d * dar2_peripheral2
    d/dt(dar2_peripheral1) <-  k12 * dar2_central - k21d * dar2_peripheral1
    d/dt(dar2_peripheral2) <-  k13 * dar2_central - k31d * dar2_peripheral2

    # DAR1: receives from DAR2 at k21, loses to DAR0 at k10
    d/dt(dar1_central)     <-  k21 * dar2_central - kel * dar1_central - k10 * dar1_central -
                               k12 * dar1_central + k21d * dar1_peripheral1 -
                               k13 * dar1_central + k31d * dar1_peripheral2
    d/dt(dar1_peripheral1) <-  k12 * dar1_central - k21d * dar1_peripheral1
    d/dt(dar1_peripheral2) <-  k13 * dar1_central - k31d * dar1_peripheral2

    # DAR0: receives from DAR1 at k10, no outgoing deconjugation (terminal)
    d/dt(dar0_central)     <-  k10 * dar1_central - kel * dar0_central -
                               k12 * dar0_central + k21d * dar0_peripheral1 -
                               k13 * dar0_central + k31d * dar0_peripheral2
    d/dt(dar0_peripheral1) <-  k12 * dar0_central - k21d * dar0_peripheral1
    d/dt(dar0_peripheral2) <-  k13 * dar0_central - k31d * dar0_peripheral2

    # Derived per-DAR concentrations (ug/mL)
    Cdar0 <- dar0_central / vc
    Cdar1 <- dar1_central / vc
    Cdar2 <- dar2_central / vc
    Cdar3 <- dar3_central / vc
    Cdar4 <- dar4_central / vc
    Cdar5 <- dar5_central / vc
    Cdar6 <- dar6_central / vc
    Cdar7 <- dar7_central / vc

    # Composite outputs: T-DM1 assay reads the conjugated fraction (DAR1..DAR7);
    # TT (total trastuzumab) assay reads every antibody backbone (DAR0..DAR7).
    Cc  <- (dar1_central + dar2_central + dar3_central + dar4_central +
            dar5_central + dar6_central + dar7_central) / vc
    Ctt <- (dar0_central + dar1_central + dar2_central + dar3_central +
            dar4_central + dar5_central + dar6_central + dar7_central) / vc

    Cc  ~ prop(CcpropSd)
    Ctt ~ prop(CttpropSd)
  })
}
