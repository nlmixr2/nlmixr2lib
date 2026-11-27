Li_2016_citrate <- function() {
  description <- "One-compartment endogenous turnover model for plasma citrate in healthy adults (Li 2016), fitted to literature mean profiles after intravenous citrate infusion and oral citrate. Citrate is produced at a constant zero-order rate that holds plasma citrate at its baseline, distributes in a single body-weight-normalised volume, and is cleared by first-order elimination; oral citrate is absorbed first-order with bioavailability F. Hepatic clearance (Eq. 4) and the hepatic influx of plasma citrate (Eq. 5) are derived outputs, and a scenario knob f_nact scales the hepatic (NaCT-mediated) clearance to reproduce the NaCT-inhibition simulations of Eqs. 6-7 and Figure 6."
  reference <- "Li Z, Erion DM, Maurer TS. Model-Based Assessment of Plasma Citrate Flux Into the Liver: Implications for NaCT as a Therapeutic Target. CPT Pharmacometrics Syst Pharmacol. 2016;5(3):132. doi:10.1002/psp4.12062"
  vignette <- "Li_2016_citrate"
  dosing <- c("depot", "central")
  units <- list(time = "min", dosing = "mmol/kg", concentration = "mmol/L")

  compartmentData <- list(
    depot = list(analyte = "citrate", units = "mmol/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "citrate", units = "mmol/kg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 62L,
    n_studies = 5L,
    disease_state = "Healthy adults (literature control groups).",
    dose_range = "Intravenous trisodium citrate infusion 0.33-0.5 mmol/kg/h for 120-220 min (total 0.68-1.21 mmol/kg); single oral citrate 0.19-0.26 mmol/kg.",
    notes = "Model fitted to published MEAN plasma citrate profiles from five literature studies (Kramer 2003 n = 16, Apsner 1997 n = 6, Bauer 2005 n = 11, Rudman 1980 n = 22, Fegan 1992 n = 7; sizes from Supporting Information II), weighted by sqrt(n)/SD (Li 2016 Methods). The NONMEM control stream (Supporting Information I) places one random effect on CL and one on baseline citrate, with one record per STUDY, so those etas are between-study variability; their final variances are not reported, so this file carries the typical values only. The residual was the literature standard error of the mean (SIGMA fixed to 1, W = SD/sqrt(n)), not an estimated error model, so no residual error is encoded. Doses and parameters are per kg body weight; the paper's g/day flux figures correspond to a 70 kg adult and citric acid MW 192.12 g/mol."
  )

  ini({
    lcl <- log(7.79); label("Plasma citrate elimination clearance CLelim (mL/min/kg)") # Table 1 'CLelim' = 7.79 mL/min/kg [5.70-10.45]
    lvc <- log(259); label("Central distribution volume Vd (mL/kg)") # Table 1 'Vd' = 259 mL/kg [242-281]
    lka <- log(0.02); label("First-order oral absorption rate constant kabs (1/min)") # Table 1 'kabs' = 0.02 per min [0.01-0.04]
    lfdepot <- log(0.37); label("Oral citrate bioavailability F (fraction)") # Table 1 'F' = 0.37 [0.20-0.64]
    lbl_citrate <- log(0.13); label("Baseline plasma citrate Cp(0) (mmol/L)") # Table 1 'Cp(0)' = 0.13 mM [0.11-0.14]

    qh <- fixed(20); label("Hepatic blood flow Qh (mL/min/kg)") # Methods, Eq. 4 text: 'Qh ... is 20 mL/min/kg' (ref 22)
    f_plasma <- fixed(0.55); label("Volume fraction of plasma in whole blood fp (fraction)") # Methods, Eq. 4 text: 'fp is 55%'
    fa <- fixed(0.97); label("Fraction of oral citrate absorbed (fraction)") # Methods, Eq. 4 text: 'citrate absorption of 97%' (ref 15, Fegan 1992)

    f_nact <- fixed(1); label("Fraction of baseline NaCT (hepatic citrate uptake) activity remaining, 1 - fiNaCT (fraction)") # Discussion Eqs. 6-7 and Figure 6 scenario knob; 1 = no inhibition
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    ka <- exp(lka)
    fdepot <- exp(lfdepot)
    bl_citrate <- exp(lbl_citrate)

    # Eq. 4: hepatic clearance of plasma citrate, CLh = Qh * fp * (0.97 - F)
    # (6.60 mL/min/kg, ~85% of CLelim, at the typical values).
    clh <- qh * f_plasma * (fa - fdepot)
    # Non-hepatic remainder of the whole-body clearance.
    clnh <- cl - clh
    # Discussion Eqs. 6-7: NaCT is assumed to carry all of the hepatic clearance,
    # so inhibition by fiNaCT = 1 - f_nact scales CLh only.
    clh_nact <- clh * f_nact
    cl_tot <- clnh + clh_nact

    # Zero-order endogenous production set by the uninhibited steady state,
    # k_p = CL * Cp(0) (Supporting Information I, $PK k_p = CL1*Bo/1000).
    # The 1/1000 converts mL/kg * mmol/L to mmol/kg.
    kprod <- cl * bl_citrate / 1000

    # Eq. 1
    d/dt(depot) <- -ka * depot
    # Eq. 2, written in amounts (mmol/kg) as in Supporting Information I $DES
    d/dt(central) <- kprod - cl_tot / vc * central + ka * depot

    central(0) <- bl_citrate * vc / 1000
    f(depot) <- fdepot

    # Plasma citrate (mmol/L = mM)
    Cc <- central / vc * 1000

    # Eq. 5: hepatic influx of plasma citrate, CLh * Cp (mmol/min/kg)
    flux_hep <- clh_nact * Cc / 1000
  })
}
