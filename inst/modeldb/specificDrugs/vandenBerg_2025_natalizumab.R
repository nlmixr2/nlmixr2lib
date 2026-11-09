vandenBerg_2025_natalizumab <- function() {
  description <- "Two-compartment population PK model with parallel linear and Michaelis-Menten elimination for intravenous natalizumab in adults with relapsing-remitting multiple sclerosis on standard or extended-interval maintenance dosing, with allometric weight effects on linear clearance and central volume (van den Berg 2025)"
  reference <- "van den Berg SPH, Toorop AA, Hooijberg F, Wolbink G, Voelkner NMF, Gelissen LMY, Killestein J, van Kempen ZLE, Dorlo TPC, Rispens T. Pharmacokinetic Model-Informed Precision Dosing of Natalizumab in Multiple Sclerosis. CPT Pharmacometrics Syst Pharmacol. 2025;14(6):1032-1041. doi:10.1002/psp4.70014"
  vignette <- "vandenBerg_2025_natalizumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline body weight; power-function effects on linear CL and central volume, both normalised to 70 kg (Supplement Section 6 NONMEM $PK, '(WT/70)**THETA(7)' and '(WT/70)**THETA(8)').",
      source_name = "WT"
    )
  )

  compartmentData <- list(
    central = list(analyte = "natalizumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "natalizumab", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 379,
    n_studies = 3,
    n_observations = 2523,
    age_range = "IQR 32-47 years",
    age_median = "38 years",
    weight_range = "IQR 64-83 kg",
    weight_median = "73 kg",
    sex_female_pct = 78.6,
    race_ethnicity = "Not reported",
    disease_state = "Relapsing-remitting multiple sclerosis in the maintenance phase of intravenous natalizumab treatment (at least 6 months of treatment in the extended-interval studies); patients with anti-drug antibodies excluded",
    dose_range = "300 mg IV every 4 weeks (standard interval) or at personalised extended intervals (dose interval range 21-64 days, median 35 days)",
    regions = "Netherlands (20 clinics)",
    notes = "Table 1 of van den Berg 2025 (median and IQR). Pooled from the STDMS biobank cohort (standard 4-week interval), the PDNMS trial (NCT03516526) and the NEXTMS trial (NCT04225312) main (10 ug/mL trough target) and low (5 ug/mL target) arms (Supplement Table S1). 2489 trough samples and 34 peak samples drawn 1 h after infusion; serum natalizumab by ELISA, LLOQ 0.025 ug/mL."
  )

  ini({
    # Final-model estimates: Table 2 (main text) and Table S4; model form from
    # the Supplement Section 6 final NONMEM control stream.
    lcl <- log(0.164); label("Linear clearance at 70 kg (L/day)") # Table 2 'Clearance (CL, L/d)' = 0.164; THETA(1)
    lvc <- log(3.48); label("Central volume of distribution at 70 kg (L)") # Table 2 'Central Volume of distribution (Vc, L)' = 3.48; THETA(2)
    lvp <- log(3.48); label("Peripheral volume of distribution (L)") # Table 2 'Peripheral Volume of distribution (Vp, L)' = 3.48; THETA(3)
    lq <- fixed(log(0.4)); label("Intercompartmental clearance (L/day)") # Table 2 'Intercompartmental clearance (Q, L/d)' = 0.4 (Fixed); Section 6 'Q = 0.4'
    lvmax <- log(2.86); label("Maximum elimination rate of the saturable pathway (mg/day)") # Table 2 'Maximum elimination rate (Vmax, mg/d)' = 2.86; THETA(4)
    lkm <- log(0.875); label("Concentration at half-maximal saturable elimination (mg/L)") # Table 2 'Concentration at half Vmax (Km, mg/L)' = 0.875; THETA(5)

    e_wt_cl <- 0.6; label("Power exponent of body weight on linear clearance (unitless)") # Table 2 'WT on CL' = 0.6; THETA(7)
    e_wt_vc <- 0.26; label("Power exponent of body weight on central volume (unitless)") # Table 2 'WT on Vc' = 0.26; THETA(8)

    # IIV reported as CV%; omega^2 = log(CV^2 + 1) (Supplement Section 2, Mould
    # and Upton). No IIV on Vc, Q or Km (Section 6 $PK).
    etalcl ~ 0.050678 # Table 2 IIV CL = 22.8 CV%
    etalvp ~ 0.094605 # Table 2 IIV Vp = 31.5 CV%
    etalvmax ~ 0.065413 # Table 2 IIV Vmax = 26 CV%

    propSd <- 0.22; label("Proportional residual error (fraction)") # Table 2 'Proportional error' = 22 CV%; Section 6 'W = SQRT(THETA(6)**2*IPRED**2)'
  })
  model({
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp + etalvp)
    q <- exp(lq)
    vmax <- exp(lvmax + etalvmax)
    km <- exp(lkm)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    Cc <- central / vc

    d/dt(central) <- -vmax * Cc / (km + Cc) - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc ~ prop(propSd)
  })
}
