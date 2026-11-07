Staudt_2023_voriconazole_rat <- function() {
  description <- paste(
    "Preclinical (rat, Wistar and Sprague-Dawley male). Four-compartment",
    "population PK model for intravenous voriconazole in healthy rats and",
    "in rats with Cryptococcus neoformans meningoencephalitis, fit jointly to",
    "total plasma, unbound brain interstitial fluid (microdialysis) and CSF",
    "concentrations pooled from three studies. The ODE states carry UNBOUND",
    "drug amounts: central (plasma) with Michaelis-Menten elimination, a",
    "peripheral compartment, the brain extracellular fluid (volume fixed to",
    "the 0.41 mL physiological value) linked to central by influx and efflux",
    "clearances, and a CSF compartment exchanging with the brain",
    "extracellular fluid. Total plasma is unbound plasma divided by the",
    "unbound fraction 0.35. Infection has its own Km and brain efflux",
    "clearance. All parameters are scaled allometrically by body weight",
    "(exponent 0.75 on Vmax and Km, 1 on the rest).",
    sep = " "
  )
  reference <- paste(
    "Staudt KJ, Dias BB, Alves IA, Lelievre B, Bouchara JP, de Araujo BV.",
    "Modeling and Simulation as a Tool to Assess Voriconazole Exposure in the",
    "Central Nervous System. Pharmaceutics. 2023;15(7):1781.",
    "doi:10.3390/pharmaceutics15071781. PMCID: PMC10384042.",
    "Parameter estimates: Table 1. Structural ODEs: Equations 4-7.",
    "Allometric scaling: Equation 1. Unbound-fraction observation equation:",
    "Results text after Equation 7. Group sizes: Supplementary Table S2.",
    "Human-to-rat dose scaling: Supplementary Equations S1-S2 and Table S1.",
    sep = " "
  )
  vignette <- "Staudt_2023_voriconazole_rat"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "voriconazole (unbound)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "voriconazole (unbound)", units = "mg", specimen = "tissue", verified = TRUE),
    brain_ecf = list(analyte = "voriconazole (unbound)", units = "mg", specimen = "brain ISF", verified = TRUE),
    brain_csf = list(analyte = "voriconazole", units = "mg", specimen = "CSF", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      source_name = "BW",
      notes = paste(
        "Equation 1: every PK parameter is scaled by (BW / BWmedian)^k with",
        "k = 0.75 for Vmax and Km and k = 1 for all other parameters, to",
        "absorb the weight difference between the Wistar and Sprague-Dawley",
        "rats of the pooled studies. The paper does not report BWmedian; the",
        "model uses 0.35 kg, the rat body mass the same paper uses for its",
        "allometric dose scaling (Supplementary Equations S1-S2). At WT =",
        "0.35 kg the parameters equal the Table 1 estimates."
      )
    ),
    DIS_INFECT_ACTIVE = list(
      description = "Experimental Cryptococcus neoformans infection: 1 = rat infected with C. neoformans (meningoencephalitis model); 0 = healthy rat.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy)",
      source_name = "Infected status",
      notes = paste(
        "Experimental rather than clinical criterion: the infected animals",
        "are the C. neoformans-infected male Wistar rats of Alves et al. 2017",
        "(Antimicrob Agents Chemother, doi:10.1128/AAC.00321-17), studied",
        "with brain microdialysis after a 5 mg/kg intravenous dose.",
        "Results: 'The parameters Km and Qout were calculated separately for",
        "the infected group', so the indicator selects the infected Km",
        "(8.135 mg) and brain efflux clearance Qout (0.388 L/h) in place of",
        "the healthy values (4.76 mg, 0.634 L/h)."
      )
    )
  )

  population <- list(
    species = "rat (Wistar and Sprague-Dawley, male)",
    n_subjects = 49L,
    n_studies = 3L,
    age_range = "Adult (age not reported)",
    weight_range = "Not reported (allometric reference 0.35 kg assumed)",
    sex_female_pct = 0,
    disease_state = paste(
      "Healthy rats, and rats infected with Cryptococcus neoformans",
      "(experimental cryptococcal meningoencephalitis; Alves et al. 2017",
      "study only)."
    ),
    dose_range = "Single 5, 10 or 30 mg/kg intravenous bolus.",
    regions = "Brazil (Federal University of Rio Grande do Sul, Porto Alegre); France (Angers)",
    notes = paste(
      "Data pooled from three preclinical studies (Section 2.1;",
      "Supplementary Table S2, animals/observations): Alves et al. 2017,",
      "Wistar, 5 mg/kg, plasma healthy 11/100 and infected 5/49, brain",
      "microdialysis healthy 11/183 and infected 5/101; Li 2008, Wistar,",
      "plasma 5 mg/kg 6/60 and 10 mg/kg 6/60; Lelievre et al. 2018,",
      "Sprague-Dawley, 30 mg/kg, 21 rats with 18 total plasma and 18 CSF",
      "observations (three animals sacrificed per time point, analysed as",
      "naive pooled data). Results text: 287 plasma, 284 microdialysate and",
      "18 CSF concentrations. Plasma protein binding 66% (Methods 2.2).",
      "NONMEM 7.4.3, FOCE-I, ADVAN13."
    )
  )

  ini({
    # Structural parameters: Table 1 (estimate, %RSE), for a rat at the
    # allometric reference weight. Units as printed in Table 1 except where
    # noted; the abstract prints Q1 'L' and Qin / Qout 'h-1', the table 'L'
    # for Q1 and 'L/h' for Qin / Qout -- all three are clearances (L/h).
    lvmax <- log(1.06); label("Maximum elimination rate Vmax (mg/h)") # Table 1: Vmax = 1.06 mg/h (RSE 13%)
    # Equation 4 writes the elimination as Vmax * A1 / (Km + A1), with the
    # unbound central AMOUNT A1, so Km is encoded as an amount (mg) although
    # Table 1 labels it mg/L; see the vignette for the reasoning.
    lkm <- log(4.76); label("Michaelis-Menten constant Km, healthy rats (mg)") # Table 1: Km = 4.76 mg/L (RSE 20%); abstract and Results text: Km = 4.76
    lkm_infected <- log(8.135); label("Michaelis-Menten constant Km, C. neoformans-infected rats (mg)") # Results after Eq 7 and Discussion: Km,infected = 8.135; Table 1: Km,infected = 8.13 (RSE 44%)
    lvc <- log(1.15); label("Central volume of distribution V1, unbound (L)") # Table 1: V1 = 1.15 L (RSE 7%)
    lq <- log(2.69); label("Intercompartmental clearance central-peripheral Q1 (L/h)") # Table 1: Q1 = 2.69 (RSE 18%), unit printed as 'L'
    lvp <- log(0.37); label("Peripheral volume of distribution V2 (L)") # Table 1: V2 = 0.37 L (RSE 18%)
    lclin <- log(0.81); label("Influx clearance central to brain extracellular fluid Qin (L/h)") # Table 1: Qin = 0.81 L/h (RSE 13%)
    lclef <- log(0.634); label("Efflux clearance brain extracellular fluid to central Qout, healthy rats (L/h)") # Results after Eq 7 and Discussion: Qout = 0.634; Table 1: Qout = 0.63 L/h (RSE 13%)
    lclef_infected <- log(0.388); label("Efflux clearance brain extracellular fluid to central Qout, C. neoformans-infected rats (L/h)") # Results after Eq 7 and Discussion: Qout,infected = 0.388; Table 1: Qout,infected = 0.38 L/h (RSE 36%)
    lvecf <- fixed(log(0.00041)); label("Brain extracellular fluid volume V3 (L)") # Table 1: V3 = 0.00041 L FIX; Results: 'the interstitial physiological value [36]'
    lqcsf <- log(0.33); label("Exchange clearance brain extracellular fluid to CSF Q4 (L/h)") # Table 1: Q4 = 0.33 L/h (RSE 19%)
    lvcsf <- log(0.13); label("CSF volume of distribution V4 (L)") # Table 1: V4 = 0.13 L (RSE 17%)
    fu <- fixed(0.35); label("Unbound fraction of voriconazole in rat plasma (fraction)") # Results after Eq 7: estimates 'divided by 0.35, the unbound fraction of the drug [22]'

    # Allometric exponents (Equation 1): 'k = 0.75 for maximum rate of
    # metabolism (Vmax) and Michaelis-Menten constant (Km), and k = 1 for all
    # the other parameters'.
    e_wt_vmax_km <- fixed(0.75); label("Allometric exponent on Vmax and Km (unitless)") # Methods 2.2, after Eq 1
    e_wt_q_vc <- fixed(1); label("Allometric exponent on every clearance and volume (unitless)") # Methods 2.2, after Eq 1

    # IIV: Table 1 reports %CV for exponential IIV; omega^2 = log(CV^2 + 1).
    # The infected Km and Qout share the healthy eta (one eta per parameter).
    etalkm ~ 0.1349 # Table 1: omega Km = 38 %CV -> log(0.38^2 + 1); shrinkage 5.8%
    etalvc ~ 0.0807 # Table 1: omega V1 = 29 %CV -> log(0.29^2 + 1); shrinkage 5.5%
    etalvp ~ 0.1980 # Table 1: omega V2 = 47 %CV -> log(0.47^2 + 1); shrinkage 26.7%
    etalclef ~ 0.1418 # Table 1: omega Qout = 39 %CV -> log(0.39^2 + 1); shrinkage 20.1%

    # Residual error: log-additive, 'described separately for plasma,
    # microdialysate and CSF data' (Methods 2.2), encoded as SDs on the log
    # scale.
    expSd <- 0.107; label("Log-additive residual SD, total plasma (log scale)") # Table 1: plasma log-additive error = 0.107 (RSE 5%)
    expSd_Cbrain_ecf <- 0.074; label("Log-additive residual SD, brain microdialysate (log scale)") # Table 1: microdialysis log-additive error = 0.074 (RSE 5%)
    expSd_Cbrain_csf <- 0.450; label("Log-additive residual SD, CSF (log scale)") # Table 1: CSF log-additive error = 0.450 (RSE 18%)
  })

  model({
    # Allometric reference weight; the paper's BWmedian is not reported.
    wt_ref <- 0.35 # kg; rat body mass of Supplementary Equations S1-S2
    allo_mm <- (WT / wt_ref)^e_wt_vmax_km
    allo_lin <- (WT / wt_ref)^e_wt_q_vc

    # Results: 'The parameters Km and Qout were calculated separately for
    # the infected group' (Equations 8-9 as printed are not usable: Eq 8
    # gives Km = 0 for a healthy rat).
    vmax <- exp(lvmax) * allo_mm
    km <- exp(lkm * (1 - DIS_INFECT_ACTIVE) + lkm_infected * DIS_INFECT_ACTIVE + etalkm) * allo_mm
    vc <- exp(lvc + etalvc) * allo_lin
    q <- exp(lq) * allo_lin
    vp <- exp(lvp + etalvp) * allo_lin
    clin <- exp(lclin) * allo_lin
    clef <- exp(lclef * (1 - DIS_INFECT_ACTIVE) + lclef_infected * DIS_INFECT_ACTIVE + etalclef) * allo_lin
    vecf <- exp(lvecf) * allo_lin
    qcsf <- exp(lqcsf) * allo_lin
    vcsf <- exp(lvcsf) * allo_lin

    # Equations 4-7 (unbound amounts). A1 = central, A2 = peripheral1,
    # A3 = brain_ecf, A4 = brain_csf.
    d/dt(central) <- -q / vc * central - vmax * central / (km + central) + q / vp * peripheral1 + clef / vecf * brain_ecf - clin / vc * central
    d/dt(peripheral1) <- -q / vp * peripheral1 + q / vc * central
    d/dt(brain_ecf) <- -clef / vecf * brain_ecf + clin / vc * central + qcsf / vcsf * brain_csf - qcsf / vecf * brain_ecf
    d/dt(brain_csf) <- -qcsf / vcsf * brain_csf + qcsf / vecf * brain_ecf

    # Unbound plasma, unbound brain extracellular fluid (microdialysis) and
    # CSF. The CSF observations were total CSF concentrations fitted
    # directly to A4 / V4.
    Cu <- central / vc
    Cbrain_ecf <- brain_ecf / vecf
    Cbrain_csf <- brain_csf / vcsf
    # Observed plasma concentrations are total: 'the estimates were divided
    # by 0.35, the unbound fraction of the drug' (Results after Eq 7).
    Cc <- Cu / fu

    Cc ~ lnorm(expSd)
    Cbrain_ecf ~ lnorm(expSd_Cbrain_ecf)
    Cbrain_csf ~ lnorm(expSd_Cbrain_csf)
  })
}
