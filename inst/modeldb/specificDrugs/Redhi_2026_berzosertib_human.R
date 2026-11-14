Redhi_2026_berzosertib_human <- function() {
  description <- paste(
    "Human antimalarial exposure-response model for intravenous berzosertib (VE-822 /",
    "VX-970 / M6620), an ATR kinase inhibitor assessed for repurposing against",
    "Plasmodium falciparum malaria. Two-compartment population PK model of the UNBOUND",
    "plasma concentration after 1-h IV infusions, fitted to concentration-time",
    "profiles simulated from the authors' PK-Sim PBPK model (20 virtual adults), with",
    "fixed allometric weight exponents (1 on volumes, 0.75 on clearances). It drives the",
    "parasitemia model of the companion mouse fit, transferred unchanged with its",
    "unbound-plasma EC50, from a 1% starting parasitemia. The PD layer carries a known",
    "source inconsistency: with the published growth and kill rate constants the",
    "maximal net parasite decline is 0.020 /h, whereas the paper's Figure 4 regimens",
    "clear two logs of parasitemia over about 140 h, implying roughly 0.03 /h; the",
    "values are carried verbatim, not tuned (see the vignette).",
    "See modellib('Redhi_2026_berzosertib_mouse') for the PD source."
  )
  reference <- paste(
    "Redhi D, Mulubwa M, Mmonwa MM, Njoroge M, Krugmann L, Chibale K.",
    "Repurposing Drugs for Malaria through a Human Dose Prediction:",
    "A Case Study with Berzosertib. Clin Pharmacol Ther. 2026;120(5):1315-1325.",
    "doi:10.1002/cpt.70303. PK from Supplementary Table S6 and Supplementary",
    "Material S1 section 1.12; PD from Table 1."
  )
  vignette <- "Redhi_2026_berzosertib"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL",
    parasites = "percent parasitemia"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling applied to CL, Vc, Q and Vp with fixed exponents",
        "(0.75 on CL and Q, 1 on Vc and Vp), centred on the 72.4 kg mean of the",
        "20 simulated adults (Supplementary Material S1 section 1.12). The Table S6",
        "estimates are reported per kg of body weight; multiplying by WT gives the",
        "absolute parameter for an individual."
      ),
      source_name = "Body weight"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      notes = "Tested and not retained in the human population PK model (Supplementary Material S1 section 1.12)."
    ),
    SEXF = list(
      description = "Sex indicator (1 = female).",
      units = "(binary)",
      type = "binary",
      notes = "Tested and not retained (Supplementary Material S1 section 1.12). The simulated cohort was 10 male and 10 female."
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "berzosertib (unbound)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "berzosertib (unbound)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    parasites = list(
      analyte = "Plasmodium falciparum asexual blood stage",
      units = "% parasitemia",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 20L,
    n_studies = 1L,
    weight_range = "48.9-95.8 kg",
    weight_median = "72.4 kg (mean)",
    sex_female_pct = 50,
    disease_state = "Simulated uncomplicated Plasmodium falciparum blood-stage malaria (exposure-response); PK fitted to PBPK-simulated cancer-population profiles.",
    dose_range = paste(
      "Simulated 1-h IV infusions. PBPK/PK development: 18, 90 and 210 mg/m^2.",
      "Antimalarial dose prediction: single 440 mg; 145 mg q48h x 3; 75 mg q24h x 6;",
      "35 mg q12h x 12. Proposed efficacious regimen 75 mg q24h x 6."
    ),
    regions = "Virtual population matched to the Middleton and Yap phase-I cancer studies; PK analysis in Monolix 2024R1 (Supplementary Material S1 section 1.12).",
    notes = paste(
      "The human PK parameters were estimated from unbound-plasma concentration-time",
      "profiles simulated by a PK-Sim PBPK model, not from measured clinical PK.",
      "PK samples at 0.5, 1, 1.5, 2, 4, 8, 24, 48 and 72 h. The unbound EC50 and MPC",
      "(0.58 and 0.79 ng/mL) and the parasite growth and kill rate constants come from",
      "the P. falciparum-infected humanized-mouse PK-PD fit (Table 1), assuming the",
      "free drug acts on the parasite identically in human and mouse blood. Exposure-",
      "response simulations started from 1% parasitemia with a 0.01% limit of detection."
    )
  )

  ini({
    # Unbound-plasma PK: Supplementary Table S6 (population PK on the PBPK-
    # simulated data). Estimates are per kg; lX is log of the per-kg value.
    lcl <- log(10.7)
    label("Unbound clearance CL per kg (L/h/kg)") # Table S6: Cl = 10.7 L/h/kg (RSE 5.9%)
    lvc <- log(26.3)
    label("Unbound central volume Vc per kg (L/kg)") # Table S6: V c = 26.3 L/kg (RSE 5.3%)
    lq <- log(29.03)
    label("Unbound intercompartmental clearance Q per kg (L/h/kg)") # Table S6: Q = 29.03 L/h/kg (RSE 7.9%)
    lvp <- log(167.3)
    label("Unbound peripheral volume Vp per kg (L/kg)") # Table S6: V p = 167.3 L/kg (RSE 9.7%)

    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of weight on CL and Q (unitless)") # S1 section 1.12: exponent 0.75 fixed on Cl and Q
    e_wt_vc <- fixed(1)
    label("Allometric exponent of weight on Vc and Vp (unitless)") # S1 section 1.12: exponent 1 fixed on V c and V p

    # Parasitemia PD transferred from the mouse fit (Table 1), driven by the
    # UNBOUND plasma concentration with the unbound-plasma EC50.
    lkgrow <- fixed(log(0.038))
    label("Parasite growth rate constant Kgrowth (1/h)") # Table 1: K growth = 0.038 /h
    lemax <- fixed(log(0.058))
    label("Maximal parasite kill rate constant Kkill (1/h)") # Table 1: K kill = 0.058 /h
    lec50 <- fixed(log(0.58))
    label("Unbound-plasma EC50 of parasite killing (ng/mL)") # Table 1: unbound plasma EC50 = 0.58 ng/mL
    lhill <- fixed(log(7))
    label("Hill coefficient H (unitless)") # Table 1: H = 7, fixed
    lrbase_parasites <- fixed(log(1))
    label("Initial parasitemia at the start of treatment (%)") # Methods 'Human exposure-response': initial parasitemia 1%

    # IIV, Table S6 (SD of the log-normal random effects; variances are squares).
    etalcl ~ 0.0625 # Table S6: omega Cl = 0.25
    etalvc ~ 0.0121 # Table S6: omega Vc = 0.11
    etalq ~ 0.0784 # Table S6: omega Q = 0.28
    etalvp ~ 0.1764 # Table S6: omega Vp = 0.42

    propSd <- 0.2
    label("Proportional residual error, unbound berzosertib (fraction)") # Table S6: proportional = 0.2 (RSE 7.0%)
  })

  model({
    # Allometric scaling on absolute parameters, centred on 72.4 kg.
    cl <- exp(lcl + etalcl) * WT * (WT / 72.4)^(e_wt_cl - 1)
    vc <- exp(lvc + etalvc) * WT * (WT / 72.4)^(e_wt_vc - 1)
    q <- exp(lq + etalq) * WT * (WT / 72.4)^(e_wt_cl - 1)
    vp <- exp(lvp + etalvp) * WT * (WT / 72.4)^(e_wt_vc - 1)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- k21 * peripheral1 - k12 * central - kel * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Unbound plasma concentration (ng/mL); dose in mg, V in L.
    Cc <- 1000 * central / vc

    kgrow <- exp(lkgrow)
    kkill <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # Equation 2: direct-effect sigmoidal kill on parasitemia (percent).
    effect <- Cc^hill / (Cc^hill + ec50^hill)
    d/dt(parasites) <- kgrow * parasites - kkill * effect * parasites
    parasites(0) <- exp(lrbase_parasites)

    # Equation 3: minimum parasiticidal concentration, the unbound EC90 (ng/mL).
    mpc <- ec50 * 9^(1 / hill)

    Cc ~ prop(propSd)
  })
}
