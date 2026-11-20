Redhi_2026_berzosertib_mouse <- function() {
  description <- paste(
    "Preclinical (mouse). Sequential population PK-PD model of oral berzosertib",
    "(VE-822 / VX-970 / M6620), an ATR kinase inhibitor being assessed for repurposing",
    "as an antimalarial, in NSG mice engrafted with human erythrocytes and infected",
    "with Plasmodium falciparum Pf3D7 0087/N9. One-compartment whole-blood PK with",
    "first-order absorption and linear elimination drives a direct-effect parasitemia",
    "model: exponential parasite growth minus a sigmoidal (Hill = 7, fixed) whole-blood",
    "concentration-dependent first-order kill. Time zero is the first dose, given on",
    "day 3 post-infection. KNOWN SOURCE INCONSISTENCY: with the published growth and",
    "kill rate constants the maximal net parasite decline is 0.020 /h, whereas the",
    "paper's own Figure 4 regimens clear two logs of parasitemia over about 140 h,",
    "implying roughly 0.03 /h; the values are carried verbatim, not tuned",
    "(see the vignette)."
  )
  reference <- paste(
    "Redhi D, Mulubwa M, Mmonwa MM, Njoroge M, Krugmann L, Chibale K.",
    "Repurposing Drugs for Malaria through a Human Dose Prediction:",
    "A Case Study with Berzosertib. Clin Pharmacol Ther. 2026;120(5):1315-1325.",
    "doi:10.1002/cpt.70303. PK from Supplementary Table S4; PD from Table 1."
  )
  vignette <- "Redhi_2026_berzosertib"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL",
    parasites = "percent parasitemia (infected human erythrocytes)"
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      notes = "Explored as a covariate on the PK and PK-PD models but not retained (Results). The mice weighed 25-35 g; doses were given in mg/kg, so the user converts each dose to mg (dose_mg = mg/kg * weight_kg) before simulating."
    ),
    HCT = list(
      description = "Hematocrit.",
      units = "%",
      type = "continuous",
      notes = "Percentage hematocrit was explored as a covariate and not retained (Results)."
    ),
    HUMAN_RBC_PCT = list(
      description = "Percentage of human (engrafted) erythrocytes.",
      units = "%",
      type = "continuous",
      notes = "Percentage human and mouse erythrocytes were explored as covariates and not retained (Methods; Results)."
    ),
    PARA_BL = list(
      description = "Initial parasitemia.",
      units = "%",
      type = "continuous",
      notes = "Explored as a covariate (Methods) and not retained in the final model."
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "berzosertib",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "berzosertib",
      units = "mg",
      specimen = "whole blood",
      verified = TRUE
    ),
    parasites = list(
      analyte = "Plasmodium falciparum Pf3D7 0087/N9 asexual blood stage",
      units = "% parasitemia",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (NOD-scid IL-2Rgamma-null [NSG], engrafted with human erythrocytes)",
    n_subjects = 26L,
    n_studies = 2L,
    age_range = "6-10 weeks",
    weight_range = "25-35 g",
    sex_female_pct = 0,
    disease_state = "Plasmodium falciparum Pf3D7 0087/N9 asexual blood-stage infection (2e7 parasites IV on day 0); treatment from day 3 post-infection.",
    dose_range = paste(
      "Oral gavage. Dose escalation: 0.5, 5, 7.5, 10, 15 and 20 mg/kg once daily",
      "on days 3-6 post-infection. Dose fractionation: 10, 20, 30, 40 and 60 mg/kg",
      "once on day 3, or 20 and 60 mg/kg total split into two equal doses on days 3",
      "and 5."
    ),
    regions = "University of Cape Town (Supplementary Material S1, section 1.10).",
    notes = paste(
      "Two mice per treatment group: 6 dose-escalation groups and 7 dose-fractionation",
      "groups (26 treated mice), plus an untreated vehicle group in each study.",
      "Mice were engrafted daily with human erythrocytes, so the hematocrit during",
      "treatment was about 77%. Whole-blood berzosertib LLOQ 2 ng/mL (censored);",
      "parasitemia limit of detection 0.01% by flow cytometry. Fitted with SAEM in",
      "Monolix 2024R1, PK first and PK-PD second with the individual PK parameters as",
      "regressors."
    )
  )

  ini({
    # Whole-blood PK: Supplementary Table S4 (oral apparent parameters; no IV
    # arm in NSG mice, so CL and V are CL/F and V/F).
    lka <- log(0.44)
    label("Absorption rate constant ka (1/h)") # Table S4: ka = 0.44 /h (RSE 11.7%)
    lcl <- log(0.25)
    label("Apparent whole-blood clearance CL/F (L/h)") # Table S4: Cl blood = 0.25 L/h (RSE 5.73%)
    lvc <- log(1.39)
    label("Apparent whole-blood volume V/F (L)") # Table S4: V d,blood = 1.39 L (RSE 10.5%)

    # Parasitemia PD: Table 1.
    lkgrow <- log(0.038)
    label("Parasite growth rate constant Kgrowth (1/h)") # Table 1: K growth = 0.038 /h (RSE 6.65%)
    lemax <- log(0.058)
    label("Maximal parasite kill rate constant Kkill (1/h)") # Table 1: K kill = 0.058 /h (RSE 8.85%)
    lec50 <- log(13.24)
    label("Whole-blood EC50 of parasite killing (ng/mL)") # Table 1: whole blood EC50,blood = 13.24 ng/mL (RSE 43.7%)
    lhill <- fixed(log(7))
    label("Hill coefficient H (unitless)") # Table 1: H = 7, fixed after a scan of 0.25-9 (Table S5)
    # Not a model parameter in the paper: the PD fit started each mouse at its own
    # observed day-3 parasitemia. 0.6% is the day-3 group mean of Figure 1a
    # (digitised; all groups lie between 0.50% and 0.75%).
    lrbase_parasites <- fixed(log(0.6))
    label("Parasitemia at the first dose, day 3 post-infection (%)") # figure-derived: Figure 1a day-3 group means, 0.50-0.75%

    # IIV: Table S4 and Table 1 report the SD of the log-normal random effects;
    # the variances below are their squares.
    etalcl ~ 0.0441 # Table S4: omega Clb = 0.21
    etalvc ~ 0.1156 # Table S4: omega Vb = 0.34
    etalkgrow ~ 0.01 # Table 1: omega Kgrow = 0.1
    etalemax ~ 0.000144 # Table 1: omega Kkill = 0.012
    etalec50 ~ 0.3844 # Table 1: omega EC50 = 0.62

    propSd <- 0.42
    label("Proportional residual error, whole-blood berzosertib (fraction)") # Table S4: proportional = 0.42 (RSE 5.71%)
    propSd_parasites <- 0.3
    label("Proportional residual error, parasitemia (fraction)") # Table 1: proportional = 0.3 (RSE 7.91%)
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    kel <- cl / vc

    kgrow <- exp(lkgrow + etalkgrow)
    kkill <- exp(lemax + etalemax)
    ec50 <- exp(lec50 + etalec50)
    hill <- exp(lhill)

    # Equations 1a-1c. Dose in mg and V in L, so x 1000 gives ng/mL.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    Cc <- 1000 * central / vc

    # Unbound plasma concentration: whole blood / (B:P) x fu,p. Table S1:
    # fu,p (mouse) = 0.09 and B:P = 2.06 at 77% hematocrit, which reproduce the
    # Table S4 unbound-plasma V 7.50 L and CL 1.35 L/h, and the Table 1 unbound
    # EC50 0.58 and MPC 0.79 ng/mL.
    fu <- 0.09
    bp <- 2.06
    Cu <- Cc * fu / bp

    # Equation 2: direct-effect sigmoidal kill driven by whole-blood Cc.
    effect <- Cc^hill / (Cc^hill + ec50^hill)
    d/dt(parasites) <- kgrow * parasites - kkill * effect * parasites
    parasites(0) <- exp(lrbase_parasites)

    # Equation 3: minimum parasiticidal concentration, the EC90 (ng/mL).
    mpc <- ec50 * 9^(1 / hill)

    Cc ~ prop(propSd)
    parasites ~ prop(propSd_parasites)
  })
}
