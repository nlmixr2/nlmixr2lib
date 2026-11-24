Krause_2016_ACT451840 <- function() {
  description <- paste(
    "Joint population PK/PD model of the antimalarial ACT-451840",
    "(Actelion-451840) in healthy male volunteers given an induced",
    "blood-stage Plasmodium falciparum (strain 3D7) infection (Krause 2016).",
    "PK: two-compartment model with first-order absorption, an absorption",
    "lag time and first-order elimination (micro-rate-constant",
    "parameterisation, apparent volume V/F), single 500 mg oral dose under",
    "fed conditions. PD: exponential parasite growth minus a sigmoid-Hill",
    "drug kill driven by the plasma concentration,",
    "dP/dt = P * (G - D * C^gamma / (C^gamma + IC50^gamma)), with the parasite",
    "density P0 at inoculation as initial condition. Model time 0 is the",
    "inoculation; the study dosed the drug 168 h later.",
    sep = " "
  )
  reference <- paste(
    "Krause A, Dingemanse J, Mathis A, Marquart L, Moehrle JJ, McCarthy JS.",
    "Pharmacokinetic/pharmacodynamic modelling of the antimalarial effect of",
    "Actelion-451840 in an induced blood stage malaria study in healthy",
    "subjects. Br J Clin Pharmacol. 2016;82(2):412-421.",
    "doi:10.1111/bcp.12962.",
    sep = " "
  )
  vignette <- "Krause_2016_ACT451840"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # The oral drug dose enters `depot`; the parasite inoculum is the initial
  # condition of `parasites` at time 0 (no dosing record).
  dosing <- "depot"

  compartmentData <- list(
    depot = list(
      analyte = "ACT-451840",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "ACT-451840",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "ACT-451840",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    parasites = list(
      analyte = "Plasmodium falciparum (3D7) asexual blood-stage parasites",
      units = "parasites/mL",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 8L,
    n_studies = 1L,
    age_range = "19-38 years",
    age_mean = "24 years",
    weight_range = "59-98 kg",
    weight_mean = "77 kg",
    sex_female_pct = 0,
    race_ethnicity = c(Caucasian = 4, Indian = 2, Asian = 1, `Native Hawaiian` = 1),
    disease_state = paste(
      "Healthy malaria-naive male volunteers inoculated intravenously with",
      "about 1800 viable P. falciparum 3D7-infected erythrocytes (induced",
      "blood-stage malaria challenge)"
    ),
    dose_range = "500 mg ACT-451840 single oral dose (suspension) under fed conditions, 168 h after inoculation",
    regions = "Australia",
    notes = paste(
      "Single-centre open-label phase 1b proof-of-concept study. Counts in",
      "race_ethnicity are subject numbers (n = 8). Demographics from the",
      "Results 'Subject demographics' paragraph; BMI 24 (20-28) kg/m^2,",
      "height 179 (167-191) cm. Rescue artemether/lumefantrine was given",
      "on recrudescence; parasite data were fitted jointly with PK in",
      "Monolix 4.2 (SAEM) with M3 handling of BLQ/below-LOD data."
    )
  )

  ini({
    # PK (Krause 2016 Table 1, population parameters). Micro-rate-constant
    # parameterisation; V is apparent (V/F) because the dose was oral.
    lka <- log(0.23); label("First-order absorption rate constant ka (1/h)")             # Table 1 row 'ka (1/h) = 0.23' (RSE 6%)
    ltlag <- log(0.45); label("Absorption lag time tlag (h)")                            # Table 1 row 'tlag (h) = 0.45' (RSE 2%)
    lkel <- log(0.31); label("First-order elimination rate constant k (1/h)")            # Table 1 row 'k (1/h) = 0.31' (RSE 12%)
    lvc <- log(1337.89); label("Apparent central volume of distribution V/F (L)")        # Table 1 row 'V/F (l) = 1337.89' (RSE 9%)
    lk12 <- log(0.11); label("Central-to-peripheral rate constant k12 (1/h)")            # Table 1 row 'k12 (1/h) = 0.11' (RSE 7%)
    lk21 <- log(0.03); label("Peripheral-to-central rate constant k21 (1/h)")            # Table 1 row 'k21 (1/h) = 0.03' (RSE 7%)

    # PD (Krause 2016 Table 1). Table 1 prints G and D as 'unitless' values
    # 6.67 and 18.09; they are rates in percent per hour, so the per-hour rate
    # constants are the tabulated values divided by 100. The scaling is
    # pinned by Figure 5B (typical log10 parasite count rises from
    # log10(0.11) = -0.96 at inoculation to about 3.9 at dosing 168 h later,
    # i.e. ln(10) * 4.87 / 168 = 0.067 per hour) and by Figure 3, whose
    # y-axis 'Drug effect (growth rate reduction)' plateaus at 18.
    lrbase_parasites <- log(0.11); label("Parasite density at inoculation P0 (parasites/mL)") # Table 1 row 'P0 (count/ml) = 0.11' (RSE 58%)
    lkgrow <- log(6.67 / 100); label("First-order parasite growth rate G (1/h)")               # Table 1 row 'G (unitless) = 6.67' (RSE 4%), percent per hour
    lemax <- log(18.09 / 100); label("Maximum drug-induced parasite kill rate D (1/h)")        # Table 1 row 'D (unitless) = 18.09' (RSE 5%), percent per hour
    lec50 <- log(2.70); label("Plasma concentration for half-maximal kill IC50 (ng/mL)")       # Table 1 row 'IC50 (ng/ml) = 2.70' (RSE 10%)
    lhill <- log(2.46); label("Hill coefficient gamma of the kill curve (unitless)")           # Table 1 row 'gamma (unitless) = 2.46' (RSE 7%)

    # IIV: Table 1 'Random effects (standard deviations)'. Monolix
    # log-normal random effects; variance = SD^2.
    etalka ~ 0.0196            # Table 1 random effect 'ka' SD = 0.14; 0.14^2
    etaltlag ~ 0.0025          # Table 1 random effect 'tlag' SD = 0.05; 0.05^2
    etalkel ~ 0.0841           # Table 1 random effect 'k' SD = 0.29; 0.29^2
    etalvc ~ 0.04              # Table 1 random effect 'V' SD = 0.20; 0.20^2
    etalk12 ~ 0.0036           # Table 1 random effect 'k12' SD = 0.06; 0.06^2
    etalk21 ~ 0.0064           # Table 1 random effect 'k21' SD = 0.08; 0.08^2
    etalrbase_parasites ~ 1.3924 # Table 1 random effect 'P0' SD = 1.18; 1.18^2
    etalkgrow ~ 0.0009         # Table 1 random effect 'G' SD = 0.03; 0.03^2
    etalemax ~ 0.0016          # Table 1 random effect 'D' SD = 0.04; 0.04^2
    etalec50 ~ 0.0625          # Table 1 random effect 'IC50' SD = 0.25; 0.25^2
    etalhill ~ 0.0064          # Table 1 random effect 'gamma' SD = 0.08; 0.08^2

    # Residual error: Monolix 4.2 combined error model, a + b * f.
    addSd <- 4.38e-9; label("Additive residual error on plasma ACT-451840 (ng/mL)")            # Table 1 row 'Additive error PK = 4.38E-09' (RSE 3.72E+06%; the paper notes it is effectively zero)
    propSd <- 0.265; label("Proportional residual error on plasma ACT-451840 (fraction)")      # Table 1 row 'Multiplicative error PK = 0.265' (RSE 7%)
    addSd_parasites <- 26.1; label("Additive residual error on parasite density (parasites/mL)") # Table 1 row 'Additive error PD = 26.1' (RSE 31%)
    propSd_parasites <- 0.558; label("Proportional residual error on parasite density (fraction)") # Table 1 row 'Multiplicative error PD = 0.558' (RSE 7%)
  })

  model({
    # 1. Individual PK parameters.
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag + etaltlag)
    kel <- exp(lkel + etalkel)
    vc <- exp(lvc + etalvc)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)

    # 2. Individual PD parameters.
    rbase <- exp(lrbase_parasites + etalrbase_parasites)
    kgrow <- exp(lkgrow + etalkgrow)
    emax <- exp(lemax + etalemax)
    ec50 <- exp(lec50 + etalec50)
    hill <- exp(lhill + etalhill)

    # 3. PK ODEs (dose in mg, volume in L).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    alag(depot) <- tlag

    # Plasma concentration in ng/mL (mg/L x 1000).
    Cc <- 1000 * central / vc

    # 4. Parasite growth and drug kill (Krause 2016 PK/PD modelling equation).
    #    The inoculum density P0 is the state at time 0 (inoculation).
    kill <- emax * Cc^hill / (Cc^hill + ec50^hill)
    d/dt(parasites) <- parasites * (kgrow - kill)
    parasites(0) <- rbase

    # 5. Observations: plasma ACT-451840 and PCR parasite density.
    Cc ~ add(addSd) + prop(propSd) + combined1()
    parasites ~ add(addSd_parasites) + prop(propSd_parasites) + combined1()
  })
}
