McCarthy_2016_artefenomel <- function() {
  description <- paste(
    "Population PK / PD model of the synthetic ozonide antimalarial artefenomel",
    "(OZ439) in healthy adult volunteers inoculated with Plasmodium falciparum",
    "blood-stage parasites (induced blood-stage malaria challenge model;",
    "McCarthy 2016). PK is two-compartment with first-order absorption after a",
    "lag time and a dose-dependent apparent clearance",
    "CL/F = 40.9 * (DOSE / 500)^-0.55 L/h, which reproduces the",
    "supra-proportional exposure seen across the 100, 200 and 500 mg single",
    "oral doses. PD is a parasite growth-and-kill model: parasitaemia grows at",
    "a fixed first-order rate G and is killed at a sigmoid-Emax rate",
    "D * C^H / (C^H + EC50^H) driven by the predicted plasma concentration,",
    "starting from the subject's parasitaemia at dosing. Residual error was not",
    "reported and is fixed to zero."
  )
  reference <- paste(
    "McCarthy JS, Baker M, O'Rourke P, Marquart L, Griffin P, Hooft van",
    "Huijsduijnen R, Moehrle JJ. Efficacy of OZ439 (artefenomel) against early",
    "Plasmodium falciparum blood-stage malaria infection in healthy volunteers.",
    "J Antimicrob Chemother. 2016;71(9):2620-2627. doi:10.1093/jac/dkw174."
  )
  vignette <- "McCarthy_2016_artefenomel"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "artefenomel", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "artefenomel", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "artefenomel", units = "mg", specimen = "tissue", verified = TRUE),
    parasites = list(
      analyte = "Plasmodium falciparum blood-stage parasites (qPCR of the 18S rRNA gene, all stages)",
      units = "parasites/mL",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  covariateData <- list(
    DOSE_ARTEFENOMEL_MG = list(
      description = "Administered single oral artefenomel dose",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed per subject (one oral dose of 100, 200 or 500 mg in the",
        "three cohorts). Enters the apparent clearance as",
        "CL/F = CL/F0 * (DOSE_ARTEFENOMEL_MG / 500)^-0.55 (Table 2 footnote c).",
        "The paper states the doses as administered amounts of artefenomel",
        "given as the mesylate salt and does not say whether they are salt or",
        "free-base masses; populate this column -- and the dose-record amt --",
        "with the same nominal doses the paper uses."
      ),
      source_name = "dose"
    ),
    PARA = list(
      description = "Plasmodium falciparum parasitaemia at the time of artefenomel dosing",
      units = "parasites/uL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed per subject. Seeds the parasite state at dosing:",
        "parasites(0) = PARA * 1000, because the model works in the paper's",
        "parasites/mL scale (1 parasite/uL = 1000 parasites/mL). The paper",
        "measured parasitaemia by qPCR of the 18S rRNA gene, which counts all",
        "parasite stages, and dosed once parasitaemia had crossed 1000",
        "parasites/mL (PARA = 1); the cohort medians at dosing in Figure 2",
        "lie between about 2500 and 6500 parasites/mL. The paper does not",
        "report how the initial parasitaemia was handled in the PK/PD fit."
      ),
      source_name = "parasitaemia at dosing (no column name printed)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 24L,
    n_studies = 1L,
    age_range = "20-34 years (mean 25.5, SD 3.7)",
    weight_range = "53.9-97.6 kg (mean 73.1, SD 12)",
    sex_female_pct = 50,
    race_ethnicity = c(Caucasian = 100),
    disease_state = paste(
      "Healthy volunteers inoculated with ~1800 P. falciparum-infected",
      "erythrocytes (induced blood-stage malaria challenge model)"
    ),
    dose_range = "Single oral dose of 100, 200 or 500 mg artefenomel (8 subjects per dose)",
    regions = "Australia (QIMR Berghofer, Brisbane)",
    notes = paste(
      "Demographics from Supplementary Table S1. Artefenomel was given as a",
      "suspension after 200 mL of full-cream milk once qPCR parasitaemia had",
      "reached 1000 parasites/mL. PK was fitted to all 24 subjects. The 100 mg",
      "cohort was rescued at 60 h and contributed no PD analysis of its own;",
      "the PK/PD parameters (Table 3) were estimated sequentially with the",
      "subject-specific predicted PK profiles driving the parasite model",
      "(Supplementary data, 'PK/PD model structure')."
    )
  )

  ini({
    # ---- PK: Table 2 (= Supplementary Table S5) ---------------------------
    lka <- log(0.18) ; label("First-order absorption rate constant ka (1/h)") # Table 2 'Ka (h-1) = 0.18 (0.023-13%)'
    ltlag <- log(0.40) ; label("Absorption lag time ALAG (h)") # Table 2 'ALAG (h) = 0.40 (0.03-8%)'
    lcl <- log(40.9) ; label("Apparent clearance CL/F0 at the 500 mg reference dose (L/h)") # Table 2 'CL/F (L/h) = 40.9 (3.8-9%)'
    lvc <- log(189) ; label("Apparent central volume Vc/F (L)") # Table 2 'Vc/F (L) = 189 (61-32%)'
    lvp <- log(1573) ; label("Apparent peripheral volume Vp/F (L)") # Table 2 'Vp/F (L) = 1573 (177-11%)'
    lq <- log(12.7) ; label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2 'Q/F (L/h) = 12.7 (2.0-16%)'
    e_dose_cl <- -0.55 ; label("Power exponent of dose on CL/F, reference 500 mg (unitless)") # Table 2 'Dose on CL/F = 0.55 (0.11-20%)'; footnote c 'CL/Fdose = CL/F0 x (dose/500)^-0.55'

    # ---- PD: Table 3 (= Supplementary Table S6) ---------------------------
    lkgrow <- fixed(log(0.06)) ; label("First-order parasite growth rate G (1/h)") # Table 3 'G (1/h) = 0.06 (fixed)'
    lemax <- log(0.20) ; label("Maximum first-order parasite kill rate D (1/h)") # Table 3 'D (1/h) = 0.20 (0.015)'
    lec50 <- log(6.19) ; label("Artefenomel concentration giving half the maximum kill rate EC50 (ng/mL)") # Table 3 'EC50 (ng/mL) = 6.19 (1.69)'
    lhill <- log(2.03) ; label("Hill coefficient H of the concentration-kill relationship (unitless)") # Table 3 'H = 2.03 (1.02)'

    # ---- IIV: variances of the random effects (Table 2 footnote a) --------
    etalvc ~ 0.269 # Table 2 'Vc/F' IIV variance 0.269 (SE 0.271)
    etalvp ~ 0.071 # Table 2 'Vp/F' IIV variance 0.071 (SE 0.08)
    etalcl ~ 0.064 # Table 2 'CL/F' IIV variance 0.064 (SE 0.04)
    etaltlag ~ 0.001 # Table 2 'ALAG' IIV variance 0.001 (SE 0.001)
    etalemax ~ 0.017 # Table 3 'D' IIV 0.017; Supplementary Table S6 'IIV (SE) = 0.017 (0.03)'
    etalec50 ~ 0.408 # Table 3 'EC50' IIV 0.408; Supplementary Table S6 'IIV (SE) = 0.408 (0.45)'

    # ---- Residual error: not reported anywhere in the paper or supplement --
    propSd <- fixed(0) ; label("Proportional residual error on artefenomel concentration (fraction; not reported in the source)") # not reported in Table 2 or Supplementary Table S5
    propSd_parasites <- fixed(0) ; label("Proportional residual error on parasitaemia (fraction; not reported in the source)") # not reported in Table 3 or Supplementary Table S6
  })
  model({
    # ---- Individual PK parameters ----------------------------------------
    ka <- exp(lka)
    tlag <- exp(ltlag + etaltlag)
    # Table 2 footnote c: CL/Fdose = CL/F0 x (dose/500)^-0.55
    cl <- exp(lcl + etalcl) * (DOSE_ARTEFENOMEL_MG / 500)^e_dose_cl
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- Individual PD parameters (Table 3) ------------------------------
    kgrow <- exp(lkgrow)
    emax <- exp(lemax + etalemax)
    ec50 <- exp(lec50 + etalec50)
    hill <- exp(lhill)

    d/dt(depot) <- -ka * depot
    alag(depot) <- tlag
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg and volume in L give mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc

    # Supplementary data Equation 1:
    #   dP/dt = P x (G - D x C(t)^H / (C(t)^H + EC50^H))
    kkill <- emax * Cc^hill / (Cc^hill + ec50^hill)
    d/dt(parasites) <- parasites * (kgrow - kkill)
    # PARA is in parasites/uL; the paper's parasitaemia scale is parasites/mL.
    parasites(0) <- PARA * 1000

    Cc ~ prop(propSd)
    parasites ~ prop(propSd_parasites)
  })
}
