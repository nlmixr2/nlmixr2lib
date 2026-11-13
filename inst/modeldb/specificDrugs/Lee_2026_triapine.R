Lee_2026_triapine <- function() {
  description <- "Two-compartment population PK model of triapine with a two-transit-compartment Erlang oral absorption chain, logit-normal oral bioavailability, allometric weight scaling and a current-smoking effect on clearance, linked through an effect compartment to a multiplicative sigmoid Emax model of blood methemoglobin (% of total hemoglobin), in adults with locally advanced cervical / vaginal cancer or gastroenteropancreatic neuroendocrine tumours receiving IV and oral triapine"
  reference <- "Lee H, Dunn A, Taylor SE, Chauhan A, Ivy SP, Gobburu J, Beumer JH. Enhancing triapine treatment: strategies for dose optimization and methemoglobin level mitigation. Cancer Chemother Pharmacol. 2026;96(1):55. doi:10.1007/s00280-026-04898-6"
  vignette <- "Lee_2026_triapine"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ug/L",
    methb = "% of total hemoglobin"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Actual body weight; allometric scaling with reference 70 kg and exponents fixed at 0.75 (CL, Q) and 1 (Vc, Vp) (Methods 'Covariate PK Model' equation and Methods S1 code: wtcl = (WEIGHT/70)^0.75, wtv = WEIGHT/70).",
      source_name = "WEIGHT"
    ),
    SMOKE = list(
      description = "Current (active) smoker indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (nonsmoker)",
      notes = "Results smoking equation: 'smoking status is a binary variable (0 = nonsmoker, 1 = active smoker)'. Multiplicative linear effect on CL only, CL * (1 + e_smoke_cl * SMOKE). 6 of 36 subjects were smokers (all in NCT02595879, Table 1).",
      source_name = "SMOKE_NUM"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "triapine", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "triapine", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "triapine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "triapine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "triapine", units = "mg", specimen = "plasma", verified = TRUE),
    effect = list(analyte = "triapine", units = "ug/L", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 36L,
    n_studies = 2L,
    age_range = "23-74 years",
    weight_range = "41.8-138.4 kg",
    sex_female_pct = 66.7,
    race_ethnicity = c(White = 91.7, Asian = 2.8, Black = 2.8, `Not reported` = 2.8),
    disease_state = "Locally advanced cervical or vaginal cancer receiving chemoradiation (NCT02595879, NCI 9892; n = 20) and well-differentiated somatostatin-receptor-positive gastroenteropancreatic neuroendocrine tumours receiving Lu-177 dotatate (NCT04234568, NCI 10388; n = 16)",
    dose_range = "Triapine 50 mg IV 2-h infusion (n = 20, day 1) and 100 or 150 mg oral once daily",
    regions = "United States",
    smokers = "6 of 36 active smokers (all NCT02595879)",
    notes = "Table 1. 573 quantifiable plasma triapine concentrations (32 BLQ at 24 h omitted, M1) from 36 patients; 252 methemoglobin observations from 33 patients. Smoking status 30 no / 6 yes. NCT04234568 weight reported only as [53.8, 138.4]; NCT02595879 weight mean 77.4 (SD 14.1) kg, range 41.8-103.9 kg."
  )

  ini({
    # PK (Table 2, 'Final Model Estimate'; clearances and volumes for a 70 kg subject)
    lktr <- log(3.11); label("Transit rate constant Ktr (1/h)") # Table 2 Ktr = 3.11 1/h (bootstrap 95% CI 2.63-3.72)
    lcl <- log(51.8); label("Clearance CL for a 70 kg nonsmoker (L/h)") # Table 2 CL = 51.8 L/h
    lvc <- log(67.9); label("Central volume of distribution Vc for 70 kg (L)") # Table 2 Vc = 67.9 L
    lq <- log(12.6); label("Intercompartmental clearance Q for 70 kg (L/h)") # Table 2 Q = 12.6 L/h
    lvp <- log(65.5); label("Peripheral volume of distribution Vp for 70 kg (L)") # Table 2 Vp = 65.5 L
    logitfdepot <- qlogis(0.76); label("Logit of oral bioavailability F (logit of fraction)") # Table 2 F = 0.76
    e_wt_cl <- fixed(0.75); label("Allometric exponent of weight on CL and Q (unitless)") # Methods 'Covariate PK Model': fixed to 0.75 for all clearance terms
    e_wt_vc <- fixed(1); label("Allometric exponent of weight on Vc and Vp (unitless)") # Methods 'Covariate PK Model': fixed to 1 for all volume terms
    e_smoke_cl <- 0.38; label("Fractional increase in CL in active smokers (unitless)") # Table 2 Effect_smoker = 0.38 (95% CI 0.19-0.62)

    # PD: effect compartment + multiplicative sigmoid Emax on methemoglobin
    lke0 <- log(0.36); label("Effect-compartment equilibration rate constant ke0 (1/h)") # Table 2 ke0 = 0.36 1/h
    le0 <- log(0.95); label("Baseline methemoglobin Ebase (% of total hemoglobin)") # Table 2 Ebase = 0.95 %
    lemax <- fixed(log(100)); label("Maximum fractional increase of methemoglobin over baseline Emax (unitless)") # Table 2 'Emax (fixed)' = 100
    lec50 <- log(2535); label("Effect-site triapine concentration at half-maximal effect EC50 (ug/L)") # Table 2 EC50 = 2535 (header says mg/L; Results text: 2535 ug/L)
    lhill <- log(1.63); label("Hill coefficient (unitless)") # Table 2 Hill = 1.63

    # IIV. Log-normal etas: omega^2 = log(1 + CV^2) from the Table 2 %CV
    # (Methods: %CV = sqrt(exp(omega^2) - 1) x 100).
    etalktr ~ 0.2049 # Table 2 omega_Ktr 47.8 %CV -> log(1 + 0.478^2)
    etalcl ~ 0.0797 # Table 2 omega_CL 28.8 %CV -> log(1 + 0.288^2)
    etalvc ~ 0.0481 # Table 2 omega_Vc 22.2 %CV -> log(1 + 0.222^2)
    # Logit-normal F. Methods: %CV = (1 - Fpop) * sqrt(omega^2) = 29.5 % with
    # Fpop = 0.76 gives omega^2 = (0.295 / 0.24)^2 = 1.511. Methods S1 code
    # passes that parameter as the SD argument of Julia Normal(mu, sigma)
    # (ebio ~ Normal(logit(tvbio), Omega2_bio)), so the fitted logit-scale SD
    # is 1.511 and the variance 1.511^2 = 2.283.
    etalogitfdepot ~ 2.283 # Table 2 omega_F 29.5 %CV with Methods S1 SD encoding
    etale0 ~ 0.1137 # Table 2 omega_Ebase 34.7 %CV -> log(1 + 0.347^2)
    etalec50 ~ 0.2891 # Table 2 omega_EC50 57.9 %CV -> log(1 + 0.579^2)

    # Residual error: combined, variances summed (Methods S1/S2 code:
    # Normal(f, sqrt(f^2 * sigma2_prop + sigma2_add)))
    addSd <- 3.94; label("Additive residual error, triapine (ug/L)") # Table 2 sigma PK_additive = 3.94 ug/L
    propSd <- 0.332; label("Proportional residual error, triapine (fraction)") # Table 2 sigma PK_proportional = 33.2 %CV
    addSd_methb <- 0.13; label("Additive residual error, methemoglobin (% of total hemoglobin)") # Table 2 sigma PD_additive = 0.13 (unit printed as ug/L; the endpoint is in %)
    propSd_methb <- 0.273; label("Proportional residual error, methemoglobin (fraction)") # Table 2 sigma PD_proportional = 27.3 %CV
  })

  model({
    # Individual PK parameters (Methods S1 @pre block)
    ktr <- exp(lktr + etalktr)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (1 + e_smoke_cl * SMOKE)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q <- exp(lq) * (WT / 70)^e_wt_cl
    vp <- exp(lvp) * (WT / 70)^e_wt_vc
    fdepot <- expit(logitfdepot + etalogitfdepot)

    # Individual PD parameters (Table 2; BSV on Ebase and EC50)
    ke0 <- exp(lke0)
    e0 <- exp(le0 + etale0)
    emax <- exp(lemax)
    ec50 <- exp(lec50 + etalec50)
    hill <- exp(lhill)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-transit Erlang absorption: depot -> transit1 -> transit2 -> central,
    # every transfer at rate ktr (Methods S1 @dynamics). IV doses go to central.
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(central) <- ktr * transit2 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    f(depot) <- fdepot

    # Dose in mg and volumes in L give mg/L; x 1000 for ug/L
    Cc <- 1000 * central / vc

    # Effect compartment (Methods ke0 equation), concentration units ug/L
    d/dt(effect) <- ke0 * (Cc - effect)

    # Methods effect equation: E = Ebase * (1 + Emax * Ceff^h / (EC50^h + Ceff^h))
    methb <- e0 * (1 + emax * effect^hill / (ec50^hill + effect^hill))

    Cc ~ add(addSd) + prop(propSd) + combined2()
    methb ~ add(addSd_methb) + prop(propSd_methb) + combined2()
  })
}
