Choules_2024_enfortumab <- function() {
  description <- "Compartmental reduction of the Simcyp ADC-module physiologically based pharmacokinetic (PBPK) model for enfortumab vedotin (conjugated antibody) and its released payload monomethyl auristatin E (MMAE) in patients with locally advanced/metastatic urothelial carcinoma and other solid tumours (Choules 2024). The conjugated antibody is a linear 2-compartment model with first-order elimination from the central compartment; every milligram of eliminated antibody releases its full MMAE payload (mean drug-antibody ratio 3.735) into a 1-compartment MMAE disposition model whose volume and clearance are the paper's published Simcyp compound-file inputs. Combined P-gp / CYP3A4 perpetrators (ketoconazole, rifampin) act as multiplicative factors on MMAE clearance. The published model itself is a Simcyp platform PBPK model whose whole-body physiology is not tabulated; the antibody parameters here were determined from the paper's published simulated summary statistics and concentration-time figures, so this file is a compartmental reduction and not a PBPK model. See the vignette for the derivation and the held-out validation gates."
  reference <- "Choules MP, Zuo P, Otsuka Y, Garg A, Tang M, Bonate P. Physiologically based pharmacokinetic model to predict drug-drug interactions with the antibody-drug conjugate enfortumab vedotin. J Pharmacokinet Pharmacodyn. 2024;51(5):417-428. doi:10.1007/s10928-023-09877-5. PMID 37624557; PMCID PMC11576838."
  vignette <- "Choules_2024_vedotin_ddi"

  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ug/mL"
  )

  # Cc (conjugated antibody) is in ug/mL == mg/L, matching the paper's ADC
  # figures and tables. Cc_mmae is carried in the SAME mg/L system for internal
  # consistency; the paper reports MMAE in ng/mL, which is 1000 * Cc_mmae.
  compartmentData <- list(
    central      = list(analyte = "enfortumab vedotin (conjugated antibody)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1  = list(analyte = "enfortumab vedotin (conjugated antibody)", units = "mg", specimen = "plasma", verified = TRUE),
    central_mmae = list(analyte = "monomethyl auristatin E (MMAE)",           units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CONMED_KETOCONAZOLE = list(
      description        = "Concomitant ketoconazole (combined P-glycoprotein and strong CYP3A4 inhibitor) coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "1 = ketoconazole 400 mg orally once daily coadministered with enfortumab vedotin (Choules 2024 Table S2 virtual trial design: enfortumab vedotin 1.25 mg/kg on day 4 of daily ketoconazole, 984 h simulation), 0 = enfortumab vedotin alone. Time-varying in principle; in the source simulation ketoconazole dosing starts before the enfortumab vedotin dose and continues throughout, so the indicator is effectively on for the whole MMAE observation window and is used here as a step indicator. The effect is applied as a multiplicative factor on MMAE clearance only; the conjugated antibody is unaffected (Choules 2024 predicted no change in ADC exposure). The magnitude is back-calculated from the published MMAE AUC geometric mean ratio, not from an inhibition constant: because this reduction has MMAE AUC(inf) = payload dose / MMAE CL exactly, a CL multiplier of 1/GMR reproduces the published AUC ratio by construction. Reference category is 'no ketoconazole'.",
      source_name        = NA_character_
    ),
    CONMED_RIF = list(
      description        = "Concomitant rifampin (rifampicin; combined P-glycoprotein and strong CYP3A4 inducer) coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "1 = rifampin 600 mg orally once daily coadministered with enfortumab vedotin (Choules 2024 Table S2 virtual trial design: enfortumab vedotin 1.25 mg/kg on day 8 of daily rifampin, 936 h simulation), 0 = enfortumab vedotin alone. The 7 days of rifampin before the enfortumab vedotin dose put CYP3A4 induction at the post-induction equilibrium, so this is a step indicator at full induction, matching the chronic-induction semantics of the canonical. The effect is applied as a multiplicative factor on MMAE clearance only. Magnitude back-calculated from the published MMAE AUC geometric mean ratio as described for CONMED_KETOCONAZOLE. Note that the source's rifampin simulation additionally applied a fold increase in the P-gp relative activity factor to MMAE to represent transporter induction, which the platform did not support natively; that mechanism is folded into the single lumped CL multiplier here. Reference category is 'no rifampin'.",
      source_name        = NA_character_
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 100L,
    n_studies      = 2L,
    age_range      = "24-83 years (1.25 mg/kg virtual trial); 48-81 years (1.0 mg/kg virtual trial); 40-90 years (multiple-dose and drug-interaction virtual trials)",
    age_mean       = "64.9 years (simulated); 65.1 years (observed)",
    weight_mean    = "72.0 kg (simulated); 80.3 kg (observed)",
    sex_female_pct = 26,
    disease_state  = "Locally advanced or metastatic Nectin-4-expressing urothelial carcinoma and other malignant solid tumours; Simcyp cancer population model with the plasma tissue-volume scaling factor changed from 1.2 to 1.0.",
    dose_range     = "1.0 and 1.25 mg/kg as a 30-minute intravenous infusion; single dose and days 1, 8 and 15 of a 28-day cycle.",
    notes          = "The simulated virtual trials used the Simcyp cancer population: 10 trials of 10 participants for the single-dose simulations and 10 trials of 15 participants for the multiple-dose verification (Choules 2024 Supplementary Methods). Observed comparator data are the phase 1 EV-101 (NCT02091999, n = 201 for the demographic comparison in Table S4) and phase 2 EV-201 (NCT03219333, n = 125) studies. Simulated versus observed demographics (mean (SD)): age 64.9 (12.2) vs 65.1 (10.1) years; weight 72.0 (14.1) vs 80.3 (20.1) kg; plasma albumin 38.2 (7.1) vs 37.5 (4.1) g/dL as printed; haematocrit 37.6 (4.65) vs 35.3 (4.76)% (Choules 2024 Table S4). All parameters in this file are calibrated at the 72.0 kg simulated mean weight; see the vignette for why no body-weight relationship is identifiable from the published outputs."
  )

  ini({
    # Conjugated-antibody disposition. The source is a Simcyp minimal-PBPK ADC
    # model whose compartment volumes are platform outputs and are not
    # tabulated, so these four values were determined from the paper's own
    # published simulated summary statistics for enfortumab vedotin: the
    # geometric-mean-adjacent arithmetic mean Cmax and AUC(d0-7) at 1.0 and
    # 1.25 mg/kg and the mean terminal half-life (Choules 2024 Table 3). The
    # plotted simulated curves in Figure 3 are arithmetic means, so the
    # arithmetic-mean column of Table 3 is the calibration target.
    lcl <- log(0.10969); label("Conjugated-antibody clearance (CL, L/h)")                                # derived: Choules 2024 Table 3 simulated mean AUC(d0-7) 32.6 ug*day/mL (1.25 mg/kg) and 26.0 (1.0 mg/kg) at 72.0 kg
    lvc <- log(3.4842);  label("Conjugated-antibody central volume (V1, L)")                             # derived: Choules 2024 Table 3 simulated mean Cmax 25.3 ug/mL (1.25 mg/kg) and 20.6 (1.0 mg/kg)
    lq  <- log(0.040504); label("Conjugated-antibody inter-compartmental clearance (Q, L/h)")            # derived: Choules 2024 Table 3 simulated mean terminal half-life 1.96 day with Figure 3a curve shape
    lvp <- log(1.6225);  label("Conjugated-antibody peripheral volume (V2, L)")                          # derived: Choules 2024 Table 3 simulated mean terminal half-life 1.96 day with Figure 3a curve shape

    # MMAE disposition. Both values are the paper's published Simcyp compound-file
    # inputs for MMAE and are NOT fitted here; they are used exactly as printed.
    lcl_mmae <- fixed(log(2.72));   label("MMAE clearance (CL, L/h)")                                    # Choules 2024 Table 2: Simcyp retrograde calculation based on apparent CLiv 2.72 L/h
    lvc_mmae <- fixed(log(218.2));  label("MMAE volume of distribution (Vss, L)")                        # Choules 2024 Table 2: Vss 3.03 L/kg x 72.0 kg simulated mean weight (Table S4)

    # Payload stoichiometry. Each milligram of eliminated conjugated antibody
    # releases dar molecules of MMAE; Frel (deconjugation) and Frel (catabolic)
    # are both 1 in Choules 2024 Table 1, so the released fraction is 1 and the
    # mass conversion is purely stoichiometric.
    dar     <- fixed(3.735);    label("Mean drug-antibody ratio (mol MMAE per mol conjugated antibody)") # Choules 2024 Table 1: mean DAR 3.735, automatically calculated by Simcyp from the discrete DAR distribution
    mw  <- fixed(146664.1); label("Enfortumab vedotin molecular weight (Da)")                        # Choules 2024 Table 1: 146,664.1 Da, from the investigational new drug application
    mw_mmae <- fixed(717.98);   label("MMAE molecular weight (Da)")                                      # Choules 2024 Table 2: 717.98 Da

    # Combined P-gp / CYP3A4 perpetrator effects on MMAE clearance. The
    # reduction has MMAE AUC(inf) = released payload dose / MMAE CL, so a CL
    # multiplier equal to the reciprocal of the published MMAE AUC geometric
    # mean ratio reproduces that ratio exactly. Choules 2024 states that
    # simulations were run long enough for the AUC(last) ratio to equal the
    # AUC(inf) ratio, so the two are interchangeable here.
    e_keto_cl_mmae <- fixed(0.72464); label("Multiplier of ketoconazole coadministration on MMAE CL: cl_mmae *= e_keto_cl_mmae^CONMED_KETOCONAZOLE") # derived: 1 / 1.38, Choules 2024 Table 5 predicted MMAE AUC ratio for enfortumab vedotin 1.25 mg/kg + ketoconazole
    e_rif_cl_mmae  <- fixed(2.12766); label("Multiplier of rifampin coadministration on MMAE CL: cl_mmae *= e_rif_cl_mmae^CONMED_RIF")               # derived: 1 / 0.47, Choules 2024 Table 5 predicted MMAE AUC ratio for enfortumab vedotin 1.25 mg/kg + rifampin

    # Residual error. Choules 2024 is a deterministic PBPK simulation study; it
    # reports no residual-error model. The CV% values in Tables 3, 4 and 6 are
    # between-subject spread of the simulated virtual population arising from
    # Simcyp physiological variability, not a residual-error estimate, so they
    # are not reused here.
    propSd      <- fixed(0); label("Proportional residual error on Cc (fraction) -- not reported by the source")           # not reported in Choules 2024; see vignette Errata
    propSd_mmae <- fixed(0); label("Proportional residual error on Cc_mmae (fraction) -- not reported by the source")      # not reported in Choules 2024; see vignette Errata
  })

  model({
    # 1. Derived terms. Mass of MMAE released per unit mass of conjugated
    # antibody eliminated, from the mean DAR and the two molecular weights.
    fpayload <- dar * mw_mmae / mw

    # 2. Individual parameters
    cl <- exp(lcl)
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp)
    cl_mmae <- exp(lcl_mmae) *
      e_keto_cl_mmae^CONMED_KETOCONAZOLE *
      e_rif_cl_mmae^CONMED_RIF
    vc_mmae <- exp(lvc_mmae)

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    kel_mmae <- cl_mmae / vc_mmae

    # 4. ODE system. All conjugated antibody leaving the central compartment
    # (catabolism, additional non-specific plasma clearance and deconjugation
    # in Choules 2024 Table 1 are lumped into cl) releases its payload.
    d/dt(central)      <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1)  <-  k12 * central - k21 * peripheral1
    d/dt(central_mmae) <-  fpayload * kel * central - kel_mmae * central_mmae

    # 5. Observations. Cc is the conjugated-antibody plasma concentration in
    # ug/mL. Cc_mmae is the unconjugated MMAE plasma concentration in the same
    # mg/L system; multiply by 1000 to obtain the ng/mL scale used by the paper.
    Cc      <- central / vc
    Cc_mmae <- central_mmae / vc_mmae
    Cc      ~ prop(propSd)
    Cc_mmae ~ prop(propSd_mmae)
  })
}
