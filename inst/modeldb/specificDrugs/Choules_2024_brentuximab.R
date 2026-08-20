Choules_2024_brentuximab <- function() {
  description <- "Compartmental reduction of the Simcyp ADC-module physiologically based pharmacokinetic (PBPK) model for brentuximab vedotin (conjugated antibody) and its released payload monomethyl auristatin E (MMAE) in patients with haematological malignancy (Choules 2024). This is the comparator model that Choules 2024 built alongside enfortumab vedotin, using the same valine-citrulline-MMAE linker and the same MMAE compound file, in order to verify the MMAE drug-drug interaction predictions against an observed clinical interaction study. The conjugated antibody is a linear 2-compartment model with first-order elimination from the central compartment; every milligram of eliminated antibody releases its full MMAE payload (mean drug-antibody ratio 4.026) into a 1-compartment MMAE disposition model whose volume and clearance are the paper's published Simcyp compound-file inputs. Combined P-gp / CYP3A4 perpetrators (ketoconazole, rifampin) act as multiplicative factors on MMAE clearance. The published model itself is a Simcyp platform PBPK model whose whole-body physiology is not tabulated; the antibody parameters here were determined from the paper's published simulated summary statistics, so this file is a compartmental reduction and not a PBPK model. See the vignette for the derivation and the held-out validation gates."
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
    central      = list(analyte = "brentuximab vedotin (conjugated antibody)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1  = list(analyte = "brentuximab vedotin (conjugated antibody)", units = "mg", specimen = "plasma", verified = TRUE),
    central_mmae = list(analyte = "monomethyl auristatin E (MMAE)",            units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CONMED_KETOCONAZOLE = list(
      description        = "Concomitant ketoconazole (combined P-glycoprotein and strong CYP3A4 inhibitor) coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "1 = ketoconazole 400 mg orally once daily coadministered with brentuximab vedotin, 0 = brentuximab vedotin alone. Used here as a step indicator at full inhibition; the source's brentuximab vedotin interaction simulations mirror the clinical study NCT01026415, in which ketoconazole dosing brackets the antibody dose. The effect is applied as a multiplicative factor on MMAE clearance only; the conjugated antibody is unaffected. The magnitude is back-calculated from the published MMAE AUC geometric mean ratio, not from an inhibition constant: because this reduction has MMAE AUC(inf) = payload dose / MMAE CL exactly, a CL multiplier of 1/GMR reproduces the published AUC ratio by construction. Reference category is 'no ketoconazole'.",
      source_name        = NA_character_
    ),
    CONMED_RIF = list(
      description        = "Concomitant rifampin (rifampicin; combined P-glycoprotein and strong CYP3A4 inducer) coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "1 = rifampin 600 mg orally once daily coadministered with brentuximab vedotin, 0 = brentuximab vedotin alone. Step indicator at the post-induction CYP3A4 equilibrium, matching the chronic-induction semantics of the canonical and the design of the clinical interaction study NCT01026415. The effect is applied as a multiplicative factor on MMAE clearance only. Magnitude back-calculated from the published MMAE AUC geometric mean ratio as described for CONMED_KETOCONAZOLE. Note that the source's rifampin simulations additionally applied a fold increase in the P-gp relative activity factor to MMAE to represent transporter induction, which the platform did not support natively; that mechanism is folded into the single lumped CL multiplier here. Reference category is 'no rifampin'.",
      source_name        = NA_character_
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 120L,
    n_studies      = 2L,
    age_range      = "22-70 years (assumed by the source when constructing the virtual trial)",
    weight_mean    = "72.0 kg (assumed the same Simcyp cancer population mean as the enfortumab vedotin simulations; see notes)",
    sex_female_pct = 40,
    disease_state  = "Haematological malignancy (CD30-positive lymphoma); unmodified Simcyp cancer population model.",
    dose_range     = "1.8 and 2.7 mg/kg as a 30-minute intravenous infusion, single dose; 1.2 and 1.8 mg/kg in the drug-interaction simulations.",
    notes          = "The simulated virtual trials used 10 trials of 12 participants each over a 504 h study duration (Choules 2024 Supplementary Methods); ages 22-70 years and 40% women were explicitly stated by the source as assumed. Observed comparator data were digitised from the brentuximab vedotin clinical pharmacology submission to the US Food and Drug Administration, and the drug-interaction comparator is the clinical study NCT01026415. The source does not print a body weight for the brentuximab vedotin virtual population; 72.0 kg (the enfortumab vedotin simulated mean from Table S4) is used here because both simulations draw from the Simcyp cancer population, and it is corroborated by the reduction reproducing the published MMAE AUC(inf) at both dose levels to within 1% (see the vignette). All parameters in this file are calibrated at that weight."
  )

  ini({
    # Conjugated-antibody disposition. The source is a Simcyp minimal-PBPK ADC
    # model whose compartment volumes are platform outputs and are not
    # tabulated, so these four values were determined from the paper's own
    # published simulated summary statistics for brentuximab vedotin: the
    # arithmetic mean Cmax and AUC(inf) at 1.8 and 2.7 mg/kg and the mean
    # terminal half-life (Choules 2024 Table S3).
    lcl <- log(0.065773);  label("Conjugated-antibody clearance (CL, L/h)")                              # derived: Choules 2024 Table S3 simulated mean AUC(inf) 82.2 ug*day/mL (1.8 mg/kg) and 123 (2.7 mg/kg) at 72.0 kg
    lvc <- log(4.0966);    label("Conjugated-antibody central volume (V1, L)")                           # derived: Choules 2024 Table S3 simulated mean Cmax 31.5 ug/mL (1.8 mg/kg) and 47.2 (2.7 mg/kg)
    lq  <- log(0.0057794); label("Conjugated-antibody inter-compartmental clearance (Q, L/h)")           # derived: Choules 2024 Table S3 simulated mean terminal half-life 5.71 day and MMAE mean Cmax 4.62 / 6.96 ng/mL
    lvp <- log(1.0127);    label("Conjugated-antibody peripheral volume (V2, L)")                        # derived: Choules 2024 Table S3 simulated mean terminal half-life 5.71 day and MMAE mean Cmax 4.62 / 6.96 ng/mL

    # MMAE disposition. Choules 2024 used the SAME MMAE compound file for
    # brentuximab vedotin and enfortumab vedotin, so these are the identical
    # published Simcyp inputs used in Choules_2024_enfortumab and are not
    # fitted here.
    lcl_mmae <- fixed(log(2.72));   label("MMAE clearance (CL, L/h)")                                    # Choules 2024 Table 2: Simcyp retrograde calculation based on apparent CLiv 2.72 L/h
    lvc_mmae <- fixed(log(218.2));  label("MMAE volume of distribution (Vss, L)")                        # Choules 2024 Table 2: Vss 3.03 L/kg x 72.0 kg assumed mean weight

    # Payload stoichiometry. Frel (deconjugation) and Frel (catabolic) are both
    # 1 in Choules 2024 Table S1, so the released fraction is 1 and the mass
    # conversion is purely stoichiometric.
    dar     <- fixed(4.026);  label("Mean drug-antibody ratio (mol MMAE per mol conjugated antibody)")   # Choules 2024 Table S1: mean DAR 4.026, automatically calculated by Simcyp from the discrete DAR distribution
    mw  <- fixed(148081); label("Brentuximab vedotin molecular weight (Da)")                         # Choules 2024 Table S1: 148,081 Da, from the Adcetris CHMP assessment report
    mw_mmae <- fixed(717.98); label("MMAE molecular weight (Da)")                                        # Choules 2024 Table 2: 717.98 Da

    # Combined P-gp / CYP3A4 perpetrator effects on MMAE clearance, from the
    # brentuximab vedotin rows of the interaction table. The reduction has
    # MMAE AUC(inf) = released payload dose / MMAE CL, so a CL multiplier equal
    # to the reciprocal of the published MMAE AUC geometric mean ratio
    # reproduces that ratio exactly.
    e_keto_cl_mmae <- fixed(0.72993); label("Multiplier of ketoconazole coadministration on MMAE CL: cl_mmae *= e_keto_cl_mmae^CONMED_KETOCONAZOLE") # derived: 1 / 1.37, Choules 2024 Table 5 predicted MMAE AUC ratio for brentuximab vedotin 1.8 mg/kg + ketoconazole
    e_rif_cl_mmae  <- fixed(2.12766); label("Multiplier of rifampin coadministration on MMAE CL: cl_mmae *= e_rif_cl_mmae^CONMED_RIF")               # derived: 1 / 0.47, Choules 2024 Table 5 predicted MMAE AUC ratio for brentuximab vedotin 1.8 mg/kg + rifampin

    # Residual error. Choules 2024 is a deterministic PBPK simulation study; it
    # reports no residual-error model. The CV% values in Table S3 are
    # between-subject spread of the simulated virtual population arising from
    # Simcyp physiological variability, not a residual-error estimate.
    propSd      <- fixed(0); label("Proportional residual error on Cc (fraction) -- not reported by the source")       # not reported in Choules 2024; see vignette Errata
    propSd_mmae <- fixed(0); label("Proportional residual error on Cc_mmae (fraction) -- not reported by the source")  # not reported in Choules 2024; see vignette Errata
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
    # in Choules 2024 Table S1 are lumped into cl) releases its payload.
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
