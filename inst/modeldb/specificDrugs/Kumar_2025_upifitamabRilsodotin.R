Kumar_2025_upifitamabRilsodotin <- function() {
  description <- "Compartmental reduction of the Simcyp physiologically based pharmacokinetic (PBPK) model for the antibody-drug conjugate upifitamab rilsodotin (NaPi2b-targeted, ester-linked auristatin F-hydroxypropylamide (AF-HPA) payload) in adults with solid tumours (Kumar 2025). Conjugated AF-HPA (the ADC dose expressed as payload equivalents) follows a 2-compartment model whose peripheral compartment is the paper's single adjusting compartment; a fixed fraction of its clearance (payload release by linker hydrolysis) forms unconjugated AF-HPA, which follows a 1-compartment model. Itraconazole (combined CYP3A4 and P-gp inhibitor) acts as a multiplicative factor on unconjugated AF-HPA clearance. The published model is a Simcyp minimal-PBPK (conjugate) plus full-PBPK (payload) model whose whole-body physiology is not tabulated; the clearances here are derived from the paper's printed renal clearances and rat excretion fractions, the conjugate central volume from the paper's simulated terminal half-life, so this file is a compartmental reduction and not a PBPK model. The Child-Pugh hepatic-impairment scenarios are not encoded (see the vignette)."
  reference <- paste(
    "Kumar N, Dixit V, Burt H, Gill KL, Jones HM, Keirstead N, Su D, Toader D,",
    "Lowinger TB. PBPK Modeling to Predict Clinical Drug-Drug Interaction and",
    "Impact of Hepatic Impairment for an ADC With the Payload Auristatin",
    "F-Hydroxypropylamide. CPT Pharmacometrics Syst Pharmacol.",
    "2025;14(10):1661-1672. doi:10.1002/psp4.70088. PMID 40742885; PMCID",
    "PMC12521058.",
    sep = " "
  )
  vignette <- "Kumar_2025_upifitamabRilsodotin"

  units <- list(
    time = "h",
    dosing = "mg (AF-HPA payload equivalents)",
    concentration = "ng/mL"
  )

  # Both species are tracked as AF-HPA mass (MW 803 g/mol for conjugated and
  # unconjugated AF-HPA, Kumar 2025 Table 1), so payload release is a 1:1 mass
  # transfer. Doses are the conjugated AF-HPA payload dose from Kumar 2025
  # Supplementary Equation 1: ADC (antibody) dose / 145000 * mean DAR * 803.
  compartmentData <- list(
    central = list(
      analyte = "conjugated AF-HPA (upifitamab rilsodotin, as payload equivalents)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "conjugated AF-HPA (single adjusting compartment)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_afhpa = list(
      analyte = "unconjugated (released) AF-HPA",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    CONMED_ITRACONAZOLE = list(
      description = "Concomitant itraconazole (combined strong CYP3A4 and P-glycoprotein inhibitor) coadministration indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "1 = itraconazole 200 mg twice daily on day 1 and 200 mg once daily thereafter, with the ADC given on day 14 of 42 days of itraconazole dosing (Kumar 2025 Methods and Table 2 footnote), 0 = upifitamab rilsodotin alone. Itraconazole is on for the whole observation window, so the indicator is used as a step indicator. The effect is a multiplicative factor on unconjugated AF-HPA clearance only; conjugated AF-HPA is unaffected (Kumar 2025 Table S6 GMR 1.00). The magnitude is back-calculated from the published unconjugated AF-HPA AUC0-inf geometric mean ratio for the base itraconazole model (CYP3A4 and P-gp both inhibited): in this reduction the unconjugated AUC0-inf equals released payload / CL, so a CL multiplier of 1/GMR reproduces that ratio exactly, and the published Cmax ratio is held out as a validation gate. The single multiplier lumps the CYP3A4 and P-gp mechanisms; it is not an inhibition constant.",
      source_name = NA_character_
    )
  )

  population <- list(
    species = "human",
    n_subjects = 15L,
    n_studies = 1L,
    age_range = "53-85 years (model-development cohort); 57-79 years (multiple-dose verification cohort)",
    sex_female_pct = 87,
    disease_state = "Adults with advanced solid tumours in the phase 1/2 study MER-XMT-1536-1 (NCT03319628); simulated with the Simcyp Sim-Cancer population.",
    dose_range = "Upifitamab rilsodotin 3-52 mg/m2 (antibody) single IV dose (dose escalation) and 36 mg/m2 or 43 mg/m2 every 4 weeks (dose expansion). 36 mg/m2 corresponds to a conjugated AF-HPA payload dose of 2.0-2.12 mg/m2 depending on the lot drug-antibody ratio.",
    notes = "Top-down model development used the 36 mg/m2 cycle-1 dose-expansion cohort (15 patients, 87% female). Verification used the dose-escalation cohort (20-43 mg/m2) and the 36 mg/m2 q4w cycle-3 cohort (8 patients aged 57-79 years, 88% female). Drug-interaction simulations used the Simcyp healthy-volunteer population (10 trials of 10 subjects aged 20-50 years, 50% female). Body weight, height and body surface area of the cohorts are not reported; the reduction is calibrated at a reference body surface area of 1.71 m2 and a reference body weight of 70 kg (see the vignette)."
  )

  ini({
    # Conjugated AF-HPA disposition (minimal PBPK with a single adjusting
    # compartment, SAC == peripheral1). kin / kout move drug MASS into and out
    # of the SAC (Kumar 2025 Figure 2 legend), so they are k12 / k21.
    lk12 <- log(0.00475); label("Conjugated AF-HPA transfer rate into the single adjusting compartment, kin (1/h)") # Kumar 2025 Table 1: kin 0.00475 1/h (Supplementary Section 10, Simcyp parameter estimation)
    lk21 <- log(0.03); label("Conjugated AF-HPA transfer rate out of the single adjusting compartment, kout (1/h)") # Kumar 2025 Table 1: kout 0.03 1/h (Supplementary Section 10: estimate 0.0150 x 2 by manual adjustment)

    # Conjugated AF-HPA clearance (payload release + biliary + renal). The
    # paper's iterative retrograde method matched the human clearance-pathway
    # fractions to the rat excretion fractions (Table S1), so total CL =
    # CL_R / fe_renal. The Table 1 CL_R input is the healthy-renal-function
    # value that Simcyp scales by 0.63 in the Sim-Cancer population
    # (Supplementary Section 12).
    lcl <- log(0.00302 * 0.63 / 0.0690); label("Conjugated AF-HPA total clearance, Sim-Cancer population (L/h)") # derived: Kumar 2025 Table 1/Table S2 CL_R 0.00302 L/h x 0.63 (Suppl. Section 12) / renal fe 0.0690 (Table S1) = 0.02757 L/h

    # Conjugated AF-HPA central volume. Simcyp's systemic volume is a platform
    # output, so it is solved from the paper's simulated geometric-mean terminal
    # half-life 107 h (Table S3) with the printed k12 / k21 and the derived CL:
    # k10 = beta * (k12 + k21 - beta) / (k21 - beta), beta = ln2 / 107 h.
    lvc <- log(3.5414); label("Conjugated AF-HPA central volume (L)") # derived: Kumar 2025 Table S3 simulated terminal t1/2 107 h; cross-check: Table 1 Vss 0.0686 L/kg x 70 kg = 4.80 L minus the liver volume

    # Payload release. The fraction of conjugated AF-HPA clearance that releases
    # unconjugated AF-HPA is the recovery-corrected rat fraction excreted by
    # 'Metabolism - other routes' for conjugated AF-HPA, which the paper assigns
    # to CL_Payload-Release (Methods; Supplementary Section 4).
    fm_afhpa <- fixed(0.525); label("Fraction of conjugated AF-HPA clearance that releases unconjugated AF-HPA (unitless)") # Kumar 2025 Table S1: conjugated AF-HPA 'Metabolism - other routes' 0.525

    # Unconjugated AF-HPA disposition. Total CL from the printed renal
    # clearance and the rat renal fraction excreted, as for the conjugate.
    lcl_afhpa <- log(0.0361 * 0.63 / 0.0263); label("Unconjugated AF-HPA total clearance, Sim-Cancer population (L/h)") # derived: Kumar 2025 Table 1/Table S2 CL_R 0.0361 L/h x 0.63 (Suppl. Section 12) / renal fe 0.0263 (Table S1) = 0.8648 L/h
    lvc_afhpa <- log(0.524 * 70); label("Unconjugated AF-HPA volume of distribution (L)") # Kumar 2025 Supplementary Section 11: full-PBPK Vss 0.524 L/kg (Kp scalar 0.022) x 70 kg reference weight

    # Itraconazole (base model, CYP3A4 and P-gp inhibited) on unconjugated
    # AF-HPA clearance: CL multiplier = 1 / published AUC0-inf GMR.
    e_itraconazole_cl_afhpa <- fixed(1 / 1.34); label("Multiplier of itraconazole coadministration on unconjugated AF-HPA CL, power form on CONMED_ITRACONAZOLE (unitless)") # derived: Kumar 2025 Table 2 base model AUC0-inf GMR 1.34

    # Residual error. Kumar 2025 is a deterministic PBPK simulation study and
    # reports no residual-error model; the CV% values in Tables S3-S5 are
    # observed between-subject spread, not residual error.
    propSd <- fixed(0); label("Proportional residual error on Cc (fraction), not reported by the source") # not reported in Kumar 2025; see vignette
    propSd_afhpa <- fixed(0); label("Proportional residual error on Cc_afhpa (fraction), not reported by the source") # not reported in Kumar 2025; see vignette
  })

  model({
    # 1. Individual parameters
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    cl <- exp(lcl)
    vc <- exp(lvc)
    cl_afhpa <- exp(lcl_afhpa) * e_itraconazole_cl_afhpa^CONMED_ITRACONAZOLE
    vc_afhpa <- exp(lvc_afhpa)

    # 2. Micro-constants
    kel <- cl / vc
    kel_afhpa <- cl_afhpa / vc_afhpa

    # 3. ODE system. Payload release, biliary and renal clearance of conjugated
    # AF-HPA all act on the systemic (central) compartment; the fraction
    # fm_afhpa of that flux appears as unconjugated AF-HPA (1:1 by mass).
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(central_afhpa) <- fm_afhpa * kel * central - kel_afhpa * central_afhpa

    # 4. Observations (mg/L x 1000 = ng/mL)
    Cc <- 1000 * central / vc
    Cc_afhpa <- 1000 * central_afhpa / vc_afhpa
    Cc ~ prop(propSd)
    Cc_afhpa ~ prop(propSd_afhpa)
  })
}
