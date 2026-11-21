Cleary_2026_methylprednisolone <- function() {
  description <- paste(
    "PBPK (GastroPlus 8.3.3) long-acting intramuscular methylprednisolone",
    "acetate paper, compartmental disposition layer. Two-compartment",
    "intravenous methylprednisolone model in a 70 kg adult with",
    "body-weight-proportional central volume and clearance, fitted with the",
    "GastroPlus PKPlus module to published mean plasma concentrations after",
    "an 80 mg intravenous bolus of methylprednisolone sodium succinate. The",
    "paper's mechanistic intramuscular depot (particle dissolution and",
    "injection-site absorption) and its simulated in vivo dissolution",
    "profiles are not included (see the vignette). Deterministic",
    "typical-value model with no between-subject or residual variability.",
    sep = " "
  )
  reference <- paste(
    "Cleary H, Fotaki N, Persoons T, D'Arcy DM. Using PBPK to Simulate Target",
    "Biopredictive Dissolution Profiles for Long-Acting Injectables - Where",
    "to Begin With Critical Bioavailability Attributes? CPT Pharmacometrics",
    "Syst Pharmacol. 2026;15(3):e70212. doi:10.1002/psp4.70212.",
    "Disposition parameters are from Table 1; the peripheral-to-central rate",
    "constant K21 is not printed and was back-solved from the simulated",
    "concentration at 0.18 h in Supplementary Table S1.",
    sep = " "
  )
  vignette <- "Cleary_2026_methylprednisolone"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  dosing <- c("central")

  compartmentData <- list(
    central = list(analyte = "methylprednisolone", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "methylprednisolone", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "GastroPlus compartmental parameters are entered per kg and scale",
        "linearly with the simulated patient's weight; the paper simulated",
        "the GastroPlus default body weight of 70 kg (Table 1)."
      ),
      source_name = "Body weight"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 1L,
    age_range = "Adults (literature mean profile; ages not reported in this paper)",
    weight_range = "70 kg reference individual (GastroPlus default, Table 1)",
    disease_state = "Not reported in this paper (literature intravenous study, Al-Habet and Rogers 1989)",
    dose_range = "80 mg intravenous bolus of methylprednisolone sodium succinate (this model); 80 mg and 120 mg intramuscular methylprednisolone acetate (paper, not included)",
    regions = "Not reported",
    notes = paste(
      "No individual data were fitted. The 1-, 2- and 3-compartment models",
      "were fitted with the GastroPlus PKPlus module to a single digitised",
      "mean plasma methylprednisolone profile after an 80 mg intravenous",
      "bolus of methylprednisolone sodium succinate (in vivo dataset A,",
      "Al-Habet and Rogers 1989; Supplementary S2 Table S2_1, 18 time points",
      "to 10 h); the 2-compartment model was selected (AIC -89.56). The",
      "intramuscular layer of the paper was evaluated against mean profiles",
      "after 80 mg methylprednisolone acetate (dataset B, FDA ANDA review)",
      "and 120 mg on days 0 and 15 in 48 patients with active rheumatoid",
      "arthritis (dataset C, Metselaar 2022).",
      sep = " "
    )
  )

  ini({
    lvc <- fixed(log(0.84106)); label("Central volume per kg body weight (L/kg)") # Table 1 'Volume of distribution (Vd)' = 0.84106 L/kg (PKPlus simulation, step 1)
    lk12 <- fixed(log(45.107)); label("Central-to-peripheral rate constant K12 (1/h)") # Table 1 'K12' = 45.107 1/h (PKPlus simulation, step 1)
    lk21 <- fixed(log(89.14)); label("Peripheral-to-central rate constant K21 (1/h)") # not printed; back-solved so that C(0.18 h) after 80 mg IV in 70 kg = Table S1 simulated 'Cmax' 842 ng/mL (rounding interval 88.98-89.31)
    lcl <- fixed(log(0.46751)); label("Systemic clearance per kg body weight (L/h/kg)") # Table 1 'Clearance (Cl)' = 0.46751 L/h/kg (PKPlus simulation, step 1)
  })

  model({
    vc <- exp(lvc) * WT
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    cl <- exp(lcl) * WT
    kel <- cl / vc

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Plasma methylprednisolone concentration, mg/L (= ug/mL).
    Cc <- central / vc
  })
}
