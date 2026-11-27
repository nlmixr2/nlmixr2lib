Xu_2020_jnj54175446_ro <- function() {
  description <- "Exposure-occupancy model for brain P2X7 receptor occupancy by JNJ-54175446 in healthy men, measured by [18F]-JNJ-64413739 PET after single oral doses (Xu 2020). Sigmoid Emax function of the average plasma JNJ-54175446 concentration over the PET scan (CAV, about 4 to 6 h after the dose), with a maximal occupancy below 100 percent. Typical-value model: the source reports no inter-individual variability and no residual error. The output ro is the occupancy in percent. CAV is supplied as a covariate, for example from the companion population PK model Xu_2020_jnj54175446."
  reference <- paste(
    "Xu Y, Miao X, Ravenstijn P, Hijzen A, Schmidt ME, Nandy P, Zhou H.",
    "Translational Model-Informed Dose Selection for a Human Positron",
    "Emission Tomography Imaging Study of JNJ-54175446, a P2X7 Receptor",
    "Antagonist. Clin Transl Sci. 2020;13(2):309-317.",
    "doi:10.1111/cts.12711.",
    sep = " "
  )
  vignette <- "Xu_2020_jnj54175446"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CAV = list(
      description = "Average plasma JNJ-54175446 concentration over the PET scan",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Average of the plasma concentrations measured at about 4, 4.75 and 6 h after the dose, immediately before tracer injection, midway through and at the end of the scan (Xu 2020 Methods 'Human study design and sample collection', Eq. 3). For simulation, take the time-average of Cc from Xu_2020_jnj54175446 over the 4-6 h window. The averaging window is the PET scan, not a dosing interval.",
      source_name = "C ave"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 8L,
    n_studies = 1L,
    n_observations = "16 brain receptor occupancy measures",
    sex_female_pct = 0,
    disease_state = "Healthy men.",
    dose_range = "Single oral doses of 5, 10, 20, 50, 150 and 300 mg JNJ-54175446 with a standard breakfast, two doses per subject at least one week apart.",
    studies = "Open-label [18F]-JNJ-64413739 PET study NCT03088644.",
    notes = "Each subject had a baseline scan and two post-dose scans (Xu 2020 Methods and Table S1). Occupancy was derived from the change in the regional distribution volume of the tracer; maximal occupancy below 100 percent reflects the non-displaceable part of the signal (Discussion). Age and weight of the PET cohort are not reported in Xu 2020."
  )

  ini({
    # Xu 2020 Table 2, 'Human estimates' column (estimate (RSE%)); Eq. 3.
    lemax <- log(78.3); label("Maximal brain P2X7 receptor occupancy (percent)") # Table 2 human 'Emax' = 78.3 (1.78)
    lec50 <- log(44.3); label("Average plasma concentration giving half-maximal occupancy (ng/mL)") # Table 2 human 'EC50,RO' = 44.3 (4.25)
    lhill <- log(1.43); label("Sigmoidicity of the exposure-occupancy relationship (unitless)") # Table 2 human 'gamma' = 1.43 (24.7)
  })

  model({
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # Eq. 3: E = Emax * Cave^gamma / (EC50,RO^gamma + Cave^gamma)
    ro <- emax * CAV^hill / (ec50^hill + CAV^hill)
  })
}
