Xu_2020_jnj54175446_ro_monkey <- function() {
  description <- "Preclinical (rhesus monkey). Exposure-occupancy model for brain P2X7 receptor occupancy by JNJ-54175446 measured by [18F]-JNJ-64413739 PET after single intravenous doses (Xu 2020). Sigmoid Emax function of the average plasma JNJ-54175446 concentration over the PET scan (CAV, about 10 to 130 min after the dose). Typical-value model fitted to five occupancy measures from two animals; the source reports no variability terms. The output ro is the occupancy in percent. Xu 2020 projected human occupancy from this model by scaling EC50 by the ratio of the monkey to the human plasma fraction unbound (17.0 / 11.2 percent, giving 65.6 ng/mL) with Emax and gamma unchanged."
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
      notes = "Average of the plasma concentrations measured at about 10, 70 and 130 min after the intravenous dose, immediately before tracer injection, midway through and at the end of the scan (Xu 2020 Methods 'Monkey PET study design and sample collection', Eq. 3). The averaging window is the PET scan, not a dosing interval.",
      source_name = "C ave"
    )
  )

  population <- list(
    species = "rhesus monkey (Macaca mulatta)",
    n_subjects = 2L,
    n_studies = 1L,
    n_observations = "5 brain receptor occupancy measures",
    weight_range = "8.7-12.1 kg",
    sex_female_pct = 50,
    disease_state = "Healthy animals.",
    dose_range = "Single intravenous doses of 0.1, 0.4, 2.5 and 5.3 mg/kg JNJ-54175446 starting 10 min before tracer injection; the male received three and the female two blockade scans, at least two weeks apart.",
    notes = "Monkey PET study performed at Molecular NeuroImaging, New Haven, CT (Xu 2020 Methods). Plasma fraction unbound 17.0 percent in monkey versus 11.2 percent in human (Methods Eq. 4)."
  )

  ini({
    # Xu 2020 Table 2, 'Monkey estimates' column (estimate (RSE%)); Eq. 3.
    lemax <- log(56.0); label("Maximal brain P2X7 receptor occupancy (percent)") # Table 2 monkey 'Emax' = 56.0 (4.41); Results 'maximal brain occupancy (Emax) was 56.0%'
    lec50 <- log(43.2); label("Average plasma concentration giving half-maximal occupancy (ng/mL)") # Table 2 monkey 'EC50,RO' = 43.2 (9.02)
    lhill <- log(3.71); label("Sigmoidicity of the exposure-occupancy relationship (unitless)") # Table 2 monkey 'gamma' = 3.71 (37.2)
  })

  model({
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # Eq. 3: E = Emax * Cave^gamma / (EC50,RO^gamma + Cave^gamma)
    ro <- emax * CAV^hill / (ec50^hill + CAV^hill)
  })
}
