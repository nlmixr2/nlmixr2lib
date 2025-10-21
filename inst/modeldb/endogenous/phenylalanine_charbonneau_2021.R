phenylalanine_charbonneau_2021 <- function() {
  description <- "Phenylalanine model for absorption and metabolism in healthy subjects and patients with PKU"
  reference <- "Charbonneau, M.R., Denney, W.S., Horvath, N.G. et al. Development of a mechanistic model to predict synthetic biotic activity in healthy volunteers and patients with phenylketonuria. Commun Biol 4, 898 (2021). https://doi.org/10.1038/s42003-021-02183-1"
  covariates <-
    list(
      WT = "Body weight in kg",
      time = "Time in hours",
      f_pah = "Fraction of healthy PAH activity (healthy = 1; PKU patient = 0 to 0.03)",
      bl_phe = "Typical values are about 0.075 mmol/L in healthy subjects and 1.18 mmol/L in patients"
    )
  # parameters come from Table 4 in paper
  ini({
    bl_phe <- 1.18; label("Baseline Phenylalanine (Phe) concentration (mmol/L)")
    bl_gut <- 0; label("Baseline Phe in the gut (mg)")

    ka_gut <- 0.25; label("Absorption rate from gut to plasma")
    v_npd <- 0.012; label("Rate of net protein breakdown ((mmol/L)/hr)")

    vmax_pah <- 0.9; label("Maximum rate of Phe breakdown by PAH in a healthy subject ((mmol/L)/hr)")
    f_pah <- 0; label("Fraction of healthy PAH activity (PKU patient = 0 to 0.02)")
    km_pah <- 0.51; label("Michaelis-Menten constant for Phe with PAH (mmol/L)")
    kact_pah <- 0.54; label("Phe activation constant for PAH")

    vmax_trans <- 0.063; label("Maximum rate of Phe breakdown by transaminase ((mmol/L)/hr)")
    km_trans <- 1.37; label("Michaelis-Menten constant for Phe with transaminase (mmol/L)")

    cl_renal <- 5.696e-4; label("Renal clearance of Phe per body weight ((L/kg)/hr)")

    vd <- 0.5; label("Body-weight normalized volume distribution of Phe (L/kg)")
  })
  model({
    # Molecular weight of Phe (g/mol)
    mw_phe <- 165.19
    # Unit conversion adjustment from Gut to Plasma concentrations (mmol/L)/mg
    f_gut_plasma <- 1/(mw_phe * vd_phe * WT)

    v_pah <- vmax_pah*f_pah / (1 + km_pah/phe + km_pah*kact_pah/(phe^2)) # units: (mmol/L)/hr
    v_trans <- vmax_trans / (1 + km_trans/phe) # units: (mmol/L)/hr
    v_renal <- phe * cl_renal * vd # units: (mmol/L)/hr

    d/dt(gut) <- -ka_gut*gut
    d/dt(phe) <- ka_gut*gut*f_gut_plasma + v_npd - v_pah - v_trans - v_renal
    gut(0) <- bl_gut
    phe(0) <- bl_phe
    phe_umol <- phe * 1000 # units: umol/L (more commonly used in clinical laboratories)

    # The following is an augmentation of the model reported in the paper.  It
    # indicates the approximate daily Phe intake (in mg) to achieve
    # steady-state.
    daily_phe_intake <- 24 * vd * (v_pah + v_trans + v_renal - v_npd) / f_gut_plasma
  })
}
