Wahl_2017_tenofovir_frt_tfvdp <- function() {
  description <- paste0(
    "Preclinical (mouse). Static landmark exposure-response model for the ",
    "fraction of BLT-humanized mice protected from vaginal HIV-1 ",
    "acquisition by systemic tenofovir disoproxil fumarate (TDF) ",
    "pre-exposure prophylaxis, driven by the female reproductive tract ",
    "(FRT) tissue concentration of the active metabolite tenofovir ",
    "diphosphate (TFVdp) (Wahl 2017). This is the FRT-TFVdp arm of a ",
    "five-model PK-PD comparison (dose, plasma TFV, FRT TFV, FRT TFVdp, ",
    "cervicovaginal lavage TFV). The protected fraction is a fractional ",
    "sigmoidal Hill (Emax) function of the FRT TFVdp concentration, E0 + ",
    "(Emax - E0) * C^H / (EC50^H + C^H) with Emax fixed to 1 and the Hill ",
    "coefficient H FIXED to 1 (Wahl 2017 Eqn 2, Table 2 'FRT TFVdp' ",
    "column). There is no PK layer and no ODE: the FRT TFVdp ",
    "concentration enters as the CONC_TFVDP_FRT_FMOLG covariate column ",
    "(on the BLT scale), and the probability of protection is emitted as ",
    "the landmark output prob_protection. Because FRT TFVdp did not ",
    "increase monotonically with dose (it was highest at 140 mg/kg), the ",
    "model predicts no further efficacy gain at 300 mg/kg. The model was ",
    "fitted to BLT-equivalent FRT TFVdp obtained by scaling the parallel ",
    "BALB/c FRT TFVdp by 0.62. Fitted by maximum likelihood in ADAPT5."
  )
  reference <- paste(
    "Wahl A, Ho PT, Denton PW, Garrett KL, Hudgens MG, Swartz G, O'Neill C,",
    "Veronese F, Kashuba AD, Garcia JV.",
    "Predicting HIV Pre-exposure Prophylaxis Efficacy for Women using a",
    "Preclinical Pharmacokinetic-Pharmacodynamic In Vivo Model.",
    "Sci Rep. 2017 Jan 31;7:41098. doi:10.1038/srep41098. PMCID PMC5286499.",
    sep = " "
  )
  vignette <- "Wahl_2017_tenofovir_prep"

  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; the FRT TFVdp concentration enters as the CONC_TFVDP_FRT_FMOLG covariate column)",
    concentration = "CONC_TFVDP_FRT_FMOLG in fmol/g of FRT tissue; prob_protection is a probability (0-1)"
  )

  covariateData <- list(
    CONC_TFVDP_FRT_FMOLG = list(
      description = "Female reproductive tract (FRT) tissue tenofovir-diphosphate (TFVdp) concentration driving the protected-fraction model, on the BLT-mouse scale, in fmol/g.",
      units = "fmol/g",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Landmark FRT (uterus + cervix + vagina) TFVdp concentration at ",
        "the time of vaginal HIV challenge. TFVdp is the intracellular ",
        "active anabolite of tenofovir. The model was FITTED on ",
        "BLT-equivalent concentrations, obtained by multiplying the ",
        "parallel BALB/c FRT TFVdp by the BLT:BALB/c median ratio 0.62 ",
        "(Wahl 2017 Results 'Comparative TFV and TFVdp concentrations', ",
        "and Fig. 4c which scales TFVdp by 0.62). This covariate carries ",
        "the BLT-scale concentration: to reproduce a dose cohort from the ",
        "BALB/c measurements in Table S4, multiply the BALB/c median by ",
        "0.62. BALB/c median FRT TFVdp by TDF dose (Table S4): 20 mg/kg ",
        "1193083, 50 mg/kg 1146214, 140 mg/kg 8167641, 300 mg/kg 6828724 ",
        "fmol/g (note the 140 mg/kg median exceeds the 300 mg/kg median, ",
        "which is why predicted efficacy does not rise from 140 to 300 ",
        "mg/kg). A value of 0 gives prob_protection = E0. Member of the ",
        "applied-concentration CONC_<DRUG>_<UNITS> family; distinct from ",
        "the registered CONC_TFVDP_FMOLMC (mucosal-mononuclear-cell TFVdp ",
        "in fmol per million cells, Jayachandran 2021) by both matrix ",
        "(whole FRT tissue) and unit basis (per gram of tissue), which is ",
        "why this matrix-and-unit-qualified sibling name is registered ",
        "rather than reusing that canonical."
      ),
      source_name = "FRT TFVdp (fmol/g)"
    )
  )

  population <- list(
    species = "mouse (BLT humanized, protection endpoint; FRT TFVdp exposures from parallel BALB/c mice, scaled to the BLT level)",
    n_subjects = "Protection from 78 BLT mice across four TDF arms; FRT TFVdp from BALB/c mice (n=8 per dose)",
    n_studies = 1L,
    n_observations = "Binary protection outcome per BLT mouse related to the BLT-scaled cohort-median FRT TFVdp at the matching TDF dose",
    age_range = "Not reported (young adult laboratory mice)",
    sex_female_pct = 100,
    race_ethnicity = "n/a (mouse)",
    disease_state = "HIV-1 pre-exposure prophylaxis challenge model. BLT humanized NSG mice (HIV-susceptible) for protection; BALB/c mice as the FRT-tissue exposure surrogate.",
    dose_range = "TDF 20, 50, 140, or 300 mg/kg intraperitoneally once daily",
    regions = "USA (University of North Carolina at Chapel Hill)",
    notes = paste0(
      "FRT TFVdp was the active-metabolite exposure; it did not increase ",
      "monotonically with dose (highest median at 140 mg/kg), so the ",
      "model predicts no efficacy increase from 140 to 300 mg/kg, with ",
      "considerable overlap in the concentrations providing 84-100% ",
      "protection (Results 'PK-PD analysis of tenofovir')."
    )
  )

  ini({
    # ================================================================
    # Wahl 2017 Eqn 2 (fractional sigmoidal Hill / Emax) fitted to the
    # BLT-mouse protection data against the BLT-scale FRT TFVdp
    # concentration, final ADAPT5 estimates from Table 2, 'FRT TFVdp
    # (fmol/g)' column. Values are mean (relative standard error [%SE]).
    # E0 natural scale; EC50 log-transformed. Emax fixed to 1 (all
    # models). The Hill coefficient H is FIXED to 1 for the tissue and
    # CVL models (Table 2 reports H = NA; estimating it gave false
    # minima, Methods 'PD model development').
    # ================================================================
    e0 <- 0.2622
    label("Baseline protected fraction at zero FRT TFVdp E0 (probability)")  # Wahl 2017 Table 2, FRT TFVdp column E0 = 0.2622 (%SE 7.405)
    lemax <- fixed(log(1))
    label("Log maximum protected fraction Emax (unitless); Wahl 2017 held Emax at 1 for all models")  # Wahl 2017 Methods 'PD model development': 'Emax was fixed to one for all models'
    lec50 <- log(6.052e5)
    label("Log FRT TFVdp concentration giving half-maximal protection EC50 (fmol/g)")  # Wahl 2017 Table 2, FRT TFVdp column EC50 = 6.052e5 fmol/g (%SE 22.75)
    lhill <- fixed(log(1))
    label("Log Hill coefficient H of the FRT TFVdp concentration-response curve (unitless)")  # Wahl 2017 Table 2 reports H = NA; Methods: H fixed to 1 for tissue and CVL models (estimating it gave false minima)
    # ================================================================
    # Placeholder additive residual for the deterministic landmark
    # probability; Wahl 2017 used a proportional error whose final
    # estimate it did not report. Not from source.
    # ================================================================
    addSd_prob_protection <- fixed(0.001)
    label("Placeholder additive residual SD on the typical-value protected fraction; Wahl 2017 used a proportional error whose final estimate it did not report")  # not from source; see vignette Assumptions and deviations
  })

  model({
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # Fractional sigmoidal Hill (Emax) model, Wahl 2017 Eqn 2:
    # E = E0 + (Emax - E0) * DV^H / (EC50^H + DV^H), Emax = 1, H = 1.
    # Driver DV is the BLT-scale FRT TFVdp concentration (fmol/g).
    prob_protection <-
      e0 + (emax - e0) * CONC_TFVDP_FRT_FMOLG^hill /
        (ec50^hill + CONC_TFVDP_FRT_FMOLG^hill)

    prob_protection ~ add(addSd_prob_protection)
  })
}
