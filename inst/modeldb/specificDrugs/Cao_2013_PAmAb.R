Cao_2013_PAmAb <- function() {
  description <- "Second-generation minimal physiologically-based PK (mPBPK) model for PAmAb in adults (Cao 2013 Model A; clearance from plasma)"
  reference <- "Cao Y, Balthasar JP, Jusko WJ. Second-generation minimal physiologically-based pharmacokinetic model for monoclonal antibodies. J Pharmacokinet Pharmacodyn. 2013 Oct;40(5):597-607. doi:10.1007/s10928-013-9332-2. PAmAb plasma data digitized from Subramanian GM et al. Clin Infect Dis. 2005;41(1):12-20 (PMID 15937757)."
  vignette <- "Cao_2013_PAmAb"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    n_subjects     = NA_integer_,
    n_studies      = 1,
    age_range      = "healthy adult volunteers (per Subramanian 2005 source study)",
    weight_range   = "70 kg reference body weight (Cao 2013 Table 2 footnote)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA,
    disease_state  = "Phase 1 single-dose study of PAmAb (a fully human monoclonal antibody against Bacillus anthracis protective antigen) in healthy volunteers (Subramanian 2005).",
    dose_range     = "1, 3, 10, 20, 40 mg/kg IV (Cao 2013 Figure 5 PAmAb panel)",
    regions        = NA,
    notes          = "Cao 2013 Table 2, Model A. Parameters fit by Cao et al. to plasma concentration profiles digitized from Subramanian GM et al. Clin Infect Dis 2005;41:12-20 (PMID 15937757). sigma1 was fixed at 0.950 in Cao 2013 (not estimated; CV reported as not applicable in Table 2)."
  )

  ini({
    sigma1 <- 0.950; label("Vascular reflection coefficient for tight tissues (unitless; fixed at 0.950 in Cao 2013)")  # Cao 2013 Table 2 (Model A): 0.950, fixed (footnote c "Not applicable")
    sigma2 <- 0.779; label("Vascular reflection coefficient for leaky tissues (unitless)")                                # Cao 2013 Table 2 (Model A): 0.779 (CV 5.24%)
    lcl   <- log(0.20808); label("Plasma clearance (CLp, L/day)")                                                        # Cao 2013 Table 2 (Model A): CLp = 0.00867 L/hr (CV 2.29%) = 0.20808 L/day
  })

  model({
    sigmal     <- 0.2
    kp         <- 0.8
    vplasma    <- 2.6
    visf       <- 15.6
    lymphflow  <- 2.9

    vtight  <- 0.65 * visf * kp
    vleaky  <- 0.35 * visf * kp
    l1      <- 0.33 * lymphflow
    l2      <- 0.67 * lymphflow
    vlymph  <- vplasma

    cl <- exp(lcl)

    cp     <- plasma / vplasma
    ctight <- tight  / vtight
    cleaky <- leaky  / vleaky
    clymph <- lymph  / vlymph

    d/dt(plasma) <- clymph * lymphflow -
                    cp * l1 * (1 - sigma1) -
                    cp * l2 * (1 - sigma2) -
                    cl * cp
    d/dt(tight)  <- l1 * (1 - sigma1) * cp -
                    l1 * (1 - sigmal) * ctight
    d/dt(leaky)  <- l2 * (1 - sigma2) * cp -
                    l2 * (1 - sigmal) * cleaky
    d/dt(lymph)  <- l1 * (1 - sigmal) * ctight +
                    l2 * (1 - sigmal) * cleaky -
                    clymph * lymphflow

    Cc <- cp
  })
}
