Cao_2013_mepolizumab <- function() {
  description <- "Second-generation minimal physiologically-based PK (mPBPK) model for mepolizumab in adults (Cao 2013 Model A; clearance from plasma)"
  reference <- "Cao Y, Balthasar JP, Jusko WJ. Second-generation minimal physiologically-based pharmacokinetic model for monoclonal antibodies. J Pharmacokinet Pharmacodyn. 2013 Oct;40(5):597-607. doi:10.1007/s10928-013-9332-2. Mepolizumab plasma data digitized from Smith DA et al. Clin Pharmacokinet. 2011;50(4):215-227 (PMID 21348536)."
  vignette <- "Cao_2013_mepolizumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    n_subjects     = NA_integer_,
    n_studies      = 1,
    age_range      = "adults (per Smith 2011 source study)",
    weight_range   = "70 kg reference body weight (Cao 2013 Table 2 footnote)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA,
    disease_state  = "Mepolizumab anti-IL-5 humanized IgG1 (population per Smith 2011 PK/PD analysis).",
    dose_range     = "IV; see Cao 2013 Figure 5 mepolizumab panel.",
    regions        = NA,
    notes          = "Cao 2013 Table 2, Model A. Parameters fit by Cao et al. to plasma concentration profiles digitized from Smith DA et al. Clin Pharmacokinet 2011;50:215-227 (PMID 21348536). sigma1 was fixed at 0.950 in Cao 2013 (not estimated; CV reported as not applicable in Table 2)."
  )

  ini({
    sigma1 <- 0.950; label("Vascular reflection coefficient for tight tissues (unitless; fixed at 0.950 in Cao 2013)")  # Cao 2013 Table 2 (Model A): 0.950, fixed (footnote c "Not applicable")
    sigma2 <- 0.750; label("Vascular reflection coefficient for leaky tissues (unitless)")                                # Cao 2013 Table 2 (Model A): 0.750 (CV 1.48%)
    lclp   <- log(0.20424); label("Plasma clearance (CLp, L/day)")                                                        # Cao 2013 Table 2 (Model A): CLp = 0.00851 L/hr (CV 1.50%) = 0.20424 L/day
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

    clp <- exp(lclp)

    cp     <- plasma / vplasma
    ctight <- tight  / vtight
    cleaky <- leaky  / vleaky
    clymph <- lymph  / vlymph

    d/dt(plasma) <- clymph * lymphflow -
                    cp * l1 * (1 - sigma1) -
                    cp * l2 * (1 - sigma2) -
                    clp * cp
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
