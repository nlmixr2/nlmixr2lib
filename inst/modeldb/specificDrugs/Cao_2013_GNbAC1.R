Cao_2013_GNbAC1 <- function() {
  description <- "Second-generation minimal physiologically-based PK (mPBPK) model for GNbAC1 in adults (Cao 2013 Model A; clearance from plasma)"
  reference <- "Cao Y, Balthasar JP, Jusko WJ. Second-generation minimal physiologically-based pharmacokinetic model for monoclonal antibodies. J Pharmacokinet Pharmacodyn. 2013 Oct;40(5):597-607. doi:10.1007/s10928-013-9332-2. GNbAC1 plasma data digitized from Curtin F et al. Clin Ther. 2012;34(11):2268-2278 (PMID 23200102)."
  vignette <- "Cao_2013_GNbAC1"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    n_subjects     = NA_integer_,
    n_studies      = 1,
    age_range      = "healthy adult volunteers (first-in-human; per Curtin 2012)",
    weight_range   = "70 kg reference body weight (Cao 2013 Table 2 footnote)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA,
    disease_state  = "First-in-human single-dose escalation of GNbAC1 (humanized IgG4 against MSRV envelope protein, evaluated as a candidate therapy for Multiple Sclerosis-associated endogenous retrovirus).",
    dose_range     = "0.15, 0.6, 2, 6 mg/kg IV (Cao 2013 Figure 5 GNbAC1 panel)",
    regions        = NA,
    notes          = "Cao 2013 Table 2, Model A. Parameters fit by Cao et al. to plasma concentration profiles digitized from Curtin F et al. Clin Ther 2012;34:2268-2278 (PMID 23200102). GNbAC1 is the only IgG4 mAb in the Cao 2013 human cohort; native IgG4 has a lower available ISF fraction (Kp = 0.4 per Cao 2013 Methods, refs 22-23) than IgG1 (Kp = 0.8). Cao 2013 Table 2 does not state the Kp used per drug; this model file uses the IgG1 default Kp = 0.8 for consistency with the rest of the Cao 2013 cohort. Users wanting an IgG4-specific representation should refit sigma1/sigma2/CLp with Kp = 0.4."
  )

  ini({
    sigma1 <- 0.915; label("Vascular reflection coefficient for tight tissues (unitless)")  # Cao 2013 Table 2 (Model A): 0.915 (CV 5.21%)
    sigma2 <- 0.831; label("Vascular reflection coefficient for leaky tissues (unitless)")  # Cao 2013 Table 2 (Model A): 0.831 (CV 5.71%)
    lclp   <- log(0.17136); label("Plasma clearance (CLp, L/day)")                          # Cao 2013 Table 2 (Model A): CLp = 0.00714 L/hr (CV 2.93%) = 0.17136 L/day
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
