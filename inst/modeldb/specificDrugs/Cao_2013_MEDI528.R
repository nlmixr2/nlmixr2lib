Cao_2013_MEDI528 <- function() {
  description <- "Second-generation minimal physiologically-based PK (mPBPK) model for MEDI-528 in adults (Cao 2013 Model A; clearance from plasma)"
  reference <- "Cao Y, Balthasar JP, Jusko WJ. Second-generation minimal physiologically-based pharmacokinetic model for monoclonal antibodies. J Pharmacokinet Pharmacodyn. 2013 Oct;40(5):597-607. doi:10.1007/s10928-013-9332-2. MEDI-528 plasma data digitized from White B et al. Clin Ther. 2009;31(4):728-740 (PMID 19446146)."
  vignette <- "Cao_2013_MEDI528"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    n_subjects     = NA_integer_,
    n_studies      = 1,
    age_range      = "healthy adult volunteers (per White 2009 source study)",
    weight_range   = "70 kg reference body weight (Cao 2013 Table 2 footnote)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA,
    disease_state  = "Two first-in-human, open-label, phase I dose-escalation safety trials of MEDI-528, a humanized monoclonal antibody against interleukin-9, in healthy adult volunteers (White 2009).",
    dose_range     = "0.3, 1, 3, 9 mg/kg IV (Cao 2013 Figure 5 MEDI-528 panel)",
    regions        = NA,
    notes          = "Cao 2013 Table 2, Model A. Parameters fit by Cao et al. to plasma concentration profiles digitized from White B et al. Clin Ther 2009;31:728-740 (PMID 19446146). The drug is referred to as MEDI-528 (with hyphen) in the source publications; the model file name uses MEDI528 because R identifiers cannot contain hyphens."
  )

  ini({
    sigma1 <- 0.987; label("Vascular reflection coefficient for tight tissues (unitless)")  # Cao 2013 Table 2 (Model A): 0.987 (CV 6.17%)
    sigma2 <- 0.754; label("Vascular reflection coefficient for leaky tissues (unitless)")  # Cao 2013 Table 2 (Model A): 0.754 (CV 7.00%)
    lclp   <- log(0.12888); label("Plasma clearance (CLp, L/day)")                          # Cao 2013 Table 2 (Model A): CLp = 0.00537 L/hr (CV 6.88%) = 0.12888 L/day
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
