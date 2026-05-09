Cao_2013_adecatumumab <- function() {
  description <- "Second-generation minimal physiologically-based PK (mPBPK) model for adecatumumab in adults (Cao 2013 Model A; clearance from plasma)"
  reference <- "Cao Y, Balthasar JP, Jusko WJ. Second-generation minimal physiologically-based pharmacokinetic model for monoclonal antibodies. J Pharmacokinet Pharmacodyn. 2013 Oct;40(5):597-607. doi:10.1007/s10928-013-9332-2. Adecatumumab plasma data digitized from Oberneder R et al. Eur J Cancer. 2006;42(15):2530-2538 (PMID 16930989)."
  vignette <- "Cao_2013_adecatumumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    n_subjects     = NA_integer_,
    n_studies      = 1,
    age_range      = "adults (per Oberneder 2006 source study)",
    weight_range   = "70 kg reference body weight (Cao 2013 Table 2 footnote)",
    sex_female_pct = 0,
    race_ethnicity = NA,
    disease_state  = "Hormone-refractory prostate cancer (adecatumumab anti-EpCAM human IgG1).",
    dose_range     = "Multiple-dose IV (per Oberneder 2006); see Cao 2013 Figure 5 adecatumumab panel.",
    regions        = NA,
    notes          = "Cao 2013 Table 2, Model A. Parameters fit by Cao et al. to plasma concentration profiles digitized from Oberneder R et al. Eur J Cancer 2006;42:2530-2538 (PMID 16930989). Adecatumumab has the highest CLp in the Cao 2013 human cohort (~0.030 L/hr), corresponding to its short ~15-day half-life. The fit used a typical-value structural mPBPK model in ADAPT 5; no IIV or covariate effects were estimated."
  )

  ini({
    sigma1 <- 0.883; label("Vascular reflection coefficient for tight tissues (unitless)")  # Cao 2013 Table 2 (Model A): 0.883 (CV 9.64%)
    sigma2 <- 0.524; label("Vascular reflection coefficient for leaky tissues (unitless)")  # Cao 2013 Table 2 (Model A): 0.524 (CV 7.37%)
    lclp   <- log(0.720); label("Plasma clearance (CLp, L/day)")                            # Cao 2013 Table 2 (Model A): CLp = 0.0300 L/hr (CV 3.19%) = 0.720 L/day
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
