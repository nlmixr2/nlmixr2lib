Cao_2013_gevokizumab <- function() {
  description <- "Second-generation minimal physiologically-based PK (mPBPK) model for gevokizumab in adults (Cao 2013 Model A; clearance from plasma)"
  reference <- "Cao Y, Balthasar JP, Jusko WJ. Second-generation minimal physiologically-based pharmacokinetic model for monoclonal antibodies. J Pharmacokinet Pharmacodyn. 2013 Oct;40(5):597-607. doi:10.1007/s10928-013-9332-2. Gevokizumab plasma data digitized from Cavelti-Weder C et al. Diabetes Care. 2012;35(8):1654-1662 (PMID 22699287)."
  vignette <- "Cao_2013_gevokizumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    n_subjects     = NA_integer_,
    n_studies      = 1,
    age_range      = "adults (per Cavelti-Weder 2012 source study)",
    weight_range   = "70 kg reference body weight (Cao 2013 Table 2 footnote)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA,
    disease_state  = "Type 2 diabetes mellitus (gevokizumab anti-IL-1beta humanized IgG2).",
    dose_range     = "0.01, 0.03, 0.1, 0.3, 1, 3 mg/kg IV (Cao 2013 Figure 5 gevokizumab panel)",
    regions        = NA,
    notes          = "Cao 2013 Table 2, Model A. Parameters fit by Cao et al. to plasma concentration profiles digitized from Cavelti-Weder C et al. Diabetes Care 2012;35:1654-1662 (PMID 22699287). The fit used a typical-value structural mPBPK model in ADAPT 5; no IIV or covariate effects were estimated."
  )

  ini({
    sigma_tight <- 0.931; label("Vascular reflection coefficient for tight tissues (unitless)")  # Cao 2013 Table 2 (Model A): 0.931 (CV 2.58%)
    sigma_leaky <- 0.837; label("Vascular reflection coefficient for leaky tissues (unitless)")  # Cao 2013 Table 2 (Model A): 0.837 (CV 2.63%)
    lcl   <- log(0.16032); label("Plasma clearance (CLp, L/day)")                          # Cao 2013 Table 2 (Model A): CLp = 0.00668 L/hr (CV 1.87%) = 0.16032 L/day
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
                    cp * l1 * (1 - sigma_tight) -
                    cp * l2 * (1 - sigma_leaky) -
                    cl * cp
    d/dt(tight)  <- l1 * (1 - sigma_tight) * cp -
                    l1 * (1 - sigmal) * ctight
    d/dt(leaky)  <- l2 * (1 - sigma_leaky) * cp -
                    l2 * (1 - sigmal) * cleaky
    d/dt(lymph)  <- l1 * (1 - sigmal) * ctight +
                    l2 * (1 - sigmal) * cleaky -
                    clymph * lymphflow

    Cc <- cp
  })
}
