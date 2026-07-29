Cao_2013_PRO95780 <- function() {
  description <- "Second-generation minimal physiologically-based PK (mPBPK) model for PRO95780 (drozitumab) in adults (Cao 2013 Model A; clearance from plasma)"
  reference <- "Cao Y, Balthasar JP, Jusko WJ. Second-generation minimal physiologically-based pharmacokinetic model for monoclonal antibodies. J Pharmacokinet Pharmacodyn. 2013 Oct;40(5):597-607. doi:10.1007/s10928-013-9332-2. PRO95780 plasma data digitized from Camidge DR et al. Clin Cancer Res. 2010;16(4):1256-1263 (PMID 20145186)."
  vignette <- "Cao_2013_PRO95780"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    n_subjects     = NA_integer_,
    n_studies      = 1,
    age_range      = "adults (per Camidge 2010 source study)",
    weight_range   = "70 kg reference body weight (Cao 2013 Table 2 footnote)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA,
    disease_state  = "Advanced malignancies (PRO95780 / drozitumab, a death receptor 5 (DR5) agonistic human IgG1 antibody).",
    dose_range     = "4, 10, 15, 20 mg/kg IV (Cao 2013 Figure 5 PRO95780 panel)",
    regions        = NA,
    notes          = "Cao 2013 Table 2, Model A. Parameters fit by Cao et al. to plasma concentration profiles digitized from Camidge DR et al. Clin Cancer Res 2010;16:1256-1263 (PMID 20145186). The drug is also known as drozitumab in the literature."
  )

  ini({
    sigma_tight <- 0.984; label("Vascular reflection coefficient for tight tissues (unitless)")  # Cao 2013 Table 2 (Model A): 0.984 (CV 10.3%)
    sigma_leaky <- 0.638; label("Vascular reflection coefficient for leaky tissues (unitless)")  # Cao 2013 Table 2 (Model A): 0.638 (CV 19.4%)
    lcl   <- log(0.2976); label("Plasma clearance (CLp, L/day)")                           # Cao 2013 Table 2 (Model A): CLp = 0.0124 L/hr (CV 18.3%) = 0.2976 L/day
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
