Cao_2013_visilizumab <- function() {
  description <- "Second-generation minimal physiologically-based PK (mPBPK) model for visilizumab in adults (Cao 2013 Model A; clearance from plasma)"
  reference <- "Cao Y, Balthasar JP, Jusko WJ. Second-generation minimal physiologically-based pharmacokinetic model for monoclonal antibodies. J Pharmacokinet Pharmacodyn. 2013 Oct;40(5):597-607. doi:10.1007/s10928-013-9332-2. Visilizumab plasma data digitized from Carpenter PA et al. Blood. 2002;99(8):2712-2719 (PMID 11929757)."
  vignette <- "Cao_2013_visilizumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    n_subjects     = NA_integer_,
    n_studies      = 1,
    age_range      = "adults (per Carpenter 2002 source study)",
    weight_range   = "70 kg reference body weight (Cao 2013 Table 2 footnote)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA,
    disease_state  = "Steroid-refractory acute graft-versus-host disease (visilizumab humanized non-FcR-binding anti-CD3 IgG2).",
    dose_range     = "3 mg/m^2 IV (Cao 2013 Figure 5 visilizumab panel)",
    regions        = NA,
    notes          = "Cao 2013 Table 2, Model A. Parameters fit by Cao et al. to plasma concentration profiles digitized from Carpenter PA et al. Blood 2002;99:2712-2719 (PMID 11929757). Doses are reported in mg/m^2 (BSA-normalized) in Cao 2013 Figure 5; this model treats dose in mg, so external scaling of mg/m^2 -> mg requires the subject body surface area."
  )

  ini({
    sigma_tight <- 0.949; label("Vascular reflection coefficient for tight tissues (unitless)")  # Cao 2013 Table 2 (Model A): 0.949 (CV 3.63%)
    sigma_leaky <- 0.834; label("Vascular reflection coefficient for leaky tissues (unitless)")  # Cao 2013 Table 2 (Model A): 0.834 (CV 5.44%)
    lcl   <- log(0.3648); label("Plasma clearance (CLp, L/day)")                           # Cao 2013 Table 2 (Model A): CLp = 0.0152 L/hr (CV 3.02%) = 0.3648 L/day
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
