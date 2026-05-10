Cao_2013_mab8C2 <- function() {
  description <- "Preclinical mPBPK model for the murine anti-topotecan IgG1 mAb 8C2 in mice (Cao 2013 Model A; clearance from plasma)"
  reference <- "Cao Y, Balthasar JP, Jusko WJ. Second-generation minimal physiologically-based pharmacokinetic model for monoclonal antibodies. J Pharmacokinet Pharmacodyn. 2013 Oct;40(5):597-607. doi:10.1007/s10928-013-9332-2. 8C2 plasma data from Abuqayyas L, Balthasar JP. Int J Pharm. 2012;439(1-2):8-16 (PMID 23018115)."
  vignette <- "Cao_2013_mab8C2"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    species        = "Mouse (Mus musculus); preclinical PK study (Cao 2013 mPBPK fit; underlying data Abuqayyas 2012).",
    n_subjects     = NA_integer_,
    n_studies      = 1,
    weight_range   = "20 g body weight (Cao 2013 Table 1 footnote 'Assumed 20 g body weight')",
    disease_state  = "Healthy mice (wild-type and FcgammaR knock-out cohorts in source data); preclinical PK characterisation of 8C2, a murine IgG1 anti-topotecan mAb used as a non-binding carrier antibody.",
    dose_range     = "0.04, 0.1, 0.4, 8 mg/kg single IV (Cao 2013 Figure 3)",
    regions        = "Preclinical (academic study, University at Buffalo)",
    scope_note     = "Preclinical-only mPBPK fit. Filed under inst/modeldb/pharmacokinetics/ rather than specificDrugs/ because nlmixr2lib's specificDrugs tier is reserved for human drugs.",
    notes          = "Cao 2013 Table 1, Model A. Parameters fit by Cao et al. to plasma concentration profiles from Abuqayyas L & Balthasar JP. Int J Pharm 2012;439:8-16 (PMID 23018115). The function name is mab8C2 because R identifiers cannot start with a digit; the antibody is referred to as 8C2 in the source publications. The 8C2 plasma data in Cao 2013 Figure 3 are reported in nM (molar units) rather than mg/L; this model returns Cc in mg/L (dose in mg, volume in L). Conversion: nM = (mg/L) / MW(mAb) * 1e6 with MW(mAb) ~ 1.5e5 g/mol."
  )

  ini({
    sigma1 <- 0.943; label("Vascular reflection coefficient for tight tissues (unitless)")  # Cao 2013 Table 1 (Model A): 0.943 (CV 30.7%)
    sigma2 <- 0.378; label("Vascular reflection coefficient for leaky tissues (unitless)")  # Cao 2013 Table 1 (Model A): 0.378 (CV 34.2%)
    lclp   <- log(1.260e-4); label("Plasma clearance (CLp, L/day)")                          # Cao 2013 Table 1 (Model A): CLp = 0.525e-5 L/hr (CV 46.5%) = 1.260e-4 L/day
  })

  model({
    # Mouse system parameters; Cao 2013 Table 1 footnote and Methods, for a 20 g BW mouse.
    # Vplasma = 0.85 mL = 0.00085 L; ISF = 4.35 mL = 0.00435 L;
    # total lymph flow = 0.12 mL/hr = 0.00288 L/day; sigmaL = 0.2; Kp = 0.8 for native IgG1.
    sigmal     <- 0.2
    kp         <- 0.8
    vplasma    <- 0.00085
    visf       <- 0.00435
    lymphflow  <- 0.00288

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
