Cao_2013_mab7E3 <- function() {
  description <- "Preclinical mPBPK model for the murine anti-platelet IgG1 mAb 7E3 in mice (Cao 2013 Model A; clearance from plasma)"
  reference <- "Cao Y, Balthasar JP, Jusko WJ. Second-generation minimal physiologically-based pharmacokinetic model for monoclonal antibodies. J Pharmacokinet Pharmacodyn. 2013 Oct;40(5):597-607. doi:10.1007/s10928-013-9332-2. 7E3 plasma data from Garg A, Balthasar JP. J Pharmacokinet Pharmacodyn. 2007;34(5):687-709 (PMID 17636457)."
  vignette <- "Cao_2013_mab7E3"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    species        = "Mouse (Mus musculus); preclinical PK study (Cao 2013 mPBPK fit; underlying data Garg 2007).",
    n_subjects     = NA_integer_,
    n_studies      = 1,
    weight_range   = "20 g body weight (Cao 2013 Table 1 footnote 'Assumed 20 g body weight')",
    disease_state  = "Healthy mice; preclinical PK characterisation of 7E3, the parent murine anti-GPIIb/IIIa IgG1 mAb (chimeric Fab fragment is approved as abciximab).",
    dose_range     = "8 mg/kg single IV (Cao 2013 Figure 2)",
    regions        = "Preclinical (academic study, University at Buffalo)",
    scope_note     = "Preclinical-only mPBPK fit. Filed under inst/modeldb/pharmacokinetics/ rather than specificDrugs/ because nlmixr2lib's specificDrugs tier is reserved for human drugs.",
    notes          = "Cao 2013 Table 1, Model A. Parameters fit by Cao et al. to plasma concentration profiles from Garg A & Balthasar JP. J Pharmacokinet Pharmacodyn 2007;34:687-709 (PMID 17636457). The function name is mab7E3 because R identifiers cannot start with a digit; the antibody is referred to as 7E3 in the source publications. sigma1 was fixed at 0.95 in Cao 2013 (footnote 'Assumed', CV reported as not applicable in Table 1)."
  )

  ini({
    sigma1 <- 0.95;  label("Vascular reflection coefficient for tight tissues (unitless; fixed at 0.95 in Cao 2013)")  # Cao 2013 Table 1 (Model A): 0.95, fixed (footnote b "Assumed")
    sigma2 <- 0.421; label("Vascular reflection coefficient for leaky tissues (unitless)")                                # Cao 2013 Table 1 (Model A): 0.421 (CV 10.4%)
    lclp   <- log(1.1976e-4); label("Plasma clearance (CLp, L/day)")                                                       # Cao 2013 Table 1 (Model A): CLp = 0.499e-5 L/hr (CV 14.1%) = 1.1976e-4 L/day
  })

  model({
    # Mouse system parameters; Cao 2013 Table 1 footnote and Methods, for a 20 g BW mouse.
    # Vplasma = 0.85 mL = 0.00085 L; ISF = 4.35 mL = 0.00435 L;
    # total lymph flow = 0.12 mL/hr = 0.00288 L/day; sigmaL = 0.2; Kp = 0.8 for native IgG1 (refs 22-23).
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
