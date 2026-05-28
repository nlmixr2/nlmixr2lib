Cao_2013_siltuximab <- function() {
  description <- "Second-generation minimal physiologically-based PK (mPBPK) model for siltuximab in adults (Cao 2013 Model A; clearance from plasma)"
  reference <- "Cao Y, Balthasar JP, Jusko WJ. Second-generation minimal physiologically-based pharmacokinetic model for monoclonal antibodies. J Pharmacokinet Pharmacodyn. 2013 Oct;40(5):597-607. doi:10.1007/s10928-013-9332-2. Siltuximab plasma data digitized from Puchalski T et al. Clin Cancer Res. 2010;16(5):1652-1661 (PMID 20179212)."
  vignette <- "Cao_2013_siltuximab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    n_subjects     = NA_integer_,
    n_studies      = 1,
    age_range      = "adults (per Puchalski 2010 source study)",
    weight_range   = "70 kg reference body weight (Cao 2013 Table 2 footnote)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA,
    disease_state  = "Metastatic renal cell carcinoma (siltuximab anti-IL-6 chimeric IgG1).",
    dose_range     = "1, 3, 6, 12 mg/kg IV (Cao 2013 Figure 5, siltuximab panel)",
    regions        = NA,
    notes          = "Cao 2013 Table 2, Model A. Parameters fit by Cao et al. to plasma concentration profiles digitized from Puchalski T et al. Clin Cancer Res 2010;16:1652-1661 (PMID 20179212). The fit used a typical-value structural mPBPK model in ADAPT 5; no IIV or covariate effects were estimated. Population demographics for the underlying Puchalski 2010 trial are not reproduced in Cao 2013 -- consult the source publication for those details."
  )

  ini({
    # mPBPK Model A parameters from Cao 2013 Table 2 (siltuximab row)
    sigma_tight <- 0.964; label("Vascular reflection coefficient for tight tissues (unitless)")  # Cao 2013 Table 2 (Model A): 0.964 (CV 10.3%)
    sigma_leaky <- 0.673; label("Vascular reflection coefficient for leaky tissues (unitless)")  # Cao 2013 Table 2 (Model A): 0.673 (CV 9.27%)
    lcl   <- log(0.276); label("Plasma clearance (CLp, L/day)")                            # Cao 2013 Table 2 (Model A): CLp = 0.0115 L/hr (CV 5.19%) = 0.276 L/day
  })

  model({
    # Physiological system parameters for a 70 kg adult; Cao 2013 Table 2 footnote and Methods.
    # Vplasma = 2.6 L; ISF = 15.6 L; total lymph flow = 2.9 L/day; sigmaL = 0.2 (Cao 2013 Methods).
    # Kp = 0.8 for native IgG1 (Cao 2013 Methods, refs 22-23). Siltuximab is a chimeric IgG1.
    sigmal     <- 0.2
    kp         <- 0.8
    vplasma    <- 2.6
    visf       <- 15.6
    lymphflow  <- 2.9

    # Lumped tissue volumes and lymph splits (Cao 2013 Eqs 6 and 7)
    vtight  <- 0.65 * visf * kp
    vleaky  <- 0.35 * visf * kp
    l1      <- 0.33 * lymphflow
    l2      <- 0.67 * lymphflow
    vlymph  <- vplasma   # Cao 2013 Methods, ref 21

    cl <- exp(lcl)

    # Concentrations (mg/L); compartments hold amounts (mg)
    cp     <- plasma / vplasma
    ctight <- tight  / vtight
    cleaky <- leaky  / vleaky
    clymph <- lymph  / vlymph

    # Cao 2013 Eqs 1-4 (Model A: clearance from plasma, CLp)
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

    # Plasma concentration (mg/L = ug/mL when dose is in mg and volume in L)
    Cc <- cp
  })
}
