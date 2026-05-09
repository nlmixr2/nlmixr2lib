Cao_2013_tefibazumab <- function() {
  description <- "Second-generation minimal physiologically-based PK (mPBPK) model for tefibazumab in adults (Cao 2013 Model A; clearance from plasma)"
  reference <- "Cao Y, Balthasar JP, Jusko WJ. Second-generation minimal physiologically-based pharmacokinetic model for monoclonal antibodies. J Pharmacokinet Pharmacodyn. 2013 Oct;40(5):597-607. doi:10.1007/s10928-013-9332-2. Tefibazumab plasma data digitized from Hetherington S et al. Antimicrob Agents Chemother. 2006;50(10):3499-3500 (PMID 17005843)."
  vignette <- "Cao_2013_tefibazumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    n_subjects     = NA_integer_,
    n_studies      = 1,
    age_range      = "adults (per Hetherington 2006 source study)",
    weight_range   = "70 kg reference body weight (Cao 2013 Table 2 footnote)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA,
    disease_state  = "End-stage renal disease requiring hemodialysis (tefibazumab humanized IgG1 anti-clumping factor A from Staphylococcus aureus).",
    dose_range     = "10, 20 mg/kg IV (Cao 2013 Figure 5 tefibazumab panel)",
    regions        = NA,
    notes          = "Cao 2013 Table 2, Model A. Parameters fit by Cao et al. to plasma concentration profiles digitized from Hetherington S et al. Antimicrob Agents Chemother 2006;50:3499-3500 (PMID 17005843). Cao 2013 notes that Model A is preferred for tefibazumab on the basis of the latent constraint sigma1 > sigma2, even though Model B has a slightly lower objective function value."
  )

  ini({
    sigma1 <- 0.902; label("Vascular reflection coefficient for tight tissues (unitless)")  # Cao 2013 Table 2 (Model A): 0.902 (CV 7.37%)
    sigma2 <- 0.815; label("Vascular reflection coefficient for leaky tissues (unitless)")  # Cao 2013 Table 2 (Model A): 0.815 (CV 7.56%)
    lclp   <- log(0.22392); label("Plasma clearance (CLp, L/day)")                          # Cao 2013 Table 2 (Model A): CLp = 0.00933 L/hr (CV 4.19%) = 0.22392 L/day
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
