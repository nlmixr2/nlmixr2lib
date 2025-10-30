igg_kim_2006 <- function() {
  description <- "Immunoglobulin G (IgG) model for nonlinear metabolism in healthy subjects"
  reference <- "Kim J, Hayton WL, Robinson JM, Anderson CL. Kinetics of FcRn-mediated recycling of IgG and albumin in human: pathophysiology and therapeutic implications using a simplified mechanism-based model. Clin Immunol. 2007 Feb;122(2):146-55. doi: 10.1016/j.clim.2006.09.001. Epub 2006 Oct 13. PMID: 17046328; PMCID: PMC2791364."
  # parameters from table 1 in paper
  ini({
    ljmax <- log(147); label("Maximal rate of FcRn-mediated recycling (Jmax); i.e., in vivo recycling capacity (mg/d/kg")
    lkm <- log(21); label("Plasma concentration at which a half Jmax reached (Km); The Michaelis constant (mg/mL)")
    kint <- 0.18; label("Fractional intrinsic catabolic rate (1/day)")
    v1 <- 42; label("Vascular compartment volume (mg/kg)")
    css <- 12.1; label("Steady-state plasma concentration (mg/mL)")
  })
  model({
    jmax <- exp(ljmax)
    km <- exp(lkm)
    # Equation 2
    krmr_0 <- jmax/(v1*(km + css))
    jrmr <- krmr_0 * v1 * css
    # Equation 3
    kcat_0 <- kint - krmr_0
    # Equation 1
    jcat_0 <- kcat_0 * v1 * css
    jpro_0 <- jcat_0

    krmr <- jmax/(v1*(km + igg))
    kcat <- kint - krmr
    d/dt(igg) <- jpro_0/v1 - kcat*igg
    igg_0 <- css
    igg(0) <- css
  })
}
