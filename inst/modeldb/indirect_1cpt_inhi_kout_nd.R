indirect_1cpt_inhi_kout_nd <- function() 
{
    description <- "One compartment indirect response model with inhibition of kout."
    ini({
        lkel <- 0.534
        label("elimination rate (1/d)")
        lIC50 <- 0.67
        label("Drug concentration producing 50% of maximum inhibition at effect site (IC50)")
        lkin <- 0.48
        label("Zero-order rate constant for production of drug response(1/d)")
        lkout <- 0.34
        label("First-order rate constant for loss of drug response")
        propSd <- c(0, 0.5)
        label("Proportional residual error (fraction)")
    })
    model({
        kel <- exp(lkel)
        IC50 <- exp(lIC50)
        kin <- exp(lkin)
        kout <- exp(lkout)
        d/dt(central) <- -kel * central
        d/dt(effect) <- kin - kout * (1 - Cc/(Cc + IC50)) * effect
        Cc <- central/vc
        Cc ~ prop(propSd)
    })
}
