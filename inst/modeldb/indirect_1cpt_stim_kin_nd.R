indirect_1cpt_stim_kin_nd <- function() 
{
    description <- "One compartment indirect response model with stimulation of kin.Parameterized using rate cosntants"
    ini({
        lkel <- 0.534
        label("elimination rate (1/d)")
        lEC50 <- 0.67
        label("Drug concentration producing 50% of maximum stimulation at effect site (EC50)")
        lEmax <- 0.85
        label("Maximum effect attributed to drug (Emax)")
        lkin <- 0.48
        label("Zero-order rate constant for production of drug response(1/d)")
        lkout <- 0.34
        label("First-order rate constant for loss of drug response")
        propSd <- c(0, 0.5)
        label("Proportional residual error (fraction)")
    })
    model({
        kel <- exp(lkel)
        EC50 <- exp(lEC50)
        Emax <- exp(lEmax)
        kin <- exp(lkin)
        kout <- exp(lkout)
        d/dt(central) <- -kel * central
        d/dt(effect) <- kin * (1 + Emax * Cc/(Cc + IC50)) - kout * 
            effect
        Cc <- central/vc
        Cc ~ prop(propSd)
    })
}
