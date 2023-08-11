indirect_2cpt_stim_kin_nd <- function() 
{
    description <- "Two compartment indirect response model with stimulation of kin."
    ini({
        lvc <- 3.45
        label("Central volume of distribution (Vc)")
        lkel <- 0.534
        label("elimination rate (1/d)")
        lk12 <- 0.48
        label("central-to-peripheral rate (1/d)")
        lk21 <- 0.34
        label("peripheral-to-central rate (1/d)")
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
        vc <- exp(lvc)
        kel <- exp(lkel)
        k12 <- exp(lk12)
        k21 <- exp(lk21)
        EC50 <- exp(lEC50)
        Emax <- exp(lEmax)
        kin <- exp(lkin)
        kout <- exp(lkout)
        d/dt(central) <- -kel * central
        d/dt(peripheral1) <- k12 * central - k21 * peripheral1
        d/dt(effect) <- kin * (1 + Emax * Cc/(Cc + IC50)) - kout * 
            effect
        Cc <- central/vc
        Cc ~ prop(propSd)
    })
}
