indirect_2cpt_inhi_kout_CLV_nd <- function() 
{
    description <- "Two compartment indirect response model with inhibition of kout."
    ini({
        lvc <- 3.45
        label("Central volume of distribution (Vc)")
        lcl <- 1
        label("Clearance (Cl)")
        lvp <- 5
        label("Peripheral volume of distribution (Vp)")
        lq <- 0.1
        label("Intercompartmental clearance (Q)")
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
        vc <- exp(lvc)
        cl <- exp(lcl)
        vp <- exp(lvp)
        q <- exp(lq)
        IC50 <- exp(lIC50)
        kin <- exp(lkin)
        kout <- exp(lkout)
        kel <- cl/vc
        k12 <- q/vc
        k21 <- q/vp
        d/dt(central) <- -kel * central
        d/dt(peripheral1) <- k12 * central - k21 * peripheral1
        d/dt(effect) <- kin - kout * (1 - Cc/(Cc + IC50)) * effect
        Cc <- central/vc
        Cc ~ prop(propSd)
    })
}
