indirect_2cpt_inhi_kin_r0rmaxcrmax_nd <- function() 
{
    description <- "Two compartment indirect response model with inhibition of kin."
    ini({
        lvc <- 3.45
        label("Central volume of distribution (Vc)")
        lcl <- 0.85
        label("Clearance (Cl)")
        lvp <- 5
        label("Peripheral volume of distribution (Vp)")
        lq <- 0.1
        label("Intercompartmental clearance (Q)")
        lr0 <- 0.2
        label("Baseline response prior to drug administration (R0)")
        lrmax <- 0.9
        label("Maximal response (CRmax)")
        ls1 <- 1
        label("Initial slope of the response versus time curve (S1)")
        limax <- 0.56
        label("Maximum inhibitory factor attributed to drug (Imax)")
        lcrmax <- 0.67
        label("Plasma concentration of drug at the time of maximal response (CRmax)")
        propSd <- c(0, 0.5)
        label("Proportional residual error (fraction)")
    })
    model({
        vc <- exp(lvc)
        cl <- exp(lcl)
        vp <- exp(lvp)
        q <- exp(lq)
        r0 <- exp(lr0)
        rmax <- exp(lrmax)
        s1 <- exp(ls1)
        imax <- exp(limax)
        crmax <- exp(lcrmax)
        kel <- cl/vc
        k12 <- q/vc
        k21 <- q/vp
        imax <- (r0 - rmax)/r0
        kin <- -s1/imax
        kout <- kin/r0
        IC50 <- crmax * (rmax - (1 - imax) * r0)/(r0 - rmax)
        d/dt(central) <- -kel * central
        d/dt(peripheral1) <- k12 * central - k21 * peripheral1
        d/dt(effect) <- kin * (1 - Cc/(Cc + IC50)) - kout * effect
        Cc <- central/vc
        Cc ~ prop(propSd)
    })
}
