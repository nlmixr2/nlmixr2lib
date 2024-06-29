oncology_sdm_lobo_2002 <- function() {
  description <- "Signal transduction model for delayed concentration effects on cancer cell growth"
  reference <- "Lobo ED, Balthasar JP. Pharmacodynamic modeling of chemotherapeutic effects: Application of a transit compartment model to characterize methotrexate effects in vitro. AAPS J. 2002;4(4):212-222. doi:10.1208/ps040442"
  depends<-"Cc"
  units<-list(time="hr")
  # Values for lkng, ltau, lec50, and kmax are for methotrexate from Lobo 2002,
  # Table 2.  propErr and addErr are added as reasonable values though not from
  # Lobo 2002 where no value is apparent in the paper.
  ini({
    lkng <- log(0.02) ; label("Cell net growth rate (growth minus death) (1/hr)")
    ltau <- log(34.1) ; label("Mean transit time of each transit compartment (hr)")
    lec50 <- log(0.1) ; label("Drug concentration reducing the cell growth by 50% (ug/mL)")
    kmax <- 0.29 ; label("Maximum drug-related reduction in cell growth (1/hr)")

    tumorVolpropSd <- c(0, 0.3) ; label("Proportional residual error (fraction)")
    tumorVoladdSd <- c(0, 50, 1000) ; label("Additive residual error (tumor volume units)")
  })
  model({
    # Cc is the drug concentration
    kng <- exp(lkng)
    tau <- exp(ltau)
    ec50 <- exp(lec50)

    drugEffectTumorVol <- kmax*Cc/(ec50 + Cc)

    tumorVol(0) <- tumorVol0
    d/dt(tumorVol) <- kng*tumorVol - transit4*tumorVol
    d/dt(transit1) <- (drugEffectTumorVol - transit1)/tau
    d/dt(transit2) <- (transit1 - transit2)/tau
    d/dt(transit3) <- (transit2 - transit3)/tau
    d/dt(transit4) <- (transit3 - transit4)/tau
    tumorVol ~ prop(tumorVolpropSd) + add(tumorVoladdSd)
  })
}
