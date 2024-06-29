oncology_xenograft_simeoni_2004 <- function() {
  description <- "Oncology tumor growth model in xenograft models"
  reference <- "Monica Simeoni, Paolo Magni, Cristiano Cammia, Giuseppe De Nicolao, Valter Croci, Enrico Pesenti, Massimiliano Germani, Italo Poggesi, Maurizio Rocchetti; Predictive Pharmacokinetic-Pharmacodynamic Modeling of Tumor Growth Kinetics in Xenograft Models after Administration of Anticancer Agents. Cancer Res 1 February 2004; 64 (3): 1094–1101. https://doi.org/10.1158/0008-5472.CAN-03-2524"
  depends<- "Cc"
  units<-list(time="day")
  # Values for damageTransit (k1), drugSlope (k2), tumorExpGrowth (lambda0),
  # tumorLinGrowth (lambda1) are from paclitaxel experiment 1 reported in Table
  # 2 from the reference (limits are not from the reference). The values from
  # Table 2 will be estimated on the log scale to ensure positive values.
  # Residual errors are not in the original reference.
  ini({
    ldamageTransit <- log(c(0.1, 0.968, 10)) ; label("Transit rate through damage (1/day)")
    ldrugSlope <- log(c(0.00001, 0.000629, 0.1)) ; label("Linear drug effect on cycling cells (1/(day*ng/mL))")
    ltumorExpGrowth <- log(c(0.001, 0.273, 2)) ; label("Tumor exponential growth rate (1/day)")
    ltumorLinGrowth <- log(c(0.01, 0.814, 5)) ; label("Tumor linear growth rate (tumor volume/day)")
    tumorVolpropSd <- 0.2 ; label("Proportional residual error (fraction)")
    tumorVoladdSd <- 30 ; label("Additive residual error (tumor volume)")
  })
  model({
    damageTransit <- exp(ldamageTransit)
    drugSlope <- exp(ldrugSlope)
    tumorExpGrowth <- exp(ltumorExpGrowth)
    tumorLinGrowth <- exp(ltumorLinGrowth)

    # tumorVol0 is provided in the data as the initial volume of the tumor.  It
    # can also be estimated.
    cyclingCells(0) <- tumorVol0
    psi <- 20 # psi is defined in the paper to cause a rapid switch between exponential and linear regimes
    tumorVol <- cyclingCells + damagedCells1 + damagedCells2 + damagedCells3
    # Cc is provided in the data (or in an appended model) as the drug
    # concentration.  Units for Cc will be apply to k2.
    drugEffectCyclingCells <- drugSlope*Cc
    d/dt(cyclingCells) <- tumorExpGrowth*cyclingCells/(1 + (tumorExpGrowth/tumorLinGrowth * tumorVol)^psi)^(1/psi) - drugEffectCyclingCells*cyclingCells
    d/dt(damagedCells1) <- drugEffectCyclingCells*cyclingCells - damageTransit*damagedCells1
    d/dt(damagedCells2) <- damageTransit*(damagedCells1 - damagedCells2)
    d/dt(damagedCells3) <- damageTransit*(damagedCells2 - damagedCells3)
    tumorVol ~ prop(tumorVolpropSd) + add(tumorVoladdSd)
  })
}

attr(oncology_xenograft_simeoni_2004, "message") <- "You can modify the number of damaged cell compartments in the model using the function updateOncologyXenograftSimeoni2004(model, ncmt)"
oncology_xenograft_simeoni_2004
