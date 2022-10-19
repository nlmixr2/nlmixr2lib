oncology_xenograft_simeoni_2004 <- function() {
  description <- "Oncology tumor growth model in xenograft models"
  reference <- "Monica Simeoni, Paolo Magni, Cristiano Cammia, Giuseppe De Nicolao, Valter Croci, Enrico Pesenti, Massimiliano Germani, Italo Poggesi, Maurizio Rocchetti; Predictive Pharmacokinetic-Pharmacodynamic Modeling of Tumor Growth Kinetics in Xenograft Models after Administration of Anticancer Agents. Cancer Res 1 February 2004; 64 (3): 1094–1101. https://doi.org/10.1158/0008-5472.CAN-03-2524"
  # Values for k1, k2, lambda0, lambda2 are from paclitaxel experiment 1
  # reported in Table 2 from the reference (limits are not from the reference).
  # The values from Table 2 will be estimated on the log scale to ensure
  # positive values.  Residual errors are not in the original reference.
  ini({
    lk1 <- log(c(0.1, 0.968, 10)) ; label("Transit rate through damage (1/day)")
    lk2 <- log(c(0.00001, 0.000629, 0.1)) ; label("Linear drug effect on cycling cells (1/(day*ng/mL))")
    llambda0 <- log(c(0.001, 0.273, 2)) ; label("Tumor exponential growth rate (1/day)")
    llambda1 <- log(c(0.01, 0.814, 5)) ; label("Tumor linear growth rate (tumor volume/day)")
    propErr <- 0.2 ; label("Proportional residual error (fraction)")
    addErr <- 30 ; label("Additive residual error (tumor volume)")
  })
  model({
    k1 <- exp(lk1)
    k2 <- exp(lk2)
    lambda0 <- exp(llambda0)
    lambda1 <- exp(llambda1)

    # tumor0 is provided in the data as the initial volume of the tumor.  It can
    # also be estimated.
    cyclingCells(0) <- tumor0
    psi <- 20 # psi is defined in the paper to cause a rapid switch between exponential and linear regimes
    tumor <- cyclingCells + damagedCells1 + damagedCells2 + damagedCells3
    # cp is provided in the data (or in an appended model) as the drug
    # concentration.  Units for cp will be apply to k2.
    drugEffect <- k2*cp
    d/dt(cyclingCells) <- lambda0*cyclingCells/(1 + (lambda0/lambda1 * tumor)^psi)^(1/psi) - drugEffect*cyclingCells
    d/dt(damagedCells1) <- k2*cp*cyclingCells - k1*damagedCells1
    d/dt(damagedCells2) <- k1*(damagedCells1 - damagedCells2)
    d/dt(damagedCells3) <- k1*(damagedCells2 - damagedCells3)
    tumor ~ prop(propErr) + add(addErr)
  })
}
