Cirincione_2017_exenatide <- function() {
  description <- "Exenatide immediate-release PK model (Cirincione 2017)"
  reference <- "Cirincione B, Mager DE. Population pharmacokinetics of exenatide. British Journal of Clinical Pharmacology. 2017;83(3):517-526. doi:10.1111/bcp.13135"
  covariateData <-
    list(
      AMT = "Dose (ug)",
      DV = "Exenatide plasma concentration (pg/mL)",
      eGFR = "Modification of Diet in Renal Disease estimate of glomerular filtration rate (mL/min/(1.73m^2))",
      WT = "Baseline body weight (kg)",
      DVID = "Was the subject in Study 1 ('study1'), Study 5 ('study5'), or another study ('otherStudy')?  Typically 'otherStudy'"
    )
  # parameters are from Table 2 in the paper
  ini({
    lcl <- log(4.58) ; label("Linear clearance rate (L/hr)")
    etalcl ~ log(1.339)
    e_cl_gfr <- 0.838; label("Effect of eGFR on clearance (unitless)")

    lq <- log(3.72); label("Intercompartmental clearance (L/hr)") # written as Cld in the model table

    lkm <- log(567); label("Michaelis-Menten constant for clearance (pg/mL)")
    etalkm ~ log(1.957)

    lvmax <- log(1.55); label("Maximum Michaelis-Menten clearance (ug/hr)")

    lvp <- log(7.04); label("Peripheral compartment volume (L)")

    lvc <- log(7.03); label("Typical central compartment clearance (L)")
    etalvc ~ log(1.805)
    e_vc_wt <- 2.67; label("Effect of body weight on central volume (unitless)")

    lkamax <- log(0.0813); label("Maximum first-order absorption rate (1/hr)")
    lkmka <- log(16.9); label("Michaelis-Menten constant for absorption (ug)")
    ttau <- fixed(1.35); label("Duration of zero-order absorption")
    fdepot <- fixed(1); label("Bioavailability (fraction)")
    logitfr <- logit(0.628); label("Fraction of dose with first-order absorption")

    expSdOther <- 0.373 ; label("Exponential residual error for all other studies")
    expSdStudy1 <- 0.39 ; label("Exponential residual error for Study 1")
    expSdStudy5 <- 0.08 ; label("Exponential residual error for Study 5")
  })
  model({
    # declare compartment order
    cmt(depot)
    cmt(central)
    cmt(peripheral1)
    # cl equation is from table 2 in the paper
    cl <- exp(lcl + etalcl)*(eGFR/80)^e_cl_gfr
    q <- exp(lq)
    km <- exp(lkm + etalkm)
    vmax <- exp(lvmax)
    vp <- exp(lvp)
    # vc equation is from table 2 in the paper
    vc  <- exp(lvc + etalvc)*(WT/84.8)^e_vc_wt
    kamax <- exp(lkamax)
    kmka <- exp(lkmka)
    fr <- expit(logitfr)

    kel <- cl/vc + vmax/(km*vc + central)
    k12 <- q/vc
    k21 <- q/vp

    # Need to turn k0 off at time > tau
    mtime(tau) <- ttau

    kzero <- (1-fr)*podo(depot)/tau
    if (tad(depot) > tau) kzero <- 0.0

    # Need to turn ka on at time > tau
    ka <- fr*kamax/(kmka + depot)
    if (tad(depot) <= tau) ka <- 0.0

    d/dt(depot) <-  -ka*depot - kzero
    d/dt(central) <- ka*depot + kzero - kel*central - k12*central + k21*peripheral1
    d/dt(peripheral1) <-                              k12*central - k21*peripheral1
    f(depot) <- fdepot

    cp <- central/vc

    cp ~ lnorm(expSdOther) | otherStudy
    cp ~ lnorm(expSdStudy1) | study1
    cp ~ lnorm(expSdStudy5) | study5
  })
}
