Soehoel_2022_tralokinumab <- function() {
  description <- "Tralokinumab PK model (Soehoel 2022)"
  reference <- "Soehoel A, Larsen MS, Timmermann S. Population Pharmacokinetics of Tralokinumab in Adult Subjects With Moderate to Severe Atopic Dermatitis. Clinical Pharmacology in Drug Development. 2022;11(8):910-921. doi:10.1002/cpdd.1113"
  units<-list(time="day",dosing="mg") 
  # From Table 2 footnotes
  covariateData <-
    list(
      nonECZTRA = "1 = any study other than ECZTRA; 0 = ECZTRA study",
      WT = "Body weight in kg",
      dilution = "Was the drug diluted as it was in study D2213C00001? 1 = yes, 0 = no (0 is typical)"
    )
  ini({
    lka <- log(0.184); label("Absorption rate (1/day)")
    lvc <- log(2.71); label("Central volume of distribution (L)")
    lcl <- log(0.149); label("Clearance (L/day)")
    lvp <- log(1.44); label("Peripheral volume of distribution (L)")
    lq <- log(0.159); label("Intercompartmental clearance (L/day)")
    lfdepot <- log(0.761); label("Subcutaneous bioavailability (fraction)")
    CcaddSd <- 0.238; label("Additive residual error (ug/mL)")
    CcpropSd <- 0.216; label("Proportional residual error (fraction)")

    e_wt_vcvp <- 0.793; label("Effect of body weight on central and peripheral volumes (unitless)")
    e_wt_clq <- 0.873; label("Effect of body weight on clearance and intercompartmental clearance (unitless)")
    e_nonECZTRA_cl <- 0.344; label("Effect of non-ECZTRA trials on clearance (unitless)")
    e_nonECZTRA_vc <- 0.258; label("Effect of non-ECZTRA trials on central volume (unitless)")
    e_f_dilution <- 0.354; label("Effect of dilution on bioavailability (unitless)")
    e_ka_dilution <- -0.519; label("Effect of dilution trials on absorption rate (unitless)")

    etavc + etacl ~ c(0.386148, 0.2683494, 0.3057157)
  })
  model({
    fdepot <- exp(lfdepot)*(1 + e_f_dilution*dilution)
    ka <- exp(lka)*(1 + e_ka_dilution*dilution)
    cl <- exp(lcl + etacl)*(WT/75)^e_wt_clq * (1 + e_nonECZTRA_cl*nonECZTRA)
    vc <- exp(lvc + etavc)*(WT/75)^e_wt_vcvp * (1 + e_nonECZTRA_vc*nonECZTRA)
    q <- exp(lq)*(WT/75)^e_wt_clq
    vp <- exp(lvp)*(WT/75)^e_wt_vcvp

    # No unit conversion is required to change mg/L (dosing amount/central
    # volume unit) to ug/mL (measurement unit)
    Cc <- linCmt()
    f(depot) <- fdepot
    Cc ~ add(CcaddSd) + prop(CcpropSd)
  })
}

# etavc, etacl, and the covariance were calculated from the Table 2 footnotes
# as:
# etavc: sqrt(log(0.313^2 + 1)) = 0.386148
# etacl: sqrt(log(0.401^2 + 1)) = 0.3057157
# cov(etavc, etacl): sqrt(0.61*0.386148*0.3057157)
