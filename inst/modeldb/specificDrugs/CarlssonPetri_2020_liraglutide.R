CarlssonPetri_2020_liraglutide <- function() {
  reference <- "Carlsson Petri KC, Hale PM, Hesse D, Rathor N, Mastrandrea LD. Liraglutide pharmacokinetics and exposure-response in adolescents with obesity. Pediatric Obesity. 2021;16(10):e12799. doi:10.1111/ijpo.12799"
  covariateData <-
    list(
      WT = "Baseline body weight in kg",
      CHILD = "Is the subject a child? (1 for yes, 0 for no)",
      ADOLESCENT = "Is the subject an adolescent? (1 for yes, 0 for no)",
      SEXM = "1 for male, 0 for female"
    )
  ini({
    lka <- fixed(log(0.0813)) ; label("Absorption rate (1/hr)")
    lcl <- log(1.01) ; label("Apparent clearance (L/h)")
    e_wt_cl <- 0.762; label("Body weight exponent on CL/F")
    e_sex_cl <- 1.12; label("Sex contrast (male/female) on CL/F")
    e_age_child_cl <- 1.11; label("Age contrast (child/adult) on CL/F")
    e_age_adolescent_cl <- 1.06; label("Age contrast (adolescent/adult) on CL/F")
    lvc  <- fixed(log(13.8)) ; label("Apparent central volume of distribution (L)")
    e_wt_vc <- 0.587; label("Body weight exponent on Vc/F")
    
    etalcl ~ log(1.312)
    etalvc ~ log(1.317)
    propSd <- 0.433 ; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    cl_wt <- (WT/100)^e_wt_cl # Equation 2 in the paper
    cl_sex <- SEXM^e_sex_cl # Equation 3 in the paper
    cl_age <- CHILD^e_age_child_cl * ADOLESCENT^e_age_adolescent_cl # Equation 4 in the paper
    cl <- exp(lcl + etalcl)*cl_wt*cl_sex*cl_age # Equation 1 in the paper
    vc_wt <- (WT/100)^e_wt_vc # Not in the paper, based on Equation 2 in the paper
    vc  <- exp(lvc + etalvc)*vc_wt # Equation 5 in the paper
    
    cp <- linCmt()
    cp ~ prop(propSd)
  })
}
