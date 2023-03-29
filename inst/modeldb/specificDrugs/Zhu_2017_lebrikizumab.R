Zhu_2017_lebrikizumab <- function() {
  reference <- "Zhu R, Zheng Y, Dirks NL, et al. Model-based clinical pharmacology profiling and exposure-response relationships of the efficacy and biomarker of lebrikizumab in patients with moderate-to-severe asthma. Pulmonary Pharmacology & Therapeutics. 2017;46:88-98. doi:10.1016/j.pupt.2017.08.010"
  covariateData <-
    list(
      WT = "Baseline body weight in kg",
      AGE = "Age in years",
      SEXF = "1 for female, 0 for male",
      FORM_NS0 = "Is the formulation NS0? 1 for yes, 0 for no (typically no)",
      FORM_CHO_PHASE2 = "Is the formulation CHO from Phase 2? 1 for yes, 0 for no (typically no)",
      ADA = "Is the subject ADA positive ever? 1 for yes, 0 for no",
      RACE_BLACK = "Is the race of the subject black or African American? 1 for yes, 0 for no",
      RACE_ASIAN = "Is the race of the subject Asian? 1 for yes, 0 for no",
      RACE_OTHER = "Is the race of the subject 'other'? 1 for yes, 0 for no"
    )
  ini({
    lcl <- log(0.156); label("Clearance (L/day)")
    lvc <- log(4.10); label("Central volume of distribution (L)")
    lvp <- log(1.45); label("Peripheral volume of distribution (L)")
    lq <- log(0.284); label("Intercompartmental clearance (L/day)")
    lka <- log(0.239); label("Absorption rate (1/day)")
    lfdepot <- log(0.856); label("Subcutaneous bioavailability (fraction)")

    e_cl_wt <- log(1.00); label("Effect of body weight on clearance (unitless)")
    e_vc_wt <- log(0.814); label("Effect of body weight on central volume (unitless)")
    e_vp_wt <- log(0.692); label("Effect of body weight on peripheral volume (unitless)")
    e_q_wt <- log(0.479); label("Effect of body weight on intercompartmentl clearance (unitless)")
    e_cl_age <- log(0.0241); label("Effect of age on clearance (unitless)")
    e_cl_sexf <- log(1.06); label("Effect of sex on clearance (unitless)")
    e_cl_race_black <- log(1.07); label("Effect of race (black or African American) on clearance (unitless)")
    e_cl_race_asian <- log(1.09); label("Effect of race (Asian) on clearance (unitless)")
    e_cl_race_other <- log(1.11); label("Effect of race (other) on clearance (unitless)")
    e_ka_form_nso <- log(0.981); label("Effect of NSO formulation on absorption rate (unitless)")
    e_ka_form_cho_phase2 <- log(0.989); label("Effect of CHO formulation used during Phase 2 on absorption rate (unitless)")
    e_f_form_nso <- log(1.00); label("Effect of NSO formulation on bioavailability (unitless)")
    e_f_form_cho_phase2 <- log(0.973); label("Effect of CHO formulation used during Phase 2 on bioavailability (unitless)")
    e_cl_ada_positive<- log(1.04); label("Effect of anti-drug antibody (ADA) positivity on clearance (unitless)")

    # converted from covariance matrix reported in Table 3
    etacl + etavc + etaka ~
      c(
        0.32403703,
        0.28844410, 0.35213634,
        0.04505552, 0.06625708, 0.39242834
      )

    cppropSd <- 0.0490; label("Proportional residual error (fraction)")
    cpaddSd <- 0.00154; label("Additive residual error (ug/mL)")
  })
  model({
    # Normalized continuous covariate values based on footnote to Table 3
    WTNORM <- WT/70
    AGENORM <- AGE/40

    cl <-
      exp(lcl + etacl) *
      WTNORM^e_cl_wt * AGENORM^e_cl_age * e_cl_sexf^SEXF *
      e_cl_race_black^RACE_BLACK * e_cl_race_asian^RACE_ASIAN * e_cl_race_other^RACE_OTHER *
      e_cl_ada_positive^ADA
    vc <- exp(lvc + etavc) * WTNORM^e_vc_wt
    vp <- exp(lvp) * WTNORM^e_vp_wt
    q <- exp(lq) * WTNORM^e_q_wt
    ka <- exp(lka + etaka) * e_ka_form_nso^FORM_NSO * e_ka_form_cho_phase2^FORM_CHO_PHASE2
    fdepot <- exp(lfdepot) * e_f_form_nso^FORM_NSO * e_f_form_cho_phase2^FORM_CHO_PHASE2
    cp <- linCmt()
    f(depot) <- fdepot
    cp ~ add(cpaddSd) + prop(cppropSd)
  })
}
