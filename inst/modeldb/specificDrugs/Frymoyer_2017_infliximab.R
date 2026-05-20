Frymoyer_2017_infliximab <- function() {
  description <- "Two-compartment population PK model of intravenous infliximab in children and adults with Crohn's disease (Frymoyer 2017; structural model and parameter values from the Fasanmade et al. REACH + ACCENT I analysis reproduced in Frymoyer 2017 Methods)"
  reference <- paste(
    "Frymoyer A, Hoekman DR, Piester TL, de Meij TG, Hummel TZ, Benninga MA,",
    "Kindermann A, Park KT. Application of population pharmacokinetic modeling",
    "for individualized infliximab dosing strategies in Crohn disease.",
    "J Pediatr Gastroenterol Nutr. 2017;65(6):639-645.",
    "doi:10.1097/MPG.0000000000001620.",
    "Structural model and parameter values were originally developed by",
    "Fasanmade et al. from 112 children in the REACH pediatric Crohn's disease",
    "trial and 580 adults in the ACCENT I adult Crohn's disease trial",
    "(reference 3 of Frymoyer 2017); Frymoyer 2017 reproduces the equations",
    "and final parameter estimates verbatim in its Methods (Model Evaluation)",
    "section and externally validates the model in 34 Dutch children with",
    "Crohn's disease.",
    sep = " "
  )
  vignette <- "Frymoyer_2017_infliximab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL, Vc, and Vp; the source equations are written on per-kilogram CL/V values with reference 65 kg. After conversion to total CL (L/day) and total V (L), the implicit allometric exponents on WT are 1 + 0.313 = 1.313 on CL, 1 + 0.233 = 1.233 on Vc, and 1 + 0.588 = 1.588 on Vp; Q is a constant per-kg value, giving a total-Q exponent of 1.0. The published exponents (0.313, 0.233, 0.588) are preserved in the parameter labels as 'e_wt_<param>'; the +1 conversion is applied in model() and called out in the in-file comments.",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Serum albumin (baseline; treated as static covariate in Frymoyer 2017)",
      units              = "g/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL; normalized as ALB/4.1 per Frymoyer 2017 Methods equation. Reference value 4.1 g/dL. Unit g/dL (US convention) -- distinct from g/L used by other infliximab popPK papers (Hanzel 2021).",
      source_name        = "ALB"
    ),
    ADA_POS = list(
      description        = "Anti-drug antibody positivity (antibodies to infliximab)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative)",
      notes              = "Power-of-coefficient effect on CL: CL_typ * 1.292^ADA_POS (i.e. +29.2% when ADA-positive). Source paper labels this covariate 'ATI' (antibodies to infliximab); renamed to canonical ADA_POS per covariate-columns.md.",
      source_name        = "ATI"
    ),
    CONMED_IMMUNOMOD = list(
      description        = "Concomitant immunomodulator therapy (any of: purine analogue (azathioprine, 6-mercaptopurine) or methotrexate)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant immunomodulator)",
      notes              = "Power-of-coefficient effect on CL: CL_typ * 0.863^CONMED_IMMUNOMOD (i.e. -13.7% when on an immunomodulator). Composite of CONMED_AZA, CONMED_MP, and CONMED_MTX -- the source paper pools all three under the single 'IMM' indicator (Table 1 footnote: 'Concomitant immunomodulation refers to purine-analogue or methotrexate.'). Source paper labels this covariate 'IMM'; renamed to canonical CONMED_IMMUNOMOD per covariate-columns.md.",
      source_name        = "IMM"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 692L,
    n_studies      = 2L,
    age_range      = "Pediatric (REACH, n = 112) and adult (ACCENT I, n = 580); REACH enrolled subjects 6-17 years and ACCENT I enrolled adults 18-65 years.",
    weight_range   = "Not reported in Frymoyer 2017; reference value used by the model is 65 kg.",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported in Frymoyer 2017 for the model-development cohort.",
    disease_state  = "Moderate-to-severe Crohn's disease (pooled REACH pediatric and ACCENT I adult trials).",
    dose_range     = "5-10 mg/kg IV infusion at induction (weeks 0, 2, 6) followed by maintenance every 8 weeks (ACCENT I and REACH design).",
    regions        = "Multi-regional pivotal trials of infliximab in Crohn's disease.",
    notes          = paste(
      "The structural PK model and parameter values were originally developed by",
      "Fasanmade et al. (reference 3 of Frymoyer 2017) by pooling 112 children",
      "from the REACH trial (Moderate-to-Severe Crohn's Disease in Pediatric",
      "Subjects) and 580 adults from the ACCENT I trial. Frymoyer 2017 reproduces",
      "the model's equations and final parameter estimates verbatim in its Methods",
      "(Model Evaluation) section. The validation cohort in Frymoyer 2017 itself",
      "(n = 34 Dutch children with Crohn's disease, median age 14.9 y, median",
      "weight 53 kg, 38% female, median albumin 4.3 g/dL, 12% ADA-positive,",
      "44% on concomitant immunomodulator) is used to assess the predictive",
      "performance of the model, not to develop or update its parameters.",
      "Reference covariate values for the typical subject: WT = 65 kg,",
      "ALB = 4.1 g/dL, ADA-negative (ATI = 0), no concomitant immunomodulator",
      "(IMM = 0).",
      sep = " "
    )
  )

  ini({
    # Structural parameters -- typical total CL (L/day) and total volumes (L)
    # for the reference subject (WT = 65 kg, ALB = 4.1 g/dL, ADA-negative, no
    # concomitant immunomodulator). The source equations are stated on a
    # per-kilogram basis:
    #   CL (mL/kg/day) = 5.42 * (WT/65)^0.313 * (ALB/4.1)^(-0.855)
    #                       * 0.863^IMM * 1.292^ATI
    #   Vc (mL/kg)     = 52.4 * (WT/65)^0.233
    #   Vp (mL/kg)     = 19.6 * (WT/65)^0.588
    #   Q  (mL/kg/day) = 2.26  (constant per-kg)
    # The typical-value parameters here are the per-kg values multiplied by the
    # reference weight 65 kg, then converted to L (or L/day):
    #   lcl = log(5.42 * 65 / 1000) = log(0.3523 L/day)
    #   lvc = log(52.4 * 65 / 1000) = log(3.406 L)
    #   lvp = log(19.6 * 65 / 1000) = log(1.274 L)
    #   lq  = log(2.26 * 65 / 1000) = log(0.1469 L/day)
    lcl <- log(0.3523); label("Clearance for the reference subject (CL, L/day)")          # Frymoyer 2017 Methods (CL equation): 5.42 mL/kg/day * 65 kg / 1000
    lvc <- log(3.406);  label("Central volume of distribution for the reference subject (Vc, L)")     # Frymoyer 2017 Methods (Vc equation): 52.4 mL/kg * 65 kg / 1000
    lvp <- log(1.274);  label("Peripheral volume of distribution for the reference subject (Vp, L)")  # Frymoyer 2017 Methods (Vp equation): 19.6 mL/kg * 65 kg / 1000
    lq  <- log(0.1469); label("Inter-compartmental clearance for the reference subject (Q, L/day)")   # Frymoyer 2017 Methods: 2.26 mL/kg/day * 65 kg / 1000

    # Covariate effect parameters. Power exponents reported by Frymoyer 2017 are
    # on the per-kg form of CL / Vc / Vp; the +1 conversion to total parameters
    # is applied in model().
    e_wt_cl  <-  0.313; label("Power exponent of body weight on per-kg CL ((WT/65)^e_wt_cl on CL/kg)")          # Frymoyer 2017 Methods (CL equation)
    e_wt_vc  <-  0.233; label("Power exponent of body weight on per-kg Vc ((WT/65)^e_wt_vc on Vc/kg)")          # Frymoyer 2017 Methods (Vc equation)
    e_wt_vp  <-  0.588; label("Power exponent of body weight on per-kg Vp ((WT/65)^e_wt_vp on Vp/kg)")          # Frymoyer 2017 Methods (Vp equation)
    e_alb_cl <- -0.855; label("Power exponent of serum albumin on CL ((ALB/4.1)^e_alb_cl)")                     # Frymoyer 2017 Methods (CL equation)
    e_imm_cl <-  0.863; label("Concomitant-immunomodulator power-of-coefficient on CL (e_imm_cl^IMM; -13.7%)")  # Frymoyer 2017 Methods (CL equation)
    e_ada_cl <-  1.292; label("ATI power-of-coefficient on CL (e_ada_cl^ADA_POS; +29.2% when ADA-positive)")    # Frymoyer 2017 Methods (CL equation)

    # Inter-individual variability (diagonal). Frymoyer 2017 Methods reports the
    # IIV as an exponential error model with %CV; the on-log-scale variance is
    # omega^2 = log(1 + CV^2). %CV values: 30.7% (CL), 12.6% (Vc), 55.3% (Vp).
    #   var_cl = log(1 + 0.307^2) = 0.09010
    #   var_vc = log(1 + 0.126^2) = 0.01575
    #   var_vp = log(1 + 0.553^2) = 0.26687
    etalcl ~ 0.09010   # Frymoyer 2017 Methods: 30.7% CV on CL
    etalvc ~ 0.01575   # Frymoyer 2017 Methods: 12.6% CV on Vc
    etalvp ~ 0.26687   # Frymoyer 2017 Methods: 55.3% CV on Vp

    # Residual error -- combined proportional (29.2% CV) + additive
    # (SD 0.371 ug/mL) per Frymoyer 2017 Methods.
    propSd <- 0.292;  label("Proportional residual error (fraction)")   # Frymoyer 2017 Methods
    addSd  <- 0.371;  label("Additive residual error (ug/mL)")          # Frymoyer 2017 Methods
  })
  model({
    # Individual PK parameters. Reference subject: WT = 65 kg, ALB = 4.1 g/dL,
    # ADA-negative, no concomitant immunomodulator. Covariate forms per
    # Frymoyer 2017 Methods (Model Evaluation):
    #   CL (mL/kg/day) = 5.42 * (WT/65)^0.313 * (ALB/4.1)^(-0.855)
    #                       * 0.863^IMM * 1.292^ATI
    #   Vc (mL/kg)     = 52.4 * (WT/65)^0.233
    #   Vp (mL/kg)     = 19.6 * (WT/65)^0.588
    #   Q  (mL/kg/day) = 2.26  (constant)
    # Converting per-kg parameters to totals by multiplying by WT (i.e. carrying
    # the implicit allometric exponent of +1 on each parameter), the model
    # becomes:
    #   CL (L/day) = CL_typ * (WT/65)^(1 + e_wt_cl) * (ALB/4.1)^e_alb_cl
    #                  * e_imm_cl^IMM * e_ada_cl^ADA_POS
    #   Vc (L)     = Vc_typ * (WT/65)^(1 + e_wt_vc)
    #   Vp (L)     = Vp_typ * (WT/65)^(1 + e_wt_vp)
    #   Q  (L/day) = Q_typ  * (WT/65)
    cl <- exp(lcl + etalcl) *
      (WT / 65)^(1 + e_wt_cl) *
      (ALB / 4.1)^e_alb_cl *
      e_imm_cl^CONMED_IMMUNOMOD *
      e_ada_cl^ADA_POS
    vc <- exp(lvc + etalvc) * (WT / 65)^(1 + e_wt_vc)
    vp <- exp(lvp + etalvp) * (WT / 65)^(1 + e_wt_vp)
    q  <- exp(lq)           * (WT / 65)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment model (IV infusion administered to the central compartment).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL
    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
