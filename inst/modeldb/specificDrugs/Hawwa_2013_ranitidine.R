Hawwa_2013_ranitidine <- function() {
  description <- paste(
    "One-compartment population PK model for ranitidine in critically ill",
    "children (n = 78, age 15 days to 15.5 years, weight 1.3 to 47 kg)",
    "receiving oral and/or intravenous bolus doses for stress-ulcer or",
    "GORD prophylaxis. First-order absorption with allometric scaling of",
    "clearance (fixed exponent 0.75) and central volume (fixed exponent",
    "1.0) to a 70 kg adult. Cardiac failure or cardiac surgery (pooled",
    "binary indicator) multiplicatively reduces clearance by 53.7%. IIVs",
    "on absorption rate constant and bioavailability were dropped during",
    "model building so the model could minimize; only CL and V carry IIV.",
    "Proportional residual error (Hawwa 2013)."
  )
  reference <- paste(
    "Hawwa AF, Westwood PM, Collier PS, Millership JS, Yakkundi S,",
    "Thurley G, Shields MD, Nunn AJ, Halliday HL, McElnay JC.",
    "Prophylactic ranitidine treatment in critically ill children -",
    "a population pharmacokinetic study.",
    "Br J Clin Pharmacol. 2013;75(5):1265-1276.",
    "doi:10.1111/j.1365-2125.2012.04473.x."
  )
  vignette <- "Hawwa_2013_ranitidine"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric size covariate; reference weight 70 kg per Hawwa 2013",
        "Methods/Results ('Both CL and V were allometrically scaled to an",
        "adult of 70 kg with power values of 0.75 and 1.0 for CL and V,",
        "respectively'). Cohort weight range 1.3 to 47 kg (Table 1)."
      ),
      source_name        = "WT"
    ),
    DIS_HF_OR_CARDSURG = list(
      description        = paste(
        "Pooled binary indicator for cardiac failure or cardiac surgery,",
        "the only non-weight covariate retained after backward stepwise",
        "elimination."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no cardiac failure and no cardiac surgery)",
      notes              = paste(
        "1 = subject has cardiac failure or has undergone (or is in the",
        "postoperative period of) cardiac surgery; 0 = neither. Hawwa 2013",
        "Table 1 reports 14 subjects with cardiac failure, 10 with cardiac",
        "surgery, and 17 in the union (7-subject overlap) out of 78.",
        "Surgical indications include patent ductus arteriosus (n = 2),",
        "atrial or ventricular septal defect (n = 3), hypoplasia of the",
        "left ventricle (n = 1), bilateral pleural effusion (n = 1), and",
        "other (n = 3). Patients who had undergone cardiac surgery",
        "received only intravenous therapy during the first week",
        "following surgery (Table 1 footnote). The covariate is treated",
        "as time-fixed per subject."
      ),
      source_name        = "HEART"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 78L,
    n_studies      = 1L,
    age_range      = "15 days to 15.51 years",
    age_median     = "4.57 years (mean; SD 4.48; Table 1)",
    weight_range   = "1.3 to 47 kg",
    weight_median  = "16.27 kg (mean; SD 12.24; Table 1)",
    sex_female_pct = 100 * 41 / 78,
    race_ethnicity = NULL,
    disease_state  = paste(
      "Critically ill children in paediatric intensive care receiving",
      "ranitidine for prophylaxis against stress ulcers and",
      "gastrointestinal bleeding, or treatment of gastro-oesophageal",
      "reflux disease. 17 of 78 subjects (21.8%) had cardiac failure or",
      "had undergone cardiac surgery; 9 of 78 (11.5%) had renal failure",
      "(not retained as a covariate in the final model)."
    ),
    dose_range     = paste(
      "Oral and intravenous bolus dosing; mean 1.18 mg/kg per dose,",
      "SD 0.43; route 20.5% IV only, 15.4% oral only, 64.1% both",
      "(Table 1)."
    ),
    regions        = paste(
      "United Kingdom (Royal Belfast Hospital for Sick Children and",
      "Alder Hey Children's Hospital, Liverpool)"
    ),
    notes          = paste(
      "Plasma assay (HPLC, Methods 'Ranitidine assay'); LOQ 25 ng/mL.",
      "Concentrations below LOQ but detectable (n = 26, 10.4% of the",
      "final 248 samples) were imputed at LOQ/2 (12.5 ng/mL) per Hing et",
      "al. 2001 [reference 15 in the source]. 248 samples in the final",
      "analysis dataset (median 2 per patient, range 1 to 13).",
      "Original cohort included 91 children; 13 were excluded because",
      "all of their concentrations were below the LOQ. NONMEM VI (level",
      "1.0) with FOCE-I (first-order conditional estimation with",
      "interaction) was used for the population fit."
    )
  )

  ini({
    # Structural parameters: Hawwa 2013 Table 2 FINAL model column
    # (n = 184 plasma concentrations from 78 subjects after pre-dose
    # exclusions; CL and V are allometrically standardised to a 70 kg
    # adult). The RSEs are point-estimate relative standard errors from
    # the NONMEM covariance step, recorded for traceability.
    lka     <- log(1.31);  label("Absorption rate constant ka (1/h)")                           # Table 2 FINAL: ka = 1.31 1/h (RSE 26.1%)
    lcl     <- log(32.1);  label("Total clearance CL standardised to 70 kg (L/h)")              # Table 2 FINAL: CL = 32.1 L/h (RSE 27.4%); Methods page 5 final equation
    lvc     <- log(285);   label("Central volume of distribution V standardised to 70 kg (L)") # Table 2 FINAL: V  = 285 L  (RSE 34.3%); Methods page 5 final equation
    lfdepot <- log(0.275); label("Oral bioavailability F1 (fraction)")                          # Table 2 FINAL: F1 = 0.275 (RSE 27.1%)

    # Allometric exponents: Hawwa 2013 Results paragraph 3 says
    # "Investigation of models where the power values were not fixed,
    # but included as additional thetas (q s), did not result in any
    # significant improvement in model fit." -> the published model
    # keeps 0.75 and 1.0 as fixed structural exponents.
    allo_cl <- fixed(0.75); label("Allometric exponent on CL (unitless, fixed)")  # Methods/Results page 5: fixed at 0.75
    allo_vc <- fixed(1);    label("Allometric exponent on V (unitless, fixed)")   # Methods/Results page 5: fixed at 1.0

    # Covariate effect of cardiac failure / cardiac surgery on CL.
    # Hawwa 2013 Methods page 5 final equation:
    #   TVCL = 32.1 * (WT/70)^0.75 * 0.463^HEART
    # i.e., CL is multiplied by 0.463 when HEART = 1. Encoded on the
    # log scale as e_dis_hf_or_cardsurg_cl = log(0.463); applied as
    # exp(e_dis_hf_or_cardsurg_cl * DIS_HF_OR_CARDSURG).
    e_dis_hf_or_cardsurg_cl <- log(0.463); label("Log-scale effect of cardiac failure or cardiac surgery on CL")  # Table 2 FINAL: theta(HEART,CL) = 0.463 (RSE 23.5%); OFV decrease 12.618, P < 0.001

    # Interindividual variability: exponential (log-normal) on CL and V
    # only. Paper Methods page 4: "the IIVs on both bioavailability and
    # absorption rate constant had to be removed for the model to
    # minimize successfully." Table 2 reports IIV as CV%; the variance
    # is omega^2 = log(1 + CV^2).
    #   IIV CL  = 60.1% -> omega^2 = log(1 + 0.601^2) = 0.30828    # Table 2 FINAL IIV CL  = 60.1% (RSE 33.0%)
    #   IIV V   = 85.0% -> omega^2 = log(1 + 0.850^2) = 0.54382    # Table 2 FINAL IIV V   = 85.0% (RSE 44.5%)
    # No correlation between etalcl and etalvc was reported in the
    # paper; off-diagonal covariance is implicit zero.
    etalcl ~ 0.30828
    etalvc ~ 0.54382

    # Residual error: proportional only. Paper Methods page 4
    # ("Simplification of the residual error model was considered
    # during model building by removing the residual variance
    # component that has a value close to zero") and Results paragraph
    # 1 confirm the additive term was removed and a proportional-only
    # model was retained. Table 2 FINAL: sigma_prop = 59.5% CV ->
    # propSd = 0.595 in nlmixr2's linear-scale parameterisation.
    propSd <- 0.595; label("Proportional residual SD on plasma Cc (fraction)")  # Table 2 FINAL: residual CV 59.5% (RSE 13.4%)
  })

  model({
    # Individual PK parameters. Allometric scaling applies the fixed
    # exponents (0.75 on CL, 1 on V) at each subject's body weight
    # against the 70 kg reference. The pooled cardiac-failure /
    # cardiac-surgery indicator multiplies CL by 0.463 (1 -> 0.463;
    # 0 -> 1).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^allo_cl *
          exp(e_dis_hf_or_cardsurg_cl * DIS_HF_OR_CARDSURG)
    vc <- exp(lvc + etalvc) * (WT / 70)^allo_vc

    # Micro-constant
    kel <- cl / vc

    # ODE system: first-order absorption from depot into central, with
    # first-order elimination from central. Equivalent to NONMEM
    # ADVAN2/TRANS2 (Hawwa 2013 Results, "implemented using NONMEM
    # subroutines ADVAN2 and TRANS2").
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Bioavailability: F1 applies to the depot compartment (oral
    # doses). Intravenous bolus doses given directly into central
    # bypass the depot and are not scaled by F1.
    f(depot) <- exp(lfdepot)

    # Observation and proportional residual error.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
