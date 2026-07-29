Sidhu_1998_artemisinin <- function() {
  description <- "One-compartment population PK model with first-order absorption for oral artemisinin in 23 paediatric (2-12 y) and 31 adult (16-45 y) Vietnamese patients with uncomplicated falciparum malaria, fit to sparse capillary plasma samples from a 5-day 10 mg/kg/day field-setting regimen. Separate population estimates for CL/F and V/F are carried for adults (per-subject) and children (per-kg body weight) via a CHILD age-group covariate. Time-dependency in artemisinin disposition is modelled as a 6.9-fold systematic decrease in oral bioavailability between the Day 1 and Day 5 doses, with the published inter-occasion variability on apparent CL/F and V/F retained per occasion. Inter-individual variability on CL/F and V/F is collectively estimated for both age groups; the published etaCL/etaVc correlation 'near unity' is encoded at 0.95 for numerical stability."
  reference <- paste(
    "Sidhu JS, Ashton M, Huong NV, Hai TN, Karlsson MO, Sy ND,",
    "Jonsson EN, Cong LD. (1998).",
    "Artemisinin population pharmacokinetics in children and adults with",
    "uncomplicated falciparum malaria.",
    "Br J Clin Pharmacol 45(4):347-354.",
    "doi:10.1046/j.1365-2125.1998.t01-1-00686.x.",
    sep = " "
  )
  vignette <- "Sidhu_1998_artemisinin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight at study inclusion.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Sidhu 1998 Table 1: paediatric median 20 kg (range 8-32), adult median 46.5 kg (range 34-56). WT scales CL/F and V/F in children only (per-kg parameterisation; see footnote b of Table 2); adult typical values are not WT-scaled.",
      source_name        = "WT"
    ),
    CHILD = list(
      description        = "Age-group indicator (1 = paediatric subject aged 2-12 y, 0 = adult subject aged 16-45 y).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (adult, 16-45 y)",
      notes              = "Sidhu 1998 Results: 'A one-compartment model with separate pharmacokinetic estimates for children and adults was found best to describe the disposition of artemisinin' (Table 2). The categorical age-group covariate selected over total body weight in adults; in children, WT was retained as an additional within-group covariate (footnote b of Table 2). Age cutoffs: paediatric = 2-12 y, adult = 16-45 y.",
      source_name        = "patient group (adults / children)"
    ),
    OCC = list(
      description        = "Integer-valued occasion / period indicator for inter-occasion-variability and the systematic time-dependent change in bioavailability.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Values 1 = Day 1 (first oral dose of the 5-day regimen) and 2 = Day 5 (last oral dose). Decomposed inside model() into binary indicators oc1 and oc2 that (a) multiplex the IOV etas on log CL/F and log V/F per Sidhu 1998 Table 2 (pi_CL/F = 53%, pi_V/F = 86%) and (b) gate the 6.9-fold systematic decrease in oral bioavailability between Day 1 and Day 5 (delta_F_Day1_to_Day5 = 6.9 in Table 2).",
      source_name        = "OCC (Day-1 vs Day-5 sampling occasion)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 54L,
    n_studies      = 1L,
    n_observations = 140L,
    age_range      = "2-45 years",
    weight_range   = "8-56 kg",
    sex_female_pct = 35,
    race_ethnicity = c(Asian = 100),
    disease_state  = "uncomplicated Plasmodium falciparum malaria",
    dose_range     = paste(
      "Oral artemisinin 10 mg/kg/day for 5 days.",
      "Paediatric subjects: 10 mg/kg single dose on Days 1 and 5",
      "(morning), 5 mg/kg twice daily (approx 07:00 and 19:00) on Days 2-4.",
      "Adult subjects: 2 x 250 mg single dose on Days 1 and 5 (morning),",
      "1 x 250 mg twice daily on Days 2-4. Hard gelatine capsules",
      "(25, 50, 100, 150 or 250 mg strengths)."
    ),
    regions        = "Vietnam (Phu Rieng rubber plantation, Song Be Province)",
    notes          = paste(
      "Sidhu 1998 Table 1: 23 paediatric subjects (16 M / 7 F, age 2-12 y,",
      "weight 8-32 kg; 10 in the 2-7 y stratum and 13 in the 8-12 y stratum)",
      "and 31 adult subjects (19 M / 12 F, age 16-45 y, weight 34-56 kg;",
      "19 in the 16-30 y stratum and 12 in the 31-45 y stratum). Total 140",
      "plasma concentrations (107 on Day 1, 33 on Day 5). Sparse-sampling",
      "scheme: 2-3 capillary blood samples per patient on Day 1 (at either",
      "2.5 h and 6 h, or 4 h and 8 h / 10 h post-dose) plus one capillary",
      "sample on Day 5 in 18 children and 15 adults. Children were more",
      "anaemic than adults (paediatric Hb median 102 g/L vs adult 118 g/L,",
      "P = 0.004); haemoglobin was not retained as a PK covariate. Gender",
      "and smoking were screened but excluded from the final model due to",
      "biased estimates and confounding with the adult-male subgroup."
    )
  )

  ini({
    # Structural population parameters - Sidhu 1998 Table 2 ('Estimates of
    # CL/F, V/F and Ka are those for Day 1.'). Adult typical values are the
    # whole-population estimates (not weight-normalised); paediatric typical
    # values are reported per kg total body weight and applied as a
    # CHILD-gated shift on log CL/F and log V/F that combines a fold-change
    # vs the adult anchor with log(WT) (see model() block).

    lcl     <- log(432)                ; label("Adult typical apparent oral clearance CL/F (L/h)")         # Sidhu 1998 Table 2: CL/F_adults = 432 L/h (rel. s.e. 19%)
    lvc     <- log(1600)               ; label("Adult typical apparent volume of distribution V/F (L)")   # Sidhu 1998 Table 2: V/F_adults = 1600 L (rel. s.e. 28%)
    lka     <- log(1.7)                ; label("First-order absorption rate constant ka (1/h)")           # Sidhu 1998 Table 2: ka = 1.7 1/h (rel. s.e. 25%); constrained ka > CL/V during fitting
    lfdepot <- fixed(log(1))           ; label("Oral bioavailability anchor on Day 1 (unitless)")         # Sidhu 1998 Methods: F is identifiable only as a relative quantity; anchored at unity on Day 1 to let the Day-5 covariate carry the time-dependent decrease

    # CHILD covariate effect on log CL/F and log V/F. The pediatric typical
    # value reported as 14.4 L/h/kg (CL/F) / 37.9 L/kg (V/F) is reconstructed
    # inside model() as (lcl_adult + e_child_cl + log(WT)) so the canonical
    # lcl / lvc parameters carry the adult reference and the e_child_*
    # parameters carry the adult->child log shift combined with body weight.
    e_child_cl <- log(14.4 / 432)      ; label("CHILD log-shift on CL/F combined with log(WT) (unitless)")   # Sidhu 1998 Table 2: CL/F_children = 14.4 L/h/kg (rel. s.e. 24%); ratio 14.4/432 = 0.0333 (log = -3.402)
    e_child_vc <- log(37.9 / 1600)     ; label("CHILD log-shift on V/F  combined with log(WT) (unitless)")   # Sidhu 1998 Table 2: V/F_children  = 37.9 L/kg   (rel. s.e. 33%); ratio 37.9/1600 = 0.02369 (log = -3.741)

    # Systematic time-dependent decrease in oral bioavailability. The paper
    # reports a single multiplicative ratio delta_F_Day1_to_Day5 = F_Day1/F_Day5
    # = 6.9 (rel. s.e. 20%); on the log scale, F_Day5 = exp(-log(6.9)) * F_Day1.
    # Gated by the Day-5 occasion indicator oc2 inside model().
    e_occ2_fdepot <- log(1 / 6.9)      ; label("Day-5 log shift on oral bioavailability F (unitless)")    # Sidhu 1998 Table 2: delta_F_Day1->Day5 = 6.9 (rel. s.e. 20%); log(1/6.9) = -1.932

    # Inter-individual variability (IIV). The paper reports omega CL/F = 45%,
    # omega V/F = 104%, omega ka = 576% (Table 2, with relative s.e. 44%,
    # 36%, 21% respectively). Variances on the log scale are translated as
    # omega^2 = log(CV^2 + 1) (standard log-normal interpretation, matching
    # the convention used in the sibling artemisinin extraction
    # Birgersson_2016_artemisinin.R):
    #   CL/F  : 45%  -> omega^2 = log(0.45^2  + 1) = 0.1843
    #   V/F   : 104% -> omega^2 = log(1.04^2  + 1) = 0.7334
    #   ka    : 576% -> omega^2 = log(5.76^2  + 1) = 3.5316
    # The paper's text records 'the correlation between eta_V/F_i and
    # eta_CL/F_i was found to be near unity' (Results, p. 350) but does not
    # publish the numeric off-diagonal. The encoded correlation is set to
    # 0.95 (cov = 0.95 * sqrt(0.1843 * 0.7334) = 0.3493) for numerical
    # stability while preserving the qualitative finding; documented as a
    # deviation in the validation vignette.
    etalcl + etalvc ~ c(0.1843,
                        0.3493, 0.7334)                                                                    # Sidhu 1998 Table 2 + Results: IIV CL/F 45% (rel.s.e. 44%), IIV V/F 104% (rel.s.e. 36%), cor near unity
    etalka          ~ 3.5316                                                                               # Sidhu 1998 Table 2: IIV ka 576% (rel.s.e. 21%); retained per paper's Results choice to estimate rather than fix

    # Inter-occasion variability (IOV) on log CL/F and log V/F. Two
    # occasions are observed in the dataset (Day 1 and Day 5). The paper
    # reports a single IOV estimate per parameter (pi_CL/F = 53%, pi_V/F =
    # 86%; rel. s.e. 32% and 36% respectively), implying the variance is
    # shared across occasions (NONMEM $OMEGA BLOCK(1) SAME idiom). nlmixr2
    # has no SAME shortcut so occasion 1 carries the estimated variance and
    # occasion 2 is fix()-ed to the shared value, mirroring the
    # Aregbe_2012_alvespimycin and Jonsson_2011_ethambutol pattern.
    #   pi_CL/F : 53% -> pi^2 = log(0.53^2 + 1) = 0.2475
    #   pi_V/F  : 86% -> pi^2 = log(0.86^2 + 1) = 0.5536
    etaiov_cl_1 ~ 0.2475                                                                                   # Sidhu 1998 Table 2: pi_CL/F = 53% (rel. s.e. 32%)
    etaiov_cl_2 ~ fix(0.2475)                                                                              # NONMEM $OMEGA BLOCK(1) SAME (occasion 2 shares occasion-1 variance)
    etaiov_vc_1 ~ 0.5536                                                                                   # Sidhu 1998 Table 2: pi_V/F  = 86% (rel. s.e. 36%)
    etaiov_vc_2 ~ fix(0.5536)                                                                              # NONMEM $OMEGA BLOCK(1) SAME (occasion 2 shares occasion-1 variance)

    # Residual variability. The paper reports a 'essentially proportional'
    # residual error model with sigma = 47% (CV; rel. s.e. 41%), expressed
    # in Methods as ln(C_obs) = ln(C_pred) + epsilon, which maps to a
    # proportional residual error on the linear (mg/L = ug/L * 1000) scale
    # in nlmixr2 (see references/verification-checklist.md D and the
    # sibling Birgersson_2016_artemisinin.R for the same convention).
    propSd <- 0.47 ; label("Proportional residual error (fraction)")                                       # Sidhu 1998 Table 2: sigma = 47% (rel. s.e. 41%)
  })

  model({
    # Decompose the integer-valued OCC column into binary occasion
    # indicators. The dataset assembler is expected to set OCC = 1 for all
    # Day-1 dosing and sampling records, and OCC = 2 for all Day-5 dosing
    # and sampling records (the paper's two sampling occasions).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)

    # IOV multiplexers on log CL/F and log V/F. Each subject's two IOV etas
    # are drawn independently; only one is active per occasion via oc1/oc2.
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2
    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2

    # Age-group-dependent typical values for log CL/F and log V/F.
    # CHILD = 0 (adult): ln_cl_typ = lcl (no WT scaling per Sidhu 1998 Results).
    # CHILD = 1 (child): ln_cl_typ = lcl + e_child_cl + log(WT), which equals
    # log(432 * 14.4/432 * WT) = log(14.4 * WT), reproducing the paper's
    # per-kg pediatric parameterisation.
    ln_cl_typ <- lcl + CHILD * (e_child_cl + log(WT))
    ln_vc_typ <- lvc + CHILD * (e_child_vc + log(WT))

    # Individual PK parameters.
    cl  <- exp(ln_cl_typ + etalcl + iov_cl)
    vc  <- exp(ln_vc_typ + etalvc + iov_vc)
    ka  <- exp(lka + etalka)
    kel <- cl / vc

    # One-compartment oral PK with first-order absorption.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Bioavailability: anchored at exp(lfdepot) = 1 on Day 1, and
    # exp(lfdepot + e_occ2_fdepot) = 1/6.9 on Day 5 (the systematic
    # time-dependent decrease in oral bioavailability the paper attributes
    # to autoinduction-like behaviour).
    f(depot) <- exp(lfdepot + e_occ2_fdepot * oc2)

    # Plasma concentration. Dose units mg, vc units L -> central/vc has
    # units of mg/L = 1000 ug/L. The multiplicative 1000 converts to the
    # ug/L units reported in Figure 1 of the paper.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
