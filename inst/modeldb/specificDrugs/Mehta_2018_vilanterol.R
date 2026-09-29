Mehta_2018_vilanterol <- function() {
  description <- "Two-compartment population PK model with first-order absorption for inhaled vilanterol in patients with COPD receiving single-inhaler fluticasone furoate/umeclidinium/vilanterol triple therapy, with body weight and age on apparent clearance, refit on the FULFIL study pooled with the historical umeclidinium/vilanterol program"
  reference <- "Mehta R, Pefani E, Beerahee M, Brealey N, Barnacle H, Birk R, Zhu CQ, Lipson DA. Population Pharmacokinetic Analysis of Fluticasone Furoate/Umeclidinium/Vilanterol via a Single Inhaler in Patients with COPD. J Clin Pharmacol. 2018;58(11):1461-1467. doi:10.1002/jcph.1253"
  vignette <- "Mehta_2018_fluticasoneFuroate_umeclidinium_vilanterol"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")
  # Unit note: doses are entered in ug and volumes are in L, so `Cc` is in
  # ug/L == ng/mL. Mehta 2018 reports vilanterol concentrations and exposures
  # in pg/mL and pg*h/mL (assay LLQ 10 pg/mL, upper limit 1000 pg/mL, Sect.
  # 'Pharmacokinetic Assessments'); multiply `Cc` by 1000 to compare against
  # the published values. The source control stream (Supplementary Table 3)
  # carries the same factor as `S2 = V2/1000`; it is NOT reproduced inside the
  # model here so that dose / volume / clearance stay mutually consistent.
  #
  # Structure note: this is the vilanterol model of the UMECLIDINIUM/VILANTEROL
  # (Anoro) program -- two compartments with FIRST-ORDER absorption. It is NOT
  # the three-compartment, zero-order-absorption vilanterol model of the
  # fluticasone furoate/vilanterol program, which is extracted separately as
  # `Siederer_2016_vilanterol`. The Methods narrative of Mehta 2018 describes
  # the latter ('A 3-compartment linear model with zero-order absorption ...
  # with covariates of effect of age (on CL/F and V1/F), body weight (on CL/F),
  # sex and smoking (on V1/F)'), but both printed sources for the model
  # actually fitted here contradict that sentence and agree with each other:
  # Table 2 lists exactly CL/F, V2/F, Q/F, V3/F and KA (a two-compartment
  # first-order-absorption parameterisation, with no zero-order duration and no
  # second peripheral compartment), and Supplementary Table 3 specifies
  # '$SUBS ADVAN4 TRANS4' -- NONMEM's two-compartment first-order-absorption
  # model -- with only body weight and age on CL/F, reading
  # 'Final_Anoro_FulFill_VI_AllDoses.csv'. The equations win over the
  # narrative; see the vignette 'Assumptions and deviations' section.

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power-model effect on CL/F normalised to 70 kg, per Supplementary Table 3 '$PK MU_1=LOG(THETA(1)) + WTEX1*LOG(WT/70) + AGEX*LOG(AGE/60)'. Mean weight in the FULFIL PK population was 81 kg (Table 3). The exponent is not re-estimated or reported for the combined dataset and is held at a structural zero -- see the ini() comments.",
      source_name = "WT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power-model effect on CL/F normalised to 60 years, per Supplementary Table 3 '$PK ... + AGEX*LOG(AGE/60)'. Mean age in the FULFIL PK population was 64 years (Table 3); FULFIL required patients to be at least 40 years old. The exponent is not re-estimated or reported for the combined dataset and is held at a structural zero.",
      source_name = "AGE"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "vilanterol", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "vilanterol", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vilanterol", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 74,
    n_studies = 1,
    age_mean = "64 years",
    weight_mean = "81 kg",
    sex_female_pct = 26,
    race_ethnicity = c(`White` = 100),
    disease_state = "symptomatic chronic obstructive pulmonary disease (mean 45% predicted FEV1)",
    dose_range = "vilanterol 25 ug once daily by oral inhalation, as the fluticasone furoate/umeclidinium/vilanterol 100/62.5/25 ug single-inhaler triple combination (Ellipta)",
    regions = "global (162 centers in 15 countries: Russian Federation, Ukraine, Mexico, Germany, Greece, Czech Republic, Romania, Bulgaria, China, Estonia, Hungary, Italy, Poland, Republic of Korea, Slovakia)",
    notes = "Demographics are the FULFIL (CTT116853, NCT02345161) PK population of 74 patients randomized to fluticasone furoate/umeclidinium/vilanterol who provided serial (n = 10) or sparse (n = 64) samples at weeks 12 and 24 (Table 3); mean BMI 28 kg/m2, mean height 171 cm. The PARAMETER ESTIMATES in this file were obtained on a COMBINED dataset that pools these FULFIL data with the historical vilanterol data of the umeclidinium/vilanterol program (source file 'Final_Anoro_FulFill_VI_AllDoses.csv') used to build the Goyal 2014 model; the size of the historical half is not restated in Mehta 2018, so n_subjects records the FULFIL contribution only. Data below the 10 pg/mL quantification limit were treated as censored using the NONMEM M3 full-likelihood approach (Ahn 2008)."
  )

  ini({
    # -----------------------------------------------------------------------
    # Structural parameters. Mehta 2018 Table 2 reports the vilanterol THETAs
    # only on the natural scale, in its 'Model Parameter Estimates With
    # Combined Dataset (RSE%)' column; the 'Historical Model Parameter
    # Estimates' column of the same table is NOT used here. Supplementary
    # Table 3 parameterises them as MU_1 = LOG(THETA(1)) + <covariates>, i.e.
    # the THETAs are natural-scale, so each is wrapped in log() below.
    #
    # The deposited $THETA block holds round-number INITIAL estimates
    # ((0, 45) for CL, (0, 200) for V2, (0, 160) for Q, (0, 100) for V3,
    # (0, 5) for KA) that match neither column of Table 2 -- note in particular
    # that the V3 initial of 100 L is more than an order of magnitude below the
    # 1280 L estimate, so these are search starting points, not results.
    # -----------------------------------------------------------------------
    lka <- log(19.6); label("Apparent first-order absorption rate constant after oral inhalation (1/h)") # Table 2 vilanterol 'KA (h-1)' combined-dataset estimate 19.6 (RSE 9.5%)
    lcl <- log(41.6); label("Apparent inhaled clearance CL/F at the covariate reference (L/h)") # Table 2 vilanterol 'CL/F (L/h)' combined-dataset estimate 41.6 (RSE 1.5%)
    lvc <- log(271); label("Apparent central volume of distribution V2/F (L)") # Table 2 vilanterol 'V2/F (L)' combined-dataset estimate 271 (RSE 2.1%)
    lq <- log(116); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2 vilanterol 'Q/F (L/h)' combined-dataset estimate 116 (RSE 4.4%)
    lvp <- log(1280); label("Apparent peripheral volume of distribution V3/F (L)") # Table 2 vilanterol 'V3 /F (L)' combined-dataset estimate 1280 (RSE 4.6%)

    # -----------------------------------------------------------------------
    # Covariate effects. Supplementary Table 3 gives the exact functional form
    #   MU_1 = LOG(THETA(1)) + WTEX1*LOG(WT/70) + AGEX*LOG(AGE/60)
    # i.e. plain power models on WT/70 and AGE/60, applied to CL/F only
    # (MU_2 through MU_5 carry no covariate term). The $PROBLEM line agrees:
    # 'Wt and AGE on CL'.
    #
    # The EXPONENTS are not reported for the combined dataset: Table 2 lists
    # only the five structural rows, and Methods states 'a covariate analysis
    # was not planned. The same covariate relationship was assumed'. The
    # stream's $THETA entries for them (0.5, -0.5) are round-number initial
    # estimates from the same block shown above to hold initials rather than
    # results. They are therefore held at a structural zero rather than
    # transcribed from an initial estimate; the functional form is preserved so
    # a downstream user can supply exponents from the historical analysis.
    # See the vignette 'Assumptions and deviations' section.
    # -----------------------------------------------------------------------
    e_wt_cl <- fixed(0); label("Power exponent on (WT/70) for CL/F (unitless)") # Supplementary Table 3 '$PK WTEX1 = THETA(8)' in 'MU_1=LOG(THETA(1)) + WTEX1*LOG(WT/70)'; no combined-dataset estimate reported in Table 2
    e_age_cl <- fixed(0); label("Power exponent on (AGE/60) for CL/F (unitless)") # Supplementary Table 3 '$PK AGEX = THETA(9)' in 'MU_1=... + AGEX*LOG(AGE/60)'; no combined-dataset estimate reported in Table 2

    # -----------------------------------------------------------------------
    # Inter-individual variability. Supplementary Table 3 places an exponential
    # ETA on all five structural parameters, and Results confirms individual
    # MAP Bayes estimates were obtained for all 74 patients. The magnitudes are
    # not reported for the combined dataset (Table 2 has no OMEGA rows) and the
    # stream's $OMEGA block (0.5, 0.7, 0.3, 0.4, 0.3) is round-number initials,
    # so each is declared at fixed(0).
    #
    # The stream also carries INTER-OCCASION variability -- three $OMEGA
    # BLOCK(1) SAME pairs on CL and on V2 across OCC = 1, 2, 3, plus an eta on
    # the proportional residual term (ETA(12)). Per nlmixr2lib convention IOV
    # is not encoded in library models; it is recorded here and in the vignette
    # so the omission is auditable.
    # -----------------------------------------------------------------------
    etalka ~ fixed(0) # Supplementary Table 3 'KA = EXP(MU_5+ETA(5))'; combined-dataset variance not reported
    etalcl ~ fixed(0) # Supplementary Table 3 'CLB = EXP(MU_1+ETA(1))'; combined-dataset variance not reported
    etalvc ~ fixed(0) # Supplementary Table 3 'V2B = EXP(MU_2+ETA(2))'; combined-dataset variance not reported
    etalq ~ fixed(0) # Supplementary Table 3 'Q = EXP(MU_3+ETA(3))'; combined-dataset variance not reported
    etalvp ~ fixed(0) # Supplementary Table 3 'V3 = EXP(MU_4+ETA(4))'; combined-dataset variance not reported

    # -----------------------------------------------------------------------
    # Residual error. Supplementary Table 3 builds a COMBINED additive plus
    # proportional standard deviation on the untransformed concentration scale,
    # identical in form to the umeclidinium stream:
    #   SIG  = THETA(6)                      (additive SD, pg/mL)
    #   SIG2 = F*THETA(7)*EXP(ETA(12))       (proportional SD)
    #   SD   = SQRT(SIG*SIG + SIG2*SIG2)
    #   Y    = F + ERR(1)*SD                 with $SIGMA 1 FIXED
    # Neither magnitude is reported for the combined dataset and the stream's
    # initials (15 and 0.2) are round numbers, so both are declared at fixed(0).
    # The eta on the proportional term is not representable in this file and is
    # documented in the vignette instead.
    # -----------------------------------------------------------------------
    addSd <- fixed(0); label("Additive residual error (ng/mL)") # Supplementary Table 3 '$ERROR SIG=THETA(6)'; combined-dataset magnitude not reported
    propSd <- fixed(0); label("Proportional residual error (fraction)") # Supplementary Table 3 '$ERROR SIG2=F*THETA(7)*EXP(ETA(12))'; combined-dataset magnitude not reported
  })

  model({
    # Power covariate models on CL/F only, exactly as written in Supplementary
    # Table 3's $PK block: the log-additive MU_ formulation is algebraically
    # identical to the multiplicative power form used here.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (AGE / 60)^e_age_cl
    vc <- exp(lvc + etalvc)
    ka <- exp(lka + etalka)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
