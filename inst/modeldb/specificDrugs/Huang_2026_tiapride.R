Huang_2026_tiapride <- function() {
  description <- paste(
    "One-compartment oral population PK model for tiapride in Chinese children and adolescents",
    "aged 5-15 years treated for tic disorders (Huang 2026), fitted to 214 opportunistic steady-state",
    "plasma samples from 38 outpatients. First-order absorption with linear elimination; the",
    "absorption rate constant (Ka = 0.219 1/h) is far smaller than the elimination rate constant",
    "(CL/F / Vd/F = 15.3 / 5.77 = 2.65 1/h), so disposition is flip-flop and the apparent terminal",
    "half-life is set by absorption (ln(2)/Ka = 3.17 h, matching the 3.23 h literature value the",
    "paper cites). Fat-free mass is the only retained covariate, power-scaled on apparent clearance",
    "with exponent 0.553 referenced to the cohort median 30.62 kg; it displaced creatinine clearance,",
    "with which it is strongly correlated. Exponential interindividual variability is carried on all",
    "three structural parameters, and residual variability is combined proportional plus additive.",
    "Monte Carlo simulation from this model supports 75 mg three times daily as the regimen keeping",
    "steady-state peak concentration inside the 560-2000 ng/mL therapeutic window.",
    "The paper's companion plasma-saliva joint model is included: a saliva compartment is driven",
    "from central by saturable Michaelis-Menten transport (Vmax = 34.7 mg/h, Km = 762 ng/mL fixed)",
    "and cleared first-order (K30 = 6.24 1/h, fat-free mass exponent 0.38), giving a salivary",
    "concentration that exceeds plasma and whose saliva:plasma ratio falls as concentration rises.",
    "That layer was fitted sequentially with the plasma parameters held fixed, and the saliva state",
    "is driven rather than mass-balance-coupled, so it does not deplete central. The apparent saliva",
    "volume converting the saliva amount to the reported ng/mL is NOT given anywhere in the paper;",
    "it is recovered from the Figure 4 visual predictive check and is the one assumed value in this",
    "file -- see the vignette Errata before using the saliva output quantitatively.",
    sep = " "
  )
  reference <- paste(
    "Huang W, Shen J, Luo X, Wu Y, Zheng Y, Zhou J, Xu B, Yin X, Wu X.",
    "Population Pharmacokinetics of Tiapride in Children and Adolescents with Tic Disorders:",
    "Leveraging Plasma and Saliva Concentration to Guide Individualized Dosing.",
    "Drug Des Devel Ther. 2026;20. doi:10.2147/DDDT.S587387",
    sep = " "
  )
  vignette <- "Huang_2026_tiapride"

  # The paper reports every concentration in ng/mL. This model works in the
  # library-standard mg / L / h set, so Cc is in mg/L and 1 mg/L = 1000 ng/mL
  # (the 560-2000 ng/mL therapeutic window is 0.56-2.0 mg/L).
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(
      analyte = "tiapride",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "tiapride",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    # Secreted saliva, the matrix a Salivette collects (Methods 'Sample
    # collection'). Figure 2 draws the central-to-saliva transport arrow
    # DASHED where every other arrow in the diagram is solid, and the joint
    # model was fitted with the plasma parameters held fixed, so this state is
    # driven by -- not drawn out of -- central. See the model() note.
    saliva = list(
      analyte = "tiapride",
      units = "mg",
      specimen = "saliva",
      verified = TRUE
    )
  )

  covariateData <- list(
    FFM = list(
      description = "Fat-free mass at baseline.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The only covariate retained in the final model, power-scaled on apparent clearance and",
        "referenced to the cohort median 30.62 kg (Huang 2026 Table 1, IQR 26.93-34.65 kg;",
        "Equation 6). Age, weight, height, BSA, FFM and CLCR were all screened on CL/F by stepwise",
        "regression and only FFM survived (dOFV = -25.1). The Discussion is explicit that FFM won",
        "BECAUSE it is strongly correlated with creatinine clearance -- tiapride is predominantly",
        "renally eliminated, so FFM is standing in for renal function here rather than for a",
        "distribution volume. The paper does NOT state which fat-free-mass equation was used to",
        "derive the column; for a 5-15 year old cohort the Al-Sallami et al. paediatric correction",
        "to the Janmahasatian adult form is the usual choice, and the vignette assumes it. Note the",
        "exponent 0.553 has a wide bootstrap 95% CI (0.277-0.771) that excludes neither the",
        "theory-based allometric 0.75 nor a linear-per-kg 1, so it is not sharply identified."
      ),
      source_name = "FFM"
    )
  )

  # Screened in the covariate analysis and NOT retained in the final model
  # (Huang 2026 Results "PopPK Model Based on Plasma Concentration" and
  # Discussion). Documented here for provenance; deliberately absent from
  # model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      notes = "Median 8 years (IQR 7-10), Huang 2026 Table 1. Screened on CL/F, not retained."
    ),
    WT = list(
      description = "Total body weight.",
      units = "kg",
      type = "continuous",
      notes = "Median 36.7 kg (IQR 31.0-42.5), Huang 2026 Table 1. Screened on CL/F, not retained; FFM won."
    ),
    HT = list(
      description = "Height.",
      units = "cm",
      type = "continuous",
      notes = "Median 135 cm (IQR 130-146.5), Huang 2026 Table 1. Screened on CL/F, not retained."
    ),
    BSA = list(
      description = "Body surface area.",
      units = "m^2",
      type = "continuous",
      notes = paste(
        "Median 1.18 (IQR 1.06-1.32), Huang 2026 Table 1. Screened on CL/F, not retained.",
        "Table 1 prints the unit as 'cm2', which cannot be right for values near 1.18 in an",
        "8-year-old; the values are m^2 and the printed unit is a typographical error."
      )
    ),
    CLCR = list(
      description = "Creatinine clearance.",
      units = "mL/min",
      type = "continuous",
      notes = paste(
        "Median 118.3 mL/min (IQR 106.77-136.79), Huang 2026 Table 1. Screened on CL/F and NOT",
        "retained even though tiapride is predominantly renally eliminated: the Discussion states",
        "FFM gave the greater dOFV reduction 'owing to its strong correlation with CLCR'. The",
        "estimating equation is not stated in the paper."
      )
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units = "(binary)",
      type = "categorical",
      notes = "31 male / 7 female, Huang 2026 Table 1. Screened and not significant (Discussion)."
    ),
    CONMED_ANY = list(
      description = "Any concomitant medication.",
      units = "(binary)",
      type = "categorical",
      notes = paste(
        "23 of 38 patients had a combined-medication case, Huang 2026 Table 1. The Discussion names",
        "aripiprazole, topiramate, clonidine, sodium valproate and traditional Chinese medicine and",
        "reports that none exhibited a statistically significant effect on tiapride PK. Recorded as",
        "a single screened-and-rejected any-comedication flag because the paper reports no",
        "per-drug effect estimates."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 38,
    n_studies = 1,
    age_range = "5-15 years; median 8 (IQR 7-10)",
    weight_median = "36.7 kg (IQR 31.0-42.5)",
    height_median = "135 cm (IQR 130-146.5)",
    ffm_median = "30.62 kg (IQR 26.93-34.65)",
    bmi_median = "19.1 (IQR 16.72-21.89)",
    bsa_median = "1.18 m^2 (IQR 1.06-1.32)",
    sex_female_pct = 18.4,
    race_ethnicity = c(Asian = 100),
    disease_state = "tic disorders diagnosed by DSM-5, without organic disease or other neuropsychiatric comorbidity",
    renal_function = "creatinine clearance median 118.3 mL/min (IQR 106.77-136.79); serum creatinine median 47 umol/L (IQR 41-52)",
    dose_range = "oral tiapride 2-10 mg/kg/day given two or three times daily; median total daily dose 215 mg/day (IQR 150-300)",
    regions = "China (single centre, Fujian)",
    notes = paste(
      "Single-centre prospective observational outpatient study at Fujian Medical University Union",
      "Hospital, April 2024 to October 2025, with 6 months of follow-up per patient. Paired plasma",
      "and saliva samples were taken before and after the final dose after at least 7 days of",
      "continuous treatment, so all data are at steady state; the post-dose sampling interval had a",
      "median of 4.88 h (IQR 2.17-13.81). Sampling was opportunistic and tied to clinic visits:",
      "45 samples (21%) fell in the absorption phase (0-2 h), 18 (8%) around Tmax (2-2.5 h) and 101",
      "(49%) in the late elimination period (> 10 h), leaving the 2.5-10 h window sparse -- which is",
      "why a two-compartment model was unstable and a one-compartment model was selected despite",
      "the biphasic disposition reported for tiapride in adults. Of 215 plasma samples collected,",
      "one was below the 2 ng/mL LLOQ and was discarded (Beal M1), leaving 214 in the analysis.",
      "Saliva was collected with Salivette cotton-swab devices; 205 saliva samples were taken",
      "(fewer than plasma because of insufficient volume and contamination) and one below the LLOQ",
      "was likewise dropped by M1, leaving 204. The LC-MS/MS calibrated range was 2-1000 ng/mL in",
      "plasma and 4-2000 ng/mL in saliva, so a substantial part of the observed saliva data in",
      "Figure 4B sits above the highest saliva calibrator.",
      "Tiapride is supplied as 100 mg tablets divisible into halves, thirds and quarters, so",
      "clinical doses are rounded to 50, 66.6 or 75 mg per administration."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters. Huang 2026 Table 2 "Final Model / Estimate"
    # column, written out as Equations 5-7. All are APPARENT (CL/F, Vd/F):
    # the study is oral-only, so bioavailability is not identifiable and no
    # F term is reported.
    # ------------------------------------------------------------------
    lka <- log(0.219); label("First-order absorption rate constant, Ka (1/h)")   # Table 2 theta_Ka = 0.219 (RSE 3%; bootstrap median 0.22, 95% CI 0.209-0.234); Equation 5
    lcl <- log(15.3); label("Apparent clearance, CL/F (L/h)")                    # Table 2 theta_CL/F = 15.3 (RSE 3%; bootstrap median 15.4, 95% CI 14.5-16.3); Equation 6
    lvc <- log(5.77); label("Apparent central volume of distribution, Vd/F (L)") # Table 2 theta_Vd/F = 5.77 (RSE 16%; bootstrap median 7.25, 95% CI 1.99-11.6); Equation 7

    # ------------------------------------------------------------------
    # Covariate effect. Equation 6 prints
    #   CL/F = 15.3 * (FFM/30.62)^0.553 * exp(eta_CL/F)
    # with 30.62 kg the Table 1 cohort median fat-free mass.
    # ------------------------------------------------------------------
    e_ffm_cl <- 0.553; label("Power exponent for fat-free mass on CL/F (unitless)") # Table 2 theta_FFM-CL/F = 0.553 (RSE 20%; bootstrap median 0.553, 95% CI 0.277-0.771); dOFV = -25.1

    # ------------------------------------------------------------------
    # Interindividual variability. Exponential, P_i = theta * exp(eta_i)
    # with eta ~ N(0, omega^2) (Methods, Equations 1-3).
    #
    # VARIANCE CONVENTION. Table 2 heads these rows 'eta_Ka (%)',
    # 'eta_CL/F (%)' and 'eta_Vd/F (%)' -- the eta itself as a percentage,
    # i.e. omega * 100, so variance = (percent/100)^2. The paper prints no
    # CV% column and no control stream that would settle it independently.
    # Read instead as a log-normal CV% (omega = sqrt(log(1 + CV^2))) the two
    # small terms barely move (17.6% -> 0.1745, 22.8% -> 0.2249) but Vd/F
    # would drop from 0.843 to 0.7444; the choice is recorded in the
    # vignette Errata.
    # ------------------------------------------------------------------
    etalka ~ 0.030976   # Table 2 eta_Ka = 17.6% -> 0.176^2 (RSE 13%, shrinkage 31%; bootstrap 16.5%, 95% CI 9.9-21.8)
    etalcl ~ 0.051984   # Table 2 eta_CL/F = 22.8% -> 0.228^2 (RSE 15%, shrinkage 22%; bootstrap 21.4%, 95% CI 12.6-29.7)
    etalvc ~ 0.710649   # Table 2 eta_Vd/F = 84.3% -> 0.843^2 (RSE 26%, shrinkage 56%; bootstrap 72.8%, 95% CI 30.2-133.4)

    # ------------------------------------------------------------------
    # Residual unexplained variability. Combined proportional plus additive
    # (Results: 'incorporating a combined error model for residual
    # variability'). The additive term is printed in ng/mL and is converted
    # to the mg/L working unit of this model.
    # ------------------------------------------------------------------
    propSd <- 0.156; label("Proportional residual SD for plasma Cc (fraction)") # Table 2 epsilon_prop = 15.6% (RSE 25%, shrinkage 43%; bootstrap 15.2%, 95% CI 9.1-21.3)
    addSd <- 0.0879; label("Additive residual SD for plasma Cc (mg/L)")         # Table 2 epsilon_add = 87.9 ng/mL = 0.0879 mg/L (RSE 57%, shrinkage 43%; bootstrap 81.4, 95% CI 15.3-294.9)

    # ==================================================================
    # SALIVA LAYER. Huang 2026 Table 3, written out as Equations 8-10.
    #
    # Fitted SEQUENTIALLY (Methods 'PopPK modeling'): the plasma model was
    # established first, its parameters were then held fixed, and the saliva
    # compartment was added with the NONMEM ADVAN6 subroutine. Every value
    # below therefore belongs to the second estimation step; nothing above
    # this divider moves when the saliva layer is present.
    # ==================================================================
    lvmax_saliva <- log(34.7); label("Maximum rate of saturable central-to-saliva transport, Vmax (mg/h)")      # Table 3 theta_Vmax = 34.7 (RSE 13%; bootstrap median 35, 95% CI 27.3-46.9); Equation 8
    lkm_saliva <- fixed(log(0.762)); label("Michaelis constant of central-to-saliva transport, Km (mg/L)")      # Table 3 theta_Km = 762 ng/mL = 0.762 mg/L, held constant because its RSE was unacceptably high (Results 'PopPK Model Based on Saliva Concentration'); Equation 9
    lkel_saliva <- log(6.24); label("First-order elimination rate constant from saliva, K30 (1/h)")             # Table 3 theta_K30 = 6.24 (RSE 4%; bootstrap median 6.22, 95% CI 5.72-6.8); Equation 10
    e_ffm_kel_saliva <- 0.38; label("Power exponent for fat-free mass on K30 (unitless)")                       # Table 3 theta_FFM-K30 = 0.38 (RSE 37%; bootstrap median 0.367, 95% CI 0.158-0.627)

    # ------------------------------------------------------------------
    # APPARENT SALIVA VOLUME -- THE ONE ASSUMED VALUE IN THIS FILE.
    #
    # Table 3 lists exactly six parameters and no volume, and Equations 8-10
    # (all rendered from the PDF's vector graphics and read) introduce none,
    # so the scale converting the saliva compartment AMOUNT to the observed
    # ng/mL saliva CONCENTRATION is not reported anywhere in the paper. There
    # is no supplement to consult: the EuropePMC supplementaryFiles endpoint
    # returns HTTP 404 for PMC13111164 and the core record reports
    # hasSuppl 'N', so this is a reporting gap, not an acquisition gap.
    #
    # PROVENANCE: recovered by digitising Figure 4 (visual predictive check)
    # at 500 dpi, via the saliva quasi-steady state. K30 = 6.24 1/h is a
    # 6.7-minute half-life, so the saliva state tracks plasma essentially
    # instantaneously and
    #     A_saliva = Vmax * Cc / ((Km + Cc) * K30).
    # At the ~2 h peak the model's median plasma prediction is ~790 ng/mL and
    # its median saliva prediction ~2100 ng/mL, giving
    #     A_saliva = 34.7 * 0.790 / ((0.762 + 0.790) * 6.24) = 2.83 mg
    #     vsaliva  = 1000 * 2.83 / 2100 = 1.35 L.
    # Operator-ratified 2026-09-21 (sidecar request-001 / response-001 q1,
    # option B) under the standing figure-digitisation policy. Rejected
    # readings, and the roughly twofold uncertainty this value carries, are
    # recorded in the vignette Errata.
    # ------------------------------------------------------------------
    lvsaliva <- fixed(log(1.35)); label("Apparent saliva volume scaling saliva amount to concentration (L)")    # NOT REPORTED by Huang 2026; recovered from the Figure 4 VPC medians -- see the note above

    etalkel_saliva ~ 0.09   # Table 3 eta_K30 = 30% -> 0.30^2 (RSE 12%, shrinkage 15%; bootstrap 29.6%, 95% CI 22.9-35.6)

    propSd_Csaliva <- 0.284; label("Proportional residual SD for saliva Csaliva (fraction)")                    # Table 3 epsilon_prop = 28.4% (RSE 15%, shrinkage 18%; bootstrap 28.3%, 95% CI 24.5-32.6); saliva residual is proportional only
  })

  model({
    # 1. Reference fat-free mass: the Table 1 cohort median, written into
    #    Equation 6 as the normalising constant.
    ffm_ref <- 30.62  # kg

    # 2. Individual parameters (Equations 5-7).
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (FFM / ffm_ref)^e_ffm_cl
    vc <- exp(lvc + etalvc)

    # 3. Micro-constant.
    kel <- cl / vc

    # 4. Saliva-layer parameters (Equations 8-10). Only K30 carries IIV and a
    #    covariate; Vmax has neither (the data could not support IIV on Vmax
    #    or Km) and Km is held constant.
    vmax_saliva <- exp(lvmax_saliva)
    km_saliva <- exp(lkm_saliva)
    kel_saliva <- exp(lkel_saliva + etalkel_saliva) * (FFM / ffm_ref)^e_ffm_kel_saliva
    vsaliva <- exp(lvsaliva)

    # 5. Plasma concentration. Defined before the ODE block because the
    #    saliva transport term reads it. Cc is in mg/L; multiply by 1000 for
    #    the paper's ng/mL.
    Cc <- central / vc

    # 6. ODE system. Plasma is one-compartment with first-order absorption;
    #    NONMEM ADVAN2 equivalent. Note kel = 15.3/5.77 = 2.65 1/h is an
    #    order of magnitude larger than ka = 0.219 1/h, so the model is
    #    flip-flop and the observed terminal slope reports absorption, not
    #    elimination.
    #
    #    SALIVA IS DRIVEN, NOT MASS-BALANCE-COUPLED. The central equation
    #    carries no -vmax_saliva*Cc/(km_saliva+Cc) loss term, so saliva does
    #    not deplete plasma. Three independent facts fix this reading:
    #      (a) Figure 2 draws this one arrow DASHED and every other arrow
    #          solid, the same notation Nguyen 2026 uses for its explicitly
    #          hypothetical saliva bio-compartment;
    #      (b) the layer was fitted with the plasma parameters held FIXED, so
    #          a term that removed mass from central would have invalidated
    #          them, yet Figure 3A's plasma goodness of fit is unchanged;
    #      (c) the printed constants are not mass-conserving. Vmax = 34.7
    #          mg/h exceeds the cohort's entire average absorption rate
    #          (~9 mg/h at the median 215 mg/day dose), and at a plasma
    #          concentration of 780 ng/mL the transport term would be
    #          34.7*0.78/(0.762+0.78) = 17.6 mg/h against a true elimination
    #          of cl*Cc = 15.3*0.78 = 11.9 mg/h -- salivary loss would
    #          outrun systemic clearance, which the data exclude.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    d/dt(saliva) <- vmax_saliva * Cc / (km_saliva + Cc) - kel_saliva * saliva

    # 7. Saliva concentration. The saliva amount is rescaled by the apparent
    #    saliva volume, which the paper does not report -- see the extended
    #    note on lvsaliva in ini() and the vignette Errata.
    Csaliva <- saliva / vsaliva

    Cc ~ prop(propSd) + add(addSd)
    Csaliva ~ prop(propSd_Csaliva)
  })
}
