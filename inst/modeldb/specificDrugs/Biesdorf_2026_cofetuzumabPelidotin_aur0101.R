Biesdorf_2026_cofetuzumabPelidotin_aur0101 <- function() {
  description <- paste0(
    "Population pharmacokinetic model of the UNCONJUGATED Aur0101 ",
    "auristatin payload (PF-06380101) released from cofetuzumab ",
    "pelidotin (Cofe-P, ABBV-647) in adults with advanced solid tumours ",
    "(Biesdorf 2026, n = 202 pooled from the first-in-human phase 1 ",
    "study B7661001 and the phase 1b NSCLC study M19-611). ",
    "One-compartment linear disposition fed by a first-order rate Kdeg, ",
    "'representing the deconjugation and release of payload from the ",
    "antibody-drug conjugate into systemic circulation', from a depot ",
    "that represents the administered conjugate. Interindividual ",
    "variability on CL, Vc and Kdeg; combined proportional plus ",
    "additive residual error. Body weight and baseline albumin act on ",
    "CL and Vc, liver metastases on CL and the phase 1b study on Vc. ",
    "Because Kdeg (0.243 /day) is far below CL/Vc (about 6.5 /day) the ",
    "payload profile is formation-rate limited and its terminal ",
    "half-life is set by deconjugation. The depot is dosed directly; ",
    "the conjugate model Biesdorf_2026_cofetuzumabPelidotin is not a ",
    "prerequisite. NOTE: the source does not report the drug-antibody ",
    "ratio, the molecular weights, or the payload fraction needed to ",
    "convert a milligram dose of conjugate into the payload mass this ",
    "depot receives, so the depot must be dosed in Aur0101-equivalent ",
    "micrograms; see the vignette Assumptions and deviations section."
  )
  reference <- paste(
    "Biesdorf C, Rinas M, Engelhardt B, Saab R, Guo C, Ferlini C, Freise KJ,",
    "Menon RM, Polepally AR. Population Pharmacokinetics and Exposure-Response",
    "Analyses to Support Cofetuzumab Pelidotin (ABBV-647) Dose Optimization in",
    "Non-Small Cell Lung Cancer Patients. Clin Pharmacol Ther. 2026.",
    "doi:10.1002/cpt.70432. PMCID PMC13473755.",
    "All parameter values are from Table 2, 'Parameter estimates for Aur0101",
    "(unconjugated payload) population pharmacokinetic model'.",
    sep = " "
  )
  vignette <- "Biesdorf_2026_cofetuzumabPelidotin"

  units <- list(
    time = "day",
    dosing = "ug (Aur0101 equivalents entering the depot, NOT milligrams of conjugate)",
    concentration = "ng/mL"
  )

  # CL is in L/day and Vc in L (Table 2). The payload is measured and
  # plotted in ng/mL (LLOQ 0.015 ng/mL; Figures S2 and S4), and the Table 2
  # additive error (variance 0.00276, SD 0.053) is only plausible on that
  # scale, so the depot is dosed in ug: ug / L == ng/mL.
  compartmentData <- list(
    depot = list(
      analyte = "cofetuzumab pelidotin conjugate, as Aur0101 equivalents awaiting deconjugation",
      units = "ug",
      specimen = "not applicable",
      verified = FALSE
    ),
    central = list(
      analyte = "unconjugated Aur0101 (PF-06380101)",
      units = "ug",
      specimen = "serum",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight at baseline.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Normalised to 70 kg inside model(). Biesdorf 2026 normalised continuous covariates to the dataset median, which is not printed; Table S3 gives the mean (SD) 70.3 (18.2) kg, and 70 kg is used as the rounded reference, as in the companion conjugate model. Estimated exponents 0.252 on CL and 0.451 on Vc (Table 2).",
      source_name = "Weight"
    ),
    ALB = list(
      description = "Serum albumin at baseline.",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Normalised to 37 g/L inside model(). The dataset median used by Biesdorf 2026 as the reference is not printed; Table S3 gives the mean (SD) 37.2 (5.16) g/L overall (36.3 in B7661001, 39.2 in M19-611), and 37 g/L is the rounded dataset mean. Exponents 1.79 on CL and 1.90 on Vc (Table 2): higher albumin raises both payload CL and Vc in nearly equal proportion, so it lowers payload exposure (AUC = dose / CL) but barely changes the payload elimination rate CL/Vc.",
      source_name = "Albumin"
    ),
    LMET = list(
      description = "Liver metastases at baseline; 1 = present, 0 = absent.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no liver metastases)",
      notes = "Multiplies payload CL by 0.696 (Table 2 'Liver metastasis effect on CL'). Table S3: 48 of 202 patients (23.8%) had liver metastases (29.2% in B7661001, 12.3% in M19-611).",
      source_name = "Liver metastasis"
    ),
    STUDY_M19611 = list(
      description = "1 = patient from the phase 1b study M19-611 (NCT04189614); 0 = patient from the first-in-human phase 1 study B7661001 (NCT02222922).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (first-in-human study B7661001)",
      notes = "Multiplies payload Vc by 0.828 (Table 2 'Study on Vc'). Table 2 does not name the reference study; it is taken to be coded as in the companion conjugate model, where the phase 1b study carries the factor (established there from the Table 3 median Cavg). See the vignette Assumptions and deviations section.",
      source_name = "Study"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 202L,
    n_studies = 2L,
    n_observations = "3530 unconjugated Aur0101 serum concentrations",
    age_range = "mean 61.1 years (SD 11.0) (Biesdorf 2026 Table S3)",
    weight_range = "mean 70.3 kg (SD 18.2) (Biesdorf 2026 Table S3)",
    sex_female_pct = 69.8,
    disease_state = "Advanced solid tumours (first-in-human B7661001, n = 137: triple-negative breast, NSCLC and ovarian cancer; phase 1b M19-611, n = 65: recurrent PTK7-expressing NSCLC).",
    dose_range = "Cofe-P 0.2-3.7 mg/kg every 3 weeks and 2.1-3.2 mg/kg every 2 weeks (B7661001); 2.8 mg/kg every 3 weeks (M19-611); 60-minute intravenous infusion.",
    liver_metastasis_pct = 23.8,
    albumin = "mean 37.2 g/L (SD 5.16) (Table S3)",
    notes = paste0(
      "Same 202 patients as the conjugate model. Aur0101 was measured ",
      "by LC-MS/MS with an LLOQ of 0.015 ng/mL; concentrations below the ",
      "LLOQ were excluded. NONMEM 7.5, FOCE with interaction."
    )
  )

  ini({
    # ==================================================================
    # Structural disposition, Biesdorf 2026 Table 2 (natural scale with
    # %RSE and 95% CI, hence log()). Typical values are for a 70 kg
    # patient with albumin 37 g/L, no liver metastases, from the
    # first-in-human study.
    #
    # lka is the deconjugation rate the source calls Kdeg ('first-order
    # rate constant representing the deconjugation and release of
    # payload from the antibody-drug conjugate into systemic
    # circulation', Table 2 footnote). Results describe the payload model
    # as 'a one-compartment model with first-order absorption and linear
    # elimination', so Kdeg is the first-order input rate from a depot
    # and is encoded as ka. ln(2)/0.243 = 2.85 days, matching the 2.95-day
    # harmonic-mean payload half-life of the phase 1b study.
    # ==================================================================
    lka <- log(0.243); label("First-order deconjugation (release) rate of Aur0101 from the conjugate, Kdeg (1/day)") # Table 2 'Kdeg (1/day)' 0.243, %RSE 1.39, 95% CI (0.237, 0.250)
    lcl <- log(127);   label("Unconjugated Aur0101 clearance (L/day)")                                                # Table 2 'CL (L/day)' 127, %RSE 3.91, 95% CI (118, 138)
    lvc <- log(19.4);  label("Unconjugated Aur0101 volume of distribution (L)")                                       # Table 2 'Vc (L)' 19.4, %RSE 5.83, 95% CI (17.3, 21.8)

    # ==================================================================
    # Covariate effects: power functions of covariate / reference for
    # continuous covariates, multiplicative factors for categorical ones
    # (Biesdorf 2026 Methods).
    # ==================================================================
    e_wt_cl           <- 0.252; label("Power exponent on (WT/70 kg) for Aur0101 CL (unitless)")                 # Table 2 'Weight effect on CL' 0.252, %RSE 52.4, 95% CI (-0.00701, 0.511)
    e_wt_vc           <- 0.451; label("Power exponent on (WT/70 kg) for Aur0101 Vc (unitless)")                 # Table 2 'Weight effect on Vc' 0.451, %RSE 39.8, 95% CI (0.0987, 0.802)
    e_alb_cl          <- 1.79;  label("Power exponent on (ALB/37 g/L) for Aur0101 CL (unitless)")               # Table 2 'Albumin effect on CL' 1.79, %RSE 12.5, 95% CI (1.36, 2.23)
    e_alb_vc          <- 1.90;  label("Power exponent on (ALB/37 g/L) for Aur0101 Vc (unitless)")               # Table 2 'Albumin effect on Vc' 1.90, %RSE 16.8, 95% CI (1.27, 2.53)
    e_lmet_cl         <- 0.696; label("Multiplicative factor on Aur0101 CL for liver metastases present versus absent (unitless)") # Table 2 'Liver metastasis effect on CL' 0.696, %RSE 7.57, 95% CI (0.600, 0.807)
    e_study_m19611_vc <- 0.828; label("Multiplicative factor on Aur0101 Vc for the phase 1b study M19-611 versus the first-in-human study (unitless)") # Table 2 'Study on Vc' 0.828, %RSE 10.1, 95% CI (0.680, 1.01)

    # ==================================================================
    # Interindividual variability: the first column is omega^2; the
    # footnote %CV = SQRT[exp(omega^2) - 1] x 100 reproduces the printed
    # column (sqrt(exp(0.199) - 1) = 47.0%, sqrt(exp(0.386) - 1) =
    # 68.7%, sqrt(exp(0.0210) - 1) = 14.6%).
    # ==================================================================
    etalcl ~ 0.199   # Table 2 'IIV on CL' 0.199, 47.0 %CV, 7.56% shrinkage
    etalvc ~ 0.386   # Table 2 'IIV on Vc' 0.386, 68.7 %CV, 7.35% shrinkage
    etalka ~ 0.0210  # Table 2 'IIV on Kdeg' 0.0210, 14.6 %CV, 24.0% shrinkage

    # ==================================================================
    # Residual error, read as NONMEM $SIGMA VARIANCES and entered as
    # square roots, on the same grounds as the conjugate model: with
    # N = 3530 observations the %RSE expected of a variance is about
    # 2.4% and of an SD about 1.2%; the printed 3.52% and 5.11% match the
    # variance reading. The additive SD is then sqrt(0.00276) = 0.053
    # ng/mL, about 3.5 x the 0.015 ng/mL LLOQ.
    # ==================================================================
    propSd <- sqrt(0.191);   label("Proportional residual error on Cc (fraction)") # Table 2 'Proportional error' 0.191 (variance), %RSE 3.52, 95% CI (0.178, 0.204)
    addSd  <- sqrt(0.00276); label("Additive residual error on Cc (ng/mL)")        # Table 2 'Additive error' 0.00276 (variance), %RSE 5.11, 95% CI (0.00248, 0.00304)
  })

  model({
    # ----- Individual parameters -----
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (ALB / 37)^e_alb_cl * e_lmet_cl^LMET
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc * (ALB / 37)^e_alb_vc * e_study_m19611_vc^STUDY_M19611
    ka <- exp(lka + etalka)

    # ----- Micro-constant -----
    kel <- cl / vc

    # ----- One-compartment disposition with first-order deconjugation input -----
    # depot holds the Aur0101 still bound to the conjugate; dose it with
    # the Aur0101-equivalent mass in ug, not with the mg dose of conjugate.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # ----- Observation -----
    # Cc is the unconjugated Aur0101 serum concentration in ng/mL.
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
