Yamamoto_2023_tominersen <- function() {
  description <- paste(
    "Six-compartment population PK model for tominersen (2'-O-methoxyethyl",
    "antisense oligonucleotide targeting huntingtin mRNA) following intrathecal",
    "lumbar-puncture administration in adults with Huntington's disease",
    "(Yamamoto 2023). Two coupled three-compartment subsystems: the intrathecal",
    "bolus enters the central CSF compartment, which exchanges with two",
    "CNS-tissue compartments and drains unidirectionally into the plasma central",
    "compartment via CL_CSF (bioavailability F1 and F2 both fixed to 1); the",
    "plasma central compartment has two peripheral compartments and first-order",
    "elimination. Baseline CSF total protein, age, and antidrug antibodies are",
    "covariates on CL_CSF; body weight is an allometric power covariate on all",
    "plasma clearances and volumes, and antidrug antibodies and female sex are",
    "covariates on plasma CL. Residual error is additive on the log scale",
    "(lognormal) for both CSF and plasma, each with interindividual variability",
    "on the residual magnitude."
  )
  reference <- paste(
    "Yamamoto Y, Sanwald Ducray P, Bjornsson M, Smart K, Grimsey P,",
    "Vatakuti S, Portron A, Massonnet B, Norris DA, Silber Baumann HE.",
    "Development of a population pharmacokinetic model to characterize the",
    "pharmacokinetics of intrathecally administered tominersen in",
    "cerebrospinal fluid and plasma.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(9):1213-1226.",
    "doi:10.1002/psp4.13001",
    sep = " "
  )
  vignette <- "Yamamoto_2023_tominersen"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  # The two peripheral CSF compartments are named for CNS tissue rather than
  # peripheral1 / peripheral2 because those two canonical names are already
  # taken by the plasma subsystem's peripheral compartments in this same model.
  # Yamamoto 2023 Table 2 / Figure 3 label them "V in the first / second
  # peripheral CSF compartment"; the Appendix S1 NONMEM control stream names
  # them COMP=(BRAIN) and COMP=(BRAIN2), and the Discussion attributes the long
  # CSF terminal half-life to "a slow rate of ASOs dissociating from the CNS
  # tissue". The numbering follows the canonical peripheral1 / peripheral2
  # pattern. Extends the single cns_tissue compartment declared by the sibling
  # intrathecal-ASO model Luu_2017_nusinersen.R, which has only one.
  paper_specific_compartments <- c("cns_tissue1", "cns_tissue2")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric power scaling on the plasma-subsystem parameters only,",
        "normalised to a 75 kg reference (Appendix S1 $PK: `(WT/75)**THETA(1)`",
        "for CL_plasma, Q2_plasma and Q3_plasma, and `(WT/75)**THETA(2)` for",
        "V1_plasma, V2_plasma and V3_plasma). Table 2 reports these two",
        "exponents as BW_CLs = 0.687 and BW_Vs = 0.866. No body-weight effect",
        "was retained on any CSF-subsystem parameter. Analysis-population range",
        "40.1 to 116 kg, median 69.1 to 75.3 kg across the contributing studies",
        "(Table S1); the 5th and 95th percentiles quoted in the Results are 52",
        "and 95.5 kg."
      ),
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Exponential (centered) covariate on CSF clearance, referenced to 49",
        "years (Appendix S1 $PK: `CLCAGE = EXP(THETA(17)*(AGE - 49))`; the",
        "49-year reference is corroborated by the Figure 5 typical-profile",
        "settings). Analysis-population range 25 to 66 years, median 46 to 50",
        "(Table S1); the 5th and 95th percentiles quoted in the Results are 31",
        "and 64 years."
      ),
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "male (SEXF = 0)",
      notes              = paste(
        "Fractional-difference covariate on plasma clearance. Appendix S1 $PK",
        "codes `IF(SEX.EQ.0) CLSEX = 1` and `IF(SEX.EQ.1) CLSEX = (1 +",
        "THETA(19))`; THETA(19) is negative (-0.186) and the Results state that",
        "CL_plasma was decreased by approximately 19% in women, so the source",
        "SEX = 1 level is female and maps directly onto the canonical SEXF",
        "coding with no value inversion. Population 350 of 750 female (46.7%,",
        "Table S1)."
      ),
      source_name        = "SEX"
    ),
    ADA_POS = list(
      description        = "Antidrug-antibody positivity, 1 = ADA-positive, 0 = ADA-negative",
      units              = "(binary)",
      type               = "binary",
      reference_category = "ADA-negative (ADA_POS = 0)",
      notes              = paste(
        "Time-varying fractional-difference covariate applied to BOTH CSF",
        "clearance and plasma clearance (Appendix S1 $PK: `IF(ADA.EQ.1) CLCADA",
        "= 1+THETA(16)` and `IF(ADA.EQ.1) CLADA = 1+THETA(18)`). ADA status was",
        "carried per concentration record, not per subject (Tables S3A and",
        "S3B); ADAs were not measured in the phase I/IIa study or in GEN-PEAK,",
        "so those records were ADA-missing and took the ADA-negative reference",
        "multiplier of 1."
      ),
      source_name        = "ADA"
    ),
    CSF_TPRO = list(
      description        = "Total protein concentration in cerebrospinal fluid",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Exponential (centered) covariate on CSF clearance, referenced to",
        "0.35 g/L (Appendix S1 $PK: `CLCTPCSF = EXP(THETA(15)*(TPCSF -",
        "0.35))`; the 0.35 g/L reference is corroborated by the Figure 5",
        "typical-profile settings). The control stream additionally disables the",
        "effect (multiplier = 1) for implausible values above 2 g/L and for the",
        "-99 missing-data sentinel; the >2 g/L guard is reproduced here, the -99",
        "sentinel is not (see the vignette's Assumptions and deviations).",
        "Analysis-population median 0.305 to 0.372 g/L across the contributing",
        "studies, overall range 0 to 15.5 g/L with the 15.5 g/L maximum in",
        "GEN-EXTEND (Table S1); the 5th and 95th percentiles quoted in the",
        "Results are 0.19 and 0.54 g/L."
      ),
      source_name        = "TPCSF"
    )
  )

  # Covariates that Yamamoto 2023 screened in the stepwise covariate model
  # (Methods, "The stepwise covariate model building procedure") but did NOT
  # retain in the final model. Documentation only -- these are deliberately not
  # referenced in model(). Only the three with canonical register entries are
  # listed here; the remaining screened covariates were Huntington's-specific
  # or per-administration measures with no canonical column name yet (CAG
  # repeat length, CAG age-product score, caudate volume, ventricle volume,
  # whole-brain volume, total protein in blood, the volume of CSF withdrawn
  # before dosing, and the lumbar site / route of intrathecal administration)
  # and are described in population$notes instead.
  covariatesDataExcluded <- list(
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = paste(
        "Screened on the CSF model parameters and not retained. Range 108 to",
        "202 cm, median 170 to 173 cm (Table S1)."
      ),
      source_name = "HT"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = paste(
        "Screened on the plasma model parameters and not retained; the",
        "Discussion notes that hepatic insufficiency is expected to have little",
        "effect on oligonucleotide clearance. Range 4.88 to 114 U/L, median",
        "13.5 to 17.0 U/L (Table S1)."
      ),
      source_name = "ALT"
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Screened on the plasma model parameters and not retained. The",
        "Discussion reports the analysed range as 51 to 150 mL/min; Table S1",
        "spans 50.6 to 168.7 mL/min with medians of 98.4 to 108 mL/min. Urine",
        "is the primary elimination route for antisense oligonucleotides, but",
        "no renal-function effect was detectable over this range."
      ),
      source_name = "CCRCLT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 750L,
    n_studies      = 5L,
    studies        = paste(
      "Pooled: phase I/IIa (NCT02519036), its open-label extension",
      "(NCT03342053), GEN-PEAK (NCT04000594), GENERATION HD1 (NCT03761849),",
      "GEN-EXTEND (NCT03842969)"
    ),
    age_range      = "25 to 66 years (5th-95th percentiles 31 to 64; reference 49)",
    weight_range   = "40.1 to 116 kg (5th-95th percentiles 52 to 95.5; reference 75)",
    sex_female_pct = 46.7,
    disease_state  = "Adults with early manifest or manifest Huntington's disease",
    dose_range     = "10 to 120 mg intrathecal bolus, Q4W / Q8W / Q16W, for up to 25 months",
    administration_routes = "Intrathecal bolus, almost all by lumbar puncture (L3-L4 in 60.7 to 67.6% of doses, Table S2); the GEN-PEAK doses were given via an indwelling intrathecal catheter",
    renal_function = "Creatinine clearance 50.6 to 168.7 mL/min (Table S1); the Discussion quotes the analysed range as 51 to 150 mL/min. Renal function was not a significant covariate",
    notes          = paste(
      "6302 CSF and 5454 plasma concentrations (Yamamoto 2023 Table 1). 4% of",
      "CSF and 19% of plasma samples were below the limit of quantification and",
      "were excluded from the fit (the M3 method was attempted for plasma but",
      "did not converge). Baseline demographics are in Table S1, per-dose",
      "administration covariates in Table S2, and time-varying ADA in Tables",
      "S3A and S3B. Estimation with NONMEM 7.4.4 FOCE-I (Appendix S1",
      "$ESTIMATION METHOD=1 INTER); covariate selection by PsN stepwise",
      "covariate modelling (SCM) with adaptive scope reduction. The CSF",
      "subsystem was fitted first and its parameters were then held fixed while",
      "the plasma subsystem was fitted (Methods, 'PopPK model development').",
      "Covariates screened but NOT retained: on the CSF model, age (retained),",
      "sex, CAG repeat length, CAG age-product score, caudate volume, ventricle",
      "volume, whole-brain volume and height; on the plasma model, age, sex",
      "(retained), alanine aminotransferase, total protein in plasma and",
      "creatinine clearance. Route of administration, lumbar site of",
      "administration and the volume of CSF withdrawn before dosing were also",
      "tested and not retained. Dose was tested as a covariate on the CSF",
      "parameters and was not significant, supporting a linear CSF model."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # CSF subsystem structural parameters.
    # The CSF subsystem was fitted first, so in the final combined control
    # stream these carry the CSF fit's final estimates as `FIX` values
    # (Appendix S1 $THETA 3-8). Table 2 reports the same numbers rounded to
    # three significant figures, each with a relative standard error, so they
    # are estimated parameters and are NOT wrapped in fixed() here. The full
    # precision below is the control stream's.
    #
    # The bioavailabilities F1 (into CSF) and F2 (CSF to plasma) were both
    # fixed to 1 (Table 2), so neither carries a parameter.
    # ------------------------------------------------------------------
    lcl_csf <- log(0.0176577)   ; label("CSF clearance CL_CSF (unidirectional CSF-to-plasma transfer; L/h)")   # Table 2 CL_CSF 0.0177 L/h (RSE 7.96%); Appendix S1 $THETA(3) 0.0176577 FIX
    lvcsf   <- log(0.0436042)   ; label("Central CSF volume V1_CSF (L)")                                       # Table 2 V1,CSF 0.0436 L (RSE 33.1%); Appendix S1 $THETA(4) 0.0436042 FIX
    lqcsf   <- log(0.00982466)  ; label("CSF inter-compartmental clearance Q2_CSF (csf <-> cns_tissue1; L/h)")  # Table 2 Q2,CSF 0.00982 L/h (RSE 27.3%); Appendix S1 $THETA(5) 0.00982466 FIX
    lvcns   <- log(0.0517384)   ; label("First peripheral CSF (CNS-tissue) volume V2_CSF (L)")                  # Table 2 V2,CSF 0.0517 L (RSE 11.8%); Appendix S1 $THETA(6) 0.0517384 FIX
    lqcsf2  <- log(1.20953e-05) ; label("CSF inter-compartmental clearance Q3_CSF (csf <-> cns_tissue2; L/h)")  # Table 2 Q3,CSF 0.0000121 L/h (RSE 15.7%); Appendix S1 $THETA(7) 1.20953E-05 FIX
    lvcns2  <- log(0.0117556)   ; label("Second peripheral CSF (CNS-tissue) volume V3_CSF (L)")                 # Table 2 V3,CSF 0.0118 L (RSE 15.8%); Appendix S1 $THETA(8) 0.0117556 FIX

    # ------------------------------------------------------------------
    # Plasma subsystem structural parameters. These were the parameters being
    # estimated in the final run, so Appendix S1 $THETA 9-14 lists only their
    # INITIAL estimates (9.95883, 41.1879, 14.2579, 273.266, 0.276129,
    # 619.497). The final values below are therefore the Table 2 point
    # estimates.
    # ------------------------------------------------------------------
    lcl  <- log(10.7) ; label("Plasma clearance CL_plasma at 75 kg (L/h)")                                  # Table 2 CL_plasma 10.7 L/h (RSE 3.30%)
    lvc  <- log(41.0) ; label("Central plasma volume V1_plasma at 75 kg (L)")                               # Table 2 V1,plasma 41.0 L (RSE 2.88%)
    lq   <- log(13.4) ; label("Plasma inter-compartmental clearance Q2_plasma at 75 kg (L/h)")              # Table 2 Q2,plasma 13.4 L/h (RSE 7.98%)
    lvp  <- log(283)  ; label("First peripheral plasma volume V2_plasma at 75 kg (L)")                      # Table 2 V2,plasma 283 L (RSE 6.15%)
    lq2  <- log(0.267); label("Second plasma inter-compartmental clearance Q3_plasma at 75 kg (L/h)")       # Table 2 Q3,plasma 0.267 L/h (RSE 4.84%)
    lvp2 <- log(613)  ; label("Second peripheral plasma volume V3_plasma at 75 kg (L)")                     # Table 2 V3,plasma 613 L (RSE 6.35%)

    # ------------------------------------------------------------------
    # Covariate effects on CSF clearance. Continuous covariates use the
    # exponential centered form and ADA a fractional difference to the most
    # common category (Yamamoto 2023 Equation 2; Appendix S1 $PK). Appendix S1
    # carries these as FIX for the same sequential-fitting reason as the CSF
    # structural parameters, and Table 2 reports an RSE for each, so they are
    # estimated and not wrapped in fixed().
    # ------------------------------------------------------------------
    e_tpcsf_cl_csf <- -0.539797 ; label("Exponential coefficient of CSF total protein on CL_CSF (per g/L, centered at 0.35 g/L)")  # Table 2 Total protein_CL_CSF -0.540 (RSE 11.2%); Appendix S1 $THETA(15) -0.539797 FIX
    e_ada_cl_csf   <- -0.100697 ; label("Fractional change in CL_CSF when ADA-positive (unitless)")                                # Table 2 ADA_CL_CSF -0.101 (RSE 6.89%); Appendix S1 $THETA(16) -0.100697 FIX
    e_age_cl_csf   <- -0.005498 ; label("Exponential coefficient of age on CL_CSF (per year, centered at 49 years)")               # Table 2 Age_CL_CSF -0.00550 (RSE 12.8%); Appendix S1 $THETA(17) -0.005498 FIX

    # ------------------------------------------------------------------
    # Covariate effects on the plasma subsystem. The two body-weight exponents
    # are shared: THETA(1) across CL_plasma, Q2_plasma and Q3_plasma, and
    # THETA(2) across V1_plasma, V2_plasma and V3_plasma. Appendix S1 lists
    # only their initial estimates (1.08182 and 0.913675), as it does for the
    # ADA and sex effects (-0.672001 and -0.001), so all four final values
    # below come from Table 2.
    # ------------------------------------------------------------------
    e_wt_cl   <- 0.687  ; label("Allometric exponent of body weight on the plasma clearances CL, Q2, Q3 (unitless)")  # Table 2 BW_CLs 0.687 (RSE 2.60%); Appendix S1 $PK `(WT/75)**THETA(1)`
    e_wt_vc   <- 0.866  ; label("Allometric exponent of body weight on the plasma volumes V1, V2, V3 (unitless)")     # Table 2 BW_Vs 0.866 (RSE 1.54%); Appendix S1 $PK `(WT/75)**THETA(2)`
    e_ada_cl  <- -0.673 ; label("Fractional change in CL_plasma when ADA-positive (unitless)")                        # Table 2 ADA_CL_plasma -0.673 (RSE 0.639%); Appendix S1 $THETA(18)
    e_sexf_cl <- -0.186 ; label("Fractional change in CL_plasma in women (unitless)")                                 # Table 2 Sex_CL_plasma -0.186 (RSE 12.1%); Appendix S1 $THETA(19)

    # ------------------------------------------------------------------
    # Interindividual variability. Exponential IIV (Yamamoto 2023 Equation 1).
    # Table 2's IIV percentages are sqrt(omega^2) * 100, verified against the
    # FIXED CSF $OMEGA values in Appendix S1: sqrt(0.0273094) = 0.1653 against
    # the reported 16.5%, sqrt(0.312055) = 0.5586 against 55.9%, and
    # sqrt(0.199008) = 0.4461 against 44.6%. They are NOT the lognormal
    # sqrt(exp(omega^2) - 1) CV. CSF variances are therefore taken at full
    # precision from Appendix S1; plasma variances are squared from the Table 2
    # percentages because Appendix S1 lists only the plasma initial estimates.
    #
    # No IIV was carried on V2_CSF, Q3_CSF, V3_CSF, V1_plasma or Q3_plasma
    # (Appendix S1 $OMEGA 5, 6, 7, 9 and 12 are all `0 FIX`, and Table 2 lists
    # no IIV row for them), so no eta is defined for those parameters.
    # ------------------------------------------------------------------
    etalcl_csf ~ 0.0273094 # Table 2 IIV on CL_CSF 16.5% (RSE 3.22%, shrinkage 6.65%); Appendix S1 $OMEGA(2) 0.0273094 FIX
    etalvcsf   ~ 0.312055  # Table 2 IIV on V1,CSF 55.9% (RSE 41.6%, shrinkage 88.4%); Appendix S1 $OMEGA(3) 0.312055 FIX
    etalqcsf   ~ 0.199008  # Table 2 IIV on Q2,CSF 44.6% (RSE 51.2%, shrinkage 88.7%); Appendix S1 $OMEGA(4) 0.199008 FIX
    etalcl     ~ 0.082944  # Table 2 IIV on CL_plasma 28.8% (RSE 2.84%, shrinkage 21.9%) -> 0.288^2; Appendix S1 $OMEGA(8) initial 0.090009
    etalq      ~ 5.7121    # Table 2 IIV on Q2,plasma 239% (RSE 6.35%, shrinkage 35.0%) -> 2.39^2; Appendix S1 $OMEGA(10) initial 5.82967
    etalvp     ~ 0.416025  # Table 2 IIV on V2,plasma 64.5% (RSE 9.33%, shrinkage 52.0%) -> 0.645^2; Appendix S1 $OMEGA(11) initial 0.425571
    etalvp2    ~ 0.367236  # Table 2 IIV on V3,plasma 60.6% (RSE 8.39%, shrinkage 44.3%) -> 0.606^2; Appendix S1 $OMEGA(13) initial 0.367345

    # ------------------------------------------------------------------
    # Interindividual variability on the residual magnitude. Appendix S1
    # $ERROR uses `Y = IPRED + EPS(n) * EXP(ETA(k))`, i.e. the log-scale
    # residual SD is scaled per subject by exp(eta). Encoded here as
    # expSd_i <- expSd * exp(etaexpSd), following the Chandasana 2024
    # etapropSd precedent.
    # ------------------------------------------------------------------
    etaexpSd_Ccsf ~ 0.214278 # Table 2 IIV RUV CSF 46.3% (RSE 2.84%, shrinkage 10.1%); Appendix S1 $OMEGA(1) 0.214278 FIX
    etaexpSd      ~ 0.0841   # Table 2 IIV RUV plasma 29.0% (RSE 3.67%, shrinkage 14.6%) -> 0.29^2; Appendix S1 $OMEGA(14) initial 0.088158

    # ------------------------------------------------------------------
    # Residual unexplained variability. The analysis was performed on
    # log-transformed concentrations with an additive residual error
    # (Yamamoto 2023 Methods, "Residual error for tominersen in CSF and
    # plasma"), which is a lognormal residual on the linear scale. Table 2
    # reports each sigma as sqrt(sigma^2) * 100, verified against the FIXED
    # CSF $SIGMA in Appendix S1: sqrt(0.0975477) = 0.3123 against the
    # reported 31.2%.
    # ------------------------------------------------------------------
    expSd_Ccsf <- sqrt(0.0975477); label("CSF log-scale residual SD (typical value before the per-subject exp(etaexpSd_Ccsf) modifier)") # Table 2 RUV CSF 31.2% (RSE 2.38%, shrinkage 1.61%); Appendix S1 $SIGMA(1) 0.0975477 FIX
    expSd      <- 0.536          ; label("Plasma log-scale residual SD (typical value before the per-subject exp(etaexpSd) modifier)")   # Table 2 RUV plasma 53.6% (RSE 1.94%, shrinkage 7.97%); Appendix S1 $SIGMA(2) initial 0.291491
  })

  model({
    # ---- 1. Covariate multipliers -------------------------------------
    # CSF clearance: exponential (centered) effects of CSF total protein and
    # age, and a fractional-difference effect of ADA positivity
    # (Appendix S1 $PK, "CSF model" block). The control stream sets the
    # total-protein multiplier to 1 for implausible CSF protein above 2 g/L
    # (Table S1 shows a GEN-EXTEND maximum of 15.5 g/L).
    cov_tpcsf_cl_csf <- ifelse(CSF_TPRO > 2, 1, exp(e_tpcsf_cl_csf * (CSF_TPRO - 0.35)))
    cov_age_cl_csf   <- exp(e_age_cl_csf * (AGE - 49))
    cov_ada_cl_csf   <- 1 + e_ada_cl_csf * ADA_POS

    # Plasma: allometric body-weight scaling on all clearances and volumes,
    # plus fractional-difference ADA and sex effects on CL only.
    wt_cl <- (WT / 75)^e_wt_cl
    wt_v  <- (WT / 75)^e_wt_vc

    # ---- 2. Individual parameters -------------------------------------
    cl_csf <- exp(lcl_csf + etalcl_csf) * cov_tpcsf_cl_csf * cov_age_cl_csf * cov_ada_cl_csf
    vcsf   <- exp(lvcsf + etalvcsf)
    qcsf   <- exp(lqcsf + etalqcsf)
    vcns   <- exp(lvcns)
    qcsf2  <- exp(lqcsf2)
    vcns2  <- exp(lvcns2)

    cl  <- exp(lcl + etalcl) * wt_cl * (1 + e_ada_cl * ADA_POS) * (1 + e_sexf_cl * SEXF)
    vc  <- exp(lvc)            * wt_v
    q   <- exp(lq + etalq)     * wt_cl
    vp  <- exp(lvp + etalvp)   * wt_v
    q2  <- exp(lq2)            * wt_cl
    vp2 <- exp(lvp2 + etalvp2) * wt_v

    # ---- 3. Micro-constants (Appendix S1 $PK) -------------------------
    # Reproduces the control stream's K14 / K12 / K21 / K13 / K31 (CSF) and
    # K40 / K45 / K54 / K46 / K64 (plasma). K14 is the one-way CSF-to-plasma
    # transfer; there is no return path from plasma to CSF.
    k14 <- cl_csf / vcsf
    k12 <- qcsf   / vcsf
    k21 <- qcsf   / vcns
    k13 <- qcsf2  / vcsf
    k31 <- qcsf2  / vcns2

    k40 <- cl / vc
    k45 <- q  / vc
    k54 <- q  / vp
    k46 <- q2 / vc
    k64 <- q2 / vp2

    # ---- 4. ODE system (Appendix S1 $DES) -----------------------------
    # Compartment order matches the control stream's A(1)..A(6):
    #   1 csf (CSF, DEFDOSE), 2 cns_tissue1 (BRAIN), 3 cns_tissue2 (BRAIN2),
    #   4 central (CENTRAL), 5 peripheral1 (PERIPH1), 6 peripheral2 (PERIPH2).
    # States carry amounts in mg, the dosing unit; F1 = F2 = 1 so there is no
    # bioavailability term.
    d/dt(csf)         <- k21 * cns_tissue1 + k31 * cns_tissue2 - (k14 + k12 + k13) * csf
    d/dt(cns_tissue1) <- k12 * csf - k21 * cns_tissue1
    d/dt(cns_tissue2) <- k13 * csf - k31 * cns_tissue2
    d/dt(central)     <- k14 * csf + k54 * peripheral1 + k64 * peripheral2 - (k40 + k45 + k46) * central
    d/dt(peripheral1) <- k45 * central - k54 * peripheral1
    d/dt(peripheral2) <- k46 * central - k64 * peripheral2

    # ---- 5. Observed concentrations -----------------------------------
    # Reproduces the control stream's `S1 = VC1/1000` and `S4 = V1/1000`
    # scaling: an amount in mg divided by a volume in L gives mg/L, and
    # dividing the volume by 1000 converts that to ug/L = ng/mL.
    Ccsf <- csf     / (vcsf / 1000)
    Cc   <- central / (vc / 1000)

    # ---- 6. Residual error -------------------------------------------
    # Additive residual error on log-transformed concentrations, with the
    # per-subject exp(eta) modifier on the residual SD from Appendix S1 $ERROR.
    expSd_Ccsf_i <- expSd_Ccsf * exp(etaexpSd_Ccsf)
    expSd_i      <- expSd      * exp(etaexpSd)

    Ccsf ~ lnorm(expSd_Ccsf_i)
    Cc   ~ lnorm(expSd_i)
  })
}
