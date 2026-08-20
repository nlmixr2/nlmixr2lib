Sinha_2026_oxcarbazepine <- function() {
  description <- "Joint parent-metabolite population PK model for oral oxcarbazepine (OXZ) and its active 10-monohydroxy derivative (MHD) in children with and without obesity, aged 44 days to 20.9 years (Sinha 2026). One-compartment OXZ plus one-compartment MHD with first-order absorption, complete conversion of OXZ to MHD (the OXZ elimination clearance IS the MHD formation clearance), molar-mass-corrected bidirectional mass transfer, and a first-order MHD-to-OXZ back-transformation. All clearances and volumes are allometrically scaled to pharmacokinetic weight (PKWT), a fat-free-mass-based body-size descriptor, referenced to 50 kg; body size was the only covariate retained."
  reference   <- "Sinha J, Zimmerman K, Balevic SJ, Hornik C, Muller WJ, Rathore M, Meyer M, Finkelstein Y, Al-Uzri A, Lakhotia A, Goldstein S, Chen JY, Anand R, Gonzalez D. Population Pharmacokinetic Modeling of Oxcarbazepine and Its Active Metabolite 10-Monohydroxy Derivative to Inform Dosing in Children with Obesity. Clin Pharmacokinet. 2026;65:329-344. doi:10.1007/s40262-025-01579-0. Correction: Clin Pharmacokinet. 2026;65:345. doi:10.1007/s40262-025-01613-1 (Eq. 10 was printed as TV_V,MHD = 50*(PKWT/50)^0.752 but should read TV_V,OXZ = 33.1*(PKWT/50)^0.752; no parameter value changed)."
  vignette    <- "Sinha_2026_oxcarbazepine"
  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass used as the body-size descriptor for allometric scaling of every clearance and volume in the model. The source paper calls this quantity 'pharmacokinetic weight' (PKWT).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference value 50 kg (Sinha 2026 ESM2 control stream: every scaling term is (PKWT/50)^theta; Discussion, 'a reference adult with 50 kg FFM').",
        "PKWT is a COMPOSITE fat-free-mass metric defined in Sinha 2026 Methods 2.3.2: FFM for children and young adults aged >= 3 years, and total body weight (WT) for children younger than 3 years.",
        "The WT fallback below 3 years exists because the Al-Sallami et al. paediatric FFM equation was developed in children >= 3 years, so the authors assumed a WT-equivalent FFM below that age.",
        "FFM itself is computed from WT, BMI, postnatal age (PNA, years) and sex:",
        "for 18 <= PNA < 21 years by Janmahasatian et al. (ESM1 Eqs. S1-S2) as 9270*WT/(6680+216*BMI) for males and 9270*WT/(8780+244*BMI) for females;",
        "for 3 <= PNA < 18 years by Al-Sallami et al. (ESM1 Eqs. S3-S4) as those adult forms multiplied by (0.88 + 0.12/(1+(PNA/13.4)^-12.7)) for males and (1.11 - 0.11/(1+(PNA/7.1)^-1.1)) for females.",
        "NOTE the sign in the female paediatric multiplier: ESM1 Eq. S4 prints '+0.11', which cannot be correct because the paediatric multiplier must converge to 1 at adolescence so that Eq. S4 reduces to the adult Eq. S2, and because the male numerator is exactly (1 - 0.88) = 0.12. The corresponding female numerator is (1 - 1.11) = -0.11. The vignette asserts this convergence numerically.",
        "Exponents applied to (FFM/50): 1 (fixed) on OXZ clearance, 0.671 on MHD clearance, and a shared 0.752 on both volumes.",
        "The model consumes FFM/PKWT directly as a data column; the derivation above is reproduced in the validation vignette."
      ),
      source_name        = "PKWT"
    )
  )

  # Covariates that Sinha 2026 screened but did not retain in the final model
  # (Table S5 of ESM1 lists the univariate AIC for each; none improved on the
  # PKWT-only model M6, whose AIC was 6698.54). Documented here so the full
  # screening set is preserved without declaring covariates the model never
  # references.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight (screened as the body-size descriptor, superseded by PKWT)",
      units       = "kg",
      type        = "continuous",
      notes       = "Sinha 2026 Table 3: the WT-based model M3 had AIC 6705.62 versus 6698.54 for the PKWT-based final model M6, an ~7-point improvement that the Discussion cites as the reason PKWT was preferred. WT is still required upstream to compute PKWT (see covariateData$FFM notes). Cohort median 27.2 kg, range 3.825-156.9 kg."
    ),
    PNA = list(
      description = "Postnatal age (screened, not retained)",
      units       = "years",
      type        = "continuous",
      notes       = "Sinha 2026 Table S5: model M8a (PKWT + PNA) had AIC 6699.51 versus 6698.54 for M6, so PNA was not retained. Cohort median 9 years, range 44 days to 20.90 years. PNA is still required upstream to compute PKWT."
    ),
    PMA = list(
      description = "Postmenstrual age, tested as a sigmoidal (Hill) maturation function on clearance (screened, not retained)",
      units       = "weeks",
      type        = "continuous",
      notes       = "Sinha 2026 Table S5: model M8b (PKWT + PMA Hill maturation, estimating TM50 and the Hill coefficient gamma) had AIC 6701.14 versus 6698.54 for M6. The Discussion attributes the absent maturation signal to the small number of infants (6% of the cohort; no neonates) and to early maturation of the glucuronidation pathway. Gestational age was imputed at 40 weeks for 96% of participants."
    ),
    SCR = list(
      description = "Serum creatinine (screened, not retained)",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Sinha 2026 Table S5: model M9a (PKWT + SCR) had AIC 6700.15 versus 6698.54 for M6. The Discussion attributes the null result to under-representation of renal impairment (SCR IQR 0.23-0.63 mg/dL). Table 2 medians: 0.42 mg/dL (AED01), 0.40 mg/dL (POP01)."
    ),
    EGFR = list(
      description = "Estimated glomerular filtration rate (screened, not retained)",
      units       = "mL/min/1.73m^2",
      type        = "continuous",
      notes       = "Sinha 2026 Table S5: model M9b (PKWT + eGFR) had AIC 6701.47 versus 6698.54 for M6. Discussion reports the cohort eGFR IQR as 95-156 mL/min/1.73 m^2, i.e. no renal impairment to detect."
    ),
    CONMED_PB = list(
      description = "Concomitant phenobarbital, the only enzyme-inducing antiepileptic drug testable in this dataset (screened, not retained)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Sinha 2026 Table S5: model M10a (PKWT + phenobarbital) had AIC 6700.34 versus 6698.54 for M6. Only 6 of 100 participants received phenobarbital; no participant received carbamazepine, and phenytoin/fosphenytoin exposure was unknown for POP01, so phenobarbital was the only EIAED assessable. The Discussion notes that previous population PK studies reported a 17-35% increase in MHD clearance with concomitant EIAEDs, which this dataset was underpowered to detect."
    ),
    CONMED_VPA = list(
      description = "Concomitant valproic acid (screened, not retained)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Sinha 2026 Table S5: model M10b (PKWT + valproic acid) had AIC 6700.34 versus 6698.54 for M6. Seven of 100 participants received valproic acid."
    ),
    CONMED_LEV = list(
      description = "Concomitant levetiracetam (screened, not retained)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Sinha 2026 Table S5: model M10c (PKWT + levetiracetam) had AIC 6700.34 versus 6698.54 for M6. Ten of 100 participants received levetiracetam."
    ),
    CONMED_TPM = list(
      description = "Concomitant topiramate (screened, not retained)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Sinha 2026 Table S5: model M10d (PKWT + topiramate) had AIC 6700.34 versus 6698.54 for M6. Eight of 100 participants received topiramate."
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "oxcarbazepine", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "oxcarbazepine", units = "mg", specimen = "plasma", verified = TRUE),
    central_mhd = list(analyte = "10-monohydroxy derivative (MHD) of oxcarbazepine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 100,
    n_studies      = 2,
    age_range      = "44 days to 20.90 years (postnatal age)",
    age_median     = "9 years",
    weight_range   = "3.825-156.9 kg",
    weight_median  = "27.2 kg",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported in the publication (race, ethnicity, and race-white/black were carried in the analysis dataset per the ESM2 $INPUT record but no summary is tabulated).",
    disease_state  = "Children and young adults receiving oxcarbazepine as standard of care, predominantly for epilepsy / partial-onset seizures. 52 of 100 participants (52%) had obesity, defined as a body-mass-index-for-age at or above the 95th percentile of the CDC growth charts in participants >= 2 years of age.",
    dose_range     = "Standard-of-care oxcarbazepine dosing (immediate-release tablet or suspension); doses were not protocol-assigned. Label-recommended maintenance target doses used in the paper's simulations are 60 mg/kg/day for ages 2-4 years, and 900, 1200, and 1800 mg/day for ages > 4 years weighing 20-29, > 29-39, and > 39 kg respectively, all twice daily (Table 1).",
    regions        = "United States and Canada (Pediatric Trials Network sites).",
    n_observations = "425 plasma concentrations: 212 oxcarbazepine and 213 MHD. Median (range) 3 (1-7) repeated observations per participant for each analyte. No concentration was below the 1 ng/mL lower limit of quantification.",
    obesity_status = "52% with obesity (BMI-for-age >= 95th CDC percentile); 25 from POP01 and 27 from AED01. AED01 enrolled only participants with obesity.",
    co_medication  = "Phenobarbital (n = 6), valproic acid (n = 7), topiramate (n = 8), and levetiracetam (n = 10) were co-administered on the dosing day or the prior day. No participant received carbamazepine. Phenytoin/fosphenytoin status was known only for AED01, where none was administered.",
    notes          = "Pooled from two multicentre prospective standard-of-care PK studies: POP01 (NICHD-2011-POP01, NCT01431326, n = 73, any patient < 21 years with or without obesity) and AED01 (NICHD-2015-AED01, NCT02993861, n = 27, ages 2-18 years, all with obesity). Demographics from Sinha 2026 Table 2 and Results 3.1. Gestational-age records existed for only four POP01 participants (37.6-41 weeks), so postmenstrual age was computed from an imputed 40-week gestational age for 96% of the pooled dataset."
  )

  ini({
    # Structural parameters - final model M6 of Sinha 2026, reported in Table 4
    # and reproduced in ESM1 Table S3 (row M6) and the ESM2 NONMEM control
    # stream. Reference body size is 50 kg PKWT for every scaling term.
    #
    # There is no bioavailability parameter in the source control stream:
    # oxcarbazepine is completely absorbed (Introduction, citing reference [5]),
    # so F is implicitly 1 and the clearances / volumes below are absolute
    # rather than apparent.
    lka          <- log(0.269);     label("First-order absorption rate constant for OXZ (1/h)")                    # Sinha 2026 Table 4: Ka = 0.269 1/h, RSE 28%; ESM2 $THETA(1)
    lcl          <- log(220);       label("OXZ elimination clearance at 50 kg PKWT (L/h)")                         # Sinha 2026 Table 4: CL_OXZ = 220 L/h, RSE 13%; ESM2 $THETA(4). This is simultaneously the MHD formation clearance.
    lvc          <- log(33.1);      label("OXZ central volume of distribution at 50 kg PKWT (L)")                   # Sinha 2026 Table 4: V_OXZ = 33.1 L, RSE 34%; ESM2 $THETA(2). The 2026 Correction (doi:10.1007/s40262-025-01613-1) restates Eq. 10 as TV_V,OXZ = 33.1*(PKWT/50)^0.752; the original printed it as a duplicate of the V_MHD equation.

    # MHD metabolite parameters (suffix "_mhd", registered as a
    # metabolite-suffix entry in inst/references/compartment-names.md;
    # founding example Rodrigues_2017_oxcarbazepine.R).
    lcl_mhd      <- log(3.05);      label("MHD elimination clearance at 50 kg PKWT (L/h)")                          # Sinha 2026 Table 4: CL_MHD = 3.05 L/h, RSE 7%; ESM2 $THETA(5)
    lvc_mhd      <- fixed(log(50)); label("MHD central volume of distribution at 50 kg PKWT (L), literature value") # Sinha 2026 Table 4: V_MHD = 50 L (fixed); ESM2 $THETA(3) "50 FIX". Fixed to the healthy-adult NCA estimate of reference [5] because the joint parent-metabolite framework leaves it structurally unidentifiable (Results 3.2.1 model M5, Discussion).

    # Back-transformation of MHD to OXZ (metabolite-to-parent conversion).
    lkbt         <- log(0.0433);    label("MHD -> OXZ back-transformation rate constant (1/h)")                     # Sinha 2026 Table 4: Kbt = 0.0433 1/h, RSE 22%; ESM2 $THETA(9)

    # Allometric exponents on (PKWT / 50). The volume exponent is a SINGLE
    # shared parameter applied to both V_OXZ and V_MHD (ESM2 $PK: TVV2 and TVV3
    # both use THETA(6)).
    e_ffm_cl     <- fixed(1);       label("Allometric exponent on OXZ clearance (unitless)")                        # Sinha 2026 Table 4: theta_CL,OXZ = 1 (fixed); ESM2 $THETA(7) "1 FIX". Estimated at 0.98 in model M3b, then fixed to 1 (Results 3.2.1).
    e_ffm_cl_mhd <- 0.671;          label("Allometric exponent on MHD clearance (unitless)")                        # Sinha 2026 Table 4: theta_CL,MHD = 0.671, RSE 10%; ESM2 $THETA(8)
    e_ffm_vc     <- 0.752;          label("Shared allometric exponent on the OXZ and MHD volumes (unitless)")       # Sinha 2026 Table 4: theta_V = 0.752, RSE 10%; ESM2 $THETA(6), shared by TVV2 and TVV3

    # Between-subject variability. Sinha 2026 Methods 2.2 defines BSV as "the
    # standard deviation of log-normally distributed random effects with a zero
    # mean", and Table 4 reports those standard deviations as percentages. The
    # reported percentage is therefore taken as omega itself, and ini() needs
    # omega^2. See the vignette "Assumptions and deviations" for the full
    # argument, including the RSE-based check that the reported quantities are
    # on the SD (not variance) scale.
    #
    # BSV was fixed to zero for Ka, Kbt, and V_OXZ (Table 4 "0 (fixed)";
    # ESM2 $OMEGA ETA(1), ETA(2), ETA(9) all "0 FIX"), so those parameters
    # carry no eta at all.
    etalcl       ~ 0.2601                                                                                          # 0.51^2; Sinha 2026 Table 4: BSV on CL_OXZ = 51% (RSE 19%, shrinkage 26%); ESM2 $OMEGA ETA(4)
    etalcl_mhd   ~ 0.1296                                                                                          # 0.36^2; Sinha 2026 Table 4: BSV on CL_MHD = 36% (RSE 12%, shrinkage 16%); ESM2 $OMEGA ETA(5)
    etalvc_mhd   ~ 0.64                                                                                            # 0.80^2; Sinha 2026 Table 4: BSV on V_MHD = 80% (RSE 20%, shrinkage 33%); ESM2 $OMEGA ETA(3)

    # Residual unexplained variability: separate proportional errors for the two
    # analytes, no additive component (Results 3.2.1 model M2 selection). For a
    # NONMEM proportional error Y = F*(1+EPS) the reported CV% is exactly
    # sqrt(sigma^2)*100, so these are direct SD fractions.
    propSd       <- 0.52;           label("OXZ proportional residual SD (fraction)")                                # Sinha 2026 Table 4: sigma_prop,OXZ = 52% CV (RSE 13%, shrinkage 14%); ESM1 Table S2 lists the same quantity as "Proportional error (CV%)"
    propSd_mhd   <- 0.24;           label("MHD proportional residual SD (fraction)")                                # Sinha 2026 Table 4: sigma_prop,MHD = 24% CV (RSE 28%, shrinkage 20%)
  })

  model({
    # Reference body size for allometric scaling: 50 kg PKWT (ESM2 $PK, every
    # scaling term is (PKWT/50)^theta).
    ref_ffm <- 50

    # Salt / molar-mass factors for the bidirectional parent <-> metabolite
    # transfer. Sinha 2026 Methods 2.3.1: "Both the forward and backward mass
    # flow between the parent and the metabolite were based on their molar
    # masses and were implemented by accounting for the salt factors SF_ft and
    # SF_bt of the respective steps. The salt factors were calculated from the
    # molecular weight of OXZ (252.27 Da) and MHD (254.28 Da)." The ESM2 control
    # stream hardcodes these as 1.00796 (MHD/OXZ, forward) and 0.992095
    # (OXZ/MHD, backward); the ratio form below is identical to 6 significant
    # figures and is exactly reciprocal, which the compartment amounts require.
    sf_ft <- 254.28 / 252.27
    sf_bt <- 252.27 / 254.28

    # Individual parameters. Ka, V_OXZ, and Kbt carry no BSV (fixed to zero in
    # the source model).
    ka     <- exp(lka)
    cl     <- exp(lcl + etalcl)             * (FFM / ref_ffm)^e_ffm_cl
    vc     <- exp(lvc)                      * (FFM / ref_ffm)^e_ffm_vc
    cl_mhd <- exp(lcl_mhd + etalcl_mhd)     * (FFM / ref_ffm)^e_ffm_cl_mhd
    vc_mhd <- exp(lvc_mhd + etalvc_mhd)     * (FFM / ref_ffm)^e_ffm_vc
    kbt    <- exp(lkbt)

    # Micro-constants. OXZ is assumed to be eliminated EXCLUSIVELY by conversion
    # to MHD (Methods 2.3.1), so kel is simultaneously the OXZ elimination rate
    # constant and the MHD formation rate constant.
    kel  <- cl     / vc
    kelm <- cl_mhd / vc_mhd

    # ODE system - ESM2 $DES (final model M6), Fig. 1 Eqs. 1-3.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central + kbt * central_mhd * sf_bt
    d/dt(central_mhd) <-  kel * central * sf_ft - kelm * central_mhd - kbt * central_mhd

    # Observations.
    Cc     <- central     / vc
    Cc_mhd <- central_mhd / vc_mhd

    Cc     ~ prop(propSd)
    Cc_mhd ~ prop(propSd_mhd)
  })
}
