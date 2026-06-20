Gupta_2015_ixazomib <- function() {
  description <- "Three-compartment population pharmacokinetic model for the oral proteasome inhibitor ixazomib (MLN9708) developed from pooled data of 226 adult patients with advanced multiple myeloma, lymphoma, or solid tumours across four phase I dose-escalation studies (Gupta 2015). Combined intravenous and oral data are described by a three-compartment model with first-order absorption and linear elimination; IV and oral data share the same disposition kinetics. Inter-individual variability is estimated on clearance, central volume V2, the second peripheral volume V4, absorption rate constant Ka, and bioavailability F; IIV on Q3, V3, and Q4 was fixed to zero. Body surface area on V4 (reference 1.90 m^2, exponent 2.3) is the only retained covariate; weight, age, gender, race, creatinine clearance, ALT, AST, albumin, and bilirubin had no clinically relevant effect on ixazomib pharmacokinetics. Residual error is additive on log-transformed concentration (NONMEM Y = LOG(F) + EPS(1)) which maps to a proportional error in linear concentration space. This analysis supported the switch from BSA-based to fixed (4 mg) dosing in subsequent ixazomib clinical studies."
  reference <- paste(
    "Gupta N, Zhao Y, Hui A-M, Esseltine D-L, Venkatakrishnan K. (2015).",
    "Switching from body surface area-based to fixed dosing for the",
    "investigational proteasome inhibitor ixazomib: a population",
    "pharmacokinetic analysis.",
    "Br J Clin Pharmacol 79(5):789-800.",
    "doi:10.1111/bcp.12542.",
    sep = " "
  )
  vignette <- "Gupta_2015_ixazomib"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    BSA = list(
      description        = "Body surface area at baseline.",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used as a power covariate on the second peripheral volume V4 (canonical name `lvp2`) with reference BSA = 1.90 m^2 and exponent 2.3 (Gupta 2015 Table 3 BSA-on-V4 entry; Appendix 1 NONMEM control stream V0 = (BSAC/1.90)**THETA(9); V4 = (THETA(6)*V0)*EXP(ETA(6))). The paper does not state the BSA formula (DuBois / Mosteller / Haycock); recorded as unspecified. The covariate-effect coefficient is encoded as `e_bsa_vp2` consistent with the e_<cov>_<param> convention. BSA was retained on V4 in preference to body weight because BSA and weight gave the same drop in OFV (35 points each), they are collinear, and BSA matched the BSA-based dosing convention used in the source phase I studies (Gupta 2015 Results, BSA versus weight comparison).",
      source_name        = "BSAC"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 226L,
    n_studies      = 4L,
    age_range      = "23-86 years (median 62.0)",
    age_median     = "62.0 years",
    weight_range   = "35.5-132.3 kg (median 78.0)",
    weight_median  = "78.0 kg",
    bsa_range      = "1.3-2.6 m^2 (median 1.9)",
    bsa_median     = "1.9 m^2",
    sex_female_pct = 42.5,
    race_ethnicity = c(Caucasian = 87, Other = 13),
    disease_state  = "Adult patients with advanced haematological or non-haematological malignancies: advanced solid tumours (n = 88), relapsed/refractory lymphoma (n = 30), relapsed/refractory multiple myeloma (n = 53 + 55 = 108).",
    dose_range     = "Body-surface-area-based oral or intravenous dosing 0.125-3.95 mg/m^2 once or twice weekly. Per-study ranges (Gupta 2015 Table 1): C16001 0.125-2.34 mg/m^2 IV twice weekly; C16002 0.125-3.11 mg/m^2 IV weekly; C16003 0.24-2.23 mg/m^2 oral twice weekly; C16004 0.24-3.95 mg/m^2 oral weekly. Oral ixazomib capsule strengths were 0.2, 0.5, and 2 mg; BSA-based doses were rounded to the available capsule strengths.",
    regions        = "Not stated (four pooled phase I oncology trials).",
    studies        = "C16001 (n = 88), C16002 (n = 30), C16003 (n = 53), C16004 (n = 55) (Gupta 2015 Table 1).",
    baseline_labs  = "Albumin 38.0 g/L (23-48); ALT 20.0 U/L (7-100); AST 24.0 U/L (9-82); bilirubin 7.0 uM (1.7-39.3); creatinine clearance 88.0 mL/min (21.9-213.7, Cockcroft-Gault). Median values with ranges per Gupta 2015 Table 2.",
    notes          = "Demographics from Gupta 2015 Table 2. The four pooled phase I studies all used BSA-based dosing; this analysis showed that BSA does not influence ixazomib clearance and supported the transition to fixed (4 mg) dosing in subsequent phase III studies (e.g., TOURMALINE-MM1 NCT01564537, NCT01850524). A later popPK analysis (Gupta 2017, modellib('Gupta_2017_ixazomib')) refit the model on a larger combined phase I/II/III dataset (n = 755) including data from TOURMALINE-MM1; that analysis retained BSA on V4 with similar exponent magnitude but updated reference (1.87 m^2 vs 1.90 m^2 here)."
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight; screened against CL, V2, and V4 in covariate analysis. WT was significant only on V4 (same OFV drop as BSA) and was rejected in favour of BSA in the final model because BSA-based dosing was used in the source phase I studies (Gupta 2015 Results).",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened; not retained in the final model."
    ),
    AGE = list(
      description = "Subject age in years; screened against CL with no clinically relevant effect over the 23-86 year range (Gupta 2015 Figure 3B; Discussion).",
      units       = "years",
      type        = "continuous",
      notes       = "Screened; not retained."
    ),
    SEXF = list(
      description = "Female-sex indicator; screened against CL, V2, and V4 with no significant effect (Gupta 2015 Figure 3E; Results).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened; not retained. Source coded gender as M/F (96 female of 226)."
    ),
    RACE = list(
      description = "Race (Caucasian vs other); screened on CL with no significant effect (Gupta 2015 Figure 3F).",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened; not retained. 196 Caucasian (87%) vs 30 other (13%)."
    ),
    CRCL = list(
      description = "Creatinine clearance (Cockcroft-Gault) over 22-213.7 mL/min; screened on CL with no clinically relevant effect (Gupta 2015 Figure 3C; Discussion).",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened; not retained. Renal elimination is a minor clearance pathway for ixazomib."
    ),
    ALT = list(
      description = "Alanine aminotransferase; screened on CL with no clinically relevant effect (Gupta 2015 Results, narrative).",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened; not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase; screened on CL with no clinically relevant effect (Gupta 2015 Results, narrative).",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened; not retained."
    ),
    ALB = list(
      description = "Serum albumin; screened on CL with no clinically relevant effect (Gupta 2015 Results, narrative).",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened; not retained."
    ),
    BILI = list(
      description = "Total bilirubin; screened on CL with no clinically relevant effect (Gupta 2015 Results, narrative).",
      units       = "uM",
      type        = "continuous",
      notes       = "Screened; not retained."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural PK parameters - final model estimates from Gupta 2015
    # Table 3. NONMEM theta names (Appendix 1 control stream) map to
    # nlmixr2 canonical names as: V2 -> central volume `lvc`;
    # V3 -> first peripheral volume `lvp`; V4 -> second peripheral
    # volume `lvp2`; Q3 -> inter-compartmental clearance central <->
    # peripheral1 `lq`; Q4 -> inter-compartmental clearance central
    # <-> peripheral2 `lq2`. The Mean column in Table 3 is the
    # back-transformed (linear-scale) point estimate; log-transform
    # for use in `ini()`.
    # ------------------------------------------------------------------
    lka      <- log(0.5)   ; label("First-order absorption rate constant (1/h)")          # Gupta 2015 Table 3: Ka = 0.5 (RSE 7.4%)
    lcl      <- log(2.0)   ; label("Systemic clearance (L/h)")                              # Gupta 2015 Table 3: CL = 2.0 (RSE 4.8%)
    lvc      <- log(14.3)  ; label("Central volume of distribution V2 (L)")                 # Gupta 2015 Table 3: V2 = 14.3 (RSE 10.0%)
    lq       <- log(9.7)   ; label("Inter-compartmental clearance Q3 (L/h)")                # Gupta 2015 Table 3: Q3 = 9.7 (RSE 15.3%)
    lvp      <- log(412.0) ; label("Volume of first peripheral compartment V3 (L)")         # Gupta 2015 Table 3: V3 = 412.0 (RSE 7.0%)
    lq2      <- log(22.3)  ; label("Inter-compartmental clearance Q4 (L/h)")                # Gupta 2015 Table 3: Q4 = 22.3 (RSE 5.7%)
    lvp2     <- log(83.4)  ; label("Reference volume of second peripheral compartment V4 (L) at BSA = 1.90 m^2") # Gupta 2015 Table 3: V4 = 83.4 (RSE 17.7%)
    lfdepot  <- log(0.6)   ; label("Absolute oral bioavailability F (fraction)")            # Gupta 2015 Table 3: F = 0.6 (RSE 6.0%)

    # ------------------------------------------------------------------
    # Covariate effect: BSA on V4 in power form with reference BSA =
    # 1.90 m^2. NONMEM Appendix 1:
    #   V0 = (BSAC/1.90)**THETA(9)
    #   V4 = (THETA(6) * V0) * EXP(ETA(6))
    # i.e. V4_i = lvp2_typical * (BSA_i / 1.90)^e_bsa_vp2 * exp(etalvp2_i)
    # ------------------------------------------------------------------
    e_bsa_vp2 <- 2.3       ; label("Power exponent for BSA effect on V4 (unitless)")        # Gupta 2015 Table 3 BSA-on-V4 = 2.3 (RSE 18.9%); Appendix 1 THETA(9)

    # ------------------------------------------------------------------
    # Inter-individual variability. The paper uses an exponential
    # (log-normal) IIV model on CL, V2, V4, Ka, F. IIV on Q3, V3, Q4
    # was fixed to zero in NONMEM ($OMEGA 0 FIXED on ETA(3), ETA(4),
    # ETA(5)) and so is omitted here. Table 3 reports %BSV* as the CV%
    # on the linear scale; the corresponding internal log-scale
    # variance is omega^2 = log(1 + CV^2).
    # ------------------------------------------------------------------
    etalcl     ~ 0.1647    # Gupta 2015 Table 3 %BSV CL 42.3% -> omega^2 = log(1 + 0.423^2) = 0.1647
    etalvc     ~ 0.6981    # Gupta 2015 Table 3 %BSV V2 100.5% -> omega^2 = log(1 + 1.005^2) = 0.6981
    etalvp2    ~ 0.1858    # Gupta 2015 Table 3 %BSV V4 45.2% -> omega^2 = log(1 + 0.452^2) = 0.1858
    etalka     ~ 0.2912    # Gupta 2015 Table 3 %BSV Ka 58.1% -> omega^2 = log(1 + 0.581^2) = 0.2912
    etalfdepot ~ 0.1660    # Gupta 2015 Table 3 %BSV F 42.5% -> omega^2 = log(1 + 0.425^2) = 0.1660

    # ------------------------------------------------------------------
    # Residual error: additive on log-transformed concentration
    # (NONMEM $ERROR Y = LOG(F) + ERR(1); $SIGMA 0.61 was the initial
    # value, final variance from Table 3 = 0.3). Per the NONMEM
    # default convention, the value 0.3 in Table 3 row "Additive
    # (Mean)" is interpreted as the variance of EPS(1) on the
    # log-concentration scale; the corresponding SD is sqrt(0.3)
    # which is passed to nlmixr2 `prop()` (additive-on-log maps to
    # proportional in linear space; see references/nonmem-translation.md
    # and the Assumptions section of the validation vignette for the
    # variance-vs-SD reading).
    # ------------------------------------------------------------------
    propSd     <- sqrt(0.3); label("Proportional residual error SD (additive on log-concentration; ~CV in linear space)") # Gupta 2015 Table 3 Additive variance = 0.3 (RSE 6.1%) -> SD = sqrt(0.3) ~ 0.5477
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual PK parameters. BSA scales V4 by (BSA / 1.90)^2.3
    # per Gupta 2015 Table 3 and Appendix 1. No other retained
    # covariate effects on any parameter in the final model. IIV is
    # applied as a multiplicative log-normal random effect on CL, V2,
    # V4 (lvp2), Ka, and F (lfdepot).
    # ------------------------------------------------------------------
    ka     <- exp(lka + etalka)
    cl     <- exp(lcl + etalcl)
    vc     <- exp(lvc + etalvc)
    q      <- exp(lq)
    vp     <- exp(lvp)
    q2     <- exp(lq2)
    vp2    <- exp(lvp2 + etalvp2) * (BSA / 1.90) ^ e_bsa_vp2
    fdepot <- exp(lfdepot + etalfdepot)

    # ------------------------------------------------------------------
    # 2. ODE system. NONMEM ADVAN12;TRANS4: first-order absorption
    # from depot -> central; central <-> peripheral1 (Q3) and central
    # <-> peripheral2 (Q4) distribution; first-order linear
    # elimination (CL) from central. IV doses land directly in
    # `central` via the cmt column of the event table; oral doses
    # land in `depot` and transit via ka.
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          (cl / vc) * central -
                          (q  / vc) * central + (q  / vp ) * peripheral1 -
                          (q2 / vc) * central + (q2 / vp2) * peripheral2
    d/dt(peripheral1) <-  (q  / vc) * central - (q  / vp ) * peripheral1
    d/dt(peripheral2) <-  (q2 / vc) * central - (q2 / vp2) * peripheral2

    # ------------------------------------------------------------------
    # 3. Bioavailability on the oral depot compartment (Gupta 2015
    # Appendix 1 F1 = THETA(8) * EXP(ETA(8))). IV doses bypass this.
    # ------------------------------------------------------------------
    f(depot) <- fdepot

    # ------------------------------------------------------------------
    # 4. Concentration in ng/mL. Dose units are mg, vc is in L, so
    # central/vc has units mg/L = 1 ug/mL = 1000 ng/mL. Multiply by
    # 1000 to express the observation in ng/mL, matching the assay
    # reporting units (Gupta 2015 Methods; LC-MS/MS LLOQ 0.5 ng/mL).
    # This matches the NONMEM scale S2 = V2/1000.
    # ------------------------------------------------------------------
    Cc <- (central / vc) * 1000

    # ------------------------------------------------------------------
    # 5. Proportional residual error. NONMEM additive-on-log-scale
    # error (Y = LOG(F) + ERR(1)) maps to nlmixr2 prop() in linear
    # concentration space (per references/nonmem-translation.md).
    # ------------------------------------------------------------------
    Cc ~ prop(propSd)
  })
}
