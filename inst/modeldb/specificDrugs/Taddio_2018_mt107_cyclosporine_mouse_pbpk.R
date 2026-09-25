Taddio_2018_mt107_cyclosporine_mouse_pbpk <- function() {
  description <- paste(
    "Preclinical (mouse). PBPK (semi-physiological whole-body), fitted to dynamic",
    "whole-body PET. Disposition of the investigational carbon-11 CD80 PET tracer",
    "[11C]MT107 in xenograft-bearing female SCID mice that received cyclosporine",
    "50 mg/kg intravenously 30 to 50 minutes before the tracer, to inhibit",
    "cyclosporine-sensitive hepatic and renal transporters. Structurally identical",
    "to the control model Taddio_2018_mt107_mouse_pbpk: nine states carry",
    "decay-corrected radioactivity as an AMOUNT in blood plasma (central), a",
    "reversible and an irreversible hepatic pool (liver_exchange, liver_deep), the",
    "combined gallbladder plus intestinal lumen (gallbladder_intestine), two renal",
    "pools (kidney_exchange, kidney_deep), two parallel peripheral-tissue pools",
    "(peripheral1, peripheral2) and urine. Cyclosporine is not a covariate here: the",
    "authors fitted each scan independently and estimated no drug-interaction",
    "coefficient, so the treatment effect is carried entirely by this file's",
    "parameter values. Relative to the control model, hepatobiliary clearance falls",
    "from about 35 to about 17 uL/min and the hepatic extraction ratio from 0.035 to",
    "0.017 (both significant), renal clearance rises about 1.7-fold, and the",
    "excretion pattern shifts from predominantly hepatobiliary (CLH/CLR 3.8) to",
    "balanced (CLH/CLR 0.9). The transintestinal route kbg is set to zero, as the",
    "authors did for [11C]MT107. Parameters are the mean of four independently",
    "fitted cyclosporine scans; the paper reports no hierarchical random effects, so",
    "the model is deterministic.",
    sep = " "
  )
  reference <- paste(
    "Taddio MF, Mu L, Keller C, Schibli R, Kramer SD. Physiologically Based",
    "Pharmacokinetic Modelling with Dynamic PET Data to Study the In Vivo Effects of",
    "Transporter Inhibition on Hepatobiliary Clearance in Mice.",
    "Contrast Media Mol Imaging. 2018;2018:5849047. doi:10.1155/2018/5849047.",
    "PMCID: PMC6008768.",
    "Model structure and compartment labels: Figure 2. Mass-transfer rate constants,",
    "peripheral-tissue volume and infusion duration: Supplementary Table 1, row",
    "'MT107 Ciclosporin / Average' (n = 4). Derived clearance, extraction-ratio and",
    "tissue-distribution equations: Materials and Methods section 2.3, equations 2,",
    "3 and 4. Group comparison and significance tests: Table 1 and Results section",
    "3.3. Organ volumes, blood volume, hematocrit, plasma flows and glomerular",
    "filtration rate: Materials and Methods sections 2.2 and 2.3, and the footer of",
    "Supplementary Table 1.",
    sep = " "
  )
  vignette <- "Taddio_2018_hepatobiliary_transport_pet_pbpk"
  units <- list(time = "min", dosing = "MBq", concentration = "MBq/mL")

  # See Taddio_2018_mt107_mouse_pbpk for the rationale: the two-pool-per-organ
  # parameterisation describes the PET time-activity curves rather than anatomy,
  # and the authors report that a parallel arrangement fitted equally well.
  paper_specific_compartments <- c(
    "liver_exchange",
    "liver_deep",
    "gallbladder_intestine",
    "kidney_exchange",
    "kidney_deep"
  )

  compartmentData <- list(
    central = list(analyte = "[11C]MT107 radioactivity", units = "MBq", specimen = "plasma", verified = TRUE),
    liver_exchange = list(analyte = "[11C]MT107 radioactivity", units = "MBq", specimen = "tissue", verified = TRUE),
    liver_deep = list(analyte = "[11C]MT107 radioactivity", units = "MBq", specimen = "tissue", verified = TRUE),
    gallbladder_intestine = list(
      analyte = "[11C]MT107 radioactivity",
      units = "MBq",
      specimen = "bile",
      verified = TRUE
    ),
    kidney_exchange = list(analyte = "[11C]MT107 radioactivity", units = "MBq", specimen = "tissue", verified = TRUE),
    kidney_deep = list(analyte = "[11C]MT107 radioactivity", units = "MBq", specimen = "tissue", verified = TRUE),
    peripheral1 = list(analyte = "[11C]MT107 radioactivity", units = "MBq", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "[11C]MT107 radioactivity", units = "MBq", specimen = "tissue", verified = TRUE),
    urine = list(analyte = "[11C]MT107 radioactivity", units = "MBq", specimen = "urine", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight of the mouse. Converts the paper's per-gram organ and plasma volumes into the absolute volumes that the concentration observables and the derived clearances need.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Does not enter any rate constant. The mass-transfer rate constants are first-order (1/min) and were fitted per scan on the amount scale, so body weight affects only the volumes used to turn amounts into concentrations and rate constants into clearances.",
        "Supplementary Table 1 reports body weight per scan; the four cyclosporine [11C]MT107 scans were 19.3, 21.2, 21.0 and 19.1 g (mean 20.15 g, i.e. 0.02015 kg). The register carries WT in kg, so the model multiplies by 1000 internally to reach the paper's per-gram basis.",
        sep = " "
      ),
      source_name = "BW"
    )
  )

  population <- list(
    species = "mouse (female C.B.17 SCID, carrying hCD80-positive Raji xenografts)",
    n_subjects = 4L,
    n_studies = 1L,
    age_range = "7-10 weeks",
    weight_range = "19.1-21.2 g",
    weight_median = "20.15 g (mean of the four cyclosporine scans)",
    sex_female_pct = 100,
    disease_state = "hCD80-positive Raji xenograft-bearing immunodeficient mice; the xenograft radioactivity fraction was negligible and was not modelled",
    dose_range = "3-14 MBq [11C]MT107 (< 20 nmol/kg) as a ~10 s intravenous injection in 100-200 uL saline with 5% ethanol",
    regions = "Switzerland (ETH Zurich)",
    co_medication = "cyclosporine 50 mg/kg intravenously 30-50 min before tracer, as a 1:1 dilution of Sandimmun (50 mg/mL cyclosporine in 26.1% ethanol / 65% PEG-35 castor oil) in water for injection",
    notes = paste(
      "Data were repurposed from PET experiments that had not been designed for PBPK modelling. Scans were acquired on a SuperArgus PET/CT, started 60 s after tracer injection and lasting 60 min, so the first minute after injection is missing from every data set.",
      "Each scan was fitted independently with a custom MATLAB script (ode45, fmincon, MultiStart with 128 random initial-parameter sets); the parameters in this model are the arithmetic mean of the four cyclosporine fits reported in Supplementary Table 1. Because the derived clearances are non-linear functions of the rate constants, evaluating equation 2 at these mean rate constants gives a hepatobiliary clearance of about 17.9 uL/min rather than the 17.1 uL/min mean of the four per-scan values reported in Table 1; the vignette reproduces the per-scan values exactly.",
      "The control experiments cannot exclude an effect of the Sandimmun adjuvants (ethanol, PEG-35 castor oil) on tracer pharmacokinetics; the authors list this as limitation (i) in the Discussion. Cyclosporine can also displace drugs from plasma protein binding, which would raise the unbound fraction and could explain part of the non-significant rises in renal clearance and tissue distribution.",
      "The between-animal standard deviations of the fitted rate constants are reported in Supplementary Table 1 but are dispersions of four independent least-squares fits, not estimated random effects, so they are not encoded as IIV.",
      sep = " "
    )
  )

  ini({
    # ---- Hepatobiliary mass-transfer rate constants (1/min) --------------
    # Suppl. Table 1, row 'MT107 Ciclosporin / Average' (n = 4). The trailing SD
    # is the between-animal dispersion of the four independent fits.
    kbh1 <- 1.1767; label("Plasma to reversible hepatic pool mass-transfer rate constant (1/min)") # Suppl. Table 1 'kBH1' 1.1767, SD 0.3832
    kh1b <- 2.4244; label("Reversible hepatic pool to plasma mass-transfer rate constant (1/min)") # Suppl. Table 1 'KH1B' 2.4244, SD 0.6461; t-test vs control p = 0.0109
    kh1h2 <- 0.0573; label("Reversible to irreversible hepatic pool mass-transfer rate constant (1/min)") # Suppl. Table 1 'kH1H2' 0.0573, SD 0.022
    kh2g <- 0.1576; label("Irreversible hepatic pool to gallbladder and intestine mass-transfer rate constant (1/min)") # Suppl. Table 1 'kH2G' 0.1576, SD 0.0706
    kgh1 <- 0.008; label("Intestinal reabsorption to reversible hepatic pool mass-transfer rate constant (1/min)") # Suppl. Table 1 'kGH1' 0.008, SD 0.0039; Results 3.3 reports 0.0080 +/- 0.0038 versus 0.0077 +/- 0.0033 in control, p = 0.90
    kbg <- fixed(0); label("Transintestinal plasma to intestine mass-transfer rate constant (1/min)") # Results 3.3: 'kBG ... was set to 0' for the [11C]MT107 scans; Figure 2 grey arrow

    # ---- Renal mass-transfer rate constants (1/min) ----------------------
    kbr1 <- 0.1576; label("Plasma to first renal pool mass-transfer rate constant (1/min)") # Suppl. Table 1 'kBR1' 0.1576, SD 0.0495
    kr1b <- 0.9242; label("First renal pool to plasma mass-transfer rate constant (1/min)") # Suppl. Table 1 'kR1B' 0.9242, SD 0.5265
    kr1r2 <- 0.0801; label("First to second renal pool mass-transfer rate constant (1/min)") # Suppl. Table 1 'kR1R2' 0.0801, SD 0.051
    kr2r1 <- 0.0479; label("Second to first renal pool mass-transfer rate constant (1/min)") # Suppl. Table 1 'kR2R1' 0.0479, SD 0.0356
    kr1u <- 0.2224; label("First renal pool to urine mass-transfer rate constant (1/min)") # Suppl. Table 1 'kR1U' 0.2224, SD 0.13

    # ---- Peripheral-tissue mass transfer ---------------------------------
    # Figure 2: plasma splits into the fast pool at fbt1 * kbt and into the slow
    # pool at (1 - fbt1) * kbt, so kbt is the total plasma-to-tissue rate
    # constant. Suppl. Table 1 tabulates '1-fBT1' rather than fbt1 itself.
    kbt <- 1.2932; label("Total plasma to peripheral tissue mass-transfer rate constant (1/min)") # Suppl. Table 1 'kBT' 1.2932, SD 0.8487
    fbt1 <- 0.8834; label("Fraction of the plasma to peripheral tissue rate constant directed to the fast pool (unitless)") # Suppl. Table 1 '1-fBT1' 0.1166, SD 0.0898, so fbt1 = 1 - 0.1166
    kt1b <- 2.8165; label("Fast peripheral pool to plasma mass-transfer rate constant (1/min)") # Suppl. Table 1 'kT1B' 2.8165, SD 2.4081
    kt2b <- 0.0335; label("Slow peripheral pool to plasma mass-transfer rate constant (1/min)") # Suppl. Table 1 'kT2B' 0.0335, SD 0.0043

    # ---- Volumes and physiological constants -----------------------------
    vtissue_bw <- 0.7919; label("Peripheral tissue volume per unit body weight (mL/g)") # Suppl. Table 1 'VTissue/BW' 0.7919, SD 0.0955; Methods 2.2 quote 0.79 +/- 0.10 cm3/g, not significantly different from control (p = 0.47). Estimated during fitting
    vplasma_bw <- fixed(0.03276); label("Plasma volume per unit body weight (mL/g)") # Suppl. Table 1 footer 'vPlasma 0.03276 ml/g'; Methods 2.3: VBlood * (1 - hematocrit) = 0.0585 * 0.56
    vliver_bw <- fixed(0.065); label("Liver volume per unit body weight (mL/g)") # Methods 2.2: '0.065 cm3 per g BW for liver', Davies and Morris
    vkidney_bw <- fixed(0.0164); label("Kidney volume per unit body weight (mL/g)") # Methods 2.2: '0.0164 cm3 per g BW for kidneys', Davies and Morris
    vblood_bw <- fixed(0.0585); label("Blood volume per unit body weight (mL/g)") # Methods 2.2: 'V Blood 0.0585 ml per g body weight'
    hct <- fixed(0.44); label("Hematocrit (unitless)") # Methods 2.2: 'The hematocrit was assumed 0.44'
    qh <- fixed(1000); label("Hepatic plasma flow (uL/min)") # Suppl. Table 1 footer 'QH 1000 ul/min'; Methods 2.3 QP,H 1.0 mL/min
    qr <- fixed(730); label("Renal plasma flow (uL/min)") # Suppl. Table 1 footer 'QR 730 ul/min'; Methods 2.3 QP,R 0.73 mL/min
    gfr <- fixed(160); label("Glomerular filtration rate (uL/min)") # Suppl. Table 1 footer 'GFR 160 ul/min'; Results 3.2 'maximal expected CL by glomerular filtration would be ~160 ul/min'

    # ---- Residual error ---------------------------------------------------
    propSd <- fixed(0); label("Proportional residual SD on the plasma radioactivity concentration (fraction); magnitude not reported") # Methods 2.3: the objective was a weighted sum of squared residuals with the first two plasma, liver and kidney residuals weighted 5-fold; no residual-error model is reported
  })

  model({
    # ---- Absolute volumes for this animal --------------------------------
    bwg <- WT * 1000
    vplasma <- vplasma_bw * bwg
    vliver <- vliver_bw * bwg
    vkidney <- vkidney_bw * bwg
    vtissue <- vtissue_bw * bwg

    # ---- Compartment system (Figure 2) -----------------------------------
    d/dt(central) <- kh1b * liver_exchange - kbh1 * central +
      kr1b * kidney_exchange - kbr1 * central +
      kt1b * peripheral1 + kt2b * peripheral2 - kbt * central -
      kbg * central
    d/dt(liver_exchange) <- kbh1 * central - kh1b * liver_exchange -
      kh1h2 * liver_exchange + kgh1 * gallbladder_intestine
    d/dt(liver_deep) <- kh1h2 * liver_exchange - kh2g * liver_deep
    d/dt(gallbladder_intestine) <- kh2g * liver_deep + kbg * central -
      kgh1 * gallbladder_intestine
    d/dt(kidney_exchange) <- kbr1 * central - kr1b * kidney_exchange -
      kr1r2 * kidney_exchange + kr2r1 * kidney_deep -
      kr1u * kidney_exchange
    d/dt(kidney_deep) <- kr1r2 * kidney_exchange - kr2r1 * kidney_deep
    d/dt(peripheral1) <- fbt1 * kbt * central - kt1b * peripheral1
    d/dt(peripheral2) <- (1 - fbt1) * kbt * central - kt2b * peripheral2
    d/dt(urine) <- kr1u * kidney_exchange

    # ---- Observables -----------------------------------------------------
    Cc <- central / vplasma
    Cblood <- Cc * (1 - hct)
    Aliver <- liver_exchange + liver_deep + vblood_bw * vliver * Cblood
    Akidney <- kidney_exchange + kidney_deep + vblood_bw * vkidney * Cblood
    Atissue <- peripheral1 + peripheral2 + vblood_bw * vtissue * Cblood
    Agallbladder <- gallbladder_intestine
    Aurine <- urine

    # ---- Derived pharmacokinetic parameters (equations 2 to 4) -----------
    clh <- kbh1 * kh1h2 / (kh1b + kh1h2) * vplasma * 1000
    clbg <- kbg * vplasma * 1000
    clr <- kbr1 * kr1u / (kr1b + kr1u) * vplasma * 1000
    cltot <- clh + clbg + clr
    eh <- clh / qh
    clr_gfr <- clr / gfr
    # Results 3.2 frames CLR against the maximal renal clearance achievable with
    # transporter-mediated secretion on top of filtration, i.e. renal plasma flow.
    clr_qr <- clr / qr
    dtissue <- (fbt1 * kbt / kt1b + (1 - fbt1) * kbt / kt2b) * vplasma / vtissue

    Cc ~ prop(propSd)
  })
}
