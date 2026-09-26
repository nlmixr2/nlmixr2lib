Taddio_2018_am7_mouse_pbpk <- function() {
  description <- paste(
    "Preclinical (mouse). PBPK (semi-physiological whole-body), fitted to dynamic",
    "whole-body PET. Disposition of the investigational carbon-11 CD80 PET tracer",
    "[11C]AM7 in a xenograft-bearing female CD1 nude mouse under control",
    "conditions. Nine states carry decay-corrected radioactivity as an AMOUNT:",
    "blood plasma (central), a reversible and an irreversible hepatic pool",
    "(liver_exchange, liver_deep), the combined gallbladder plus intestinal lumen",
    "(gallbladder_intestine), two renal pools (kidney_exchange, kidney_deep), two",
    "parallel peripheral-tissue pools (peripheral1, peripheral2) and urine. Unlike",
    "its structurally modified successor [11C]MT107, [11C]AM7 carries a fitted",
    "transintestinal plasma-to-intestine route kbg, which contributes 18.3 percent",
    "of the total hepatobiliary clearance reported in Table 1, and is cleared",
    "renally at essentially the glomerular filtration rate (CLR/GFR 1.0),",
    "consistent with its low lipophilicity (logD 0.1) and correspondingly little",
    "tubular reabsorption. The urinary bladder does not empty and there is no",
    "defecation during the 60-minute scan, so gallbladder_intestine and urine are",
    "terminal accumulating sinks and the system conserves mass exactly. The model",
    "also returns the paper's derived hepatobiliary and renal clearances, hepatic",
    "extraction ratio and tissue distribution coefficient (equations 2 to 4).",
    "Parameters are those of the single reference scan of Figure 6(a), the only scan",
    "in the study for which the whole urinary bladder was inside the field of view",
    "and dissection data were available, so that the modelled urinary curve could be",
    "validated against data that were not fitted. The paper reports no hierarchical",
    "random effects, so the model is deterministic. Companion models:",
    "Taddio_2018_mt107_mouse_pbpk and Taddio_2018_mt107_cyclosporine_mouse_pbpk.",
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
    "'AM7 Control / Fig. 6A'. Fitted and experimental time-activity curves for this",
    "scan, including the urinary curve that was predicted rather than fitted: Figure",
    "6(a). Derived clearance, extraction-ratio and tissue-distribution equations:",
    "Materials and Methods section 2.3, equations 2, 3 and 4. Summary clearances and",
    "the transintestinal contribution: Table 1 and its footnote c. Organ volumes,",
    "blood volume, hematocrit, plasma flows and glomerular filtration rate:",
    "Materials and Methods sections 2.2 and 2.3, and the footer of Supplementary",
    "Table 1.",
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
    central = list(analyte = "[11C]AM7 radioactivity", units = "MBq", specimen = "plasma", verified = TRUE),
    liver_exchange = list(analyte = "[11C]AM7 radioactivity", units = "MBq", specimen = "tissue", verified = TRUE),
    liver_deep = list(analyte = "[11C]AM7 radioactivity", units = "MBq", specimen = "tissue", verified = TRUE),
    gallbladder_intestine = list(analyte = "[11C]AM7 radioactivity", units = "MBq", specimen = "bile", verified = TRUE),
    kidney_exchange = list(analyte = "[11C]AM7 radioactivity", units = "MBq", specimen = "tissue", verified = TRUE),
    kidney_deep = list(analyte = "[11C]AM7 radioactivity", units = "MBq", specimen = "tissue", verified = TRUE),
    peripheral1 = list(analyte = "[11C]AM7 radioactivity", units = "MBq", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "[11C]AM7 radioactivity", units = "MBq", specimen = "tissue", verified = TRUE),
    urine = list(analyte = "[11C]AM7 radioactivity", units = "MBq", specimen = "urine", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight of the mouse. Converts the paper's per-gram organ and plasma volumes into the absolute volumes that the concentration observables and the derived clearances need.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Does not enter any rate constant. The mass-transfer rate constants are first-order (1/min) and were fitted per scan on the amount scale, so body weight affects only the volumes used to turn amounts into concentrations and rate constants into clearances.",
        "Supplementary Table 1 gives 24.4 g (0.0244 kg) for this scan; the second [11C]AM7 control scan, Figure 6(b), was a 19.7 g animal. The register carries WT in kg, so the model multiplies by 1000 internally to reach the paper's per-gram basis.",
        sep = " "
      ),
      source_name = "BW"
    )
  )

  population <- list(
    species = "mouse (female CD1 nude, carrying hCD80-positive Raji xenografts)",
    n_subjects = 1L,
    n_studies = 1L,
    age_range = "7-10 weeks",
    weight_range = "24.4 g (this scan); 19.7 g for the second control scan",
    sex_female_pct = 100,
    disease_state = "hCD80-positive Raji xenograft-bearing immunodeficient mice; the xenograft radioactivity fraction was negligible and was not modelled",
    dose_range = "13.9 MBq [11C]AM7 (< 20 nmol/kg) as a ~10 s intravenous injection in 100-200 uL saline with 5% ethanol",
    regions = "Switzerland (ETH Zurich)",
    co_medication = "vehicle only (13% ethanol, 2 mL/kg)",
    notes = paste(
      "Two [11C]AM7 control scans were modelled. This file carries the Figure 6(a) scan, which the authors used to evaluate the whole modelling approach because it was the only scan that both contained the complete urinary bladder in the field of view and had terminal dissection data; the modelled urinary curve was therefore predicted rather than fitted and could be checked against the image-derived urine curve. The second control scan (Figure 6(b), 19.7 g) is reproduced in the vignette by overriding ini().",
      "The single [11C]AM7 scan after cyclosporine (Supplementary Figure 3, 25.3 g) is deliberately NOT extracted: the authors state that its kidney curve was poorly defined, that the data 'were not suitable for modelling' and that the fitting results are unreliable, and Table 1 records every modelled [11C]AM7 cyclosporine quantity as not determined. Its rate constants nonetheless appear in Supplementary Table 1 and are listed in the vignette Errata.",
      "Data were repurposed from PET experiments that had not been designed for PBPK modelling. The scan was started 60 s after tracer injection and lasted 60 min, so the first minute after injection is missing. Compared with dissection at 68 min, the image-derived liver and combined gallbladder-plus-intestine radioactivities were underestimated and the peripheral-tissue radioactivity (shoulder volume of interest) was higher than the dissected muscle.",
      "The infusion duration used in the fits was T = 0.15 min (Supplementary Table 1), slightly shorter than the ~10 s stated in Materials and Methods section 2.1.",
      sep = " "
    )
  )

  ini({
    # ---- Hepatobiliary mass-transfer rate constants (1/min) --------------
    # Suppl. Table 1, row 'AM7 Control / Fig. 6A'. Unlike the [11C]MT107 scans,
    # kbg was fitted rather than set to zero for [11C]AM7.
    kbh1 <- 0.2595; label("Plasma to reversible hepatic pool mass-transfer rate constant (1/min)") # Suppl. Table 1 'kBH1' 0.2595 (Fig. 6B scan: 0.5099)
    kh1b <- 0.2306; label("Reversible hepatic pool to plasma mass-transfer rate constant (1/min)") # Suppl. Table 1 'KH1B' 0.2306 (Fig. 6B scan: 0.3992)
    kh1h2 <- 0.1569; label("Reversible to irreversible hepatic pool mass-transfer rate constant (1/min)") # Suppl. Table 1 'kH1H2' 0.1569 (Fig. 6B scan: 0.208)
    kh2g <- 0.0193; label("Irreversible hepatic pool to gallbladder and intestine mass-transfer rate constant (1/min)") # Suppl. Table 1 'kH2G' 0.0193 (Fig. 6B scan: 0.0145)
    kgh1 <- 0.0397; label("Intestinal reabsorption to reversible hepatic pool mass-transfer rate constant (1/min)") # Suppl. Table 1 'kGH1' 0.0397 (Fig. 6B scan: 0.0226)
    kbg <- 0.0236; label("Transintestinal plasma to intestine mass-transfer rate constant (1/min)") # Suppl. Table 1 'kBG' 0.0236 (Fig. 6B scan: 0.0107). Table 1 footnote c: contributes 18.3% of the reported CLH for this scan

    # ---- Renal mass-transfer rate constants (1/min) ----------------------
    kbr1 <- 0.3819; label("Plasma to first renal pool mass-transfer rate constant (1/min)") # Suppl. Table 1 'kBR1' 0.3819 (Fig. 6B scan: 0.4691)
    kr1b <- 0.2497; label("First renal pool to plasma mass-transfer rate constant (1/min)") # Suppl. Table 1 'kR1B' 0.2497 (Fig. 6B scan: 0.2939)
    kr1r2 <- fixed(0); label("First to second renal pool mass-transfer rate constant (1/min)") # Suppl. Table 1 'kR1R2' 0 for this scan; Methods 2.3 set a lower bound of 0 on every k, so the fit returned the bound and the second renal pool stays empty (Fig. 6B scan: 0.0035)
    kr2r1 <- 0.2676; label("Second to first renal pool mass-transfer rate constant (1/min)") # Suppl. Table 1 'kR2R1' 0.2676 (Fig. 6B scan: 0.0655). Unidentifiable while kr1r2 is 0
    kr1u <- 0.2796; label("First renal pool to urine mass-transfer rate constant (1/min)") # Suppl. Table 1 'kR1U' 0.2796 (Fig. 6B scan: 0.4241)

    # ---- Peripheral-tissue mass transfer ---------------------------------
    # Figure 2: plasma splits into the fast pool at fbt1 * kbt and into the slow
    # pool at (1 - fbt1) * kbt, so kbt is the total plasma-to-tissue rate
    # constant. Suppl. Table 1 tabulates '1-fBT1' rather than fbt1 itself.
    kbt <- 0.8149; label("Total plasma to peripheral tissue mass-transfer rate constant (1/min)") # Suppl. Table 1 'kBT' 0.8149 (Fig. 6B scan: 0.4235)
    fbt1 <- 0.6962; label("Fraction of the plasma to peripheral tissue rate constant directed to the fast pool (unitless)") # Suppl. Table 1 '1-fBT1' 0.3038, so fbt1 = 1 - 0.3038 (Fig. 6B scan: 1-fBT1 = 1, i.e. fbt1 = 0)
    kt1b <- 163.382; label("Fast peripheral pool to plasma mass-transfer rate constant (1/min)") # Suppl. Table 1 'kT1B' 163.382 (Fig. 6B scan: 110.3334). Very fast, so this pool is effectively in instantaneous equilibrium with plasma
    kt2b <- 0.1222; label("Slow peripheral pool to plasma mass-transfer rate constant (1/min)") # Suppl. Table 1 'kT2B' 0.1222 (Fig. 6B scan: 0.1948)

    # ---- Volumes and physiological constants -----------------------------
    vtissue_bw <- 0.6001; label("Peripheral tissue volume per unit body weight (mL/g)") # Suppl. Table 1 'VTissue/BW' 0.6001 (Fig. 6B scan: 0.6). Estimated during fitting
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
    # Table 1 reports CLH for [11C]AM7 as the sum of the equation-2 clearance and
    # the transintestinal clearance kbg * vplasma (footnote c), so cltot below
    # matches the Table 1 'CL total' row.
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
