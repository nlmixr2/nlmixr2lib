# nlmixr2lib

# development version

* Add Novakovic 2017 cladribine ([doi:10.1208/s12248-016-9977-z](https://doi.org/10.1208/s12248-016-9977-z)) [DDMODEL00000223] — adults with relapsing-remitting multiple sclerosis (8-item EDSS Item Response Theory disease-progression model with FREM covariates).
* Add Schindler 2016 sunitinib ([doi:10.1002/psp4.12057](https://doi.org/10.1002/psp4.12057)) [DDMODEL00000221] — adults with imatinib-resistant or imatinib-intolerant advanced gastrointestinal stromal tumor.
* Add Jager 2011 gemtuzumab ([doi:10.1371/journal.pone.0024265](https://doi.org/10.1371/journal.pone.0024265)) [DDMODEL00000229] — patients with acute myeloid leukemia (mechanism-based PKPD with explicit CD33 binding and leukemic-blast cell killing; DDMORE Monolix re-fit with added peripheral PK compartment).
* Add Khan 2015 ciprofloxacin ([doi:10.1093/jac/dkv233](https://doi.org/10.1093/jac/dkv233)) [DDMODEL00000225] — in vitro time-kill experiments on E. coli K-12 wild-type and quinolone-resistant single-step mutants.
* Add Ter Heine 2014 tamoxifen ([doi:10.1111/bcp.12388](https://doi.org/10.1111/bcp.12388)) [DDMODEL00000212] — adults with breast cancer on chronic 20 mg PO QD tamoxifen at steady state (joint parent-metabolite popPK with CYP2D6 and CYP3A4/5 covariates on endoxifen formation).
* Add Mohamed 2016 colistin + meropenem ([doi:10.1093/jac/dkv488](https://doi.org/10.1093/jac/dkv488)) [DDMODEL00000173] — in vitro time-kill PK/PD against P. aeruginosa wild-type ATCC 27853 and meropenem-resistant ARU552.
* Add Bizzotto 2016 glucose ([doi:10.1152/ajpendo.00045.2016](https://doi.org/10.1152/ajpendo.00045.2016)) [DDMODEL00000227] — adults across the glucose-tolerance spectrum (mechanistic glucose-tracer kinetics simulator).
* Add Netterberg 2017 docetaxel ([doi:10.1007/s00280-017-3366-x](https://doi.org/10.1007/s00280-017-3366-x)) [DDMODEL00000224] — adults receiving docetaxel chemotherapy (Friberg-style myelosuppression PD model with Kloft 2006 parameter set).
* Add Zecchin 2016 survival ([doi:10.1111/bcp.12994](https://doi.org/10.1111/bcp.12994)) [DDMODEL00000218] — adults with advanced epithelial ovarian cancer.
* Add Hansson 2013 sunitinib fatigue Markov + proportional-odds model ([doi:10.1038/psp.2013.62](https://doi.org/10.1038/psp.2013.62)) [DDMODEL00000222] — adults with imatinib-resistant gastrointestinal stromal tumors.
* Add Svensson 2016 bedaquiline ([doi:10.1002/psp4.12147](https://doi.org/10.1002/psp4.12147)) [DDMODEL00000219] — adults with multidrug-resistant tuberculosis (parent + N-desmethyl M2 metabolite, time-varying weight and albumin).
* Add Jonsson 2011 ethambutol ([doi:10.1128/AAC.00274-11](https://doi.org/10.1128/AAC.00274-11)) [DDMODEL00000220] — adult South African pulmonary tuberculosis patients.
* Add Zecchin 2016 tumorovarian ([doi:10.1111/bcp.12994](https://doi.org/10.1111/bcp.12994)) [DDMODEL00000217] — adults with advanced epithelial ovarian cancer receiving carboplatin monotherapy or carboplatin + gemcitabine combination chemotherapy.
* Add Girard 2012 pimasertib ([www.page-meeting.org/?abstract=2458](https://www.page-meeting.org/?abstract=2458)) [DDMODEL00000215] — adults with advanced solid tumours and hematological malignancies (joint K-PD / cumulative-logit Markov ocular-AE-grade and Weibull-TTE dropout model).
* Add Lestini 2015 TGF-β inhibitor ([doi:10.1007/s11095-015-1693-3](https://doi.org/10.1007/s11095-015-1693-3)) [DDMODEL00000192] — simulated 50-subject oncology cohort (one-compartment first-order absorption PK + indirect-response biomarker turnover; first nlmixr2lib `inst/modeldb/ddmore/` entry).
* Add Li 2006 meropenem ([doi:10.1177/0091270006291035](https://doi.org/10.1177/0091270006291035)) [DDMODEL00000213] — adult patients receiving meropenem.
* Add Friberg 2002 paclitaxel ([doi:10.1200/JCO.2002.02.140](https://doi.org/10.1200/JCO.2002.02.140)) [DDMODEL00000186] — adult cancer patients receiving paclitaxel chemotherapy (semi-mechanistic myelosuppression PK/PD with leukocyte output).
* Add Hansson 2013 sunitinib ([doi:10.1038/psp.2013.61](https://doi.org/10.1038/psp.2013.61)) [DDMODEL00000197] — adults with imatinib-resistant gastrointestinal stromal tumours (biomarker PD model for VEGF, sVEGFR-2, sVEGFR-3, sKIT).
* Add Henin 2009 capecitabine ([doi:10.1038/clpt.2008.220](https://doi.org/10.1038/clpt.2008.220)) [DDMODEL00000214] — adults with metastatic colorectal or advanced/metastatic breast cancer (Markov-proportional-odds model for hand-and-foot syndrome grades 0-2 with K-PD capecitabine exposure).
* Add Plan 2012 pain ([doi:10.1038/clpt.2011.301](https://doi.org/10.1038/clpt.2011.301)) [DDMODEL00000194] — adults in placebo arm of three Phase III neuropathic-pain trials (Markov Integer Model for daily 11-point Likert pain scores).
* Add Kloft 2004 sibrotuzumab ([doi:10.1023/B:DRUG.0000006173.72210.1c](https://doi.org/10.1023/B:DRUG.0000006173.72210.1c)) [DDMODEL00000195] — adults with metastatic FAP-positive solid tumors.
* Change: Standardize Xie 2019 agomelatine parameter, compartment, IOV-eta, and residual-error naming. Compartments rename `DEPOT1` / `DEPOT2` / `LIVER` / `CENTPRNT` / `ALMTPERI` / `METB3OH` / `METB7DM` / `METB7DMPERI` to canonical `depot` / `depot2` / `liver` / `central` / `peripheral1` / `central_3oh` / `central_7dm` / `peripheral1_7dm`; the per-occasion IOV etas (formerly `e.IOV1`–`e.IOV5` and `eta17`–`eta35`) become descriptive `etaiov_<param>_<occ>` names paired with the structural parameter they apply to (`k13`, `alag2`, `k23`, `clint`, `fpop`); residual-error parameters `sdalmt` / `sd3oh` / `sd7dm` become canonical `addSd_lcalmt` / `addSd_lc3oh` / `addSd_lc7dm`; primary/secondary `1` suffixes are dropped where redundant (`fDepot1` -> `fDepot`, `alag1` -> `alag`, `ltvalag1` -> `ltvalag`, `etaltvalag1` -> `etaltvalag`, `F1` -> `fpop`, `etaF1` -> `etafpop`); per-metabolite molecular-weight ratios are renamed descriptively (`mpr1` / `mpr2` -> `mpr_3oh` / `mpr_7dm`). Convention infrastructure (`R/conventions.R`, `R/checkModelConventions.R`) extended to register the `liver` compartment, the `depot[0-9]+` numbered-depot pattern, the `3oh` / `7dm` metabolite suffixes, and the `etaiov_<param>_<occ>` IOV-eta naming pattern.
* Change: Standardize parameter, compartment, covariate-effect, and residual-error naming across the model library. Volumes use `vc` / `vp` / `vp2`; Michaelis-Menten Vmax uses `lvmax` / `vmax`; multi-component CL uses `lcl_ss` / `lcl_time`; shared exponents use `e_<cov>_<param1>_<param2>`; reversed-order covariate effects are flipped to `e_<cov>_<param>`. Parent-metabolite ADC models drop the parent `_adc` suffix and add canonical `<canonical>_<metab>` compartments (`central_mmae`, `central_dxd`, `central_sn38`, `central_tab`, `central_nab`, `central_dm4`, `central_medm4`); metabolite outputs become `Cc_<metab>`. Residual-error parameters now use the parameter-first form `propSd_<X>` / `addSd_<X>` (e.g. `propSd_dxd`, `addSd_tumorSize`) for every non-parent output, while the parent observation `Cc` keeps bare `propSd` / `addSd`. Convention infrastructure (`R/checkModelConventions.R`) extended to enforce the new canonical forms and flag deprecated names.
* Fix Robbie 2012 palivizumab population metadata after revisit against the full-text PDF: pediatric `n` corrected from 1,767 (derived) to 1,684 (paper PK dataset), `sex_female_pct` denominator from 1,660 to 1,684, and added a vignette note clarifying that the 2012 erratum's equation 10 correction (6,900 → 16,900) is a NONMEM `$PRIOR` variance for adult Vp, not a final-model parameter. All `ini()` parameter values continue to match Table 2 exactly.
* Add de Vries Schultink 2020 zenocutuzumab (MCLA-128) ([doi:10.1007/s40262-020-00858-2](https://doi.org/10.1007/s40262-020-00858-2)) — adults with advanced solid tumors (first bispecific antibody in nlmixr2lib).
* Add Almquist 2022 anifrolumab ([doi:10.1002/jcph.2055](https://doi.org/10.1002/jcph.2055)) — adults with moderate-to-severe systemic lupus erythematosus and healthy volunteers (QSS-TMDD with dynamic IFNAR1 receptor pool and time-varying linear clearance).
* Fix Clegg 2024 nirsevimab vignette Figure 4: pass disjoint `id_offset` per cohort and carry `trial` via `rxSolve(keep = )` so the four panels no longer collapse into a single Frankenstein-subject simulation (predictions had been ~3× too high).
* Add Lu 2022 patritumab deruxtecan ([doi:10.1002/jcph.2137](https://doi.org/10.1002/jcph.2137)) — adults with HER3-expressing advanced or metastatic solid tumors (NSCLC, breast cancer, colorectal cancer).
* Add Cheng 2026 immunoglobulin ([doi:10.1002/bcp.70420](https://doi.org/10.1002/bcp.70420)) — children with primary immunodeficiency or secondary antibody deficiency receiving intravenous immunoglobulin replacement therapy.
* Add Frey 2010 tocilizumab ([doi:10.1177/0091270009350623](https://doi.org/10.1177/0091270009350623)) — adults with moderate-to-severe rheumatoid arthritis.
* Add Hayashi 2007 omalizumab ([doi:10.1111/j.1365-2125.2006.02803.x](https://doi.org/10.1111/j.1365-2125.2006.02803.x)) — Japanese adults with atopic asthma or seasonal allergic rhinitis.
* Add Hu 2014 bapineuzumab ([doi:10.1002/jcph.393](https://doi.org/10.1002/jcph.393)) — adults with mild-to-moderate Alzheimer's disease.
* Add Huang 2017 VRC01 ([doi:10.1080/19420862.2017.1311435](https://doi.org/10.1080/19420862.2017.1311435)) — HIV-uninfected healthy adults receiving IV or SC HIV-1 broadly neutralizing monoclonal antibody.
* Add Ide 2020 elotuzumab ([doi:10.1002/jcph.1698](https://doi.org/10.1002/jcph.1698)) — Japanese and non-Japanese adults with multiple myeloma.
* Add Marquez-Megias 2023 adalimumab ([doi:10.3390/biomedicines11102822](https://doi.org/10.3390/biomedicines11102822)) — adults with inflammatory bowel disease.
* Add Nestorov 2014 factor VIII Fc fusion protein ([doi:10.1002/cpdd.167](https://doi.org/10.1002/cpdd.167)) — previously treated patients with severe hemophilia A.
* Add Pérez-Ruixo 2025 posdinemab ([doi:10.1002/cpt.70173](https://doi.org/10.1002/cpt.70173)) — healthy adults and adults with Alzheimer's disease.
* Add Petrov 2024 romiplostim ([doi:10.1002/cpdd.1494](https://doi.org/10.1002/cpdd.1494)) — adults with chronic immune thrombocytopenia (ITP).
* Add Pouzin 2022 tusamitamab ravtansine ([doi:10.1007/s10928-021-09799-0](https://doi.org/10.1007/s10928-021-09799-0)) — adults with advanced solid tumors expressing CEACAM5 (multi-analyte ADC + NAB + DM4 + MeDM4 model).
* Add Takeuchi 2023 ozoralizumab ([doi:10.1002/jcph.2380](https://doi.org/10.1002/jcph.2380)) — Japanese adults with rheumatoid arthritis.
* Add Toukam 2025 BIIB107 ([doi:10.1002/jcph.70109](https://doi.org/10.1002/jcph.70109)) — healthy adult volunteers (first-in-human study supporting MS dose optimization).
* Add Yao 2018 guselkumab ([doi:10.1002/jcph.1063](https://doi.org/10.1002/jcph.1063)) — adults with moderate-to-severe plaque psoriasis.
* Add Frey 2013 tocilizumab DAS28 exposure-response ([doi:10.1177/0091270012437585](https://doi.org/10.1177/0091270012437585)) — adults with moderate-to-severe rheumatoid arthritis (OPTION + TOWARD phase III pool).
* Add Hwang 2023 monalizumab ([doi:10.1002/jcph.2220](https://doi.org/10.1002/jcph.2220)) — adults with advanced solid tumors or squamous cell carcinoma of the head and neck.
* Add Suri 2018 brentuximab vedotin ([doi:10.1002/cpt.1037](https://doi.org/10.1002/cpt.1037)) — adults with CD30-positive Hodgkin lymphoma, systemic anaplastic large-cell lymphoma, mycosis fungoides, or primary cutaneous ALCL pooled across six clinical studies (including the phase III ALCANZA trial in CTCL).
* Add Zhong 2026 abatacept ([doi:10.1002/jcph.70156](https://doi.org/10.1002/jcph.70156)) — pooled adults with rheumatoid arthritis, patients aged 2-17 years with polyarticular juvenile idiopathic arthritis, and patients aged 6+ years with hematologic malignancies receiving HLA-matched unrelated-donor HSCT.
* Add Mukai 2019 mogamulizumab ([doi:10.1002/jcph.1564](https://doi.org/10.1002/jcph.1564)) — adults with cutaneous T-cell lymphoma or adult T-cell lymphoma.
* Add Lin 2024 casirivimab ([doi:10.1007/s11095-024-03764-5](https://doi.org/10.1007/s11095-024-03764-5)) — pediatric and adult subjects, non-infected or with SARS-CoV-2 infection or household contacts.
* Add Faelens 2021 infliximab ([doi:10.3390/pharmaceutics13101623](https://doi.org/10.3390/pharmaceutics13101623)) — adults with moderate-to-severe ulcerative colitis.
* Add Sathe 2024 sacituzumab govitecan ([doi:10.1007/s40262-024-01366-3](https://doi.org/10.1007/s40262-024-01366-3)) — adults with metastatic triple-negative breast cancer and other solid tumors.
* Add Takahashi 2023 abatacept ([doi:10.1182/blood.2023020035](https://doi.org/10.1182/blood.2023020035)) — pooled adult/pediatric rheumatoid-arthritis / pJIA patients and adult/pediatric hematopoietic-cell-transplant recipients in the ABA2 trial.
* Add Hong 2025 datopotamab deruxtecan ([doi:10.1002/psp4.70118](https://doi.org/10.1002/psp4.70118)) — adults with advanced or metastatic NSCLC and breast cancer.
* Add Fau 2020 isatuximab ([doi:10.1002/psp4.12561](https://doi.org/10.1002/psp4.12561)) — adults with relapsed/refractory multiple myeloma.
* Add Lu 2019 polatuzumab vedotin ([doi:10.1002/psp4.12482](https://doi.org/10.1002/psp4.12482)) — adults with relapsed/refractory or previously untreated B-cell non-Hodgkin lymphoma.
* Add Okada 2025 rocatinlimab ([doi:10.1002/psp4.70069](https://doi.org/10.1002/psp4.70069)) — adults with moderate-to-severe atopic dermatitis (pooled with plaque psoriasis, ulcerative colitis, and healthy volunteers).
* Add Zhou 2025 brentuximab vedotin ([doi:10.1002/cpt.3629](https://doi.org/10.1002/cpt.3629)) — pediatric patients (5-18 years) with newly diagnosed Hodgkin lymphoma or relapsed/refractory systemic anaplastic large-cell lymphoma.
* Add Lin 2024 pozelimab ([doi:10.1007/s10928-024-09941-8](https://doi.org/10.1007/s10928-024-09941-8)) — healthy adult volunteers, adults with paroxysmal nocturnal hemoglobinuria, and pediatric and adult patients with CHAPLE disease.
* Add Yin 2021 trastuzumab deruxtecan ([doi:10.1002/cpt.2096](https://doi.org/10.1002/cpt.2096)) — adults with HER2-positive breast cancer or other HER2-expressing solid tumors.
* Add Hwang 2022 tremelimumab ([doi:10.1111/bcp.15622](https://doi.org/10.1111/bcp.15622)) — adults with advanced solid tumours.
* Add Papathanasiou 2025 belantamab mafodotin ([doi:10.1007/s40262-025-01508-1](https://doi.org/10.1007/s40262-025-01508-1)) — adults with relapsed/refractory multiple myeloma (ADC moiety only).
* Add Kuchimanchi 2024 dostarlimab ([doi:10.1111/bcp.16325](https://doi.org/10.1111/bcp.16325)) — adults with primary advanced or recurrent endometrial cancer (RUBY Part 1) and advanced solid tumours (GARNET).
* Add Diao 2016 daclizumab CD25 occupancy ([doi:10.1111/bcp.13051](https://doi.org/10.1111/bcp.13051)) — adults with relapsing-remitting multiple sclerosis (PD model with Othman 2014 PK backbone).
* Add Diao 2016 daclizumab CD56 bright NK expansion ([doi:10.1111/bcp.13051](https://doi.org/10.1111/bcp.13051)) — adults with relapsing-remitting multiple sclerosis (PD model with Othman 2014 PK backbone).
* Add Diao 2016 daclizumab Treg reduction ([doi:10.1111/bcp.13051](https://doi.org/10.1111/bcp.13051)) — adults with relapsing-remitting multiple sclerosis (PD model with Othman 2014 PK backbone).
* Add Fiedler-Kelly 2020 fremanezumab exposure-response ([doi:10.1111/head.13855](https://doi.org/10.1111/head.13855)) — adults with episodic migraine and adults with chronic migraine (two PD-only models: `FiedlerKelly_2020_fremanezumab_em` and `FiedlerKelly_2020_fremanezumab_cm`).
* Add Koopman 2023 factor IX-Fc ([doi:10.1111/bcp.15881](https://doi.org/10.1111/bcp.15881)) — children, adolescents, and adults with haemophilia B (real-world data including patients aged < 12 years).
* Add Le Tilly 2021 trastuzumab ([doi:10.1002/cpt.2188](https://doi.org/10.1002/cpt.2188)) — adults with HER2+ breast cancer leptomeningeal carcinomatosis receiving intrathecal and intravenous trastuzumab.
* Add Wang 2024 sugemalimab ([doi:10.1111/bcp.16276](https://doi.org/10.1111/bcp.16276)) — adults with advanced solid tumours or lymphomas across nine Phase I-III sugemalimab trials.
* Add Yang 2024 axatilimab ([doi:10.1002/cpt.3503](https://doi.org/10.1002/cpt.3503)) — pooled healthy adults, adults with advanced solid tumors, and adults / children with chronic graft-versus-host disease.
* Add Hood 2021 MEDI7836 ([doi:10.3390/pharmaceutics13050613](https://doi.org/10.3390/pharmaceutics13050613)) — healthy adult males in a first-in-human single-ascending-dose trial.
* Amend Castro-Suarez 2020 nimotuzumab: V1 decreased 53% for 50 mg cohort per Figure 4 visual inspection (direction not stated in paper text; corresponding author contacted).
* Add Castro-Suárez 2020 nimotuzumab ([doi:10.3390/pharmaceutics12121147](https://doi.org/10.3390/pharmaceutics12121147)) — adults with autosomal dominant polycystic kidney disease.
* Add Yang 2021 cemiplimab ([doi:10.1007/s10928-021-09739-y](https://doi.org/10.1007/s10928-021-09739-y)) — adults with advanced solid tumors including cutaneous squamous cell carcinoma.
* Add Papachristos 2020 bevacizumab ([doi:10.3390/ijms21113753](https://doi.org/10.3390/ijms21113753)) — adults with metastatic colorectal cancer (three co-equal final models: descriptive PK, binding QSS TMDD, and PK/PD Imax inhibition of free VEGF-A).
* Add Ngo 2020 HL2351 ([doi:10.1002/psp4.12552](https://doi.org/10.1002/psp4.12552)) — healthy adult Korean men.
* Add Wojciechowski 2022 domagrozumab ([doi:10.1111/cts.13418](https://doi.org/10.1111/cts.13418)) — healthy adult volunteers and pediatric patients with Duchenne muscular dystrophy.
* Add Yu 2022 ofatumumab ([doi:10.1007/s40263-021-00895-w](https://doi.org/10.1007/s40263-021-00895-w)) — adults with relapsing multiple sclerosis.
* Add Melhem 2022 dostarlimab ([doi:10.1111/bcp.15339](https://doi.org/10.1111/bcp.15339)) — adults with advanced solid tumours.
* Add Brillac 2025 isatuximab ([doi:10.1007/s00280-025-04832-2](https://doi.org/10.1007/s00280-025-04832-2)) — pediatric and adult patients with relapsed/refractory acute leukemias.
* Add Wu 2024 inotuzumab ozogamicin ([doi:10.1007/s40262-024-01386-z](https://doi.org/10.1007/s40262-024-01386-z)) — pediatric and adult patients with relapsed/refractory B-cell precursor acute lymphoblastic leukemia and adults with B-cell non-Hodgkin's lymphoma.
* Fix vignettes: derive concentration units in `labs()` from model `$units` metadata; replace inline trapezoidal NCA with PKNCA; add PKNCA sections to nalmefene and clesrovimab vignettes; add PKNCA treatment grouping to benralizumab and durvalumab vignettes.
* Add Chen 2022 guselkumab ([doi:10.1111/cts.13197](https://doi.org/10.1111/cts.13197)) — adults with active psoriatic arthritis (DISCOVER-1 and DISCOVER-2 phase 3 trials).
* Add `concentration` field to `units` metadata for `Wang_2017_benralizumab` and `Ogasawara_2020_durvalumab` models.
* Add Chen 2020 luspatercept ([doi:10.1002/psp4.12515](https://doi.org/10.1002/psp4.12515)) — adults with anemia due to lower-risk myelodysplastic syndromes.
* Add Martinez 2019 alirocumab ([doi:10.1007/s40262-018-0669-y](https://doi.org/10.1007/s40262-018-0669-y)) — healthy volunteers and adults with hypercholesterolemia.
* Add Li 2017 brentuximab vedotin ([doi:10.1002/jcph.920](https://doi.org/10.1002/jcph.920)) — adults with relapsed/refractory CD30-expressing hematologic malignancies (Hodgkin lymphoma and systemic anaplastic large cell lymphoma).
* Add Quartino 2019 trastuzumab ([doi:10.1007/s00280-018-3728-z](https://doi.org/10.1007/s00280-018-3728-z)) — adults with metastatic breast cancer, early breast cancer, advanced gastric cancer, or other solid tumors.
* Add Suleiman 2019 risankizumab ([doi:10.1007/s40262-019-00759-z](https://doi.org/10.1007/s40262-019-00759-z)) — healthy volunteers and adults with moderate-to-severe plaque psoriasis.
* Add Berends 2019 infliximab ([doi:10.1007/s10928-019-09652-5](https://doi.org/10.1007/s10928-019-09652-5)) — adults with moderate-to-severe ulcerative colitis.
* Add Nikanjam 2019 siltuximab ([doi:10.1007/s00280-019-03939-7](https://doi.org/10.1007/s00280-019-03939-7)) — adults pooled across healthy volunteers and oncology cohorts including Castleman's disease and smoldering multiple myeloma.
* Add Sanghavi 2020 ipilimumab ([doi:10.1002/psp4.12477](https://doi.org/10.1002/psp4.12477)) — adults with advanced solid tumors receiving ipilimumab alone or in combination with nivolumab.
* Add Zhang 2019 nivolumab ([doi:10.1002/psp4.12476](https://doi.org/10.1002/psp4.12476)) — adults with advanced solid tumors, alone or in combination with ipilimumab or chemotherapy.
* Add Wang 2020 ontamalimab ([doi:10.1002/jcph.1590](https://doi.org/10.1002/jcph.1590)) — adults with moderate-to-severe ulcerative colitis or Crohn's disease.
* Add Hanzel 2021 infliximab CT-P13 ([doi:10.1111/apt.16609](https://doi.org/10.1111/apt.16609)) — adults with Crohn's disease or ulcerative colitis.
* Add Zhou 2021 belimumab ([doi:10.1007/s40268-021-00363-2](https://doi.org/10.1007/s40268-021-00363-2)) — adult and pediatric patients with systemic lupus erythematosus (Chinese and non-Chinese).
* Add Aguiar 2021 ustekinumab ([doi:10.3390/pharmaceutics13101587](https://doi.org/10.3390/pharmaceutics13101587)) — adults with active Crohn's disease.
* Add Li 2019 abatacept ([doi:10.1002/jcph.1308](https://doi.org/10.1002/jcph.1308)) — adults with rheumatoid arthritis.
* Add Gandhi 2021 abatacept ([doi:10.1002/jcph.1781](https://doi.org/10.1002/jcph.1781)) — pooled adults with rheumatoid arthritis and patients aged 2-17 years with polyarticular juvenile idiopathic arthritis.
* Add Mulyukov 2018 ranibizumab ([doi:10.1002/psp4.12322](https://doi.org/10.1002/psp4.12322)) — anti-VEGF-naive adults with neovascular age-related macular degeneration.
* Add Bajaj 2017 nivolumab ([doi:10.1002/psp4.12143](https://doi.org/10.1002/psp4.12143)) — patients with advanced solid tumors (melanoma, NSCLC, RCC, other).
* Add Kielbasa 2020 galcanezumab ([doi:10.1002/jcph.1511](https://doi.org/10.1002/jcph.1511)) — healthy adults and adults with episodic or chronic migraine.
* Add Othman 2014 daclizumab ([doi:10.1007/s40262-014-0159-9](https://doi.org/10.1007/s40262-014-0159-9)) — healthy adult volunteers.
* Add Timmermann 2019 brodalumab ([doi:10.1111/bcpt.13202](https://doi.org/10.1111/bcpt.13202)) — adults with moderate-to-severe plaque psoriasis.
* Add Kuchimanchi 2018 evolocumab ([doi:10.1007/s10928-018-9592-y](https://doi.org/10.1007/s10928-018-9592-y)) — healthy adults and patients with hypercholesterolemia (including heterozygous familial hypercholesterolemia and statin-intolerant cohorts).
* Add Mo 2018 olaratumab ([doi:10.1007/s40262-017-0562-0](https://doi.org/10.1007/s40262-017-0562-0)) — patients with advanced or metastatic cancer (soft tissue sarcoma, NSCLC, GIST, glioblastoma multiforme).
* Add Kovalenko 2016 dupilumab ([doi:10.1002/psp4.12136](https://doi.org/10.1002/psp4.12136)) — healthy volunteers and adults with moderate-to-severe atopic dermatitis.
* Add Narwal 2013 sifalimumab ([doi:10.1007/s40262-013-0085-2](https://doi.org/10.1007/s40262-013-0085-2)) — adults with moderately active systemic lupus erythematosus.
* Add Valenzuela 2025 nipocalimab ([doi:10.1002/psp4.70109](https://doi.org/10.1002/psp4.70109)) — healthy adults and adults with generalized myasthenia gravis.
* Add Hua 2015 anrukinzumab ([doi:10.1111/bcp.12589](https://doi.org/10.1111/bcp.12589)) — healthy volunteers, asthma, and ulcerative colitis patients.
* Add Robbie 2012 palivizumab ([doi:10.1128/AAC.06446-11](https://doi.org/10.1128/AAC.06446-11)) — adults and pediatric patients at high risk of RSV lower respiratory tract disease.
* Add Bender 2014 trastuzumab emtansine (T-DM1, ADC) ([doi:10.1208/s12248-014-9618-3](https://doi.org/10.1208/s12248-014-9618-3)) — preclinical: rats + cynomolgus monkeys; two model variants (reduced + mechanistic DAR chain).
* Add Rosario 2015 vedolizumab ([doi:10.1111/apt.13243](https://doi.org/10.1111/apt.13243)) — adults with moderately-to-severely active ulcerative colitis or Crohn's disease.
* Add Long 2017 necitumumab ([doi:10.1007/s40262-016-0452-x](https://doi.org/10.1007/s40262-016-0452-x)) — cancer patients with advanced solid tumors.
* Add Wade 2015 certolizumab ([doi:10.1002/jcph.491](https://doi.org/10.1002/jcph.491)) — adults with Crohn's disease.
* Add Lon 2013 abatacept ([doi:10.1007/s10928-013-9341-1](https://doi.org/10.1007/s10928-013-9341-1)) — male Lewis rats with collagen-induced arthritis (preclinical).
* Add Lu 2014 trastuzumab emtansine ([doi:10.1007/s00280-014-2500-2](https://doi.org/10.1007/s00280-014-2500-2)) — adults with HER2-positive locally advanced or metastatic breast cancer.
* Add Gupta 2016 amatuximab ([doi:10.1007/s00280-016-2984-z](https://doi.org/10.1007/s00280-016-2984-z)) — adults with advanced mesothelin-expressing cancers including malignant pleural mesothelioma.
* Add Zheng 2016 sifalimumab ([doi:10.1111/bcp.12864](https://doi.org/10.1111/bcp.12864)) — adults with systemic lupus erythematosus.
* Add Mould 2007 alemtuzumab ([doi:10.1111/j.1365-2125.2007.02914.x](https://doi.org/10.1111/j.1365-2125.2007.02914.x)) — adults with B-cell chronic lymphocytic leukaemia.
* Add Farrell 2012 farletuzumab ([doi:10.1007/s00280-012-1959-y](https://doi.org/10.1007/s00280-012-1959-y)) — women with advanced epithelial ovarian cancer.
* Add Xu 2011 sirukumab ([doi:10.1111/j.1365-2125.2011.03964.x](https://doi.org/10.1111/j.1365-2125.2011.03964.x)) — healthy adult volunteers in a first-in-human study.
* Add Yamada 2025 zolbetuximab ([doi:10.1111/cts.70280](https://doi.org/10.1111/cts.70280)) — patients with locally advanced unresectable or metastatic gastric/gastroesophageal junction adenocarcinoma.
* Add Jackson 2022 ixekizumab ([doi:10.1111/bcp.15034](https://doi.org/10.1111/bcp.15034)) — paediatric patients with moderate-to-severe plaque psoriasis.
* Add Moein 2022 etrolizumab ([doi:10.1002/psp4.12846](https://doi.org/10.1002/psp4.12846)) — patients with moderately-to-severely active ulcerative colitis.
* Add Chua 2025 mirikizumab ([doi:10.1111/cts.70320](https://doi.org/10.1111/cts.70320)) — patients with moderately-to-severely active Crohn's disease.
* Add Ma 2020 sarilumab ANC ([doi:10.1007/s40262-020-00899-7](https://doi.org/10.1007/s40262-020-00899-7)) — adults with rheumatoid arthritis.
* Add Ma 2020 sarilumab DAS28-CRP ([doi:10.1007/s40262-020-00899-7](https://doi.org/10.1007/s40262-020-00899-7)) — adults with rheumatoid arthritis.
* Add Tiraboschi 2025 amlitelimab ([doi:10.1002/psp4.70121](https://doi.org/10.1002/psp4.70121)) — adults with moderate-to-severe atopic dermatitis.
* Add Budha 2023 tislelizumab ([doi:10.1002/psp4.12880](https://doi.org/10.1002/psp4.12880)) — patients with advanced tumors.
* Add Masters 2022 avelumab ([doi:10.1002/psp4.12771](https://doi.org/10.1002/psp4.12771)) — patients with advanced solid tumors.
* Add Xu 2019 sarilumab ([doi:10.1007/s40262-019-00765-1](https://doi.org/10.1007/s40262-019-00765-1)) — adults with rheumatoid arthritis.
* `checkModelConventions()` — new function that reports deviations from the package's parameter-naming, covariate, compartment, and metadata conventions for a single model or the entire `modeldb`. Called automatically during `buildModelDb()` so convention drift surfaces at package-build time (existing grandfathered deviations continue to build). Canonical standards (PK parameter prefixes, `eta`-prefixed IIV, `propSd`/`addSd` residual error, canonical covariate column register with aliases, canonical compartment vocabulary) are codified in an internal `.nlmixr2libConventions` list that mirrors the `extract-literature-model` skill references. Addresses issue #39.
* Added canonical TMDD archetype models under `inst/modeldb/pharmacokinetics/`: `PK_1cmt_tmdd_full` (Mager & Jusko 2001), `PK_1cmt_tmdd_qss`, `PK_1cmt_tmdd_mm`, `PK_2cmt_tmdd_qss`, and `PK_2cmt_tmdd_mm` (Gibiansky et al. 2008 QSS and MM approximations). Built to the `extract-literature-model` skill conventions with `reference` / `units` / `population` metadata, per-parameter source-trace comments, and CL/V parameterization with `kel` derived inside `model()` so later transforms can re-parameterize. Replaces the 38 draft models from PR #60 (#60).
* Registered `target`, `complex`, and `total_target` as canonical TMDD compartment names in the `extract-literature-model` skill naming conventions (per @iamstein's proposal on PR #60).
* Added `vignettes/tmdd_archetypes.Rmd` comparing the five TMDD archetypes with typical-value trajectories of drug, free target, and complex, and a regime-convergence check showing QSS/MM collapsing onto the full Mager & Jusko 2001 model in the fast-binding regime of Gibiansky 2008.
* Add Thakre 2022 risankizumab ([doi:10.1007/s40744-022-00495-0](https://doi.org/10.1007/s40744-022-00495-0)) — patients with active psoriatic arthritis.
* `addEta()`, `addResErr()`, `addDepot()`, `removeDepot()`, `addTransit()`, and `removeTransit()` now accept `model` as a deprecated alias for `ui` (issue #84). Passing `model = ...` emits a deprecation warning; passing both `ui` and `model` is an error.
* `addDepot()` and `addTransit()` now work correctly when `d/dt(central)` or `d/dt(depot)` appears at the beginning or end of the model block, or when transit-compartment ODEs and residual-error (`~`) specs are interleaved with assignment lines. The newly introduced helper and ODE lines are inserted immediately adjacent to the modified ODE so that the relative order of every pre-existing model line is preserved (#77, #78).
* Markov modeling creation functions including `createMarkovModel()` were added
* Add Kotani 2022 astegolimab ([doi:10.1002/jcph.2021](https://doi.org/10.1002/jcph.2021)) — adults with severe asthma.
* Add Fasanmade 2009 infliximab ([doi:10.1007/s00228-009-0718-4](https://doi.org/10.1007/s00228-009-0718-4)) — adults with moderately-to-severely active ulcerative colitis.
* Add Fiedler-Kelly 2019 fremanezumab ([doi:10.1111/bcp.14096](https://doi.org/10.1111/bcp.14096)) — adults with chronic or episodic migraine.
* Add Hu 2026 clesrovimab ([doi:10.1002/cpt.70199](https://doi.org/10.1002/cpt.70199)) — preterm and full-term infants.
* Add Clegg 2024 nirsevimab ([doi:10.1002/jcph.2401](https://doi.org/10.1002/jcph.2401)) — preterm and term infants.
* Verified all published-literature specific-drug and mAb-consensus models against their source papers and fixed several parameter-encoding bugs that had been latent in the package since their original addition:
  - **CarlssonPetri 2021 liraglutide**: fixed categorical covariate encoding that was zeroing individual clearance for subjects not in the indexed group. `(1 - SEXF)^e_sex_cl` → `e_sex_cl^(1 - SEXF)` (previously evaluated `0^1.12 = 0` for females); `CHILD^e_age_child_cl * ADOLESCENT^e_age_adolescent_cl` → `e_age_child_cl^CHILD * e_age_adolescent_cl^ADOLESCENT` (previously evaluated `0^1.11 * 0^1.06 = 0` for adults). IIV rewritten as `omega^2 = log(1 + CV^2)` per Table 3's explicit `%CV = sqrt(exp(omega^2) - 1) * 100` footnote.
  - **Zhu 2017 lebrikizumab**: fixed IIV variance-covariance block that was storing `sqrt(variance)` (SDs) instead of variances/covariances from Table 3.
  - **Soehoel 2022 tralokinumab**: fixed IIV block (SDs → variance-covariance matrix from Table 2 footnote `IIV = sqrt(exp(omega^2) - 1)`, correlation 0.61 applied on variances not SDs); corrected body-weight exponent on V2/V3 from `0.793` to Table 2's `0.783`.
  - **Kovalenko 2020 dupilumab**: squared the five IIV values so they store variances. Paper Methods explicitly defines `omega` as "the standard deviation [SD] of between-subject variability", and nlmixr2 `etaX ~ value` stores the variance (omega^2), so SDs needed squaring.
  - **Davda 2014 mAb consensus (PK_2cmt_mAb_Davda_2014)**: all parameters verified; no changes required.
* Filled in previously-TODO `population` metadata blocks for the five verified models above with demographics from each paper's Table 1.
* Added validation vignettes for the five verified models above (`CarlssonPetri_2021_liraglutide.Rmd`, `Zhu_2017_lebrikizumab.Rmd`, `Soehoel_2022_tralokinumab.Rmd`, `Kovalenko_2020_dupilumab.Rmd`, `PK_2cmt_mAb_Davda_2014.Rmd`). Each vignette follows the `extract-literature-model` skill conventions with a population description, source trace, virtual cohort, simulation, figure replication, and PKNCA-based NCA validation.
* Retrofit Cirincione 2017 exenatide model to the `extract-literature-model` skill conventions and fix parameter encoding bugs: `ka_max` corrected from `0.0813` to paper value `12.8` /hr, `Km` rescaled to ng/mL so units are consistent with `Cc = central / vc`, and IIV variances rewritten as `log(1 + CV^2)` rather than the `log(1 + CV)` shortcut. Replaces the character-valued `DVID` covariate with `STUDY1` / `STUDY5` binary indicators and adds a companion validation vignette with PKNCA checks against the paper's Figure 5 typical values.
* Fix endogenous-model bugs surfaced during a parameter audit:
  - **Charbonneau 2021 phenylalanine**: fixed `vd_phe` typo on the `f_gut_plasma` line (the undefined symbol made `rxSolve()` fail) — should read `vd`; removed a stray `vd *` factor from the `daily_phe_intake` augmentation output so the reported value is in mg/day rather than (L/kg)·mg/day. All Table 4 parameter values verified against the authors' Zenodo companion notebook.
  - **Kim 2006 IgG**: corrected `V1` units label from `(mg/kg)` to `(mL/kg)` (paper Table 1); closed a missing `)` in the `ljmax` label; removed a redundant `igg_0 <- css` line. Parameter values verified against Table 1.
* Extended the `extract-literature-model` skill to cover endogenous, mechanistic, and turnover models: added naming conventions for Vmax/Km/baseline/rate parameters, a new `references/endogenous-validation.md` reference covering steady-state / perturbation-recovery / mass-balance / dimensional-analysis patterns, and an explicit dimensional-analysis item in `references/verification-checklist.md` (driven by the Kim units-label and Charbonneau `daily_phe_intake` bugs, which were only catchable by dimensional analysis).
* Retrofit the remaining specific-drug models and the Davda 2014 mAb 2-compartment model to the `extract-literature-model` skill conventions (no parameter-value changes). Structured `covariateData` with `description` / `units` / `type` / `reference_category` / `notes` / `source_name`, canonical `units` list (time/dosing/concentration), and `population` metadata blocks added (with TODO placeholders where not yet sourced). IIV etas renamed to `eta` + transformed-parameter name (e.g., `etacl` → `etalcl`, `iiv_lka` → `etalka`, `bsv_fpla_*` → `etalfpla_*`). Race columns renamed to the `RACE_` prefix (`BLACK` → `RACE_BLACK`, `ASIAN` → `RACE_ASIAN`, `MULTIRACIAL` → `RACE_MULTI`, `BLACK_OTH` → `RACE_BLACK_OTH`, `ASIAN_AMIND_MULTI` → `RACE_ASIAN_AMIND_MULTI`); `ADA` → `ADA_POS`; `SEXM` → `SEXF` (value inverted to keep the original effect magnitude). Touched vignettes updated to use the canonical column names. Models touched: `CarlssonPetri_2021_liraglutide`, `Clegg_2024_nirsevimab`, `Grimm_2023_gantenerumab`, `Grimm_2023_trontinemab`, `Hu_2026_clesrovimab`, `Kovalenko_2020_dupilumab`, `Kyhl_2016_nalmefene`, `PK_2cmt_mAb_Davda_2014`, `Soehoel_2022_tralokinumab`, `Xie_2019_agomelatine`, `Zhu_2017_lebrikizumab`.

# Version 0.3.2

* Add Kim 2006 model for IgG metabolism
* Add Xie 2019 agomelatine PK model
* Drop `qs` since it will be archived
* Update tumor growth inhibition models
* `addResErr()` now works with multiple-endpoint models
* Additional testing

# Version 0.3.1

* Bug fix for replacement of multiplicative expressions in `nlmixr2lib`
* phenylalanine_charbonneau_2021 had its net protein breakdown parameter corrected
* Kyhl_2016_nalmefene model was added

# Version 0.3.0

* Added ability to choose style type when modifying models.  Currently
  supported styles are: "camel" for `variablesLikeThis`, "snake" for
  `variables_like_this`, "dot" for `variables.like.this` and "blank"
  for `variableslikethis`.  This can be selected with
  `setCombineType()`.

* With the new combination style, you can change how `eta` variables
  are constructed with the `option(nlmixr2lib.etaCombineType="camel")`
  or whatever you wish it to the variable style to be.

* Added new model building framework for building models

  - **PK model building functions**

     - `addTransit()`/`removeTransit()` which were present before, but now modified and
       made a bit more robust, more closely matching literature method
       of transit compartments.

     - `addDepot()`/`removeDepot()` which were present before, but
       modified to be a bit more robust.

     - `addWeibullAbs()` which adds a Weibull absorption to a PK model

     - `convertMM()` converts linear elimination to Michaelis-Menten elimination

     - `transPK()` converts the `cl` style parameter transformations
       to various other PK transformations like `k`, `aob`, `alpha`,
       `k12`

  - **PD model building functions**

   - `addIndirectLin()` -- this adds an indirect effect model to a PK
     model that has a concentration `Cc` in the model.  This purposely
     uses a simple linear effect of `Cc*Ek` or `Cc*Ik` so it will be
     easy to parse and turn into other functional forms (like `Emax`
     or `Hill`).  If the PK model is not present it will use `Cc` as a
     covariate in a purely PD models.

   - `addIndirect()` -- this builds on `addIndirectLin()` and adds
     `Emax` or `Hill` models to a PK model. You can also set `imax=1`
     or `emax=1` to drop these parameters from being estimated in the
     model.  Additionally `hill=TRUE` will add a Hill coefficient to
     the sigmoid model.

   - `addEffectCmtLin()` -- this adds an effect compartment based on
     the `Cc` in the model.  The linear effect can be modified into
     other function forms.

   - `addDirectLin()` -- this adds a direct effect model based on the
     `Cc` in the model.

   - **Changing functional forms of Effect models**

     - `convertEmax()` changes linear effect models to Emax models

     - `convertEmaxHill()` changes linear effect models to Hill models

     - `convertQuad()` changes linear effect models to quadratic models

     - `convertLogLin()` changes linear effect models to log-linear models

   - **Changing functional forms of Baselines in non-indirect response models**

     - `addBaselineConst()` changes the zero baseline to a estimated
       constant

     - `addBaselineLin()` changes the zero baseline to a estimated
       constant and a linear constant with respect to `time`.

     - `addBaselineExp()` changes the zero baseline to a exponential
       decay with respect to time

     - `addBaseline1exp()` -- the baseline effect is changed from zero
       to to an exponential approaching to a constant (with respect to
       time).

   - **Changing model properties** (all use `addCmtProp()`)

      - `addBioavailability()` adds bioavailability property to a
        compartment

      - `addRate()` adds a modeled rate to a compartment

      - `addDur()` adds modeled duration to a compartment

      - `addIni()` adds an initial value to a compartment

      - `addLag()` adds a lag time to the a compartment

* Add Carlsson Petri (2021) liraglutide PK model
* Add Cirincione (2017) exenatide immediate-release PK model
* Add a variety of indirect response models
* Add a variety of tumor growth inhibition models and move all oncology models
  into a new model database directory
* Add a variety of double-absorption PK models
* `cp` and related `cpddSd` and `cppropSd` were renamed to `Cc`, `CcAddSd` and
  `CcPropSd` (fix #70).
* Multiple-endpoint models will have the `DV` column in the modeldb separated by
  commas.

# Version 0.2.0

* Work with the new `rxode2` version 2.0.12 `model()` and `ini()` assignment
  methods.
* Therapeutic-area specific models have begun being added.
* Models can now give the user some additional information load via the
  `message` meta-data.
* Models can now be in different directories.  The change is for ease of
  maintaining the library, it is not a change that affects users.
* A regression where `addEta()` did not change the parameter, related to a
  change in `rxode2`, was fixed.
* `addEta()` detects where to add etas more robustly when covariates are on the
  parameter.

## Models added

* Add Davda (2014) mAb consensus model
* Add Liu (2017) time-dependent clearance model based on nivolumab
* Add Kovalenko (2020) dupilumab PK model
* Add Soehoel (2022) tralokinumab PK model
* Add Zhu (2017) lebrikizumab PK model

# Version 0.1.0

* Initial version
