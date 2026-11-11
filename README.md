# Ground-truth simulation of EEG beta-burst source localization

A simulation framework that tests whether transient beta bursts can be **detected and localized** from
32-channel scalp EEG with a template head model. Bursts are placed at known cortical and deep locations,
projected to the scalp through a realistic head model, buried in realistic background activity, then
detected and localized by several standard pipelines. Every map is scored against the ground truth.

The question is practical: **before a burst source map from real EEG is interpreted, how accurate can that
map be?**

In this setup, the answer is: **at 0 dB SNR or better, correct pipelines are about 1 cm accurate and identify
M1 reliably; below that, burst detection (not the inverse) fails. Premotor and SMA sources 14 mm apart are
partly confused by dSPM, a deep source always appears on the cortex, noise alone still produces a
plausible-looking peak, and a common shortcut (one threshold pooled across all channels) does not localize
at all.**

![Example: one simulated recording, three detectors and seven localization pipelines](figures/fig1_example.png)

> **Version 2 (Nov 2026).** Electrode positions now come from the 10-05 template, which fits fsaverage.
> MNE's `biosemi32` template positions sit about 2 cm posterior on fsaverage; in version 1 this produced an
> apparent M1 → premotor bias that is now gone. All numbers below are from the corrected run.

> **Scope.** One template head (fsaverage), one montage (BioSemi-32 channel names at 10-05 positions), simulated bursts. The framework shows
> what such a pipeline can and cannot resolve; it does not validate any particular real-data map. All
> numbers below come from the reference run (20 seeds per cell, 120 s recordings).

---

## Contents

| File | What it is |
|---|---|
| `beta_burst_sim.ipynb` | The full framework: simulation, detection, localization, scoring, experiments, figures |
| `requirements.txt` | Python dependencies |
| `figures/` | Figures from the reference run |
| `results/` | Per-recording result tables (CSV) from the reference run |
| `LICENSE` | MIT |
| `real_data/` | Application to three open Parkinson's disease EEG cohorts, with all sensitivity analyses (see below) |

## Quick start

```bash
git clone https://github.com/tanvir-hasib/beta-burst-localization-simulation-model.git
cd beta-burst-localization-simulation-model
pip install -r requirements.txt
jupyter notebook beta_burst_sim.ipynb      # then Run All
```

The first run downloads the fsaverage template (about 200 MB, via `mne.datasets.fetch_fsaverage`) and
builds two 3-layer BEMs. Run time for the reference run is about 40 min on a desktop CPU.

| Switch (top of notebook) | Default | Effect |
|---|---|---|
| `HEAD_MODEL` | `"fsaverage"` | `"offline"` uses a no-download approximation (fsaverage5 cortex from nilearn, MNE's bundled inner skull). Use it only for testing. |
| `QUICK_MODE` | `False` | `True` gives a ~5 min smoke test (3 seeds, 2 SNRs, 60 s recordings) |
| `BURST_MODEL` | `"event"` | `"jr"` generates bursts from a noise-driven Jansen–Rit column instead of windowed oscillations |
| `REAL_BACKGROUND_FILE` | `None` | A resting EEG file to use as background instead of simulated activity |
| `CFG["run_matched_control"]` | `True` | Run the matched-head-model control (section 7) |

All parameters live in the `CFG` dictionary. Seeds are fixed, so a run reproduces the reference numbers.

## How the simulation works

```
burst generator ─► source patch (10 mm) ─► generating head model ─► + background ─► scalp EEG (avg ref)
                                                                                      │
          scores ◄─ map vs truth ◄─ inverse (inverse head model) ◄─ trough-aligned events ◄─ detection
```

| Stage | Implementation |
|---|---|
| **Bursts** | Hann-windowed 16–24 Hz events with log-normal duration and amplitude, sharpened troughs (non-sinusoidal beta). Optional Jansen–Rit column driven across its Hopf point. Ground truth: onset, offset and deepest trough of every burst. |
| **Sources** | Seven locations (fsaverage ≈ MNI): left and right M1 hand area, left and right dorsal premotor (PMd), left SMA, a left occipital control, and a deep left thalamic dipole. Cortical sources are 10 mm patches with fixed normal orientation. |
| **Head models** | fsaverage, 3-layer BEM, BioSemi-32 channels at 10-05 template positions (`colin27_1005`). Data are **generated** with skull conductivity 0.0042 S/m, ico5 sources (20 484) and electrodes jittered by 5 mm, and **inverted** with 0.006 S/m, ico4 sources (5 124) and template electrodes. This avoids the "inverse crime". |
| **Background** | 300 random cortical sources with 1/f activity, posterior alpha, weak ongoing beta, white sensor noise. SNR = burst RMS ÷ background beta-band RMS at the strongest channel. |
| **Detection** | One beta envelope, 75th-percentile threshold, ≥ 2 cycles, gaps < 50 ms merged, one event per burst at its deepest trough. |
| **Localization** | Noise covariance from burst-free data, `nave` = number of bursts, depth weighting 0.8. Maps from trough-aligned averages (±25 ms RMS), except LCMV, which uses beta-band power in bursts ÷ burst-free power. |

### Pipelines compared

| Detector | Detection signal |
|---|---|
| `pooled_any` | One percentile pooled over all channels; a sample counts if **any** channel exceeds it. Every flagged sample is a trigger. Included as a negative control. |
| `channel` | One fixed channel (C3) |
| `ged` | Leading generalized-eigendecomposition component, beta band vs broadband |
| `oracle` | The true burst troughs (upper bound) |

Localization methods: **dSPM**, **sLORETA**, **eLORETA**, **LCMV** beamformer, **single-dipole scan**,
each implemented directly on the gain matrix (whitened, average-reference null space removed).

### Scores

| Score | Definition |
|---|---|
| Peak error | Distance (mm) from the map maximum to the true source centroid |
| Spatial dispersion | √(Σ m²d² / Σ m²); for LCMV the weights are max(ratio − 1, 0)² |
| AUC | Vertices inside the true patch vs all others |
| Region accuracy | Map maximum closer to the true region than to any other region |
| Hemisphere accuracy | Map maximum in the true hemisphere |
| Detection F1 | Detected vs true troughs matched within ±100 ms |

## Results (reference run)

### 1. The pooled "any channel" threshold does not localize

The pooled rule flags **94 %** of all samples as burst (detection F1 = 0.003). Its maps peak **68–74 mm** from
the source and pick the right region in **13–20 %** of recordings (about chance among six cortical regions)
at every SNR. Its error does not improve when the bursts get stronger.

### 2. With a correct pipeline, detection limits accuracy

| SNR | GED F1 | dSPM, GED events: error / region acc. | dSPM, true events: error / region acc. | LCMV, GED events: error / region acc. |
|---|---|---|---|---|
| −10 dB | 0.18 | 69 mm / 18 % | 14 mm / 86 % | 72 mm / 20 % |
| −5 dB | 0.25 | 40 mm / 37 % | 12 mm / 88 % | 33 mm / 44 % |
| 0 dB | 0.66 | 14 mm / 77 % | 12 mm / 89 % | 11 mm / 79 % |
| +5 dB | 0.84 | 12 mm / 91 % | 12 mm / 91 % | 11 mm / 88 % |

Peak error pooled over all seven locations. With the true burst times, the error barely changes with SNR
(14 → 12 mm): at −10 and −5 dB almost all of the error comes from detection. A single fixed channel (C3)
never exceeds F1 = 0.48. At +5 dB, sLORETA, eLORETA and the dipole scan reach region accuracies of
0.98–0.99 (dSPM 0.91, LCMV 0.88).

![Localization error vs SNR, per source location](figures/fig2_error_vs_snr.png)
![Detection accuracy](figures/fig4_detection.png)

### 3. M1 is resolved; premotor cortex is partly confused with SMA by dSPM

At +5 dB, with dSPM:

| True source | Peak error | Region accuracy (GED events) | Region accuracy (true events) |
|---|---|---|---|
| L M1 hand | 10.9 mm | 100 % | 100 % |
| R M1 hand | 9.8 mm | 100 % | 100 % |
| L PMd | 4.4 mm | **70 %** (30 % → SMA) | 65 % |
| R PMd | 13.9 mm | **75 %** (25 % → SMA) | 80 % |
| L SMA | 11.1 mm | 100 % | 100 % |
| L occipital | 9.2 mm | 100 % | 100 % |

M1 and PMd centroids are 23 mm (left) and 19 mm (right) apart and are separated reliably. PMd and SMA are
only 14 mm apart; dSPM's depth weighting pulls some PMd maxima medially onto SMA, the same with true and
detected events, so this is a resolution limit, not a detection error. sLORETA, eLORETA and the dipole scan
separate them almost perfectly at this SNR.

![Region confusion, dSPM at +5 dB](figures/fig3_confusion_dSPM.png)

### 4. Head-model mismatch is not the main source of error

The control generates data with the inverse head model itself (same conductivities and electrode positions;
only the source-space density differs). Median peak error:

| Source | dSPM, true events: mismatched | matched | LCMV: mismatched | matched |
|---|---|---|---|---|
| L M1 hand | 10.9 mm | 12.1 mm | 10.9 mm | 25.7 mm |
| R M1 hand | 11.0 mm | 8.6 mm | 8.6 mm | 8.6 mm |
| L PMd | 4.4 mm | 4.4 mm | 4.4 mm | 4.4 mm |
| R PMd | 13.9 mm | 13.2 mm | 17.0 mm | 11.7 mm |
| L SMA | 8.1 mm | 5.2 mm | 14.9 mm | 5.2 mm |

A perfect head model does not consistently improve accuracy; in some cells it is worse (LCMV, left M1). The
remaining ~1 cm error reflects 32-channel spatial sampling, the ~6 mm inverse grid and the patch extent, not
the forward model. Single cells vary with the seed and should not be over-interpreted.

![Matched-head-model control](figures/fig7_matched_control.png)

### 5. A deep source always appears on the cortex

A thalamic source (MNI −12, −18, 8) is mislocalized by 25–29 mm with every method. dSPM assigns it to
**SMA in 50 %**, M1 in 30 % and occipital cortex in 20 % of recordings. A "motor cortex" peak in EEG does not
rule out a deep generator.

### 6. Noise-only recordings still produce a peak

With no bursts at all, every pipeline still produces a maximum. dSPM and LCMV place it mostly in frontal and
central midline cortex, overlapping the motor regions of interest. Maps need statistics against a null (e.g.,
phase-randomised surrogates or permutation of burst/non-burst labels), not just a peak location.

![Noise-only maxima](figures/fig5_noise_only.png)

### 7. OFF vs ON: differences from the same source

Fifteen simulated patients each have a fixed left-M1 source (jittered by 6 mm between patients). OFF and
ON differ only in burst rate, duration and amplitude, at −5 dB. Every pair is analysed twice: with a
threshold set per recording, and with one GED filter and threshold shared by OFF and ON.

| Same true source in both conditions | Threshold per recording | Shared threshold |
|---|---|---|
| True burst rate, OFF vs ON | 0.76 vs 0.66 /s (p = 0.0006) | same |
| Detected burst rate | 0.91 vs 0.88 /s (p = 0.083) | 0.93 vs 0.87 /s (p = 0.012) |
| dSPM peak error, OFF − ON | −18.7 mm (p = 0.087) | 1.9 mm (p = 0.068) |
| dSPM dispersion, OFF − ON | **−7.4 mm (p = 0.0003)** | −0.8 mm (p = 0.45) |
| dSPM dispersion with true events, OFF − ON | −3.7 mm (p = 0.003) | same |

Wilcoxon signed-rank tests, n = 15.

- With a per-recording threshold, a **"significant" difference in map spread appears with no change in
  location**. A shared threshold removes it and recovers the true rate difference.
- **Detected rates are inflated by background beta.** At −5 dB about three quarters of GED-detected events
  are not true bursts (precision 0.24 in the sweep), so the detected rate difference is much smaller than
  the true one.
- **Dispersion still differs with perfect detection**, because weaker ON bursts give broader maps. Map
  spread can only be compared between conditions if burst SNR and the number of bursts averaged are matched.
- At −5 dB, dSPM peak error per patient is bimodal (some near 10–20 mm, others 60–120 mm): at realistic SNR
  only part of a cohort is localized well.

![OFF vs ON with the same true source](figures/fig6_group_offon.png)

## Real-data analysis (`real_data/`)

The `real_data/` folder applies the validated pipeline to three open resting-state EEG cohorts from OpenNeuro.
It holds every analysis reported in the paper *Simulation-based validation of scalp EEG beta-burst imaging in Parkinson's disease*.

### Data

Download the three datasets from OpenNeuro (BIDS format):

| Cohort | OpenNeuro | Version | Design |
|---|---|---|---|
| San Diego | [ds002778](https://openneuro.org/datasets/ds002778) | 1.0.5 | 15 PD (OFF/ON) + 16 HC, 32-ch BioSemi |
| New Mexico | [ds003490](https://openneuro.org/datasets/ds003490) | 1.1.0 | 25 PD (OFF/ON) + 25 HC, 64-ch BrainVision |
| Iowa | [ds004584](https://openneuro.org/datasets/ds004584) | 1.0.0 | 100 PD (ON) + 49 HC, 64-ch BrainVision |

Then set the three `"root"` paths in the `DATASETS` dictionary in the settings cell of `pd_beta_three_datasets.ipynb` to where you saved them.

### Notebooks and run order

Run the notebooks from inside `real_data/`, in this order. All outputs are written to `real_data/pd_beta_three/`. Every expensive step is cached, so a run that stops can simply be restarted.

| # | Notebook | What it does | Paper | Runtime* |
|---|---|---|---|---|
| 1 | `pd_beta_three_datasets.ipynb` | Preprocessing (filtering, bad channels, ICA/ICLabel), GED detector with quality gate, common aperiodic-normalized threshold, burst metrics, dSPM/eLORETA/LCMV maps, surrogate tests, split-half reliability, SNR calibration, group statistics and meta-analysis, figures | Sections 2.6–2.9, 3.2–3.6; Tables 2, 10–12; Figs 8–14 | 3–4 h |
| 2 | `iowa_confound_check.ipynb` | Iowa covariate models M0–M3 (ICA cleaning, aperiodic slope, recording length) | Table 13 | < 1 min |
| 3 | `revision_analyses.ipynb` | Sensitivity analyses: software versions and source-space overlap; simulation with the real-data detector; region-specific SNR calibration; model sensitivity (patch size, frequency, two generators, STN, EMG); quality-gate sensitivity; amplitude decomposition; extended split-half reliability; 200/100 surrogates; random-effects table, bootstrap CIs, minimum detectable effects, TOST | Sections 2.10, 3.1.7–3.1.8; Tables 8, S1–S6 | 4–5 h |
| 4 | `harmonic_check.ipynb` | Why the real-data detector fails at high SNR: the second harmonic of non-sinusoidal bursts biases the gate's aperiodic fit | Table 9 | ~20 min |
| 5 | `waveform_check.ipynb` | Burst amplitude under four aperiodic models (incl. FOOOF), gate re-run with a harmonic-safe fit, waveform shape (sharpness, rise–decay asymmetry, harmonic excess) | Section 3.5.4; Table 14 | ~10 min |

\*On a desktop PC (Intel i7, 64 GB RAM). Most of the time in notebook 3 goes to the map surrogates. Set `N_SURR_MAPS = 0` to skip them.

Notebooks 3–5 do not repeat the main analysis. They import the functions of the notebook before them and read its cached outputs. Notebook 1 loads the simulation framework from `../beta_burst_sim.ipynb`.

### Outputs

```
real_data/pd_beta_three/
    preprocessed/   cleaned recordings + preprocessing_qc.csv   (not in the repository; regenerate with notebook 1)
    cache/          head models and calibration simulations     (not in the repository)
    results/        per-recording tables, statistics, meta-analysis     → copied to real_data/results/
    revision/       R0–R7, H1 and W1 tables of the sensitivity analyses → copied to real_data/results/revision/
    figures/        all real-data figures                               → copied to real_data/figures/
```

The CSV files in `real_data/results/` and the figures in `real_data/figures/` are the exact outputs behind the numbers in the paper.

### Additional requirements

The real-data notebooks need `mne-icalabel`, `onnxruntime` and `python-picard` (preprocessing), and `fooof` (notebook 5). All are listed in `requirements.txt`.

## Practical recommendations for real data

1. Detect bursts on **one** spatial signal (a GED component or a source time course), never as "any channel
   above a pooled threshold".
2. Align events to the burst trough, set `nave` to the number of bursts, and estimate the noise covariance
   from burst-free data.
3. For condition comparisons, use one filter and threshold per subject shared across conditions, and match
   the number of bursts and burst amplitude before comparing map dispersion.
4. With 32 channels and a template head, M1 can be identified at adequate SNR, but do not separate premotor
   cortex from SMA with dSPM alone (check with sLORETA/eLORETA), and do not rule out deep sources.
5. Test maps against a noise-only or label-permuted null.
6. Run this framework with your own montage, SNR and burst statistics before interpreting a real map.

## Limitations

- **One template head.** Individual anatomy, a denser montage or MEG would change the numbers, especially
  for PMd vs SMA. Accuracy also depends on electrode positions matching the head: the version 1 template
  offset alone created a systematic M1 → PMd error.
- **Simulated bursts and background.** Real bursts vary in waveform, spatial extent and orientation; real
  background contains real bursts. `REAL_BACKGROUND_FILE` allows a real background, but then its ground
  truth is impure.
- **Single-source recordings.** Each recording has one burst source; simultaneous bilateral or
  cortico-subcortical bursts are not simulated.
- **Location-dependent results.** Errors differ between hemispheres and between neighbouring locations
  because of local folding, so the tables describe these seeds, not the motor cortex in general.
- **Coarse region score.** Region accuracy assigns the maximum to the nearest of six seeds; it measures
  separability among these regions only.

## Reproducing the figures

The notebook writes all figures and tables to `burst_sim_results/`. The files in `figures/` and `results/`
are copies from the reference run.

| Figure | Notebook output |
|---|---|
| Example recording and maps | `fig1_example.png` |
| Error vs SNR per location | `fig2_error_vs_snr.png` |
| Region confusion (dSPM, LCMV) | `fig3_confusion_dSPM.png`, `fig3b_confusion_LCMV.png` |
| Detection accuracy | `fig4_detection.png` |
| Noise-only maxima | `fig5_noise_only.png` |
| OFF vs ON | `fig6_group_offon.png` |
| Matched-head-model control | `fig7_matched_control.png` |

## Citation

If you use this code, please cite this repository:

```bibtex
@software{hasib_beta_burst_sim_2026,
  author = {Hasib, Tanvir},
  title  = {Ground-truth simulation of {EEG} beta-burst source localization},
  year   = {2026},
  url    = {https://github.com/tanvir-hasib/beta-burst-localization-simulation-model}
}
```

Key method references: Dale et al. (2000), *Neuron* 26:55–67 (dSPM) · Pascual-Marqui (2002), *Methods Find.
Exp. Clin. Pharmacol.* 24D:5–12 (sLORETA) · Pascual-Marqui (2007), arXiv:0710.3341 (eLORETA) · Van Veen et
al. (1997), *IEEE Trans. Biomed. Eng.* 44:867–880 (LCMV) · Gramfort et al. (2013), *Front. Neurosci.* 7:267
(MNE-Python) · Cohen (2022), *NeuroImage* 247:118809 (GED) · Shin et al. (2017), *eLife* 6:e29086 (beta
events).

## License

MIT — see [`LICENSE`](LICENSE).

## Contact

Tanvir Hasib · Multimedia University, Malaysia · open an [issue](../../issues) for questions or bug reports.
