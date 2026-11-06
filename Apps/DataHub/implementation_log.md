# Implementation Log

## [2026-11-06] — AreaComparison: FEF vs SC timing (which area's value signal rises first)
**Files changed:** `Batch_ObjLocValue_AreaComparison.m` (timing is now the main analysis; overlap kept behind `Params.RunOverlap`), `ObjLocSingleNeuron/` (new: ComputeValueLatency, ComputePopulationOnset, ComputeAreaLatencyData, ComputeLatencyAreaComparison, PlotLatencySlide).
**What:**
- Single-neuron onset: rates use a causal 30 ms boxcar, so no backward leak. The onset is the first bin starting ≥5 consecutive 10 ms bins with high-vs-low Mann-Whitney p < .05 and the same sign as the 100–200 ms AUC.
- Population onset: half-max of the sign-aligned AUC−0.5, with a 2000× bootstrap over neurons.
- Variants: main = each trial's response minus its own −300 to 0 ms baseline, on all value-coding neurons; control = raw responses, only neurons without anticipation.
- Reference: visual-response onset of the same neurons (signrank of rate minus baseline), and the per-neuron value−visual lag.
- FEF vs SC tests: bootstrap onset difference and ranksum, with Holm correction within each signal family.
- Outputs: `Timing_*_<Subgroup>.png/.fig` and `Timing_<Subgroup>.csv` in `Results/ObjectLocationPSTH/AreaComparison/`.

**Result (All):**
- Visual onset: SC leads by about 27 ms in both tasks (population 50 vs 76 ms; single-neuron medians 45 vs 85–105 ms; p < .002).
- Value onset (corrected): SC is earlier at the population level, object −13 ms (Holm .34) and location −43 ms (Holm .20, wide SC CI 42–91 ms). Single-neuron medians are equal. Not significant.
- Value−visual lag: SC is longer (object +15 ms, Holm .038). The SC visual lead is not carried over into an earlier value signal.

**Why:** The user's FEF vs SC comparison focuses on timing.

---

## [2026-11-06] — Batch_ObjLocValue_AreaComparison: FEF vs SC comparison from the DataHub GUI
**Files changed:** `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_AreaComparison.m` (new), `ObjLocSingleNeuron/ComputeBatchRowArea.m` (new), `ObjLocSingleNeuron/PlotSlideScatterVenn.m` (Venn circle without duplicate vertex: removes the polyshape warning), both switchyards.
**What:** A GUI batch protocol. Load `FEFSC_BothRA_StableValue.txt` and each row is assigned to FEF or SC by session plus channel range. The match uses both because simultaneous sessions (Adams042624/060624FRACSAC1) are listed in both areas with different channels. If only one area is loaded, the other area's batch file is loaded automatically. For All / Fast / Others it computes overlap vs chance (all neurons, with and without anticipation, baseline-corrected), coding fractions, and FEF vs SC tests (Fisher, log-OR z, Holm). Figures are made for `PlotSubgroups` (default All): scatter+Venn per area, overlap vs chance, anticipation split, baseline-corrected, area proportions, area overlap OR. Outputs go to `Results/ObjectLocationPSTH/AreaComparison/`. Results for All match ObjLocValue_SlideDeck exactly (FEF 302, SC 215).
**Why:** The per-area batch scripts take one batch file each, so the GUI had no way to compare FEF and SC.

---

## [2026-11-05] — ObjLocValue_SlideDeck: talk-ordered figures (scatter+Venn → heatmap → 4-type AUC → overlap vs chance → anticipation split → baseline-corrected → FEF vs SC)
**Files changed:** `AnalysisProtocol/BatchAnalysis/ObjLocValue_SlideDeck.m` (new), `ObjLocSingleNeuron/` (new: ComputeOverlapVsChance, ComputeSlideDeckData, ComputeSlideDeckComparison, PlotSlideScatterVenn, PlotSlideOverlapVsChance; PlotPresentationEffectDots gains `LegendPos`).
**What:** Overlap vs chance: E[both] = nObj·nLoc/N with a hypergeometric 95% range, and a Fisher test. The verdict follows the Fisher p. Slides 1–5 use the protocols' saved AUCs; the anticipation split (oneDR pre-target ranksum) and slide 6 use trial-level ROC. FEF vs SC tests are Holm-corrected over 7 tests. Headlines are generated from results. Outputs go to `Results/ObjectLocationPSTH/SlideDeck/`.
**Why:** Presentation in the user's requested order.

---

## [2026-11-05] — ObjLocValue_SingleUnitSummary: cross-area summary, FEF vs SC tests, slide figures
**Files changed:** `AnalysisProtocol/BatchAnalysis/ObjLocValue_SingleUnitSummary.m` (new), `ObjLocSingleNeuron/` (new: BuildObjLocBatchData, ComputeSplitHalfReliability, ComputeSingleUnitSummary, ComputeAreaComparison, ComputePresentationData, FormatSingleUnitSummaryText, PresentationStyle, Slide* helpers, PlotPresentationFractionBars, PlotPresentationEffectDots, PlotPresentationAUCScatter).
**What:** Runs FEF and SC together and collects every single-unit metric: coding fractions, pre-target differences, overlap OR, Spearman rho, concordance, and a split-half noise-corrected r (Spearman-Brown, 100 splits, bootstrap CI over neurons). Formal FEF vs SC tests use Fisher, ranksum, log-OR z, Fisher z and bootstrap, Holm-corrected over 16 tests. Produces 5 slide figures whose headlines are worded from the results. Area colors are CVD-validated (#2a78d6 / #eb6834). Outputs go to `Results/ObjectLocationPSTH/SingleUnitSummary/`.
**Why:** Gives one reproducible summary of the single-unit observations and the area comparison for the presentation.

---

## [2026-11-05] — Baseline-corrected ROC + split by oneDR baseline difference
**Files changed:** `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_BaselineCorrectedROC.m` (new), `Batch_ObjLocValue_BaselineSplit.m` (new), `ObjLocSingleNeuron/` (new: MannWhitneyAUC, ComputeTrialWindowAUC, ComputeAUCVariantTable, ComputeOverlapStats, ComputeBaselineCorrectionResults, ComputeBaselineSplitResults, PlotRawVsCorrectedAUC, PlotOverlapStatsComparison, PlotObjLocAUCScatterBySubset; LoadObjLocNeuronTable gains `LoadTrials`), both switchyards.
**What:** Trial-level ROC from `psth_targeton_raw` (RF trials): raw response [100 200] ms, response minus the same trial's [-300 0] ms baseline, and baseline alone. Uses exact Mann-Whitney AUC + ranksum. Raw AUC vs saved TgtAUC_RF: Spearman rho ≈ 0.97 (saved uses smoothed FR). Script 1 compares object vs location coding raw vs corrected. Script 2 splits neurons by the oneDR baseline difference (p < .05) and compares the object–location relationship in each subset. Time courses use a 50 ms boxcar (`TCSmoothMs`).
**Why:** The oneDR block design creates a pre-target value offset. These scripts test whether the object–location relationship depends on it.

---

## [2026-11-05] — Trial-count check applied to all AUC analyses; location-sorted heatmap; heatmap range
**Files changed:** `Batch_ObjLocValue_Organize.m`, `Batch_ObjLocValue_SingleNeuron.m`, `ObjLocSingleNeuron/{LoadObjLocNeuronTable,ComputeObjLocSelection,ComputeObjLocSingleNeuronResults,PlotAUCHeatmapPair}.m`
**What:** `Sel_All` now requires >= 10 trials in each RF condition (high/low x stable/oneDR), the check that was previously applied only to the saved PopPSTH (`Sel_Raw`). The saved PopPSTH neuron set is unchanged. The new script uses the same rule (`Params.MinTrials`). Neuron sets match Organize exactly: FEF 302/45/204, SC 215/33/139 (All/Fast/Others). Also added a heatmap sorted by location AUC, set the heatmap x-range to [-200 400] ms, and both scripts now add the helper folder to the path themselves (fixes GUI runs).
**Why:** Keeps the AUC figures and the decoder input on the same neurons, and drops noisy low-trial AUCs.

---

## [2026-11-05] — Batch_ObjLocValue_SingleNeuron: object vs location value at the single-neuron level
**Files changed:** `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_SingleNeuron.m` (new), `AnalysisProtocol/BatchAnalysis/ObjLocSingleNeuron/*` (new: LoadObjLocNeuronTable, ComputeObjLocSelection, ComputeValueGroups, ComputeSignAlignedAUC, ComputeGroupAUCTimecourse, ComputeSignConcordance, ComputePreTargetAUC, ComputeOtherValueEquivalence, ComputeObjLocSingleNeuronResults, ComputeAUCSortOrder, Plot*), `Batch_ObjLocValue_Organize.m` (heatmap helpers moved out to shared files), `RunBatchObjLoc.m` (BatchFunc / flags / SingleNeuronParams overridable), both switchyards.
**What:** Groups neurons as ObjOnly/LocOnly/Both/None/Ambiguous/Unclassified using the 100–200 ms ranksum p. The "no coding" rule is an option (`NoCodingCriterion` = pvalue | effect | both). Diagnostics: sign concordance (binomial), pre-target AUC (own-sign aligned), other-value equivalence (sign-aligned mean, bootstrap CI, TOST), sign-aligned group AUC time courses, and AUC heatmaps sorted by object AUC. FEF/SC are run separately, each as All/Fast/Others.
**Why:** To test whether the "other" value AUC in ObjOnly/LocOnly groups is real sub-threshold coding or a confound.

---

## [2026-11-05] — Fix OneDR NaN fallback + add object-sorted AUC heatmaps
**Files changed:** `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Organize.m`
**What:** (1) Removed NaN assignments in the OneDR `else` branch that overwrote stable-task variables (`p_Vis_RF`, `VisDiff_RF`, `p_pureSac_RF`, `hPSTH_TgtOn_Value_RF`, …); the `_oneDR` fallbacks below are unchanged. (2) Added side-by-side AUC-vs-time heatmaps (object | location) for Sel_All, both sorted by object-task TgtAUC_RF (`HeatmapSortMode` = 'signed' or 'abs'). New local helpers: `ComputeAUCSortOrder` (pure) and `PlotAUCHeatmapPair` (render only).
**Why:** The bug could corrupt stable-task fields for neurons without OneDR data. The heatmap shows whether location value coding follows object value coding neuron by neuron.

---

## [2026-11-05] — Gaze-duration histogram range fixed to 350-1500 ms

**Files changed:** `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Behavior.m`
**What:** New local `MetricDisplayRange` ([350 1500] for GazeDur, [] otherwise). In `PlotFeatureHistogram` (fig 205, all four monkey × task panels) it sets the bin edges (350:25:1500) and the x-limits. Other figures keep automatic ranges.
**Why:** User request; consistent range across the four histogram panels.

---

## [2026-11-05] — Fig 208 back to error rate (selectable)

**Files changed:** `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Behavior.m`
**What:** New setting `FailureMeasure` at the top of the batch script ('Error' default, or 'BreakFix'). It is passed to `ComputeBatchValueBehavior`, which reads `Summary.<Measure>Rate_High/_Low` and `<Measure>N_*` generically; fig 208 labels and the console summary follow it. Both measures are still exported by ValueTaskPerformance; the per-session figure still shows the break-fixation rate.
**Why:** User: the error-rate difference is more visible.

---

## [2026-11-05] — Fix: fixation breaks in the object task were never classified as breaks

**Files changed:** `AnalysisProtocol/ComputeValueTaskPerformance.m`
**What:** Error stage is now taken from the outcome codes and their timing: DelayBreak = BROKEFIX logged before any FIXOFF; NoSaccade = FAILED_SAC (or FIXOFF without SAC_OCCURRED); otherwise GazeBreak. Previously DelayBreak meant "no FIXOFF", but in the object task the task logs FIXOFF as cleanup after a break, so every object-task break was classified as NoSaccade/GazeBreak and the break-fixation rate was 0 in all sessions (the user's batch: 0 in 35 Robin and 23 Adams sessions). OneDR trials are unchanged (no cleanup FIXOFF).
Check: Robin062923FRACSAC1 breaks H/L 5/7 (was 0/0); Adams102125FRACSAC1 25/44, p = 0.036; OneDR sessions identical to before.
Also found: all 58 OneDR_Performance exports in Results predate the BreakFix field, so the location panels of fig 208 were empty (n = 0). ValueTaskPerformance must be re-run on the OneDR batch.
**Why:** User saw zero object-task and empty location-task panels in fig 208.

---

## [2026-11-05] — Break-fixation / error rate only over trials with the target shown; recover uncoded values

**Files changed:** `AnalysisProtocol/ComputeValueTaskPerformance.m`, `AnalysisProtocol/ValueTaskPerformance.m`, `AnalysisProtocol/PlotValueTaskPerformance.m`
**What:** Rates now use only All Data trials with TGTCD (target shown) and a known value. The GOODOBJ/BADOBJ code is sent at a task-dependent time (object: usually 50–320 ms before target onset; OneDR: 200–400 ms after), so in OneDR a break between target onset and the code left the trial uncoded and it was dropped. New local `InferMissingValue` recovers the value:
- Location tasks: target side (LEFTHEMI 80 / RIGHTHEMI 81; if missing, a target-ID code in 1000–1999 mapped to side when that map is ≥95% consistent) compared with the high-value side of the neighbouring coded trials; left unknown if they disagree (a reversal in between).
- Object task: target-ID code → value from coded trials, validated leave-one-out.
- A rule is applied only if it reproduces the coded trials' values with accuracy ≥ `MinInferenceAccuracy` (0.95).

New outputs: TrialTable `ValueInferred`, `Perf.NumValueInferred`, `Perf.ValueInference` (Method / Accuracy / NInferred), and `Summary.NumValueInferred`. The per-session figure note reports the recovered and unknown counts.
Tests: Robin110221OneDR4 recovered 26 trials (validation 100%), Robin021324OneDR1 recovered 9 (100%); Adams and object sessions needed none. The object ID rule failed validation in Adams102125 (59%) and is not applied there.
**Why:** User request — break-fixation rate only among trials after the target was shown.

---

## [2026-11-05] — Break-fixation rate replaces error rate; full-colour bars in fig 208

**Files changed:** `AnalysisProtocol/ComputeValueTaskPerformance.m`, `AnalysisProtocol/ValueTaskPerformance.m`, `AnalysisProtocol/PlotValueTaskPerformance.m`, `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Behavior.m`
**What:**
- New `Perf.BreakFix`: fraction of value-coded trials (same denominator as the error rate) with fixation broken before the go signal (ErrorStage 'DelayBreak'); Fisher test high vs low. Exported as `Summary.BreakFixRate_High/_Low`, `BreakFixN_*`, `BreakFix_P`. The error rate is still computed and exported.
- Per-session figure: the left panel is now the break-fixation rate (value colours) instead of the stacked error stages.
- Batch fig 208 and the console summary use the break-fixation rate. Sessions exported before this change (no `BreakFixRate_*`) are left out of that panel until ValueTaskPerformance is re-run.
- Fig 208 bars are full colour, with white session dots (black edge) on top, dark grey pairing lines, and the sem error bars drawn last.
**Why:** User preference: break-fixation rate isolates the value-dependent failure (holding fixation during the delay) better than the overall error rate.

---

## [2026-11-05] — Value colour convention: object red/blue, location purple/green

**Files changed:** `ToolBox/ValueColors.m` (new), `AnalysisProtocol/OneDR/OneDRAnalysis.m`, `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Behavior.m`, `AnalysisProtocol/PlotValueSaccadeBehavior.m`, `AnalysisProtocol/PlotValueTaskPerformance.m`, `AnalysisProtocol/ComputeValueSaccadeBehavior.m`
**What:** `ValueColors(Context)` returns [high low] RGB: object = red [0.85 0.15 0.15] / blue [0.15 0.30 0.85]; location = purple [0.49 0.18 0.56] / green [0.13 0.69 0.47]. It accepts 'Obj'/'Object'/'Loc'/'Location' or a task code (ONEDR_TASK, OBJ_LOCATIONVALUE_TASK → location).
- OneDRAnalysis: all high/low PSTH and raster colours are now purple/green (6+6 shadedErrorBar, 4 SpikeRaterPlot, 4 trial-order raster plots, annotation text). Excluded switch trials stay grey.
- Batch figs 201 and 203–208: colours are set per task panel.
- Per-session plots: colours chosen from `Behavior.TaskCode` (new field, set in ComputeValueSaccadeBehavior).
- PlotValueTaskPerformance: the GazeBreak stage colour changed from purple to brown so it can't be confused with location high value; the "no value code" note moved to the bottom.
**Why:** Match the lab figure convention (object value red/blue, location value purple/green).

---

## [2026-11-05] — OneDRAnalysis: target eccentricity/position for the 3-D layout

**Files changed:** `AnalysisProtocol/OneDR/OneDRAnalysis.m`
**What:** In the 3-D layout (Adams: trials × [ID target1 target2] × [x y] or [angle ecc]), `Data.TargetPolarLoc(Selection,2:end)` mixed target 2 into the columns, so `TargetEcc` was target 2's angle (e.g. 70/250 instead of 20 deg) and `UniqueEcc` had 2 values. Target 1 is now read explicitly as `[TPL(:,2,1), TPL(:,2,2)]`, and the same for `TargetLoc`. The 2-D layout (Robin) is unchanged. The main-figure titles now print `Angle: 70, Ecc: 20.0` instead of `7.0e+01`.
**Why:** Wrong eccentricity in the titles and in the layout logic that depends on `length(UniqueEcc)`. Robin's output still matches the baseline, apart from the stochastic permutation p-values.

---

## [2026-11-05] — OneDRAnalysis: reversal detection fix + excluded switch trials shown in grey; fig 208 bars + significance

**Files changed:** `AnalysisProtocol/OneDR/OneDRAnalysis.m`, `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Behavior.m`
**What:**
- OneDRAnalysis: value switches are now detected as a change in the high-value side, not `diff(TargetIndex(:,2))`, which found no switches in Adams's 3-D layout, so his first 5 post-switch trials were never excluded. Analyses still run only on kept trials (exclusion unchanged).
- Before the exclusion, all-trial copies (`AllTrials`) are saved. With figures on, the local function `BuildAllTrialRasters` builds target- and saccade-aligned rasters for all trials using the same alignment as the main pipeline. The trial-order rasters in fig 110 draw every trial: red = high, blue = low, grey = excluded switch trial. An annotation gives the excluded count. The exclusion message now prints correctly (the old `disp` concatenated a number into a char).
- Regression check: Robin110221OneDR4 (same switches under both rules) gives OutputData identical to the pre-change baseline, except the two permutation-test p-values (`p_AUC_Tgt_*`), which vary run to run.
- Fig 208 (`PlotPerformance`): mean ± sem bars, jittered session dots joined by grey lines, and a significance bracket with stars from the paired signrank (`SigStars`).
**Why:** Apply the switch exclusion correctly to both monkeys and make the excluded trials visible; make the error/anticipation figure presentation-ready.

---

## [2026-11-05] — Batch error/anticipation summary; Adams TargetLoc fix; reward-based gaze check

**Files changed:** `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Behavior.m`, `AnalysisProtocol/ComputeValueSaccadeBehavior.m`
**What:**
- Batch now also reads `StableValuePerformance` / `OneDR_Performance` (optional per session). Fig 208, `PlotPerformance`: rows = monkey; columns = error rate (object, location) and anticipation rate (object, location); session lines high→low, mean ± sem, paired signrank across sessions, and Fisher on summed counts. Also prints a console summary.
- Fix: in Adams's 3-D `TargetLoc` layout (trials × [ID target1 target2] × [x y]), the target centre is `[TargetLoc(:,2,1), TargetLoc(:,2,2)]`. It was read as 2-D, which put the gaze window at the wrong place, so every Adams OneDR trial was LandedOutside and gaze duration was all NaN.
- `GazeTooShort` now means the eye left the target before reward onset (RWDCD), replacing the fixed 400 ms floor. The task's hold timer starts at window entry, mid-saccade, so valid gaze measured from saccade end can be slightly under 400 ms.
- Default `TargetWindowRadius` raised from 5 to 6 deg (Adams's landing error reaches about 5 deg).
**Why:** The user asked for batch summaries of error rate and anticipation; the Adams format bug and the gaze check showed up while testing them.

---

## [2026-11-05] — ValueTaskPerformance: per-session error rate + anticipatory saccades; reversal detection fix

**Files changed:** `AnalysisProtocol/ValueTaskPerformance.m` (new), `AnalysisProtocol/ComputeValueTaskPerformance.m` (new), `AnalysisProtocol/PlotValueTaskPerformance.m` (new), `AnalysisProtocol/SaveProtocolEntry.m` (new), `AnalysisProtocol/SaveValueBehaviorOutput.m`, `AnalysisProtocol/ComputeValueSaccadeBehavior.m`, `DataAnalysisListSwitchyard.m`, `DataAnalysisProtocolSwitchyard.m`
**What:** New protocol `ValueTaskPerformance`, registered for OBJ_DIR_SAC_TASK, ONEDR_TASK and OBJ_LOCATIONVALUE_TASK. It loads Good and All Data from the packed file (`Data.Path`) and maps the Good-trial range onto All Data (Good = All trials with CORRECTCD, in order).
- Error rate = fraction of value-coded All trials not completed, split by stage (DelayBreak / NoSaccade / GazeBreak); trials where the target was shown but the value is unknown are counted separately.
- Anticipation = fraction of completed trials with RT < 80 ms (also RT < 0).
- High vs low compared with a Fisher exact test.
- Exports `OutputData.{StableValue|OneDR|ObjLocationValue|FlexibleValue}Performance` with `.Summary` (rates, n, p) and `.Perf` (`TrialTable`, by-stage counts).
- New `SaveProtocolEntry` handles the per-channel file writes; `SaveValueBehaviorOutput` now uses it.
- Reversal fix: `ComputeValueSaccadeBehavior` now detects OneDR reversals by a change in the high-value side, not by `diff(TargetIndex(:,2))`. In Adams's 3-D TargetIndex layout that column is all zeros, so no reversals were found and post-reversal trials were never excluded. The new method matches TargetIndex exactly for Robin.
**Why:** RT does not differ between high and low value in this delayed-saccade design; error rate and anticipation do.

---

## [2026-11-05] — Batch_ObjLocValue_Behavior: inset bar plots in histogram figures

**Files changed:** `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Behavior.m`
**What:** Each panel of the histogram figures (203 RT, 204 PeakV, 205 GazeDur) now has an inset bar plot (top right) of high vs low mean ± sem, computed from the same trials. The main axis gets 60% headroom with no ticks in it; the legend moved to the northwest. Insets are drawn after `sgtitle` so panel positions are final. The violin plot (207) is unchanged.
**Why:** Show the mean difference alongside each distribution.

---

## [2026-11-05] — Batch_ObjLocValue_Behavior: violin plot of saccade parameters

**Files changed:** `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Behavior.m`
**What:** Fig 207, `PlotViolin`: rows = monkey, columns = RT / PeakV / GazeDur. Each panel shows four violins (object high/low, location high/low) of trials pooled across sessions, drawn with the built-in `violinplot` (R2024b+), plus a median dot, an IQR bar, and rank-sum p for high vs low in each task.
**Why:** Show the shape of the high vs low value difference for each saccade parameter in both value tasks.

---

## [2026-11-05] — JLab-style target saccade detection + gaze hold check

**Files changed:** `AnalysisProtocol/DetectTargetSaccade.m` (new), `AnalysisProtocol/ComputeValueSaccadeBehavior.m`
**What:** Replaced SelectSaccade-based detection with `DetectTargetSaccade`, ported from JLab `CalculateRT/detectSaccade`: Savitzky-Golay (7-sample) smoothed 2-D speed; noise >3000 deg/s; baseline speed over [−250, −100] ms re: go; trigger = baseline + 30 deg/s held 5 samples. Onset is searched back from the `SAC_OCCURRED` marker (eye left the fixation window, used in place of JLab's Fixation_exited) to −100 ms, with JLab's strict→released criteria, then back-traced. Offset comes from the speed trace within 200 ms (raising the threshold if needed); JLab's Choicetime anchor was dropped because there is no saccade-off marker. Differences from JLab: peak velocity is taken over onset→offset rather than from the trigger (JLab missed the peak after anticipatory saccades), and the baseline window ends at −100 ms because some onsets are as early as −67 ms. Gaze shorter than the task's required hold (ParaLib `GazeHoldForTarget`, 400 ms) is flagged `GazeTooShort` and set to NaN.
**Why:** To follow the lab's validated saccade detection. Rewarded trials need ≥400 ms on target, so a shorter gaze duration is a measurement error.

---

## [2026-11-05] — Batch_ObjLocValue_Behavior: histograms for all three metrics; per-channel save fix

**Files changed:** `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Behavior.m`, `AnalysisProtocol/SaveValueBehaviorOutput.m`, `AnalysisProtocol/ObjDirectedSaccadeTask/RTGoodBadAnalysis.m`, `AnalysisProtocol/OneDR/OneDRRT.m`
**What:** Replaced `PlotRTHistogram` with `PlotFeatureHistogram`, which plots high vs low value distributions for RT (fig 203), PeakV (fig 204) and GazeDur (fig 205); the 3-D plot moved to fig 206. RT keeps fixed bins (−150:20:500 ms); PeakV/GazeDur use 25-unit bins over each panel's range. Rows = monkey, columns = task. SaveValueBehaviorOutput now writes one file per channel in `Data.SpikeChannel`.
**Why:** The user asked for histograms of all three metrics. Bug: GUI batch mode passes a channel range (e.g. 4:46), and `sprintf` on the vector produced junk filenames (`…_C4.mat5_6_N7…`). When a junk name already existed from April, `load` errored and the GUI batch stopped silently at row 20 of FEFSC_BothRA_OneDR.

---

## [2026-11-05] — Batch_ObjLocValue_Behavior: RT histograms + 3-D feature space per monkey

**Files changed:** `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Behavior.m`
**What:** Compute now builds `Stats.Trials.(Obj|Loc)`, a pooled trial table (Monkey, Session, IsGood, RT, PeakV, GazeDur; trials with any NaN dropped). Monkey comes from the session name ('A…' = Adams, otherwise Robin). Fig 203: RT histograms of high vs low value, 20-ms bins, dashed lines at means, rank-sum p; rows = monkey, columns = task. Fig 204: 3-D scatter of trials in (RT, PeakV, GazeDur), with the high and low means joined by a line; rows = monkey, columns = task.
**Why:** Show trial-level value separation in each task, separately for each monkey.

---

## [2026-11-05] — Batch_ObjLocValue_Behavior: population high vs low value behavior, object vs location

**Files changed:** `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Behavior.m` (new), `DataAnalysisListSwitchyard.m`, `DataAnalysisProtocolSwitchyard.m`
**What:** Batch protocol, run on a stable value batch file. Reads `StableValueBehavior` and `OneDR_Behavior` (TrialTable/SessionMean format) from the per-neuron result files. Unit = session: the first channel file of each batch row that has both tasks exported; rows without both are skipped and listed. Stats: per-session high/low means; paired signrank across sessions for each task × metric; object-vs-location value effect (High−Low) compared with signrank and Spearman correlation. Fig 201: 2×3 paired session plots (filled marker = within-session p<0.05). Fig 202: 1×3 scatter of the object vs location value effect. Saves `Sessions` and `Stats` to `Results/ObjLocBehavior/<Area>/ObjLocBehavior_<BatchFile>.mat`.
**Why:** Population-level behavioral validation that high vs low value changes saccades in both value contexts (ProjectOverview, behavior quantification).

---

## [2026-11-05] — Value behavior: target-window gaze rule + per-trial/session export

**Files changed:** `AnalysisProtocol/ComputeValueSaccadeBehavior.m`, `AnalysisProtocol/PlotValueSaccadeBehavior.m`, `AnalysisProtocol/SaveValueBehaviorOutput.m`
**What:** GazeDur now runs from the end of the saccade until the eye leaves the target window. The window is a circle of `TargetWindowRadius` = 5° around `Data.TargetLoc`. Leaving means being outside the window, or missing samples, for ≥20 ms. The eye must enter the window within 100 ms of saccade end; otherwise the trial is flagged `LandedOutside`. Removed GazePostTgtOff. The export now also includes `TrialTable` (per-trial RT/PeakV/GazeDur plus flags) and `SessionMean` (All/Good/Bad mean, sem, n). New DataStamp keys: `RT`/`PeakV`/`GazeDur` for all trials, `_Good`/`_Bad` splits, and `_Mean`/`_Good_Mean`/`_Bad_Mean`; the legacy `RT_*` keys are kept.
**Why:** To match the task's definition of gaze (eye inside the target window) and to support session-level aggregation.

---

## [2026-11-05] — Per-session value behavior: RT, peak velocity, gaze duration (object + location value)

**Files changed:** `AnalysisProtocol/ComputeValueSaccadeBehavior.m` (new), `AnalysisProtocol/PlotValueSaccadeBehavior.m` (new), `AnalysisProtocol/SaveValueBehaviorOutput.m` (new), `AnalysisProtocol/ObjDirectedSaccadeTask/RTGoodBadAnalysis.m` (rewritten), `AnalysisProtocol/OneDR/OneDRRT.m` (rewritten), `DataAnalysisProtocolSwitchyard.m` (OneDR_RT edit target)
**What:** The two duplicated per-session RT protocols now share one pure compute function, one plot function and one save helper. Metrics for high vs low value: RT (re: FIXOFF), 2-D peak velocity, gaze duration (saccade end to next large saccade/blink), and gaze leave time re: target offset. Target saccade = first saccade >3° starting >50 ms after target onset; it is rejected as an artifact if amplitude is off target eccentricity by >40%, or eye samples are missing within ±10 ms of it. OneDR drops the first 5 trials after each reversal (matches OneDRAnalysis) and also summarizes by target side. Summary per metric: mean/sem/n, rank-sum p, AUC. Saves to `OutputData.StableValueBehavior` / `OneDR_Behavior` with the legacy `RT_*` DataStamp keys kept, plus new keys and the full `Behavior` struct. OneDR now saves to channel `Data.SpikeChannel` (was hardcoded C1).
**Why:** Behavioral validation of object vs location value (ProjectOverview). Legacy code took the first saccade after FIXOFF, which missed anticipatory saccades (they start ~20–30 ms before the FIXOFF code), and it kept blink artifacts (4000 deg/s, 35°). GazeDur is partly task-imposed: the target stays on through the reward, and the high-value reward is longer (252 vs 82 ms), so GazePostTgtOff was added as a confound-free gaze measure.

---

## [2026-03-30] — PopulationAnalysis: oneDR decoder + cross-context decoders (Stable↔OneDR)

**Files changed:**
- `PopulationAnalysisSeparate/Batch_ObjLoc_PopPSTH_FEF.m`

**What:** Four additions to the FEF population decoding pipeline:
1. **OneDR preprocessing** — smoothed (Gaussian 30 ms) and trimmed to [-500, 300] ms, with per-neuron `n_min_oneDR`.
2. **OneDR within-context decoder** — identical pseudo-population structure (70/30 split, 100 repeats, 100 permutations, cluster-corrected significance, cross-temporal matrix). Saves to `Decoder_OneDR_RF_FEF_Others.mat`.
3. **Cross-decoder: Train Stable → Test OneDR** — trains on stable (70% of trials, resampled), tests on all oneDR trials; z-score from stable training stats applied to oneDR test set. Saves diagonal accuracy + cross-temporal matrix to `Decoder_Cross_Stable2OneDR_RF_FEF.mat`.
4. **Cross-decoder: Train OneDR → Test Stable** — trains on oneDR (70%), tests on all stable trials; z-score from oneDR training stats. Saves to `Decoder_Cross_OneDR2Stable_RF_FEF.mat`.

**Why:** Assess value-coding generalization across task contexts (stable value task ↔ one-direction task). The cross-decoders reveal whether the population geometry for value is shared between the two tasks.

## [2026-03-30] — Paper Figure 1: scatter tgtAUC stable vs oneDR with marginal histograms

**Files changed:**
- `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Organize.m`

**What:** Added two paper figure sections (raw and rectified AUC scatter) plus a `plotAUCScatter` local helper function. Also added `TgtAUC_oneDR_RF_rect` scalar rectification (previously missing — only the PSTH time-series version existed).

**Why:** Paper Figure 1 requires a scatter of each neuron's value-discrimination AUC in the stable task (x) vs oneDR task (y), with marginal histograms and significance-coded markers (4 groups: sig both = filled circle, sig stable only = filled square, sig oneDR only = filled triangle, neither = open circle).

## [2026-03-23 ~session] — Project setup & peakAUC bug fix

**Files changed:**
- `CLAUDE.md` (created)
- `implementation_log.md` (created)
- `AnalysisProtocol/OneDR/OneDRAnalysis.m`
- `AnalysisProtocol/ObjDirectedSaccadeTask/ObjectResponsesDifferentLocation.m`

**What:** Created CLAUDE.md with project context, scientific questions, and embedded skills. Fixed peakAUC/peakTime bug in both single-neuron analysis files: replaced `max(auc_sm)` with `max(abs(auc_sm - 0.5))` so neurons that prefer "bad" (AUC < 0.5) are correctly identified.

**Why:** `max(AUC)` only finds the highest AUC value, missing neurons whose peak discriminability lies below 0.5. The permutation test and latency logic were already correct but were receiving a wrong peakTime as input due to this bug.

---

## [2026-03-23] — Batch bug fix + population PSTH save

**Files changed:**
- `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Organize.m`

**What:** Fixed missing `(index)` on `peakTime_AUC_OneDR_RF` and `LatencyTimeAUC_OneDR_RF` in the OneDR-absent else block (was overwriting entire arrays with scalar NaN). Removed redundant duplicate assignments. Added population data assembly and save block after the loading loop, writing `psth_object_high`, `psth_object_low`, `psth_location_high`, `psth_location_low`, `auc_psth_object`, `auc_psth_location`, `time_object`, `time_location` to `Results/PopData_ObjLoc_<BatchFileName>.mat`.

**Why:** Bug caused `peakTime`/`LatencyTime` arrays to be corrupted (reset to scalar) whenever a neuron lacked the OneDR task. Save block added to support downstream population-level analysis.

---

## [2026-03-23] — Save PSTH by opto group (Fast / Others)

**Files changed:**
- `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Organize.m`

**What:** Moved save to after `Sel_Fast`/`Sel_Others` are defined. Now saves two separate files — `PopData_ObjLoc_<BatchFileName>_Fast.mat` and `PopData_ObjLoc_<BatchFileName>_Others.mat` — each containing `psth_object_high`, `psth_object_low`, `psth_location_high`, `psth_location_low`, `auc_psth_object`, `auc_psth_location`, `time_object`, `time_location` subsetted to the respective group.

**Why:** User requested separate files for Fast Excitation and No-Effect opto groups to enable group-specific downstream population analysis.

---

## [2026-03-23] — Save FEF and SC into separate subfolders

**Files changed:**
- `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Organize.m`

**What:** Brain area (FEF or SC) is parsed from `BatchFileName` prefix. Files are now saved to `Results/ObjectLocationPSTH/FEF/` or `Results/ObjectLocationPSTH/SC/` (created automatically if absent). File names include area label: `PopData_ObjLoc_<Area>_<Group>.mat`.

**Why:** FEF and SC are run as separate batches (`FEF_BothRA_StableValue.txt` / `SC_BothRA_StableValue.txt`), so outputs need to be kept in separate folders for downstream area-level comparisons.

---

## [2026-03-23] — Switch save format to -v7.3 for Python compatibility

**Files changed:**
- `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Organize.m`

**What:** Added `-v7.3` flag to `save()` call so output `.mat` files use HDF5 format, loadable in Python via `h5py`.

**Why:** User plans to switch to Python for population analysis. `-v7.3` removes the 2 GB size limit and is compatible with both `h5py` (Python) and MATLAB.

---

## [2026-03-23] — Headless batch runner + SC FreeViewNeural try-catch fix

**Files changed:**
- `AnalysisProtocol/BatchAnalysis/RunBatchObjLoc.m` (created)
- `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Organize.m`

**What:** Created `RunBatchObjLoc.m` to run FEF and SC batches headless via `matlab -batch`. Commented out `keyboard` to prevent hang. Fixed SC crash: wrapped FreeViewNeural data extraction in try-catch so neurons with incomplete DataStamp (missing keys like `PSTH_Stim_Mean_norm`) fall through to NaN defaults with `OptoIndex=0` instead of crashing. FEF completed successfully (77 Fast, 407 Others neurons).

**Why:** Some older SC recordings have a `FreeViewNeural` field but with an older/incomplete DataStamp that lacks newer keys, causing a `containers.Map` key-not-found error at runtime. After user re-ran the FreeViewNeural analysis for Robin060722FRACSAC8, SC rerun completed cleanly: 41 Fast, 233 Others neurons saved with no warnings.

---

## [2026-03-23] — Add per-trial decoder variables to single-neuron DataStamp

**Files changed:**
- `AnalysisProtocol/OneDR/OneDRAnalysis.m`
- `AnalysisProtocol/ObjDirectedSaccadeTask/ObjectResponsesDifferentLocation.m`

**What:** Added 8 new keys to the DataStamp export in both single-neuron analysis files: `psth_targeton_raw` (all trials × T, raw firing rate, no z-score), `target_index` (1=high value, 0=low value, all trials), `RF_index` (1=into RF, 0=out of RF, all trials), and trial counts `n_high_RF`, `n_low_RF`, `n_high_nonRF`, `n_low_nonRF`. Time selection uses the existing `sel`/`SelectTime` mask ([-400, 400] ms). RF selection uses `ismember(TargetAngle, UniqueAngle(max_idx))` (OneDR) and `ismember(TargetAngle, UniqueAngle(RF==1))` (ObjDir).

**Why:** Provides all necessary inputs for a single-neuron decoder: a raw PSTH matrix with trial labels for value and spatial condition, plus trial counts per condition for balanced sampling or cross-validation.

---

## [2026-03-23] — Remove redundant nonRF peak stat assignments in ObjectResponsesDifferentLocation

**Files changed:**
- `AnalysisProtocol/ObjDirectedSaccadeTask/ObjectResponsesDifferentLocation.m`

**What:** Removed two sets of redundant nonRF peak stat assignments. In the single-location block (`length(UniqueAngle)==1`), the correct NaN assignments at lines 1055–1058 were immediately overwritten by a buggy block that tried to index into empty arrays (`find(RF==0,1,'first')` returns `[]` when only one location exists) — removed the overwrite. In the multi-location `else` block, a verbatim duplicate assignment of `peakTime_AUC_nonRF` etc. appeared after the correct one — removed the duplicate. The loop-level peak finding already uses `max(abs(auc_sm_sel - 0.5))`, so nonRF values in the computed arrays are correct.

**Why:** The single-location overwrite was a latent bug (indexing with empty index). The multi-location duplicate was dead code.

---

## [2026-03-23] — Fix latency search direction in both single-neuron analysis files

**Files changed:**
- `AnalysisProtocol/OneDR/OneDRAnalysis.m`
- `AnalysisProtocol/ObjDirectedSaccadeTask/ObjectResponsesDifferentLocation.m`

**What:** Replaced `dirMask = AUC > 0.5` with a direction-aware mask based on `peakAUC_avg(i)`: if peak > 0.5, `dirMask = AUC > 0.5` (good-preferring); if peak < 0.5, `dirMask = AUC < 0.5` (bad-preferring). The loop condition simplified from `(all(segDir) || all(~segDir))` to `all(segDir)`. Also updated the comment in ObjDir to drop the misleading "if auc > 0.5" description.

**Why:** The old `all(~segDir)` branch (all bins ≤ 0.5) did not constrain the search to the same direction as the neuron's actual peak — it could trigger on any cluster of non-significant baseline bins that happened to be ≤ 0.5. Tying the direction mask to `peakAUC_avg` ensures latency is only detected in the direction the neuron actually responds.

## [2026-03-23 00:00] — TDR_ValueContext.m: Targeted Dimensionality Reduction script
**Files changed:** `TDR_ValueContext.m` (new)
**What:** Implements the full Targeted Dimensionality Reduction (TDR) pipeline from Tafazoli et al. (2026) adapted to the value × context experimental design. Includes pseudo-PSTH generation for testing, GLM regression per timepoint, PCA denoising, QR-based orthogonal axis extraction, projection of all 8 conditions onto TDR axes, value encoding distance computation, cross-context value axis angle, and a 12-panel figure.
**Why:** Enables population-level analysis of whether value is encoded in a shared subspace across object and location contexts, directly addressing the core scientific question of shared vs. distinct coding populations.

## [2026-03-23 00:01] — CrossDecoder_ValueContext.py + RunCrossDecoder.m
**Files changed:** `CrossDecoder_ValueContext.py` (new), `RunCrossDecoder.m` (new)
**What:** Python cross-context and cross-temporal value decoder with MATLAB wrapper. Implements same-time cross-context decoding (train obj→test loc and vice versa), within-context cross-temporal matrices, and cross-context cross-temporal matrix. Key optimisation: vectorised prediction at all test times via single matrix multiply (O(T) fits not O(T²)). Parallelised over t_train with joblib.
**Why:** Addresses whether value is decoded by the same rule across object and location contexts (cross-context generalisation), and whether value representations are temporally stable (cross-temporal generalisation).

## [2026-03-24 00:00] — CrossDecoder_ValueContext.m + PreprocessPSTH.m + PlotDecoderResults.m (MATLAB)
**Files changed:** `PreprocessPSTH.m` (new), `CrossDecoder_ValueContext.m` (new), `PlotDecoderResults.m` (new)
**What:** Full MATLAB implementation of the cross-context value decoder per DecoderPlan.md. PreprocessPSTH applies reflect-pad, 100ms Gaussian smooth, and 2x downsample. CrossDecoder_ValueContext runs 4 decodings (within same-time, cross same-time, within cross-temporal, cross cross-temporal) with 10-fold stratified CV, z-scoring inside folds, parfor parallelism, per-t_train checkpointing, full weight saving, progress logging, and quick-test mode. PlotDecoderResults produces 6-panel figure.
**Why:** Implements planned decoder pipeline to test whether value is encoded in a shared subspace across object and location contexts.

## [2026-03-24 12:35] — Load & save new PSTH fields in Batch_ObjLocValue_Organize
**Files changed:** `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Organize.m`, `plans/v1.md`
**What:** Added loading of 7 new per-neuron fields (`psth_targeton_raw`, `target_index`, `RF_index`, `n_high/low_RF/nonRF`) from both StableValue and OneDR DataStamps; built nonRF neurons×time PSTH matrices; updated save block to include all new fields, nonRF matrices, and ModuType for both Sel_Fast and Sel_Others groups.
**Why:** New fields were saved by upstream analysis but not yet collected or exported to the ObjectLocationPSTH folder.

## [2026-03-24 00:01] — TDR_ValueContext_RF.m
**Files changed:** `PopulationAnalysis/TDR_ValueContext_RF.m` (new)
**What:** RF-only variant of TDR_ValueContext.m. Uses 4 conditions (high/low value × object/location context, RF trials only). Drops the RF position predictor, retaining 3 axes: value, context, value×context. Otherwise identical pipeline: z-score, PCA denoise, GLM, QR orthogonalisation, projection, value encoding distance, cross-context angle.
**Why:** Matches the RF-only decision made for the cross-decoder, enabling consistent comparison between TDR geometry and decoding results.

## [2026-03-25] — Split per-trial PSTH into RF and nonRF trials

**Files changed:**
- `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Organize.m`

**What:** Added 4 loops (one per group × task combination) that use `RF_index_*_{Fast,Others}` to split each neuron's `psth_targeton_raw` matrix by row into RF trials and nonRF trials. Produces 8 new cell arrays named `psth_targeton_raw_{stable,oneDR}_{Fast,Others}_{RF,nonRF}`.

**Why:** Downstream decoders need RF and nonRF conditions separated to compare value coding in and out of the receptive field.

---

## [2026-03-25] — Trial-count filter + Fast/Others split for per-trial PSTH arrays

**Files changed:**
- `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Organize.m`

**What:** Added `TrialCount_OK` mask (n_high_RF ≥ 10 AND n_low_RF ≥ 10 for both stable and oneDR tasks). Combined with existing `Sel_Fast`/`Sel_Others` to form `Sel_Raw_Fast` and `Sel_Raw_Others`. Applied both masks to `psth_targeton_raw`, `target_index`, `RF_index` (cell arrays) and `n_high_RF`, `n_low_RF` (numeric arrays) for both tasks, producing 8 psth cell arrays, 8 index cell arrays, and 8 trial-count vectors named with `_Fast`/`_Others` suffix.

**Why:** Downstream population decoders require sufficient trials per condition (≥10) and need labels split by opto modulation group.

---

## [2026-03-25] — Load per-trial decoder variables in Batch_ObjLocValue_Organize

**Files changed:**
- `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Organize.m`

**What:** Added loading of `psth_targeton_raw`, `target_index`, `RF_index`, `n_high_RF`, `n_low_RF`, `n_high_nonRF`, `n_low_nonRF` for both StableValue (suffix `_stable`) and OneDR (suffix `_oneDR`) tasks. NaN fills added to the OneDR else-branch for neurons that lack that task.

**Why:** These keys were already saved in both single-neuron DataStamps but were not collected in the batch organizer, making them unavailable for downstream population decoding.

---

## [2026-03-25] — Save 4 RF PSTH matrices to FEF/SC subfolders
**Files changed:** `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Organize.m`
**What:** Added save block after the selection assignments. Parses `BatchFileName` with `contains()` to set `AreaLabel` = 'FEF' or 'SC'. Creates `Results/ObjectLocationPSTH/<Area>/` if absent. Saves two `.mat` files (Fast and Others), each containing the 4 RF mean PSTHs: `PSTH_TgtOn_Good/Bad_RF_Stable_<Group>`, `PSTH_TgtOn_Good/Bad_RF_OneDR_<Group>`, plus `PSTH_Time_TgtOn_stable` and `PSTH_Time_TgtOn_OneDR`. Format: `-v7.3`.
**Why:** Population RF PSTHs needed for downstream analysis; separate files per opto group and brain area keep results organized.

---

## [2026-03-25] — Panel 3: 4-circle cyclic Venn diagram
**Files changed:** `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Organize.m`
**What:** Replaced 2-circle Venn + outer ring with a 4-circle cyclic Venn. Circles: A=top (stable value, green), B=right (oneDR value, purple), C=bottom (stable non-value, pale green), D=left (oneDR non-value, pale purple). Adjacent pairs overlap (A∩B=both, B∩C=oneDR only, C∩D=neither, D∩A=stable only); opposite pairs (A-C, B-D) are geometrically disjoint by design (h=0.36, r=0.30). Every neuron falls in exactly one of the 4 intersection regions.
**Why:** Original 2-circle Venn couldn't show non-value neurons as a circle. The cyclic arrangement is the only topology where all 4 groups form valid pairwise intersections with no spurious 3-way or 4-way overlaps.

---

## [2026-03-25] — Panels 1–2 → pie charts; Panel 3 Venn add outer "neither" circle
**Files changed:** `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Organize.m`
**What:** Replaced panels 1 & 2 circle diagrams with `pie()` charts (3 slices: High pref blue, Low pref orange, No pref gray; counts in labels). Fixed panel 3 Venn by adding a large outer gray circle (all neurons, drawn first opaque) with two semi-transparent overlapping inner circles (stable green, oneDR purple) on top; the visible gray ring represents "neither" neurons, with count labeled below.
**Why:** Pie charts cleanly show proportions for mutually exclusive groups. The Venn outer circle was missing, so the "neither" region had no visual representation.

---

## [2026-03-25] — Add normalized AUC scatter plot
**Files changed:** `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Organize.m`
**What:** Added `x_norm = abs(x_sc - 0.5)` and `y_norm = abs(y_sc - 0.5)` and a second scatter+histogram figure identical in layout to the raw-AUC scatter. Axes range [0, 0.52]; same 4 significance color groups and filled/open marker convention; marginal histograms show norm-AUC distributions with sig neurons filled.
**Why:** Normalized AUC collapses high-prefer and low-prefer into a single unsigned discriminability axis, making overall value coding strength directly comparable between tasks.

---

## [2026-03-25] — Fix scatter/Venn selection + add outer circle to Venn panels 1–2
**Files changed:** `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Organize.m`
**What:** Bug 1: replaced `Sel_All` (required OptoIndex==1, effectively showing only opto-recorded neurons) with `Sel_AUC = (OneDRIndex'==1) & ResponseValid & VisualResponsive'` so scatter/Venn include all neurons with both tasks. Bug 2: Venn panels 1 and 2 now draw a large gray outer circle (all neurons) with two opaque colored inner circles (high/low prefer) on top; the gray ring visible around the inner circles represents no-preference neurons.
**Why:** Plots were only showing Sel_Others-sized data due to opto filter; Venn had no visual region for non-value-coding neurons.

---

## [2026-03-25] — Venn diagrams: value preference and cross-task overlap
**Files changed:** `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Organize.m`
**What:** Added a 3-panel Venn figure. Panels 1 & 2 show high-prefer (AUC>0.5, sig) and low-prefer (AUC<0.5, sig) as two non-overlapping circles, with no-preference count outside, for stable and oneDR tasks respectively. Panel 3 is a classic overlapping Venn showing cross-task value coding: stable-only, both, oneDR-only, and neither. Colors: blue=high-prefer, orange=low-prefer, green=stable, purple=oneDR.
**Why:** Quantifies how many neurons prefer high vs low value in each task, and how much the value-coding populations overlap across tasks.

---

## [2026-03-25] — Scatter plot: peakAUC stable vs oneDR with marginal histograms
**Files changed:** `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Organize.m`
**What:** Added scatter plot (X = peakAUC_Stable_RF, Y = peakAUC_OneDR_RF, all Sel_All neurons) with 4 color groups (oneDR-only orange filled, stable-only blue filled, neither gray open, both green filled). Top and right marginal histograms overlaid with open bars (all neurons) and filled bars (significant neurons); top uses blue for sig_stable, right uses orange for sig_oneDR. Axes are linked for synchronized zoom. Reference dashed lines at AUC=0.5.
**Why:** Visual inspection of how many neurons encode value in each task context and whether they overlap.

---

## [2026-03-24 00:00] — TDR_ValueContext_Separate: independent TDR per context
**Files changed:** `PopulationAnalysisSeparate/TDR_ValueContext_Separate.m` (new)
**What:** Implements Targeted Dimensionality Reduction run separately for the object value context (psth_high/low_value_object_RF) and the location value context (psth_high/low_value_location_RF). Each context uses a single predictor (value: +1/−1). A shared PCA denoising matrix (fit on all 4 conditions) keeps both value axes in a common reference frame. The function outputs own-axis and cross-axis projections, value encoding distances, and the angle between the two value axes.
**Why:** Allows direct comparison of the context-specific value subspace geometry without the shared-predictor assumptions of TDR_ValueContext_RF.m.

---

## [2026-03-25] — Fix save block: replace old PSTH matrices with high/low per-trial cell arrays
**Files changed:** `AnalysisProtocol/BatchAnalysis/Batch_ObjLocValue_Organize.m`
**What:** Replaced the 4 old population-average matrices (`PSTH_TgtOn_Good/Bad_RF_Stable/OneDR_Fast/Others`) in both save calls with the 8 per-trial high/low cell arrays produced by the `%% High/low value split` section: `psth_targeton_raw_{stable,oneDR}_{Fast,Others}_{RF,nonRF}_{high,low}`, plus time vectors.
**Why:** The save block was saving the wrong variables — the old population-average matrices instead of the new per-trial raw PSTH cell arrays needed for downstream decoding.
