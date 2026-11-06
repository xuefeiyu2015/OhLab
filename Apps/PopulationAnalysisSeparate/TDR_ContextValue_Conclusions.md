# Context-Dependent Value Modulation along the FEF–SC Pathway

## Overview

We applied Targeted Dimensionality Reduction (TDR; Mante et al., 2013) to population activity in FEF and SC recorded during two value-based tasks: a stable object-value task and a one-directional reward (oneDR) spatial-value task. The TDR framework extracts three orthogonal axes from the population response — a context axis (stable vs. oneDR), a value main-effect axis (high vs. low reward), and a context × value interaction axis — and projects all four task conditions (stable-high, stable-low, oneDR-high, oneDR-low) onto this low-dimensional space.

---

## Main Findings

### 1. Object value and spatial value are encoded in orthogonal neural subspaces

Within-context value axes were computed separately for each task (d_val_stable from stable-high vs. stable-low; d_val_oneDR from oneDR-high vs. oneDR-low). The angular separation between these axes was near-orthogonal in both areas:

- FEF: |cos(θ)| = 0.093 (~85°)
- SC:  |cos(θ)| = 0.214 (~78°)

A cos(θ) near zero means the neural population uses fundamentally different linear combinations of neurons to encode value in the two tasks. Neither area maintains a single value dimension that applies across task contexts.

### 2. Context-dependent value coding is as strong as context-independent coding

The near-zero cos(θ) between within-context value axes has a direct mathematical consequence. In the balanced 2×2 factorial design:

```
d_val_stable ∝ β_value + β_interaction
d_val_oneDR  ∝ β_value − β_interaction
```

For their dot product to equal zero: ||β_interaction|| ≈ ||β_value||. The interaction (context × value) term is comparable in magnitude to the value main effect. There is no dominant context-independent value signal — value coding in FEF and SC is inherently contextualized.

### 3. FEF and SC show similar value geometry after removing the common response

The apparent difference in value separation between FEF and SC in the raw TDR projections is a scale artifact. SC has a large non-selective response to target onset that projects onto the value axis, inflating the denominator without contributing to discrimination. After de-meaning across conditions, value discrimination strength is comparable between areas. The difference between areas lies in the magnitude of the shared (non-discriminative) response, not in value coding per se.

### 4. FEF and SC do not support cross-context value generalization

A linear decoder trained on value in one task context would perform at chance when tested on the other, because the value axes are orthogonal. This distinguishes FEF and SC from areas thought to carry abstract economic value representations (e.g., OFC, vmPFC), where value axes are expected to generalize across task contexts.

---

## Complementary Evidence: Population Decoding

Population decoding provides a functional counterpart to the geometric findings from TDR, testing the same question — does value information generalize across task contexts — from a different angle.

### Within-context decoding succeeds in both areas

Linear SVM decoders trained and tested within the same task context (stable-only or oneDR-only, cross-validated) achieved well above-chance accuracy in both FEF and SC. This confirms that value information is reliably encoded in the population response of each area, and that the near-orthogonal value subspaces identified by TDR are not a result of weak or noisy value coding — the value signal within each context is strong and decodable.

### Cross-context decoding fails

When decoders were trained on one task context and tested on the other (train stable → test oneDR, or vice versa), accuracy dropped to near chance in both FEF and SC. A decoder that learned to distinguish high from low value in the object-value task could not transfer that knowledge to the spatial-value task, and vice versa.

### Convergence between TDR and decoding

The two analyses arrive at the same conclusion through independent approaches:

| Analysis | Method | Finding |
|---|---|---|
| TDR | Geometry of population subspaces | Object-value and spatial-value axes are nearly orthogonal (cos ≈ 0.09–0.21) |
| Cross-context decoding | Functional generalization of a linear classifier | Decoding accuracy drops to chance when context changes |

This convergence strengthens the conclusion: the failure of cross-context decoding is not a limitation of the decoder or a noise issue, but a direct consequence of the population geometry revealed by TDR. The value axis in one context is orthogonal to the value axis in the other, so any linear readout optimized for one is uninformative for the other.

---

## Opto-Group Analysis: FEF-Recipient Neurons in SC Show Partial Value Generalization

ChR2-expressing neurons were classified as Fast (direct opto excitation) and compared to Others. In FEF, Fast neurons (SC-projecting) are indistinguishable from the broader population in their value geometry. In SC, Fast neurons (FEF-recipient) show notably higher cross-context value consistency:

| Group | |cos(θ)| | Angle |
|---|---|---|
| FEF Fast (SC-projecting) | 0.10 | ~84° |
| FEF Others | 0.09 | ~85° |
| **SC Fast (FEF-recipient)** | **0.40** | **~66°** |
| SC Others | 0.16 | ~81° |

### Interpretation

FEF projection neurons encode value in a highly context-specific way, identical to the rest of FEF. Yet their postsynaptic targets in SC show substantially higher cross-context value consistency. **Context generalization of value emerges at or after the FEF→SC synapse, not in FEF itself.**

One likely mechanism is convergence: each FEF-recipient SC neuron integrates input from many FEF projection neurons. Averaging across multiple context-specific vectors with partially different orientations amplifies the shared (context-general) component while context-specific components partially cancel. A filtering or gain-normalization effect at the FEF→SC synapse may contribute as well.

The non-FEF-recipient SC population (Others, cos = 0.16) remains largely context-specific, indicating that the more context-generalized value coding in SC is concentrated in the FEF-driven subpopulation rather than being a global property of SC.

---

## Summary and Interpretation

FEF and SC encode value in a task-contextualized manner: object value (stable task) and spatial value (oneDR task) occupy largely orthogonal subspaces in the population state space of both areas. Neither area maintains an abstract value signal that generalizes across task contexts.

This is consistent with the roles of FEF and SC as late-stage sensorimotor areas that compute value in the context of a specific action plan. Value is represented as "what this object is worth in this task" or "where this spatial target is worth looking in this task," rather than as a context-free economic quantity.
obj
The FEF→SC pathway introduces a partial transformation: FEF sends a context-specific value signal, and its recipient neurons in SC express modestly more cross-context value consistency. This may reflect the requirement that SC, as the final arbiter of saccade target selection, integrates value information across different task demands to generate a single spatial decision.

### The complete story in four points

- **1. Heterogeneous object and location value representation among FEF and SC populations.**
  *FEF and SC contain distinct subpopulations selectively tuned to object value, spatial value, or both, revealing task-specific value channels at the single-neuron level.* Within both FEF and SC, individual neurons are not uniformly tuned: some respond selectively to object value (stable task only), others to spatial/location value (oneDR task only), some to both, and some to neither. This heterogeneity is the cellular substrate for everything that follows — task-specific value coding at the population level is not an emergent surprise but a direct consequence of neurons being tuned to different value dimensions.

- **2. Object and spatial value occupy orthogonal subspaces in the FEF and SC population state space.**
  *Segregated object-value and spatial-value neuron populations produce orthogonal coding subspaces in the FEF and SC population state space.* Object-value-only neurons drive the value axis in the stable context; location-value-only neurons drive a different axis in the oneDR context. Because these populations are largely non-overlapping, the two within-context value axes end up nearly orthogonal (cos ≈ 0.09–0.21 across areas and groups). The population geometry is a direct reflection of the single-neuron tuning landscape: the more segregated the two value populations, the more orthogonal the subspaces.

- **3. Partial value integration across contexts emerges selectively in FEF-recipient neurons in SC.**
  *FEF-recipient neurons in SC exhibit greater cross-context value consistency than FEF projection neurons or other SC neurons, indicating that partial value integration emerges at the FEF→SC synapse.* Despite both areas having similarly orthogonal value subspaces overall, FEF-recipient neurons in SC (Fast SC, cos ≈ 0.40) stand out from the rest of SC (cos ≈ 0.16) and from FEF projection neurons themselves (cos ≈ 0.10). Neurons that respond to both object and location value — the "both" population — are the most likely candidates to drive this higher cross-context consistency in FEF-recipient SC neurons, as they carry value information that partially bridges the two task contexts. The FEF→SC synapse may act as a convergence point where the object-value and location-value streams begin to be combined.

- **4. Context generalization of value emerges along the FEF–SC pathway, not in FEF itself.**
  *The FEF–SC pathway implements a hierarchical transformation from task-specific value encoding in FEF to a partially generalized value signal in SC that supports context-flexible saccade target selection.* FEF encodes value in a task-specific, segregated manner. As this signal propagates to SC via direct FEF axons, the FEF-recipient SC subpopulation expresses a modestly more generalized value representation. SC as a whole, however, remains largely context-specific — the generalization is not global but confined to the FEF-driven subpopulation. This architecture allows SC to retain task-specific value signals for precise saccade programming while the FEF-recipient neurons provide a partially integrated value estimate that could support more flexible, context-tolerant target selection.

---

## Caveats

- The context axis separation in TDR reflects all differences between the two task states (sensory, motor, cognitive), not pure context representation. Claiming "context coding" requires ruling out sensory/motor confounds.
- Near-orthogonal value axes are partly expected given that object identity and spatial location are fundamentally different value-binding dimensions. The geometry does not prove active maintenance of separate representations.
- The cos(θ) advantage of SC Fast neurons should be verified against potential SNR confounds: if FEF-recipient neurons have higher firing rates or stronger target responses, this alone could produce higher apparent cos(θ).

---

*Analysis: AUC_ContextDiscrimination.m*
*Date: 2026-04-02*
