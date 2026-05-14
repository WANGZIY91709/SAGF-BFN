# SAGF-BFN

## Overview

This repository provides the implementation of **SAGF-BFN** (*Structure-Aware Graph Fusion for Brain Functional Networks*), a multi-view brain functional network (BFN) fusion framework for EEG-based depression analysis.

The proposed framework jointly models:

- **View-specific topological structures** through structure-aware graph embedding learning;
- **Cross-view shared neural interaction patterns** through tensor-based shared structure compression;
- **Unified fused brain network representation** via adaptive multi-view fusion.

The repository mainly includes the optimization procedure of SAGF-BFN. The statistical analysis and machine learning classifiers used in the paper can be directly implemented using standard MATLAB built-in functions.

---

# Repository Structure

```text
SAGF-BFN/
│
├── solver.m                 % Main optimization procedure of SAGF-BFN
├── Gshrink.m                % Auxiliary function for shrinkage operation
├── L2_distance_1.m          % Pairwise Euclidean distance computation
├── normalize_columns.m      % Column normalization function
├── README.md
