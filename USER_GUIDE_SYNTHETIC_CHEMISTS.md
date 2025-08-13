# 🧪 User Guide for Synthetic Chemists

This document provides a practical, step-by-step guide for synthetic chemists to evaluate reaction barriers and transition states using ConfoTS.

## 1. Purpose of ConfoTS in Synthetic Planning

ConfoTS can be used to:

- Predict **relative energy barriers** for competing reaction pathways.
- Compare **diastereoselectivity** or **regioselectivity** under the Curtin–Hammett scenario.
- Validate proposed mechanistic hypotheses before committing to multi-step synthesis.

## 2. Minimum Input Required

You only need:

- A **reaction scheme** in `.rxn` format (can be drawn in ChemDraw or similar).
- **Basic knowledge of your reaction** (reactants, products, stereochemistry).
- No need to manually prepare 3D geometries or atom mappings—ConfoTS does this automatically.

## 3. Quick Start Workflow

### 1. Prepare Your Reaction File

- Draw your reaction in ChemDraw.
- Save/export as **MDL RXN file** (`.rxn`) — in ChemDraw, this appears as **“ISIS Reaction V2000 (\*.rxn)”** in the save dialog.
- Include correct stereochemistry and charges.
- If you wish to evaluate **multiple selectivity outcomes** (e.g., regioisomers or diastereomers), draw the **products side-by-side** in the same reaction scheme.
  - _Note:_ In ChemDraw, products must be placed **horizontally** to ensure correct recognition in the `.rxn` file.

### 2. Set Up the Working Directory

- Create a new folder for your reaction.
- Copy the provided `confots_step*.py` scripts and `template/` folder into it.
- Place your `.rxn` file inside.

### 3. Choose Your Mode of Execution

You can run the workflow in **two ways**:

#### Option 1: Single-step run (recommended for most users) — executes all steps automatically.

1. Open a terminal or Command Prompt ().
2. Move into your working directory:
   ```bash
   cd path/to/your/working_directory
   ```
3. Run the workflow:
   ```bash
   python3 run_confots.py
   ```
   - On **HPC clusters** with job scheduling (e.g., PBS):
     ```bash
     qsub run_confots.sh
     ```

#### Option 2: Step-by-step run (advanced) — runs each stage individually for more control.

1. Move into your working directory (as above).
2. Run the scripts in order:
   ```bash
   python3 confots_step0.py
   python3 confots_step1.py
   ...
   python3 confots_step5.py
   ```

_Tip:_ On most systems, you can right-click inside your working directory and choose **“Open Terminal”** or **“Open Command Prompt”** to begin.

### 4. Check the Output

- After completion, view the **energy diagrams** in:
  - `130_energy_diagram_xtb/` (fast, semiempirical estimate)
  - `310_energy_diagram_dft/` (more accurate DFT results)
- Transition state structures are stored in:
  - `126_emp_ts_xyz/` — all validated TS structures from semiempirical calculations (for all conformers)
  - `312_dft_min_xyz/` — DFT-optimized structures of the lowest-energy TS for each pathway
- In both the diagrams and output file names, compounds are automatically assigned **lowercase letter indices** (`a`, `b`, `c`, …) in the order they appear from **left to right** in the original reaction scheme.

## 4. Examples of Potential Use Cases for Synthetic Chemists

The following are representative reaction types where ConfoTS **may** be applicable, provided the reaction mechanism is consistent with a single-transition-state model:

- **Diastereoselective cycloadditions** — Compare TS energies for exo vs. endo approaches.
- **1,3-dipolar cycloadditions** — Assess regio- and diastereoselectivity via competing concerted TSs.
- **[3,3]-sigmatropic rearrangements (Claisen/Cope)** — Evaluate substituent and conformational effects on a single concerted TS.
- **SN2 substitutions** — Compare competing nucleophiles, leaving groups, or stereoelectronic approaches.

## 5. Practical Tips

- Start with default settings in `config.ini`—only adjust if needed.
- If your molecule is very large (> 500 Da), consider increasing `max_mw` in `config.ini`.
- For rapid screening, rely on the `xtb`-based diagram (`130_energy_diagram_xtb`) before committing to DFT.
- Always confirm key results with **DFT-level** calculations before publication.

## 6. Scope and Limitations

- Designed for **single-transition-state** processes (e.g., pericyclic reactions).
- Not intended for **complex multistep mechanisms** or **catalytic cycles**.
- Evaluates only the specific reaction outcomes defined in the input file.
- The program performs highly complex automated procedures and cannot account perfectly for every possible reaction type.
- Depending on the substrate structure and the type of reaction, conformational sampling for transition states may still miss certain possibilities due to inherent workflow limitations.
- **Always verify not only the energy diagrams but also the intermediate files, calculation steps, and log files to ensure that the results are chemically meaningful.**
