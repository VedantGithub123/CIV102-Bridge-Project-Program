# CIV102 Bridge Project

This repository contains the design, structural analysis, and SolidWorks files for a matboard bridge. The Python model analyzes a simply supported bridge under a moving three-car train load and compares the resulting shear-force and bending-moment envelopes with several material, glue, and buckling failure limits.

All calculations use Newtons and millimetres. The current design models a 1200 mm bridge with 1.27 mm matboard.

## Getting Started

The analysis requires Python 3 and the following packages:

- NumPy
- Matplotlib

Clone the repository and enter its directory:

```text
git clone <repository-url>
cd CIV102-Bridge-Project-Program
```

Install the dependencies with pip:

```text
python -m pip install numpy matplotlib
```

## Repository Structure

```text
CIV102-Bridge-Project-Program/
├── archive/       Previous design iterations and the physical-design model
├── build/         SolidWorks assembly, parts, drawings, and cross-sections
├── config/        Current material, geometry, load, and failure parameters
├── media/         Photos documenting the construction process
├── outputs/       Generated or exported shear-force and bending-moment plots
├── reports/       Design report, calculations, and assembly documentation
├── src/           Python analysis code
└── README.md
```

## Running the Analysis

Run the simulation from the repository root:

```text
python src/simulation.py
```

The script opens Matplotlib figures for the bridge, selected cross-sections, shear-force diagrams, and bending-moment diagrams. It also prints the maximum allowable load and the factor of safety for each modeled failure mode. The calculation is performed when the file is run; there is currently no separate command-line interface.

The script imports the design parameters from `config/design.py`. After changing the geometry or material values, rerun the simulation to recalculate the results.

## Design Configuration

`config/design.py` contains:

- Matboard tensile, compressive, shear strength, Young's modulus, and Poisson's ratio
- Glue shear strength
- Bridge length and matboard thickness
- Cross-section key points along the span
- Glue-tab locations
- Diaphragm locations
- Locations and dimensions used for the four buckling failure models
- The optional `case1` load-case switch

Cross-sections and related geometry are linearly interpolated between the configured key positions. The current load model distributes the train weight across six wheel loads and evaluates the train at successive positions along the bridge.

## Analysis

The main module, `src/simulation.py`, calculates:

- Cross-section centroid and second moment of area
- Support reactions, shear-force diagrams, and bending-moment diagrams
- Maximum shear and bending-moment envelopes for all train positions
- Tensile and compressive matboard failure limits
- Matboard and glue shear failure limits
- Four buckling failure limits
- Overall maximum load and failure-mode-specific factors of safety

The plotting functions compare the load envelopes against the corresponding allowable values. The simulation neglects bridge self-weight and treats diaphragms as infinitely thin. It also uses simplified assumptions for glue-tab behavior and side-plate shear buckling; see the comments in `src/simulation.py` for the full list.

## Outputs and Project Files

The `outputs/` directory currently contains exported examples of:

- Shear-force diagrams and envelopes
- Bending-moment diagrams and envelopes

The simulation itself displays figures but does not save them automatically. Save figures manually from the Matplotlib windows if new plots are required.

The `build/` directory contains the CAD assembly, individual matboard parts, cross-section parts, and engineering drawings. The `reports/` directory contains the design report, calculations, and assembly documentation. `archive/` preserves earlier design iterations for comparison.

## Design Process

The bridge design uses a variable cross-section, diaphragms, glue tabs, and local reinforcement to balance strength, stiffness, and material use. The geometry was developed in SolidWorks and represented in Python as interpolated rectangular regions so that section properties and failure limits can be evaluated along the entire span.

The analysis workflow is:

1. Define material properties and bridge geometry in `config/design.py`.
2. Calculate the structural response for all relevant train positions.
3. Calculate allowable loads for material, glue, and buckling failure modes.
4. Determine the governing maximum load and factors of safety.
5. Compare the analysis with the CAD drawings and the submitted reports.