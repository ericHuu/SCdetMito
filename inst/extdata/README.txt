Bundled example data
====================

demo_pbmc_multidonor_seurat.rds
-------------------------------

`demo_pbmc_multidonor_seurat.rds` is the bundled SCdetMito demonstration
object. It is derived from public 10x Genomics GEM-X Human PBMC Donor 1-4
datasets generated with GEM-X Single Cell 3' chemistry and Cell Ranger 9.0.0.
The source datasets are distributed by 10x Genomics under the CC BY 4.0
license.

This object is intended for SCdetMito examples, software tests, and
multi-sample or group-aware QC demonstrations. It is not intended for
biological inference.

The bundled object contains four donor-level samples:

- donor1: male, age group 18_35
- donor2: male, age group 18_35
- donor3: female, age group 36_50
- donor4: female, age group 36_50

Metadata columns include sample, donor, sex, age_group, group, condition,
tissue, species, dataset, data_source, and demo_note.

Raw counts and mitochondrial features are retained. `mitoRatio` was calculated
with `ensure_mito_ratio(species = "human")` and is stored as a fraction in
0-1 scale.

The maintainer build script downsampled each donor to 2,500 cells using fixed
seed 1401, yielding 10,000 cells total. The maintainer-side construction script
is kept in the GitHub repository under `data-raw/`, is excluded from the source
package build, and is not run during package installation, tests, vignettes, or
R CMD check.

No private, unpublished, or patient-identifiable data are included.
