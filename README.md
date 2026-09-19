# Oxygen Subduction in the Pacific Ocean

MATLAB code used to calculate and decompose oxygen-subduction fluxes and to reproduce the figures associated with the accompanying study.

## Repository structure

```text
oxygen-subduction-ecco/
├── cal_script/     Calculation and preprocessing scripts
├── fig_script/     Main-text and supplementary figure scripts
├── function/       Project functions and bundled third-party MATLAB tools
├── input_data/     User-supplied input data (not tracked by Git)
├── result_data/    Calculated results (not tracked by Git)
└── figs/           Generated figures (not tracked by Git)
```

## Data

The scientific input and result files are not included because several individual files exceed GitHub's regular 100 MiB file limit and the full working directory is several gigabytes. The excluded prepared input and result files will be available from Zenodo (link: https://zenodo.org/uploads/22155114). Small `.mat` resources required by functions and bundled toolboxes remain in `function/`.

See [DATA_MANIFEST.md] for the required file names, expected locations, producing scripts, and known gaps. Before running a calculation or figure script, place the corresponding data files in `input_data/` or `result_data/`.

`PhscData_monthly.mat` integrates physical variables from monthly ECCO NetCDF files and includes vertical averages within the mixed layer. `oxydata_monthly.mat` integrates GOBAI-O2 and includes oxygen averages within the mixed layer and within the 10-m water column immediately below it. `Ot_section_monclim.mat` and `UncO2_Hmax.mat` are secondary products derived from GOBAI-O2.

## Usage

1. Clone or download this repository.
2. Supply the files listed in `DATA_MANIFEST.md`.
3. Open MATLAB.
4. Run a script directly from `cal_script/` or `fig_script/`. The script adds `function/` and its subdirectories to the MATLAB path automatically.

For example:

```matlab
run(fullfile('cal_script','cal_oxygen_flux.m'))
run(fullfile('fig_script','FIG1_SoxSpatial.m'))
```

Calculation outputs are written to `result_data/`, and figures are written to `figs/`.

## Main calculation workflow

```text
Physical and oxygen input fields
              |
              v
   cal_oxygen_flux.m
              |
              v
OSresults_monthly.mat
              |
              +--> cal_seas_corr_maxR.m
              +--> cal_eddy_regions.m
              +--> figure scripts
```

`run_calculate_mxl_monthly.m` calculates monthly mixed-layer depth (`Hml`), maximum mixed-layer depth (`Hmax`), and seasonal thermocline thickness (`Hst`). It directly reads monthly ECCO files named `THETA_YYYY_MM.nc` and `SALT_YYYY_MM.nc` from `input_data/raw/`.

## MATLAB dependencies

- MATLAB with support for `matfile`, `ncread`, `tiledlayout`, and `exportgraphics`
- Signal Processing Toolbox for filtering used by `FIGS8_ENSO.m`
- Statistics and Machine Learning Toolbox for statistical routines used by selected scripts
- GSW Oceanographic Toolbox, bundled in `function/gsw_matlab_v3_06_16_1/`
- M_Map, bundled in `function/m_map/`

Additional bundled functions retain their original notices where supplied. See [THIRD_PARTY_NOTICES.md].

## Citation and license

When the associated article and repository DOI are available, their complete citation will add here and cite them when reusing this code.

This repository currently does not declare a project-level software license. Before public release, the authors will select a license consistent with their intended reuse terms and verify compatibility with the bundled third-party components.
