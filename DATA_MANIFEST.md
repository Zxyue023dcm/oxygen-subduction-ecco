# Data manifest

Large scientific data are intentionally excluded from Git tracking. The prepared input and result files listed below will be archived on Zenodo (link: https://zenodo.org/uploads/22155114). Directory placeholders are included so that the expected layout is preserved after cloning.

Original ECCO NetCDF files should be downloaded from the official PO.DAAC collection and placed in `input_data/raw/`. They are not redistributed in this code repository.

## Input data

Place these files in `input_data/` unless another location is stated.

| File | Approximate working-copy size | Source and processing | Status in this code release |

| `PhscData_monthly.mat` | 1312.21 MiB | Integrates physical variables from monthly ECCO raw files, including vertical averages within the mixed layer and the seasonal thermocline. `run_calculate_mxl_monthly.m` appends `Hml`, `Hmax`, and `Hst`. | Excluded; Zenodo link |
| `oxydata_monthly.mat` | 294.85 MiB | Integrates GOBAI-O2, including oxygen averages within the maximum mixed layer and within the 10-m water column immediately below it. | Excluded; Zenodo link |
| `nino34.mat` | 0.01 MiB | Niño 3.4 index from NOAA PSL, restricted to 2004–2017 for this analysis. | Excluded; Zenodo link |
| `GOBAI-O2-v2.2.nc` | Not present in the prepared working directory | Source GOBAI-O2 field used by `cal_o2_uncer.m` and `cal_Ot_section.m`. | Excluded; GOBAI-O2 official source |

## Raw ECCO files

`run_calculate_mxl_monthly.m` directly reads monthly potential-temperature and salinity files from `input_data/raw/` using this naming convention:

```text
input_data/raw/
├── THETA_2004_01.nc
├── SALT_2004_01.nc
├── ...
├── THETA_2017_12.nc
└── SALT_2017_12.nc
```

The script expects one `THETA_YYYY_MM.nc` file and one `SALT_YYYY_MM.nc` file for every month from January 2004 through December 2017, downloaded from https://ecco.jpl.nasa.gov/drive/files/Version4/Release4/interp_monthly. Each file must use the ECCO 0.5-degree product variables `THETA` or `SALT` together with `longitude`, `latitude`, and `Z`. Only `input_data/raw/.gitkeep` is tracked; the NetCDF files themselves are ignored by Git.

## Result data

Generated or precomputed result files belong in `result_data/`.

| File | Approximate working-copy size | Producer or provenance | Status in this code release |

| `OSresults_monthly.mat` | 769.90 MiB | Calculated by `cal_oxygen_flux.m` from the prepared ECCO and GOBAI-O2 inputs. | Excluded; Zenodo link |
| `Wek.mat` | 388.79 MiB | Calculated by `cal_ekman.m`. | Excluded; Zenodo link |
| `seas_corr_maxR.mat` | 7.43 MiB | Calculated by `cal_seas_corr_maxR.m`. | Excluded; Zenodo link |
| `BASINed_maxR.mat` | 0.24 MiB | Calculated by `cal_eddy_regions.m`. | Excluded; Zenodo link |
| `UncO2_Hmax.mat` | 71.39 MiB | Secondary product derived from GOBAI-O2 by `cal_o2_uncer.m`. | Excluded; Zenodo link |
| `Ot_section_monclim.mat` | 1.18 MiB | Secondary climatological section product derived from GOBAI-O2. | Excluded; Zenodo link |
| `OSresults_monthly_ecco.mat` | 765.21 MiB | Precomputed ECCO comparison product. | Excluded; Zenodo link |
| `Straj_particle.mat` | 39.84 MiB | Precomputed particle-trajectory product. | Excluded; Zenodo link |

## Small support resources retained in the repository

The following `.mat` files are not research outputs. They are small resources required by plotting functions, masks, or bundled third-party tools and are therefore tracked:

- `function/Areaall.mat`
- `function/colorData.mat`
- `function/PO_subbasins.mat`
- `function/slanCM_Data.mat`
- `function/gsw_matlab_v3_06_16_1/library/gsw_data_v3_0.mat`
- `function/m_map/private/igrf.mat`
- `function/m_map/private/m_coasts.mat`
- `function/m_map/private/m_topo.mat`

## Data references

1. ECCO Consortium, Fukumori, I., Wang, O., Fenty, I., Forget, G., Heimbach, P., & Ponte, R. M. (2021). *ECCO Ocean Temperature and Salinity—Monthly Mean 0.5 Degree (Version 4 Release 4)* (V4r4) [Data set]. PO.DAAC. https://doi.org/10.5067/ECG5M-OTS44
2. Sharp, J. D., Fassbender, A. J., Carter, B. R., Johnson, G. C., Schultz, C., & Dunne, J. P. (2022). *GOBAI-O2: A Global Gridded Monthly Dataset of Ocean Interior Dissolved Oxygen Concentrations Based on Shipboard and Autonomous Observations* (NCEI Accession 0259304) [Data set]. NOAA National Centers for Environmental Information. https://doi.org/10.25921/Z72M-YZ67
3. NOAA Physical Sciences Laboratory. (2025). *Niño 3.4 SST index from NOAA ERSSTv5, 2004–2017* [Data set]. https://psl.noaa.gov/data/timeseries/month/Nino34_CPC/


