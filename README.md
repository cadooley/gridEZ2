
# gridEZ2

gridEZ generates a sampling frame consisting of gridded enumeration
zones (EZs). EZs are constructed based on a target population count per
EZ whilst being restricted by a maximum geographic size.

## Installation

You can install gridEZ2 like so:

``` r
# Install devtools if it's not already installed
install.packages("devtools")

# Install the package from GitHub
devtools::install_github("cadooley/gridEZ2")
```

## Executing the function

You will need to prepare three raster input files before running the
function. These are: - population raster - stratification raster 1 -
stratification raster 2 The three rasters need to be in the same
coordinate reference system (crs), e.g. WGS84, and be of the same
spatial resolution. When running the gridEZ function, you provide the
file path for each of the rasters: population_raster_path (for the
population raster), strata1_raster_path (for one of your stratification
raster) and strata2_raster_path (for your second stratification raster).
For household or population surveys, the stratification rasters are
commonly settlement type (urban/rural) and administrative units.
However, there are many possible ways a survey could be stratified,
e.g. distance to a point of interest or topological feature, climate
zones, etc. and you can use any desired stratification rasters you wish
for gridEZ. Each stratum is treated as having uncrossable boundaries and
therefore EZs will never include grid cells of different strata.

You can execute the gridEZ function using one of the following code
examples after editing your raster and output file paths as well as
“ncores”:

``` r
library(parallel)

# gridEZ utilises parallel processing to generate EZs for different strata. We recommend using your total number of computer cores minus 1. This will tell you the total number of cores available on your computer:
parallel::detectCores() 
# we advise that you specify the following value for the "ncores" parameter in the gridEZ function:
parallel::detectCores()  - 1

library(gridEZ)

?gridEZ

# run gridEZ using default settings. This will produce EZs based on predefined parameter for "medium" sized EZs
gridEZ(
  population_raster_path = "example/filepathway/pop_raster.tif",
  strata1_raster_path = "example/filepathway/strata1_raster.tif",
  strata2_raster_path = "example/filepathway/strata2_raster.tif",
  output_path = "example/filepathway/outputfiles/",
  run_ID = "_run1",
  ncores = 3
)

# run gridEZ using predefined parameters for "small" EZs. You can also specify "large" EZs using this code format
gridEZ(
  population_raster_path = "example/filepathway/pop_raster.tif",
  strata1_raster_path = "example/filepathway/strata1_raster.tif",
  strata2_raster_path = "example/filepathway/strata2_raster.tif",
  predefined_EZ_size = TRUE,
  EZ_target_size = "small",
  output_path = "example/filepathway/outputfiles/",
  run_ID = "_run2",
  ncores = 3
)

# run gridEZ using your own parameters by setting predefined_EZ_size as FALSE and then providing numeric values for target_pop_per_EZ and max_cells_per_EZ
gridEZ(
  population_raster_path = "example/filepathway/pop_raster.tif",
  strata1_raster_path = "example/filepathway/strata1_raster.tif",
  strata2_raster_path = "example/filepathway/strata2_raster.tif",
  predefined_EZ_size = FALSE,
  target_pop_per_EZ = 600,
  max_cells_per_EZ = 2000,
  output_path = "example/filepathway/outputfiles/",
  run_ID = "_run2",
  ncores = 3
)
```

## Example output

This is an example section of a gridEZ sampling frame for Samoa. The
target population count per EZ is 500 people, i.e. the gridEZ function’s
“medium” size EZs.

<img src="man/figures/WSM_example_output_medEZ.png" width=500>

Map data: Google, Maxar Technologies ©2025
