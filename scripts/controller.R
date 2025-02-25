# run_controller.R
# Controller script to update state and clinic population densities and render dashboard

# Print start message with timestamp
cat("Starting controller script:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# 1. Run the population density script for states
cat("Running population_density_state.R...\n")
source("scripts/population_density_state.R")

# 2. Run the population density script for clinics
cat("Running population_density_clinics.R...\n")
source("scripts/population_density_clinics.R")

# 3. Render the flexdashboard to generate index.html in the docs folder
cat("Rendering state_dashboard.Rmd to docs/index.html...\n")
rmarkdown::render("state_dashboard.Rmd",
                  output_file = "index.html",
                  output_dir = "docs")

cat("Controller script finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
