# PACKAGES & SETUP -----------------------------------------------------------
pacman::p_load(
  rio,
  here,
  sf,
  janitor,
  tidyverse,
  reactable,
  htmltools
)

# IMPORT DATA ---------------------------------------------------------------
# 1) Shapefile for Malaysia map (district-level)
map_malaysia <- suppressMessages(suppressWarnings(
  st_read("https://raw.githubusercontent.com/dosm-malaysia/data-open/main/datasets/geodata/administrative_2_district.geojson", quiet = TRUE)
))
map_malaysia <- st_make_valid(map_malaysia)

# 2) Facility data (assumes facility code, clinic name, lat/long, state, etc.)
facility_data <- import(here("data", "Facility.csv"))

# 3) Vaccine Slots dataset (facility-level vaccination data)
vaccine_slots <- import(here("data", "full_data.csv"))

# 4) Facebook Population Density dataset (for elderly 60+ counts)
facebook_density <- import(here("data", "mys_elderly_60_plus_2020.csv"))

# 5) Updated Prevalence dataset (with DM, HPT, and Obesity prevalence)
prevalence <- import(here("data", "prevalence.csv")) %>%
  clean_names() %>%
  # Convert prevalence values to percentages (if originally proportions)
  mutate(
    diabetes_prevalence_60 = diabetes_prevalence_60 * 100,
    hpt_prevalence_60      = hpt_prevalence_60_est * 100,
    obesity_prevalence_60  = obesity_prevalence_60_est * 100
  )

# DATA CLEANING -------------------------------------------------------------
# Clean column names for vaccine slots and facility data, aggregate by facility
vaccine_slots_clean <- vaccine_slots %>% 
  clean_names() %>% 
  rename(facility_code = hf_code,
         facility_name = hf_name) %>% 
  mutate(total_stok_vaksin = as.numeric(total_stok_vaksin),
         state = str_to_sentence(state),
         state = recode(state,
                        "Negeri sembilan" = "Negeri Sembilan",
                        "Pulau pinang" = "Pulau Pinang",
                        "W.p. Kuala lumpur" = "Wilayah Persekutuan Kuala Lumpur",
                        "W.p. Putrajaya" = "Wilayah Persekutuan Putrajaya",
                        "W.p. Labuan" = "Wilayah Persekutuan Labuan")) %>% 
  group_by(facility_code, state) %>%
  summarise(
    total_slots     = sum(total_slots, na.rm = TRUE),
    total_available = sum(available_slots, na.rm = TRUE),
    total_stock     = sum(total_stok_vaksin, na.rm = TRUE),
    total_booked    = sum(booked, na.rm = TRUE)
  ) %>%
  mutate(total_remaining = total_slots - total_booked,
         perc_uptake = (total_booked/total_slots)*100
  ) %>%  
  filter(facility_code != "DB-9927")

# Clean up facility data
facility_data_clean <- facility_data %>% clean_names() 

# Merge vaccine_slots with facility data (to obtain latitude, longitude, clinic name, etc.)
vaccine_slots_joined <- vaccine_slots_clean %>%
  left_join(facility_data_clean %>% 
              select(facility_code, latitude, longitude, facility_name),
            by = "facility_code")

# Convert to an sf object (using WGS84)
vaccine_slots_sf <- vaccine_slots_joined %>%
  filter(!is.na(latitude) & !is.na(longitude)) %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# TRANSFORM SPATIAL OBJECTS --------------------------------------------------
# Transform map and vaccine_slots to UTM (EPSG:32647)
map_malaysia_utm <- st_transform(map_malaysia, crs = 32647)
vaccine_slots_sf_utm <- st_transform(vaccine_slots_sf, crs = 32647)

# BUFFERING FOR POPULATION DATA ----------------------------------------------
# Convert Facebook density data to an sf object and transform to UTM
population_sf <- st_as_sf(facebook_density, coords = c("longitude", "latitude"), crs = 4326)
population_sf_utm <- st_transform(population_sf, crs = 32647)

# Create buffers around each facility: 5 km, 10 km, 20 km
buffers_5km  <- st_buffer(vaccine_slots_sf_utm, dist = 5000)
buffers_10km <- st_buffer(vaccine_slots_sf_utm, dist = 10000)
buffers_20km <- st_buffer(vaccine_slots_sf_utm, dist = 20000)

# Function to calculate total elderly population within a buffer
calc_pop_within_buffer <- function(buffer_sf, pop_sf, buffer_label) {
  buffer_join <- st_join(buffer_sf, pop_sf, join = st_intersects)
  pop_summary <- buffer_join %>%
    group_by(facility_code) %>%
    summarise(!!sym(buffer_label) := sum(mys_elderly_60_plus_2020, na.rm = TRUE))
  return(pop_summary)
}

# Calculate population totals for each buffer distance
pop_5km  <- calc_pop_within_buffer(buffers_5km,  population_sf_utm, "pop_5km")
pop_10km <- calc_pop_within_buffer(buffers_10km, population_sf_utm, "pop_10km")
pop_20km <- calc_pop_within_buffer(buffers_20km, population_sf_utm, "pop_20km")

# --- MERGE POPULATION SUMMARIES INTO THE FACILITY-LEVEL SF OBJECT -----------
vaccine_slots_sf_utm <- vaccine_slots_sf_utm %>%
  left_join(st_drop_geometry(pop_5km),  by = "facility_code") %>%
  left_join(st_drop_geometry(pop_10km), by = "facility_code") %>%
  left_join(st_drop_geometry(pop_20km), by = "facility_code")

# --- MERGE IN PREVALENCE DATA -----------------------------------------------
vaccine_slots_sf_utm <- vaccine_slots_sf_utm %>%
  left_join(prevalence %>% 
              select(state, diabetes_prevalence_60, hpt_prevalence_60, obesity_prevalence_60),
            by = "state") %>% 
  mutate(diabetes_cases_60_5km = pop_5km * (diabetes_prevalence_60/100),
         diabetes_cases_60_10km = pop_10km * (diabetes_prevalence_60/100),
         diabetes_cases_60_20km = pop_20km * (diabetes_prevalence_60/100),
         hpt_cases_60_5km = pop_5km * (hpt_prevalence_60/100),
         hpt_cases_60_10km = pop_10km * (hpt_prevalence_60/100),
         hpt_cases_60_20km = pop_20km * (hpt_prevalence_60/100),
         obesity_cases_60_5km = pop_5km * (obesity_prevalence_60/100),
         obesity_cases_60_10km = pop_10km * (obesity_prevalence_60/100),
         obesity_cases_60_20km = pop_20km * (obesity_prevalence_60/100)
  )

# --- PCA FOR DIMENSIONALITY REDUCTION ---------------------------------------
# Select columns for PCA: the 5km, 10km, 20km population and the derived case estimates
pca_data <- vaccine_slots_sf_utm %>%
  ungroup() %>% 
  st_drop_geometry() %>%
  select(pop_5km, pop_10km, pop_20km, 
         diabetes_cases_60_5km, diabetes_cases_60_10km, diabetes_cases_60_20km,
         hpt_cases_60_5km, hpt_cases_60_10km, hpt_cases_60_20km,
         obesity_cases_60_5km, obesity_cases_60_10km, obesity_cases_60_20km)

pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)

# Extract the first principal component as a "need score"
vaccine_slots_sf_utm <- vaccine_slots_sf_utm %>%
  ungroup() %>%
  mutate(need_score = pca_result$x[, 1])

# --- PREPARE FACILITY-LEVEL DATA FOR MODELING -------------------------------
# Drop geometry and ensure total_booked is numeric
facility_model_data <- vaccine_slots_sf_utm %>%
  st_drop_geometry() %>%
  mutate(
    total_booked = as.numeric(total_booked),
    pop_5km = if_else(is.na(pop_5km), 0, pop_5km)
  )

# --- LINEAR MODEL USING THE NEED SCORE FROM PCA ----------------------------
lm_model <- lm(log(total_booked + 1) ~ need_score, data = facility_model_data)
summary(lm_model)

# Predict demand on the log-scale then back-transform
facility_model_data <- facility_model_data %>%
  mutate(
    predicted_demand = exp(predict(lm_model, newdata = .)) - 1,
    predicted_demand = if_else(predicted_demand < 0, 1, predicted_demand)
  )

# Allocate a total of 170,000 vaccines proportionally based on predicted demand
total_predicted_demand <- sum(facility_model_data$predicted_demand, na.rm = TRUE)
facility_model_data <- facility_model_data %>%
  mutate(
    predicted_vaccines = round((predicted_demand / total_predicted_demand) * 170000, 0)
  )

# --- BUILD FACILITY (CLINIC)-LEVEL TABLE ------------------------------------
# Reorder, rename, and calculate additional metrics; remove facility_code
final_facility_data <- facility_model_data %>%
  mutate(
    percent_slots_booked = if_else(total_slots > 0, as.numeric(total_booked) / as.numeric(total_slots), NA_real_),
    supply_diff = predicted_vaccines - as.numeric(total_stock),
    percent_supply_diff = if_else(as.numeric(total_stock) > 0, supply_diff / as.numeric(total_stock), 1),
    adjustment = case_when(
      supply_diff < 0 ~ paste("Pull out", abs(supply_diff), "vaccines"),
      supply_diff > 0 ~ paste("Supply", supply_diff, "more vaccines"),
      TRUE ~ "No adjustment"
    )
  ) %>%
  select(
    facility_name,         # Clinic name
    total_slots,           # Slots
    total_booked,          # Booked
    percent_slots_booked,  # % Booked
    total_stock,           # Stock
    predicted_vaccines,    # Predicted Vaccines
    supply_diff,           # Supply Difference
    percent_supply_diff,   # % Supply Difference
    adjustment,
    pop_5km, pop_10km, pop_20km,
    diabetes_prevalence_60,
    hpt_prevalence_60,
    obesity_prevalence_60,
    diabetes_cases_60_5km, diabetes_cases_60_10km, diabetes_cases_60_20km,
    hpt_cases_60_5km, hpt_cases_60_10km, hpt_cases_60_20km,
    obesity_cases_60_5km, obesity_cases_60_10km, obesity_cases_60_20km
  ) %>%
  rename(
    Clinic = facility_name,
    Slots = total_slots,
    Booked = total_booked,
    `Population > 60 (5km radius)` = pop_5km,
    `Population > 60 (10km radius)` = pop_10km,
    `Population > 60 (20km radius)` = pop_20km,
    `Dm Prevalence` = diabetes_prevalence_60,
    `Hpt Prevalence` = hpt_prevalence_60,
    `Obesity Prevalence` = obesity_prevalence_60,
    `Dm Cases > 60 (5km)` = diabetes_cases_60_5km,
    `Dm Cases > 60 (10km)` = diabetes_cases_60_10km,
    `Dm Cases > 60 (20km)` = diabetes_cases_60_20km,
    `Hpt Cases > 60 (5km)` = hpt_cases_60_5km,
    `Hpt Cases > 60 (10km)` = hpt_cases_60_10km,
    `Hpt Cases > 60 (20km)` = hpt_cases_60_20km,
    `Obesity Cases > 60 (5km)` = obesity_cases_60_5km,
    `Obesity Cases > 60 (10km)` = obesity_cases_60_10km,
    `Obesity Cases > 60 (20km)` = obesity_cases_60_20km
  )

# --- BUILD INTERACTIVE TABLE WITH REACTABLE (WITH SEARCH BARS) -------------
facility_table <- reactable(
  final_facility_data,
  pagination = FALSE,       # Show all rows on one page
  filterable = TRUE,        # Enable a search bar on every column
  defaultSorted = "Clinic",
  defaultColDef = colDef(
    filterable = TRUE,      # Ensure each column is searchable
    style = list(fontFamily = "Calibri, Arial, sans-serif", fontSize = "14px")
  ),
  columns = list(
    Clinic = colDef(
      name = "Clinic",
      minWidth = 200
    ),
    Slots = colDef(
      name = "Slots",
      align = "center",
      format = colFormat(separators = TRUE, digits = 0)
    ),
    Booked = colDef(
      name = "Booked",
      align = "center",
      format = colFormat(separators = TRUE, digits = 0)
    ),
    percent_slots_booked = colDef(
      name = "% Booked",
      align = "center",
      format = colFormat(percent = TRUE, digits = 1),
      style = function(value) {
        if (!is.na(value)) {
          if (value < 0.5) {
            list(background = "#ffcccc", color = "#a10000")
          } else if (value < 0.8) {
            list(background = "#ffffcc", color = "#666600")
          } else {
            list(background = "#ccffcc", color = "#006600")
          }
        }
      }
    ),
    total_stock = colDef(
      name = "Stock",
      align = "center",
      format = colFormat(separators = TRUE, digits = 0)
    ),
    predicted_vaccines = colDef(
      name = "Predicted Vaccines",
      align = "center",
      format = colFormat(separators = TRUE, digits = 0)
    ),
    supply_diff = colDef(
      name = "Supply Difference",
      align = "center",
      format = colFormat(separators = TRUE, digits = 0),
      style = function(value) {
        if (!is.na(value)) {
          if (value < 0) {
            list(background = "#ffcccc", color = "#a10000")
          } else if (value > 0) {
            list(background = "#ccffcc", color = "#006600")
          } else {
            list()
          }
        }
      }
    ),
    percent_supply_diff = colDef(
      name = "% Supply Difference",
      align = "center",
      format = colFormat(percent = TRUE, digits = 1),
      style = function(value) {
        if (!is.na(value)) {
          if (value > 0) {
            list(background = "#ffcccc", color = "#a10000")
          } else if (value < 0) {
            list(background = "#ccffcc", color = "#006600")
          } else {
            list()
          }
        }
      }
    ),
    adjustment = colDef(
      name = "Adjustment",
      minWidth = 180,
      align = "center",
      cell = function(value) {
        if (grepl("Pull out", value)) {
          div(style = "color: #a10000; font-weight: bold;", value)
        } else if (grepl("Supply", value)) {
          div(style = "color: #006600; font-weight: bold;", value)
        } else {
          value
        }
      }
    ),
    # Additional columns moved to the back:
    `Population > 60 (5km radius)` = colDef(
      align = "center",
      format = colFormat(separators = TRUE, digits = 0)
    ),
    `Population > 60 (10km radius)` = colDef(
      align = "center",
      format = colFormat(separators = TRUE, digits = 0)
    ),
    `Population > 60 (20km radius)` = colDef(
      align = "center",
      format = colFormat(separators = TRUE, digits = 0)
    ),
    `Dm Prevalence` = colDef(
      align = "center",
      format = colFormat(digits = 1)
    ),
    `Hpt Prevalence` = colDef(
      align = "center",
      format = colFormat(digits = 1)
    ),
    `Obesity Prevalence` = colDef(
      align = "center",
      format = colFormat(digits = 1)
    ),
    `Dm Cases > 60 (5km)` = colDef(
      align = "center",
      format = colFormat(separators = TRUE, digits = 0)
    ),
    `Dm Cases > 60 (10km)` = colDef(
      align = "center",
      format = colFormat(separators = TRUE, digits = 0)
    ),
    `Dm Cases > 60 (20km)` = colDef(
      align = "center",
      format = colFormat(separators = TRUE, digits = 0)
    ),
    `Hpt Cases > 60 (5km)` = colDef(
      align = "center",
      format = colFormat(separators = TRUE, digits = 0)
    ),
    `Hpt Cases > 60 (10km)` = colDef(
      align = "center",
      format = colFormat(separators = TRUE, digits = 0)
    ),
    `Hpt Cases > 60 (20km)` = colDef(
      align = "center",
      format = colFormat(separators = TRUE, digits = 0)
    ),
    `Obesity Cases > 60 (5km)` = colDef(
      align = "center",
      format = colFormat(separators = TRUE, digits = 0)
    ),
    `Obesity Cases > 60 (10km)` = colDef(
      align = "center",
      format = colFormat(separators = TRUE, digits = 0)
    ),
    `Obesity Cases > 60 (20km)` = colDef(
      align = "center",
      format = colFormat(separators = TRUE, digits = 0)
    )
  ),
  theme = reactableTheme(
    style = list(fontFamily = "Calibri, Arial, sans-serif", fontSize = "14px")
  )
)

# Display the interactive facility (clinic-level) table
facility_table

# --- SAVE THE REACTABLE OBJECT AS RDS ---------------------------------------
saveRDS(facility_table, file = here("output", "clinic_vaccine.reactable.rds"))
