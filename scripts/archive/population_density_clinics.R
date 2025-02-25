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
# 1) Shapefile (for Malaysia map)
map_malaysia <- suppressMessages(suppressWarnings(
  st_read("https://raw.githubusercontent.com/dosm-malaysia/data-open/main/datasets/geodata/administrative_2_district.geojson", quiet = TRUE)
))
map_malaysia <- st_make_valid(map_malaysia)

# 2) Facility data
facility_data <- import(here("data", "Facility.csv"))

# 3) Vaccine Slots dataset
vaccine_slots <- import(here("data", "vaccine_slots.csv"))

# 4) Facebook Population Density dataset
facebook_density <- import(here("data", "mys_elderly_60_plus_2020.csv"))

# 5) DM Prevalence by state
dm_prevalence <- import(here("data", "dm_prevalence.csv"))

# DATA CLEANING -------------------------------------------------------------
vaccine_slots_clean <- vaccine_slots %>% clean_names()
facility_data_clean <- facility_data %>% clean_names()
dm_prevalence_clean <- dm_prevalence %>% clean_names()  
# (Expected: state, diabetes_prevalence_60, etc.)

# Merge vaccine_slots with facility data (to obtain lat/long)
vaccine_slots_joined <- vaccine_slots_clean %>%
  left_join(facility_data_clean %>% select(facility_code, latitude, longitude),
            by = "facility_code")

# Convert to an sf object (WGS84)
vaccine_slots_sf <- vaccine_slots_joined %>%
  filter(!is.na(latitude) & !is.na(longitude)) %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# DATA TRANSFORM -------------------------------------------------------------
# Transform map_malaysia and vaccine_slots_sf to UTM (EPSG:32647)
map_malaysia_utm <- st_transform(map_malaysia, crs = 32647)
vaccine_slots_sf_utm <- st_transform(vaccine_slots_sf, crs = 32647)

# BUFFERING FOR POPULATION DATA ----------------------------------------------
population_sf <- st_as_sf(facebook_density, coords = c("longitude", "latitude"), crs = 4326)
population_sf_utm <- st_transform(population_sf, crs = 32647)

# Create buffers: 5 km, 10 km, 20 km
buffers_5km  <- st_buffer(vaccine_slots_sf_utm, dist = 5000)
buffers_10km <- st_buffer(vaccine_slots_sf_utm, dist = 10000)
buffers_20km <- st_buffer(vaccine_slots_sf_utm, dist = 20000)

# Function: calculate total population within a buffer
calc_pop_within_buffer <- function(buffer_sf, pop_sf, buffer_label) {
  buffer_join <- st_join(buffer_sf, pop_sf, join = st_intersects)
  pop_summary <- buffer_join %>%
    group_by(facility_code) %>%
    summarise(!!sym(buffer_label) := sum(mys_elderly_60_plus_2020, na.rm = TRUE))
  return(pop_summary)
}

# Calculate population totals for each buffer
pop_5km  <- calc_pop_within_buffer(buffers_5km,  population_sf_utm, "pop_5km")
pop_10km <- calc_pop_within_buffer(buffers_10km, population_sf_utm, "pop_10km")
pop_20km <- calc_pop_within_buffer(buffers_20km, population_sf_utm, "pop_20km")

# Merge population summaries back into facility sf object
vaccine_slots_sf_utm <- vaccine_slots_sf_utm %>%
  left_join(st_drop_geometry(pop_5km), by = "facility_code") %>%
  left_join(st_drop_geometry(pop_10km), by = "facility_code") %>%
  left_join(st_drop_geometry(pop_20km), by = "facility_code")

# CALCULATE AT-RISK POPULATION -----------------------------------------------
# Merge DM prevalence by state and calculate estimated 60+ with DM for each buffer
vaccine_slots_sf_utm <- vaccine_slots_sf_utm %>%
  left_join(dm_prevalence_clean %>% select(state, diabetes_prevalence_60),
            by = "state") %>%
  mutate(
    pop_5km_dm_est  = pop_5km  * (diabetes_prevalence_60 / 100),
    pop_10km_dm_est = pop_10km * (diabetes_prevalence_60 / 100),
    pop_20km_dm_est = pop_20km * (diabetes_prevalence_60 / 100)
  )

# MODEL 1: PCA with Total 60+ (All Buffers) -------------------------------
pop_matrix_all <- vaccine_slots_sf_utm %>% 
  st_drop_geometry() %>%
  select(pop_5km, pop_10km, pop_20km)
pca_all <- prcomp(pop_matrix_all, center = TRUE, scale. = TRUE)
vaccine_slots_sf_utm <- vaccine_slots_sf_utm %>%
  mutate(pca_all_score = pca_all$x[,1])
if(min(vaccine_slots_sf_utm$pca_all_score, na.rm = TRUE) < 0) {
  vaccine_slots_sf_utm <- vaccine_slots_sf_utm %>%
    mutate(pca_all_score = pca_all_score - min(pca_all_score, na.rm = TRUE))
}
total_pca_all <- sum(vaccine_slots_sf_utm$pca_all_score, na.rm = TRUE)
vaccine_slots_sf_utm <- vaccine_slots_sf_utm %>%
  mutate(predicted_vaccines_model1 = (pca_all_score / total_pca_all) * 170000)

# MODEL 2: PCA with 60+ with DM (All Buffers) ------------------------------
pop_matrix_dmstate <- vaccine_slots_sf_utm %>% 
  st_drop_geometry() %>%
  select(pop_5km_dm_est, pop_10km_dm_est, pop_20km_dm_est)
pca_dmstate <- prcomp(pop_matrix_dmstate, center = TRUE, scale. = TRUE)
vaccine_slots_sf_utm <- vaccine_slots_sf_utm %>%
  mutate(pca_dmstate_score = pca_dmstate$x[,1])
if(min(vaccine_slots_sf_utm$pca_dmstate_score, na.rm = TRUE) < 0) {
  vaccine_slots_sf_utm <- vaccine_slots_sf_utm %>%
    mutate(pca_dmstate_score = pca_dmstate_score - min(pca_dmstate_score, na.rm = TRUE))
}
total_pca_dmstate <- sum(vaccine_slots_sf_utm$pca_dmstate_score, na.rm = TRUE)
vaccine_slots_sf_utm <- vaccine_slots_sf_utm %>%
  mutate(predicted_vaccines_model2_new = (pca_dmstate_score / total_pca_dmstate) * 170000)

# COMBINE MODELS WITH A LINEAR MODEL ---------------------------------------
combined_df <- vaccine_slots_sf_utm %>%
  st_drop_geometry() %>%
  select(facility_code, predicted_vaccines_model1, predicted_vaccines_model2_new) %>%
  mutate(final_vaccines = (predicted_vaccines_model1 + predicted_vaccines_model2_new) / 2)
lm_final <- lm(final_vaccines ~ predicted_vaccines_model1 + predicted_vaccines_model2_new, data = combined_df)
combined_df <- combined_df %>%
  mutate(final_predicted_vaccines = predict(lm_final, newdata = combined_df))
vaccine_slots_sf_utm <- vaccine_slots_sf_utm %>%
  left_join(combined_df %>% select(facility_code, final_predicted_vaccines), by = "facility_code") %>%
  mutate(final_predicted_vaccines = round(final_predicted_vaccines, 0))

# FINAL STATE DATASET (for state-level table) -------------------------------
final_state_data <- state_merged %>%
  select(
    state,
    total_slots,
    slots = total_available,
    booked = total_booked,
    stock = total_stock,
    remainder = total_remaining,
    pop_60_plus,
    diabetes_prevalence_60,
    predicted_vaccines
  ) %>% 
  mutate(predicted_vaccines = round(predicted_vaccines, 0)) %>%
  mutate(
    percent_slots_booked = if_else(slots > 0, booked / slots, NA_real_),
    supply_diff = predicted_vaccines - stock,
    percent_supply_diff = if_else(stock > 0, supply_diff / stock, 1),
    adjustment = case_when(
      supply_diff < 0 ~ paste("Pull out", abs(supply_diff), "vaccines"),
      supply_diff > 0 ~ paste("Supply", supply_diff, "more vaccines"),
      TRUE ~ "No adjustment"
    ),
    diabetes_prevalence_60 = diabetes_prevalence_60 / 100
  ) %>% 
  select(-c(total_slots, remainder)) %>% 
  select(state, slots, booked, percent_slots_booked, stock, predicted_vaccines, 
         supply_diff, percent_supply_diff, adjustment, pop_60_plus, diabetes_prevalence_60)

# BUILD STATE TABLE WITH REACTABLE -----------------------------------------
final_table <- reactable(
  final_state_data,
  pagination = FALSE,
  filterable = FALSE,
  defaultSorted = "slots",
  defaultColDef = colDef(
    style = list(fontFamily = "Calibri, Arial, sans-serif", fontSize = "14px")
  ),
  columns = list(
    state = colDef(
      name = "State",
      minWidth = 150,
      cell = function(value) {
        img_src <- knitr::image_uri(sprintf("images/%s.png", value))
        tagList(
          img(src = img_src, height = "32px", style = list(marginRight = "8px")),
          div(value, style = list(display = "inline-block"))
        )
      }
    ),
    slots = colDef(
      name = "Slots",
      align = "center",
      format = colFormat(separators = TRUE, digits = 0)
    ),
    booked = colDef(
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
          if (value < 50) {
            list(background = "#ffcccc", color = "#a10000")
          } else if (value < 80) {
            list(background = "#ffffcc", color = "#666600")
          } else {
            list(background = "#ccffcc", color = "#006600")
          }
        }
      }
    ),
    stock = colDef(
      name = "Stock",
      align = "center",
      format = colFormat(separators = TRUE, digits = 0)
    ),
    predicted_vaccines = colDef(
      name = "Predicted Vaccines",
      align = "center",
      format = colFormat(separators = TRUE, digits = 0)
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
    pop_60_plus = colDef(
      name = "Pop 60 Plus",
      align = "center",
      format = colFormat(separators = TRUE, digits = 0)
    ),
    diabetes_prevalence_60 = colDef(
      name = "Dm Prevalence",
      align = "center",
      format = colFormat(percent = TRUE, digits = 1)
    )
  ),
  theme = reactableTheme(
    style = list(fontFamily = "Calibri, Arial, sans-serif", fontSize = "14px")
  )
)

# Wrap state table with caption/footnote for state-level table
state_table_widget <- tagList(
  final_table,
  div(
    style = "font-size: 12px; text-align: center; margin-top: 10px;",
    "Predicted vaccines are based on population estimates (DOSM) and NCD estimates (NHMS 2023).",
    br(),
    "Digital Health Division"
  )
)

# BUILD FACILITY TABLE WITH REACTABLE --------------------------------------
# Create a facility-level dataset: facility name (in sentence case), slots, booked, stock, final predicted vaccines.
final_facility_data <- vaccine_slots_sf_utm %>%
  st_drop_geometry() %>%
  select(facility_code, facility = klinik_kesihatan, slots = count_of_available_slots, 
         booking = booked_slots, stock = total_stok_vaksin, final_predicted_vaccines) %>%
  mutate(facility = str_to_sentence(facility),
         slots = round(slots, 0),
         booking = round(booking, 0),
         stock = round(stock, 0),
         final_predicted_vaccines = round(final_predicted_vaccines, 0))

facility_table <- reactable(
  final_facility_data,
  pagination = FALSE,
  filterable = FALSE,
  defaultSorted = "facility",
  defaultColDef = colDef(
    style = list(fontFamily = "Calibri, Arial, sans-serif", fontSize = "14px")
  ),
  columns = list(
    facility = colDef(
      name = "Facility",
      minWidth = 200
    ),
    slots = colDef(
      name = "Slots",
      align = "center",
      format = colFormat(separators = TRUE, digits = 0)
    ),
    booking = colDef(
      name = "Booked",
      align = "center",
      format = colFormat(separators = TRUE, digits = 0)
    ),
    stock = colDef(
      name = "Stock",
      align = "center",
      format = colFormat(separators = TRUE, digits = 0)
    ),
    final_predicted_vaccines = colDef(
      name = "Predicted Vaccines",
      align = "center",
      format = colFormat(separators = TRUE, digits = 0)
    )
  ),
  theme = reactableTheme(
    style = list(fontFamily = "Calibri, Arial, sans-serif", fontSize = "14px")
  )
)

# Wrap facility table with caption/footnote for facility-level table
facility_table_widget <- tagList(
  facility_table,
  div(
    style = "font-size: 12px; text-align: center; margin-top: 10px;",
    "Predicted vaccines are based on population estimates (DOSM) and NCD estimates (NHMS 2023).",
    br(),
    "Digital Health Division"
  )
)

# DISPLAY TABLES -----------------------------------------------------------
# You can display both tables in your RStudio Viewer or browser.
state_table_widget
facility_table_widget

# SAVE TABLES AS SELF-CONTAINED HTML FILES ---------------------------------
library(htmlwidgets)
saveWidget(state_table_widget, "final_state_table.html", selfcontained = TRUE)
saveWidget(facility_table_widget, "final_facility_table.html", selfcontained = TRUE)
