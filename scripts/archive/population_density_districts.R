# PACKAGES & SETUP -----------------------------------------------------------
pacman::p_load(
  rio,
  here,
  sf,
  tidyverse,
  janitor
)

# IMPORT & CLEAN DATA ------------------------------------------------------
# 1) Shapefile
map_malaysia <- suppressMessages(suppressWarnings(
  st_read("https://raw.githubusercontent.com/dosm-malaysia/data-open/main/datasets/geodata/administrative_2_district.geojson", quiet = TRUE)
)) %>%
  st_make_valid() %>%
  rename(district_code_dosm = code_district) %>%
  # NEW: Convert district_code_dosm to character to match joining dataset
  mutate(district_code_dosm = as.character(district_code_dosm))

# 2) District-level population (replace old Facebook data)
district_pop <- import(here("data", "district_population2023.xlsx")) %>%
  clean_names() %>%
  # Convert only numeric columns (adjust these column names as needed)
  mutate(across(
    c(pop, 
      male60_64, male65_69, male70_74, male75_79, male80_84, male_above85,
      female60_64, female65_69, female70_74, female75_79, female80_84, female_above85),
    ~ as.numeric(gsub(",", "", .))
  )) %>%
  mutate(
    pop_60_plus = rowSums(across(
      c(male60_64, male65_69, male70_74, male75_79, male80_84, male_above85,
        female60_64, female65_69, female70_74, female75_79, female80_84, female_above85)
    ), na.rm = TRUE)
  )

# 3) DM Prevalence by state
dm_prevalence_clean <- import(here("data", "dm_prevalence.csv")) %>%
  clean_names()  
# (Expected columns include: state, diabetes_prevalence_60, etc.)

# 4) Vaccine Slots data (which already has district information)
vaccine_slots <- import(here("data", "vaccine_slots.csv")) %>%
  clean_names()

# AGGREGATION TO DISTRICT -------------------------------------------------
# (We assume vaccine_slots already has district-level information, so we group by state & district)
district_vax_data <- vaccine_slots %>%
  group_by(state, district) %>%
  summarise(
    total_slots     = sum(sum_of_total_slots, na.rm = TRUE),
    total_available = sum(count_of_available_slots, na.rm = TRUE),
    total_stock     = sum(total_stok_vaksin, na.rm = TRUE),
    total_booked    = sum(booked_slots, na.rm = TRUE),
    total_remaining = sum(remaining_slots, na.rm = TRUE),
    overall_uptake  = if_else(
      sum(sum_of_total_slots, na.rm = TRUE) > 0,
      sum(booked_slots, na.rm = TRUE) / sum(sum_of_total_slots, na.rm = TRUE),
      NA_real_
    )
  ) %>%
  ungroup()

# MERGE DISTRICT VAX DATA WITH POPULATION ---------------------------------
# We join on "district" (adjust if you have a district code to join on)
district_merged <- district_vax_data %>%
  left_join(district_pop, by = "district")

# NEW: Merge DM prevalence data based on state
district_merged <- district_merged %>%
  left_join(dm_prevalence_clean %>% select(state, diabetes_prevalence_60),
            by = "state") %>%
  mutate(
    # Calculate the estimated number of individuals 60+ with DM for the district.
    # Assumes diabetes_prevalence_60 is a percentage.
    pop_60_dm_est = pop_60_plus * (diabetes_prevalence_60 / 100)
  )

# MERGE WITH SHAPEFILE -------------------------------------------------
# Join the aggregated district data with the district shapefile (using district & district_code_dosm if available)
map_malaysia_district <- map_malaysia %>%
  select(state, district_code_dosm, district, geometry)

# For the join, ensure that the common fields have matching types.
# If district in district_merged is character, then this should work:
district_merged_sf <- map_malaysia_district %>%
  left_join(district_merged, by = "district")

# MODEL: PCA WITH POPULATION 60+ & POPULATION 60+ WITH DM ------------------
# We want to combine two columns: pop_60_plus and pop_60_dm_est
pca_data <- district_merged_sf %>%
  st_drop_geometry() %>%
  select(pop_60_plus, pop_60_dm_est) %>% 
  na.omit()

pca_res <- prcomp(pca_data, center = TRUE, scale. = TRUE)

district_merged_sf <- district_merged_sf %>%
  mutate(pca_score = pca_res$x[,1])

# Shift scores if negative
min_score <- min(district_merged_sf$pca_score, na.rm = TRUE)
if (min_score < 0) {
  district_merged_sf <- district_merged_sf %>%
    mutate(pca_score = pca_score - min_score)
}

total_score <- sum(district_merged_sf$pca_score, na.rm = TRUE)
district_merged_sf <- district_merged_sf %>%
  mutate(predicted_vaccines_model2 = (pca_score / total_score) * 170000)

# FINAL DATASET -----------------------------------------------------------
final_district_data <- district_merged_sf %>%
  st_drop_geometry() %>%
  select(
    district_code_dosm,
    district,
    total_slots,
    total_available,
    total_booked,
    total_remaining,
    overall_uptake,
    pop,
    pop_60_plus,
    diabetes_prevalence_60,
    pop_60_dm_est,
    pca_score,
    predicted_vaccines_model2
  )

print(final_district_data, n = Inf)
