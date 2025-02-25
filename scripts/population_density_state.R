# PACKAGES & SETUP -----------------------------------------------------------
pacman::p_load(
  rio,
  here,
  sf,
  reactable,
  htmltools,
  tidyverse,
  knitr,
  htmlwidgets,
  webshot2,
  janitor
)

# IMPORT & CLEAN DATA -----------------------------------------------------------
# 1) Shapefile (district-level, to be aggregated to state)
map_malaysia <- suppressMessages(suppressWarnings(
  st_read("https://raw.githubusercontent.com/dosm-malaysia/data-open/main/datasets/geodata/administrative_2_district.geojson", quiet = TRUE)
)) %>%
  st_make_valid() %>%
  rename(district_code_dosm = code_district) %>%
  mutate(district_code_dosm = as.character(district_code_dosm))

# Aggregate district-level shapefile to state-level boundaries
map_malaysia_state <- map_malaysia %>%
  group_by(state) %>%
  summarise(geometry = st_union(geometry)) %>%
  ungroup()

# 2) District-level population (replace old Facebook data)
district_pop <- import(here("data", "district_population2023.xlsx")) %>%
  clean_names() %>%
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
  ) %>% 
  filter(level == "state") %>% 
  mutate(district = str_to_sentence(district),
         district = recode(district,
                           "Wp kuala lumpur" = "Wilayah Persekutuan Kuala Lumpur",
                           "Wp putrajaya"   = "Wilayah Persekutuan Putrajaya",
                           "Wp labuan" = "Wilayah Persekutuan Labuan",
                           "Negerisembilan" = "Negeri Sembilan",
                           "Pulau pinang" = "Pulau Pinang"
         )) %>% 
  group_by(district) %>%
  mutate(pop_60_plus = if_else(district == "Wilayah Persekutuan Kuala Lumpur & Putrajaya",
                               sum(pop_60_plus, na.rm = TRUE),
                               pop_60_plus)) %>% 
  select(state = district,
         pop_60_plus) %>% 
  distinct()

# 3) DM Prevalence and additional prevalence data by state
prevalence <- import(here("data", "prevalence.csv")) %>%
  clean_names() %>%
  # Convert percentages into proportions (if they are given as percentages)
  mutate(
    diabetes_prevalence_60 = diabetes_prevalence_60*100,
    hpt_prevalence_60      = hpt_prevalence_60_est*100,
    obesity_prevalence_60  = obesity_prevalence_60_est*100
  )
# Expected columns: state, diabetes_prevalence_60, hpt_prevalence_60, obesity_prevalence_60

# 4) Vaccine Slots data (already contains state-level info)
vaccine_slots <- import(here("data", "full_data.csv")) %>%
  clean_names()

# AGGREGATION TO STATE -----------------------------------------------------------
# Aggregate facility-level (vaccine slots) data by state.
state_vax_data <- vaccine_slots %>%
  filter(hf_code != "DB-9927") %>% 
  mutate(total_stok_vaksin = as.numeric(total_stok_vaksin),
         state = str_to_sentence(state),
         state = recode(state,
                        "Negeri sembilan" = "Negeri Sembilan",
                        "Pulau pinang" = "Pulau Pinang",
                        "W.p. Kuala lumpur" = "Wilayah Persekutuan Kuala Lumpur",
                        "W.p. Putrajaya" = "Wilayah Persekutuan Putrajaya",
                        "W.p. Labuan" = "Wilayah Persekutuan Labuan")) %>% 
  group_by(state) %>%
  summarise(
    total_slots     = sum(total_slots, na.rm = TRUE),
    total_available = sum(available_slots, na.rm = TRUE),
    total_stock     = sum(total_stok_vaksin, na.rm = TRUE),
    total_booked    = sum(booked, na.rm = TRUE)) %>%
  ungroup() %>% 
  mutate(total_remaining = total_slots - total_booked,
         perc_uptake = (total_booked/total_slots)*100
  ) 

# MERGE STATE VAX DATA WITH POPULATION & PREVALENCE -------------------------
state_merged <- state_vax_data %>%
  left_join(district_pop, by = "state") %>%
  left_join(prevalence, by = "state")

# MODEL: LM WITH POPULATION 60+ & DM, HPT, OBESITY ---------------------
lm_model <- lm(log(total_booked + 1) ~ pop_60_plus + diabetes_cases_60 + 
                 hpt_cases_60_est + obesity_cases_60_est, 
               data = state_merged)
summary(lm_model)

state_merged <- state_merged %>%
  mutate(predicted_demand = exp(predict(lm_model, newdata = .)) - 1,
         predicted_demand = if_else(predicted_demand < 0, 1, predicted_demand))

total_predicted_demand <- sum(state_merged$predicted_demand, na.rm = TRUE)

state_merged <- state_merged %>%
  mutate(predicted_vaccines = (predicted_demand / total_predicted_demand) * 170000)


# Final dataset -----------------------------------------------------------
final_state_data <- state_merged %>%
  select(
    state,
    total_slots,
    slots = total_available,
    booked = total_booked,
    stock = total_stock,
    predicted_vaccines,
    pop_60_plus,
    diabetes_prevalence_60,
    diabetes_cases_60,
    hpt_prevalence_60,
    hpt_cases_60_est,
    obesity_prevalence_60,
    obesity_cases_60_est
  ) %>% 
  mutate(predicted_vaccines = round(predicted_vaccines, 0)) %>%
  # Calculate additional metrics
  mutate(
    percent_slots_booked = if_else(total_slots > 0, booked / slots, NA_real_),
    supply_diff = predicted_vaccines - stock,
    percent_supply_diff = if_else(stock > 0, supply_diff / stock, 1),
    adjustment = case_when(
      supply_diff < 0 ~ paste("Pull out", abs(supply_diff), "vaccines"),
      supply_diff > 0 ~ paste("Supply", supply_diff, "more vaccines"),
      TRUE ~ "No adjustment"
    )
  ) %>%
  # Combine cases with prevalence: cases followed by prevalence in brackets
  mutate(
    DM = paste0(diabetes_cases_60, " (", diabetes_prevalence_60, "%)"),
    HPT = paste0(hpt_cases_60_est, " (", hpt_prevalence_60, "%)"),
    Obesity = paste0(obesity_cases_60_est, " (", obesity_prevalence_60, "%)")
  ) %>%
  # Select/reorder final columns
  select(
    state,
    slots = total_slots,
    booked,
    percent_slots_booked,
    stock,
    predicted_vaccines,
    supply_diff,
    percent_supply_diff,
    adjustment,
    pop_60_plus,
    DM,
    HPT,
    Obesity
  )

# BUILD INTERACTIVE TABLE WITH REACTABLE ------------------------------------------------
final_table <- reactable(
  final_state_data,
  pagination = FALSE,      # Fit all rows on one page
  filterable = FALSE,      # No search bar at the top (change to TRUE if needed)
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
    pop_60_plus = colDef(
      name = "Pop 60 Plus",
      align = "center",
      format = colFormat(separators = TRUE, digits = 0)
    ),
    DM = colDef(
      name = "DM Cases (Prevalence)",
      align = "center"
    ),
    HPT = colDef(
      name = "HPT Cases (Prevalence)",
      align = "center"
    ),
    Obesity = colDef(
      name = "Obesity Cases (Prevalence)",
      align = "center"
    )
  ),
  theme = reactableTheme(
    style = list(fontFamily = "Calibri, Arial, sans-serif", fontSize = "14px")
  )
)

# Display the interactive table in RStudio
final_table

# SAVE THE REACTABLE OBJECT AS RDS ---------------------------------------------
saveRDS(final_table, file = here("output", "state_vaccine.reactable.rds"))
