
# Load packages -----------------------------------------------------------
pacman::p_load(
  rio,         # for import
  here,        # for file paths
  janitor,     # for clean_names
  tidyverse    # for dplyr, etc.
)


# Import the new dataset ----------------------------------------------
#import new slots
updated_slots <- import(here("data", "updated_slots.csv"))

# Import and Clean Names
stocks <- import(here("data", "stock.csv")) %>%
  clean_names()


# Group by facility code & date, then sum columns of interest ----------
# We'll rename 'hf_code' to 'facility_code' to match the old data's "FACILITY CODE".
# Adjust the columns you want to sum as needed (e.g. total_slots, booked, cancelled, noshow).
updated_slots <- updated_slots %>%
  clean_names() %>% 
# Columns now include: appointment_booking_summary_id, date, hf_code, hf_name, state,
# total_slots, available_slots, booked, cancelled, noshow, uptake, remaining, stock, etc.
  rename(facility_code = hf_code,
         klinik_kesihatan = hf_name) %>%
  group_by(facility_code,
           klinik_kesihatan,
           state) %>%
  summarise(
    total_slots_new  = sum(total_slots, na.rm = TRUE),
    total_booked_new = sum(booked, na.rm = TRUE),
    total_cancelled  = sum(cancelled, na.rm = TRUE),
    total_noshow     = sum(noshow, na.rm = TRUE),
  ) %>% 
  mutate(klinik_kesihatan = str_to_upper(klinik_kesihatan))

#    Also keep negeri, pkd if you want to preserve them.
stocks_selected <- stocks %>%
  select(
    negeri,
    pkd,
    clinic_1 = klinik_kesihtan,
    stock_1 = bilangan_stok,
    clinic_2 = klinik_kesihatan_belum_buka_slot,
    stock_2 = bilangan_stok_2
  )

# 3) Create Two Data Frames and Stack Them
df1 <- stocks_selected %>%
  select(negeri, pkd, clinic = clinic_1, stock = stock_1)

df2 <- stocks_selected %>%
  select(negeri, pkd, clinic = clinic_2, stock = stock_2)

# Bind rows to stack them
stocks_stacked <- bind_rows(df1, df2) %>%
  # Remove rows where 'clinic' is NA or blank
  filter(!is.na(clinic) & clinic != "") %>%
  # 4) Convert clinic names to uppercase
  mutate(clinic = str_to_upper(clinic),
         negeri = str_to_upper(negeri)) %>% 
  rename(klinik_kesihatan = clinic,
         stok = stock) %>% 
  mutate_if(is.character, ~ gsub("[0-9]", "", .)) %>% 
  mutate(klinik_kesihatan = na_if(klinik_kesihatan, "")) %>% 
  drop_na(klinik_kesihatan) %>% 
  mutate(klinik_kesihatan = str_squish(klinik_kesihatan),
         negeri = str_squish(negeri),
         klinik_kesihatan = str_replace_all(klinik_kesihatan, "KLINIK KESIHATAN IBU DAN ANAK", "KKIA"),
         negeri = na_if(negeri, "")) %>%  
  fill(negeri, .direction = "down") 
  

# Now 'stocks_stacked' has columns:
#   negeri, pkd, clinic, stock
# with clinic names in uppercase.

# 4) Merge with the original vaccine slots data ---------------------------


# If the old data doesn't have date, you can omit 'date' from the join.
vaccine_slots <- updated_slots %>%
  rename(negeri = state) %>% 
  mutate(klinik_kesihatan = str_squish(klinik_kesihatan)) %>% 
  left_join(stocks_stacked, by = c("klinik_kesihatan", "negeri"))


vaccine_slots_check <- vaccine_slots %>% 
  filter(is.na(stok))


anti_join(updated_slots, stocks_stacked, by = "klinik_kesihatan")

x <- updated_slots %>%  select(klinik_kesihatan, negeri = state)
y <- stocks_stacked %>%  select(klinik_kesihatan, negeri)

write.csv(x, "base.csv")
write.csv(y, "right.csv")
