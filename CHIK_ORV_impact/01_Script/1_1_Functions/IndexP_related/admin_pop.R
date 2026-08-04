#https://gadm.org/download_country.html

## extract population in brazil-------------------------------------------------------

setwd("C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine")
foi_df <- readRDS("CHIK_mapping/MainData/foi_shrink_df.RDS")
gadm_bra <- st_read("CHIK_vaccine_impact/1_Data/world_admin_boundary/gadm41_BRA_1.shp")
foi_sf <- st_as_sf(foi_df, coords = c("x", "y"), crs = 4326)
brazil_sf <- foi_sf %>% filter(country == "Brazil")
brazil_sf_with_admin <- st_join(brazil_sf, gadm_bra, join = st_intersects)
unmatched_points <- brazil_sf_with_admin %>%
  filter(is.na(NAME_1))
# Filter unmatched points
unmatched_points <- brazil_sf_with_admin %>% filter(is.na(NAME_1))

# Find the nearest polygon for each unmatched point
nearest_indices <- st_nearest_feature(unmatched_points, gadm_bra)

# Extract the corresponding NAME_1 values
nearest_names <- gadm_bra$NAME_1[nearest_indices]

# Assign the nearest NAME_1 to unmatched points
unmatched_points$NAME_1 <- nearest_names
brazil_sf_with_admin$ID <- seq_len(nrow(brazil_sf_with_admin))
# Update the original sf object with corrected NAME_1
brazil_sf_with_admin <- brazil_sf_with_admin %>%
  mutate(NAME_1 = ifelse(is.na(NAME_1), nearest_names[ID], NAME_1))
bra_rio <- brazil_sf_with_admin %>% filter(NAME_1 == "Rio de Janeiro")
# Remove geometry and sum columns 8 to 25, ignoring NA values
rio_pop <- colSums(st_drop_geometry(bra_rio[, 8:25]), na.rm = TRUE)

brazil_pop_all <-   brazil_sf_with_admin %>%
  group_by(NAME_1) %>%
  summarise(across(c("1":"18"), sum, na.rm = TRUE)) %>%
  ungroup()
brazil_pop_all <- st_drop_geometry(brazil_pop_all)
brazil_pop_all <- brazil_pop_all[-28, ]
brazil_pop_transposed <- as.data.frame(t(brazil_pop_all))
colnames(brazil_pop_transposed) <- brazil_pop_transposed[1, ]
brazil_pop_transposed <- brazil_pop_transposed[-1, ]
brazil_pop_transposed[] <- as.data.frame(
  lapply(brazil_pop_transposed, function(x) as.numeric(unlist(x)))
)
brazil_pop_transposed <- brazil_pop_transposed %>% mutate(tot_pop = rowSums(select(., everything())))

save(brazil_pop_transposed, file = "00_Data/0_2_Processed/brazil_pop_transposed.RData")

## extract population in costa rica-------------------------------------------------------
gadm_cri <- st_read("CHIK_vaccine_impact/1_Data/world_admin_boundary/gadm41_CRI_1.shp")
foi_sf <- st_as_sf(foi_df, coords = c("x", "y"), crs = 4326)
cri_sf <- foi_sf %>% filter(country == "Costa Rica")
cri_sf_with_admin <- st_join(cri_sf, gadm_cri, join = st_intersects)
unmatched_points <- cri_sf_with_admin %>%
  filter(is.na(NAME_1))
# Filter unmatched points
unmatched_points <- cri_sf_with_admin %>% filter(is.na(NAME_1))

# Find the nearest polygon for each unmatched point
nearest_indices <- st_nearest_feature(unmatched_points, gadm_cri)

# Extract the corresponding NAME_1 values
nearest_names <- gadm_cri$NAME_1[nearest_indices]

# Assign the nearest NAME_1 to unmatched points
unmatched_points$NAME_1 <- nearest_names
cri_sf_with_admin$ID <- seq_len(nrow(cri_sf_with_admin))
# Update the original sf object with corrected NAME_1
cri_sf_with_admin <- cri_sf_with_admin %>%
  mutate(NAME_1 = ifelse(is.na(NAME_1), nearest_names[ID], NAME_1))

cri_pop_all <-   cri_sf_with_admin %>%
  group_by(NAME_1) %>%
  summarise(across(c("1":"18"), sum, na.rm = TRUE)) %>%
  ungroup()
cri_pop_all <- st_drop_geometry(cri_pop_all)
cri_pop_all <- cri_pop_all[-8, ]
cri_pop_transposed <- as.data.frame(t(cri_pop_all))
colnames(cri_pop_transposed) <- cri_pop_transposed[1, ]
cri_pop_transposed <- cri_pop_transposed[-1, ]
cri_pop_transposed[] <- as.data.frame(
  lapply(cri_pop_transposed, function(x) as.numeric(unlist(x)))
)

cri_pop_transposed <- cri_pop_transposed %>% mutate(tot_pop = rowSums(select(., everything())))

save(cri_pop_transposed, file = "00_Data/0_2_Processed/cri_pop_transposed.RData")

## extract population in dom-------------------------------------------------------

foi_df <- readRDS("00_Data/0_1_Raw/foi_shrink_df.RDS")
gadm_dom <- st_read("00_Data/0_1_Raw/world_admin_boundary/gadm41_DOM_1.shp")
foi_sf <- st_as_sf(foi_df, coords = c("x", "y"), crs = 4326)
dom_sf <- foi_sf %>% filter(country == "Dominican Republic")
dom_sf_with_admin <- st_join(dom_sf, gadm_dom, join = st_intersects)
unmatched_points <- dom_sf_with_admin %>%
  filter(is.na(NAME_1))
# Filter unmatched points
unmatched_points <- dom_sf_with_admin %>% filter(is.na(NAME_1))

# Find the nearest polygon for each unmatched point
nearest_indices <- st_nearest_feature(unmatched_points, gadm_dom)

# Extract the corresponding NAME_1 values
nearest_names <- gadm_dom$NAME_1[nearest_indices]

# Assign the nearest NAME_1 to unmatched points
unmatched_points$NAME_1 <- nearest_names
dom_sf_with_admin$ID <- seq_len(nrow(dom_sf_with_admin))
# Update the original sf object with corrected NAME_1
dom_sf_with_admin <- dom_sf_with_admin %>%
  mutate(NAME_1 = ifelse(is.na(NAME_1), nearest_names[ID], NAME_1))

dom_pop_all <- dom_sf_with_admin %>%
  group_by(NAME_1) %>%
  summarise(across(c("1":"18"), sum, na.rm = TRUE)) %>%
  ungroup()

dom_pop_all <- st_drop_geometry(dom_pop_all)
dom_pop_all <- dom_pop_all[-27, ]
dom_pop_transposed <- as.data.frame(t(dom_pop_all))
colnames(dom_pop_transposed) <- dom_pop_transposed[1, ]
dom_pop_transposed <- dom_pop_transposed[-1, ]
dom_pop_transposed[] <- as.data.frame(
  lapply(dom_pop_transposed, function(x) as.numeric(unlist(x)))
)

dom_pop_transposed <- dom_pop_transposed %>% mutate(tot_pop = rowSums(select(., everything())))

save(dom_pop_transposed, file = "00_Data/0_2_Processed/dom_pop_transposed.RData")









