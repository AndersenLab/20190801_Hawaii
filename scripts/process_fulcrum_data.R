#!/usr/bin/env Rscript
#Load necessary packages
library(tidyverse)
library(lubridate)
library(geosphere)
library(googlesheets)
library(geonames)
library(rebus)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

###################################################################
### Define functions                                            ###
###################################################################

filter_box <- function(longitude, latitude, coords) {
  between(longitude, coords[1], coords[3]) &
    between(latitude, coords[2], coords[4]) &
    !is.na(longitude)
}

FtoC <- function(F) {
  (F - 32)*(5/9)
}

###################################################################
### 1: Read field collection data (C-labels)                    ###
###################################################################
collection <- readr::read_csv("data/fulcrum/nematode_field_sampling.csv") %>%
  dplyr::mutate(c_label = stringr::str_to_upper(c_label)) %>%
  # name created_by to specify who picked up the sample
  dplyr::rename(collection_by = created_by) %>%
  dplyr::select(-updated_at,
                -system_created_at,
                -system_updated_at,
                -date) %>%
  # choose one sample photo only. This takes the first sample photo and warns if additional photos are discarded
  tidyr::separate(col = sample_photo, into = "sample_photo", sep = ",", extra = "warn") %>%
  # this is UTC time (very important if you want to convert to HST time)
  dplyr::mutate(collection_datetime_UTC = lubridate::ymd_hms(created_at, tz = "UTC")) %>% 
  # again this is UTC date (very important if you want to convert to HST date)
  dplyr::mutate(collection_date_UTC = lubridate::date(created_at)) %>% 
  dplyr::select(-created_at) %>%
  # Fix Fahrenheit observations
  dplyr::mutate(substrate_temperature = ifelse(substrate_temperature > 35,
                                               FtoC(substrate_temperature),
                                               substrate_temperature)) %>% 
  # Fix ambient temp F to C
  dplyr::mutate(ambient_temperature = ifelse(ambient_temperature_c > 50,
                                             FtoC(ambient_temperature_c),
                                             ambient_temperature_c)) %>%
  # force ambient temp to numeric
  dplyr::mutate(ambient_temperature = as.numeric(ambient_temperature)) %>%
  # drop ambient temp c
  dplyr::select(-ambient_temperature_c) %>%
  # add flags for runs of temperature data
  dplyr::arrange(collection_datetime_UTC) %>%
  dplyr::mutate(flag_ambient_temperature_run = (ambient_humidity == dplyr::lag(ambient_humidity)) &
                  (ambient_temperature == dplyr::lag(ambient_temperature))
                & (gridsect == "no"))

# OPTIONAL - remove change set variables that we don't need.Should not have exported change set information. And remove collections that are not part of Stefan's collection.
collection <-  collection %>%
  dplyr::select(-version, -change_type, -changeset_id) %>%
  dplyr::filter(collection_by == "stefanzdraljevic2018@u.northwestern.edu")

###################################################################
### 2: Read isolation data (S-labels)                           ###
###################################################################

# Read in S-labels
isolation <- readr::read_csv("data/fulcrum/nematode_isolation.csv") %>%
  dplyr::select(c_label_id = c_label,
                isolation_id = fulcrum_id,
                isolation_datetime_UTC = system_created_at,
                isolation_by = created_by,
                worms_on_sample,
                approximate_number_of_worms,
                isolation_date_UTC = date,
                isolation_local_time = time,
                isolation_latitude = latitude,
                isolation_longitude = longitude)

# OPTIONAL - filter out isolations that are not really part of this project
isolation <- isolation %>%
  dplyr::filter(isolation_by %in% c("claire.buchanan@northwestern.edu", "nicoleroberto2018@u.northwestern.edu"))

#############################################################################
### 3: Use exiftool to extract lat, long elevation. ONLY NEED TO RUN ONCE ###
#############################################################################

# # Read in data from photos. Need to install using ‘brew install exiftool’ in terminal.
# comm <- paste0("exiftool -coordFormat '%+.6f' -csv -ext jpg ",
#                 getwd(),
#                 "/photos/*")
# 
# # Exif Data
#  exif <- readr::read_csv(pipe(comm)) %>%
#    dplyr::mutate(SourceFile = stringr::str_replace(basename(SourceFile), ".jpg", "")) %>%
#    dplyr::select(sample_photo = SourceFile,
#                  altitude = GPSAltitude,
#                  latitude = GPSLatitude,
#                  longitude = GPSLongitude,
#                  ExposureTime,
#                  Artist,
#                  Aperture,
#                  BrightnessValue,
#                  FOV) %>%
#    dplyr::mutate(altitude =  as.numeric(stringr::str_replace(altitude, " m", ""))) %>%
#    dplyr::mutate(FOV =  as.numeric(stringr::str_replace(FOV, " deg", ""))) %>%
#    dplyr::group_by(sample_photo) %>%
#    # Only retain data from one sample photo.
#    dplyr::distinct(.keep_all=T)
# save(file = "data/fulcrum/exif.Rda", exif)

# load data from images already processed by Exif
load("data/fulcrum/exif.Rda")

###################################################################
### 4: Join collection, isolation, and location data            ###
###################################################################

#prevent scientific notation
options(scipen=999)

# join collection, isolation, and location data
df1 <- dplyr::full_join(isolation, collection, by = c("c_label_id" = "fulcrum_id")) %>%
  #rename the lat and long from fulcrum to collection_fulcrum_latitude and collection_fulcrum_longitude so that we can specify lat and long from exif tool
  dplyr::rename(collection_fulcrum_latitude = latitude, collection_fulcrum_longitude = longitude) %>%
  dplyr::select(c_label,
                everything(),
                -c_label_id,
                -sample_photo_url) %>%
  # Join position data from exif by sample_photo. in some cases there is not position data from the photos
  dplyr::left_join(exif) %>%
  # Create flag to track if lat and long come from record or photo
  dplyr::mutate(collection_lat_long_method = ifelse(is.na(latitude), "fulcrum", "photo")) %>%
  # In cases where lat/lon are not available from photo set to collection_fulcrum_latitude and collection_fulcrum_longitude 
  dplyr::mutate(latitude = ifelse(is.na(latitude), collection_fulcrum_latitude, latitude)) %>%
  dplyr::mutate(longitude = ifelse(is.na(longitude), collection_fulcrum_longitude, longitude)) %>%
  #dplyr::mutate_at(.vars = vars(dplyr::starts_with("gps")),
  #                 .funs = funs(as.numeric)) %>%
  dplyr::rename(fulcrum_altitude = gps_altitude) %>%
  dplyr::mutate(worms_on_sample = ifelse(is.na(worms_on_sample), "?", worms_on_sample)) %>%
  dplyr::filter(!is.na(c_label)) %>%
  dplyr::select(-assigned_to,
                -status,
                -Artist) %>%
  # Calculate the Haversine distance between fulcrum record_latitude and record_longitue and photo latitude and longitude
  dplyr::rowwise() %>%
  dplyr::mutate(collection_lat_long_method_diff = geosphere::distHaversine(c(longitude, latitude),
                                                   c(collection_fulcrum_longitude, collection_fulcrum_latitude)),
                # adjust collection_lat_long_method_diff to NA if there is only a fulcrum GPS postion
                collection_lat_long_method_diff = ifelse(collection_lat_long_method == "fulcrum", NA, collection_lat_long_method_diff)) %>%
  # rename the latitude and longitude to include 'collection_' prefix
  dplyr::ungroup() %>%
  dplyr::rename(collection_latitude = latitude,
                collection_longitude = longitude,
                collection_local_time = time)

####################################################################
###      (OPTIONAL) CORRECT DUPLICATE ISOLATIONS (OPTIONAL)      ###
####################################################################

###################################################################
### 5: Joining C_lables with S_labels                           ###
###################################################################

df2 <- readr::read_csv("data/fulcrum/nematode_isolation_s_labeled_plates.csv") %>%
  # OPTIONAL - select claire and nicole records only b/c project was not specified and we had to export data based on record dates.
  dplyr::filter(created_by %in% c("claire.buchanan@northwestern.edu", "nicoleroberto2018@u.northwestern.edu")) %>%
  dplyr::select(fulcrum_parent_id, s_label) %>%
  dplyr::full_join(df1, by = c("fulcrum_parent_id" = "isolation_id")) %>% # this used to be a left join
  dplyr::select(-fulcrum_parent_id, -updated_by, -geometry, -gps_horizontal_accuracy,
                -gps_vertical_accuracy, -gps_speed, -gps_course, -ExposureTime, 
                -Aperture, -BrightnessValue, -FOV, -altitude) %>%
  # set S-labels to NA if isolation entry is 'no plates' there is 1 instance of this.
  dplyr::mutate(s_label = ifelse(s_label == "no plates", NA, s_label))

# OPTIONAL: remove duplicated s_label(s)
df2 <- df2 %>%
  # add a count of row number to grouping variable to remove duplicate s_label (S-0298).
  dplyr::group_by(s_label) %>%
  dplyr::mutate(n = row_number()) %>%
  dplyr::mutate(n = ifelse(is.na(s_label), NA, n)) %>%
  dplyr::filter(is.na(n) | n == "1") %>%
  dplyr::select(-n)

###################################################################
### 6: Handle Substrates                                        ###
###################################################################

# adjust substrate categories if needed

###################################################################
### 7: Get alititudes from ggs locations                        ###
###################################################################
# # only need to run once to get altitudes.This can take a long time ~10 minutes for 1,000 collections.
# options(geonamesUsername="katiesevans")
# altitudes <- df1 %>%
#   dplyr::ungroup() %>%
#   dplyr::select(c_label, collection_latitude, collection_longitude) %>%
#   dplyr::rowwise() %>%
#   # Use collection_latitidue and collection_longitude to find altitudes. Note, these lat and longs should be spot checked to ensure proper collection locations.
#   dplyr::mutate(geonames_altitude = geonames::GNsrtm3(collection_latitude, collection_longitude)$srtm3) %>%
#   dplyr::ungroup()
# 
# save(altitudes, file = "data/fulcrum/altitude.Rda")

# join geonames altitude data to record altitude data
load("data/fulcrum/altitude.Rda")

df3 <- df2 %>%
  dplyr::ungroup() %>%
  dplyr::inner_join(altitudes) %>%
  # make altitude variable and altitude_method variables to track which altitude is being used
  dplyr::mutate(altitude_method_diff = geonames_altitude - fulcrum_altitude,
                altitude = ifelse(collection_lat_long_method == "photo", geonames_altitude,
                                  ifelse(is.na(fulcrum_altitude), geonames_altitude, fulcrum_altitude)),
                altitude_method = ifelse(collection_lat_long_method == "photo", "geonames",
                                         ifelse(is.na(fulcrum_altitude), "geonames", "fulcrum")))

###################################################################
### 8: Assign Islands and trails                                ###
###################################################################

# Create Island Column
df3$collection_island <- "?"
df3[filter_box(df3$collection_longitude, df3$collection_latitude, c(-158.3617,21.1968,-157.5117,21.7931)), "collection_island"] <- "Oahu"
df3[filter_box(df3$collection_longitude, df3$collection_latitude, c(-159.9362, 21.6523, -159.1782, 22.472)), "collection_island"] <- "Kauai"
df3[filter_box(df3$collection_longitude, df3$collection_latitude, c(-157.327, 21.0328, -156.685, 21.2574)), "collection_island"] <- "Molokai"
df3[filter_box(df3$collection_longitude, df3$collection_latitude, c(-156.7061, 20.4712, -155.9289, 21.0743)), "collection_island"] <- "Maui"
df3[filter_box(df3$collection_longitude, df3$collection_latitude, c(-156.1346, 18.6619, -154.6985, 20.4492)), "collection_island"] <- "Big Island"

# Create Trail Column
df3$collection_location <- NA
df3[filter_box(df3$collection_longitude, df3$collection_latitude, c(-157.72537,21.303309,-157.71919,21.32122)), "collection_location"] <- "Kuliouou Ridge Trail"
df3[filter_box(df3$collection_longitude, df3$collection_latitude, c(-158.0192352613,21.5014265529,-158.0145925283,21.5041245046)), "collection_location"] <- "Wahiawa Botanical Garden"
df3[filter_box(df3$collection_longitude, df3$collection_latitude, c(-157.8598800302,21.3149311581,-157.855797708,21.3182194587)), "collection_location"] <- "Foster Community Garden"
df3[filter_box(df3$collection_longitude, df3$collection_latitude, c(-157.7829487403,21.3569863645,-157.7752268314,21.3655295525)), "collection_location"] <- "Maunawili Demonstration Trail"
df3[filter_box(df3$collection_longitude, df3$collection_latitude, c(-157.8014534712,21.3322593,-157.798127532,21.3427719396)), "collection_location"] <- "Manoa Falls Trail"
df3[filter_box(df3$collection_longitude, df3$collection_latitude, c(-157.8135502338,21.3779082884,-157.7915561199,21.3970691079)), "collection_location"] <- "Ho'omaluhia Botanical Garden"
df3[filter_box(df3$collection_longitude, df3$collection_latitude, c(-159.613624,22.167098,-159.575601,22.226422)), "collection_location"] <- "Na Pali Coast State Wilderness Park"

###################################################################
### 9: Add photo URLs                                           ###
###################################################################

# df <-df %>% dplyr::rowwise() %>%
#   dplyr::group_by(c_label) %>%
#   dplyr::mutate(photo = paste0(c_label,
#                                ".",
#                                stringr::str_to_lower(stringr::str_replace_all(substrate, "[^[:alnum:]]", "_")),
#                                ".1.jpg"),
#                 photo_url = paste0("https://storage.googleapis.com/elegansvariation.org/photos/hawaii2017/",
#                                    c_label,
#                                    ".jpg"),
#                 photo_url_thumb =  paste0("https://storage.googleapis.com/elegansvariation.org/photos/hawaii2017/",
#                                           c_label,
#                                           ".thumb.jpg")) %>%
#   dplyr::ungroup()
# 
# photo_comms <- df %>% dplyr::mutate(sample_photo = str_split(sample_photo, ",")) %>%
#   dplyr::select(-s_label) %>%
#   dplyr::select(c_label, sample_photo, substrate) %>%
#   tidyr::unnest() %>%
#   dplyr::group_by(c_label) %>%
#   dplyr::mutate(comm = paste0("cp ../data/photos/id/",
#                               sample_photo,
#                               ".jpg",
#                               " ",
#                               "../data/photos/c/",
#                               c_label,
#                               ".",
#                               stringr::str_to_lower(str_replace_all(substrate, "[^[:alnum:]]", "_")),
#                               ".",
#                               dplyr::row_number(c_label),
#                               ".jpg")) %>%
# 
#   dplyr::select(-c_label, comm)
# 
# writeLines(photo_comms$comm, con = file("scripts/rename_photos.sh"))

###################################################################
### 10: Merge automated blast data                              ###
###################################################################

# # Merge in blast data; Take top hit
# blast_results <- readr::read_tsv("data/sanger/blast_results.tsv") %>%
#   dplyr::group_by(s_plate) %>%
#   dplyr::filter(row_number() == 1)

###################################################################
### 11: Merge manual blast results                              ###
###################################################################

# Each collection should have unique gs_key. Your Google Sheet key can be found in your PUBLISHED Google sheets URL.
# Select the string of data found between the slashes after spreadsheets/d in your Google Sheet URL.
# use stringr to find s_labels from s_label column

genotyping_sheet_raw <- googlesheets::gs_key("14bZVya6eWSOZq_Ldth1wovt_fVmiOzWQW74L_l24Gtg") %>%
  googlesheets::gs_read("Sheet1", na = c("#N/A", "NA", ""),
                        by = c("c_label", "s_label")) %>%
  dplyr::filter(!is.na(s_label)) %>%
  # remove c_label variable (this column was hand typed and contains at least 2 errors)
  dplyr::select(s_label, pcr_date, lysis_date, ITS2_pcr_product, notes, species_id) %>%
  # OPTIONAL - add 0s to ITS2_pcr_product variable where NA
  dplyr::mutate(ITS2_pcr_product = ifelse(is.na(ITS2_pcr_product), 0, ITS2_pcr_product))

# find s_labels in genotyping sheet
slabels <- str_subset(genotyping_sheet_raw$s_label, pattern = "S-")

# filter genotyping sheet by s_labels matching "S-" pattern
genotyping_sheet <- genotyping_sheet_raw %>%
  dplyr::filter(s_label %in% slabels)

# Remove any duplicated s_labels. We removed duplicated s_labels in genotyping sheet (S-0374, S-0382, S-0381).

# Join genotyping sheet with collection and isolation data
fulcrum_dat <- df3 %>% 
  dplyr::full_join(genotyping_sheet) %>%
  # Rename variables
  dplyr::rename(project_id = project,
                collection_id = c_label,
                isolation_id = s_label) %>%
  # Fill project_id variable incase there are NAs introduced
  tidyr::fill(project_id) %>%
  # OPTIONAL - change project_id to `HawaiiAugust2019` b/c this was made without a unique fulcrum project_id and Stefan entered them under the HawaiiOctober2018 project.
  dplyr::mutate(project_id = ifelse(!is.na(project_id), "HawaiiAugust2019", "HawaiiAugust2019")) %>%
  # Reorder variables
  dplyr::select(project_id,
                collection_id,
                isolation_id,
                species_id,
                collection_by,
                collection_datetime_UTC,
                collection_date_UTC,
                collection_local_time,
                collection_island,
                collection_location,
                collection_latitude,
                collection_longitude,
                collection_fulcrum_latitude,
                collection_fulcrum_longitude,
                collection_lat_long_method,
                collection_lat_long_method_diff,
                ambient_temperature,
                flag_ambient_temperature_run,
                ambient_humidity,
                substrate_temperature,
                substrate_moisture,
                fulcrum_altitude,
                geonames_altitude,
                altitude,
                altitude_method,
                altitude_method_diff,
                landscape,
                sky_view,
                substrate,
                substrate_other,
                substrate_notes,
                sample_photo,
                sample_photo_caption,
                gridsect,
                gridsect_index,
                grid_sect_direction,
                gridsect_radius,
                isolation_by,
                isolation_datetime_UTC,
                isolation_date_UTC,
                isolation_local_time,
                isolation_latitude,
                isolation_longitude,
                worms_on_sample,
                approximate_number_of_worms,
                lysis_date,
                pcr_date,
                ITS2_pcr_product,
                notes)

# export R dataframe
save(file = "data/fulcrum/fulcrum_dat.Rda", fulcrum_dat)

###################################################################
### 12: project specific report                                 ###
###################################################################

######################
### Issue 1        ###
######################
# Stefan did not have a project ID so he put all his collections in under HawaiiOctober2018 project.
# SOLUTION: exporting records that were created in august and september of 2019 then filtereing out the collections not made by Stefan.

#####################
### Issue 2       ###
#####################
# Because there was no project_ID Nicole did not enter a project ID when helping Claire isolate nematode. 
# SOLUTION: export isolation records made from august to september then filter to Claire and Nicole as isolators. Loraina, Gaotian, and Emily all siolated during this time as well so they were filtered out.
