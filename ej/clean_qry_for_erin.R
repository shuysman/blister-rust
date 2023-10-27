library(tidyverse)
library(lubridate)

data <- read_csv("./Qry_forErin.csv")

data %>%
    mutate(Start_Date = strptime(Start_Date, format = "%m/%d/%Y %H:%M:%S")) %>%
    transmute(
        network = "fivepks",
        transect_id = Loc_Name,
        year = year(Start_Date),
        aspect_deg = NA,
        slope_deg = NA,
        elevation_m = NA,
        tree_id = TreeID_Number,
        tree_status = Tree_Status,
        tree_species = Species_Code,
        dbh_cm = TreeDBH_cm,
        br_status = if_else(
            if_any(
                ends_with(
                    c(
                        "A_Upper_YN",
                        "A_Mid_YN",
                        "A_Lower_YN",
                        "I_Upper_YN",
                        "I_Mid_YN",
                        "I_Lower_YN"
                    )
                ), ~ .x == "Y"
            ),
            1, 0
        ),
        mpb_status = NA,
        notes = ""
    ) %>%
    view()
