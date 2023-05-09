Prepare_User_Extrapolation_Data_Fn<- function (input_grid, strata.limits = NULL, projargs = NA, zone = NA,  flip_around_dateline = TRUE, ...) 
{

    if (FALSE) {
        input_grid = extrap_df
        strata.limits = data.frame("STRATA" = c("All", "NMFS", "DFO"))
        projargs = NA
        zone = NA
        flip_around_dateline = TRUE
    }
    
    if (is.null(strata.limits)) {
        strata.limits = data.frame(STRATA = "All_areas")
    }
    message("Using strata ", strata.limits)
    Data_Extrap <- input_grid
    Area_km2_x = Data_Extrap[, "Area_km2"]
    Tmp = cbind(BEST_LAT_DD = Data_Extrap[, "Lat"], BEST_LON_DD = Data_Extrap[,  "Lon"])
    if ("Depth" %in% colnames(Data_Extrap)) {
        Tmp = cbind(Tmp, BEST_DEPTH_M = Data_Extrap[, "Depth"])
    }
    if("STRATA" %in% colnames(Data_Extrap)){
        Tmp = cbind(Tmp, BEST_STRATA = as.character(Data_Extrap[, "STRATA"]))
    }
    a_el = as.data.frame(matrix(NA, nrow = nrow(Data_Extrap), ncol = nrow(strata.limits), dimnames = list(NULL, strata.limits[, "STRATA"])))
    for (l in 1:ncol(a_el)) {
        a_el[, l] = apply(Tmp, MARGIN = 1, FUN = match_strata_fn, 
            strata_dataframe = strata.limits[l, , drop = FALSE])
        a_el[, l] = ifelse(is.na(a_el[, l]), 0, Area_km2_x)
    }
    tmpUTM = project_coordinates(X = Data_Extrap[, "Lon"], Y = Data_Extrap[, 
        "Lat"], projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline)
    Data_Extrap = cbind(Data_Extrap, Include = 1)
    if (all(c("E_km", "N_km") %in% colnames(Data_Extrap))) {
        Data_Extrap[, c("E_km", "N_km")] = tmpUTM[, c("X", "Y")]
    } else {
        Data_Extrap = cbind(Data_Extrap, E_km = tmpUTM[, "X"], N_km = tmpUTM[, "Y"])
    }
    Return = list(a_el = a_el, Data_Extrap = Data_Extrap, zone = attr(tmpUTM, 
        "zone"), projargs = attr(tmpUTM, "projargs"), flip_around_dateline = flip_around_dateline, 
        Area_km2_x = Area_km2_x)
    return(Return)
}

match_strata_fn<-  function (x, strata_dataframe)  {
    match_latitude_TF = match_longitude_TF = match_depth_TF = match_strata_TF = rep(TRUE, nrow(strata_dataframe))
    if (all(c("south_border", "north_border") %in% names(strata_dataframe))) {
        match_latitude_TF = as.numeric(x["BEST_LAT_DD"]) > strata_dataframe[, "south_border"] & as.numeric(x["BEST_LAT_DD"]) <=  strata_dataframe[, "north_border"]
    }
    if (all(c("west_border", "east_border") %in% names(strata_dataframe))) {
        match_longitude_TF = as.numeric(x["BEST_LON_DD"]) > strata_dataframe[, "west_border"] & as.numeric(x["BEST_LON_DD"]) <=  strata_dataframe[, "east_border"]
    }
    if (all(c("shallow_border", "deep_border") %in% names(strata_dataframe))) {
        match_depth_TF = as.numeric(x["BEST_DEPTH_M"]) > strata_dataframe[, "shallow_border"] & as.numeric(x["BEST_DEPTH_M"]) <= strata_dataframe[, "deep_border"]
    }
    if(names(strata_dataframe) == "STRATA"){
        match_strata_TF = as.character(x["BEST_STRATA"]) == strata_dataframe[, "STRATA"]
    }
    Char = as.character(strata_dataframe[match_latitude_TF & match_longitude_TF &  match_depth_TF & match_strata_TF, "STRATA"])
    return(ifelse(length(Char) == 0, NA, Char))
}

cog_from_dens<- function(vast_fit, proj_dens, proj_ind){

    if (FALSE) {
        # vast_fit = vast_fit
        # proj_dens = dens_df
        # proj_ind = ind_df
    }
    
    # Check units
    units(proj_dens$D_gct) <- units(vast_fit$Report$D_gct)
    units(proj_ind$Index)<- units(vast_fit$Report$Index_ctl)
    
    ## Calculate Index_gctl (D_gct * a_gl)
    # Get extrapolation info
    if (vast_fit$data_list$n_l > 1) {
        extrap_dat <- data.frame(vast_fit$extrapolation_list$Data_Extrap, vast_fit$extrapolation_list$a_el) %>%
            filter(., Region == "All") %>%
            rename(., "a_gl" = "All")
        proj_ind <- proj_ind %>%
            filter(., Region == "All")
    } else {
        extrap_dat <- data.frame(vast_fit$extrapolation_list$Data_Extrap, "a_gl" = vast_fit$extrapolation_list$a_el)
    }
    
    # Join to area and multiply --- this is different!
    if(FALSE){
        check <- data.frame(vast_fit$Report$Index_gctl[, , 1, ]) 
        index_gctl <- proj_dens %>%
            left_join(., extrap_dat) %>%
            drop_na() %>%
            mutate(.,
                "Index_gctl" = D_gct * a_gl
            )
            mutate("RowID" = seq(from = 1, to = nrow(.))) %>%
            pivot_longer(., !RowID, names_to = "Strata", values_to = "Index_gctl") %>%
            filter(., Strata == "Stratum_1")
        str(check)
    }

    index_gctl <- proj_dens %>%
        left_join(., extrap_dat) %>%
        drop_na() %>%
        mutate(.,
            "Index_gctl" = D_gct * a_gl
        )
    
    # Join to overall index ctl
    index_gctl <- index_gctl %>%
        left_join(., proj_ind) %>%
        mutate(., "Index_Prop" = Index_gctl / Index) %>%
        distinct(., Lat, Lon, Time, Sim_Scenario, D_gct, E_km, N_km, Index_gctl, Index, Index_Prop) %>%
        mutate(.,
            "Lon_wt" = E_km * Index_Prop,
            "Lat_wt" = N_km * Index_Prop
        )
    
    # Now get COG measures by multiplying Index_Prob by lon/lat
    cog_out <- index_gctl %>%
        group_by(Time, Sim_Scenario) %>%
        summarize(.,
            "Mean_Lon" = sum(Lon_wt),
            "Mean_Lat" = sum(Lat_wt)
        )
    return(cog_out)
}


eff_area_from_dens<- function(vast_fit, proj_dens, proj_ind){

    if (FALSE) {
        # vast_fit = vast_fit
        # proj_dens = dens_df
        # proj_ind = ind_df
    }
    
    # Check units
    proj_dens <- proj_dens %>%
        drop_units()
    proj_ind <- proj_ind %>%
        drop_units
    
    # One region?
    if(vast_fit$data_list$n_l == 1){
        # Getting mean density
        # Extrapolation info
        extrap_dat <- data.frame(vast_fit$extrapolation_list$Data_Extrap, "a_gl" = vast_fit$extrapolation_list$a_el)
        
        # First, need index_gctl - density multiplied by area of knot
        index_gctl <- proj_dens %>%
            left_join(., extrap_dat, by = c("Lon" = "Lon", "Lat" = "Lat")) %>%
            drop_na() %>%
            mutate(.,
                "Index_gctl" = D_gct * a_gl
            )
        
        # Next, need information on the total index within the area
        index_gctl <- index_gctl %>%
            left_join(., proj_ind)
        
        # Finally, get the mean density and then use that to get the effective area occupied Summarize across knots and get effective area
        eff_area <- index_gctl %>%
            mutate(., "mean_d_ctl_temp" = D_gct * Index_gctl / Index) %>%
            group_by(., Time, Sim_Scenario, Region) %>%
            summarize(.,
                "mean_d_ctl" = sum(mean_d_ctl_temp, na.rm = TRUE)
            ) %>%
            left_join(., proj_ind) %>%
                mutate(., "Eff_Area" = Index / mean_d_ctl) %>%
                dplyr::select(., Time, Sim_Scenario, Region, Eff_Area)
    } else {
        # Multiple regions...
        regs_vec <- as.vector(unlist(vast_fit$settings$strata.limits))
        
        for(i in seq_along(regs_vec)){
            # Extrapolation info
            extrap_dat <- data.frame(vast_fit$extrapolation_list$Data_Extrap, vast_fit$extrapolation_list$a_el) %>%
                filter(., Region == regs_vec[i]) %>%
                rename(., "a_gl" = regs_vec[i])
            proj_ind_use2 <- proj_ind %>%
                filter(., Region == regs_vec[i])
            
            # First, need index_gctl - density multiplied by area of knot
            index_gctl <- proj_dens %>%
                left_join(., extrap_dat, by = c("Lon" = "Lon", "Lat" = "Lat")) %>%
                drop_na() %>%
                mutate(.,
                    "Index_gctl" = D_gct * a_gl
                )
            
            # Next, need information on the total index within the area
            index_gctl <- index_gctl %>%
                left_join(., proj_ind_use2)
                   
            # Finally, get the mean density and then use that to get the effective area occupied Summarize across knots and get effective area
            eff_area_temp <- index_gctl %>%
                mutate(., "mean_d_ctl_temp" = D_gct * Index_gctl / Index) %>%
                group_by(., Time, Sim_Scenario, Region) %>%
                summarize(.,
                    "mean_d_ctl" = sum(mean_d_ctl_temp, na.rm = TRUE)
                ) %>%
                left_join(., proj_ind_use2) %>%
                    mutate(., "Eff_Area" = Index / mean_d_ctl) %>%
                    dplyr::select(., Time, Sim_Scenario, Region, Eff_Area)
                
            if(i == 1){
                eff_area<- eff_area_temp
            } else {
                eff_area<- bind_rows(eff_area, eff_area_temp)
            }
        }
    }
    return(eff_area)
}

summary.sim_results<- function (vast_fit, sim_obj, resp_scale, nice_times = NULL, out_t_scale = NULL, nice_category_names = nice_category_names, climate_scenario = climate_scenario, 
    out_dir) {
    if (FALSE) {
        # tar_load(vast_fit)
        # tar_load(vast_projections)
        # sim_obj <- vast_projections
        # what <- "Index_ctl"
        # nice_times <- nice_times
        # out_t_scale <- "annual"
        # probs <- c(0.1, 0.5, 0.9)
        # mean_instead <- FALSE
        # nice_category_names <- nice_category_names
        # climate_scenario <- climate_scenario
        # out_dir <- paste0(res_root, "prediction_df")
        
        # Capelin
        date_dir<- here::here("2023-02-17/Capelin_BC/")
        vast_fit = fit_full
        sim_obj = uncert_res_full
        resp_scale = "raw"
        nice_times <- nice_times
        out_t_scale = NULL
        probs = c(0.1, 0.5, 0.9)
        mean_instead = FALSE
        nice_category_names = "Capelin_Random"
        climate_scenario = climate_scenario = paste0("gfdl", "_full")
        out_dir = date_dir

        # ## Cod -- COG is off...
        # date_dir <- "~/GitHub/mixedmod_projections/2022-10-25/Cod_BC/"
        # vast_fit = readRDS(paste0(date_dir, "SpST_mod_fit.rds"))
        # sim_obj = readRDS(file = paste0(date_dir, "SpST", "_random_ProjectionsList.rds"))
        # resp_scale = "raw"
        # nice_times <- as.Date(c(paste0(seq(from = 1985, to = 2100, by = 1), "-03-16"), paste0(seq(from = 1985, to = 2100, by = 1), "-07-16"), paste0(seq(from = 1985, to = 2100, by = 1), "-10-16")))
        # nice_times <- nice_times[order(nice_times)]
        # out_t_scale = NULL
        # nice_category_names = paste0("Cod_", "Base_None")
        # climate_scenario = paste0("EnvOnly_Base_5thpercentile", "_SSP5_85")
        # out_dir = date_dir
    }
    
    time_ind <- seq(from = 1, to = length(nice_times))
    time_labels <- nice_times
    index_regions_ind <- seq(from = 1, to = vast_fit$data_list$n_l)
    index_regions <- vast_fit$settings$strata.limits$STRATA[index_regions_ind]
    categories_ind <- seq(from = 1, to = vast_fit$data_list$n_c)
    grid_ind <- seq(from = 1, to = vast_fit$data_list$n_g)

    for (i in seq_along(sim_obj)) {

        if(FALSE){
            # Checking sim_obj 
            summary(sim_obj[[100]][["Index_ctl"]])
        }

        # Going to want the Index...
        ind_array <- array(unlist(sim_obj[[i]][which(names(sim_obj[[i]]) ==
            "Index_ctl")]),
            dim = c(unlist(vast_fit$data_list[c("n_c")]), n_t = length(nice_times), unlist(vast_fit$data_list[c("n_l")])),
            dimnames = list(
                categories_ind, time_labels,
                index_regions
            )  
        )
        ind_df <- data.frame(aperm(ind_array, c(2, 3, 1)))
        colnames(ind_df) <- gsub(".1", "", colnames(ind_df))
        ind_df$Sim_Scenario <- rep(paste0("Sim_", i), nrow(ind_df))
        ind_df$Time <- nice_times
        # ind_df[, 1] <- ifelse(ind_df[, 1] == 0, NA, ind_df[, 1])
        # ind_df <- na.omit(ind_df)
        ind_df <- ind_df %>%
            pivot_longer(., !c(
                Sim_Scenario,
                Time
            ), names_to = "Region", values_to = "Index")

        # Check
        if(FALSE){
            true <- vast_fit$Report$Index_ctl
            str(true)
            str(ind_df)
        }
        
        # Density
        dens_array <- array(unlist(sim_obj[[i]][which(names(sim_obj[[i]]) ==
            "D_gct")]),
            dim = c(unlist(vast_fit$data_list[c("n_g")]), (vast_fit$data_list[c("n_c")]), n_t = length(nice_times)),
            dimnames = list(grid_ind, categories_ind, time_labels)
        )
        dens_df <- data.frame(aperm(dens_array, c(1, 3, 2)))
        colnames(dens_df) <- nice_times
        dens_df$Lat <- vast_fit$spatial_list$latlon_g[, "Lat"]
        dens_df$Lon <- vast_fit$spatial_list$latlon_g[, "Lon"]
        dens_df <- dens_df %>%
            pivot_longer(.,
                !c(Lat, Lon),
                names_to = "Time", values_to = "D_gct"
            ) %>%
            arrange(Time, Lat, Lon)
        dens_df$Time <- as.Date(dens_df$Time)
        dens_df$Sim_Scenario <- paste0("Sim_", i)

        # Check
        if(FALSE){
            true <- vast_fit$Report$D_gct
            str(true)
            str(dens_df)
        }

        # Center of gravity
        cog_array <- array(unlist(sim_obj[[i]][which(names(sim_obj[[i]]) ==
            "mean_Z_ctm")]),
        dim = c(unlist(vast_fit$data_list[c("n_c")]), n_t = length(nice_times), 2),
        dimnames = list(
            categories_ind, time_labels,
            c("Lon", "Lat")
        )
        )
        cog_true_df <- data.frame(aperm(cog_array, c(2, 3, 1)))
        names(cog_true_df)<- c("Eastings", "Northings")
        cog_true_df$Time<- nice_times
        cog_true_df$Time <- as.Date(cog_true_df$Time)
        cog_true_df$Sim_Scenario <- paste0("Sim_", i)

        cog_df <- cog_from_dens(vast_fit = vast_fit, proj_dens = dens_df, proj_ind = ind_df)

        # Check
        if(FALSE){
            true_lon <- vast_fit$Report$mean_Z_ctm[,,1]
            str(true_lon)
            true_lat<- vast_fit$Report$mean_Z_ctm[,,2]
            str(true_lat)
            str(cog_df)
        }

        # Effective area
        eff_area_array <- array(unlist(sim_obj[[i]][which(names(sim_obj[[i]]) ==
            "effective_area_ctl")]),
            dim = c(unlist(vast_fit$data_list[c("n_c")]), n_t = length(nice_times), unlist(vast_fit$data_list[c("n_l")])),
            dimnames = list(
                categories_ind, time_labels,
                index_regions
            )  
        )
        eff_area_true_df <- data.frame(aperm(eff_area_array, c(2, 3, 1)))
        colnames(eff_area_true_df) <- gsub(".1", "", colnames(eff_area_true_df))
        eff_area_true_df$Sim_Scenario <- rep(paste0("Sim_", i), nrow(eff_area_true_df))
        eff_area_true_df$Time <- nice_times
        # ind_df[, 1] <- ifelse(ind_df[, 1] == 0, NA, ind_df[, 1])
        # ind_df <- na.omit(ind_df)
        eff_area_true_df <- eff_area_true_df %>%
            pivot_longer(., !c(
                Sim_Scenario,
                Time
            ), names_to = "Region", values_to = "Eff_Area")

        
        eff_area_df <- eff_area_from_dens(vast_fit = vast_fit, proj_dens = dens_df, proj_ind = ind_df)

        # Check
        if (FALSE) {
            vast_fit$Report$effective_area_ctl
            str(true_eff_area)
            str(eff_area_df)
        }
        
        if(resp_scale == "log"){
            ind_df$Index <- log(ind_df$Index)
            dens_df$D_gct <- log(dens_df$D_gct)
        }
        
        if (i == 1) {
            res_out_ind <- ind_df
            res_out_dens <- dens_df
            res_out_cog <- cog_df
            res_out_cog_true<- cog_true_df
            res_out_eff_area <- eff_area_df
            res_out_eff_area_true<- eff_area_true_df
        } else {
            res_out_ind <- bind_rows(res_out_ind, ind_df)
            res_out_dens <- bind_rows(res_out_dens, dens_df)
            res_out_cog <- bind_rows(res_out_cog, cog_df)
            res_out_cog_true <- bind_rows(res_out_cog_true, cog_true_df)
            res_out_eff_area <- bind_rows(res_out_eff_area, eff_area_df)
            res_out_eff_area_true<- bind_rows(res_out_eff_area_true, eff_area_true_df)
        }
    }

    # Calculate summaries across all runs
    res_out_ind <- res_out_ind %>%
        group_by(., Time, Region) %>%
        summarise(
            Prob_0.5 = quantile(Index, probs = 0.5, na.rm = TRUE),
            Prob_0.1 = quantile(Index, probs = 0.1, na.rm = TRUE),
            Prob_0.9 = quantile(Index, probs = 0.9, na.rm = TRUE)
        )
    
    res_out_dens <- res_out_dens %>%
        group_by(., Lat, Lon, Time) %>%
        summarise(
            Prob_0.5 = quantile(D_gct, probs = 0.5, na.rm = TRUE),
            Prob_0.1 = quantile(D_gct, probs = 0.1, na.rm = TRUE),
            Prob_0.9 = quantile(D_gct, probs = 0.9, na.rm = TRUE)
        )
    
    res_out_cog <- res_out_cog %>%
        group_by(., Time) %>%
        summarise(
            Lon_Prob_0.5 = quantile(Mean_Lon, probs = 0.5, na.rm = TRUE),
            Lon_Prob_0.1 = quantile(Mean_Lon, probs = 0.1, na.rm = TRUE),
            Lon_Prob_0.9 = quantile(Mean_Lon, probs = 0.9, na.rm = TRUE),
            Lat_Prob_0.5 = quantile(Mean_Lat, probs = 0.5, na.rm = TRUE),
            Lat_Prob_0.1 = quantile(Mean_Lat, probs = 0.1, na.rm = TRUE),
            Lat_Prob_0.9 = quantile(Mean_Lat, probs = 0.9, na.rm = TRUE)
        )
    
    res_out_cog_true<- res_out_cog_true %>%
        group_by(., Time) %>%
        summarise(
            Lon_Prob_0.5 = quantile(Eastings, probs = 0.5, na.rm = TRUE),
            Lon_Prob_0.1 = quantile(Eastings, probs = 0.1, na.rm = TRUE),
            Lon_Prob_0.9 = quantile(Eastings, probs = 0.9, na.rm = TRUE),
            Lat_Prob_0.5 = quantile(Northings, probs = 0.5, na.rm = TRUE),
            Lat_Prob_0.1 = quantile(Northings, probs = 0.1, na.rm = TRUE),
            Lat_Prob_0.9 = quantile(Northings, probs = 0.9, na.rm = TRUE)
        )
    
    res_out_eff_area <- res_out_eff_area %>%
        group_by(., Time, Region) %>%
        summarise(
            Prob_0.5 = quantile(Eff_Area, probs = 0.5, na.rm = TRUE),
            Prob_0.1 = quantile(Eff_Area, probs = 0.1, na.rm = TRUE),
            Prob_0.9 = quantile(Eff_Area, probs = 0.9, na.rm = TRUE)
        )
    
    res_out_eff_area_true <- res_out_eff_area_true %>%
        group_by(., Time, Region) %>%
        summarise(
            Prob_0.5 = quantile(Eff_Area, probs = 0.5, na.rm = TRUE),
            Prob_0.1 = quantile(Eff_Area, probs = 0.1, na.rm = TRUE),
            Prob_0.9 = quantile(Eff_Area, probs = 0.9, na.rm = TRUE)
        )
    
    res_out<- list("Index" = res_out_ind, "Dens" = res_out_dens, "COG" = res_out_cog, "COG_True" = res_out_cog_true, "EffArea" = res_out_eff_area, "EffArea_True" = res_out_eff_area_true)
 
    saveRDS(res_out, file = paste0(out_dir, "/", nice_category_names, 
        "_", climate_scenario, ".rds"))
    return(res_out)
}
    