#' gridEZ: generate a gridded sampling frame
#'
#' @description
#' gridEZ generates a sampling frame consisting of gridded enumeration zones (EZs). EZs are constructed based on a target population count per EZ whilst being restricted by a maximum geographic size.
#' gridEZ uses gridded population counts and two stratification datasets that must all be rasters in the same coordinate reference system (CRS) with the same spatial extents.
#' For household or population surveys, the stratification datasets could be settlement types and administrative units, for example. If one or no stratification is needed, the "strata1" and/or "strata2" rasters (specified in the function parameters) can simply contain one value everywhere. Note that large areas for a single stratum is not advised as it will require a high level of computation and likely lead to very slow running time and/or a memory error.
#'
#' @param population_raster_path Character string. Full file path to a population raster file. Each pixel/cell should contain a population count. NAs allowed. Cells containing NAs are excluded from EZs.
#' @param strata1_raster_path Character string. Full file path to a strata raster file. Each pixel/cell should contain an ID corresponding to a stratum, e.g. administrative unit. NAs allowed. Cells containing NAs are excluded from EZs.
#' @param strata2_raster_path Character string. Full file path to a strata raster file. Each pixel/cell should contain an ID corresponding to a stratum. Strata 2 will usually be different to strata 1 depending on the stratification requirements of a survey, e.g. settlement type classification, climate zone, categories of distance to a geographical feature. If only one form of stratification is required, simply repeat the file path for Strata 1 here. NAs allowed. Cells containing NAs are excluded from EZs.
#' @param output_path Character string. File path for the location you would like output files to be saved.
#' @param run_ID Character string. Default is "_run1". Text to be added to the end of output file names.
#' @param predefined_EZ_size Logical. Default is TRUE. If TRUE, gridEZ will use pre-set parameters for target_pop_per_EZ & max_cells_per_EZ to construct EZs. Set to FALSE if you want to specify your own target_pop_per_EZ & max_cells_per_EZ.
#' @param EZ_target_size Character sting. Considered only if predefined_EZ_size = TRUE. Default is "medium". Can be either "small", "medium" (or "med") or "large". For "small" target_pop_per_EZ = 75 and max_cells_per_EZ = 100; for "medium" target_pop_per_EZ = 500 and max_cells_per_EZ = 900; for "large" target_pop_per_EZ = 1200 and max_cells_per_EZ = 2500
#' @param target_pop_per_EZ Numeric, single number. Considered only if predefined_EZ_size = FALSE. Default is 500. This should be specified if you would like to define your own target population per EZ.
#' @param max_cells_per_EZ Numeric, single number. Considered only if predefined_EZ_size = FALSE. Default is 900 raster cells. This should be specified if you would like to define your own maximum cells per EZ. Note, this is not a hard maximum to allow flexibility in EZs reaching the target_pop_per_EZ.
#' @param ncores Numeric, single number. Default is 1. The number of computer cores you would like for processing. Recommend using your total number of cores minus 1. Use parallel::detectCores() to find the total number of cores available on your computer.
#' @returns Raster and csv saved to file. Raster contain the EZ ID for each grid cell. The csv contains columns of EZ IDs ("EZ_IDs"), population per EZ ("total_pop") and number of grid cells per EZ ("N").
#' @export
gridEZ <- function(population_raster_path, strata1_raster_path, strata2_raster_path,
                   predefined_EZ_size = TRUE, EZ_target_size = "medium",
                   target_pop_per_EZ = 500, max_cells_per_EZ = 900,
                   output_path, run_ID = "run1",
                   ncores=1){

  ############################################
  # initial raster and specifications checks #
  ############################################

  population_raster <- normalizePath(paste(population_raster_path))
  strata1_raster <- normalizePath(paste(strata1_raster_path))
  strata2_raster <- normalizePath(paste(strata2_raster_path))

  if(file.exists(paste(output_path,"/EZ_raster_master", run_ID, ".tif",sep=""))){
    stop(paste("Desired output filename already exists; delete file or change run_ID so that output file is given a different name. File location and name is:", output_path,"/EZ_raster_master", run_ID, ".tif",sep=""))
  }
  if(file.exists(paste(output_path,"/EZ_pop_raster_master", run_ID, ".tif",sep=""))){
    stop(paste("Desired output filename already exists; delete file or change run_ID so that output file is given a different name. File location and name is:", output_path,"/EZ_pop_raster_master", run_ID, ".tif",sep=""))
  }
  if(terra::crs(terra::rast(population_raster)) != terra::crs(terra::rast(strata2_raster)) | terra::crs(terra::rast(population_raster)) != terra::crs(terra::rast(strata1_raster))) {
    stop(paste("Projection system of rasters differs. You can use terra::crs(terra::rast()) to confirm raster projections."))
  }
  if(terra::ext(terra::rast(population_raster)) != terra::ext(terra::rast(strata2_raster)) | terra::ext(terra::rast(population_raster)) != terra::ext(terra::rast(strata1_raster))) {
    stop(paste("Extent of rasters differs. You can use terra::ext(terra::rast()) to confirm raster extents."))
  }
  if(predefined_EZ_size == TRUE){
    if(EZ_target_size == "small"){
      target_pop_per_EZ <- 75
      max_cells_per_EZ <- 100
    }
    if(EZ_target_size == "medium"){
      target_pop_per_EZ <- 500
      max_cells_per_EZ <- 900
    }
    if(EZ_target_size == "med"){
      target_pop_per_EZ <- 500
      max_cells_per_EZ <- 900
    }
    if(EZ_target_size == "large"){
      target_pop_per_EZ <- 1200
      max_cells_per_EZ <- 2500
    }
    print(paste("Creating EZs based on predefined EZ size criteria:"))
    print(paste("target population per EZ = ", target_pop_per_EZ, " & max number of cells per EZ = ", max_cells_per_EZ))
  }else{
    if(is.numeric(target_pop_per_EZ) == FALSE | !length(max_cells_per_EZ) == 1){
      stop("target_pop_per_EZ needs to be a single number")
    }
    print("EZs will be created based on user's specified population per EZ")
    if(is.numeric(max_cells_per_EZ) == FALSE | !length(max_cells_per_EZ) == 1){
      stop("max_cells_per_EZ needs to be a single number")
    }
    print(paste("Creating EZs based on user defined EZ size criteria:"))
    print(paste("target population per EZ = ", target_pop_per_EZ, " & max number of cells per EZ = ", max_cells_per_EZ))
    }

  ###########################
  # create temporary folder #
  ###########################

  print(paste("Creating temporary folder", " /temp_folder_", run_ID, sep=""))
  dir.create(paste(output_path,"/temp_folder_", run_ID, sep=""))

  ####################################################################
  # split study are into to clumps defined by uncrossable boundaries #
  ####################################################################

  r_crs <- terra::crs(terra::rast(strata1_raster))
  r_xmn <- terra::ext(terra::rast(strata1_raster))[1]
  r_xmx <- terra::ext(terra::rast(strata1_raster))[2]
  r_ymn <- terra::ext(terra::rast(strata1_raster))[3]
  r_ymx <- terra::ext(terra::rast(strata1_raster))[4]
  strata1_ids <- unique(terra::values(terra::rast(strata1_raster)))
  strata1_ids <- stats::na.omit(strata1_ids)
  strata2_ids <- unique(terra::values(terra::rast(strata2_raster)))
  strata2_ids <- stats::na.omit(strata2_ids)
  ids <- cbind(rep(strata1_ids,each=length(strata2_ids)),rep(strata2_ids,length(strata1_ids)))
  clump_extents <- lapply(split(ids, seq(nrow(ids))),function(x){
    i <- x[1]; j <- x[2]
    target_cells <- which(terra::values(terra::rast(strata1_raster)) == i &
                            terra::values(terra::crop(terra::rast(strata2_raster),terra::rast(strata1_raster),extend=TRUE)) == j)
    if(length(target_cells)>0){ # possible to be zero as not all strata2 types in each stratum
      value_extent <- terra::ext(terra::rast(strata1_raster), cells = target_cells)
      list(x,c(as.numeric(value_extent[1]),
                 as.numeric(value_extent[2]),
                 as.numeric(value_extent[3]),
                 as.numeric(value_extent[4])))
    }
  })
  clump_extents <- Filter(Negate(is.null), clump_extents)

  rm(strata1_ids,strata2_ids,ids)
  gc()

  ###################################
  # specify EZ generation function  #
  ###################################

  EZ_gen_fn <- function(clump_x, target_pop_per_EZ, max_cells_per_EZ, output_path, run_ID, population_raster, strata2_raster, strata1_raster){

    # save a temporary file for clump currently being processed
    fake <- list(1)
    tempfileID <- paste("strata1_",(clump_x[[1]][1]), "_strata2_", (clump_x[[1]][2]), sep="")
    save(fake, file = paste(output_path,"/temp_folder_", run_ID, "/temp_", tempfileID, ".RData",sep=""))

    ###################################
    # prep & check clump information  #
    ###################################

    min_pop_per_EZ <- target_pop_per_EZ * 0.66666  # 66% of target pop per EZ
    max_pop_per_EZ <- target_pop_per_EZ * 1.33333  # 33% of target pop per EZ
    subpop_rast <- terra::crop(terra::rast(population_raster),
                               terra::ext(c(clump_x[[2]][1], clump_x[[2]][2],
                                            clump_x[[2]][3], clump_x[[2]][4])))
    subpop_rast <- terra::mask(subpop_rast,
                               terra::crop(terra::rast(strata1_raster),
                                           terra::ext(c(clump_x[[2]][1], clump_x[[2]][2],
                                                        clump_x[[2]][3], clump_x[[2]][4]))),
                               inverse=TRUE,maskvalues=clump_x[[1]][1])
    subpop_rast <- terra::mask(subpop_rast,
                               terra::crop(terra::rast(strata2_raster),
                                           terra::ext(c(clump_x[[2]][1], clump_x[[2]][2],
                                                        clump_x[[2]][3], clump_x[[2]][4]))),
                               inverse=TRUE,maskvalues=clump_x[[1]][2])
    clump_pop_size <- sum(terra::values(subpop_rast), na.rm = TRUE)

    ###################################################################################
    # if clump contains small pop and small number of cells, skip to contiguity check #
    ###################################################################################

    if(((clump_pop_size <= max_pop_per_EZ) & (length(subpop_rast[!is.na(subpop_rast)]) <= max_cells_per_EZ * 1.5)) | length(subpop_rast[!is.na(subpop_rast)])==1){
      EZ_mat <- subpop_mat <- as.matrix(subpop_rast,wide=TRUE)
      EZ_ID_vec <- EZ_mat[!is.na(EZ_mat)] <- 1

    }else{

    #####################################################################################
    # if clump contains pop and cells above min requirements, apply EZ generation steps #
    #####################################################################################

      #################################################
      # create initial blocked mesh of potential EZs  #
      #################################################

      subpop_mat <- as.matrix(subpop_rast,wide=TRUE)
      max_cell_pop_size <- max(subpop_mat[!(is.na(subpop_mat))])
      med_cell_pop_size <- stats::median(subpop_mat[!(is.na(subpop_mat))])
      min_stratum_row <- min_stratum_col <- 1
      max_stratum_row <- nrow(subpop_mat)
      max_stratum_col <- ncol(subpop_mat)
      if((max_cell_pop_size * max_cells_per_EZ) < min_pop_per_EZ){
        block_horiz <- block_vert <- round(sqrt(max_cells_per_EZ))
      }else{
        val1 <- round(sqrt(target_pop_per_EZ/med_cell_pop_size)) # var1 = est of cells horiz (or vert)
        if(val1 <= 0){val1 <- 1}
        if(val1 > round(sqrt(max_cells_per_EZ))){
          block_horiz <- block_vert <- round(sqrt(max_cells_per_EZ))
        }else{
          if(val1 %in% 1:2){
            block_vert <- block_horiz <- val1
          }else{
            if(val1 < 5){val2 <- val1 - 1}
            if(val1 >= 5){val2 <- val1 - 2}
            row_range <- max_stratum_row - min_stratum_row
            col_range <- max_stratum_col - min_stratum_col
            if(row_range >= col_range){
              block_vert <- val1
              block_horiz <- val2
            }else{
              block_vert <- val2
              block_horiz <- val1
            }
          }
        }
      }
      # ID location of max pop cell and use this as the centre of a block. Build block frame from there
      ID_max_pop_cell <- (which(subpop_mat == max_cell_pop_size, arr.ind = TRUE))[1,]
      y_start <- min((ID_max_pop_cell[1] + floor(block_vert/2)), max_stratum_row) #/2 so that max cell is middle of block
      x_start <- min((ID_max_pop_cell[2] + floor(block_horiz/2)), max_stratum_col)
      seq_vert <- y_start; seq_horiz <- x_start
      # build grid for blocks. If block edge is close to the edge of the raster section, use alternative initial block mesh location
      if((y_start - block_vert) > min_stratum_row){seq_vert <- rev(seq(y_start, min_stratum_row, -block_vert))}
      if((y_start + block_vert) < max_stratum_row){seq_vert <- c(seq_vert, seq((y_start + block_vert), max_stratum_row, block_vert))}
      if(!seq_vert[1] == min_stratum_row){seq_vert <- c(min_stratum_row, seq_vert)}
      if(!seq_vert[length(seq_vert)] == max_stratum_row){seq_vert <- c(seq_vert, max_stratum_row)}
      if((x_start - block_horiz) > min_stratum_col){seq_horiz <- rev(seq(x_start, min_stratum_col, -block_horiz))}
      if((x_start + block_horiz) < max_stratum_col){seq_horiz <- c(seq_horiz, seq((x_start + block_horiz), max_stratum_col, block_horiz))}
      if(!seq_horiz[1] == min_stratum_col){seq_horiz <- c(min_stratum_col, seq_horiz)}
      if(!seq_horiz[length(seq_horiz)] == max_stratum_col){seq_horiz <- c(seq_horiz, max_stratum_col)}
      if(block_vert == 1 | length(seq_vert) == 1){block_ref_vert <- cbind(seq_vert,seq_vert)}else{
        block_ref_vert <- cbind(seq_vert[1:(length(seq_vert)-1)], c((seq_vert[2:(length(seq_vert)-1)]-1), max_stratum_row))
      }
      if(block_horiz == 1 | length(seq_horiz) == 1){block_ref_horiz <- cbind(seq_horiz,seq_horiz)}else{
        block_ref_horiz <- cbind(seq_horiz[1:(length(seq_horiz)-1)], c((seq_horiz[2:(length(seq_horiz)-1)]-1), max_stratum_col))
      }
      block_ref_df <- cbind(block_ref_vert[rep(seq_len(nrow(block_ref_vert)), nrow(block_ref_horiz)),], block_ref_horiz[rep(seq_len(nrow(block_ref_horiz)), each = nrow(block_ref_vert)),])
      # assign temporary EZ IDs to all blocks within extent of clump
      EZ_mat <- matrix(NA, nrow = max_stratum_row, ncol = max_stratum_col)
      for(block_ID in 1:nrow(block_ref_df)){    # loop instead of function or apply because want vals in one matrix to be updated
        EZ_mat[unique(block_ref_df[block_ID,1]:block_ref_df[block_ID,2]), unique(block_ref_df[block_ID,3]:block_ref_df[block_ID,4])] <- block_ID
      }
      EZ_mat[which(is.na(subpop_mat))] <- NA
      dt <- data.table::data.table(EZ_ID = as.vector(EZ_mat), pop = as.vector(subpop_mat))
      dt <- dt[!is.na(dt[["EZ_ID"]]), ]
      dt <- data.table::as.data.table(stats::aggregate(pop ~ EZ_ID, data = dt, FUN = function(x) sum(x, na.rm = TRUE)))

      rm(seq_horiz,seq_vert,x_start,y_start,val1,block_horiz,block_vert,block_ID,block_ref_df,block_ref_horiz,block_ref_vert,ID_max_pop_cell)
      gc()

      ###############################################################################
      # add blocks (EZs) with pop below min pop per EZ to neighbouring blocks (EZs) #
      ###############################################################################

      small_pop_EZs <- dt[dt[["pop"]] < min_pop_per_EZ, ]
      if(nrow(small_pop_EZs)>0){
        small_pop_EZs <- small_pop_EZs[["EZ_ID"]]
        # of these small pop EZs, exclude those that are close to or above max geog size
        cells_count_test <- table(EZ_mat[EZ_mat %in% small_pop_EZs]) > (max_cells_per_EZ * 0.75)
        if("TRUE" %in% cells_count_test){
          max_geog_in_small_pop_EZs <- as.numeric(names(cells_count_test[cells_count_test == "TRUE"]))
          small_pop_EZs <- small_pop_EZs[!small_pop_EZs %in% max_geog_in_small_pop_EZs]
        }
        small_pop_EZs <- as.numeric(names(sort(table(EZ_mat[EZ_mat %in% small_pop_EZs])))) # re-order small_pop_EZs by geog size (smallest first)
      }else{
        small_pop_EZs <- NULL
      }
      # start adding small pop EZs to suitable neighbouring EZs
      if(length(small_pop_EZs)>0){
        rm(dt)
        repeat{
          EZ_ID <- small_pop_EZs[1]
          # check out neighbour EZs
          EZ_cell_IDs <- which(EZ_mat == EZ_ID, arr.ind = TRUE)
          nb_EZs <- cell_nbs <- NULL
          for(df_row in 1:nrow(EZ_cell_IDs)){
            loc_vert <- EZ_cell_IDs[df_row,1]
            loc_horiz <- EZ_cell_IDs[df_row,2]
            adj_cells <- cbind(c((loc_vert + 1), (loc_vert - 1), loc_vert, loc_vert), c(loc_horiz, loc_horiz, (loc_horiz + 1), (loc_horiz - 1)))
            if(TRUE %in% (!adj_cells[,1] %in% 1:max_stratum_row)){adj_cells <- adj_cells[-which((!adj_cells[,1] %in% 1:max_stratum_row)),]}
            if(TRUE %in% (!adj_cells[,2] %in% 1:max_stratum_col)){adj_cells <- adj_cells[-which((!adj_cells[,2] %in% 1:max_stratum_col)),]}
            cell_nbs <- EZ_mat[adj_cells]
            nb_EZs <- c(nb_EZs, cell_nbs)
          }
          nb_EZs <- unique(nb_EZs)
          nb_EZs <- nb_EZs[!is.na(nb_EZs) & !nb_EZs == EZ_ID]
          if(length(nb_EZs)==0){ # if no neighbours, omit EZ ID from small_pop_EZs list and move to next
            small_pop_EZs <- small_pop_EZs[-which(small_pop_EZs == EZ_ID)]
          }else{
            # measure compactness for each of the possible combinations of blocks/EZs
            compact_vals <- sapply(nb_EZs,function(new_EZ){
            test_EZ_locs <- which(EZ_mat == EZ_ID | EZ_mat == new_EZ ,arr.ind = TRUE)
            test_mat <- matrix(NA,nrow=dim(EZ_mat)[1],ncol=dim(EZ_mat)[2])
            test_mat[test_EZ_locs] <- 1
            test_poly <- terra::as.polygons(terra::rast(test_mat,ext = terra::ext(subpop_rast),crs = terra::crs(subpop_rast)), dissolve = TRUE)
            poly_area <- terra::expanse(test_poly, unit = "m")
            edge_peri <- terra::perim(test_poly)
            (4 * pi * poly_area / (edge_peri^2))
          })
            EZ_ID_to_add <- (nb_EZs[compact_vals == max(compact_vals)])[1]
            EZ_mat[EZ_mat == EZ_ID_to_add] <- EZ_ID
            if(EZ_ID_to_add %in% small_pop_EZs){small_pop_EZs <- small_pop_EZs[-which(small_pop_EZs == EZ_ID_to_add)]}
            if(sum(subpop_mat[which(EZ_mat == EZ_ID)]) > min_pop_per_EZ){
              small_pop_EZs <- small_pop_EZs[-which(small_pop_EZs == EZ_ID)]
            }else{
            if(length(which(EZ_mat[] == EZ_ID)) > (max_cells_per_EZ * 0.8)){
              small_pop_EZs <- small_pop_EZs[-which(small_pop_EZs == EZ_ID)]
            }
          }
          }
          if(!length(small_pop_EZs)>0) break
          small_pop_EZs <- as.numeric(names(sort(table(EZ_mat[EZ_mat %in% small_pop_EZs]))))
        }
        dt <- data.table::data.table(EZ_ID = as.vector(EZ_mat), pop = as.vector(subpop_mat))
        dt <- dt[!is.na(dt[["EZ_ID"]]), ]
        dt <- data.table::as.data.table(stats::aggregate(pop ~ EZ_ID, data = dt, FUN = function(x) sum(x, na.rm = TRUE)))

      }

      ###############################################################
      # split blocks (EZs) with pop greater than 2 x max pop per EZ #
      ###############################################################

      vbig_pop_EZs <- dt[dt[["pop"]] >= (2 * target_pop_per_EZ),]
      if(nrow(vbig_pop_EZs)>0){
        vbig_pop_EZs <- vbig_pop_EZs[["EZ_ID"]]
        # of these very big pop EZs, exclude any single cell EZs, as these can't be split
        cells_count_test <- table(EZ_mat[EZ_mat %in% vbig_pop_EZs]) >= 2
        if("TRUE" %in% cells_count_test){
          vbig_pop_EZs <- as.numeric(names(cells_count_test[cells_count_test == "TRUE"]))
        }
      }else{
        vbig_pop_EZs <- NULL
      }
      # start splitting very big pop EZs in two, along axis with most cells
      if(length(vbig_pop_EZs)>0){
        EZ_IDs_new <- (max(dt[["EZ_ID"]]) + 1):(max(dt[["EZ_ID"]]) + 2)
        rm(dt)
        repeat{
          EZ_ID <- vbig_pop_EZs[1]
          EZ_cell_IDs <- which(EZ_mat == EZ_ID, arr.ind = TRUE)
          tot_pop_size <- sum(subpop_mat[EZ_cell_IDs])
          # create list, an element for each potential split, element contains matrix of cell row & col ids one EZ of the potential split, do this because can't assume EZs are only four edged shapes
          row_IDs <- sort(unique(EZ_cell_IDs[,1]))
          col_IDs <- sort(unique(EZ_cell_IDs[,2]))
          if(length(row_IDs) >= length(col_IDs)){
            side_IDs <- row_IDs
            min_side_IDs <- min(row_IDs)
            row_or_col <- 1
          }else{
            side_IDs <- col_IDs
            min_side_IDs <- min(col_IDs)
            row_or_col <- 2
          }
          cell_IDs_ls <- lapply(side_IDs[1:(length(side_IDs) - 1)], function(x){matrix(EZ_cell_IDs[EZ_cell_IDs[,row_or_col] %in% (min_side_IDs:x),], ncol=2)})
          pop_sizes <- sapply(cell_IDs_ls, function(x){sum(subpop_mat[x])})
          pop_sizes_second <- tot_pop_size - pop_sizes
          keepers1 <- which(pop_sizes >= min_pop_per_EZ)
          keepers2 <- which(pop_sizes_second >= min_pop_per_EZ)
          test_splits <- keepers1[which(keepers1 %in% keepers2)] # so only keep split options that lead to both new EZs having pop larger than min_pop_per_EZ value
          if(!length(test_splits) > 0){
            vbig_pop_EZs <- vbig_pop_EZs[-(which(vbig_pop_EZs == EZ_ID))]
            if(!length(vbig_pop_EZs)>0){
              break
            }else{
              next    # don't do anything if splitting will result in new EZs with pop too small
            }
          }
          if(length(test_splits) > 10){
            test_splits <- test_splits[seq(1,length(test_splits),ceiling(length(test_splits)/10))]
          }
          compact_vals <-  lapply(test_splits,function(test_EZ_split){
            compacts <- NULL
            for(x in 1:2){
              if(x==1){test_EZ_locs <- cell_IDs_ls[[test_EZ_split]]}
              if(x==2){
                d1 <- duplicated(rbind(cell_IDs_ls[[test_EZ_split]],EZ_cell_IDs))
                d1 <- d1[(nrow(cell_IDs_ls[[test_EZ_split]])+1):length(d1)]
                test_EZ_locs <- EZ_cell_IDs[which(d1 == FALSE),]
              }
              if(!class(test_EZ_locs)[1] == "matrix"){test_EZ_locs <- matrix(test_EZ_locs, nrow=1)}
              colnames(test_EZ_locs) <- c("row","col")
              test_mat <- matrix(NA,nrow=dim(EZ_mat)[1],ncol=dim(EZ_mat)[2])
              test_mat[test_EZ_locs] <- 1
              test_poly <- terra::as.polygons(terra::rast(test_mat,ext = terra::ext(subpop_rast),crs = terra::crs(subpop_rast)), dissolve = TRUE)
              poly_area <- terra::expanse(test_poly, unit = "m")
              edge_peri <- terra::perim(test_poly)
              compacts <- c(compacts, (4 * pi * poly_area / (edge_peri^2)))
            }
            compacts
          })
          compact_vals <- do.call(rbind, compact_vals)
          compact_sum <- round((compact_vals[,1] + compact_vals[,2]),3)
          compact_threshold <- max(compact_sum) # take the most compact EZ option as all splits being tested will result in EZs with pop larger than min pop per EZ. EZs with too big pop will be split in later steps
          selected_split_ID <- test_splits[which(compact_sum == compact_threshold)][1]
          EZ_mat[cell_IDs_ls[[selected_split_ID]]] <- EZ_IDs_new[1]
          EZ_ID_pop <- sum(subpop_mat[which(EZ_mat == EZ_ID)])
          if(EZ_ID_pop < (2 * target_pop_per_EZ)){vbig_pop_EZs <- vbig_pop_EZs[-(which(vbig_pop_EZs == EZ_ID))]}
          EZ_ID_new_pop <- sum(subpop_mat[which(EZ_mat == EZ_IDs_new[1])])
          if(EZ_ID_new_pop >= (2 * target_pop_per_EZ)){vbig_pop_EZs <- c(vbig_pop_EZs, EZ_IDs_new[1])}
          EZ_IDs_new <- EZ_IDs_new[-1]; EZ_IDs_new <- c(EZ_IDs_new, (max(EZ_IDs_new) + 1))
          if(!length(vbig_pop_EZs)>0) break
        }
        dt <- data.table::data.table(EZ_ID = as.vector(EZ_mat), pop = as.vector(subpop_mat))
        dt <- dt[!is.na(dt[["EZ_ID"]]), ]
        dt <- data.table::as.data.table(stats::aggregate(pop ~ EZ_ID, data = dt, FUN = function(x) sum(x, na.rm = TRUE)))
      }

      #################################################################################
      # add blocks (EZs) with pop greater than max pop per EZ, then split into 2 or 3 #
      #################################################################################

      big_pop_EZs <- dt[dt[["pop"]] > max_pop_per_EZ,]
      if(nrow(big_pop_EZs)>0){
        big_pop_EZs <- big_pop_EZs[["EZ_ID"]]
        # of these big pop EZs, exclude any single cell EZs, as these can't get any smaller
        cells_count_test <- table(EZ_mat[EZ_mat %in% big_pop_EZs]) >= 2
        if("TRUE" %in% cells_count_test){
          big_pop_EZs <- as.numeric(names(cells_count_test[cells_count_test == "TRUE"]))
        }
      }else{
        big_pop_EZs <- NULL
      }
      # start processing big pop EZs
      if(length(big_pop_EZs)>0){
        EZ_IDs_new <- (max(dt[["EZ_ID"]]) + 1):(max(dt[["EZ_ID"]]) + 3)
        rm(dt)
        repeat{
          EZ_ID <- big_pop_EZs[1]
          EZ_cell_IDs <- which(EZ_mat == EZ_ID, arr.ind = TRUE)
          #neighbour EZs
          nb_EZs <- cell_nbs <- NULL
          for(df_row in 1:nrow(EZ_cell_IDs)){
            loc_vert <- EZ_cell_IDs[df_row,1]
            loc_horiz <- EZ_cell_IDs[df_row,2]
            adj_cells <- cbind(c((loc_vert + 1), (loc_vert - 1), loc_vert, loc_vert), c(loc_horiz, loc_horiz, (loc_horiz + 1), (loc_horiz - 1)))
            if(TRUE %in% (!adj_cells[,1] %in% 1:max_stratum_row)){adj_cells <- adj_cells[-which((!adj_cells[,1] %in% 1:max_stratum_row)),]}
            if(TRUE %in% (!adj_cells[,2] %in% 1:max_stratum_col)){adj_cells <- adj_cells[-which((!adj_cells[,2] %in% 1:max_stratum_col)),]}
            cell_nbs <- EZ_mat[adj_cells]
            nb_EZs <- c(nb_EZs, cell_nbs)
          }
          nb_EZs <- unique(nb_EZs)
          nb_EZs <- nb_EZs[!is.na(nb_EZs) & !nb_EZs == EZ_ID]
          if(length(nb_EZs) > 0){
            EZ_pop_size <- sum(subpop_mat[EZ_cell_IDs])
            nb_EZ_pop_sizes <- NULL
            for(nb_EZ in nb_EZs){
              nb_EZ_cell_IDs <- which(EZ_mat == nb_EZ, arr.ind = TRUE)
              nb_EZ_pop_size <- sum(subpop_mat[nb_EZ_cell_IDs])
              nb_EZ_pop_sizes <- c(nb_EZ_pop_sizes, nb_EZ_pop_size)
            }
            new_pop_sizes <- (nb_EZ_pop_sizes + EZ_pop_size)/2
            if(max(new_pop_sizes) < min_pop_per_EZ){
              big_pop_EZs <- big_pop_EZs[-which(big_pop_EZs == EZ_ID)]
              if(!length(big_pop_EZs)>0){
                break
              }else{
                next   # don't do anything if splitting will result in new EZs with too small pop
              }
            }
            # combine with appropriate neighbouring EZ
            nb_EZ_ID <- which(new_pop_sizes >= min_pop_per_EZ)
            if(length(nb_EZ_ID) > 1){
              nb_EZs <- nb_EZs[nb_EZ_ID]
              compact_vals <- sapply(nb_EZs,function(new_EZ){
                test_EZ_locs <- which(EZ_mat == EZ_ID | EZ_mat == new_EZ ,arr.ind = TRUE)
                test_mat <- matrix(NA,nrow=dim(EZ_mat)[1],ncol=dim(EZ_mat)[2])
                test_mat[test_EZ_locs] <- 1
                test_poly <- terra::as.polygons(terra::rast(test_mat,ext = terra::ext(subpop_rast),crs = terra::crs(subpop_rast)), dissolve = TRUE)
                poly_area <- terra::expanse(test_poly, unit = "m")
                edge_peri <- terra::perim(test_poly)
                (4 * pi * poly_area / (edge_peri^2))
              })
              nb_EZ <- (nb_EZs[compact_vals == max(compact_vals)])[1]
            }else{
              nb_EZ <- nb_EZs[nb_EZ_ID]
            }
            # create single big EZ to then split into three or two
            EZ_mat[EZ_mat == nb_EZ] <- EZ_ID
            if(nb_EZ %in% big_pop_EZs){big_pop_EZs <- big_pop_EZs[-which(big_pop_EZs == nb_EZ)]}
          }
          EZ_cell_IDs <- which(EZ_mat == EZ_ID, arr.ind = TRUE)
          tot_pop_size <- sum(subpop_mat[EZ_cell_IDs])
          if(((tot_pop_size/3) > min_pop_per_EZ) & (nrow(EZ_cell_IDs) >= 3)){
            split_types <- c(1,2)     # splitting twice
          }else{
            split_types <- 2          # splitting once (this split is same as second split when splitting twice)
          }
          #split shape into most compact shapes
          for(split_type in split_types){
            if(split_type == 2 & length(split_types) == 2){     # this updates cell IDs if doing 2nd split
              EZ_cell_IDs <- which(EZ_mat == EZ_ID, arr.ind = TRUE)
              tot_pop_size <- sum(subpop_mat[EZ_cell_IDs])
              if(nrow(EZ_cell_IDs) == 1) break
            }
            row_IDs <- sort(unique(EZ_cell_IDs[,1]))
            col_IDs <- sort(unique(EZ_cell_IDs[,2]))
            if(length(row_IDs) >= length(col_IDs)){
              side_IDs <- row_IDs
              min_side_IDs <- min(row_IDs)
              row_or_col <- 1
            }else{
              side_IDs <- col_IDs
              min_side_IDs <- min(col_IDs)
              row_or_col <- 2
            }
            cell_IDs_ls <- lapply(side_IDs[1:(length(side_IDs) - 1)], function(x){matrix(EZ_cell_IDs[EZ_cell_IDs[,row_or_col] %in% (min_side_IDs:x),], ncol=2)})
            pop_sizes <- sapply(cell_IDs_ls, function(x){sum(subpop_mat[x])})
            if(split_type == 1){
              diff_to_ideal_pop <- abs((tot_pop_size/3) - pop_sizes)
            }
            if(split_type == 2){
              diff_to_ideal_pop <- abs(tot_pop_size - (2*pop_sizes))  # A = B + C, this gives |C - B| i.e. diff in the two potential EZs' pop sizes
            }
            test_splits <- order(diff_to_ideal_pop)
            if(length(test_splits) >= 10){test_splits <- test_splits[1:10]}
            compact_vals <-  lapply(test_splits,function(test_EZ_split){
              compacts <- NULL
              for(x in 1:2){
                if(x==1){test_EZ_locs <- cell_IDs_ls[[test_EZ_split]]}
                if(x==2){
                  d1 <- duplicated(rbind(cell_IDs_ls[[test_EZ_split]],EZ_cell_IDs))
                  d1 <- d1[(nrow(cell_IDs_ls[[test_EZ_split]])+1):length(d1)]
                  test_EZ_locs <- EZ_cell_IDs[which(d1 == FALSE),]
                }
                if(!class(test_EZ_locs)[1] == "matrix"){test_EZ_locs <- matrix(test_EZ_locs, nrow=1)}
                colnames(test_EZ_locs) <- c("row","col")
                test_mat <- matrix(NA,nrow=dim(EZ_mat)[1],ncol=dim(EZ_mat)[2])
                test_mat[test_EZ_locs] <- 1
                test_poly <- terra::as.polygons(terra::rast(test_mat,ext = terra::ext(subpop_rast),crs = terra::crs(subpop_rast)), dissolve = TRUE)
                poly_area <- terra::expanse(test_poly, unit = "m")
                edge_peri <- terra::perim(test_poly)
                compacts <- c(compacts, (4 * pi * poly_area / (edge_peri^2)))
              }
              compacts
            })
            compact_vals <- do.call(rbind, compact_vals)
            compact_sum <- round((compact_vals[,1] + compact_vals[,2]),3)
            compact_threshold <- max(compact_sum) * 0.75 # this gives a threshold value of relatively desirable compactness given the options
            selected_split_ID <- test_splits[which(compact_sum >= compact_threshold)[1]] # as this list is ordered in terms of preferred pop sizes, we want the first option that is within the "relative" desirable compactness range
            EZ_mat[cell_IDs_ls[[selected_split_ID]]] <- EZ_IDs_new[1]
            EZ_IDs_new <- EZ_IDs_new[-1]; EZ_IDs_new <- c(EZ_IDs_new, (max(EZ_IDs_new) + 1))
          }
          big_pop_EZs <- big_pop_EZs[-which(big_pop_EZs == EZ_ID)]
          if(!length(big_pop_EZs)>0) break
        }
      }else{
        rm(dt)
      }

      #########################################################################
      # split blocks (EZs) with geog. size greater than max geog. size per EZ #
      #########################################################################

      dt <- data.table::data.table(EZ_ID = as.vector(EZ_mat))
      dt <- dt[!is.na(dt[["EZ_ID"]]), ]
      dt <- table(dt)
      dt <- data.table::data.table(EZ_ID=as.numeric(names(dt)),N=as.vector(dt))
      big_geog_EZs <- dt[dt[["N"]] > (1.5 * max_cells_per_EZ), ]
      if(nrow(big_geog_EZs)>0){
        big_geog_EZs <- big_geog_EZs[["EZ_ID"]]
      }else{
        big_geog_EZs <- NULL
      }
      EZ_ID_vec <- dt[["EZ_ID"]]
      rm(dt)
      if(length(big_geog_EZs)>0){
        EZ_IDs_new <- (max(EZ_ID_vec) + 1):(max(EZ_ID_vec) + length(big_geog_EZs))
        for(EZ_ID in big_geog_EZs){
          EZ_cell_IDs <- which(EZ_mat == EZ_ID, arr.ind = TRUE)
          row_IDs <- sort(unique(EZ_cell_IDs[,1]))
          col_IDs <- sort(unique(EZ_cell_IDs[,2]))
          if(length(row_IDs) >= length(col_IDs)){
            side_IDs <- row_IDs
            min_side_IDs <- min(row_IDs)
            row_or_col <- 1
          }else{
            side_IDs <- col_IDs
            min_side_IDs <- min(col_IDs)
            row_or_col <- 2
          }
          cell_IDs_ls <- lapply(side_IDs[1:(length(side_IDs) - 1)], function(x){matrix(EZ_cell_IDs[EZ_cell_IDs[,row_or_col] %in% (min_side_IDs:x),], ncol=2)})
          cell_counts <- sapply(cell_IDs_ls, function(x){nrow(x)})
          target_geog <- nrow(EZ_cell_IDs)/2
          selected_split_ID <- which(abs(target_geog - cell_counts) == min(abs(target_geog - cell_counts)))[1]
          EZ_mat[cell_IDs_ls[[selected_split_ID]]] <- EZ_IDs_new[1]
          EZ_IDs_new <- EZ_IDs_new[-1]
        }
        EZ_ID_vec <- unique(EZ_mat[which(!is.na(EZ_mat))])
      }

    }

    ################################
    # check all EZs are contiguous #
    ################################

    EZ_ID_max <- max(EZ_ID_vec)
    for(EZ_ID in EZ_ID_vec){
        cell_IDs <- which(EZ_mat[] == EZ_ID, arr.ind = TRUE)
        min_EZ_row <- min(cell_IDs[,1])
        max_EZ_row <- max(cell_IDs[,1])
        min_EZ_col <- min(cell_IDs[,2])
        max_EZ_col <- max(cell_IDs[,2])
        subEZ_IDs_mat <- EZ_mat[min_EZ_row:max_EZ_row, min_EZ_col:max_EZ_col]
        subEZ_IDs_mat[which(!subEZ_IDs_mat[] == EZ_ID)] <- NA
        if(class(subEZ_IDs_mat)[1] == "matrix"){
          EZ_clump_raster <- terra::patches(terra::rast(subEZ_IDs_mat), directions = 4, allowGaps = FALSE)
          clump_ls <- sort(unique(terra::values(EZ_clump_raster)))
          if(length(clump_ls) > 1){
            new_EZ_IDs_vec <- c(EZ_ID, seq((EZ_ID_max + 1), (EZ_ID_max + length(clump_ls) - 1), length.out = (length(clump_ls) - 1)))
            EZ_clump_mat <- as.matrix(EZ_clump_raster,wide=TRUE)
            rm(EZ_clump_raster)
            for(new_clump in clump_ls[-1]){
              new_clump_sub_IDs <- which(EZ_clump_mat[] == new_clump, arr.ind = TRUE)
              new_clump_IDs <- cbind(new_clump_sub_IDs[,1] + min_EZ_row - 1, new_clump_sub_IDs[,2] + min_EZ_col - 1)
              EZ_mat[new_clump_IDs] <- new_EZ_IDs_vec[which(clump_ls == new_clump)]
            }
            EZ_ID_max <- max(new_EZ_IDs_vec) + 1
            new_EZ_pop <- sapply(new_EZ_IDs_vec, function(EZ_index){sum(subpop_mat[which(EZ_mat == EZ_index)])})
            new_EZ_size <- sapply(new_EZ_IDs_vec, function(EZ_index){length(which(EZ_mat == EZ_index))})
            small_new_EZs <- new_EZ_IDs_vec[which(new_EZ_size < (0.75 * max_cells_per_EZ) & (new_EZ_pop < min_pop_per_EZ))]
            for(small_EZ in small_new_EZs){
              # neighbour EZs
              EZ_cell_IDs <- which(EZ_mat == small_EZ, arr.ind = TRUE)
              nb_EZs <- cell_nbs <- NULL
              for(df_row in 1:nrow(EZ_cell_IDs)){
                loc_vert <- EZ_cell_IDs[df_row,1]
                loc_horiz <- EZ_cell_IDs[df_row,2]
                adj_cells <- cbind(c((loc_vert + 1), (loc_vert - 1), loc_vert, loc_vert), c(loc_horiz, loc_horiz, (loc_horiz + 1), (loc_horiz - 1)))
                if(TRUE %in% (!adj_cells[,1] %in% 1:nrow(subpop_mat))){adj_cells <- adj_cells[-which((!adj_cells[,1] %in% 1:nrow(subpop_mat))),]}
                if(TRUE %in% (!adj_cells[,2] %in% 1:ncol(subpop_mat))){adj_cells <- adj_cells[-which((!adj_cells[,2] %in% 1:ncol(subpop_mat))),]}
                cell_nbs <- EZ_mat[adj_cells]
                nb_EZs <- c(nb_EZs, cell_nbs)
              }
              nb_EZs <- unique(nb_EZs)
              nb_EZs <- nb_EZs[!is.na(nb_EZs) & !nb_EZs == small_EZ]
              if(length(nb_EZs) > 0){
                # measure compactness for each of the possible combinations of blocks/EZs
                compact_vals <- sapply(nb_EZs,function(new_EZ){
                  test_EZ_locs <- which(EZ_mat == EZ_ID | EZ_mat == new_EZ ,arr.ind = TRUE)
                  test_mat <- matrix(NA,nrow=dim(EZ_mat)[1],ncol=dim(EZ_mat)[2])
                  test_mat[test_EZ_locs] <- 1
                  test_poly <- terra::as.polygons(terra::rast(test_mat,ext = terra::ext(subpop_rast),crs = terra::crs(subpop_rast)), dissolve = TRUE)
                  poly_area <- terra::expanse(test_poly, unit = "m")
                  edge_peri <- terra::perim(test_poly)
                  (4 * pi * poly_area / (edge_peri^2))
                })
                EZ_ID_to_add <- (nb_EZs[compact_vals == max(compact_vals)])[1]   # compactness is likely most important when adding on smaller chunks, so don't consider pop size here
                EZ_mat[EZ_mat == small_EZ] <- EZ_ID_to_add
              }
            }
          }
       }
      }

    ####################
    # create EZ raster #
    ####################

    EZs_raster <- subpop_rast
    names(EZs_raster) <- "EZ_ID"
    terra::values(EZs_raster) <- NA
    terra::values(EZs_raster) <- as.vector(t(EZ_mat))

    ###################
    # save EZs raster #
    ###################

    unlink(paste(output_path,"/temp_folder_", run_ID, "/temp_", tempfileID, ".RData",sep=""), recursive = TRUE)
    terra::writeRaster(EZs_raster,paste(output_path,"/temp_folder_", run_ID, "/EZs_", tempfileID, ".tif",sep=""))
  }
  # end of EZ generation function

  ###############################################################
  # run EZ generation function across specified number of cores #
  ###############################################################

  cl <- parallel::makeCluster(ncores, type="PSOCK")
  print("Starting generation of EZs")
  parallel::clusterEvalQ(cl, {
    library(terra)
    library(data.table)
  })
  invisible(parallel::parLapply(cl, clump_extents, EZ_gen_fn,
                      target_pop_per_EZ=target_pop_per_EZ,
                      max_cells_per_EZ=max_cells_per_EZ,
                      output_path = output_path, run_ID = run_ID,
                      population_raster = population_raster,
                      strata2_raster = strata2_raster,
                      strata1_raster = strata1_raster))
  parallel::stopCluster(cl)
  print("EZ generation complete. Compiling output files")

  ###############################################################
  #  load EZ files, and create EZ ID and EZ population rasters  #
  ###############################################################

  # list of individual EZ rasters for each strata1-strata2 combination
  current_list <- list.files(path = paste(output_path,"/temp_folder_", run_ID, sep=""), pattern =".tif$", full.names=TRUE)
  rasters <- lapply(current_list, terra::rast)

  # unique IDs across all EZ rasters
  ids <- lapply(rasters, function(r) sort(unique(stats::na.omit(terra::values(r))))) # EZ ids in each raster
  EZ_counts <- sapply(ids,length) # number of unique values to get new non-overlapping ID sets per raster
  EZ_counts_fin <- cumsum(EZ_counts) # vector of final EZ id values for each raster
  EZ_counts_init <- c(1,EZ_counts_fin[-length(EZ_counts_fin)]+1) # vector of initial EZ id values for each raster
  rasters_updated_ids <- lapply(1:length(rasters), function(i){
    terra::subst(rasters[[i]], from = ids[[i]], to = EZ_counts_init[i]:EZ_counts_fin[i])
  })
  rm(rasters,current_list,ids,EZ_counts,EZ_counts_fin,EZ_counts_init)

  # merge and save master raster
  merged_raster <- do.call(terra::merge, rasters_updated_ids)
  terra::writeRaster(merged_raster, paste(output_path,"/EZ_IDs_", run_ID, ".tif",sep=""))

  # create pop per EZ table and save as csv
  pop_rast <- terra::crop(terra::rast(population_raster),merged_raster,extend=TRUE)
  dt <- data.table::data.table(EZ_ID = as.vector(merged_raster), pop = as.vector(pop_rast))
  dt <- dt[!is.na(dt[["EZ_ID"]]), ]
  cell_count_dt <- data.table::as.data.table(table(dt[["EZ_ID"]]))
  names(cell_count_dt) <- c("EZ_ID", "N")
  cell_count_dt[["EZ_ID"]] <- as.numeric(cell_count_dt[["EZ_ID"]])
  pop_count_dt <- data.table::as.data.table(stats::aggregate(pop ~ EZ_ID, data = dt, FUN = function(x) sum(x, na.rm = TRUE)))
  dt <- merge(pop_count_dt, cell_count_dt, by='EZ_ID')
  utils::write.csv(dt,paste(output_path,"/EZ_Pop_Ncells_",run_ID,".csv",sep=""),row.names = FALSE)

  unlink(paste(output_path,"/temp_folder_", run_ID, sep=""), recursive = TRUE)
  print("Output files saved and temporary folder deleted. gridEZ has successfully finished")
  invisible(NA)
}
