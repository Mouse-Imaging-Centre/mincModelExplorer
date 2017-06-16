
#' Create Minc Model Explorer App
#' 
#' Read in pydpiper files and directories and produce an appropriate app directory
#'
#' @param title The name of the app
#' @param out_dir Where to create the app directory
#' @param metadata A file containing the data needed to fit models of
#' absolute and relative jacobians. Can include a labels column indicating
#' MAGeT labelled data files
#' @param defs A definitions file
#' @param model_formula The basic model to run
#' @param blur_level For filtering pydpiper output
#' @param relevels A named character vector of column_name = base_factor_level
#' defaults to c("genotype" = "WT")
#' @param plot_choices values to offer the use for plot choices
#' @param label_data A file containing paths to MAGeT labelled data files
#' @param study_average Where to find the final average, defaults to the nlin-3 average
#' from pydpiper
#' @param mask A mask to use for the analysis, defaults to the nlin-3 mask
#' @param clobber whether to overwrite files in an existing app directory
#' @return The final merged data frame invisibly.
#' @export
make_explorer_app <-
    function(title
           , pydpiper_dir
           , out_dir
           , metadata
           , defs
           , model_formula = Filenames ~ genotype
           , config_file = system.file("config.R", package = "mincModelExplorer")
           , blur_level = 0.2
           , input_file_col = "Filenames"
           , relevels = c("genotype" = "WT")
           , plot_choices = "genotype"
           , label_data = NULL
           , study_average = NULL
           , mask = NULL
           , clobber = FALSE){

        # Fix strings as factors
        old_opts <- options()
        on.exit({ options(old_opts) })
        options(stringsAsFactors = FALSE)

        #Setup output directory structure
        if(file.exists(out_dir) & !clobber)
            stop("The output directory exists already, specify clobber to overwrite")

        if(!file.exists(out_dir))
            dir.create(out_dir)

        if(!dir.exists(file.path(out_dir, "volumes")))
            dir.create(file.path(out_dir, "volumes"))

        # Copy over app files
        lapply(c("global.R", "screenplot.R", "server.R", "ui.R")
             , function(f) file.copy(system.file(file.path("app_files", f)
                                               , package = "minc_model_explorer")
                                     , out_dir))

        # Global options
        globalOptions <- list(title = title, plotChoices = plot_choices)

        # Setup data frames
        dets_file <- file.path(pydpiper_dir, "determinants.csv")
        xfms_file <- file.path(pydpiper_dir, "transforms.csv")

        if(!file.exists(dets_file) | !file.exists(xfms_file))
            stop("Unable to read determinants or transforms file from your pydpiper directory\n"
               , "are you sure you have the right path")

        if(!file.exists(metadata))
            stop("Unable to read your metadata file")

        if(!is.null(label_data) && !file.exists(label_data))
            stop("Unable to find your maget data file")

        # Read in data frames
        dets <- read.csv(dets_file) %>% filter(fwhm == blur_level)
        xfms <- read.csv(xfms_file)
        metadata <- read.csv(metadata)

        # Check xfms has the native file, otherwise search for it
        if(is.null(xfms$native_file)){
            link_err <- paste0("Unable to link determinants to input files \n"
                       , "consider adding a `native_file` column to "
                       , xfms_file)
            
            native_files <- sub("\\.mnc", "", basename(metadata[[input_file_col]]))
            nf_matches <- vapply(native_files, function(f){
                matches <- grep(f, dets$log_full_det)
                if(length(matches) > 1) #possible dups in asym pipeline for example
                    matches <- grep(paste0(f, "_N_I"), dets$log_full_det)

                if(length(matches) != 1)
                    stop(link_err)

                matches
            }, numeric(1))

            if(any(duplicated(nf_matches)))
               stop(link_err)

            dets <- dets[nf_matches, ] ##removes indivs not in metadata
            dets$native_file <- metadata[[input_file_col]]
        } else {
            dets <- inner_join(dets, xfms, by = c("inv_xfm" = "lsq12_nlin_xfm"))
        }

        # Make native files absolute
        # check the first file, if it doesn't start with /, make all absolute
        if(!grepl("^/", dets$native_file[1])) 
            dets$native_file <- file.path(pydpiper_dir, dets$native_file)

        if(!grepl("^/", metadata[[input_file_col]][1]))
            metadata[[input_file_col]] <- file.path(pydpiper_dir, metadata[[input_file_col]]) 

        # Merge the data frames
        gfs <- inner_join(dets, metadata                        
                        , by = c("native_file" = input_file_col))

        # Check if we can run anatomy
        if(is.null(label_data) & is.null(gfs$labels))
            message("Running without structural analysis, if this is not what you want, supply a "
                    , "`label_data` argument, or add a `labels` column to `metadata`")

        # Locate and validate in mask and study average
        if(is.null(study_average))
            study_average <- Sys.glob(file.path(pydpiper_dir, "*_nlin/*-nlin-3.mnc"))

        if(length(study_average) != 1)
            stop("Trouble finding study average, please pass on the command line")

        if(is.null(mask))
            mask <- Sys.glob(file.path(pydpiper_dir, "*_nlin/*-nlin-3_mask.mnc"))

        if(length(study_average) != 1)
            stop("Trouble finding mask, please pass on the command line")

        
        #Read in average volume
        anatVol <- mincArray(mincGetVolume(study_average))
        d <- dim(anatVol)
        
        # Prepare formulae
        abs_model <- update(model_formula, log_full_det ~ .)
        rel_model <- update(model_formula, log_nlin_det ~ .)
        summary_model <- update(model_formula, voxel ~ .)

        ## Relevel according to the user supplied "dictionary" (really char vector)
        gfs <-
            Reduce(function(base, update){
                base[[update]] <- relevel(factor(base[[update]]), relevels[[update]])
                base
            }
          , names(relevels), gfs)

        # Move volume files
        abs_files <- gfs$log_full_det
        new_abs_files <- file.path("volumes", basename(abs_files))
        if(!grepl("^/", abs_files[1])) 
            abs_files <- file.path(pydpiper_dir, abs_files)

        rel_files <- gfs$log_nlin_det
        new_rel_files <- file.path("volumes", basename(rel_files)) 
        if(!grepl("^/", rel_files[1])) 
            rel_files <- file.path(pydpiper_dir, rel_files)

        mapply(file.copy, from = abs_files, to = file.path(out_dir, new_abs_files))
        mapply(file.copy, from = rel_files, to = file.path(out_dir, new_rel_files))

        # Set files for running locally
        gfs$log_full_det <- abs_files
        gfs$log_nlin_det <- rel_files
            
        # Run minclms
        vs <- mincLm(abs_model, gfs, mask=mask)
        vsrel <- mincLm(rel_model, gfs, mask=mask)

        # Compute FDR
        vsFDR<-mincFDR(vs)
        vsrelFDR<-mincFDR(vsrel)

        # Fix filenames for use in the app
        gfs$log_full_det <- new_abs_files
        gfs$log_nlin_det <- new_rel_files

        # Setup specific voxel model function
        ## stupid lugging of environments causes this function to explode in size
        ## and ruin everything. 
        modelfunc <- function(data) { summary(lm(summary_model, data)) }
        environment(modelfunc) <- globalenv()
        environment(modelfunc)$summary_model <- summary_model

        # Bring in config (I don't reaaallly like this, but the config is so ugly)
        source(config_file, local = TRUE)

        # Set up objects to be exported 
        export_file_args <-
            alist(file = file.path(out_dir, "data.rda")
                , globalOptions
                , anatVol
                , d
                , statsList
                , gfs)
        
        ## Anatomy, if we can
        if(!is.null(label_data) | !is.null(gfs$labels)){
            # Run anat gets
            allvols <- anatGetAll(gfs$labels, gfs$labels[1], method="labels"
                            , defs=defs)
            gfs$vols <- anatCombineStructures(allvols, defs=defs)

            # Compute brain volume
            brainvolumes <- rowSums(gfs$vols) 

            # Compute relative vols
            gfs$rel_vols <- gfs$vols/brainvolumes * 100

            # Run anatLms
            avs <- anatLm(~genotype, gfs, gfs$vols)
            avsrel <- anatLm(~genotype, gfs, gfs$rel_vols)

            # Run anatFDRs
            qavs <- anatFDR(avs)
            qavsrel <- anatFDR(avsrel)

            export_file_args <-
                c(export_file_args, alist(avs, avsrel, qavs, qavsrel))
        }

         

        # Save files required for the shiny app
        do.call(save, export_file_args)
        
        invisible(gfs)
    }



