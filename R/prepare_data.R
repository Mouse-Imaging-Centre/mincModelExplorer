
#' Create Minc Model Explorer App
#' 
#' Read in pydpiper files and directories and produce an appropriate app directory
#' @export
make_explorer_app <-
    function(title
           , pydpiper_dir
           , out_dir
           , metadata
           , defs
           , model_formula = Filenames ~ genotype
           , config = system.file("config.R", package = "minc_model_explorer")
           , blur_level = 0.2
           , input_file_col = "mouse"
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
                                               , package = "minc_model_explorer")))

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

        if(!is.null(maget_data) & !file.exists(maget_data))
            stop("Unable to find your maget data file")

        # Read in data frames
        dets <- read.csv(dets_file) %>% filter(fwhm == blur_level)
        xfms <- read.csv(xfms_file)
        metadata <- read.csv(metadata)

        # Make native files absolute
        # check the first file, if it doesn't start with /, make all absolute
        if(!grepl("^/", xfms$native_file[1])) 
            xfms$native_file <- file.path(pydpiper_dir, xfms$native_file)

        if(!grepl("^/", metadata[[input_col_names]][1]))
            metadata[[input_col_names]] <- file.path(pydpiper_dir, metadata[[input_col_names]]) 

        # Merge the data frames
        gfs <- inner_join(metadata
                        , inner_join(dets, xfms, by = c("inv_xfm" = "lsq12_nlin_xfm"))
                          #small peice of cleverness here, `by` expects a character vector
                          #where the name is the LHS column, setNames avoids the uninterpolatable
                          # c(lname = "rname") problem. bquote and substitute *won't work*
                        , setNames("native_file", input_file_col))

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
                base[[update]] <- relevel(factor(base[[update]]), rlvls[[update]])
                base
            }
          , names(relevels), gfs)

        # Move volume files
        abs_files <- gfs$log_full_det
        new_abs_files <- file.path(out_dir, "volumes", basename(abs_files))
        if(!grepl("^/", abs_files[1])) 
            abs_files <- file.path(pydpiper_dir, abs_files)

        rel_files <- gfs$log_nlin_det
        new_rel_files <- file.path(out_dir, "volumes", basename(rel_files)) 
        if(!grepl("^/", rel_files[1])) 
            rel_files <- file.path(pydpiper_dir, rel_files)

        mapply(file.copy, from = abs_files, to = new_abs_files)
        mapply(file.copy, from = rel_files, to = rel_abs_files)

        gfs$log_full_det <- new_abs_files
        gfs$log_nlin_det <- new_rel_files
            
        # Run minclms
        vs <- mincLm(abs_model, gfs, mask=mask)
        vsrel <- mincLm(rel_model, gfs, mask=mask)

        # Compute FDR
        vsFDR<-mincFDR(vs)
        vsrelFDR<-mincFDR(vsrel)

        # Setup specific voxel model function
        modelfunc <- function(data) { summary(lm(summary_model, data)) }

        # Bring in config (I don't reaaallly like this, but the config is so ugly)
        source(config_file)

        # Set up objects to be exported 
        export_file_args <-
            alist(file = file.path(out_dir, "data.rds")
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



