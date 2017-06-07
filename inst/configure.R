statsList <-
    list(AbsVol=
             list(data = mincArray(vs, "tvalue-genotypeHT")
                , xvar = "genotype"
                , colour=FALSE
                , fill="genotype"
                , symmetric=TRUE
                , modelfunc=modelfunc
                ##### change this! 
                , description=paste("Absolute Volume Differences"
                                  , " - PTEN 167 (HT vs WT). "
                                  , "1% FDR - 3.4 and at 5% FDR - 2.7")
                , filenames=gfs$AbsJacobians
                , legend="T-Stats")
         
      ,  RelVol =
             list(data = mincArray(vsrel, "tvalue-genotypeHT")
                , xvar = "genotype"
                , colour=FALSE
                , fill="genotype"
                , symmetric=TRUE
                , modelfunc=modelfunc
                ######change this!
                , description=paste0("Relative Volume Differences"
                                   , " - PTEN 167 (HT vs WT). "
                                   , "1% FDR - 4.2 and at 5% FDR - 3.4")
               ,  filenames=gfs$RelJacobians,
                 legend="T-Stats"))
