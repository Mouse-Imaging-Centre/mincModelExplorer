
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)
library(RMINC)
library(plotrix)
library(ggplot2)
library(reshape2)
library(gridBase)
library(grid)

source("screenplot.R")

cat(names(statsList))

shinyServer(function(input, output, clientData, session) {
    source("server_utils.R", local = TRUE)
    
    # Setup drop-down with selectable stats
    updateSelectInput(session, "statistic", choices=names(statsList))
  
    # Initialize the description text
    descriptiveText <-
        reactiveValues(description="")
  
    # create a local copy of input data
    rvoxel <- reactive({ getVoxel() })
    rlocation2 <- reactive({ getLocation2() })
    rlocation1 <- reactive({ getLocation1() })
    rlocation3 <- reactive({ getLocation3() })

    # Determine the starting location for voxel references
    v <- reactiveValues(
        loc1 = ceiling(d[3]/2),
        loc2 = ceiling(d[2]/2),
        loc3 = ceiling(d[1]/2),
        voxel = NULL
    )

    # Check for plot clicks on a coronal slice
    observeEvent(input$plot_click, {
        location <- rlocation2()
        v$loc1 <- location[1]
        v$loc2 <- location[2]
        v$loc3 <- location[3]
        if (input$updatePlot) { v$voxel <- rvoxel() }
    })

    # Check for a clicks on a sagittal slice
    observeEvent(input$click_sagittal, {
        location <- getLocation1()
        v$loc1 <- location[1]
        v$loc2 <- location[2]
        v$loc3 <- location[3]
        if(input$updatePlot) { v$voxel <- rvoxel() }
    })

    # Check for clicks on an axial slcie
    observeEvent(input$click_axial, {
        location <- rlocation3()
        v$loc1 <- location[1]
        v$loc2 <- location[2]
        v$loc3 <- location[3]
        if(input$updatePlot) { v$voxel <- rvoxel() }
    })

    # If the slice series is clicked, print where
    observeEvent(input$sliceSeriesClick, {
        cat("SS click:", input$sliceSeriesClick$x, input$sliceSeriesClick$y, "\n")
    })
    
    # Adjust slider increment by current dimension size
    observe({
        dval <- as.integer(input$dimension)
        diff20 <- round(d[dval] * 0.2)
        cat("in observe slice range:", dval, d[dval], diff20, "\n")
        updateSliderInput(session, "sliceRange", max=d[dval], value=c(diff20, d[dval]-diff20))
    })

    # Adjust current stat maximum and description according to the
    # selected statistic
    observe({
        currentStat <- input$statistic
        updateTextInput(session, "ssTitle", value=input$statistic)
        maxstat <- round(max(abs(range(statsList[[currentStat]]$data))), 1)
        cat("in observe ", maxstat, "\n")
        if (!is.infinite(maxstat)) {
            updateSliderInput(session, "range", max=maxstat, value=c(2, maxstat))
            descriptiveText$description = statsList[[currentStat]]$description

        }
    })

    # Output the descriptive text
    output$statsDescription <- renderUI({
        p(descriptiveText$description)
    })

    # Render the slice series
    output$seriesPlot <- renderPlot({
        shinySliceSeries(input$statistic)
    }, bg="black")

    # Handle slice series downloader input
    output$ssDownload <- downloadHandler(
        filename = function() { paste("sliceSeries", tolower(input$ssFileType), sep=".") },
        content = function(file) {
            if (input$ssFileType == "PDF") {
                pdf(file=file, height=input$ssHeight, width=input$ssWidth, bg="black")#, units="in", res=72)
            }
            else if (input$ssFileType == "PNG") {
                png(file=file, height=input$ssHeight, width=input$ssWidth, bg="black", units="in", res=input$ssRes)
            }
            else if (input$ssFileType == "TIFF") {
                tiff(file=file, height=input$ssHeight, width=input$ssWidth, bg="black", units="in", res=input$ssRes)
            }
            shinySliceSeries(input$ssTitle)
            dev.off()
        }
        
    )
    
    # Output a plot for the current voxel data with base size 12
    output$graphPlot <- renderPlot({
        graphVoxelPlot(12)
    })

    # Handle triplanar slice plot rendering and download
    output$spDownload <- downloadHandler(
        filename = function() { paste("sliceAndPlot", tolower(input$spFileType), sep=".") },
        content = function(file) {
            if (input$spFileType == "PDF") {
                cairo_pdf(file=file, height=input$spHeight, width=input$spWidth, bg="black")#, units="in", res=72)
            }
            else if (input$spFileType == "PNG") {
                png(file=file, height=input$spHeight, width=input$spWidth, bg="black", units="in", res=input$spRes)
            }
            else if (input$spFileType == "TIFF") {
                tiff(file=file, height=input$spHeight, width=input$spWidth, bg="black", units="in", res=input$spRes)
            }
            
            # use split.screen for a complicated design
            lmatrix <- screenplotlayout()
            close.screen(all.screens=T)
            split.screen(lmatrix)
            # first screen: axial slice
            screen(1)
            par(mar=c(0,0,0,0), bg="black")
            plotAxialSlice()
            # second screen: coronal slice
            screen(2)
            par(mar=c(0,0,0,0))
            plotCoronalSlice()
            # third screen: sagittal slice
            screen(3)
            par(mar=c(0,0,0,0))
            plotSagittalSlice()
            # fourth screen: colour legend.
            # NOTE: still copied too much code from RMINC
            screen(4)
            par(mar=c(0,0,0,0))
            
            rcol <- colorRampPalette(c("blue", "turquoise1"))(255)
            high <- input$range[2]
            low <- input$range[1]
            legend <- statsList[[input$statistic]]$legend
            if (statsList[[input$statistic]]$symmetric==TRUE) {
                col <- colorRampPalette(c("red", "yellow"))(255)
                color.legend(0, 0.05, 0.2, 0.45, c(high*-1, low*-1), rev(rcol), gradient="y", align="rb", col="white")
                color.legend(0, 0.55, 0.2, 0.95, c(low, high), col, gradient="y", align="rb", col="white")
                text(0.8, 0.5, labels=legend, srt=90, col="white")
            }
            else {
                col <- rainbow(255)
                color.legend(0, 0.25, 0.2, 0.75, c(low, high), col, gradient="y", align="rb", col="white")
                text(0.8, 0.5, labels=legend, srt=90, col="white")
            }
            
            # fifth screen: the plot
            # uses baseViewports to get the viewport to plot into.
            screen(5)
            plot.new()
            vps <- baseViewports()
            pushViewport(vps$figure)
            vp1 <- plotViewport(c(0,0,0,0))
            p <- graphVoxelPlot(8)
            print(p, vp=vp1)
            close.screen(all.screens=T)
            dev.off()
            #close.screen(all.screens=T)
        })

    # Render the coronal slice
    output$coronalPlot <- renderPlot({
        plotCoronalSlice()
    }, bg="black")

    # Render the saggital slice
    output$sagittalPlot <- renderPlot({
        par(mar=c(0,0,0,0))
        plotSagittalSlice()  
    }, bg="black")

    # Render the axial slice
    output$axialPlot <- renderPlot({
        par(mar=c(0,0,0,0))
        plotAxialSlice()
    }, bg="black")


    # Output the results of printing modelfunc applied to gfs
    # and the current voxel data
    output$summaryText <- renderPrint({
        gfs$voxel <- v$voxel 
        #func <- statsList[[input$statistic]]$modelfunc
        #environment(func) <- 
        statsList[[input$statistic]]$modelfunc(gfs)
    })

    # Render a data table of anatomy results
    output$volumesTable <- DT::renderDataTable({
        
        usableNames <- unique(statsList[[input$statistic]][c("xvar", "colour", "fill")])
        isFalse <- usableNames %in% FALSE
        usableNames <- usableNames[!isFalse]
        gfs$newGrouping <- ""
        for (i in 1:length(usableNames)) {
            gfs$newGrouping <- paste(gfs$newGrouping, gfs[,usableNames[[i]]])
            cat("usableNames ", i, ": ", usableNames[[i]], "\n", sep = "")
        }
        cat("\n")
        
        gfs$newGrouping <- factor(gfs$newGrouping)
        cat(gfs$newGrouping[1])
        
        dt <- as.data.frame(t(apply(gfs$vols, 2, function(x) { tapply(x, gfs[,"newGrouping"], mean)})))
        dtsd <- as.data.frame(t(apply(gfs$vols, 2, function(x) { tapply(x, gfs[,"newGrouping"], sd)})))
        l <- levels(gfs[,"newGrouping"])
        for (i in 2:length(l)) {
            cat("pre dt\n")
            dt[,paste("effect size:", l[i])] <- (dt[,l[i]] - dt[,l[1]]) / dtsd[,l[1]]
            cat("post dt\n")
        }
        
        colnames(qavs) <- paste("q-value", colnames(qavs))
        dt <- cbind(qavs, dt)
        dt
        
    })

    # Render the volume/anatomy plot
    output$volumesPlot <- renderPlot({
        selectedStructures <- colnames(gfs$vols)[input$volumesTable_rows_selected]
        usableNames <- unique(statsList[[input$statistic]][c("xvar", "colour", "fill")])

        isFalse <- usableNames %in% FALSE
        usableNames <- usableNames[!isFalse]

        v <- as.data.frame.matrix(gfs$vols)
        for (i in 1:length(usableNames)) {
            v[,usableNames[[i]]] = gfs[,usableNames[[i]]]
        }
        usableNames <- unlist(usableNames)

        cat(" ", usableNames[1], " \n")
        m <- melt(v, id.vars=usableNames) 

        m <- subset(m, variable %in% selectedStructures)

        cat(names(m), "\n")
        cat(nrow(m), "\n")

        
        update_geom_defaults("point", list(colour = "black"))
        update_geom_defaults("errorbar", list(colour = "black"))
        update_geom_defaults("boxplot", list(colour = "black", fill="transparent"))
        graphData(m$value, "Volume", m) + theme_classic() + facet_wrap(~variable, scales="free_y") 
    })

})
