
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

# theme_black, stolen from: https://jonlefcheck.net/2013/03/11/black-theme-for-ggplot2-2/
# with a few minor modifications to make it look more like theme_classic
theme_black=function(base_size=12,base_family="") {
  theme_grey(base_size=base_size,base_family=base_family) %+replace%
    theme(
      # Specify axis options
      axis.line=element_line(colour="white"), #element_blank(), 
      axis.text.x=element_text(size=base_size*0.8,color="white",
                               lineheight=0.9,vjust=1), 
      axis.text.y=element_text(size=base_size*0.8,color="white",
                               lineheight=0.9,hjust=1), 
      axis.ticks=element_line(color="white",size = 0.2), 
      axis.title.x=element_text(size=base_size,color="white",vjust=1, margin=margin(5,0,0,0)), 
      axis.title.y=element_text(size=base_size,color="white",angle=90,
                                vjust=0.5, margin=margin(0,5,0,0)), 
      axis.ticks.length=unit(0.3,"lines"), 
      axis.ticks.margin=unit(0.5,"lines"),
      # Specify legend options
      legend.background=element_rect(color=NA,fill="black"), 
      legend.key=element_rect(color=NA, fill="black"), 
      legend.key.size=unit(1.2,"lines"), 
      legend.key.height=NULL, 
      legend.key.width=NULL,     
      legend.text=element_text(size=base_size*0.8,color="white"), 
      legend.title=element_text(size=base_size*0.8,face="bold",hjust=0,
                                color="white"), 
      legend.position="right", 
      legend.text.align=NULL, 
      legend.title.align=NULL, 
      legend.direction="vertical", 
      legend.box=NULL,
      # Specify panel options
      panel.background=element_rect(fill="black",color = NA), 
      panel.border=element_blank(), #element_rect(fill=NA,color="white"), 
      panel.grid.major=element_blank(), 
      panel.grid.minor=element_blank(), 
      panel.margin=unit(0.25,"lines"),  
      # Specify facetting options
      strip.background=element_rect(fill="grey30",color="grey10"), 
      strip.text.x=element_text(size=base_size*0.8,color="white"), 
      strip.text.y=element_text(size=base_size*0.8,color="white",
                                angle=-90), 
      # Specify plot options
      plot.background=element_rect(color="black",fill="black"), 
      plot.title=element_text(size=base_size*1.2,color="white"), 
      plot.margin=unit(c(1,1,0.5,0.5),"lines")
    )
}



cat(names(statsList))

shinyServer(function(input, output, clientData, session) {
  
  # update some elements based on the data being accessed
  updateSelectInput(session, "statistic", choices=names(statsList))
  
  # a function to plot either voxels or volumes based on the input data, colours, fill, etc.
  graphData <- function(ydata, ycaption, inputdata) {
    inputdata$data = ydata
    inputdata$xvar <- gfs[,input$xvar]#inputdata[,statsList[[input$statistic]]$xvar]
    if (input$colour != "None") {
      inputdata$colour <- gfs[,input$colour]#inputdata[,statsList[[input$statistic]]$colour]
    }
    if (input$fill != "None") {
      inputdata$fill <- gfs[,input$fill]
    }    
#    if (statsList[[input$statistic]]$fill == FALSE) {
#      p <- qplot(xvar, ydata, data=inputdata, geom=input$graphType, colour=colour) + 
#        xlab(input$xvar) +
#        ylab(ycaption) +
#        scale_colour_brewer(statsList[[input$statistic]]$colour, palette="Set1")
#    }
#    else {
#      inputdata$fill <- inputdata[,statsList[[input$statistic]]$fill]

    mywidth=0.5
    myposition <- position_dodge()
    if (input$graphType == "point") { myposition <- position_jitterdodge(dodge.width=mywidth) }
    
    p <- do.call(qplot, c(list(inputdata$xvar
                             , inputdata$data
                             , data=inputdata
                             , geom=input$graphType
                             , position=myposition),
                          list(colour=inputdata$colour)[input$colour!="None"], 
                          list(fill=inputdata$fill)[input$fill!="None"])) +
      
      xlab(input$xvar) +
      ylab(ycaption) +
      scale_colour_brewer(input$colour, palette="Set1") +
      scale_fill_grey(input$fill)
    if (input$graphType == "point") {
      cat("is point\n")
      p <- p + stat_summary(fun.data=mean_cl_boot, geom="errorbar", width=0.3, size=2, position=position_dodge(width=mywidth))
      
    }
    else {
      cat("not point\n", input$graphType, "\n")
    }
#    }
    return(p)
  }   
  
  # a function to plot the slice series
  shinySliceSeries <- function(plottitle) {
    mincPlotSliceSeries(anatVol, statsList[[input$statistic]]$data,
                        anatLow=700, anatHigh=1400, low=input$range[1], high=input$range[2], 
                        begin=input$sliceRange[1], end=input$sliceRange[2], plottitle = plottitle, 
                        dim=as.integer(input$dimension), symmetric=statsList[[input$statistic]]$symmetric,
                        legend=statsList[[input$statistic]]$legend, mfrow=c(input$rows, input$columns))
  }  
  
    
  getLocation3 <- function() {
    location <- c(input$click_axial$x, input$click_axial$y)
    #location <- ceiling(location * c(d[1], d[2]))
    location <- c(v$loc1, location[2], location[1])
    return(location)
  }
  
  getLocation2 <- function() {
    location <- c(input$plot_click$x, input$plot_click$y)
    #location <- ceiling(location * c(d[1], d[3]))
    location <- c(location[2], v$loc2,location[1])
    return(location)
  }
  
  getLocation1 <- function(){
    location <- c(input$click_sagittal$x, input$click_sagittal$y)
    #location <- ceiling(location * c(d[2], d[3]))
    location <- c(location[2],location[1],v$loc3)
    return(location)
  }
  getVoxel <- function() {
    cat(file=stderr(), "in getVoxel\n")
    voxel <- mincGetVoxel(statsList[[input$statistic]]$filenames, v$loc1, v$loc2, v$loc3)
    cat(file=stderr(), "File1:", statsList[[input$statistic]]$filenames[1], "\n")
    cat(file=stderr(), "File2:", statsList[[input$statistic]]$filenames[2], "\n")
    return(voxel)
  }

  descriptiveText <- reactiveValues(
    description="some description")
  
  # create a local copy of input data
  rvoxel <- reactive({ getVoxel() })
  rlocation2 <- reactive({ getLocation2() })
  rlocation1 <- reactive({ getLocation1() })
  rlocation3 <- reactive({ getLocation3() })
  
  v <- reactiveValues(
    loc1 = ceiling(d[3]/2),
    loc2 = ceiling(d[2]/2),
    loc3 = ceiling(d[1]/2),
    voxel = NULL
  )
  
  observeEvent(input$plot_click, {
    location <- rlocation2()
    v$loc1 <- location[1]
    v$loc2 <- location[2]
    v$loc3 <- location[3]
    if (input$updatePlot) { v$voxel <- rvoxel() }
  })
  observeEvent(input$click_sagittal, {
    location <- getLocation1()
    v$loc1 <- location[1]
    v$loc2 <- location[2]
    v$loc3 <- location[3]
    if(input$updatePlot) { v$voxel <- rvoxel() }
  })
  observeEvent(input$click_axial, {
    location <- rlocation3()
    v$loc1 <- location[1]
    v$loc2 <- location[2]
    v$loc3 <- location[3]
    if(input$updatePlot) { v$voxel <- rvoxel() }
  })
  
  observeEvent(input$sliceSeriesClick, {
    cat("SS click:", input$sliceSeriesClick$x, input$sliceSeriesClick$y, "\n")
  })
  
  # update user interface elements based on other selections
  observe({
    # update min and max slice based on the dimension being displaced
    dval <- as.integer(input$dimension)
    diff20 <- round(d[dval] * 0.2)
    cat("in observe slice range:", dval, d[dval], diff20, "\n")
    updateSliderInput(session, "sliceRange", max=d[dval], value=c(diff20, d[dval]-diff20))
    #updateSliderInput(session, "end", min=-d[dval])
  })
  
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

  output$statsDescription <- renderUI({
    p(descriptiveText$description)
  })
  
  output$seriesPlot <- renderPlot({
    shinySliceSeries(input$statistic)
  }, bg="black")
  
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
  
  plotCoronalSlice <- function() {
    mincPlotAnatAndStatsSlice(anatVol,
                              statsList[[input$statistic]]$data,
                              anatLow=700, anatHigh=1400, 
                              low=input$range[1], high=input$range[2], 
                              slice=v$loc2, symmetric=statsList[[input$statistic]]$symmetric,
                              dim=2)    
    abline(h=v$loc1)
    abline(v=v$loc3) 
  }
  
  output$coronalPlot <- renderPlot({
    plotCoronalSlice()

  }, bg="black")
  
  plotSagittalSlice <- function() {
    mincPlotAnatAndStatsSlice(anatVol,
                              statsList[[input$statistic]]$data,
                              anatLow=700, anatHigh=1400, 
                              low=input$range[1], high=input$range[2], 
                              slice=v$loc3, symmetric=statsList[[input$statistic]]$symmetric,
                              dim=1, legend=statsList[[input$statistic]]$legend)    
    abline(h=v$loc1)
    abline(v=v$loc2)
  }
  
  output$sagittalPlot <- renderPlot({
    par(mar=c(0,0,0,0))
    plotSagittalSlice()  
  }, bg="black")
  
  plotAxialSlice <- function() {
    par(mar=c(0,0,0,0))
    mincPlotAnatAndStatsSlice(anatVol,
                              statsList[[input$statistic]]$data,
                              anatLow=700, anatHigh=1400, 
                              low=input$range[1], high=input$range[2], 
                              slice=v$loc1, symmetric=statsList[[input$statistic]]$symmetric,
                              dim=3)    
    abline(v=v$loc3)
    abline(h=v$loc2)
  }
  
  output$axialPlot <- renderPlot({
    par(mar=c(0,0,0,0))
    plotAxialSlice()
  }, bg="black")
    
  graphVoxelPlot <- function(baseSize=12) {
    update_geom_defaults("point", list(colour = "white"))
    update_geom_defaults("errorbar", list(colour = "white"))
    update_geom_defaults("boxplot", list(colour = "white", fill="transparent"))
    
    graphData(exp(v$voxel), "jacobians", gfs) + theme_black(baseSize) + theme(legend.position="top") + 
      guides(colour = guide_legend(nrow=1), fill=guide_legend(nrow=1))
  }
  
  output$graphPlot <- renderPlot({
    graphVoxelPlot(12)
  })

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
    
  output$summaryText <- renderPrint({
    #location <- c(input$plot_click$x, input$plot_click$y)
    #location <- ceiling(location * c(d[1], d[3]))
    gfs$voxel <- v$voxel #mincGetVoxel(gfs$reljacobians02, location[2], input$slice, location[1])
    
    #anova(lm(voxel ~ mouse.gender + Neonatal, gfs))
    statsList[[input$statistic]]$modelfunc(gfs)
  })
  
  output$volumesTable <- DT::renderDataTable({
    
    usableNames <- unique(statsList[[input$statistic]][c("xvar", "colour", "fill")])
    isFalse <- usableNames %in% FALSE
    usableNames <- usableNames[!isFalse]
    gfs$newGrouping <- ""
    for (i in 1:length(usableNames)) {
      gfs$newGrouping <- paste(gfs$newGrouping, gfs[,usableNames[[i]]])
      cat("usableNames i: ", usableNames[[i]], " ")
    }
    cat("\n")
    
    #gfs$newGrouping <- paste(gfs[,as.vector(usableNames)], sep="::")
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
    #dt$'effect size' <- (dt[,"male"] - dt[,"female"]) / dtsd[,"male"]
    cat("Huh", colnames(qavs), "\n")
    colnames(qavs) <- paste("q-value", colnames(qavs))
    dt <- cbind(qavs, dt)
    dt
    
  })
  output$volumesPlot <- renderPlot({
    selectedStructures <- input$volumesTable_rows_selected
    
    usableNames <- unique(statsList[[input$statistic]][c("xvar", "colour", "fill")])
    isFalse <- usableNames %in% FALSE
    usableNames <- usableNames[!isFalse]
    #cat(1)
    v <- as.data.frame(gfs$vols)
    for (i in 1:length(usableNames)) {
      v[,usableNames[[i]]] = gfs[,usableNames[[i]]]
    }
    usableNames <- unlist(usableNames)
    #cat(2)
    #v[,statsList[[input$statistic]]$xvar] = gfs[,statsList[[input$statistic]]$xvar]
    cat(" ", usableNames[1], " ")
    m <- melt(v, id.vars=usableNames) #statsList[[input$statistic]]$xvar)
    #cat(3)
    m <- subset(m, variable %in% selectedStructures)
    #cat(4)
    cat(names(m), "\n")
    cat(nrow(m), "\n")
    cat("5\n")
    
    update_geom_defaults("point", list(colour = "black"))
    update_geom_defaults("errorbar", list(colour = "black"))
    update_geom_defaults("boxplot", list(colour = "black", fill="transparent"))
    graphData(m$value, "Volume", m) + theme_classic() + facet_wrap(~variable, scales="free_y") 
    #qplot(m[,statsList[[input$statistic]]$xvar], value, colour=m[,statsList[[input$statistic]]$xvar],
    #      geom=input$graphType, data=m, ylab="Volume (mm3)", xlab=input$statistic) +
    #  scale_colour_brewer(input$statistic, palette="Set1") +
    #  facet_wrap(~variable, scales="free_y") 
    })
  #output$volumesPlot2 <- renderPlot({
  #  selectedStructure <- input$volumesTable_rows_selected[length(input$volumesTable_rows_selected)-1]
  #  qplot(mouse.gender, vols[,selectedStructure], 
  #        geom="boxplot", data=gfs, main=selectedStructure, ylab="Volume (mm3)")
  #})
})
