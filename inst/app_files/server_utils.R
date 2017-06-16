graphData <- function(ydata, ycaption, inputdata) {
    inputdata$data = ydata
    inputdata$xvar <- gfs[,input$xvar]
    if (input$colour != "None") {
        inputdata$colour <- gfs[,input$colour]
    }
    if (input$fill != "None") {
        inputdata$fill <- gfs[,input$fill]
    }    

    mywidth=0.5
    myposition <- position_dodge()
    
    if (input$graphType == "point")
        myposition <- position_jitterdodge(dodge.width=mywidth) 
    
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
        p <- p + stat_summary(fun.data=mean_cl_boot
                            , geom="errorbar"
                            , width=0.3
                            , size=2
                            , position=position_dodge(width=mywidth))
        
    }
    else {
        cat("not point\n", input$graphType, "\n")
    }
    
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
    location <- c(v$loc1, location[2], location[1])
    return(location)
}

getLocation2 <- function() {
    location <- c(input$plot_click$x, input$plot_click$y)
    location <- c(location[2], v$loc2,location[1])
    return(location)
}

getLocation1 <- function(){
    location <- c(input$click_sagittal$x, input$click_sagittal$y)
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

graphVoxelPlot <- function(baseSize=12) {
    update_geom_defaults("point", list(colour = "white"))
    update_geom_defaults("errorbar", list(colour = "white"))
    update_geom_defaults("boxplot", list(colour = "white", fill="transparent"))
    
    graphData(exp(v$voxel), "jacobians", gfs) + theme_black(baseSize) + theme(legend.position="top") + 
        guides(colour = guide_legend(nrow=1), fill=guide_legend(nrow=1))
}

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
