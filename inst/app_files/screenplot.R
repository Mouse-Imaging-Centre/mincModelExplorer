# uses split.screen to lay out a plot with 5 elements (3 cardinal slices,
# a plot of an indicated voxel, and the colour-bar)

# hardcode the layout for now
screenplotlayout <- function() {
  lmatrix <- rbind(c(0, 315, 0, 478), # axial
                   c(315, 315+315, 0, 241), # coronal
                   c(315, 315+478, 241, 241+241), # sagittal
                   c(315+458, 315+478+40, 241, 241+241), # legend
                   c(315+315, 315+478+40, 0, 241)) # plot
  lmatrixn <- lmatrix
  lmatrixn[,1:2] <- lmatrix[,1:2] / max(lmatrix[,1:2])
  lmatrixn[,3:4] <- lmatrix[,3:4] / max(lmatrix[,3:4])
  return(lmatrixn)
}
