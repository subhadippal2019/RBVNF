


#' @export
load_packages<-function(){
#################################
library(gsl)
library(Bessel)
library(ghyp)
library(ggplot2)
library(Rfast)

require(cowplot)
library(grid)
library(gridExtra)
  library(mvtnorm) #install.packages("pgdraw")
  library(glmnet)
 # MCMCpack
  #Rfast::rvmf(n=1, mu=c(1,0,0), k=10)
 # library("rfast") #library("pgdraw")
  #library(nimble)
}
#library()
#install.packages("multicore")


#' @export
load_additional_packages<-function(){
  library(plotly)
  #library(emdbook)
  #Here, I load a map with the raster package.
  library(ggplot2)
  require(dplyr)
  library(tidyverse)
  #library(rgdal)
  #library(rayshader)
  library(ggfx)
}
