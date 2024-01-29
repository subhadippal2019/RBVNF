


#' plot_Acf_violiin_Slice
#' @param data :  A dataframe of size n*lag X number of different kappa seq. ACF_simNum_d_samSize[,,,1]
#' @examples
#' load("/Users/subhadippal/Documents/Google Drive/MyRPackages/MC_SIM_Reg_Dir_Data/SimulationSummary_ACF.Rdata")
#' plot_Acf_violiin_RBVNF(ACF_simNum_d_samSize[,,,1], Ylimit=c(.95, 1))
#' plot_Acf_violiin_RBVNF(ACF_simNum_d_samSize[,,,1], size=.2, ifGray=TRUE)
#' plot_Acf_violiin_RBVNF(ACF_simNum_d_samSize[,,,1], size=.2, ifGray=TRUE, GridLines="Vertical")
#' plot_Acf_violiin_RBVNF(ACF_simNum_d_samSize[,,,1],size=.2, ifGray=TRUE, GridLines="none")
#' @export

plot_Acf_violiin_RBVNF<-function(data, MAXLAG=15,scale_type="width",size=.5, Ylimit=NULL, ifGray=TRUE, GridLines="Vertical"){
  #MAXLAG=MAXLAG+1

  #MAXLAG=15
  lagsInd=(as.vector(replicate(100, 1:(MAXLAG))))
  temp_dd<-data[2:(MAXLAG+1),,]
  dd<-as.data.frame(matrix(temp_dd, nrow=prod(dim(temp_dd)[-3])))
  dd$lagsInd= (lagsInd)
  names(dd)=c("d_02", "d_03", "d_05", "d_10", "Lag")
  var_names=names(dd);  var_names1=gsub("d_", "d=", var_names)
  n_var=length(var_names)
  p <- ggplot(dd, aes(factor(get(var_names[5])), get(var_names[1])))
  #browser()
  #if(n_var>=3){
    code_str=paste0("p<-p+geom_violin(aes(y =", (var_names[1:(n_var-1)]),", fill = '",var_names1[1:(n_var-1)],"'),size=size, alpha = 0.7,scale =  scale_type)")
    eval(parse(text=code_str))
  #}
  p=p + theme( panel.grid.minor = element_blank())+ theme(legend.title=element_blank())+theme(legend.position="bottom")
  p=p+xlab("Lag Numbers")+ylab("Autocorrelation")
  if( GridLines=="none"){p=p+theme(panel.grid.major = element_blank())}
  if( GridLines=="Vertical"){p=p+theme(panel.grid.major.y = element_blank())}

  if(!is.null(Ylimit)){
    p=p+ylim(Ylimit[1], Ylimit[2])
  }
  if(ifGray){p=p+scale_fill_grey(start = 1, end = 0)}
  return( plt=p)
}




