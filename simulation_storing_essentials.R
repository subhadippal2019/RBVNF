




## Codes for storing the pivotal quantities

Run_Time=Sys.time()-Start_Time
#Sys_info=get_sys_details()


###################### put this chunk inside a function
function_name=as.character(match.call()[1])
function_def0<-  capture.output(print(get( function_name)));
function_def0[1]=paste0( function_name,"<-",function_def0[1])
function_def=paste( function_def0 , collapse = "\n")
#function_def1=dput(get(function_name))
#########

Sys_info=Sys.info()

# MC<-list(Mc_Beta=Store_Beta,
# Mc_sigma_sq=Store_sigma_sq,
# Mc_xi=Store_xi,
# Mc_eta=Store_eta,
# Mc_Nu=Store_Nu)


MC<-list(         Mc_Beta=  Store_Beta,
                  Mc_sigma_sq=Store_sigma_sq,
                  Mc_xi=Store_xi,
                  Mc_lambda=Store_lambda,
                  Mc_Nu=Store_Nu)



lst<-list(MC=MC,
          Data=list(y=y, X=X),
          Methodtype="Regularized",
          Run_Time=Run_Time,
          System_info=Sys_info,
          call=match.call(),
          function_def=function_def)



