args <- commandArgs(TRUE)
if(length(args)==0){
    print("No arguments supplied.")
    ##supply default values
}else{
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}
library(carx)
source("./carx_zeger.R")
source("./simulation_setting.R")
source("./simSingle.R")

cwd <- getwd()
setwd(ind_result_dir)
print(iRep)
ret <- simSingle(iRep)
setwd(cwd)
