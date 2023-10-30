check_pak <- function(pak_list){
install_fun <- function(need_install_pak=pak_list){
  for (f in need_install_pak){
    if(!require(f,character.only = TRUE)){
      install.packages(f)
      library(f,character.only = TRUE)
      } 
      else{
      check <- paste("loading",f)
      print(check)
      library(f,character.only = TRUE)
      }
  }
}
    tryCatch({
      install_fun(need_install_pak=pak_list)
    }, warning = function(war) {
      print(war$message)
    },error = function(err) {
      print(err$message)
     stop(err)
    })
}
