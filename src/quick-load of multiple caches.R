##quickload caches

multilazyLoad <- function(path){
      for(f in unique(gsub("\\..+$", "", grep("RData", list.files(path, full.names = TRUE), value = T)))){
            print(f)
            lazyLoad(f)
      }
}
