################################################################
##get k-clique subGraph
getKcSubGraph<-function(k=4,graphList=getDefaultUndirectedGraph()){
      if(k<1)
            stop("k can't <1 ")
      graphList<-graphList[sapply(graphList,function(x) length(x)>0)]
      subGraphIndex<-0
      subGraph<-list()
      subNames<-character()
      for(i in 1:length(graphList)){
         kc<-kCliques(graphList[[i]])
         if(length(kc)>0){
            if(k<=length(kc)){
                  kkc<-kc[[k]]
            }
            else{
                  kkc<-kc[[length(kc)]]
            }

            for(j in 1:length(kkc)){
                  subGraphIndex<-subGraphIndex+1
                  subGraph[subGraphIndex]<-subGraph(kkc[[j]],graphList[[i]])
                  subNames[subGraphIndex]<-paste(names(graphList)[i],j,sep="_")
            }
          }
      }
      names(subGraph)<-subNames
      return(subGraph)
}
#subGraphList<-getKnnSubGraph()
