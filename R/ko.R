#######################################################
isKOMPath<-function(top){
   name<-names(top)
   for(i in 1:length(name)){
       if(name[i]=="reaction"){
           return (TRUE)
       }
   }
   return (FALSE)
}


################################################################################
##
getKOPathGraph<-function(top){
#top<-xmlRoot(xmlTreeParse(path)) 
  g1<-list()
  nodes<-getKONodes(top)
 if(length(nodes)>2){
  g1<-new("graphNEL",nodes,edgemode="undirected")
  if(length(top)>0&&(!is.na(top))){ 
      idKO<-list()   
      k<-0
      for(i in 1:length(top)){
       if(xmlName(top[[i]])=="entry"){
            attrs<-xmlAttrs(top[[i]])
            if(!is.na(attrs[["name"]])&&!is.na(attrs[["type"]])){
               if(attrs[["type"]]=="ortholog"){
                  tmp1<-unlist(strsplit(attrs[["name"]]," "))
                  if(length(tmp1)>0){ 
                     k<-k+1 
                     
                     idKO[[k]]<-tmp1
                     names(idKO)[k]<-attrs[["id"]]
                  }
               }
            }
       }
      }
      for(i in 1:length(top)){
       if(xmlName(top[[i]])=="relation"){
            attrs<-xmlAttrs(top[[i]])
            if(!is.na(attrs[["entry1"]])&&!is.na(attrs[["entry2"]])){

                one<-idKO[[attrs[["entry1"]]]]
                two<-idKO[[attrs[["entry2"]]]]
                if(length(one)>0&length(two)>0){
                 for(j in 1:length(one)){
                  for(k in 1:length(two)){
                   if(!isAdjacent(g1,one[j],two[k])){ 
                   g2<-addEdge(one[j],two[k],g1,1)
                   g1<-g2
                   }
                  }
                 }
                 if(length(one)>1){
                 for(j in 1:(length(one)-1)){
                  for(k in (j+1):length(one)){
                   if(!isAdjacent(g1,one[j],one[k])){ 
                   
                   g2<-addEdge(one[j],one[k],g1,1)
                   g1<-g2
                   }
                  }
                 }
                 }
                 if(length(two)>1){
                 for(j in 1:(length(two)-1)){
                  for(k in (j+1):length(two)){
                   if(!isAdjacent(g1,two[j],two[k])){ 
                   g2<-addEdge(two[j],two[k],g1,1)
                   g1<-g2
                   }
                  }
                 }
                 }
                }
                
            }
       }
      }     

  }
 }
 return (g1)
}
#top<-xmlRoot(xmlTreeParse("D:/r/kgml/ko04011.xml")) 
#g1<-getKOPathGraph("e:/ko/ko04011.xml")

#plot(g1,"neato")

#######################

#########################################################
##
getKONodes<-function(top){

      KOList<-character()
      KOListIndex<-0
      for(i in 1:length(top)){
            if(xmlName(top[[i]])=="entry"){
                  if(!is.na(xmlAttrs(top[[i]])["type"])){
                        if(!is.na(xmlAttrs(top[[i]])["name"])){
                              if(xmlAttrs(top[[i]])["type"]=="ortholog"){
                              tmp1<-unlist(strsplit(xmlAttrs(top[[i]])["name"]," "))
                              if(length(tmp1)>0){
                              for(k in 1:length(tmp1)){
                              KOListIndex<-KOListIndex+1
                              KOList[KOListIndex]<-tmp1[k]
                              }
                              }
                              }
                        }                       
                  }
            }
      }
      KOList<-intersect(KOList,KOList)
      return(KOList)
}

#top<-xmlRoot(xmlTreeParse("D:/r/kgml/ko04010.xml"))
#getKONodes(top)

#top<-xmlRoot(xmlTreeParse("D:/r/kgml/ko00010.xml"))
#getKONodes(top)


#########################################################
##get all pathway graph
getAllKOPathGraph<-function(fileList,path,verbose=TRUE){
    if(length(fileList)<1&&length(path)<1){
      print("fileList or path is error")
    }
    else{
      gList<-list()
      if(verbose==TRUE)
         print("Note that the programming is time consumming!!!!!!!!!!!!!!")
         print(paste("begin to build graph!!!!!!!!!!!!!!!!!!!!!!!!!"))
      for(i in 1:length(fileList)){
         if(verbose==TRUE)
            print(paste("building No.",i," graph. it is constructed from file ",fileList[i]))
         path1<-paste(path,fileList[i],sep="/")
         top<-xmlRoot(xmlTreeParse(path1)) 
         if(isKOMPath(top)==TRUE){
            if(verbose==TRUE)
              print("This is  a metabolic pathway")
            gList[[i]]<-getKOMPathGraph(top)
         }
         else{
            if(verbose==TRUE)
              print("This is  other kinds of pathway")
            gList[[i]]<-getKOPathGraph(top)  
         }
         if(verbose==TRUE){          
            if(length(gList[[i]])==0){
                  print("the graph is empty!")
            }
            else{
                  print(paste("finished. Nodes number:",numNodes(gList[[i]]),"Edges number:",numEdges(gList[[i]])))
            }
         }
      }
      if(verbose==TRUE)
         print(paste("finished!!!!!!!!!!!!!!!!!!!!!!!"))
      
      return(gList)
    }
}
########################################
#KOList<-list.files("e:/ko")
#g2<-getAllKOPathGraph(KOList[175:length(KOList)],"e:/ko")
#g2<-getAllKOPathGraph(KOList,"e:/ko")
#g2<-getAllKOPathGraph(KOList[1],"e:/ko")




###############################################################
##get pathway graph
getKOMPathGraph<-function(top){
      #top<-xmlRoot(xmlTreeParse(path))      
      nodes<-getKONodes(top)
      g1<-list()
      if(length(nodes)>2){
            g1<-new("graphNEL",nodes,edgemode="undirected")
            for(i in 1:length(nodes)){
                   adj<-getAdjKOMID(nodes[i],top)
                   if(length(adj)>0){
                        for(j in 1:length(adj)){    
                              if(!isAdjacent(g1,nodes[i],adj[j])){                
                              g2<-addEdge(nodes[i],adj[j],g1,1)
                              g1<-g2    
                              }                 
                        }
                   }
            }
      }
      return(g1)
}
#top<-xmlRoot(xmlTreeParse("e:/ko/ko00010.xml"))
#g3<-getKOMPathGraph(top)






################################################################################
##get adjacent KOs of a KO
getAdjKOMID<-function(KOID,top){
  reactionKOIndex<-0
  reactionKO<-character()
  if(length(KOID)>0&&(!is.na(KOID))&&length(top)>0&&(!is.na(top))){      
      KOReaction<-character()
      KOReactionIndex<-0
      tmp<-character()
      for(i in 1:length(top)){
            attrs<-xmlAttrs(top[[i]])
            if(!is.na(attrs["name"])){
                  tmp1<-unlist(strsplit(attrs["name"]," "))
                  if(length(tmp1)>0){  
                  for(k in 1:length(tmp1)){
                  if(tmp1[k]==KOID){                        
                        if(!is.na(attrs["reaction"])){
                              tmp<-unlist(strsplit(attrs["reaction"]," "))  
                              if(length(tmp)>0){
                                    for(j in 1:length(tmp)){
                                          KOReactionIndex<-KOReactionIndex+1
                                          KOReaction[KOReactionIndex]<-tmp[j]
                                    }   
                              }    
                        }                  
                  }
                  }
                  }
            }
      }
     
      KOReaction<-intersect(KOReaction,KOReaction) 
      #print(KOReaction[1])
      productIndex<-0
      product<-character()
      if(length(KOReaction)>0){
            #print(KOReaction[1])
            for(j in 1:length(KOReaction)){
                  for(i in 1:length(top)){
                        attrs<-xmlAttrs(top[[i]])
                        if(!is.na(attrs["name"])){  
                               if(attrs["name"]==KOReaction[j]){
                                     reactionType<-attrs["type"]
                                     
                                     topReaction<-top[[i]]
                                     if(length(topReaction)>0){
                                           #print(topReaction)
                                          for(k in 1:length(topReaction)){
                                                
                                               
                                                             productIndex<-productIndex+1
                                                             product[productIndex]<-xmlAttrs(topReaction[[k]])["name"]
                                                         
                                         
                                          }
                                     }
                                           
                               }
                        }
                 }
            }                 
      }
      
      product<-intersect(product,product)
      #print(product)
      substrateReaction<-character()
      substrateReactionIndex<-0
      if(length(product)>0){
            for(j in 1:length(product)){
                  for(i in 1:length(top)){
                        if(xmlName(top[[i]])=="reaction"){
                              reactionType<-xmlAttrs(top[[i]])["type"]
                              for(k in 1:length(top[[i]])){
                                    
                                    
                                          if(xmlAttrs(top[[i]][[k]])["name"]==product[j]){
                                                substrateReactionIndex<-substrateReactionIndex+1
                                                substrateReaction[substrateReactionIndex]<-xmlAttrs(top[[i]])["name"]
                                                break
                                          }
                                    
                              }
                        }           
                  }
            }
      }
      substrateReaction<-intersect(substrateReaction,substrateReaction)
      #print(substrateReaction)

      if(length(substrateReaction)>0){
            for(j in 1:length(substrateReaction)){
                  for(i in 1:length(top)){
                        if(!is.na(xmlAttrs(top[[i]])["reaction"])){
                              if(xmlAttrs(top[[i]])["type"]=="ortholog"){
                              splitReaction<-unlist(strsplit(xmlAttrs(top[[i]])["reaction"]," "))
                              for(k in 1:length(splitReaction)){
                                    if(splitReaction[k]==substrateReaction[j]){
                                          tmp1<-unlist(strsplit(xmlAttrs(top[[i]])["name"]," "))
                                          if(length(tmp1)>0){
                                          for(c in 1:length(tmp1)){
                                          if(tmp1[c]!=KOID){                                                
                                                reactionKOIndex<-reactionKOIndex+1
                                                reactionKO[reactionKOIndex]<-tmp1[c]
                                          } 
                                          }
                                          }
                                    }
                              }
                              }
                        }
                  }
            }
      }
   }
   else{ 
       print("getAdjKOID function error") 
   }    
      return(intersect(reactionKO,reactionKO))
}

#top<-xmlRoot(xmlTreeParse("e:/ko/ko00010.xml"))
#KO2<-getAdjKOMID("ko:K00169",top)


####################################################################
##get KO sub-pathway annotation
getKOAnn<-function(geneList,background=getDefaultBackground(),
   order="pvalue",decreasing=FALSE,graphList){
      if(typeof(geneList)!="character"){
	  print("warning: your geneList must be 'character' vector. Because the type of your current geneList is not correct, it has been conveted arbitrarily using the function as.character().")
	  as.character(geneList)
	  }
      if(!exists("ke2g")) initialize_ke2g()
      #subGraph<-getKcSubGraph(k,graphList)
      graphList<-graphList[sapply(graphList,function(x) length(x)>0)]     
      keggpathid2name<-get("keggpathid2name",envir=ke2g)
      
      annList<-list()
      for(i in 1:length(graphList)){
            ann<-list(pathwayName="not known",annGeneList=character(),annGeneNumber=0,
                      annBgGeneList=character(),annBgNumber=0,geneNumber=0,bgNumber=0,pvalue=1,qvalue=1)
            
                  graphGeneList<-getGeneFromKO(nodes(graphList[[i]]))         
                      
            annotatedGeneList<-intersect(graphGeneList,geneList)
            annotatedBackgroundList<-intersect(graphGeneList,background)
            
            pathwayName<-keggpathid2name[[substring(names(graphList)[i],6,10)]]
            #pathwayName<-keggpathid2name[names(graphList)[i]]
            if(length(pathwayName)!=0)
                ann$pathwayName<-pathwayName
            ann$annGeneList<-annotatedGeneList   
         
            ann$annGeneNumber<-length(annotatedGeneList)
			  ann$annBgGeneList<-annotatedBackgroundList
            ann$annBgNumber<-length(annotatedBackgroundList)

            ann$geneNumber<-length(geneList)
            ann$bgNumber<-length(background)

            ann$pvalue<-1-phyper(ann$annGeneNumber-1,ann$annBgNumber,
                 ann$bgNumber-ann$annBgNumber,ann$geneNumber)
            
            annList[[i]]<-ann
      } 
      qvalueList<-fdrtool(sapply(annList,function(x) x$pvalue), 
                   statistic="pvalue",plot=FALSE,verbose=FALSE)$qval
      for(i in 1:length(annList)){
            annList[[i]]$qvalue<-qvalueList[i]
      }
      names(annList)<-names(graphList)
      annList<-annList[sapply(annList,function(x) x$annGeneNumber>0)]
      annList<-annList[order(sapply(annList,function(x) x[[order]]),decreasing=decreasing)]

      return(annList)
}
