################################################################################
##get adjacent enzymes of an enzyme
getAdjECID<-function(ECID,top,direction){
  reactionECIndex<-0
  reactionEC<-character()
  if(length(ECID)>0&&(!is.na(ECID))&&length(top)>0&&(!is.na(top))){      
      ECReaction<-character()
      ECReactionIndex<-0
      tmp<-character()
      for(i in 1:length(top)){
            attrs<-xmlAttrs(top[[i]])
            if(!is.na(attrs["name"])){
                  tmp1<-unlist(strsplit(attrs["name"]," "))
                  if(length(tmp1)>0){  
                  for(k in 1:length(tmp1)){
                  if(tmp1[k]==ECID){                        
                        if(!is.na(attrs["reaction"])){
                              tmp<-unlist(strsplit(attrs["reaction"]," "))  
                              if(length(tmp)>0){
                                    for(j in 1:length(tmp)){
                                          ECReactionIndex<-ECReactionIndex+1
                                          ECReaction[ECReactionIndex]<-tmp[j]
                                    }   
                              }    
                        }                  
                  }
                  }
                  }
            }
      }
     
      ECReaction<-intersect(ECReaction,ECReaction) 
      productIndex<-0
      product<-character()
      if(length(ECReaction)>0){
            #print(ECReaction[1])
            for(j in 1:length(ECReaction)){
                  for(i in 1:length(top)){
                        attrs<-xmlAttrs(top[[i]])
                        if(!is.na(attrs["name"])){  
                               if(attrs["name"]==ECReaction[j]){
                                     reactionType<-attrs["type"]
                                     
                                     topReaction<-top[[i]]
                                     if(length(topReaction)>0){
                                           #print(topReaction)
                                          for(k in 1:length(topReaction)){
                                                if(direction=="directed"){
                                                if(reactionType=="reversible"){
                                                      productIndex<-productIndex+1
                                                      product[productIndex]<-xmlAttrs(topReaction[[k]])["name"]
                                                      
                                                }      
                                                if(reactionType=="irreversible"){
                                                      if(xmlName(topReaction[[k]])=="product"){
                                                            productIndex<-productIndex+1
                                                            product[productIndex]<-xmlAttrs(topReaction[[k]])["name"]
                                                            
                                                      }
                                                }
                                                }
                                                else if(direction=="undirected"){
                                                             productIndex<-productIndex+1
                                                             product[productIndex]<-xmlAttrs(topReaction[[k]])["name"]
                                                }            
                                         
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
                                    if(direction=="directed"){
                                    if(reactionType=="reversible"){
                                          if(xmlAttrs(top[[i]][[k]])["name"]==product[j]){
                                                substrateReactionIndex<-substrateReactionIndex+1
                                                substrateReaction[substrateReactionIndex]<-xmlAttrs(top[[i]])["name"]
                                                break
                                          }
                                    }
                                    if(reactionType=="irreversible"){
                                          if(xmlName(top[[i]][[k]])=="substrate"){
                                                if(xmlAttrs(top[[i]][[k]])["name"]==product[j]){
                                                      substrateReactionIndex<-substrateReactionIndex+1
                                                      substrateReaction[substrateReactionIndex]<-xmlAttrs(top[[i]])["name"]
                                                      break
                                                }
                                          }
                                    }
                                    }
                                    else if(direction=="undirected"){
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
      }
      substrateReaction<-intersect(substrateReaction,substrateReaction)
      #print(substrateReaction)

      if(length(substrateReaction)>0){
            for(j in 1:length(substrateReaction)){
                  for(i in 1:length(top)){
                        if(!is.na(xmlAttrs(top[[i]])["reaction"])){
                              if(xmlAttrs(top[[i]])["type"]=="enzyme"){
                              splitReaction<-unlist(strsplit(xmlAttrs(top[[i]])["reaction"]," "))
                              for(k in 1:length(splitReaction)){
                                    if(splitReaction[k]==substrateReaction[j]){
                                          tmp1<-unlist(strsplit(xmlAttrs(top[[i]])["name"]," "))
                                          if(length(tmp1)>0){
                                          for(c in 1:length(tmp1)){
                                          if(tmp1[c]!=ECID){                                                
                                                reactionECIndex<-reactionECIndex+1
                                                reactionEC[reactionECIndex]<-tmp1[c]
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
       print("getAdjECID function error") 
   }    
      return(intersect(reactionEC,reactionEC))
}
##EC2<-getAdjECID()


#########################################################
##
getNodes<-function(top){
      ECList<-character()
      ECListIndex<-0
      for(i in 1:length(top)){
            if(xmlName(top[[i]])=="entry"){
                  if(!is.na(xmlAttrs(top[[i]])["reaction"])){#should revise
                        if(!is.na(xmlAttrs(top[[i]])["name"])){
                              if(xmlAttrs(top[[i]])["type"]=="enzyme"){
                              tmp1<-unlist(strsplit(xmlAttrs(top[[i]])["name"]," "))
                              if(length(tmp1)>0){
                              for(k in 1:length(tmp1)){
                              ECListIndex<-ECListIndex+1
                              ECList[ECListIndex]<-tmp1[k]
                              }
                              }
                              }
                        }                       
                  }
            }
      }
      ECList<-intersect(ECList,ECList)
      return(ECList)
}
###############################################################
##get pathway graph
getPathGraph<-function(path,direction){
      top<-xmlRoot(xmlTreeParse(path))      
      nodes<-getNodes(top)
      g1<-list()
      if(length(nodes)>2){
            g1<-new("graphNEL",nodes,edgemode=direction)
            for(i in 1:length(nodes)){
                   adj<-getAdjECID(nodes[i],top,direction)
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

#########################################################
##get all pathway graph
getAllPathGraph<-function(fileList,path,direction,verbose=TRUE){
      gList<-list()
      if(verbose==TRUE)
         print("Note that the programming is time consumming!!!!!!!!!!!!!!")
         print(paste("begin to build ",direction," graph!!!!!!!!!!!!!!!!!!!!!!!!!"))
      for(i in 1:length(fileList)){
         if(verbose==TRUE)
            print(paste("building No.",i," graph. it is constructed from file ",fileList[i]))
         gList[[i]]<-getPathGraph(paste(path,fileList[i],sep="/"),direction)  
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