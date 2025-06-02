predfunction = function(ourtree,x,ptimes,splitrule){
  #splitrule==0 means use Z %*% v < 0 for splitting rule, splitrule==1 means use 
  # Z %*% v <= 0 for splitting rule. 
  z = c(1,x)
  node = ourtree
  v = node$optv
  while (isNotLeaf(node)){
    if (splitrule == 0)
    {
      splittest=(v%*%z<0)
    }
    else {
      splittest = (v%*%z<=0)
    }
    if (splittest){
      node=node$children[[1]]
      v = node$optv
    } else {
      node=node$children[[2]]
      v= node$optv 
    } 
  }
  outnode=node$name
  pred=summary(node$KMest,times=ptimes,extend=TRUE)[[6]]
  return(list(outnode,pred))
}