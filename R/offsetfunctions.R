offset.tree <- function(tree, data, s, response, exposure){
  off.tree <- predict(tree$tree.b, type = "tree",newdata = data,
                       s = s, exposure = exposure, response = response)

  if(length(off.tree)==1){off.tree <- rep(0,nrow(data))}
  off.tree
}
