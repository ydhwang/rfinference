

#' finding the previous conditions in the tree
#'
#' finding the "parents"
#' using this function recursively, we can find the family ancestry.
#'
#' @param tree a single tree object obtained from using \code{"randomForest"} and \code{"getTree"}
#' @param i index of the tree
#'
#' @return a list containing \code{"cond"} and \code{"id"})
#'
#' @examples
#' # prevCond(tree, i)
#'
#' @export

prevCond <- function(tree, i) {
  id <- ifelse(i %in% tree$right_daughter,
               which(tree$right_daughter == i),
               which(tree$left_daughter == i))
  cond <- ifelse(i %in% tree$right_daughter,
                 paste(tree$split_var[id], ">", tree$split_point[id]),
                 paste(tree$split_var[id], "<", tree$split_point[id]))
  list(cond = cond, id = id)
}


#' id_root
#'
#' finding the "parents"; recovering the sequence of node ids from the given leaf id.
#' using prevCond function recursively, we can find the family ancestry
#'
#' @param tree a single tree object obtained from using \code{"randomForest"} and \code{"getTree"}
#' @param id index of the tree
#'
#' @return a vector containing the entire "ancestry"
#'
#' @examples
#' # id_root(tree, id)
#'
#' @export

id_root <- function(tree, id){
  prevConds <- prevCond(tree, id)
  id_star <- prevConds$id[1]
  if (id_star == 1){
    return(id_star)
  }else{
    return( c(id_star, id_root(tree, id_star)) )
  }
}


#' cond_root
#'
#' finding the "parents"
#' using this function recursively, we can find the family ancestry.
#'
#'
#' @param tree a single tree object obtained from using \code{"randomForest"} and \code{"getTree"}
#' @param id index of the tree
#'
#' @return a statement defining the leaf (associated with "id"), by collecting the entire "ancestry"
#'
#' @examples
#' # cond_root(tree, id)
#'
#' @export

cond_root <- function(tree, id){
  prevConds <- prevCond(tree, id)
  id_star <- prevConds$id
  cond_star <- prevConds$cond[1]
  if (id_star == 1){
    return(cond_star)
  }else{
    return(paste(cond_star, " & ", cond_root(tree, id_star)))
  }
}

#' getConds
#'
#'
#' give the partition statements for each leaf; a wrapper using cond_root and other functions
#' included for backward compatibility
#'
#' using this function recursively, we can find the family ancestry.
#'
#'
#' @param tree a single tree object obtained from using \code{"randomForest"} and \code{"getTree"}
#'
#' @return a statement defining the leaf (associated with "id"), by collecting the entire "ancestry"
#'
#' @examples
#' # getConds(tree)
#'
#' @export

getConds <- function(tree){
  # backward compatibility
  id.leafs <- id.leafs <- which(tree$status == -1)
  tibble::tibble(condition = sapply(id.leafs, cond_root, tree = tree_b, simplify = TRUE), prediction = tree_b[tree_b$status == -1, "prediction"])
}

#' borehole
#'
#' a test function
#'
#' @param xx a design matrix of 8 dimension
#'
#' @return a vector containing evaluated values of borehole function.
#'
#' @examples
#' # borehole(xx)
#'
#' @export


borehole <- function(xx)
{
  ##########################################################################
  #
  # BOREHOLE FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it
  # and/or modify it under the terms of the GNU General Public License as
  # published by the Free Software Foundation; version 2.0 of the License.
  # Accordingly, this program is distributed in the hope that it will be
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # OUTPUT AND INPUT:
  #
  # y  = water flow rate
  # xx = c(rw, r, Tu, Hu, Tl, Hl, L, Kw)
  #
  ##########################################################################
  #   [0.05, 0.15]	radius of borehole (m)
  # r ∈ [100, 50000]	radius of influence (m)
  # Tu ∈ [63070, 115600]   	transmissivity of upper aquifer (m2/yr)
  # Hu ∈ [990, 1110]	potentiometric head of upper aquifer (m)
  # Tl ∈ [63.1, 116]	transmissivity of lower aquifer (m2/yr)
  # Hl ∈ [700, 820]	potentiometric head of lower aquifer (m)
  # L ∈ [1120, 1680]	length of borehole (m)
  # Kw ∈ [9855, 12045]	hydraulic conductivity of borehole (m/yr)
  rw <- xx[1] * 0.1 + 0.05
  r  <- xx[2] * 49900 + 100
  Tu <- xx[3] * 52530 + 63070
  Hu <- xx[4] * 120 + 990
  Tl <- xx[5] * 52.9 + 63.1
  Hl <- xx[6] * 120 + 700
  L  <- xx[7] * 560 + 1120
  Kw <- xx[8] * 2190 + 9855

  frac1 <- 2 * pi * Tu * (Hu-Hl)

  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r/rw) * (1+frac2a+frac2b)

  y <- frac1 / frac2
  return(y)
}
