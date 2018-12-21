#Global settings
MAX_CLUSTERS = 10
AXIS_TITLE_SIZE = 19
AXIS_TEXT_SIZE = 8
PCH_TEXT_SIZE = 2
AXIS_FONT = "italic"
COLOR_SAMPLES_DEF = "#000099"

# Creates a circle
circleFun = function(center = c(0, 0), diameter = 2, npoints = 100) {

  r = diameter/2
  tt = seq(0, 2 * pi, length.out = npoints)
  xx = center[1] + r * cos(tt)
  yy = center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

#' Print the variance of a component
#'
#' Prints the percent of explained variance for a component of a block (by default, the superblock or the last one) analysed by R/SGCCA
#'
#' @param rgcca A list giving the results of a R/SGCCA
#' @param n An integer giving the index of the analysis component
#' @param i An integer giving the index of a list of blocks
#' @seealso \code{\link[RGCCA]{rgcca}}, \code{\link[RGCCA]{sgcca}}
#' @examples
#' AVE = list(c(0.6, 0.5), c(0.7, 0.45))
#' rgcca.res = list(AVE = list(AVE_X = AVE))
#' # For the superblock (or the last block)
#' printAxis(rgcca.res, 1)
#' # "Axis 1 (70%)"
#' # For the first block
#' printAxis(rgcca.res, 2, 1)
#' # "Axis 2 (50%)"
#' @export printAxis
printAxis = function (rgcca, n, i = NULL){

  # by default, take the last block
  if ( is.null(i) )
    i = length(rgcca$AVE$AVE_X)

  paste("Axis ", n, " (", round(rgcca$AVE$AVE_X[[i]][n] * 100 , 1),"%)", sep="")
}

#' Default font for plots
theme_perso = function() {

  theme(
    legend.text = element_text(size = 13),
    legend.title = element_text(face="bold.italic", size=16),
    plot.title = element_text(size = 25, face = "bold", hjust=0.5, margin = margin(0,0,20,0))
  )
}

#' Plot of samples space
#'
#' Plots samples on two components of a RGCCA
#'
#' @param rgcca A list giving the results of a R/SGCCA
#' @param resp A vector of characters corresponding either to a qualitative variable with levels or a continuous variable
#' @param comp_x An integer giving the index of the analysis component used for the x-axis
#' @param comp_y An integer giving the index of the analysis component used for the y-axis
#' @param i_block An integer giving the index of a list of blocks
#' @examples
#' coord = lapply(1:3, function(x) matrix(runif(15 * 2, min = -1), 15, 2))
#' AVE_X = lapply(1:3, function(x) runif(2))
#' rgcca.res = list(Y = coord, AVE = list(AVE_X = AVE_X))
#' # Using a superblock
#' plotSamplesSpace(rgcca.res, rep(LETTERS[1:3], each = 5))
#' # Using the first block
#' plotSamplesSpace(rgcca.res, runif(15, min=-15, max = 15), 1, 2, 1)
#' @export plotSamplesSpace
plotSamplesSpace = function (rgcca, resp, comp_x = 1, comp_y = 2, i_block = NULL){
  # resp : color the points with a vector

  if ( is.null(i_block) )
    i_block = length(rgcca$Y)

  df2 = data.frame(rgcca$Y[[i_block]])

  # if the resp is numeric
  if ( ! is.null(resp) ){
    if( ! unique(isCharacter(as.vector(resp)))){
      # add some transparency
      p = ggplot(df2, aes(df2[,comp_x], df2[,comp_y], alpha = (resp - min(resp)) / max(resp - min(resp)))) +
        # get a color scale by quantile
            scale_alpha_continuous(
            name = "resp",
            breaks = seq(0, 1, .25),
            labels = round(quantile(resp))
        ) + geom_text(color = COLOR_SAMPLES_DEF, aes(label = rownames(df2)), size = PCH_TEXT_SIZE)
        #+ geom_text_repel(color=COLOR_SAMPLES_DEF, aes(label= rownames(df2)), size = PCH_TEXT_SIZE, force=2)
    }else
      p = NULL
  }else
    p = NULL

  p = plotSpace(rgcca, df2, "Samples", resp, "resp", comp_x, comp_y, i_block, p)

  # remove legend if missing
  if (is.null(resp)){
    p + theme(legend.position = "none")
  }else
    p
}

#' Get the blocs of each variables
#'
#' Get a vector of block name for each corresponding variable. The last block is considered as the superblock and ignored.
#'
#' @param rgcca A list giving the results of a R/SGCCA
#' @return A vector of character giving block name for each corresponding variable.
#' @seealso \code{\link[RGCCA]{rgcca}}, \code{\link[RGCCA]{sgcca}}
#' @examples
#' rgcca.res = list(a = rep(NA, 4))
#' names(rgcca.res$a) = LETTERS[1:4]
#' getBlocsVariables(rgcca.res)
#' # a, b, c
#' @export getBlocsVariables
getBlocsVariables = function(rgcca){

  rep( names(rgcca$a)[-length(rgcca$a)],
       sapply(rgcca$a[1:(length(rgcca$a)-1)], NROW))
}

#' Plot of variables space
#'
#' Correlation circle highlighting the contribution of each variables to the construction of the RGCCA components
#' @param rgcca A list giving the results of a R/SGCCA
#' @param blocks A list of matrix
#' @param comp_x An integer giving the index of the analysis component used for the x-axis
#' @param comp_y An integer giving the index of the analysis component used for the y-axis
#' @param superblock A boolean giving the presence (TRUE) / absence (FALSE) of a superblock
#' @param i_block An integer giving the index of a list of blocks
#' @examples
#' setMatrix = function(nrow, ncol, iter = 3) lapply(1:iter, function(x) matrix(runif(nrow * ncol), nrow, ncol))
#' blocks = setMatrix(10, 5)
#' blocks[[4]] = Reduce(cbind, blocks)
#' for (i in 1:4)
#'     colnames(blocks[[i]]) = paste( LETTERS[i], as.character(1:NCOL(blocks[[i]])), sep="" )
#' coord = setMatrix(10, 2, 4)
#' a = setMatrix(5, 2)
#' a[[4]] = matrix(runif(15 * 2), 15, 2)
#' AVE_X = lapply(1:4, function(x) runif(2))
#' rgcca.res = list(Y = coord, a = a, AVE = list(AVE_X = AVE_X))
#' names(rgcca.res$a) = LETTERS[1:4]
#' # Using a superblock
#' plotVariablesSpace(rgcca.res, blocks, 1, 2, TRUE)
#' # Using the first block
#' plotVariablesSpace(rgcca.res, blocks, 1, 2, FALSE, 1)
#' @export plotVariablesSpace
plotVariablesSpace = function(rgcca, blocks, comp_x = 1, comp_y = 2, superblock = TRUE, i_block = NULL){

  x = y = NULL

  if ( is.null(i_block) )
    i_block = length(blocks)

  df =  data.frame(
    #correlation matrix within a block for each variables and each component selected
    sapply ( c(comp_x, comp_y), function(x) cor( blocks[[i_block]], rgcca$Y[[i_block]][, x] ) ) ,
    row.names = colnames(blocks[[i_block]])
  )

  # if superblock is selected, color by blocks
  if ( superblock & ( i_block == length(blocks)) )
    color = getBlocsVariables(rgcca)
  else
    color = rep(1, NROW(df))

  df = data.frame(df, color)

  p = plotSpace(rgcca, df, "Variables", color, "Blocks", comp_x, comp_y, i_block) +
    geom_path(aes(x, y), data = circleFun(), col = "grey", size = 1) +
    geom_path(aes(x, y), data = circleFun()/2, col = "grey", size = 1, lty = 2)

  # remove legend if not on superblock
  if ( !superblock || !( i_block == length(blocks) ) )
    p + theme(legend.position = "none")
  else
    p

}

#' Plot of components space
#'
#' Plots RGCCA components in a bi-dimensional space
#'
#' @param rgcca A list giving the results of a R/SGCCA
#' @param df A dataframe
#' @param title A character with the name of the space (either "Variables" or "Samples")
#' @param group A vector of character with levels used to color the points
#' @param name_group A character giving the type of groups (either "Blocs"  or "Response")
#' @param comp_x An integer giving the index of the analysis component used for the x-axis
#' @param comp_y An integer giving the index of the analysis component used for the y-axis
#' @param i_block An integer giving the index of a list of blocks
#' @param p A ggplot object
#' @examples
#' df = as.data.frame(matrix(runif(20*2, min = -1), 20, 2))
#' AVE =  lapply(1:4, function(x) runif(2))
#' rgcca.res = list(AVE = list(AVE_X = AVE))
#' plotSpace(rgcca.res, df, "Samples", rep(c("a","b"), each=10), "Response")
#' @export plotSpace
plotSpace = function (rgcca, df, title, group, name_group, comp_x = 1, comp_y = 2, i_block = 1, p = NULL){

  #if (comp_x > NB_COMP) comp_x = 1
  #if (comp_y > NB_COMP) comp_y = 2

  if (is.null(p)){
    if (name_group == "Blocks"){
      # For variablesPlot
      x = 1; y = 2
    }else{
      x = comp_x; y = comp_y
    }
    p = ggplot(df, aes(df[,x], df[,y], colour = group)) +
      geom_text(aes(label = rownames(df)), size = PCH_TEXT_SIZE)
      #geom_text_repel(aes(label= rownames(df)), size = PCH_TEXT_SIZE, force=2)
  }

  p + theme_classic() +
    geom_vline(xintercept = 0, col = "grey", linetype = "dashed", size = 1) +
    geom_hline(yintercept = 0, col = "grey", linetype = "dashed", size = 1) +
    labs ( title = paste(title, "space"),
           x = printAxis(rgcca, comp_x, i_block),
           y = printAxis(rgcca, comp_y, i_block),
           color = name_group) +
    scale_y_continuous(breaks = NULL) +
    scale_x_continuous(breaks = NULL) +
    theme_perso() +
    theme(
      axis.text = element_blank(),
      axis.title.y = element_text(face = AXIS_FONT, margin = margin(0,20,0,0), size = AXIS_TITLE_SIZE),
      axis.title.x = element_text(face = AXIS_FONT, margin = margin(20,0,0,0), size =AXIS_TITLE_SIZE)
    )
  #+ stat_ellipse()
  #TODO: if NB_VAR > X
}

#' Histogram of a fingerprint
#'
#' Histogram of the higher outer weight vectors for a component of a block (by default, the superblock or the last one) analysed by R/SGCCA
#'
#' @param rgcca A list giving the results of a R/SGCCA
#' @param comp An integer giving the index of the analysis component
#' @param n_mark An integer giving the number of best fingerprint to select
#' @param superblock A boolean giving the presence (TRUE) / absence (FALSE) of a superblock
#' @param i_block An integer giving the index of a list of blocks
#' @seealso \code{\link[RGCCA]{rgcca}}, \code{\link[RGCCA]{sgcca}}
#' @examples
#' weights = lapply(1:3, function(x) matrix(runif(10*2), 10, 2))
#' weights[[4]] = Reduce(rbind, weights)
#' rgcca.res = list(a = weights)
#' names(rgcca.res$a) = LETTERS[1:4]
#' # With the 1rst component of the superblock
#' plotFingerprint(rgcca.res, 1, TRUE)
#' # With the 2nd component of the 1rst block by selecting the 5 higher weights
#' plotFingerprint(rgcca.res, 2, FALSE, 10, 1)
#' @export plotFingerprint
plotFingerprint = function(rgcca, comp = 1, superblock = TRUE, n_mark = 100, i_block = NULL){

  color = NULL

  # if no specific block is selected, by default, the superblock is selected (or the last one)
  if ( is.null(i_block) )
    i_block = length(rgcca$a)

  # select the weights
  df = data.frame(rgcca$a[[i_block]])

  # Get a qualitative variable with which block is associated with each variables
  if (  superblock & ( i_block == length(rgcca$a) ) )
    df = data.frame(df, color = getBlocsVariables(rgcca) )

  # sort in decreasing order
  df = data.frame(df[order(abs(df[,comp]), decreasing = TRUE),], order = nrow(df):1)

  # if the superblock is selected, color the text of the y-axis according to their belonging to each blocks
  #TODO: change this with a booleean with/without superblock
  if (  superblock & ( i_block == length(rgcca$a) ) ){
    color2 = df$color; levels(color2) = hue_pal()(length(rgcca$a)-1)
  }else{
    color2 = "black"
  }

  # max threshold for n
  if(NROW(df) >= n_mark) df = df[1:n_mark,]

  if (  superblock & i_block == length(rgcca$a) ){
    p = ggplot(df, aes(order, df[,comp], fill = color))
  }else{
    p = ggplot(df, aes(order, df[,comp]))
  }

    p = plotHistogram(p, df, "Variable weights", as.character(color2))
    p + labs (subtitle = printAxis(rgcca, comp, i_block))
}

#' Histogram of Average Variance Explained
#'
#' Histogram of the model quality (base on Average Variance Explained) for each blocks and sort in decreasing order
#'
#' @param rgcca A list giving the results of a R/SGCCA
#' @param comp An integer giving the index of the analysis component
#' @seealso \code{\link[RGCCA]{rgcca}}, \code{\link[RGCCA]{sgcca}}
#' @examples
#' random_val = function() lapply(1:4, function(x) runif(1))
#' rgcca.res = list(AVE = list(AVE_X = random_val()), a = random_val())
#' names(rgcca.res$a) = LETTERS[1:4]
#' plotAVE(rgcca.res, 1)
#' @export plotAVE
plotAVE = function(rgcca, comp = 1){

  df = as.matrix ( sapply(rgcca$AVE$AVE_X, function(x) x[comp]) )

  row.names(df) = names(rgcca$a)

  # order by decreasing
  #TODO: catch : Error in data.frame: row names contain missing values : the length of the header is not the same of the row number
  df = data.frame(df[order(abs(df), decreasing = TRUE),], order = nrow(df):1)

  p = ggplot(df, aes(order, df[,1]))
  plotHistogram(p, df, "Average Variance Explained")
}


#' Histogram settings
#'
#' Default font for a vertical barplot.
#'
#' @param p A ggplot object.
#' @param df A dataframe with a column named "order"
#' @param title A character string giving a graphic title
#' @param color A vector of character giving the colors for the rows
#' @examples
#' df = data.frame(x = runif(30), order = 30:1)
#' library("ggplot2")
#' p = ggplot(df, aes(order, x))
#' plotHistogram(p, df, "This is my title", "red")
#' # Add colors per levels of a variable
#' df$color = rep(c(1,2,3), each=10)
#' p = ggplot(df, aes(order, x, fill = color))
#' plotHistogram(p, df, "Histogram", as.character(df$color))
#' @export plotHistogram
plotHistogram = function(p, df, title = "", color = "black"){

    p +
    #TODO: if NB_ROW > X, uncomment this
    #geom_hline(yintercept = c(-.5,.5), col="grey", linetype="dotted", size=1) +
    geom_hline(yintercept = 0, col = "grey", size = 1) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_x_continuous(breaks = df$order, labels = rownames(df)) +
    labs(
      title = title,
      x = "", y = "",
      fill = "Cluster") +
    theme_classic() +
    theme_perso() +
    theme(
      axis.text.y = element_text(size = AXIS_TEXT_SIZE, face = AXIS_FONT, color = color),
      axis.text.x = element_text(size = AXIS_TEXT_SIZE, face = AXIS_FONT, color = "darkgrey"),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      plot.subtitle = element_text(hjust = 0.5, size = 16, face = "italic"))
}

plotDiscriminantVariables = function(t, n, cl, d, m){

  if(m > ncol(d))
    m = ncol(d)

  ctr = 100 * getCtrVar(t, n, cl, d)
  max_ctr= apply(ctr, 2, max)
  which_max_ctr= apply(ctr, 2, which.max)
  df = data.frame(max_ctr[order(max_ctr, decreasing = TRUE)], color = as.character(which_max_ctr), order = length(max_ctr):1)[0:(m), ]
  color2 = as.factor(which_max_ctr); levels(color2) = hue_pal()(length(max_ctr))
  p = ggplot(df, aes(order, df[,1], fill = color))

  plotHistogram(p, df, "Main discriminant variables")
}