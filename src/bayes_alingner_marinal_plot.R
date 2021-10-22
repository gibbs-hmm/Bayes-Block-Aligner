#' bayes_alingner_marinal_plot - plot marginal alignment probability generated bt BayesAligner.
#'
#' @param filename - a marginal probability file *.marginal.*
#' @param title - optional title string
#'
#' @return a data frame with columns Position and Probability
#' 
bayes_alingner_marinal_plot <- function(filename,
                                        title = NULL) {
  require(tidyverse)
  
  df <- read.table(filename, 
                   skip = 6,
                   col.names = c('Position', 'Probability'))
  
  p <- ggplot(df) +
    geom_linerange(aes(x = Position, ymin = 0, ymax = Probability)) +
    ylab('Marginal Probability of Alignment') +
    ylim(c(0,1))
  
  if(! is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  print(p)
  
  return(df)
}