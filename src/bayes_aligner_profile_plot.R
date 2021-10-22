#' bayes_alingner_marinal_plot - plot 2d plot of alignment probabilities generated bt BayesAligner.
#'
#' @param filename - a marginal probability file *.profile.*
#' @param title - optional title string
#'
#' @return a data frame with columns Query_Position, Data_position, and Probability
#' 
bayes_aligner_profile_plot <- function(filename,
                                       title = NULL) {
  require(tidyverse)
  
  df <- read.table(filename, 
                   skip = 12,
                   col.names = c('Query_Position', 'Data_Position', 'Probability'))
  
  p <- ggplot(df) +
    stat_summary_2d(aes(x = Query_Position, y = Data_Position, z = Probability),
                    bins = c(df[nrow(df), 'Query_Position'], df[nrow(df), 'Data_Position'])) +
    xlab('Query Position') +
    ylab('Data Position') +
    guides(fill = guide_legend(title = 'Probability')) +
    scale_fill_gradientn(colors = terrain.colors(10))

  if(! is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  print(p)
  
  return(df)
}