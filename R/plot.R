#' Default font for plots
theme_perso <- function() {
    theme(
        legend.text = element_text(size = 13),
        legend.title = element_text(face = "bold.italic", size = 16),
        plot.title = element_text(
            size = 25,
            face = "bold",
            hjust = 0.5,
            margin = margin(0, 0, 20, 0)
        )
    )
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
#' df <- data.frame(x = runif(30), order = 30:1)
#' library("ggplot2")
#' p <- ggplot(df, aes(order, x))
#' plotHistogram(p, df, "This is my title", "red")
#' # Add colors per levels of a variable
#' df$color <- rep(c(1, 2, 3), each = 10)
#' p <- ggplot(df, aes(order, x, fill = color))
#' plotHistogram(p, df, "Histogram", as.character(df$color))
#' @export plotHistogram
plotHistogram <- function(p, df, title = "", color = "black") {
    p +
        # TODO: if NB_ROW > X, uncomment this
        # geom_hline(yintercept = c(-.5,.5), col="grey", linetype="dotted", size=1) +
        geom_hline(yintercept = 0, col = "grey", size = 1) +
        geom_bar(stat = "identity") +
        coord_flip() +
        scale_x_continuous(breaks = df$order, labels = rownames(df)) +
        labs(
            title = title,
            x = "",
            y = "",
            fill = "Cluster"
        ) +
        theme_classic() +
        theme_perso() +
        theme(
            axis.text.y = element_text(
                size = 12,
                face = "italic",
                color = color
            ),
            axis.text.x = element_text(
                size = 12,
                face = "italic",
                color = "darkgrey"
            ),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            plot.subtitle = element_text(
                hjust = 0.5,
                size = 16,
                face = "italic"
            )
        )
}

getDiscriminantVariables <- function(t, n, cl, d, m) {
    if (m > ncol(d)) {
        m <- ncol(d)
    }

    ctr <- 100 * getCtrVar(t, n, cl, d)
    max_ctr <- apply(ctr, 2, sum)
    # which_max_ctr= apply(ctr, 2, which.max)
    # color = as.character(which_max_ctr),
    res <- data.frame(discr_var = max_ctr[order(max_ctr, decreasing = TRUE)], order = length(max_ctr):1)[0:(m), ]
    # color2 = as.factor(which_max_ctr); levels(color2) = hue_pal()(length(max_ctr))
    # fill = color
}

plotDiscriminantVariables <- function(discr) {
    p <- ggplot(discr, aes(order, discr[, 1]))
    plotHistogram(p, discr, "Main discriminant variables")
}
