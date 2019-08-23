#' readGenes
#'
#' This function reads gene data from an excel file.
#' Por ahora tiene que tener las columnas: Well, Raton, Tratamiento, N0
#' @param excel Path of the excel file to read.
#' @return A dataframe with data read from the excel file.
#' @importFrom magrittr %>%
#' @export
readGenes <- function(excel){
  # Lee todos los genes de un archivo de excel, un gen por hoja.
  # Nombre de hoja = Nombre gen
  Genes <- readxl::excel_sheets(excel)
  names(Genes) <- Genes
  Genes %>%
    purrr::map_dfr(function(x) readxl::read_excel(path = excel, sheet = x,
                                   col_types = c("text",
                                                 "text", "text", "numeric")
    ), .id = "Gen")
}

#' filterDuplicates
#'
#' This function takes a dataframe produced by readGenes.
#' Removes duplicates if ratio is greater than a given threshold.
#' @param geneDataFrame Dataframe produced by readGenes.
#' @param threshold Threshold for removing bad duplicates.
#' @return geneDataFrame without bad duplicates.
#' @importFrom magrittr %>%
#' @export
filterDuplicates <- function(geneDataFrame, threshold = 2){
  # Chequea los duplicados de cada raton para cada gen
  # TODO: si el gen que falla es un HK, avisar y bajar todos los genes de ese
  #       raton. Luego sacar los na.omit de las otras funciones
  checkedDF <- dplyr::left_join(
    geneDataFrame,
    geneDataFrame %>%
      dplyr::group_by(Gen, Raton) %>%
      dplyr::summarise(
        duplicateRatio = max(N0) / min(N0),
        duplicateOk = duplicateRatio <= threshold
        )
    )

  if (sum(!checkedDF$duplicateOk)) {
    print("Failed duplicates:")
    print(
      checkedDF %>%
        dplyr::filter(duplicateOk == F)
      )
  }

  checkedDF %>%
    dplyr::filter(duplicateOk == T) %>%
    dplyr::select(-duplicateOk)
}

#' averageDuplicates
#'
#' This function takes a dataframe produced by readGenes.
#' @param geneDataFrame Dataframe produced by readGenes.
#' @return geneDataFrame with duplicates averaged.
#' @importFrom magrittr %>%
#' @export
averageDuplicates <- function(geneDataFrame){
  # calcula el promedio
  geneDataFrame %>%
    dplyr::group_by(Gen, Raton, Tratamiento) %>%
    dplyr::summarize(N0 = mean(N0))
}

#' plotHouseKeepings
#'
#' This function plots average N0 vs Sample for housekeepings diagnostics.
#' @param geneDataFrame Dataframe produced by readGenes with averaged duplicates.
#' @param houseKeepingNames Vector with names of genes to include in the plot.
#' @importFrom magrittr %>%
#' @export
plotHouseKeepings <- function(geneDataFrame, houseKeepingNames){
  # Hace el clasico grafico de housekeepings
  geneDataFrame %>%
    dplyr::filter(Gen %in% houseKeepingNames) %>%
    ggplot2::ggplot(ggplot2::aes(x = Raton, y = N0, color = Gen, group = Gen)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_y_log10() +
    ggplot2::ylab("Log N0")
}

#' relativizeGeneExpresion
#'
#' This function relativizes averaged gene expression to the geometric mean of a set of housekeepings.
#' @param geneDataFrame Dataframe produced by readGenes with averaged duplicates.
#' @param houseKeepingNames Vector with names of housekeepings that will be used for relativization.
#' @return A dataframe with relativized gene expression levels
#' @importFrom magrittr %>%
#' @export
relativizeGeneExpresion <- function(geneDataFrame, houseKeepingNames){
  HK <- geneDataFrame %>%
    dplyr::filter(Gen %in% houseKeepingNames) %>%
    tidyr::spread(Gen, N0) %>%
    na.omit() %>% # Descarta los que para algun HK fueron filtrados
    tidyr::gather(Gen, N0, -Raton, -Tratamiento) %>%
    dplyr::group_by(Raton, Tratamiento) %>%
    dplyr::summarize(HKGeomMean = exp(mean(log(N0))))

  dplyr::left_join(
      geneDataFrame %>%
        dplyr::filter(!Gen %in% houseKeepingNames),
      HK
    ) %>%
    na.omit() %>% # Descarta los genes que tienen HK pero fallan
    dplyr::mutate(relativeExpression = N0 / HKGeomMean) %>%
    dplyr::select(-N0, -HKGeomMean)
}

#' boxplotRelativeGeneExpression
#'
#' This function plots relative expression by group
#' @param geneDataFrame Dataframe produced by readGenes with averaged duplicates.
#' @importFrom magrittr %>%
#' @export
boxplotRelativeGeneExpression <- function(geneDataFrame){
  ggplot2::ggplot(
      geneDataFrame,
      ggplot2::aes(x = Tratamiento, y = relativeExpression)
    ) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_dotplot(binaxis='y', stackdir='center') +
    ggplot2::scale_x_discrete(limits = c("NP", "LP")) +
    ggplot2::facet_wrap(~Gen, ncol = 4, scales = "free_y")
}

#' barplotRelativeGeneExpression
#'
#' This function plots relative expression by group
#' @param geneDataFrame Dataframe produced by readGenes with averaged duplicates.
#' @importFrom magrittr %>%
#' @export
barplotRelativeGeneExpression <- function(geneDataFrame){
  geneDataFrame %>%
    dplyr::group_by(Gen, Tratamiento) %>%
    dplyr::summarize(
      avgRelativeExpression = mean(relativeExpression),
      semRelativeExpression = sd(relativeExpression) / sqrt(dplyr::n())
    ) %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = Tratamiento,
        y = avgRelativeExpression,
        fill = Tratamiento)
      ) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_dotplot(
      data = geneDataFrame,
      binaxis = "y",
      stackdir = "center",
      color = "black",
      ggplot2::aes(x = Tratamiento, y = relativeExpression)
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = avgRelativeExpression - semRelativeExpression,
        ymax = avgRelativeExpression + semRelativeExpression
      ),
      width = 0.2
    ) +
    ggplot2::scale_x_discrete(limits = c("NP", "LP")) +
    ggplot2::facet_wrap(~Gen, ncol = 5, scales = "free_y") +
    ggplot2::theme(legend.position = "none")
}

#' barplotControlRelativizedGeneExpression
#'
#' This function plots relative expression by group relativized to the average
#' value of the control group.
#' @param geneDataFrame Dataframe produced by readGenes with averaged duplicates.
#' @param labels Whether to include sample labels.
#' @importFrom magrittr %>%
#' @export
barplotControlRelativizedGeneExpression <- function(geneDataFrame, labels = F){
  controlMeanExpression <- geneDataFrame %>%
    dplyr::filter(Tratamiento == "NP") %>%
    dplyr::group_by(Gen) %>%
    dplyr::summarize(controlMeanExpression = mean(relativeExpression))

  controlRelativizedGeneDataFrame <- dplyr::left_join(
      geneDataFrame,
      controlMeanExpression,
      by = "Gen"
    ) %>%
    dplyr::mutate(
      controlRelativizedExpression = relativeExpression / controlMeanExpression
    )

  if (labels) {
    controlRelativizedGeneDataFrame %>%
      dplyr::group_by(Gen, Tratamiento) %>%
      dplyr::summarize(
        avgRelativeExpression = mean(controlRelativizedExpression),
        semRelativeExpression = sd(controlRelativizedExpression) / sqrt(dplyr::n())
      ) %>%
      ggplot2::ggplot(
        ggplot2::aes(
          x = Tratamiento,
          y = avgRelativeExpression,
          fill = Tratamiento
        )
      ) +
      ggplot2::geom_bar(stat = "identity", width = 0.5) +
      ggplot2::geom_text(
        ggplot2::aes(
          label = round(avgRelativeExpression, 2)),
        color = "red",
        nudge_x = -0.4
      ) +
      ggplot2::geom_dotplot(
        data = controlRelativizedGeneDataFrame,
        binaxis = "y",
        stackdir = "center",
        color = "black",
        ggplot2::aes(x = Tratamiento, y = controlRelativizedExpression)
      ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = avgRelativeExpression - semRelativeExpression,
                     ymax = avgRelativeExpression + semRelativeExpression),
        width = 0.2
      ) +
      ggplot2::geom_text(
        data = controlRelativizedGeneDataFrame,
        ggplot2::aes(
          x = Tratamiento,
          y = controlRelativizedExpression,
          label = Raton
        ),
        nudge_x = 0.4
      ) +
      ggplot2::scale_x_discrete(limits = c("NP", "LP")) +
      ggplot2::facet_wrap(~Gen, ncol = 4) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::ylab("Relative Expression") +
      ggplot2::xlab("Treatment")
  } else {
    controlRelativizedGeneDataFrame %>%
      dplyr::group_by(Gen, Tratamiento) %>%
      dplyr::summarize(
        avgRelativeExpression = mean(controlRelativizedExpression),
        semRelativeExpression = sd(controlRelativizedExpression) / sqrt(dplyr::n())
      ) %>%
      ggplot2::ggplot(
        ggplot2::aes(
          x = Tratamiento,
          y = avgRelativeExpression,
          fill = Tratamiento
        )
      ) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::geom_dotplot(
        data = controlRelativizedGeneDataFrame,
        binaxis = "y",
        stackdir = "center",
        color = "black",
        ggplot2::aes(x = Tratamiento, y = controlRelativizedExpression)
      ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(
          ymin = avgRelativeExpression - semRelativeExpression,
          ymax = avgRelativeExpression + semRelativeExpression),
        width = 0.2
      ) +
      ggplot2::scale_x_discrete(limits = c("NP", "LP")) +
      ggplot2::facet_wrap(~Gen, ncol = 4) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::ylab("Relative Expression") +
      ggplot2::xlab("Treatment")
  }

}


#' testGene
#'
#' This function test for significant differences in expression between
#' two groups for a single gene within a relativeExpression dataframe.
#' @param relativeExpDF Dataframe with gene relative expression information.
#' @param groupCol Grouping column to be used for difference testing.
#' @param geneCol Column containing gene names.
#' @param expressionCol Column containing relative expression.
#' @param gene Gene to test.
#' @param alternative Alternative hypothesis for the test
#' ("two.sided", "greater" or "less").
#' @importFrom magrittr %>%
#' @export
testGene <- function(relativeExpDF, groupCol, geneCol, expressionCol, gene, alternative = "two.sided") {

  groups <- unique(
    relativeExpDF %>%
      dplyr::filter({{geneCol}} == gene) %>%
      dplyr::select({{groupCol}}) %>%
      dplyr::pull()
  )


  if (length(groups) > 2) {
    return("more than two groups")
  } else{
    n1 <- relativeExpDF %>%
      dplyr::filter({{groupCol}} == groups[1], {{geneCol}} == gene) %>%
      dplyr::count() %>%
      dplyr::pull()

    n2 <- relativeExpDF %>%
      dplyr::filter({{groupCol}} == groups[2], {{geneCol}} == gene) %>%
      dplyr::count() %>%
      dplyr::pull()

    s1 <- relativeExpDF %>%
      dplyr::filter({{groupCol}} == groups[1], {{geneCol}} == gene) %>%
      dplyr::select({{expressionCol}}) %>%
      dplyr::pull() %>%
      stats::shapiro.test()

    s2 <- relativeExpDF %>%
      dplyr::filter({{groupCol}} == groups[2], {{geneCol}} == gene) %>%
      dplyr::select({{expressionCol}}) %>%
      dplyr::pull() %>%
      stats::shapiro.test()


    if (s1$p.value >= 0.5 & s2$p.value >= 0.5){
      # Data meets normality
      var <- stats::var.test(
        relativeExpDF %>%
          dplyr::filter({{groupCol}} == groups[1], {{geneCol}} == gene) %>%
          dplyr::select({{expressionCol}}) %>%
          dplyr::pull(),
        relativeExpDF %>%
          dplyr::filter({{groupCol}} == groups[2], {{geneCol}} == gene) %>%
          dplyr::select({{expressionCol}}) %>%
          dplyr::pull()
      )

      test <- stats::t.test(
        relativeExpDF %>%
          dplyr::filter({{groupCol}} == groups[1], {{geneCol}} == gene) %>%
          dplyr::select({{expressionCol}}) %>%
          dplyr::pull(),
        relativeExpDF %>%
          dplyr::filter({{groupCol}} == groups[2], {{geneCol}} == gene) %>%
          dplyr::select({{expressionCol}}) %>%
          dplyr::pull(),
        var.equal = (var$p.value>=0.5),
        alternative = alternative
        )

    } else {
      # Data doesn't meet normality
      var <- c(method = "", p.value = NA)

      test <- stats::wilcox.test(
        relativeExpDF %>%
          dplyr::filter({{groupCol}} == groups[1], {{geneCol}} == gene) %>%
          dplyr::select({{expressionCol}}) %>%
          dplyr::pull(),
        relativeExpDF %>%
          dplyr::filter({{groupCol}} == groups[2], {{geneCol}} == gene) %>%
          dplyr::select({{expressionCol}}) %>%
          dplyr::pull(),
        alternative = alternative
      )
    }

    return(
      c(
        g1 = groups[1],
        n1 = n1,
        g2 = groups[2],
        n2 = n2,
        norMethod = s1$method,
        nor1 = s1$p.value,
        nor2 = s2$p.value,
        varMethod = var$method,
        var = var$p.value,
        testMethod = test$method,
        test = test$p.value)
    )
  }
}
