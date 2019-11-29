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

#' filterOutliers
#'
#' This function removes outliers in two groups for a single gene within a
#' relativeExpression dataframe.
#' @param relativeExpDF Dataframe with gene relative expression information.
#' @param groupCol Grouping column to be used for difference testing.
#' @param geneCol Column containing gene names.
#' @param expressionCol Column containing relative expression.
#' @param gene Gene to test.
#' @param controlGroup Name of the control group.
#' @param treatedGroup Name of the treated group.
#' ("two.sided", "greater" or "less").
#' @importFrom magrittr %>%
#' @export
filterOutliers <- function(relativeExpDF, groupCol, geneCol, expressionCol, gene, controlGroup, treatedGroup) {
  tControlMax <- relativeExpDF %>%
    dplyr::filter({{geneCol}} == gene, {{groupCol}} == controlGroup) %>%
    dplyr::select({{expressionCol}}) %>%
    dplyr::pull() %>%
    outliers::grubbs.test(opposite = F)

  tControlMin <- relativeExpDF %>%
    dplyr::filter({{geneCol}} == gene, {{groupCol}} == controlGroup) %>%
    dplyr::select({{expressionCol}}) %>%
    dplyr::pull() %>%
    outliers::grubbs.test(opposite = T)

  if(tControlMax$p.value < tControlMin$p.value & tControlMax$p.value < 0.05) {
    # Filtro max
    outlierC <- relativeExpDF %>%
      dplyr::filter({{geneCol}} == gene, {{groupCol}} == controlGroup) %>%
      dplyr::arrange(desc({{expressionCol}})) %>%
      head(1)

    max <- relativeExpDF %>%
      dplyr::filter({{geneCol}} == gene, {{groupCol}} == controlGroup) %>%
      dplyr::arrange(desc({{expressionCol}})) %>%
      dplyr::select({{expressionCol}}) %>%
      utils::head(1) %>%
      dplyr::pull()

    relativeExpDF <- relativeExpDF %>%
      dplyr::filter(
        {{geneCol}} != gene |
        {{groupCol}} != controlGroup |
        {{expressionCol}} != max
    )
  } else if(tControlMin$p.value < tControlMax$p.value & tControlMin$p.value < 0.05) {
    # Filtro min
    outlierC <- relativeExpDF %>%
      dplyr::filter({{geneCol}} == gene, {{groupCol}} == controlGroup) %>%
      dplyr::arrange({{expressionCol}}) %>%
      head(1)

    min <- relativeExpDF %>%
      dplyr::filter({{geneCol}} == gene, {{groupCol}} == controlGroup) %>%
      dplyr::arrange({{expressionCol}}) %>%
      dplyr::select({{expressionCol}}) %>%
      utils::head(1) %>%
      dplyr::pull()

    relativeExpDF <- relativeExpDF %>%
      dplyr::filter(
        {{geneCol}} != gene |
        {{groupCol}} != controlGroup |
        {{expressionCol}} != min
      )
  }

  tTreatedMax <- relativeExpDF %>%
    dplyr::filter({{geneCol}} == gene, {{groupCol}} == treatedGroup) %>%
    dplyr::select({{expressionCol}}) %>%
    dplyr::pull() %>%
    outliers::grubbs.test(opposite = F)

  tTreatedMin <- relativeExpDF %>%
    dplyr::filter({{geneCol}} == gene, {{groupCol}} == treatedGroup) %>%
    dplyr::select({{expressionCol}}) %>%
    dplyr::pull() %>%
    outliers::grubbs.test(opposite = T)

  if(tTreatedMax$p.value < tTreatedMin$p.value & tTreatedMax$p.value < 0.05) {
    # Filter max
    outlierT <- relativeExpDF %>%
      dplyr::filter({{geneCol}} == gene, {{groupCol}} == treatedGroup) %>%
      dplyr::arrange(desc({{expressionCol}})) %>%
      head(1)

    max <- relativeExpDF %>%
      dplyr::filter({{geneCol}} == gene, {{groupCol}} == treatedGroup) %>%
      dplyr::arrange(desc({{expressionCol}})) %>%
      dplyr::select({{expressionCol}}) %>%
      utils::head(1) %>%
      dplyr::pull()

    relativeExpDF <- relativeExpDF %>%
      dplyr::filter(
        {{geneCol}} != gene |
        {{groupCol}} != treatedGroup |
        {{expressionCol}} != max
      )
  } else if(tTreatedMin$p.value < tTreatedMax$p.value & tTreatedMin$p.value < 0.05) {
    # Filter min
    outlierT <- relativeExpDF %>%
      dplyr::filter({{geneCol}} == gene, {{groupCol}} == treatedGroup) %>%
      dplyr::arrange({{expressionCol}}) %>%
      head(1)

    min <- relativeExpDF %>%
      dplyr::filter({{geneCol}} == gene, {{groupCol}} == treatedGroup) %>%
      dplyr::arrange({{expressionCol}}) %>%
      dplyr::select({{expressionCol}}) %>%
      utils::head(1) %>%
      dplyr::pull()

    relativeExpDF <- relativeExpDF %>%
      dplyr::filter(
        {{geneCol}} != gene |
        {{groupCol}} != treatedGroup |
        {{expressionCol}} != min
      )
  }

  # print outliers
  if(base::exists("outlierC") | base::exists("outlierT")) {
    print("Outliers:")
    if (base::exists("outlierC")) {
      print(outlierC)
    }
    if (base::exists("outlierT")) {
      print(outlierT)
    }
  }

  # returns filtered DF
  return(relativeExpDF)
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
#' @param filterOutliers Whether to remove outliers or not. Only one outlier per
#' group based on Grubb's test.
#' @param controlGroup Name of the control group.
#' @param treatedGroup Name of the treated group.
#' ("two.sided", "greater" or "less").
#' @importFrom magrittr %>%
#' @export
testGene <- function(relativeExpDF, groupCol, geneCol, expressionCol, gene, alternative = "two.sided", filterOutliers = T, controlGroup, treatedGroup) {

  groups <- unique(
    relativeExpDF %>%
      dplyr::filter({{geneCol}} == gene) %>%
      dplyr::select({{groupCol}}) %>%
      dplyr::pull()
  )


  if (length(groups) > 2) {
    return("more than two groups")
  }

  if (filterOutliers) {
    relativeExpDF <- filterOutliers(
      relativeExpDF = relativeExpDF,
      groupCol = {{groupCol}},
      geneCol = {{geneCol}},
      expressionCol = {{expressionCol}},
      gene = gene,
      controlGroup = controlGroup,
      treatedGroup = treatedGroup
    )
  }

  n1 <- relativeExpDF %>%
    dplyr::filter({{groupCol}} == controlGroup, {{geneCol}} == gene) %>%
    dplyr::count() %>%
    dplyr::pull()

  n2 <- relativeExpDF %>%
    dplyr::filter({{groupCol}} == treatedGroup, {{geneCol}} == gene) %>%
    dplyr::count() %>%
    dplyr::pull()

  s1 <- relativeExpDF %>%
    dplyr::filter({{groupCol}} == controlGroup, {{geneCol}} == gene) %>%
    dplyr::select({{expressionCol}}) %>%
    dplyr::pull() %>%
    stats::shapiro.test()

  s2 <- relativeExpDF %>%
    dplyr::filter({{groupCol}} == treatedGroup, {{geneCol}} == gene) %>%
    dplyr::select({{expressionCol}}) %>%
    dplyr::pull() %>%
    stats::shapiro.test()


  if (s1$p.value >= 0.05 & s2$p.value >= 0.05){
    # Data meets normality
    var <- stats::var.test(
      relativeExpDF %>%
        dplyr::filter({{groupCol}} == controlGroup, {{geneCol}} == gene) %>%
        dplyr::select({{expressionCol}}) %>%
        dplyr::pull(),
      relativeExpDF %>%
        dplyr::filter({{groupCol}} == treatedGroup, {{geneCol}} == gene) %>%
        dplyr::select({{expressionCol}}) %>%
        dplyr::pull()
    )

    test <- stats::t.test(
      relativeExpDF %>%
        dplyr::filter({{groupCol}} == treatedGroup, {{geneCol}} == gene) %>%
        dplyr::select({{expressionCol}}) %>%
        dplyr::pull(),
      relativeExpDF %>%
        dplyr::filter({{groupCol}} == controlGroup, {{geneCol}} == gene) %>%
        dplyr::select({{expressionCol}}) %>%
        dplyr::pull(),
      var.equal = (var$p.value>=0.05),
      alternative = alternative
      )
  } else {

    test <- stats::wilcox.test(
      relativeExpDF %>%
        dplyr::filter({{groupCol}} == treatedGroup, {{geneCol}} == gene) %>%
        dplyr::select({{expressionCol}}) %>%
        dplyr::pull(),
      relativeExpDF %>%
        dplyr::filter({{groupCol}} == controlGroup, {{geneCol}} == gene) %>%
        dplyr::select({{expressionCol}}) %>%
        dplyr::pull(),
      alternative = alternative
    )
    var <- list(p.value = NA, method = "")
  }

  if(substr(test$method,1,1) == " ") {
    test$method = substr(test$method,2,nchar(test$method))
  }

  return(
    list(
      g1 = controlGroup,
      n1 = n1,
      g2 = treatedGroup,
      n2 = n2,
      norMethod = s1$method,
      nor1 = s1$p.value,
      nor2 = s2$p.value,
      varMethod = var$method,
      var = var$p.value,
      testMethod = test$method,
      test = test$p.value,
      alternative = alternative
    )
  )
}

#' plotGene
#'
#' This function plots a bar and dots graph for a single gene
#' comparing two groups using the testGene function.
#' @param relativeExpDF Dataframe with gene relative expression information.
#' @param groupCol Grouping column to be used for difference testing.
#' @param geneCol Column containing gene names.
#' @param expressionCol Column containing relative expression.
#' @param gene Gene to test.
#' @param alternative Alternative hypothesis for the test
#' ("two.sided", "greater" or "less").
#' @param filterOutliers Whether to remove outliers or not. Only one outlier per
#' group based on Grubb's test.
#' @param controlRelativized Display expression relativized to control group
#' (i.e. control group mean expression = 1)
#' @param controlGroup Name of the control group.
#' @param treatedGroup Name of the treated group.
#' @param maxY Upper limit of y axis.
#' @importFrom magrittr %>%
#' @importFrom rlang :=
#' @export
plotGene <- function(relativeExpDF, groupCol, geneCol, expressionCol, gene, alternative = "two.sided", filterOutliers = T, controlRelativized = F, controlGroup, treatedGroup, maxY = 2) {
  groups <- unique(
    relativeExpDF %>%
      dplyr::filter({{geneCol}} == gene) %>%
      dplyr::select({{groupCol}}) %>%
      dplyr::pull()
  )

  if (length(groups) > 2) {
    return("more than two groups")
  }

  relativeExpDF <- relativeExpDF %>%
    dplyr::filter({{geneCol}} == gene) %>%
    dplyr::ungroup()

  if (filterOutliers) {
    relativeExpDF <- filterOutliers(
      relativeExpDF = relativeExpDF,
      groupCol = {{groupCol}},
      geneCol = {{geneCol}},
      expressionCol = {{expressionCol}},
      gene = gene,
      controlGroup = controlGroup,
      treatedGroup = treatedGroup
    )
  }

  testResult <- testGene(
    relativeExpDF = relativeExpDF,
    groupCol = {{groupCol}},
    geneCol = {{geneCol}},
    expressionCol = {{expressionCol}},
    gene = gene,
    alternative = alternative,
    filterOutliers = F,
    controlGroup = controlGroup,
    treatedGroup = treatedGroup
  )

  if (controlRelativized) {
    controlMeanExpression <- relativeExpDF %>%
      dplyr::filter({{groupCol}} == controlGroup) %>%
      dplyr::select({{expressionCol}}) %>%
      dplyr::pull() %>%
      base::mean()

    relativeExpDF <- relativeExpDF %>%
      dplyr::mutate(
        {{expressionCol}} := {{expressionCol}} / controlMeanExpression
      )
  }

  plot <- relativeExpDF %>%
    dplyr::group_by({{groupCol}}) %>%
    dplyr::summarize(
      avgRelativeExpression = mean({{expressionCol}}),
      semRelativeExpression = stats::sd({{expressionCol}}) / base::sqrt(dplyr::n())
    ) %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = {{groupCol}},
        y = avgRelativeExpression
      )
    ) +
    ggplot2::geom_bar(stat = "identity", color = "black", alpha = 0, width = .75) +
    ggplot2::geom_dotplot(
      data = relativeExpDF,
      method = "histodot",
      binaxis = "y",
      stackdir = "center",
      binwidth = 0.05,
      dotsize = 1,
      alpha = 0.75,
      ggplot2::aes(
        x = {{groupCol}},
        y = {{expressionCol}},
        fill = {{groupCol}},
        color = {{groupCol}}
      )
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = avgRelativeExpression - semRelativeExpression,
        ymax = avgRelativeExpression + semRelativeExpression),
      width = 0.2
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expand_scale(add = c(0, .1)),
      limits = c(0, maxY)
    ) +
    ggplot2::scale_x_discrete(
      limits = c(controlGroup, treatedGroup),
      labels = c(
        base::paste0(controlGroup, "\n(n = ", testResult$n1, ")"),
        base::paste0(treatedGroup, "\n(n = ", testResult$n2, ")")
      )
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      axis.line = ggplot2::element_line()
    ) +
    ggplot2::ylab("Relative Expression") +
    ggplot2::xlab("Treatment")

  if(testResult$alternative == "greater"){
    alternative = base::paste(treatedGroup, "greater than", controlGroup)
  } else if(testResult$alternative == "less") {
    alternative = base::paste(treatedGroup, "less than", controlGroup)
  } else {
    alternative = base::paste(treatedGroup, "different than", controlGroup)
  }

  maxExp <- relativeExpDF %>%
    dplyr::select({{expressionCol}}) %>%
    dplyr::pull() %>%
    base::max()

  if (base::signif(testResult$test < 0.001, 2)) {
    plot <- plot +
      ggplot2::ggtitle(
        label = gene,
        subtitle = base::paste0(
          testResult$testMethod,
          ": p < 0.001\nH\u2090 : ",
          alternative
        )
      ) +
      ggplot2::annotate("text",x = 1.5,y = maxExp * 1.1, label="***", size = 8)
  } else if (base::signif(testResult$test < 0.01, 2)) {
    plot <- plot +
      ggplot2::ggtitle(
        label = gene,
        subtitle = base::paste0(
          testResult$testMethod,
          ": p < 0.01\nH\u2090 : ",
          alternative
        )
      ) +
      ggplot2::annotate("text",x = 1.5,y = maxExp * 1.1, label="**", size = 8)
  } else if (base::signif(testResult$test < 0.05, 2)) {
    plot <- plot +
      ggplot2::ggtitle(
        label = gene,
        subtitle = base::paste0(
          testResult$testMethod,
          ": p = ",
          base::signif(testResult$test, 2),
          "\nH\u2090 : ",
          alternative
        )
      ) +
      ggplot2::annotate("text",x = 1.5,y = maxExp * 1.1, label="*", size = 8)
  } else {
    plot <- plot +
      ggplot2::ggtitle(
        label = gene,
        subtitle = base::paste0(
          testResult$testMethod,
          ": p = ",
          base::signif(testResult$test, 2),
          "\nH\u2090 : ",
          alternative
        )
      )
  }
  return(plot)
}
