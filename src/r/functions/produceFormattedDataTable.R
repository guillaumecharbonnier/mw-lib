produceDataTableWithBgCol <- function(df){
    table <- datatable(df,
                       escape = FALSE,
                       rownames=FALSE,
                       options = list(scrollX = TRUE,
                                      columnDefs = list(list(className = 'dt-center', 
                                                             targets = "_all")))
    )

    pal <- brewer.pal(n=12, name='Set3')

    for (factor_col_name in names(Filter(is.factor, table$x$data))){
        factor_levels <- levels(table$x$data[,factor_col_name])

        table <- formatStyle(table,
                             factor_col_name,
                             backgroundColor = styleEqual(levels = factor_levels,
                                                          values = pal[1:length(factor_levels)]))
    }
    return(table)
}

# This function take a data object (either a matrix or a data frame)
# and produce a datatable with my favourite default settings for display.
produceDataTableWithButtons <- function(df){
    table <- datatable(df,
                              extensions = 'Buttons',
                              options = list(dom = 'Bfrtip',
                                             buttons = c('csv', 'excel'),
                                             scrollX = TRUE,
                                             scrollY = TRUE),
                              class = 'display nowrap',
                              escape = FALSE,
                              rownames=FALSE)
    return(table)
}

