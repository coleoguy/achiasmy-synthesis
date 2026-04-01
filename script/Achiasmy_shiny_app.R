library(shiny)
library(ggplot2)
library(plotly)

# =====================================
#  VECTORIZED DATA PREP
# =====================================
achiasmy <- read.csv("Achiasmy_full_data.csv", stringsAsFactors = FALSE)

cn <- names(achiasmy)
achiasmy <- achiasmy[apply(achiasmy, 1, function(r) sum(r %in% cn, na.rm = TRUE) < 3), ]

# Helper: strictly clean unique values
uniq_clean <- function(x) {
  if (is.null(x)) return(character(0))
  x <- x[!is.na(x) & x != ""]
  sort(unique(x))
}

# Pretty display names (convert periods to spaces in UI/legend only)
pretty_name <- function(x) gsub("\\.+", " ", x)

# 1. Vectorized Cleaning
char_cols <- sapply(achiasmy, is.character)
achiasmy[char_cols] <- lapply(achiasmy[char_cols], function(x) {
  x <- trimws(x)
  is.na(x) <- x == ""
  x
})

# 2. Vectorized Type Conversion
if ("Diploid.Number" %in% names(achiasmy)) {
  achiasmy$Diploid.Number <- suppressWarnings(as.numeric(achiasmy$Diploid.Number))
}

# 2.1. Normalize "Distance Pairing" Title (robust to exact column naming)
mt_col <- names(achiasmy)[tolower(trimws(names(achiasmy))) %in% c("meiosis type", "meiosis.type")]
if (length(mt_col) == 1) {
  x <- achiasmy[[mt_col]]
  x <- trimws(x)
  x <- gsub("\\s+", " ", x)
  x[tolower(x) == "distance pairing"] <- "Distance Pairing"
  achiasmy[[mt_col]] <- x
}

# 2.2 Normalize Achiasmatic Sex labels (merge "both male and female?" -> "Both")
as_col <- names(achiasmy)[tolower(trimws(names(achiasmy))) %in% c("achiasmatic sex", "achiasmatic.sex")]
if (length(as_col) == 1) {
  x <- achiasmy[[as_col]]
  x <- trimws(x)
  x <- gsub("\\s+", " ", x)
  x[tolower(x) == "both male and female?"] <- "Both"
  achiasmy[[as_col]] <- x
}

# 3. Vectorized String Paste
achiasmy$GenusSpecies <- paste(achiasmy$Genus, achiasmy$Species)
achiasmy$GenusSpecies[is.na(achiasmy$Genus) | is.na(achiasmy$Species)] <- NA

# Constants
TAXA_LEVELS <- c("Kingdom", "Class", "Order", "Family")
AXIS_VARS   <- setdiff(names(achiasmy), c("Citation", "GenusSpecies"))

# =====================================
#  UI (Compact Generation)
# =====================================
ui <- fluidPage(
  tags$head(tags$style("table, th, td { padding: 2px !important; font-size: 12px !important; }")),
  titlePanel("Achiasmy Karyotype Database"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3, h4("Filters"),
      lapply(TAXA_LEVELS, function(L) {
        selectInput(L, paste0(L, ":"), choices = c("All", uniq_clean(achiasmy[[L]])))
      }),
      
      selectizeInput("Genus", "Genus:", choices = NULL, options = list(placeholder = "Search Genus")),
      selectizeInput("Species", "Species:", choices = NULL, options = list(placeholder = "Search Species")),
      hr(), h4("Search"),
      selectizeInput("genus_species", "Genus + species:", choices = NULL, options = list(placeholder = "Type name..."))
    ),
    
    mainPanel(
      tabsetPanel(
        # ---- SUMMARY TAB ----
        tabPanel(
          "Summary",
          div(
            HTML("
              <br><br>
              A karyotype is a highly variable and complex trait that can reveal
              changes in genome organization, uncover phylogenetic history, and
              distinguish cryptic species. This database compiles achiasmatic
              karyotypes from the literature and allows you to explore patterns
              of recombination and chromosome variation across taxa.
              <br><br>
              Use the filters on the left to subset the dataset by Kingdom, class,
              order, family, genus, or species, or search directly by
              genus–species name. Once you have selected your dataset, you can:
              <ul>
                <li>View a customized table of records in the <b>Data table</b> tab,
                including only the columns you select and downloading the results as a CSV.</li>
                <li>Visualize relationships among variables in the <b>Plot</b> tab,
                choosing any combination of x-axis, y-axis, and color-by variables
                (including categorical traits such as meiosis type or achiasmatic sex).</li>
              </ul>
              <br>
              <b>Submitting data:</b> If you are aware of additional records that should
              be added to the database, please email us and we will incorporate the
              missing data:
              <br><b>Email:</b> blackmon@tamu.edu
              <br><br>
              Data taken from the database must not be reproduced in published lists,
              online databases, or other formats, nor redistributed without permission.
              The information in this database is provided solely for personal and academic
              use, and must not be used for the purposes of financial gain.
              <br><br>
              <b>The database should be cited as follows:
              Alfieri, blah blah blah, &amp; H. Blackmon. Achiasmy Synthesis.
              XXXXX. XX:XX. XX–XX.</b>
              <br><br>
              Current version of the database is 0.1 (last updated 5 May 2022).
            "),
            style = "font-size:100%"
          )
        ),
        
        # ---- DATA TABLE TAB ----
        tabPanel(
          "Data table",
          checkboxGroupInput(
            "columns", "Columns:",
            choices  = setNames(names(achiasmy), pretty_name(names(achiasmy))),
            selected = names(achiasmy),
            inline   = TRUE
          ),
          downloadButton("downloadData", "Download CSV"),
          tableOutput("table")
        ),
        
        # ---- PLOT TAB ----
        tabPanel(
          "Plot", br(),
          fluidRow(
            column(
              4,
              selectInput(
                "xvar", "X:",
                choices  = setNames(AXIS_VARS, pretty_name(AXIS_VARS)),
                selected = AXIS_VARS[1]
              )
            ),
            column(
              4,
              selectInput(
                "yvar", "Y:",
                choices  = setNames(AXIS_VARS, pretty_name(AXIS_VARS)),
                selected = AXIS_VARS[2]
              )
            ),
            column(
              4,
              selectInput(
                "colorvar", "Color:",
                choices = c("None" = "None", setNames(AXIS_VARS, pretty_name(AXIS_VARS)))
              )
            )
          ),
          checkboxInput("flip", "Flip Axes", FALSE),
          plotlyOutput("plot", height = "600px")
        )
      )
    )
  )
)

# =====================================
#  SERVER
# =====================================
server <- function(input, output, session) {
  
  # --- 1. Filter Logic (The Mask) ---
  get_base_mask <- reactive({
    mask <- rep(TRUE, nrow(achiasmy))
    
    if (nzchar(input$genus_species)) {
      return(achiasmy$GenusSpecies == input$genus_species)
    }
    
    if (input$Kingdom != "All") mask <- mask & (achiasmy$Kingdom == input$Kingdom)
    if (input$Class   != "All") mask <- mask & (achiasmy$Class   == input$Class)
    if (input$Order   != "All") mask <- mask & (achiasmy$Order   == input$Order)
    if (input$Family  != "All") mask <- mask & (achiasmy$Family  == input$Family)
    
    mask
  })
  
  # --- 2. Cascading Updates (FIXED) ---
  update_taxa_inputs <- function(dat, current_level) {
    if (current_level %in% c("Kingdom", "Class", "Order"))
      updateSelectInput(session, "Family", choices = c("All", uniq_clean(dat$Family)))
    
    if (current_level != "Genus")
      updateSelectizeInput(session, "Genus", choices = c("All", uniq_clean(dat$Genus)), server = TRUE)
    
    updateSelectizeInput(session, "Species", choices = c("All", uniq_clean(dat$Species)), server = TRUE)
  }
  
  observeEvent(input$Kingdom, {
    sub <- achiasmy
    if (input$Kingdom != "All") sub <- sub[sub$Kingdom == input$Kingdom, ]
    updateSelectInput(session, "Class", choices = c("All", uniq_clean(sub$Class)))
    updateSelectInput(session, "Order", choices = c("All", uniq_clean(sub$Order)))
    update_taxa_inputs(sub, "Kingdom")
  })
  
  observeEvent(input$Class, {
    sub <- achiasmy
    if (input$Kingdom != "All") sub <- sub[sub$Kingdom == input$Kingdom, ]
    if (input$Class != "All") sub <- sub[sub$Class == input$Class, ]
    updateSelectInput(session, "Order", choices = c("All", uniq_clean(sub$Order)))
    update_taxa_inputs(sub, "Class")
  })
  
  observeEvent(input$Order, {
    sub <- achiasmy[get_base_mask(), ]
    update_taxa_inputs(sub, "Order")
  })
  
  observeEvent(input$Family, {
    sub <- achiasmy[get_base_mask(), ]
    update_taxa_inputs(sub, "Family")
  })
  
  observeEvent(input$Genus, {
    sub <- achiasmy[get_base_mask(), ]
    if (!is.null(input$Genus) && input$Genus != "All") {
      sub <- sub[sub$Genus == input$Genus, ]
    }
    updateSelectizeInput(session, "Species", choices = c("All", uniq_clean(sub$Species)), server = TRUE)
  })
  
  updateSelectizeInput(session, "genus_species", choices = c("", uniq_clean(achiasmy$GenusSpecies)), server = TRUE)
  
  # --- 3. Final Data Construction (DEBOUNCED) ---
  raw_filtered_data <- reactive({
    mask <- get_base_mask()
    if (!is.null(input$Genus)   && input$Genus   != "All") mask <- mask & (achiasmy$Genus   == input$Genus)
    if (!is.null(input$Species) && input$Species != "All") mask <- mask & (achiasmy$Species == input$Species)
    achiasmy[mask, ]
  })
  
  final_data <- raw_filtered_data |> debounce(500)
  
  # --- 4. Outputs ---
  output$table <- renderTable({
    req(input$columns)
    final_data()[, input$columns, drop = FALSE]
  })
  
  output$downloadData <- downloadHandler(
    filename = function() paste0("achiasmy_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(final_data()[, input$columns, drop = FALSE], file, row.names = FALSE)
  )
  
  # --- 5. Interactive Plot Output ---
  output$plot <- renderPlotly({
    dat <- final_data()
    req(nrow(dat) > 0, input$xvar, input$yvar)
    
    vars <- c(input$xvar, input$yvar)
    is_num <- sapply(dat[vars], function(x) sum(!is.na(suppressWarnings(as.numeric(x)))) > 0)
    
    if (is_num[1]) dat[[vars[1]]] <- as.numeric(dat[[vars[1]]])
    if (is_num[2]) dat[[vars[2]]] <- as.numeric(dat[[vars[2]]])
    
    cols_to_keep <- unique(c(vars, if (input$colorvar != "None") input$colorvar, "GenusSpecies"))
    dat <- na.omit(dat[, cols_to_keep])
    
    p <- ggplot(dat, aes_string(x = vars[1], y = vars[2], text = "GenusSpecies"))
    
    if (all(is_num)) {
      p <- p + geom_point(aes_string(color = if (input$colorvar != "None") input$colorvar), size = 3, alpha = 0.7)
    } else {
      p <- p + geom_jitter(
        aes_string(color = if (input$colorvar != "None") input$colorvar),
        width = 0.2, height = if (is_num[1]) 0.2 else 0, alpha = 0.7
      )
      if (sum(is_num) == 1) p <- p + geom_boxplot(outlier.shape = NA, alpha = 0.4, fill = "grey90")
    }
    
    p <- p + theme_bw() + labs(title = paste(vars[2], "vs", vars[1]))
    
    # Fix legend title (remove periods) when color is mapped
    if (!is.null(input$colorvar) && input$colorvar != "None") {
      p <- p + labs(color = pretty_name(input$colorvar))
    }
    
    if (!is_num[1]) p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    if (input$flip) p <- p + coord_flip()
    
    ggplotly(p)
  })
}

shinyApp(ui, server)
