#
# App to calculate reference points and compare results of various assumptions
#
################################################################################

# Load necessary libraries
library(shiny); library(tidyr); library(ggplot2); 
library(patchwork); library(shinybusy); library(scales); library(DT)
source("0_ref_pt_helper_functions_231221.R")

# Text for the instructions page
instructions_text <- "<p>
This Shiny app enables users to input age-based data for a stock and investigate 
how time-varying productivity and productivity assumptions affect reference point 
estimates. In particular, this app calculates the equilibrium-based reference points 
φ<sub>0</sub>, SSB<sub>0</sub>, F<sub>MSY</sub>, SSB<sub>MSY</sub>, F<sub>X%SPR</sub>, 
and SSB<sub>X%SPR</sub>.</p>

<p>The data must be submitted in the form of a CSV file, with each row of the file 
corresponding to a specific year-age combination (e.g., a row for year 2000 and age 1). 
To perform the necessary calculations, the app requires the following columns:
<ul>
  <li>Year</li>
  <li>Age</li>
  <li>Weight-at-age (in kg)</li>
  <li>Maturity-at-age</li>
  <li>Selectivity-at-age</li>
</ul>
Natural mortality can be provided as a constant (e.g., M = 0.2) or a column. Missing 
values are not allowed in any of these columns. Additionally, the application requires 
parameters to describe the stock-recruit relationship. You can specify this relationship 
using either a Beverton-Holt or a Ricker model, which involve alpha and beta parameters, 
or as a constant value if no stock-recruit relationship is assumed.</p>

<p>For support and other queries, please contact 
<a href='mailto:Nathan.Hebert@dfo-mpo.gc.ca'>Nathan.Hebert@dfo-mpo.gc.ca</a>.</p>"

# Set font size for plots
font_size <- 16

# Vectors used by the plots
palette <- hue_pal()(5)
labels <- c("All four","Natural mortality", "Maturity-at-age","Selectivity-at-age",
            "Weight-at-age")

#############################HELPER FUNCTIONS#############################

# Function to render an empty/placeholder plot with an optional message
empty_plot <- function(xlab, ylab, message = "Click 'Calculate'")
{
  ggplot() + labs(x = xlab, y = ylab) + theme_classic() + 
    theme(text = element_text(size = font_size),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank()) +
    geom_text(aes(x = 0.5, y = 0.5, label = message), vjust = 0.5, hjust = 0.5, size = 6)
}

# Function to generate one of the time series plots
ts_plot <- function(data, ylab, palette, labels)
{
  ggplot() + geom_path(data = data, size = 1.08, aes(x = yr, y = value, 
                                                               color = variable)) + 
    theme_classic() + scale_y_continuous(ylab, expand = c(0,0), 
                                         limits = c(0, NA)) +
    labs(x="Year") + scale_color_manual(values = palette, labels = labels) +
    labs(col = "Time-varying\ncomponent(s)") + theme(text = element_text(size = font_size))
}

# Function to generate one of the plots comparing time periods
RP_time_period_plot <- function(RP, xlab, ref_points_df, fill = palette[1])
{
  labels <- c("Time Period 3", "Time Period 2","Time Period 1")
  ggplot(ref_points_df, aes_string(y = "factor(`Time Period`, levels = labels)", 
                                   x = paste0("`", RP, "`"))) +
           geom_bar(stat = "identity", position = "dodge", fill = fill) +
           labs(y = "", x = xlab) + theme_classic() + 
           theme(text = element_text(size = font_size)) +
           geom_text(aes(label = format(round(get(RP), 2), nsmall = 2), 
                         group = `Time Period`, x = get(RP) / 2), 
                     position = position_dodge(width = 0.9), size = 6, color = "black")
}

###############################DEFINE UI##################################

SPR_default <- 40 # Start with 40%SPR as default for calcs

ui <- fluidPage(
  tags$head(tags$style(HTML(".irs-grid-text { font-size: 8.5pt; }"))), # Slider font size
  titlePanel("Equilibrium Reference Points Calculator"),
  
  mainPanel(width = 12,
  tabsetPanel(
    tabPanel("Instructions", class = "tab-content",
             add_busy_spinner(spin = "fading-circle"),
             fluidRow(
               HTML(instructions_text)
             )
    ),
    tabPanel("Setup",
             tabsetPanel(
    tabPanel("Data Upload", class = "tab-content", style = "height: 70vh;",
             add_busy_spinner(spin = "fading-circle"),
             fluidRow(
               column(width = 2,
                      wellPanel(
                      fileInput("dataFile", "Upload Data File (CSV):", placeholder = "No file"),
                      textInput("sep", "Seperator:", ","),
                      radioButtons("header", "Header:",
                                   choices = c("True" = "TRUE", 
                                               "False" = "FALSE")), class = "custom-well"),
                      wellPanel(
                        selectInput("years", "Year Column:", ""),
                        selectInput("ages", "Age Column:", ""), class = "custom-well")
                      ),
               column(width = 10,
                      wellPanel(
                        style = "height: 70vh%; width: 98%;",
                        DTOutput("table"), class = "custom-well"
                      ))
             )
    ),
    
    tabPanel("Growth, Maturity, and Selectivity", 
             class = "tab-content", style = "height: 70vh;",
             add_busy_spinner(spin = "fading-circle"),
             fluidRow(
               column(width = 2,
                      wellPanel(
                          selectInput("waa", "Weight-At-Age Column (in kg):", ""), class = "custom-well"),
                      wellPanel(
                          selectInput("mat", "Maturity-At-Age Column:", ""),
                          selectInput("selectivity", "Selectivity-At-Age Column:", ""), 
                          class = "custom-well")),
               column(width = 10,
                      plotOutput("growth_selmat_plot", height = "70vh")),
             absolutePanel(
               top = "70%", left = "90%",
               width = 120, height = 60,
               draggable = TRUE,
               selectInput("year_selmat", "Year to Plot At-Age Maturity and/or Selectivity:", ""),
               style = "transform: translate(-50%, -50%);")
             )
    ),
    
    tabPanel("Natural Mortality", class = "tab-content", style = "height: 70vh;",
             add_busy_spinner(spin = "fading-circle"),
             fluidRow(
               column(width = 2,
                      wellPanel(
                        radioButtons("mType", "M-At-Age Data Type:", 
                                     choices = c("Constant" = "constant", 
                                                 "Column" = "column")),
                        conditionalPanel(
                          condition = 'input.mType == "constant"',  # Show when "Constant" selected
                          sliderInput("mConstant", "Constant Value:", value = 0.2, 
                                      step = 0.05, min = 0, max = 2)
                        ),
                        conditionalPanel(
                          condition = 'input.mType != "constant"',  # Hide when "Constant" selected
                          selectInput("m_column", "Column:", "")
                        ), class = "custom-well")
               ),
               column(width = 10,
                      plotOutput("survivorshipPlot", height = "70vh")
               )
             ),
             absolutePanel(
               top = "25%", left = "90%",
               width = 120, height = 60,
               draggable = TRUE,
               selectInput("year_survivor", "Year to Plot Unfished Survivorship:", ""),
               style = "transform: translate(-50%, -50%);"
             )
    ),
    
    tabPanel("Stock-Recruit Relationship", class = "tab-content",
             add_busy_spinner(spin = "fading-circle"), style = "height: 70vh;",
             fluidRow(
               column(width = 2,
                      wellPanel(
                      radioButtons("srr", "Model:", 
                                   choices = c("Beverton-Holt" = "beverton-holt", 
                                               "Ricker" = "ricker", 
                                               "Constant (No Relationship)" = "constant")), 
                      conditionalPanel(
                        condition = 'input.srr == "constant"',  # Show when "Constant" selected
                        numericInput("recruitment_mean", "Constant Value:", value = 1, 
                                     step = 1, min = 0)
                      ),
                      conditionalPanel(
                        condition = 'input.srr != "constant"',  # Hide when "Constant" selected
                        numericInput("alpha", "Alpha Value:", value = 0.09, 
                                     step = 0.01, min = 0),
                        numericInput("beta", "Beta Value:", value = 0.03, 
                                     step = 0.01, min = 0)), class = "custom-well")),
                column(width = 10,
                      plotOutput("stockRecruitPlot", height = "70vh")),
               absolutePanel(
                 top = "25%", left = "90%",
                 width = 120, height = 120,
                 draggable = TRUE,
                 sliderInput("max_y", "Y-Axis Maximum for Plot", 5, min = 0.0001, max = 10000, step = 0.01),
                 sliderInput("max_x", "X-Axis Maximum for Plot", 400, min = 10, max = 2000, step = 10),
                 style = "transform: translate(-50%, -50%);"
               )
             )
    ))),
    
    tabPanel("Results",
             tabsetPanel(
    tabPanel("Reference Points Over Time", class = "tab-content", style = "height: 70vh;", 
             add_busy_spinner(spin = "fading-circle"),
             fluidRow(
               column(width = 2, 
                      wellPanel(
                      sliderInput("SPR1", "% SPR", value = SPR_default, min = 10, max = 100, 
                                  step = 1), class = "custom-well"),
                      wellPanel(
                      actionButton("calculateBtn", "Calculate"))),
               column(width = 10, 
                      plotOutput("combined_time_series_plot", height = "70vh"))
             )
    ),
    
    tabPanel("Compare Time Periods for Aggregation", class = "tab-content", style = "height: 70vh;",
             add_busy_spinner(spin = "fading-circle"),
             fluidRow(
               column(width = 2,
                      wellPanel(
                      sliderInput("SPR2", "% SPR", value = SPR_default, min = 10, max = 100, 
                                  step = 1), class = "custom-well"),
                      wellPanel(
                      sliderInput("yearRange1", "Time Period 1:", 
                                  min = 1970, max = 2022, value = c(2000, 2010),
                                  step = 1, sep = ""),
                      sliderInput("yearRange2", "Time Period 2:", 
                                  min = 1970, max = 2022, value = c(2000, 2005),
                                  step = 1, sep = ""),
                      sliderInput("yearRange3", "Time Period 3:", 
                                  min = 1970, max = 2022, value = c(2005, 2005),
                                  step = 1, sep = ""), class = "custom-well"),
                      wellPanel(
                      actionButton("calculateBtn2", "Calculate"))),
               column(
                 width = 10,
                 plotOutput("combinedPlot", height = "70vh"))
             )))))
  ),
    
  # Define CSS styles to control margins, etc.
  tags$style(HTML('
    .tab-content {
      margin-left: 10px; margin-top: 15px;
    }
  ')),
  tags$style(HTML(
    ".custom-well { padding-bottom: 2px; }"
    ))
)

############################DEFINE SERVER LOGIC##################################

server <- function(input, output, session) {
  
  #####USER CONTROL CODE####
  
  # Let user load in csv files
  data <- reactive({
    req(input$dataFile)
    read.csv(input$dataFile$datapath, sep = input$sep, header = as.logical(input$header))
  })
  
  # Display placeholder for dataframe initially
  output$table <- renderDT({return(data.frame(Placeholder = character(0)))})
  # Display csv when uploaded
  observeEvent(input$dataFile, {
  output$table <- renderDT({
    req(data())
    datatable(data(), options = list(
      scrollY = "400px", pageLength = 10, ordering = FALSE, 
      searching = FALSE, scrollX = T))})
  })
  
  # Once data uploaded...
  observeEvent(input$dataFile, {
    # Update dropdowns with columns from data file
    choices <- colnames(data())
    updateSelectInput(session, "years", choices = choices, selected = "")
    updateSelectInput(session, "ages", choices = choices, selected = "")
    updateSelectInput(session, "selectivity", choices = choices, selected = "")
    updateSelectInput(session, "mat", choices = choices, selected = "")
    updateSelectInput(session, "waa", choices = choices, selected = "")
    updateSelectInput(session, "m_column", choices = choices, selected = "")
    
    # Clear year-based inputs
    updateSelectInput(session, "year_selmat", choices = "", selected = "")
    updateSelectInput(session, "year_survivor", choices = "", selected = "")
    updateSliderInput(session, "yearRange1", 
                      min = 1970, max = 2022, value = c(2000, 2010))
    updateSliderInput(session, "yearRange2",
                      min = 1970, max = 2022, value = c(2000, 2005))
    updateSliderInput(session, "yearRange3",
                      min = 1970, max = 2022, value = c(2005, 2005))
    
    # Clear reference point time series plots
    output$combined_time_series_plot <- renderPlot({
      combined_plot <- empty_plot("Year", expression("\u03A6"["0"]~"(kg/recruit)")) +
        empty_plot("Year", expression("SSB"["0"]~"(kt)")) +
        empty_plot("Year", expression("F"["MSY"])) +
        empty_plot("Year", expression("SSB"["MSY"]~"(kt)")) + 
        empty_plot("Year", bquote(F[.(input$SPR1)*"%SPR"])) +
        empty_plot("Year", bquote(SSB[.(input$SPR1)*"%SPR"]~(kt))) +
        plot_layout(nrow = 3, ncol = 2)
      print(combined_plot)
    })
    
    # Clear reference point barplots
    output$combinedPlot <- renderPlot({ 
      combined_plot <- empty_plot(expression("\u03A6"["0"]~"(kg/recruit)"),"Time Period") +
        empty_plot(expression("SSB"["0"]~"(kt)"),"Time Period") +
        empty_plot(expression("F"["MSY"]),"Time Period") +
        empty_plot(expression("SSB"["MSY"]~"(kt)"),"Time Period") + 
        empty_plot(bquote(F[.(input$SPR2)*"%SPR"]),"Time Period") +
        empty_plot(bquote(SSB[.(input$SPR2)*"%SPR"]~(kt)),"Time Period") +
        plot_layout(nrow = 3, ncol = 2)
      print(combined_plot)
    })
  })
  
  # Update year inputs to allow time interval covered by the data when get years data
  observeEvent(input$years, {
    req(input$years)
    year_values <- unique(data()[[input$years]])
    updateSliderInput(session, "yearRange1", min = min(year_values), 
                      max = max(year_values))
    updateSliderInput(session, "yearRange2", min = min(year_values), 
                      max = max(year_values))
    updateSliderInput(session, "yearRange3", min = min(year_values), 
                      max = max(year_values))
    updateSelectInput(session, "year_selmat", 
                      choices = year_values)
    updateSelectInput(session, "year_survivor", 
                      choices = year_values)
  })

  observeEvent(input$SPR1, {
    # When SPR1 changes, update SPR2 to match
    updateSliderInput(session, "SPR2", value = input$SPR1)
  })
  observeEvent(input$SPR2, {
    # When SPR2 changes, update SPR1 to match
    updateSliderInput(session, "SPR1", value = input$SPR2)
  })
  
  #####PLOT CODE####
  
  # Plot weight-at-age + maturity and selectivity
  observe({
    # Plot weight-at-age when have the necessary inputs
    if(input$waa!=""&& 
       input$years!=""&&input$ages!="")
    {
      # Fill in empirical weights if using empirical weight-at-age
        plot1 <- ggplot() + geom_path(size = 1.08, aes(y = data()[[input$waa]], 
                                                x = data()[[input$years]], 
                                                colour = as.factor(data()[[input$ages]]))) + 
            theme_classic() + scale_y_continuous("Weight (kg)") + 
            labs(x="Year",colour="Age (years)") + theme(text = element_text(size = font_size),
                                                        legend.text = element_text(size = font_size),
                                                        legend.position = "bottom")
    }
    else
    {
      # Render an empty weight-at-age plot
        plot1 <- empty_plot("Year", "Weight (kg)", "Further information required")
    }
    
    # If necessary inputs, plot maturity and/or selectivity
    if (input$years != "" && input$ages != "" && input$year_selmat %in% data()[[input$years]]
        &&(input$selectivity!=""|input$mat!="")) {
      n_ages <- length(unique(data()[[input$ages]]))
      
      # Initialize an empty ggplot
      plot2 <- ggplot() + labs(x = "Age (years)", y = "Proportion") + theme_classic()
      # If selectivity is available, add it to the plot
      if (input$selectivity != "") {
        sel_df <- data.frame(age = 1:n_ages, value = data()[[input$selectivity]]
                             [which(data()[[input$years]] == input$year_selmat)])
        plot2 <- plot2 +
          geom_path(data = sel_df, aes(x = age, y = value, color = "Selectivity"), size = 1.08) +
          geom_point(data = sel_df, aes(x = age, y = value, color = "Selectivity"), size = 3)
      }
      # If maturity is available, add it to the plot
      if (input$mat != "") {
        mat_df <- data.frame(age = 1:n_ages, value = data()[[input$mat]]
                             [which(data()[[input$years]] == input$year_selmat)])
        plot2 <- plot2 +
          geom_path(data = mat_df, aes(x = age, y = value, color = "Maturity"), size = 1.08) +
          geom_point(data = mat_df, aes(x = age, y = value, color = "Maturity"), size = 3)
      }
      # Customize the legend and x-axis
      plot2 <- plot2 +
        scale_color_manual(values = c("Selectivity" = palette[1], "Maturity" = palette[2]), 
                           name = "") + theme(legend.text = element_text(size = font_size),
                                              legend.position = "bottom",
                                              text = element_text(size = font_size)) +
        scale_x_continuous(breaks = 1:n_ages, labels = unique(data()[[input$ages]]))
    }
    else
    {
      # Render an empty selectivity/maturity plot
        plot2 <- empty_plot("Age (years)", "Proportion", "Further information required")
    }
    
    # Combine selectivity/maturity plot with weight-at-age plot
    output$growth_selmat_plot <- renderPlot({plot1+plot2})
  })
  
  # Calculate unfished survivorship-at-age if have the necessary info
  survivorship_values <- reactive({
    if(input$ages!="" && input$years!="" && 
       ((!is.na(input$mConstant) && input$mType == "constant") || 
        (!input$m_column == "" && input$mType != "constant")))
    {
      # Use either column or constant depending on user choice...
      # Return NULL if get error that arises when switching datasets
      result <- tryCatch({
        if (input$mType == "constant") 
        {
          m_value <- input$mConstant
        } 
        else 
        {
          m_value <- data()[[input$m_column]][which(data()[[input$years]] == 
                                                      input$year_survivor)]
        }
        survivorship(n_ages = length(unique(data()[[input$ages]])),
                     m = m_value,
                     selectivity = rep(0, length(unique(data()[[input$ages]]))),
                     f = 0)
      }, error = function(e) {NULL})
    }
    # Return null if missing information
    else
    {
      return(NULL)
    }
  })
  # Plot unfished survivorship-at-age if have the necessary info
  observe({
    if (!is.null(survivorship_values()))
    {
      n_ages <- length(unique(data()[[input$ages]]))
      output$survivorshipPlot <- renderPlot({
        ggplot() + 
          geom_path(size = 1.08, aes(x = 1:n_ages, 
                                     y = survivorship_values()), col = palette[1]) + 
          geom_point(size = 3, aes(x = 1:n_ages, 
                                   y = survivorship_values()), col = palette[1]) + 
          xlab("Age (years)") + theme_classic() + 
          scale_y_continuous(expression("l"["0"]~"-at-age"), limits = c(0,NA),expand=c(0,0)) +
          theme(text = element_text(size = font_size)) +
          scale_x_continuous(breaks = 1:n_ages,
                             labels = unique(data()[[input$ages]]))
      })
    }
    # If there isn't valid values, render empty plot
    else
    {
      output$survivorshipPlot <- renderPlot({empty_plot("Age (years)", 
                                                        expression("l"["0"]~"-at-age"),
                                                        "Further information required")})
    }
  })
  
  # Plot the SRR...
  observe({
    # ...when have the necessary values... either beverton-holt, ricker, or constant
    if((!is.na(input$recruitment_mean)&&input$srr == "constant")||
       (!is.na(input$alpha)&&!is.na(input$beta)&&input$srr != "constant"))
    {
      output$stockRecruitPlot <- renderPlot({
        if(input$srr == "beverton-holt")
        {
          srr_funct <- function(x) (input$alpha * x / (1 + input$beta * x))
        }
        else if(input$srr == "ricker")
        {
          srr_funct <- function(x) (input$alpha * x * exp(-input$beta*x))
        }
        else
        {
          srr_funct <- function(x) (input$recruitment_mean)
        }
        ggplot() +
          theme_classic() +
          labs(x = 'SSB (kt)', y = 'Recruitment (millions)') +
          scale_x_continuous(expand = c(0, 0), limits = c(0, input$max_x)) +
          scale_y_continuous(expand = c(0, 0), limits = c(0, input$max_y)) +
          geom_function(fun = srr_funct, 
                        colour = palette[1], size = 1.2) +
          theme(text = element_text(size = font_size))
      }) 
    }
    # Otherwise, render an empty plot
    else
    {
      output$stockRecruitPlot <- renderPlot({empty_plot("SSB (kt)", "Recruitment (billions)",
                                                        "Further information required")})
    }
  })
  
  # Render empty time series plots initially using default 40%SPR
  output$combined_time_series_plot <- renderPlot({
    combined_plot <- empty_plot("Year", expression("\u03A6"["0"]~"(kg/recruit)")) +
      empty_plot("Year", expression("SSB"["0"]~"(kt)")) +
      empty_plot("Year", expression("F"["MSY"])) +
      empty_plot("Year", expression("SSB"["MSY"]~"(kt)")) + 
      empty_plot("Year", bquote(F[.(SPR_default)*"%SPR"])) +
      empty_plot("Year", bquote(SSB[.(SPR_default)*"%SPR"]~(kt))) +
      plot_layout(nrow = 3, ncol = 2)
    print(combined_plot)
  })
  # Update combined time series plot once have values calculated
  observeEvent(input$calculateBtn, {
    if(!is.null(time_series_df()))
    {
      DF_RP <- time_series_df()
      
      # Generate plot
      combined_plot <- ts_plot(DF_RP[grepl("phi_0", DF_RP$variable), ], 
                               expression("\u03A6"["0"]~"(kg/recruit)"), 
                               palette, labels) +
        ts_plot(DF_RP[grepl("SSB_0", DF_RP$variable), ], 
                expression("SSB"["0"]~"(kt)"), 
                palette, labels) +
        ts_plot(DF_RP[grepl("F_MSY", DF_RP$variable), ], 
                expression("F"["MSY"]), 
                palette, labels) +
        ts_plot(DF_RP[grepl("SSB_MSY", DF_RP$variable), ], 
                expression("SSB"["MSY"]~"(kt)"), 
                palette, labels) +
        ts_plot(DF_RP[grepl("F_XSPR", DF_RP$variable)&!grepl("SSB", DF_RP$variable), ], 
                bquote(F[.(input$SPR1)*"%SPR"]), 
                palette, labels) +
        ts_plot(DF_RP[grepl("SSB_XSPR", DF_RP$variable), ], 
                bquote(SSB[.(input$SPR1)*"%SPR"]~(kt)), 
                palette, labels) + 
        plot_layout(nrow = 3, ncol = 2, guides = 'collect')
    }
    # If there isn't valid values after plot initialization, display blank plot again
    else
    {
      combined_plot <- empty_plot("Year", expression("\u03A6"["0"]~"(kg/recruit)")) +
        empty_plot("Year", expression("SSB"["0"]~"(kt)")) +
        empty_plot("Year", expression("F"["MSY"])) +
        empty_plot("Year", expression("SSB"["MSY"]~"(kt)")) + 
        empty_plot("Year", bquote(F[.(input$SPR1)*"%SPR"])) +
        empty_plot("Year", bquote(SSB[.(input$SPR1)*"%SPR"]~(kt))) +
        plot_layout(nrow = 3, ncol = 2)
    }
    output$combined_time_series_plot <- renderPlot({print(combined_plot)})
  })
  
  # Render an empty series of reference point barplots initially (using default 40%SPR)
  output$combinedPlot <- renderPlot({ 
    combined_plot <- empty_plot(expression("\u03A6"["0"]~"(kg/recruit)"),"Time Period") +
      empty_plot(expression("SSB"["0"]~"(kt)"),"Time Period") +
      empty_plot(expression("F"["MSY"]),"Time Period") +
      empty_plot(expression("SSB"["MSY"]~"(kt)"),"Time Period") + 
      empty_plot(bquote(F[.(SPR_default)*"%SPR"]),"Time Period") +
      empty_plot(bquote(SSB[.(SPR_default)*"%SPR"]~(kt)),"Time Period") +
      plot_layout(nrow = 3, ncol = 2)
    print(combined_plot)
  })
  # Update combined plot of reference points once have valid values calculated
  observeEvent(input$calculateBtn2, {
    if (!is.null(ref_points_df())) 
      {
        combined_plot <- RP_time_period_plot("phi_0", expression("\u03A6"["0"]~"(kg/recruit)"), 
                                          ref_points_df()) +
          RP_time_period_plot("SSB_0", expression("SSB"["0"]~"(kt)"), ref_points_df()) +
          RP_time_period_plot("F_MSY", expression("F"["MSY"]), ref_points_df()) +
          RP_time_period_plot("SSB_MSY", expression("SSB"["MSY"]~"(kt)"), ref_points_df()) +
          RP_time_period_plot("F_X%SPR", bquote(F[.(input$SPR2)*"%SPR"]), ref_points_df()) +
          RP_time_period_plot("SSB_X%SPR", bquote(SSB[.(input$SPR2)*"%SPR"]~(kt)), ref_points_df()) +
          plot_layout(nrow = 3, ncol = 2)
      }
      # If there isn't valid values after plot initialization, display blank plot again
      else
      {
        combined_plot <- empty_plot(expression("\u03A6"["0"]~"(kg/recruit)"),"Time Period") +
          empty_plot(expression("SSB"["0"]~"(kt)"),"Time Period") +
          empty_plot(expression("F"["MSY"]),"Time Period") +
          empty_plot(expression("SSB"["MSY"]~"(kt)"),"Time Period") + 
          empty_plot(bquote(F[.(input$SPR2)*"%SPR"]),"Time Period") +
          empty_plot(bquote(SSB[.(input$SPR2)*"%SPR"]~(kt)),"Time Period") +
          plot_layout(nrow = 3, ncol = 2)
      }
      output$combinedPlot <- renderPlot({print(combined_plot)}) 
    })

  #####REFERENCE POINT CALCULATIONS CODE####
    
  # Calculate time series of reference points... allow one component to change at
  # a time (and then allow them all to change)... must hit calculate button first
  # (and have necessary info)
  time_series_df <- eventReactive(input$calculateBtn, {
    
    # Initialize a vector to track missing information
    missing_items <- character(0) 
    # Check conditions and add missing items to the vector
    missing_items <- c(
      if (input$years == "") "years",
      if (input$ages == "") "age classes",
      if (input$waa == "") "weight-at-age",
      if (input$selectivity == "") "selectivity-at-age",
      if (input$mat == "") "maturity-at-age",
      if ((is.na(input$recruitment_mean) && input$srr == "constant") ||
          (is.na(input$alpha) && input$srr != "constant") ||
          (is.na(input$beta) && input$srr != "constant")) "recruitment",
      if ((is.na(input$mConstant) && input$mType == "constant") ||
          (input$m_column == "" && input$mType != "constant")) "natural mortality"
    )
    # Remove empty strings from the vector
    missing_items <- missing_items[missing_items != ""]
    
    # Move forward with calcs if have all info needed
    if (length(missing_items) == 0)
    {
      yr1 <- min(data()[[input$years]]) # Initial year
      yr2 <- max(data()[[input$years]]) # Terminal year
      
      # Count number of age classes
      n_ages <- length(unique(data()[[input$ages]]))
      
      # Define a dataframe to hold calculated reference point values 
      # for each year... loop through and fill it in
      DF_RP <- data.frame(yr = yr1:yr2)
      for (variables in c("M","selectivity", "waa", "mat", "all four"))
      {
        for (i in 1:length(yr1:yr2))
        {
          # Initialize variables to their initial yr values
          selectivity <- data()[[input$selectivity]][which(data()[[input$years]] 
                                                           == yr1)]
          waa <- data()[[input$waa]][which(data()[[input$years]] 
                                           == yr1)]
          mat <- data()[[input$mat]][which(data()[[input$years]] 
                                           == yr1)]
          if(input$mType == "constant")
          {
            m <- input$mConstant
          }
          else
          {
            m <- data()[[input$m_column]][which(data()[[input$years]] 
                                                == yr1)]
          }
          
          # Depending on the case, update the variables
          if (variables == "selectivity" || variables == "all four") {
            selectivity <- data()[[input$selectivity]][which(data()[[input$years]] 
                                                             == (yr1-1 + i))]
          }
          if (variables == "waa" || variables == "all four") {
            waa <- data()[[input$waa]][which(data()[[input$years]] 
                                             == (yr1-1 + i))]
          }
          if (variables == "mat" || variables == "all four") {
            mat <- data()[[input$mat]][which(data()[[input$years]] 
                                             == (yr1-1 + i))]
          }
          if (variables == "M" || variables == "all four") {
            if(input$mType == "constant")
            {
              m <- input$mConstant
            }
            else
            {
              m <- data()[[input$m_column]][which(data()[[input$years]] 
                                                  == yr1-1 + i)]
            }
          }
          
          # Compute and store reference point estimates for that year
          ref_pts <- equilibrium_RPs(n_ages = n_ages, 
            m = m, 
            selectivity = selectivity,
            waa = waa, mat = mat,
            alpha = input$alpha, 
            beta = input$beta,
            srr = if (input$srr == "constant") NULL else input$srr,
            recruitment_mean = recruitment_mean,
            SPR = input$SPR1/100)
          
          # Store the results in the dataframe
          DF_RP[i, paste(variables, "SSB_0", sep = "_")] <- unlist(ref_pts["SSB_0"])
          DF_RP[i, paste(variables, "SSB_MSY", sep = "_")] <- unlist(ref_pts["SSB_MSY"])
          DF_RP[i, paste(variables, "F_XSPR", sep = "_")] <- unlist(ref_pts["F_X%SPR"])
          DF_RP[i, paste(variables, "SSB_XSPR", sep = "_")] <- unlist(ref_pts["SSB_X%SPR"])
          DF_RP[i, paste(variables, "F_MSY", sep = "_")] <- unlist(ref_pts["F_MSY"])
          DF_RP[i, paste(variables, "phi_0", sep = "_")] <- unlist(ref_pts["phi_0"])
        }
      }
      # Get in long format to make plotting easier
      DF_RP <- DF_RP %>%
        gather(key = "variable", value = "value", -yr)
      
      # Hide selectivity from phi0 and SSB_0 plots
      DF_RP[DF_RP$variable %in% c("selectivity_SSB_0", "selectivity_phi_0"), "value"]  <- -1 
      
      return(DF_RP)
    }
    
    # If missing information
    else
    {
      # Show a notification to the user indicating which items are missing
      missing_message <- paste("Missing information on:", paste(missing_items, 
                                                                 collapse = ", "))
      showNotification(missing_message, type = "warning")
      return(NULL)
    }
  })

  # Calculate reference points using user-defined intervals of time
  ref_points_df <- eventReactive(input$calculateBtn2, {
    
    # Initialize a vector to track missing information
    missing_items <- character(0) 
    # Check conditions and add missing items to the vector
    missing_items <- c(
      if (input$years == "") "years",
      if (input$ages == "") "age classes",
      if (input$waa == "") "weight-at-age",
      if (input$selectivity == "") "selectivity-at-age",
      if (input$mat == "") "maturity-at-age",
      if ((is.na(input$recruitment_mean) && input$srr == "constant") ||
          (is.na(input$alpha) && input$srr != "constant") ||
          (is.na(input$beta) && input$srr != "constant")) "recruitment",
      if ((is.na(input$mConstant) && input$mType == "constant") ||
          (input$m_column == "" && input$mType != "constant")) "natural mortality"
    )
    # Remove empty strings from the vector
    missing_items <- missing_items[missing_items != ""]

    # Move forward with calcs if have all info needed
    if (length(missing_items) == 0)
    {
      ref_points <- average_time_period_RPs(
        DATA = data(),
        years = input$years,
        average_yrs_list = list(c(input$yearRange1[1]:input$yearRange1[2]),
                                c(input$yearRange2[1]:input$yearRange2[2]),
                                c(input$yearRange3[1]:input$yearRange3[2])),
        ages = input$ages,
        m = if (input$mType == "constant") input$mConstant else input$m_column,
        selectivity = input$selectivity,
        waa = input$waa,
        mat = input$mat,
        srr = if (input$srr == "constant") NULL else input$srr,
        alpha = input$alpha,
        beta = input$beta,
        recruitment_mean = input$recruitment_mean,
        SPR = input$SPR2/100)
      return(ref_points)
    }
    # If missing information
    else
    {
      # Show a notification to the user indicating which items are missing
      missing_message <- paste("Missing information on:", paste(missing_items, 
                                                                collapse = ", "))
      showNotification(missing_message, type = "warning")
      return(NULL)
    }
  })
}

# Run the Shiny app
shinyApp(ui, server)