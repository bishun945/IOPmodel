library(shiny)
library(shinyBS)
library(IOPmodel)
library(magrittr)

library(patchwork)
library(ggplot2)
library(data.table)
# library(ggvis)

#### pre-functions ####
inline_numericInput <- function(inputId, label, value, ...) {
  list(
    div(style = "display: inline-block; vertical-align: sub; width: 40%",
        strong(paste0(label, ": "))),
    div(style = "display: inline-block; vertical-align: sub; width: 50%",
        numericInput(inputId, label = NULL, value = value, ...)),
    br()
  )
}


inline_sliderInput <- function(inputId, label, value, ...) {
  list(
    div(style = "display: inline-block; text-align: left; vertical-align: middle; padding: 0px; width: 20%",
        strong(paste0(label, ": "))),
    div(style = "display: inline-block; text-align: left; vertical-align: middle; padding: 0px; width: 70%",
        sliderInput(
          inputId, label = NULL, value = value, ...
        )),
    br()
  )
}

inline_updateSliderInput <- function(session, inputId, label, value, ...) {
  list(
    div(style = "display: inline-block; text-align: left; vertical-align: middle; padding: 0px; width: 20%",
        strong(paste0(label, ": "))),
    div(style = "display: inline-block; text-align: left; vertical-align: middle; padding: 0px; width: 70%",
        updateSliderInput(
          session, inputId, label = NULL, value = value, ...
        )),
    br()
  )
}

show_spec <- function(tmp) {

  tmp <- tmp[between(wavelen,  400, 750)]

  Rrs_scale <- max(c(tmp$bp, tmp$aph * 5, tmp$ad * 5,
                     tmp$acdom * 5)) * 0.7 / max(c(tmp$Rrs))

  dt_label <- data.table(
    x = 650,
    y = max(c(tmp$bp, tmp$aph * 5, tmp$ad * 5, tmp$acdom * 5)) * 1.3,
    label = sprintf(paste0(c("[Chl] = %.2f mg/m<sup>3</sup><br>",
                             "[ISM] = %.2f g/m<sup>3</sup><br>",
                             "*a*<sub>*g*</sub>(440) = %.4f m<sup>-1</sup>"),
                           collapse = ""),
                    tmp$Chl[1], tmp$ISM[1], tmp$ag440[1])
  )

  y_title_left <- "IOP (m<sup>-1</sup>)"
  y_title_right <- "*R*<sub>*rs*</sub> (sr<sup>-1</sup>)"

  tmp %>%
    ggplot(aes(x = wavelen)) +
    geom_path(aes(y = bp, col = "bp")) +
    geom_path(aes(y = bph, col = "bph")) +
    geom_path(aes(y = aph * 5, col = "aph")) +
    geom_path(aes(y = ad * 5, col = "ad")) +
    geom_path(aes(y = acdom * 5, col = "ag")) +
    geom_path(aes(y = Rrs * Rrs_scale, col = "Rrs")) +
    ggtext::geom_richtext(inherit.aes = FALSE, data = dt_label,
                          aes(x = 400, y = y, label = label),
                          hjust = 0, vjust = 1,
                          label.padding = unit(rep(0, 4), "lines"),
                          label.color = NA,
                          fill = alpha("white", 0.7)) +
    scale_x_continuous(name =  "Wavelength (nm)",
                       limits = c(400, 750),
                       breaks = seq(300, 900, 100)) +
    scale_y_continuous(name = y_title_left,
                       sec.axis =
                         sec_axis(~./Rrs_scale,
                                  name = y_title_right)) +
    scale_color_manual(
      "Component",
      values = c(
        "aph" = "orange",
        "ad" = "red",
        "ag" = "cyan",
        "bp" = "black",
        "bph" = "green",
        "Rrs" = "blue"
      ),
      labels = c(
        "aph" = "*a*<sub>*ph*</sub> x 5",
        "ad" = "*a*<sub>*d*</sub> x 5",
        "ag" = "*a*<sub>*g*</sub> x 5",
        "bp" = "*b*<sub>*p*</sub>",
        "bph" = "*b*<sub>*ph*</sub>",
        "Rrs" = "*R*<sub>*rs*</sub>"
      )
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.title.y.left = ggtext::element_markdown(),
      axis.title.y.right = ggtext::element_markdown(color = "blue"),
      axis.text.y.right = ggtext::element_markdown(color = "blue"),
      legend.text = ggtext::element_markdown(size = 12),
      legend.title = ggtext::element_markdown(),
      legend.position = "bottom"
    )

}

can_convert_to_numeric <- function(x) {
  all(grepl('^(?=.)([+-]?([0-9]*)(\\.([0-9]+))?)$', x, perl = TRUE))
}

#### ui ####

ui <- fluidPage(

  shinyjs::useShinyjs(),

  titlePanel("IOP model"),

  tabsetPanel(

    tabPanel(

      "Run",

      fluidRow(
        column(
          3,
          wellPanel(

            list(
              div(style = "display: inline-block; vertical-align: sub",
                  h4("Water constitute")),
              div(style = "display: inline-block; vertical-align: sub",
                  actionButton("oac_dice", NULL, icon = icon("dice"))),
              br()
            ),

            # h4("Water constitute"),

            # actionButton("oac_dice", NULL, icon = icon("dice")),

            # list(
            #   div(style = "display: inline-block; vertical-align: sub",
            #       actionButton("oac_dice", NULL, icon = icon("dice"))),
            #   div(style = "display: inline-block; vertical-align: sub",
            #       checkboxInput("logOAC", "log10 scale?", value = FALSE)),
            #   br()
            # ),

            sliderInput("Chl", "Chlorophyll a concentration [mg/m^3]",
                        min = 0.00, max = 1000, value = 0.3, step = 0.001),
            textInput("Chl_text", value = 0.3, label = NULL),


            sliderInput("ISM", "Inorganic suspended matter concentration [g/m^3]",
                        min = 0.000, max = 2000, value = 0.5, step = 0.001),
            textInput("ISM_text", value = 0.5, label = NULL),


            sliderInput("ag440", "CDOM absorption coefficient at 440 nm [m^-1]",
                        min = 0.0000, max = 50, value = 0.015, step = 0.0001),
            textInput("ag440_text", value = 0.015, label = NULL),


            bsTooltip("Chl", "Concentration of chlorophyll a, Chl", "right"),
            bsTooltip("ISM", "Concentration of inorganic suspended matter, ISM", "right"),
            bsTooltip("ag440", "CDOM absorption coefficient at 440 nm [m^-1], ag440", "right"),

          ),

          wellPanel(
            h4("Pure water"),

            list(
              div(style = "display: inline-block; vertical-align: sub",
                  "aw version:"),
              div(style = "display: inline-block; vertical-align: sub",
                  radioButtons("aw_version", NULL,
                               c("1" = 1, "2" = 2, "3" = 3),
                               selected = 3,
                               inline = TRUE)),
              br()
            ),
            bsTooltip("aw_version", "Version of pure water absorption:<br>(1) 300-420 nm by Ed Fry et al; <br>(2) 300-420 nm by Morel et al. 2007;<br>(3) 300-510 nm by Mason et al. 2016", "right"),
            sliderInput("Temp", "Water temperature [degC]",
                        min = 0, max = 30, value = 15, step = 0.1),
            sliderInput("Sal", "Water salinity [PSU]",
                        min = 0, max = 30, value = 10, step = 0.1)
          )

        ),
        column(
          3,

          wellPanel(

            h4("Phytoplankton group fraction"),

            fluidPage(

              column(
                6,
                radioButtons("rand_frac_case", "Class option",
                             c("Case-1" = 1,
                               "Case-2" = 2,
                               "Custom" = 3),
                             inline = FALSE)
              ),

              column(
                6,
                h5("Roll the dice for random fraction values:"),
                actionButton("frac_dice", NULL, icon = icon("dice"))
              ),

            ),


            br(),

            inline_sliderInput("Brown", "Brown", min = 0, max = 1, value = 1/7, step = 0.01),
            inline_sliderInput("Green", "Green", min = 0, max = 1, value = 1/7, step = 0.01),
            inline_sliderInput("Crypt", "Crypt", min = 0, max = 1, value = 1/7, step = 0.01),
            inline_sliderInput("CyanB", "CyanB", min = 0, max = 1, value = 1/7, step = 0.01),
            inline_sliderInput("CyanR", "CyanR", min = 0, max = 1, value = 1/7, step = 0.01),
            inline_sliderInput("Cocco", "Cocco", min = 0, max = 1, value = 1/7, step = 0.01),
            inline_sliderInput("PhyC1", "PhyC1", min = 0, max = 1, value = 1/7, step = 0.01),

            bsTooltip("Brown", "Brown-colored phytoplankton group (Heterokontophyta [including diatoms], Dinophyta, and Haptophyta)", "right"),
            bsTooltip("Green", "Green-colored phytoplankton group (Chlorophyta)", "right"),
            bsTooltip("Crypt", "Cryptophytes phytoplankton group (Cryptophyceae)", "right"),
            bsTooltip("CyanB", "Blue-green-colored Cyanobacteria phytoplankton group (Cyanobacteria)", "right"),
            bsTooltip("CyanR", "Red-colored Cyanobacteria phytoplankton group (Cyanobacteria)", "right"),
            bsTooltip("Cocco", "Coccolithophores phytoplankton group (Coccolithus huxleyi)", "right"),
            bsTooltip("PhyC1", "Phytoplankton group for oligotrophic Case-1 water (Synechococcus and Prochlorococcus)", "right"),

            h5("If you selected 'Custom' and changed the sliders, please click 'Rescale':"),
            actionButton("frac_rescale", "Rescale", icon = icon("scale-balanced"))

            # place for pie plot
            # textOutput("test")

          )

        ),
        column(
          5,
          wellPanel(

            h4("Run option"),

            list(
              # div(style = "display: inline-block; vertical-align: sub; width: 10%",
              div(style = "display: inline-block; vertical-align: sub",
                  "Mode:"),
              # div(style = "display: inline-block; vertical-align: sub; width: 80%",
              div(style = "display: inline-block; vertical-align: sub",
                  radioButtons("which_term", NULL,
                               c("Four-term" = 4, "Two-term" = 2),
                               inline = TRUE)),
              br()
            ),

            actionButton("Run", "Run IOP model", icon = icon("play")),
          ),
          wellPanel(
            h4("Result visualization"),
            # textOutput("test")
            plotOutput("model_plot", height = "500px")
          )
        )
      )

    ),

    tabPanel(

      "Table",

      div(
        style = "display: inline-block",
        strong("You can download the model results for further analysis:   "),
        downloadButton("downloadData", NULL)
        ),
      tableOutput("model_table")

    ),

    tabPanel(

      "Settings",

      column(
        3,
        wellPanel(
          h4("Water constitute range"),
          inline_numericInput("Chl_min", "Chl min", value = 0),
          inline_numericInput("Chl_max", "Chl max", value = 1e3),
          inline_numericInput("Chl_step", "Chl step", value = 0.001),

          inline_numericInput("ISM_min", "ISM min", value = 0),
          inline_numericInput("ISM_max", "ISM max", value = 2e3),
          inline_numericInput("ISM_step", "ISM step", value = 0.001),

          inline_numericInput("ag440_min", "ag440 min", value = 0),
          inline_numericInput("ag440_max", "ag440 max", value = 50),
          inline_numericInput("ag440_step", "ag440 step", value = 0.0001),
        )
      ),

      column(
        3,
        wellPanel(
          h4("Model parameter", actionButton("model_parm_dice", NULL, icon = icon("dice"))),
          # actionButton("model_parm_dice", NULL, icon = icon("dice")),
          sliderInput("A_d", "Detritus single scattering albedo",
                      min = 0, max = 1, value = 0.9542, step = 0.0001),
          sliderInput("G_d", "Power law exponent of attenuation of detritus",
                      min = 0, max = 1, value = 0.3835, step = 0.0001),
          sliderInput("qt_bd", "Quantile values of biogenic detritus absorption",
                      min = 0, max = 1, value = 0.5, step = 0.1),
          sliderInput("qt_md", "Quantile values of minerogenic detritus absorption",
                      min = 0, max = 1, value = 0.5, step = 0.1),
          checkboxInput("varyA", "Vary A of aph676=Ax[Chl]^E?", value = FALSE),
          numericInput("ag_seed", "Seed of aCDOM library (Four-term):", 8888),
          numericInput("Two_term_seed", "Seed of Two-term model:", 8888),
          # inline_numericInput("ag_seed", "Seed number for aCDOM library (Four-term)", value = 8888),
          # inline_numericInput("Two_term_seed", "Seed number for Two-term model", value = 8888),

          bsTooltip("varyA", "The default A is 0.0237. If checked, value varis between 0.0112 and 0.0501", "right"),

        )
      ),

      column(
        3,
        wellPanel(
          h4("Lifecycle of coccolithohpores"),
          sliderInput("cc_frac", "Fraction of Cocco absorption",
                      min = 0, max = 1, value = 1, step = 0.001),
          textInput("cc_frac_text", value = 1, label = NULL),
          p("Change default parameters to", actionButton("cc_parm", "Cocco")),
          p(em("Let me know if you have any preference of default settings"))
        )
      )

    ),

    tabPanel(

      "HelpDoc",

      h4("Here are some quick intro for this tool"),

      p("Three columns in 'Run' panel: ",
        tags$ul(
          tags$li("the first controls water constitue concentration and model parameters."),
          tags$li("the second is for phytoplankton group fraction and pure water settings."),
          tags$li("the third part inlcudes model running setups and result visualization")),
        ),

      p("In the 'Table' panel, you can check and download the model results after running."),

      br(),

      h5("Class option for different water cases."),
      p("In the model developing, we assume some biological constraints for the Case-1 water, PhyC1 is dominant when Chl is low and Brown is then mixed up when Chl is increasing a bit."),
      p("But for the Case-2 water, no constraints for all groups."),
      p("In 'Custom' option, you can also modify slider bars for different groups and click 'Rescale' to re-normalize them."),

      br(),

      h5("'Four-term' vs 'Two-term'"),
      p("We define 'Four-term' as the IOP model using four components, i.e., pure water, phytoplankton, detritus, and CDOM and 'Two-term' as one using two components, i.e., pure water and phytoplankton."),
      p("Components in the Four-term model are not covariant but all constitues in the Two-term model (except for pure water) are covaried with phytoplankton concentration."),
      p("In principle, one can simulate the similar results by the Two-term model in the Four-term."),
      p("The reason we develop this Two-term model is to effectively simulate so-called Case-1 waters which is FAR more than Case-2 waters."),

      br(),
      h5("For people who are lazy to choose parameters like me"),
      p("There are some dices ", icon("dice"), " provided for several chunks. Random values will be generated when you click them. "),
      p("These generator are based on uniform and normal distributions, parameters of which are trained from many published data sets."),

      br(),
      h5("Updates for Model Parameter 2024-11-17"),
      p("Thanks Yulun for pointing out a bug that sometimes backscattering of total particles (bp) will randomly distribute for the constant OAC values."),
      p("This is simply because the improper parameters leading to negatibe bp values and the IOP model will regenerate them until bp is positive."),
      p("Now the app will raise a warning and update the used model parameters once this error is catched."),

      h5("To be finished..."),

      br(),
      p("If you want ot know more details about the model, please refere to the first in the reference list."),

      p("If you have any questions or find any bugs, please email me via", code("Shun.Bi@outlook.com"),
        "or open an issue in ", a("the GitHub page", href = "https://github.com/bishun945/IOPmodel/issues"),
        "with a 'shiny' label"),

      p("Have fun!"),

      p(em(sprintf("Last updated on %s by", "2024-11-17"),
        tags$a("Shun", href = "https://bishun945.github.io/CV/")))

    ),

    tabPanel(

      "Workflow",

      tags$img(src = "workflow.jpg", width = 800)

    ),

    tabPanel(

      "Reference",

      tags$div(
        tags$h4("Some references related to the model:"),
        tags$ul(
          tags$li(a("Bi, S., M. Hieronymi, and R. Röttgers. 2023. Bio-geo-optical modelling of natural waters. Front. Mar. Sci. 10: 1196352. doi:10.3389/fmars.2023.1196352",
                    href = "https://doi.org/10.3389/fmars.2023.1196352")),
          br(),
          tags$li(a("Bi, S., and M. Hieronymi. 2024. Holistic optical water type classification for ocean, coastal, and inland waters. Limnology & Oceanography lno.12606. doi:10.1002/lno.12606",
                    href = "https://aslopubs.onlinelibrary.wiley.com/doi/10.1002/lno.12606")),
          br(),
          tags$li(a("R package 'IOPmodel'",
                    href = "https://github.com/bishun945/IOPmodel")),
          br(),
          tags$li(a("Python package 'pyOWT'",
                    href = "https://github.com/bishun945/pyOWT")),
          br(),
          tags$li(a("WOPP processor: Röttgers, R, R Doerffer, D McKee, and W Schönfeld. “The Water Optical Properties Processor (WOPP): Pure Water Spectral Absorption, Scattering and Real Part of Refractive Index Model.” Technical Report No WOPP-ATBD/WRD6, 2016. https://calvalportal.ceos.org/tools.",
                    href = "https://calvalportal.ceos.org/tools")),
          br(),
          tags$li(a("Data set used in this model: Röttgers, Rüdiger, Shun Bi, Henning Burmester, Kerstin Heymann, Martin Hieronymi, Hajo Krasemann, and Wolfgang Schönfeld. “A Data Set of Water Inherent Optical Properties and Concentrations of Water Constituents from the German Bight and Adjacent Regions.” PANGAEA, 2022.",
                    href = "https://doi.pangaea.de/10.1594/PANGAEA.950774"))
        )
      )

    )

  )

)


# Define server logic required to draw a histogram
server <- function(input, output, session) {

  #### GUI setup ####
  observeEvent(input$Chl_min | input$Chl_max | input$Chl_step, {
    updateSliderInput(session = session,
                      inputId = "Chl",
                      label = "Chlorophyll a concentration [mg/m^3]",
                      value = input$Chl,
                      min = input$Chl_min, max = input$Chl_max, step = input$Chl_step)
  })

  observeEvent(input$ISM_min | input$ISM_max | input$ISM_step, {
    updateSliderInput(session = session,
                      inputId = "ISM",
                      label = "Inorganic suspended matter concentration [g/m^3]",
                      value = input$ISM,
                      min = input$ISM_min, max = input$ISM_max, step = input$ISM_step)
  })

  observeEvent(input$ag440_min | input$ag440_max | input$ag440_step, {
    updateSliderInput(session = session,
                      inputId = "ag440",
                      label = "CDOM absorption coefficient at 440 nm [m^-1]",
                      value = input$ag440,
                      min = input$ag440_min, max = input$ag440_max, step = input$ag440_step)
  })


  #### Synchronous slider and text inputs ####
  # Chl
  observeEvent(input$Chl, {
    if(can_convert_to_numeric(input$Chl_text)) {
      if(input$Chl != as.numeric(input$Chl_text)) {
        updateTextInput(session = session,
                        inputId = "Chl_text",
                        value = input$Chl)
      }
    }
  })
  observeEvent(input$Chl_text, {
    if(can_convert_to_numeric(input$Chl_text)) {
      if(input$Chl != as.numeric(input$Chl_text)) {
        updateSliderInput(session = session,
                          inputId = "Chl",
                          label = "Chlorophyll a concentration [mg/m^3]",
                          value = input$Chl_text,
                          min = input$Chl_min, max = input$Chl_max, step = input$Chl_step)
      }
    }
  })

  # ISM
  observeEvent(input$ISM, {
    if(can_convert_to_numeric(input$ISM_text)) {
      if(input$ISM != as.numeric(input$ISM_text)) {
        updateTextInput(session = session,
                        inputId = "ISM_text",
                        value = input$ISM)
      }
    }
  })
  observeEvent(input$ISM_text, {
    if(can_convert_to_numeric(input$ISM_text)) {
      if(input$ISM != as.numeric(input$ISM_text)) {
        updateSliderInput(session = session,
                          inputId = "ISM",
                          label = "Inorganic suspended matter concentration [g/m^3]",
                          value = input$ISM_text,
                          min = input$ISM_min, max = input$ISM_max, step = input$ISM_step)
      }
    }
  })

  # ag440
  observeEvent(input$ag440, {
    if(can_convert_to_numeric(input$ag440_text)) {
      if(input$ag440 != as.numeric(input$ag440_text)) {
        updateTextInput(session = session,
                        inputId = "ag440_text",
                        value = input$ag440)
      }
    }
  })
  observeEvent(input$ag440_text, {
    if(can_convert_to_numeric(input$ag440_text)) {
      if(input$ag440 != as.numeric(input$ag440_text)) {
        updateSliderInput(session = session,
                          inputId = "ag440",
                          label = "CDOM absorption coefficient at 440 nm [m^-1]",
                          value = input$ag440_text,
                          min = input$ag440_min, max = input$ag440_max, step = input$ag440_step)
      }
    }
  })

  # cc_frac
  observeEvent(input$cc_frac, {
    if(can_convert_to_numeric(input$cc_frac_text)) {
      if(input$cc_frac != as.numeric(input$cc_frac_text)) {
        updateTextInput(session = session,
                        inputId = "cc_frac_text",
                        value = input$cc_frac)
      }
    }
  })
  observeEvent(input$cc_frac_text, {
    if(can_convert_to_numeric(input$cc_frac_text)) {
      if(input$cc_frac != as.numeric(input$cc_frac_text)) {
        updateSliderInput(session = session,
                          inputId = "cc_frac",
                          label = "Fraction of Cocco absorption",
                          value = input$cc_frac_text,
                          min = 0, max = 1, step = 0.001)
      }
    }
  })

  #### Cocco ####
  # change OAC ranges to Cocco theme
  observeEvent(input$cc_parm, {

    updateTextInput(session, "Chl_max", NULL, value = 30)
    updateTextInput(session, "ISM_max", NULL, value = 10)
    updateTextInput(session, "ag440_max", NULL, value = 1)

    updateTextInput(session, "Chl_text", NULL, value = 10)
    updateTextInput(session, "ISM_text", NULL, value = 0)
    updateTextInput(session, "ag440_text", NULL, value = 0.0015)

    updateTextInput(session, "cc_frac_text", NULL, value = 0.01)

    group_name <- c("Brown", "Green", "Crypt", "CyanB", "CyanR",
                    "Cocco", "PhyC1")
    group_frac <- round(as.numeric(IOPmodel:::def_frac("Cocco", 1)), 2)

    # print(group_frac)

    for(i in 1:length(group_name)) {
      inline_updateSliderInput(session = session,
                               inputId = group_name[i],
                               label = group_name[i],
                               value = group_frac[i], min = 0, max = 1, step = 0.01)
    }


  })


  # Add more OAC range themes here... ####


  #### Rescale fraction for any inputs ####
  observeEvent(input$frac_rescale, {

    frac_value <- c(input$Brown, input$Green, input$Crypt, input$CyanB,
                    input$CyanR, input$Cocco, input$PhyC1)

    frac_value <- round(frac_value / sum(frac_value), 2)

    group_name <- c("Brown", "Green", "Crypt", "CyanB", "CyanR",
                    "Cocco", "PhyC1")

    for(i in 1:length(group_name)) {
      inline_updateSliderInput(session = session, inputId = group_name[i],
                               label = group_name[i],
                               value = frac_value[i], min = 0, max = 1, step = 0.01)
    }

    # output$test <- renderText({paste(frac_value, collapse = ", ")})

  })

  #### Dice for Phytoplankton group fraction ####
  observeEvent(input$frac_dice, {

    case = input$rand_frac_case
    Chl  = input$Chl

    if(case == 1 | case == 2) {

      frac_value <- as.numeric(round(IOPmodel:::rand_frac(case, Chl), 2))

      group_name <- c("Brown", "Green", "Crypt", "CyanB", "CyanR",
                      "Cocco", "PhyC1")

      for(i in 1:length(group_name)) {
        inline_updateSliderInput(session = session, inputId = group_name[i],
                                 label = group_name[i],
                                 value = frac_value[i], min = 0, max = 1, step = 0.01)
      }

      # output$test <- renderText({paste(frac_value, collapse = ", ")})

    } else {

      showNotification("You selected 'Custom' so the dice is not working!", type = "warning")

    }

  })

  observeEvent(input$rand_frac_case, {

    if(input$rand_frac_case == 3) {

      showNotification("'Custom' selected! Remember to click 'Resacle'", type = "message")

    }

  })


  #### Dice for model parameters ####
  observeEvent(input$model_parm_dice, {

    qt_bd <- round(runif(1, 0, 1), 1)
    qt_md <- round(runif(1, 0, 1), 1)

    A_mean <- -1.3390
    A_sd <- 0.0618
    G_mean <-  0.3835
    G_sd <- 0.1277

    A_d <- rnorm(1, mean = A_mean, sd = A_sd)

    G_d <- IOPmodel:::rnorm_bound(
      1,
      mean = G_mean,
      sd = G_sd,
      lo = 0,
      up = NA
    )

    A_d <- 1 - 10^A_d

    updateSliderInput(session, "A_d", "Detritus single scattering albedo",
                      min = 0, max = 1, value = A_d, step = 0.0001)
    updateSliderInput(session, "G_d", "Power law exponent of attenuation of detritus",
                      min = 0, max = 1, value = G_d, step = 0.0001)
    updateSliderInput(session, "qt_bd", "Quantile values of biogenic detritus absorption",
                      min = 0, max = 1, value = qt_bd, step = 0.1)
    updateSliderInput(session, "qt_md", "Quantile values of minerogenic detritus absorption",
                      min = 0, max = 1, value = qt_md, step = 0.1)


  })

  #### Dice for OAC concentrations ####
  observeEvent(input$oac_dice, {

    Chl <- round(10^IOPmodel:::rnorm_bound(1, 0.121, 0.858, log10(0.02), log10(1000)), 2)
    ISM <- round(10^IOPmodel:::rnorm_bound(1, 0.386, 1.637, log10(0.001), log10(2000)), 3)
    ag440 <- round(10^IOPmodel:::rnorm_bound(1, -0.671, 0.969, log10(0.0001), log10(50)), 4)

    updateSliderInput(session, "Chl", "Chlorophyll a concentration [mg/m^3]",
                      min = 0.00, max = 1000, value = Chl, step = 0.01)
    updateSliderInput(session, "ISM", "Inorganic suspended matter concentration [g/m^3]",
                      min = 0.000, max = 2000, value = ISM, step = 0.001)
    updateSliderInput(session, "ag440", "CDOM absorption coefficient at 440 nm [m^-1]",
                      min = 0.0000, max = 50, value = ag440, step = 0.0001)

    # updateSliderInput(session, "Chl", "Chlorophyll a concentration [mg/m^3]",
    #                   min = input$Chl_min, max = input$Chl_max,
    #                   value = Chl, step = input$Chl_step)
    # updateSliderInput(session, "ISM", "Inorganic suspended matter concentration [g/m^3]",
    #                   min = input$ISM_min, max = input$ISM_max,
    #                   value = input$ISM, step = input$ISM_step)
    # updateSliderInput(session, "ag440", "CDOM absorption coefficient at 440 nm [m^-1]",
    #                   min = input$ag440_min, max = input$ag440_max,
    #                   value = input$ag440, step = input$ag440_step)

    updateTextInput(session, "Chl_text", value = Chl, label = NULL)
    updateTextInput(session, "ISM_text", value = ISM, label = NULL)
    updateTextInput(session, "ag440_text", value = ag440, label = NULL)


  })

  #### Switch Four-term or Two-term model ####
  observeEvent(input$which_term, {

    which_term = input$which_term

    if(which_term == 4) {

      showNotification("Four-term model selected! Components are not co-varied!", type = "message")

      shinyjs::enable("ISM")
      shinyjs::enable("ISM_text")

      shinyjs::enable("ag440")
      shinyjs::enable("ag440_text")

      shinyjs::enable("A_d")
      shinyjs::enable("G_d")
      shinyjs::enable("qt_bd")
      shinyjs::enable("qt_md")

      shinyjs::enable("model_parm_dice")
      # shinyjs::enable("oac_dice")


    } else if (which_term == 2) {

      showNotification("Two-term model selected! Detritus and CDOM are co-varying with Chl!\nDetritus and CDOM related widgets are disabled!", type = "message")

      shinyjs::disable("ISM")
      shinyjs::disable("ISM_text")

      shinyjs::disable("ag440")
      shinyjs::disable("ag440_text")

      shinyjs::disable("A_d")
      shinyjs::disable("G_d")
      shinyjs::disable("qt_bd")
      shinyjs::disable("qt_md")

      shinyjs::disable("model_parm_dice")
      # shinyjs::disable("oac_dice")

    }

  })


  #### Click-run ####
  rv <- reactiveValues()

  # model_run <- eventReactive(input$Run, {
  observeEvent(input$Run, {

    showNotification("Model is running!", type = "message", duration = 5)

    frac_phyto <-
      c(input$Brown, input$Green, input$Crypt, input$CyanB,
        input$CyanR, input$Cocco, input$PhyC1)

    if(sum(frac_phyto) < 0.98) {
      showNotification("The sum of fraction less than 1! Did you click 'Rescale'?",
                       type = "warning", duration = 5)
    }

    a_frac <- IOPmodel:::def_frac("Cocco", input$cc_frac)
    a_frac[-6] <- 1

    if(input$which_term == 4) {

      res <-
        IOP_four_comp(
          Chl   = input$Chl,
          ag440 = input$ag440,
          ISM   = input$ISM,
          Temp  = input$Temp,
          Sal   = input$Sal,
          qt_bd = input$qt_bd,
          qt_md = input$qt_md,
          frac_phyto = frac_phyto,
          A_d   = input$A_d,
          G_d   = input$G_d,
          varyA = input$varyA,
          # ag_seed = NULL,
          ag_seed = input$ag_seed,
          aw_version = input$aw_version,
          a_frac = a_frac
        )

      attr(res, "term") <- 4

      if (attr(res, "warn_neg_bd")) {
        showNotification("The calculated bd values based on inputs have negatives. Random 'Model Parameters' of 'Settings' have been updated to keep bd positve.",
                         type = "warning", duration = 5)

        if(res$parm$A_d < 0) {
          A_d = 1 - 10^res$parm$A_d
        } else {
          A_d = res$parm$A_d
        }

        # They should be updated correspondingly...
        updateSliderInput(session, "A_d", "Detritus single scattering albedo",
                          min = 0, max = 1, value = A_d, step = 0.0001)
        updateSliderInput(session, "G_d", "Power law exponent of attenuation of detritus",
                          min = 0, max = 1, value = res$parm$G_d, step = 0.0001)
        updateSliderInput(session, "qt_bd", "Quantile values of biogenic detritus absorption",
                          min = 0, max = 1, value = res$parm$qt_bd, step = 0.1)
        updateSliderInput(session, "qt_md", "Quantile values of minerogenic detritus absorption",
                          min = 0, max = 1, value = res$parm$qt_md, step = 0.1)
      }

    } else if(input$which_term == 2) {

      res <-
        IOP_two_comp(
          Chl   = input$Chl,
          Temp  = input$Temp,
          Sal   = input$Sal,
          frac_phyto = frac_phyto,
          varyA = input$varyA,
          aw_version = input$aw_version,
          # seed = NULL,
          seed = input$Two_term_seed,
          a_frac = a_frac
        )

      attr(res, "term") <- 2

    }

    rv$res <- res

    return(res)

  })

  # output$test <- renderText({paste(res$parm, collapse = ", ")})


  #### Click-plot ####
  output$model_plot <- renderPlot({

    if(!is.null(rv$res)) {

      res <- rv$res

      if(attr(res, "term") == 4) {

        spec <- as.data.table(res$spec) %>%
          dplyr::mutate(
            Rrs = Rrs_L11,
            Chl = res$parm$Chl,
            ag440 = res$parm$ag440,
            ISM = res$parm$ISM
          )

      } else if(attr(res, "term") == 2) {

        spec <- as.data.table(res$spec) %>%
          dplyr::mutate(
            Rrs = Rrs_L11,
            Chl = res$parm$Chl,
            ag440 = as.data.table(res$spec)[wavelen == 440, acdom],
            ISM = -99
          )

      }

      p <- show_spec(spec)

      showNotification("Done!", type = "message")

      print(p)

    }

  })

  #### Click-table ####
  output$model_table <- renderTable({

    if(!is.null(rv$res)) {
      res <- rv$res
      as.data.table(res$spec) %>%
        dplyr::mutate(wavelen = as.integer(wavelen))
    }

  }, striped = TRUE, hover = TRUE, spacing = "xs", digits = 6, width = "500px")

  observeEvent(input$Run, {
    if(is.null(input$Run)) {
      shinyjs::disable("downloadData")
    } else {
      shinyjs::enable("downloadData")
    }
  })

  output$downloadData <- downloadHandler(

    filename = "IOP_spectra_outputs.csv",

    content = function(file) {
      if (!is.null(rv$res)) {
        res <- rv$res
        write.csv(as.data.table(res$spec), file, row.names = FALSE)
      }
    }

  )

}


# Run the application
shinyApp(ui = ui, server = server)









