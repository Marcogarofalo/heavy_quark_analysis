
library(shiny)
library(Rose)
library(ggplot2)
library(plotly)

library(shiny)
library(bslib)

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      selectInput(inputId = 'ensemble', label='ensemble', 
                  c("cB211.072.64", "cB211.072.96", "cC211.06.80",
                    "cD211.054.96", "cE211.044.112" )
      ),
      selectInput(inputId = 'Z', label='Z', 
                  c("Z0", "Z1", "Z2","Z0_part")
      ),
      selectInput(inputId = 'sep', label='sep', 
                  c("t56_44", "t78_62", "t65_51", "t91_72")
      ),
      selectInput(inputId = 'th', label='th', 
                  c("1","2","3","4","5","6","7","8","9","9.5")
      ),
      selectInput(inputId = 'sigma', label='sigma', 
                  c("0.005000",  "0.010000",  "0.020000",  "0.030000",  "0.050000",
                    "0.070000",  "0.090000",  "0.110000",  "0.130000",  "0.150000",
                    "0.170000",  "0.190000",  "0.210000",  "0.230000",
                    "0.250000",  "0.270000",  "0.290000",  "0.310000",  "0.330000",
                    "0.350000",  "0.400000",
                    "0.500000")
      ),
      selectInput(inputId = 'lambda', label='lambda', 
                  2^c(5:25)
      ),
      uiOutput('variables'),
      uiOutput('variables_diff')
      
    ),
    # Show a plot of the generated distribution
    mainPanel(
      # withMathJax(),
      plotlyOutput(outputId = "fit_plot", height = "600px"),
      plotlyOutput(outputId = "plot_reconstruction", height = "600px"),
    )
  )
)

server <- function(input, output, session) {
  output$fit_plot <- renderPlotly({
    Z<-input$Z
    
    df<-NULL
    fileA<-paste0("/home/garofalo/analysis/heavy_quarks/data_Inc/out/",input$ensemble,"_th",input$th,"_",input$sep,".dat_HLT_AoverB")
    if(file.exists(fileA)){
      mtA<-read_df(fileA)
      all_obsA<- get_all_corr(mtA)
      dfA<-get_block(string =paste0("HLT_",Z,"-sig",input$sigma,"-alpha-1.99"),
                     all_obs=all_obsA,mt=mtA,df = NULL,
                     log = FALSE, number = NULL,nudge = 0,
                     print_res = FALSE,logx=10,ix=11,
                     iy=6,ierr=7,ifit=9,ierrfit=10,
                     iplateau = 1
      )
      
      dfA<-get_block(string =paste0("HLT_",Z,"-sig",input$sigma,"-alpha0.00"),
                     all_obs=all_obsA,mt=mtA,df = dfA,
                     log = FALSE, number = NULL,nudge = 0,
                     print_res = FALSE,logx=10,ix=11,
                     iy=6,ierr=7,ifit=9,ierrfit=10,
                     iplateau = 1
      )
      dfA<-get_block(string =paste0("HLT_",Z,"-sig",input$sigma,"-alpha2.00"),
                     all_obs=all_obsA,mt=mtA,df = dfA,
                     log = FALSE, number = NULL,nudge = 0,
                     print_res = FALSE,logx=10,ix=11,
                     iy=6,ierr=7,ifit=9,ierrfit=10,
                     iplateau = 1
      )
      
      
      dfA$x<-dfA$x-log10(dfA[,12])# A/A0_ref
      dfA$xfit<-dfA$xfit-log10(dfA[,12])# A/A0_ref
      dfA$tmin<-dfA$tmin-log10(dfA[,12])# A/A0_ref
      dfA$tmax<-dfA$tmax-log10(dfA[,12])# A/A0_ref
      gg<- plot_df_corr_ggplot(dfA,width = 0.01,stroke = 1,alpha = 0.2,size_error=0.6)
      
      avep<-(max(dfA$fit)+min(dfA$fit))/2.0
      sdp<- max(dfA$fit)-min(dfA$fit)
     fig<- myplotly(gg,"","log10(A/A0)_ref", paste0("rho-th",input$th), to_print=FALSE, 
                   output = "HTML", to_webgl = FALSE,
                   yrange=c(avep-4*sdp,avep+4*sdp), legend_position = c(0.2,1))
    
    }
    else{
      gg<-myggplot()
      fig<- myplotly(gg," DATA NOT FOUD ")
    }
    return(fig)
  })
  output$plot_reconstruction <- renderPlotly({
    Z<-input$Z
    sigma<-input$sigma
    
    df<-NULL
    fileA<-paste0("/home/garofalo/analysis/heavy_quarks/data_Inc/out/",
                  input$ensemble,"_th",input$th,"_",input$sep,".dat_HLT_kernel")
    if(file.exists(fileA)){
      mtA<-read_df(fileA)
      all_obsA<- get_all_corr(mtA)
      name<-paste0("HLT_",Z,"-sig",sigma,"-alpha-1.99_lam",
                   input$lambda)
      fit<-get_full_res(string =name,
                   all_obs=all_obsA,mt=mtA)
      df1<-add_corr_to_df(string =name,
                          all_obs=all_obsA,mt=mtA,df = NULL,
                          log = FALSE, number = NULL,nudge = 0,print_res = FALSE,
                          rename = paste0("alpha-1.99 A/A0=",log10(fit[1,1]/fit[1,2]) )
      )
      name<-paste0("HLT_",Z,"-sig",sigma,"-alpha0.00_lam",
                   input$lambda)
      fit<-get_full_res(string =name,
                        all_obs=all_obsA,mt=mtA)
      df1<-add_corr_to_df(string =name,
                          all_obs=all_obsA,mt=mtA,df = df1,
                          log = FALSE, number = NULL,nudge = 0,print_res = FALSE,
                          rename = paste0("alpha0.00 A/A0=",log10(fit[1,1]/fit[1,2]) )
      )
       
      
      name<-paste0("HLT_",Z,"-sig",sigma,"-alpha2.00_lam",
                   input$lambda)
      fit<-get_full_res(string =name,
                        all_obs=all_obsA,mt=mtA)                   
      df1<-add_corr_to_df(string =name,
                          all_obs=all_obsA,mt=mtA,df = df1,
                          log = FALSE, number = NULL,nudge = 0,print_res = FALSE,
                          rename = paste0("alpha2.00 A/A0=",log10(fit[1,1]/fit[1,2]) )
      )
      
      
  
      tmp<-data.frame("omega"=df1$x, "K"=df1$y, "Kb"=df1$fit, 
                      "label"=paste0(df1$label)
      )
      df<-rbind(df,tmp)
      
      ################ plot
      
      gg<-myggplot(shape = FALSE,fill = FALSE)
      gg<-gg + geom_line(aes(x=df$omega  , y=df$Kb,
                             color=paste0("K-",df$label)
      ),linetype="dashed" )
      gg<-gg + geom_line(aes(x=df$omega  , y=df$K,
                             color="exact"), linetype="solid")
      # gg<-gg + geom_line(aes(x=df$omega  , y=df$K-df$Kb,
      #                        color=paste0("diff-",df$label)
      #                        ),linetype="dotted")
      fig<- myplotly(gg,"","omega", "K", to_print=FALSE ,
                     output = "HTML", to_webgl = FALSE,
                     legend_position = c(0.5,0.95),legend_title = "")
      return(fig)
    }
  })
  
}

# Create Shiny app ----
shinyApp(ui = ui, server = server) 
