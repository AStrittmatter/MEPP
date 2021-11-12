
# Load packages
library(shiny)
library(AER) # Used for IV regression
library(ggplot2) # Used for density plot
library(dplyr)


ui <- fluidPage(titlePanel(div("MEPP-weak instrument on Artifcial Data",style="color:#6E6E6E",align="center")) ,
         sidebarLayout(
           sidebarPanel(tags$head(
             tags$style("body {background-color: #EFFBFB; }"),
             tags$style(type='text/css', "#title1 { height: 25px; }"),
             tags$style(type='text/css', "#title2 { height: 25px; }"),
             tags$style(type="text/css", "#loadmessage { padding: 5px 0px 5px 0px; text-align: center; font-weight: bold; font-size: 100%; color: #000000; background-color: #ff8533; z-index: 105; }")
           ),
           div(align="center",h4("Parameters for DGP")),
           numericInput("N", "Select sample size", 50,step = 1),
           numericInput("beta", "Specify size of effect of D on Y", 0.5,step = 0.1),
           numericInput("bias_ols", "Specify the size of the OLS bias", 0.2,step = 0.1),
           numericInput("Fd", "Specify first stage F-statistic", 2,step = 1),
           numericInput("replications", "Select number of replications for Monte Carlo", 1000,step = 1),
           numericInput("sigma_Z","Specify noise of the DGP",1),
           div(align="center",submitButton()),
           width=3,style="background-color: #9b9b9b;padding: 4px"),
  
           mainPanel(tabsetPanel(
             tabPanel("Description", 
                 'In this ShinyApps, we use a Monte-Carlo simulation to illustrate 
                 the bias of OLS under violations of the exclusion restriction and 
                 the vulnerability of IV when the first stage is weak. 
                 Because the data is artificial, we have control over the data 
                 generating process (DGP). Accordingly, 
                 we know the true underlying data model.'),
             
              tabPanel("First Stage Regression", 'We generate one random draw from our DGP 
                       and test the power of the first stage.',
                       br(),
                       verbatimTextOutput('FirstRegression'),br(),
                       verbatimTextOutput('F_emp_text'),br(),
                       verbatimTextOutput('F_true_text')),
                       
              tabPanel("OLS Regression", br(),
                       verbatimTextOutput('OlsRegression')),
              tabPanel("IV Regression", br(),
                       verbatimTextOutput('IVRegression')),
              tabPanel("Monte Carlo", 'Up to now, we generated only one
                        random draw of the DGP and estimated the parameters of interst in this sample. 
                         But are the biases that we found systematic or did they occur just by chance? 
                         To answer this, we take many independent random draws of the DGP and estimate
                         the parameters of interest in each sample. Afterwards, we can compare the 
                         distributions of the parameters, which gives us information about 
                       the finite sample properties of the estimators.', br(),
                       plotOutput("MonteCarloPlot")),
              tabPanel("Data",  
                       br(),
                       DT::dataTableOutput("DF"))
            ))
                )
    )

server <- function(input, output){ 
  userchoise <- reactive(
    list(input$N,
         input$beta,
         input$bias_ols,
         input$Fd,
         input$replications,
         input$sigma_Z)
  )
  observeEvent(userchoise(),{
    set.seed(12345689) # Set starting values
    
    # Generate Instrument
    Z <- reactive({
      rnorm(n=isolate(input$N), mean = 0, sd =isolate(input$sigma_Z))  
    }
    )
    
    # Generate Treatment
    sigma_D <- reactive({ 
      sqrt((isolate(input$N)-2)*isolate(input$sigma_Z)^2/isolate(input$Fd)) # Specify Variance of v
    }
    )
    
    v <- reactive({  
      rnorm(n=isolate(input$N), mean = 0, sd =sigma_D())
    }
    )
    
    D <-reactive({
      Z() + v()
    }
    )
    
    # Generate Outcome
    frac  <- reactive({
      isolate(input$bias_ols)*(isolate(input$sigma_Z)^2+sigma_D()^2)/sigma_D()^2
    }
    )
    Y <- reactive({
     isolate(input$beta)*D() + frac()*v() + rnorm(n=isolate(input$N))
    }) 
    dataMEPP <- reactive(
      as.data.frame(cbind(Y(),D(),Z(),v()))
    )
    output$DF <- DT::renderDataTable({
      df <- dataMEPP()
      names(df) <- c('Y','D','Z','v')
      df 
    },
    options = list("orderClasses" = TRUE, "responsive" = TRUE, "pageLength" = 10), rownames=FALSE
    )
    
    ## First Stage Regression ##
    #OLS regression
    
    first_stage <- reactive({
        lm(D() ~ Z(),data = dataMEPP())
    }
    )
    sumFS <- reactive({
      summary(first_stage())
    })
    
    # Empirical F-statistic
    F_emp <- reactive({
    (isolate(input$N)-2)*(first_stage()$coefficient[2]^2*var(Z())/(var(D())))/(1-(first_stage()$coefficient[2]^2*var(Z())/var(D())))
    })  
    
    output$F_emp_text <- renderText({ 
      paste("Empirical F-statistic: ", round(F_emp(),2))
    })
    
    # Population F-statistic
    F_true <- reactive({
      (isolate(input$N)-2)*(isolate(input$sigma_Z)^2/(sigma_D()^2+isolate(input$sigma_Z)^2))/(1-(isolate(input$sigma_Z)^2/(sigma_D()^2+isolate(input$sigma_Z)^2)))
    }) 
    
    output$F_true_text <- renderText({ 
      paste("True F-statistic: ", round(F_true(),2))
    })
    
    output$FirstRegression <-  renderPrint({
      sumFS()
    }
    )
    
    sumFS <- reactive({
      summary(first_stage())
    })
    
    ## OLS ##
    OLS <- reactive({
      lm(Y() ~ D(),data = dataMEPP())
    }
    )
    
    sumOLS <- reactive({
      summary(OLS())
    })
    
    output$OlsRegression <-  renderPrint({
      sumOLS()
    }
    )
    
    ## IV ##
    IV <- reactive({
      ivreg(formula = Y() ~ D() | Z(),data = dataMEPP())
    }
    )
    
    sumIV <- reactive({
      summary(IV())
    })
    
    output$IVRegression <-  renderPrint({
      sumIV()
    }
    )
    
    ####### Monte Carlo Simulation #######
    
    datE <- reactive({
      # Generate matrix to store results
      df = matrix(NA, nrow = isolate(input$replications), ncol = 2)
      colnames(df) <- c("ols", "iv")
      
    
      for(r in c(1:isolate(input$replications))){
        
        ## Data Generating Process ##
        
        # Generate Instrument
        Z <- rnorm(n=isolate(input$N), mean = 0, sd =isolate(input$sigma_Z))  
        
        
        # Generate Treatment
        sigma_D <- sqrt((isolate(input$N)-2)*isolate(input$sigma_Z)^2/isolate(input$Fd)) # Specify Variance of v
        
        v <- rnorm(n=isolate(input$N), mean = 0, sd =sigma_D)
        
        D <- Z + v
      
        
        # Generate Outcome
        frac  <- isolate(input$bias_ols)*(isolate(input$sigma_Z)^2+sigma_D^2)/sigma_D^2
    
        Y <- isolate(input$beta)*D + frac*v + rnorm(n=isolate(input$N))
        
        ## Estimation ##
        
        # OLS
        ols <- lm(Y ~ D, data = as.data.frame(cbind(Y,D)))
        df[r,1] <- ols$coefficients[2]
        
        # IV
        iv <- ivreg(formula = Y ~ D | Z, data = as.data.frame(cbind(Y,D,Z)))
        df[r,2] <-iv$coefficients[2]
        
      }
      df
    })
    
    # Density plot
    
    dat <- reactive({
      data.frame(estimator = factor(rep(c("ols","iv"), each=isolate(input$replications))),
       effect = rbind(as.matrix(datE()[,1]),as.matrix(datE()[,2])))
    }
    )
    
    output$MonteCarloPlot <- renderPlot({
      ggplot(dat(), aes(x=effect, colour=estimator)) + geom_density() + xlim((isolate(input$beta)-1.5),(isolate(input$beta)+1.5)) +
        geom_vline(xintercept = isolate(input$beta), color = "green") +
        geom_vline(xintercept = (isolate(input$beta)+isolate(input$bias_ols)), color = "blue")
    }
    ) 
    })
    }

## Data Generating Process ##

shinyApp(ui = ui,server)


