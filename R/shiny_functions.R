
`%!in%` = Negate(`%in%`)


m_default_gen <- function(){ #Might add arguments to this function
  
  m_default <- matrix(nrow = 3, ncol = 2)
  colnames(m_default) <- c("Cum Prob", "Expert_1")
  #rownames(m_default) <- rep("Expert_1",3)
  m_default[,1] <- c(0.025, 0.5, 0.975)
  return(m_default)
}

return_pooled_info <- function(input_mat, St_indic = 1,dist = "best", mode =NULL){
  #dist_considered <- c("normal","t","gamma", "lognormal", "beta") 
  
  if(St_indic == 1){
    lower_bound = 0
    upper_bound = 1
  }else{
    lower_bound = -Inf
    upper_bound = Inf
  }
  
  fit.eval <- fitdist_mod(input_mat[,2:ncol(input_mat), drop = F],
                          probs = input_mat[,1], upper = upper_bound, lower = lower_bound, 
                          expertnames = paste0("Expert_",1:(ncol(input_mat)-1)),
                          mode = mode)
  
  plts_pool <- makePoolPlot(fit= fit.eval,
                            xl =lower_bound,
                            xu =upper_bound,
                            d = dist,
                            w = 1,
                            lwd =1,
                            xlab = "x",
                            ylab =expression(f[X](x)),
                            legend_full = TRUE,
                            ql = NULL,
                            qu = NULL,
                            nx = 200,
                            addquantile = FALSE,
                            fs = 12,
                            expertnames = paste0("Expert_",1:(ncol(input_mat)-1)),
                            St_indic =St_indic)
  
  dfs_pool <-  plts_pool[["data"]]
  if(dist == "best"){
    selc_fit <- fit.eval$best.fitting[,"best.fit"]
  }else{
    selc_fit <- rep(dist, length(fit.eval$best.fitting[,"best.fit"]))
  }
  selc_fit_loc <- sapply(selc_fit, function(x){which(x  == names(fit.eval$ssq))})
  
  pool.df_output <- matrix(nrow = length(selc_fit),ncol = 3)
  colnames(pool.df_output) <- c("param1", "param2", "param3")
  
  for(j in 1:length(selc_fit_loc)){
    pool.df_output[j,1:length(fit.eval[[selc_fit_loc[j]]][j,])] <-  as.numeric(as.vector(fit.eval[[selc_fit_loc[j]]][j,]))
  }
  dfs_expert <- data.frame(dist = names(selc_fit_loc), wi = 1/nrow(pool.df_output), pool.df_output)
  
  return(list(dfs_expert, plts_pool))
}


elicit_surv <- function (){
  
  #Fixes required for expertsurv implemented in this file!!
  #Fixes - fix error_mod_normal
  #Fix DLSsurvspline was accessed in the wrong order with test.deriv
  #Fix this mess in runMLE where id_trt and id_St don't work -- Fixed - pargSurv/expert wasn't reading the correct ids.
  #Fix the lik_lno function the likelihood was incorrect for the compute ICS_stan
  
  ## Load Packages
  list.of.packages <- need<-c("expertsurv", "dplyr", "shiny", "ggplot2", "ggfortify", "shinyWidgets", "survminer","shinycssloaders", "shinyjs", "shinyMatrix", "readxl") #needed libraries
  res <- lapply(list.of.packages, require, character.only = TRUE)
  not.loaded <-   list.of.packages[which(sapply(res, unlist) ==F)]
  not.loaded2 <- lapply(not.loaded, require, character.only = TRUE)
  not.installed <-   not.loaded2[which(sapply(not.loaded2, unlist) ==F)]
  #load the packages
  if(length(not.installed)) install.packages(not.loaded)
  
  options(spinner.color="#0275D8", spinner.color.background="#ffffff", spinner.size=2)
  
  
  runApp(list(ui = fluidPage(useShinyjs(),
                             # tags$style('.container-fluid {
                             #               background-color: #7b8cde;
                             #}'),
                             titlePanel("ShinyExpertsurv"),
                             sidebarPanel(
                               wellPanel(
                                 #fluidRow(column(3, downloadButton("report", "Download report")),
                                 #column(2, offset = 1, actionButton("exit", "Quit"))),
                                 fileInput('df_upload', 'Choose Excel data file to upload',
                                           accept = c(".xlsx")),
                                 p("Data should have the following columns: time and status. If your data has two treatment arms please include an arm column with numeric values indicating the treatments."),
                                 numericInput("n_expert", "Number of Experts", value = 1, min = 1),
                                 numericInput("n_timepoint", "Number of Timepoints", value = 1,min = 1,  max = 2),
                                 numericInput("scale1", "Scale for density", value = 1),
                                 numericInput("xlim", "Limit of x-axis on Kaplan-Meier curve", value = round(10,#max(df$time)*2,
                                                                                                             digits = 0)),
                                 selectInput(inputId ="pool_type_eval", label = "Pooling approach for experts", 
                                             choices = c("Linear Pool" = "linear pool",
                                                         "Logarithmic Pool"= "log pool"), 
                                             selected = "linear pool"),
                                 selectInput(inputId ="dist_select", label = "Select the best fitting distribution for Expert Pooling", 
                                             choices = c("Best Fitting" = "best",
                                                         "Normal"= "normal",
                                                         "T-distribution" = "t",
                                                         "Gamma" = "gamma",
                                                         "Log-Normal" = "lognormal",
                                                         "Beta" = "beta"), 
                                             selected = "best"),
                                 actionButton(paste0('update_expert'), "Plot/Update Expert Opinions")
                                 
                               ),
                               
                               hr(),
                               
                               tabsetPanel(id = "Timepoints",
                                           tabPanel("Timepoints1",
                                                    numericInput(paste0("time1"), label= "Timepoint", value= 1),
                                                    textInput('quant_vec1', 'Enter a Vector of Quantiles', "0.025,0.5,0.975"),
                                                    matrixInput(
                                                      inputId = "matrix1",
                                                      value = m_default_gen(),
                                                      class = "numeric",
                                                      cols = list(names = TRUE,
                                                                  editableNames = FALSE),
                                                      rows = list(names = FALSE,
                                                                  editableNames = FALSE)),
                                                    helpText("Enter the judgements in the table below,
                            one column per expert. Enter quantile values 
                        corresponding to the cumulative probabilities."),
                                                    
                                                    plotOutput(paste0("expert_plot1"))),
                                           
                                           tabPanel("Timepoints2",
                                                    numericInput(paste0("time2"), label= "Timepoint", value= 1),
                                                    textInput('quant_vec2', 'Enter a Vector of Quantiles', "0.025,0.5,0.975"),
                                                    matrixInput(
                                                      inputId = "matrix2",
                                                      value = m_default_gen(),
                                                      class = "numeric",
                                                      cols = list(names = TRUE,
                                                                  editableNames = FALSE),
                                                      rows = list(names = FALSE,
                                                                  editableNames = FALSE)),
                                                    helpText("Enter the judgements in the table below,
                            one column per expert. Enter quantile values 
                        corresponding to the cumulative probabilities."),
                                                    plotOutput(paste0("expert_plot2")))
                               )),
                             mainPanel(
                               #withSpinner(tableOutput('tb'), type = 2),
                               h3("Kaplan-Meier Survival Plot"),
                               plotOutput(paste0("plot_km_expert1")),
                               selectInput("opinion_type", label = "Choose opinion type", 
                                           choices = c("Survival at timepoint(s)" = "survival",
                                                       "Mean difference between survival"= "mean",
                                                       "No expert opinion" = "no_expert"), 
                                           selected = "survival"),
                               
                               selectInput("stat_type", label = "Choose statistical approach", 
                                           choices = c("Frequentist" = "mle","Bayesian" = "hmc"), 
                                           selected = "mle"),
                               numericInput("id_trt", label = "Select treatment ID corresponding to expert opinion", value = 1),
                               
                               
                               pickerInput(
                                 inputId = "param_mod", 
                                 label = "Choose models:", 
                                 choices = c("Exponential" = "exp",
                                             "Weibull" = "wei",
                                             "Gompertz" = "gomp",
                                             "Log-Logistic"= "llo",
                                             "Log-normal" = "lno",
                                             "Generalized-Gamma" = "gga",
                                             "Royston-Parmar" = "rps"), 
                                 options = list(
                                   `actions-box` = TRUE, 
                                   size = 10,
                                   `selected-text-format` = "count > 3"
                                 ), 
                                 multiple = TRUE,
                                 selected  = c("exp", "wei")
                               ),
                               selectInput("incl_psa", label = "Include Statistical Uncertainty in Plots", 
                                           choices = c("Yes" = "yes",
                                                       "No"= "no"), 
                                           selected = "no"),  
                               
                               actionButton("run_analysis", "Run Analysis"),
                               #br(),
                               #h4("Output of the Statistical Analysis"),
                               selectInput("gof_type", label = "Choose goodness of fit measure", 
                                           choices = c("AIC" = "aic","BIC" = "bic"), 
                                           selected = "AIC"),
                               plotOutput("plot_gof"),
                               # selectInput("outFormat", label = "Report format",
                               #             choices = list('html' = "html_document",
                               #                            'pdf' = "pdf_document",
                               #                            'Word' = "word_document"))    #textOutput("txt"),
                               
                               textInput('file_name', 'Enter a File name to save output', "Output-File"),
                               actionButton("save_output", "Save current files")
                             )
                             
  ),
  server = function(input, output, session) {
    
    value <- reactiveValues(
      m_default = m_default_gen(),
      n_expert_prev = 1,
      quant_vec2 = NULL,
      id_trt = NULL) #Up to max timepoints
    
    
    observeEvent(input$stat_type,{
      if(input$stat_type == "mle"){
        updateSelectInput(session,"gof_type",choices = c("AIC" = "aic", "BIC" = "bic"))
      }else{
        updateSelectInput(session,"gof_type",choices = c("WAIC" = "waic", "PML" = "pml"))
      }
      
    })
    
    observeEvent(input$df_upload,{
      inFile <- input$df_upload
      df_upload <- readxl::read_excel(inFile$datapath)
      
      if(is.null(df_upload[["arm"]])){
        #browser()
        result.km <- survfit(Surv(time, status) ~ 1, data = df_upload, conf.type="log-log")
        km.data <- data.frame(cbind(result.km[[c("time")]],
                                    result.km[[c("surv")]],
                                    result.km[[c("upper")]],
                                    result.km[[c("lower")]],
                                    arm = 1))
        
        updateSelectInput(session,"opinion_type",choices = c("Survival at timepoint(s)" = "survival",
                                                             "No expert opinion" = "no_expert"))
        hide("id_trt") #hide id_trt panel
        value$id_trt <- NULL
        
        df_upload$arm <- 1
      }else{
        #browser()
        show("id_trt") #hide id_trt panel
        updateNumericInput(inputId = "id_trt", min =min(df_upload$arm),max = max(df_upload$arm), value = max(df_upload$arm))
        km.data <- NULL
        for(i in unique(df_upload$arm)){
          df_temp <- df_upload %>% filter(arm == i)
          result.km_temp <- survfit(Surv(time, status) ~ 1, data = df_temp, conf.type="log-log")
          km.data_temp <- data.frame(cbind(result.km_temp[[c("time")]],
                                           result.km_temp[[c("surv")]],
                                           result.km_temp[[c("upper")]],
                                           result.km_temp[[c("lower")]],
                                           arm = i))
          
          km.data <-  rbind(km.data,km.data_temp)
        }
        
        updateSelectInput(session,"opinion_type",
                          choices = c("Survival at timepoint(s)" = "survival",
                                      "Mean difference between survival"= "mean",
                                      "No expert opinion" = "no_expert"),
                          selected = "survival")
        
      }
      
      colnames(km.data) <- c("Time", "Survival", "upper", "lower", "arm")
      
      value$km.data <- km.data
      value$df_upload <- df_upload
      value$id_trt <- input$id_trt
      #Need to adjust for arm
      
      plot_fit <- ggplot(value$km.data, aes(x = Time,y =Survival, col = factor(arm)))+
        geom_step()+
        ylim(0,1)+
        xlim(0, input$xlim)+
        geom_step(aes(x  = Time, y =upper, col = factor(arm)))+
        geom_step(aes(x  = Time, y =lower, col = factor(arm)))+
        theme_light()#+
      #scale_x_continuous(expand = c(0, 0))+#, breaks=seq(0, 30, 2)) + 
      #scale_y_continuous(expand = c(0, 0))#, breaks=seq(0, 1, 0.05))
      
      
      
      output$plot_km_expert1<- renderPlot(plot_fit)
      
    })
    
    observeEvent(input$n_timepoint, {
      
      if(input$n_timepoint > 1){
        showTab(inputId = "Timepoints", target = "Timepoints1")
        showTab(inputId = "Timepoints", target = "Timepoints2")
      }
      if(input$n_timepoint == 1){
        showTab(inputId = "Timepoints", target = "Timepoints1")
        hideTab(inputId = "Timepoints", target = "Timepoints2")
      }
    })
    
    observeEvent(input$opinion_type,{
      if(input$opinion_type == "survival"){
        
        if(input$n_timepoint > 1){
          showTab(inputId = "Timepoints", target = "Timepoints1")
          showTab(inputId = "Timepoints", target = "Timepoints2")
        }
        if(input$n_timepoint == 1){
          showTab(inputId = "Timepoints", target = "Timepoints1")
          hideTab(inputId = "Timepoints", target = "Timepoints2")
        }
        
      }
      
      if(input$opinion_type == "mean"){
        hideTab(inputId = "Timepoints", target = "Timepoints2")
      }
      
      if(input$opinion_type == "no_expert"){
        hideTab(inputId = "Timepoints", target = "Timepoints1")
        hideTab(inputId = "Timepoints", target = "Timepoints2")
      }
      
      if(input$opinion_type == "survival" | input$opinion_type == "mean"){
        show("n_expert")
        show("n_timepoint")
        show("scale1")
        show("pool_type_eval")
        show("dist_select")
        if(is.null(value$id_trt)){
          hide("id_trt")
        }else{
          hide("id_trt")
        }
        
        
      }else{ #No Expert Opinion
        hide("n_expert")
        hide("n_timepoint")
        hide("scale1")
        hide("pool_type_eval")
        hide("dist_select")
        hide("id_trt")
      }
      
    })
    
    observeEvent(input$id_trt,{
      value$id_trt <- input$id_trt
    })
    
    observeEvent(input$n_expert, {
      # browser()
      if(input$n_expert == 1){
        shinyjs::hideElement(id = "pool_type_eval")
        hide("pool_type_eval")
        hide("dist_select")
      }else{
        shinyjs::showElement(id = "pool_type_eval")
        show("pool_type_eval")
        show("dist_select")
      }
      
      for(i in 1:2){ #Modify this force it to me 2 which is the max number of timepoints
        mat_exist <- input[[paste0("matrix",i)]]
        if(input$n_expert > value$n_expert_prev){
          extra_cols <- input$n_expert - value$n_expert_prev 
          mat_bind <- matrix(nrow = nrow(mat_exist), ncol = extra_cols)
          mat_exist <- cbind(mat_exist,mat_bind)
        }else if(input$n_expert == value$n_expert_prev){
        } else{
          mat_exist <- mat_exist[,1:(input$n_expert+1),drop = F]
        }
        colnames(mat_exist) <- c("Cum Prob", paste0("Expert_", 1:input$n_expert))
        updateMatrixInput(session, paste0("matrix",i), value=mat_exist )
      }
      value$n_expert_prev <- input$n_expert
      
    })
    
    
    toListen <- reactive({
      list(input$quant_vec1,input$quant_vec2)
    })
    
    observeEvent(toListen(),{
      
      for(i in 1:input$n_timepoint){#max number of quant_vec
        #browser()
        if(!is.null(input[[paste0("quant_vec",i)]])){
          quant_vec_temp <- input[[paste0("quant_vec",i)]]
          quant_num <- as.numeric(unlist(strsplit(quant_vec_temp,",")))
          if(length(quant_num)==0){
            new_mat <-   matrix(ncol = input$n_expert +1, nrow = 1) # Handle case when nothing is entered
            colnames(new_mat) <- c("Cum Prob", paste0("Expert_",1:input$n_expert ))
            
          }else{
            mat_exist <- input[[paste0("matrix",i)]]
            retain_quant_index <-which(mat_exist[,1] %in% quant_num)
            retain_quant <- mat_exist[retain_quant_index,1]
            change_quant_index <- which(quant_num %!in% retain_quant)
            new_mat <- matrix(ncol = input$n_expert +1, nrow = length(quant_num))
            colnames(new_mat) <- c("Cum Prob", paste0("Expert_",1:input$n_expert ))
            
            if(length(retain_quant_index)>0){
              new_mat[1:length(retain_quant_index),] <-mat_exist[retain_quant_index,]
              
            }
            if(length(change_quant_index)>0){
              new_mat[(length(retain_quant_index)+1):nrow(new_mat),1] <- quant_num[change_quant_index]
            }
          }
          updateMatrixInput(session, paste0("matrix",i), value=new_mat)
          
        }
      }})
    
    
    observeEvent(input$update_expert,{
      #pool_type_eval <- "linear pool"
      times_expert_vec <- c()
      df.linear_all <- NULL
      param_expert <- list()
      
      
      plot_fit <- ggplot(value$km.data, aes(x = Time,y =Survival, col = factor(arm)))+
        geom_step()+
        ylim(0,1)+
        xlim(0, input$xlim)+
        geom_step(aes(x  = Time, y =upper, col = factor(arm)))+
        geom_step(aes(x  = Time, y =lower, col = factor(arm)))+
        theme_light()#+
      #scale_x_continuous(expand = c(0, 0))+#, breaks=seq(0, 30, 2)) + 
      #scale_y_continuous(expand = c(0, 0))
      
      if(!any(is.na(input[["matrix1"]][,2]))){ # If Expert opinions are not NA values
        
        
        for(i in 1:input$n_timepoint){ #Update
          
          #want to allow for zeros
          output_pool <- return_pooled_info(input[[paste0("matrix",i)]], St_indic = 0,dist = input$dist_select, mode =NULL)
          if(input$opinion_type == "survival"){
            output_pool[[2]] <- output_pool[[2]]+ xlim(c(0,1)) #If survival we want to truncate.
          }
          
          
          output[[paste0("expert_plot",i)]] <- renderPlot(output_pool[[2]])
          
          times_expert = input[[paste0("time",i)]]
          times_expert_vec <- c(times_expert_vec, times_expert)
          
          df.linear <- subset(output_pool[[2]]$data, ftype == input$pool_type_eval) %>% rename(y = x) %>% 
            mutate(x = times_expert + fx*input[[paste0("scale1")]], 
                   times_expert = times_expert)
          df.linear_all <- rbind(df.linear_all, df.linear)
          
          output_pool[[1]][,"dist"] <- gsub("normal", "norm", output_pool[[1]][,"dist"])
          
          param_expert[[i]] <- output_pool[[1]]
        }
        
        value$param_expert <- param_expert
        value$timepoint_expert <- times_expert_vec
        value$df.linear_all <- df.linear_all
        
        if(input$opinion_type == "survival"){
          plot_fit <- plot_fit+
            geom_ribbon(data = df.linear_all, aes(x = x, y = y, xmin= x, xmax =times_expert, group=times_expert), 
                        fill = "sky blue", alpha = 0.5, colour = "grey")
        }
        
        
      }
      
      output$plot_km_expert1<- renderPlot(plot_fit)
      
      
      
    })
    
    observeEvent(input$run_analysis, {
      
      #browser()
      if(length(unique(value$df_upload[["arm"]]))==1){
        formula_text <- "Surv(time,status)~1"
      }else{
        formula_text <- "Surv(time,status)~factor(arm)"
      }
      
      if(!is.null(value$param_expert)& input$opinion_type != "no_expert"){
        
        mod_fit  <- fit.models.expert(formula=as.formula(formula_text),data=value$df_upload,
                                      distr=input$param_mod,
                                      method=input$stat_type,
                                      pool_type = input$pool_type_eval,#"log pool", 
                                      opinion_type = input$opinion_type,
                                      times_expert = value$timepoint_expert, 
                                      param_expert = value$param_expert,
                                      id_trt = input$id_trt,
                                      id_St  = input$id_trt,
                                      k = 1)
        
        value$mod_fit <- mod_fit
      }
      
      if(input$opinion_type == "no_expert"){
        
        
        param_expert_vague <- list()
        param_expert_vague[[1]] <- data.frame(dist = "beta", wi = 1, param1 = 1, param2 = 1, param2 = NA)
        
        mod_fit  <- fit.models.expert(formula=as.formula(formula_text),data=value$df_upload,
                                      distr=input$param_mod,
                                      method=input$stat_type,
                                      pool_type = input$pool_type_eval,#"log pool", 
                                      opinion_type = "survival",
                                      times_expert = 2, 
                                      param_expert = param_expert_vague,
                                      k = 1,
                                      id_St  = 1)
        
        value$mod_fit <- mod_fit
      }
      
    })
    
    value$plot_km_expert1 <- eventReactive(value$mod_fit,{
      
      if(input$incl_psa == "yes"){
        #browser()
        models <- names(value$mod_fit$models)
        psa_outuput <- list()
        
        for(i in 1:length(value$mod_fit$models)){
          psa <- make.surv(fit = value$mod_fit,mod = i, nsim = 1000, t = seq(0,input$xlim,length.out = 1000))
          df_temp  <- t(apply(psa$mat[[1]], 1,quantile, probs = c(0.025, 0.5,.975))) %>% data.frame()
          df_temp$time <- seq(0,input$xlim,length.out = 1000)
          mod_name <- names(value$mod_fit$models)[i]
          psa_outuput[[mod_name]] <- df_temp %>% mutate(model = mod_name)
        }
        
        df_final_plot <- do.call(rbind.data.frame, psa_outuput)
        df_final_plot$models <- factor(df_final_plot$model, levels = unique(df_final_plot$model))
        #browser()
        if(input$opinion_type == "survival"){
          ggplot(data = df_final_plot, aes(y = X50., x = time, group = models, colour = models))+
            geom_line()+
            geom_line(data = df_final_plot ,aes(y = X97.5., x = time),linetype="dotdash")+
            geom_line(data = df_final_plot ,aes(y = X2.5., x = time), linetype="dotdash")+
            geom_step(data =value$km.data, mapping = aes(x = Time,y =Survival, col = factor(arm)),inherit.aes = FALSE)+
            geom_ribbon(data = value$df.linear_all, aes(x = x, y = y, xmin= x, xmax =times_expert, group=times_expert), 
                        fill = "sky blue", alpha = 0.5, colour = "grey")
        }else{
          ggplot(data = df_final_plot, aes(y = X50., x = time, group = models, colour = models))+
            geom_line()+
            geom_line(data = df_final_plot ,aes(y = X97.5., x = time),linetype="dotdash")+
            geom_line(data = df_final_plot ,aes(y = X2.5., x = time), linetype="dotdash")+
            geom_step(value$km.data, aes(x = Time,y =Survival, col = factor(arm)),inherit.aes = FALSE)+
            geom_step(value$km.data, aes(x = Time,y =lower, col = factor(arm)),inherit.aes = FALSE)+
            geom_step(value$km.data, aes(x = Time,y =upper, col = factor(arm)),inherit.aes = FALSE)
        }
        
      }else{
        
        if(input$opinion_type == "survival"){
          plot(value$mod_fit, add.km = TRUE,t = seq(0,input$xlim,length.out = 1000))+
            geom_ribbon(data = value$df.linear_all, aes(x = x, y = y, xmin= x, xmax =times_expert, group=times_expert), 
                        fill = "sky blue", alpha = 0.5, colour = "grey")
          
        }else{
          plot(value$mod_fit, add.km = TRUE,t = seq(0,input$xlim,length.out = 1000))
          
        } 
        
      }  
    })
    
    value$plot_gof <- eventReactive(value$mod_fit,{
      #browser()
      model.fit.plot(value$mod_fit,type = input$gof_type)
    })
    
    observeEvent(input$update_expert, {
      # browser()
      hide("plot_gof")
    })
    
    
    observeEvent(input$run_analysis, {
      
      output$plot_km_expert1 <- renderPlot(
        value$plot_km_expert1())
      
      output$plot_gof <- renderPlot(
        value$plot_gof())
      
      show("plot_gof")
      
    })
    
    observeEvent(input$save_output,{
      #browser()
      saveRDS(list(model = value$mod_fit, surv_plt = value$plot_km_expert1(), gof_plt = value$plot_gof()),
              file = paste0(input$file_name,".rds"))
      
      #readRDS(file = paste0(input$file_name,".rds"))
    })
    
    
  }),launch.browser = TRUE)
  
  
}
