# load packages
library(tidyverse)
library(haven)
library(survey)
library(broom)
library(ggplot2)
library(purrr)
library(pROC)
library(car) 
library(survey)
library(surveyROC)
library(svyROC)
library(srvyr)
library(vroom)
library(caret)
library(misty)
library(mvnmle)
library(naniar)
library(gridExtra)
library(mice)
library(effsize)
library(kableExtra)
library(DiagrammeR)
library(DiagrammeRsvg)
library(webshot2)
library(rsvg)
library(khsmisc)
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())


# general functions
# Function to format p-values
format_p <- function(p) {
  ifelse(p < 0.001, "<0.001", round(p, digits = 3))
}

# function to print and examine ROC curve, AUC, and other metrics
model_eval_fun <- function(formula = "HIPR ~ HIC15_log", data = madymo_nhtsa_df_listwise) {
  # Extract required components from survey design
  design <- svydesign(
    id = ~PSU,
    strata = ~PSUSTRAT,
    weights = ~RATWGT,
    data = data
  )
  
  model <- svyglm(
    as.formula(formula),
    design = design,
    family = quasibinomial
  )
  
  roc_data <- data.frame(
    y = design$variables$HIPR,
    phat = fitted(model),
    weights = design$variables$RATWGT
  )
  
  
  # Create weighted ROC object
  invisible(mycurve <- wroc(
    response.var = "y",
    phat.var = "phat",
    weights.var = "weights",
    data = roc_data,
    tag.event = levels(roc_data$y)[2],
    tag.nonevent = levels(roc_data$y)[1],
    cutoff.method = "Youden"
  ))
  
  # Calculate prevalence of the event
  event_level <- levels(roc_data$y)[2]
  prevalence <- sum(roc_data$weights[roc_data$y == event_level]) / sum(roc_data$weights)
  
  # browser()
  
  # Extract performance metrics
  metrics_list <- list(
    AUC = mycurve$wauc,
    Optimal_Cutoff = mycurve$optimal.cutoff$cutoff.value,
    Sensitivity = mycurve$optimal.cutoff$Sew,
    Specificity = mycurve$optimal.cutoff$Spw,
    Correlation = cor(as.numeric(roc_data$y), roc_data$phat, use = "pairwise.complete.obs"),
    # get p-value and standard error for model (only HIC15_log)
    beta = summary(model)$coefficients[2, 1],
    p_value = summary(model)$coefficients[2, 4],
    se = summary(model)$coefficients[2, 2],
    Cases = mycurve$basics$n.event,
    Controls = mycurve$basics$n.nonevent
    # PPV - excluding because low prevalence makes PPV low
    # PPV = 100 * ((mycurve$optimal.cutoff$Sew * prevalence) / 
    #                ((mycurve$optimal.cutoff$Sew * prevalence) + 
    #                   ((1 - mycurve$optimal.cutoff$Spw) * (1 - prevalence))))
  )
  
  
  # Generate plot
  roc_plot <- wroc.plot(
    x = mycurve,
    print.auc = TRUE,
    print.cutoff = TRUE
  )
  
  # Return both plot and metrics
  return(list(
    mycurve = mycurve,
    metrics = metrics_list,
    model = model
  ))
}

# function to compare models in dataframe using loop
compare_df_fun <- function(models) {
  comparison_df <- data.frame()
  for (model in models) {
    # Extract the AUC and other metrics
    # β <- get(model)$model$coefficients
    auc <- get(model)$metrics$AUC
    # optimal_cutoff <- get(model)$metrics$Optimal_Cutoff
    sensitivity <- get(model)$metrics$Sensitivity
    specificity <- get(model)$metrics$Specificity
    # correlation <- get(model)$metrics$Correlation
    # get n.event from mycurve$basics 
    # num_cases <- get(model)$metrics$Cases
    # count number of HIPR == 0
    # num_controls <- get(model)$metrics$Controls
    p_value <- get(model)$metrics$p_value
    se <- get(model)$metrics$se
    beta <- get(model)$metrics$beta
    # PPV <- get(model)$metrics$PPV
    
    # browser()
    
    # Create a new row for the comparison dataframe
    new_row <- data.frame(
      Model = model,
      
      AUC = auc,
      # Optimal_Cutoff = optimal_cutoff,
      # how to interpret results as cutoff of injury for HIC15
      # Optimal_Cutoff_log = (optimal_cutoff),
      Sensitivity = sensitivity,
      Specificity = specificity,
      # Correlation = correlation,
      # Cases = num_cases,
      # Controls = num_controls,
      HIC_Beta = beta,
      HIC_SE = se,
      P_value = format_p(p_value)
    )
    
    # Append the new row to the comparison dataframe
    comparison_df <- rbind(comparison_df, new_row)
  }
  return(comparison_df)
}

roc_grid_prep_fun <- function(models_obj) {
  
  model_names <- names(models_obj)
  
  # Extract data from each model object
  roc_df <- imap_dfr(models_obj, function(model, model_name) {
    # Extract ROC curve info from wroc object
    wroc_data <- model$mycurve$wroc.curve
    
    # Extract correlation and other values for display on side
    beta <- model$metrics$beta
    # r_value <- model$metrics$Correlation
    AUC <- model$metrics$AUC
    SE <- model$metrics$se
    P_value <- model$metrics$p_value
    cases <- model$metrics$Cases
    Controls <- model$metrics$Controls
    
    
    tibble(
      FPR = 1 - wroc_data$Spw.values,  # 1 - Specificity
      TPR = wroc_data$Sew.values,      # Sensitivity
      Model = model_name,
      AUC = AUC,
      # Optimal_Cutoff = optimal_cutoff,
      Beta = beta,
      SE = SE,
      P_value = P_value,
      Cases = cases,
      Controls = Controls
      # Correlation = r_value
    )
  })
  
  
  # preserve order of models (don't get alphabetized by summarise
  roc_df <- roc_df %>%
    mutate(Model = factor(Model, levels = model_names))
}


roc_compare_grid_fun <- function(models_obj, col_num = 3) {
  # Generate ROC grid with adjusted columns
  
  # browser()
  
  roc_df <- models_obj
  
  
  # Get unique correlation per model for annotation, next line p-value
  label_df <- roc_df %>%
    group_by(Model) %>%
    # summarize to get unique values for Correlation and AUC
    summarise(
      # Correlation = first(Correlation),
      Cases = first(Cases),
      Beta = first(Beta),
      AUC = first(AUC),
      # Optimal_Cutoff = first(Optimal_Cutoff),
      SE = first(SE),
      P_value = first(P_value),
      Cases = first(Cases),
      Controls = first(Controls)
    ) %>%
    mutate(label = paste0("Cases = ", Cases,
                          "\nControls = ", Controls, 
                          "\nAUC = ", round(AUC, 2),
                          "\nBeta = ", round(Beta, 2),
                          "\nSE = ", round(SE, 2), 
                          "\nP-value = ", format_p(P_value)
                          # "\nOptimal Cutoff = ", round(Optimal_Cutoff, 2),
                          # "\n\nr = ", round(Correlation, 2)
    )
    )
  
  # Plot curves
  ggplot(roc_df, aes(x = FPR, y = TPR)) +
    geom_line(color = "blue", linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
    facet_wrap(~ Model, ncol = col_num) +
    geom_text(
      data = label_df,
      aes(x = 0.65, y = 0.25, label = label),  # Adjust position if needed
      inherit.aes = FALSE,
      size = 3,
      hjust = 0
    ) +
    theme_minimal() +
    labs(
      title = "ROC Curves by Model",
      x = paste0("1 - Specificity (False Positive Rate)"),
      y = paste0("Sensitivity (True Positive Rate)")
    )
}

# creates grid of ROC curves that are split by quantiles (defined in argument) for 2 continuous variables to compare.
# var_x is variable that appears on x-axis, var y is variable on y-axis
roc_grid_interaction_fun <- function(
    df = madymo_nhtsa_df_listwise, 
    var_x = DVEST_MPH, 
    var_y = BMI, 
    var_x_lab = "Delta V (MPH)",
    var_y_lab = "BMI (kg/m^2)",
    num_levelsx = 3, 
    num_levelsy = 3, 
    outcome_var = "HIPR", 
    predictor = "HIC15_log",
    split_x_by_quantile = TRUE,
    split_x_breaks = NULL,
    split_y_by_quantile = TRUE,
    split_y_breaks = NULL) {
  
  # Capture variable names
  var_x_name <- deparse(substitute(var_x))
  var_y_name <- deparse(substitute(var_y))
  
  var_x_quantile_name <- paste0(var_x_name, "_quantile")
  var_y_quantile_name <- paste0(var_y_name, "_quantile")
  
  # Validate custom breaks if required
  if (!split_x_by_quantile) {
    if (is.null(split_x_breaks)) {
      stop("split_x_breaks must be provided when split_x_by_quantile is FALSE")
    }
    if (!is.numeric(split_x_breaks) || length(split_x_breaks) < 2) {
      stop("split_x_breaks must be a numeric vector of length >= 2")
    }
  }
  
  if (!split_y_by_quantile) {
    if (is.null(split_y_breaks)) {
      stop("split_y_breaks must be provided when split_y_by_quantile is FALSE")
    }
    if (!is.numeric(split_y_breaks) || length(split_y_breaks) < 2) {
      stop("split_y_breaks must be a numeric vector of length >= 2")
    }
  }
  
  # Create temporary df with quantiles or custom splits
  temp_df <- df
  
  # Handle var_x splitting
  if (split_x_by_quantile) {
    temp_df <- temp_df %>%
      mutate(!!var_x_quantile_name := as.factor(ntile({{var_x}}, num_levelsx)))
  } else {
    temp_df <- temp_df %>%
      mutate(
        !!var_x_quantile_name := as.numeric(cut({{var_x}}, 
                                                breaks = split_x_breaks, 
                                                include.lowest = TRUE, 
                                                right = TRUE))
      ) %>%
      mutate(!!var_x_quantile_name := as.factor(!!sym(var_x_quantile_name)))
  }
  
  # Handle var_y splitting
  if (split_y_by_quantile) {
    temp_df <- temp_df %>%
      mutate(!!var_y_quantile_name := as.factor(ntile({{var_y}}, num_levelsy)))
  } else {
    temp_df <- temp_df %>%
      mutate(
        !!var_y_quantile_name := as.numeric(cut({{var_y}}, 
                                                breaks = split_y_breaks, 
                                                include.lowest = TRUE, 
                                                right = TRUE))
      ) %>%
      mutate(!!var_y_quantile_name := as.factor(!!sym(var_y_quantile_name)))
  }
  
  # Determine number of levels for each variable
  n_level_x <- if (split_x_by_quantile) num_levelsx else (length(split_x_breaks) - 1)
  n_level_y <- if (split_y_by_quantile) num_levelsy else (length(split_y_breaks) - 1)
  
  # Get min/max for each quantile or split
  # For var_x
  if (split_x_by_quantile) {
    var_x_min_max <- temp_df %>%
      group_by(!!sym(var_x_quantile_name)) %>%
      summarise(
        min = min(round({{var_x}}), na.rm = TRUE),
        max = max(round({{var_x}}), na.rm = TRUE)
      )
  } else {
    n_x <- length(split_x_breaks) - 1
    var_x_min_max <- data.frame(
      level = 1:n_x,
      min = round(split_x_breaks[1:n_x]),
      max = round(split_x_breaks[2:(n_x + 1)])
    )
    colnames(var_x_min_max)[1] <- var_x_quantile_name
  }
  
  # For var_y
  if (split_y_by_quantile) {
    var_y_min_max <- temp_df %>%
      group_by(!!sym(var_y_quantile_name)) %>%
      summarise(
        min = min(round({{var_y}}), na.rm = TRUE),
        max = max(round({{var_y}}), na.rm = TRUE)
      )
  } else {
    n_y <- length(split_y_breaks) - 1
    var_y_min_max <- data.frame(
      level = 1:n_y,
      min = round(split_y_breaks[1:n_y]),
      max = round(split_y_breaks[2:(n_y + 1)])
    )
    colnames(var_y_min_max)[1] <- var_y_quantile_name
  }
  
  # Create all combinations of subset rules
  subset_rules <- list()
  
  # 1. Individual cells (var_x levels × var_y levels)
  for (j in 1:n_level_y) {
    for (i in 1:n_level_x) {
      name <- paste0(
        var_x_lab, " ", var_x_min_max$min[i], "-", var_x_min_max$max[i], 
        ", ", var_y_lab, " ", var_y_min_max$min[j], "-", var_y_min_max$max[j]
      )
      rule <- paste0(var_x_quantile_name, " == ", i, " & ", var_y_quantile_name, " == ", j)
      subset_rules[[name]] <- rule
    }
    
    # 2. Row summaries (all var_x, fixed var_y level)
    name <- paste0(
      "All ", var_x_lab, 
      ", ", var_y_lab, " ", var_y_min_max$min[j], "-", var_y_min_max$max[j]
    )
    rule <- paste0(var_y_quantile_name, " == ", j)
    subset_rules[[name]] <- rule
  }
  
  # 3. Column summaries (fixed var_x level, all var_y)
  for (i in 1:n_level_x) {
    name <- paste0(
      var_x_lab, " ", var_x_min_max$min[i], "-", var_x_min_max$max[i],
      ", All ", var_y_lab
    )
    rule <- paste0(var_x_quantile_name, " == ", i)
    subset_rules[[name]] <- rule
  }
  
  # 4. Full model (all var_x, all var_y)
  name <- paste0("All ", var_x_lab, ", All ", var_y_lab)
  rule <- paste0(var_x_quantile_name, " %in% 1:", n_level_x, " & ", 
                 var_y_quantile_name, " %in% 1:", n_level_y)
  subset_rules[[name]] <- rule
  
  # Create models for each subset
  models <- lapply(names(subset_rules), function(name) {
    rule <- subset_rules[[name]]
    formula <- paste0(outcome_var, " ~ ", predictor)
    model_eval_fun(
      formula = formula,
      data = subset(temp_df, eval(parse(text = rule)))
    )
  })
  names(models) <- names(subset_rules)
  
  # Generate ROC grid with adjusted columns
  col_num = n_level_x + 1
  roc_df <- roc_grid_prep_fun(models)
  
  # Get unique correlation per model for annotation, next line p-value
  label_df <- roc_df %>%
    group_by(Model) %>%
    summarise(
      Cases = first(Cases),
      Beta = first(Beta),
      AUC = first(AUC),
      SE = first(SE),
      P_value = first(P_value),
      Cases = first(Cases),
      Controls = first(Controls)
    ) %>%
    mutate(label = paste0("Cases = ", Cases,
                          "\nControls = ", Controls, 
                          "\nAUC = ", round(AUC, 2),
                          "\nBeta = ", round(Beta, 2),
                          "\nSE = ", round(SE, 2), 
                          "\nP-value = ", format_p(P_value)))
  
  # Return plot components
  return(list(
    roc_df = roc_df,
    label_df = label_df,
    col_num = col_num,
    var_x_lab = var_x_lab,
    var_y_lab = var_y_lab
  ))
}

roc_grid_plot <- function(roc_data) {
  ggplot(roc_data$roc_df, aes(x = FPR, y = TPR)) +
    geom_line(color = "blue", linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
    facet_wrap(~ Model, ncol = roc_data$col_num) +
    geom_text(
      data = roc_data$label_df,
      aes(x = 0.65, y = 0.3, label = label),
      inherit.aes = FALSE,
      size = 1.5,
      hjust = 0
    ) +
    theme_minimal() +
    labs(
      title = "ROC Curves by Model",
      x = paste0(roc_data$var_y_lab, "\n1 - Specificity (False Positive Rate)"),
      y = paste0(roc_data$var_x_lab, "\nSensitivity (True Positive Rate)")
    ) +
    theme(text=element_text(size=6)
          )
}


# for continuous variables, t.test
# for categorical variables, chisq.test (if less than 5, fisher.test)
baseline_desc_fun <- function(df) {
  # Helper function to format p-values
  format_p <- function(p) {
    if (is.na(p)) {
      return("")
    } else if (p < 0.001) {
      return("<0.001")
    } else {
      return(sprintf("%.3f", p))
    }
  }
  
  # Create cases and controls data frames
  cases_df <- df %>% filter(HIPR == 1)
  controls_df <- df %>% filter(HIPR == 0)
  
  # Define variables with labels and types
  baseline_vars <- list(
    list(var_name = "DV_MPH", label = "Delta V (mph)", type = "continuous"),
    list(var_name = "SEATPOS", label = "Seat Position", type = "categorical"),
    list(var_name = "BELTUSE_BIN", label = "Belt Use (%)", type = "categorical"),
    list(var_name = "BAGDEPLOY_BIN", label = "Airbag Deployment (%)", type = "categorical"),
    list(var_name = "AGE", label = "Age", type = "continuous"),
    list(var_name = "BMI", label = "BMI in kg/m²", type = "continuous"),
    list(var_name = "SEX_BIN", label = "Sex (%)", type = "categorical"),
    list(var_name = "ALCINV", label = "Alcohol Involvement (%)", type = "categorical"),
    list(var_name = "WEATHER", label = "Weather Conditions (%)", type = "categorical"),
    list(var_name = "VEHWGT", label = "Vehicle Weight (lbs)", type = "continuous"),
    list(var_name = "MODELYR", label = "Vehicle Model Year", type = "continuous"),
    list(var_name = "TIME_cat", label = "Time of Day", type = "categorical"),
    list(var_name = "VEHTYPE", label = "Vehicle Type", type = "categorical")
  )
  
  # Initialize empty dataframe for results
  baseline_descriptives <- data.frame(
    Characteristic = character(),
    Value = character(),
    Control = character(),
    Case = character(),
    OR = character(),
    P_Value = character(),
    stringsAsFactors = FALSE
  )
  
  # Process each variable
  for (var in baseline_vars) {
    if (var$type == "continuous") {
      # Continuous variables: Logistic regression for OR and p-value
      formula <- as.formula(paste("HIPR ~", var$var_name))
      
      # create survey design object
      svy_design_temp <- svydesign(
        ids = ~ PSU,
        strata = ~ PSUSTRAT,
        weights = ~ RATWGT,
        data = df,
        nest = TRUE
      )
      
      # Fit the model
      model <- svyglm(
        formula,
        family = quasibinomial,
        design = svy_design_temp
      )
      
      coef_table <- summary(model)$coefficients
      
      if (nrow(coef_table) > 1) {
        or_value <- exp(coef_table[2, 1])
        p_value <- coef_table[2, 4]
      } else {
        or_value <- NA
        p_value <- NA
      }
      
      # Calculate means and SDs
      case_mean <- mean(cases_df[[var$var_name]], na.rm = TRUE)
      case_sd <- sd(cases_df[[var$var_name]], na.rm = TRUE)
      control_mean <- mean(controls_df[[var$var_name]], na.rm = TRUE)
      control_sd <- sd(controls_df[[var$var_name]], na.rm = TRUE)
      
      baseline_descriptives <- rbind(
        baseline_descriptives,
        data.frame(
          Characteristic = var$label,
          Value = "Mean (SD)",
          Control = sprintf("%.1f (%.1f)", control_mean, control_sd),
          Case = sprintf("%.1f (%.1f)", case_mean, case_sd),
          OR = ifelse(is.na(or_value), "NA", sprintf("%.2f", or_value)),
          P_Value = format_p(p_value),
          stringsAsFactors = FALSE
        )
      )
      
    } else {
      # Categorical variables: Process each level
      current_var <- df[[var$var_name]]
      has_na <- anyNA(current_var)
      
      # Convert to factor if needed (preserve existing levels if already factor)
      if (!is.factor(current_var)) {
        # If not a factor, convert to factor using existing unique values in their order of appearance
        current_var <- factor(current_var, levels = unique(current_var))
      } else {
        # If already a factor, preserve the existing levels
        current_var <- factor(current_var, levels = levels(current_var))
      }
      
      # Add NA as a level if missing values exist
      if (has_na) {
        current_var <- addNA(current_var)
      }
      all_levels <- levels(current_var)
      
      # Get counts for cases and controls
      case_factor <- factor(cases_df[[var$var_name]], levels = all_levels)
      control_factor <- factor(controls_df[[var$var_name]], levels = all_levels)
      if (has_na) {
        case_factor <- addNA(case_factor)
        control_factor <- addNA(control_factor)
      }
      
      case_counts <- table(case_factor, useNA = "ifany")
      control_counts <- table(control_factor, useNA = "ifany")
      
      # Identify levels with zero counts in both cases and controls
      zero_levels <- all_levels[(case_counts == 0) & (control_counts == 0)]
      
      # Remove zero-count levels (unless it's the only level left)
      if (length(zero_levels) > 0 && length(all_levels) > length(zero_levels)) {
        current_var <- factor(current_var, 
                              levels = setdiff(all_levels, zero_levels))
        all_levels <- levels(current_var)
        
        # Recalculate counts with updated levels
        case_factor <- factor(cases_df[[var$var_name]], levels = all_levels)
        control_factor <- factor(controls_df[[var$var_name]], levels = all_levels)
        if (has_na) {
          case_factor <- addNA(case_factor)
          control_factor <- addNA(control_factor)
        }
        
        case_counts <- table(case_factor, useNA = "ifany")
        control_counts <- table(control_factor, useNA = "ifany")
      }
      
      # Calculate percentages
      case_pct <- prop.table(case_counts) * 100
      control_pct <- prop.table(control_counts) * 100
      
      # Create temporary factor in df for modeling
      temp_var_name <- paste0("temp_", var$var_name)
      df[[temp_var_name]] <- current_var
      
      # Create survey design object with the temporary variable
      svy_design_temp <- svydesign(
        ids = ~ PSU,
        strata = ~ PSUSTRAT,
        weights = ~ RATWGT,
        data = df,
        nest = TRUE
      )
      
      # Run survey logistic regression (handle errors)
      model_formula <- as.formula(
        paste("HIPR ~", temp_var_name))
      model <- svyglm(
        model_formula, 
        design = svy_design_temp, 
        family = quasibinomial())
      
      # Precompute ORs and p-values for each level
      or_values <- rep(NA, length(all_levels))
      p_values <- rep(NA, length(all_levels))
      or_values[1] <- 1.0  # Reference level OR = 1
      
      if (!is.null(model)) {
        coef_table <- summary(model)$coefficients
        n_coef <- nrow(coef_table)
        for (i in 2:length(all_levels)) {
          coef_name <- paste0(temp_var_name, all_levels[i])
          if (coef_name %in% rownames(coef_table)) {
            or_values[i] <- exp(coef_table[coef_name, "Estimate"])
            p_values[i] <- coef_table[coef_name, "Pr(>|t|)"]
          }
        }
      }
      
      # Remove temporary variable
      df[[temp_var_name]] <- NULL
      
      # Add rows for each level (reference level first)
      for (idx in 1:length(all_levels)) {
        level <- all_levels[idx]
        level_label <- ifelse(is.na(level), "Missing", as.character(level))
        
        control_count <- ifelse(is.na(control_counts[level]), 0, control_counts[level])
        control_pct_val <- ifelse(is.na(control_pct[level]), 0, control_pct[level])
        case_count <- ifelse(is.na(case_counts[level]), 0, case_counts[level])
        case_pct_val <- ifelse(is.na(case_pct[level]), 0, case_pct[level])
        
        baseline_descriptives <- rbind(
          baseline_descriptives,
          data.frame(
            Characteristic = if (idx == 1) var$label else "",
            Value = level_label,
            Control = sprintf("%d (%.1f%%)", control_count, control_pct_val),
            Case = sprintf("%d (%.1f%%)", case_count, case_pct_val),
            OR = ifelse(is.na(or_values[idx]), "NA", sprintf("%.2f", or_values[idx])),
            P_Value = format_p(p_values[idx]),
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
  
  # View the final table using kable
  knitr::kable(baseline_descriptives, format = "latex")
}












# load madymo_nhtsa_df_listwise from csv
madymo_nhtsa_df_listwise <- readRDS(file.path("data", "madymo_nhtsa_df_listwise.rds"))

# select only relevant columns
columns_final <- c("ID", "HIPR", "BAGDEPLOY_BIN", "BELTUSE_BIN", "BMI", 
                   "SEATPOS", "SEX_BIN", "YEAR", "AGE", "DV_MPH", "DVEST_MPH", 
                   "PSU", "PSUSTRAT", "RATWGT", "ALCINV", "DRGINV", "MAKE", 
                   "MODELYR", "VEHWGT", "WEATHER", "TIME", "TIME_cat", 
                   "VEHTYPE", "HIC15", "HIC15_log")
madymo_nhtsa_df_listwise <- madymo_nhtsa_df_listwise %>%
  select(
    all_of(columns_final))


# define occupant variables
occ_vars <- c("SEX_BIN", "AGE", "BMI")
# define environment vars
env_vars <- c("WEATHER", "MAKE", "MODELYR", "VEHWGT", "TIME", "VEHTYPE")
# define simulation variables
sim_vars <- c("DV_MPH", "SEATPOS", "BELTUSE_BIN", "BAGDEPLOY_BIN")


# ≥20 years old,\ndriver or front passenger seat positions,\nsingle-event frontal crashes\nbetween 5 and 40 miles per hour delta V

criteria_diagram_png_path <- file.path("r_docs", "images", "criteria_diagram.png")
design <- tibble::tribble(
  ~left,                           ~n_left, ~right,              ~n_right,
  "Base population",        113579,    "Demographic and crash\nexclusion criteria applied",       98487,
  "Study population",              15092,     "Participants with\nmissing exposure data", 1610,
  "Complete-case set", 13482,     "",                  NA_integer_
  )

# Plot
exclusion_flowchart(design, width = 2) %>%
  export_svg() %>%
  charToRaw() %>%
  rsvg_png(file = criteria_diagram_png_path,
           width = 1800)




# criteria_diagram_png_long_path <- file.path("r_docs", "images", "criteria_diagram.png")
# grViz(
#   "digraph flowchart {
#       # node definitions with substituted label text
#       node [fontname = Helvetica, shape = rectangle]        
#       tab1 [label = '@@1']
#       tab2 [label = '@@2']
#       tab3 [label = '@@3']
#       tab4 [label = '@@4']
#       tab5 [label = '@@5']
#       tab6 [label = '@@6']
#       tab7 [label = '@@7']
#       tab8 [label = '@@8']
#       tab9 [label = '@@9']
# 
#       # edge definitions with the node IDs
#       tab1 -> tab2 -> tab3 -> tab4 -> tab5 -> tab6 -> tab7 -> tab8 -> tab9;
#       }
# 
#       [1]: 'Count before any filters: 113579'
#       [2]: 'Count after ≥20 years old: 84268'
#       [3]: 'Count after only driver or front passenger seat positions: 73209'
#       [4]: 'Count after single-event crashes: 34663'
#       [5]: 'Count after only frontal crashes: 21349'
#       [6]: 'Count after between 5 and 40 miles per hour delta V: 15092'
#       [7]: 'Count after no seatbelt NA: 13490'
#       [8]: 'Count after no airbag NA: 13482'
#       [9]: 'Final count: 13482'
#       "
# ) %>%
#   export_svg() %>%
#   charToRaw() %>%
#   rsvg_png(file = criteria_diagram_png_long_path,
#            width = 1800)



# load models from madymo_nhtsa_models.RData
load(file.path("data", "madymo_nhtsa_models.RData"))

# load in- exclusion criteria statements 
load(file.path("data", "criteria_list.RData"))

# compare performance of models
comp <- c(
  "univariate_model_listwise",
  "multivariate_model_listwise_final",
  "multivariate_model_listwise_final_interaction")

models_listwise_comparison_df <- compare_df_fun(comp)
models_listwise_comparison_df$Model <- c(
  "Univariate", 
  "Final", 
  "Interaction")
print(models_listwise_comparison_df)


# present models on single grid
models <- list(
 "Univariate" = univariate_model_listwise,
 "Final" = multivariate_model_listwise_final
)
models_prep <- roc_grid_prep_fun(models)
univariate_final_roc_grid_viz <- roc_compare_grid_fun(models_prep, col_num = 2)



# compare relevant terms in interaction grid
load(
  file = file.path("data", "interaction_roc_grid_obj.RData"
  ))
interaction_grid_viz <- roc_grid_plot(interaction_roc_grid_obj)

load(
  file = file.path("data", "interaction_roc_grid_obj_new_bmi_split.RData"
  ))
interaction_grid_viz_new_bmi_split <- roc_grid_plot(interaction_roc_grid_obj)

table_1 <- baseline_desc_fun(madymo_nhtsa_df_listwise)
