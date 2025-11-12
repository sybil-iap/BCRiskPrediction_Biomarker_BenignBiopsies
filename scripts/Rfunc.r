library(glmnet)
library(caret)
library(ROCR);
library(plotROC)
library(pROC)

library(ggplot2)
library(ggrepel)

#########
#### correlation plot
corr_heat <- function(RNAseq_log2fc_sig = RNAseq_log2fc_sig_ERneg,
                      filepath = "./covar_mat.pdf")
{
    cor(RNAseq_log2fc_sig, use = "pairwise.complete.obs") %>%
    reshape2::melt() %>%
    ggplot(aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Pearson\nCorrelation") +
    theme(text = element_text(size = 24, face = "bold"),
          axis.text.x = element_text(angle = 90, vjust =0.3,
                                     hjust = 1, face = "bold"), 
          legend.position = "right") +
    labs(x = NULL, y = NULL)
ggsave(filepath, width = 14, height = 10)
}

######### Volcano plot function
volcano_plot <- function(DGE_pval =  DGE_pval, log2fc_threshold = 0.6,
                         pvalue_threshold = 0.05, filepath = "", file_title = ""){
    DGE_pval$DE_pos <- "NO"
    DGE_pval$DE_pos[DGE_pval$log2FoldChange > 0.6 & DGE_pval$pvalue < 0.001] <- "UP"
    DGE_pval$DE_pos[DGE_pval$log2FoldChange < -0.6 & DGE_pval$pvalue < 0.001] <- "DOWN"
    
    UP_gene_list=DGE_pval[DGE_pval$DE_pos=="UP", ]
    DOWN_gene_list=DGE_pval[DGE_pval$DE_pos=="DOWN", ]
    
    up_top_gene = UP_gene_list[order(UP_gene_list$pvalue, decreasing=FALSE)[1:5], "gene_name"]
    down_top_gene = DOWN_gene_list[order(DOWN_gene_list$pvalue,decreasing=FALSE)[1:5], "gene_name"]
    up_down_list = unlist(c(up_top_gene,down_top_gene))
    
    DGE_pval$delabel <- NA
    DGE_pval$delabel[DGE_pval$gene_name %in% up_down_list] <- unlist(DGE_pval[DGE_pval$gene_name %in% up_down_list, "gene_name"])
    
    DGE_pval_label <- dplyr::filter(DGE_pval, !is.na(delabel))
    
    # x:log2fodchange, y:pvalue
    options(ggrepel.max.overlaps = Inf)
    ggplot(data=DGE_pval, aes(x=log2FoldChange, y=-log10(pvalue),
                              col=DE_pos)) +
        geom_point() +
        geom_vline(xintercept=c(-log2fc_threshold, log2fc_threshold), col="red") +
        geom_hline(yintercept=-log10(pvalue_threshold), col="red") +
        labs(x = "Log2FC", y = "-log10(p-value)", title = file_title, colour = "DE") +
        geom_label_repel(data = DGE_pval_label,
                         aes(x=log2FoldChange, y=-log10(pvalue),
                             label = delabel),
                         position = position_jitter(width = 0.5, seed = 2),
                         arrow = arrow(length = unit(0.01, "npc")),
                         size = 5,  fontface="bold", angle = 20,
                         segment.curvature = -0.1,
                         segment.ncp = 3, segment.angle = 20,
                         point.padding = 0.1, box.padding = 0.5) +
        theme(text = element_text(size = 24, face = "bold"),
              legend.position="bottom")
    ggsave(filepath, width = 8, height = 6)
}



######
getAUC <- function(D, pred){
    return( roc(D, pred, quiet = TRUE)$auc )
}

#var_list_x = cov_list
#data_train = glm_dt
## My EN_CV function
# var_list = colnames(dataset_fit_lm_noMCI_sub4)

myENcv <- function(data_train, var_list_x, cv_fold = 10, tune_length = 10, 
                   alpha_vec = seq(0, 1, length.out = 10)){
    set.seed(2019)            
    trainIndex_folds = groupKFold(data_train$ID, k = cv_fold)
    lambda_vec = exp( seq(-3, 6, length.out = 100) )
    best_lambda_vec <- best_auc <- rep(NA, tune_length)
    
    for(i in 1:tune_length) {
        Test = NULL
        Test_pred = NULL
        for(k in 1:cv_fold){
            temp_train = data_train[trainIndex_folds[[k]], ]
            temp_test = data_train[-trainIndex_folds[[k]], ]
            Test = c(Test, temp_test[, "DeveloppedCancer"])
            fit_temp = glmnet(x = as.matrix(temp_train[, eval(var_list_x)]), 
                              y = temp_train[, "DeveloppedCancer"],
                              family = "binomial", alpha = alpha_vec[i], lambda = lambda_vec)
            pred_temp = predict(fit_temp, newx = as.matrix(temp_test[, eval(var_list_x)]), type = "response")
            Test_pred = rbind(Test_pred, pred_temp)
        }
        AUC_temp = apply(Test_pred, 2, getAUC, D = Test)
        best_lambda_vec[i] = fit_temp$lambda[which.max(AUC_temp)]
        best_auc[i] = max(AUC_temp)
    }
    best_alpha = alpha_vec[which.max(best_auc)]
    best_lambda = best_lambda_vec[which.max(best_auc)]
    print("Best alpha, best lambda, max auc:")
    print(c(best_alpha, best_lambda, max(best_auc)))
    
    best_fit = glmnet(x = as.matrix(data_train[, eval(var_list_x)]),
                      y = data_train[, "DeveloppedCancer"],
                      family = "binomial", alpha = best_alpha, lambda = best_lambda)
    return(list(best_alpha = best_alpha, best_lambda = best_lambda,
                best_auc = max(best_auc), best_fit = best_fit) )
}

my_logistic <- function(data_train, var_list_x, cv_fold = 10){
  set.seed(2025)
  # groupKFold: splits the data based on a grouping factor
  trainIndex_fold = groupKFold(data_train$ID, k = cv_fold)
  auc_vec <- rep(NA, cv_fold)
  
  for (k in 1:cv_fold){
    temp_train <- data_train[trainIndex_fold[[k]],]
    temp_test <- data_train[-trainIndex_fold[[k]],]
    
    temp_fit <- glm(DeveloppedCancer ~ ., # y ~ . means all variables in the data except the responce variable
                    data = temp_train[, c(var_list_x, "DeveloppedCancer")],
                    family = "binomial")
    temp_pred <- predict(temp_fit, 
                         newdata = temp_test[, c(var_list_x, "DeveloppedCancer")],
                         type = "response")
    # Calculate AUC using the roc function
    auc_vec[k] <- roc(response = temp_test[, "DeveloppedCancer"], predictor = temp_pred, quiet = TRUE)$auc
  }
  
  final_fit <- glm(DeveloppedCancer ~ ., 
                   data = data_train[, c(var_list_x, "DeveloppedCancer")],
                   family = "binomial")
  mean_auc <- mean(auc_vec, na.rm = TRUE)
  return(list(model = final_fit, mean_auc = mean_auc))
}

plot_LRcoef <- function(fit = fit_temp$best_fit, title = "Trained Model", 
                        pdf_file = "beta.pdf", add_legend = TRUE,
                        beta_abs_cut = 0.1, 
                        width = 5, height = 6){
  
  coef <-  fit %>%  coefficients() %>% as.matrix()
  coef <-  coef[-1, ]
  coef <- coef %>% as_tibble() %>%
    dplyr::mutate(Var = names(coef), beta = coef) %>%
    dplyr::filter( abs(beta) >  beta_abs_cut) %>% dplyr::arrange(beta) %>% as.data.frame()
  
  plot <- coef %>% ggplot(aes(x = Var, y = beta, fill=beta>0)) +
    geom_bar(stat = "identity", width = 0.5, show.legend = add_legend) +
    coord_flip() +
    scale_x_discrete(limits = (coef$Var) ) +
    labs(x = NULL, y = NULL, title = title) +
    theme(text = element_text(size = 20, face = "bold"),
          axis.text.x = element_text(angle = -90, vjust = -0.1))
  
  ggsave(pdf_file, width = width, height = height)
  return(plot)
}

# Plot EN coefficients
plot_ENcoef <- function(fit = fit_temp$best_fit, title = "Trained Model", 
                        pdf_file = "beta.pdf", add_legend = TRUE,
                        beta_abs_cut = 0.1, 
                        width = 5, height = 6){
    
    coef <-  fit %>%  coefficients() %>% as.matrix()
    coef <-  coef[-1, ]
    coef <- coef %>% as_tibble() %>%
        dplyr::mutate(Var = names(coef), beta = coef) %>%
        dplyr::filter( abs(beta) >  beta_abs_cut) %>% dplyr::arrange(beta) %>% as.data.frame()
    
    
    plot <- coef %>% ggplot(aes(x = Var, y = beta, fill=beta>0)) +
        geom_bar(stat = "identity", width = 0.5, show.legend = add_legend) +
        coord_flip() +
        scale_x_discrete(limits = (coef$Var) ) +
        labs(x = NULL, y = NULL, title = title) +
        theme(text = element_text(size = 20, face = "bold"),
              axis.text.x = element_text(angle = -90, vjust = -0.1))

    ggsave(pdf_file, width = width, height = height)
    return(plot)
}



