#' @title process_study_scores_kidney OR kidney_om_lb_mi_tox_score_list
#'
#' @description
#' working: gives exact accurate representation as Amin's Liver f6 function works, but doesn't work for all studyids, use the function for further harmonization.
#' This function processes kidney organ toxicity scores: orchestrates BW, Kidney-to-BW, lb, and MI score calculations for the kidney domain
#' External function dependencies: Assumes external functions loaded in the package:
#'                        get_compile_data(), get_bw_score(), kidneytobw_score(), kidney_lb_score(), kidney_mi_score()
#' for a set of studies or XPT files. It can output individual scores, z-scores by USUBJID, or averaged scores
#' for multiple studies, and handles errors (like no data presents, but function still continues and load data from available studyids) during the processing steps.
#'
#' @param studyid_or_studyids A character vector or a single study ID to process.
#' If multiple studies (for machine learning) are provided, the function processes each study sequentially. (Mandatory)
#'
#' @param path_db A character string specifying the path to the database or directory containing the data files.
#' (Mandatory)
#'
#' @param fake_study A boolean flag indicating if the study data is simulated (`TRUE`) or real (`FALSE`). Default is `FALSE`. (Optional)
#'
#' @param use_xpt_file A boolean flag indicating whether to use an XPT file for the study data. Default is `FALSE`. (Optional)
#'
#' @param output_individual_scores A boolean flag indicating whether individual scores should be returned (`TRUE`) or averaged scores (`FALSE`). Default is `FALSE`. (Optional)
#'
#' @param output_zscore_by_USUBJID A boolean flag indicating whether to output z-scores by `USUBJID` (`TRUE`) or averaged scores (`FALSE`). Default is `FALSE`. (Optional)
#'
#' @return A data frame containing the calculated scores for each study. The type of result depends on the flags passed:
#' - If `output_individual_scores` is `TRUE`, a data frame with individual scores for each study is returned.
#' - If `output_zscore_by_USUBJID` is `TRUE`, a data frame with z-scores by `USUBJID` for each study is returned.
#' - If neither flag is set, the function returns a data frame with averaged scores for each study.
#'
#' @examples
#' \dontrun{
#' # Get averaged scores for a single study
#' result_default <- kidney_mi_score(studyid, path_db)
#'
#' # Get individual scores for multiple studies
#' result_individual_multiple_score <- kidney_mi_score(
#'                                   studyid = "123456",
#'                                   path_db = "path/to/database.db",
#'                                   master_compiledata = NULL,
#'                                   return_individual_scores = TRUE,
#'                                   return_zscore_by_USUBJID = FALSE)
#'
#' result_zscore <- kidney_mi_score("123456", path_db,
#'                                   master_compiledata = NULL,
#'                                   return_individual_scores = FALSE,
#'                                   return_zscore_by_USUBJID = TRUE)
#' }
#'
#' @export
process_kidney_scores_limited <- function(studyid_or_studyids, path_db, fake_study = FALSE,
                                  use_xpt_file = FALSE, output_individual_scores = FALSE,
                                  output_zscore_by_USUBJID = FALSE) {
    
    # Embedded safe_execute
    
    # calculate_kidneyToBW_zscore
    calculate_kidneyToBW_zscore <- function(studyid, path_db, fake_study, use_xpt_file, master_compiledata,
                                            output_individual_scores, output_zscore_by_USUBJID,
                                            bwzscore_BW, master_kidneyToBW, FOUR_Kidney_Score_avg, master_error_df) {
      # If using XPT files, no study ID is needed (NULL); otherwise, use studyid
      studyid_local <- if (use_xpt_file) NULL else studyid
      
      if (output_individual_scores) {
        HD_kidney_zscore_df <- kidneytobw_score(
          studyid = studyid_local,
          path_db = path_db,
          fake_study = fake_study,
          use_xpt_file = use_xpt_file,
          master_compiledata = master_compiledata,
          bwzscore_BW = bwzscore_BW,
          return_individual_scores = TRUE,
          return_zscore_by_USUBJID = FALSE
        )
        master_kidneyToBW <- rbind(master_kidneyToBW, as.data.frame(HD_kidney_zscore_df))
        
      } else if (output_zscore_by_USUBJID) {
        kidneyTOBW_zscore_by_USUBJID_HD <- kidneytobw_score(
          studyid = studyid_local,
          path_db = path_db,
          fake_study = fake_study,
          use_xpt_file = use_xpt_file,
          master_compiledata = master_compiledata,
          bwzscore_BW = bwzscore_BW,
          return_individual_scores = FALSE,
          return_zscore_by_USUBJID = TRUE
        )
        kidneyTOBW_study_identifier <- unique(kidneyTOBW_zscore_by_USUBJID_HD$STUDYID)
        master_kidneyToBW[[as.character(kidneyTOBW_study_identifier)]] <- as.data.frame(kidneyTOBW_zscore_by_USUBJID_HD)
        
      } else {
        calculated_kidneyToBW_value <- safe_execute({
          averaged_kidneyToBW_df <- kidneytobw_score(
            studyid = studyid_local,
            path_db = path_db,
            fake_study = fake_study,
            use_xpt_file = use_xpt_file,
            master_compiledata = master_compiledata,
            bwzscore_BW = bwzscore_BW,
            return_individual_scores = FALSE,
            return_zscore_by_USUBJID = FALSE
          )
          kidneyToBW_df <- dplyr::rename(averaged_kidneyToBW_df, kidneyToBW_avg = avg_kidneyToBW_zscore)
          kidneyToBW_df$kidneyToBW_avg[1]
        }, studyid, "KidneyToBW", use_xpt_file, path_db, master_error_df)
        
        study_id <- as.character(unique(master_compiledata$STUDYID))
        row_index <- which(as.character(FOUR_Kidney_Score_avg$STUDYID) == study_id)
        if (length(row_index) > 0 && !is.null(calculated_kidneyToBW_value)) {
          FOUR_Kidney_Score_avg$kidneyToBW_avg[row_index] <- calculated_kidneyToBW_value
        }
      }
      
      return(list(updated_kidneyToBW = master_kidneyToBW, updated_score = FOUR_Kidney_Score_avg))
    }
    
    safe_execute <- function(expr, studyid, block, use_xpt_file, path_db, master_error_df) {
      tryCatch(
        expr,
        error = function(e) {
          error_record <- data.frame(
            STUDYID = if (use_xpt_file) path_db else studyid,
            Block = block,
            ErrorMessage = e$message,
            stringsAsFactors = FALSE
          )
          master_error_df <<- rbind(master_error_df, error_record)
          return(NULL)
        }
      )
    }
    
    # Flag guard
    if (output_individual_scores && output_zscore_by_USUBJID) {
      stop("Error: Both 'output_individual_scores' and 'output_zscore_by_USUBJID' cannot be TRUE at the same time.")
    }
    
    # Initialize containers
    master_kidneyToBW <- if (output_individual_scores) data.frame() else list()
    master_lb_score_six <- data.frame()
    master_lb_score <- list()
    master_mi_df <- data.frame()
    master_mi_score <- list()
    FOUR_Kidney_Score_avg <- data.frame(STUDYID = character(), BWZSCORE_avg = numeric(),
                                        kidneyToBW_avg = numeric(), LB_score_avg = numeric(), MI_score_avg = numeric())
    master_error_df <- data.frame(STUDYID = character(), Block = character(),
                                  ErrorMessage = character(), stringsAsFactors = FALSE)
    
    for (studyid in studyid_or_studyids) {
      if (use_xpt_file) path_db <- studyid
      print(paste("Processing study:", studyid))
      
      master_compiledata <- safe_execute({
        get_compile_data(
          studyid = if (use_xpt_file) NULL else studyid,
          path_db = path_db,
          fake_study = fake_study,
          use_xpt_file = use_xpt_file
        )
      }, studyid, "compiledata", use_xpt_file, path_db, master_error_df)
      
      if (is.null(master_compiledata)) next
      
      if (!output_individual_scores && !output_zscore_by_USUBJID) {
        FOUR_Kidney_Score_avg <- rbind(FOUR_Kidney_Score_avg, data.frame(
          STUDYID = unique(master_compiledata$STUDYID),
          BWZSCORE_avg = NA,
          kidneyToBW_avg = NA,
          LB_score_avg = NA,
          MI_score_avg = NA
        ))
      }
      
      res_bw <- {
        studyid_local <- if (use_xpt_file) NULL else studyid
        result <- safe_execute({
          if (output_individual_scores) {
            get_bw_score(studyid_local, path_db, fake_study, use_xpt_file, master_compiledata,
                         TRUE, FALSE)
          } else if (output_zscore_by_USUBJID) {
            as.data.frame(get_bw_score(studyid_local, path_db, fake_study, use_xpt_file, master_compiledata,
                                       FALSE, TRUE))
          } else {
            averaged_HD_BWzScore <- get_bw_score(studyid_local, path_db, fake_study, use_xpt_file,
                                                 master_compiledata, FALSE, FALSE)
            calculated_BWzScore_value <- averaged_HD_BWzScore$BWZSCORE_avg[1]
            FOUR_Kidney_Score_avg$BWZSCORE_avg[FOUR_Kidney_Score_avg$STUDYID == unique(master_compiledata$STUDYID)] <- calculated_BWzScore_value
          }
        }, studyid, "BWZscore", use_xpt_file, path_db, master_error_df)
        list(result = result, updated_score = FOUR_Kidney_Score_avg)
      }
      FOUR_Kidney_Score_avg <- res_bw$updated_score
      
      res_lb <- {
        studyid_local <- if (use_xpt_file) NULL else studyid
        if (output_individual_scores) {
          master_lb_scores <- kidney_lb_score(
            studyid_local, path_db, fake_study, use_xpt_file,
            master_compiledata, TRUE, FALSE
          )
          master_lb_score_six <- rbind(master_lb_score_six, master_lb_scores)
        } else if (output_zscore_by_USUBJID) {
          LB_zscore_by_USUBJID_HD <- kidney_lb_score(
            studyid_local, path_db, fake_study, use_xpt_file,
            master_compiledata, FALSE, TRUE
          )
          lb_study_identifier <- unique(LB_zscore_by_USUBJID_HD$STUDYID)
          master_lb_score[[as.character(lb_study_identifier)]] <- LB_zscore_by_USUBJID_HD
        } else {
          calculated_LB_value <- safe_execute({
            averaged_LB_score <- kidney_lb_score(
              studyid_local, path_db, fake_study, use_xpt_file,
              master_compiledata, FALSE, FALSE
            )
            averaged_LB_score$LB_score_avg[1]
          }, studyid, "LB", use_xpt_file, path_db, master_error_df)
          FOUR_Kidney_Score_avg$LB_score_avg[FOUR_Kidney_Score_avg$STUDYID == unique(master_compiledata$STUDYID)] <- calculated_LB_value
        }
      }
      
      res_mi <- {
        studyid_local <- if (use_xpt_file) NULL else studyid
        if (output_individual_scores) {
          mi_score_final_list_df <- kidney_mi_score(
            studyid_local, path_db, fake_study, use_xpt_file,
            master_compiledata, TRUE, FALSE
          )
          master_mi_df <- dplyr::bind_rows(master_mi_df, mi_score_final_list_df)
        } else if (output_zscore_by_USUBJID) {
          MI_score_by_USUBJID_HD <- kidney_mi_score(
            studyid_local, path_db, fake_study, use_xpt_file,
            master_compiledata, FALSE, TRUE
          )
          mi_study_identifier <- unique(MI_score_by_USUBJID_HD$STUDYID)
          master_mi_score[[as.character(mi_study_identifier)]] <- MI_score_by_USUBJID_HD
        } else {
          calculated_MI_value <- safe_execute({
            averaged_MI_score <- kidney_mi_score(
              studyid_local, path_db, fake_study, use_xpt_file,
              master_compiledata, FALSE, FALSE
            )
            averaged_MI_score$MI_score_avg[1]
          }, studyid, "MI", use_xpt_file, path_db, master_error_df)
          FOUR_Kidney_Score_avg$MI_score_avg[FOUR_Kidney_Score_avg$STUDYID == unique(master_compiledata$STUDYID)] <- calculated_MI_value
        }
      }
      
      res_kidneyToBW <- calculate_kidneyToBW_zscore(
        studyid, path_db, fake_study, use_xpt_file, master_compiledata,
        output_individual_scores, output_zscore_by_USUBJID,
        bwzscore_BW = NULL, master_kidneyToBW, FOUR_Kidney_Score_avg, master_error_df
      )
      master_kidneyToBW <- res_kidneyToBW$updated_kidneyToBW
      FOUR_Kidney_Score_avg <- res_kidneyToBW$updated_score
      
    }
    
    # Return based on mode
    if (output_individual_scores) {
      return(dplyr::full_join(master_kidneyToBW, master_lb_score_six, by = "STUDYID") %>%
               dplyr::full_join(master_mi_df, by = "STUDYID"))
    } else if (output_zscore_by_USUBJID) {
      combined_df <- dplyr::bind_rows(master_kidneyToBW) %>%
        dplyr::full_join(dplyr::bind_rows(master_lb_score), by = c("STUDYID", "USUBJID")) %>%
        dplyr::full_join(dplyr::bind_rows(master_mi_score), by = c("STUDYID", "USUBJID"))
      return(combined_df)
    } else {
      FOUR_Kidney_Score_avg[, 2:ncol(FOUR_Kidney_Score_avg)] <- round(FOUR_Kidney_Score_avg[, 2:ncol(FOUR_Kidney_Score_avg)], 2)
      return(FOUR_Kidney_Score_avg)
    }
  }
  