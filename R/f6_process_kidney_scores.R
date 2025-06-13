#' @title process_study_scores_kidney OR kidney_om_lb_mi_tox_score_list
#'
#' @description
#' Working: giving extra cols such as USUBJID, may need to remove extra cols if needed.
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

process_kidney_scores <- function(studyid_or_studyids, path_db, fake_study = FALSE,
                                  use_xpt_file = FALSE, output_individual_scores = FALSE,
                                  output_zscore_by_USUBJID = FALSE) {

  # Helper to catch errors and suppress warnings
  ##safe_execute: wrap any expression (expr) in a tryCatch so that if an error occurs,
  ##it records that error in master_error_df and returns NULL instead of halting everything.
  safe_execute <- function(expr, studyid, block, use_xpt_file, path_db, master_error_df) {
    tryCatch(
      suppressWarnings(expr),       # expr - actual R expression to evaluate calling function such as get_compile_data
      error = function(e) {
        error_record <- tibble(
          STUDYID = if (use_xpt_file) path_db else studyid,
          Block = block,           #a string such as "BWZscore", "LB" etc
          ErrorMessage = e$message
        )
        master_error_df <<- bind_rows(master_error_df, error_record)
        return(NULL)
      }
    )
  }

  # Nested function (1): calculate Kidney:BW ratio or return per-USUBJID when requested
  calculate_kidneyToBW_zscore <- function(studyid, path_db, fake_study, use_xpt_file, master_compiledata,
                                          output_individual_scores, output_zscore_by_USUBJID,
                                          bwzscore_BW, master_error_df) {
# Note: All three names—studyid, studyid_local, and sid within each nested functions refer to “the same thing,”
# namely “the current study ID (or NULL if we’re using an XPT file).”
    studyid_local <- if (use_xpt_file) NULL else studyid
    if (output_zscore_by_USUBJID) {
      df_kbw <- safe_execute(
        kidneytobw_score(
          studyid = studyid_local,
          path_db = path_db,
          fake_study = fake_study,
          use_xpt_file = use_xpt_file,
          master_compiledata = master_compiledata,
          bwzscore_BW = bwzscore_BW,
          return_individual_scores = FALSE,
          return_zscore_by_USUBJID = TRUE
        ),
        studyid, 'KidneyToBW', use_xpt_file, path_db, master_error_df
      )
# Ensures that the returned data frame always has a USUBJID column.
# If the helper returned no USUBJID column, we insert a dummy column of NA
      if (!is.null(df_kbw) && !('USUBJID' %in% names(df_kbw))) {
        df_kbw <- df_kbw %>% mutate(USUBJID = NA_character_)
      }
      return(df_kbw)
#This branch handles “study‐level summary” when neither flag is true. In other words, we want a single numeric kidneyToBW_avg for the entire study.
    } else if (!output_individual_scores) {
      val <- safe_execute({
        avg_df <- kidneytobw_score(
          studyid = studyid_local,
          path_db = path_db,
          fake_study = fake_study,
          use_xpt_file = use_xpt_file,
          master_compiledata = master_compiledata,
          bwzscore_BW = bwzscore_BW,
          return_individual_scores = FALSE,
          return_zscore_by_USUBJID = FALSE
        )
        dplyr::rename(avg_df, kidneyToBW_avg = avg_kidneyToBW_zscore)$kidneyToBW_avg[1]
      }, studyid, 'KidneyToBW', use_xpt_file, path_db, master_error_df)
      return(val)
    } else {
      # Individual mode does not call this nested function for per-USUBJID
      return(NULL)
    }
  }

  # Prevent conflicting flags
  if (output_individual_scores && output_zscore_by_USUBJID) stop('Cannot request both individual and USUBJID-level scores.')

  # Initialize containers
  if (output_individual_scores) {
# if TRUE, create below empty tibbles, each will hold that study’s single row (per study) for BW, LB, and MI, respectively.
    master_bw <- tibble()
    master_lb <- tibble()
    master_mi <- tibble()
  }
  if (output_zscore_by_USUBJID) {
    master_usubjid_bw <- tibble()
    master_usubjid_lb <- tibble()
    master_usubjid_mi <- tibble()
  }
  FOUR_Kidney_Score_avg <- tibble(
# FOUR_Kidney_Score_avg—an empty tibble with columns (STUDYID, BWZSCORE_avg, kidneyToBW_avg, LB_score_avg, MI_score_avg)
# that will collect the study‐level summary when both flags are FALSE.
    STUDYID = character(),
    BWZSCORE_avg = numeric(),
    kidneyToBW_avg = numeric(),
    LB_score_avg = numeric(),
    MI_score_avg = numeric()
  )
  master_error_df <- tibble(STUDYID = character(), Block = character(), ErrorMessage = character())

  # Loop over each study ID
  for (study in studyid_or_studyids) { #loops over each provided study identifier or, if use_xpt_file=TRUE, each provided file path.
    if (use_xpt_file) path_db <- study  #If we’re in “XPT mode,” treat study itself as the path_db (the downstream helper functions expect path_db to be a file path).
    message('Processing study: ', study) #Prints a progress message to the console so the user knows which study is currently being handled.

    # Nested function (2): Retrieve compiled data
#returns a data frame (mc) containing all raw data from the SEND dataset for that study (BW, LB, MI domains, etc.).
    mc <- safe_execute(
      get_compile_data(
        studyid = if (use_xpt_file) NULL else study,
        path_db = path_db,
        fake_study = fake_study,
        use_xpt_file = use_xpt_file
      ),
      study, 'compiledata', use_xpt_file, path_db, master_error_df
    )
#If we failed to fetch or process the compiled data, skip to the next iteration (i.e. move on to the next study ID without doing any further work here)
    if (is.null(mc)) next

    # For summary mode: add one placeholder row per study
    if (!output_individual_scores && !output_zscore_by_USUBJID) {
      FOUR_Kidney_Score_avg <- bind_rows(
        FOUR_Kidney_Score_avg,
        tibble(
          STUDYID = as.character(unique(mc$STUDYID)),
          BWZSCORE_avg = NA_real_,
          kidneyToBW_avg = NA_real_,
          LB_score_avg = NA_real_,
          MI_score_avg = NA_real_
        )
      )
    }

    # Nested function (3)----- BW z-scores -----
    bw_res <- safe_execute({
      sid <- if (use_xpt_file) NULL else study
      get_bw_score(
        studyid = sid,
        path_db = path_db,
        fake_study = fake_study,
        use_xpt_file = use_xpt_file,
        master_compiledata = mc,
        return_individual_scores = output_individual_scores,
        return_zscore_by_USUBJID = output_zscore_by_USUBJID
      )
    }, study, 'BWZscore', use_xpt_file, path_db, master_error_df)

    if (output_individual_scores && !is.null(bw_res)) {
      master_bw <- bind_rows(master_bw, bw_res)
    }
    if (output_zscore_by_USUBJID && !is.null(bw_res)) {
      if (!('USUBJID' %in% names(bw_res))) bw_res <- bw_res %>% mutate(USUBJID = NA_character_)
      master_usubjid_bw <- bind_rows(master_usubjid_bw, bw_res %>% mutate(USUBJID = as.character(USUBJID)))
    }
    if (!output_individual_scores && !output_zscore_by_USUBJID && !is.null(bw_res)) {
      FOUR_Kidney_Score_avg$BWZSCORE_avg[FOUR_Kidney_Score_avg$STUDYID == unique(mc$STUDYID)] <- bw_res$BWZSCORE_avg[1]
    }

    # Nested function (4)----- LB scores -----
    lb_res <- safe_execute({
      sid <- if (use_xpt_file) NULL else study
      kidney_lb_score(
        studyid = sid,
        path_db = path_db,
        fake_study = fake_study,
        use_xpt_file = use_xpt_file,
        master_compiledata = mc,
        return_individual_scores = output_individual_scores,
        return_zscore_by_USUBJID = output_zscore_by_USUBJID
      )
    }, study, 'LB', use_xpt_file, path_db, master_error_df)

    if (output_individual_scores && !is.null(lb_res)) {
      master_lb <- bind_rows(master_lb, lb_res)
    }
    if (output_zscore_by_USUBJID && !is.null(lb_res)) {
      if (!('USUBJID' %in% names(lb_res))) lb_res <- lb_res %>% mutate(USUBJID = NA_character_)
      master_usubjid_lb <- bind_rows(master_usubjid_lb, lb_res %>% mutate(USUBJID = as.character(USUBJID)))
    }
    if (!output_individual_scores && !output_zscore_by_USUBJID && !is.null(lb_res)) {
      FOUR_Kidney_Score_avg$LB_score_avg[FOUR_Kidney_Score_avg$STUDYID == unique(mc$STUDYID)] <- lb_res$LB_score_avg[1]
    }

    # Nested function (5)----- MI scores -----
    mi_res <- safe_execute({
      sid <- if (use_xpt_file) NULL else study
      kidney_mi_score(
        studyid = sid,
        path_db = path_db,
        fake_study = fake_study,
        use_xpt_file = use_xpt_file,
        master_compiledata = mc,
        return_individual_scores = output_individual_scores,
        return_zscore_by_USUBJID = output_zscore_by_USUBJID
      )
    }, study, 'MI', use_xpt_file, path_db, master_error_df)

    if (output_individual_scores && !is.null(mi_res)) {
      master_mi <- bind_rows(master_mi, mi_res)
    }
    if (output_zscore_by_USUBJID && !is.null(mi_res)) {
      if (!('USUBJID' %in% names(mi_res))) mi_res <- mi_res %>% mutate(USUBJID = NA_character_)
      master_usubjid_mi <- bind_rows(master_usubjid_mi, mi_res %>% mutate(USUBJID = as.character(USUBJID)))
    }
    if (!output_individual_scores && !output_zscore_by_USUBJID && !is.null(mi_res)) {
      FOUR_Kidney_Score_avg$MI_score_avg[FOUR_Kidney_Score_avg$STUDYID == unique(mc$STUDYID)] <- mi_res$MI_score_avg[1]
    }

    # ----- Kidney:BW ratio for summary mode -----
    if (!output_individual_scores && !output_zscore_by_USUBJID) {
      ktbw_val <- calculate_kidneyToBW_zscore(
        study, path_db, fake_study, use_xpt_file, mc,
        output_individual_scores = FALSE,
        output_zscore_by_USUBJID = FALSE,
        bwzscore_BW = NULL,
        master_error_df = master_error_df
      )
      FOUR_Kidney_Score_avg$kidneyToBW_avg[FOUR_Kidney_Score_avg$STUDYID == unique(mc$STUDYID)] <- ktbw_val
    }
  }

  # Return results based on requested mode
  if (output_individual_scores) {
    # Take the first row per study from each of BW, LB, MI and join by STUDYID to ensure one row per study
    master_bw1 <- master_bw %>% group_by(STUDYID) %>% slice(1) %>% ungroup()
    master_lb1 <- master_lb %>% group_by(STUDYID) %>% slice(1) %>% ungroup()
    master_mi1 <- master_mi %>% group_by(STUDYID) %>% slice(1) %>% ungroup()
    result_individual <- master_bw1 %>%
      full_join(master_lb1, by = "STUDYID") %>%
      full_join(master_mi1, by = "STUDYID")
    return(result_individual)
  } else if (output_zscore_by_USUBJID) {
    combined_df <- master_usubjid_bw %>%
      full_join(master_usubjid_lb, by = c("STUDYID", "USUBJID")) %>%
      full_join(master_usubjid_mi, by = c("STUDYID", "USUBJID"))
    return(combined_df)
  } else {
    FOUR_Kidney_Score_avg <- FOUR_Kidney_Score_avg %>%
      mutate(across(-STUDYID, ~ round(.x, 2)))
    return(FOUR_Kidney_Score_avg)
  }
}
