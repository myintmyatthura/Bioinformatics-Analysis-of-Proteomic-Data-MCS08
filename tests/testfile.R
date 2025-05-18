# =============================================================
# Shiny‑App Test Suite
# =============================================================
# This single file contains both low‑level unit tests (white‑box)
# and high‑level end‑to‑end UI tests (black‑box) for the entire
# proteomics‑analysis Shiny application.
#
# ▸ To run everything:  devtools::test()   (in the app root)
# ▸ To view coverage:   covr::file_coverage("main.R", "tests/testfile.R")
# =============================================================

library(testthat)
library(shinytest2)    # UI automation
library(DBI)           # mock DB
library(digest)        # password hashing

# ---- 1. Mock database layer -------------------------------------
#   Isolate DB‑related code by connecting to an in‑memory SQLite
#   database that mirrors the schema of the production Postgres DB.
# -----------------------------------------------------------------

mock_pool <- dbConnect(RSQLite::SQLite(), ":memory:")

dbExecute(mock_pool, "
CREATE TABLE users (
  user_id INTEGER PRIMARY KEY AUTOINCREMENT,
  username TEXT UNIQUE,
  password_hash TEXT,
  email TEXT UNIQUE,
  created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
")

dbExecute(mock_pool, "
CREATE TABLE analysis_history (
  analysis_id INTEGER PRIMARY KEY AUTOINCREMENT,
  user_id     INTEGER,
  filename    TEXT,
  s3_key      TEXT,
  ci_columns  TEXT,
  he_columns  TEXT,
  log2_threshold FLOAT,
  pval_threshold FLOAT,
  created_at  TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
")

# ---- 2. Copies of server‑side functions ----------------
#   re‑source the pure functions we need for unit testing.
# -----------------------------------------------------------------

registerUser <- function(username, password, email, pool = mock_pool) {
  user_exists <- dbGetQuery(
    pool,
    sprintf("SELECT 1 FROM users WHERE username ='%s' OR email ='%s' LIMIT 1", username, email)
  )
  if (nrow(user_exists) > 0) return(FALSE)
  password_hash <- digest(password, algo = "sha256")
  dbExecute(pool, sprintf(
    "INSERT INTO users (username, password_hash, email) VALUES ('%s', '%s', '%s')",
    username, password_hash, email))
  TRUE
}

verifyLogin <- function(username, password, pool = mock_pool) {
  password_hash <- digest(password, algo = "sha256")
  res <- dbGetQuery(pool, sprintf(
    "SELECT user_id FROM users WHERE username = '%s' AND password_hash = '%s'",
    username, password_hash))
  if (nrow(res) == 1) res$user_id else NULL
}

impute_mean <- function(group_data) {
  for (i in 1:nrow(group_data)) {
    row <- group_data[i, ]
    row_mean <- mean(row, na.rm = TRUE)
    row[is.na(row)] <- row_mean
    group_data[i, ] <- row
  }
  return(group_data)
}

label_significance <- function(result, log2_threshold, pval_threshold) {
  result$Significance <- ifelse(
    result$adj.P.Val < pval_threshold &
      abs(result$logFC) > log2_threshold,
    ifelse(
      result$logFC > log2_threshold,
      "Upregulated",
      "Downregulated"
    ),
    "Not Significant"
  )
  return(result)
}
# -----------------------------------------------------------------
#                      WHITE‑BOX UNIT TESTS
# -----------------------------------------------------------------

test_that("registerUser() works for new user", {
  expect_true(registerUser("testing", "secret", "testing@example.com"))
})

test_that("registerUser() blocks duplicates", {
  expect_false(registerUser("testing", "another", "testing@example.com"))
})

test_that("verifyLogin() fails for wrong pwd", {
  expect_null(verifyLogin("testing", "wrong"))
})

test_that("verifyLogin() succeeds for correct pwd", {
  expect_true(!is.null(verifyLogin("testing", "secret")))
})

test_that("impute_mean fills NAs by row mean", {
  m <- matrix(c(1, NA, 3,
                4, 5, NA), nrow = 2, byrow = TRUE)
  imp <- impute_mean(m)
  expect_equal(imp[1,2], mean(c(1,3))) 
  expect_equal(imp[2,3], mean(c(4,5))) 
})

test_that("label_significance categorizes correctly", {
  df <- data.frame(
    logFC = c(2, -2, 0.5),
    adj.P.Val = c(0.01, 0.04, 0.2)
  )
  result <- label_significance(df, log2_threshold = 1, pval_threshold = 0.05)
  expect_equal(result$Significance, c("Upregulated", "Downregulated", "Not Significant"))
})

# -----------------------------------------------------------------
#                BLACK‑BOX / END‑TO‑END UI TESTS
# -----------------------------------------------------------------
#   These tests the entire Shiny app exactly as a user
#   would: upload file → configure inputs → run analysis → toggle
#   outputs → download files with shinytest2
# -----------------------------------------------------------------

# Path to the Shiny app root
app_dir <- normalizePath("..", mustWork = TRUE)  # tests/ -> root

# A small sample CSV (bundled with the app in inst/ or tests/testthat/)
#   12‑row, 12‑column anonymised dataset to keep CI fast.
#   Place file at tests/testthat/data/ProteinquantitationHE&CI.csv
sample_csv <- file.path(app_dir, "tests", "testthat", "data", "ProteinquantitationHE&CI.csv")

# Skip UI tests if sample file missing
if (file.exists(sample_csv)) {
  
  test_that("FULL WORKFLOW: upload → analysis → outputs", {
    app <- AppDriver$new(app_dir,
                         name   = "e2e_full_workflow",
                         height = 1600, width = 1200, seed = 987)
    
    # 1. LANDING PAGE  --------------------------------------------------
    app$click("goToLogin")
    app$expect_ui("#username")
    
    # 2. REGISTER + LOGIN  ---------------------------------------------
    app$click("show_register_form")
    app$set_inputs(reg_username         = "user")
    app$set_inputs(reg_email            = "user@example.com")
    app$set_inputs(reg_password         = "pass123")
    app$set_inputs(reg_confirm_password = "pass123")
    app$click("register")
    
    # back to login
    app$set_inputs(username = "user", password = "pass123")
    app$click("login")
    
    # 3. UPLOAD FILE + PARAMETERS  --------------------------------------
    app$expect_ui("#file")
    app$upload_file(file = sample_csv)
    
    app$set_inputs(protein_columns = "Accession")
    app$set_inputs(gene_columns    = "Gene Name")
    app$set_inputs(ci_columns      = "CI-1,CI-2,CI-3,CI-4")
    app$set_inputs(he_columns      = "HE-1,HE-2,HE-3,HE-4")
    app$set_inputs(log2_threshold  = "1")
    app$set_inputs(pval_threshold  = "0.05")
    
    app$click("start_analysis")
    
    # Wait until status message reads "Analysis complete!"
    app$wait_for(function() {
      app$get_value(output = "status_message") == "Analysis complete!"
    }, timeout = 10000)
    
    # 5. TOGGLE VOLCANO & HEATMAP  --------------------------------------
    app$click("toggle_volcano")
    app$expect_ui("#volcanoPlot")
    
    app$click("toggle_heatmap")
    app$expect_ui("#heatmapPlot")
    
    # 5. DOWNLOAD PNG AND CHECK FILE  -----------------------------------
    tmp_png <- tempfile(fileext = ".png")
    app$download("download_volcano", dest = tmp_png)
    expect_true(file.info(tmp_png)$size > 0)
    
    # 6. GENERATE & DOWNLOAD PDF REPORT  --------------------------------
    tmp_pdf <- tempfile(fileext = ".pdf")
    app$download("download_report", dest = tmp_pdf)
    expect_true(file.info(tmp_pdf)$size > 0)
    
    app$stop()
  })
  
}

# -----------------------------------------------------------------
#                    CODE COVERAGE HINT (optional)                 
# -----------------------------------------------------------------
# (Uncomment locally)  covr::report(covr::package_coverage())
