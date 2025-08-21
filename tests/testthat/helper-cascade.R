# Helper file for Cascade tests
# Use this space for shared utilities or global skips if needed.

# Muffle only warnings whose message matches a regex pattern,
# keeping other warnings visible so we don't hide real issues.
muffle_warnings_matching <- function(expr, patterns) {
  pats <- as.character(patterns)
  withCallingHandlers(expr, warning = function(w) {
    msg <- conditionMessage(w)
    if (any(vapply(pats, function(p) grepl(p, msg, ignore.case = TRUE), logical(1)))) {
      invokeRestart("muffleWarning")
    }
  })
}
