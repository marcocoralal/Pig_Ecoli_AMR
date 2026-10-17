# Antimicrobial non-susceptibility in porcine E. coli
# Multivariable logistic regression and sensitivity analyses
# Input: data/AMR_deidentified_analysis.csv
# Output: results/
# Run from the repository root: source("R/00_master_analysis.R")
options(stringsAsFactors = FALSE)
packages <- c("dplyr", "car", "logistf", "geepack", "igraph")
missing_pkgs <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Install required packages first: ", paste(missing_pkgs, collapse = ", "))
input_file <- "data/AMR_deidentified_analysis.csv"
if (!file.exists(input_file)) stop("Run from repository root; missing: ", input_file)
dir.create("results", showWarnings = FALSE, recursive = TRUE)
antibiotics <- c("sxt", "gm", "cip", "caz", "tet", "s", "chl", "cf", "ak", "azm", "ff", "amc")
analysis_vars <- paste0("ns_", antibiotics)
amr_data <- utils::read.csv(input_file, check.names=FALSE)
stopifnot(identical(names(amr_data), c("week_cluster", analysis_vars)),
          nrow(amr_data)==344L, ncol(amr_data)==13L,
          !anyNA(amr_data),
          all(vapply(amr_data[analysis_vars], function(x) is.numeric(x) && all(x %in% 0:1), logical(1))),
          is.numeric(amr_data$week_cluster),
          all(amr_data$week_cluster == as.integer(amr_data$week_cluster)),
          all(amr_data$week_cluster >= 1))
amr_data <- amr_data[order(amr_data$week_cluster), , drop=FALSE]
# Laboratory interpretation (I+R = 1, S = 0) must be confirmed independently.
utils::write.csv(data.frame(metric=c("ast_complete_rows", "week_clusters"),
  value=c(nrow(amr_data), length(unique(amr_data$week_cluster)))),
  "results/data_flow_audit.csv", row.names=FALSE)
utils::write.csv(do.call(rbind, lapply(analysis_vars, function(v)
  data.frame(variable=v, n_0=sum(amr_data[[v]]==0), n_1=sum(amr_data[[v]]==1)))),
  "results/phenotype_counts.csv", row.names=FALSE)

# Fit 12 outcome-specific GLMs (11 predictors each).
models <- list(); diagnostics <- list(); vif_rows <- list(); glm_rows <- list()
for (outcome in analysis_vars) {
  predictors <- setdiff(analysis_vars, outcome)
  form <- stats::reformulate(predictors, response = outcome)
  fit <- stats::glm(form, data = amr_data, family = stats::binomial(link = "logit"))
  null <- stats::glm(stats::reformulate("1", response = outcome),
                     data = amr_data, family = stats::binomial(link = "logit"))
  if (!isTRUE(fit$converged)) warning("GLM not converged: ", outcome)
  co <- summary(fit)$coefficients[predictors, , drop = FALSE]
  vifs <- car::vif(fit)
  events <- sum(amr_data[[outcome]] == 1L); nonevents <- sum(amr_data[[outcome]] == 0L)
  models[[outcome]] <- fit
  diagnostics[[outcome]] <- data.frame(outcome=outcome, n=stats::nobs(fit),
    events_1=events, events_0=nonevents, minority_events_per_predictor=min(events, nonevents)/length(predictors),
    converged=fit$converged, mcfadden_r2=1-as.numeric(stats::logLik(fit))/as.numeric(stats::logLik(null)),
    aic=stats::AIC(fit), max_vif=max(vifs))
  vif_rows[[outcome]] <- data.frame(outcome=outcome, predictor=names(vifs), vif=as.numeric(vifs))
  glm_rows[[outcome]] <- data.frame(outcome=outcome, predictor=predictors,
    beta=co[, "Estimate"], se=co[, "Std. Error"], p=co[, "Pr(>|z|)"], row.names=NULL)
}
model_diagnostics <- dplyr::bind_rows(diagnostics)
vif <- dplyr::bind_rows(vif_rows)
glm_results <- dplyr::bind_rows(glm_rows)
stopifnot(nrow(glm_results)==132L, !anyNA(glm_results$p))
glm_results$or <- exp(glm_results$beta)
glm_results$ci_low <- exp(glm_results$beta - stats::qnorm(.975)*glm_results$se)
glm_results$ci_high <- exp(glm_results$beta + stats::qnorm(.975)*glm_results$se)
glm_results$q <- stats::p.adjust(glm_results$p, method="BH")
glm_results$association <- ifelse(glm_results$or > 1, "Positive", "Negative")
glm_results$significant <- glm_results$q < .05
significant <- glm_results[glm_results$significant, , drop=FALSE]
utils::write.csv(model_diagnostics, "results/model_diagnostics.csv", row.names=FALSE)
utils::write.csv(vif, "results/vif_by_model.csv", row.names=FALSE)
utils::write.csv(glm_results, "results/GLM_132_FDR.csv", row.names=FALSE)

# Firth penalized logistic regression for ALL 12 models and all 132 tests.
# Profile penalized-likelihood CI and p-values; BH correction across 132 tests.
firth_rows <- list()
for (outcome in analysis_vars) {
  predictors <- setdiff(analysis_vars, outcome)
  fit <- logistf::logistf(stats::reformulate(predictors, response=outcome),
    data=amr_data, firth=TRUE, pl=TRUE,
    control=logistf::logistf.control(maxit=1000))
  firth_rows[[outcome]] <- data.frame(outcome=outcome, predictor=predictors,
    or_firth=exp(fit$coefficients[predictors]),
    ci_low_firth=exp(fit$ci.lower[predictors]),
    ci_high_firth=exp(fit$ci.upper[predictors]),
    p_firth=unname(fit$prob[predictors]), row.names=NULL)
}
firth <- dplyr::bind_rows(firth_rows)
stopifnot(nrow(firth)==132L, !anyNA(firth$p_firth))
firth$q_firth <- stats::p.adjust(firth$p_firth, method="BH")
firth$significant_firth <- firth$q_firth < .05
utils::write.csv(firth, "results/Firth_132_FDR.csv", row.names=FALSE)

# GEE: week clusters, exchangeable correlation, robust standard errors.
# geeglm requires contiguous observations within each cluster.
gee_rows <- list()
for (outcome in analysis_vars) {
  predictors <- setdiff(analysis_vars, outcome)
  fit <- geepack::geeglm(stats::reformulate(predictors, response=outcome),
    data=amr_data, id=week_cluster, family=stats::binomial(link="logit"),
    corstr="exchangeable")
  co <- summary(fit)$coefficients[predictors, , drop=FALSE]
  gee_rows[[outcome]] <- data.frame(outcome=outcome, predictor=predictors,
    beta_gee=co[, "Estimate"], se_gee=co[, "Std.err"], p_gee=co[, "Pr(>|W|)"], row.names=NULL)
}
gee_results <- dplyr::bind_rows(gee_rows)
stopifnot(nrow(gee_results)==132L, !anyNA(gee_results$p_gee))
gee_results$or_gee <- exp(gee_results$beta_gee)
gee_results$ci_low_gee <- exp(gee_results$beta_gee - stats::qnorm(.975)*gee_results$se_gee)
gee_results$ci_high_gee <- exp(gee_results$beta_gee + stats::qnorm(.975)*gee_results$se_gee)
gee_results$q_gee <- stats::p.adjust(gee_results$p_gee, method="BH")
gee_results$significant_gee <- gee_results$q_gee < .05
utils::write.csv(gee_results, "results/GEE_week_132_FDR.csv", row.names=FALSE)

# Cross-method comparison.
comparison <- dplyr::left_join(glm_results, firth, by=c("outcome", "predictor")) |>
  dplyr::left_join(gee_results, by=c("outcome", "predictor"))
comparison$significant_all_three <- comparison$significant & comparison$significant_firth & comparison$significant_gee
utils::write.csv(comparison, "results/GLM_Firth_GEE_comparison.csv", row.names=FALSE)

# Table 4: all and only BH-significant GLM associations.
labels <- c(sxt="Trimethoprim-sulfamethoxazole", gm="Gentamicin", cip="Ciprofloxacin",
  caz="Ceftazidime", tet="Tetracycline", s="Streptomycin", chl="Chloramphenicol",
  cf="Cefoxitin", ak="Amikacin", azm="Azithromycin", ff="Florfenicol",
  amc="Amoxicillin-clavulanate")
label_var <- function(v) unname(labels[sub("^ns_", "", v)])
table4 <- dplyr::left_join(significant, model_diagnostics[, c("outcome", "mcfadden_r2")], by="outcome")
table4 <- table4[order(table4$outcome, table4$q), , drop=FALSE]
table4$outcome_name <- label_var(table4$outcome)
table4$predictor_name <- label_var(table4$predictor)
utils::write.csv(table4, "results/Table4_numeric.csv", row.names=FALSE)
format_p <- function(x) ifelse(x < .001, "<0.001", sprintf("%.3f", x))
table4_print <- data.frame(
  `Non-susceptibility outcome`=table4$outcome_name,
  `Non-susceptibility predictor`=table4$predictor_name,
  `Adjusted OR (95% CI)`=sprintf("%.2f (%.2f–%.2f)", table4$or, table4$ci_low, table4$ci_high),
  `p-value`=format_p(table4$p), `FDR-adjusted q-value`=format_p(table4$q),
  `McFadden pseudo-R2`=sprintf("%.3f", table4$mcfadden_r2),
  `Association`=table4$association, check.names=FALSE)
utils::write.csv(table4_print, "results/Table4_publication.csv", row.names=FALSE, fileEncoding="UTF-8")

# Figure 1: directed edges are predictor -> outcome in SEPARATE models, not causality.
edges <- data.frame(from=sub("^ns_", "", significant$predictor),
  to=sub("^ns_", "", significant$outcome),
  association=significant$association, or=significant$or, q=significant$q)
vertices <- data.frame(name=antibiotics)
g <- igraph::graph_from_data_frame(edges, directed=TRUE, vertices=vertices)
stopifnot(igraph::vcount(g)==12L, igraph::ecount(g)==nrow(significant))
ends_mat <- igraph::ends(g, igraph::E(g), names=TRUE)
curvature <- vapply(seq_len(nrow(ends_mat)), function(i) {
  if (any(ends_mat[,1]==ends_mat[i,2] & ends_mat[,2]==ends_mat[i,1])) .22 else 0
}, numeric(1))
order_circle <- c("sxt", "chl", "tet", "cip", "gm", "ak", "s", "ff", "caz", "azm", "amc", "cf")
a <- pi/2 - 2*pi*(seq_along(order_circle)-1)/length(order_circle)
xy <- cbind(2.2*cos(a), 2.2*sin(a)); rownames(xy) <- order_circle
xy <- xy[igraph::V(g)$name, , drop=FALSE]
color <- ifelse(igraph::E(g)$association=="Positive", "#0072B2", "#D55E00")
utils::write.csv(edges, "results/Figure1_edges.csv", row.names=FALSE)

# Quantify both directed total degree and distinct undirected neighbors.
neighbor_counts <- vapply(igraph::V(g), function(v) length(igraph::neighbors(g, v, mode="all")), integer(1))
node_metrics <- data.frame(antimicrobial=igraph::V(g)$name,
  total_degree=igraph::degree(g, mode="all"),
  distinct_neighbors=neighbor_counts)
utils::write.csv(node_metrics, "results/Figure1_node_connectivity.csv", row.names=FALSE)

grDevices::pdf("results/Figure1_network.pdf", width=12, height=11, useDingbats=FALSE)
layout(matrix(c(1,2), ncol=1), heights=c(8.5,1.5))
par(mar=c(1,2,3,2))
plot(g, layout=xy, vertex.size=27, vertex.color="white", vertex.frame.color="#333333",
  vertex.label=toupper(igraph::V(g)$name), vertex.label.cex=.85, vertex.label.font=2,
  edge.color=color, edge.width=1.8, edge.curved=curvature, edge.arrow.size=.5,
  edge.lty=ifelse(igraph::E(g)$association=="Positive",1,2), margin=.08,
  main="Antimicrobial non-susceptibility associations")
par(mar=c(0,0,0,0)); plot.new()
legend("center", legend=c("Positive association (aOR > 1)","Negative association (aOR < 1)",
  "Arrow: predictor to outcome (not causality)"), col=c("#0072B2","#D55E00","#444444"),
  lty=c(1,2,1), lwd=c(2.5,2.5,1.5), bty="o", bg="white", cex=.9)
grDevices::dev.off()

# Validation of analytical outputs against reference counts.
checks <- c(n_complete=nrow(amr_data), tests=nrow(glm_results), nominal=sum(glm_results$p<.05),
  fdr_glm=sum(glm_results$significant), positive=sum(significant$association=="Positive"),
  negative=sum(significant$association=="Negative"),
  fdr_gee=sum(gee_results$significant_gee), glm_gee_overlap=sum(comparison$significant & comparison$significant_gee),
  glm_firth_overlap=sum(comparison$significant & comparison$significant_firth),
  all_three=sum(comparison$significant_all_three))
expected <- c(n_complete=344, tests=132, nominal=40, fdr_glm=25, positive=22,
  negative=3, fdr_gee=19, glm_gee_overlap=17, glm_firth_overlap=25, all_three=17)
utils::write.csv(data.frame(metric=names(checks), observed=unname(checks), expected=unname(expected)),
  "results/manuscript_validation.csv", row.names=FALSE)
print(checks)
if (!identical(unname(as.numeric(checks)), unname(as.numeric(expected)))) {
  stop("Validation failed. See results/manuscript_validation.csv.")
}
capture.output(sessionInfo(), file="results/sessionInfo.txt")
cat("Analysis completed; validation checks passed.\n")
