benchmark_write_runtime_plot <- function(summary_results, output_dir,
                                         metric = "mean_fit_time_sec") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("Package 'ggplot2' is not installed; skipping runtime plot.")
    return(character())
  }
  if (!metric %in% names(summary_results)) {
    stop("Metric '", metric, "' is not present in summary_results.")
  }

  plot_dt <- data.table::copy(summary_results[!is.na(get(metric))])
  if (!nrow(plot_dt)) {
    warning("No successful benchmark rows available for plotting.")
    return(character())
  }
  plot_dt[, x_value := dimension_value]
  plot_dt[, plot_metric := pmax(get(metric), .Machine$double.eps)]
  plot_dt[, n_points := .N, by = .(sweep, package)]

  p <- ggplot2::ggplot(
    plot_dt,
    ggplot2::aes(x = x_value, y = plot_metric, color = package, group = package)
  ) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_wrap(~ sweep, scales = "free", labeller = ggplot2::as_labeller(c(
      N = "Scaling N (fixed J)",
      J = "Scaling J (fixed N)",
      full_grid = "Full N x J grid"
    ))) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = "Scaling dimension",
      y = "Average runtime (seconds, log10 scale)",
      color = "Package"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )
  if (any(plot_dt$n_points > 1L)) {
    p <- p + ggplot2::geom_line(data = plot_dt[n_points > 1L], linewidth = 0.7)
  }

  png_path <- file.path(output_dir, "runtime_scaling.png")
  pdf_path <- file.path(output_dir, "runtime_scaling.pdf")
  ggplot2::ggsave(png_path, p, width = 8, height = 5, dpi = 160)
  ggplot2::ggsave(pdf_path, p, width = 8, height = 5)
  c(png_path, pdf_path)
}
