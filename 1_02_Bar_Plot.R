out_dir <- 'figures/figures/ab_fig_2_cMCL_and_nnMCL_10.12.2025' # Changed output directory to differentiate
load(file.path(out_dir, 'ev.100k_multi_samples_cMCL_nnMCL_10.21.2025.rda'))

all_ev_lists <- list(
    RBL = EV.rbl,
    LCL = EV.lcl,
    GCBC = EV.gcbc,
    MBC = EV.mbc,
    NBC = EV.nbc,
    PC = EV.pc,
    cMCL = EV.cMCL,
    nnMCL = EV.nnMCL
)

sample_names <- names(all_ev_lists)


# generate compartment a and b figures
chroms <- paste0('chr',c(1:22,'X'))

for(chr in chroms){
  # Adjust height to accommodate all 6 samples
  pdf(file = file.path(out_dir, paste0(chr, '.tad.ab.pdf')), width=10, height=14)
  par(mfrow=c(8,1), font.lab=2, cex.lab=1.2, lty = 0, mar=c(0,4,0,0), oma = c(0,0,0,0), mgp = c(2, 1, 0), xaxs='i', yaxs='i')
  barplot(EV.rbl[[chr]], col = ifelse(EV.rbl[[chr]] > 0,"red","blue"),space = 0,ylab='RBL EV',ylim=c(-0.12,0.12))
  barplot(EV.lcl[[chr]], col = ifelse(EV.lcl[[chr]] > 0,"red","blue"),space = 0,ylab='LCL EV',ylim=c(-0.12,0.12))
  barplot(EV.gcbc[[chr]], col = ifelse(EV.gcbc[[chr]] > 0,"red","blue"),space = 0,ylab='GCBC EV',ylim=c(-0.12,0.12))
  barplot(EV.mbc[[chr]], col = ifelse(EV.mbc[[chr]] > 0,"red","blue"),space = 0,ylab='MBC EV',ylim=c(-0.12,0.12))
  barplot(EV.nbc[[chr]], col = ifelse(EV.nbc[[chr]] > 0,"red","blue"),space = 0,ylab='NBC EV',ylim=c(-0.12,0.12))
  barplot(EV.pc[[chr]], col = ifelse(EV.pc[[chr]] > 0,"red","blue"),space = 0,ylab='PC EV',ylim=c(-0.12,0.12))
  barplot(EV.cMCL[[chr]], col = ifelse(EV.cMCL[[chr]] > 0,"red","blue"),space = 0,ylab='cMCL EV',ylim=c(-0.12,0.12))
  barplot(EV.nnMCL[[chr]], col = ifelse(EV.nnMCL[[chr]] > 0,"red","blue"),space = 0,ylab='nnMCL EV',ylim=c(-0.12,0.12))
  
  
  dev.off()
}

# 1.0 check range of EV
out_file <- file.path(out_dir, "EV_range_summary.txt")
# Open a connection for writing
sink(out_file)
cat("EV Range Summary (min and max values per sample)\n")
cat("------------------------------------------------\n\n")
for (sample_name in names(all_ev_lists)) {
  cat(sample_name, "\n")
  for (chr in names(all_ev_lists[[sample_name]])) {
    ev <- all_ev_lists[[sample_name]][[chr]]
    cat("  ", chr, ": ",
        round(min(ev, na.rm = TRUE), 4), " to ",
        round(max(ev, na.rm = TRUE), 4), "\n")
  }
  cat("\n")
}
# Close the connection
sink()
cat(" EV ranges written to:", out_file, "\n")

# 2.1 EV histogram (for each chromosome in sample individually)
pdf(file = file.path(out_dir, 'ev.hist_all_samples.pdf'), width = 9, height = length(chroms) * 2) # Adjusted height
par(mfrow=c(length(chroms), 6), font.lab=2, cex.lab=1.2) # Adjusted columns to 6
x_limit <- c(-0.15, 0.15)
pdf(file = file.path(out_dir, 'ev.hist_all_samples_fixed.pdf'),
    width = 16, height = length(chroms) * 8)
par(mfrow = c(length(chroms), 8), font.lab = 2, cex.lab = 1.2,
    mar = c(3,3,2,1), oma = c(0,0,0,0))
for (i in seq_along(EV.rbl)) {
  hist(EV.rbl[[i]], n = 50, main = paste0(names(EV.rbl)[i], ': RBL'),
       xlab = "EV Value", xlim = x_limit)
  hist(EV.lcl[[i]], n = 50, main = paste0(names(EV.lcl)[i], ': LCL'),
       xlab = "EV Value", xlim = x_limit)
  hist(EV.gcbc[[i]], n = 50, main = paste0(names(EV.gcbc)[i], ': GCBC'),
       xlab = "EV Value", xlim = x_limit)
  hist(EV.mbc[[i]], n = 50, main = paste0(names(EV.mbc)[i], ': MBC'),
       xlab = "EV Value", xlim = x_limit)
  hist(EV.nbc[[i]], n = 50, main = paste0(names(EV.nbc)[i], ': NBC'),
       xlab = "EV Value", xlim = x_limit)
  hist(EV.pc[[i]], n = 50, main = paste0(names(EV.pc)[i], ': PC'),
       xlab = "EV Value", xlim = x_limit)
  hist(EV.cMCL[[i]], n = 50, main = paste0(names(EV.cMCL)[i], ': cMCL'),
       xlab = "EV Value", xlim = x_limit)
  hist(EV.nnMCL[[i]], n = 50, main = paste0(names(EV.nnMCL)[i], ': nnMCL'),
       xlab = "EV Value", xlim = x_limit)
}
dev.off()

#2.2 EV histogram single sample (for each sample containing all chromosomes)
# Fixed axes
x_limit <- c(-0.15, 0.15)
breaks  <- seq(x_limit[1], x_limit[2], length.out = 51)  # 50 bins
# Compute a global y-limit across ALL samples (all chromosomes)
y_max <- 0
for (sn in sample_names) {
  ev_vals <- unlist(all_ev_lists[[sn]], use.names = FALSE)
  ev_vals <- ev_vals[!is.na(ev_vals)]
  if (length(ev_vals) > 0) {
    h <- hist(ev_vals, breaks = breaks, plot = FALSE)
    y_max <- max(y_max, h$counts, na.rm = TRUE)
  }
}
ylim_fixed <- c(0, ceiling(y_max * 1.05))  # add 5% headroom
# Plot 8 samples as a 2x4 grid (4 per row) on one page
pdf(file = file.path(out_dir, 'ev.hist_all_samples_fixed_4perrow_50bins.pdf'),
    width = 14, height = 7.5)
par(mfrow = c(2, 4), mar = c(4,4,3,1), oma = c(0,0,0,0), font.lab = 2, cex.lab = 1.2)
for (sn in sample_names) {
  ev_vals <- unlist(all_ev_lists[[sn]], use.names = FALSE)
  ev_vals <- ev_vals[!is.na(ev_vals)]
  if (length(ev_vals) == 0) {
    plot(NA, xlim = x_limit, ylim = ylim_fixed,
         main = paste(sn, "(no data)"), xlab = "EV Value", ylab = "Frequency")
  } else {
    hist(ev_vals, breaks = breaks, xlim = x_limit, ylim = ylim_fixed,
         main = paste(sn, "EV Distribution"), xlab = "EV Value", ylab = "Frequency",
         col = "gray70", border = "gray30")
    abline(v = 0, lty = 2, lwd = 2, col = "red")
  }
}
dev.off()
cat(" Saved:", file.path(out_dir, 'ev.hist_all_samples_fixed_4perrow_50bins.pdf'), "\n")

# 3.0 EV Histograms: Pairwise Comparison (where signs differ, one PDF per chromosome) ---
xlim_fixed <- c(-0.15, 0.15)
breaks     <- seq(xlim_fixed[1], xlim_fixed[2], length.out = 51)  # 50 bins
# 3-per-row layout
cols_per_page  <- 3          # ev1, ev2, difference
rows_per_page  <- 10         # 10 pairs per page => 30 panels per page
page_width_in  <- 15         # ~5 in per panel horizontally
page_height_in <- rows_per_page * 3.5
# 2) Compute ONE global y-limit across ALL chromosomes & ALL pairs
global_ymax <- 0
have_any    <- FALSE
for (chr_name in chroms) {
  for (s1_idx in 1:(length(sample_names) - 1)) {
    for (s2_idx in (s1_idx + 1):length(sample_names)) {
      s1 <- sample_names[s1_idx]; s2 <- sample_names[s2_idx]
      ev1 <- all_ev_lists[[s1]][[chr_name]]
      ev2 <- all_ev_lists[[s2]][[chr_name]]
      if (is.null(ev1) || is.null(ev2)) next
      ok <- which(is.finite(ev1) & is.finite(ev2))
      if (length(ok) == 0) next
      ds <- sign(ev1[ok]) != sign(ev2[ok])
      if (!any(ds)) next
      have_any <- TRUE
      v1 <- ev1[ok][ds]
      v2 <- ev2[ok][ds]
      vd <- v2 - v1
      # IMPORTANT: counts are computed using the same fixed breaks
      h1 <- hist(v1, breaks = breaks, plot = FALSE)
      h2 <- hist(v2, breaks = breaks, plot = FALSE)
      hD <- hist(vd, breaks = breaks, plot = FALSE)
      global_ymax <- max(global_ymax, h1$counts, h2$counts, hD$counts, na.rm = TRUE)
    }
  }
}
if (!have_any) {
  message("No diff-sign bins found anywhere. Writing small info PDFs per chromosome...")
  for (chr_name in chroms) {
    pdf(file = file.path(out_dir, paste0(chr_name, '.ev.hist_pairwise_diff_signs.pdf')),
        width = 8, height = 3)
    par(mfrow = c(1,1), mar = c(2,2,2,1), oma = c(0,0,0,0))
    plot.new(); title(main = paste0(chr_name, ": No diff-sign bins across pairs"))
    dev.off()
  }
} else {
  ylim_fixed <- c(0, ceiling(global_ymax * 1.05))  # headroom
  # 3) Plot per chromosome, using uniform xlim/ylim everywhere
  for (chr_name in chroms) {
    pdf(file = file.path(out_dir, paste0(chr_name, '.ev.hist_pairwise_diff_signs.pdf')),
        width = page_width_in, height = page_height_in)
    par(mfrow = c(rows_per_page, cols_per_page),
        font.lab = 2, cex.lab = 1.0,
        mar = c(4,4,3,1), oma = c(0,0,0,0))
    any_pair_plotted <- FALSE
    for (s1_idx in 1:(length(sample_names) - 1)) {
      for (s2_idx in (s1_idx + 1):length(sample_names)) {
        s1 <- sample_names[s1_idx]; s2 <- sample_names[s2_idx]
        ev1 <- all_ev_lists[[s1]][[chr_name]]
        ev2 <- all_ev_lists[[s2]][[chr_name]]
        if (is.null(ev1) || is.null(ev2)) next
        ok <- which(is.finite(ev1) & is.finite(ev2))
        if (length(ok) == 0) next
        ds <- sign(ev1[ok]) != sign(ev2[ok])
        if (any(ds)) {
          any_pair_plotted <- TRUE
          v1 <- ev1[ok][ds]
          v2 <- ev2[ok][ds]
          vd <- v2 - v1
          # ev1
          hist(v1, breaks = breaks, xlim = xlim_fixed, ylim = ylim_fixed,
               main = paste0(chr_name, ': ', s1, ' (vs ', s2, ')'),
               xlab = "EV Value", ylab = "Frequency",
               col = "gray70", border = "gray30")
          abline(v = 0, col = "red", lty = 2, lwd = 2)
          # ev2
          hist(v2, breaks = breaks, xlim = xlim_fixed, ylim = ylim_fixed,
               main = paste0(chr_name, ': ', s2, ' (vs ', s1, ')'),
               xlab = "EV Value", ylab = "Frequency",
               col = "gray70", border = "gray30")
          abline(v = 0, col = "red", lty = 2, lwd = 2)
          # difference
          hist(vd, breaks = breaks, xlim = xlim_fixed, ylim = ylim_fixed,
               main = paste0(chr_name, ': ', s2, ' - ', s1),
               xlab = "EV Difference", ylab = "Frequency",
               col = "gray70", border = "gray30")
          abline(v = 0, col = "red", lty = 2, lwd = 2)
        } else {
          # consume a full row (3 panels) for pagination consistency
          plot.new(); title(paste0(chr_name, ': ', s1, ' - ', s2, ' (No diff signs)'))
          plot.new(); plot.new()
        }
      }
    }
    if (!any_pair_plotted) {
      # In case this chromosome specifically has no diff-sign pairs, show a single notice
      par(mfrow = c(1,1))
      plot.new(); title(main = paste0(chr_name, ": No diff-sign bins across pairs"))
    }
    dev.off()
  }
}


# 4.0 MA plot: pairwise comparison across all samples
# This will now generate a separate PDF for each chromosome.
# Each PDF will contain (N * (N-1) / 2) plots for N samples.
# For 6 samples, this is 15 plots per chromosome PDF.
for(chr_name in chroms){ # Loop by chromosome name for consistent access
  # Adjust PDF dimensions to make individual plots roughly square in a 5x3 grid
  # A width of 15 inches for 3 columns gives 5 inches per column.
  # A height of 25 inches for 5 rows gives 5 inches per row.
  # This results in roughly 5x5 inch plots, which are square.
  pdf(file = file.path(out_dir, paste0(chr_name, '.ev.ma_pairwise_plots.pdf')),
      width = 15, # 3 columns * ~5 inches/column
      height = 25) # 5 rows * ~5 inches/row
  
  # Set mfrow to arrange 5 rows and 3 columns of plots on the page
  par(mfrow=c(rows_per_page, plots_per_row),
      font.lab=2, cex.lab=1.2,
      mar=c(4, 4, 3, 1), # Added margins (bottom, left, top, right)
      oma=c(0, 0, 0, 0), # Outer margins
      mgp = c(2.5, 0.8, 0), # Adjusted axis title and label positions
      xaxs='i', yaxs='i' # Ensures plotting region extends to axis limits
  ) # <--- ADDED THIS CLOSING PARENTHESIS FOR par()
  
  # The rest of your code block follows below, unchanged
  for (s1_idx in 1:(length(sample_names) - 1)) {
    for (s2_idx in (s1_idx + 1):length(sample_names)) {
      sample1_name <- sample_names[s1_idx]
      sample2_name <- sample_names[s2_idx]
      
      ev1 <- all_ev_lists[[sample1_name]][[chr_name]]
      ev2 <- all_ev_lists[[sample2_name]][[chr_name]]
      
      # Ensure both EV vectors are of the same length and contain valid numbers
      common_idx <- which(!is.na(ev1) & !is.na(ev2))
      
      if (length(common_idx) > 0) {
        EV_diff <- ev2[common_idx] - ev1[common_idx]
        EV_mean <- (ev1[common_idx] + ev2[common_idx]) / 2
        
        # Determine appropriate x and y limits for the plot
        x_limit <- c(-0.1, 0.1) # Mean EV ranges from -0.1 to 0.1
        y_limit <- c(-0.2, 0.2) # Diff EV ranges from -0.2 to 0.2
        
        plot(EV_mean, EV_diff,
             pch=16, cex=0.25, # Smaller points for dense plots
             main=paste0(chr_name, ": ", sample2_name, " vs ", sample1_name),
             xlab = "Mean EV", ylab = "EV Difference",
             xlim = x_limit, # Set consistent x-axis limits
             ylim = y_limit, # Set consistent y-axis limits
             cex.main=1.0, cex.lab=1.0, cex.axis=0.8) # Adjust font sizes within plot
        abline(h=0,lty="dashed", col="red") # Make horizontal line red for visibility
      } else {
        # Plot placeholder if no common data points exist
        plot(NA, xlim=c(min(x_limit), max(x_limit)), ylim=c(min(y_limit), max(y_limit)), # Use consistent limits for empty plots too
             type="n", xlab="", ylab="",
             main=paste0(chr_name, ": ", sample2_name, " vs ", sample1_name, " (No common data)"),
             cex.main=1.0)
      }
    }
  }
  dev.off() # Close the PDF for the current chromosome
}

# sign flipped: all samples (density plot of differences for sign-flipped bins)
# This will now include differences for ALL pairwise sample combinations where signs flip.
EV.flip_all_pairwise <- c()
for (s1_idx in 1:(length(sample_names) - 1)) {
  for (s2_idx in (s1_idx + 1):length(sample_names)) {
    sample1_name <- sample_names[s1_idx]
    sample2_name <- sample_names[s2_idx]
    
    for(i in seq_along(chroms)){
      chr <- chroms[i]
      ev1 <- all_ev_lists[[sample1_name]][[chr]]
      ev2 <- all_ev_lists[[sample2_name]][[chr]]
      
      idx <- sign(ev1) != sign(ev2)
      EV_diff_flipped <- ev2[idx] - ev1[idx]
      EV.flip_all_pairwise <- c(EV.flip_all_pairwise, EV_diff_flipped)
    }
  }
}

# The density plot will aggregate all these pairwise differences
ev_flip_density_pairwise <- density(EV.flip_all_pairwise, na.rm = TRUE)
pdf(file.path(out_dir, 'ev.flip_density_pairwise.pdf'),6,6)
plot(ev_flip_density_pairwise, main = "Density of EV Differences for Pairwise Sign-Flipped Bins (All Samples)")
abline(v=0,lty="dashed")
dev.off()


# 5.0 chr17 染色体示例绘图
chr = 'chr17'
idx = as.integer((3.35e+7 / 1e+5 + 1):(4.9e+7 / 1e+5))
EV.rbl.chr = -EV.rbl[[chr]][idx]
EV.lcl.chr = -EV.lcl[[chr]][idx]
EV.gcbc.chr = -EV.gcbc[[chr]][idx]
EV.mbc.chr = -EV.mbc[[chr]][idx]
EV.nbc.chr = -EV.nbc[[chr]][idx]
EV.pc.chr = -EV.pc[[chr]][idx]
EV.cMCL.chr = -EV.cMCL[[chr]][idx]
EV.nnMCL.chr = -EV.nnMCL[[chr]][idx]


pdf(file = file.path(out_dir, paste0(chr, '.tad.ab.examp.pdf')), width=6, height=12)  # 高度加大以容纳6个样本
par(mfrow = c(6,1), font.lab=2, cex.lab=1.2, lty=0, mar=c(0,4,0,0), oma=c(0,0,0,0), mgp=c(2,1,0), xaxs='i', yaxs='i')
barplot(EV.rbl.chr, col=ifelse(EV.rbl.chr > 0, "red", "blue"), space=0, ylab='RBL EV', ylim=c(-0.05, 0.05))
barplot(EV.lcl.chr, col=ifelse(EV.lcl.chr > 0, "red", "blue"), space=0, ylab='LCL EV', ylim=c(-0.05, 0.05))
barplot(EV.gcbc.chr, col=ifelse(EV.gcbc.chr > 0, "red", "blue"), space=0, ylab='GCBC EV', ylim=c(-0.05, 0.05))
barplot(EV.mbc.chr, col=ifelse(EV.mbc.chr > 0, "red", "blue"), space=0, ylab='MBC EV', ylim=c(-0.05, 0.05))
barplot(EV.nbc.chr, col=ifelse(EV.nbc.chr > 0, "red", "blue"), space=0, ylab='NBC EV', ylim=c(-0.05, 0.05))
barplot(EV.pc.chr,  col=ifelse(EV.pc.chr  > 0, "red", "blue"), space=0, ylab='PC EV',  ylim=c(-0.05, 0.05))
barplot(EV.cMCL.chr,  col=ifelse(EV.cMCL.chr  > 0, "red", "blue"), space=0, ylab='cMCL EV',  ylim=c(-0.05, 0.05))
barplot(EV.nnMCL.chr,  col=ifelse(EV.nnMCL.chr  > 0, "red", "blue"), space=0, ylab='nnMCL EV',  ylim=c(-0.05, 0.05))

dev.off()

# 6.0 chr3 染色体示例绘图
chr = 'chr3'
idx = as.integer((1.2e+8 / 1e+5 + 1):(1.44e+8 / 1e+5))
EV.rbl.chr = -EV.rbl[[chr]][idx]
EV.lcl.chr = -EV.lcl[[chr]][idx]
EV.gcbc.chr = -EV.gcbc[[chr]][idx]
EV.mbc.chr = -EV.mbc[[chr]][idx]
EV.nbc.chr = -EV.nbc[[chr]][idx]
EV.pc.chr = -EV.pc[[chr]][idx]
EV.cMCL.chr = -EV.cMCL[[chr]][idx]
EV.nnMCL.chr = -EV.nnMCL[[chr]][idx]


pdf(file = file.path(out_dir, paste0(chr, '.tad.ab.examp.pdf')), width=10, height=12)
par(mfrow = c(6,1), font.lab=2, cex.lab=1.2, lty=0, mar=c(0,4,0,0), oma=c(0,0,0,0), mgp=c(2,1,0), xaxs='i', yaxs='i')
barplot(EV.rbl.chr, col=ifelse(EV.rbl.chr > 0, "red", "blue"), space=0, ylab='RBL EV', ylim=c(-0.05, 0.05))
barplot(EV.lcl.chr, col=ifelse(EV.lcl.chr > 0, "red", "blue"), space=0, ylab='LCL EV', ylim=c(-0.05, 0.05))
barplot(EV.gcbc.chr, col=ifelse(EV.gcbc.chr > 0, "red", "blue"), space=0, ylab='GCBC EV', ylim=c(-0.05, 0.05))
barplot(EV.mbc.chr, col=ifelse(EV.mbc.chr > 0, "red", "blue"), space=0, ylab='MBC EV', ylim=c(-0.05, 0.05))
barplot(EV.nbc.chr, col=ifelse(EV.nbc.chr > 0, "red", "blue"), space=0, ylab='NBC EV', ylim=c(-0.05, 0.05))
barplot(EV.pc.chr,  col=ifelse(EV.pc.chr  > 0, "red", "blue"), space=0, ylab='PC EV',  ylim=c(-0.05, 0.05))
barplot(EV.cMCL.chr,  col=ifelse(EV.cMCL.chr  > 0, "red", "blue"), space=0, ylab='cMCL EV',  ylim=c(-0.05, 0.05))
barplot(EV.nnMCL.chr,  col=ifelse(EV.nnMCL.chr  > 0, "red", "blue"), space=0, ylab='nnMCL EV',  ylim=c(-0.05, 0.05))

dev.off()
