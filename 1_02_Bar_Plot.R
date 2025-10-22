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

# EV histogram (for each sample individually)
pdf(file = file.path(out_dir, 'ev.hist_all_samples.pdf'), width = 9, height = length(chroms) * 2) # Adjusted height
par(mfrow=c(length(chroms), 6), font.lab=2, cex.lab=1.2) # Adjusted columns to 6
for(i in seq_along(EV.rbl)){
  hist(EV.rbl[[i]],n=50,main=paste0(names(EV.rbl)[i],': RBL'))
  hist(EV.lcl[[i]],n=50,main=paste0(names(EV.lcl)[i],': LCL'))
  hist(EV.gcbc[[i]],n=50,main=paste0(names(EV.gcbc)[i],': GCBC'))
  hist(EV.mbc[[i]],n=50,main=paste0(names(EV.mbc)[i],': MBC'))
  hist(EV.nbc[[i]],n=50,main=paste0(names(EV.nbc)[i],': NBC'))
  hist(EV.pc[[i]],n=50,main=paste0(names(EV.pc)[i],': PC'))
  hist(EV.cMCL[[i]],n=50,main=paste0(names(EV.cMCL)[i],': cMCL'))
  hist(EV.nnMCL[[i]],n=50,main=paste0(names(EV.nnMCL)[i],': nnMCL'))
}
dev.off()


# --- EV Histograms: Pairwise Comparison (where signs differ, one PDF per chromosome) ---
num_pairs <- length(sample_names) * (length(sample_names) - 1) / 2
for(chr_name in chroms){ # Loop by chromosome name for consistent access
  pdf(file = file.path(out_dir, paste0(chr_name, '.ev.hist_pairwise_diff_signs.pdf')), width = 18, height = 24)
  par(mfrow=c(num_pairs, 3), font.lab=2, cex.lab=1.2, mar=c(4,4,3,1), oma = c(0,0,0,0)) # Explicit margins
  
  for (s1_idx in 1:(length(sample_names) - 1)) {
    for (s2_idx in (s1_idx + 1):length(sample_names)) {
      sample1_name <- sample_names[s1_idx]
      sample2_name <- sample_names[s2_idx]
      
      ev1 <- all_ev_lists[[sample1_name]][[chr_name]] # Use chr_name for access
      ev2 <- all_ev_lists[[sample2_name]][[chr_name]] # Use chr_name for access
      
      idx_diff_sign <- sign(ev1) != sign(ev2)
      
      if (sum(idx_diff_sign, na.rm = TRUE) > 0) {
        hist(ev1[idx_diff_sign], n=50, main=paste0(chr_name, ': ', sample1_name, ' (VS ', sample2_name, ' Diff Sign)'), xlab="EV Value")
        hist(ev2[idx_diff_sign], n=50, main=paste0(chr_name, ': ', sample2_name, ' (VS ', sample1_name, ' Diff Sign)'), xlab="EV Value")
        hist(ev2[idx_diff_sign] - ev1[idx_diff_sign], n=50, main=paste0(chr_name, ': ', sample2_name, '-', sample1_name, ' (Diff Bins)'), xlab="EV Difference")
      } else {
        plot(NA, xlim=c(0,1), ylim=c(0,1), type="n", xlab="", ylab="", main=paste0(chr_name, ': ', sample1_name, '-', sample2_name, ' (No diff signs)'))
        plot(NA, xlim=c(0,1), ylim=c(0,1), type="n", xlab="", ylab="", main="")
        plot(NA, xlim=c(0,1), ylim=c(0,1), type="n", xlab="", ylab="", main="") # <--- This was the misplaced line
      }
    }
  }
  dev.off()
}

# EV histogram updated 10.22.2025 #
# --- EV Histograms: Pairwise Comparison (where signs differ, one PDF per chromosome) ---
num_pairs <- length(sample_names) * (length(sample_names) - 1) / 2
rows_per_page <- 10  # number of pairs (rows) per page; device will auto-advance pages

for (chr_name in chroms) {
  pdf(file = file.path(out_dir, paste0(chr_name, '.ev.hist_pairwise_diff_signs.pdf')),
      width = 18, height = 25)  # comfortable page size
  par(mfrow = c(rows_per_page, 3), font.lab = 2, cex.lab = 1.0,
      mar = c(4,4,3,1), oma = c(0,0,0,0))

  pair_count <- 0
  for (s1_idx in 1:(length(sample_names) - 1)) {
    for (s2_idx in (s1_idx + 1):length(sample_names)) {
      pair_count <- pair_count + 1

      sample1_name <- sample_names[s1_idx]
      sample2_name <- sample_names[s2_idx]

      ev1 <- all_ev_lists[[sample1_name]][[chr_name]]
      ev2 <- all_ev_lists[[sample2_name]][[chr_name]]

      # robust sign compare: only where both are non-NA
      idx <- which(!is.na(ev1) & !is.na(ev2))
      has_data <- length(idx) > 0

      if (has_data) {
        ds <- sign(ev1[idx]) != sign(ev2[idx])
        if (any(ds, na.rm = TRUE)) {
          hist(ev1[idx][ds], n = 50,
               main = paste0(chr_name, ': ', sample1_name, ' (vs ', sample2_name, ' diff sign)'),
               xlab = "EV Value")
          hist(ev2[idx][ds], n = 50,
               main = paste0(chr_name, ': ', sample2_name, ' (vs ', sample1_name, ' diff sign)'),
               xlab = "EV Value")
          hist(ev2[idx][ds] - ev1[idx][ds], n = 50,
               main = paste0(chr_name, ': ', sample2_name, ' - ', sample1_name, ' (diff bins)'),
               xlab = "EV Difference")
        } else {
          plot.new(); title(paste0(chr_name, ': ', sample1_name, ' - ', sample2_name, ' (No diff signs)'))
          plot.new(); plot.new()
        }
      } else {
        plot.new(); title(paste0(chr_name, ': ', sample2_name, ' vs ', sample1_name, ' (No common data)'))
        plot.new(); plot.new()
      }
      # no need to manually start a new page; when mfrow grid fills, PDF auto-advances
    }
  }
  dev.off()
}



# MA plot: pairwise comparison across all samples
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


# chr17 染色体示例绘图
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

# chr3 染色体示例绘图
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
