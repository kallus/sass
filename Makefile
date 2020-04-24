fig/report.pdf: figures src/report.tex
	(cd src && texcompile) # change texcompile to pdflatex report.tex
	#(cd src && pdflatex report.tex)
	mv src/report.pdf fig/

figures: fig/clusters.png fig/hierarchic.png fig/hist.pdf fig/improvement.pdf fig/roc.pdf fig/sim.pdf fig/pathway.pdf fig/auc.pdf fig/sim_cluster.pdf fig/sim_scale-free.pdf fig/overlapping.pdf fig/adjdiff.pdf fig/pcor.pdf fig/rocs.pdf

fig/clusters.png: src/fig_clusters.R cache/contender_output/contenders_run
	Rscript src/fig_clusters.R

fig/hierarchic.png: src/fig_hierarchic.R cache/contender_output/contenders_run
	Rscript src/fig_hierarchic.R

fig/hist.pdf: src/fig_hist.R cache/contender_output/contenders_run
	Rscript src/fig_hist.R

fig/improvement.pdf: src/fig_improvement.R cache/contender_output/contenders_run
	Rscript src/fig_improvement.R

fig/roc.pdf: src/fig_roc.R cache/contender_output/contenders_run
	Rscript src/fig_roc.R

fig/sim.pdf: src/fig_sim.R cache/contender_output/contenders_run
	Rscript src/fig_sim.R

fig/sim_cluster.pdf: src/fig_sim_cluster.R cache/contender_output/contenders_run
	Rscript src/fig_sim_cluster.R

fig/sim_scale-free.pdf: src/fig_sim_scale-free.R cache/contender_output/contenders_run
	Rscript src/fig_sim_scale-free.R

fig/overlapping.pdf: src/fig_overlapping.R cache/contender_output/contenders_run
	Rscript src/fig_overlapping.R

fig/pathway.pdf: src/fig_pathway.R cache/contender_output/contenders_run
	Rscript src/fig_pathway.R

fig/auc.pdf: src/fig_auc.R cache/contender_output/contenders_run
	Rscript src/fig_auc.R

fig/adjdiff.pdf: src/fig_adjdiff.R cache/contender_output/contenders_run
	Rscript src/fig_adjdiff.R

fig/pcor.pdf: src/fig_pcor.R cache/glasso_sassa_cmp.tsv
	Rscript src/fig_pcor.R

fig/rocs.pdf: src/fig_rocs.R cache/contender_output/contenders_run
	Rscript src/fig_rocs.R

cache/glasso_sassa_cmp.tsv: src/run_adaptive.R
	Rscript src/run_adaptive.R

cache/contender_output/contenders_run: src/run_contenders.R
	Rscript src/prepare_dirs.R
	Rscript src/run_contenders.R

# Run GRAB manually:
#   cd src/grab
#   source .venv/bin/activate
#   python estimate.py
#   deactivate
# do this after running run_contenders.R and then run run_contenders.R again
