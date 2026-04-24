# Run pipeline
sbatch launch.sh -with-singularity singularity/dgarrimar-ml\@sha256-94a298fc496e09fc6e1c5886661b3b9856741eefc88654bde4e76df4c7f2109c.img

# Downstream analyses
for i in result/*; do p=$(cat $i/perf.txt); echo -e "$i\t$p" | sed 's/result\///'; done | sed 's/_/\t/' | sed '1i tc\tte\trmse\tr' > summary.txt
