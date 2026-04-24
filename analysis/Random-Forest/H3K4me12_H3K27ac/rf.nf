/*
 * Train RF with variable importance (Gini, SHAP) estimation
 */

// Input params
params.chromDir = "$baseDir/data/marks"
params.marks = "$baseDir/data/marks.txt"
params.expr = "$baseDir/data/gene.expression.tsv"
params.genes = "$baseDir/data/subset/upregulated.bed"
params.tps = "$baseDir/data/timepoints.txt"
params.outDir = "$baseDir/result"
params.random = 1
params.cv = 3
params.niter = 10

// Timepoints
tps = params.tps ? file(params.tps) : System.in
Channel.from(tps.readLines()).into {tc_ch; te_ch}

// Run
 
process rf {
 
  tag {"TC: $tc, TE: $te"}

  publishDir "${params.outDir}/${tc}_${te}", mode: 'copy'

  input:
  each tc from tc_ch
  val te from te_ch
  
  output:
  set file("*.txt"), file("*.pickle") into out_ch 
  
  script:
  """
  rf.py -c ${params.chromDir} -m ${params.marks} -e ${params.expr} -g ${params.genes} -x $tc -y $te -r ${params.random} --cv ${params.cv} --n_iter ${params.niter}
  """
}
