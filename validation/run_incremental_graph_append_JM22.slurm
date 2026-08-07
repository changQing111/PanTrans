#!/bin/bash -l
#SBATCH --job-name=pantrans-append-JM22
#SBATCH --partition=m256Partition
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=3-00:00:00
#SBATCH --output=/data/changq/PanTrans/.worktrees/append-construct-flow/validation/incremental_graph_append_JM22/slurm-%j.out
#SBATCH --error=/data/changq/PanTrans/.worktrees/append-construct-flow/validation/incremental_graph_append_JM22/slurm-%j.err

set -euo pipefail

worktree=/data/changq/PanTrans/.worktrees/append-construct-flow
output_dir=${worktree}/validation/incremental_graph_append_JM22
history_dir=${output_dir}/history_graph
history_construct_dir=${history_dir}/construct
append_dir=${output_dir}/append
comparison_dir=${output_dir}/comparison
python_bin=/data/changq/miniforge3/bin/python

mkdir -p "${history_construct_dir}" "${append_dir}" "${comparison_dir}"
export PYTHONPATH=${worktree}/src
export PATH=/home/changq/my_data/bio_soft/cactus-bin-v2.9.8/bin:/data/changq/miniforge3/bin:${PATH}
export LC_ALL=C

history_prefix=five_varieties
history_manifest=${history_dir}/five_varieties.graph.json
history_gtf=${history_construct_dir}/${history_prefix}_unrenamed.gtf
history_cdna=${history_construct_dir}/${history_prefix}_unrenamed_cdna.refer.fasta

if [[ ! -s "${history_gtf}" || ! -s "${history_cdna}" ]]; then
    "${python_bin}" -m pantrans.main construct \
        --name /data/changq/PanTrans/test/tmp_variety.txt \
        --cdna /data/changq/PanTrans/test/tmp_cdna.fasta \
        --gdna /data/changq/PanTrans/test/tmp_gdna.fasta \
        --bed /data/changq/PanTrans/test/tmp.bed \
        --bam /data/changq/PanTrans/test/pantrans_construct_laste/test_pantrans_cdna_align_gdna.filtered.bam \
        --main-chroms /data/changq/PanTrans/test/tmp_chrom.txt \
        --reference CS \
        --prefix "${history_prefix}" \
        --threads "${SLURM_CPUS_PER_TASK}" \
        --output "${history_construct_dir}"
fi

if [[ ! -s "${history_manifest}" ]]; then
    "${python_bin}" "${worktree}/scripts/export_graph_package.py" \
        --filtered-bam /data/changq/PanTrans/test/pantrans_construct_laste/test_pantrans_cdna_align_gdna.filtered.bam \
        --cdna /data/changq/PanTrans/test/tmp_cdna.fasta \
        --gdna /data/changq/PanTrans/test/tmp_gdna.fasta \
        --bed /data/changq/PanTrans/test/tmp.bed \
        --name /data/changq/PanTrans/test/tmp_variety.txt \
        --reference CS \
        --main-chroms /data/changq/PanTrans/test/tmp_chrom.txt \
        --coverage-min 0.80 \
        --identity-min 0.90 \
        --soft-clip-max 0.10 \
        --filter-logic-id pantrans-filter-v1 \
        --output "${history_manifest}"
fi

/usr/bin/time -v -o "${output_dir}/append.time.txt" \
    "${python_bin}" -m pantrans.main append \
        --name JM22 \
        --cdna "${history_cdna}" \
        --history-gtf "${history_gtf}" \
        --query-cdna /data/changq/PanTrans/test/JM22_cdna.fasta \
        --gdna /data/changq/PanTrans/test/JM22.gdna.fasta \
        --bed /data/changq/PanTrans/test/test_pantrans_pre.refer_append_JM22.bed \
        --history-graph "${history_manifest}" \
        --prefix incremental_append_JM22 \
        --threads "${SLURM_CPUS_PER_TASK}" \
        --output "${append_dir}"

"${python_bin}" "${worktree}/scripts/compare_append_to_construct.py" \
    --incremental-edges "${append_dir}/incremental_append_JM22.graph.edges.tsv" \
    --construct-bam /data/changq/PanTrans/test2/pantrans_construct_out/test_pantrans_cdna_align_gdna.filtered.bam \
    --incremental-pre "${append_dir}/incremental_append_JM22_pre.tmp.cluster" \
    --construct-pre /data/changq/PanTrans/test2/pantrans_construct_out/test_pantrans_pre.tmp.cluster \
    --incremental-last "${append_dir}/incremental_append_JM22_last.tmp.cluster" \
    --construct-last /data/changq/PanTrans/test2/pantrans_construct_out/test_pantrans_last.tmp.cluster \
    --output-dir "${comparison_dir}"
