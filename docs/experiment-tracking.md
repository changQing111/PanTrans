# PKU-IAAS experiment tracking

Use `/data/changq/miniforge3/bin/python /data/changq/Project/PanTrans/scripts/run_experiment.py`.
The runner requires Linux, Python 3.9+, Slurm, GNU time and the existing PanTrans dependencies.

## Prepare (does not submit)

```bash
python scripts/run_experiment.py prepare construct5 --purpose 'hub-direct baseline'
python scripts/run_experiment.py prepare append5plus5 --purpose 'hub-direct baseline'
python scripts/run_experiment.py prepare construct10 --purpose 'hub-direct baseline'
```

Run these commands using the absolute Python/repository paths above on PKU-IAAS.
Each prints a unique run directory under `/data/changq/Project/pantrans_test/experiment_runs`.
An actual source snapshot, Git commit/status/diff, SHA256 manifest, resolved CLI
parameters including defaults, input identities, commands and Slurm scripts are saved.
Uncommitted and untracked Python modules under `src/` are included in the snapshot.
The job uses this snapshot through PYTHONPATH, not the live editable installation.
Dependencies are recorded at execution time, not frozen into a container.

Presets reproduce the inspected scripts:

| Preset | Input directory | Important settings |
|---|---|---|
| construct5 | `/data/changq/Project/pantrans_test` | `five_variety_*`; reuse `five_variety_construct_out2/five_variety_cdna_align_gdna.bam`; rescue enabled |
| append5plus5 | same base, `five_append/` | history `five_variety_construct_out8/`; rescue enabled |
| construct10 | same base, `ten_variety/` | `ten_variety_*`; full alignment; rescue disabled |

For construct5, `--full` removes BAM reuse. Use `--rescue` or `--no-rescue` to
explicitly control rescue. Use `--history-dir /path/to/previous/run/results`
to append from a newly completed tracked construct5 run; its run ID is linked.
Prepare append only after historical files exist. Existing legacy history is
explicitly marked `legacy_external`, never misattributed to the current version.
The runner records graph sidecar files and original historical sequence inputs.

Defaults: 32 CPUs, fatPartition, 240G requested memory, two-day time limit.
Override with `--cpus`, `--partition`, `--memory`, `--time-limit`. These are requests,
not measured consumption. Original result folders are never overwritten.

## Submit

```bash
/data/changq/miniforge3/bin/python /data/changq/Project/PanTrans/scripts/run_experiment.py submit /absolute/run/directory
```

Use the runner, not bare `sbatch`: it holds the main job, records its ID, submits
an `afterany` collector, records that ID and releases the main job. If collector
submission fails, the main job stays held for repair. Do not resubmit the same run;
prepare another directory. Snapshot hashes and input metadata are checked before
analysis; input SHA256 hashes are computed before analysis timing starts.

## Records and resource accounting

`run.json` is the authoritative record; `runs.tsv` is atomically rebuilt with an
exclusive file lock. `environment.json` records Python/packages and executable
versions/paths, plus TransDecoder file hashes. `logs/time.txt` records GNU time
wall/CPU time and maximum RSS; stdout and stderr are separate. The analysis timer
excludes input hashing and environment capture. Slurm elapsed time includes them.
The post-job collector saves all Slurm job/step rows to `logs/sacct.txt` and JSON,
including CPU time, start/end, exit code, elapsed time and MaxRSS; it updates the
index even after OOM/cancellation where the application cannot finalize its record.

`gnu_time_max_rss_kib` and `slurm_step_max_rss_kib` retain their measurement scope.
Neither is asserted to be the simultaneous aggregate peak of all child processes.
A cgroup total peak is not currently collected. Missing accounting is null, not 0.
Slurm accounting availability and sampling depend on cluster configuration.
Accounting is retried for up to approximately two minutes; if still delayed or the
collector itself was cancelled, run:

```bash
python scripts/run_experiment.py collect /absolute/run/directory --retry
```

Four automated tests exercise execution records/indexing, modified-snapshot
rejection, OOM finalization and held-job/collector submission order without
submitting real Slurm jobs. Production datasets are validated during preparation;
full biological runs have not been launched by installation of this runner.
