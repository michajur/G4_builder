# Current input contract

There are three distinct inputs:

1. **Simulation configuration:** a `SimulationConfig`, or JSON containing its fields.
2. **Execution arguments:** CLI flags or Python method arguments for a particular action.
3. **Prepared run metadata:** `run.json`, written by the pipeline and read when resuming.

The package code defines these conventions. `dev_notebook.ipynb` is an independent
testing workspace and does not define the input format or defaults.

## CLI

| Mode | Configuration source | Additional inputs |
| --- | --- | --- |
| `prepare_nogmx` | Required `-i/--input` JSON | Required `-w/--working-dir`; optional `--overwrite` |
| `prepare` | Required `-i/--input` JSON | Required `-w/--working-dir`; optional `--overwrite` |
| `run_sim` | `run.json` inside the working directory | Required `-w/--working-dir` |
| `extract` | `run.json` inside the working directory | Required `-w/--working-dir`; optional extraction path flags below |

`run_sim` and `extract` do not accept `--input`: they use the preparation's saved
settings and reference labels rather than rebuilding the model from a new JSON.

Default filenames come from
[g4_builder/data/file_names.json](g4_builder/data/file_names.json). For MD stages,
only TPR names are defined: `em_tpr` defaults to `em.tpr`, and `tpr` to `smd.tpr`.
The `-deffnm` prefixes and MD output filenames are derived from these TPR names
in code. Neither prefixes nor derived MD output names are configurable entries.
`RunPaths(working_dir)` also reads the metadata filename (default `run.json`).
Keep those definitions unchanged throughout the stages of an existing run.
With the supplied defaults, simulation files are lowercase `smd.mdp`, `smd.tpr`,
`smd.pdb`, and `smd.xtc`.

```bash
g4-builder prepare_nogmx -i examples/input.json -w testing/new_run
g4-builder prepare -i examples/input.json -w testing/new_run
g4-builder run_sim -w testing/new_run

# Use the paths recorded by a completed run_sim:
g4-builder extract -w testing/new_run

# Or specify paths to an externally run simulation:
g4-builder extract -w testing/new_run \
    -trj testing/new_run/trjs/trj.part0001.xtc \
    -s smd.tpr -colvar COLVAR

# Supply multiple trajectory parts in their intended order:
g4-builder extract -w testing/new_run \
    -trj testing/new_run/trjs/trj.part0001.xtc testing/new_run/trjs/trj.part0002.xtc \
    -s smd.tpr -colvar COLVAR
```

The long flags are `--trajectories`, `--structure` and `--colvar`.
`-trj` takes one or more paths; `-s` and `-colvar` each take one.
Explicit flags override recorded paths. Without `-trj`, recorded simulation
outputs must exist. Without `-s`, use the recorded structure or prepared TPR;
without `-colvar`, use the recorded COLVAR or `working_dir/COLVAR`.

### Path resolution

- `--input` and `--working-dir` are relative to the terminal's current directory.
- Bare extraction filenames such as `smd.xtc` resolve from `--working-dir`.
- Relative extraction paths containing directories, such as
  `testing/new_run/smd.xtc`, resolve from the terminal's current directory.
- Absolute paths work in either case; `~` is expanded.
- `gromacs_executable` is a command on PATH or an executable path. A relative
  executable path containing `/` is resolved from the caller's current directory.
- Output frames are written into the working directory.

For example, with `-w /tmp/run -trj part1.xtc`, the trajectory is
`/tmp/run/part1.xtc`. With `-trj trjs/part1.xtc`, it is `trjs/part1.xtc` relative
to the terminal's current directory. To point inside a subdirectory of the run,
provide its absolute path or its complete path relative to the current directory.
These rules apply equally to trajectory, structure and COLVAR inputs and do not
depend on whether candidate files happen to exist.

### What belongs in the JSON

Only simulation configuration belongs in the JSON. There are no `preparation`,
`run`, `extract`, or `logging` blocks. In particular, trajectory/structure/COLVAR
paths belong in extraction flags, not JSON.

The full example is [examples/input.json](examples/input.json). Its mandatory
fields are:

```json
{
  "sequence": "GGGTAGGGCGGCGGGGCAGGGT",
  "GT_resids": {
    "GT_1": [11, 1, 2],
    "GT_2": [8, 7, 6],
    "GT_3": [19, 20, 21],
    "GT_4": [15, 14, 13]
  },
  "polarity_pattern": "ssa",
  "helicity": "RH",
  "force_field": 4
}
```

`GT_resids` uses ordered, 1-based residue indices. There must be four equal-length
tracts, each with at least two guanines and matching the polarity-pattern length.
Residues must be unique across tracts. Tract dictionary order determines template
correspondence. `force_field` accepts a short name or a positive integer menu
number. Optional fields use defaults:

| JSON field | Default |
| --- | --- |
| `gromacs_executable` | `"gmx"` |
| `dielectric_const` | `15` |
| `saving_stride` | `5000` steps |
| `after_folding_steps` | `2500000` steps |
| `max_bulge_length` | `3` |
| `n_frames` | `5` |
| `schedule` | `ScheduleConfig` defaults, as shown in the full example |

Partial `schedule` objects use defaults for omitted settings. The only schedule
inputs are `gt_initial_steps` (default 400000), `core_start` (1000000),
`core_steps` (3000000) and `tetrads_steps` (3000000), all positive integers.
`core_start` must exceed `gt_initial_steps`, `core_steps` must be divisible by
three, and `tetrads_steps` must exceed `core_steps / 3`.

Stage arrays, targets, force constants and `core_spacing` are derived rather
than configurable. All RMSD targets are zero. The GT force constant ramps from
0 to 12000 by `gt_initial_steps`. Each core block ramps to 20000 at one third
of `core_steps`, holds until two thirds, then falls to zero at `core_steps`.
Core blocks start two thirds of `core_steps` apart. The tetrad block ramps to
20000 at `core_steps / 3`, then ends at `tetrads_steps` with a force constant
of 5000. GT holding continues until the tetrad block starts and GT release
ends at its second stage. The default model timestep is 0.002 ps; protocol
steps are converted to COLVAR time using that timestep for extraction.

The CLI fixes `dt_ps` at 0.002 ps and rejects that field in CLI JSON. It also
uses no run/extraction timeout, a 60-second pdb2gmx timeout, a 2-second discovery
timeout, and INFO console messages plus `gromacs.log`.

Both API and CLI always use `-maxwarn 1` for minimization preparation,
`-maxwarn 2` for simulation preparation, `-cpi` when running the simulation,
and the `System` group for frame extraction. These command settings are not
configuration fields or method arguments.

Unknown fields, duplicate keys, execution blocks and unsupported names are
errors. No keys are silently dropped, renamed or migrated.

### What comes from run.json

Preparation writes the validated config, full-core RMSD label, protocol step
count, total simulation steps, reference metadata and completed stages. Running
the simulation adds output trajectory/structure/COLVAR paths. Extraction records
the selected times/RMSDs and output frame paths.

Run metadata also records the working directory so paths can be relocated when
moving a run. This internal manifest is separate from the input JSON format;
the pipeline separates the directory from the settings when restoring a config.

The client restores these facts from `run.json`. It does not read or persist a
second set of execution-option blocks. Missing or unsupported config names in
metadata raise an error; there is no compatibility conversion.

`prepare` can continue matching offline preparation. Repeating preparation may
update only the executable and `n_frames` without changing the model; other
differences require a new run or explicit `--overwrite`. This flag regenerates
inputs but does not delete previous MD output files.

## Python API

Use exactly the same simulation field names as CLI JSON. Three equivalent
construction methods are available:

```python
import g4_builder as g4

# 1. JSON simulation fields plus an explicit working directory:
config = g4.SimulationConfig.from_json(
    "examples/input.json", working_dir="testing/api_run"
)

# 2. Simulation settings plus a separate working directory:
values = config.to_dict()
config = g4.SimulationConfig.from_dict(values, working_dir="testing/api_run")

# 3. Normal constructor keyword arguments:
config = g4.SimulationConfig(
    sequence="GGGTAGGGCGGCGGGGCAGGGT",
    GT_resids={"GT_1": [11, 1, 2], "GT_2": [8, 7, 6],
               "GT_3": [19, 20, 21], "GT_4": [15, 14, 13]},
    polarity_pattern="ssa", helicity="RH", force_field=4,
    working_dir="testing/api_run",
)
```

The API also allows `dt_ps`. `from_dict(values, working_dir=...)` and
`from_json(path, working_dir=...)` both require the directory as a separate
argument. Neither the input dictionary nor the JSON may contain `working_dir`,
even if it matches the separate argument. `config.to_dict()` omits that field
so its output can be reused as input settings. The normal constructor still
accepts `working_dir` directly. A relative working directory resolves from
Python's current directory. API JSON accepts only simulation settings, with no
special execution blocks or automatic name conversion.

Pass execution/file options directly to methods:

```python
pipeline = g4.FoldingPipeline(config)
g4.configure_logging(pipeline.paths.log)
pipeline.prepare()
pipeline.minimize(pdb2gmx_timeout=60, discovery_timeout=2)
pipeline.prepare_simulation()
pipeline.run(timeout=None)
pipeline.extract(
    trajectories=["testing/api_run/trjs/part1.xtc", "testing/api_run/trjs/part2.xtc"],
    structure="smd.tpr", colvar="COLVAR", n=3,
    timeout=None,
)
```

Extraction paths in the API use the same rules as CLI flags: bare filenames use
the pipeline working directory, directory-containing relative paths use Python's
current directory, and absolute paths retain their location. Omitted paths use recorded outputs where
available. API callers can change execution arguments and logging explicitly;
the CLI uses their fixed defaults. `FoldingPipeline.from_directory(directory)`
restores a run for further API actions using the same `run.json` as the CLI.

The lower-level extractor takes only the runner at construction. Selection
parameters belong to `find_best_frames`; file inputs for extraction belong to
`extract_frames`:

```python
extractor = g4.TrajectoryFrameExtractor(runner)
extractor.find_best_frames(
    rmsd_label=protocol.full_ref_label, protocol_steps=protocol.final_steps,
    colvar="COLVAR", n=5, dt_ps=0.002,
)
extractor.extract_frames(trajectories="smd.xtc", structure="smd.tpr")
```

`colvar=None` uses `runner.paths.colvar`. `find_best_frames` returns and stores
the selected times/RMSDs; every new selection clears the previous results.
`extract_frames` requires a successful selection and accepts one trajectory path
or an ordered iterable of paths, plus optional `timeout`. Saved paths are printed
and stored in `extractor.saved_paths`.

All reference PDB files, the PLUMED protocol file and the simulation `.tpr`
must share the simulation working directory. Protocol generation uses reference
filenames without a directory prefix. The pipeline places these files together;
when using lower-level components, save them in the same directory yourself.

Reference preparation is a function call returning a `ReferenceSet`:

```python
references = g4.build_references(
    polarity_pattern="ssa", helicity="RH", strand=strand,
    GT_resids=config.GT_resids, max_bulge_length=3,
)
references.save_full_ref(paths.full_ref)
references.save_partial_refs(paths.working_dir)

if references.needs_serial_remap(processed_structure):
    references.remap_serials(processed_structure)
    references.save_full_ref(paths.full_ref)
    references.save_partial_refs(paths.working_dir)
```

The function accepts the strand DataFrame and ordered G-tract residue mapping
directly; it does not read simulation input JSON. Optional `template_dir`
overrides the packaged reference-template root. The returned object owns the
reference tables, metadata, saving and atom-serial remapping.

`build_references(..., build_template=False)` loads a stored template by default.
Pass `build_template=True` to generate one using `build_reference_template()`.
When overriding `template_dir`, building reads the flat PDBs from its `layers/`
subdirectory; loading uses its `<helicity>_<n>tet/<pattern>.pdb` file.
This switch is a function argument only: it is not a `SimulationConfig` field,
JSON setting, CLI option or pipeline option. For a template DataFrame alone,
use `g4.build_reference_template(pattern, helicity, layer_dir=None)`.

For lower-level use, create `runner = g4.GromacsRunner(working_dir)` and get
paths through `runner.paths`. Construct preparation with
`g4.GromacsPreparation(runner)`; it and
`g4.TrajectoryFrameExtractor(runner)` share that paths object. There is no
separate paths argument to preparation. Pass `dt_ps` to `prep_sim_file()` and
`find_best_frames()`. It must be the same when preparing
the simulation MDP and selecting post-protocol frames; the pipeline passes the
configuration's value to both automatically.
