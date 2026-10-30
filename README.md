# G4 builder library

`g4_builder` provides a Python API and terminal client for strand construction,
reference preparation, restrained GROMACS simulations and frame extraction.
See [INPUTS.md](INPUTS.md) for the complete current input contract.

## Terminal client

After installation use `g4-builder`; from the source directory you can also use
`conda run -n jupyter_env python -m g4_builder` with the same arguments.

```bash
# JSON is used only for preparation; working_dir is always a separate flag.
g4-builder prepare_nogmx -i examples/input.json -w testing/new_run
g4-builder prepare       -i examples/input.json -w testing/new_run

# Later modes read the prepared run.json; they do not accept --input.
g4-builder run_sim -w testing/new_run
g4-builder extract -w testing/new_run

# Explicit extraction paths override the paths recorded by run_sim.
g4-builder extract -w testing/new_run \
    -trj testing/new_run/trjs/trj.part0001.xtc testing/new_run/trjs/trj.part0002.xtc \
    -s smd.tpr -colvar COLVAR
```

Preparation reads simulation settings from [examples/input.json](examples/input.json).
Mandatory fields are `sequence`, `GT_resids`, `polarity_pattern`, `helicity`, and
`force_field`. Omitted optional fields use `SimulationConfig` defaults. All
input names are identical to the Python API. JSON contains only simulation
fields: no `preparation`, `run`, `extract` or `logging` blocks. Unknown keys,
duplicate keys, and invalid values raise errors. No aliases or name migrations
are performed, including when loading saved metadata.

`prepare_nogmx` writes strand, references, PLUMED protocol and run metadata.
`prepare` additionally performs topology generation, boxing, minimization and
simulation TPR preparation. It can continue a matching offline preparation.
Completed stages are reused. `--overwrite` explicitly regenerates inputs and
does not delete old MD outputs; a fresh directory is preferable for a new run.
Existing preparation settings must match the JSON, except the executable and
number of extracted frames may be updated by repeating preparation.

`run_sim` and `extract` load `run.json` for configuration, protocol duration,
full-core RMSD label and recorded output paths. Extraction accepts:

| Flag | Meaning | If omitted |
| --- | --- | --- |
| `-trj`, `--trajectories` | One XTC or an ordered list of parts | Recorded simulation trajectories; error if none are recorded |
| `-s`, `--structure` | Structure/TPR input | Recorded structure or prepared TPR |
| `-colvar`, `--colvar` | COLVAR input | Recorded COLVAR or working_dir/COLVAR |

For extraction, bare filenames such as `smd.part0001.xtc` resolve from
`--working-dir`. Relative paths containing directories, such as
`testing/new_run/smd.part0001.xtc`, resolve from the terminal's current directory.
Absolute paths are also accepted, and `~` is expanded. The JSON filename and the
working-directory flag are relative to the terminal's current directory.
Multiple trajectory parts are joined with `trjcat` before extraction.
The saved `n_frames` determines the number of extracted structures.

CLI-only policy fixes the timestep at 0.002 ps, pdb2gmx timeout at 60 s, discovery
timeout at 2 s, and run/extraction timeout unlimited. These values remain
configurable in the Python API. `dt_ps` must not appear in CLI JSON, even with
its default value. CLI logging writes full output to gromacs.log with brief
console messages; API callers may configure logging themselves.

For both API and CLI, the commands always use `-maxwarn 1` for minimization
preparation, `-maxwarn 2` for simulation preparation, `-cpi` for simulation
execution, and the `System` group for frame extraction. These are not options.

Use `g4-builder --help` or `g4-builder extract --help` to inspect available flags.
Runtime failures return exit code 1, parser errors 2 and handled interruption 130.

## Installation and environment

From the project directory, imports work immediately in `jupyter_env`. To
make them available when Python starts elsewhere, install the package:

```bash
conda run -n jupyter_env python -m pip install -e . --no-build-isolation --no-deps
```

NumPy and pandas must already be installed with that command. Normal
`pip install .` also installs the declared Python dependencies. GROMACS and
PLUMED are external programs, not Python dependencies.

GROMACS is only needed when running command stages. Activate its environment
before starting Python/Jupyter, for example `spack load /34tpcvn`. Running the
SMD stage requires a GROMACS build with PLUMED support. Offline preparation
does not look for GROMACS or start any subprocesses.

## Default filenames and derived GROMACS outputs

[g4_builder/data/file_names.json](g4_builder/data/file_names.json) defines package-level
filenames, including the two TPR inputs: `"em_tpr": "em.tpr"` and
`"tpr": "smd.tpr"`. It does not contain prefix settings or MD output names.
`RunPaths(working_dir)` reads this file to construct all run paths. There is no
`prefix` field in `SimulationConfig` or input JSON, and no prefix argument to
`RunPaths`.

`RunPaths` takes each TPR filename's stem as `mdrun -deffnm`. Thus `em.tpr`
produces prefix `em`, and `smd.tpr` produces prefix `smd`. Simulation MDP, PDB
and trajectory-part names are derived in code, without duplicate definitions
in JSON. For example, changing `tpr` to `fold.tpr` yields `fold.mdp`, `fold.pdb`,
and `fold.part*.xtc`. PDB output is explicitly requested with `-c`;
other MD outputs follow GROMACS's `-deffnm` convention.

Simulation execution always uses `-noappend` and discovers only numbered XTC
parts such as `smd.part0001.xtc`. There is no plain `smd.xtc` default path or
fallback. Extraction still accepts any explicitly supplied XTC file path.

The file retains independent filenames for references, COLVAR, application
logging, metadata and analysis outputs. Patterns use `{name}` for reference
identifiers and `{index}` for frame numbers. The packaged MDP template source
filenames are also defined there and are independent of the generated names.

Change this package data file before preparing a new run if you want different
defaults. Keep the same definitions throughout preparation, execution and
resumption; do not rename outputs of an existing run by changing the defaults
midway. Installed packages include this JSON resource.

## Python API usage

See [the pipeline example notebook](notebooks/api_pipeline.ipynb) for dictionary
and JSON configuration options, staged preparation and optional local simulation.
See [the step-by-step example notebook](notebooks/api_step_by_step.ipynb) for
the individual API calls without the pipeline, including instructions for
running the simulation on another machine and extracting frames after transfer.

```python
from g4_builder import SimulationConfig, FoldingPipeline, configure_logging

config = SimulationConfig(
    sequence="GGGTAGGGCGGCGGGGCAGGGT",
    GT_resids={
        "GT_1": [11, 1, 2],
        "GT_2": [8, 7, 6],
        "GT_3": [19, 20, 21],
        "GT_4": [15, 14, 13],
    },
    polarity_pattern="saa",
    helicity="RH",
    working_dir="testing/new_run",
    force_field="amber99bsc1",  # Use the short -ff name on your installation.
    gromacs_executable="gmx",
    dielectric_const=80,
    saving_stride=5000,
    after_folding_steps=2_500_000,
    n_frames=5,
)

pipeline = FoldingPipeline(config)
configure_logging(pipeline.paths.log)  # Optional; full output in gromacs.log.

# Offline: creates the directory, strand, references, PLUMED file, run.json.
result = pipeline.prepare()

# Inspect these before running external commands.
print(result.protocol_steps, result.simulation_steps)
print(result.full_ref_label)
print(pipeline.references.ref_info)
print(pipeline.protocol.text)

# Explicit external-command stages:
pipeline.minimize()
pipeline.prepare_simulation()

# Run only with a PLUMED-enabled GROMACS installation:
pipeline.run()
result = pipeline.extract()

print(result.best_frames)  # {time_in_ps: rmsd}, ordered from lowest RMSD.
print(result.frame_paths)  # {1: Path(.../frame_1.pdb), ...}
```

Choose a new directory for a new simulation. `prepare()` refuses to overwrite
existing run inputs by default. `prepare(overwrite=True)` explicitly replaces
inputs, but does not clean out old trajectory files; use it only when that is
intended. Command stages also use GROMACS's usual output/backup behavior.

`configure_logging()` installs brief console messages and a detailed append-only
log. Importing or constructing library objects does not install handlers or
open log files. You may instead configure the `g4_builder` Python logger in
your application. GROMACS stdout/stderr are streamed to logging and are never
returned by command methods. Command errors raise exceptions and log their
diagnostics. Set the console level to `logging.DEBUG` to display full output.
Reference and extracted-frame save methods still print their names and paths.

## What each stage does

| Stage | Behavior |
| --- | --- |
| `prepare()` | Validates geometry/template correspondence; builds strand and full/partial references; generates PLUMED text and metadata. |
| `minimize()` | Runs `pdb2gmx`, checks reference atom serials against processed atom identifiers, prints a notice and remaps/saves references only if they changed, updates PLUMED, adds the box, prepares and runs minimization. |
| `prepare_simulation()` | Uses protocol steps plus post-folding steps to render the MDP, then creates the package-defined simulation TPR. |
| `run()` | Runs PLUMED-enabled `mdrun` with `-cpi` and `-noappend`; records actual XTC part filenames. |
| `extract()` | Selects post-protocol COLVAR minima using the full-core RMSD label and writes ranked PDB frames. |

Stages return the same inspectable `RunResult`. GROMACS wrapper methods return
`None`. All paths are available via `pipeline.paths`, and reference paths via
`result.reference_paths`. Input residue lists are 1-based and ordered; tract
dictionary order defines template correspondence. Four tracts of equal length
are required, with at least two residues per tract. Each selected residue must
be guanine, and tract length must match the polarity-pattern length.

### Files required in the simulation directory

**All partial reference PDB files, the PLUMED protocol file and the GROMACS
simulation input (`.tpr`) must be in the same directory when running the
simulation. Run GROMACS with that directory as its working directory.**
`build_protocol()` writes bare reference filenames such as `ref_1.pdb`; it
does not accept a reference-directory option. The pipeline places these files
in `pipeline.paths.working_dir` and runs GROMACS there automatically.

By default, template combinations come from the stored package data. A missing
polarity/helicity template raises `FileNotFoundError`. The lower-level
`build_references(..., build_template=True)` call can instead generate a template
from the packaged flat layers, as described below.

## Resume or analyze a run from another session

After every successful pipeline stage, `run.json` records the validated config,
protocol duration and label, reference metadata, completed stages and known
outputs. It does not auto-skip stages or restart failed commands.

```python
pipeline = FoldingPipeline.from_directory("testing/new_run")
configure_logging(pipeline.paths.log)

# For example, if minimization was the last completed stage:
pipeline.prepare_simulation()
# pipeline.run()
```

If you launch MD externally, supply its trajectory path explicitly:

```python
result = pipeline.extract(trajectories="smd.part0001.xtc")

# Multiple restart parts are joined with trjcat before frame extraction:
result = pipeline.extract(
    trajectories=["smd.part0001.xtc", "smd.part0002.xtc"],
    structure="smd.tpr",
    colvar="COLVAR",
    n=3,
)
```

Bare filenames resolve from the working directory; relative paths containing
directories resolve from Python's current directory. Absolute paths are accepted.
These rules apply to trajectories, the structure input and COLVAR. For an
internally executed run, `extract()` uses the output paths recorded by `run()`.
The prepared TPR is used as the structure input so extraction uses authoritative
topology atom ordering. Metadata can be restored after moving the whole run
directory; paths inside the old directory are relocated. Externally referenced
absolute paths are preserved and may need explicit replacements.

COLVAR and the selected trajectories must cover the same run. Selection is
stable by RMSD, uses one entry per rounded integer time, and requires at least
`n` distinct post-protocol times. Times are picoseconds. The protocol boundary
is `protocol_steps * dt_ps`, and that same `dt_ps` is written into the MDP
(default 0.002 ps). `trjconv -dump` selects the nearest saved trajectory frame;
matching COLVAR/trajectory output strides is recommended. If COLVAR extends
beyond available trajectory files, supply the matching files/data rather than
expecting COLVAR alone to establish the XTC duration.

For lower-level extraction, configure only the shared runner on the object and
pass each operation its own inputs:

```python
import g4_builder as g4

extractor = g4.TrajectoryFrameExtractor(runner)
extractor.find_best_frames(
    rmsd_label=protocol.full_ref_label,
    protocol_steps=protocol.final_steps,
    colvar="COLVAR", n=5, dt_ps=0.002,
)
extractor.extract_frames(
    trajectories=["smd.part0001.xtc", "smd.part0002.xtc"],
    structure="smd.tpr",
)
```

Selection reads COLVAR only; it needs no trajectory, structure or GROMACS
execution. The selected `{time_ps: RMSD}` mapping is returned and stored in
`extractor.best_frames`. Extraction uses that selection, prints output paths and
stores them in `extractor.saved_paths`. Starting a new selection clears both
results, including if selection fails. The same object can be reused for another
selection or a different set of trajectory files.

## Lower-level API

The individual components remain usable without the pipeline:

```python
from g4_builder import (
    DNAStrandBuilder, build_references, build_protocol, ScheduleConfig,
    GromacsRunner, GromacsPreparation, read_pdb,
)

# The runner owns the shared paths and controls subprocess execution.
# Constructing it does not launch GROMACS or create the working directory.
runner = GromacsRunner(config.working_dir, config.gromacs_executable)
paths = runner.paths
gromacs = GromacsPreparation(runner)

strand_builder = DNAStrandBuilder()
strand = strand_builder.build(config.sequence)

references = build_references(
    config.polarity_pattern, config.helicity, strand, config.GT_resids,
    max_bulge_length=config.max_bulge_length,
)
# references is a fully built ReferenceSet: full reference, partials and metadata.
paths.working_dir.mkdir(parents=True, exist_ok=True)
strand_builder.save_pdb(paths.strand)
references.save_full_ref(paths.full_ref)
protocol = build_protocol(references, strand.serial, schedule=config.schedule,
                          print_stride=config.saving_stride, print_file=paths.colvar.name)
protocol.write(paths.plumed)
references.save_partial_refs(paths.working_dir)  # Prints; paths in references.saved_paths.

# After pdb2gmx, check atom identities and update changed serial numbers:
# processed = read_pdb(paths.conf)
# if references.needs_serial_remap(processed):
#     references.remap_serials(processed)
#     references.save_full_ref(paths.full_ref)
#     references.save_partial_refs(paths.working_dir)

# Force-field discovery is explicit:
# print(gromacs.discover_force_fields())
```

`GromacsPreparation(runner)` and `TrajectoryFrameExtractor(runner)` use
`runner.paths` directly. Supply the timestep when calling
`gromacs.prep_sim_file(..., dt_ps=0.002)` and
`extractor.find_best_frames(..., dt_ps=0.002)`.
`FoldingPipeline` also shares that object, so the working directory and filename
defaults are defined once per runner. Use `runner.paths` to inspect or save inputs;
`runner.working_dir` is a property exposing its working directory. A runner passed
to the pipeline must use the configuration's working directory.
`RunPaths(directory)` remains available for inspecting filenames independently.
Custom runners should subclass `GromacsRunner` and override `run()` / `probe()`
as needed, retaining the shared paths initialized by its constructor.

`build_references(polarity_pattern, helicity, strand, GT_resids,
max_bulge_length=3, template_dir=None, build_template=False)` returns a fully built
`ReferenceSet`. It reads the template, validates atom correspondence and creates the full and
partial reference tables without writing files or changing the inputs.
There is no separate builder object or intermediate build state.

To generate the template from planar layers instead of loading a stored stack:

```python
references = build_references(
    config.polarity_pattern, config.helicity, strand, config.GT_resids,
    build_template=True,
)
```

This option exists only on `build_references`; configs, JSON, CLI and pipeline
defaults continue loading stored templates. With `template_dir` supplied, loading
reads `<template_dir>/<helicity>_<n>tet/<pattern>.pdb`, while building reads
`<template_dir>/layers/a.pdb` and `s.pdb`.

`g4_builder.build_reference_template(pattern, helicity, layer_dir=None)` is also
available directly and returns the generated template DataFrame. It uses the
packaged `data/reference_templates/layers/` files by default. Rise and signed
twist values are defined in
[g4_builder/data/stack_properties.json](g4_builder/data/stack_properties.json),
which is loaded on each build. They accumulate for adjacent layer pairs.
Residues are grouped across layers before advancing to the next original residue. The top-level
`ref_from_scratch.py` imports this implementation rather than maintaining a copy.
Editing top-level layer files does not change the packaged copies; use an
explicit directory override when working with experimental layers.

The resulting object provides `full_ref`, `partial_refs`, `ref_info`, `entries`
and `full_ref_label`. Use `save_full_ref(path)` and `save_partial_refs(directory)`
to write references without TER/END records; partial saving prints paths and
records them in `references.saved_paths`. Use `needs_serial_remap(structure)`
and `remap_serials(structure)` after topology preparation. The pipeline uses
these same methods, and `RunResult.reference_paths` holds its partial file paths.

Both the pipeline configuration and lower-level wrapper accept short force-field
names or positive integer menu numbers, for example `force_field=4`. With a
number, the menu is discovered automatically when `pdb2gmx` runs, then the
selection is validated and sent through stdin. Construction and offline
preparation still start no processes. Menu numbers depend on the installation;
short names are more reproducible across installations. TIP3P remains fixed
as in the existing implementation.

`gromacs.pdb2gmx(force_field, strand="inputs/custom.pdb")` accepts an optional
strand PDB path. Relative paths use the runner's working directory. If omitted,
it uses `runner.paths.strand`. Automatic force-field discovery uses the same
supplied strand, so the default strand file is not required for this override.

`ReferenceSet.entries` holds each reference's name, category (`GT`, `core`, or
`tetrads`), ordered residues, RMSD label and atom DataFrame. `is_full_core`
identifies the final cumulative core explicitly; equal DataFrame sizes cannot
make the label ambiguous. Block splitting and schedule timing are defined in the package configuration. `ScheduleConfig` separates those timings from formatting:

```python
from dataclasses import replace
modified = replace(config, schedule=ScheduleConfig(core_spacing=3_000_000))
```

`config.to_dict()` exports JSON-compatible input settings without `working_dir`.
Supply the directory separately when loading a mapping:

```python
config = SimulationConfig.from_dict(config.to_dict(), working_dir="testing/new_run")
```

Load an input file directly through the Python API:

```python
config = SimulationConfig.from_json("examples/input.json", working_dir="testing/new_run")
pipeline = FoldingPipeline(config)
```

This uses the same parameter names as the Python constructor, with defaults for
omitted optional settings. API-only settings such as `dt_ps`
are accepted here; the CLI restrictions do not apply. Both `from_dict()` and
`from_json()` require a separate `working_dir` argument and reject that key in
the input dictionary or JSON. Relative working directories resolve from Python's
current working directory. The normal `SimulationConfig(...)` constructor still
takes `working_dir` as an ordinary argument.
Only simulation fields are accepted. Execution/file settings are passed to
the corresponding methods, such as `pipeline.extract(trajectories=...)`,
`pipeline.run(timeout=3600)` or `configure_logging(...)`. Unsupported blocks
and unknown names raise errors; nothing is silently discarded or renamed.

## Code layout and tests

`g4_builder/config.py`, `paths.py`, `pdb_io.py`, `strand.py`, `references.py`,
`protocol.py`, `runner.py`, `gromacs.py`, `extraction.py`, `pipeline.py` and
`logging.py` separate the responsibilities. `g4_builder/data/` contains copies
of the 12 nucleotide templates, supported reference templates, flat a/s layers
and two MDPs.
They are included in installed wheels. Editing the old top-level template
files does not modify these package copies.

Run the independent package tests using the project environment:

```bash
conda run -n jupyter_env python -m unittest discover -s tests_new -v
```

They exercise offline preparation with the real templates, config validation,
reference metadata and remapping, schedules, full pipeline orchestration with
a fake GROMACS runner, metadata restoration, extraction and real subprocess
logging/failure handling. They do not run molecular dynamics or require PLUMED.
