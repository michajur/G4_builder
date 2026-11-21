# G4 Builder

G4 Builder prepares DNA strands and G-quadruplex reference structures, generates
PLUMED folding restraints, prepares and runs GROMACS simulations, and extracts
trajectory frames with the lowest post-protocol RMSD. It provides a command-line
interface for the complete workflow and a Python API for running individual
steps or composing a pipeline. Molecular structures are represented as pandas
DataFrames and can be saved as PDB files.

## Installation

### Requirements

| Dependency | Requirement |
| --- | --- |
| Python | 3.9 or newer |
| NumPy | 1.23 or newer |
| pandas | 1.5 or newer |
| GROMACS | Required for topology generation, minimization, simulation-input preparation and frame extraction |
| PLUMED-enabled GROMACS | Required for running folding simulations |

Strand construction, reference building and PLUMED protocol generation can run
without GROMACS. GROMACS and PLUMED are external software; installing this
Python package does not install them.

### Install from a checkout

In your chosen Python environment, run the following from the directory
containing `pyproject.toml`:

```bash
python -m pip install .
```

For development, use an editable installation so imports use the source tree:

```bash
python -m pip install -e .
```

If the environment already contains the declared dependencies and setuptools
build tools, dependency installation and build isolation can be skipped:

```bash
python -m pip install -e . --no-build-isolation --no-deps
```

The installable distribution and terminal command are named `g4-builder`;
the Python import is `g4_builder`. Check the installation with:

```bash
g4-builder --help
python -m g4_builder --help
```

### Make GROMACS available

Activate the environment that provides GROMACS before starting Python, Jupyter
or the terminal client. On a cluster this may involve an environment module or
Spack. Set the executable using `gromacs_executable` in the simulation
settings, or `executable` when constructing a `GromacsRunner`.

The executable may be a command on `PATH`, such as `gmx`, or an absolute path.
A relative executable path containing `/` is resolved from the caller's current
directory. The package does not activate external software environments itself.
For notebooks, select a kernel in the environment where the package is installed.

## Command-line interface

### Basic workflow

Edit [examples/input.json](examples/input.json) for your sequence, G-tracts and
available force field. Choose a new run directory, then run:

```bash
# Prepare strand, references and PLUMED protocol without calling GROMACS.
g4-builder prepare_nogmx -i examples/input.json -w runs/example

# Continue preparation through topology generation, minimization and the TPR.
g4-builder prepare -i examples/input.json -w runs/example

# Run with a PLUMED-enabled GROMACS installation.
g4-builder run_sim -w runs/example

# Extract the lowest-RMSD frames from the supplied trajectory.
g4-builder extract -w runs/example -trj smd.xtc
```

`prepare` can also be used directly in a new directory; it performs all the
preparation steps, including those in `prepare_nogmx`.

| Mode | Action | Configuration source |
| --- | --- | --- |
| `prepare_nogmx` | Build and save the strand, references, PLUMED protocol and run metadata | Required input JSON |
| `prepare` | Perform the above, generate topology, add a box, minimize and prepare the simulation TPR | Required input JSON |
| `run_sim` | Run the prepared simulation | Saved `run.json` |
| `extract` | Select and extract low-RMSD post-protocol frames | Saved `run.json`, with optional file-path overrides |

All modes require `-w/--working-dir`. Preparation modes additionally require
`-i/--input` and accept `--overwrite`. `run_sim` and `extract` do not accept
an input JSON: they use the configuration and protocol information recorded
during preparation and stored inside the working directory.

### Input JSON

All simulation settings for CLI preparation are provided through the input
JSON file. Supply the working directory separately with `-w/--working-dir`;
do not include `working_dir` in the JSON. Extraction file paths are supplied
through their own CLI flags. Missing optional entries use their defaults.
Unknown fields and duplicate JSON keys raise errors.

**Required entries**

#### `sequence`

DNA sequence written in the 5′ → 3′ direction, using `A`, `T`, `G` and `C`.
At least two bases are required.

#### `GT_resids`

An object containing four lists of 1-based residue indices identifying
guanosines in `sequence`. The lists define the four G-columns in the G4
structure, in counterclockwise order. Their order in the JSON object determines
their correspondence to the reference template; the example names them
`GT_1` through `GT_4`.

Within each list, residues are ordered from the bottom tetrad to the top.
For example, `"GT_1": [11, 1, 2]` places G11 in the bottom tetrad, G1 in the
middle tetrad and G2 in the top tetrad. All four lists must have the same
length, which defines the number of tetrads in the resulting G4 structure.

Each list must contain at least two residues. All indices must refer to
guanines in the sequence and must be unique across all four lists. Preserve
the intended order of both the lists and their indices rather than sorting them.

> **Figure placeholder:** G-column order and the mapping of residue indices
> to tetrads.

#### `polarity_pattern`

A string of `a` and `s` characters defining the orientation of guanine bases
in each tetrad, from bottom to top. There must be one character per tetrad,
so the string's length must match the length of each list in `GT_resids`.

> **Figure placeholder:** Guanine orientations corresponding to `a` and `s`
> and their order in a tetrad stack.

#### `helicity`

Either `"RH"` or `"LH"`, selecting right-handed or left-handed twist between
tetrads in the G4 structure.

> **Figure placeholder:** Right-handed and left-handed tetrad stacks.

#### `force_field`

A GROMACS short force-field name or a menu number. Run
`gmx pdb2gmx` to inspect the available force fields and their menu numbers.

**Force-field menu numbers depend on the GROMACS installation.** The number `4`
in the example is an example selection, not a universal force-field identifier.
Choose the appropriate value for your installation. The lower-level API can
display the available menu before selection.

**Optional entries**

#### `gromacs_executable`

**Default:** `"gmx"`.

GROMACS command or executable path.

#### `dielectric_const`

**Default:** `15`.

Dielectric constant written to the simulation molecular dynamics parameter
(MDP) file as `epsilon_r`. Must be positive and finite.

#### `saving_stride`

**Default:** `5000`.

Number of MD steps between compressed trajectory frames and PLUMED output
records. Must be a positive integer.

#### `after_folding_steps`

**Default:** `2500000`.

Number of additional MD steps after the folding protocol. Tetrad restraints
remain active during this period. Must be a nonnegative integer.

#### `max_bulge_length`

**Default:** `3`.

Maximum number of intervening residues (bulge size) allowed between successive
residues in a G-tract reference block. Larger gaps split the G-tract into
separate reference blocks. Must be a nonnegative integer.

#### `n_frames`

**Default:** `5`.

Number of structures to extract from the post-protocol trajectory. Frames
are ranked by RMSD to the full-core reference: frame 1 has the lowest RMSD,
frame 2 has the second-lowest RMSD, and so on. Must be a positive integer.

#### `schedule` — restraint schedule

**Default:** The schedule shown in the example below.

An object containing four configurable timings, all in MD steps. Omitted
fields use their defaults:

- `gt_initial_steps` (default `400000`): Duration of the initial G-tract
  ramp from a force constant of 0 to 12000.
- `core_start` (default `1000000`): Starting step for the first cumulative
  core restraint. Must be greater than `gt_initial_steps`.
- `core_steps` (default `3000000`): Total duration of each core restraint
  block. Must be divisible by three so that the intermediate stages fall
  on integer MD steps.
- `tetrads_steps` (default `3000000`): Total duration of the tetrad restraint
  block. Must be greater than `core_steps / 3`.

All four values must be positive integers. Targets and force constants are
fixed; the package generates the following `[step, target, kappa]` stages:

```text
gt_initial:
  [0,                  0.0,     0]
  [gt_initial_steps,   0.0, 12000]

core:
  [0,                  0.0,     0]
  [core_steps / 3,     0.0, 20000]
  [2 * core_steps / 3, 0.0, 20000]
  [core_steps,         0.0,     0]

tetrads:
  [0,                  0.0,     0]
  [core_steps / 3,     0.0, 20000]
  [tetrads_steps,      0.0, 10000]
```

G-tract stage steps are measured from the simulation start. Core and tetrad
stage steps are offsets from the corresponding restraint block's start.
Successive core blocks start `core_spacing = 2 * core_steps / 3` steps apart,
so they overlap for the final third of each block. Stage arrays and
`core_spacing` are derived values, not input fields.

For `n_core` cumulative core references, tetrad restraints start at:

```text
tetrad_start = core_start + n_core * core_spacing
```

G-tract restraints hold their force constant at 12000 until `tetrad_start`,
then reach zero at `tetrad_start + core_steps / 3`, the second tetrad stage.
The protocol length is the largest step used by any restraint, and the total
simulation length adds `after_folding_steps`.

The generated protocol uses upper-wall `MOVINGRESTRAINT` restraints
(`VERSE=U`); PLUMED interpolates targets and force constants between stage
points. See the [PLUMED MOVINGRESTRAINT documentation](https://www.plumed.org/doc-v2.10/user-doc/html/_m_o_v_i_n_g_r_e_s_t_r_a_i_n_t.html).
The writer does not change PLUMED's default units: RMSD targets are in nm,
time output is in ps, and RMSD restraint force constants use kJ/mol/nm².
PDB coordinates remain in Å. See [PLUMED units](https://www.plumed.org/doc-v2.10/user-doc/html/_u_n_i_t_s.html).

> **Figure placeholder:** Restraint schedule showing G-tract, overlapping
> core and tetrad stages, followed by the post-folding simulation period.

#### Example input JSON

This complete example uses the default optional settings and restraint schedule.
It is also available as [examples/input.json](examples/input.json).

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
  "force_field": 4,
  "gromacs_executable": "gmx",
  "dielectric_const": 15,
  "saving_stride": 5000,
  "after_folding_steps": 2500000,
  "max_bulge_length": 3,
  "n_frames": 5,
  "schedule": {
    "gt_initial_steps": 400000,
    "core_start": 1000000,
    "core_steps": 3000000,
    "tetrads_steps": 3000000
  }
}
```

### Extract from an explicit trajectory file

Use these flags for a simulation run outside the client, or to override the
recorded output paths:

```bash
# Use the trajectory in the run directory.
g4-builder extract -w runs/example -trj smd.xtc

# Specify the structure and COLVAR as well.
g4-builder extract -w runs/example \
    -trj smd.xtc \
    -s smd.tpr -colvar COLVAR
```

| Flag | Input | Default when omitted |
| --- | --- | --- |
| `-trj`, `--trajectories` | XTC trajectory path | Trajectory recorded by `run_sim`; required explicitly if no outputs are recorded |
| `-s`, `--structure` | Matching structure or TPR file | Recorded structure, otherwise the prepared simulation TPR |
| `-colvar`, `--colvar` | COLVAR file | Recorded COLVAR, otherwise the run directory's default COLVAR |

Extraction uses the saved full-core RMSD label and protocol duration. The
number of frames comes from `n_frames` in the saved configuration. All input
files must belong to the same simulation, with matching atom order.

### Path rules

`--input` and `--working-dir` resolve relative to the terminal's current
directory. For extraction inputs, in both the CLI and Python API:

| Supplied path | Resolution |
| --- | --- |
| `smd.xtc` | A bare filename: relative to the run's working directory |
| `results/smd.xtc` | A path containing directories: relative to the caller's current directory |
| `/data/run/smd.xtc` | Absolute path: used directly |

Output frames are always written to the working directory.

### Reuse prepared runs

Each successful stage updates `run.json` with the validated settings, reference
metadata, protocol duration, completed stages and available results. It is
generated run metadata, distinct from the user input JSON.

Repeating `prepare` with matching settings reuses completed preparation stages.
The executable and `n_frames` may be updated without rebuilding the model;
other configuration changes require a new directory or explicit `--overwrite`.
`--overwrite` regenerates input files and run metadata but does not delete
existing MD output files. Use a fresh directory for a different simulation.

`run_sim` uses `-cpi`. GROMACS uses an available checkpoint for continuation
and appends to the existing output files by default. Keep the checkpoint and
output files together when continuing a simulation. See the
[GROMACS mdrun documentation](https://manual.gromacs.org/current/onlinehelp/gmx-mdrun.html).

Use `g4-builder --help` or `g4-builder extract --help` for parser help. Exit
codes are 0 for success, 1 for handled runtime errors, 2 for argument errors,
and 130 for an interrupted command.

## Python API

```python
import g4_builder as g4
```

The pipeline coordinates the same components that can be called individually.
Executable examples are available in the [pipeline notebook](notebooks/api_pipeline.ipynb)
and [step-by-step notebook](notebooks/api_step_by_step.ipynb).

### Configure a pipeline

Choose one of the following configuration methods. Both loaders require
`working_dir` as a separate argument and reject that key inside the dictionary
or JSON file.

**From JSON:**

```python
config = g4.SimulationConfig.from_json(
    "examples/input.json",
    working_dir="runs/api_example",
)
```

**From a dictionary:**

```python
settings = {
    "sequence": "GGGTAGGGCGGCGGGGCAGGGT",
    "GT_resids": {
        "GT_1": [11, 1, 2],
        "GT_2": [8, 7, 6],
        "GT_3": [19, 20, 21],
        "GT_4": [15, 14, 13],
    },
    "polarity_pattern": "ssa",
    "helicity": "RH",
    "force_field": 4,
}
config = g4.SimulationConfig.from_dict(settings, working_dir="runs/api_example")
```

**Using constructor arguments:**

```python
config = g4.SimulationConfig(
    sequence="GGGTAGGGCGGCGGGGCAGGGT",
    GT_resids={
        "GT_1": [11, 1, 2],
        "GT_2": [8, 7, 6],
        "GT_3": [19, 20, 21],
        "GT_4": [15, 14, 13],
    },
    polarity_pattern="ssa",
    helicity="RH",
    force_field=4,
    working_dir="runs/api_example",
)
```

Constructing a config validates its values; it does not create files or run
GROMACS. Configs are immutable; use `dataclasses.replace()` to derive changes:

```python
from dataclasses import replace

config = replace(config, schedule=g4.ScheduleConfig(core_steps=4_500_000))
```

Dictionary and JSON loaders accept a nested schedule dictionary. The direct
constructor expects a `ScheduleConfig` object. `config.to_dict()` exports
JSON-compatible settings without `working_dir`. Use the input JSON example
above when preparing settings for the CLI.

### Run the pipeline

```python
pipeline = g4.FoldingPipeline(config)
g4.configure_logging(pipeline.paths.log)

result = pipeline.prepare()       # Strand, references, protocol and run.json.
pipeline.minimize()               # Topology, reference remapping, box and EM.
pipeline.prepare_simulation()     # Simulation MDP and TPR.

# Run these after preparation with a PLUMED-enabled GROMACS executable.
pipeline.run()
result = pipeline.extract()

print(result.protocol_steps, result.simulation_steps)
print(result.full_ref_label)
print(result.best_frames)         # {time_ps: RMSD}, ranked by RMSD.
print(result.frame_paths)         # {1: Path(...), 2: Path(...), ...}
```

| Method | Behavior |
| --- | --- |
| `prepare(overwrite=False)` | Build and save the strand, full/partial references and PLUMED protocol without calling GROMACS. Refuse to replace existing run inputs unless explicitly allowed. |
| `minimize(pdb2gmx_timeout=60, discovery_timeout=2)` | Generate topology, detect and remap changed reference atom serials, regenerate the protocol, add a box, prepare and run energy minimization. |
| `prepare_simulation()` | Create the simulation MDP and TPR using protocol length plus post-folding steps. |
| `run(timeout=None)` | Run the PLUMED simulation and record the XTC output path. |
| `extract(trajectories=None, structure=None, colvar=None, n=None, timeout=None)` | Select and extract frames using recorded outputs or explicit file overrides. `n` overrides the configured `n_frames`. |

Stages return the same `RunResult` object. It contains paths, protocol and
simulation step counts, the full-core label, reference paths, completed stages,
simulation outputs and extraction results. Intermediate objects are available
as `pipeline.strand`, `pipeline.references`, `pipeline.protocol` and, after
extraction, `pipeline.extractor`. `pipeline.paths.to_dict()` lists all planned
paths and patterns, including files not yet created.

### Restore a pipeline or analyze an external run

```python
pipeline = g4.FoldingPipeline.from_directory("runs/api_example")
g4.configure_logging(pipeline.paths.log)

# Call the next required stage, for example after minimization:
pipeline.prepare_simulation()
```

Restoration reads `run.json` and prepared files without starting GROMACS.
Paths within a run directory are relocated when the entire directory is moved;
external absolute paths retain their original locations. API stages run when
called; unlike CLI preparation, they do not automatically skip completed stages.

For a simulation run separately, supply its trajectory explicitly:

```python
result = pipeline.extract(
    trajectories="smd.xtc",
    structure="smd.tpr",
    colvar="COLVAR",
    n=3,
)
```

### Lower-level API: prepare files step by step

This example uses the individual components without `FoldingPipeline` or an
input configuration file. Choose the force field for your installation before
running the topology step.

```python
import g4_builder as g4

runner = g4.GromacsRunner("runs/manual_example", executable="gmx")
paths = runner.paths
paths.working_dir.mkdir(parents=True, exist_ok=False)
g4.configure_logging(paths.log)

# Build and save the strand.
strand_builder = g4.DNAStrandBuilder()
strand = strand_builder.build("GGGTAGGGCGGCGGGGCAGGGT")
strand_builder.save_pdb(paths.strand)

# Build and save references using ordered, 1-based residue indices.
GT_resids = {
    "GT_1": [11, 1, 2],
    "GT_2": [8, 7, 6],
    "GT_3": [19, 20, 21],
    "GT_4": [15, 14, 13],
}
references = g4.build_references("ssa", "RH", strand, GT_resids)
references.save_full_ref(paths.full_ref)
references.save_partial_refs(paths.working_dir)

# Generate topology and inspect any changes to reference atom numbering.
preparation = g4.GromacsPreparation(runner)
print(preparation.discover_force_fields())
preparation.pdb2gmx(force_field=4)
processed = g4.read_pdb(paths.conf)
if references.needs_serial_remap(processed):
    print("Reference atom serials changed; remapping references.")
    references.remap_serials(processed)
    references.save_full_ref(paths.full_ref)
    references.save_partial_refs(paths.working_dir)

# Add the box and minimize.
preparation.add_box()
preparation.prep_enmin_file()
preparation.run_enmin()

# Generate the protocol using the processed structure's atom serials.
protocol = g4.build_protocol(
    references,
    atom_serials=processed["serial"],
    schedule=g4.ScheduleConfig(),
    print_stride=5000,
    print_file=paths.colvar.name,
)
protocol.write(paths.plumed)

# Render the MDP and prepare the simulation TPR.
preparation.prep_sim_file(
    n_steps=protocol.final_steps + 2_500_000,
    saving_stride=5000,
    dielectric_const=15,
)
```

Each GROMACS preparation method returns `None` and records successful outputs
in `preparation.saved_paths`. `pdb2gmx()` accepts an optional `strand` PDB path;
without it, the input is `paths.strand`. Relative strand overrides resolve from
the runner's working directory. `discover_force_fields(timeout=2, strand=None)`
stores and returns a menu-number-to-name dictionary.

The wrappers use fixed command choices: TIP3P in `pdb2gmx`, a centered
dodecahedron box with `editconf -d 1.4`, `grompp -maxwarn 1` for minimization,
and `grompp -maxwarn 2` for the simulation. Simulation execution always uses
`-cpi`, `-v` and `-plumed`; frame extraction selects `System`.
These are not configuration entries or method options. The packaged MDP files
provide the remaining simulation settings.

### Lower-level API: run and extract

With PLUMED-enabled GROMACS, continue the manual preparation example with:

```python
preparation.run_sim(timeout=None)

extractor = g4.TrajectoryFrameExtractor(runner)
best_frames = extractor.find_best_frames(
    rmsd_label=protocol.full_ref_label,
    protocol_steps=protocol.final_steps,
    colvar=paths.colvar,
    n=5,
)
extractor.extract_frames(
    trajectories=paths.trajectory,
    structure=preparation.outputs.structure,
    timeout=None,
)

print(best_frames)
print(extractor.saved_paths)
```

`TrajectoryFrameExtractor` takes only a runner at construction. Selection reads
COLVAR and does not require a trajectory, structure or GROMACS execution.
`colvar=None` uses the default path. Extraction takes a trajectory path and a
matching structure/TPR path, then runs `trjconv` for each selected time.

`find_best_frames()` returns and stores `{time_ps: RMSD}` in `best_frames`.
Every new selection clears the previous selection and saved paths, including
when selection fails. `extract_frames()` requires a successful selection,
returns `None`, prints saved paths and records them in `saved_paths`.

Manual component calls do not create pipeline `run.json` metadata. Retain the
protocol's `final_steps` and `full_ref_label` if selecting frames
in a later session.

### Lower-level API: structures and references

`DNAStrandBuilder(nucleotide_dir=None)` loads 12 packaged nucleotide templates
for 5′, middle and 3′ positions. `build(sequence, chain_id="A")` aligns
consecutive templates using their phosphate anchors and returns the strand
DataFrame. At least two bases are required by the terminal template set.
Custom nucleotide templates can be supplied through `nucleotide_dir`.

`build_references(polarity_pattern, helicity, strand, GT_resids,
max_bulge_length=3, layer_dir=None)` returns a
`ReferenceSet` containing `full_ref`, `entries`, `partial_refs`, `ref_info` and
`full_ref_label`. Partial references cover G-tract blocks, cumulative cores and
individual tetrads. Template atom names must match the strand's selected base
atoms row for row before coordinates are copied.

Each `ReferenceSet.entries` record contains its name, category, ordered residue
indices, RMSD label and atom DataFrame. Its `is_full_core` flag identifies the
final cumulative core used for frame selection.

Reference templates are always generated by `build_reference_template()` from
the flat `a.pdb` and `s.pdb` layers in
`g4_builder/data/reference_templates/layers/`. This applies to the CLI,
pipeline and lower-level API. Both RH and LH stacks support any a/s pattern
whose length matches the supplied G-tracts; simulation configurations require
at least two tetrads.

```python
references = g4.build_references(
    "ssa", "RH", strand, GT_resids,
)

# Generate a template DataFrame directly, without mapping it to a strand.
template = g4.build_reference_template("sass", "LH")
```

Both `build_references()` and `build_reference_template()` accept `layer_dir`
to supply a custom directory containing the flat `a.pdb` and `s.pdb` files.
This is a lower-level API option; CLI and pipeline preparation use the
packaged layers.

The flat layers must contain residues 1–4 with finite coordinates in the Z=0
plane. Rise and signed twist values are read on every template build from
[g4_builder/data/stack_properties.json](g4_builder/data/stack_properties.json).
Rise is in Å and twist in degrees; adjacent-pair transformations accumulate
relative to the preceding layer. Output residues are grouped by original
residue number across layers before advancing to the next original residue.

`ReferenceSet.save_full_ref(path)` and `save_partial_refs(directory)` save
reference PDBs without TER/END records. Prepared references have
`temp_factor=1.0`, with blank chain and element fields. Partial saving prints
paths and stores them in `saved_paths`. Atom-serial remapping matches
`(res_seq, atom_name)`, preserves coordinates and ordering, and rejects missing
or ambiguous matches.

Standalone PDB I/O is also available:

```python
atoms = g4.read_pdb("structure.pdb")
g4.save_pdb(atoms, "copy.pdb", add_ter_end=False)
```

`read_pdb()` reads ATOM/HETATM records into columns `record_name`, `serial`,
`atom_name`, `alt_loc`, `res_name`, `chain_id`, `res_seq`, `i_code`, `x`, `y`,
`z`, `occupancy`, `temp_factor`, `seg_id`, `element` and `charge`. Coordinates
are in Å. Blank optional numeric fields are represented by NaN.
`save_pdb()` preserves row order and atom/residue numbers; its default
`add_ter_end=True` adds TER and END records.

## Simulation files and frame selection

### Keep simulation inputs together

**All partial reference PDBs, the PLUMED protocol and simulation TPR must be in
the same directory when running the simulation.** Run GROMACS from that
directory. `build_protocol()` writes bare reference filenames, and the pipeline
places the required files together automatically.

COLVAR must have its column labels on the first line, for example:

```text
#! FIELDS time rm1 rm2 rm3
```

The first column is time in ps; the requested RMSD label identifies the
selection column. Selection filters times at or after the end of the folding
protocol, ignores non-finite time/RMSD values, rounds times to
integer ps, and ranks records by increasing RMSD. Duplicate rounded times are
removed, keeping the lowest RMSD. At least `n` distinct eligible times are
required. Frame 1 corresponds to the lowest selected RMSD.

COLVAR and the trajectory must cover the same run.
`trjconv -dump` extracts the saved frame nearest each requested time; matching
trajectory and COLVAR saving strides helps their samples correspond. Selection
does not inspect XTC duration to detect mismatched COLVAR data.

### Default output files

With the supplied [file_names.json](g4_builder/data/file_names.json), the main
artifacts are:

| File | Purpose |
| --- | --- |
| `strand.pdb` | Constructed DNA strand |
| `full_ref.pdb`, `ref_1.pdb`, … | Full and partial reference structures |
| `plumed_protocol.dat` | Generated restraint protocol |
| `conf.pdb`, `topol.top`, `posre.itp` | Topology-preparation outputs |
| `conf_box.pdb` | Boxed structure |
| `em.tpr`, `em.pdb` | Minimization input and minimized structure |
| `smd.mdp`, `smd.tpr` | Simulation parameter file and run input |
| `smd.xtc` | Simulation trajectory |
| `COLVAR` | PLUMED time series used for frame selection |
| `frame_1.pdb`, … | Extracted frames ranked by RMSD |
| `gromacs.log` | Application log when configured by the CLI or API |
| `run.json` | Pipeline state for later stages and restoration |

Only the two MD TPR names are defined for stage naming: `em_tpr` and `tpr`.
Their stems determine `mdrun -deffnm`; related names are derived in code.
For example, changing `tpr` to `fold.tpr` produces a `fold.mdp` input and
`fold.xtc` trajectory. With the default `smd.tpr`, the final structure is
saved as `smd.pdb`.

There is no prefix setting in the input config. Filename patterns use `{name}`
for references and `{index}` for frames.
`RunPaths.to_dict()` returns absolute `Path` values and retains those patterns;
it does not scan for existing files. Use `paths.reference("ref_1")` and
`paths.frame(1)` for specific generated paths.

The package data directory also contains the nucleotide/reference PDBs and
`mdp/em.mdp` and `mdp/SMD.mdp` source templates. Edit maintained package data
in an editable checkout to change package-wide defaults. Keep filename defaults
consistent across all stages and machines for an existing run. Editing the
top-level experimental template directories does not update packaged copies.

## Run a simulation on another machine

1. Complete preparation through the simulation TPR, using CLI `prepare` or
   `pipeline.prepare_simulation()` after minimization.
2. Copy the entire working directory to the simulation machine, keeping the
   TPR, PLUMED protocol, references, run metadata and any checkpoints together.
3. Activate a PLUMED-enabled GROMACS installation and run from inside that
   directory. With the default filenames:

   ```bash
   gmx mdrun -deffnm smd -cpi -v -plumed plumed_protocol.dat -c smd.pdb
   ```

   Use the site's job-submission system when required. The Python package
   is not needed remotely for this direct GROMACS command.
4. Copy `COLVAR` and `smd.xtc` back, preserving the TPR
   used for simulation. Analyze the results locally with explicit paths:

   ```bash
   g4-builder extract -w runs/example \
       -trj smd.xtc \
       -s smd.tpr -colvar COLVAR
   ```

Running GROMACS directly does not record the output trajectory in the pipeline,
so an explicit trajectory path is required for this route. In a fresh Python
session, restore the pipeline with `FoldingPipeline.from_directory()` and use
`pipeline.extract(trajectories=...)`. Extraction still needs GROMACS locally,
but does not need PLUMED support. If GROMACS is unavailable on the local machine,
perform the GROMACS preparation stages on the remote machine as well.

## Logging and troubleshooting

The CLI configures INFO console messages and a detailed `gromacs.log` file.
For the Python API, logging is opt-in:

```python
import logging

g4.configure_logging(
    log_file="runs/api_example/gromacs.log",
    console=True,
    console_level=logging.INFO,
)
```

GROMACS stdout and stderr are streamed to the logger at DEBUG level rather
than returned by command methods. Output follows a labelled block header,
with no timestamps, level labels or stream prefixes added to individual lines.
Indentation, blank lines and trailing whitespace are preserved. File logging
retains the full streams; set `console_level=logging.DEBUG` to show them in
the console too. Repeated
configuration replaces package-owned handlers. Applications can configure
the `g4_builder` logger themselves instead. Without logging configuration,
normal command messages are silent; failures still raise exceptions.
Reference and extracted-frame save messages are printed directly.

| Problem | What to check |
| --- | --- |
| GROMACS executable cannot be launched | Activate its environment before starting Python/Jupyter, or supply the executable path. |
| Simulation rejects `-plumed` | Use a GROMACS executable with PLUMED support, or run the prepared simulation on another machine. |
| Force-field menu number is unavailable | Choose a number from this installation's menu or use an available short name. |
| Reference generation cannot find layer files | Check that the packaged layers are installed, or that the API's `layer_dir` contains both `a.pdb` and `s.pdb`. |
| Extraction cannot find an input | Check the bare-filename versus directory-containing path rules, the run directory and supplied files. |
| Too few eligible frames | Check the RMSD label, protocol duration and COLVAR coverage, or reduce `n_frames`/API `n`. |
| Existing preparation conflicts with new settings | Use the matching settings, a new run directory or explicit preparation `--overwrite`. |
| Notebook reports a config type mismatch after source edits | Restart the kernel, rerun imports and recreate the config; module reloads can leave different class versions in memory. |

API timeouts are supplied to operations, for example
`pipeline.minimize(pdb2gmx_timeout=120, discovery_timeout=5)` or
`pipeline.run(timeout=3600)`. The CLI uses a 60-second topology timeout,
2-second force-field discovery timeout and no simulation/extraction timeout.
Timeout and logging choices are not stored as simulation settings.

## Development

The maintained implementation is in `g4_builder/`, package assets in
`g4_builder/data/`, JSON examples in `examples/`, and example notebooks in
`notebooks/`. `pyproject.toml` defines dependencies, package data and the
`g4-builder` console entry point. Importing `g4_builder` exposes the Python API;
`python -m g4_builder` runs the same client as the installed command.

`build/` and `g4_builder.egg-info/` are generated packaging artifacts. Make
changes in the maintained package source rather than `build/lib/`. Reinstall
after changing installation metadata; restart running kernels when needed
after source changes.

Run the package tests from the project root:

```bash
python -m unittest discover -s tests_new -v
```

The tests cover real template geometry, reference mapping, protocol generation,
configuration and CLI inputs, metadata restoration, extraction and subprocess
logging. GROMACS workflow tests use a fake runner; they do not require a
PLUMED-enabled executable or run molecular dynamics.
