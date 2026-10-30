# Installation, package structure, API and CLI

This guide summarizes how this project becomes an importable Python library
and a terminal application. Simulation usage is described in [README.md](README.md).
The current JSON fields, extraction flags and API argument rules are documented
in [INPUTS.md](INPUTS.md).

## 1. The three names

| Role | Name | Defined by | Usage |
| --- | --- | --- | --- |
| Distribution (installable project) | `g4-builder` | `[project].name` in `pyproject.toml` | Installation tools and package metadata |
| Python import package | `g4_builder` | The `g4_builder/` directory and package discovery configuration | `import g4_builder as g4` |
| Terminal command | `g4-builder` | `[project.scripts]` in `pyproject.toml` | `g4-builder prepare ...` |

These names are related by our choices, but they serve different purposes.
Changing the distribution name does not automatically rename the import
directory or the terminal command. The terminal command name could also differ
from the distribution name.

The spelling `g4_builder` is appropriate for a Python identifier. The hyphen in
`g4-builder` is suitable for a distribution or terminal command, but not for an
ordinary Python import statement.

## 2. What you edit and what tools generate

```text
G4_builder/                         Project root; not itself the import package
├── pyproject.toml                  Build/install configuration
├── README.md                       Usage documentation
├── g4_builder/                     Maintained package source
│   ├── __init__.py                 Public Python API
│   ├── __main__.py                 Entry for python -m g4_builder
│   ├── cli.py                      Terminal argument parsing and dispatch
│   ├── client_input.py             CLI JSON settings and fixed-option checks
│   ├── config.py                   Validated simulation/schedule settings
│   ├── pipeline.py                 Workflow orchestration and run metadata
│   ├── runner.py                   Shared GROMACS subprocess execution
│   ├── template_builder.py         Build reference templates from flat layers
│   ├── ...                        Strand, reference, protocol, I/O modules
│   └── data/                      Maintained, distributable templates
├── examples/input.json             Complete example client input
├── tests_new/                      Tests for this package
├── build/                         Generated build workspace
│   ├── lib/g4_builder/             Generated copies of package code and data
│   └── bdist.linux-x86_64/         Temporary distribution staging directory
└── g4_builder.egg-info/            Generated setuptools metadata
```

Edit `g4_builder/`, not `build/lib/g4_builder/`. Setuptools copies the selected
modules and templates into a build workspace while assembling an installable
distribution. Those generated copies can be stale and can be overwritten by
the next build. They are not continuously synchronized source files.

The package's `.py` files remain Python source; this build does not compile
GROMACS or turn the application into a native executable. A wheel such as
`g4_builder-0.1.0-py3-none-any.whl` is an installable archive containing package
files and metadata. `build/lib` is a staging location, not the environment's
`site-packages` directory.

The templates in `g4_builder/data/` are different: they are maintained source
assets included in the package. They were copied from the original project's
data directories so the new implementation could be independent while leaving
the old files intact. Editing a top-level legacy template does not update its
package copy. For the new implementation, maintain the files in
`g4_builder/data/`.

`build/` and `*.egg-info/` are normally excluded from version control. They
contain no simulation results and should not be edited manually.

## 3. What pyproject.toml declares

The current [pyproject.toml](pyproject.toml) contains:

```toml
[build-system]
requires = ["setuptools>=61"]
build-backend = "setuptools.build_meta"

[project]
name = "g4-builder"
version = "0.1.0"
description = "DNA strand, G-quadruplex reference and GROMACS folding preparation"
requires-python = ">=3.9"
dependencies = ["numpy>=1.23", "pandas>=1.5"]
readme = "README.md"

[project.scripts]
g4-builder = "g4_builder.cli:main"

[tool.setuptools.packages.find]
include = ["g4_builder*"]
namespaces = false

[tool.setuptools.package-data]
g4_builder = ["data/*.json", "data/nucleotides/*.pdb", "data/reference_templates/*/*.pdb", "data/mdp/*.mdp"]
```

- **`build-system`** selects setuptools as the build backend. Pip is a frontend
  that asks the backend to build the package. These requirements are build
  tools, separate from the application's runtime dependencies.
- **`project`** declares the distribution's identity, version, Python
  compatibility, runtime Python dependencies and description source.
- **`project.scripts`** tells the installer to create a launcher named
  `g4-builder`. `g4_builder.cli:main` means “call `main` in `g4_builder.cli`”.
- **`packages.find`** selects the intended Python packages for distribution.
  `namespaces = false` uses ordinary package discovery rather than implicit
  namespace-package discovery.
- **`package-data`** includes non-Python templates needed after installation.
  Without them, installed Python code could be present but unable to find its
  PDB/MDP inputs.

This file describes the application and its installation. A simulation input
JSON describes a particular run. They are separate configurations.

The dependency entries specify minimum versions, not a fully pinned environment.
GROMACS and PLUMED are external software and are not installed by these Python
dependency declarations. They must be available separately, for example through
your Spack environment.

## 4. Regular versus editable installation

Run installation commands from the project root using the intended environment.

### Regular installation

```bash
conda run -n jupyter_env python -m pip install .
```

Pip builds an installable distribution and installs its contents into the
environment. Later source edits do not automatically update that installed
copy; reinstall to update it.

### Editable installation for development

```bash
conda run -n jupyter_env python -m pip install -e . --no-build-isolation --no-deps
```

- **`-e .`** installs the current project in editable mode. Imports are connected
  to the source tree, so ordinary Python source edits generally do not require
  reinstalling.
- **`--no-build-isolation`** uses build tools already installed in the environment.
- **`--no-deps`** skips dependency installation; NumPy and pandas must already be
  available.

A modern editable installation normally creates a directory such as
`g4_builder-0.1.0.dist-info/` in `site-packages`, plus an import mechanism such
as a `.pth` file and possibly an import helper. It need not copy the full
`g4_builder/` source directory into `site-packages`. Exact helper names and
mechanisms depend on the backend and version.

Both installation modes create the declared terminal launcher in the
environment's executable directory, normally `jupyter_env/bin/g4-builder` on
this system. Activate the environment to make that command available on PATH:

```bash
conda activate jupyter_env
g4-builder --help
```

After changing entry points, dependencies or other installation metadata,
reinstall even when using editable mode. After changing Python code, a running
notebook may need module reloading or a kernel restart because it retains
objects already imported into memory.

During the earlier packaging verification, wheels were installed into temporary
directories under `/tmp`; an editable installation into `jupyter_env` was not
performed by the assistant. The generated `build/` and `.egg-info/` directories
alone do not prove that the package was installed into that environment.

## 5. Why imports can work before installation

Python searches locations listed in `sys.path`. Your project root can be one
of those locations simply because of how Python or the notebook kernel started:

| Startup method | Typical initial search location |
| --- | --- |
| Interactive Python or `python -c ...` | Current directory, represented by `''` |
| `python script.py` | Directory containing the script |
| `python -m g4_builder` | Current directory |
| Notebook | Kernel/IPython configuration, commonly including the notebook directory |

If `/home/michal/Symulacje/G4_builder` is searchable, Python finds its
`g4_builder/` child without an installation. This does not create a terminal
launcher or make the source discoverable from every other working directory.

Inspect the actual notebook/process rather than assuming which copy it imports:

```python
import os
import sys
import g4_builder

print(sys.executable)
print(os.getcwd())
print(sys.path)
print(g4_builder.__file__)
```

An empty `sys.path` entry means the current working directory. Local source
directories can also shadow installed copies, so the reported `__file__` path
is especially useful when debugging imports.

## 6. Python API: __init__.py

When you run:

```python
import g4_builder as g4
```

Python initializes the package by executing `g4_builder/__init__.py`, normally
once per interpreter session unless it is reloaded. Our file exposes the public
API through imports such as:

```python
from .config import SimulationConfig, ScheduleConfig
from .pipeline import FoldingPipeline, RunResult
```

This enables both styles:

```python
import g4_builder as g4
pipeline = g4.FoldingPipeline(config)

# Or import directly from the implementation module:
from g4_builder.pipeline import FoldingPipeline
```

The leading `.` in `.config` means an import relative to this package.
`__all__` lists the names exported by `from g4_builder import *` and documents
the intended public API. It does not prevent access to other modules.

The initialization file also contains:

```python
import logging as _logging
_logging.getLogger(__name__).addHandler(_logging.NullHandler())
```

This attaches a silent handler to the package logger. Importing the library
does not choose a log file or configure console output; callers can use
`configure_logging()`. Importing the package does not start GROMACS, parse CLI
arguments or create simulation files.

## 7. Terminal client: __main__.py and cli.py

`g4_builder/__main__.py` contains:

```python
from .cli import main
raise SystemExit(main())
```

It is the entry for:

```bash
conda run -n jupyter_env python -m g4_builder --help
```

The installed `g4-builder` launcher calls the same `cli.main()` function because
of `[project.scripts]`. The two paths share one parser and implementation:

```text
Python API:        import g4_builder → __init__.py → library classes/functions
Module execution:  python -m g4_builder → __main__.py → cli.main()
Installed client:  g4-builder → generated launcher → cli.main()
```

`cli.py` handles arguments, JSON input, mode selection and errors. Business
logic stays in the pipeline and component modules. `SystemExit(main())` passes
the returned status to the operating system: 0 for success, 1 for handled
runtime failures, 2 for argument-parser errors, and 130 for handled keyboard
interruptions.

Use `python -m g4_builder`, rather than executing `g4_builder/__main__.py` as a
file, so its relative imports have the correct package context.

Example client usage after installation:

```bash
g4-builder prepare_nogmx -i examples/input.json -w testing/new_run
g4-builder prepare       -i examples/input.json -w testing/new_run
g4-builder run_sim                            -w testing/new_run
g4-builder extract                            -w testing/new_run
```

`run.json` holds per-simulation state used by later modes. It is unrelated to
installation metadata.

## 8. What .egg-info contains

| Generated file | Purpose |
| --- | --- |
| `PKG-INFO` | Name, version, requirements and description, including README text |
| `SOURCES.txt` | Inventory of files selected for the source distribution |
| `requires.txt` | Runtime dependency declarations |
| `entry_points.txt` | Terminal command mapping: `g4-builder = g4_builder.cli:main` |
| `top_level.txt` | Top-level import package name: `g4_builder` |
| `dependency_links.txt` | Legacy dependency-location metadata; currently effectively empty |

Change `pyproject.toml` or the maintained files, rather than editing generated
metadata. The name `.egg-info` is a setuptools convention and does not mean we
are distributing an old-style egg. Modern installed wheels use `.dist-info`
metadata in the environment.

## 9. Standard logging versus package logging

With normal package imports:

```python
import logging                 # Absolute import: standard library
from . import logging          # Relative import: g4_builder/logging.py
import g4_builder.logging      # Explicit absolute import of our module
```

`import logging as _logging` merely gives the imported standard module a local
alias; the alias does not change import resolution. Absolute imports search
`sys.path`, rather than automatically looking in the importing file's package
directory.

If the inner `g4_builder/` directory itself is placed directly on `sys.path`,
its `logging.py` can shadow the standard library. Normal package imports and
module execution avoid that particular setup. Renaming our module to
`logging_config.py` remains a possible clarity improvement, but has not been
done.

## Further reading

- [Python Packaging User Guide: pyproject.toml](https://packaging.python.org/en/latest/guides/writing-pyproject-toml/)
- [Pip: regular and editable local installs](https://pip.pypa.io/en/stable/topics/local-project-installs/)
- [Python: __main__](https://docs.python.org/3/library/__main__.html)
- [Setuptools: pyproject.toml configuration](https://setuptools.pypa.io/en/latest/userguide/pyproject_config.html)
