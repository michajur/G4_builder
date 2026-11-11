"""Support python -m g4_builder as well as the installed entry point."""

from .cli import main

raise SystemExit(main())
