# Non-code data census

This is an initial repository-wide census to support data provenance and
licensing review. It is not a determination that any particular file is
licensed or unlicensed. A file's presence in the repository, an EnergyPlus
license header, or a source-code license check does not establish the rights
to redistribute data embedded in that file.

## Repeatable inventory

Run `python scripts/dev/non_code_data_census.py` from any directory in the
checkout to refresh [`non-code-data-inventory.csv`](non-code-data-inventory.csv).
The CSV records path, suffix, byte count, and SHA-256 for non-source files
tracked by Git. Untracked working-tree files are excluded. Source file suffixes
are intentionally excluded. The inventory catches committed additions,
removals, renames, and content changes when reviewed as a diff; hashes do not
record provenance, permitted uses, or the file's role. It is an inventory aid,
not a substitute for the data record described below.

## Initial findings

- The repository contains roughly 11,700 tracked files. The initial scan found
  more than 4,400 non-source files, spanning text, images, spreadsheets,
  weather files, compiled libraries, FMUs, and other binary assets.
- `testfiles/` contains over 1,000 files and is broadly installed as
  `ExampleFiles`; CMake excludes names beginning with `_`. These examples
  include input data and may embed third-party curves, schedules, material
  properties, geometry, or other reference data.
- `datasets/` has 93 files/subtree entries, many explicitly installed into the
  release's `DataSets` or `MacroDataSets` directories. The install list includes
  performance curves, reference properties, escalation factors, photovoltaic
  data, and packaged FMUs.
- `weather/` has 97 files/subtree entries. Only selected weather triplets are
  installed by the visible CMake list; the rest may still be used by tests,
  examples, or other workflows.
- `third_party/`, `bin/`, `design/`, `doc/`, `src/`, and `tst/` also contain
  non-source assets. The inventory intentionally includes these rather than
  assuming that only files in `datasets/` are data or that all non-code assets
  are shipped in the same way.

Counts above are a point-in-time orientation, not policy boundaries. The CSV
is the exact file-level snapshot for this branch. The packaging observations
come from `cmake/Install.cmake`; other packaging paths and generated release
contents should be checked before treating that file as the complete release
surface.

## Recommended provenance record

For data that is externally sourced or contains factual/empirical data, record
at least:

| Field | Purpose |
| --- | --- |
| Repository path (or a narrow path pattern) | Identifies the files covered |
| Data role and consumers | Explains what the data represents and how code, tests, examples, or releases use it |
| Origin/source and version or retrieval date | Makes the upstream source traceable |
| Creator/attribution | Preserves required credit |
| License, terms, or public-domain basis | Records the redistribution basis and any conditions |
| Modifications and provenance notes | Explains transformations, derivation, and unresolved questions |
| Release status | States whether the file is included in a package or only used in development/tests |

Keep this human-reviewed provenance record separate from the generated hash
inventory. Group rows only when the same source, terms, transformation, and
release status genuinely apply to every covered file. Record "unknown" and
assign an owner for follow-up rather than inferring a license from file type,
comments, or location.

## Suggested review order

1. Release-distributed files: `testfiles/`, explicitly installed datasets and
   weather, plus assets copied by other release scripts. Trace file lists from
   `cmake/Install.cmake` and `.github/workflows/release_*.yml`.
2. Data-rich inputs: dataset IDFs and example/test IDFs that contain curves,
   performance tables, materials, schedules, or measured reference values.
   Check in-file comments and historical/source documentation, then verify
   attribution and redistribution terms at the upstream source.
3. Weather and geographic data: EPW/DDY/STAT files, maps, and associated
   metadata. Record the dataset provider and terms; filenames alone are not
   adequate provenance.
4. Other binaries and media: FMUs, libraries, images, spreadsheets, and
   archives. Review their own notices and packaging consumers.

An eventual automated check should require a provenance entry for new or
changed files in governed data areas and flag files lacking a reviewed status.
It should not claim to decide whether a license is legally sufficient. The
existing `scripts/dev/license-check.py` validates license headers for C/C++
and Python source; it does not scan data files or identify their contents.
