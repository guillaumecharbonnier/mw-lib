MW Snakemake layout and conventions

Purpose

This document describes the conventions and organization used for Snakemake "mw-*" metaworkflows in this workspace so future agents and humans can find and modify workflow components reliably.

Overview

- Each workflow lives under an mw-<name> directory. Many projects include a local mw-lib/ and there is a shared top-level mw-lib/ used as the canonical metaworkflow entry point.
- The canonical Snakemake entry points are Snakefile files (often at mw-lib/Snakefile or mw-<name>/.../Snakefile) that include smaller .smk/.rules/.snake fragments under src/snakemake.

Typical layout (per project)

mw-<name>/
  src/snakemake/
    rules/           # .smk or .rules rule fragments included by the top-level Snakefile
    functions/       # python/.snake helpers to map files/wildcards
    imports.py       # python imports used by Snakefiles
    variables.py     # lists of ids, parameters, sometimes variables.snake
    tables/          # TSVs (e.g. *ids.tsv) used to populate mwconf['ids']
    hooks/
      onstart/
      onsuccess/
      onerror/
    envs/            # conda/env YAMLs (R/Python) referenced by rules
    scripts/         # helper scripts called by rules (R/Python/SH)
  mw-lib/            # per-project library (optional)
  out/               # produced outputs (reports, aggregates) - workflow targets live here
  inp/ or data/      # (convention) input/raw data may be referenced here or provided via -data repo

Shared metaworkflow (mw-lib/Snakefile)

- The top-level Snakefile uses glob() to include components from ../mw*/src/snakemake: imports.py, wildcard_constraints.smk, variables.{py,snake}, functions/*, rules/*.smk, and tables/*ids.tsv.
- It constructs an mwconf dict (mwconf['ids'], mwconf['targets']) used by rules and scripts; this avoids passing huge config objects into script: blocks.
- Hooks (onstart/onsuccess/onerror) are discovered by glob and executed by the top-level Snakefile.
- The top-level Snakefile assumes workflows and data are organized as sibling mw-* directories so relative globs (../mw*/...) resolve consistently.

inp / out folders: principle and usage

- Principle: separate raw inputs from generated outputs.
  - inp/ (or data/) contains raw inputs, imports, symlinks to data repos, or small reference files. It should not contain derived files.
  - out/ contains all generated outputs: intermediate files, final reports, aggregated tables, and rendered HTML (examples seen: out/aggregate_coverage_ratio/..., out/quarto/render/...).
- Rules and the top-level rule all frequently reference out/... as final targets. Keeping outputs under out/ makes it easy to clean, snapshot, or publish results.
- Recommended practice: do not write raw inputs into inp/ from workflows. Keep inp/ immutable for reproducibility; write everything generated to out/.
- Many Snakefiles and rules assume these conventions; search for "out/" targets when adding or debugging rules.

Public mw-lib vs per-project mw-lib

- Public/shared mw-lib (repository-root mw-lib or top-level mw-lib/Snakefile) acts as the canonical metaworkflow driver that collects rule fragments and helper modules from all sibling mw-* projects. It is the usual entry point for running multi-project analyses.
- Per-project mw-lib (mw-<name>/mw-lib/) holds reusable rules, scripts, and envs specific to that project. This allows project-level customization while keeping shared behavior centralized.
- Logic: shared mw-lib holds the orchestration and includes fragments from per-project src/snakemake directories. Per-project modules expose rules, variables, and functions that the shared mw-lib can include via relative globs.

Data repos (mw-<name>-data) and symlink organization

- Common pattern: separate code and large data into sibling repositories/directories.
  - Code: mw-<name>
  - Data: mw-<name>-data
- This separation keeps code repos small and avoids committing large binary/FASTQ/BAM files.

Symlink conventions and how they make workflows seamless

- To make the metaworkflow and per-project Snakefiles find data transparently, create sibling directories and symlinks so expected relative paths resolve. Typical setups:
  - Workspace contains: mw-foo/ and mw-foo-data/ as siblings. Inside mw-foo, create a symlink named inp or data pointing to ../mw-foo-data (e.g. mw-foo/inp -> ../mw-foo-data).
  - Alternatively, create explicit symlinks for specific subfolders (e.g. mw-foo/src/snakemake/tables -> ../mw-foo-data/tables).
- Because the top-level Snakefile uses relative globs like ../mw*/src/snakemake, having mw-foo and mw-foo-data as siblings ensures includes and table-loading work without modifications.
- ID mapping tables (src/snakemake/tables/*ids.tsv) commonly contain absolute or repository-relative paths that point into the -data repo. Keeping these paths consistent (or using symlinks) ensures mwconf['ids'] resolves correctly.

Practical recommendations

- Keep raw data in mw-<name>-data and commit only small index files or path tables (tables/*ids.tsv) into the code repo. Use tables/*ids.tsv to map logical ids to file paths in the data repo.
- Create a simple script or README describing the symlink setup required for each code repo; e.g.: ln -s ../mw-foo-data mw-foo/inp
- Ensure the working directory when running the top-level Snakefile matches expected relative paths (often repository root or mw-lib/ depending on invocation). Document the recommended run location in project README.
- Use out/ for all generated files so downstream report consumers always know where to look and cleaning is straightforward: rm -rf out/

Troubleshooting

- If includes fail, check the relative glob patterns in mw-lib/Snakefile and confirm sibling mw-* directories exist and are accessible.
- If an id in mwconf['ids'] points to a non-existing path, check tables/*ids.tsv and the symlink structure to the -data repo.
- If rules try to write into inp/ or the -data repo, rework the rule to write into out/ instead to preserve immutability of inputs.

If anything in this structure changes, update this MW.md so agents can continue to rely on these conventions.
