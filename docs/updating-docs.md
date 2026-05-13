# Keeping docs up to date

## When you add or change a parameter

Follow these steps whenever a parameter is added, removed, or renamed in `nextflow.config`.

### 1. Update `nextflow.config`

Add or modify the entry in the `params {}` block.

### 2. Rebuild the schema

Run the nf-core schema builder from the repo root. This is an interactive TUI — it reads `nextflow.config`, diffs against `nextflow_schema.json`, and prompts you to accept additions, removals, and description edits:

```bash
nf-core pipelines schema build
```

Work through the prompts:
- Accept new parameters and fill in their descriptions.
- Confirm removal of parameters that no longer exist.
- Save when done - this writes `nextflow_schema.json`.

### 3. Regenerate the README parameters table

Generate updated markdown from the schema and write it to a temp file:

```bash
nf-core pipelines schema docs --output docs/params.md
```

Open `README.md` and replace the entire `## Parameters` section (from `## Parameters` down to the next `##` heading) with the contents of `docs/params.md`. Then delete the temp file:

```bash
rm docs/params.md
```
