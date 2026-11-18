# EMC perturbation proteomics

Every published proteome in which a subunit of the ER membrane protein complex (EMC) was knocked out or knocked down, re-analysed on common terms, with a membrane-protein filter and links to human orthologues.

Site: https://xiaozzxi.github.io/emc-atlas/

## What is here

- `index.html`, `dataset.html`, `protein.html`, `access.html`, `methods.html`: the pages (plain HTML, CSS and JavaScript; no build step).
- `assets/`: styles and scripts.
- `data/`: the JSON the pages read, written by `downloads/code/export_site.py`.
- `downloads/`: all values as CSV (`emc_long.csv`, `conditions.csv`, one CSV per dataset in `reanalysis/`), the original supplementary files in `original/`, and the scripts in `code/`.

## Rebuilding the data

From a folder holding the original files and the scripts:

```
python build.py <folder with the original files> <output folder>
python export_site.py <output folder> <folder with the original files> <this site folder>
```

`build.py` writes `emc_long.csv` and `conditions.csv`; `export_site.py` writes `data/` and `downloads/` for the site. Dataset page text lives in `datasets_meta.py`.

## Previewing locally

The pages load their data with `fetch`, which browsers block for files opened directly from disk. From this folder run

```
python -m http.server 8000
```

and open http://localhost:8000.

## Licence

Code: MIT. Our re-analysed tables (`downloads/*.csv`, `downloads/reanalysis/`, `data/`): CC BY 4.0. Files in `downloads/original/` are copied unchanged from the papers' supplements and repositories and remain under their publishers' terms. See `LICENSE`.
