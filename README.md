# Anchor-based ortholog LCR database

Build one FASTA per human LCR and its anchor-mapped ortholog regions, a strict-source subset, and one XLSX with essential UniProt/OMA identifiers. Everything runs on local files. Python 3.10 or newer is required. Uniprot Release 2026_01 and OMA Release All.Jul2024 (31st release) were used to generate the supplied snapshot.

## Run

From this directory in PowerShell or a terminal:

```powershell
python -m pip install -r requirements.txt
python run_pipeline.py --overwrite
```

The default run uses the supplied anchor mappings and writes exactly these outputs:

```text
outputs/
  anchor_fastas_default/    # seg_mode=default; 59,284 groups, 59,284 FASTAs
  anchor_fastas_strict/     # seg_mode=strict; 2,385 groups, 2,385 FASTAs
  metadata.xlsx
```

Strict refers to the human source's SEG detection mode. It does not impose a stricter anchor-coverage threshold. SEG default retains all anchored groups, including strict sources. Both sets contain extracted, unaligned LCR sequences. They are not multiple-sequence alignments.

Each file is named `family_id_source_lcr_id.fasta`. The first record is the human source, followed by mapped orthologs sorted by LCR ID. Identifiers and coordinates remain in the original header format.

The archive already includes verified outputs, so the command above replaces them. Use `--output PATH` to write elsewhere. Existing output directories are refused unless `--overwrite` is supplied. Outputs are staged before replacement so an unsuccessful run preserves previous results. `--family 851` selects one local family; repeat the option for several families.

## Inputs

Four tables under `inputs/tables/` are required. They are snapshots of the original human LCR response and OMA bulk files:

| File | Essential contents |
| --- | --- |
| `protein.csv` | Protein IDs, separate UniProt/OMA IDs, species, available taxids, full sequences |
| `lcr_seq.csv` | Human seed LCRs and anchor-derived ortholog intervals, with explicit coordinate systems |
| `lcr_family.csv` | Local family ID, human source accession, OMA-derived member protein IDs |
| `phylogeny_map.csv` | Family ID, human source LCR ID, mapped ortholog LCR ID |

The data includes 139,257 protein records, 67,712 human seed LCRs, 314,202 mapped ortholog LCRs and 14,820 families. Human seeds without successful mappings remain available as inputs for rebuilding; no orphan-only FASTA files are produced. The 12,814 independent non-human SEG calls are omitted because the anchor mapper does not use them.

The loader validates IDs, family membership, sequences and coordinates before writing outputs. No database, online request, annotation service or raw OMA download is required for the supplied snapshot.

## Essential metadata

`metadata.xlsx` has two flat, filterable sheets:

- **Groups:** local family ID, human source LCR ID/accession/mode, sequence/species counts, FASTA filename.
- **Proteins:** protein ID, UniProt accession, OMA ID, species, available taxid and protein length in amino acids. Only proteins represented in the output FASTAs are included. Join to each FASTA using `pid` in its headers.

Non-human identifiers are in `oma_id`; unavailable UniProt accessions stay blank. Human OMA IDs, missing taxids, gene names, protein names and original OMA group IDs that could not be recovered stay blank. `family_id` is the existing local database identifier, not an OMA group identifier.

If a local OMA/UniProt crosswalk becomes available, fill the identifier gaps without changing FASTA contents:

```powershell
python run_pipeline.py --oma-uniprot "C:\data\oma-uniprot.txt.gz" --overwrite
```

The crosswalk must have OMA protein ID and UniProt accession as its first two whitespace-separated columns; `#` comments are allowed. Multiple matches are retained in sorted, semicolon-separated cells. This is an optional local input, not a network lookup.

## Rebuild anchors

The anchor-alignment implementation is in `pipeline/anchors.py`. To recompute mappings from the cached human LCRs, full ortholog proteins and OMA-derived family memberships:

```powershell
python run_pipeline.py --rebuild-anchors --output outputs_rebuilt
```

The retained method uses 100-residue human flanks, local alignment with BLOSUM62 for mammals/BLOSUM45 for other species, gap scores -11/-1, and minimum query-flank coverage 0.50. Candidate intervals must be 1–500 residues. Two accepted anchors bound the interval. With only one accepted anchor, the fallback spans up to `2 * source_length + 1` residues, clipped at the protein boundary.

An anchor query or target containing a residue unsupported by those matrices, such as selenocysteine (`U`), produces no hit. Biological sequences are retained unchanged; residues are not substituted to force an alignment.

## Code and checks

```text
run_pipeline.py
inputs/
  pipeline/
    oma/                      # OMA bulk files and extracted target-species members of the human-seed-containing groups
    anchor_map.py             # Align human LCR flanks against ortholog proteins and record accepted intervals
    config.py                 # Hard-coded parameters for anchor mapping and FASTA export
    seed_lcrs.py              # Save the LCRAnnotationsDB human LCR response and its exact request receipt
    seed_orthologs_bulk.py    # Save the OMA bulk files and extract the target species members of the human-seed-containing groups
    upstream.py               # Build the four local input tables from saved human LCRs and OMA bulk files
  tables/
    inputs.zip                # Compact tables; unzip before running
  fetch_human_lcrs.py   # Save the LCRAnnotationsDB human LCR response and its exact request receipt
  prepare_inputs.py     # Build human seeds and OMA protein/family CSVs from saved upstream files

pipeline/
  data.py           # Local inputs, integrity checks, optional ID crosswalk
  anchors.py        # Anchor mapping
  export.py         # FASTAs and metadata XLSX
tests/test_pipeline.py
requirements.txt
inputs.zip          # Compact inputs; unzip to `inputs/` before running
outputs/            # Generated; ignored by git
```

Run the synthetic checks with:

```powershell
python -m unittest discover -s tests -v
```

## Input provenance

The supplied input CSV files are snapshots. Their preparation involved three stages:

### 1. Human proteins and LCRs

The original `seed_lcrs.py` queried the LCRAnnotationsDB endpoint `https://lcrannotdb.lcr-lab.org/api/proteins/`, selecting `Homo sapiens`. The response supplied UniProt accessions, full protein sequences, LCR sequences, coordinates and identification methods.

The identification-method field was converted into `seg_mode`, including `default` and `strict`. These human LCRs are precomputed annotations; the current pipeline does not rerun SEG to obtain them.

### 2. OMA ortholog proteins and families

The original `seed_orthologs_bulk.py` used three OMA bulk files:

- `oma-uniprot.txt.gz`: mapped human UniProt accessions to OMA protein IDs.
- `oma-groups.txt.gz`: identified OMA groups containing the human seeds and selected target-species members.
- `oma-seqs.fa.gz`: supplied full sequences for the selected ortholog proteins.

The target species were chimpanzee, macaque, mouse, rat, cow, dog, pig, opossum, chicken, Xenopus and zebrafish. The ingestion script retained the first encountered member per target species within each selected group.

Protein records and family memberships were stored using local database IDs. These local family IDs are not original OMA group IDs.

### 3. Anchor-derived regions and table exports

The `anchor_map.py` script aligned human LCR flanks against the ortholog proteins. Accepted intervals were stored as ortholog LCRs, with their human-source and family relationships recorded in `phylogeny_map`.

The resulting `protein`, `lcr_seq`, `lcr_family` and `phylogeny_map` tables were exported to CSV. The supplied inputs retain the fields needed for FASTA generation, essential metadata and anchor rebuilding.

### Reproducibility boundary

The current pipeline reproduces FASTA collections from these saved tables. `--rebuild-anchors` recomputes anchor mappings using the saved human seeds, ortholog sequences and family memberships; it does not retrieve the upstream databases again.

The original raw API responses and OMA bulk files are not included in this repository. Upstream release versions, retrieval dates and file checksums are not encoded in the CSV exports. Consequently, the supplied tables support reproduction from the saved snapshot, while reconstruction from the original database releases requires those additional inputs pipeline run.
