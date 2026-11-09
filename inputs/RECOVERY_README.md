# ChaMM recovered data

The source archive was recovered from the official QMOF Figshare article 13147324, version 18, file 59573735. The published archive checksum was verified before recovery.

## Included

- `data/chamm_atom_charges_recovered.csv.gz`: all atoms in the original three-charge cohort, with `qmof_id`, `name`, zero-based `atom_index`, `element`, `ddec6`, `cm5`, `bader`.
- `data/chamm_common_cohort.csv`: structure IDs, names and atom counts.
- `data/recovery_report.json`: observed counts, source version/checksum and comparison with the retained benchmark.
- `inputs/expanded_benchmark_records.json`: complete retained structures for the 34-structure connectivity benchmark.
- `inputs/selection.json`: benchmark selection history.
- `inputs/figshare_article_v18.json`: publisher metadata snapshot.
- `code/13_recover_charge_data.py`: the script used to rebuild the table.
- `ChaMM_restart_package.zip`: the earlier saved code and initial connectivity audit.

The charge table is a minimal reconstruction, not a byte-for-byte restoration of the lost CSV schema. It retains source site order and the complete three-scheme cohort. Pairwise differences, ranges and standard deviations can be recalculated from its three numeric charge columns.

## Read the recovered table

Pandas automatically decompresses it:

```python
import pandas as pd
for chunk in pd.read_csv("data/chamm_atom_charges_recovered.csv.gz", chunksize=200_000):
    print(len(chunk))
```

No DFT, connectivity inference, descriptor construction, modeling or scientific association test was performed by this recovery step. The final graph and descriptor table still need validation. The lost local descriptor CSV was provisional and should not be reused as a final result.

## Original public archive

Official landing page: https://figshare.com/articles/dataset/QMOF_Database/13147324

Verified file URL: https://ndownloader.figshare.com/files/59573735

File: `qmof_database.zip`; size: 392,088,304 bytes; MD5: `0d89aaf66f2c306e86e47fd91cb1346e`.

Inside it, the source member is `qmof_database/qmof_structure_data.json` (3,153,444,757 bytes uncompressed). This archive is not embedded in the recovery package because it is large and publicly retrievable. The geometry/charge member was streamed directly from the ZIP during recovery, without extracting the 3.15 GB JSON to disk.

The public source includes full structure information for rebuilding local descriptors later. The 34-record benchmark alone is not a replacement for that full structure data.

## Reproduce recovery

From this directory after downloading the public archive to it:

```powershell
python -m pip install ijson
python code\13_recover_charge_data.py --archive qmof_database.zip --metadata inputs\figshare_article_v18.json --benchmark inputs\expanded_benchmark_records.json --outdir recovered_again
```

The output directory must not already exist. A cohort or benchmark mismatch prevents finalization of the reconstructed table; no full-cohort changes are accepted silently.

## Attribution

QMOF source: Andrew S. Rosen and collaborators, official QMOF database (Figshare DOI 10.6084/m9.figshare.13147324, version 18). Official source documentation: https://github.com/Andrew-S-Rosen/QMOF . The authors make the database available under CC BY 4.0. Changes here: selection of structures with all three finite charge schemes and export of a minimal atom-level charge table; no charge values are fitted or recalculated. Cite the original database papers when using these data scientifically.
