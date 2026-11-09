# Claim-specific primary-source review

Publisher records, publisher-deposited DOI metadata, and author-hosted
originals were checked during release preparation. Raw metadata are in
`reference_metadata.json` and `figshare_live_metadata.json`. Some publisher
full texts returned HTTP 403; this is not recorded as a successful full-text
check. All cited bibliography keys resolve.

| Claim/reference | Primary evidence and outcome |
|---|---|
| Prior QMOF charge comparison | [Publisher](https://www.nature.com/articles/s41524-022-00796-6), [author-hosted manuscript](https://perssongroup.lbl.gov/papers/rosen_high_throughput_2022.pdf): scheme differences can exceed the tested functional changes. |
| Adsorption motivation | [Publisher article](https://link.springer.com/article/10.1007/s10450-019-00187-2) supports framework-dependent charge sensitivity, particularly at open metal sites. Its [correction](https://link.springer.com/article/10.1007/s10450-021-00301-3) is now cited; no individual adsorption number is imported. |
| DDEC6 | Publisher metadata for 10.1039/C6RA04656H confirm identity. The [primary correction PDF](https://pubs.rsc.org/en/content/articlepdf/2022/ra/d2ra90050e) corrects equation 25 and is now cited. This study inherits charges rather than implementing DDEC6. |
| CM5 and Bader | Publisher-deposited records for 10.1021/ct200866d and 10.1088/0953-8984/21/8/084204 confirm identity. Citations support definitions, not new calculations. |
| PACMAN | [Publisher-deposited metadata](https://api.crossref.org/works/10.1021/acs.jctc.4c00434) give **5368–5380**, correcting supplied **5352–5367**. Preserve published “Predicter.” The [authors' software record](https://github.com/Chung-Research-Group/chung-research-group.github.io/blob/main/Software%20%26%20Data.dc.html) confirms structural charge prediction. The correction DOI 10.1021/acs.jctc.4c00976 resolves, but its deposited record lacks final volume/pages; unsupported supplied page 7384 is replaced by DOI and publication date. No PACMAN accuracy number or correction-specific technical claim is used. |
| pymatgen, ASE, local environments, MOF identifiers | Publisher DOI records match titles, authors, volumes and pages. The [author-hosted local-environment paper](https://escholarship.org/uc/item/4qr765x8) provides the methods reference. Exact algorithm parameters come from saved code and decisions. |
| XGBoost | [Author preprint](https://arxiv.org/abs/1603.02754) and KDD DOI 10.1145/2939672.2939785 confirm the algorithm citation. |
| Source release | Live Figshare version-18 metadata confirm file 59573735, size 392088304, MD5 and CC BY 4.0. The local source archive matches the published checksum. |

A targeted overlap search for conditional charge-definition disagreement in
MOFs found established charge comparisons and predictors; this is not proof
of exhaustive novelty. No first-ever claim is made. No citation for an
unperformed GNN, DFT, redox or causal analysis was introduced.

This ledger supersedes the supplied reference notes' inaccurate PACMAN
pagination and unverified correction pagination. The source package is
preserved unchanged outside this release workspace.
