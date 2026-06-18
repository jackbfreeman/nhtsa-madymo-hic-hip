# MADYMO HIC₁₅ vs. real-world head injury (NASS-CDS)

Analysis code and manuscript source for a study validating MADYMO-simulated
Head Injury Criterion (HIC₁₅) scores against real-world head-injury outcomes in
frontal motor-vehicle crashes, using the U.S. NASS-CDS database (2000–2015) and
a complex-survey logistic-regression framework.

## Repository layout

```
.
├── r_docs/
│   ├── Final paper cleaned.Rmd   # the analysis + manuscript (knit this)
│   ├── references.bib            # bibliography
│   ├── elsevier-vancouver.csl    # citation style
│   ├── import.sty, title.sty     # LaTeX includes for the PDF build
│   └── Figure_1.png … Figure_3.png, roc_small.png   # figures used by the manuscript
├── data/
│   └── MADYMO_df.csv             # validated MADYMO HIC results (small; included)
├── nhtsa-madymo-hic-hip.Rproj    # RStudio project
├── LICENSE
└── README.md
```

## Reproducing the analysis

1. Open `nhtsa-madymo-hic-hip.Rproj` in RStudio (analysis run under **R 4.4.2**).
2. Obtain the large data files (see **Data availability** below) and place them
   under `data/`.
3. Open `r_docs/Final paper cleaned.Rmd`. Near the top, the **reproduction
   switches** control how much is rebuilt from scratch:

   | switch | default | effect when `TRUE` |
   |---|---|---|
   | `REBUILD_RAW` | `FALSE` | re-derive the cached data from the raw NASS-CDS + MADYMO files |
   | `REBUILD_MODELS` | `FALSE` | re-fit the survey models from the cached data |
   | `REBUILD_MICE` | `FALSE` | re-run the multiple-imputation sensitivity analysis (slow) |

   With all switches `FALSE` the document loads cached results and knits quickly.
4. Knit to PDF (uses `bookdown::pdf_document2`).

> Note: chunk `eval` is only honored when **knitting**. If you run chunks
> interactively in RStudio, the rebuild chunks will execute regardless of the
> switches.

## Data availability

The large inputs are not stored in this repository. They are openly available on
the Open Science Framework: **https://doi.org/10.17605/OSF.IO/VHK2P**

- raw NASS-CDS databases (`data/NHTSA_databases/`; also available directly from NHTSA)
- raw MADYMO simulation output (`data/minimal_simulations/`)
- cached extracts and fitted-model objects (`*.csv`, `*.rds`, `*.RData`)

Download the data bundle from OSF and unpack it into `data/` before knitting,
or set the relevant `REBUILD_*` switch to regenerate from the raw inputs.

## License

See [LICENSE](LICENSE).
