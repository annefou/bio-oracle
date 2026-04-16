# Nanopublication Plan — SST Fish Community Replication

Create these 7 nanopubs on [Science Live](https://platform.sciencelive4all.org/create) in order.

---

## Nanopub 1: Paper Quotation

**Template:** Annotate Quotation (`ANNOTATE_QUOTATION`)

- **Paper DOI:** https://doi.org/10.1111/gcb.16633
- **Quotation:** "In conclusion, our study clearly identifies sea temperature as the primary driver of fish community composition across the Northeast Atlantic continental shelf."
- **Comment:** Rutterford et al. analysed standardised abundance data for 198 marine fish species across 23 fisheries trawl surveys (31,502 sampling events, 2005-2018) from southern Spain to northern Norway. Their PCoA analysis showed SST was the environmental variable most closely associated with the first axis of variation, explaining 29% of variance (R²=0.890). Bio-ORACLE (Assis et al., 2018) was used for environmental layers.

---

## Nanopub 2: AIDA Sentence

**Template:** AIDA Sentence (`AIDA_SENTENCE`)

- **Sentence:** According to Rutterford et al. (2023), analysis of standardised abundance data for 198 fish species across the Northeast Atlantic (2005-2018) identifies sea surface temperature as the primary driver of fish community structure, explaining 29% of community variation (R²=0.890), followed by salinity (21%) and depth (12%).
- **Paper DOI:** https://doi.org/10.1111/gcb.16633

---

## Nanopub 3: FORRT Replication Claim

**Template:** FORRT Claim (`FORRT_CLAIM`)

- **Original study DOI:** https://doi.org/10.1111/gcb.16633
- **Claim:** Sea surface temperature is the primary driver of fish community structure across Northeast Atlantic shelf seas.
- **Domain:** Marine ecology / Fish community ecology
- **Evidence type:** Observational (standardised trawl survey analysis + PCoA + regression)

---

## Nanopub 4: FORRT Replication Study — Reproduction (DATRAS)

**Template:** FORRT KL Replication (`FORRT_KL_REPLICATION`)

- **Original study DOI:** https://doi.org/10.1111/gcb.16633
- **Replication type:** Reproduction (same data source, same method)
- **Replication study author:** Anne Fouilloux (ORCID: 0000-0002-1784-2920)
- **Description:** Reproduction using the same ICES DATRAS trawl survey data (22 NE Atlantic surveys, 2005-2018, 258 fish species, 247 grid cells). Filters applied: shelf depth ≤200m, fish species only (WoRMS), Baltic surveys excluded.
- **Data sources:**
  - ICES DATRAS: https://datras.ices.dk/
  - Bio-ORACLE: https://bio-oracle.org/
- **Code repository:** https://github.com/annefou/bio-oracle
- **Analysis type:** PCoA on Bray-Curtis dissimilarity + regression

---

## Nanopub 5: FORRT Replication Outcome — Reproduction (DATRAS)

**Template:** FORRT KL Replication Outcome (`FORRT_KL_REPLICATION_OUTCOME`)

- **Links to:** Nanopub 4 (the replication study)
- **Outcome:** Confirmed
- **Direction:** Consistent with original claim
- **Summary:** PCoA analysis of 258 fish species across 247 NE Atlantic grid cells (DATRAS, 2005-2018) confirms SST as the primary driver of community structure. SST was most strongly associated with PCoA Axis 1 (29.2% variance, r=-0.757, R²=0.573, p<0.001). Full model R²=0.879. Variance explained (29.2%) closely matches Rutterford et al.'s 29%.
- **Scope note:** Reproduction uses the same data source (ICES DATRAS) and methodology. Slightly more species (258 vs 198) and grid cells (247 vs 193) due to different filtering thresholds.

---

## Nanopub 6: FORRT Replication Study — Replication (ClimateFish + MEDITS)

**Template:** FORRT KL Replication (`FORRT_KL_REPLICATION`)

- **Original study DOI:** https://doi.org/10.1111/gcb.16633
- **Replication type:** Replication (different region, different data sources)
- **Replication study author:** Anne Fouilloux (ORCID: 0000-0002-1784-2920)
- **Description:** Independent replication in the Mediterranean Sea using two datasets: (1) ClimateFish visual census (15 species, 28 sites, 2009-2021), and (2) MEDITS demersal trawl surveys filtered to shelf depth ≤200m (230 species, 115 grid cells, 2005-2024). Same analytical approach: PCoA on Bray-Curtis dissimilarity + regression against SST, salinity, depth.
- **Data sources:**
  - ClimateFish: https://www.seanoe.org/data/00756/86784/ (CC-BY 4.0)
  - JRC MEDITS: https://data.jrc.ec.europa.eu/dataset/ef36af5d-eb4e-4a9f-9513-82a877fe71fe
  - Bio-ORACLE: https://bio-oracle.org/
- **Code repository:** https://github.com/annefou/bio-oracle
- **Analysis type:** PCoA on Bray-Curtis dissimilarity + regression

---

## Nanopub 7: FORRT Replication Outcome — Replication (ClimateFish + MEDITS)

**Template:** FORRT KL Replication Outcome (`FORRT_KL_REPLICATION_OUTCOME`)

- **Links to:** Nanopub 6 (the replication study)
- **Outcome:** Confirmed
- **Direction:** Consistent with original claim
- **Summary:** Both Mediterranean datasets confirm SST as the primary driver of shelf fish community structure. ClimateFish: SST most strongly associated with PCoA Axis 1 (45.2% variance, r=-0.784, R²=0.615, p<0.001). MEDITS shelf: SST most strongly associated with PCoA Axis 1 (18.9% variance, r=-0.605, R²=0.367, p<0.001). The weaker effect in MEDITS reflects the narrower Mediterranean SST range (5.6°C vs 10°C in NE Atlantic). Without the shelf depth filter, depth dominates — confirming the claim is specific to continental shelf communities.
- **Scope note:** Replication uses different region (Mediterranean vs NE Atlantic), different data sources (visual census + demersal trawl vs standardised trawl surveys), and different species sets (15+230 vs 198). Both confirm the same conclusion.

---

## Creation Order

1. Paper Quotation → get URI
2. AIDA Sentence → get URI
3. FORRT Claim → get URI (references nanopub 1)
4. FORRT Replication Study (DATRAS) → get URI (references nanopub 3)
5. FORRT Replication Outcome (DATRAS) → get URI (references nanopub 4)
6. FORRT Replication Study (Med) → get URI (references nanopub 3)
7. FORRT Replication Outcome (Med) → get URI (references nanopub 6)

## Embedding in Jupyter Book

After creating each nanopub, note the URI and embed with:

```markdown
<iframe src="https://platform.sciencelive4all.org/embed/view?uri=NANOPUB_URI&showShare=false" 
        width="100%" height="500" frameborder="0"></iframe>
```
