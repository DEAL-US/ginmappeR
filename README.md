# ginmappeR

ginmappeR provides functionalities to translate gene or protein identifiers between state-of-art biological databases: [CARD](https://card.mcmaster.ca/), [NCBI](https://www.ncbi.nlm.nih.gov/) Protein, Nucleotide and Gene, [UniProt](https://www.uniprot.org/) and [KEGG](https://www.kegg.jp). Also offers complementary functionality like NCBI identical proteins or UniProt similar genes clusters retrieval.

## Launch as API

``` r
plumber::plumb_api("ginmappeR", "ginmappeR-API")$run(port=8000, swagger=FALSE)
```
