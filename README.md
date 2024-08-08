# ginmappeR

ginmappeR provides functionalities to translate gene or protein identifiers between state-of-art biological sequence databases: [CARD](https://card.mcmaster.ca/), [NCBI](https://www.ncbi.nlm.nih.gov/) Protein, Nucleotide and Gene, [UniProt](https://www.uniprot.org/) and [KEGG](https://www.kegg.jp). Also offers complementary functionality like NCBI identical proteins or UniProt similar genes clusters retrieval.

ginmappeR is available at [Bioconductor repository](https://bioconductor.org/packages/ginmappeR/).


### FAQ
---

- CARD Database manual download

If the user has limited access to its machine or CARD automatic download fails, it can be manually downloaded from the CARD website (https://card.mcmaster.ca/download) and extracted. The extracted files can be placed in any folder as with function changeCARDPath() the user can specify the path to them.