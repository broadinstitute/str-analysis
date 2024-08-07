This folder contains catalogs of known disease-associated TR loci in a format that can be used with [ExpansionHunter](https://github.com/Illumina/ExpansionHunter) or [TRGT](https://github.com/PacificBiosciences/trgt).
The catalogs also include 9 loci (TMEM185A, DMD, etc.) that are not currently considered causal for any Mendelian disease, but that were historically listed as candidate disease loci and so were often included in publications that cataloged disease-associated STRs.

**TRGT catalogs:**

There are two versions of the TRGT catalog: one with, and one without adjacent repeats (such as the CCG repeat next to the main CAG repeat at the HTT locus). 


**ExpansionHunter catalogs:**

There are four versions of the ExpansionHunter catalog: for GRCh38 and GRCh37, with and without off-target regions.

[ExpansionHunter docs](https://github.com/Illumina/ExpansionHunter/blob/master/docs/04_VariantCatalogFiles.md) describe the `LocusId`, `LocusStructure`, `ReferenceRegion`, `VariantType`, `VariantId` and `OfftargetRegions` fields in more detail. 
These fields modify ExpansionHunter behavior for each locus. 

ExpansionHunter ignores extra fields in the locus definition, so I added reference information to the catalog in the fields listed below:

* `Gene` - the name of the gene that contains the STR locus  
* `GeneId` - the Ensembl ESNG gene id
* `GeneRegion` - what part of the gene contains the STR locus 
* `DiscoveryYear` - year when the disease association was first established [source: [Depienne 2021](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8205997/)]
* `DiscoveryMethod` - how the disease association was originally discovered [source: [Depienne 2021](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8205997/)]. 
    * `cl` is `cloning` 
    * `CGA` is `candidate gene analysis`
    * `ExScr` is `expansion screening`
    * `L` is `linkage`
    * `RED` is `repeat expansion detection`
* `Diseases` - list of 1 or more diseases associated with this locus
  * `Symbol` - disease symbol (ie. `DM1`)
  * `Name` - disease name (ie. `Myotonic dystrophy 1`)
  * `Inheritance` - the inheritance mode(s) of the disease(s) associated with this locus
    * `AR` is `autosomal recessive`
    * `AD` is `autosomal dominant`
    * `XR` is `x-linked recessive`
    * `XD` is `x-linked dominant`
  * `OMIM` - OMIM phenotype MIM number of the disease (ie. [`160900`](https://omim.org/entry/160900?search=160900&highlight=160900))
  * `NormalMax` - the maximum number of repeats at this locus that do not lead to a disease phenotype in the individual harboring this number of repeats
  * `PathogenicMin` - the smallest number of repeats that lead to fully-penetrant disease. This threshold is based on sources such as [[Depienne 2021](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8205997/)], [STRipy](https://stripy.org/database), [OMIM](https://www.omim.org/), and [GeneReviews](https://www.ncbi.nlm.nih.gov/books/NBK1116/). These thresholds should be considered approximate and subject to revision. 
  * `IntermediateRange` - some of the disease-associated STR loci also have an intermediate range below the pathogenic threshold that is associated with milder disease or reduced penetrance. These ranges should be considered approximate and subject to revision. 
  * `Note` - any additional notes about this disease association
* `RepeatUnit` - the main repeat motif 
* `PathogenicMotif` - for loci such as RFC1 or DAB1 at which the reference genome motif is not the same as the pathogenic motif, this field reports the subset of  motif(s) that are known to be pathogenic when expanded
* `MainReferenceRegion` - some loci have adjacent repeat regions specified to improve genotyping accuracy. For example, the FXN locus definition includes an adjacent poly-A repeat. The MainReferenceRegion field provides the reference coordinates of the main disease-associated region. 

These variant catalogs were used to generate the [gnomAD STR callset](https://gnomad.broadinstitute.org/short-tandem-repeats?dataset=gnomad_r3)
described in more detail in this [blog post](https://gnomad.broadinstitute.org/news/2022-01-the-addition-of-short-tandem-repeat-calls-to-gnomad/).
