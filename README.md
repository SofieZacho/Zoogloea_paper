# *Zoogloea* study

**Title**: Genome-resolved taxonomy, global distribution, and EPS potential of *Zoogloea*, the canonical floc-former in wastewater treatment systems

**Authors**: Sofie Zacho Vestergaard, Lei Liu, Morten Kam Dahl Dueholm, Per Halkjær Nielsen

Center for Microbial Communities, Department of Chemistry and Bioscience, Aalborg University, Aalborg, Denmark.


## Intro to study

*Zoogloea* is a globally abundant wastewater bacterium involved in floc formation and nutrient removal, yet its taxonomy and floc-forming mechanisms remain poorly resolved. Using 47 *Zoogloea* genomes from the MiDAS global genome catalogue and 3 DSMZ isolate recovered *Zoogloea* genomes together with global abundance data, we clarified the genus’ taxonomy, distribution, and metabolic diversity. We further investigated a recently described floc-forming polysaccharide biosynthesis gene cluster and found it distributed across multiple genera within the Burkholderiales. The cluster was highly conserved and consistently co-occurred with floc-regulating genes, suggesting a widespread and shared genetic basis for floc formation in activated sludge ecosystems.


## Repo structure


| Folder | Content |
| --- | --- |
| ├── 00_Figures/ | The code used to create figures in the manuscript named after which figure they. |
| ├── README.md | A brief overview of the project and workflow |
| └── LICENSE | The [license](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/licensing-a-repository) for the repo. |


## Overview of main analyses in study
* Cultivated three *Zoogloea* strains, generated complete genome assemblies using long-read sequencing, and curated assemblies by polishing and contaminant removal.
* Integrated newly generated genomes with publicly available *Zoogloea* genomes and the MiDAS global genome catalogue for comparative analyses.
* Assessed the global distribution and relative abundance of *Zoogloea* species using 16S rRNA amplicon data across 480 activated sludge wastewater treatment plants worldwide.
* Performed genome annotation and metabolic reconstruction to investigate functional potential relevant to activated sludge nutrient cycling.
* Reconstructed genome-based and 16S rRNA phylogenies and calculated average nucleotide identity (ANI) to resolve taxonomic and evolutionary relationships within *Zoogloea*.
* Identified, curated, and compared EPS (extracellular polymeric substance) biosynthesis gene clusters across *Zoogloea* genomes.
* Built HMM-based screening pipelines to map the phylogenetic distribution of the *Zoogloea* EPS gene cluster across thousands of genome catalogues.

## Data availability
The raw and assembled sequencing data of the 3 DSMZ *Zoogloea* isolates have been deposited in the European Nucleotide Archive (ENA) at EMBL-EBI under accession number `PRJEB115673`.
16S rRNA gene amplicon data are available through the [MiDAS 4 project](https://www.nature.com/articles/s41467-022-29438-7). MAGs from the MiDAS global genome catalogue are deposited in ENA under BioProject `PRJEB83983`.
