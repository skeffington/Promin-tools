## Promin Tools

This is the project page for Promin Tools: a collection of bioinformatic tools for the analysis of protein seqeunces thought to be involved in biomineralization. 

So far there are two tools:

1. **Protein Motif Finder** which finds motifs overrepresented in one set of protein sequences compared to another.
2. **Sequence Properties Analyzer** which performs statistical analyses of one set of proteins compared to another.

These tools are stil being updated, so if you want new features or have suggestions for improvements please get in touch.

##### Table of Contents  

- [Quick start guidevia the discovery environment](#quick-start-guide-via-the-discovery-environment)
- [Quick start guide via the docker image](#quick-start-guide-via-the-docker-image)
- [Protein Motif Finder](#protein-motif-finder)
  * [Description of inputs](#description-of-inputs)
  * [Description of Outputs](#description-of-outputs)
  * [Description of the program](#description-of-the-program)
  * [The motif-x software](#the-motif-x-software)
- [Sequence properties analyzer](#sequence-properties-analyzer)
  * [Description of inputs](#description-of-inputs)
  * [Description of Outputs](#description-of-outputs)
- [Contact](#contact)

### Quick start guide via the discovery environment

The tools are hosted on Cyverse: and online cyberinfrastructure funded by the National Science Foundation’s Directorate for Biological Sciences.

* Go to https://cyverse.org/
* Click **Launch** and then **Discovery Environment**
* Set yourself up with a Cyverse account
* Login
* Open the **Data** window
* Upload your data as two fasta files
* Open the **Apps** window
* Search for **Protein Motif Finder** and select
* Select your data, choose the parameters an **Launch Analysis**
* You will get a notification message when your analysis is complete

### Quick start guide via the docker image

If you know your way round Docker, the Protein Motif Finder Docker image is available here:

https://hub.docker.com/repository/docker/biologistatsea/promofi/

run with the command:

```
perl -S motif_finder.pl [fasta input of foreground sequences (> 3 sequences)] [fast input of background sequnces] [window size (must be an odd number)] [p-value cutoff] [output prefix]
```
For example:
```
perl -S motif_finder.pl foreground.fasta background.fasta 9 1e-6 my_output
```

For further details about the inputs and outputs see below.

The Sequence Properties Analyzer is coming soon.

```
perl -S seq_properties.pl [foreground fasta] [background fasta] [flps p-value (eg 1e-4)] [seg window size (eg 12)] [seg locut (eg 2.2)] [seg hicut (eg 2.5)] [output prefix]

```

### Protein Motif Finder

#### Description of inputs:

*Fasta inputs:* Foreground and background sequences should be in fasta format and protein sequences should be represented with uppercase letters using the single-letter amino acid code. There should not be a * or any other non amino-acid character at the end of, or elsewhere in, the sequence. 

*Foreground sequences:* These are the protein sequences of interest. There should be a minimum of three sequences. If there are less than three seqeunces the program will abort. They may be entire protein sequences or partial protein sequnces (eg a particular domain extracted from a number of proteins). This set of proteins does not have to conatin every amino acid. 

*Background sequences:* This is the set of sequences to which the foreground seqeunces are compared when identifying statistically overrepresented motifs. Typically this is the pridicted proteome for your organism of interest, however it could also be another set of proteins depending on the biological question which you are asking. The background set may include the foreground proteins or exclude them. In the former case, a motif enrichment factor of one means the motif is only presnt in the foreground sequences and an enrichment factor less than one is not possible. In the latter place motifs that are present in the foreground set but abscent in the background set are reported as infinitely enriched and enrichment factors between zero and one are possible. Enrichment factors less than one mean that the motif is found more often in the foreground sequences than in the background sequences. It is assumed in the program that the background sequence set is large enough that all amino acids are represented.

*Window size:* This is effectively the maximum length of a motif that can be found by the program. Motifs shorter than the the window size can be found. It must be an odd number. The larger the window size, the longer it will take the program to run. In practice, window sizes of 7, 9, or 11 are often a good starting point. 

*p-value cutoff:* A value of 1e-6 is a good starting point. The family-wise error rate represented by the threshold will vary with the window size. The following information is taken from  Wagih et al. 2011:

"The significance parameter corresponds to the binomial probability threshold necessary to “fix” each motif position during the motif-building phase of the algorithm. It is critical to note that this value does not take into account a correction for multiple hypotheses (such as the Bonferroni correction). On any given motif-x search step there are (number of possible characters at each position) * (number of nonfixed positions) hypotheses being tested. For example, in an S-centered analysis of width 15, there would be (20) * (14) = 80 hypotheses tested. To ensure an alpha-value of at least 0.05 by the Bonferroni method, one would need to divide the desired alphavalue by the total number of hypotheses tested (i.e., 0.05/280 = 0.00018). Thus, for the previous example inputting a significance value of 0.00018 into motif-x would in fact correspond to a p-value of 0.05. The use of a motif-x significance threshold 15 greater than 0.0005 is not suggested as it may result in the extraction of motifs that are not statistically significant. The default value of this parameter is 0.000001, which corresponds to an actual alpha-value of approximately 0.0003 for a protein motif analysis of width 15 after Bonferroni correction. A significance threshold of 0.000001 was chosen for the present example."

*Output prefix:* This is the text which will be at the beginning of each of the output file names 

#### Description of Outputs:

* _bgmotifs.txt : counts of each of the overrepresented motifs in each of the proteins in the background sequence set
* _fgmotifs.txt : counts of each of the overrepresented motifs in each of the proteins in the foreground sequence set
* _fgenrich.txt : the enrichment of each of the overrepresented motifs in each of the proteins in the foreground sequence set with repspect to the background sequence set  
*_motifsummary.txt : for each overrepresented motif, the enrichemnt in the forgraound sequnces relative to the background; the count of proteins in which it appears and the median count of the motif per  protein  
* _motifs.html : A html summary of the analysis with some helpful plots  
* _Rscript.R : The R script that generated the html output. This can be used as a basis to produce your own custom figures and plots in R.  
* _wordclouds.svg : An svg file of the wordclouds from the html report. This can be imported into inkscape or illustrator to make publication ready figures.  

These filenames will be prefixed by the user specified file prefix.

#### Description of the program:

The motif-x motif finding engine accepts aligned input sequences centred on a particular residue. Thus the first part of the program processes the fasta inputs and generates the input files necessary to run motif-x. Motif-x is then run via the R module, rmotifx, centered on each of the amino acids in turn and the results are combined. Some of the significant motifs identified are likely to be redundant. For example, if the central residue was K and the window size was 7, the motifs ...K..K and K..K... may both be identified separately. The program removes any leading or trailing wildcards from the motifs and collapses this redundancy, so that in this example only K..K would be reported in the final list of motifs. It should be noted that this procedure is conservative with respect to the p-values computed by motif-x (ie all motifs in the final list will also meet the p-value threshold defined by the user). Next, all the found motifs are enumerated in the forground and the background data sets, enrichments calculated and the output tables printed. An R script then reads the ouput tables and produces a grphical output with explanation in the html document.


#### The motif-x software:

Wagih O, Sugiyama N, Ishihama Y, Beltrao P. (2015) Uncovering phosphorylation-based specificities through functional interaction networks (2015). Mol. Cell. Proteomics 
Chou MF & Schwartz D (2011). Biological sequence motif discovery using motif-x. Curr Protoc Bioinformatics. Chapter 13:Unit 13.15-24. doi:10.1002/0471250953.bi1315s35. 

The motif-x software may be used freely by users from academic and non-profit organizations. Users from the commercial sector should contact Daniel Schwartz (daniel.schwartz(at)uconn.edu).

The software in this package is provided under the GPL-3 licence.


### Sequence properties analyzer

#### Description of inputs

As with Protein Motif Finder, the main inputs are:

*Foreground fasta input*

*Background fasta input*

The other parameters are for the Seg program and should be left as defaut except for advanced users.

#### Description of Outputs

The main outputs of the program are the data tables. For convenience these are processed to an html document with
some visualizations of the data. However one analysis pipeline can never be appropriate for every dataset,
so the R script used to generate these diagrams is also an output of the program. This is to allow you, or an
infromatically minded colleague to easily tweak the plots or to develop these analyses further in a way 
appropriate for your data. 

The main tabular outputs of the program are listed below:  

* _AAabundance.txt: The data displayed in table 2.
* _AAenrich.txt: Enrichment (fold change) for each amino acid in the POI with respect to the background (displayed in table 1). The second numeric column is the data for figure 1.
* _bg.comp: Relative frequencies of amino acids in the background protein set
* _fg.comp: Relative frequencies of amino acids in the POI protein set
* _chargedclus_POI_summary.txt: For each protein in the POI set, the number, length and percentage of the sequence covered by positive and negative clusters of amino acids. Data used in figure 5.
* _chargedclus_POI.txt: The list of POI proteins containing positive or negative clusters of amino acids, along with the positions of the clusters.
* _chargedclus_Proteome_summary.txt: For each protein in the background proteome set, the number, length and percentage of the sequence covered by positive and negative clusters of amino acids. Data used in figure 5.
* _chargedclus_Proteome.txt: The list of background proteome proteins containing positive or negative clusters of amino acids, along with the positions of the clusters.
* _chargedclus_seqtbl.txt: The data displayed in table 4.
* _complexityPOI_summary.txt: The data displayed in table 3. Data used in figure 3.
* _complexityProteome_summary.txt: Similar to the data displayed in table 3, but for the background proteome protein set. Data used in figure 3.
* _POI_disorder_summary.txt: For each POI protein, gives the length and number of disordered regions and the percentage of the sequence predicted to be disordered.
* _POI_disorder.txt: Gives the position in the POI proteins of each disordered region. Data used in figure 4.
* _Proteome_disorder_summary.txt: For each background proteome protein, gives the length and number of disordered regions and the percentage of the sequence predicted to be disordered.
* _Proteome_disorder.txt: Gives the position in the background proteome proteins of each disordered region. Data used in figure 4.
* _positions.txt: Gives the positions of all compositionally biased regions in the POI protein set, including the residues that are biased and the inverse of the binomial p-value for the bias.
* _protcount.txt: For each type of compositional bias (amino acid or group of amino acids), gives the count of proteins in the POI set containing this bias. Data used in figure 2.
* _Rscript.R : The R script that generated the html output. This can be used as a basis to produce your own custom figures and plots in R  


### Contact: 
This software was written by Alastair Skeffington of the Max Planck Institute for Molecular Plant Physiology, Potsdam, Germany
