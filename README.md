<h1> <i>Quercus havardii - Forests</i> Manuscript </h1>

Raw data and R scripts for analyzing <i>Quercus havardii</i> genetic, morphological, and environmental analyses. 
This manuscript is in preparation to be submitted to the open access journal <i>Forests</i> in a Special Issue titled “<i>Quercus</i> Genetics: Insights into the Past, Present, and Future of Oaks.

<h2>Overview</h2>
The genetic diversity patterns of widespread species can be influenced by large scale and local factors. Environmental heterogeneity, as well as degree of connectivity, across the distribution of a species can play a critical role in the genetic and morphological divergence of populations and can act as a speciation mechanism. We examine these relationships within the genus <i>Quercus</i>, a model clade for exploring evolutionary processes due to its diversity, dominance, and importance in ecosystems. We studied <i>Quercus havardii</i>, a uniquely adapted desert oak with a disjunct distribution. We used population-level analyses of genetic, morphological, and environmental variables over various geographic scales to quantify differentiation and understand the forces that influence the divergence of populations in the species.


<h2>Development</h2>

We welcome contributions from any individual, whether code, documentation, or issue tracking.
All participants are expected to follow the code of conduct for this project.


<h3>Authors and Contributions:</h3>

<b>Bethany A. Zumwalde (BAZ) <br>
Ross A. McCauley (RM) <br>
Drew Duckett (DD) <br>
Ian J. Fullinwider (IF) <br>
Emma Suzuki Spence (ESS) <br>
Sean Hoban (SH)</b> <br>

BAZ carried out lab work, analyzed genetic and environmental data, and wrote the paper.  <br>
RM carried out fieldwork, collected and analyzed morphological data, and assisted with writing the paper.  <br>
DD carried out fieldwork, labwork, analyses and assisted with writing the manuscript.  <br>
IF collected morphological data and reviewed the manuscript.  <br>
ESS assisted with lab work, analyses, and assisted with writing the manuscript.  <br>
SH conceived the project, carried out fieldwork, supervised lab work, assisted with genetic and environmental analyses, and co‐authored the manuscript. <br>


<h2>Data Collection</h2>
Seed, leaves, soil, and herbarium vouchers were sampled from havardii in two separate trips referred to as East and West. The East sampling was performed by Sean and Drew in Texas, New Mexico, and Oklahoma from 8/12/16-8/16/16. Samples were also sent from other researchers and institutions. The West sampling was performed by Sean and Ross McCauley (Fort Lewis College) from 8/25/16-8/30/16. Samples were also sent from other researchers and institutions.
<br>
A total of 667 samples from 26 primary populations and 10 auxiliary populations of georeferenced locations were used in this study.   Populations were chosen by contacting land managers of private and public land in the region, by consulting GBIF and SEINet, and via suggestions from the International Oak Society.  The objective was to ensure populations were sampled throughout the geographic range. 

<h2> Delimiting Populations</h2>
It was decided that populations must be roughly 5 km (~3 miles) apart from each other.  This distance was based off of pollen flow limits.  Emma Spence used GIS to determine population edge proximities to each other and total area of each population which can be found here. For populations with individual GPS coordinates also had calculations for individual-individual distances (mean, min, etc.) within populations. Debates were had of which populations may or may not need to be split or recombined. To summarize the issue briefly, most of our populations are sampled on the scale of a few hundred meters to about 1 km among individuals sampled per population.



<h2>Available Data and Description of Workflow</h2>

Below is a list of main items and short descriptions for separate analyses used in the paper

<ol>
<b><li>Genetic Dataset</li></b>
  <ol>
    <b><li>Qhavardii</li></b>
      <ol>
         <b> <li>7indiv_11loci_RedefinedPopsMay2019</b> (Primary folder containing raw files and outputs used to calculate genetic statistics)</li>
         <b> <li>code_havardii_analysis_final.R </b>(Final script used to execute analyses) </li>
         <b> <li>Older </b>(This contains raw data and code for preliminary analyses)</li>
         </ol>
  </ol>
  </li>
<b><li>Morphometric Dataset</li></b>
  <ol>
      <b><li>Raw Data</li></b>
      <ol>
        <b> <li>QH_total_ANOVA.csv </b>(Used for ANOVAs)</li>
        <b> <li>QH_total.csv </b>(Used for PCA analyses and boxplots)</li>
         </ol>
     <b><li>Scripts</b></li>
         <ol>
        <b> <li>Final harvardii morphology scripts.R</b> (R script for all morphological analyses)</li>
         </ol>
  </ol>
  </li>
<b><li>Environmental Dataset</li></b>
  <ol>
    <b><li>Raw Data</b></li>
      <ol>
      <b>  <li>QH_Pops_Final.csv</b> (Geographic localities used for bioclimatic analyses)</li>
        </ol>
    <b><li>Scripts</b></li>
         <ol>
         <b><li>QH_Environ_Final.R</b> (Used to make environmental PCA)</li>
        <b> <li>QH_Environ_Normality_Correlations.R</b> (Used to check for correlated bioclimatic variables)</li>
         </ol>
  </ol>
  </li>
<b><li>Regressions</li></b>
  <ol>
    <b> <li>Hoban_work</b> (Contains all raw code, scrips, and outputs for regressions of genetic and morphological traits with environmental variables)</li>
  </ol>
  </li>
<b><li>Figures</li></b>
 <ol>
    <b><li>PCA Panel Figure</b> (Scripts for making panel of genetic, morphological, and environmental PCA's) </li>
    <b><li>Preliminary Figures</b> (Archival preliminary figures) </li>
         </ol>
</ol>
