# Bayes-Block-Aligner
The detailed theory of Bayesian Aligner can be found in<br>
<br>
   (1) Bayesian Adaptive Alignment and Inference<br>
       &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Jun Zhu, Jun S. Liu and Charles E. Lawrence<br>
       &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ISMB-97 pp.358-368.<br>
	   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;https://www.aaai.org/Papers/ISMB/1997/ISMB97-055.pdf<br>
   (2) Bayesian adaptive sequence alignment algorithm<br>
       &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Jun Zhu, Jun S. Liu and Charles E. Lawrence<br>
       &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Bioinformatics (in press) 2/98.<br>
	   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;https://academic.oup.com/bioinformatics/article/14/1/25/267263<br>

## Installation ##
Clone this repo. There is a 64-bit binary compiled under Ubuntu 20.04 LTS in the bin/ directory.<br>
To compile a new version cd to teh src/ directory and run: make. \<br>
<br>

## Options ##

<pre><code>
COMMAND LINE OPTIONS:
    -QUERY=: filename of sequence to align.
    -DATA=: The filename of the 2nd sequence to align.
    -NON_PROTEIN: switch between protein/DNA (Protein is the default).
    -MAX_NUM_BLOCKS=: The number of Sankoff blocks to align
    -DNA_MATRIX,STEP: DNA relation matrix number: [1-500] PAM,step (ie use a range).
    -MATRIX_NUM: protein relation matrix number:[0-7] BLOSUM30-100, [8] ALL BLOSUM, 
                  [9-35] PAM40-300, step 10,[36] ALL PAM.

     Advanced options: (These options should rarely be used, please send us email 
                         (steve@wadsworth.org) if you plan to use them)
         -DEBUG=: gives you detail information.
         -OUTFILE=: The output filename (Default is just the screen).
         -BEST_ALIGN: a straight forward Sankoff method.
         -BASE_BLOCKSIZE=: the base size of a block used in calculation of the
                           number of Sankoff blocks.(Advanced).
         -EXACT_BLOCK: Don't backsample proportional (Default is false).
         -EXACT_POST: Calculate the exact posterior probability (Default is false).
         -OUT_SEQ:output matching sequence if best alignments are calculated.
         -NO_PROFILE: Don't create a full(AKA big) .profile file .
         -FULL_PROFILE:  Create a full(AKA big) .profile file (Default is a sparse file) .
         -BACK_S_SIZE=: The number of backsamples.
         -BACK_S_OUTSEQ: Output all alignments
         -BACK_CUTOFF: find all alignments above a cutoff value.
         -BACK_C_VALUE=: Cutoff value
         -BACK_S_BLOCKS=: find the [N] most likely alignments (Default is 20)
</code></pre>

## Example ##
<pre><code>
bin/Bayesaligner -DATA=data/1GKY.fa -QUERY=data/2AK3A.fa  -MATRIX_NUM=1
</code></pre>

This code will produce 3 output files:<br>
&nbsp;&nbsp;&nbsp;gi1065285pdb-gi494043pdb.alignment.0 - lists the aligned blocks<br>
&nbsp;&nbsp;&nbsp;gi1065285pdb-gi494043pdb.marginal.0 - marginal probability of each of the query positions being aligned<br>
&nbsp;&nbsp;&nbsp;gi1065285pdb-gi494043pdb.profile.0 - base pair alkignment probabilities<br>

## Understanding the Results ##

The Bayesaligner is a "statistical" procedure that returns the entire<br>
alignment space (not just the maximum). The products of this program are<br>
the posterior probabilities of alignment variables and their visualizations.<br>
&nbsp;&nbsp;&nbsp;P(k|R)---posterior probability of k blocks given this pair of sequences.<br>
&nbsp;&nbsp;&nbsp;P(matrix|R)---posterior probability that this relation matrix correctly relates these sequences.<br>
<br>
The posterior probability of each pair of bases aligning can be obtained using
     backsampling or by exact calculation.  Backsampling is the default method
     because it uses less memory and less compute time . Backsampling will generate
     three files ```****.profile.*```, ```***.marginal.*```, and ```****.alignment.*```. The posterior 
     probability of each base aligning (```***.profile```) can be viewed using the supplied R code.
     The "peaks" represent the locations
     in each sequence where there is a strong alignment. The marginal posterior probability
     (```***.marginal.*```) of a particular base in the query sequence aligning with any base 
     in the data sequence can be view in R.
     The marginal plot is useful in finding the conserved regions of the query sequence
     in the data sequence(Phylogenetic footprinting). The alignment (```****.alignment.*```) file
     displays the most likely aligned blocks between the query and data sequences.<br>
	 
## Utilities  ##
The src/ directory contains two R files for plotting alignment probabilities.<br>

