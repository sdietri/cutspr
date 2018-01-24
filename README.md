# CutSPR
CutSPR is a tool that assists the construction of a mutagenesis vector based on the pJoe8999 plasmid system.

It uses genomic references in the form of GenBank or Fasta formatted sequence files to generate primers for the construction of a specific CRISPR-Cas9 based deletion or insertion vector.


CutSPR is written for Python3 and requires the following additional packages:
<ul>
<li>PyQT5</li>
<li>biopython</li>
</ul>

In addition, the NCBI Blast software package needs to be installed for CutSPR to function which can be obtained here:
https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download


On Ubuntu 16.10, the following steps should resolve all necessary dependencies:
<ul>
  <li>sudo apt-get update</li>
  <li>sudo apt-get install python3.6</li>
  <li>sudo apt-get install python3-pyqt5</li>
  <li>sudo apt-get install python3-biopython</li>
  <li>sudo apt-get install ncbi-blast+</li>
</ul>

If you’re using another version of Ubuntu (e.g. the latest LTS release), we recommend using the deadsnakes PPA to install Python 3.6:
<ul>
  <li>sudo apt-get install software-properties-common</li>
  <li>sudo add-apt-repository ppa:deadsnakes/ppa</li>
  <li>sudo apt-get update</li>
</ul>

Start CutSPR with the flowing command line:<br>
$ python3 cutspr.py
  



