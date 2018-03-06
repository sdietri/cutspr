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

<b>How to start CutSPR</b><br>
If you use the python version from GitHub, ensure that the dependencies are resolved and run cutspr.py:<br>
$ python3 cutspr.py

Alternatively, prepackaged portable versions are available for the 64bit versions of Windows, MacOSX and Linux at
http://appmibio.uni-goettingen.de/index.php?sec=sw

These versions only require the NCBI blast suite to be installed on your system and should require no additional setup.

CutSPR has been tested successfully on Ubuntu 16.04 LTS, Linux Mint 18.3, Windows 7 and 10 and MacOSX El Capitan, Sierra and High Sierra.
Biolinux8 14.04 was found to be incompatible.

<b>Feedback and Support</b><br>
We aimed to make this application as compatible as possible to the commonly used operating systems. However it is impossible to cover everything and foresee all potential obstacles. So if you have faced difficulties using CutSPR, do not hesitate to contact us so we can find a solution.

Mail to:<br>
Robert Hertel: rhertel(at)gwdg.de<br>
Sascha Dietrich: sdietri(at)gwdg.de

