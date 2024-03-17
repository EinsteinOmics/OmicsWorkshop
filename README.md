# Intro to Omics Workshop Vignettes

An asynchronous learning workshop introducing key topics in Omics. For schedule, location, and other details, see our website: https://einsteinomics.github.io/omics-workshop2024/



<br />

## General instructions

The Omics Workshop sessions will feature code examples with real data, which we call "vignettes". Each vignette will be available in three ways:

1. Jupyter Notebook files (.ipynb), for those who want to follow along interactively using the HPC Jupyter workflow, where all software has been pre-configured
2. local RMarkdown or Jupyter files, for those who want to follow along interactively using R or Python on their own computers ("locally"). Note that you must install relevant software packages for each vignette
3. View-only HTML files, for those who want to follow along on their computer without any interactive coding (scroll to bottom for this)


<br />
<br />

## About this guide

This guide will walk you through:
- Getting yourself onto the HPC and https://cluster.einsteinmed.edu/
- Getting the code vignettes onto your account on HPC
- Using pre-installed software packages on cluster.einsteinmed.edu
- At the bottom are instructions specific to each vignette, including non-HPC options


<br />
<br />


## Getting yourself on the HPC

To follow along on Einstein's High Performance Computer (HPC), it is reqired that you sign up for an HPC account, and access https://cluster.einsteinmed.edu/. Optionally, is also recommended to familiarize yourself with the remote access option (GlobalProtect) for HPC use outside of the EinsteinSecure wifi.

- [HPC account setup](https://montefioreorg.sharepoint.com/sites/Einstein-IT-HPC/SitePages/HPC3.0-UQuick-Start.aspx) - login required. It has all the relevant guides including how to set up HPC account.
- [cluster.einsteinmed.edu usage](https://montefioreorg.sharepoint.com/sites/Einstein-IT-HPC/Shared%20Documents/Forms/AllItems.aspx?id=%2Fsites%2FEinstein%2DIT%2DHPC%2FShared%20Documents%2FGeneral%2FHPC3%2E0%20docs%2FIntroduction%20to%20Einstein%20HPC%20Portal%2Dv5%2Epdf&parent=%2Fsites%2FEinstein%2DIT%2DHPC%2FShared%20Documents%2FGeneral%2FHPC3%2E0%20docs) - try to follow this and activate the Jupyter Notebook workflow. Note that it may take a few minutes to start running.
- [A guide for opening Jupyter Notebooks on cluster.einsteinmed.edu](https://drive.google.com/file/d/1yZ9yI0QSOlWYwSjmssMvgqt0up1-9wa1/view?usp=sharing), prepared by PhD Candidate and Omics Club board member Marliette Rodriguez Matos.
- [Remote access to HPC via GlobalProtect](https://montefioreorg.sharepoint.com/sites/Einstein-IT-HPC/Shared%20Documents/Forms/AllItems.aspx?id=%2Fsites%2FEinstein%2DIT%2DHPC%2FShared%20Documents%2FGeneral%2FHPC3%2E0%20docs%2FIT%2DREF%2D2023%2D094%20Einstein%20Academic%20Research%20Systems%20Portal%282%29%2Epdf&parent=%2Fsites%2FEinstein%2DIT%2DHPC%2FShared%20Documents%2FGeneral%2FHPC3%2E0%20docs)

If you run into errors, you can reach out to us or submit an IT ticket.

<br />
<br />

## Getting the code vignettes onto your account on HPC

For HPC usage: first, try to activate the Jupyter workflow on https://cluster.einsteinmed.edu/. This will create a folder in your HPC home directory called "pw". It is recommended to store all vignette code in there to be easily accessible to Jupyter.


Then there are two options for downloading the vignettes to your account: a terminal based option with Git, and a non-terminal option of folder mounting similar to data.einsteinmed.edu.

### 1. Unix shell (bash / zsh / terminal) with `git clone`. 

- First, open a Unix shell, such as Mac terminal. On Windows, there are several options including [Git Bash](https://git-scm.com/download/win), [Cygwin](https://www.cygwin.com/), and [WSL](https://learn.microsoft.com/en-us/windows/wsl/about). 
- Next, use `ssh` to get to the HPC. Follow the instructions above, and make sure to use GlobalProtect if off campus. In short, you must use `ssh username@username.hpc.einsteinmed.edu` (replace "username" with your active domain).
- Once there, use the "cd" command to get to the pw folder in your home directory: `cd ~/pw/`.
- Finally, you can download the scripts: `git clone https://github.com/EinsteinOmics/OmicsWorkshop.git`
- Note that later vignettes will be added. To update the folder once those are added, you can cd into the folder and use `git pull`.
- If `git clone` does not work when updating the folders, you can also try `git fetch --all` followed by `git reset --hard origin/main`

<br />

### 2. Folder Mounting to HPC just like data.einsteinmed.edu

This makes use of the folder mounting, identical to how data.einsteinmed.edu is used. (HPC and data.einsteinmed.edu have access to the same files)
- [Folder mounting guide](https://it.einsteinmed.edu/documentation/how-to-mount-the-hpc-file-system/) - Note that if "einsteinmed.**org**" as the guide says does not work, you can try "einsteinmed.**edu**" instead.

First, download the vigenttes from Github. This can be done by clicking the big green button that says "Code" at the top of the screen. Then click "Download Zip". Make sure to unzip the file. Then, use the guide above to move the folder to HPC. You should place it in the "pw" folder in your home directory.

Note that later vignettes will be added. To download only a certain folder, you can use the steps above, unzip, and transfer only the folder(s) not yet uploaded to your account. You can repeat this once the later vignettes are released.

<br />
<br />



## Using pre-installed software packages on cluster.einsteinmed.edu

Once you have logged on to https://cluster.einsteinmed.edu/, and opened the relevant vignette files ending in .ipynb, you should load the software packages which have been pre-installed before this workshop. These are implemented as Conda environments configured as Jupyter kernels. This means it is basically one click to use the right software.

When you first open the Jupyter notebook, you may be prompted to select a kernel. Alternatively you can always switch kernal at the top right of the notebook. You should select either "omics-workshop-R" or "omics-workshop-python" based on the vignette.

Use "omics-workshop-R" for the following vignettes:
- 03_BulkRNAATAC_KithXiang
- 05_scATAC_LindsayAnthony
- 06_scRNA_DeyouAlex

Use "omics-workshop-python" for the following vignettes:
- 04_AncestryAnalysis_SriDavid
- 07_Proteomics_SimoneJuan


<br />
<br />

# Further Guides for Specific Vignettes


### 03 Dr. Kith Pradhan and Xiang Yu Zheng → bulk RNA/ATAC seq analysis

Dr. Pradhan has made available a Dropbox link with options for running the vignette locally (on your own computer) in R, looking along with a non-interactive .HTML file, the main code itself as a .R file, and a guide for installing the key R software packages such as Deseq2.
- [Dropbox link to vignette](https://www.dropbox.com/scl/fo/uyo4mtqp9aze1u5ckphe2/h?rlkey=hni2pwjl9p3tiwewpkpljjl82&dl=0)

For local use, the relevant files are called "tutorialCode.R" for interactive code, or the .HTML files for viewing only. Note that required software packages are pre-installed and configured on HPC, while for local interactive use, you must install the software before the session. For HPC use, see the guide from the top of this page.

<br />


### 04 Dr. Srilakshmi Raj and David Yang → Ancestry analysis

Vignettes are now available. 

For HPC option: If you used git before, you should be able to navigate to the folder and just use `git pull`. If you used the folder mounting method (ie the data.einsteinmed.edu) method, please follow instructions on the [vignette website](https://github.com/EinsteinOmics/OmicsWorkshop).

For local use, navigate to [this folder](https://github.com/EinsteinOmics/OmicsWorkshop/tree/main/04_AncestryAnalysis_SriDavid), click on the .html files, click the Download button top right of screen, which looks like an arrow pointing downward. Or, click the button with three dots on top right and click download.



<br />

### 05 Dr. Lindsay LaFave and Anthony Griffen → scATAC-seq

There will only be a view-only .HTML vignette for this session. Please navigate to [this folder](https://github.com/EinsteinOmics/OmicsWorkshop/tree/main/05_scATAC_LindsayAnthony), click on the .html file, click the Download button top right of screen, which looks like an arrow pointing downward. Or, click the button with three dots on top right and click download.


<br />

### 06 Dr. Deyou Zheng and Alexander Ferrena → scRNAseq

- For view only: navigate to [this folder](https://github.com/EinsteinOmics/OmicsWorkshop/tree/main/06_scRNA_DeyouAlex), and download the .HTML file: click on the .html file, then click the Download button top right of screen, which looks like an arrow pointing downward, or click the button with three dots on top right and click download.
- For interactive use on HPC: see instructions above for details to either use git clone / git pull, or download the folder and transfer to HPC using samba.


<br />

### 07 Dr. Simone Sidoli and Juan Sepulveda → Proteomics

- For view only: navigate to [this folder](https://github.com/EinsteinOmics/OmicsWorkshop/tree/main/07_Proteomics_SimoneJuan), and download the .HTML file: click on the .html file, then click the Download button top right of screen, which looks like an arrow pointing downward, or click the button with three dots on top right and click download.
- For interactive use on HPC: see instructions above for details to either use git clone / git pull, or download the folder and transfer to HPC using samba.
