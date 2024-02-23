# Intro to Omics Workshop Vignettes

An asynchronous learning workshop introducing key topics in Omics.

<br />

## General instructions

Each vignette will be made avaialable in three ways:
- Jupyter Notebook files (.ipynb), for those who want to follow along interactively using the HPC Jupyter workflow, where all software has been pre-configured
- local RMarkdown or Jupyter files, for those who want to follow along interactively using on their local systems
- pre-run HTML files, for those who want to follow along locally without any code

<br />
<br />


## HPC guide

To follow along on HPC, it is absolutely necessary that you sign up for an HPC account, and access https://cluster.einsteinmed.edu/. It is also recommended to familiarize yourself with the remote access option (GlobalProtect) for HPC use outside of the EinsteinSecure wifi.

- [HPC account setup](https://montefioreorg.sharepoint.com/sites/Einstein-IT-HPC/SitePages/HPC3.0-UQuick-Start.aspx) - login required. It has all the relevant guides including how to set up HPC account.
- [cluster.einsteinmed.edu usage](https://montefioreorg.sharepoint.com/sites/Einstein-IT-HPC/Shared%20Documents/Forms/AllItems.aspx?id=%2Fsites%2FEinstein%2DIT%2DHPC%2FShared%20Documents%2FGeneral%2FHPC3%2E0%20docs%2FIntroduction%20to%20Einstein%20HPC%20Portal%2Dv5%2Epdf&parent=%2Fsites%2FEinstein%2DIT%2DHPC%2FShared%20Documents%2FGeneral%2FHPC3%2E0%20docs) - try to follow this and activate the Jupyter Notebook workflow.
- [Remote access to HPC via GlobalProtect](https://montefioreorg.sharepoint.com/sites/Einstein-IT-HPC/Shared%20Documents/Forms/AllItems.aspx?id=%2Fsites%2FEinstein%2DIT%2DHPC%2FShared%20Documents%2FGeneral%2FHPC3%2E0%20docs%2FIT%2DREF%2D2023%2D094%20Einstein%20Academic%20Research%20Systems%20Portal%282%29%2Epdf&parent=%2Fsites%2FEinstein%2DIT%2DHPC%2FShared%20Documents%2FGeneral%2FHPC3%2E0%20docs)

If you run into errors, you can reach out to us or submit an IT ticket.

<br />
<br />

## Instructions for downloading code vignettes

For HPC usage: first, try to activate the Jupyter workflow on https://cluster.einsteinmed.edu/. This will create a folder in your HPC home directory called "pw". It is recommended to store all vignette code in there to be easily accessible to Jupyter.

### 1. Unix shell (bash / zsh / terminal) with `git clone`. 

- First, open a Unix shell, such as Mac terminal. On Windows, there are several options including [Git Bash](https://git-scm.com/download/win), [Cygwin](https://www.cygwin.com/), and [WSL](https://learn.microsoft.com/en-us/windows/wsl/about). 
- Next, use `ssh` to get to the HPC. Follow the instructions above, and make sure to use GlobalProtect if off campus. In short, you must use `ssh username@username.hpc.einsteinmed.edu` (replace "username" with your active domain).
- Once there, use the "cd" command to get to the pw folder in your home directory: `cd ~/pw/`.
- Finally, you can download the scripts: `git clone https://github.com/EinsteinOmics/OmicsWorkshop.git`

<br />

### 2. Folder Mounting to HPC just like data.einsteinmed.edu

This makes use of the folder mounting, identical to how data.einsteinmed.edu is used.
- [Folder mounting guide](https://it.einsteinmed.edu/documentation/how-to-mount-the-hpc-file-system/) - Note that "einsteinmed.**edu**" should be used in all instances instead of "einsteinmed.**org**" in this guide.

First, download the vigenttes from Github. This can be done by clicking the big green button that says "Code" at the top of the screen. Then click "Download Zip". Make sure to unzip the file. Then, use the guide above to move the folder to HPC. You should place it in the "pw" folder in your home directory.

<br />
<br />


# Further Guides for Specific Vignettes


### 03 Bulk RNA-seq and bulk ATAC-seq with Dr. Kith Pradhan (assisted by Xiang Yu Zheng)

Dr. Pradhan has made available a Dropbox link with options for locally the vignette locally, looking along with a non-interactive .HTML file, the main code itself as a .R file, and a guide for installing the key R software packages such as Deseq2.
- [Dropbox link to vignette](https://www.dropbox.com/scl/fo/uyo4mtqp9aze1u5ckphe2/h?rlkey=hni2pwjl9p3tiwewpkpljjl82&dl=0)

