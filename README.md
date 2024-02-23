# Intro to Omics Workshop Vignettes

An asynchronous learning workshop introducing key topics in Omics.

## General instructions

Each vignette will be made avaialable in three ways:
- Jupyter Notebook files (.ipynb), for those who want to follow along interactively using the HPC Jupyter workflow, where all software has been pre-configured
- local RMarkdown or Jupyter files, for those who want to follow along interactively using on their local systems
- pre-run HTML files, for those who want to follow along locally without any code


## HPC guide

To follow along on HPC, it is absolutely necessary that you sign up for an HPC account, and access https://cluster.einsteinmed.edu/. It is also recommended to familiarize yourself with the remote access option (GlobalProtect) for HPC use outside of the EinsteinSecure wifi.

- [HPC account setup](https://montefioreorg.sharepoint.com/sites/Einstein-IT-HPC/SitePages/HPC3.0-UQuick-Start.aspx) - login required. It has all the relevant guides including how to set up HPC account.
- [cluster.einsteinmed.edu usage](https://montefioreorg.sharepoint.com/sites/Einstein-IT-HPC/Shared%20Documents/Forms/AllItems.aspx?id=%2Fsites%2FEinstein%2DIT%2DHPC%2FShared%20Documents%2FGeneral%2FHPC3%2E0%20docs%2FIntroduction%20to%20Einstein%20HPC%20Portal%2Dv5%2Epdf&parent=%2Fsites%2FEinstein%2DIT%2DHPC%2FShared%20Documents%2FGeneral%2FHPC3%2E0%20docs) - try to follow this and activate the Jupyter Notebook workflow.
- [Remote access to HPC via GlobalProtect](https://montefioreorg.sharepoint.com/sites/Einstein-IT-HPC/Shared%20Documents/Forms/AllItems.aspx?id=%2Fsites%2FEinstein%2DIT%2DHPC%2FShared%20Documents%2FGeneral%2FHPC3%2E0%20docs%2FIT%2DREF%2D2023%2D094%20Einstein%20Academic%20Research%20Systems%20Portal%282%29%2Epdf&parent=%2Fsites%2FEinstein%2DIT%2DHPC%2FShared%20Documents%2FGeneral%2FHPC3%2E0%20docs)

If you run into errors, you can reach out to us or submit an IT ticket.


## Instructions for downloading code vignettes

For HPC usage: first, try to activate the Jupyter workflow on https://cluster.einsteinmed.edu/. This will create a folder in your HPC home directory called "pw". It is recommended to store all folders in there.

1. Unix shell (bash / zsh / terminal) with `git clone`. 

First, open a Unix shell, such as Mac terminal. On Windows, there are several options including [Git Bash](https://git-scm.com/download/win), [Cygwin](https://www.cygwin.com/), and [WSL](https://learn.microsoft.com/en-us/windows/wsl/about)https://learn.microsoft.com/en-us/windows/wsl/about. 

Next, use `ssh` to get to the HPC. Follow the instructions above, and make sure to use GlobalProtect if off campus. In short, you must use `ssh <username>@<username>.hpc.einsteinmed.edu` (replace <username> with your active domain).

Once there, use the `cd` command to get to the pw folder in your home directory: `cd ~/pw/`

Finally, you can download the scripts: `git clone `

