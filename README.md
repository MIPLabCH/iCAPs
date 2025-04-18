**License**
This toolbox is released under an [[ ssh://git@c4science.ch/source/iCAPs.git | Apache 2.0 license ]].

**Installation**
The following toolboxes need to be installed and added to your path before using iCAPs:
* [[ https://www.fil.ion.ucl.ac.uk/spm/software/spm12/ | SPM12 ]]
* [[ https://www.nitrc.org/projects/cbinifti/ | cbiNifti ]]
* [[ https://epfl-lts2.github.io/gspbox-html/ | GSPBox ]]

**Link to zip**
A nightly build (zip) can be downloaded at [[ https://miplab.epfl.ch/index.php/software/total-activation | https://miplab.epfl.ch/index.php/software/total-activation ]]

**Description**
Toolbox with functions and scripts to obtain innovation-driven coactivation patterns (iCAPs) from fMRI data. Main steps of the pipeline are:
* Total activation for regularized deconvolution of fMRI data
* Thresholding of innovation frames based on surrogate data
* Clustering of innovation frames to extract spatial brain activity patterns
* Spatio-temporal regression for time course recovery of networks

For a smooth start, we suggest to download example data for two subjects from [[ https://drive.google.com/drive/folders/17XhZ_9eJ9X65imFtaHg48WjAZRLHkRUJ?usp=sharing]]. Look at the example main and input scripts to see how to run the complete pipeline on the example fMRI data of two subjects. Note that the results are not stable for this few data, the example only exists to demonstrate how input and output data is saved.

**Installation and updates**
For installation, clone this repository using git:
```
git clone --depth=1 https://c4science.ch/source/iCAPs
```
Cloning may take a while, because the example data takes up some space (~5GB). To get the latest update of the toolbox, use 
```
git pull
```

**References**
For details about the algorithms please refer to:
* Using iCAPs to reveal functional networks
  - Karahanoglu, F. I. and Van De Ville, D. (2015). [[ https://www.nature.com/articles/ncomms8751 | Transient brain activity disentangles fMRI resting-state dynamics in terms of spatially and temporally overlapping networks. ]] Nat. Commun., 6:7751.
 - Van De Ville, D. and Karahanoglu, F. I. (2016). [[ http://spie.org/newsroom/6521-resting-state-neuroimaging-unravels-functional-organization-in-the-brain?highlight=x2416&ArticleID=x119771&SSO=1 | Resting-State Neuroimaging Unravels Functional Organization in the Brain. ]] SPIE Newsroom, August 15.

* Total activation deconvolution:
  - Karahanoglu, F. I., Bayram, I., and Van De Ville, D. (2011).[[ https://miplab.epfl.ch/pub/karahanoglu1101.pdf |  A Signal Processing Approach to Generalized 1-D Total Variation. ]] IEEE Trans. Signal Process., 59(11):5265{5274.
  - Karahanoglu, F. I., Caballero-Gaudes, C., Lazeyras, F., and Van De Ville, D. (2013). [[ https://miplab.epfl.ch/pub/karahanoglu1302.pdf | Total activation: FMRI deconvolution through spatio-temporal regularization. ]] Neuroimage, 73:121{134.
  - Farouj, Y., Karahanoglu, F. I., and Van De Ville, D. (2017). [[ https://miplab.epfl.ch/pub/farouj1701.pdf | Regularized Spatiotemporal Deconvolution of fMRI Data Using Gray-Matter Constrained Total Variation. ]] Proc. 14th IEEE Int. Symp. Biomed. Imaging From Nano to Macro, pages 472-475.

* Transient-informed regression:
  - Zöller, D. M., Bolton, T. A. W., Karahanoglu, F. I., Eliez, S., Schaer, M., and Van De Ville, D. V. D. (in press). Robust recovery of temporal overlap between network activity using transient-informed spatiotemporal regression. IEEE Transactions on Medical Imaging

* More background about dynamic functional connectivity:
 - Preti, M. G., Bolton, T., Van De Ville, D. (2017). [[ https://miplab.epfl.ch/pub/preti1701.pdf | The Dynamic Functional Connectome: State-of-the-Art and Perspectives. ]] NeuroImage 160, pages 41-54.
 - Karahanoglu, F.I., Van De Ville, D. (2017). [[ https://www.sciencedirect.com/science/article/pii/S2468451117300417?via%3Dihub | Dynamics of Large-Scale fMRI Networks: Deconstruct Brain Activity to Build Better Models of Brain Function. ]] Current Opinion in Biomedical Engineering 3, pages 28-36.

**Contact**
If you have any questions, error reports or if you simply want to chat about iCAPs, feel free to contact: flavia.petruso@epfl.ch

**History**
Initial algorithm & code development: Isik Karahanoglu
Toolbox v1.0 implementation: Thomas Bolton
Implementation of public toolbox: Daniela Zöller
Current maintenance & contact: Petruso Flavia