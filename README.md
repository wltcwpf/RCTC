# RCTC
The Rotated Combination of Two-Component Ground Motions
For the details, please refer to the PEER report at the link: https://peer.berkeley.edu/sites/default/files/2017_09_stewart_9.10.18.pdf


# Installation
1, make sure "Rcpp" and "pracma" are installed as they are dependencies;
        install.packages(c('Rcpp', 'pracma'))

2, For Mac users: you may need to download Xcode for command line tools or (by running the following command in Terminal):
		xcode-select --install
   For Windows users: you may need to download and install Rtools from the following link (note remember to add path into the system environment path when installing):
        https://cran.r- project.org/bin/windows/Rtools/

3, install "devtools":
        install.packages("devtools")

4, lastly, install "RCTC" from github:
        library(devtools)
        install_github('wltcwpf/RCTC')
	
# Citation
If you use RCTC (directly or as a dependency of another package) for work resulting in an academic publication, we would appreciate if the report is cited:
Wang, P., Stewart, J. P., Bozorgina, Y., Boore, D. M., Kishida, T. (2017). “R” Package for Computation of Earthquake Ground Motion Response Spectra. Report No. 2017/09. PEER, UC Berkeley.
