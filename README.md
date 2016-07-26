Algorithms for the Compound Latent Dirichlet Allocation Model
=============================================================

This **R** package, **clda**, implements the following Markov chain Monte Carlo 
(MCMC) and variational methods for the compound latent Dirichlet allocation 
(cLDA) model. 
* Auxiliary variable update within Gibbs sampler (AGS)
* Metropolis adjusted Langevin algorithm within Gibbs sampler (MGS)
* Variational Expectation Maximization (VEM)


For package documentation run 

``` help("clda") ```

in an R console. All major functions and datasets are documented and linked to 
the package index. 

To see all demo R scripts available in this package, run 

``` demo(package="clda") ```

in an R console. These demo scripts require commandline arguments for execution. 
Please see the documentation provided in each script before execution.    

Authors
----------------------------
* [Clint P. George](https://clintpgeorge.wordpress.com) (Please contact for questions and comments)
* Wei Xia
* George Michailidis 

Dependencies
----------------------------

This package uses the following R packages, which are already included in this R package.   
* **Rcpp**
* **RcppArmadillo** based on the **Armadillo** C++ package 
* **lattice**

Installation Guide 
------------------

* Download the package source from [Git Download Link](https://github.com/clintpgeorge/clda/archive/master.zip)
* Unzip the dowloaded file and rename the folder **clda-master** to **clda** 
* To install **clda** run ```R CMD INSTALL clda``` on the commandline 
* To uninstall **clda** run ```R CMD REMOVE clda``` on the commandline 

