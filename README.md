KineticEval - An R package for Kinetic Evaluations
===========

This package is developed based on the R package [`mkin`](https://github.com/jranke/mkin). It is used by the **KinGUII v2.1**, which is a successor of the first version of KinGUI. 

The package is used to determine kinetic parameters from results of environmental fate studies, e.g. aerobic soil degradation, by fitting respective mathematical models to the observed data.

The package allows the choice between different optimization algorithms. In particular the estimation of parameter confidence intervals is much improved/corrected as compared to other commonly used evaluation softwares by providing the methods Iteratively Reweighted Least Squares (IRLS) and Markov Chain Monte Carlo (MCMC).

## Installation

* **Officially released version**

    Official version comes together with a GUI which can be obtained from [here](http://kinguii.vrbka.net/KinGUIIv2.1.zip).


    All one need to do is to follow the link and unzip the downloaded folder to a local drive. After starting KinGUII v2.1 for the first time, you also need to set the path to R and to the Working Directory (File -> Preferences -> Paths). Please see the manual for some more details.

* **development version**

    
    ```r
    require(devtools)
    install_github("KineticEval", "zhenglei-gao")
    ```



## License

* Under the [GNU General Public License (GPL)](http://www.gnu.org/licenses/gpl.html)

## Wiki

* [Some Incomplete Information](https://github.com/zhenglei-gao/KineticEval/wiki)

## Questions, Bug reports, and new feature requests

* [Github Issues](https://github.com/zhenglei-gao/KineticEval/issues?page=1&state=open)

## Contributions

See [here](https://github.com/jranke/mkin) for a nice description of the histories and credits.


## Changes

See ChangeLog

## Disclaimer

**Warning: Still under development**, use it at your own risk. 







