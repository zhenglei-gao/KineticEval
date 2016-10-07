\name{NEWS}
\title{News for Package 'KineticEval'}
\encoding{UTF-8}
\section{Changes in KineticEval version 1.0-46 (2016-10-07)}{
  \subsection{fixed}{
    \itemize{
      \item #39, Graphing Issues by Kari from SCC.
    }
  }
}
\section{Changes in KineticEval version 1.0-45 (2016-02-19)}{
  \subsection{fixed}{
    \itemize{
      \item #39, Brigitte's DFOP parent was split when there is back transformation.
    }
  }
}
\section{Changes in KineticEval version 1.0-44 (2016-01-14)}{
  \subsection{fixed}{
    \itemize{
      \item #38, Carolin's DFOP not working because of DFOP DT90 calculations where
        uniroot function cannot find the root due to wrongly used min and max values
    }
  }
}
\section{Changes in KineticEval version 1.0-43 (2015-12-18)}{
  \subsection{fixed}{
    \itemize{
      \item #37, Carolin's DFOP not working because of DFOP DT50 DT90 calculations
    }
  }
}
\section{Changes in KineticEval version 1.0-42 (2015-10-16)}{
  \subsection{fixed}{
    \itemize{
      \item #36, No summary files after successful fit
    }
  }
}
\section{Changes in KineticEval version 1.0-41 (2014-12-11)}{
  \subsection{fixed}{
    \itemize{
      \item #32, KinGui to KinGUII.
    }
  }
}
\section{Changes in KineticEval version 1.0-40 (2014-05-12)}{
  \subsection{fixed}{
    \itemize{
      \item #31, first observed data not at time 0.
    }
  }
}
\section{Changes in KineticEval version 1.0-39 (2014-04-24)}{
  \subsection{changed}{
    \itemize{
      \item scientic format is set to scipen=10 for DT50 printing
    }
  }
}
\section{Changes in KineticEval version 1.0-38 (2014-04-24)}{
  \subsection{fixed}{
    \itemize{
      \item Water-Sediment setting returns degradation rates instead of dissipation
        rates now
    }
  }
  \subsection{changed}{
    \itemize{
      \item scientic format is set to scipen=10 for DT50 printing
    }
  }
}
\section{Changes in KineticEval version 1.0-37 (2014-04-23)}{
  \subsection{fixed}{
    \itemize{
      \item ODE numerical solutions with NaNs.
      \item Water-Sediment setting returns degradation rates instead of dissipation
        rates now
    }
  }
}
\section{Changes in KineticEval version 1.0-36 (2014-03-31)}{
  \subsection{added}{
    \itemize{
      \item More Demos: water-sediment and MCMC
    }
  }
}
\section{Changes in KineticEval version 1.0-35 (2014-03-31)}{
  \subsection{added}{
    \itemize{
      \item option to use Jacobian or hessian to calculate covariance matrix.
      \item Index for Demo files
    }
  }
  \subsection{fixed}{
    \itemize{
      \item MCMC formation fraction calculation not updated
    }
  }
}
\section{Changes in KineticEval version 1.0-33 (2014-03-31)}{
  \subsection{changed}{
    \itemize{
      \item warning text.
    }
  }
  \subsection{fixed}{
    \itemize{
      \item water-sediment bug.
    }
  }
}
\section{Changes in KineticEval version 1.0-32 (2014-03-27)}{
  \subsection{fixed}{
    \itemize{
      \item output incorrect due to no single formation fractions.
    }
  }
}
\section{Changes in KineticEval version 1.0-31 (2014-03-24)}{
  \subsection{fixed}{
    \itemize{
      \item stop message show up in the comments, not necessary.
      \item output incorrect for single Parent due to no transformed formation
        fractions.
    }
  }
}
\section{Changes in KineticEval version 1.0-30 (2014-03-06)}{
  \subsection{fixed}{
    \itemize{
      \item stop message show up in the comments, not necessary.
    }
  }
}
\section{Changes in KineticEval version 1.0-29 (2014-03-06)}{
  \subsection{added}{
    \itemize{
      \item Function MCstd to calculate the standard error of formation fractions
    }
  }
  \subsection{fixed}{
    \itemize{
      \item #28, formation fraction standard error calculation by faked posterior
        distributions
    }
  }
}
\section{Changes in KineticEval version 1.0-28 (2014-03-06)}{
  \subsection{fixed}{
    \itemize{
      \item #28, formation fraction standard error calculation by faked posterior
        distributions
    }
  }
}
\section{Changes in KineticEval version 1.0-27 (2014-03-05)}{
  \subsection{added}{
    \itemize{
      \item Covariance calculation default by Jacobian
    }
  }
  \subsection{fixed}{
    \itemize{
      \item #27, not trustable CI
      \item added model message
    }
  }
}
\section{Changes in KineticEval version 1.0-25 (2014-02-10)}{
  \subsection{changed}{
    \itemize{
      \item mkin move to imports
    }
  }
}
\section{Changes in KineticEval version 1.0-24 (2014-02-10)}{
  \subsection{added}{
    \itemize{
      \item ChangeLog
    }
  }
  \subsection{changed}{
    \itemize{
      \item NEWS
    }
  }
}
\section{Changes in KineticEval version 1.0-23 (2014-02-10)}{
  \subsection{changed}{
    \itemize{
      \item using roxyPackage instead for building pachages
    }
  }
  \subsection{fixed}{
    \itemize{
      \item Metabolite as a SINK
      \item Chi2Error being Inf or NaN
    }
  }
}
\section{Changes in KineticEval version 1.0-22 (2014-02-10)}{
  \subsection{added}{
    \itemize{
      \item Metabolite as a SINK
    }
  }
}
\section{Changes in KineticEval version 1.0-21 (2014-01-24)}{
  \subsection{fixed}{
    \itemize{
      \item Multiple Metabolite and no SINK
    }
  }
}
\section{Changes in KineticEval version 1.0-20 (2013-10-27)}{
  \subsection{fixed}{
    \itemize{
      \item Chi2 error NA when df <=0.
    }
  }
}
\section{Changes in KineticEval version 1.0-19 (2013-10-27)}{
  \subsection{added}{
    \itemize{
      \item A fit plot in pdf format is generated for reporting.
    }
  }
  \subsection{fixed}{
    \itemize{
      \item Output KineticEval Effectively fixed: formation fractions.
    }
  }
}
\section{Changes in KineticEval version 1.0-18 (2013-10-27)}{
  \subsection{changed}{
    \itemize{
      \item Add eval Methods in output
    }
  }
}
\section{Changes in KineticEval version 1.0-17 (2013-10-27)}{
  \subsection{added}{
    \itemize{
      \item Add an option to compute the exact Hessian using `numDeriv`, produce
        p-value 0.5 instead of NA's.
    }
  }
  \subsection{changed}{
    \itemize{
      \item Warning messages modified.
    }
  }
}
\section{Changes in KineticEval version 1.0-16 (2013-10-27)}{
  \subsection{added}{
    \itemize{
      \item Issue an warning when changing the initial values.
    }
  }
  \subsection{fixed}{
    \itemize{
      \item M0 not fixed: at 100 for parent substance.
    }
  }
}
\section{Changes in KineticEval version 1.0-15 (2013-10-27)}{
  \subsection{added}{
    \itemize{
      \item No added:
    }
  }
  \subsection{fixed}{
    \itemize{
      \item CI different from different optimization, NA or not NA?.
    }
  }
}
\section{Changes in KineticEval version 1.0-14 (2013-10-27)}{
  \subsection{added}{
    \itemize{
      \item No added:
    }
  }
  \subsection{fixed}{
    \itemize{
      \item changed: for `KinReport` function.
    }
  }
}
\section{Changes in KineticEval version 1.0-13 (2013-10-27)}{
  \subsection{added}{
    \itemize{
      \item No added:
    }
  }
  \subsection{fixed}{
    \itemize{
      \item #15, MCMC bug for incomplete data.
    }
  }
}
\section{Changes in KineticEval version 1.0-12 (2013-10-27)}{
  \subsection{added}{
    \itemize{
      \item No added:
    }
  }
  \subsection{fixed}{
    \itemize{
      \item internal function `norm` to `KineticEval:::norm`
      \item #15, test cases for `ginv`
    }
  }
}
\section{Changes in KineticEval version 1.0-11 (2013-10-27)}{
  \subsection{changed}{
    \itemize{
      \item Add routines generating logging files. #14
      \item Add routines generating report files, .kgg, .kgo when the model is not
        correctly set up or when the optimization could not be finished. #14
    }
  }
}
\section{Changes in KineticEval version 1.0-10 (2013-10-27)}{
  \subsection{added}{
    \itemize{
      \item Add Comment lines in the output summary file.
    }
  }
}
\section{Changes in KineticEval version 1.0-9 (2013-10-27)}{
  \subsection{fixed}{
    \itemize{
      \item #10. using generalized inverse to calculate inverse of the hessian.
      \item #16, with Generalized inverse, output NA where is appropriate.
    }
  }
}
\section{Changes in KineticEval version 1.0-8 (2013-10-27)}{
  \subsection{fixed}{
    \itemize{
      \item #11. TRR problematic case
      \item #6. Ghost compartment is correctly handled now.
      \item #5. 3 DFOP models.
    }
  }
}
\section{Changes in KineticEval version 1.0-7 (2013-10-27)}{
  \subsection{added}{
    \itemize{
      \item No added:
    }
  }
  \subsection{fixed}{
    \itemize{
      \item #1, FOMC not calculated.
    }
  }
}
\section{Changes in KineticEval version 1.0-5 (2013-10-27)}{
  \subsection{added}{
    \itemize{
      \item No added:
    }
  }
  \subsection{fixed}{
    \itemize{
      \item fixed: starting value problem. starting value reported by lab people need
        to be changed:.
    }
  }
}
\section{Changes in KineticEval version 1.0-4 (2013-10-27)}{
  \subsection{added}{
    \itemize{
      \item added: a function `summay.kingui` to summarize a model
      \item added: a few example data sets
      \item added: a demo named `Complex_Model`
    }
  }
}
\section{Changes in KineticEval version 1.0-3 (2013-10-27)}{
  \subsection{fixed}{
    \itemize{
      \item #1 FOMC model result is not calculated.
      \item #2 DFOP DT90 is not correctly calculated.
    }
  }
}
\section{Changes in KineticEval version 1.0-2 (2013-10-27)}{
  \subsection{fixed}{
    \itemize{
      \item #0, DFOP model for metabolite is re-established in the correct form.
    }
  }
}
\section{Changes in KineticEval version 1.0-1 (2013-10-27)}{
  \subsection{changed}{
    \itemize{
      \item in this NEWS file, #n means the issue number on GitHub, e.g. #2 is
        https://github.com/zhenglei-gao/KineticEval/issues/2
    }
  }
}
