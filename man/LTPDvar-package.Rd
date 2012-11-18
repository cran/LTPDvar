\name{LTPDvar-package}
\alias{LTPDvar-package}
\alias{LTPDvar}
\docType{package}
\title{
LTPD and AOQL single sampling plans for inspection by variables
}
\description{Calculation and evaluation of LTPD and AOQL plans for sampling inspection by variables which minimize mean inspection cost per lot of process average quality 
}
\author{
Nikola Kasprikova

Maintainer: Nikola Kasprikova <data@tulipany.cz>
}

\section{Introduction}{
Assume that measurements  of a single
quality characteristic \eqn{X} are independent, identically distributed
normal random variables with unknown parameters \eqn{\mu}  and \eqn{\sigma^2}.
 For the quality characteristic \eqn{X}  either  an upper specification
limit \eqn{U} is given (the  item is defective (nonconforming) if its  measurement exceeds \eqn{U}), or
a lower  specification  limit \eqn{L}  is given (the  item  is defective if its
measurement  is  smaller  than  \eqn{L}).  It  is  further assumed that
the unknown  parameter  \eqn{\sigma}  is  estimated  using  the  sample  standard
deviation \eqn{s}.

The inspection procedure is as follows:

 Draw a random sample of \eqn{n} items and compute
 \eqn{\bar{x}} and \eqn{s}.

 Accept the lot if \deqn{          {{U - \bar{x}} \over s } \ge k} or
         \deqn{   {{\bar{x} - L}\over s} \ge k. }

The operating  characteristic (see \code{\link{OC}}) is 
\deqn{                                   
L(p;n,k) = \int_{k\sqrt n}^\infty g(t;n-1,u_{1-p}\sqrt n) \,dt,}

where \eqn{                              
g(t;n-1,u_{1-p}\sqrt n) }
is probability density function of non-central \eqn{t} distribution with \eqn{(n-1)}
degrees of freedom and noncentrality parameter \eqn{\lambda=u_{1-p}\sqrt n}.

If case that we do not use exact formula for \code{\link{OC}} and we use the normal distribution
as an approximation of the non-central  \eqn{t}  distribution instead, we have
\deqn{
L(p;n,k) = \Phi \left({u_{1-p}-k \over A} \right),}
where
\deqn{
A = \sqrt{{1 \over n} + {k^2 \over 2(n-1)}}  .  }

The function \eqn{\Phi} is a standard normal distribution function
and \eqn{u_{1-p}} is a quantile of order \eqn{1-p}.
 
The task to be solved is determination of the sample size \eqn{n} and the critical
 value \eqn{k}.}    
\section{LTPD plans for acceptance sampling inspection by variables}{     
In case  of acceptance sampling by attributes (each  inspected item  is classified  as either  good or defective), there exist a procedure (Dodge and Romig, 1998) for finding sampling plans
which minimize  the  mean  number of items inspected per lot of process average quality
\deqn{I_s = N - (N-n)\cdot L(\bar{p};n,c) }
under the condition which protects   the  consumer   against  the
acceptance of a bad lot --  the probability
of accepting a submitted lot of tolerance quality \eqn{p_t} (consumer's
risk) shall be  0.10,
         \deqn{L(p_t;n,c) = 0.10  }
(LTPD single sampling plans),  where the given parameters are \eqn{N}, \eqn{\bar{p}}, \eqn{p_t}. 
\eqn{N} is the number of items in the lot,
\eqn{\bar{p}} is the process average fraction defective,
\eqn{p_t} is the lot tolerance fraction defective  (\eqn{P_t=100p_t} is the lot tolerance per cent defective -- denoted LTPD),
\eqn{n} is the number of items in the sample \eqn{(n<N)},
\eqn{c} is the acceptance number (the lot is rejected  when the number
     of defective items in the sample is greater than \eqn{c}),
\eqn{L(p)} is the operating characteristic
(the probability of accepting  a submitted lot
with fraction defective \eqn{p}).

LTPD  plans for inspection  by variables and attributes have been introduced in (Klufa, 1994).  Under the same protection of consumer, LTPD plan for inspection by variables and attributes  is  in many situations more economical with respect to inspection cost than the corresponding Dodge-Romig LTPD attribute sampling plan. 

For LTPD  plans for inspection  by variables and attributes (all items from
    the  sample are  inspected  by  variables, but  the remainder
    of rejected lots is inspected only by attributes), new parameter \eqn{c_m} is introduced, as 
the cost of inspection of one item by
variables divided by the cost of inspection of one item by attributes (usually is \eqn{c_m > 1}). Then the mean  inspection cost  per lot of process  average  quality is
\eqn{I_{ms}*c_a}, where \eqn{c_a} is the cost of inspection of one item by attributes and 
\deqn{I_{ms} = n\cdot c_m  +  (N-n)\cdot [1 - L(\bar{p};n,k)].} 
(see \code{\link{Ims}}). So we search for  the  acceptance  plan  \eqn{(n,k)}  minimizing
the mean inspection  cost per lot of process average quality (or equivalently minimizing \eqn{I_{ms}})
under the condition \eqn{L(p_t;n,k) = 0.10}. 

 If we use approximative formula for operating  characteristic, the equation \eqn{L(p_t;n,k)=0.10} has 
one and only one solution (see (Klufa, 1994))
 \deqn{k = {u_{1-p_t} - u_{0.10}.h \over g},}
where
\deqn{g = 1 - {u_{0.10}^2 \over 2(n-1)},}
\deqn{h = \sqrt{ {g \over n} + {u_{1-p_t}^2 \over 2(n-1)} }. }

If we use exact formula for \code{\link{OC}}, then exact solution of the equation \eqn{L(p_t;n,k)=0.10}
 is (see \code{\link{kL}})
\deqn{
k=\frac{t_{0.9}(n-1,u_{1-p_t}\sqrt n)}{\sqrt n},    }

where \eqn{t_{0.9}(n-1,u_{1-p_t}\sqrt n)} is a quantile of order 0.9  
of non-central \eqn{t} distribution with \eqn{(n-1)}
degrees of freedom and noncentrality parameter \eqn{\lambda=u_{1-p_t}\sqrt n}.

Then \eqn{I_{ms}} may be expressed as a function of one variable \eqn{n}
\deqn{
I_{ms}(n)=n\cdot c_m+(N-n)\cdot \alpha(n),     }
where \eqn{\alpha(n)} is the producer's risk  
(the probability of rejecting a lot of process average quality).

Function \code{\link{planLTPD}} searches for the sample size \eqn{n} minimizing \eqn{I_{ms}(n)} and gives plan with resulting \eqn{n} and corresponding \eqn{k} (see \code{\link{kL}}) as output.  In  \code{\link{planLTPD}} if \code{method="napprox"}, approximate \code{\link{OC}} is used and the solution is obtained using procedure described in (Klufa, 1994). If \code{method="exact"} (default), the optimization procedure searches for \eqn{n} in interval with centre at \eqn{n} resulting from \code{planLTPD(..., method = "napprox")}.

}
\section{AOQL plans for acceptance sampling inspection by variables}{
Under the assumption that each inspected item  is classified  as either  good or defective
(acceptance sampling by attributes) Dodge and Romig (1998) introduced sampling plans \eqn{(n, c)}
which minimize  the  mean  number of items inspected per lot of process average quality, assuming  that the remainder of rejected lots is inspected
\deqn{I_s = N - (N-n)\!\cdot\!L(\bar p;n,c)                       }
under the condition
\deqn{
\max_{0<p<1} AOQ(p) = p_L,                                    }
where \eqn{p_L} is the average outgoing quality limit (the given parameter) and AOQ is the average outgoing quality, i. e. the mean fraction defective  after  inspection (assuming that each defective item found is replaced by good one) when   the  fraction  defective  before  inspection was \eqn{p}.
Sampling plans for inspection by variables, which in comparison with sampling plans for inspection by attributes in many situations bring considerable savings in inspection cost, were then introduced in (Klufa, 1997).
 Function \code{\link{planAOQL}} searches for plan minimizing \eqn{I_{ms}(n)} under the condition that \code{\link{AOQ}} does not exceed the given value \eqn{p_L}.  In  \code{\link{planAOQL}} if \code{method="napprox"}, approximate \code{\link{OC}} is used and the solution is obtained using procedure described in (Klufa, 1997). If \code{method="exact"} (default), the optimization procedure searches for \eqn{n} in interval with centre at \eqn{n} resulting from \code{planAOQL(..., method = "napprox")}.

}
\references{
Dodge, H. F. - Romig, H. G.: Sampling  Inspection  Tables: Single and 
Double Sampling. John Wiley, 1998.

Klufa, J.:  Acceptance Sampling by Variables when the Remainder
  of Rejected Lots is Inspected. Statistical Papers, Vol.35, 337 - 349, 1994.

Klufa, J.: Exact calculation of the Dodge-Romig LTPD single sampling plans for inspection by variables. Statistical Papers, Vol. 51(2), 297-305, 2010.

Klufa J,: Dodge-Romig AOQL single sampling plans for inspection by variables. Statistical Papers 38: 111 - 119, 1997.
}
\keyword{ package }
\seealso{
\code{\link{planLTPD}}, \code{\link{planAOQL}}, \code{\link{kL}}, \code{\link{OC}}, \code{\link{AOQ}}, \code{\link{Ims}}
}
\examples{# calculation of LTPD plan
zz=planLTPD(N=1000,pt=0.1,pbar=0.001);zz
plot(zz);
# find critical value for given sample size
## and lot tolerance fraction defective 
kL(n=16,pt=0.1)
# create another plan
zz2=new("ACSPlan", n=16, k=2.71)
plot(zz2,xl=0.001, xu=0.15, xlabm="fraction defective",
ylabm="probability of acceptance",typem="l")
# calculation of AOQL plan
planAOQL(N=1000,pbar=0.005,pL=0.01, method="napprox", cm=1.5)
}
