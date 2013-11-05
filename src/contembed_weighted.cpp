
//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


//modified from https://gist.github.com/kevinushey/4561281
NumericVector row_sums( NumericMatrix& X ) {
int nRows = X.nrow();
NumericVector out = no_init(nRows);
for( int i=0; i < nRows; i++ ) {
NumericMatrix::Row tmp = X(i, _);
out[i] = sum( tmp );
}
return out;
}





void printvec(NumericVector vec)
{
  for (int i =0; i < vec.size(); i++)
    {
      std::cout << vec[i] << ",";
    }
  std::cout << "\n";
}

//Projection of point xi on the sphere of diameter dij centered on xj 
// [[Rcpp::export]]
NumericVector ProjectVec(NumericVector xi,NumericVector xj,double dij)
{
  NumericVector delta = xi-xj;
  double vnorm = sqrt(sum(pow(delta,2)));
  //  printvec(xi);printvec(xj);std::cout << "\n";
  if (vnorm == 0)
    {
      ProjectVec(rnorm(xi.size()),xj,dij);
    }
  else
    {
      //      std::cout << vnorm << "\n";
      return xj + delta * dij / vnorm;
    }
}

// Compute all surrogate points at a given timepoint
// [[Rcpp::export]]
NumericMatrix ProjectC(NumericMatrix Dtau,NumericMatrix Xtau, int i)
{
  int n=Xtau.ncol(),d=Xtau.nrow();
  NumericVector x = Xtau(_,i);
  NumericMatrix out(d,n);
  for (int j = 0; j < n; j++)
    {
      //std::cout << "j: " << j  << "n: " << n << "\n";
      if (j != (i))
	{
	  out(_,j) = ProjectVec(x,Xtau(_,j),Dtau(i,j));
	}
    }
  return out;
}


// Compute all surrogate curves
// [[Rcpp::export]]
List MajoriseC(List Dlist,List Xlist, int i,List params)
{
  int n = params["N"],d=params["D"],t=params["T"];
  List XauxList(t);
  for (int ind =0;ind < t; ind++)
    {
      XauxList(ind)=ProjectC(Dlist(ind),Xlist(ind),i);
    }

  return XauxList;
}

// [[Rcpp::export]]
NumericMatrix MinimiseC(List XauxList, int i, List params)
{
  int n = params["N"],d=params["D"],t=params["T"];
  //Step 1: compute time-wise mean of the surrogate curves
  arma::mat MeanAux(d,t);
  NumericVector x(d),m(d);
  
  arma::mat Regress = as<arma::mat>(params["Regress"]);
  arma::mat nx(d,t);
  for (int ind_t =0;ind_t < t; ind_t++)
    {
      NumericMatrix Xa = XauxList[ind_t];
      m = (row_sums(Xa))/(n-1);
      MeanAux.col(ind_t) = as<arma::colvec>(m);
    }

  //Step 2: Compute dimension-wise spline regression
  for (int ind_d =0;ind_d < d; ind_d++)
    {
      nx.row(ind_d) = MeanAux.row(ind_d)* Regress;
    }  

  return wrap(nx);
}

// [[Rcpp::export]]
void CenterConfiguration(List XList,List params)
{
  int n = params["N"],d=params["D"],t=params["T"],ind_row,ind_t;
  for (ind_t  =0;ind_t < t; ind_t++)
    {
      for (ind_row = 0;ind_row < d; ind_row++)
	{
	  NumericMatrix X = XList[ind_t];
	  X(ind_row,_) = X(ind_row,_) - sum(X(ind_row,_))/n;
	}      
    }
}

// [[Rcpp::export]]
NumericMatrix MinimiseWeighted(List XauxList, List WList, int i, List params)
{
  int n = params["N"],d=params["D"],t=params["T"],ind_row;
  //Step 1: compute time-wise mean of the surrogate curves
  arma::mat MeanAux(d,t);
  NumericVector x(d),m(d);
  arma::mat Regress = as<arma::mat>(params["Regress"]);
  arma::mat nx(d,t);
  for (int ind_t =0;ind_t < t; ind_t++)
    {
      NumericMatrix Xa = XauxList[ind_t];
      NumericMatrix W=WList[ind_t];
      for (ind_row = 0;ind_row < d; ind_row++)
	{
	  //	  
	  m[ind_row] = sum(Xa(ind_row,_)*W(i,_));
	  //	  std::cout<< "ind_row " << ind_row << " m[ind_row] " << m[ind_row] << "\n";
	}
      MeanAux.col(ind_t) = as<arma::colvec>(m);
    }

  //Step 2: Compute dimension-wise spline regression
  for (int ind_d =0;ind_d < d; ind_d++)
    {
      nx.row(ind_d) = MeanAux.row(ind_d)* Regress; //Matrix product, not element-wise
    }  
  return wrap(nx);
}



// Inner-loop: optimise a given curve x_i(t)
// [[Rcpp::export]]
int OptimiseCurve(List Dlist,List Xlist, int i, List params)
{
  int n = params["N"],d=params["D"],t=params["T"],weighted=params["weighted"];
  int nIter = 1,  maxIter = 12;
  double eps = params["eps"], delta = 2*eps;
  NumericMatrix x(d,t),old_x(d,t);
  List WList;
  if (weighted)
    {
      WList = params["WL"];
    }

  while ((delta > eps) & (nIter < maxIter))
    {
      if (weighted)
	{
	  x = MinimiseWeighted(MajoriseC( Dlist,Xlist,i,params),WList,i,params);	}
      else
	{
	  x = MinimiseC(MajoriseC( Dlist,Xlist,i,params),i,params);
	}

      
      if (nIter > 1)
	{
	  delta = mean(abs(x-old_x))/mean(abs(x));
	}
      else
	{
	  delta = 2.0* (double) params["eps"];
	}
      old_x = x;
      for (int ind_t =0;ind_t < t; ind_t++)
	{
	  NumericMatrix M = Xlist[ind_t];
	  NumericMatrix::Column xc=M(_,i);
	  xc = x(_,ind_t);
	}	  
      nIter++;
      //      std::cout << "delta: " << delta << " params[eps] " << (double) params["eps"] <<  " x[0] "<< x[0] << "\n";
    }

  return nIter;
}

// [[Rcpp::export]]
void OptimiseAllCurves(List Dlist,List Xlist, List params)
{
  int n = params["N"],d=params["D"],t=params["T"];
  //  int nIter = 1,  maxIter = 12;
  for (int ind_curve = 0;ind_curve < n;ind_curve++)
    {
      OptimiseCurve(Dlist,Xlist, ind_curve, params);
    }

}
// [[Rcpp::export]]
int OuterLoop(List Dlist,List Xlist, List params,int maxIter)
{
  int n = params["N"],d=params["D"],t=params["T"];
  int nIter = 0;
  while ((nIter < maxIter))
    {
      OptimiseAllCurves(Dlist,Xlist,params);
      CenterConfiguration(Xlist,params);
      nIter++;
    }
}
