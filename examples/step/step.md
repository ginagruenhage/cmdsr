#Weird effect using stepfunctions

Let's build a distance matrix based constant functions with even spacing:


```r
DL1 <- replicate(5, a.m(dist(c(0, 1, 2, 3, 4))))
DL1 <- alply(DL1, 3, function(d) d)
```

As expected, cmds can perfectly embed these distances in k=1 dimensions:

```r
plot.cmds(cmds(DL1, k = 1, l = 0))
```

```
## Total cost C:  0 
## The algorithm converged. (delta =  0 )
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 


Let's glue a second set of constant functions to the first set, but with double the spacing.


```r
DL2 <- replicate(5, a.m(dist(c(0, 2, 4, 6, 8))))
DL2 <- alply(DL2, 3, function(d) d)
plot.cmds(cmds(DL2, k = 1, l = 0))
```

```
## Total cost C:  0 
## The algorithm converged. (delta =  0 )
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-31.png) 

```r
plot.cmds(cmds(c(DL1, DL2), k = 1, l = 0))
```

```
## Total cost C:  0 
## The algorithm converged. (delta =  0 )
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-32.png) 

That also works fine.

Now, let's take functions with random y intercept:

```r
set.seed(3)
vec <- rnorm(5, sd = 2)
DL3 <- replicate(5, a.m(dist(vec)))
DL3 <- alply(DL3, 3, function(d) d)
plot.cmds(cmds(DL3, k = 1, l = 0))
```

```
## Total cost C:  1.233e-31 
## The algorithm converged. (delta =  0 )
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 

That works. However, if we combine this distance list with one of the previous one, the algorithm fails in the second case:

```r
plot.cmds(cmds(c(DL1, DL3), k = 1, l = 0), shepard = TRUE)
```

```
## Total cost C:  1.233e-31 
## The algorithm converged. (delta =  0 )
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-51.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-52.png) 

```r
res <- cmds(c(DL2, DL3), k = 1, l = 0, v = T)
```

```
## Initialization: Total cost C:  278.3 
## Iteration:  1 , total cost C:  37.09 
## Iteration:  2 , total cost C:  31.02 
## Iteration:  3 , total cost C:  30.99 
## Iteration:  4 , total cost C:  30.99 
## Iteration:  5 , total cost C:  30.99 
## Iteration:  10 , total cost C:  30.99 
## Iteration:  15 , total cost C:  30.99 
## Iteration:  20 , total cost C:  30.99 
## Iteration:  25 , total cost C:  30.99 
## Iteration:  30 , total cost C:  30.99 
## Iteration:  35 , total cost C:  30.99 
## Iteration:  40 , total cost C:  30.99 
## Iteration:  45 , total cost C:  30.99 
## Iteration:  50 , total cost C:  30.99 
## Total cost C:  30.99 
## The algorithm converged. (delta =  0 )
```

```r
plot.cmds(res, shepard = TRUE)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-53.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-54.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-55.png) ![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-56.png) 

To check wether this result is a local minimum of the cost function, I wrote three auxiliary functions:

```r
# transform the embedding list to a numeric vector
XL2num <- function(XL) {
    a.n(laply(XL, function(x) as.numeric(x)))
}

# wrapper for C.L (which computes the cost in the cmds algorithm), that
# takes a numeric vector as first argument
Cost <- function(x, DL, params) {
    XL <- x2list(x, params)
    C.L(DL, XL, params)
}

# calculates the numerical gradient and checks whether it's equal to zero
# with numerical tolerance
check.local.minimum <- function(res) {
    x <- XL2num(res$XL)
    n <- length(x)
    g <- grad(Cost, x, DL = res$DL, params = res$params)
    
    if (all.equal(rep(0, n), g)) {
        cat("The result is a local minimum of the cost function.\n")
    } else {
        cat("The automatic test for local minimum failed. \n You are now in a browser to check the actual value of the numerical gradient, g.\n")
        browser()
    }
}
```

The weird result is really a local minimum:

```r
check.local.minimum(res)
```

```
## The result is a local minimum of the cost function.
```


Iit seems to depend on the values in the random vector wether it works or not, for other random vectors I didn't get this effect. For example:

```r
set.seed(20)
vec <- rnorm(5, sd = 2)
DL3 <- replicate(5, a.m(dist(vec)))
DL3 <- alply(DL3, 3, function(d) d)
plot.cmds(cmds(c(DL1, DL3), k = 1, l = 0), shepard = TRUE)
```

```
## Total cost C:  8.135e-31 
## The algorithm converged. (delta =  0 )
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-81.png) ![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-82.png) 

```r
plot.cmds(cmds(c(DL2, DL3), k = 1, l = 0), shepard = TRUE)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-83.png) 

```
## Total cost C:  8.135e-31 
## The algorithm converged. (delta =  0 )
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-84.png) ![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-85.png) ![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-86.png) 

