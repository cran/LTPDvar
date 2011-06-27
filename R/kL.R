kL <-
function (n, pt, b = 0.1) 
{
    lambda <- function(n, p) qnorm(1 - p, mean = 0, sd = 1) * 
        n^0.5
    qt(p = 1 - b, df = n - 1, ncp = lambda(n, pt))/n^0.5
}
