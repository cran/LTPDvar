OC <-
function (p, n, k, type = c("t", "napprox")) 
{
    type = match.arg(type)
    Lnapprox <- function(p, n, k) pnorm((qnorm(1 - p) - k)/((1/n) + 
        ((k^2)/(2 * n - 2)))^0.5)
    Lt <- function(p, n, k) 1 - pt(q = k * n^0.5, df = n - 1, 
        ncp = qnorm(1 - p) * n^0.5)
    if (type == "t") 
        return(Lt(p, n, k))
    if (type == "napprox") 
        Lnapprox(p, n, k)
}
