planLTPD <-
function (N, pt, pbar, b = 0.1, cm = 1, method = c("exact", "napprox"), 
    intdif = 20) 
{
    intdif = floor(intdif)
    method = match.arg(method)
    nAPPROXL <- function(N, pt, b, pbar, cm) {
        if (N < 5) 
            stop("lot size < 5")
        k0L <- function(n, pt, b) {
            g = function(n, b) 1 - qnorm(b, mean = 0, sd = 1)^2/(2 * 
                n - 2)
            h = function(n, pt, b) (g(n, b)/n + qnorm(1 - pt, 
                mean = 0, sd = 1)^2/(2 * n - 2))^0.5
            return((qnorm(1 - pt, mean = 0, sd = 1) - qnorm(b, 
                mean = 0, sd = 1) * h(n, pt, b))/g(n, b))
        }
        alpha0 <- function(n, pt, b, pbar) pnorm((k0L(n, pt, 
            b) - qnorm(1 - pbar))/((1/n) + (k0L(n, pt, b)^2/(2 * 
            n - 2)))^0.5)
        fF = function(n, pt, b, pbar, cm) (cm - alpha0(n + 1, 
            pt, b, pbar))/(alpha0(n, pt, b, pbar) - alpha0(n + 
            1, pt, b, pbar)) + n
        n0 = (ceiling(uniroot(function(n) fF(n, pt, b, pbar, 
            cm) - N, c(5, N))$root))
        return(new("ACSPlan", n = n0, k = k0L(n0, pt, b)))
    }
    plan0 = nAPPROXL(N, pt, b, pbar, cm)
    if (method == "napprox") 
        return(plan0)
    if (method == "exact") {
        fOptimn <- function(N, pt, b, pbar, cm) {
            Imsf = function(n) Ims(n, N, pt, pbar, cm, b)
            fMinSearch = function(nl_, nu_) ifelse(nl_ == nu_, 
                nl_, ifelse(Imsf(nl_ + floor((nu_ - nl_)/2)) <= 
                  Imsf(nl_ + floor((nu_ - nl_)/2) + 1), fMinSearch(nl_, 
                  nl_ + floor((nu_ - nl_)/2)), fMinSearch(nl_ + 
                  floor((nu_ - nl_)/2) + 1, nu_)))
            init_ = n(nAPPROXL(N, pt, b, pbar, cm))
            nl_init_ = max(c(init_ - intdif, 1))
            nu_init_ = min(c(N, init_ + intdif))
            outfOptimn = fMinSearch(nl_init_, nu_init_)
            if (outfOptimn == nu_init_) 
                stop(" upper search interval limit reached, increase intdif parameter")
            if (outfOptimn == nl_init_) 
                stop(" lower search interval limit reached, increase intdif parameter")
            return(outfOptimn)
        }
        n = fOptimn(N, pt, b, pbar, cm)
        return(new("ACSPlan", n = n, k = kL(n, pt, b)))
    }
}
