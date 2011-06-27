Ims <-
function (n, N, pt, pbar, cm = 1, b = 0.1) 
n * cm + (N - n) * (1 - OC(pbar, n, kL(n, pt, b), type = "t"))
