
# JIT provides significant speedup of the simpleDLA function
# require(compiler)
# enableJIT(3)

require(dqrng)

simpleDLA <- function(nsam, tiny=1e-9, seed=NULL) {
  # compute a diffusion limited aggregation of discs of radius 1/2

  # 
  set.seed(seed)
  
  Hugetol <- 1.02
  
  v <- matrix(NaN, nsam, 1)
  p <- matrix(NaN, nsam, 1)
  rho <- matrix(0, nsam, 1)

  v[1] <- 0;

  # we standardize the second disc to be centered at 1:
  v[2] <- 1; p[2] <- 0; rho[2] <- 1;
  
  for (k in 3:nsam) {
    
    Huge <- 2 + rho[k-1] * Hugetol
    
    x <- exp(21 * pi * dqrunif(1)) * Huge
    
    y <- Inf; z <- Inf
    
    while (TRUE) {
      if (abs(x) > Huge * Hugetol) {
        q <- exp(2i * pi * dqrunif(1))
        x <- Huge * x * (1 - Huge * q / Conj(x)) / (Huge - x * q)
      }
      
      d <- abs(x - v[1:k-1])
      iy <- which.min(d); dy <- d[iy]; y <- v[iy]
      
      if (dy < 1 + tiny) {
        break
      }
      
      rad <- dy - 1
      x <- x + rad * exp(2i * pi * dqrunif(1))
      
    }
    
    v[k] <- x; p[k] <- y; rho[k] <- max(abs(x), rho[k-1])
    
    if (k %% 5000 == 0) {
      plot(1.06 * Huge * c(1,-1), 1.06 * Huge * c(1,-1), type="n")
      # paradoxically if you eliminate this line the function is slower?
      symbols(c(0,0), c(0,0), Huge*c(1,Hugetol), inches = FALSE, fg="cyan", add=TRUE)
      segments(Re(v[2:k]), Im(v[2:k]), Re(p[2:k]), Im(p[2:k]),
               lwd=2, col = topo.colors(k))
      title(sub=toString(k))
    }
    
  }
  
  list(v=v, p=p)
  
}
