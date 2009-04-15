
wishart.spike.par <- function( spike, n.df=NA, p.dim=NA, var=1 ) {
    ratio  <- p.dim/n.df
    above  <- spike > sqrt( ratio )*var 
    center <- ifelse( !above, NA,
                  ( spike + var )*( 1 + ratio*( var/spike ) ) )
    scale  <- ifelse( !above, NA,
                  ( ( spike + var )
                    * sqrt( 2*( 1 - ratio*( var/spike )^2 ) ) 
                    / sqrt( n.df ) ) )
    
    list( center=center, scale=scale )
}

wishart.spike.vec.par <- function( spike, n.df=NA, p.dim=NA, var=1 ) {
    ratio <- p.dim/n.df
    above <- spike > sqrt( ratio )*var 
    cor2  <- ifelse( !above, NA,
                ( ( spike^2 - ratio*var^4 )
                  / ( spike*( spike + ratio*var ) ) 
                  / sqrt( n.df ) ) )
    
    cor2
}

dwishart.spike <- function( x, spike, n.df=NA, p.dim=NA, var=1, 
                            log = FALSE ) {
    params <- wishart.spike.par( spike, n.df, p.dim, var )
    d      <- dnorm( x, mean=params$center, sd=params$scale, log=log )
    c
}

pwishart.spike <- function( q, spike, n.df=NA, p.dim=NA, var=1, 
                            lower.tail = TRUE, log.p = FALSE ) {
    params <- wishart.spike.par( spike, n.df, p.dim, var )
    p      <- pnorm( q, mean=params$center, sd=params$scale, lower.tail, log.p )
    p
}                                

qwishart.spike <- function( p, spike, n.df=NA, p.dim=NA, var=1, 
                            lower.tail = TRUE, log.p = FALSE ) {
    params <- wishart.spike.par( spike, n.df, p.dim, var )
    q      <- qnorm( q, mean=params$center, sd=params$scale, lower.tail, log.p )
    q
}                                

rwishart.spike <- function( n, spike, n.df=NA, p.dim=NA, var=1 ) {
    params <- wishart.spike.par( spike, n.df, p.dim, var )
    x      <- qnorm( n, mean=params$center, sd=params$scale )
    x
}
