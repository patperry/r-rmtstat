
lbe.spike.eig <- function( spike, nrow=NA, ncol=NA, var=1 ) {
    ratio  <- ncol/nrow
    above  <- spike > sqrt( ratio )*var 
    center <- ifelse( !above, NA,
                  ( spike + var )*( 1 + ratio*( var/spike ) ) )
    scale  <- ifelse( !above, NA,
                  ( ( spike + var )
                    * sqrt( 2*( 1 - ratio*( var/spike )^2 ) ) 
                    / sqrt( nrow ) ) )
                         
    
    list( center=center, scale=scale )
}

lbe.spike.vec <- function( spike, nrow=NA, ncol=NA, var=1 ) {
    ratio <- ncol/nrow
    above <- spike > sqrt( ratio )*var 
    cor2  <- ifelse( !above, NA,
                ( ( spike^2 - ratio*var^4 )
                  / ( spike*( spike + ratio*var ) ) 
                  / sqrt( nrow ) ) )
    
    cor2
}

dlbe.spike.eig <- function( x, spike, nrow=NA, ncol=NA, var=1, 
                            log = FALSE ) {
    params <- lbe.spike.eig( spike, nrow, ncol, var )
    d      <- dnorm( x, mean=params$center, sd=params$scale, log=log )
    c
}

plbe.spike.eig <- function( q, spike, nrow=NA, ncol=NA, var=1, 
                            lower.tail = TRUE, log.p = FALSE ) {
    params <- lbe.spike.eig( spike, nrow, ncol, var )
    p      <- pnorm( q, mean=params$center, sd=params$scale, lower.tail, log.p )
    p
}                                

qlbe.spike.eig <- function( p, spike, nrow=NA, ncol=NA, var=1, 
                            lower.tail = TRUE, log.p = FALSE ) {
    params <- lbe.spike.eig( spike, nrow, ncol, var )
    q      <- qnorm( q, mean=params$center, sd=params$scale, lower.tail, log.p )
    q
}                                

rlbe.spike.eig <- function( n, spike, nrow=NA, ncol=NA, var=1 ) {
    params <- lbe.spike.eig( spike, nrow, ncol, var )
    x      <- qnorm( n, mean=params$center, sd=params$scale )
    x
}
