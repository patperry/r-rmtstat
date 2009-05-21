
dmp <- function( x, ndf=NA, pdim=NA, var=1, svr=ndf/pdim, log = FALSE ) {
    gamma          <- svr
    
    inv.gamma.sqrt <- sqrt( 1/gamma )
    a              <- var*( 1 - inv.gamma.sqrt )^2
    b              <- var*( 1 + inv.gamma.sqrt )^2
    
    if( !log ) {
        # we have to handle +/- zero carefully when gamma=1
        d <- ifelse( gamma == 1 & x == 0 & 1/x > 0, Inf, 
                 ifelse( x <= a | x >= b, 0,
                     suppressWarnings(
                         gamma/( 2*pi*var*x )
                         *
                         sqrt( ( x-a )*( b-x ) ) 
                         ) ) )
    } else {
        d <- ifelse( gamma == 1 & x == 0 & 1/x > 0, Inf,
                 ifelse( x <= a | x >= b, -Inf,
                     suppressWarnings(
                         log( gamma ) 
                         - 
                         ( log( 2 ) + log( pi ) + log( var ) + log( x ) )
                         +
                         0.5*log( x-a ) 
                         + 
                         0.5*log( b-x )
                         ) ) )
    }
    d
}
