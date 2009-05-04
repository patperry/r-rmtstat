
dmp <- function( x, n.df=NA, p.dim=NA, var=1, concentration=n.df/p.dim, log = FALSE ) {
    gamma          <- concentration
    inv.gamma.sqrt <- sqrt( 1/gamma )
    a              <- var*( 1 - inv.gamma.sqrt )^2
    b              <- var*( 1 + inv.gamma.sqrt )^2
    
    if( !log ) {
        d <- ifelse( x <= a || x >= b, 0,
                gamma/( 2*pi*var*x )
                *
                sqrt( ( x-a )*( b-x ) ) 
                )
    } else {
        d <- ifelse( x <= a || x >= b, -Inf,
                 log( gamma ) 
                 - 
                 ( log( 2 ) + log( pi ) + log( var ) + log( x ) )
                 +
                 0.5*log( x-a ) 
                 + 
                 0.5*log( b-x )
                 )
    }
    d
}
