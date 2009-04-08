
lbe.max <- (function() {
    mu <- function( n,p ) {
        n.sqrt <- sqrt( n )
        p.sqrt <- sqrt( p )
        res    <- ( n.sqrt + p.sqrt )^2
        res
    }
    
    sigma <- function( n,p ) {
        n.sqrt <- sqrt( n )
        p.sqrt <- sqrt( p )
        res    <- ( n.sqrt + p.sqrt )*( 1/n.sqrt + 1/p.sqrt )^( 1/3 )
        res
    }

    mu.real <- function( n,p ) {
        mu( n-1/2,p-1/2 )
    }
    
    sigma.real <- function( n,p ) {
        sigma( n-1/2,p-1/2 )
    }
    
    alpha <- function( n,p ) {
        1/( 1 + ( mu( n-1/2,p+1/2 )/mu( n+1/2,p-1/2 ) )
                * sqrt( sigma( n-1/2,p+1/2 )/sigma( n+1/2,p-1/2 ) ) )
    }
    
    mu.cplx <- function( n,p ) {
        a   <- alpha( n,p )
        res <- mu( n-1/2,p+1/2 )*a + mu( n+1/2,p-1/2 )*( 1-a )
        res
    }
    
    sigma.cplx <- function( n,p ) {
        a   <- alpha( n,p )
        res <- sigma( n-1/2,p+1/2 )*a + sigma( n+1/2,p-1/2 )*( 1-a )
        res
    }
    
    function( nrow, ncol, var=1, beta=1 ) {
        n <- nrow
        p <- ncol
        
        if( beta == 1 ) {
            m <- mu.real( n,p )
            s <- sigma.real( n,p )
        } else if( beta == 2 ) {
            m <- mu.cplx( n,p )
            s <- sigma.cplx( n,p )
        } else {
            stop( "`beta' must be 1 or 2, not `", beta, "'")
        }
            
        center <- var*( m/n )
        scale  <- var*( s/n )
    
        list( center=center, scale=scale )
    }
})()

dlbe.max <- function( x, nrow, ncol, var=1, beta=1, log = FALSE ) {
    params <- lbe.max( nrow, ncol, var, beta )
    x.tw   <- ( x - params$center )/( params$scale )
    d.tw   <- dtw( x, beta, log )
    d      <- ifelse( log, d.tw - log( params$scale ),
                  d.tw / ( params$scale ) )
    d
}

plbe.max <- function( q, nrow, ncol, var=1, beta=1, 
                      lower.tail = TRUE, log.p = FALSE ) {
    params <- lbe.max( nrow, ncol, var, beta )
    q.tw   <- ( q - params$center )/( params$scale )
    p      <- ptw( q.tw, beta, lower.tail, log.p )
    p
}

qlbe.max <- function( p, nrow, ncol, var=1, beta=1,
                      lower.tail = TRUE, log.p = FALSE ) {
    params <- lbe.max( nrow, ncol, var, beta )
    q.tw   <- qtw( p, beta, lower.tail, log.p )
    q      <- params$center + q.tw*( params$scale )
    q
}

rlbe.max <- function( n, nrow, ncol, var=1, beta=1 ) {
    params <- lbe.max( nrow, ncol, var, beta )
    x.tw   <- rtw( n, beta )
    x      <- params$center + x.tw*( params$scale )
    x
}
