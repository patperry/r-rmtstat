
wishart.max.par <- (function() {
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
    
    function( n.df, p.dim, var=1, beta=1 ) {
        n <- n.df
        p <- p.dim
        
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

dwishart.max <- function( x, n.df, p.dim, var=1, beta=1, log = FALSE ) {
    params <- wishart.max.par( n.df, p.dim, var, beta )
    x.tw   <- ( x - params$center )/( params$scale )
    d.tw   <- dtw( x, beta, log )
    d      <- ifelse( log, d.tw - log( params$scale ),
                  d.tw / ( params$scale ) )
    d
}

pwishart.max <- function( q, n.df, p.dim, var=1, beta=1, 
                      lower.tail = TRUE, log.p = FALSE ) {
    params <- wishart.max.par( n.df, p.dim, var, beta )
    q.tw   <- ( q - params$center )/( params$scale )
    p      <- ptw( q.tw, beta, lower.tail, log.p )
    p
}

qwishart.max <- function( p, n.df, p.dim, var=1, beta=1,
                      lower.tail = TRUE, log.p = FALSE ) {
    params <- wishart.max.par( n.df, p.dim, var, beta )
    q.tw   <- qtw( p, beta, lower.tail, log.p )
    q      <- params$center + q.tw*( params$scale )
    q
}

rwishart.max <- function( n, n.df, p.dim, var=1, beta=1 ) {
    params <- wishart.max.par( n.df, p.dim, var, beta )
    x.tw   <- rtw( n, beta )
    x      <- params$center + x.tw*( params$scale )
    x
}
