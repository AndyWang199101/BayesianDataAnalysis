#########################################################################################
library(ggplot2)
library(latex2exp)
multiplot <- function(plots, plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    #plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}
######################################################################################


## test for binomial trails
setwd("~/Documents/Learning/BayesianDataAnalysis/chapter2/")
y = 437; n = 980 ## num. of successes; num of trails

beta_stat <- function( alpha,beta,plot=T,fig='' ){
    if(plot){
        setEPS()
        postscript(fig,width=10,height=8)
    }
    m = alpha / (alpha + beta)
    mo = (alpha-1)/ (alpha+beta-2)
    sd = sqrt( (alpha*beta)/((alpha+beta)^2*(alpha+beta+1)) )
    
    ## .95 central percentiles by Beta distribution
    upper.bound.1 <- qbeta( 0.05/2,shape1 = alpha,shape2 = beta )
    lower.bound.1 <- qbeta( 1-(0.05)/2,shape1 = alpha,shape2 = beta )
    
    ### by sampling
    posterior_samples <- rbeta(1000,shape1 = alpha,shape2 = beta)
    if( plot ){
        data.plot1 <- data.frame( x=posterior_samples )
        p <- ggplot( data=data.plot1,aes(x=x) )
        q1 <- p + geom_histogram(bins = 100) +theme_bw()+labs( x=TeX("Posterior samples of $\\theta$") )+
            theme(plot.title=element_text(size=20),
            axis.title.y=element_text(size = 25, vjust=+0.2),
            axis.title.x=element_text(size = 25, vjust=-0.2),
            axis.text.y=element_text(size = 25),
            axis.text.x=element_text(size = 25),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),legend.position = "none")
    }
    sorted_samples <- sort( posterior_samples,decreasing = F )
    upper.bound.2 <- sorted_samples[1000*0.05/2]
    lower.bound.2 <- sorted_samples[1000*(1-0.05/2)]
    
    
    ### by normal approximation
    upper.bound.3 <- qnorm( 0.05/2,mean = m,sd = sd )
    lower.bound.3 <- qnorm( 1-0.05/2,mean = m,sd = sd )
    
    ### by transform and normal approximation
    posterior_samples_transform <- log( posterior_samples/(1 - posterior_samples) )
    if( plot ){
        data.plot1 <- data.frame( x=posterior_samples_transform )
        p <- ggplot( data=data.plot1,aes(x=x) )
        q2 <- p + geom_histogram(bins = 100) +theme_bw()+labs( x=TeX("Posterior samples of $\\frac{1-\\theta}{\\theta}$") )+
            theme(plot.title=element_text(size=20),
                  axis.title.y=element_text(size = 25, vjust=+0.2),
                  axis.title.x=element_text(size = 25, vjust=-0.2),
                  axis.text.y=element_text(size = 25),
                  axis.text.x=element_text(size = 25),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),legend.position = "none")
    }
    
    ## estimate norm mean and sd
    norm_m <- mean(posterior_samples_transform)
    norm_v <- sd(posterior_samples_transform) * (1000/999)
    message(norm_m,norm_v)
    upper.bound.4 <- 1/(1+exp(-qnorm( 0.05/2,mean = norm_m,sd = norm_v )))
    lower.bound.4 <- 1/(1+exp(-qnorm( 1-0.05/2,mean = norm_m,sd = norm_v )))
    
    
    if(plot){
        multiplot( list(q1,q2) )
        dev.off()
    }
    
    return( list(mean=m,mode=mo,sd=sd,quan=c(upper.bound.1,lower.bound.1),
                 samples_quan = c(upper.bound.2,lower.bound.2),
                 norm_quan=c(upper.bound.3,lower.bound.3),
                 transform_quan = c(upper.bound.4,lower.bound.4)
                 ) 
            )
}

##########uniform prior##############
res <- beta_stat(438,544,fig = 'uniform_A.eps')
res2 <- beta_stat(13,19,fig = 'uniform_B.eps')
