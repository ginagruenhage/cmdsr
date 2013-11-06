
##' Display cMDS results using the googleVis motion widget
##'
##' Display a 2D embedding as a motion chart. See the "Personality over time" vignette for more information
##' @title googleVis.cmds
##' @param res: cMDS output (with embedding dimension k=2)
##' @param df: a data.frame giving information on the embedded points
##' @param time.var: name of the "time" variable in the data.frame, corresponding to the different distance matrices
##' @param id.var: name of the variable identifying the different points
##' @return a googleVis chart. Use plot(googleVis.cmds(res)) to view in browser. 
##' @author Simon Barthelmé
##' @export
googleVis.cmds <- function(res,df,time.var="time",id.var="id")
    {

        if (nrow(res$XL[[1]])!=2)
            {
                stop("Embedding must be two-dimensional")
            }
        nT <- length(res$DL)
        ids <- rownames(res$DL[[1]])
        if (length(names(res$DL))>0)
            {
                timeval <- as.numeric(names(res$DL))
            }
        else
            {
                timeval <- 1:nT
            }
        embed <- ldply(1:nT,function(ind)
                       {
                           df <- as.data.frame(t(res$XL[[ind]]))
                           names(df) <- c("cmds.x1","cmds.x2")
                           df[,time.var] <- timeval[ind]
                           df[,id.var] <- ids
                           df
                       })
        
        if (!missing(df))
            {
                embed <- merge(embed,df)
            }
        gvisMotionChart(embed,id=id.var,time=time.var)
    }
