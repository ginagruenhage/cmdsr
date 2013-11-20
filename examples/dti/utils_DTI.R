
add <- function(x) Reduce("+", x)

## helper function, output is dataframe with source.index and the according region, region index, hemisphere and lobe
get.names <- function(dt)
  {
    
    dt$nodes$hemi <- factor(str_detect(dt$nodes$dn_label,'_RH'),labels=c('LH','RH'))
    nn <- dt$nodes[,c('index','dn_freesurfer_structname','hemi')]
    
    nn$region <- paste0(nn$hemi,"_",str_split_fixed(nn$dn_free," ",n=2)[,1])
    
    names(nn) <- c("source.index","source.struct","source.hemi","source.region")
    
    lobes <- read.csv("region_scale33.csv")
    lobes$source.region <- paste0(lobes$hemisphere,"_",lobes$label)
    lobes <- subset(lobes,select=c(source.region,lobe))
    names(lobes) <- c("source.region","source.lobe")
    
    ## assign indeces to region names, 1:33 is LH, 34:66 is RH, names are alphabetical
    regs <- data.frame(source.reg.ind = 1:66, source.region= sort(unique(nn$source.region)) )
    regs <- merge(regs,lobes)
    
    ## add source.reg.ind and source.lobe to nn (common column is source.region)
    nn <- merge(nn,regs)
    nn
  }

## add names and indeces to dat.all$links
get.region.names.indeces <- function(dt)
  {
    nn <- get.names(dt)
    
    ## add nn infos to dt$links (common column is source.index)
    dt$links <- merge(dt$links,nn)

    names(nn) <- c("target.region","target.index","target.struct","target.hemi","target.reg.ind","target.lobe")

    ## add nn infos to dt$links (common column is target.index)
    dt$links <- merge(dt$links,nn)
    dt$links <- mutate(dt$links,intra.hemi= (target.hemi==source.hemi) )
    
    dt
  }

connectivity.matrix <- function(dt,stat="de_strength"){
  nN <- nrow(dt$nodes)
  A <- sparseMatrix(dt$links$source.index,dt$links$target.index,x=dt$links[,stat],dims=c(nN,nN))
  A
}

graph.group <- function(G,gvec){  
  if (class(gvec) == "factor"){
    gvec <- a.n(gvec)
  }
  
  A <- get.adjacency(G, type = "both", attr = "weight")
  N <- length(unique(gvec))
  gr <- a.m(expand.grid(1:N, 1:N))

  B <- aaply(gr, 1, function(p) mean(A[gvec == p[1], gvec == p[2]], na.rm = TRUE))

  B <- matrix(B,N,N)
  
  graph.adjacency(B, weighted = TRUE)
}

graph.average <- function(GL){
  AL <- llply(GL, get.adjacency, type = "both", attr = "weight")
  A <- add(AL) / length(AL)
  graph.adjacency(A, weighted = TRUE)
}
  
average.subjectA <- function(L)
  {
    B <- list()
    B$A <- add(L[c(1,2)])/2
    B <- c(B,L[c(3,4,5,6)])
    B
  }

plot.dti <- function(GL, XL, regions, hemi, plotname = "NULL", plot.betweenness = FALSE, plot.closeness = FALSE, ...)
  {
    
    colorpalette = c(rgb(31,119, 180,max=255), rgb( 255,127,14,max=255), rgb( 44, 160, 44,max=255),rgb(214, 39,40,max=255),rgb(148,103,189,max=255),rgb( 140,86,75,max=255),rgb(227,119,194,max=255),rgb(127,127,127,max=255),rgb( 188,189,34,max=255),rgb(23,190,207,max=255))

    cols = c(rgb(31,119, 180,max=255),rgb(148,103,189,max=255),rgb( 188,189,34,max=255), rgb( 44, 160, 44,max=255), rgb( 255,127,14,max=255),rgb(214, 39,40,max=255))

    regions <- subset(regions, source.hemi == hemi)
    nR <- dim(regions)[1]
    names(regions) <- c("label","index","hemi","lobe")

    df.label.alpha <- data.frame(time=as.factor(1:5),alpha=c(0.2,0.4,0.6,0.8,1))
    
    ed <- ldply(GL, function(G) data.frame(intergraph:::as.matrix.igraph(G,"edgelist")))
    names(ed)[1] <- "time"

    deg <- ldply(GL, function(G) data.frame(val=igraph::degree(G)))
    names(deg)[1] <- "time"
    closeness <- ldply(GL, function(G) data.frame(val=igraph::closeness(G)))
    names(closeness)[1] <- "time"
    betweenness <- ldply(GL, function(G) data.frame(val=igraph::betweenness(G)))
    names(betweenness)[1] <- "time"
    
    
    if (hemi=="RH")
      {
        coords <- adply(seq(1,length(XL)),1, function(i) data.frame(index=(nR+1):(2*nR),t(XL[[i]])))
      }
    else
      {
        coords <- adply(seq(1,length(XL)),1, function(i) data.frame(index=1:nR,t(XL[[i]])))
      }

    names(coords) <- c("time","index","X1","X2")
    
    edges <- adply(seq(1,length(XL)),1, function(t) data.frame(time=as.factor(t), X1 = subset(coords,time==t)[subset(ed,time==t)$X1,"X1"] ,Y1=subset(coords,time==t)[subset(ed,time==t)$X1,"X2"],X2=subset(coords,time==t)[subset(ed,time==t)$X2,"X1"],Y2=subset(coords,time==t)[subset(ed,time==t)$X2,"X2"]))
    
    coords <- cbind(coords, deg=deg$val, betweenness = betweenness$val, closeness = closeness$val)
    coords <- merge(coords,regions)
    
    ## short labels, for two stable, two moving regions
    df.label <- subset(coords,label %in% c("LH_superiorfrontal","LH_paracentral","LH_parahippocampal","LH_fusiform") )
    df.label <- mutate(df.label, label_short = rep(c("FUS","PARC","PARH","SF"),each=5))
    
    ## all plots in one figure
    if (plot.closeness){
      p <- ggplot(coords) + facet_wrap(~time,nrow=2,scales="fixed") + geom_segment(data=edges,aes(x=X1,y=Y1,xend=X2,yend=Y2),size=0.2,color=alpha(colorpalette[1],1/5)) + geom_point(aes(x=X1,y=X2,color=lobe,size=closeness),alpha=1) + scale_size_continuous("Closeness",range=c(1,4))
    } else if (plot.betweenness){
      p <- ggplot(coords) + facet_wrap(~time,nrow=2,scales="fixed") + geom_segment(data=edges,aes(x=X1,y=Y1,xend=X2,yend=Y2),size=0.2,color=alpha(colorpalette[1],1/5)) + geom_point(aes(x=X1,y=X2,color=lobe,size=betweenness),alpha=1)+ scale_size_continuous("Betweenness",range=c(1,4))
    } else {
      p <- ggplot(coords) + facet_wrap(~time,nrow=2,scales="fixed") + geom_segment(data=edges,aes(x=X1,y=Y1,xend=X2,yend=Y2),size=0.2,color=alpha(colorpalette[1],1/5)) + geom_point(aes(x=X1,y=X2,color=lobe,size=deg),alpha=1)+ scale_size_continuous("Degree",range=c(1,4))
    }
    
    p <- p + scale_color_manual("Lobe",values=cols)+ labs(x="Dimension 1",y="Dimension 2") + scale_x_continuous(limits=c(-2,3)) + scale_y_continuous(limits=c(-2,1.5)) 
        
    p <- p +  theme_classic() + theme(aspect.ratio=1, text = element_text(family="serif",size=7), axis.line = element_blank(),panel.border = element_rect(linetype = "solid",fill=NA, colour = "black"),legend.position="bottom",legend.key=element_rect(fill="white",colour="white"),legend.margin = unit(0,"cm"),plot.margin=unit(c(0,0,0,0),"cm"),strip.background=element_blank(),strip.text.x=element_blank())

    p <- p + geom_text(data=df.label,aes(x=(X1+0.15),y=X2,label=label_short,color=lobe),size=3,family="serif",hjust=0,vjust=0.5)
    
    p <- p + geom_text(aes(x=-2,y=-1.95,label=paste("alpha ==", alpha)),parse=TRUE,hjust=0,size=2.5,family="serif",data=df.label.alpha)

    print(p)
    
    if (hemi=="NULL")
      {
        ggsave(filename = paste("./figure/Agreement_Embedding",".pdf",sep=""),width=5)
        ggsave(filename = paste("./figure/Agreement_Embedding",".svg",sep=""),width=5)
      }
    else
      {
        if (plotname == "NULL"){
          ggsave(filename = paste("./figure/Agreement_Embedding",hemi,".pdf",sep=""),width=5,height=5.8)
          ggsave(filename = paste("./figure/Agreement_Embedding",hemi,".svg",sep=""),width=5,height=5.8)
        } else {
          ggsave(filename = paste("./figure/",plotname,hemi,".pdf",sep=""),width=5,height=5.8)
          ggsave(filename = paste("./figure/",plotname,hemi,".svg",sep=""),width=5,height=5.8)
        }
      }
    
  }
