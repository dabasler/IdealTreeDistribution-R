
# ITD Model
# DW Purves et al., 'Crown plasticity and competition for canopy space: a new spatially implicit model parameterised for 250 North American tree species'.

#R code by D.Basler 2017

### Canopy
#EXPECTED INPUT 
# 1. dbh, species
# 2. plot area
#
# Processing:
# 1. Match Species
# 2. Calcuate Height
# 3. Calculate Rmax
# 4. Add crown parameters (ratio,curvature,Vbias) to dataset...
# Redy for optimization
#
#    1. select initial zq=10
#    2. cacalulate atot(z):
#       if (a<1)decrease zq
#       else increase zq
#       stop at improvement<0.01
#       z=zq
#    
#    3.Calculate: 
#       Crown Radius
#       Crown shape
#       ... Further stuff to DO



import_ITD_parameters<-function(filename){
    ### GET THE ITS PRECALCULATED PARAMETERS
    ITD_par<-read.table(filename,sep=",",skip=1,header=TRUE,na.strings = "NA",stringsAsFactors = FALSE)
    ITD_par$species<-paste(ITD_par$genus,ITD_par$species_epithet)
    # Data from tabl S2 (doi:10.1371/journal.pone.0000870.s005)
    S2<-data.frame(C0=c(0.503,0.5,0.701,0.196,0.95,2.551),C1=c(3.126,10.0,3.955,0.511,0.95,4.106),fixed=c(0,1,0,0,1,0))
    row.names(S2)<-c("R0","R40","Rus","B","M","Vus")
    #ITD_par$Tj                                                             # Trait score (0-1) of species j (used in the single axis parameter estimation scheme).
    #ITD_par$Vbias                                                          # Distance above Z* (m) of the base of the crown of trees of species j.
    ITD_par$R0  = (1-ITD_par$Tj)*S2["R0","C0"]  +ITD_par$Tj*S2["R0","C1"]   # Maximum potential crown radius (m) of a tree with dbh 0 cm
    ITD_par$R40 = (1-ITD_par$Tj)*S2["R40","C0"] +ITD_par$Tj*S2["R40","C1"]  # Maximum potential crown radius (m) of a tree with dbh 40 cm
    ITD_par$Rus = (1-ITD_par$Tj)*S2["Rus","C0"] +ITD_par$Tj*S2["Rus","C1"]  # Crown radius (m) of understory tree.
    ITD_par$B   = (1-ITD_par$Tj)*S2["B","C0"]   +ITD_par$Tj*S2["B","C1"]    # crown curvature
    ITD_par$M   = (1-ITD_par$Tj)*S2["M","C0"]   +ITD_par$Tj*S2["M","C1"]    # Crown ratio at which the maximum radius is realized.
    ITD_par$Vus = (1-ITD_par$Tj)*S2["Vus","C0"] +ITD_par$Tj*S2["Vus","C1"]  # Crown depth (m) of understory tree
    ITD_par$id=c(1:nrow(ITD_par))
    return (ITD_par)
}


# Basic calculation formulas
calc_height<-function(dbh,a,b) return (10^(a+b*log10(dbh)))       #tree height from dbh
calc_rmax <- function(R0,R40,dbh) return(R0 +(R40-R0)*(dbh/40))   #tree potential crown diameter from dbh

# Crown Shape formulas: 
calc_r_fromtop     <- function(y,height,rmax,M,B)  return (rmax*(pmin(y,height*M)/ (height*M))^B)
calc_r_aboveground <- function(y,height,rmax,M,B)  return (rmax*(pmin((height-y),height*M)/ (height*M))^B)
calc_height_of_r   <- function(r,height,rmax,M,B)  return (height-((r/Rmax)^(1/B))*(height*M))


# ITD_Model:# Calculate estimated crown area (this function is used for optimisation of Z)
calc_az<-function(zq,height,Rmax,M,B,V)
{
    z_tilde=pmax(pmax(zq+V,0),zq)
    r=rep(0,length(height))
    rc=Rmax*(pmin(height-z_tilde, height*M)/(height*M))^B
    r[height>zq]<-rc[height>zq]
    a=pi*r^2
    return(sum(a)) # canopy area at zq
}


# Extracts specific parameters from ITD dataset
prepare_ITD<-function(dbh,sp,ITD_par)
{ 
  height<-calc_height(dbh,ITD_par$a_dbh[sp],ITD_par$b[sp])
  rmax<-calc_rmax(ITD_par$R0[sp],ITD_par$R40[sp],dbh) 
  M<-ITD_par$M[sp]
  B<-ITD_par$B[sp]
  V<-ITD_par$Vbias[sp]
  Rus<-ITD_par$Rus[sp]
  Vus<-ITD_par$Vus[sp]
  out<-data.frame(height=height,rmax=rmax,M=M,B=B,V=V,Rus=Rus,Vus=Vus)
  return(out)
}

# Optimized Z (canopy height), in order to match plot area
optimize_z<-function(a_plot,height,Rmax,M,B,V)
{
  if (calc_az(0,height,Rmax,M,B,V)<a_plot){print("ERROR: plot contains gaps, ITD not possible");return (0)}
  zq=0 # large enough starting value
  step<-100
  while ((2*step)>0.001) {
    #print (step)
    ta=calc_az(zq,height,Rmax,M,B,V)
    if (ta<a_plot){zq<-zq-step}
    if (ta>a_plot){zq<-zq+step}
    step=step/2
  }
  return (zq)
}

# Calculates individual tree parameters (canopy status, crown radius, exposed crown, crown density)
calc_results<-function(z,height,Rmax,M,B,V,Rus,Vus){
  z_hat<-pmax(z+V,0) # effective canopy height
  z_tilde<-pmax(z_hat,z)
  canopy_status<-height>z
  r_crown<-Rmax*(pmin(height-z_hat, height*M)/(height*M))^B
  r_exposed<-Rmax*(pmin(height-z_tilde, height*M)/(height*M))^B
  crown_depth<-pmin(pmax(height-z,0),height)
  # set understory parameters
  r_crown[canopy_status==0]<-Rus[canopy_status==0]
  r_exposed[canopy_status==0]<-0
  crown_depth[canopy_status==0]<-pmin(Vus,height)[canopy_status==0]
  df<-data.frame(canopy_status=canopy_status,r_crown=r_crown,r_exposed=r_exposed,crown_depth=crown_depth) 
  return(df)
}

#match species 
match_species<-function(species,ITD_par){
  sp<-rep(NA,length(species))
  for (s in unique(species)) sp[species==s]<-match(s,ITD_par$species)
  if (NA %in% sp) {print ("ERROR: Could not match all species, please check manually"); return (NULL)} else  {return (sp)} 
}



# ITD<-function(dbh,species,a_plot,ITD_par){
#   ds<-dataframe(dbh=dbh,sp=spid)# spid row number of speciesin ITD tab
#   ds<-cbind(ds,prepare_ITD(ds$dbh,ds$sp,ITD_par))
#   z<-optimize_z(a_plot,ds$height,ds$rmax,ds$M,ds$B,ds$V)
#   if (z==0) return(NULL)
#   ds<-cbind(ds,calc_results(z,ds$height,ds$rmax,ds$M,ds$B,ds$V,ds$Rus,ds$Vus))
#   return(list(z,ds))
# }
# 



ITD<-function(ds,a_plot,ITD_par){
  sp<-match_species(ds$species,ITD_par)
  if (!is.null(sp)){
    ds$sp<-sp  
    ds<-cbind(ds,prepare_ITD(ds$dbh,ds$sp,ITD_par))
    z<-optimize_z(a_plot,ds$height,ds$rmax,ds$M,ds$B,ds$V)
    if (z==0) return(NULL)
    ds<-cbind(ds,calc_results(z,ds$height,ds$rmax,ds$M,ds$B,ds$V,ds$Rus,ds$Vus))
    # Final result:
    return(list(z,ds))
  }
  return(NULL)
}



plot_tree<-function(result,i){
  z<-result[[1]]
  ds<-result[[2]]
  pl<-ceiling(max(ds$r_crown*2,max(ds$height)))
  plot(NA,xlim=c(-pl,pl),ylim=c(0,pl),asp=1,main=sprintf("Tree%s (%s)",i,ds$species[i]))
  segments(ds$dbh [i]/100,0,ds$dbh [i]/100,ds$height [i]-ds$crown_depth [i])
  segments(-ds$dbh [i]/100,0,-ds$dbh [i]/100,ds$height [i]-ds$crown_depth [i])
  ch<-c(seq(ds$height [i]-ds$crown_depth [i],ds$height [i],0.1),ds$height [i])
  rz<-calc_r_aboveground(ch,ds$height[i], ds$rmax[i],ds$M [i],ds$B [i])
  lines(rz,ch)
  lines(-rz,ch)
  
  ch<-seq(z,z+ds$V[i],0.1*ds$V[i]/abs(ds$V[i]))
  rz<-calc_r_aboveground(ch,ds$height [i],ds$rmax[i],ds$M [i],ds$B [i])
  lines(rz,ch)
  lines(-rz,ch)
  
  
  col<-"green"
  if (ds$canopy_status[i]==0) col<-"red"
  lines(c(-ds$r_crown[i],ds$r_crown[i]),c(0,0),lty=1,lwd=3,col=col)
  
  lines(c(-ds$r_exposed[i],ds$r_exposed[i]),c(z,z),lty=1,lwd=3)
  
  lines(c(-100,100),c(z,z),lty=2)
  lines(c(-100,100),c(0,0))
  
}


plot_trees<-function(result){
 
   z<-result[[1]]
  ds<-result[[2]]

  pl<-ceiling(max(ds$r_crown*2,max(ds$height)))
  plot(NA,xlim=c(0,sum(ds$r_exposed*2)),ylim=c(0,pl))
  treepos<-cumsum(ds$r_exposed*2)-ds$r_exposed
  understory<-which(ds$canopy_status==0)
  treepos[understory]<-sample( cumsum(ds$r_exposed[ds$r_exposed>0]*2)[-length(which(ds$canopy_status==1))] ,length(understory))
  for (i in 1:nrow(ds)){  
      segments(treepos[i]+ds$dbh [i]/100,0,treepos[i]+ds$dbh [i]/100,ds$height [i]-ds$crown_depth [i])
      segments(treepos[i]-ds$dbh [i]/100,0,treepos[i]-ds$dbh [i]/100,ds$height [i]-ds$crown_depth [i])
      ch<-c(seq(ds$height [i]-ds$crown_depth [i],ds$height [i],0.1),ds$height [i])
      rz<-calc_r_aboveground(ch,ds$height[i], ds$rmax[i],ds$M [i],ds$B [i])
      lines(treepos[i]+rz,ch)
      lines(treepos[i]-rz,ch)
      
      ch<-seq(z,z+ds$V[i],0.1*ds$V[i]/abs(ds$V[i]))
      rz<-calc_r_aboveground(ch,ds$height [i],ds$rmax[i],ds$M [i],ds$B [i])
      lines(treepos[i]+rz,ch)
      lines(treepos[i]-rz,ch)
      
      
      if (ds$canopy_status[i]==0) {col<-"red"
      lines(c(treepos[i]-ds$r_crown[i],treepos[i]+ds$r_crown[i]),c(0.2,0.2),lty=1,lwd=3,col=col)
      }else{
        col<-"green"
        lines(c(treepos[i]-ds$r_crown[i],treepos[i]+ds$r_crown[i]),c(0,0),lty=1,lwd=3,col=col)
        
        
      }
      
      lines(c(treepos[i]-ds$r_exposed[i],treepos[i]+ds$r_exposed[i]),c(z,z),lty=1,lwd=3)
      
      lines(c(-100,100),c(z,z),lty=2)
      lines(c(-100,100),c(0,0))
  }
}



plot_treelineup<-function(result){
  
  z<-result[[1]]
  ds<-result[[2]]
  
  pl<-ceiling(max(ds$r_crown*2,max(ds$height)))
  plot(NA,xlim=c(0,sum(ds$r_crown*2)),ylim=c(0,pl))
  treepos<-cumsum(ds$r_crown*2)-ds$r_crown
  for (i in 1:nrow(ds)){  
    segments(treepos[i]+ds$dbh [i]/100,0,treepos[i]+ds$dbh [i]/100,ds$height [i]-ds$crown_depth [i])
    segments(treepos[i]-ds$dbh [i]/100,0,treepos[i]-ds$dbh [i]/100,ds$height [i]-ds$crown_depth [i])
    ch<-c(seq(ds$height [i]-ds$crown_depth [i],ds$height [i],0.1),ds$height [i])
    rz<-calc_r_aboveground(ch,ds$height[i], ds$rmax[i],ds$M [i],ds$B [i])
    lines(treepos[i]+rz,ch)
    lines(treepos[i]-rz,ch)
    
    ch<-seq(z,z+ds$V[i],0.1*ds$V[i]/abs(ds$V[i]))
    rz<-calc_r_aboveground(ch,ds$height [i],ds$rmax[i],ds$M [i],ds$B [i])
    lines(treepos[i]+rz,ch)
    lines(treepos[i]-rz,ch)
    
    
    if (ds$canopy_status[i]==0) {col<-"red"}else{
      col<-"green"}
      lines(c(treepos[i]-ds$r_crown[i],treepos[i]+ds$r_crown[i]),c(0,0),lty=1,lwd=3,col=col)

    lines(c(treepos[i]-ds$r_exposed[i],treepos[i]+ds$r_exposed[i]),c(z,z),lty=1,lwd=3)
    
    lines(c(-100,100),c(z,z),lty=2)
    lines(c(-100,100),c(0,0))
  }
}


#### 
