##########load packages
library(conicfit)
library(gpclib)
library(plotrix)
##################set font for graphs
windowsFonts(A = windowsFont("Arial"))

####read in data file that contains the parameters used in each run
##this is the file that can be modified according to species-specific and 
#target characteristics
pf<-read.csv("parameters_to_modify.csv")

#####from here onward the script should NOT be modified
###################################################################################################

######read data file that contains the information about the eye model
of<-read.csv("eye_model_design_do_not_modify.csv")
str(of)
attach(of)


######################load functions######
######area of overlap of two circles,radii r1 and r2, distance between centres d#####
ict2<-function(r1,r2,d)
{
  rs<-c(r1,r2)
  if(rs[2]>rs[1])
  {r1<-rs[2]
  r2<-rs[1]}
  if(d<=(r1-r2))
  {pi*r2^2}
  else
    if(d>=(r1+r2))
    {0}
  else
  {d1<-(r1^2-r2^2+d^2)/(2*d)
  d2<-d-d1
  r1^2*acos(d1/r1)-d1*sqrt(r1^2-d1^2)+r2^2*acos(d2/r2)-d2*sqrt(r2^2-d2^2)}
}
##################
##############################
####area of overlap between and ellipse and a circle
#####the ellipse stays with horizontal long axis, and is centred at (0,0)
#a and b are the long and short half dimensions
###the circle is centred at (x,y) and has radius "radius"
ict3<-function(a,b,radius,x,y)
{
  circ<-calculateCircle(x,y,radius)
  elipse<-calculateEllipse(0,0,a,b)
  cc<-as(circ,"gpc.poly")
  ee<-as(elipse,"gpc.poly")
  ii<-intersect(cc,ee)
  area.poly(ii)/area.poly(cc)}


##set parameters: target diameter (size), altitude (alt), and horizontal distance (dist) in cm
size<-pf[3,2]
dist<-pf[4,2]
alt<-pf[5,2]
################eye data in degrees
#######inter-ommatidial angle (interomm) and receptive field (rf)
interomm<-pf[1,2]
rf<-pf[2,2]
if(!is.na(pf[5,2]>0)){


#############################################################################
##### PLOTS FOR HORIZONTAL TARGETS ############################################### 
  
#(dist is the horizontal distance from the centre of the target)
angle<-atan(alt/dist)*57.296/10
dist2<-sqrt(dist^2+alt^2)
near<-dist-size/2
far<-dist+size/2
sdn<-sqrt(near^2+alt^2)
sdf<-sqrt(far^2+alt^2)
suba<-2*atan((size/2)/dist2)*180/pi
subb<-acos((sdn^2+sdf^2-size^2)/(2*sdn*sdf))*180/pi
a<-suba/20
b<-subb/20

###### PLOT OF OMMATIDIA WITH TARGET SUPERIMPOSED #######
dev.new(noRStudioGD = TRUE, width=10.4, height=6, unit="in")
ppp<-par(xpd=NA,mai=c(.8,.8,.3,1.5),pin=c(4,4))
plot(x,y,xlim=c(-6.5,6.5),ylim=c(-6.5,6.5),xaxt="n",yaxt="n",xlab="Degrees off axis (horizontal)",
ylab="Degrees off axis (vertical)",family="A")
for (i in 1:37)
{
overlap<-as.integer((ict3(a,b,rf/20,x[i],y[i])*100))
of$olap[i]<-100-overlap
cl<-sprintf("grey%i",100-overlap)
draw.ellipse(x[i],y[i],1,1,col=cl)
}
draw.ellipse(0,0,a,b,border="red")
lbls<-c(-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60)
axis(1,at=lbls/10,labels=lbls,las=3,family="A")
axis(2,at=lbls/10,labels=lbls,las=1,family="A")
####add legend
for(i in 1:100)
{
cll<-sprintf("grey%i",i)
rect(7.5,-2+i*.04,8.5,-2+i*.04+.04,col=cll,border=NA)
}
leg<-c("100","75","50","25","0")
ys<-c(-2,-1,0,1,2)
for (i in 1:5)
{
text(10,ys[i]+.1,leg[i],pos=2,cex=0.9,family="A")
arrows(8.5,ys[i]+.1,8.7,ys[i]+.1,length=0)
}
text(9,3,"% Stimulation",cex=0.9,family="A")


############  POLAR EYE PLOT ###############
dev.new(noRStudioGD = TRUE, width=6, height=6, unit="in")
ppp<-par(pin=c(4,4))
plot(x/3,y/3-angle,xlim=c(-5,5),ylim=c(-9,0),xaxt="n",yaxt="n",
xlab=paste("Distance from target",dist,"cm"),ylab="",type="n",family="A")
for (i in 1:37)
{draw.ellipse(0,0,i,i,border="black")
}
angs<-seq(0,180,20)
ngr<-angs/57.296
for (i in 1:10)
{
lines(c(0,9*cos(ngr[i])),c(0,-9*sin(ngr[i])),col="black")
}
for (j in 1:9)
{
text(j*cos(ngr[4]),-j*sin(ngr[4]),j*10,col="black",family="A")
}
draw.ellipse(0,-angle,a,b,col=adjustcolor("red",alpha=.2))
par(ppp)

############ MAXIMUM DETECTION DISTANCE #################################
###set individual ommatidia threshold (tom) and total eye threshold (teye)
tom<-pf[6,2]
teye<-pf[7,2]

#################################################################

ict3<-function(a,b,radius,x,y)
{
  circ<-calculateCircle(x,y,radius)
  elipse<-calculateEllipse(0,0,a,b)
  cc<-as(circ,"gpc.poly")
  ee<-as(elipse,"gpc.poly")
  ii<-intersect(cc,ee)
  area.poly(ii)/area.poly(cc)}

################################################################
stimulus<-function(dist)
{
  sto<-vector()
  angle<-atan(alt/dist)*57.296/10
  dist2<-sqrt(dist^2+alt^2)
  near<-dist-size/2
  far<-dist+size/2
  sdn<-sqrt(near^2+alt^2)
  sdf<-sqrt(far^2+alt^2)
  suba<-2*atan((size/2)/dist2)*180/pi
  subb<-acos((sdn^2+sdf^2-size^2)/(2*sdn*sdf))*180/pi
  a<-suba/20
  b<-subb/20
  for (i in 1:37)
  {
    sto[i]<-as.integer((ict3(a,b,rf/20,x[i],y[i])*100))
  }
  sum(sto[sto>tom])}
################## moves target in from 500 cm ################################
lop<-vector()
for (i in 1:500)
{
  lop[i]<-stimulus(i)
}
############## PLOT STIMULUS VS DISTANCE ########################################

#######finds threshold distance
td<-max(which(lop>=teye))
print(paste("threshold detection distance for your target =",td, "cm"))

######### FOR TARGETS FROM 20 cm to 100 cm ###################################

ttt<-vector()
for (k in 1:9)
{
  size<-k*10+10
  lop<-vector()
  for (i in 1:500)
  {
    lop[i]<-stimulus(i)
  }
  td<-max(which(lop>=100))
  ttt[k]<-td
  ttt
}

####### PLOT DETECTION DISTANCE VS TARGET SIZE ###############
targsize<-c(20,30,40,50,60,70,80,90,100)
dev.new(noRStudioGD = TRUE, width=6, height=6, unit="in")
ppp<-par(pin=c(4,4))
plot(targsize,ttt,xlab="Target size",ylab="Detection distance", family="A", ylim=c(0,350))
m1<-lm(ttt~targsize)
abline(m1)
############finds coefficients for regression equation
cc<-coef(m1)
############ PRINTS EQUATION
cc1<-round(cc[1],4)
cc2<- abs(cc1)
if (cc1>0) {sgn<-"+"} else {sgn<-"-"}
print (paste("Equation for detection distance = ", "(", round(cc[2],4),"*target size)",
             sgn, cc2, sep = ""))

############################################################################
}else {
  
############################################################################
############# CALCULATIONS FOR VERTICAL TARGETS  #########################
  
subt<-atan(size/dist)*180/pi
rd<-subt/20

dev.new(noRStudioGD = TRUE, width=10.4, height=6, unit="in")
ppp<-par(xpd=NA,mai=c(.8,.8,.3,1.5),pin=c(4,4))
plot(x,y,xlim=c(-6.5,6.5),ylim=c(-6.5,6.5),xaxt="n",yaxt="n",xlab="Degrees off axis (horizontal)",
ylab="Degrees off axis (vertical)",type="n",family="A")
for (i in 1:37)
{
overlap<-as.integer((ict2(rf/20,rd,r[i]*interomm/10)/(((rf/20)^2)*pi))*100)
of$olap[i]<-100-overlap
cl<-sprintf("grey%i",100-overlap)
draw.ellipse(x[i],y[i],1,1,col=cl)
}
draw.ellipse(0,0,rd,rd,border="red")
lbls<-c(-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60)
axis(1,at=lbls/10,labels=lbls,las=3,family="A",cex=5)
axis(2,at=lbls/10,labels=lbls,las=1,family="A",cex=5)
#####add legend
for(i in 1:100)
{
cll<-sprintf("grey%i",i)
rect(7.5,-2+i*.04,8.5,-2+i*.04+.04,col=cll,border=NA)
}
leg<-c("100","75","50","25","0")
ys<-c(-2,-1,0,1,2)
for (i in 1:5)
{
text(10,ys[i]+.1,leg[i],pos=2,cex=0.9,family="A")
arrows(8.5,ys[i]+.1,8.7,ys[i]+.1,length=0)
}
text(9,3,"% Stimulation",cex=0.9, family="A")

par(ppp)

###############################################
############ MAXIMUM DETECTION DISTANCE #################################

###set individual ommatidia threshold (tom) and total eye threshold (teye)
tom<-pf[6,2]
teye<-pf[7,2]

######area of overlap of two circles,radii r1 and r2, distance between centres #####
ict2<-function(r1,r2,d)
{
  rs<-c(r1,r2)
  if(rs[2]>rs[1])
  {r1<-rs[2]
  r2<-rs[1]}
  if(d<=(r1-r2))
  {pi*r2^2}
  else
    if(d>=(r1+r2))
    {0}
  else
  {d1<-(r1^2-r2^2+d^2)/(2*d)
  d2<-d-d1
  r1^2*acos(d1/r1)-d1*sqrt(r1^2-d1^2)+r2^2*acos(d2/r2)-d2*sqrt(r2^2-d2^2)}
}

################################################################
stimulus<-function(dist)
{
  sto<-vector()
  subt<-atan(size/dist)*180/pi
  rd<-subt/20
  for (i in 1:37)
  {
    sto[i]<-as.integer((ict2(rf/20,rd,r[i]*interomm/10)/(((rf/20)^2)*pi))*100)
  }
  sum(sto[sto>tom])}
################## moves target in from 500 cm #################################
lop<-vector()
for (i in 1:500)
{
  lop[i]<-stimulus(i)
}


#######finds threshold distance
td<-max(which(lop>=teye))
print(paste("threshold detection distance for your target =",td, "cm"))

########### FOR TARGETS FROM 10cm to 100cm###########
ttt<-vector()
for (k in 1:10)
{
  size<-k*10
  lop<-vector()
  for (i in 1:500)
  {
    lop[i]<-stimulus(i)
  }
  td<-max(which(lop>=100))
  ttt[k]<-td
  ttt
}
ttt

####### PLOT DETECTION DISTANCE VS TARGET SIZE ###############
targsize<-c(10,20,30,40,50,60,70,80,90,100)
dev.new(noRStudioGD = TRUE, width=6, height=6, unit="in")
ppp<-par(pin=c(4,4))
plot(targsize,ttt,xlab="Target size",ylab="Detection distance", family="A", ylim=c(0,350))
m1<-lm(ttt~targsize)
abline(m1)
############finds coefficients for regression equation
cc<-coef(m1)
cc1<-round(cc[1],4)
cc2<- abs(cc1)
if (cc1>0) {sgn<-"+"} else {sgn<-"-"}
print (paste("Equation for detection distance = ", "(", round(cc[2],4),"*target size)",
            sgn, cc2, sep = ""))

}

