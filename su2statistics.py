import numpy as np
import matplotlib.pyplot as plt

#====================================================================#
# FUNCTIONS FOR THE FITS
# function that fits to a curve like a*x1+b*x2+c
def fit_3parameters(x,y,sigmay,x1,x2):
    sumx1=sumx2=sumy=sumsigma2=sumx1x2=sumx1y=sumx2y=sumx12=sumx22=0
    n=len(x)
    if(sigmay==-1):
        sigmay=list()
        for i in range(0,n+1):
            sigmay.append(1)

            
    for i in range(0,n):
        sumx1+=x1(x[i])/(sigmay[i])**2
        sumx2+=x2(x[i])/(sigmay[i])**2
        sumy+=y[i]/(sigmay[i])**2
        
        sumx1x2+=x1(x[i])*x2(x[i])/(sigmay[i])**2
        sumx1y+=x1(x[i])*y[i]/(sigmay[i])**2
        sumx2y+=x2(x[i])*y[i]/(sigmay[i])**2

        sumx12+=(x1(x[i]))**2/(sigmay[i])**2
        sumx22+=(x2(x[i]))**2/(sigmay[i])**2
        sumsigma2+=1/(sigmay[i])**2

    summatrix=[[sumx12,sumx1x2,sumx1],[sumx1x2,sumx22,sumx2],[sumx1,sumx2,sumsigma2]]
    coeffmatrix=[[sumx1y],[sumx2y],[sumy]]

    summatrix=np.array(summatrix)
    coeffmatrix=np.array(coeffmatrix)

    inverse=np.linalg.inv(summatrix)
    coeffmatrix=np.matmul(inverse,coeffmatrix)
    return coeffmatrix[0][0],coeffmatrix[1][0],coeffmatrix[2][0]


# function that fits to a 2 parameter curve a*x1+b
def fit_2parameters(x,y,sigmay,x1):
    sumx1=sumy=sumx1y=sumx12=sumsigma2=0
    n=len(x)
    if(sigmay==-1):
        sigmay=list()
        for i in range(0,n):
            sigmay.append(1)
    
    for i in range(0,n):
        sumx1+=x1(x[i])/(sigmay[i])**2
        sumy+=y[i]/(sigmay[i])**2
        sumx1y+=x1(x[i])*y[i]/(sigmay[i])**2
        sumx12+=(x1(x[i]))**2/(sigmay[i])**2
        sumsigma2+=1/(sigmay[i])**2

    summatrix=[[sumx12,sumx1],[sumx1,sumsigma2]]
    coeffmatrix=[[sumx1y],[sumy]]

    summatrix=np.array(summatrix)
    coeffmatrix=np.array(coeffmatrix)

    summatrix=np.linalg.inv(summatrix)
    coeffmatrix=np.matmul(summatrix,coeffmatrix)

    return coeffmatrix[0][0],coeffmatrix[1][0]


# x to the minus 1 power f(x)=1/x
def xminus1(x):
    return 1/x


# x to the fist power f(x)=x
def xfirst(x):
    return x


#=============================================================+++++++#
#   STATISTICS FUNCTIONS
#=============================================================+++++++#
def mean(data):
    n=len(data)
    databar=0
    for i in range(0,n):
        databar+=data[i]
    return databar/n


# varianve of the random variable
def varianceX(data):
    n=len(data)
    databar=mean(data)
    sigma=0
    # make the mean value
    for i in range(0,n):
        sigma+=(data[i]-databar)**2
    return sigma/(n-1)


# variance of the mean estimator
def varianceXhat(data):
    n=len(data)
    databar=mean(data)
    sigma=0
    # make the mean value
    for i in range(0,n):
        sigma+=(data[i]-databar)**2
    return sigma/(n*(n-1))


# function that given a correlation number, make the bining of a set of data
def bining(filename):
    file=open(filename,'r')
    data=list()
    sigmak=list()
    lenk=list()
    # read the data from the file
    for line in file.readlines():
        l=[float(x) for x in line.split()]
        data.append(l[1])
    file.close()

    # number of correlated points
    nc=len(data)

    # now we construct the blocks
    for k in range(2,nc-1):
        if(nc%k==0):
            blockeddata=[]
            # construct the blocks
            block=0
            nk=0
            for i in range(0,nc-nc%k):
                block+=data[i]
                nk+=1
                if(nk==k): # we have a block
                    blockeddata.append(block/k)
                    block=0
                    nk=0
            # we have the data on blocks, that the standart deviation
            if(len(blockeddata)!=1): 
                sigmak.append(np.sqrt(varianceXhat(blockeddata)))
                lenk.append(k)
    
    # now we try to fit a a curve a/k+b to this set
    # the value kcorr for with the data is good, is our nof correlation
    file=open('bining.dat','w')
    for ns in range(0,len(lenk)-1):
        file.write(' '.join([str(x) for x in [lenk[ns],sigmak[ns],'\n']]))        
        
    file.close()
    plt.scatter(lenk,sigmak)
    plt.show()


# function that given a set of uncorrelated date, compute the error by the jackkinfe
def jackkinfe(data):
    thetahat=mean(data)
    thetan=list()
    n=len(data)
    sthetatilde=0
    # construct the sets if the n-th element removed
    for i in range(0,n):
        l=[]
        l=data[::]
        del l[i]
        thetan.append(mean(l))

    for i in range(1,n):
        sthetatilde+=(thetan[i]-thetahat)**2
    return (n-1)*sthetatilde/n


# jackknife for a plot data
def jackkinfe_ploteddata(datax,datay,np,x1,x2):
    thetan1=list()
    thetan2=list()
    n=len(datax)
    sthetatilde1=0
    sthetatilde2=0

    if np==2: # 2 parameter fit
        thetahat1,thetahat2=fit_2parameters(datax,datay,-1,x1)
        
    if np==3: # 2 parameter fit
        thetahat1,thetahat2, thetahat3=fit_3parameters(datax,datay,-1,x1,x2)
        thetan3=list()
        sthetatilde3=0

    # construct the sets if the n-th element removed
    for i in range(0,n):
        x=[]
        y=[]
        x=datax[::]
        y=datay[::]
        del x[i]
        del y[i]
        
        if np==2:
            a,b=fit_2parameters(x,y,-1,x1)
            thetan1.append(a)
            thetan2.append(b)
        if np==3:
            a,b,c=fit_3parameters(x,y,-1,x1,x2)
            thetan1.append(a)
            thetan2.append(b)
            thetan3.append(c)

    # compute the variations
    for i in range(1,n):
        sthetatilde1+=(thetan1[i]-thetahat1)**2
        sthetatilde2+=(thetan2[i]-thetahat2)**2
        if(np==3): sthetatilde3+=(thetan3[i]-thetahat3)**2

    sthetatilde1*=(n-1)/n
    sthetatilde2*=(n-1)/n
    if(np==3): sthetatilde3*=(n-1)/n

    if(np==2): return sthetatilde1, sthetatilde2
    if(np==3): return sthetatilde1, sthetatilde2, sthetatilde3

