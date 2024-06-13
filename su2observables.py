import numpy as np
import matplotlib.pyplot as plt
import su2statistics
import os
#=============================================================+++++++#
# function that, given a number of correlation, put the data of a file in blocks
# the output is a file with name unc'name of the file'
# with the blocked, uncorrelated data
def block_data(filein,ncorr):
    data=list()
    file=open(filein,'r')
    # read the data from a file
    n=0
    ff=0
    for line in file.readlines():
        ff+=1
        if(ff==5001): break
        l=[float(x) for x in line.split()]

        # we will put the data of the same column in the same array
        if (n==0):
            for i in range(1,len(l)):
                # we take the element from that column and append on data
                col=[]
                col.append(l[i])
                data.append(col)
        else:
            # we append the data in the respect column of data array
            for i in range(0,len(l)-1):
                data[i].append(l[i+1])
        n+=1
    file.close()
    
    # we have read the data, now we make blocks with if
    databloked=list()
    for i in range(0,len(data)): # we block each data column
        block=0
        blockcol=[]
        nblock=0
        for j in range(0,n-n%ncorr): # we block this column
            block+=data[i][j]
            nblock+=1
            if(nblock==ncorr):
                blockcol.append(block/ncorr)
                block=0
                nblock=0
        databloked.append(blockcol)
    data=[]
    # now we make the mean and the standar deviation of each new set of uncorrelated points
    for i in range(0,len(databloked)):
        data.append([su2statistics.mean(databloked[i]),np.sqrt(su2statistics.varianceXhat(databloked[i]))])
    
    
    # we get out with the uncorrelated means
    file1=open("unc"+filein,'w')
    file2=open("mean"+filein,'w')
    #for i in range(0,len(databloked)):
    #    file1.write(' '.join([str(x) for x in databloked[i]]))
    #    file1.write('\n')

    for i in range(0,len(data)):
        file2.write(' '.join([str(x) for x in data[i]]))

        file2.write('\n')
    file1.close()
    file2.close()


# function that block the data for the wilson loops
def blockwilson(a,p,ncorr):
    # make the mean value of the wilson loops
    for i in range(1,a+1):
        # open the file
        if(p=='full'): file=f"fort.{100+i}" # full loops
        if(p=='proj'): file=f"fort.{200+i}" # vortex projected loops  
        if(p=='rem'): file=f"fort.{300+i}" # vortex remov loops
        if(p=='even'): file=f"fort.{500+i}" # vortex remov loops
        if(p=='odd'): file=f"fort.{600+i}" # vortex remov loops
        block_data(file,ncorr)
        #os.cmd(f'rm unc')

# function that computes the potential
def pot(a,b,p,nt,afm,amev):    
    w=list()
    sw=list()
    v=list()
    sv=list()
    t=list()
    r=list()
    for i in range(0,b):
        t.append(i+1)

    if(p=='full'): fileout=open(f"potqq.dat",'w') # full loops
    if(p=='proj'): fileout=open(f"potqq-proj.dat",'w') # vortex projected loops  
    if(p=='rem'): fileout=open(f"potqq-rem.dat",'w') # vortex remov loops

    # read data from file meanfort.10+a
    for i in range(1,a+1):
        w=[]
        sw=[]
        r.append(i)

        # open the file
        if(p=='full'): file=open(f"meanfort.{100+i}",'r') # full loops
        if(p=='proj'): file=open(f"meanfort.{200+i}",'r') # vortex projected loops  
        if(p=='rem'): file=open(f"meanfort.{300+i}",'r') # vortex remov loops

        # read the expetation values for the wilson loop of size rxt, with t=1,...b]]
        for line in file.readlines():
            l=[]
            l=[float(x) for x in line.split()]
            w.append(np.log(abs(l[0])))
            sw.append(abs(l[1])/l[0])
        file.close()
        # having the values we must make a linear regression to comput the potential
        # fit w=aexp(-v t)=:log(w)=log(a)-v*t
        av, bv =su2statistics.fit_2parameters(t[nt::],w[nt::],sw[nt::],su2statistics.xfirst)
        da, db=su2statistics.jackkinfe_ploteddata(t[nt::],w[nt::],2,su2statistics.xfirst,su2statistics.xfirst)
        v.append(-av)
        sv.append(np.sqrt(da))
    
        # get out with the data
        fileout.write(' '.join([str(x) for x in [i*afm,-av/(amev*1e3),da/(amev*1e3)]]))
        fileout.write('\n')
    #lt.errorbar(r,v,yerr=sv,label=p,fmt='D')

    fileout.close()
    if p=='full':
        # fit the curve v(x)=sigma*r+k/r+v0 to the potential
        sigma, k, v0=su2statistics.fit_3parameters(r,v,sv,su2statistics.xfirst,su2statistics.xminus1)
        dsigma, dk, dv0=su2statistics.jackkinfe_ploteddata(r,v,3,su2statistics.xfirst,su2statistics.xminus1)
        return sigma, dsigma, k, dk, v0, dv0
    if p=='proj':
        # fit the curve v(x)=sigma*r+k/r+v0 to the potential
        sigma,  v0=su2statistics.fit_2parameters(r,v,sv,su2statistics.xfirst)
        dsigma, dv0=su2statistics.jackkinfe_ploteddata(r,v,2,su2statistics.xfirst,su2statistics.xminus1)
        return sigma, dsigma, v0, dv0
    if p=='rem':
        # fit the curve v(x)=sigma*r+k/r+v0 to the potential
        sigma,  v0=su2statistics.fit_2parameters(r,v,sv,su2statistics.xminus1)
        dsigma, dv0=su2statistics.jackkinfe_ploteddata(r,v,2,su2statistics.xminus1,-1)
        return sigma, dsigma, v0, dv0


# function that computes the creutz rations
def creutz(a,p):
    chi=list()
    dchi=list()
    w=list()
    sw=list()
    for i in range(1,a+1):
        wi=[]
        swi=[]
        # open the file
        if(p=='full'): file=open(f"meanfort.{100+i}",'r') # full loops
        if(p=='proj'): file=open(f"meanfort.{200+i}",'r') # vortex projected loops  
        if(p=='rem'): file=open(f"meanfort.{300+i}",'r') # vortex remov loops

        # read the data
        wi=[]
        swi=[]
        for line in file.readlines():
            l=[float(x) for x in line.split()]
            wi.append(l[0])
            swi.append(l[1])
        w.append(wi)
        sw.append(swi)
        file.close() 
    # once we have the data, we compute the rations
    for i in range(0,a):
        if i==0:
            chi.append(-np.log(w[i][i]))
            dchi.append(sw[i][i]/w[i][i])
        else:
            chi.append(-np.log((w[i][i]*w[i-1][i-1])/(w[i-1][i]*w[i][i-1])))
            dchi.append(sw[i][i]/w[i][i]+sw[i-1][i-1]/w[i-1][i-1]+sw[i-1][i]/w[i-1][i]+sw[i][i-1]/w[i][i-1])

    # get out with the data
    fileout=open(p+'creutz.dat','w')
    for i in range(0,a):
        fileout.write(' '.join([str(x) for x in [i+1,chi[i],dchi[i],'\n']]))
    fileout.close()


def wilsoncomponens(a,b):
    even=list()
    odd=list()
    full=list()
    area=list()
    for i in range(1,a+1):
        for j in range(1,b+1):
            file=open(f'fort.{500+i*j}','r')
            eveni=[]
            oddi=[]
            fulli=[]
            g=0
            for line in file.readlines():
                g+=1
                if(g==5001):break
                l=[float(x) for x in line.split()]
                fulli.append(l[1])
                eveni.append(l[2])
                oddi.append(l[3])
            file.close()
            full.append(np.mean(fulli))
            even.append(np.mean(eveni))
            odd.append(np.mean(oddi))
            area.append(i*j)

    file=open('wilson-vortexcomp.dat','w')
    for i in range(0,len(area)):
        file.write(' '.join([str(x) for x in [area[i],full[i],even[i],odd[i],'\n']]))
    file.close()

# function that sepate
# function that deals with the projected-vortices (p-vortices) density
def pvorticesprob(a,b):
    even=list()
    odd=list()
    area=list()
    for i in range(1,a+1):
        for j in range(1,b+1):
            file=open(f'fort.{400+i*j}','r')
            eveni=[]
            oddi=[]
            a=0
            for line in file.readlines():
                a+=1
                if(a==5001): break
                l=[float(x) for x in line.split()]
                eveni.append(l[1])
                oddi.append(l[2])
            file.close()
            even.append(np.mean(eveni))
            odd.append(np.mean(oddi))
            area.append(i*j)

    file=open('pvortex-probability.dat','w')
    for i in range(0,len(area)):
        file.write(' '.join([str(x) for x in [area[i],even[i],odd[i],'\n']]))
    file.close()


def vortexlimitedwilson(a,b):
    w0=list()
    w1=list()
    w2=list()
    area=list()
    for i in range(1,a+1):
        for j in range(1,b+1):
            file=open(f'fort.{600+i*j}','r')
            w0i=[]
            w1i=[]
            w2i=[]
            f=0
            for line in file.readlines():
                f+=1
                if(f==5001): break
                l=[float(x) for x in line.split()]
                w0i.append(l[1])
                w1i.append(l[2])
                w2i.append(l[3])
            file.close()
            w0.append(np.mean(w0i))
            w1.append(np.mean(w1i))
            w2.append(np.mean(w2i))
            area.append(i*j)

    file=open('voterxlimied-wilson.dat','w')
    for i in range(0,len(area)):
        file.write(' '.join([str(x) for x in [area[i],w0[i],w1[i],w2[i],'\n']]))
    file.close()


def gluon(nr,infile,outfile,ap):
    file=open(infile,'r')
    dd=list()
    d=list()
    derr=list()
    p=list()
    #==================================================
    for line in file.readlines():
        l=[]
        l=[float(x) for x in line.split()]
        
        if(len(dd)==0):
            for i in range(0,len(l)):
                ll=[]
                ll.append(l[i])
                dd.append(ll)
                p.append(2*np.sin(np.pi*(i)/nr))
        else:
            for i in range(0,len(l)):
                dd[i].append(l[i])
    #==================================================

    for i in range(0,len(dd)):
        d.append(su2statistics.mean(dd[i]))
        derr.append(su2statistics.varianceXhat(dd[i]))

    file.close()
    file=open(outfile,'w')
    for i in range(0,len(d)):
        file.write(' '.join([str(x) for x in [p[i]/(ap*1e3),d[i]*(ap*1e3)**2,derr[i],'\n']]))
    file.close()
#=============================================================================================================#