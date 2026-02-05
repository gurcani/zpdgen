from numpy import reshape,shape,transpose
from inmzpd import inmweid,gmweid
from epszpd import epsweid,sigweid
rank = lambda x : len(shape(x))

def Inm(za,zb,b,n,m):
    res=inmweid(za,zb,b,n,m)
    if (shape(za) ==()):
        res=res[0]
        return res
    res=reshape(res,shape(transpose(za)))
    res=transpose(res)
    return res

def epsitg(om,pars):
    res=epsweid(om,pars)
    if (shape(om) ==()):
        res=res[0]
        return res
    res=reshape(res,shape(transpose(om)))
    res=transpose(res)
    return res

def sigmazpd(za,zb,b,anm):
    if (rank(anm)==2):
        numn,numm=shape(anm)
    elif (rank(anm)==1):
        numn=len(anm)
        numm=0
    elif (rank(anm)==0):
        numn=0
        numm=0
    else:
        print("error the rank of anm is not 0,1 or 2")
        return 0
    res=sigweid(za,zb,b,anm,numn,numm)
    if (shape(za) ==()):
        res=res[0]
        return res
    res=reshape(res,shape(transpose(za)))
    res=transpose(res)
    return res

def Gm(z1,z2,m):
    if(shape(z1)!=shape(z2)):
        print("error: shape(z1)!=shape(z2)!\n");
        res=()
        return res;
    res=gmweid(z1,z2,m)
    if (shape(z1) ==() and shape(z2)==()):
        res=res[0]
        return res
    res=reshape(res,shape(transpose(z1)))
    res=transpose(res)
    return res
