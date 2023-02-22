import numpy as np

def generation(n):
    #Génération d'un processus caché composé de n états
    #0 -> Sain et 1 -> Fievreux
    n=n+1
    prbX=np.random.choice(10,n)
    initial=prbX[0]
    etatsX=[]
    etatsXnum=[]
    if initial<6:
        etatsX.append("Sain")
        etatsXnum.append(0)
    else:
        etatsX.append("Fievreux")
        etatsXnum.append(1)        
    for i in range(1,n-1):
        if (etatsX[i-1]=="Sain"):
            if prbX[i]<7:
                etatsX.append("Sain")
                etatsXnum.append(0)
            else:
                etatsX.append("Fievreux")
                etatsXnum.append(1)                 
        else:
            if prbX[i]<6:
                etatsX.append("Sain")
                etatsXnum.append(0)
                    
            else:
                etatsX.append("Fievreux")
                etatsXnum.append(1)
    #Génération du processus observable composé de n états
    #0-> Etourdi, 1 -> Rhume, 2->Normal
    prbY=np.random.choice(10,n)
    etatsY=[]
    etatsYnum=[]
    for i in range(n-1):
        if (etatsX[i]=="Sain"):
            if prbY[i]<1:
                etatsY.append("Etourdi")
                etatsYnum.append(0)
            elif prbY[i]<5:
                etatsY.append("Rhume")
                etatsYnum.append(1)
            else:
                etatsY.append("Normal")
                etatsYnum.append(2)
        else: 
            if prbY[i]<1:
                etatsY.append("Normal")
                etatsYnum.append(2)
            elif prbY[i]<4:
                etatsY.append("Rhume")
                etatsYnum.append(1)
            else:
                etatsY.append("Etourdi")
                etatsYnum.append(0)
    return(prbX,etatsX,etatsXnum,prbY,etatsY,etatsYnum)

generation(5)

Q=np.array([[0.7,0.3],[0.4,0.6]])
psi=np.array([[0.1,0.4,0.5],[0.6,0.3,0.1]])
rho=[0.6,0.4]

k=5

def algoForward(k,gene):
    #On a seulement besoin des états observables Y de 0 à k
    etatsYnum=gene[5][0:k+1]
    #on a besoin d'une matrice de (k+1) lignes car on doit calculer les probabilités de 0 à k
    alpha=np.eye(k+1,2)
    def algoForwardRec(k,alpha):
        #n -> nombre d'états de la chaîne
        #k -> étape de l'algo (0<k<n)        
        if (k>0):
            tmp=algoForwardRec(k-1,alpha)[k-1]
            alpha[k,0]=round(psi[0,etatsYnum[k]]*(Q[0,0]+Q[1,0])*tmp[0],5)
            alpha[k,1]=round(psi[1,etatsYnum[k]]*(Q[0,1]+Q[1,1])*tmp[1],5)
            return(alpha)
        else:     
            #lorque k=0
            alpha[0,0]=round(rho[0]*psi[0,etatsYnum[0]],5)
            alpha[0,1]=round(rho[1]*psi[1,etatsYnum[0]],5)
            return(alpha)        
    alpha=algoForwardRec(k,alpha)
    return(alpha)
    
def filtrage(k,gene):
    # retourne l'état caché xk, à partir des états observables de l'instant 0 à l'instant k
    alpha=algoForward(k,gene)
    #initialisation de la matrice
    filtrage=alpha
    denom=np.eye(k,2)
    print(filtrage)
    for i in range (k):
        denom[i,]=filtrage[i,0]+filtrage[i,1]
    filtrage=filtrage/denom   
    return (filtrage)



def algoBackward(l,k,gene):
    #On veut connaître les probabilités des états Y de l+1 à k
    #On a seulement besoin des états observables Y de 0 à k
    etatsYnum=gene[5][0:k+1] 
    delta=k-l
    beta=np.eye(k-l,2)    
    def algoBackwardRec(l,beta,k):
        #n -> nombre d'états de la chaîne
        if (l<k-1):
            print(l)
            # on récupère la ligne d'au dessus
            tmp=algoBackwardRec(l+1,beta,k)[(l-1),]
            for x in range(2):
                #xprime = 0
                b0=Q[x,0]*psi[0,etatsYnum[l+1]]*tmp[0]
                #xprime = 1
                b1=Q[x,1]*psi[1,etatsYnum[l+1]]*tmp[1]
                beta[l-2,x]=b0+b1
            print(l,beta[l-2,])
            return(beta)        
        else:
            #lorsque l=k-1
            print(l,k)
            beta[delta-1,]=[1,1]
            return(beta)
    beta=algoBackwardRec(l,beta,k)
    return(beta)

 #faire lissage  


def viterbi(n,k,l):
    #Initialisation
    gene=generation(n)
    #On a seulement besoin des états observables Y de 0 à k
    etatsYnum=gene[5][0:k+1] 
    rho=[0.6,0.4]
    w=np.eye(k+1,2)
    
    def ProgDynaRec(n,k,w):
        #permet de calculer w(0:k)
        #n -> nombre d'états de la chaîne
        #k -> étape de l'algo (0<k<n-1)        
        if (k>0):
            tmp=ProgDynaRec(n,k-1,w)[k-1]
            for x in range(2):
                m=max(Q[0,x]*tmp[0],(Q[1,x]*tmp[1]))
                w[k,x]=m*psi[x,etatsYnum[k+1]]
            return(w)
        else:       
            w[0,0]=round(rho[0]*psi[0,etatsYnum[0]],5)
            w[0,1]=round(rho[1]*psi[1,etatsYnum[0]],5)
            return(w)        
    w=ProgDynaRec(n,k,w)
    xchapeauk=np.eye(k+1,1)
    for i in range(k+1):
        xchapeauk[i]=np.argmax(w[i,])
    xchapeaul=np.eye(k+1,k+1)
    def viterbiRec(n,l,xchapeaul):
        if (l<k):
            tmp=viterbiRec(n,l+1,xchapeaul)[:,l+1]
            #on obtient la matrice de taille k*k, et on récupère la colonne l+1
            for x in range(2):
                #à revoir
                m=max(Q[x,tmp]*w[l,x])
                xchapeaul[:,l]=m
            return(xchapeaul)
        else:
            #on s'arrête quand on a l=k
            xchapeaul[:,k]=xchapeauk
            return(xchapeaul)
    xchapeaul=viterbiRec(n,l,xchapeaul)
    print(xchapeaul)
    return(w)    

