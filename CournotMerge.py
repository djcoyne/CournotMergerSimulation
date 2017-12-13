"""
Updated: 12/13/17

@author: David Coyne

============================================================

This code is a modified version of Dan Greenfield & Nathan Wilson's work.
It generates predicted effects of a merger/acquisition by leveraging the following
assumptions:

a. Cournot competition
b. Linear demand, potentially in multiple markets
c. quadratic costs

Each application requires the following inputs:
1. Matrix of production quantities/capacities at the firm level
2. Market demand elasticity
3. Pre-merger price
4. Any merger-specific marginal cost efficiency
"""
import pandas as pd
import numpy as np
from numpy.linalg import solve
from IPython.display import display
import scipy.integrate as integrate
from scipy.optimize import minimize

# Define a merger object class
class cournotMerge:
    def __init__(self, F={}, markets=[], e=[], p=[], MF = [0,1], E=0):
        # Reformat the input data for calculations
        u = MF[0]
        v = MF[1]
        self.markets=markets
        ns = len(F)
        self.ns=ns
        ms = len(markets)
        self.ms=ms
        df = pd.DataFrame(columns =["Firm"]+F[u].mnames)

        i = 0
        for index,f in enumerate(F):
            df.loc[i]=[F[f].name,*F[f].q]
            i+=1
        self.df = df
            
        # Some simple calculations
        sums = df.sum(axis=0)
        Q = np.zeros(ms)
        for i in range(0,ms):
            Q[i]=sums[i+1]
        self.Q = Q
        shares=df.iloc[:,1:]/Q
        margins = shares/e
        mc = ((1-margins)*p)
        t = np.array(mc)-np.array(mc.min(axis=1))[:,None]
        k = df.sum(axis=1, numeric_only=True)/mc.min(axis=1)
        b = (np.ones(ms))/e*(p/Q)
        a = p+b*Q
        self.a = a
        self.b = b
        
        # Solve for profits
        rev = np.matmul(np.matrix(p),df.iloc[:,1:].T)
        costs = (df.iloc[:,1:].sum(axis=1))**2/(2*k).T + (t*df.iloc[:,1:]).sum(axis=1)
        prof = np.diag(rev - np.matrix(costs).T)
        self.prof = prof
        
        # Calculate market HHI values
        hhi = np.zeros(ms)
        for i in range(0,ms):
            ss = []
            for x in shares.iloc[:,i]:
                ss.append(np.square(x))
            hhi[i] = 10000*(sum(ss))
        self.hhi = hhi
        
        # Run the merger simulation
        post_F = {}
        y = list(F.keys()).index(u)
        z = list(F.keys()).index(v)
        self.mp = [F[u].name,F[v].name]
        post_F[1] = F[u].merge(F[v])
        NM = [x for i,x in enumerate(F) if x!=u and x!=v]
    
        for i in range(0,ns-2): 
            l = NM[i]
            post_F[i+2] = F[l]
        post_ns = len(post_F)
        self.post_ns = post_ns

        # Define k_post and t_post (weighted average transport cost)
        k_post = np.zeros(post_ns)
        t_post = np.zeros((post_ns,ms))
        k_post[0] = k[y]+k[z]
        t_post[0,:] = np.array(t[y]*F[u].q/post_F[1].q) + np.array(t[z]*F[v].q/post_F[1].q)
        for i in range(0,ns-2): 
            l = NM[i]
            li = list(F.keys()).index(l)
            k_post[i+1] = k[li]
            t_post[i+1] = t[li]
    
        #Set up the system of equations by market
        M_part={}
        for i in range(0,len(b)):
            M_part[i*(len(b)+1)] = b[i]*np.ones((post_ns,post_ns))-b[i]*np.diag(np.ones(post_ns))+((np.ones(post_ns)+2*b[i]*k_post)/k_post)*np.diag(np.ones(post_ns))
        l = 0
        for i in range(0,len(b)):
            for j in range(1,len(b)): 
                if l == i*(len(b)+1):
                    l+=1
                M_part[l]=(np.ones(post_ns)/k_post)*np.diag(np.ones(post_ns))
                l+=1
    
        M = np.zeros([1,ms*post_ns+1])
        for i in range(0,ms):
            Mr = np.zeros([post_ns,1])
            for j in range(0,ms):
                Mr = np.hstack((Mr,M_part[i*ms+j]))
            M = np.vstack((M,Mr))
        M = M[1:,1:]
        V = np.zeros(post_ns*ms)
        Vpart = np.vstack(np.hsplit(a-t_post,ms))
        for i in range(0,post_ns*ms):
            V[i]=Vpart[i][0]
    
        # Account for efficiencies
        if E!=0:
            for i in range(0,ms):
                for j in range(0,ms):
                    M[i*post_ns,j*post_ns] += (-E/100)/k_post[0] 

        # Solve the model
        q_post = np.linalg.solve(M,V)
        pos = np.where(q_post>=0)[0]
        if len(pos) < len(V):
            print("*WARNING*: Interior solution infeasible. Program will report a corner solution... User should check other corners\n\n")
        
        while len(pos) < len(V):
            Mp = M[pos][:,pos]
            Vp = V[pos]
            q_postp = np.linalg.solve(Mp, Vp)
            q_post = np.zeros(len(V))
            q_post[pos] = q_postp
            pos = np.where(q_post>=0.)[0]
        
        # Fill a dataframe with post-merger quantities
        df_post = pd.DataFrame(columns =["Firm"]+post_F[1].mnames)
        i=0
        for i,f in enumerate(post_F):
            postprod=[]
            for j in range(0,ms):
                l = i+j*post_ns
                postprod.append(q_post[l])
            df_post.loc[i] = [post_F[f].name, *postprod]
            i+=1
        self.df_post = df_post

        # Calculate aggregate quantities
        sums = df_post.sum(axis=0)    
        Q_post = np.zeros(ms)
        for i in range(0,ms):
            Q_post[i]=sums[i+1]
        self.Q_post = Q_post

        Q_change = 100*(Q_post-Q)/Q
        self.Q_change = Q_change
        shares_post = df_post.iloc[:,1:]/Q_post

        # Calculate post-merger prices
        p_post = a - b*Q_post
        self.p_post = p_post
        p_change = 100*(p_post-p)/p
        self.p_change = p_change
        
        # Solve for post-merger profits
        rev_post = np.matmul(np.matrix(p_post),df_post.iloc[:,1:].T)
        costs_post = (df_post.iloc[:,1:].sum(axis=1))**2/(2*k_post).T + (t_post*df_post.iloc[:,1:]).sum(axis=1)
        prof_post = np.diag(rev_post - np.matrix(costs_post).T)
        self.prof_post = prof_post
        prof_change = np.zeros(post_ns)
        prof_change[0] = 100*(prof_post[0] - (prof[y]+prof[z]))/(prof[y]+prof[z]) 
        for i in range(0,ns-2): 
            l = NM[i]
            li = list(F.keys()).index(l)
            prof_change[i+1] = 100*(prof_post[i+1] - prof[li])/prof[li]
        self.prof_change = prof_change
        
        # Calculate market HHI values
        hhi_post = np.zeros(ms)
        for i in range(0,ms):
            ss = []
            for x in shares_post.iloc[:,i]:
                ss.append(np.square(x))
            hhi_post[i] = 10000*(sum(ss))  
        self.hhi_post = hhi_post
        hhi_change = hhi_post - hhi
        self.hhi_change = hhi_change
        
        # Calculate the harm to consumers
        CS = np.zeros(ms)
        CS_post = np.zeros(ms)
        for i in range(0,ms):
            def integrand(x):
                return a[i]-b[i]*x-p[i]
            CS[i], err = integrate.quad(integrand,0,Q[i])
        for i in range(0,ms):
            def integrand(x):
                return a[i]-b[i]*x-p_post[i]
            CS_post[i], err = integrate.quad(integrand,0,Q_post[i])
    
        harm = CS - CS_post
        self.harm = harm
    
    def describe(self):
        # Describe the merger object
        ds = "This object simulates a merger between "+self.mp[0]+" and "+self.mp[1]+".\n"
        if self.ns>2:
            ds = ds + "\nAdditional competitors include:\n   "
            for i in range(1,self.post_ns):
                ds = ds+"%s\n   " %(self.df_post.iloc[i,0])
        ds = ds + "\nCompetitors compete in "+str(self.ms)+" markets:\n   "
        for i in range(0,self.ms):
            ds = ds+"%s\n   " %(self.markets[i])
        return ds
    
    def summarize(self):
        # Display inverse demand functions
        print("Demand Functions\n===================================")
        for i in range(0,self.ms):
            print("Inverse demand in %s is given by: P = %.3f - %.3f * Q" %(self.markets[i],self.a[i],self.b[i]))
        print("\n")
    
        # Display current production
        print("Pre-merger Production\n===================================")
        display(self.df)
        for i in range(0,self.ms):
            print("Total production in %s: %.2f" % (self.markets[i],self.Q[i]))
        print("\n")

        # Display pre-merger variable profits
        print("Pre-merger Profits\n===================================")
        for i in range(0,self.ns):
            print("Variable profits for %s: %.2f" %(self.df.iloc[i,0],self.prof[i]))
        print("\n")    

        # Display pre-merger HHI values
        print("Pre-merger Concentration\n===================================")
        for i in range(0,self.ms):
            print("HHI for %s is: %.0f" %(self.markets[i],self.hhi[i]))
        print("\n")  
    
        # Display the post-merger prices
        print("Post-merger Prices\n===================================")
        for i in range(0,self.ms):
            print("%s: %.2f (a change of %.2f percent)" % (self.markets[i],self.p_post[i],self.p_change[i]))
        print("\n")

        # Display the post-merger production
        print("Post-merger Production\n===================================")
        display(self.df_post)
        for i in range(0,self.ms):
            print("Total production in %s: %.2f (a change of %.2f percent)" % (self.markets[i],self.Q_post[i],self.Q_change[i]))
        print("\n")
   
        # Display post-merger variable profits
        print("Post-merger Profits\n===================================")
        for i in range(0,self.post_ns):
            print("Variable profits for %s: %.2f (a change of %.2f percent)" %(self.df_post.iloc[i,0],self.prof_post[i],self.prof_change[i]))
        print("\n")
   
        # Display post-merger HHI values
        print("Post-merger Concentration\n===================================")
        for i in range(0,self.ms):
            print("Post-merger HHI for %s is: %.0f (a change of %.0f)" %(self.markets[i],self.hhi_post[i], self.hhi_change[i]))
        print("\n")

        # Report the estimated harm to consumers
        print("Consumer Harm\n===================================")
        for i in range(0,self.ms):
            print("Estimated harm to consumers in %s: %.3f" % (self.markets[i], self.harm[i]))
        print("\n")
        
    def production(self):
        # Display the post-merger production
        ps = ""
        for i in range(0,self.ms):
            ps = ps+"\nTotal production in %s: %.2f (a change of %.2f percent)" % (self.markets[i],self.Q_post[i],self.Q_change[i])
        return ps
        
    def prices(self):
        # Display the post-merger prices
        ps = ""
        for i in range(0,self.ms):
            ps = ps+"\nPrice in %s: %.2f (a change of %.2f percent)" % (self.markets[i],self.p_post[i],self.p_change[i])
        return ps     
        
    def concentration(self):
        # Display post-merger HHI values
        ps = ""
        for i in range(0,self.ms):
            ps = ps +"\nPost-merger HHI for %s is: %.0f (a change of %.0f)" %(self.markets[i],self.hhi_post[i], self.hhi_change[i])
        return ps
        
    def profits(self):
        # Display post-merger variable profits
        ps = ""
        for i in range(0,self.post_ns):
            ps = ps+"\nVariable profits for %s: %.2f (a change of %.2f percent)" %(self.df_post.iloc[i,0],self.prof_post[i],self.prof_change[i])
        return ps
        
    def demand(self):
        # Display inverse demand functions
        ps = ""
        for i in range(0,self.ms):
            ps = ps + "\nInverse demand in %s is given by: P = %.3f - %.3f * Q" %(self.markets[i],self.a[i],self.b[i])
        return ps
        
    def harm(self):
        # Report the estimated harm to consumers
        ps = ""
        for i in range(0,self.ms):
            ps = ps + "\nEstimated harm to consumers in %s: %.3f" % (self.markets[i], self.harm[i])
        return ps

# Define a Cournot firm object class "firm"
class firm:
    def __init__(self, name="", q=[0], markets=[""]):
        self.name = name
        self.q = np.array(q)
        self.mnames = markets
    
    def describe(self):
        ds = "Production for %s \n===================================" % (self.name)
        for q,name in zip(self.q, self.mnames):
            ds = ds + "\n %.2f in %s" % (q, name)
        return ds
    
    def merge(self, partner):
        name = self.name + "/" + partner.name
        q = np.array(self.q) + np.array(partner.q)
        mf = firm(name,q,self.mnames)
        return mf
