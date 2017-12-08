"""
Created on Thu Nov 30 09:52:54 2017

@author: dcoyne

============================================================

This code is a modified version of Dan Greenfield & Nathan Wilson's work.
It generates predicted effects of a merge by leveraging the following
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

# Define a firm object class
class firm:
    def __init__(self, name="", q=[0], mnames=[""]):
        self.name = name
        self.q = np.array(q)
        self.mnames = mnames
    
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

# Define the markets (use the same order for elasticities, prices, and firm quantities below)
markets = ["North America", "Europe", "Other"]

# Define the market elasticities
e = [0.35, 0.75, 1.1]

# Define the market prices
p = [1.05, 0.95, 0.94]

# Define the Firms (Name, quantities)
F={}
F[1] = firm("Firm 1", [3.6, 3.5, 2.1], markets)
F[2] = firm("Firm 2", [2.8, 3.2, 2.6], markets)
F[3] = firm("Firm 3", [5.2, 4.4, 3.1], markets)
F[4] = firm("Firm 4", [0.6, 2.5, 1.6], markets)

# Identify the indices for the two parties that are merging
MF = [1,2]

# Identify any MC efficiencies for the merging parties (E% reduction in MC)
E = 10

# Reformat the input data for calculations
ns = len(F)
ms = len(markets)
df = pd.DataFrame(columns =["Firm"]+F[1].mnames)

i = 0
for index,f in enumerate(F):
    df.loc[i]=[F[f].name,*F[f].q]
    i+=1

# Some simple calculations
sums = df.sum(axis=0)
Q = np.zeros(ms)
for i in range(0,ms):
    Q[i]=sums[i+1]
shares=df.iloc[:,1:]/Q
margins = shares/e
mc = ((1-margins)*p)
t = np.array(mc)-np.array(mc.min(axis=1))[:,None]
k = df.sum(axis=1, numeric_only=True)/mc.min(axis=1)
b = (np.ones(ms))/e*(p/Q)
a = p+b*Q

# Display inverse demand functions
print("Demand Functions\n===================================")
for i in range(0,ms):
    print("Inverse demand in %s is given by: P = %.2f - %.2f * Q" %(markets[i],a[i],b[i]))
print("\n")

# Display current production
print("Pre-merger Production\n===================================")
display(df)
for i in range(0,ms):
    print("Total production in %s: %.2f (a change of %.2f percent)" % (markets[i],Q_post[i],Q_change[i]))
print("\n")
    
# Calculate market HHI values
hhi = np.zeros(ms)
for i in range(0,ms):
    ss = []
    for x in shares.iloc[:,i]:
        ss.append(np.square(x))
    hhi[i] = 10000*(sum(ss))    

# Display HHI values
print("Pre-merger Concentration\n===================================")
for i in range(0,ms):
    print("HHI for %s is: %.0f" %(markets[i],hhi[i]))
print("\n")

# Run the merger simulation
post_F = {}
u = MF[0]
v = MF[1]
post_F[1] = F[u].merge(F[v])
NM = [x for i,x in enumerate(F) if x!=u and x!=v]

for i in range(0,ns-2): 
    l = NM[i]
    post_F[i+2] = F[l]
post_ns = len(post_F)

# Define k_post and t_post (weighted average transport cost)
k_post = np.zeros(post_ns)
t_post = np.zeros((post_ns,ms))

k_post[0] = k[u-1]+k[v-1]
t_post[0,:] = np.array(t[u-1]*F[u].q/post_F[1].q) + np.array(t[v-1]*F[v].q/post_F[1].q)
for i in range(0,ns-2): 
    l = NM[i]
    k_post[i+1] = k[l-1]
    t_post[i+1] = t[l-1]
    
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
V = np.vstack(np.hsplit(a-t_post,ms))

# Account for efficiencies
if E!=0:
    for i in range(0,ms):
        for j in range(0,ms):
            M[i*post_ns,j*post_ns] += (-E/100)/k_post[0] 

# Solve the model
q_post = solve(M, V)

# Fill a dataframe with post-merger quantities
df_post = pd.DataFrame(columns =["Firm"]+post_F[1].mnames)
i=0
for i,f in enumerate(post_F):
    postprod=[]
    for j in range(0,ms):
        l = i+j*post_ns
        postprod.append(q_post[l][0])
    df_post.loc[i] = [post_F[f].name, *postprod]
    i+=1

# Calculate aggregate quantities
sums = df_post.sum(axis=0)    
Q_post = np.zeros(ms)
for i in range(0,ms):
    Q_post[i]=sums[i+1]

Q_change = 100*(Q_post-Q)/Q
shares_post = df_post.iloc[:,1:]/Q

# Display the post-merger production
print("Post-merger Production\n===================================")
display(df_post)
for i in range(0,ms):
    print("Total production in %s: %.2f (a change of %.2f percent)" % (markets[i],Q_post[i],Q_change[i]))
print("\n")

# Calculate market HHI values
hhi_post = np.zeros(ms)
for i in range(0,ms):
    ss = []
    for x in shares_post.iloc[:,i]:
        ss.append(np.square(x))
    hhi_post[i] = 10000*(sum(ss))    
del_hhi = hhi_post - hhi
    
# Display HHI values
print("Post-merger Concentration\n===================================")
for i in range(0,ms):
    print("Post-merger HHI for %s is: %.0f (a change of %.0f)" %(markets[i],hhi_post[i], del_hhi[i]))
print("\n")

# Calculate post-merger prices
p_post = a - b*Q_post
p_change = 100*(p_post-p)/p

# Display the post-merger prices
print("Post-merger Prices\n===================================")
for i in range(0,ms):
    print("%s: %.2f (a change of %.2f percent)" % (markets[i],p_post[i],p_change[i]))
print("\n")

# Calculate the impact on profits
prof_post = np.matmul(df_post.iloc[:,1:],p_post) - np.square(df_post.sum(axis=1))/(2*k_post) - t_post * df_post.iloc[:,1:]
print(prof_post)

# Calculate the harm to consumers
CS = np.zeros(ms)
CS_post = np.zeros(ms)
for i in range(0,ms):
    def integrand(x):
        return a[i]-b[i]*x
    CS[i], err = integrate.quad(integrand,0,Q[i])
    CS_post[i], err = integrate.quad(integrand,0,Q_post[i])
    
harm = CS - CS_post

# Report the estimated harm to consumers
print("Consumer Harm\n===================================")
for i in range(0,ms):
    print("Estimated harm to consumers in %s: %.3f" % (markets[i], harm[i]))
print("\n")
