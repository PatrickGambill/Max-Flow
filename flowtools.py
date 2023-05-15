# -*- coding: utf-8 -*-
"""
Created on Sun Dec  4 14:42:55 2022

@author: gambi
"""
import math
import matplotlib.pyplot as plt

def genrandom(n,m,U):
    
    '''
    Inputs: 
    n: number of nodes
    m: number of arcs
    U: max arc weight
    
    Outputs: 
    pointer: The pointer vector
    arcs: The arc set 
    
    This function will randomly generate the forward star representation of a random network on 
    n nodes with m arcs and a max arc weight of U.
    '''
    
    import random
    
    #There are n(n-1) possible arcs on a random network. If we remove the arcs that start with n, we get
    # n^2 - 2n + 1 possible arcs. To remove the arcs that end at 0, we remove the nums == 0 mod(n-1), except for 0.
    nums = [i for i in range(n**2 - 2*n + 1) if not(i%(n-1)==0) or i==0]
    #Sample picks some number of arcs from the arcs corresponding to nums
    sample = random.sample(nums, m)
    sample.sort()
    
    
    #Now that we have created the sample of nums, we can convert this to the arc set
    arcs = []
    for i in sample:
        arc = []
        arc.append(i//(n-1)) 
        
        if i%(n-1) < i//(n-1):
            arc.append(i%(n-1))
        else:
            arc.append(i%(n-1) + 1)
        
        arc.append(random.randint(1,U))
        arcs.append(arc)
    
    #Create the pointer vector
    pointer = [0]
    
    for i in range(n-1):
        #This helps if we run out of arcs early
        f = True
        for m in range(len(arcs)):
            if arcs[m][0] > i:
                f = False
                pointer.append(m)
                break
        if f:
            pointer.append(len(arcs))
            
    pointer.append(len(arcs))
    
    return pointer, arcs

def genrandoms(n,m,U, numnetworks, tries = 1000):
    
    '''
    Inputs: 
    n: number of nodes
    m: number of arcs
    U: max arc weight
    numnetworks: The number of random networks to create
    
    Outputs: 
    networks: A list of tuples, each containing a forward star representation of a random network
    
    This function will randomly generate a list containing the forward star representation of random networks on 
    n nodes with m arcs and a max arc weight of U.
    '''
    
    networks = []
    count = 0
    it = 0
    
    while count < numnetworks and it < tries:
        pointer, arcs = genrandom(n,m,U)
        
        #Check to make sure a path exists
        pred, order = bfs(pointer, arcs, 0)
        
        if not pred[n-1] == float('inf'):
            networks.append((pointer,arcs))
            count += 1
        
        it += 1
        
    if it == tries:
        print("Only " + str(len(networks)) + " networks were found in " + str(tries) + " tries.")
        
    return networks

def bfs(pointer, arcs, s):
    '''This function will perform a BFS on a network represented by the forward star representation. 
    pointer is the pointer vector represented as a list
    arcs is the list of of arcs (and possibly weights)
    s is the node to start the BFS at
    
    This function will return the pred and order vectors, each represented as a list.
    returns: (pred,order)
    
    Note, Python starts indexing at 0, so this algorithm will as well.'''
    
    #0 means unmarked node, 1 means marked node
    n = len(pointer)-1
    marks = [0 for i in range(n)]
    marks[s] = 1
    
    #pred[i]=-1 means the node is first. pred[i]=-2 means the node hasn't been visited yet
    pred = [float('inf') for i in range(n)]
    pred[s] = -1
    
    next = 0
    
    #-1 represents an unvisited node
    order = [-1 for i in range(n)]
    order[s] = next
    
    LIST = [s]
    
    while len(LIST) > 0:
        #FIFO for BFS
        i = LIST[0]
        #This skips the case with outdegree 0
        if pointer[i] != pointer[i+1]:
            #Check each arc
            for x in range(pointer[i],pointer[i+1]):
                try:
                    j = arcs[x][1]
                except:
                    break
                #This takes only admissable arcs
                try:
                    if marks[j] == 0:
                        marks[j] = 1
                        pred[j] = i
                        next += 1
                        order[j] = next
                        LIST.append(j)
                except:
                    break
        LIST.remove(i)
    
    return pred, order

def SAP(pointer, arcs, s, t):
    '''This function will find the max s-t flow using SAPMaxFlow 
    
    inputs:
        
    pointer: The pointer vector represented as a list
    arcs: List of of arcs with capacities. The order should be [i,j,u_ij] for the arc from i to j with capacity u_ij
    s: source node
    t: sink node
    
    outputs:
    (v,x)
    v: flow value
    x: max flow vector
    

    
    Note, Python starts indexing at 0, so this algorithm will as well.'''
    
    #Create the node set
    n = len(pointer)-1
    N = [i for i in range(n)]
    
    #Initialize the x vector and v
    v = 0
    x = [0 for i in arcs]
    
    #residual vector
    r = [arc[2] for arc in arcs]
    #Create the initial pred vector
    pred = [float('inf') for i in range(n)]
    pred[s] = -1
    
    #Get the exact distance labels
    
    d = [float('inf') for i in range(n)]
    d[t] = 0
    
    back = [t]
    front = []
    while not len(back) == 0:
        
        for arc in arcs:

            if arc[1] in back:
                if d[arc[0]] == float('inf'):
                    d[arc[0]] = d[arc[1]] + 1
                    front.append(arc[0])
        
        back = front
        front = []
    
    
    currentNode = s
    
    #Main Loop
    while d[s]<n:
        
        #Check for admissible arcs
        admissibleArc = False
        
        for i in range(len(arcs)):
            arc = arcs[i]
            
            #Advance along admissible arc
            if d[arc[0]] == d[arc[1]]+1 and arc[0] == currentNode and r[i]>0:

                admissibleArc = True
                pred[arc[1]] = currentNode
                currentNode = arc[1]
                
                #Augment if t is reached
                if arc[1] == t:
                    
                    P = [t]
                    while not currentNode == s:
                        nextNode = pred[currentNode]
                        P.append(nextNode)
                        currentNode = nextNode
                    #currentNode will be s here, so we need not repeat this line after augmenting
                    
                    #find delta
                    aP = []
                    #First find the augmenting path in terms of arcs
                    for a in range(len(arcs)):
                        for j in range(len(P)-1):
                            if arcs[a][0] == P[j+1] and arcs[a][1] == P[j]:
                                aP.append(a)
                    
                    delta = min(arcs[a][2] for a in aP)
                    
                    #Augment the flow along P
                    for a in aP:
                        r[a] -= delta
                        x[a] += delta
                    v += delta
        
        #retreat if no admissible arc is found
        if not admissibleArc:
            
            minDist = float('inf')
            for i in range(pointer[currentNode], pointer[currentNode + 1]):
                if d[arcs[i][1]] + 1 < minDist and r[i] > 0:
                    minDist = d[arcs[i][1]] + 1
            d[currentNode] = minDist
            
            if not currentNode == s:
                currentNode = pred[currentNode]
        
        
                    
        
    return v, x

def FIFOPreflow(pointer, arcs, s, t):
    '''This function will find the max s-t flow using FIFO preflow push 
    
    inputs:
        
    pointer: The pointer vector represented as a list
    arcs: List of of arcs with capacities. The order should be [i,j,u_ij] for the arc from i to j with capacity u_ij
    s: source node
    t: sink node
    
    outputs:
    
    x: max flow vector
    v: flow value
    sat: number of saturating pushes
    nonsat: number of nonsaturating pushes
    
    Note, Python starts indexing at 0, so this algorithm will as well.'''
    
    #Preprocess
    n = len(pointer)-1
    
    v = 0
    sat = 0
    nonsat = 0
    
    #Preprocess
    
    
    #Flow vector
    x = [0 for i in arcs]
    #Residual vector
    r = [i[2] for i in arcs]
    #Excess flow 
    e = [0 for i in range(n)]
    
    #Get distances from t using BFS
    
    revarcs = [[arc[1],arc[0],arc[2]] for arc in arcs]
    revarcs.sort()
    
    #Get the reverse pointer
    #Create the pointer vector
    revpointer = [0]
    
    for i in range(n-1):
        #This helps if we run out of arcs early
        f = True
        for m in range(len(arcs)):
            if revarcs[m][0] > i:
                f = False
                revpointer.append(m)
                break
        if f:
            revpointer.append(len(arcs))
            
    revpointer.append(len(arcs))
    
    pred, order = bfs(revpointer, revarcs, t)
    
    #Get the distances from pred
    d = [n for i in range(n)]
    d[t] = 0
    for i in range(n-1):
        count = 0
        p=i
        if pred[p] == float('inf'):
            d[i] = float('inf')
            continue
        for j in range(n-1):
            p = pred[p]
            count += 1
            if p == t:
                break
        d[i] = count
    
    #active nodes
    
    active = []
    for i in range(pointer[s],pointer[s+1]):
        
        x[i] = arcs[i][2]
        r[i] = 0
        e[arcs[i][1]] = arcs[i][2]
        active.append(arcs[i][1])
    
    d[s] = n
    
    
    #While we have active nodes
    while not len(active) == 0:
        
        currentNode = active[0]
        if e[currentNode] == 0:
            active.pop(0)
            continue
        
        #Check for admissible arcs
        
        admissibleArc = False
        for i in range(len(arcs)):
            arc = arcs[i]
            
            if e[t] >0:
                v += e[t]
                e[t] = 0
            #Push along arcs
            if d[arc[0]] == d[arc[1]]+1 and arc[0] == currentNode and not d[arc[0]]==float('inf') and not r[i] == 0:
                delta = min(e[currentNode],r[i])
                

                if not delta <= 0:
                    
                    
                    #Check if the flow is saturating
                    if delta == r[i]:
                        sat +=1
                    else:
                        nonsat +=1
                    
                    #push the flow
                    admissibleArc = True
                    e[currentNode] -= delta
                    r[i] -= delta
                    x[i] += delta
                    if (not arcs[i][0] == s) and (not arcs[i][1] == t):
                        e[arcs[i][1]] += delta
                    
                    if arc[1] == t:
                        v += delta
                
                else:
                    if not len(active) == 0:
                        if e[currentNode] == 0:
                            active.pop(0)
                    
                    
            #Check to send flow back along arcs
            elif d[arc[1]] == d[arc[0]]+1 and arc[1] == currentNode and not d[arc[0]]==float('inf') and not r[i] == 0:
                
                delta = min(e[currentNode],x[i] - r[i])
                
                if not delta <= 0:
                    #Check if the flow is saturating
                    if delta == x[i] - r[i]:
                        sat +=1
                    else:
                        nonsat +=1
                    #push the flow
                    admissibleArc = True
                    e[currentNode] -= delta
                    r[i] += delta
                    x[i] -= delta
                    if (not arcs[i][0] == s) and (not arcs[i][1] == t):
                        e[arcs[i][0]] += delta
                
                else:
                    minDist = float('inf')
                    for i in range(len(arcs)):
                        if d[arcs[i][1]] + 1 < minDist and r[i] > 0 and arcs[i][0] == currentNode:
                            minDist = d[arcs[i][1]] + 1
                        if d[arcs[i][0]] + 1 < minDist and arcs[i][2] - r[i] > 0 and arcs[i][1] == currentNode:
                            minDist = d[arcs[i][0]] + 1
            
                    d[currentNode] = minDist
                        

        if not admissibleArc:
            minDist = float('inf')
            for i in range(len(arcs)):
                if d[arcs[i][1]] + 1 < minDist and r[i] > 0 and arcs[i][0] == currentNode:
                    minDist = d[arcs[i][1]] + 1
                if d[arcs[i][0]] + 1 < minDist and arcs[i][2] - r[i] > 0 and arcs[i][1] == currentNode:
                    minDist = d[arcs[i][0]] + 1
                    
            d[currentNode] = minDist
            
            active.pop(0)
            continue

        for i in range(n):
            if (not e[i] == 0) and (not i in active):
                active.append(i)
    
    return x,v,sat,nonsat

def CapScale(pointer, arcs, s, t):
    '''This function will find the max s-t flow using the Capacity Scaling Algorithm 
    
    inputs:
        
    pointer: The pointer vector represented as a list
    arcs: List of of arcs with capacities. The order should be [i,j,u_ij] for the arc from i to j with capacity u_ij
    s: source node
    t: sink node
    
    outputs:
    
    x: max flow vector
    v: flow value
    
    Note, Python starts indexing at 0, so this algorithm will as well.'''
    
    import math
    
    #Preprocess
    n = len(pointer)-1
    m = len(arcs)
    v = 0
    
    #Flow vector
    x = [0 for i in arcs]
    #Residual vector
    r = [i[2] for i in arcs]

    #Get delta
    U = max([arc[2] for arc in arcs])
    delta = 2**(math.floor(math.log2(U)))
    
    while delta >= 1:
        
        #Get G(x,delta)
        xG = [0 for i in arcs]
        rG = []
        for i in range(len(r)):
            if r[i] >= delta:
                rG.append(r[i])
            else:
                rG.append(0)
                
        arcsG = []
        for i in range(len(arcs)):
            if rG[i] > 0:
                arcsG.append(arcs[i])
        
        pointerG = [0]
        if len(arcsG) == 0:
            pointerG += [len(arcs) for i in range(n-1)]
        
        else:
            for i in range(n-1):
                f = True
                for mG in range(len(arcsG)):
                    if arcsG[mG][0] > i:
                        f = False
                        pointerG.append(mG)
                        break
                if f:
                    pointerG.append(len(arcsG))
            
        pointerG.append(len(arcsG))
        
        
        #Use BFS to check for paths
        predG, orderG = bfs(pointerG, arcsG, s)
        
        #Augment as long as a path exists from s to t
        
        while not(predG[t] == float('inf')):
            
            #Get a path
            path = [t]
            p = predG[t]
            while p >= -1:
                path.append(p)
                p = predG[p]
                if p == -1:
                    break   
            path.reverse()
            
            
            #Get the arcs in the path
            patharcs = [path[i:i+2] for i in range(len(path)-1)]
            #Augment along the path
            nocap = [arc[0:2] for arc in arcs]
            arcind = [nocap.index(arc) for arc in patharcs]
            delt = min([r[ind] for ind in arcind])
            
            

            for i in arcind:
               x[i] += delt
               r[i] -= delt
               rG[i] -= delt
            v += delt
            
            #Update G(x,delta)
            xG = [0 for i in arcs]
            rG = []
            for i in range(len(r)):
                if r[i] >= delta:
                    rG.append(r[i])
                else:
                    rG.append(0)
                    
            arcsG = []
            for i in range(len(arcs)):
                if rG[i] > 0:
                    arcsG.append(arcs[i])
            
            pointerG = [0]
            if len(arcsG) == 0:
                pointerG += [len(arcs) for i in range(n-1)]
            
            else:
                for i in range(n-1):
                    f = True
                    for mG in range(len(arcsG)):
                        if arcsG[mG][0] > i:
                            f = False
                            pointerG.append(mG)
                            break
                    if f:
                        pointerG.append(len(arcsG))
                
            pointerG.append(len(arcsG))
            
            predG, orderG = bfs(pointerG, arcsG, s)
            
            
        delta //= 2
    
    return x,v

def GenConditions():
    '''This function will generate the conditions to be used for this project. 
    
    inputs: None
    
    outputs:
    
    Conditions: A list containing tuples of the form (n,m,U) where 
                n - Number of nodes
                m - Number of arcs
                U - Max arc capacity
                
    '''
    
    #Create the conditions to test
    conditions = []
    
    #Allowed n values
    for n in [10*n for n in range(1,5)]:
        possible_arcs = (n - 1) * (n - 2)
        
        #Allowed m values will be a percentage of the total possible arcs
        for m in [math.floor(possible_arcs*i) for i in [0.2,0.4,0.6,0.8]]:
            
            #The max U values will be powers of 2 from 2^4 to 2^8
            for U in [2**n for n in range(4,9)]:
                
                conditions.append((n,m,U))
    
    return conditions


def GenNetworks(conditions, seed = None, repeats = 50):
    
    '''This function will generate a list of networks, each in forward star representation. 
    Each network will have n nodes, a percentage o the possible arcs, and a max capacity. 
    Every five networks will have the same random conditions to start with.
    
    inputs:
    
    conditions: (n,m,U) A list of tuples representing the conditions for the networks to have. 
                  n - Number of nodes
                  m - Number of arcs
                  U - Max arc capacity
    seed: A random seed for reproducable results. The default value is None (which will use current system time)
    
    outputs:
    
    networks: A list of networks, each represented as a tuple. The first element of the tuple is a list
              representing the pointer. The second element is a list of arcs with capacity
    
    '''
    import random
    random.seed()
    
    #The set of networks to test
    networks = []
    #Generate 5 networks with each condition
    for cond in conditions:
        networks += genrandoms(cond[0],cond[1],cond[2],repeats)
    
    return networks
    
def GetTimes(networks, repeats = 1):
    '''
    This function will return the average runtime of the SAP, FIFOPreflow, and CapScale algorithms implemented above.
    These runtimes will be returned in the same order as the conditions used to generate the networks.
    
    inputs:
    
    networks: The list networks generated under previous conditions
    repeats: The number of networks generated for each condition
    
    outputs:
    
    times: This is a list of tuples, where each element in the list corresponds to a condition, and each 
           tuple returns the average runtime for SAP, FIFOPreflow, and CapScale in that order.
    
    '''

    import time
    
    times = []
    
    num = len(networks)
    for i in range(0,len(networks),repeats):
        print("Testing network " + str(i) + " of " + str(num))
        try:
            #SAP
            SAPtime = 0
            #FIFOPreflow
            FIFOPreflowtime = 0
            #CapScale
            CapScaletime = 0
            for j in range(repeats):  
                network = networks[i+j]
                pointer,arcs = network[0],network[1]
                t = len(pointer) - 2
                
                #SAP
                SAPtime -= time.time()
                SAP(pointer,arcs,0, t)
                SAPtime += time.time()
                
                #FIFOPreflow
                FIFOPreflowtime -= time.time()
                FIFOPreflow(pointer,arcs,0, t)
                FIFOPreflowtime += time.time()
            
                #CapScale
                CapScaletime -= time.time()
                CapScale(pointer,arcs,0, t)
                CapScaletime += time.time()
        
        except:
            
            pass
        
        
        times.append((SAPtime,FIFOPreflowtime,CapScaletime))
        
    return(times)

def plot(x,SAPtimes,FIFOPreflowtimes,CapScaletimes,SAPworst,FIFOPreflowworst,CapScaleworst,indices = False,title = '',xlabel='',ylabel='runtimes'):
    
    '''
    This function will plot the runtime of each algorithm against the given dataset x
    
    inputs:
    
    x: The full dataset of x values to conpare
    SAPtimes: The runtimes using SAP on each network
    FIFOPreflowtimes: The runtimes using FIFOPreflow on each network
    CapScaletimes: The runtimes using CapScale on each network
    indices: This will plot only the values at the indices given here. If False, this will plot all data points
    title: Title of plot
    xlabel: Label for the x-axis
    ylabel: label for the y-axis
    
    
    outputs:
    
    none
    
    '''
    
    if indices:
        n = len(x)
        x = [x[i] for i in range(n) if i in indices]
        SAPtimes = [SAPtimes[i] for i in range(n) if i in indices]
        FIFOPreflowtimes = [FIFOPreflowtimes[i] for i in range(n) if i in indices]
        CapScaletimes = [CapScaletimes[i] for i in range(n) if i in indices]
        SAPworst = [SAPworst[i] for i in range(n) if i in indices]
        FIFOPreflowworst = [FIFOPreflowworst[i] for i in range(n) if i in indices]
        CapScaleworst = [CapScaleworst[i] for i in range(n) if i in indices]
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        
        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_color('none')
        ax.spines['left'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
        
        ax1.scatter(x, SAPtimes, c='b', marker="s", label='SAP')
        ax1.scatter(x, FIFOPreflowtimes, c='r', marker="o", label='FIFO')
        ax1.scatter(x,CapScaletimes, c='g', marker="d", label='CapScale')
        ax1.legend(loc='upper left')
        ax1.set_ylabel(ylabel)
        
        ax2.scatter(x, SAPworst, c='b', marker="s", label='SAP')
        ax2.scatter(x, FIFOPreflowworst, c='r', marker="o", label='FIFO')
        ax2.scatter(x,CapScaleworst, c='g', marker="d", label='CapScale')
        ax2.legend(loc='upper left')
        ax2.set_ylabel('Worst Case Runtimes')
        
        ax.set_xlabel(xlabel)
        fig.suptitle(title)
        plt.show()
        
        return
    
    else:     
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        
        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_color('none')
        ax.spines['left'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
        
        ax1.scatter(x, SAPtimes, c='b', marker="s", label='SAP')
        ax1.scatter(x, FIFOPreflowtimes, c='r', marker="o", label='FIFO')
        ax1.scatter(x,CapScaletimes, c='g', marker="d", label='CapScale')
        ax1.legend(loc='upper left')
        ax1.set_ylabel(ylabel)
        
        ax2.scatter(x, SAPworst, c='b', marker="s", label='SAP')
        ax2.scatter(x, FIFOPreflowworst, c='r', marker="o", label='FIFO')
        ax2.scatter(x,CapScaleworst, c='g', marker="d", label='CapScale')
        ax2.legend(loc='upper left')
        ax2.set_ylabel('Worst Case Runtimes')
        
        ax.set_xlabel(xlabel)
        fig.suptitle(title)
        plt.show()
        
        return
        
    
#This step may take a bit of time, be patient
conditions = GenConditions()
networks = GenNetworks(conditions, seed = 'oof')
times = GetTimes(networks, repeats = 1)

#Turn the times into average times
avgtimes = []
for i in range(0, len(times), 50):
    t0, t1, t2 = 0,0,0
    for j in range(50):
        time = times[i+j]
        t0 += time[0]
        t1 += time[1]
        t2 += time[2]
    t0 /= 50
    t1 /= 50
    t2 /= 50
    avgtimes.append((t0,t1,t2))
        
#Split the times by algorithm
SAPtimes = [time[0] for time in avgtimes]
FIFOPreflowtimes = [time[1] for time in avgtimes]
CapScaletimes = [time[2] for time in avgtimes]

#Split the conditions into a list of n,m,U
ns = [condition[0] for condition in conditions]
ms = [condition[1] for condition in conditions]
Us = [condition[2] for condition in conditions]

#Get the number of steps in the worst cases
SAPworst = [con[0]**2 * con[1] for con in conditions]
FIFOPreflowworst = [con[0]**3 for con in conditions]
CapScaleworst = [con[1] ** 2 * math.log2(con[2]) for con in conditions]

#Plot the runtimes

#All data plots
plot(ns,SAPtimes,FIFOPreflowtimes,CapScaletimes,SAPworst,FIFOPreflowworst,CapScaleworst, indices = False, title = 'n vs runtimes', xlabel='n', ylabel='runtimes')
plot(ms,SAPtimes,FIFOPreflowtimes,CapScaletimes,SAPworst,FIFOPreflowworst,CapScaleworst, indices = False, title = 'm vs runtimes', xlabel='m', ylabel='runtimes')
plot(Us,SAPtimes,FIFOPreflowtimes,CapScaletimes,SAPworst,FIFOPreflowworst,CapScaleworst, indices = False, title = 'U vs runtimes', xlabel='U', ylabel='runtimes')

#Restrict to one variable at a time.



#Get back the possible values for each condition
npossible = [10*n for n in range(1,5)]
mpossible = [0.2,0.4,0.6,0.8]
Upossible = [2**n for n in range(4,9)]

#For fixed n,m plot the Us
for n in npossible:
    possible_arcs = (n - 1) * (n - 2)
    mn = [math.floor(possible_arcs*i) for i in mpossible]
    for m in mn:
        ind = [i for i in range(len(conditions)) if conditions[i][0] ==n and conditions[i][1] ==m]
        title = 'U vs runtimes with n = ' + str(n) + " and m = " + str(m)
        plot(Us,SAPtimes,FIFOPreflowtimes,CapScaletimes,SAPworst,FIFOPreflowworst,CapScaleworst, indices = ind, title = title, xlabel='U', ylabel='runtimes')

#For fixed n,U plot the ms
for n in npossible:
    for U in Upossible:
        ind = [i for i in range(len(conditions)) if conditions[i][0] ==n and conditions[i][2] ==U]
        title = 'm vs runtimes with n = ' + str(n) + " and U = " + str(U)
        plot(ms,SAPtimes,FIFOPreflowtimes,CapScaletimes,SAPworst,FIFOPreflowworst,CapScaleworst, indices = ind, title = title, xlabel='m', ylabel='runtimes')


#For fixed n,U plot the ms
for U in Upossible:
    for mper in mpossible:
        ind = [i for i in range(len(conditions)) if conditions[i][1] == math.floor(((conditions[i][0] - 1) * (conditions[i][0] - 2))*mper) and conditions[i][2] ==U]
        title = 'n vs runtimes with m = ' + str(mper) + "% of the maximum possible arcs and U =" + str(U)
        plot(ms,SAPtimes,FIFOPreflowtimes,CapScaletimes,SAPworst,FIFOPreflowworst,CapScaleworst, indices = ind, title = title, xlabel='m', ylabel='runtimes')