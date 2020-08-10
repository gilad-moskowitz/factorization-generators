#Basic Functions for working with Factorization Sets in Numerical Semigroups

import math
import time

def SetAdd(L1, L2):
    L3 = []
    for a in L1:
        for b in L2:
            L3.append(a+b)
    L3 = list(set(L3))
    L3.sort()
    L3 = tuple(L3)
    return L3
    
    
    
    
def LSGensWithTime(S):
    start = time.time()
    c = lcm(S.gens)
    print '\n'
    print S.gens
    E = min(S.gens)*c
    W = int(math.floor(E/(max(S.gens))))
    S.LengthSetsUpToElement(E)
    LSUTEtime = time.time()
    
    Q = []
    M = []
    for N in [0..W]:
        M.append([])
    for N in [0..E]:
        if N in S:
            L = S.LengthSet(N)
            m = min(L)
            if m <= W:
                M[m].append(tuple(L))


    for N in [0..W]:
        M[N] = list(set(M[N]))
    AGtime = time.time()

    for N in [1..W]:
        m = int(math.floor(N/2))
        for L in M[N]:
            isGen = True
            for X in range(1, m+1):
                for G in M[X]:
                    h = N - X
                    if N-X == 0:
                        h = 1
                    for H in M[h]:
                        if SetAdd(G, H) == L:
                            isGen = False
                            break
                    if isGen == False:
                        break
                if isGen == False:
                    break


            if isGen == True:
                print L
                Q.append(L)
                
    print '\n'     
    print 'Total time: ', time.time() - start
    print 'Time for LengthSetsUpToElement: ', LSUTEtime - start
    print 'Time for generating array M: ', AGtime - LSUTEtime
    print 'Time for finding generators: ', time.time() - AGtime
    return Q
    
    
    
    
def LSGens(S):
    c = lcm(S.gens)
    print '\n'
    print S.gens
    E = min(S.gens)*c
    W = int(math.floor(E/(max(S.gens))))
    S.LengthSetsUpToElement(E)
    
    Q = []
    M = []
    for N in [0..W]:
        M.append([])
    for N in [0..E]:
        if N in S:
            L = S.LengthSet(N)
            m = min(L)
            if m <= W:
                M[m].append(tuple(L))


    for N in [0..W]:
        M[N] = list(set(M[N]))
    
    for N in [1..W]:
        m = int(math.floor(N/2))
        for L in M[N]:
            isGen = True
            for X in range(1, m+1):
                for G in M[X]:
                    h = N - X
                    if N-X == 0:
                        h = 1
                    for H in M[h]:
                        if SetAdd(G, H) == L:
                            isGen = False
                            break
                    if isGen == False:
                        break
                if isGen == False:
                    break


            if isGen == True:
                print L
                Q.append(L)
                
    return Q
    
    
    
    
    def LSGenNumsWithTime(S, L):
    start = time.time()
    N = [] 
    print '\n'
    print S.gens
    c = lcm(S.gens)
    E = min(S.gens)*c
    S.LengthSetsUpToElement(E)
    LSUTEtime = time.time()
    for A in L:
        for B in [1..E]:
            if tuple(S.LengthSet(B)) == A:
                print B
                N.append(B)
                break
    N = tuple(N)
    print '\n'     
    print 'Total time: ', time.time() - start
    print 'Time for LengthSetsUpToElement: ', LSUTEtime - start
    return N
    
    
    
    
    def LSGenNums(S, L):
    N = [] 
    print '\n'
    print S.gens
    c = lcm(S.gens)
    E = min(S.gens)*c
    S.LengthSetsUpToElement(E)
    for A in L:
        for B in [1..E]:
            if tuple(S.LengthSet(B)) == A:
                print B
                N.append(B)
                break
    N = tuple(N)
    return N
    
    
    
    
def FactAdd2(L1, L2):
    L3 = []
    b = 0
    for a in L1:
        b += 1
    for c in range(0, b):
        L3.append(L1[c] + L2[c])
    L3 = tuple(L3)
    return L3
    
    
    
    
def FactAdd3(L1, L2):
    L3 = []
    for a in L1:
        for b in L2:
            L3.append(FactAdd2(a, b))
    L3 = list(set(L3))
    return L3
    
    
    
    
def FactGensSoheilWithTime(S):
    c = lcm(S.gens)
    start = time.time()
    print '\n'
    print S.gens
    E = c*3
    S.FactorizationsUpToElement(E)
    LSUTEtime = time.time()
    
 
    Q = [min(S.gens)]
    print min(S.gens)
    for i in [1..E]:
        R = S.Factorizations(i)
        u = len(R)
        if i in S:
            isGen = True
            for m in Q:
                k = i - m
                if k in S:  
                    #R.sort()
                    X = S.Factorizations(m)
                    Y = S.Factorizations(k)
                    Z = FactAdd3(X, Y)
                    if len(Z) == u:
                        isGen = False
                        break
                else:
                    continue
            if isGen == True:
                if i == c:
                    print i, "LCM"
                else:
                    print i
                Q.append(i)
            
    print '\n'     
    print 'Total time: ', time.time() - start
    print 'Time for FactorizationsUpToElement: ', LSUTEtime - start
    print 'Time for finding generators: ', time.time() - LSUTEtime
    return Q
    
    
    
    
def FactGensSoheil(S):
    c = lcm(S.gens)
    print '\n'
    print S.gens
    E = c*3
    S.FactorizationsUpToElement(E)
    Q = [min(S.gens)]
    print min(S.gens)
    for i in [1..E]:
        R = S.Factorizations(i)
        u = len(R)
        if i in S:
            isGen = True
            for m in Q:
                k = i - m
                if k in S:  
                    #R.sort()
                    X = S.Factorizations(m)
                    Y = S.Factorizations(k)
                    Z = FactAdd3(X, Y)
                    if len(Z) == u:
                        isGen = False
                        break
                else:
                    continue
            if isGen == True:
                if i == c:
                    print i, "LCM"
                else:
                    print i
                Q.append(i)
    return Q
    
    
    
    
def LSGensNP(S):
    c = lcm(S.gens)
    E = min(S.gens)*c
    W = int(math.floor(E/(max(S.gens))))
    S.LengthSetsUpToElement(E)
    
    Q = []
    M = []
    for N in [0..W]:
        M.append([])
    for N in [0..E]:
        if N in S:
            L = S.LengthSet(N)
            m = min(L)
            if m <= W:
                M[m].append(tuple(L))


    for N in [0..W]:
        M[N] = list(set(M[N]))

    for N in [1..W]:
        m = int(math.floor(N/2))
        for L in M[N]:
            isGen = True
            for X in range(1, m+1):
                for G in M[X]:
                    h = N - X
                    if N-X == 0:
                        h = 1
                    for H in M[h]:
                        if SetAdd(G, H) == L:
                            isGen = False
                            break
                    if isGen == False:
                        break
                if isGen == False:
                    break
            if isGen == True:
                Q.append(L)
    return Q
    
    
    
    
def LSGenNumsNP(S, L):
    N = [] 
    c = lcm(S.gens)
    E = min(S.gens)*c
    S.LengthSetsUpToElement(E)
    for A in L:
        for B in [1..E]:
            if tuple(S.LengthSet(B)) == A:
                N.append(B)
                break
    N = tuple(N)
    return N
    
    
    
    
def FactGensSoheilNP(S):
    c = lcm(S.gens)
    E = c*3
    S.FactorizationsUpToElement(E)
   
    Q = [min(S.gens)]
    for i in [1..E]:
        R = S.Factorizations(i)
        u = len(R)
        if i in S:
            isGen = True
            for m in Q:
                k = i - m
                if k in S:  
                    #R.sort()
                    X = S.Factorizations(m)
                    Y = S.Factorizations(k)
                    Z = FactAdd3(X, Y)
                    if len(Z) == u:
                        isGen = False
                        break
                else:
                    continue
            if isGen == True:
                Q.append(i)
    return Q
    
    
    
    
def TestingLCM(S, n):
    c = lcm(S.gens)*2
    S.FactorizationsUpToElement(c + n + 1)
    for i in range(c, c+n):
        if (i - c) not in S:
            continue
        else:
            R = S.Factorizations(i)
            X = S.Factorizations(c)
            Y = S.Factorizations(i - c)
            Z = FactAdd3(X, Y)
            if len(R) == len(Z):
                continue
            else:
                BreakFactors(S, i)
                break
    return i
    
    
    
    
def BreakFactors(S, n):
    R = S.Factorizations(n)
    for m in range(min(S.gens), int(math.floor(n/2)) + 1):
        k = n - m
        if m in S:
            if k in S:
                X = S.Factorizations(m)
                Y = S.Factorizations(k)
                Z = FactAdd3(X, Y)
                if len(R) == len(Z):
                    print "(", m, ", ", k, ")"
    return n
    
    
    
    
def FactPrimeFactorization(S, n):
    P = []
    m = n
    if n not in S:
        return P
    gens = FactGensSoheilNP(S)
    if n in gens:
        P.append(n)
    else:
        while m > 0:
            R = S.Factorizations(m)
            if m in gens:
                P.append(m)
                break
            for a in gens:
                if m - a in S:
                    X = S.Factorizations(a)
                    Y = S.Factorizations(m - a)
                    Z = FactAdd3(X, Y)
                    if len(R) == len(Z):
                        P.append(a)
                        m = m - a
                        break
    return P
    
    
    
    
def TestingAssociativity(S, n):
    Ass = True
    if n not in S:
        return Ass
    P = FactPrimeFactorization(S, n)
    G = P[:]
    G = list(set(G))
    if len(G) == 1:
        return Ass
    for i in range(len(P) - 1):
        if Ass == False:
            return Ass
        for j in range(i + 1, len(P)):
            m = P[i] + P[j]
            if m not in S:
                continue
            R = S.Factorizations(m)
            X = S.Factorizations(P[i])
            Y = S.Factorizations(P[j])
            Z = FactAdd3(X, Y)
            if len(R) != len(Z):
                Ass = False
                break
    return Ass
    
    
    
    
def TestingAssociativityDyno(S, nmax):
    D = {}
    for n in range(nmax + 1):
        D[n] = False
    g = FactGensSoheilNP(S)
    for n in range(nmax + 1):
        if n not in S:
            D[n] = True
            continue
        if n in g:
            D[n] = True
            continue
        P = FactPrimeFactorization(S, n)
        Ass = True
        for i in range(len(P) - 1):
            if Ass == False:
                D[n] = Ass
                break
            for j in range(i + 1, len(P)):
                m = P[i] + P[j]
                if m not in S:
                    continue
                if D[m]:
                    continue
                R = S.Factorizations(m)
                X = S.Factorizations(P[i])
                Y = S.Factorizations(P[j])
                Z = FactAdd3(X, Y)
                if len(R) != len(Z):
                    Ass = False
                    break
                else:
                    D[m] = True
                    continue
        if Ass:
            D[n] = True
    Asso = True
    for i in range(nmax + 1):
        if D[i]:
            continue
        else:
            print i
    return Asso
    
    
    
    
def TestingAssociativityWP(S, n):
    Ass = True
    if n not in S:
        return Ass
    P = FactPrimeFactorization(S, n)
    if len(P) == 1:
        return Ass
    for i in range(len(P) - 1):
        if Ass == False:
            return Ass
        for j in range(i + 1, len(P)):
            m = P[i] + P[j]
            print i, j, m
            if m not in S:
                continue
            R = S.Factorizations(m)
            X = S.Factorizations(P[i])
            Y = S.Factorizations(P[j])
            Z = FactAdd3(X, Y)
            if len(R) != len(Z):
                Ass = False
                break
    return Ass
    
    
    
    
def isSubset(A, B): 
    for a in A: 
        bol = True
        for b in B: 
            if a == b: 
                bol = False
                break 
        if (bol): 
            return 0
    return 1
    
    
    
    
def isExtremalAddition2(S, a, b):
    isExtremal = True
    
    n = a + b
    S.FactorizationsUpToElement(2*n)
    
    aFacts = S.Factorizations(a)
    bFacts = S.Factorizations(b)
    Minkowski = FactAdd3(aFacts, bFacts)
    nFacts = S.Factorizations(n)
    Remainder = []
    for a in nFacts:
        if a not in Minkowski:
            Remainder.append(a)
    
    MinkowskiPoly = Polyhedron(vertices = Minkowski)
    ActualPoly = Polyhedron(vertices = nFacts)
    
    q1 = [tuple(v) for v in MinkowskiPoly.vertex_generator()]
    q2 = [tuple(v) for v in ActualPoly.vertex_generator()]
    
    extras = []
    for a in q2:
        if a not in q1:
            extras.append(a)
    
    ModifiedMinkowski = Minkowski[:]
    for a in Remainder:
        if a not in extras:
            ModifiedMinkowski.append(a)
    
    modification = []
    for b in ModifiedMinkowski:
        if b not in Minkowski:
            modification.append(b)
    
    ModifiedMinkowskiPoly = Polyhedron(vertices = ModifiedMinkowski)
    q3 = [tuple(v) for v in ModifiedMinkowskiPoly.vertex_generator()]
    
    remaining = []
    for c in modification:
        if c in q3:
            continue
        else:
            remaining.append(c)
    
    while (len(remaining) > 0):
        test = remaining[:]
        ModifiedMinkowski = Minkowski[:]
        for a in remaining:
            ModifiedMinkowski.append(a)
        remaining = modificationStep(Minkowski, ModifiedMinkowski)
        print remaining
        test = sorted(test)
        remaining = sorted(remaining)
        if test == remaining:
            print test
            isExtremal = False
            break
    
#    if isExtremal:
#        print "Success on All Cases"
#    else:
#        print "Failure"
        
    return isExtremal
    
    
    
    
def modificationStep (Minkowski, ModifiedMinkowski):
    modification = []
    for b in ModifiedMinkowski:
        if b not in Minkowski:
            modification.append(b)
    
    ModifiedMinkowskiPoly = Polyhedron(vertices = ModifiedMinkowski)
    q3 = [tuple(v) for v in ModifiedMinkowskiPoly.vertex_generator()]
    
    remainMod = []
    for c in modification:
        if c in q3:
            continue
        else:
            remainMod.append(c)
    return remainMod
    
    
    
    
def isExtremalAddition(S, a, b):
    isExtremal = True
    
    n = a + b
    S.FactorizationsUpToElement(2*n)
    
    aFacts = S.Factorizations(a)
    bFacts = S.Factorizations(b)
    Minkowski = FactAdd3(aFacts, bFacts)
    nFacts = S.Factorizations(n)
    Remainder = []
    for a in nFacts:
        if a not in Minkowski:
            Remainder.append(a)
    
    MinkowskiPoly = Polyhedron(vertices = Minkowski)
    ActualPoly = Polyhedron(vertices = nFacts)
    
    q1 = [tuple(v) for v in MinkowskiPoly.vertex_generator()]
    q2 = [tuple(v) for v in ActualPoly.vertex_generator()]
    
    extras = []
    for a in q2:
        if a not in q1:
            extras.append(a)
    
    ModifiedMinkowski = Minkowski[:]
    for a in Remainder:
        if a not in extras:
            ModifiedMinkowski.append(a)
    
    remaining = modificationStep(Minkowski, ModifiedMinkowski)
    
    while (len(remaining) > 0):
        test = remaining[:]
        ModifiedMinkowski = Minkowski[:]
        for a in remaining:
            ModifiedMinkowski.append(a)
        remaining = modificationStep(Minkowski, ModifiedMinkowski)
        test = sorted(test)
        remaining = sorted(remaining)
        if test == remaining:
            print test
            isExtremal = False
            break
    
#    if isExtremal:
#        print "Success on All Cases"
#    else:
#        print "Failure"
        
    return isExtremal
    
    
    
    
def isOneExtremal(S, a, b):
    isOneExtremal = False
    
    n = a + b
    S.FactorizationsUpToElement(2*n)
    
    aFacts = S.Factorizations(a)
    bFacts = S.Factorizations(b)
    Minkowski = FactAdd3(aFacts, bFacts)
    nFacts = S.Factorizations(n)
    Remainder = []
    for a in nFacts:
        if a not in Minkowski:
            Remainder.append(a)
    if (len(Remainder) == 0):
        return True
    
    MinkowskiPoly = Polyhedron(vertices = Minkowski)
    ActualPoly = Polyhedron(vertices = nFacts)
    
    q1 = [tuple(v) for v in MinkowskiPoly.vertex_generator()]
    q2 = [tuple(v) for v in ActualPoly.vertex_generator()]
    
    extras = []
    for a in q2:
        if a not in q1:
            extras.append(a)
    
    if (len(extras) > 0):
        return True
    
#    ModifiedMinkowski = Minkowski[:]
#    for a in Remainder:
#        if a not in extras:
#            ModifiedMinkowski.append(a)
    
#    modification = []
#    for b in ModifiedMinkowski:
#        if b not in Minkowski:
#            modification.append(b)
    
#    ModifiedMinkowskiPoly = Polyhedron(vertices = ModifiedMinkowski)
#    q3 = [tuple(v) for v in ModifiedMinkowskiPoly.vertex_generator()]
    
#    for c in modification:
#        if c in q3:
#            isOneExtremal = True
#            break
#        else:
#            continue
        
    return isOneExtremal
