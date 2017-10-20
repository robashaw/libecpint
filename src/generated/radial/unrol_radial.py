MAX_UNROL_AM = 4

# DO NOT EDIT BELOW HERE

from sympy import *

class Qijk:
    def __init__(self, Ival = 0, Jval = 0, Kval = 0):
        self.i = Ival
        self.j = Jval
        self.k = Kval
        self.size = 2*self.i + self.j + 1
        self.start = self.k-self.i-self.j
        self.end = self.k+self.i
        self.subq = []
        self.terms = []
        self.bases = []
        self.f = []
        self.ga = []
        self.gb = []
        self.h = []
        for i in range(self.size):
            self.bases.append([])
    
    def print(self):
        print("Q", self.i, self.j, self.k)
        for i in range(len(self.bases)):
            print(i+self.start, ":")
            for j in range(len(self.bases[i])):
                print(self.bases[i][j])
    
    def print_simple(self):
         print("Q", self.i, self.j, self.k)
         for i in range(len(self.bases)):
            print(i+self.start, ":", simplify(self.bases[i]))
            
    def print_fgh(self):
        print("Q", self.i, self.j, self.k)
        for i in range(len(self.f)):
            ix = 2*i+self.start
            if ix > 0:
                print("F", ix, ":", simplify(self.f[i]))
        for i in range(len(self.gb)):
            ix = 2*i + self.start + 1
            if ix > 0:
                print("GB", ix, ":", simplify(self.gb[i]))
        for i in range(len(self.ga)):
            ix = 2*i + self.start + 1
            if ix > 0:
                print("GA", ix, ":", simplify(self.ga[i]))
        for i in range(len(self.h)):
            ix = 2*i+self.start + 2
            if ix > 0:
                print("H", ix, ":", simplify(self.h[i]))
    
    def write_code(self, f):
        print("\t\t\t\t\t\t\t\t\tcase", self.i*10000+self.j*100+self.k, ": {", file=f)

        for i in range(len(self.f)):
            ix = 2*i+self.start
            if ix == 2:
                simp = simplify(self.f[i])
                if simp!= 0:
                    print("\t\t\t\t\t\t\t\t\t\tresult = (", simp, ") * values[0];", file=f)
            elif ix > 2:
                simp = simplify(self.f[i])
                if simp!= 0:
                    print("\t\t\t\t\t\t\t\t\t\tresult += (", simp, ") * values[", ix-2, "];", file=f)
                
        for i in range(len(self.gb)):
            ix = 2*i + self.start + 1
            if ix == 1:
                simp = simplify(self.gb[i])
                if simp!= 0:
                    print("\t\t\t\t\t\t\t\t\t\tresult += (", simp, ") * G1B;", file=f)
            elif ix > 1:
                simp = simplify(self.gb[i])
                if simp!= 0:
                    print("\t\t\t\t\t\t\t\t\t\tresult += (", simp, ") * values[", ix-2, "];", file=f)

        for i in range(len(self.ga)):
            ix = 2*i + self.start + 1
            if ix == 1:
                simp = simplify(self.ga[i])
                if simp != 0:
                    print("\t\t\t\t\t\t\t\t\t\tresult += (", simp, ") * G1A;", file=f)

        for i in range(len(self.h)):
            ix = 2*i+self.start + 2
            if ix == 2:
                simp = simplify(self.h[i])
                if simp!= 0:
                    print("\t\t\t\t\t\t\t\t\t\tresult += (", simp, ") * H2;", file=f)
                    
        print("\t\t\t\t\t\t\t\t\t\tbreak;", file=f)
        print("\t\t\t\t\t\t\t\t\t}", file=f)

    def simplify(self):
        simple_bases = []
        for i in range(len(self.bases)):
            x = Symbol('x')
            y = Symbol('y')
            z = Symbol('z')
            z = 0
            for j in range(len(self.bases[i])):
                z = z + parse(self.bases[i][j])
            simple_bases.append(z)
        self.bases = simple_bases
        
    def sort(self):
        for i in range(len(self.bases)):
            if i % 2 == 0:
                self.f.append(self.bases[i])
            else:
                self.gb.append(self.bases[i])
                
    def eliminate(self):
        x = Symbol('x')
        y = Symbol('y')
        z = Symbol('z')
        p = Symbol('p')
        z = 0
        
        if self.start < 1:
            if self.end < 1:
                self.bases.append(z)
            if self.end < 2:
                self.bases.append(z)
                
            w = Symbol('w')
            
            N = self.start
            ix = 0
            gaix = -1
            hix  = -1
            
            while (N < 1):
                z = self.bases[ix]
                w = self.bases[ix+2]
                w = w + (2 * p / (N-1))*z
                self.bases[ix+2] = w
            
                w = self.bases[ix+1]
                w = w - (2 * y / (N-1))*z
                self.bases[ix+1] = w
                
                if ix % 2 == 0:
                    if gaix > -1:
                        w = self.ga[gaix]
                        w = w - (2*x / (N-1))*z
                        self.ga[gaix] = w
                    else:
                        w = -(2*x / (N-1))*z
                        self.ga.append(w)
                        gaix += 1
                        
                    if hix > -1:
                        z = self.h[hix]
                        w = (2 * p / (N-1))*z
                        self.h.append(w)
                        hix += 1
                        
                        w = self.ga[gaix]
                        w = w - (2 * y / (N-1))*z
                        self.ga[gaix] = w
                        
                        w = self.bases[ix+1]
                        w = w - (2 * x / (N-1))*z
                        self.bases[ix+1] = w
                else:
                    if hix > -1:
                        w = self.h[hix]
                        w = w - (2*x / (N-1))*z
                        self.h[hix] = w
                    else:
                        w = - (2*x / (N-1))*z
                        self.h.append(w)
                        hix += 1
                    
                    if gaix > -1:
                        z = self.ga[gaix]
                        w = (2 * p /(N-1)) * z
                        self.ga.append(w)
                        gaix += 1
                        
                        w = self.h[hix]
                        w = w - (2 * y / (N-1))*z
                        self.h[hix] = w
                        
                        w = self.bases[ix+1]
                        w = w - (2*x / (N-1))*z
                        self.bases[ix+1] = w
                    
                N += 1
                ix += 1

def parse(term):
    x = Symbol('x')
    y = Symbol('y')
    p = Symbol('p')
    bits = term.split(',')
    z = Symbol('z')
    z = 1
    for bit in bits:
        bi = bit[:2]
        if bi == "mu":
            ix = 2
            i = 0
            j = 0
            k = 0

            I = bit[ix]
            if I == "-":
                ix += 1
                i = -int(bit[ix])
            else:
                i = int(I)
                
            ix += 1
            J = bit[ix]
            if J == "-":
                ix += 1
                j = -int(bit[ix])
            else:
                j = int(J)

            ix += 1
            K = bit[ix]
            if K == "-":
                ix += 1
                k = -int(bit[ix])
            else:
                k = int(K)

            z = z * (2 + j - i - k)/(2*x)
        elif bi == "nu":
            z = z * (-y/x)
        elif bi == "xi":
            z = z * p/x
        elif bi == "rh":
            j = int(bit[3])
            z = z * (1 - 2*j)/(2*y)
        elif bi == "om":
            z = z * -1 / (2*y)
    return z
            
def unrol(q):
    if (q.i == 0 and q.j == 0):
        return
    elif (q.i > 0):
        q1 = Qijk(Ival = q.i-1, Jval = q.j, Kval = q.k-1)
        q.subq.append(q1)
        q.terms.append("mu" + str(q.i) + str(q.j) + str(q.k))

        q2 = Qijk(Ival = q.i-1, Jval = q.j-1, Kval = q.k)
        q.subq.append(q2)
        q.terms.append("nu")
    
        q3 = Qijk(Ival = q.i-1, Jval = q.j, Kval = q.k+1)
        q.subq.append(q3)
        q.terms.append("xi")
    elif(q.j > 1):
        q1 = Qijk(Ival = 0, Jval = q.j-2, Kval = q.k)
        q.subq.append(q1)
        q.terms.append("sigma")
        
        q2 = Qijk(Ival = 0, Jval = q.j-1, Kval = q.k-1)
        q.subq.append(q2)
        q.terms.append("rho" + str(q.j))
    else:
        q1 = Qijk(Ival = 0, Jval = 0, Kval = q.k)
        q.subq.append(q1)
        q.terms.append("ups")
        
        q2 = Qijk(Ival = 0, Jval = 0, Kval = q.k-1)
        q.subq.append(q2)
        q.terms.append("om")

    for i in range(len(q.subq)):
        unrol(q.subq[i])
    
    return

def collect(q, Q, term):
    if (q.i == 0 and q.j == 0):
        Q.bases[q.k-Q.start].append(term)
    else:
        for i in range(len(q.subq)):
            collect(q.subq[i], Q, term + q.terms[i] + ",")

def algebraic_unrol(i, j, k):
    q = Qijk(Ival = i, Jval = j, Kval = k)
    unrol(q)
    collect(q, q, "")
    q.simplify()
    return q


f = open('radial_gen.part2', 'w')
print("", file=f)
for j in range(MAX_UNROL_AM+1):
    for i in range(j+1):
        for k in range(1, 3*MAX_UNROL_AM+1-i-j):
            if (i + j + k) % 2 == 0: 
                q = algebraic_unrol(i, j, k)
                q.eliminate()
                q.sort()
                q.write_code(f)
                print("", file=f)
