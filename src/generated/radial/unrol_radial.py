# Script to turn the recurrence relations described in:
# [Shaw2017] = Shaw and Hill, JCP 147, 074108 (2017), doi: 10.1063/1.4986887
# for the radial integrals, into reduced form, and then into code

MAX_UNROL_AM = 4 # The maximum angular momentum to unrol up to

# DO NOT EDIT BELOW HERE

from sympy import *

class Qijk:
    """Representation of a term in the type2 ECP integral
       between gaussians on centres A and B
       Q_{ijk} = int_0^{infty} dr r^k exp(-pr^2) M_i(2Ar) M_j(2Br)
       and its reduction to the base integrals, F_N, G_N, H_N
       (see eqs 30 - 32 in [Shaw2017])
    
       Members:
       i, j, k - the integral indices
       size - number of terms in recurrence, 2*i+j-1
       start, end - the lowest, highest values of N in the recurrences
       subq - Qijk terms in the recurrence expansion
       terms, bases - coefficients and base integrals used
       f, ga, gb, h - all ocurrences of base integrals in recurrence
       
    """
    def __init__(self, Ival = 0, Jval = 0, Kval = 0):
        """Creates a term Qijk with i=Ival, j=Jval, k=Kval"""
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
        """Prints a list of the bases attached to the Qijk"""
        print("Q", self.i, self.j, self.k)
        for i in range(len(self.bases)):
            print(i+self.start, ":")
            for j in range(len(self.bases[i])):
                print(self.bases[i][j])
    
    def print_simple(self):
         """Prints algebraically simplified reps of the bases in Qijk"""
         print("Q", self.i, self.j, self.k)
         for i in range(len(self.bases)):
            print(i+self.start, ":", simplify(self.bases[i]))
            
    def print_fgh(self):
        """Prints simplified versions of the coefficients for each base
           integral in Qijk
        """
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
        """Writes C++ code in the appropriate format for this term
           to the file, f, given
        """
        # index in the switch-case block
        print("\t\t\t\t\t\t\t\t\tcase", self.i*10000+self.j*100+self.k, ": {", file=f)

        # print out F integral terms
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
                
        # Gb integral terms
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

        # Ga integral terms
        for i in range(len(self.ga)):
            ix = 2*i + self.start + 1
            if ix == 1:
                simp = simplify(self.ga[i])
                if simp != 0:
                    print("\t\t\t\t\t\t\t\t\t\tresult += (", simp, ") * G1A;", file=f)

        # H integrals terms
        for i in range(len(self.h)):
            ix = 2*i+self.start + 2
            if ix == 2:
                simp = simplify(self.h[i])
                if simp!= 0:
                    print("\t\t\t\t\t\t\t\t\t\tresult += (", simp, ") * H2;", file=f)
                    
        # close the case statement
        print("\t\t\t\t\t\t\t\t\t\tbreak;", file=f)
        print("\t\t\t\t\t\t\t\t\t}", file=f)

    def simplify(self):
        """Simplifies the algebraic expressions for each term 
           in the recurrence relation, using sympy
        """
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
        """Sorts all the bases found into base integral types
           F if even, Gb if odd
        """
        for i in range(len(self.bases)):
            if i % 2 == 0:
                self.f.append(self.bases[i])
            else:
                self.gb.append(self.bases[i])
                
    def eliminate(self):
        """Eliminates terms in the expanded recurrence expressions
           where N < 1, as described by eqns 37-39 in [Shaw2017]
        """
        x = Symbol('x') # k_A/2
        y = Symbol('y') # k_B/2
        z = Symbol('z') # placeholder 
        p = Symbol('p') # p
        z = 0
        
        # only need to eliminate base integrals with N < 1
        if self.start < 1:
            # append a zero term if max N < 1 
            if self.end < 1:
                self.bases.append(z)
            # and an additional zero if it is exactly 1
            if self.end == 1:
                self.bases.append(z)
                
            w = Symbol('w')  # placeholder
            N = self.start
            ix = 0    # start at lowest N and recur up
            gaix = -1 # counts where we're up to in the Ga expansion (eq 38)
            hix  = -1 # and the H expansion (eq 39)
            while (N < 1):
                # all base integrals share the first term
                # X_N: (2p/N-1) * X_{N+2}
                # so grab X_N and X_{N+2}
                z = self.bases[ix] 
                w = self.bases[ix+2]
                # then add the new coefficient to X_{N+2}
                w = w + (2 * p / (N-1))*z
                self.bases[ix+2] = w 
            
                # we also have X_{N+1} gaining a k_B/(N-1)
                w = self.bases[ix+1]
                w = w - (2 * y / (N-1))*z
                self.bases[ix+1] = w
                
                # if even, we add a Ga then H
                if ix % 2 == 0:
                    # is it the first Ga term?
                    if gaix > -1:
                        w = self.ga[gaix]
                        w = w - (2*x / (N-1))*z
                        self.ga[gaix] = w
                    else:
                        w = -(2*x / (N-1))*z
                        self.ga.append(w)
                        gaix += 1
                    
                    # is it the first H term?    
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
                    # otherwise H then Ga
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
                
                # recur upwards
                N += 1
                ix += 1

def parse(term):
    """Parses a string representation of the form of
       eqn. 28 in [Shaw2017] into the symbolic representation
       used by sympy
    """
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
    """Recursively unrols the recurrence relations described in
       [Shaw2017] equations 28 - 33, for the given Qijk objects
    """
    if (q.i == 0 and q.j == 0):
        # already all base integrals
        return
    elif (q.i > 0):
        # need to reduce i using eqn 28 of [Shaw2017]
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
        # need to reduce j using eqn 29 of [Shaw2017]
        q1 = Qijk(Ival = 0, Jval = q.j-2, Kval = q.k)
        q.subq.append(q1)
        q.terms.append("sigma")
        
        q2 = Qijk(Ival = 0, Jval = q.j-1, Kval = q.k-1)
        q.subq.append(q2)
        q.terms.append("rho" + str(q.j))
    else:
        # equation 33 of [Shaw2017]
        q1 = Qijk(Ival = 0, Jval = 0, Kval = q.k)
        q.subq.append(q1)
        q.terms.append("ups")
        
        q2 = Qijk(Ival = 0, Jval = 0, Kval = q.k-1)
        q.subq.append(q2)
        q.terms.append("om")
        
    # we then recurseively unrol all the Qijks just added
    # to this one
    for i in range(len(q.subq)):
        unrol(q.subq[i])
    
    return

def collect(q, Q, term):
    """Recursively joins all the coefficients from unrol(q)
       into Q, forming a string of the algebraic expression.
    """
    if (q.i == 0 and q.j == 0):
        # everything is in terms of base integrals
        Q.bases[q.k-Q.start].append(term)
    else:
        # we need to call again for each Qijk contained in
        # the original Qijk
        for i in range(len(q.subq)):
            collect(q.subq[i], Q, term + q.terms[i] + ",")

def algebraic_unrol(i, j, k):
    """Unrols the algebraic expression for the integral Q_ijk
       and returns the simplified Qijk object
    """
    q = Qijk(Ival = i, Jval = j, Kval = k)
    unrol(q)
    collect(q, q, "")
    q.simplify()
    return q


if __name__ == "__main__":
    f = open('radial_gen.part2', 'w')
    print("", file=f)
    # unrol over all desired ijk
    for j in range(MAX_UNROL_AM+1):
        for i in range(j+1):
            for k in range(1, 3*MAX_UNROL_AM+1-i-j):
                if (i + j + k) % 2 == 0: 
                    # generate the algebraic expression
                    q = algebraic_unrol(i, j, k)
                    # eliminate any non-base integrals
                    q.eliminate()
                    # sort terms by N
                    q.sort()
                    # write code to file
                    q.write_code(f)
                    # blank spacer before next switch case
                    print("", file=f)
