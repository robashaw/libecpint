
.. _program_listing_file__Users_robertshaw_devfiles_libecpint_src_generated_radial_unrol_radial.py:

Program Listing for File unrol_radial.py
========================================

|exhale_lsh| :ref:`Return to documentation for file <file__Users_robertshaw_devfiles_libecpint_src_generated_radial_unrol_radial.py>` (``/Users/robertshaw/devfiles/libecpint/src/generated/radial/unrol_radial.py``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: py

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
           self.ii = Ival
           self.jj = Jval
           self.kk = Kval
           self.sizesize = 2*self.ii + self.jj + 1 
           self.startstart = self.kk-self.ii-self.jj
           self.endend = self.kk+self.ii
           self.subqsubq = []
           self.termsterms = []
           self.basesbases = []
           self.ff = []
           self.gaga = []
           self.gbgb = []
           self.hh = []
           for i in range(self.sizesize):
               self.basesbases.append([])
       
       def print(self):
           """Prints a list of the bases attached to the Qijk"""
           print("Q", self.ii, self.jj, self.kk)
           for i in range(len(self.basesbases)):
               print(i+self.startstart, ":")
               for j in range(len(self.basesbases[i])):
                   print(self.basesbases[i][j])
       
       def print_simple(self):
            """Prints algebraically simplified reps of the bases in Qijk"""
            print("Q", self.ii, self.jj, self.kk)
            for i in range(len(self.basesbases)):
               print(i+self.startstart, ":", simplify(self.basesbases[i]))
               
       def print_fgh(self):
           """Prints simplified versions of the coefficients for each base
              integral in Qijk
           """
           print("Q", self.ii, self.jj, self.kk)
           for i in range(len(self.ff)):
               ix = 2*i+self.startstart
               if ix > 0:
                   print("F", ix, ":", simplify(self.ff[i]))
           for i in range(len(self.gbgb)):
               ix = 2*i + self.startstart + 1
               if ix > 0:
                   print("GB", ix, ":", simplify(self.gbgb[i]))
           for i in range(len(self.gaga)):
               ix = 2*i + self.startstart + 1
               if ix > 0:
                   print("GA", ix, ":", simplify(self.gaga[i]))
           for i in range(len(self.hh)):
               ix = 2*i+self.startstart + 2
               if ix > 0:
                   print("H", ix, ":", simplify(self.hh[i]))
       
       def write_code(self, f):
           """Writes C++ code in the appropriate format for this term
              to the file, f, given
           """
           # index in the switch-case block
           print("\t\t\t\t\t\t\t\t\tcase", self.ii*10000+self.jj*100+self.kk, ": {", file=f)
   
           # print out F integral terms
           for i in range(len(self.ff)):
               ix = 2*i+self.startstart
               if ix == 2:
                   simp = simplify(self.ff[i])
                   if simp!= 0:
                       print("\t\t\t\t\t\t\t\t\t\tresult = (", simp, ") * values[0];", file=f)
               elif ix > 2:
                   simp = simplify(self.ff[i])
                   if simp!= 0:
                       print("\t\t\t\t\t\t\t\t\t\tresult += (", simp, ") * values[", ix-2, "];", file=f)
                   
           # Gb integral terms
           for i in range(len(self.gbgb)):
               ix = 2*i + self.startstart + 1
               if ix == 1:
                   simp = simplify(self.gbgb[i])
                   if simp!= 0:
                       print("\t\t\t\t\t\t\t\t\t\tresult += (", simp, ") * G1B;", file=f)
               elif ix > 1:
                   simp = simplify(self.gbgb[i])
                   if simp!= 0:
                       print("\t\t\t\t\t\t\t\t\t\tresult += (", simp, ") * values[", ix-2, "];", file=f)
   
           # Ga integral terms
           for i in range(len(self.gaga)):
               ix = 2*i + self.startstart + 1
               if ix == 1:
                   simp = simplify(self.gaga[i])
                   if simp != 0:
                       print("\t\t\t\t\t\t\t\t\t\tresult += (", simp, ") * G1A;", file=f)
   
           # H integrals terms
           for i in range(len(self.hh)):
               ix = 2*i+self.startstart + 2
               if ix == 2:
                   simp = simplify(self.hh[i])
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
           for i in range(len(self.basesbases)):
               x = Symbol('x')
               y = Symbol('y')
               z = Symbol('z')
               z = 0
               for j in range(len(self.basesbases[i])):
                   z = z + parse(self.basesbases[i][j])
               simple_bases.append(z)
           self.basesbases = simple_bases
           
       def sort(self):
           """Sorts all the bases found into base integral types
              F if even, Gb if odd
           """
           for i in range(len(self.basesbases)):
               if i % 2 == 0:
                   self.ff.append(self.basesbases[i])
               else:
                   self.gbgb.append(self.basesbases[i])
                   
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
           if self.startstart < 1:
               # append a zero term if max N < 1 
               if self.endend < 1:
                   self.basesbases.append(z)
               # and an additional zero if it is exactly 1
               if self.endend == 1:
                   self.basesbases.append(z)
                   
               w = Symbol('w')  # placeholder
               N = self.startstart
               ix = 0    # start at lowest N and recur up
               gaix = -1 # counts where we're up to in the Ga expansion (eq 38)
               hix  = -1 # and the H expansion (eq 39)
               while (N < 1):
                   # all base integrals share the first term
                   # X_N: (2p/N-1) * X_{N+2}
                   # so grab X_N and X_{N+2}
                   z = self.basesbases[ix] 
                   w = self.basesbases[ix+2]
                   # then add the new coefficient to X_{N+2}
                   w = w + (2 * p / (N-1))*z
                   self.basesbases[ix+2] = w 
               
                   # we also have X_{N+1} gaining a k_B/(N-1)
                   w = self.basesbases[ix+1]
                   w = w - (2 * y / (N-1))*z
                   self.basesbases[ix+1] = w
                   
                   # if even, we add a Ga then H
                   if ix % 2 == 0:
                       # is it the first Ga term?
                       if gaix > -1:
                           w = self.gaga[gaix]
                           w = w - (2*x / (N-1))*z
                           self.gaga[gaix] = w
                       else:
                           w = -(2*x / (N-1))*z
                           self.gaga.append(w)
                           gaix += 1
                       
                       # is it the first H term?    
                       if hix > -1:
                           z = self.hh[hix]
                           w = (2 * p / (N-1))*z
                           self.hh.append(w)
                           hix += 1
                           
                           w = self.gaga[gaix]
                           w = w - (2 * y / (N-1))*z
                           self.gaga[gaix] = w
                           
                           w = self.basesbases[ix+1]
                           w = w - (2 * x / (N-1))*z
                           self.basesbases[ix+1] = w
                   else:
                       # otherwise H then Ga
                       if hix > -1:
                           w = self.hh[hix]
                           w = w - (2*x / (N-1))*z
                           self.hh[hix] = w
                       else:
                           w = - (2*x / (N-1))*z
                           self.hh.append(w)
                           hix += 1
                       
                       if gaix > -1:
                           z = self.gaga[gaix]
                           w = (2 * p /(N-1)) * z
                           self.gaga.append(w)
                           gaix += 1
                           
                           w = self.hh[hix]
                           w = w - (2 * y / (N-1))*z
                           self.hh[hix] = w
                           
                           w = self.basesbases[ix+1]
                           w = w - (2*x / (N-1))*z
                           self.basesbases[ix+1] = w
                   
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
