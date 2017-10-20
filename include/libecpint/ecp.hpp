/* 
 *      Copyright (c) 2017 Robert Shaw
 *		This file is a part of Libecpint.
 *
 *      Permission is hereby granted, free of charge, to any person obtaining
 *      a copy of this software and associated documentation files (the
 *      "Software"), to deal in the Software without restriction, including
 *      without limitation the rights to use, copy, modify, merge, publish,
 *      distribute, sublicense, and/or sell copies of the Software, and to
 *      permit persons to whom the Software is furnished to do so, subject to
 *      the following conditions:
 *
 *      The above copyright notice and this permission notice shall be
 *      included in all copies or substantial portions of the Software.
 *
 *      THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *      EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *      MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *      NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 *      LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *      OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 *      WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef ECP_HEAD
#define ECP_HEAD

#include <vector>
#include <array>
#include <map>

namespace libecpint {
	
	/** \struct GaussianECP
	  * \brief  Describes a Gaussian of angular momentum l of the form d r^n e^{-ax^2}
	  */
	struct GaussianECP {
		int n; ///< Power of r
		int l; ///< Angular momentum
		double a; ///< Exponent
		double d; ///< Coefficient
	
		/// Default constructor, sets n=l=d=0, a=1
		GaussianECP(); 
		
		/** 
		  * Constructs a new GaussianECP.
		  * @param n - power of r
		  * @param l - angular momentum
		  * @param a - exponent
		  * @param d - coefficient
		  */ 
		GaussianECP(int n, int l, double a, double d);
		
		/// Copy constructor
		GaussianECP(const GaussianECP& other);
	};

 	/** \struct ECP
	  * \brief  Stores the details of an ECP expanded in terms of Gaussians and spherical harmonics.
	  */
	struct ECP {
		std::vector<GaussianECP> gaussians; ///< All the primitives in the ECP expansion
		int N; ///< Number of Gaussians
		int L; ///< Maximum angular momentum
		
		std::array<double, 3> center_; ///< xyz coordinates of the atom on which the ECP is located
	
		/// Constructs an empty ECP (N = L = 0, center_ = {0, 0, 0})
		ECP();
		
		/**
		  * Constructs an ECP at the given position
		  * @param _center - xyz coordinates of the ECP
		  */ 
		ECP(const double *_center);
		
		/// Copy constructor
		ECP(const ECP &other);
	
		/** 
		  * Adds a new GaussianECP to the ECP
		  * @param n - power of r
		  * @param l - angular momentum
		  * @param a - exponent
		  * @param d - coefficient
		  * @param needSort - true = the GaussianECPs are sorted (if done once at the end, speeds up evaluation)
		  */ 
		void addPrimitive(int n, int l, double a, double d, bool needSort = true);
		
		/// @return the xyz coordinates of the ECP
		const double* center() const { return &center_[0]; }
		
		//// Sets the xyz coordinates of the ECP
		void setPos(double x, double y, double z);
		
		/// Sort primitives according to angular momentum
		void sort(); 
		
		/**
		  * @param i - the index of GaussianECP required
		  * @return a reference to the ith GaussianECP
		  */ 
		GaussianECP& getGaussian(int i) { return gaussians[i]; }
		
		/// @return the number of primitives in ECP
		int getN() const { return N; }
		
		/// @return true if the highest angular momentum functions have zero coefficients (e.g. Stuttgart-Dresden ECPs)
	    bool noType1() const; 
	
		/** 
		  * Evaluates the ECP at a given distance for a given angular momentum shell.
		  * @param r - the radius at which to evaluate
		  * @param l - the angular momentum shell to evaluate over
		  * @return the value of the l-th angular momentum shell of the ECP at radius r
		  */ 
		double evaluate(double r, int l);
  
		/// @return the maximum angular momentum in the ECP
		int getL() const { return L; }
	
	};

	/** \struct ECPBasis
	  * \brief A lightweight container for a basis set of ECP objects
	  */
	class ECPBasis {
	private:
		std::vector<ECP> basis;    ///< Vector of ECPs
		std::vector<int> atomList; ///< List of atom indices corresponding to the ECPs in basis
		int N; ///< Number of ECPs in basis
		int maxL; ///< Maximum angular momentum in the ECP basis
	
	public:
		/// Constructs an empty ECPBasis (N = maxL = 0)
		ECPBasis(); 
		
		/// A map of atomic number to the number of electrons in core of ECP for that atom
		std::map<int, int> core_electrons;
		
		/** 
		  * Adds an ECP to basis.
		  * @param U - the ECP to be added
		  * @param atom - the index of that atom on which U is located
		  */
		void addECP(ECP &U, int atom);
		
		/** 
		  * @param i - index of ECP required
		  * @return a reference to the ith ECP in basis
		  */ 
		ECP& getECP(int i);
		
		/**
		  * @param q - an atomic number 
		  * @return the number of electrons in core of ECP for the atom with atomic number q, if defined, otherwise zero
		  */
		int getECPCore(int q); 
		
		/**
		  * @param i - the index of ECP of interest
		  * @return the index of the atom on which the ith ECP is located
		  */
		int getAtom(int i) { return atomList[i]; }
		
		/// @return the maximum angular momentum GaussianECP in the entire ECP basis
		int getMaxL() const { return maxL; }
		
		/// @return the number of ECPs in basis
		int getN() const { return N; }
	};

}

#endif
