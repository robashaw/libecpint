---
title: 'libecpint: A C++ library for  the efficient evaluation of integrals over effective core potentials'
tags:
  - C++
  - computational chemistry
authors:
  - name: Robert A. Shaw
    orcid: 0000-0002-9977-0835
    affiliation: 1
  - name: J. Grant Hill
    orcid: 0000-0002-6457-5837
    affiliation: 1
affiliations:
 - name: Department of Chemistry, University of Sheffield, Sheffield, S3 7HF
   index: 1
date: January 27th, 2020
bibliography: paper.bib
---

# Summary
Effective core potentials (ECPs) are widely-used in computational chemistry both to reduce the computational cost of calculations,[@Dolg2000] and include relevant physics that would not otherwise be present.[@Dolg2002] In particular, for heavy main-group atoms [@Wadt1985] and transition metals,[@Hay1985] the number of core electrons greatly outnumbers the number of valence electrons. It is generally considered that these will not play a significant role in chemical reactivity, and thus can be `frozen`. Moreover, these electrons show significant relativistic character.[@Dolg2002, @Dolg2012] Both of these issues can be resolved with the introduction of an `effective core`, represented as a fixed electronic potential. This potential is typically represented as a linear combination of gaussians of varying angular momenta.[@Dolg2000]

The introduction of an ECP results in an additional term in the core  Hamiltonian, over which new electronic integrals must be computed. These three-center integrals are far from trivial, and they cannot in general be treated the same way as other electronic integrals.[@McMurchie1981, @FloresMoreno2006] Several widely used computational chemistry codes lack the ability to calculate these integrals due to the difficulty involved in their computation. The present library, `libecpint`, provides an open-source solution to this. It is a standalone library written in modern C++ capable of the highly efficient computation of integrals over ECPs with gaussian orbitals of arbitrary angular momentum, along with their first and second geometric derivatives. The methods implemented are based on novel algorithms that use automatic code generation and symbolic simplification of recursive expressions, along with highly optimised Gauss-Chebyshev quadrature.

# Statement of need

Effective core potentials are an essential part of modern computational chemistry. However, existing implementations are typically unavailable or inaccessible for free use by the open source community. Commonly used proprietary software, such as Gaussian [@Gaussian16] or Molpro,[@MOLPRO] do not make details of their implementations available, while the few open-source computational chemistry packages either do not include ECP functionality or use outdated implementations that would not be compatible with modern codebases. A notable example of this is the widely-used Psi4 package,[@Psi4] in which a rudimentary version of `libecpint` was originally implemented. Prior to this, the inclusion of ECPs was one of the most requested features by the user base.

Additionally, there has been a recent reconnaissance in the development of efficient algorithms for evaluating ECP integrals. In particular, multiple research groups have outlined new approaches to prescreening integrals,[@Song2015, @Shaw2017, @McKenzie2018] greatly reducing the computational expense. The `libecpint` library implements many of these new algorithms, combining the recursive methods and fine-grained screening of Shaw et al. [@Shaw2017] with the higher-level screening of other recent work.[@Song2015, @McKenzie2018] The only known implementations of the latter papers are otherwise only available in proprietary software. Therefore `libecpint` represents a necessary contribution to the wider open-source computational chemistry community. It has already been adopted by multiple packages, including Entos QCore and QCSerenity, and will be part of a future release of Psi4.

# Functionalities

The core functionality of `libecpint` is the evaluation of both type 1 and type 2 integrals over ECP integrals parametrised in terms of contracted sets of primitive gaussians, as described in Shaw et al.[@Shaw2017] The component parts divide into the following functionalities:
- a built-in library of parametrised ECPs, with generic containers for Gaussian-type ECPs;
- a highly-optimised Bessel function evaluation routine;
- screening of ECP integrals over shell pairs of orbital basis functions,[@Song2015, @McKenzie2018] across all ECPs in a system;
- fine-grained screening of the individual type 2 integrals over primitive gaussians; [@Shaw2017]
- recursive, automatically-generated radial integral code; [@Shaw2017]
- adaptive quadrature for integrals not covered by the recursive routines; [@FloresMoreno2006]
- first- and second-order geometric derivatives of ECP integrals over shell pairs of gaussians.

These features can be accessed via two levels of API:
- a high-level interface, where the user provides a molecular geometry and basis set, then calls routines that return the full tensor of integrals (or integral derivatives) across all shell pairs;
- a low-level interface, where the user provides parameters and deals with the primitive integrals directly.
This allows a great deal of flexibility for different use cases, potentially allowing for users to further develop or adapt the routines themselves.

All primitive integral routines have been designed to be thread safe, allowing users to readily parallelise their calculations. In addition, there is a built in testing and benchmarking suite, allowing for efficiency comparisons both with other codes, and when testing new algorithmic developments.  

# Acknowledgements

Thank you to Moritz Bensberg, Peter Bygraves, Thomas Dresselhaus, Christopher Junghans, Peter Kraus, Jan Unsleber, and Jens Wehner, for finding bugs and suggesting improvements to the initial release of `libecpint`, and often providing helpful solutions.

# References
