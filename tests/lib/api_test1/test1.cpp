#include "api.hpp"
#include "testutil.hpp"
#include <iostream>

double R = 2.66;
double g_coords[51]  = {0.0, 0.0, 0.0, 
					  0.0, 0.0, 0.0,
					  0.0, 0.0, 0.0,
					  0.0, 0.0, 0.0,
					  0.0, 0.0, 0.0,
					  0.0, 0.0, 0.0,
					  0.0, 0.0, 0.0,
					  0.0, 0.0, 0.0,
					  0.0, 0.0, 0.0,
					  0.0, 0.0, 0.0,
					  0.0, 0.0, 0.0,
					  0.0, 0.0, 0.0,
					  0.0, 0.0, R,
					  0.0, 0.0, R,
					  0.0, 0.0, R,
					  0.0, 0.0, R,
					  0.0, 0.0, R};
				  	
double g_exps[52]    = {2.808600E+03,4.211800E+02,5.034570E+01,1.791330E+01,3.805310E+00,1.749680E+00,4.485550E-01,1.644980E-01,
					  2.808600E+03,4.211800E+02,5.034570E+01,1.791330E+01,3.805310E+00,1.749680E+00,4.485550E-01,1.644980E-01,
				  	  4.485550E-01,1.644980E-01,5.020000E-02,
				  	  1.057520E+02,2.763680E+01,6.596560E+00,2.785220E+00,1.078120E+00,3.935370E-01,1.274690E-01,
				  	  1.057520E+02,2.763680E+01,6.596560E+00,2.785220E+00,1.078120E+00,3.935370E-01,1.274690E-01, 
				      1.274690E-01,3.940000E-02,
				  	  1.438650E+02,4.611630E+01,1.736940E+01,6.951070E+00,2.756070E+00,1.011780E+00,4.291000E-01, 
				  	  4.291000E-01,1.548000E-01,
					  1.301000E+01,1.962000E+00,4.446000E-01,1.220000E-01,
					  1.220000E-01,2.974000E-02,
					  7.270000E-01,1.410000E-01};
					  
double g_coefs[52]   = {1.606000E-03,8.393000E-03,6.957800E-02,-3.899080E-01,6.944970E-01,4.913540E-01,2.263700E-02,-3.723000E-03, 
					  -6.350000E-04,-3.492000E-03,-2.519500E-02,1.501130E-01,-3.662260E-01,-3.834220E-01,7.144680E-01,5.352530E-01, 
				  	  1.0, 1.0, 1.0, 
				  	  5.341000E-03,-8.308400E-02,4.477660E-01,5.506170E-01,1.235000E-01,-3.771000E-03,2.278000E-03, 
				      -1.308000E-03,2.292100E-02,-1.450290E-01,-2.090370E-01,9.373000E-02,6.050210E-01,4.571230E-01, 
				      1.0, 1.0, 
				      1.023700E-02,7.608300E-02,2.298070E-01,4.033470E-01,4.097280E-01,1.627900E-01,6.991000E-03, 
					  1.0, 1.0, 
					  1.968500E-02,1.379770E-01,4.781480E-01,5.012400E-01,
					  1.0, 1.0,
					  1.0, 1.0};
					  
int    g_ams[17]     = {0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 0, 0, 0, 1, 1};
int    g_lengths[17] = {8, 8, 1, 1, 1, 7, 7, 1, 1, 7, 1, 1, 4, 1, 1, 1, 1};

double u_coords[3]  = {0.0, 0.0, 0.0};
					  
double u_exps[16]    = {1.0,
					  70.024257, 31.178412, 7.156593,
				  	  46.773471, 46.184120, 21.713858, 20.941792,
				  	  50.698839, 50.644764, 15.447509, 15.500259, 2.800391, 1.077480,
				      14.465606, 21.234065};
double u_coefs[16]   = {0.0,
					  49.962834, 370.014205, 10.241439,
				      99.112244, 198.253046, 28.261740, 56.623366,
				      -18.605853, -27.923280, -0.379693, -0.780583, 0.035968, 0.094397,
				      -1.091269, -2.887691};
int    u_ams[16]     = {4, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3};
int    u_ns[16]      = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
int    u_lengths[1] = {16};

int main(int argc, char* argv[]) {
	using namespace libecpint; 
	
	ECPIntegrator factory;
	factory.set_gaussian_basis(17, g_coords, g_exps, g_coefs, g_ams, g_lengths);
	factory.set_ecp_basis(1, u_coords, u_exps, u_coefs, u_ams, u_ns, u_lengths);
	factory.init();
	factory.compute_integrals();
	
	std::vector<double> flat_result;
	std::shared_ptr<std::vector<double>> ints = factory.get_integrals();
	for (int i = 0; i < ints->size(); i++) 
		flat_result.push_back((*ints)[i]);
	
	//for (auto& v : flat_result) std::cout << std::setprecision(15) << v << std::endl;
	
	return check_file<double>("api_test1.output", flat_result, 1e-5, 1e-8);
}
