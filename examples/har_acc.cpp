// VolEsti

// Copyright (c) 2012-2017 Vissarion Fisikopoulos

// VolEsti is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// VolEsti is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with HeaDDaCHe,
// see <http://www.gnu.org/licenses/>.

#include <cstdlib>
#include <fstream>
#include <vol_rand.h>
#include <rounding.h>
#include <cmath>
#include <string>
#include <boost/program_options.hpp>


#include <chrono>

namespace funcs2{
class Timer2 {
public:
	Timer2() { start_time = std::chrono::high_resolution_clock::now(); }

	void start() {
		start_time = std::chrono::high_resolution_clock::now();
	}

	double end() {
		return elapsed_seconds();
	}

	double elapsed_seconds() {
		auto end_time = std::chrono::high_resolution_clock::now();
		auto elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
		return elapsed.count();
	}
private:
	decltype(std::chrono::high_resolution_clock::now()) start_time;
};
}

namespace po = boost::program_options; 

void completeTests(stdHPolytope<double>* P, vars& var, int nqp, int k, int l, int algoType) {
	var.verbose = true;
	Point* internalPoint = new Point(P->dimension(), CGAL::ORIGIN);
	Point chebPoint = P->create_point_representation(var, internalPoint);
	if (algoType==USE_LSH) {
		P->create_lsh_ds(k, l);
	} 
	else if (algoType==USE_LSHBOX) {
		P->create_lshbox(k, l);
	}

	CGAL::Random_points_on_sphere_d<Point> rps(P->dimension(), 1);
	funcs2::Timer2 timer;
	double timeit = 0;
	var.walk_steps = 200;
	for (int i=0; i<nqp; ++i) {
	    std::list<Point> randPoints;
	    rand_point_generator(*P, chebPoint, 1, var.walk_steps, randPoints, var);
		Vector direction = (*rps) - CGAL::ORIGIN;	

		auto it = randPoints.begin();
		Ray r((*it), direction);
		int numberOfSteps = 0;
		bool succeeded = false;
		timer.start();

		auto actualPoints = P->line_intersect(r.source(), r.direction().vector(), false);

		timer.start();
		Point appxPoint = P->compute_boundary_intersection(r, &numberOfSteps, &succeeded, 0.1, algoType, var, var.walk_steps, l);
		timeit += timer.elapsed_seconds();

		direction *= -1;
		timer.start();
		Point appxPoint2 = P->compute_boundary_intersection(r, &numberOfSteps, &succeeded, 0.1, algoType, var, var.walk_steps, l);

    	RNGType &rng = var.rng;
    	boost::random::uniform_real_distribution<> &urdist = var.urdist;
    	double lambda = urdist(rng);
    	Point newPoint = CGAL::Origin() + (lambda*(appxPoint-CGAL::ORIGIN) + (1-lambda)*(appxPoint2-CGAL::ORIGIN));

		++rps;
	}
	timeit /= nqp;
	//std::cout << timeit << "s/q" << std::endl;

	delete internalPoint;
}

void simpleTests(stdHPolytope<double>* P, vars& var, int nqp, bool exact=true, int k=-1, int l=-1) {
	Point* internalPoint = new Point(P->dimension(), CGAL::ORIGIN);
	Point chebPoint = P->create_point_representation(var, internalPoint);
	if (!exact) {
		P->create_lsh_ds(k, l);
	}

	CGAL::Random_points_on_sphere_d<Point> rps(P->dimension(), 1);
	// sample nqp points with CDHR

	int totalSteps = 0;
	int maxSteps = -1;
	int minSteps = 101;
	int failed = 0;
	auto algo = exact?USE_EXACT:USE_LSH;
	for (int i=0; i<nqp; ++i) {
	    std::list<Point> randPoints;
	    rand_point_generator(*P, chebPoint, 1, var.walk_steps, randPoints, var);
		Vector direction = (*rps) - CGAL::ORIGIN;	

		auto it = randPoints.begin();
		Ray r((*it), direction);
		int numberOfSteps = 0;
		bool succeeded = false;
		
		P->compute_boundary_intersection(r, &numberOfSteps, &succeeded, 0.1, algo, var, var.walk_steps, l);
	}
	delete internalPoint;
}

int main(int argc, char* argv[]) {
	// parse args
	po::options_description desc{"Params"};
	int n;
	int d; 
	int diameter;
	int nqp;
	bool exact = true;
	int k;
	int l;
	int probes;
	float w;
	int algoType = USE_EXACT;
	bool computeVol;
	std::string filename;

	desc.add_options()
		("help,h", "Help message")
		(",n", po::value<int>(&n), "Number of points")
		(",d", po::value<int>(&d), "Dimension")
		("diameter", po::value<int>(&diameter)->default_value(1000), "Diameter of sphere")
		("nqp", po::value<int>(&nqp)->default_value(100), "Number of query points")
		(",k", po::value<int>(&k)->default_value(-1), "k for LSH")
		(",L", po::value<int>(&l)->default_value(-1), "L for LSH")
		("probes", po::value<int>(&probes)->default_value(-1), "probes for LSH")
		(",w", po::value<float>(&w)->default_value(4.0f), "Window of euclidian LSH hash functions")
		("algo,a", po::value<int>(&algoType)->default_value(USE_EXACT), "0 -- FALCONN\n2 -- LSHBOX\n3 -- EXACT")
		("f1", po::value<std::string>(&filename), "Filename to use instead of creating random")
		("volume,v", po::bool_switch(&computeVol)->default_value(false), "Compute volume");

	po::variables_map vm; 
	po::store(po::parse_command_line(argc, argv, desc),  vm);
	po::notify(vm);
	//std::cout << "N = " << n << "\nd = " << d << "\ndiam = " << diameter << "\nnqp = " << nqp << "\nk = " << k << "\nL = " << l << std::endl;
	if (vm.count("help")) {
		std::cout << desc << std::endl;
		return 0;
	} 
	if (k>-1) {
		exact = false;
	}
	// init vars
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	RNGType rng(seed);
	boost::normal_distribution<> rdist(0,1);
	boost::variate_generator< RNGType, boost::normal_distribution<> > get_snd_rand(rng, rdist);
	boost::random::uniform_real_distribution<>(urdist);
	boost::random::uniform_real_distribution<> urdist1(-1,1);

	// create polytope with internal repr
	stdHPolytope<double>* P; //randomPolytope<double>(n, d, diameter);
	if (vm.count("f1")) {
		std::ifstream inp;
		std::vector<std::vector<double> > Pin;
		inp.open(filename.c_str(),std::ifstream::in);
		read_pointset(inp,Pin);
		std::cout<<"Read input from file..."<<std::endl;
		P = new stdHPolytope<double>();
		P->init(Pin);
		n = P->num_of_hyperplanes();
		d = P->dimension();
	}
	else {
		P = randomPolytope<double>(n, d, diameter);
	}

   	//int rnum = std::pow(1,-2) * 400 * n * std::log(n);
   	int rnum = std::pow(1,-2) * 40 * n * std::log(n);

    vars var(rnum,d,40,1,0.0000001,0,0,0,0,rng,get_snd_rand,
                 urdist,urdist1,true,false,false,false,false,false,algoType,0.1);
	var.verbose = true;
	double Chebtime;
	
	if (!computeVol) {
		//results["original"] = completeTests(P, var, nqp, k, l, algoType);
	}
	else {
		funcs::Timer2 timer;
		
		Point* internalPoint = new Point(P->dimension(), CGAL::ORIGIN);
		Point chebPoint = P->create_point_representation(var, internalPoint);
		if (algoType==USE_LSH) {
			P->create_lsh_ds(k, l);
		} 
		else if (algoType==USE_LSHBOX) {
			P->create_lshbox(k, l, w);
		}
		std::cout << "Using line_intersect from har_acc" << std::endl;
		var.algoType=USE_EXACT;
		timer.start();
        volume1_reuse2(*P,var,var,Chebtime);
		//results["original"]["time"] = timer.elapsed_seconds();
		timer.elapsed_seconds();
		var.algoType=algoType;
		std::cout << "Using oracle from har_acc" << std::endl;
		timer.start();
        volume1_reuse2(*P,var,var,Chebtime);
		timer.elapsed_seconds();
		//results["original"]["oracle_time"] = timer.elapsed_seconds();
	}

	randomTransformation<stdHPolytope<double> >(P);	

	if (!computeVol) {
		completeTests(P, var, nqp, k, l, algoType);
	}
	else {
		funcs::Timer2 timer;
		Point* internalPoint = new Point(P->dimension(), CGAL::ORIGIN);
		Point chebPoint = P->create_point_representation(var, internalPoint);
		if (algoType==USE_LSH) {
			P->create_lsh_ds(k, l);
		} 
		else if (algoType==USE_LSHBOX) {
			P->create_lshbox(k, l, w);
		}
		var.algoType=USE_EXACT;
		std::cout << "Using line_intersect from har_acc" << std::endl;
		timer.start();
        volume1_reuse2(*P,var,var,Chebtime);
		timer.elapsed_seconds();
		var.algoType=algoType;
		std::cout << "Using oracle from har_acc" << std::endl;
		timer.start();
        volume1_reuse2(*P,var,var,Chebtime);
		timer.elapsed_seconds();
	}

	delete P;
	//std::cout << double(counter)*100/nqp << "\n" << double(counter2)*100/nqp << std::endl;
}
