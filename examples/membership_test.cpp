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
#include <vol_rand.h>
#include <rounding.h>
#include <string>
#include <boost/program_options.hpp>
#include <random>


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

void create_random_external_points(stdHPolytope<double>* P, std::vector<Point>& points, std::vector<std::pair<Point, std::pair<int, double> > >& randomPoints, double epsilon) {
	//unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
	auto seed = rd();
    std::mt19937 gen(seed); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(epsilon, 50.0);
    std::uniform_int_distribution<> urnd(0, P->num_of_hyperplanes()-1);
	for (int i=0; i<points.size(); i++) {
		Point pz;
		int randomFacet;
		Hyperplane facet;
		while (true) {
			randomFacet = urnd(gen);
			facet = P->get_hyperplane(randomFacet);
			Point tmp = P->project(points[i], randomFacet);

			if (P->is_in(tmp)==-1) {
				pz = tmp;
				break;
			}
		}

		double dist = dis(gen);
		Vector dir = facet.orthogonal_vector();
		dir *= dist/std::sqrt(dir.squared_length());
		Point newPoint = CGAL::ORIGIN + ((pz-CGAL::ORIGIN) + dir);
		randomPoints.push_back(std::make_pair(newPoint, 
				std::make_pair(randomFacet, dist)));
	}
}

double pointsDist(Point const& a, Point const& b, bool squared=true) {
	double dist = ((a-CGAL::ORIGIN)-(b-CGAL::ORIGIN)).squared_length();
	return squared?dist:std::sqrt(dist);
}	

void completeTests(stdHPolytope<double>* P, vars& var, int nqp, int k, int l, int box_k, int box_l, int algoType, float w, int numProbes) {
	//var.verbose = true;
	Point* internalPoint = new Point(P->dimension(), CGAL::ORIGIN);
	Point chebPoint = P->create_point_representation(var, internalPoint);
	P->create_lsh_ds(k, l);
	P->create_lshbox(box_k, box_l, w);

	std::vector<Point> points;
	var.algoType = USE_EXACT;
	var.coordinate = true;
	rand_point_generator(*P, *internalPoint, nqp, P->dimension(), points, var);
	var.algoType = algoType;

	funcs::Timer2 timer;
	std::vector<double> isInTime;
	std::vector<bool> isInSucceeded;
	std::vector<double> exactTime;
	std::vector<double> falconnTime;
	std::vector<double> boxTime;

	std::vector<double> exactDists;
	std::vector<double> falconnDists;
	std::vector<double> boxDists;

	std::vector<bool> exactSucceeded;
	std::vector<bool> falconnSucceeded;
	std::vector<bool> boxSucceeded;

	std::vector<double> exactANNs;
	std::vector<double> falconnANNs;
	std::vector<double> boxANNs;

	int nnIndex;
	for (int i=0; i<nqp; i++) {
		Point p = points[i];
		double minDist = pointsDist(P->project(p, 0), p);
		for (int j=1; j<P->num_of_hyperplanes(); j++) {
			double tmpDist = pointsDist(P->project(p, j), points[i]);
			if (tmpDist<minDist) {
				minDist=tmpDist;
			}
		}

		timer.start();
		int asd= P->is_in(p);
		isInSucceeded.push_back(asd==-1);
		double isIn = timer.elapsed_seconds();
		isInTime.push_back(isIn);

		timer.start();
		bool lshContains = P->contains_point_lsh(p, numProbes, &nnIndex);
		falconnTime.push_back(timer.elapsed_seconds());
		falconnSucceeded.push_back(lshContains);
		falconnDists.push_back(lshContains?0.0:minDist);
		if (lshContains) {
			falconnANNs.push_back(0.0);
		}
		else {
			Point nn = P->get_site(P->num_of_hyperplanes()-1);
			Point ann = P->get_site(nnIndex);
			double annEpsilon = pointsDist(p, ann, false)/pointsDist(p, nn, false)-1;
			falconnANNs.push_back(annEpsilon);
		}

		timer.start();
		bool boxContains = P->contains_point_lshbox(p, &nnIndex);
		boxTime.push_back(timer.elapsed_seconds());
		boxSucceeded.push_back(boxContains);
		boxDists.push_back(boxContains?0.0:minDist);
		if (boxContains) {
			boxANNs.push_back(0.0);
		}
		else {
			Point nn = P->get_site(P->num_of_hyperplanes()-1);
			Point ann = P->get_site(nnIndex);
			double annEpsilon = pointsDist(points[i], ann, false)/pointsDist(p, nn, false)-1;
			boxANNs.push_back(annEpsilon);
		}

		timer.start();
		bool exactContains = P->contains_point_exact_nn(p, 0, &nnIndex);
		exactTime.push_back(timer.elapsed_seconds());
		exactSucceeded.push_back(exactContains);
		exactDists.push_back(exactContains?0.0:minDist);
		if (exactContains) {
			exactANNs.push_back(0.0);
		}
		else {
			Point nn = P->get_site(P->num_of_hyperplanes()-1);
			Point ann = P->get_site(nnIndex);
			double annEpsilon = pointsDist(points[i], ann, false)/pointsDist(p, nn, false)-1;
			exactANNs.push_back(annEpsilon);
		}
	}

	std::vector<std::pair<Point, std::pair<int, double> > > externalPoints;
	create_random_external_points(P, points, externalPoints, 0.1);

	std::vector<double> isInTime2;
	std::vector<bool> isInSucceeded2;
	std::vector<double> exactTime2;
	std::vector<double> falconnTime2;
	std::vector<double> boxTime2;

	std::vector<double> exactDists2;
	std::vector<double> falconnDists2;
	std::vector<double> boxDists2;

	std::vector<bool> exactSucceeded2;
	std::vector<bool> falconnSucceeded2;
	std::vector<bool> boxSucceeded2;

	std::vector<double> exactANNs2;
	std::vector<double> falconnANNs2;
	std::vector<double> boxANNs2;

	for (int i=0; i<nqp; i++) {
		Point p = externalPoints[i].first;
		int minFacet = externalPoints[i].second.first;
		double minDist = externalPoints[i].second.second;

		timer.start();
		isInSucceeded2.push_back(P->is_in(p)!=-1);
		isInTime2.push_back(timer.elapsed_seconds());

		timer.start();
		bool lshContains = P->contains_point_lsh(p, numProbes, &nnIndex);
		falconnTime2.push_back(timer.elapsed_seconds());
		falconnSucceeded2.push_back(lshContains==false);
		falconnDists2.push_back(lshContains?minDist:0.0);
		if (!lshContains) {
			falconnANNs2.push_back(0.0);
		}
		else {
			Point ann = P->get_site(P->num_of_hyperplanes()-1);
			Point nn = P->get_site(minFacet);
			double annEpsilon = pointsDist(p, ann, false)/pointsDist(p, nn, false)-1;
			falconnANNs2.push_back(annEpsilon);
		}

		timer.start();
		bool boxContains = P->contains_point_lshbox(p, &nnIndex);
		boxTime2.push_back(timer.elapsed_seconds());
		boxSucceeded2.push_back(boxContains==false);
		boxDists2.push_back(boxContains?minDist:0.0);
		if (!boxContains) {
			boxANNs2.push_back(0.0);
		}
		else {
			Point ann = P->get_site(P->num_of_hyperplanes()-1);
			Point nn = P->get_site(minFacet);
			double annEpsilon = pointsDist(points[i], ann, false)/pointsDist(p, nn, false)-1;
			boxANNs2.push_back(annEpsilon);
		}

		timer.start();
		bool exactContains = P->contains_point_exact_nn(p, 0, &nnIndex);
		exactTime2.push_back(timer.elapsed_seconds());
		exactSucceeded2.push_back(exactContains==false);
		exactDists2.push_back(exactContains?minDist:0.0);
		if (!exactContains) {
			exactANNs2.push_back(0.0);
		}
		else {
			Point ann = P->get_site(P->num_of_hyperplanes()-1);
			Point nn = P->get_site(minFacet);
			double annEpsilon = pointsDist(points[i], ann, false)/pointsDist(p, nn, false)-1;
			exactANNs2.push_back(annEpsilon);
		}
	}

	delete internalPoint;
	return response;
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
	int algoType=USE_EXACT;
	int probes;
	int box_k;
	int box_l;
	float w;
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
		("boxk", po::value<int>(&box_k), "#buckets for lshbox")
		("boxl", po::value<int>(&box_l), "#hashtables for lshbox");

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
   	int rnum = std::pow(0.000001,-2) * 400 * n * std::log(n);
    vars var(rnum,d,40,1,0.0000001,0,0,0,0,rng,get_snd_rand,
                 urdist,urdist1,true,false,false,false,false,true,USE_LSHBOX,0.1);

	// create polytope with internal repr
	stdHPolytope<double>* P = randomPolytope<double>(n, d, diameter);
	
	randomTransformation<stdHPolytope<double> >(P);	

	delete P;
	std::cout << results.dump();

}
