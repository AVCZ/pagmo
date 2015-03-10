/*****************************************************************************
 *   Copyright (C) 2004-2013 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

// Test code for the transformation of the optima of meta-problems from the original problems.

#include <iostream>
#include <vector>
#include <cmath>
#include "../src/pagmo.h"

using namespace pagmo;

const double EPS = 10e-9;

bool is_eq(const decision_vector & x1, const decision_vector & x2){
	if (x1.size() != x2.size()) return false;
	for (unsigned int i = 0; i < x1.size(); ++i){
		if (fabs(x1[i] - x2[i]) > EPS) return false;
	}
	return true;
}

int test_optima_transformations()
{
	// Going to test the code on the Himmelblau function,
	// see http://en.wikipedia.org/wiki/Himmelblau%27s_function
	// The four minima of the function in the range [-6,6]x[-6,6] are:
	// (x0,y0) = (3.0, 2.0)
	// (x1,y1) = (-2.805118,3.131312)
	// (x2,y2) = (-3.779310,-3.283186)
	// (x3,y3) = (3.584428,-1.848126)

	problem::himmelblau prob;
	// Himmelblau function minima in the bounding box [-6,6]:
	const std::vector<decision_vector> best_x = prob.get_best_x();

	// ------------------------- Shifting --------------------------
	std::vector<double> shift_vec1 = { 1., 1. };
	problem::shifted shifted_prob1(prob, shift_vec1);
	// Previous function minima shifted by vector (1,1):
	const std::vector<decision_vector> shif1_best_x = shifted_prob1.get_best_x();

	std::vector<double> shift_vec2 = { -2., 0. };
	problem::shifted shifted_prob2(shifted_prob1, shift_vec2);
	// Previous function minima shifted by vector (-2,0):
	const std::vector<decision_vector> shif2_best_x = shifted_prob2.get_best_x();

	std::vector<double> shift_vec3 = { 1., -1. };
	problem::shifted shifted_prob3(shifted_prob2, shift_vec3);
	// "Previous function minima shifted by vector (1,-1) (which is the original problem):
	const std::vector<decision_vector> shif3_best_x = shifted_prob3.get_best_x();

	for (unsigned int i = 0; i < best_x.size(); ++i) {
		if (!is_eq(best_x[i],shif3_best_x[i])) return 1;
	}
	// -------------------------------------------------------------

	// ------------------------- Normalization ---------------------
	problem::normalized normalized_prob(prob);
	// Normalizing bounds (to [-1,1]x[-1,1]) of the Himmelblau function and the corresponding transformed minima are:
	const std::vector<decision_vector> norm_best_x = normalized_prob.get_best_x();

	const std::vector<decision_vector> expected_normalized_minima = {
		{best_x[0][0] / 6, best_x[0][1] / 6},
		{best_x[1][0] / 6, best_x[1][1] / 6},
		{best_x[2][0] / 6, best_x[2][1] / 6},
		{best_x[3][0] / 6, best_x[3][1] / 6},
	};

	for (unsigned int i = 0; i < norm_best_x.size(); ++i) {
		if (!is_eq(expected_normalized_minima[i], norm_best_x[i])) return 1;
	}
	// -------------------------------------------------------------

	// ------------------------- Rotation --------------------------
	std::vector<std::vector<double>> rot_id_mat = {{1.,0.},{0.,1.}}; // identity matrix
	problem::rotated rotated_id_prob(normalized_prob,rot_id_mat);
	// Normalized Himmelblau rotated by the identity matrix and the corresponding minima are:
	const std::vector<decision_vector> rot_id_best_x = rotated_id_prob.get_best_x();

	for (unsigned int i = 0; i < norm_best_x.size(); ++i) {
		if (!is_eq(expected_normalized_minima[i], rot_id_best_x[i])) return 1;
	}

	const double pi = std::acos(-1);
	const double phi = pi / 2;
	const std::vector<std::vector<double>> rot_mat90 = {{std::cos(phi), -std::sin(phi)},{std::sin(phi), std::cos(phi)}};

	// Take the Himmelblau function and rotate by 90 degrees (plus normalize)
	problem::rotated rotated_prob1(prob,rot_mat90);

	// Rotate again by 90 degrees (plus normalize)
	problem::rotated rotated_prob2(rotated_prob1,rot_mat90);

	// Rotate again by 90 degrees (plus normalize)
	problem::rotated rotated_prob3(rotated_prob2,rot_mat90);

	// Rotate again by 90 degrees (plus normalize). We get back to the original problem (normalized_prob),
	// but there is a twist here, because along the way the bounds were normalized several
	// times (exactly four time). Therefore, all the vectors are scaled by the coef (1 / sqrt(2)) ** 4.
	// See src/problem/rotated.cpp for more details (function configure_new_bounds() in particular).
	problem::rotated rotated_prob4(rotated_prob3,rot_mat90);
	// Previous function rotated by 90 degrees counterclockwise and the minima are (after another normalization)
	// basically the same as in the original problem, but normalized by the factor of (1 / sqrt(2)) ** 3:
	const std::vector<decision_vector> rot4_best_x = rotated_prob4.get_best_x();

	const double sqrt2 = std::sqrt(2);
	// First normalization goes from the bounding box of [-6,6]x[-6,6] to [-1,1]x[-1,1]
	// and every other normalization goes from [-sqrt(2),sqrt(2)]x[-sqrt(2),sqrt(2)] to [-1,1]x[-1,1].
	const double coef = (1 / 6.) * (1 / sqrt2) * (1 / sqrt2) * (1 / sqrt2);

	const std::vector<decision_vector> expected_shrunk_minima = {
		{best_x[0][0] * coef, best_x[0][1] * coef},
		{best_x[1][0] * coef, best_x[1][1] * coef},
		{best_x[2][0] * coef, best_x[2][1] * coef},
		{best_x[3][0] * coef, best_x[3][1] * coef},
	};

	for (unsigned int i = 0; i < expected_shrunk_minima.size(); ++i) {
		if (!is_eq(expected_shrunk_minima[i], rot4_best_x[i])) return 1;
	}
	// -------------------------------------------------------------

	return 0;
}

int main()
{
	return test_optima_transformations();
}
