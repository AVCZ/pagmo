/*****************************************************************************
 *   Copyright (C) 2004-2009 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 3 of the License, or       *
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

#include <climits>
#include <iostream>
#include <vector>
#include <list>

#include "src/algorithms.h"
#include "src/archipelago.h"
#include "src/island.h"
#include "src/problems.h"
#include "src/topologies.h"
#include "src/archipelago.h"
#include "src/problem/base.h"
#include"src/keplerian_toolbox/keplerian_toolbox.h"

using namespace pagmo;

int main()
{
	boost::array<double,6> elem = {1.45815287 * ASTRO_AU,0.222828423,10.8289894 * ASTRO_DEG2RAD,178.7579447 * ASTRO_DEG2RAD , 304.3704767 * ASTRO_DEG2RAD,55.6339120 * ASTRO_DEG2RAD};
	::kep_toolbox::epoch epoch(55400,::kep_toolbox::epoch::MJD);
	::kep_toolbox::planet target(epoch,elem,ASTRO_MU_SUN,200.0,100.0);
	pagmo::problem::sample_return prob(target);

	pagmo::algorithm::sa_corana algo1(10000,1,0.01);
	pagmo::algorithm::de algo2(500,0.8,0.8,3);
	pagmo::algorithm::nlopt_sbplx algo3(500,1e-4);

	pagmo::archipelago a = pagmo::archipelago(pagmo::topology::rim());
	a.push_back(pagmo::island(prob,algo3,20));
	a.push_back(pagmo::island(prob,algo1,20));
	a.push_back(pagmo::island(prob,algo2,20));
	a.push_back(pagmo::island(prob,algo1,20));
	a.push_back(pagmo::island(prob,algo2,20));
	std::cout << prob << std::endl;

	for (int i=0;i<50;++i)
	{
		std::cout << a.get_island(0).get_population().champion().f << std::endl;
		a.evolve();
	}
	std::cout << prob.pretty(a.get_island(0).get_population().champion().x) <<std::endl;
	std::cout << a.get_island(0).get_population().champion().x <<std::endl;
	std::cout << target << std::endl;
	return 0;
}
