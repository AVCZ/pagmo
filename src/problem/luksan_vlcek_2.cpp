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

#include "base.h"
#include "luksan_vlcek_2.h"

namespace pagmo { namespace problem {

/// Constructor.
/**
 * Construct the problem from its dimension.
 * Setting cub=clb=0 creates an instance of the original Luksan Vlcek equality constrained problem, example  5.2
 * Using clb<cub allows the obtain a problem formulation with inequality constraints
 *
 * @param[in] N Problem dimension
 * @param[in] clb lower bounds for the constraints.
 * @param[in] cub upper bounds for the constraints.
 * @throws value_error if N is smaller than 13 and is odd, cub < clb
 *
 * @see L.Luksan and J.Vlcek, "Sparse and Parially Separable Test Problems for Unconstrained and Equality Constrained Optimization"
 */
luksan_vlcek_2::luksan_vlcek_2(size_t N, const double clb, const double cub):base(N+2,0,1,2*(N-7),2*(N-7))
{
	if (N<=13 || 2*(N/2)!=N)
	{
		pagmo_throw(value_error,"Problem dimension needs to be at least 3");
	}
	if (clb >cub)
	{
		pagmo_throw(value_error,"N needs to be at least 14 and even.");
	}
	set_lb(-5);
	set_ub(5);
	m_clb = std::vector<double>(N-7,clb);
	m_cub = std::vector<double>(N-7,cub);
}

/// Clone method.
base_ptr luksan_vlcek_2::clone() const
{
	return base_ptr(new luksan_vlcek_2(*this));
}

/// Implementation of the objective function.
void luksan_vlcek_2::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	f[0] = 0.;
	for (size_t i=0; i<(x.size()-2)/2; i++)
	{
		double a1 = x[2*i]*x[2*i] - x[2*i+1];
		double a2 = x[2*i] - 1.;
		double a3 = x[2*i+2]*x[2*i+2] - x[2*i+3];
		double a4 = x[2*i+2] - 1.;
		double a5 = x[2*i+1] + x[2*i+3] - 2.;
		double a6 = x[2*i+1] - x[2*i+3];
		f[0] += 100.*a1*a1 + a2*a2 + 90.*a3*a3 + a4*a4 + 10.*a5*a5 + .1*a6*a6;
	}
}

/// Implementation of the constraint function.
void luksan_vlcek_2::compute_constraints_impl(constraint_vector &c, const decision_vector &x) const
{
	for (size_t i=0; i<(x.size()-2)-7; i++)
	{
		c[2*i] = (2.+5.*x[i+5]*x[i+5])*x[i+5] + 1.;
		for (size_t k=std::max((size_t)0,i-5); k<=i+1; k++)
		{
			c[2*i] += x[k]*(x[k]+1.);
		}
		c[2*i] = c[2*i] - m_cub[i];
		c[2*i+1] = (2.+5.*x[i+5]*x[i+5])*x[i+5] + 1.;
		for (size_t k=std::max((size_t)0,i-5); k<=i+1; k++)
		{
			c[2*i+1] += x[k]*(x[k]+1.);
		}
		c[2*i+1] = m_clb[i] - c[2*i+1];
	}
}

/// Implementation of the sparsity structure: automated detection
void luksan_vlcek_2::set_sparsity(int& lenG, std::vector<int>& iGfun, std::vector<int>& jGvar) const
{
	//Initial point
	decision_vector x0(get_dimension(),1);
	//Numerical procedure
	estimate_sparsity(x0, lenG, iGfun, jGvar);
}

} }