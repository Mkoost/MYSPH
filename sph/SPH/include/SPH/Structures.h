#pragma once

#include<array>

#ifndef MYSPH_STRUCTURES
#define MYSPH_STRUCTURES

namespace MYSPH{
	using dvec3 = std::array<double, 3>;

	struct MYSmoothParticle
	{
		dvec3 r;
		dvec3 velocity;
		dvec3 acceleration;
		double h;
		double density;
		double pressure;
		double energy;
		double denergy;
	};



}

#endif
