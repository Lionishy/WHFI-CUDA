#pragma once
#ifndef DispersionRelation_CH
#define DispersionRelation_CH

#include "PhysicalParameters.h"

namespace iki {	namespace whfi {
	template <typename T, typename F_t>
	struct DispersionRelation {
		T operator()(T omega, T k) {
			return T(1. / 1836.) + k * k
				- p.nc * (omega / (k * p.betta_root_c) - p.bulk_to_term_c) * Z((omega - T(1.)) / (k * p.betta_root_c) - p.bulk_to_term_c)
				- p.nh * (omega / (k * p.betta_root_h) - p.bulk_to_term_h) * Z((omega - T(1.)) / (k * p.betta_root_h) - p.bulk_to_term_h);
		}

		DispersionRelation(F_t Z, PhysicalParamenters<T> p) : Z(Z), p(p) { }

	private:
		F_t Z;
		PhysicalParamenters<T> p;
	};

	template <typename T, typename F_t>
	struct DispersionRelationDerivative {
		T operator()(T omega, T k) {
			T arg_c = (omega - T(1.)) / (k * p.betta_root_c) - p.bulk_to_term_c;
			T arg_h = (omega - T(1.)) / (k * p.betta_root_h) - p.bulk_to_term_h;
			T Zc = Z(arg_c), Zh = Z(arg_h);
			return p.nc / (k * p.betta_root_c) * (-Zc + (omega / (k * p.betta_root_c) - p.bulk_to_term_c) * (Zc * arg_c + T(1.)))
				+ p.nh / (k * p.betta_root_h) * (-Zh + (omega / (k * p.betta_root_h) - p.bulk_to_term_h) * (Zh * arg_h + T(1.)));
		}

		DispersionRelationDerivative(F_t Z, PhysicalParamenters<T> p) : Z(Z), p(p) { }

	private:
		F_t Z;
		PhysicalParamenters<T> p;
	};
} /* whfi */ } /* iki */

#endif