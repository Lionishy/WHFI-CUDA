#pragma once
#ifndef SimpleTableIO_H
#define SimpleTableIO_H

#include "SimpleTable.h"

#include <iostream>
#include <limits>

template <typename T, size_t Dim>
std::ostream &operator<<(std::ostream &ascii_os, iki::Argument<T, Dim> const &argument) {
	for (size_t arg_idx = 0u; arg_idx != Dim; ++arg_idx)
		ascii_os << argument.components[arg_idx] << ' ';
	return ascii_os;
}

template <typename T, size_t Dim, size_t Scale>
std::ostream& operator<<(std::ostream &ascii_os, iki::UniformSimpleTable<T, Dim, Scale> const &table) {
	using namespace std;
	using namespace iki;
	
	Index<Dim> expanded_index = first_index<Dim>();
	for (auto it = table.data, end = table.data + collapsed_size<Dim>(&table.bounds); it != end;) {
		for (size_t arg_idx = 0u; arg_idx != Dim; ++arg_idx)
			ascii_os << uniform_argument<T,Dim>(&expanded_index,&table.space) << ' ';
		for (size_t scale_counter = 0u; scale_counter != Dim; ++scale_counter)
			ascii_os << *it++ << ' ';
		ascii_os << '\n';
		next_index<Dim>(&expanded_index, &table.bounds);
	}

	return ascii_os << flush;
}

template <typename T, size_t Dim>
std::istream &operator>>(std::istream &ascii_is, iki::UniformSpace<T, Dim> &space) {
	for (size_t axis_idx = 0u; axis_idx != Dim; ++axis_idx)
		ascii_is >> space.axes[axis_idx].begin >> space.axes[axis_idx].step;
	return ascii_is;
}

template <typename T, size_t Dim, size_t Scale>
std::istream &operator>>(std::istream &ascii_is, iki::UniformSimpleTable<T, Dim, Scale> &table) {
	for (size_t collapsed_idx = 0u, end = iki::collapsed_size(table.bounds); collapsed_idx != end; ++collapsed_idx) {
		for (size_t arg_idx = 0u; arg_idx != Dim; ++arg_idx)
			ascii_is.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
		for (size_t scl_idx = 0u; scl_idx != Scale; ++scl_idx)
			ascii_is >> table.data[scl_idx+collapsed_idx*Scale];
	}
	return ascii_is;
}

#endif /* SimpleTableIO_H */