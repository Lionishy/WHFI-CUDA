#pragma once
#ifndef SimpleTableIO_H
#define SimpleTableIO_H

#include "SimpleTable.h"

#include <iostream>

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

#endif /* SimpleTableIO_H */