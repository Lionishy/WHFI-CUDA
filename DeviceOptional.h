#pragma once
#ifndef DeviceOptional_H
#define DeviceOptioanl_H

namespace iki {	namespace math { namespace device {
	template <typename T>
	struct Optional {
		bool is_present = false;
		union {
			char stub;
			T data;
		} value;
	};
} /* device */ } /* math */ } /* iki */

#endif