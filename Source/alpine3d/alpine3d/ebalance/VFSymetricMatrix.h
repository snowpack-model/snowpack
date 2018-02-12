/***********************************************************************************/
/*  Copyright 2009-2015 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of Alpine3D.
    Alpine3D is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Alpine3D is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with Alpine3D.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef VFSYMETRICMATRIX_H
#define VFSYMETRICMATRIX_H

#if defined(__clang__)
	#include <unordered_map>
#elif defined(__GNUC__)
	#if __GNUC__ < 4
		#include <ext/hash_map>
	#else
		#include <tr1/unordered_map>
	#endif
#else
	#include <unordered_map>
#endif
#include <map>
#include <meteoio/MeteoIO.h>

/*
This data structure store the symetric matrix of view factor.
T represents the type of stored value
U represents the type of passed and returned value of the structure
*/
template<class T, class U> class VFSymetricMatrix
{
	public:
		VFSymetricMatrix();
		VFSymetricMatrix(const int anx, const int any);
		~VFSymetricMatrix();

		void resize(const int nx, const int ny);

		void GetSize(int &nx, int &ny);

		void Destroy();

		U getElement(const unsigned int x, const unsigned int y);

		const U operator()(const unsigned int& x, const unsigned int& y) const;

		void setElement(const unsigned int x, const unsigned int y, U val);

		int size();

		VFSymetricMatrix<T, U>& operator=(VFSymetricMatrix<T, U>& val);

	private:
#if defined(__clang__)
		typedef std::unordered_map< int, T > my_map;
#elif defined(__GNUC__)
	#if __GNUC__ < 4
		typedef __gnu_cxx::hash_map< int, T > my_map;
	#else
		typedef std::tr1::unordered_map< int, T > my_map;
	#endif
#else
		typedef std::tr1::unordered_map< int, T > my_map;
#endif
		//typedef map< int, T > my_map;

		my_map mapData;	//for ordered map, declare <, = and index operators in constructor
		unsigned int nx;
		unsigned int ny;

};

template<class T, class U> int VFSymetricMatrix<T, U>::size()
{
	return (int) mapData.size();
}

template<class T, class U> U VFSymetricMatrix<T, U>::getElement(const unsigned int x, const unsigned int y)
{
	#ifndef NOSAFECHECKS
	if ((x >= nx) || (y >= ny))
		throw mio::IndexOutOfBoundsException(std::string(), AT);
	#endif
	int idx = x*ny+y;
	if (x > y) {
		idx = y*ny+x;
	}
	typename my_map::iterator j = mapData.find(idx);
	if ( j == mapData.end() ) {
		return static_cast<U>(0.);
	} else {
		return static_cast<U>((*j).second);
	}
}

template<class T, class U> const U VFSymetricMatrix<T, U>::operator()(const unsigned int& x, const unsigned int& y) const
{
	#ifndef NOSAFECHECKS
	if ((x >= nx) || (y >= ny))
		throw mio::IndexOutOfBoundsException(std::string(), AT);
	#endif
	int idx = x*ny+y;
	if (x > y) {
		idx = y*ny+x;
	}
	typename my_map::const_iterator j = mapData.find(idx);
	if ( j == mapData.end() ) {
		return static_cast<U>(0.);
	} else {
		return static_cast<U>((*j).second);
	}
}

template<class T, class U> void VFSymetricMatrix<T, U>::setElement(unsigned int x, unsigned int y, U val)
{
	#ifndef NOSAFECHECKS
	if ((x >= nx) || (y >= ny))
		throw mio::IndexOutOfBoundsException(std::string(), AT);
	#endif
	T tval = static_cast<T>(val);
	if (tval != 0.) {
		int idx = x*ny+y;
		if (x > y) {
			idx = y*ny+x;
		}
		mapData[idx] = tval;
	}
}

template<class T, class U> VFSymetricMatrix<T, U>::VFSymetricMatrix()
{
	nx = ny = 0;
}

template<class T, class U> VFSymetricMatrix<T, U>::VFSymetricMatrix(const int anx, const int any)
{
	nx = ny = 0;
	resize(anx,any);
}

template<class T, class U> VFSymetricMatrix<T, U>::~VFSymetricMatrix()
{
	Destroy();
}

template<class T, class U> void VFSymetricMatrix<T, U>::resize(const int anx, const int any)
{
	Destroy();
	if ((anx > 0) && (any > 0)){
		nx = anx;
		ny = any;
	} else {
		throw mio::IndexOutOfBoundsException(std::string(), AT);
	}
}

template<class T, class U> void VFSymetricMatrix<T, U>::GetSize(int &anx, int &any)
{
	anx=nx;
	any=ny;
}

template<class T, class U> void VFSymetricMatrix<T, U>::Destroy()
{
	mapData.clear();
	nx=ny=0;
}

template<class T, class U> VFSymetricMatrix<T, U>& VFSymetricMatrix<T, U>::operator=(VFSymetricMatrix<T, U>& val)
{
	int anx,any;
	val.GetSize(anx,any);

	mapData = val.mapData;
	nx = anx;
	ny = any;

	return *this;
}

#endif
