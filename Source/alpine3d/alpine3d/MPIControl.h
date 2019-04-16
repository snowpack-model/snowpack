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
#ifndef MPIWRAPPER_H
#define MPIWRAPPER_H

#ifdef ENABLE_MPI
	#include <mpi.h>
#endif

#ifdef _OPENMP
	#include <omp.h>
#endif

#include <meteoio/MeteoIO.h>
#include <cstdio>

/**
 * @class MPIControl
 * @brief A singleton class that deals with all aspects of parallelization within Alpine3D (MPI, OpenMP, PETSc)
 * Any class that wishes to utilize the MPI/OpenMP/PETSc capabilities should include this header and make
 * sure that the respective cmake compilation options (called MPI, OPENMP, PETSC) are activated. The methods
 * are designed in a transparent way also returning meaningful values if the compilation options were not
 * activated, e. g. bool master() will return true if MPI was not activated and if it was activated then
 * it will return false on all processes that are not master and true solely for the master process.
 * @author Thomas Egger
 * @date   2014-07-15
 */
class MPIControl
{
	public:
		/**
		 *  Returns the single instance of this class. If an instance wasn't initialized yet, it is instantiated.
		 * @return a reference to the single instance
		 * @note MPI_Init and PetscInitialize are called entirely with NULL arguments
		 */
		static MPIControl& instance();

		// The following methods refer to OpenMP
		/**
		 * Returns whether or not OpenMP has been activated
		 * @return true if _OPENMP has been defined, false otherwise
		 */
		bool openmp() const;

		/**
		 * Returns the current OpenMP thread number or 0 if OpenMP has not been activated
		 * @return current thread number
		 */
		size_t thread() const;

		/**
		 * Returns the total number of OpenMP threads or 1 if OpenMP has not been activated
		 * @return number of threads
		 */
		size_t max_threads() const;

		// The following methods refer to MPI
		/**
		 * Returns whether the calling process is the master process or not
		 * @return true if caller is master, false otherwise
		 */
		bool master() const;

		/**
		 * Returns the rank of the master process
		 * @return rank of the master process
		 */
		size_t master_rank() const;

		/**
		 * Returns the rank of the calling process or 0 if MPI was not activated
		 * @return rank number of the calling process
		 */
		size_t rank() const;

		/**
		 * Returns the size of the world in usage, 1 if MPI was not activated
		 * @return the size of the MPI world
		 */
		size_t size() const;

		/**
		 * Returns the name of the processor of the calling process or "local" if MPI was not activated
		 * @return string name of the processor
		 */
		std::string name() const;

		/**
		 * This method allows to synchronize all MPI processes. It blocks the calling process until all
		 * processes have called barrier()
		 */
		void barrier() const;

		/**
		 * Gathers the integer send_values from all processes into a vector of integers on the root node.
		 * Index 5 of the receive_vector represents the send_value of process #5.
		 * @param[in] send_value The int value a process sends
		 * @param[out] receive_vector The buffer for the root process to hold all the send_values
		 * @param[in] root The process rank that will gather the send_values
		 */
		void gather(const int& send_value, std::vector<int>& receive_vector, const size_t& root = 0);

		//@{
		/**
		 * Combines values from all processes and distributes the result back to all processes.
		 * - allreduce_max distributes the maximum of all processes
		 * - allreduce_min distributes the minimum of all processes
		 * - allreduce_sum distributes the sum of all processes
		 * @param[in,out] value The value that is used to perform the reduction and to hold the result
		 */
		void allreduce_max(double& value);
		void allreduce_min(double& value);
		void allreduce_sum(double& value);
		void allreduce_sum(int& value);
		//@}

		/**
		 * This method is used when deserializing a class T from a void* representing a char*,
		 * instantiating an object from a string
		 * @param[in] in The void*, which is actually a char*
		 * @param[in] len The length of the string
		 * @param[out] obj The newly instantiated object
		 * @note class T must support the extraction and insertion operators, >> and <<
		 */
		template<class T> static void deserialize(const void* in, const size_t& len, T& obj)
		{
			std::stringstream obj_stream;
			obj_stream.write(reinterpret_cast<const char*>(in), len);
			obj_stream >> obj;
		}

		/**
		 * This method is used when serializing a class T, turning an object into a string
		 * @param[out] out The void**, which is actually a char**
		 * @param[in] obj The object to serialize
		 * @param[in] alloc allocate memory?
		 * @return The length of the string
		 * @note class T must support the extraction and insertion operators, >> and <<
		 */
		template<class T> static size_t serialize(void** out, const T& obj, const bool alloc=false)
		{
			std::stringstream obj_stream;
			obj_stream << obj;
			const size_t len = obj_stream.str().size();
			if (alloc) {
				*out = new char[len];
			}
			obj_stream.read(reinterpret_cast<char*>(*out), len);
			return len;
		}

		#ifdef ENABLE_MPI
		/**
		 * For custom objects of class T a custom sum function is required
		 * MPI_Op_create expects exactly this interface, thus it cannot be changed
		 * @param[in] in The void* representing a summand of the sum operation
		 * @param[out] out The void* representing the result of the sum operation
		 * @param[in] datatype A pointer to the custom datatype previously committed
		 * @note class T must support the extraction and insertion operators, \>\> and \<\<
		 *       and furthermore support the operator+
		 */
		template <class T> static void op_sum_func(void* in, void* out, int* /*len*/, MPI_Datatype* datatype)
		{
			int tmp_size;
			MPI_Type_size(*datatype, &tmp_size);

			T in_obj, out_obj;

			deserialize(in, static_cast<size_t>(tmp_size), in_obj);
			deserialize(out, static_cast<size_t>(tmp_size), out_obj);

			out_obj += in_obj;

			serialize(&out, out_obj);
		}
		#endif

		#ifdef ENABLE_MPI
		/** @brief Adds up the values of type T from all processes and distributes the sum back to all processes.
		 * @param obj The obj that is a summand of the global sum and will hold the result
		 * @note class T must support the extraction and insertion operators, \>\> and \<\<
		 *       and furthermore support the operator+
		 */
		template <class T> void allreduce_sum(T& obj)
		{
			if (size_ <= 1) return;

			MPI_Op op;
			MPI_Op_create(op_sum_func<T>, true, &op);

			void* in_obj_char=NULL;
			void* out_obj_char=NULL;

			const size_t len = serialize(&in_obj_char, obj, true);
			out_obj_char = new char[len];

			MPI_Datatype chars_datatype;
			MPI_Type_contiguous(static_cast<int>(len), MPI_CHAR, &chars_datatype);
			MPI_Type_commit(&chars_datatype);

			MPI_Allreduce(in_obj_char, out_obj_char, 1, chars_datatype, op, MPI_COMM_WORLD);
			MPI_Op_free(&op);

			deserialize(out_obj_char, len, obj);

			delete[] (char*)out_obj_char;
			delete[] (char*)in_obj_char;
		}
		#else
		template <class T> void allreduce_sum(T& /*obj*/) {}
		#endif

		#ifdef ENABLE_MPI
		/**
		 * @brief Broadcast an object via MPI or in case MPI is not activated don't do anything.
		 * @param obj that is broadcasted from process root to all processes
		 * @param[in] root The process rank that will commit the broadcast value, all others receive only
		 * @note Class T needs to have the serialize and deseralize operator << and >> implemented
		 */
		template <class T> void broadcast(T& obj, const size_t& root = 0)
		{
			if (size_ <= 1) return;

			std::string obj_string;

			if (rank_ == root) { // Build the string that represents the object
				std::stringstream obj_stream;
				obj_stream << obj;
				obj_string.insert(0, obj_stream.str());
			}

			broadcast(obj_string, root);

			if (rank_ != root) { // Deserialize the object broadcasted
				std::stringstream obj_stream;
				obj_stream << obj_string;
				obj_stream >> obj;
			}
		}
		#else
		template <class T> void broadcast(T& /*obj*/, const size_t& root=0) {(void)root;}
		#endif

		#ifdef ENABLE_MPI
		/**
		 * @brief Broadcast a vector of class T objects via MPI or in case MPI is not activated don't do anything.
		 * @param vec_obj A vector of class T objects that shall be broadcasted,
		 *                 or if not root the vector that will hold the pointers to the objects broadcasted
		 * @param[in] root The process rank that will commit the broadcast value, all others receive only
		 * @note Class T needs to have the serialize and deseralize operator \<\< and \>\> implemented
		 */
		template <class T> void broadcast(std::vector<T>& vec_obj, const size_t& root = 0)
		{
			if (size_ <= 1) return;

			std::string obj_string;
			size_t vec_size;

			if (rank_ == root) { // Build the string that represents the objects
				std::stringstream objs_stream;
				for (size_t ii=0; ii<vec_obj.size(); ii++)
					objs_stream << vec_obj[ii];
				obj_string.insert(0, objs_stream.str());
				vec_size = vec_obj.size();
			}

			broadcast(vec_size, root);
			broadcast(obj_string, root);

			if (rank_ != root) { // Deserialize the objects broadcasted
				vec_obj.clear();

				std::stringstream objs_stream;
				objs_stream << obj_string;

				T tmp;
				for (size_t ii=0; ii<vec_size; ii++) {
					objs_stream >> tmp;
					vec_obj.push_back(tmp);
				}
			}
		}
		#else
		template <class T> void broadcast(std::vector<T>& /*vec_obj*/, const size_t& root=0) {(void)root;}
		#endif

		#ifdef ENABLE_MPI
		/**
		 * @brief	Scatter the objects pointed to by vector<T*> by slices to all preocesses.
		 *        Internally no MPI_Scatterv or the like is used, because the size of the
		 *        strings may become extremely large. Thusly internally blocking send and receive calls
		 *        are utilized. In case MPI is not activated vec_local is not changed in any way.
		 * @param[in, out] vec_local A vector of T* pointers to objects that shall be scattered;
		 *                 if root then this vector will also hold the pointers to the scattered objects
		 * @param[in] root The process rank that will scatter the values, all others receive only
		 * @note Class T needs to have the serialize and deseralize operator << and >> implemented
		 */
		template <class T> void scatter(std::vector<T*>& vec_local, const size_t& root = 0)
		{
			if (size_ <= 1) return;

			if (rank_ == root) {
				for (size_t ii=1; ii<size_; ii++) { // HACK: Assuming root is master
					size_t startx, deltax;
					getArraySliceParams(vec_local.size(), ii, startx, deltax);

					send((int)deltax, ii);
					for (size_t jj=startx; jj<(startx+deltax); jj++) {
						std::string obj_string;
						std::stringstream objs_stream;

						objs_stream << *(vec_local[jj]);
						obj_string.insert(0, objs_stream.str());
						send(obj_string, ii);
						delete vec_local[jj];
					}
				}

				size_t startx, deltax;
				getArraySliceParams(vec_local.size(), startx, deltax);
				vec_local.resize(deltax);
			} else {
				std::string obj_string;
				int deltax;
				receive(deltax, root);
				vec_local.resize((size_t)deltax);

				T obj; // temporary object
				for (size_t ii=0; ii<vec_local.size(); ii++) {
					std::stringstream objs_stream;
					receive(obj_string, root);

					objs_stream << obj_string;

					objs_stream >> obj;
					vec_local[ii] = new T(obj);
				}
			}
		}
		#else
		template <class T> void scatter(std::vector<T*>& /*vec_local*/, const size_t& root=0) {(void)root;}
		#endif

		#ifdef ENABLE_MPI
		/**
		 * @brief	Send the objects pointed to by vector<T*> to process \#destination
		 * @param[in] vec_local A vector of T* pointers to objects that shall be sent
		 * @param[in] destination The process rank that will receive the values
		 * @param[in] tag Arbitrary non-negative integer assigned to uniquely identify a message
		 * @note Class T needs to have the serialize and deseralize operator << and >> implemented
		 */
		template <class T> void send(const std::vector<T*>& vec_local, const size_t& destination, const int& tag=0);
		template <class T> void send(const std::vector<T>& vec_local, const size_t& destination, const int& tag=0);
		#else
		template <class T> void send(const std::vector<T*>& /*vec_local*/, const size_t& /*destination*/, const int& tag=0) {(void)tag;}
		template <class T> void send(const std::vector<T>& /*vec_local*/, const size_t& /*destination*/, const int& tag=0) {(void)tag;}
		#endif

		#ifdef ENABLE_MPI
		/**
		 * @brief	Receive vector of objects from process \#source
		 * @param[in] vec_local A vector of T* pointers to receive the object pointers
		 * @param[in] source The process rank that will send the objects
		 * @param[in] tag Arbitrary non-negative integer assigned to uniquely identify a message
		 * @note Class T needs to have the serialize and deseralize operator << and >> implemented
		 */
		template <class T> void receive(std::vector<T*>& vec_local, const size_t& source, const int& tag=0);
		template <class T> void receive(std::vector<T>& vec_local, const size_t& source, const int& tag=0);
		#else
		template <class T> void receive(std::vector<T*>& /*vec_local*/, const size_t& /*source*/, const int& tag=0) {(void)tag;}
		template <class T> void receive(std::vector<T>& /*vec_local*/, const size_t& /*source*/, const int& tag=0) {(void)tag;}
		#endif

		#ifdef ENABLE_MPI
		/**
		 * @brief	Gathers the objects pointed to by vector<T*> from all processes into a vector<T*> on
		 *        the root node. Internally no MPI_Gatherv or the like is used, because the size of the
		 *        strings may become extremely large. Thusly internally blocking send and receive calls
		 *        are utilized. In case MPI is not activated vec_local is not changed in any way.
		 * @param[in, out] vec_local A vector of T* pointers to objects that shall be transmitted to root;
		 *                 if root then this vector will also hold the pointers to the gathered objects
		 * @param[in] root The process rank that will commit the broadcast value, all others receive only
		 * @note Class T needs to have the serialize and deseralize operator \<\< and \>\> implemented
		 */
		template <class T> void gather(std::vector<T*>& vec_local, const size_t& root = 0)
		{
			if (size_ <= 1) return;

			std::vector<int> vec_sizes;
			const size_t sum = vec_local.size();

			gather(vec_local.size(), vec_sizes); //global offsets
			allreduce_sum(sum); //global amount of T*

			if (rank_ == root) {
				std::string obj_string;
				size_t local_sum = vec_local.size();

				vec_local.resize(sum);

				T obj; //temporary object
				for (size_t ii=1; ii<vec_sizes.size(); ii++) {//HACK: Assuming master is always rank 0
					for (size_t jj=0; jj<vec_sizes[ii]; jj++) {
						std::stringstream objs_stream;
						receive(obj_string, ii);

						objs_stream << obj_string;
						objs_stream >> obj;
						vec_local[local_sum+jj] = new T(obj);
					}
					local_sum += vec_sizes[ii];
				}

			} else {
				for (size_t ii=0; ii<vec_local.size(); ii++) {
					std::string obj_string;
					std::stringstream objs_stream;

					objs_stream << *(vec_local[ii]);
					obj_string.insert(0, objs_stream.str());
					send(obj_string, root);
				}
			}
		}
		#else
		template <class T> void gather(std::vector<T*>& /*vec_local*/, const size_t& root=0) {(void)root;}
		#endif

		/** @brief Helper to split an array for domain decomposition.
		 * Given the overall size of an array calculate a start index and block length
		 * for the subarray on the basis of the current MPI rank
		 * @param[in] dimx The overall size of the array to split up
		 * @param[out] startx_sub The start of the subarray for the given processor rank
		 * @param[out] nx_sub The block length of the subarray for the given processor rank
		 */
		void getArraySliceParams(const size_t& dimx, size_t& startx_sub, size_t& nx_sub) const {getArraySliceParams(dimx, size_, rank_, startx_sub, nx_sub);}
		void getArraySliceParams(const size_t& dimx, const size_t& idx_wk, size_t& startx_sub, size_t& nx_sub) const {getArraySliceParams(dimx, size_, idx_wk, startx_sub, nx_sub);}

		void getArraySliceParamsOptim(const size_t& dimx, size_t& startx_sub, size_t& nx_sub,const mio::DEMObject& dem, const mio::Grid2DObject& landuse)
		     {getArraySliceParamsOptim(dimx, rank_, startx_sub, nx_sub, dem, landuse);}

		void getArraySliceParamsOptim(const size_t& dimx, const size_t& idx_wk, size_t& startx_sub, size_t& nx_sub,const mio::DEMObject& dem, const mio::Grid2DObject& landuse);

		/**
		* @brief Returns the parameters for splitting an array in several, balanced sub-arrays.
		* This is mostly usefull for parallel calculations, where an array will be split and sent to different
		* workers.
		* @param[in] dimx number of cells in the desired dimension
		* @param[in] nbworkers total number of slices
		* @param[in] idx_wk current slice index (starting at 0)
		* @param[out] startx_sub calculated start index for the current slice
		* @param[out] nx_sub calculated number of cells (in the desired dimension) of the current slice
		*/
		static void getArraySliceParams(const size_t& dimx, const size_t& nbworkers, const size_t& idx_wk, size_t& startx_sub, size_t& nx_sub);

		/**
		* @brief Returns the parameters for splitting an array in several, balanced sub-arrays.
		* This is mostly usefull for parallel calculations, where an array will be split and sent to different
		* workers.
		* @param[in] dimx number of cells in the desired dimension
		* @param[in] nbworkers total number of slices
		* @param[out] offset calculated start index for all slices
		* @param[out] nx calculated number of cells (in the desired dimension) of all slice
		*/
		static void getArraySliceParams(const size_t& dimx, const size_t& nbworkers, std::vector<size_t>& offset, std::vector<size_t>& nx);

	private:
		MPIControl();
		~MPIControl();

		MPIControl(MPIControl const&);     //No implementation allowed - singleton class
		void operator=(MPIControl const&); //No implementation allowed - singleton class
		void broadcast(std::string& message, const size_t& root);

		void send(std::string& message, const size_t& recipient, const int& tag=0);
		void receive(std::string& message, const size_t& source, const int& tag=0);
		void send(const int& value, const size_t& recipient, const int& tag=0);
		void receive(int& value, const size_t& source, const int& tag=0);

		static void checkSuccess(const int& ierr);

		size_t rank_;          // the rank of this process
		size_t size_;          // the number of all processes
		std::string name_;  // the name of the node this process is running on
};

#endif
