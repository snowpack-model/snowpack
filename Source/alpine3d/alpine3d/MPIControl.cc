/***********************************************************************************/
/*	Copyright 2009-2015 WSL Institute for Snow and Avalanche Research   SLF-DAVOS  */
/***********************************************************************************/
/* This file is part of Alpine3D.
    Alpine3D is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Alpine3D is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with Alpine3D. If not, see <http://www.gnu.org/licenses/>.
*/
#include <alpine3d/MPIControl.h>
#include <alpine3d/SnowpackInterfaceWorker.h>

#if (defined _WIN32 || defined __MINGW32__) && ! defined __CYGWIN__
	#include <winsock.h>
#else
	#include <unistd.h>
#endif
#ifdef ENABLE_PETSC
	#include <petscksp.h>
#endif

using namespace std;
using namespace mio;

// Local oprators << and >> on pairs, nedded to serialize pairs in MPI communication
std::istream& operator>>(std::istream& is, std::pair<size_t,size_t>& data)
{
	is.read(reinterpret_cast<char*>(&data.first), sizeof(data.first));
	is.read(reinterpret_cast<char*>(&data.second), sizeof(data.second));
	return is;
}

std::ostream& operator<<(std::ostream& os, const std::pair<size_t,size_t>& data)
{
	os.write(reinterpret_cast<const char*>(&data.first), sizeof(data.first));
	os.write(reinterpret_cast<const char*>(&data.second), sizeof(data.second));
	return os;
}

bool MPIControl::openmp() const
{
	#ifdef _OPENMP
		return true;
	#else
		return false;
	#endif
}

size_t MPIControl::thread() const
{
	#ifdef _OPENMP
		return static_cast<size_t>( omp_get_thread_num() );
	#else
		return 0;
	#endif
}

size_t MPIControl::max_threads() const
{
	#ifdef _OPENMP
		return static_cast<size_t>( omp_get_max_threads() );
	#else
		return 1;
	#endif
}

size_t MPIControl::rank() const
{
	return rank_;
}

size_t MPIControl::master_rank() const
{
	return 0; //HACK: this assumes the master is always 0
}

size_t MPIControl::size() const
{
	return size_;
}

std::string MPIControl::name() const
{
	return name_;
}

bool MPIControl::master() const
{
	return (rank_ == 0); //HACK: this assumes the master is always 0
}

#ifdef ENABLE_MPI
MPIControl::MPIControl()
{
	MPI_Init(NULL, NULL);

	#ifdef ENABLE_PETSC
	PetscInitialize(NULL, NULL, NULL, NULL);
	#endif

	int o_rank, o_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &o_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &o_size);

	if (o_rank<0 || o_size<0)
		throw mio::IOException("Invalid rank or size returned by MPI", AT);

	rank_ = static_cast<size_t>(o_rank);
	size_ = static_cast<size_t>(o_size);

	//MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN); // return info about errors

	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int name_len;
	MPI_Get_processor_name(processor_name, &name_len);
	name_ = std::string(processor_name);

	std::cout << "[i] Init of MPI on '" << name_ << "' with a total world size of " << size_ << " (I'm rank #" << rank_ << ")\n";
}

MPIControl::~MPIControl()
{
	MPI_Finalize();
}

void MPIControl::checkSuccess(const int& ierr)
{
	if (ierr != MPI_SUCCESS) {
		int err_len = 0;
		char err_buffer[MPI_MAX_ERROR_STRING];

		MPI_Error_string(ierr, err_buffer, &err_len);
		throw mio::IOException(string(err_buffer), AT);
	}
}

void MPIControl::broadcast(std::string& message, const size_t& root)
{
	int ierr;
	unsigned int msg_len = (unsigned int) message.size();

	//Now broadcast the size of the object and then the object itself
	ierr = MPI_Bcast(&msg_len, 1, MPI_UNSIGNED, root, MPI_COMM_WORLD);
	checkSuccess(ierr);

	if (rank_ != root) message.resize(msg_len);
	ierr = MPI_Bcast(const_cast<char*>(message.c_str()), msg_len, MPI_CHAR, root, MPI_COMM_WORLD);
	checkSuccess(ierr);
}

void MPIControl::send(std::string& message, const size_t& recipient, const int& tag)
{
	int ierr;
	unsigned int msg_len = (unsigned int) message.size();

	ierr = MPI_Send(&msg_len, 1, MPI_UNSIGNED, static_cast<int>(recipient), tag, MPI_COMM_WORLD);
	checkSuccess(ierr);

	ierr = MPI_Send(const_cast<char*>(message.c_str()), msg_len, MPI_CHAR, static_cast<int>(recipient), tag, MPI_COMM_WORLD);
	checkSuccess(ierr);
}

void MPIControl::receive(std::string& message, const size_t& source, const int& tag)
{
	int ierr;
	unsigned int msg_len;

	ierr = MPI_Recv(&msg_len, 1, MPI_UNSIGNED, static_cast<int>(source), tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	checkSuccess(ierr);

	message.resize(msg_len);
	ierr = MPI_Recv(const_cast<char*>(message.c_str()), msg_len, MPI_CHAR, static_cast<int>(source), tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	checkSuccess(ierr);
}

void MPIControl::send(const int& value, const size_t& recipient, const int& tag)
{
	int buffer = value;
	MPI_Send(&buffer, 1, MPI_INT, static_cast<int>(recipient), tag, MPI_COMM_WORLD);
}

void MPIControl::receive(int& value, const size_t& source, const int& tag)
{
	MPI_Recv(&value, 1, MPI_INT, static_cast<int>(source), tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void MPIControl::gather(const int& send_value, std::vector<int>& receive_vector, const size_t& root)
{
	int value = send_value; // make a copy to preserve const notion
	int* recv_data(NULL);
	if (rank_ == root) {
		receive_vector.resize(size_);
		recv_data = &receive_vector[0];
	}

	MPI_Gather(&value, 1, MPI_INT, recv_data, 1, MPI_INT, static_cast<int>(root), MPI_COMM_WORLD);
}

void MPIControl::allreduce_max(double& value)
{
	double buffer;
	MPI_Allreduce(&value, &buffer, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	value = buffer;
}

void MPIControl::allreduce_min(double& value)
{
	double buffer;
	MPI_Allreduce(&value, &buffer, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	value = buffer;
}

void MPIControl::allreduce_sum(double& value)
{
	double buffer;
	MPI_Allreduce(&value, &buffer, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	value = buffer;
}

void MPIControl::allreduce_sum(int& value)
{
	int buffer;
	MPI_Allreduce(&value, &buffer, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	value = buffer;
}

void MPIControl::barrier() const
{
	MPI_Barrier(MPI_COMM_WORLD);
}

#else
std::string getHostName() {
	static const size_t len = 4096;

	#if (defined _WIN32 || defined __MINGW32__) && ! defined __CYGWIN__
		TCHAR infoBuf[len];
		DWORD bufCharCount = len;
		if ( !GetComputerName( infoBuf, &bufCharCount ) )
			return std::string("N/A");

		return std::string(infoBuf);
	  #else
		char name[len];
		if (gethostname(name, len) != 0) {
			return std::string("localhost");
		}

		if (name[0] == '\0') return std::string("localhost");
		return std::string(name);
	#endif
}

MPIControl::MPIControl() : rank_(0), size_(1), name_( getHostName() )
{
	#ifdef _OPENMP
		std::cout << "[i] Init of OpenMP on '" << name_ << "' with a pool of " << omp_get_max_threads() << " threads\n";
	#endif
}

MPIControl::~MPIControl() {}
void MPIControl::barrier() const {}
void MPIControl::allreduce_max(double&) {}
void MPIControl::allreduce_min(double&) {}
void MPIControl::allreduce_sum(double&) {}
void MPIControl::allreduce_sum(int&) {}
void MPIControl::gather(const int& val, std::vector<int>& vec, const size_t&) { vec.resize(1, val); }
#endif


#ifdef ENABLE_MPI
/**
 * @brief	Send the objects pointed to by vector<T*> to process \#destination
 * @param[in] vec_local A vector of T* pointers to objects that shall be sent
 * @param[in] destination The process rank that will receive the values
 * @param[in] tag Arbitrary non-negative integer assigned to uniquely identify a message
 * @note Class T needs to have the serialize and deseralize operator << and >> implemented
 */
template <class T> void MPIControl::send(const std::vector<T*>& vec_local, const size_t& destination, const int& tag)
{
	if ((size_ <= 1) || (rank_ == destination)) return;

	const size_t v_size = vec_local.size();

	send((int)v_size, destination, tag); // first send the size of the vector
	for (size_t ii=0; ii<v_size; ii++) {
		std::string obj_string;
		std::stringstream objs_stream;

		objs_stream << *(vec_local[ii]);
		obj_string.insert(0, objs_stream.str());
		send(obj_string, destination, tag);
	}
}
/**
 * @brief	Receive vector of objects from process \#source
 * @param[in] vec_local A vector of T* pointers to receive the object pointers
 * @param[in] source The process rank that will send the objects
 * @param[in] tag Arbitrary non-negative integer assigned to uniquely identify a message
 * @note Class T needs to have the serialize and deseralize operator << and >> implemented
 */
template <class T> void MPIControl::receive(std::vector<T*>& vec_local, const size_t& source, const int& tag)
{
	if ((size_ <= 1) || (rank_ == source)) return;

	if (!vec_local.empty())
		throw mio::IOException("The vector to receive pointers has to be empty (please properly free the vector)", AT);

	std::string obj_string;
	int v_size;
	receive(v_size, source, tag);
	vec_local.resize((size_t)v_size);

	T obj; // temporary object
	for (size_t ii=0; ii<vec_local.size(); ii++) {
		std::stringstream objs_stream;
		receive(obj_string, source, tag);

		objs_stream << obj_string;

		objs_stream >> obj;
		vec_local[ii] = new T(obj);
	}
}

/**
 * @brief	Send the objects pointed to by vector<T*> to process \#destination
 * @param[in] vec_local A vector of T* pointers to objects that shall be sent
 * @param[in] destination The process rank that will receive the values
 * @param[in] tag Arbitrary non-negative integer assigned to uniquely identify a message
 * @note Class T needs to have the serialize and deseralize operator << and >> implemented
 */
template <class T> void MPIControl::send(const std::vector<T>& vec_local, const size_t& destination, const int& tag)
{
	if ((size_ <= 1) || (rank_ == destination)) return;

	const size_t v_size = vec_local.size();

	send((int)v_size, destination, tag); // first send the size of the vector
	for (size_t ii=0; ii<v_size; ii++) {
		std::string obj_string;
		std::stringstream objs_stream;

		objs_stream << vec_local[ii];
		obj_string.insert(0, objs_stream.str());
		send(obj_string, destination, tag);
	}
}
/**
 * @brief	Receive vector of objects from process \#source
 * @param[in] vec_local A vector of T* pointers to receive the object pointers
 * @param[in] source The process rank that will send the objects
 * @param[in] tag Arbitrary non-negative integer assigned to uniquely identify a message
 * @note Class T needs to have the serialize and deseralize operator << and >> implemented
 */
template <class T> void MPIControl::receive(std::vector<T>& vec_local, const size_t& source, const int& tag)
{
	if ((size_ <= 1) || (rank_ == source)) return;

	if (!vec_local.empty())
		throw mio::IOException("The vector to receive pointers has to be empty (please properly free the vector)", AT);

	std::string obj_string;
	int v_size;
	receive(v_size, source, tag);
	vec_local.resize((size_t)v_size);

	T obj; // temporary object
	for (size_t ii=0; ii<vec_local.size(); ii++) {
		std::stringstream objs_stream;
		receive(obj_string, source, tag);

		objs_stream << obj_string;

		objs_stream >> obj;
		vec_local[ii] = T(obj);
	}
}

// Since template is in cc file (not possible to template on h, becuse if definition is is
// h file, the custom >> and << declaration for std::pair should be in h file, which creates
// conflict with other h files).
template void MPIControl::receive<SnowStation>(std::vector<SnowStation*>&, const size_t&, const int&);
template void MPIControl::send<SnowStation>(const std::vector<SnowStation*>&, const size_t&, const int&);

template void MPIControl::receive<CurrentMeteo>(std::vector<CurrentMeteo*>&, const size_t&, const int&);
template void MPIControl::send<CurrentMeteo>(const std::vector<CurrentMeteo*>&, const size_t&, const int&);

template void MPIControl::receive<SurfaceFluxes>(std::vector<SurfaceFluxes*>&, const size_t&, const int&);
template void MPIControl::send<SurfaceFluxes>(const std::vector<SurfaceFluxes*>&, const size_t&, const int&);

template void MPIControl::receive< std::pair<unsigned long, unsigned long> >(std::vector<std::pair<unsigned long, unsigned long> >&, const size_t&, const int&);
template void MPIControl::send< std::pair<unsigned long, unsigned long> >(const std::vector<std::pair<unsigned long, unsigned long> >&, const size_t&, const int&);

#endif

MPIControl& MPIControl::instance()
{
	static MPIControl instance_;
	return instance_;
}

void MPIControl::getArraySliceParams(const size_t& dimx, const size_t& nbworkers, const size_t& idx_wk, size_t& startx_sub, size_t& nx_sub)
{
	if (nbworkers < dimx) {
		const size_t nx_avg = static_cast<size_t>( floor(static_cast<double>(dimx) / static_cast<double>(nbworkers)) );
		const size_t remainder = dimx % nbworkers;

		if (idx_wk < remainder) { //distribute remainder, 1 extra column per worker, on first workers
			startx_sub = idx_wk * nx_avg + idx_wk;
			nx_sub = nx_avg + 1;
		} else { //all remainder has been distributed, we now attribute a normal number of columns
			startx_sub = idx_wk * nx_avg + remainder;
			nx_sub = nx_avg;
		}
	} else {
		if (idx_wk < dimx) {
			startx_sub = idx_wk;
			nx_sub = 1;
		} else {
			nx_sub = 0;
			startx_sub = 0;
		}
	}
}

/**
 * @brief	Split the domain for MPI, balance the laod between MPI instances taking
 * into account cells where no comutation are required
 * @param[in] dimx size in x of the DEM (not necessary anymore, but retained to have simmilar call than getArraySliceParams)
 * @param[in] idx_wk MPI id of the current instance
 * @param[in] startx_sub Where the number of columns to compute for this instance will be stored
 * @param[in] nx_sub Where the starting point for this instance will be stored
 * @param[in] dem DEM used in the model
 * @param[in] landuse Landuse grid used in the model, nodata in landuse should indicate where no computaiton is required
 */
void MPIControl::getArraySliceParamsOptim(const size_t& dimx, const size_t& idx_wk, size_t& startx_sub, size_t& nx_sub,const mio::DEMObject& dem, const mio::Grid2DObject& landuse)
{
	//Nothing to do is size_=1
	if (size_==1){
		nx_sub = dimx;
		startx_sub = 0;
		return;
	}

	//Check for each col how many cells to compute
	const size_t dimy = dem.getNy();
	size_t n_skip_cell = 0;
	std::vector<size_t> cells_per_col(dimx,0);
	for (size_t ix = 0; ix < dimx; ix++) {
		for (size_t iy = 0; iy < dimy; iy++) {
			if (SnowpackInterfaceWorker::skipThisCell(landuse(ix,iy), dem(ix,iy))) {
				//skip nodata cells as well as water bodies, etc
				n_skip_cell++;
			} else {
				cells_per_col.at(ix)++;
			}
		}
	}

	const size_t num_cells_to_compute = dimx*dimy - n_skip_cell;
	const size_t mean_num_cells_per_mpi = num_cells_to_compute / size_;

	std::vector<size_t> startx(size_,0);
	std::vector<size_t> nx(size_,0);
	std::vector<size_t> n_cells(size_,0);
	size_t current_x=0;
	size_t current_num_cells=0;
	size_t total_num_cells=0;
	for (size_t i=0; i<size_-1;++i) {
		startx.at(i) = current_x;
		current_num_cells = 0;
		//do-while ot force to always add one column
		do {
			current_num_cells += cells_per_col.at(current_x);
			total_num_cells+=cells_per_col.at(current_x);
			current_x++;
		}
		while(current_x < dimx
					&& total_num_cells + cells_per_col.at(current_x) < (i+1.2)*mean_num_cells_per_mpi
					&& (dimx-current_x) > (size_-i));
		//The i+1.2 test is to allow to add a bit more, otherwise last column get all the remainders and is huge
		//The last test is to be sure that at least 1 column per MPI remains
		n_cells.at(i) = current_num_cells;
		//No -1 required because curent_x already incremented once in the "while" loop
		nx.at(i) = current_x-startx.at(i);
	}
	startx.at(size_-1) = current_x;
	current_num_cells = 0;
	//Finish to fill the last slice with remaining cells
	while( current_x < dimx) {
		current_num_cells+=cells_per_col.at(current_x);
		current_x++;
	}
	n_cells.at(size_-1) = current_num_cells;
	nx.at(size_-1) = current_x-startx.at(size_-1);

	startx_sub = startx.at(idx_wk);
	nx_sub = nx.at(idx_wk);
}
