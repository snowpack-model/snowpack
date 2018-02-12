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
#include <alpine3d/MPIControl.h>

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
