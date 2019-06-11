#include <alpine3d/OMPControl.h>
#include <alpine3d/SnowpackInterfaceWorker.h>


namespace OMPControl
{
  void getArraySliceParams(const size_t& dimx, const size_t& nbworkers, const size_t& idx_wk, size_t& startx_sub, size_t& nx_sub)
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
  void getArraySliceParamsOptim(const size_t& nbworkers, const std::vector<SnowStation*>& snow_station, const mio::DEMObject& dem,
                                const mio::Grid2DObject& landuse, std::vector<std::vector<size_t> >& omp_snow_stations_ind)
  {
    size_t dimx = dem.getNx();
    size_t dimy = dem.getNy();

    size_t n_skip_cell=0;
    size_t n_skip_cell_2=0;

    for (size_t ix = 0; ix < dimx; ix++) {
      for (size_t iy = 0; iy < dimy; iy++) {
        if (SnowpackInterfaceWorker::skipThisCell(landuse(ix,iy), dem(ix,iy)))
        {n_skip_cell++;}
      }
    }
    for(size_t i=0;i<snow_station.size();++i)
    {
      if(snow_station.at(i)==NULL)
      {n_skip_cell_2++;}
    }

    size_t n_snow_station_compute = (snow_station.size()-n_skip_cell);
    size_t n_snow_station_worker=n_snow_station_compute/nbworkers;
    size_t remainders = n_snow_station_compute%nbworkers;

    omp_snow_stations_ind.resize(nbworkers);

    size_t worker_i=0;
    size_t worker_n=0;
    size_t remain = remainders>0?1:0;

    for(size_t i=0;i<snow_station.size();++i)
    {
      omp_snow_stations_ind.at(worker_i).push_back(i);
      if(snow_station.at(i)!=NULL)
      {
        worker_n++;
      }
      if(worker_n == n_snow_station_worker+remain)
      {
        remainders=remainders>0?remainders-1:0;
        remain=remainders>0?1:0;
        worker_i=worker_i<(nbworkers-1)?worker_i+1:worker_i;
        worker_n=0;
      }
    }
  }
}
