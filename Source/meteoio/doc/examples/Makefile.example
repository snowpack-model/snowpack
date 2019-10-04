CC      = g++ -DGNU
CFLAGS  = -Wall -Wextra
DEBUG   = -g #-O3 # -DDEBUG -ggdb
METEOIODIR = ../../

LIBS      = -rdynamic -lstdc++ -ldl -L../../lib -lmeteoio -ldl

INCLUDE=-I. -I$(METEOIODIR) -I$(METEOIODIR)/include
INCLUDE_POPC=-I. -I$(METEOIODIR) -I$(METEOIODIR)/include

#####################
#	RULES
#####################
.cc.o: $*.cc $*.h
	$(CC) $(DEBUG) $(CFLAGS) -c $< $(INCLUDE)

%.o: %.cc
	$(CC) $(DEBUG) $(CFLAGS) -c $< $(INCLUDE)

#####################
#	TARGETS
#####################
all: meteo_reading dem_reading grid2d_reading coordinates 2D_interpolations time data_converter matrix sun random_numbers statistical_filters

meteo_reading: meteo_reading.o
	$(CC) $(DEBUG) -o $@ meteo_reading.o ${LIBS}
dem_reading: dem_reading.o
	$(CC) $(DEBUG) -o $@ dem_reading.o ${LIBS}
grid2d_reading: grid2d_reading.o
	$(CC) $(DEBUG) -o $@ grid2d_reading.o ${LIBS}
coordinates: coordinates.o
	$(CC) $(DEBUG) -o $@ coordinates.o ${LIBS}
2D_interpolations: 2D_interpolations.o
	$(CC) $(DEBUG) -o $@ 2D_interpolations.o ${LIBS}
time: time.o
	$(CC) $(DEBUG) -o $@ time.o ${LIBS}
data_converter: data_converter.o
	$(CC) $(DEBUG) -o $@ data_converter.o ${LIBS}
matrix: matrix.o
	$(CC) $(DEBUG) -o $@ matrix.o ${LIBS}
sun: sun.o
	$(CC) $(DEBUG) -o $@ sun.o ${LIBS}
random_numbers: random_numbers.o
	$(CC) $(DEBUG) -o $@ random_numbers.o ${LIBS}
statistical_filters: statistical_filters.o
	$(CC) $(DEBUG) -o $@ statistical_filters.o ${LIBS}

clean:
	rm -rf *~ *.o  *.out meteo_reading dem_reading grid2d_reading coordinates 2D_interpolations time data_converter matrix sun random_numbers statistical_filters

