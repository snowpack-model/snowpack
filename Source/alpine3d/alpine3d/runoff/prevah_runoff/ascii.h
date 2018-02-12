/*--------------------------------------------------------------------
 * $Id: ascii.h,v 1.1 2004/07/20 14:47:26 lehning Exp $
 * 
 * This file is part of libRadtran.
 * Copyright (c) 1997-2001 by Arve Kylling and Bernhard Mayer.
 * 
 * ######### Contact info: http://www.libradtran.org #########
 *
 * This program is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License   
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.        
 * 
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of  
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the   
 * GNU General Public License for more details.                    
 * 
 * You should have received a copy of the GNU General Public License          
 * along with this program; if not, write to the Free Software                
 * Foundation, Inc., 59 Temple Place - Suite 330, 
 * Boston, MA 02111-1307, USA.
 *--------------------------------------------------------------------*/

#ifndef ascii_h
#define ascii_h

#if defined (__cplusplus)
extern "C" {
#endif

/* errorcodes of ASCII_functions */

#define ASCIIFILE_NOT_FOUND      -1
#define ASCII_NO_MEMORY          -2
#define LESS_THAN_TWO_COLUMNS    -3
#define LESS_THAN_THREE_COLUMNS  -4



/* prototypes */

int ASCII_checkfile     (char *filename, int *rows, 
			 int *min_columns, int *max_columns, int *max_length);
int ASCII_calloc_string (char ****string, int rows, int columns, int length);
int ASCII_free_string   (char ***string, int rows, int columns);
int ASCII_calloc_int    (int ***value, int rows, int columns);
int ASCII_free_int      (int **value, int rows);
int ASCII_calloc_double (double ***value, int rows, int columns);
int ASCII_free_double   (double **value, int rows);
int ASCII_readfile      (char *filename, char ***array);
int ASCII_string2double (double **value, char *** string, int rows, int columns);
int ASCII_file2double   (char *filename, int *rows, 
			 int *max_columns, int *min_columns, double ***value);
int ASCII_calloc_float  (float ***value, int rows, int columns);
int ASCII_calloc_float_3D(float ****value, int rows, int columns, int length);
int ASCII_calloc_double_3D(double ****value, int rows, int columns, int length);
int ASCII_calloc_double_3D_arylen(double ****value, int rows, int columns, int *length);
int ASCII_calloc_float_4D(float *****value, int rows, int columns, int length, int fourth_dimension);
int ASCII_calloc_double_4D(double *****value, int rows, int columns, int length, int fourth_dimension);
int ASCII_calloc_float_5D(float ******value, int rows, int columns, int length, 
			  int fourth_dimension, int fifth_dimension);
int ASCII_free_float    (float **value, int rows);
int ASCII_free_float_3D (float ***value, int rows, int columns);
int ASCII_free_double_3D (double ***value, int rows, int columns);
int ASCII_free_float_4D (float ****value, int rows, int columns, int length);
int ASCII_free_double_4D(double ****value, int rows, int columns, int length);
int ASCII_free_float_5D (float *****value, int rows, int columns, int length, int fourth_dimension);
int ASCII_file2float    (char *filename, int *rows, 
			 int *max_columns, int *min_columns, float ***value);
int ASCII_parse         (char *string, char *separator, char ***array, int *number);
int ASCII_parsestring   (char *string, char ***array, int *number);
double *ASCII_column    (double **value, int rows, int column);
float  *ASCII_column_float (float **value, int rows, int column);
double *ASCII_row       (double **value, int columns, int row);
int read_1c_file        (char *filename, 
			 double **first, int *n);
int read_2c_file        (char *filename, 
			 double **first, double **second, int *n);
int read_3c_file        (char *filename, 
			 double **first, double **second, double **third, 
			 int *n);
int read_4c_file        (char *filename, 
			 double **first, double **second, double **third, double **fourth,
			 int *n);
int read_5c_file        (char *filename, 
			 double **first, double **second, double **third, 
			 double **fourth, double **fifth, int *n);
int read_6c_file        (char *filename, 
			 double **first, double **second, double **third, 
			 double **fourth, double **fifth, double **sixth, int *n);

int read_1c_file_float (char *filename, 
			float **first, int *n);
int read_2c_file_float (char *filename, 
			float **first, float **second, int *n);
int read_3c_file_float (char *filename, 
			float **first, float **second, float **third, int *n);
int read_4c_file_float (char *filename, 
			float **first, float **second, float **third, float **fourth, int *n);

char *substr            (char *buffer, char *string, int start, int length);

int ASCII_sortarray       (double **value, int rows, int columns, int index, char down);
int ASCII_sortarray_float (float **value,  int rows, int columns, int index, char down);


#if defined (__cplusplus)
}
#endif

#endif
