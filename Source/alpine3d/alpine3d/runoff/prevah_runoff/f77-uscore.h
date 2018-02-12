/*  f77-uscore.h                                          


Copyright (C) 1995 Arve Kylling.

This file is part of liblibRadtran.

libRadtran is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2, or (at your option) any
later version.

libRadtran is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with libRadtran; see the file COPYING.  If not, write to the Free
Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.

*/

#ifndef f77_uscore_h
#define f77_uscore_h 1

#define F77_APPEND_UNDERSCORE 1

#if defined (F77_APPEND_UNDERSCORE)
#define F77_FCN(f) f##_
#else
#define F77_FCN(f) f
#endif

#endif

/*
;;; Local Variables: ***
;;; mode: C ***
;;; End: ***
*/
