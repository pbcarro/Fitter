import numpy as np
from ctypes import *



#### GSL structs

###Level 0 structures
class GSL_block(Structure):
	_fields_=	[	
		("size",c_size_t),
		("data",POINTER(c_double)),
		]
				
class GSL_vector(Structure):
	_fields_=	[	
		("size",	c_size_t),
		("stride",	c_size_t),
		("data",	POINTER(c_double)),
		("block",	POINTER(GSL_block)),
		("owner",	c_int)
		]

class GSL_matrix(Structure):
	_fields_=	[	
		("size1",	c_size_t),
		("size2",	c_size_t),
		("tda",		c_size_t),
		("data",	POINTER(c_double)),
		("block",	POINTER(GSL_block)),
		("owner",	c_int)
		]

###Level 1 structures

###Function pointer types for Level 1 gsl structures
free_type = CFUNCTYPE(	c_void_p,
						c_void_p
					)

f_type = CFUNCTYPE(	c_int,
					POINTER(GSL_vector),
					c_void_p,
					POINTER(GSL_vector)
					)
fvv_type = CFUNCTYPE(	c_int, 
						POINTER(GSL_vector), 
						POINTER(GSL_vector), 
						POINTER(c_void_p),
						POINTER(GSL_vector)
					)
trs_alloc_type = CFUNCTYPE	(	c_void_p,
							c_void_p,
							c_size_t,
							c_size_t
						)

trs_init_type = CFUNCTYPE(	c_int,
							c_void_p,
							c_void_p
						)

preloop_type = CFUNCTYPE(	c_int,
							c_void_p,
							c_void_p
						)

step_type = CFUNCTYPE	(	c_int,
							c_void_p,
							c_double,POINTER(GSL_vector),
							c_void_p
						)

preduction_type = CFUNCTYPE(c_int,
							c_void_p,
							POINTER(GSL_vector),
							c_double,
							c_void_p
							)

scale_init_type = CFUNCTYPE(c_int,
							POINTER(GSL_matrix),
							POINTER(GSL_vector)
							)

scale_update_type = CFUNCTYPE(	c_int,
								POINTER(GSL_matrix),
								POINTER(GSL_vector)
							)

solver_alloc_type = CFUNCTYPE	(	c_void_p,
							c_void_p,
							c_size_t,
							c_size_t
						)

solver_init_type = CFUNCTYPE(	c_int,
								c_void_p,
								c_void_p
							)

solver_presolve_type = CFUNCTYPE(	c_int,
									c_double,
									c_void_p,
									c_void_p
								)

solver_rcond_type = CFUNCTYPE(	c_int,
								c_double,
								c_void_p
							)

class GSL_fdf(Structure):
	_fields_=	[	
		("f",		f_type),
		("df",		f_type),
		("fvv",		fvv_type),
		("n",		c_size_t),
		("p",		c_size_t),
		("params",	c_void_p),
		("nevalf",	c_size_t),
		("nevaldf",	c_size_t),
		("nevalfvv",c_size_t)
 		]

class GSL_nlinear_trs(Structure):
	_fields_=	[	
		("name",		POINTER(c_char)),
		("alloc",		trs_alloc_type),
		("init",		trs_init_type),
		("preloop",		preloop_type),
		("step",		step_type),
		("preduction",	preduction_type),
		("free",		free_type)
 		]

class GSL_nlinear_scale(Structure):
	_fields_=	[	
		("name",	POINTER(c_char)),
		("init",	scale_init_type),
		("update",	scale_update_type)
 		]

class GSL_nlinear_solver(Structure):
	_fields_=	[	
		("name",	POINTER(c_char)),
		("alloc",	solver_alloc_type),
		("init",	solver_init_type),
		("presolve",solver_presolve_type),
		("rcond",	solver_rcond_type),
		("free",	free_type),
 		]

###Level 2 structures
class GSL_multifit_nlinear_parameters(Structure):
	_fields_=	[	
		("trs",			POINTER(GSL_nlinear_trs)),
		("scale",		POINTER(GSL_nlinear_scale)),
		("solver",		POINTER(GSL_nlinear_solver)),
		("fdtype",		c_int),
		("factor_up",	c_double),
		("factor_down",	c_double),
		("factor_avmax",c_double),
		("h_df",		c_double),
		("h_fvv",		c_double),
 		]

###Function pointer types for Level 2 gsl structures
nlinear_alloc_type = CFUNCTYPE	(	c_void_p,
									POINTER(GSL_multifit_nlinear_parameters),
									c_size_t,
									c_size_t
								)

nlinear_init_type = CFUNCTYPE	(	c_int,
									c_void_p,
									POINTER(GSL_vector),
									POINTER(GSL_fdf),
									POINTER(GSL_vector),
									POINTER(GSL_vector),
									POINTER(GSL_matrix),
									POINTER(GSL_vector)
								)

nlinear_iterate_type = CFUNCTYPE(	c_int,
									c_void_p,
									POINTER(GSL_vector),
									POINTER(GSL_fdf),
									POINTER(GSL_vector),
									POINTER(GSL_vector),
									POINTER(GSL_matrix),
									POINTER(GSL_vector),
									POINTER(GSL_vector)
								)

nlinear_type_rcond_type = CFUNCTYPE(c_int,
									c_double,
									c_void_p
									)

nlinear_avratio_type = CFUNCTYPE(	c_double,
									c_void_p
								)

class GSL_multifit_nlinear_type(Structure):
	_fields_=	[	
		("name",	POINTER(c_char)),
		("alloc",	nlinear_alloc_type),
		("init",	nlinear_init_type),
		("iterate",	nlinear_iterate_type),
		("rcond",	nlinear_type_rcond_type),
		("avratio",	nlinear_avratio_type),
		("free",	free_type)
 		]


class GSL_multifit_nlinear_workspace(Structure):
	_fields_=	[	
		("type",			POINTER(GSL_multifit_nlinear_type)),
		("fdf",				GSL_fdf),
		("x",				POINTER(GSL_vector)),
		("f",				POINTER(GSL_vector)),
		("dx",				POINTER(GSL_vector)),
		("g",				POINTER(GSL_vector)),
		("J",				POINTER(GSL_matrix)),
		("sqrt_wts_work",	POINTER(GSL_vector)),
		("sqrt_wts",		POINTER(GSL_vector)),
		("niter",			c_size_t),
		("params",			GSL_multifit_nlinear_parameters),
		("state",			c_void_p)
 		]	











