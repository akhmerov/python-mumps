cdef extern from "mumps_c_types.h":
    ctypedef int MUMPS_INT

    ctypedef float SMUMPS_REAL
    ctypedef float SMUMPS_COMPLEX

    ctypedef double DMUMPS_REAL
    ctypedef double DMUMPS_COMPLEX

    ctypedef struct mumps_complex:
        pass

    ctypedef struct mumps_double_complex:
        pass

    ctypedef float CMUMPS_REAL
    ctypedef mumps_complex CMUMPS_COMPLEX

    ctypedef double ZMUMPS_REAL
    ctypedef mumps_double_complex ZMUMPS_COMPLEX


cdef extern from "zmumps_c.h":
    ctypedef struct ZMUMPS_STRUC_C:
        MUMPS_INT sym, par, job
        MUMPS_INT comm_fortran
        MUMPS_INT icntl[40]
        ZMUMPS_REAL cntl[15]

        MUMPS_INT n

        MUMPS_INT nz
        MUMPS_INT *irn
        MUMPS_INT *jcn
        ZMUMPS_COMPLEX *a

        MUMPS_INT nrhs, lrhs
        ZMUMPS_COMPLEX *rhs

        MUMPS_INT info[40]
        MUMPS_INT infog[40]
        ZMUMPS_REAL rinfo[40]
        ZMUMPS_REAL rinfog[40]

        MUMPS_INT nz_rhs
        ZMUMPS_COMPLEX *rhs_sparse
        MUMPS_INT *irhs_sparse
        MUMPS_INT *irhs_ptr

        MUMPS_INT size_schur
        MUMPS_INT *listvar_schur
        ZMUMPS_COMPLEX *schur

    cdef void zmumps_c(ZMUMPS_STRUC_C *)
