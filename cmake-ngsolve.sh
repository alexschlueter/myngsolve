cmake \
    ../../src/ngsolve \
    -DCMAKE_BUILD_TYPE:STRING="RelWithDebInfo" \
    -DINSTALL_DIR:PATH="/home/alex/code/ngsolve/inst-rel" \
    -DUSE_MKL:BOOL="1" \
    -DUSE_MUMPS:BOOL="1" \
    -DUSE_MPI:BOOL="1" \
    -DUSE_CCACHE:BOOL="1" \
    -DUSE_UMFPACK:BOOL="1" \
    -DNETGEN_SOURCE_DIR:PATH="/home/alex/code/ngsolve/src/netgen"
    -DUMFPACK_INCLUDE_DIR:PATH="/usr/include/suitesparse" \
    -DUMFPACK_CHOLMOD_LIB:FILEPATH="/usr/lib/x86_64-linux-gnu/libcholmod.so" \
    -DLIB_MPIF90:FILEPATH="/usr/lib/libmpi_mpifh.so" \
    -DUMFPACK_CONFIG_LIB:FILEPATH="/usr/lib/x86_64-linux-gnu/libsuitesparseconfig.so" \
    -DUMFPACK_LIBRARIES:STRING="/usr/lib/x86_64-linux-gnu/libumfpack.so;/usr/lib/x86_64-linux-gnu/libamd.so;/usr/lib/x86_64-linux-gnu/libcolamd.so;/usr/lib/x86_64-linux-gnu/libcholmod.so;/usr/lib/x86_64-linux-gnu/libsuitesparseconfig.so" \
    -DUMFPACK_LIB:FILEPATH="/usr/lib/x86_64-linux-gnu/libumfpack.so" \
    -DUMFPACK_COLAMD_LIB:FILEPATH="/usr/lib/x86_64-linux-gnu/libcolamd.so" \
    -DMKL_ROOT:PATH="/home/alex/intel/mkl" \
    -DUMFPACK_AMD_LIB:FILEPATH="/usr/lib/x86_64-linux-gnu/libamd.so"
