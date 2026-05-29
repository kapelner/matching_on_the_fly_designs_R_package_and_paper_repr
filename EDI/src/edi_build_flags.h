#ifndef EDI_BUILD_FLAGS_H
#define EDI_BUILD_FLAGS_H

#define EDI_BUILD_CAPTURE_METHOD "configure-generated header compiled into EDI.so"
#define EDI_BUILD_TIMESTAMP "2026-05-29 01:28:40 EDT"
#define EDI_BUILD_HOST "LAPTOP-J2T9TGGB"
#define EDI_BUILD_R_HOME "/usr/local/lib/R"
#define EDI_BUILD_R_VERSION "R Under development (unstable) (2026-04-23 r89955) -- \"Unsuffered Consequences\""
#define EDI_BUILD_R_CXX20 "g++ -O3 -march=native -flto -fno-math-errno"
#define EDI_BUILD_R_CXX20STD "-std=gnu++20"
#define EDI_BUILD_R_CXX20FLAGS "-g -O2 -UNDEBUG -Wall -pedantic -g -O0"
#define EDI_BUILD_R_SHLIB_OPENMP_CXXFLAGS "unavailable"
#define EDI_BUILD_ENV_EDI_PORTABLE "0"
#define EDI_BUILD_ENV_EDI_DISABLE_VECTORIZATION "0"
#define EDI_BUILD_ENV_EDI_NATIVE_SPEED "1"
#define EDI_BUILD_ENV_EDI_NATIVE_LTO "0"
#define EDI_BUILD_PKG_CPPFLAGS "-I../inst/include"
#define EDI_BUILD_PKG_CXXFLAGS "$(SHLIB_OPENMP_CXXFLAGS) -DNDEBUG -DEIGEN_NO_DEBUG -Wno-ignored-attributes -march=native -mtune=native -fno-lto; override CXXFLAGS+=-O3"
#define EDI_BUILD_PKG_LIBS "$(SHLIB_OPENMP_CXXFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) -ltbb12 -fstack-protector"

#endif
