DEST = .
CC = gcc
CFLAGS = $(INCLUDE) $(LIB) 

JNI_BLAS_SRC = .

INCLUDE = -I../../include -I/usr/include/gsl

LIB = -L/usr/local/lib #-L/usr/local/lib64 -L/usr/lib 

C_SRC_BLAS = $(JNI_BLAS_SRC)/c/Matrix.c
C_SRC_LU = $(JNI_BLAS_SRC)/c/LUDecomposition.c
C_SRC_CHOL = $(JNI_BLAS_SRC)/c/CholeskyDecomposition.c
C_SRC_EIG = $(JNI_BLAS_SRC)/c/EigenvalueDecomposition.c
C_SRC_QR = $(JNI_BLAS_SRC)/c/QRDecomposition.c
C_SRC_SVD = $(JNI_BLAS_SRC)/c/SingularValueDecomposition.c
C_SRC_RNG = $(JNI_BLAS_SRC)/c/JniGslRng.c
DYLIB_BLAS = ../lib/libgsl_blas.dylib
DYLIB_LU = ../lib/libgsl_lu.dylib
DYLIB_CHOL = ../lib/libgsl_chol.dylib
DYLIB_EIG = ../lib/libgsl_eig.dylib
DYLIB_QR = ../lib/libgsl_qr.dylib
DYLIB_SVD = ../lib/libgsl_svd.dylib
DYLIB_RNG = ../lib/libgsl_random.dylib

JC = javac
JFLAGS = -d
#TODO: update with other package specifications. Don't forget to rename the files in the java to point to the new DYLIB specifications
JAVA_SRC = $(JNI_BLAS_SRC)/java/Matrix.java $(JNI_BLAS_SRC)/java/LUDecomposition.java $(JNI_BLAS_SRC)/java/CholeskyDecomposition.java $(JNI_BLAS_SRC)/java/EigenvalueDecomposition.java $(JNI_BLAS_SRC)/java/QRDecomposition.java $(JNI_BLAS_SRC)/java/SingularValueDecomposition.java $(JNI_BLAS_SRC)/java/JniGslRng.java $(JNI_BLAS_SRC)/java/BayesNIG.java $(JNI_BLAS_SRC)/java/BayesNIW.java $(JNI_BLAS_SRC)/java/fileinput.java
PACKAGE = ../class/JGSL

all: $(PACKAGE) $(DYLIB_BLAS) $(DYLIB_LU) $(DYLIB_CHOL) $(DYLIB_EIG) $(DYLIB_QR) $(DYLIB_SVD) $(DYLIB_RNG) $(DYLIB_BNIG) $(DYLIB_BNIW)

$(PACKAGE): $(JAVA_SRC)
	@mkdir -p $(@D)
	$(JC) $(JFLAGS) ../class $(JAVA_SRC)
	
$(DYLIB_BLAS): $(C_SRC_BLAS) 
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -shared -fPIC $(C_SRC_BLAS) -lgsl -lgslcblas -lm -o $(DYLIB_BLAS)

$(DYLIB_LU): $(C_SRC_LU) 
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -shared -fPIC $(C_SRC_LU) -lgsl -lgslcblas -lm -o $(DYLIB_LU)

$(DYLIB_CHOL): $(C_SRC_CHOL) 
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -shared -fPIC $(C_SRC_CHOL) -lgsl -lgslcblas -lm -o $(DYLIB_CHOL)

$(DYLIB_EIG): $(C_SRC_EIG) 
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -shared -fPIC $(C_SRC_EIG) -lgsl -lgslcblas -lm -o $(DYLIB_EIG)

$(DYLIB_QR): $(C_SRC_QR) 
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -shared -fPIC $(C_SRC_QR) -lgsl -lgslcblas -lm -o $(DYLIB_QR)

$(DYLIB_SVD): $(C_SRC_SVD) 
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -shared -fPIC $(C_SRC_SVD) -lgsl -lgslcblas -lm -o $(DYLIB_SVD)

$(DYLIB_RNG): $(C_SRC_RNG)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -shared -fPIC $(C_SRC_RNG) -lgsl -lgslcblas -lm -o $(DYLIB_RNG)

clean: 
	rm -fr ../class
	rm -fr ../lib
