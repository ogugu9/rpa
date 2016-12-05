FC = gfortran
#FC=ifort

# 日付
TODAY=`date +%y-%m-%d_%Hh%Mm`

# デバッグオプション(y or n)
CHECK = y

#OpenMP による並列化（y or n）
OMP = n

#コンパイルするファイル その１（以下の拡張子をfとしたもの）
OBJS = 	Lindhard90.o src/LoadHopping.o src/LoadKpath.o \
	Utils/MenuUtils.o Utils/Utilities.o\
	src/InitConfig.o src/MeanField.o src/CalcChemical.o\
	src/CalcEigenvalue.o src/Hamiltonian.o\
	src/OutEigenvalue.o\

#コンパイルするファイル その２（以下の拡張子をfとしたもの）
OBJ2 = Common.o

SRCS=${OBJS:.o=.f}
SRC2=${OBJ2:.o=.f}
MOD=${OBJ2:.o=.mod}

#### for gfortran
ifeq ($(FC),gfortran)
  #FFLAGS= -O3 -mcmodel=large
  ifeq ($(CHECK),y)
    FFLAGS+= -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow
  endif
  LIBS= ~/lib/LAPACK/liblapack.a ~/lib/LAPACK/librefblas.a
endif

# コンパイル
PROG = Lindhard90
all:	$(PROG)

${OBJ2} : ${SRC2}
${MOD} : ${SRC2} ${OBJ2}
${OBJS} : ${OBJ2} ${MOD}

$(PROG): $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) $(OBJ2) $(LIBS)  -o $(PROG).exe

# 不要なファイルの削除
#clean:
#	rm -f $(OBJS) $(OBJ2) *% core

# 全ファイルを一つにまとめて保存
#save:
#	cat  $(SRCS) $(SRC2) > allcode$(TODAY).f
