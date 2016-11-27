FC = gfortran
#FC=ifort

# 日付
TODAY=`date +%y-%m-%d_%Hh%Mm`

# デバッグオプション(y or n)
CHECK = y

#OpenMP による並列化（y or n）
OMP = n

#コンパイルするファイル その１（以下の拡張子をfとしたもの）
OBJS = 	Lindhard90.o src/LoadHopping.o Utils/MenuUtils.o Utils/Utilities.o\
#        zheev.o dsyev.o driver.o zgesv.o

#コンパイルするファイル その２（以下の拡張子をfとしたもの）
OBJ2 = Common.o

SRCS=${OBJS:.o=.f}
SRC2=${OBJ2:.o=.f}

# コンパイルオプション設定：環境に合わせて適宜変更する
#### for gfortran
ifeq ($(FC),gfortran)
  FFLAGS= -O3 -mcmodel=large
  ifeq ($(CHECK),y)
    FFLAGS+= -fbounds-check
  endif
  LIBS= ~/lib/LAPACK/liblapack.a ~/lib/LAPACK/librefblas.a
endif

#### for ifort
ifeq ($(FC),ifort)
  FFLAGS= -O3 -mcmodel=large -i-dynamic
  ifeq ($(OMP),y)
    FFLAGS+= -openmp
  endif
  ifeq ($(CHECK),y)
    FFLAGS+= -par_report3
#    FFLAGS+= -par_threshold0
    FFLAGS+= -g
  endif
  LIBS=-mkl
endif


# コンパイル
PROG = Lindhard90
all:	$(PROG)

$(PROG): $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) $(OBJ2) $(LIBS)  -o $(PROG).exe


${OBJ2} : ${SRC2}
${OBJS} : ${OBJ2}

# 不要なファイルの削除
#clean:
#	rm -f $(OBJS) $(OBJ2) *% core

# 全ファイルを一つにまとめて保存
#save:
#	cat  $(SRCS) $(SRC2) > allcode$(TODAY).f
