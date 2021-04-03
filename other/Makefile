SRC_KYNRLINE=xsKYNrline.c
SRC_KYNCONV=xsKYNconv.c
SRC_KYNRLPLI=xsKYNrlpli.c
SRC_KYNCLP=xsKYNclp.c
SRC_KYNBBPHEN=xsKYNBBphen.c
SRC_KYNBB=xsKYNBB.c
SRC_KYNLPCR=xscutoffpl.cxx xsKYNlpcr.c
SRC_KYNREFIONX=xscutoffpl.cxx xsKYNrefionx.c
SRC_KYNXILLVER=xscutoffpl.cxx xsKYNxillver.c
SRC_KYNHREFL=xsKYNhrefl.c

SRC=xside_threads.c xside.c
SRC1=libxspec.c

CC=gcc
INCLUDE=-I/usr/local/share/heasoft-6.24/x86_64-pc-linux-gnu-libc2.31/include/
INCLUDE1=-I/usr/local/share/heasoft/Xspec/src/
LIBRARY_PATH=-L/usr/local/share/heasoft-6.24/x86_64-pc-linux-gnu-libc2.31/lib
LIBRARY= -lcfitsio -lpthread
LIBRARY1= -lXSFunctions -lXSModel -lXSUtil -lXS  -lCCfits_2.5 -lwcs-5.16 -lhdsp_3.0 -lstdc++ -lgfortran

CFLAGS=-fPIC -O3 -Wall -DOUTSIDE_XSPEC -lm

kynrline: $(SRC) $(SRC1) $(SRC_KYNRLINE)
		 $(CC) $(INCLUDE) -o $@ $(SRC) $(SRC1) $(SRC_KYNRLINE) $(LIBRARY_PATH) $(LIBRARY) $(CFLAGS)

kynconv: $(SRC) $(SRC1) $(SRC_KYNCONV)
		$(CC) $(INCLUDE) -o $@ $(SRC) $(SRC1) $(SRC_KYNCONV) $(LIBRARY_PATH) $(LIBRARY) $(CFLAGS)

kynrlpli: $(SRC) $(SRC1) $(SRC_KYNRLPLI)
		 $(CC) $(INCLUDE) -o $@ $(SRC) $(SRC1) $(SRC_KYNRLPLI) $(LIBRARY_PATH) $(LIBRARY) $(CFLAGS)

kynclp: $(SRC) $(SRC1) $(SRC_KYNCLP)
		  $(CC) $(INCLUDE) -o $@ $(SRC) $(SRC1) $(SRC_KYNCLP) $(LIBRARY_PATH) $(LIBRARY) $(CFLAGS)

kynbbphen: $(SRC) $(SRC1) $(SRC_KYNBBPHEN)
		  $(CC) $(INCLUDE) -o $@ $(SRC) $(SRC1) $(SRC_KYNBBPHEN) $(LIBRARY_PATH) $(LIBRARY) $(CFLAGS)

kynbb: $(SRC) $(SRC1) $(SRC_KYNBB)
	      $(CC) $(INCLUDE) -o $@ $(SRC) $(SRC1) $(SRC_KYNBB) $(LIBRARY_PATH) $(LIBRARY) $(CFLAGS)

kynlpcr: $(SRC) $(SRC1) $(SRC_KYNLPCR)
		$(CC) $(INCLUDE) $(INCLUDE1) -o $@ $(SRC) $(SRC1) $(SRC_KYNLPCR) $(LIBRARY_PATH) $(LIBRARY1) $(LIBRARY) $(CFLAGS)

kynreflionx: $(SRC) $(SRC1) $(SRC_KYNREFIONX)
		$(CC) $(INCLUDE) $(INCLUDE1) -o $@ $(SRC) $(SRC1) $(SRC_KYNREFIONX) $(LIBRARY_PATH) $(LIBRARY1) $(LIBRARY) $(CFLAGS)

kynxillver: $(SRC) $(SRC1) $(SRC_KYNXILLVER)
		$(CC) $(INCLUDE) $(INCLUDE1) -o $@ $(SRC) $(SRC1) $(SRC_KYNXILLVER) $(LIBRARY_PATH) $(LIBRARY1) $(LIBRARY) $(CFLAGS)

kynhrefl: $(SRC) $(SRC1) $(SRC_KYNHREFL)
		 $(CC) $(INCLUDE) -o $@ $(SRC) $(SRC1) $(SRC_KYNHREFL) $(LIBRARY_PATH) $(LIBRARY) $(CFLAGS)
