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
INCLUDE=-I/usr/local/share/heasoft/x86_64-unknown-linux-gnu-libc2.19/include
INCLUDE1=-I/usr/local/share/heasoft/Xspec/src
LIBRARY_PATH=-L/usr/local/share/heasoft/x86_64-unknown-linux-gnu-libc2.19/lib
LIBRARY= -lcfitsio_3.37 -lpthread -lm
LIBRARY1=-lXSFunctions -lXSModel -lXSUtil -lXS  -lCCfits_2.4 -lwcs-4.20 -lstdc++ -lgfortran

CFLAGS=-O3 -Wall -DOUTSIDE_XSPEC

kynrline: $(SRC) $(SRC1) $(SRC_KYNRLINE)
		 $(CC) $(CFLAGS) $(INCLUDE) -o $@ $(SRC) $(SRC1) $(SRC_KYNRLINE) $(LIBRARY_PATH) $(LIBRARY)

kynconv: $(SRC) $(SRC1) $(SRC_KYNCONV)
		$(CC) $(CFLAGS) $(INCLUDE) -o $@ $(SRC) $(SRC1) $(SRC_KYNCONV) $(LIBRARY_PATH) $(LIBRARY)

kynrlpli: $(SRC) $(SRC1) $(SRC_KYNRLPLI)
		 $(CC) $(CFLAGS) $(INCLUDE) -o $@ $(SRC) $(SRC1) $(SRC_KYNRLPLI) $(LIBRARY_PATH) $(LIBRARY)

kynclp: $(SRC) $(SRC1) $(SRC_KYNCLP)
		  $(CC) $(CFLAGS) $(INCLUDE) -o $@ $(SRC) $(SRC1) $(SRC_KYNCLP) $(LIBRARY_PATH) $(LIBRARY)

kynbbphen: $(SRC) $(SRC1) $(SRC_KYNBBPHEN)
		  $(CC) $(CFLAGS) $(INCLUDE) -o $@ $(SRC) $(SRC1) $(SRC_KYNBBPHEN) $(LIBRARY_PATH) $(LIBRARY)

kynbb: $(SRC) $(SRC1) $(SRC_KYNBB)
	      $(CC) $(CFLAGS) $(INCLUDE) -o $@ $(SRC) $(SRC1) $(SRC_KYNBB) $(LIBRARY_PATH) $(LIBRARY)

kynlpcr: $(SRC) $(SRC_KYNLPCR)
		$(CC) $(CFLAGS) $(INCLUDE) $(INCLUDE1) -o $@ $(SRC) $(SRC_KYNLPCR) $(LIBRARY_PATH) $(LIBRARY1) $(LIBRARY)

kynreflionx: $(SRC) $(SRC_KYNREFIONX)
		$(CC) $(CFLAGS) $(INCLUDE) $(INCLUDE1) -o $@ $(SRC) $(SRC_KYNREFIONX) $(LIBRARY_PATH) $(LIBRARY1) $(LIBRARY)

kynxillver: $(SRC) $(SRC_KYNXILLVER)
		$(CC) $(CFLAGS) $(INCLUDE) $(INCLUDE1) -o $@ $(SRC) $(SRC_KYNXILLVER) $(LIBRARY_PATH) $(LIBRARY1) $(LIBRARY)

kynhrefl: $(SRC) $(SRC1) $(SRC_KYNHREFL)
		 $(CC) $(CFLAGS) $(INCLUDE) -o $@ $(SRC) $(SRC1) $(SRC_KYNHREFL) $(LIBRARY_PATH) $(LIBRARY)
