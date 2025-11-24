lib_LIBRARIES = libgpdf.a
libgpdf_a_SOURCES = gpdf_inmzpd.c wofzwh.f wofzwh2.f cpolyev.f cpolyev2.f utils.f inmzpd.f epszpd.f cbj0.f
libgpdf_a_CFLAGS = -fPIC -I./
libgpdf_a_FFLAGS = -I./ -fPIC -fimplicit-none -fbounds-check
libgpdf_a_INCLUDES=$(all_includes)
EXTRA_DIST=gpdf.py
include_HEADERS=gpdf.h
noinst_HEADERS=doxy.cfg
SUBDIRS = quadp
ACLOCAL_AMFLAGS = -I m4

noinst_PROGRAMS = ftest ctest test

ftest_SOURCES = ftest.f
ftest_LDADD = -L./ -L./quadp -lgpdf -lquadp -lgfortran
ftest_CFLAGS = -I./
ftest_FFLAGS = -I./ -fPIC -fimplicit-none -fbounds-check
ftest_DEPENDENCIES = ./libgpdf.a

test_SOURCES = test.f
test_LDADD = -L./ -L./quadp -lgpdf -lquadp -lgfortran
test_CFLAGS = -I./
test_FFLAGS = -I./ -fPIC -fimplicit-none -fbounds-check
test_DEPENDENCIES = ./libgpdf.a

ctest_SOURCES = ctest.c
ctest_LDADD = -L./ -L./quadp -lgpdf -lquadp -lgfortran
ctest_CFLAGS = -I./
ctest_FFLAGS = -I./ -fPIC -fimplicit-none -fbounds-check
ctest_DEPENDENCIES = ./libgpdf.a

if ZPDGEN_PYTHON_ON
bin_SCRIPTS=inmzpd.pyf epszpd.pyf
PYLD=-L$(abs_srcdir) -L$(abs_srcdir)/quadp -lgpdf -lquadp -lgfortran
endif
CLEANFILES=$(bin_SCRIPTS) *.so

inmzpd.pyf: inmzpd.f Makefile.am
	$(F2PY) --overwrite-signature -m inmzpd -h inmzpd.pyf only : inmweid gmweid : inmzpd.f
	sed 's/:: res/,intent(out) :: res/' <inmzpd.pyf | sed 's/:: za/, intent(in) :: za/'| sed 's/:: b/, intent(in) :: b/'| sed 's/:: zb/, intent(in) :: zb/'| sed 's/:: n/, intent(in) :: n/'| sed 's/:: m/, intent(in) :: m/'| sed 's/:: z1/, intent(in) :: z1/'| sed 's/:: z2/, intent(in) :: z2/'>tmpout
	mv tmpout inmzpd.pyf
	$(F2PY) -c inmzpd.pyf $(PYLD)

epszpd.pyf: epszpd.f Makefile.am
	$(F2PY) --overwrite-signature -m epszpd -h epszpd.pyf only : epsweid sigweid : epszpd.f
	sed 's/:: res/,intent(out) :: res/' <epszpd.pyf | sed 's/:: om/, intent(in) :: om/'| sed 's/:: pars/, intent(in) :: pars/'| sed 's/:: za/, intent(in) :: za/'| sed 's/:: b/, intent(in) :: b/'| sed 's/:: zb/, intent(in) :: zb/'>tmpout
	mv tmpout epszpd.pyf
	$(F2PY) -c epszpd.pyf $(PYLD)
