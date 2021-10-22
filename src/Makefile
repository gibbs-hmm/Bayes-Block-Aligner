# 
# ------ MAKEFILE FILE FOR bayesaligner (on Sun Sparc station) ----------------
# 

CC     =cc 

CFLAGS = -g

PGM    = ../bin/Bayesaligner

LIB    = $(PGM).a

.PRECIOUS :   $(LIB) $(PGM)
 
ELEMENTS =  \
 $(LIB)(Bayesaligner.o)  \
 $(LIB)(Gibbs_NON_GUI.o)  \
 $(LIB)(fill_Sank_rel_matrix.o)\
 $(LIB)(Sankoff.o)\
 $(LIB)(Sankoff_Bayesian_align.o)\
 $(LIB)(Sankoff_get_string.o)\
 $(LIB)(common.o)\
 $(LIB)(mem_mgmt.o)

$(PGM) : $(LIB) common.h
	 $(CC) $(LIB) -lm -o ../bin/Bayesaligner

$(LIB) :  $(ELEMENTS)
