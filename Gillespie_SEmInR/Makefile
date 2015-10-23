##################################################################
####
####   MAKEFILE FOR GILLESPIE SEmInR
####
##################################################################



# ==================================
#     === Compiling Variables ===
# ==================================

C_CLANG := clang++

C_CLANG_OPT := -Wall -O3 -Wno-predefined-identifier-outside-function
# -Wno-deprecated

GCC := g++
GCC_FLAGS := -Wall -O3 -std=c++11 -Wno-sign-compare

# ==================
# === FILE LISTS ===
# ==================

SOURCE_LIST := mc.cpp simulator.cpp individual.cpp dcTools.cpp RV.cpp globalVar.cpp dcMatrix.cpp

OBJ_LIST_U := main_unit.cpp $(SOURCE_LIST)
OBJ_LIST_SF := launch_serialfarming.cpp $(SOURCE_LIST)
OBJ_LIST_SF_GCC := launch_serialfarming_gcc.cpp $(SOURCE_LIST)


# ======================
# ==== COMPILE ONLY ====
# ======================

PROG_NAME := gilSEIR
PROG_NAME_U := gilSEIR_U
PROG_NAME_SF := gilSEIR_SF

PROG_NAME_GCC := gilSEIR_GCC
PROG_NAME_U_GCC := gilSEIR_U_GCC
PROG_NAME_SF_GCC := gilSEIR_SF_GCC


# ==== CLANG COMPLIER ====

$(PROG_NAME_U): $(OBJ_LIST_U)
	$(C_CLANG) $(C_CLANG_OPT) $(OBJ_LIST_U) -o $@

$(PROG_NAME_SF): $(OBJ_LIST_SF)
	$(C_CLANG) $(C_CLANG_OPT) $(OBJ_LIST_SF) -o $@



# ==== GCC COMPLIER ====

$(PROG_NAME_U_GCC): $(OBJ_LIST_U)
	$(GCC) $(GCC_FLAGS) $(OBJ_LIST_U) -o $@

$(PROG_NAME_SF_GCC): $(OBJ_LIST_SF_GCC)
	$(GCC) $(GCC_FLAGS) $(OBJ_LIST_SF_GCC) -o $@



# ===================
# === EXECUTABLES ===
# ===================

runSimul: $(PROG_NAME)
	./clean
	./$< &

launchsf: $(PROG_NAME_SF) $(PROG_NAME_U)
	./clean
	./$< &

launchsf_gcc: $(PROG_NAME_SF_GCC) $(PROG_NAME_U_GCC)
	./clean
	./$< &



PRM_LIST := simul_param.csv simul_param_R0.csv simul_param_DOI.csv simul_param_DOL.csv

simul_param.out: generate_param_files.R $(PRM_LIST) 
	Rscript $<

simul.out: $(PROG_NAME_SF) $(PROG_NAME_U) simul_param.out
	./clean
	./$(PROG_NAME_SF)

simul_gcc.out: $(PROG_NAME_SF_GCC) $(PROG_NAME_U_GCC) simul_param.out
	./clean
	./$(PROG_NAME_SF_GCC)


# --------------------------------------------------
# --- Specific for "Generation Interval" project ---

# Execute: make simul.out        (<--runs the simulations)
# --------------------------------------------------
# --------------------------------------------------
