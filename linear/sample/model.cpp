#include <vector>
#include <tuple>
#include <random>
#include "program.h"

using namespace std;

void populate(Program& model) {
	model.obj_coeff = { };
	model.ineq_coeff = { { 5.0 } };
	model.ineq_rhs = { 1.0 };
	model.eq_coeff = { { 1.0 } };
	model.eq_rhs = { 0.0 };
	model.cont_coeff = { { 1.0 }, { }, { -1.0 } };
	model.cont_rhsconst = { 0.0, -10.0, 0.5 };
	model.bv = { { make_tuple('c', 0, 0), make_tuple('i', 0, 0), make_tuple('e', 0, 0), make_tuple('c', 2, 0) } };
}

void set_initial_contvar(Walker& w, const Program& model) {
	w.contvar = { 0.0 };
}