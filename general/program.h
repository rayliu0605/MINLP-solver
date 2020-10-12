#include <functional>
#include <vector>
#include <tuple>

using namespace std;

struct Program {
	vector<function<double(vector<bool>)>> obj_term;
	vector<vector<function<double(vector<bool>)>>> ineq_term;
	vector<double> ineq_rhs;
	vector<double> ineq_penalty_coeff;
	vector<vector<function<double(vector<bool>)>>> eq_term;
	vector<double> eq_rhs;
	vector<double> eq_penalty_coeff;
	vector<vector<function<double(vector<bool>)>>> cont_term;
	vector<double> cont_rhsconst;
	vector<double> cont_penalty_coeff;
	vector<vector<tuple<char, int, int>>> bv;
};

struct Walker {
	vector<bool> bits;
	vector<double> contvar;
	double obj_bv = 0.0;
	double penalty_bv = 0.0;
	double energy_cont = 0.0;
	bool cont_feasible = true;
	vector<double> penalty_ineq;
	vector<double> penalty_eq;
	vector<double> ineq_lhs;
	vector<double> eq_lhs;
	vector<double> cont_rhs;
	double energy = 0.0;
};

void populate(Program& model);

void set_initial_contvar(Walker& w, const Program& model);

tuple<double, bool, bool> solve_subproblem(const vector<bool>& b, const vector<double>& const_rhs, const vector<double>& cont_penalty_coeff, vector<double>& contvar);