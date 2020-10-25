#include <iostream> // cout
#include <fstream> // ofstream
#include <random> // default_random_engine, %distribution
#include <functional> // bind, ref
#include <algorithm> // generate, max, min, abs
#include <vector>
#include <tuple> // tuple, get
#include <thread>
#include <time.h>
#include "program.h"



/*global variables related to random number generator*/
int seed = 5;
std::default_random_engine generator(seed);

/*bernoulli random number generator*/
std::bernoulli_distribution bool_distribution(0.5);
auto randombool = std::bind(bool_distribution, generator);

/*uniform random number generator*/
std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);
auto randomuni = std::bind(uniform_distribution, generator);



int move_rule(Walker& w, const Program& model) {
	uniform_int_distribution<> uniindex(0, w.bits.size() - 1);
	return uniindex(generator);
}

void move_walker(Walker& w, const Program& model) {
	/*set the next move according to the rule*/
	int k = move_rule(w, model);

	/*update the variables according to the move*/
	int j;
	for (const auto& coeff : model.bv[k]) {
		switch (std::get<0>(coeff)) {
		case 'i':
			j = std::get<1>(coeff);
			w.penalty_bv -= w.penalty_ineq[j];
			w.ineq_lhs[j] -= model.ineq_term[j][std::get<2>(coeff)](w.bits);
			break;
		case 'e':
			j = std::get<1>(coeff);
			w.penalty_bv -= w.penalty_eq[j];
			w.eq_lhs[j] -= model.eq_term[j][std::get<2>(coeff)](w.bits);
			break;
		case 'c':
			j = std::get<1>(coeff);
			w.cont_rhs[j] -= model.cont_term[j][std::get<2>(coeff)](w.bits);
			break;
		case 'o':
			w.obj_bv -= model.obj_term[std::get<2>(coeff)](w.bits);
		}
	}

	/*flip the bit and update variables*/
	w.bits[k] = !w.bits[k];

	for (const auto& coeff : model.bv[k]) {
		switch (std::get<0>(coeff)) {
		case 'i':
			j = std::get<1>(coeff);
			w.ineq_lhs[j] += model.ineq_term[j][std::get<2>(coeff)](w.bits);
			w.penalty_ineq[j] = model.ineq_penalty_coeff[j] * std::max(0.0, w.ineq_lhs[j] - model.ineq_rhs[j]);
			w.penalty_bv += w.penalty_ineq[j];
			break;
		case 'e':
			j = std::get<1>(coeff);
			w.eq_lhs[j] += model.eq_term[j][std::get<2>(coeff)](w.bits);
			w.penalty_eq[j] = model.eq_penalty_coeff[j] * std::abs(w.eq_lhs[j] - model.eq_rhs[j]);
			w.penalty_bv += w.penalty_eq[j];
			break;
		case 'c':
			j = std::get<1>(coeff);
			w.cont_rhs[j] += model.cont_term[j][std::get<2>(coeff)](w.bits);
			break;
		case 'o':
			w.obj_bv += model.obj_term[std::get<2>(coeff)](w.bits);
		}
	}


	/*calculate the energy*/
	if (w.contvar.size() == 0) {
		w.energy = w.obj_bv + w.penalty_bv;
	}
	else {
		std::tuple<double, bool, bool> result = solve_subproblem(w.bits, w.cont_rhs, model.cont_penalty_coeff, w.contvar);
		w.energy_cont = std::get<0>(result);
		w.cont_feasible = std::get<1>(result);
		w.energy = w.obj_bv + w.penalty_bv + w.energy_cont;
	}
}



void initialize_walker(Walker& w, const Program& model) {
	/*randomly initialize the bits of walker*/
	const int bitsize = model.bv.size();
	w.bits.resize(bitsize);
	for (int k = 0; k < bitsize; ++k)
		w.bits[k] = randombool();

	w.ineq_lhs.resize(model.ineq_rhs.size());
	w.penalty_ineq.resize(model.ineq_rhs.size());
	w.eq_lhs.resize(model.eq_rhs.size());
	w.penalty_eq.resize(model.eq_rhs.size());
	w.cont_rhs = model.cont_rhsconst;

	/*calculate the left hand side of constraints and objective functions containing only binary variables*/
	/*also calculate the right hand side of constraints containing continuous variables*/
	for (const auto& term : model.obj_term)
		w.obj_bv += term(w.bits);

	for (int j = 0; j < model.ineq_term.size(); ++j) {
		for (const auto& term : model.ineq_term[j])
			w.ineq_lhs[j] += term(w.bits);
	}

	for (int j = 0; j < model.eq_term.size(); ++j) {
		for (const auto& term : model.eq_term[j])
			w.eq_lhs[j] += term(w.bits);
	}

	for (int j = 0; j < model.cont_term.size(); ++j) {
		for (const auto& term : model.cont_term[j])
			w.cont_rhs[j] += term(w.bits);
	}


	/*populate initial penalties corresponding to constraints with only binary variables*/
	for (int j = 0; j < model.ineq_rhs.size(); ++j) {
		w.penalty_ineq[j] = model.ineq_penalty_coeff[j] * std::max(0.0, w.ineq_lhs[j] - model.ineq_rhs[j]);
		w.penalty_bv += w.penalty_ineq[j];
	}
	for (int j = 0; j < model.eq_rhs.size(); ++j) {
		w.penalty_eq[j] = model.eq_penalty_coeff[j] * std::abs(w.eq_lhs[j] - model.eq_rhs[j]);
		w.penalty_bv += w.penalty_eq[j];
	}

	/*assign initial feasible values to continuous variables, if any*/
	set_initial_contvar(w, model);

	/*calculate energy*/
	if (w.contvar.size() == 0) {
		w.energy = w.obj_bv + w.penalty_bv;
	}
	else {
		std::tuple<double, bool, bool> result = solve_subproblem(w.bits, w.cont_rhs, model.cont_penalty_coeff, w.contvar);
		w.energy_cont = std::get<0>(result);
		w.cont_feasible = std::get<1>(result);
		w.energy = w.obj_bv + w.penalty_bv + w.energy_cont;
	}
}



int main() {
	/*populate the optimization model*/
	Program model;
	populate(model);
	if (model.ineq_term.size()) {
		model.ineq_penalty_coeff.resize(model.ineq_term.size());
		std::fill(model.ineq_penalty_coeff.begin(), model.ineq_penalty_coeff.end(), 50.0);
	}
	if (model.eq_term.size()) {
		model.eq_penalty_coeff.resize(model.eq_term.size());
		std::fill(model.eq_penalty_coeff.begin(), model.eq_penalty_coeff.end(), 50.0);
	}
	if (model.cont_term.size()) {
		model.cont_penalty_coeff.resize(model.cont_term.size());
		std::fill(model.cont_penalty_coeff.begin(), model.cont_penalty_coeff.end(), 10000.0);
	}

	const int bitsize = model.bv.size();



	/*initialize walkers*/
	Walker w;
	Walker w_neighbor;
	initialize_walker(w, model);

	/*variables related to recording and outputing the best solution discovered so far*/
	int best_index = -1;
	int record_interval = 1;
	double best_obj = 1000000000000000.0;
	time_t previous_time, current_time;
	previous_time = current_time = time(NULL);
	std::vector<bool> best_solution(bitsize);
	std::ofstream file;

	/*run the simulated annealing*/
	double T{ exp(0.032 * bitsize + 9.3) };
	bool moved = false;
	while (T > 0.1) {
		if (!moved)
			w_neighbor = w;

		move_walker(w_neighbor, model);

		if (exp((w.energy - w_neighbor.energy) / T) >= randomuni()) {
			w = w_neighbor;

			/*record the new best solutions found, if any*/
			if ((w.energy < best_obj) && (w.penalty_bv < 0.1) && w.cont_feasible) {
				best_obj = w.energy;
				std::cout << "The current objective value is " << best_obj << std::endl;

				if (current_time - previous_time >= record_interval) {
					previous_time = current_time = time(NULL);
					file.open("solution.csv");
					file << "objective value," << w.energy << "\n";
					for (int k = 0; k < bitsize; ++k)
						file << "b[" << k << "]," << w.bits[k] << "\n";
					for (int l = 0; l < w.contvar.size(); ++l)
						file << "x[" << l << "]," << w.contvar[l] << "\n";
					file.close();
				}
			}
		}

		/*cooling*/
		T *= 0.99;
	}
}