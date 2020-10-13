#include <cppad/cppad.hpp>
#include <cppad_ipopt_nlp.hpp>
#include <vector>
#include <algorithm>
#include <tuple>

using namespace std;

namespace {
	using namespace cppad_ipopt;

	class FG_info : public cppad_ipopt_fg_info
	{
	private:
		bool retape_{ false };
		vector<bool> b;
		vector<double> mu;
	public:
		// derived class part of constructor
		FG_info(vector<bool> bits, vector<double> cont_penalty_coeff) : b(bits), mu(cont_penalty_coeff) {}

		// evaluation of the objective f(x), and constraints g(x)
		// using an algorithmic differentiation (AD) class
		ADVector eval_r(size_t k, const ADVector& x_i){
			ADVector fg(4);

			vector<ADNumber> x;
			for (size_t i = 0; i < 5; ++i)
				x.push_back(x_i[i]);

			fg[0] = 2.0 * x[0] + mu[0] * x[1] + mu[1] * x[2] + mu[2] * x[3] + mu[1] * x[4];

			fg[1] = x[0] - x[1];
			fg[2] = x[0] - x[2] + x[4];
			fg[3] = b[0] * x[0] - x[3];

			return fg;
		}

		bool retape(size_t k){
			return retape_;
		}
	};
}

tuple<double, bool, bool> solve_subproblem(const vector<bool>& b, const vector<double>& cont_rhs, const vector<double>& mu, vector<double>& contvar) {
	bool ok{ true };
	size_t i;
	double temp;


	// number of original variables
	size_t n_x = 1;
	// number of constraints
	size_t m = 3;
	// number of equality constraints
	size_t me = 1;
	// number of original and auxiliary variables
	size_t n = m + me + n_x;


	// initial value of variables
	NumberVector x(n);
	for (i = 0; i < n_x; ++i)
		x[i] = contvar[i];

	x[1] = max(0.0, x[0] - cont_rhs[0]);
	temp = x[0] - cont_rhs[1];
	x[2] = max(0.0, temp);
	x[4] = max(0.0, -temp);
	x[3] = max(0.0, b[0] * x[0] - cont_rhs[2]);

	// lower and upper limits for x
	NumberVector x_l(n);
	NumberVector x_u(n);
	x_l[0] = -20.0;	 x_u[0] = 50.0;
	for (i = n_x; i < n; ++i) {
		x_l[i] = 0.0;	x_u[i] = 1.0e19;
	}


	// lower and upper limits for g
	NumberVector g_l(m);
	NumberVector g_u(m);
	g_l[0] = -1.0e19;	g_u[0] = cont_rhs[0];
	g_l[1] = cont_rhs[1];	g_u[1] = cont_rhs[1];
	g_l[2] = -1.0e19;	g_u[2] = cont_rhs[2];


	// create the Ipopt interface
	FG_info fg_info(b, mu);

	// create the Ipopt interface
	cppad_ipopt_solution s;
	Ipopt::SmartPtr<Ipopt::TNLP> cppad_nlp = new cppad_ipopt_nlp(
		n, m, x, x_l, x_u, g_l, g_u, &fg_info, &s
	);

	// create an instance of the IpoptApplication
	using Ipopt::IpoptApplication;
	Ipopt::SmartPtr<IpoptApplication> app = new IpoptApplication();

	// Ipopt options
	app->Options()->SetIntegerValue("print_level", 0);
	app->Options()->SetStringValue("sb", "yes");
	app->Options()->SetNumericValue("tol", 1e-9);

	// Initialize the IpoptApplication and process the options
	Ipopt::ApplicationReturnStatus status = app->Initialize();
	ok    &= status == Ipopt::Solve_Succeeded;
	// Run the IpoptApplication
	status = app->OptimizeTNLP(cppad_nlp);
	ok    &= status == Ipopt::Solve_Succeeded;
	ok    &= s.status == cppad_ipopt_solution::success;

	for (i = 0; i < n_x; ++i)
		contvar[i] = s.x[i];

	bool feasible{ true };
	for (i = n_x; i < n; ++i) {
		if (s.x[i] >= 0.1) {
			feasible = false;
			break;
		}
	}

	double energy = 2.0 * s.x[0] + mu[0] * s.x[1] + mu[1] * s.x[2] + mu[2] * s.x[3] + mu[1] * s.x[4];
	return make_tuple(energy, feasible, ok);
}