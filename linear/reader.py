import os

def categorize_binary_terms(bv, line, symbol, index):
    j = 0
    negative = False
    term_count = 0
    coeffs = ''
    end = False
    while not end:
        i = j
        while line[i] != '+' and line[i] != '-' and line[i] != '\n':
            i += 1

        term = line[j:i].split('*')
        if len(term) == 1 and term[0].strip()[:2] == 'b[':
            for k in range(len(term)):
                term[k] = term[k].strip()

            if negative:
                coeffs += ' -1.0,'
            else:
                coeffs += ' 1.0,'

            bv[int(term[0][2:-1])] += " make_tuple('" + symbol + "', " + index + ', ' + str(term_count) + '),'
            term_count += 1

        elif len(term) == 2:
            for k in range(len(term)):
                term[k] = term[k].strip()

            if negative:
                coeffs += ' -' + term[0] + ','
            else:
                coeffs += ' ' + term[0] + ','

            bv[int(term[1][2:-1])] += " make_tuple('" + symbol + "', " + index + ', ' + str(term_count) + '),'
            term_count += 1

        if line[i] == '+':
            negative = False
            j = i + 1
        elif line[i] == '-':
            negative = True
            j = i + 1
        else:
            end = True

    return coeffs[:-1] + ' '










model_obj_coeff = '\tmodel.obj_coeff = {'
model_ineq_coeff = '\tmodel.ineq_coeff = {'
model_ineq_rhs = '\tmodel.ineq_rhs = {'
model_eq_coeff = '\tmodel.eq_coeff = {'
model_eq_rhs = '\tmodel.eq_rhs = {'
model_cont_coeff = '\tmodel.cont_coeff = {'
model_cont_rhsconst = '\tmodel.cont_rhsconst = {'
model_bv = '\tmodel.bv = {'

with open('input') as f:
    reading_bounds = False
    #ineq_count
    #eq_count
    #cont_obj

    # read the number of variables
    line = f.readline()
    bivar_count, contvar_count = line.split()[-2:]
    bivar_count = int(bivar_count)
    contvar_count = int(contvar_count)

    bv = []
    for i in range(bivar_count):
        bv.append('')

    # read terms involving only binary variables in the objective function
    line = f.readline()
    coeffs = categorize_binary_terms(bv, line[2:], 'o', '0')
    model_obj_coeff += coeffs

    
    # read the remaining part of the objective function
    line = f.readline()
    cont_obj = line[2:].strip()

    cont_constraints = []
    cont_constraint_count = 0
    cont_constraint_bounds = []
    cont_eqconstraint_index = []
    ineq_count = 0
    eq_count = 0
    contvar_bounds = []
    for line in f:
        # read the bounds of continuous variables
        if reading_bounds:
            term = line.split()[-2:]
            if term[0].lower() == '-inf':
                term[0] = '-1.0e19'

            if term[1].lower() == 'inf':
                term[1] = '1.0e19'

            contvar_bounds.append(term)

        else:
            # read the constraints with only binary variables
            if line[:2] == 'bc':
                j = line.find('=')
                if line[j-1] == '<':
                    model_ineq_rhs += ' ' + line[j+1:].strip() + ','
                    coeffs = categorize_binary_terms(bv, line[2:j-1] + '\n', 'i', str(ineq_count))
                    model_ineq_coeff += ' {' + coeffs + '},'
                    ineq_count += 1
                else:
                    model_eq_rhs += ' ' + line[j+1:].strip() + ','
                    coeffs = categorize_binary_terms(bv, line[2:j] + '\n', 'e', str(eq_count))
                    model_eq_coeff += ' {' + coeffs + '},'
                    eq_count += 1

            # read the constraints with also continuous variables
            elif line[:2] == 'cc':
                j = line.find('=')
                index = str(cont_constraint_count)
                if line[j-1] == '<':
                    cont_constraints.append(line[2:j-1].strip())
                    cont_constraint_bounds.append(['-1.0e19', 'cont_rhs[' + index + ']'])
                else:
                    cont_constraints.append(line[2:j].strip())
                    cont_constraint_bounds.append(['cont_rhs[' + index + ']', 'cont_rhs[' + index + ']'])
                    cont_eqconstraint_index.append(cont_constraint_count)


                if 'b' in line[j+1:]:
                    i = line.index('b', j + 1)
                    while line[i] != '+' and line[i] != '-' and line[i] != '\n' and i > j + 1:
                        i -= 1

                    if line[j+1:i].strip() == '':
                        model_cont_rhsconst += ' 0.0,'
                    else:
                        model_cont_rhsconst += ' ' + line[j+1:i].strip() + ','

                    coeffs = categorize_binary_terms(bv, line[i:], 'c', index)
                    model_cont_coeff += ' {' + coeffs + '},'

                else:
                    model_cont_coeff += ' { },'
                    model_cont_rhsconst += ' ' + line[j+1:-1].strip() + ','

                cont_constraint_count += 1


            # change mode if it reaches the section of bounds on continuous variables
            elif line[:2] == 'cb':
                reading_bounds = True

# populate the information related to binary variables
for terms in bv:
    model_bv += ' {' + terms[:-1] + ' },'


# the top part in model.cpp
model = ''
model += '#include <vector>\n'
model += '#include <tuple>\n'
model += '#include <random>\n'
model += '#include "program.h"\n'
model += '\n'
model += 'using namespace std;\n'
model += '\n'


# the part in model.cpp that populates the Program class
model += 'void populate(Program& model) {\n'

if model_obj_coeff[-1] != '{':
    model += model_obj_coeff + '};\n'

if model_ineq_coeff[-1] != '{':
    model += model_ineq_coeff[:-1] + ' };\n'
    model += model_ineq_rhs[:-1] + ' };\n'

if model_eq_coeff[-1] != '{':
    model += model_eq_coeff[:-1] + ' };\n'
    model += model_eq_rhs[:-1] + ' };\n'

if model_cont_coeff[-1] != '{':
    model += model_cont_coeff[:-1] + ' };\n'
    model += model_cont_rhsconst[:-1] + ' };\n'

model += model_bv[:-1] + ' };\n'
model += '}\n'
model += '\n'


# the part in model.cpp that assigns the initial values of continuous variables
model += 'void set_initial_contvar(Walker& w, const Program& model) {\n'
contvar_initial = '\tw.contvar = {'
for term in contvar_bounds:
    contvar_initial += ' ' + str(min(max(0.0, float(term[0])), float(term[1]))) + ','
model += contvar_initial[:-1] + ' };\n'
model += '}'


# the part in model.cpp that determines how the walker should move (i.e., the topology of the unit hypercube)





subproblem = ''
subproblem += '#include <cppad/cppad.hpp>\n'
subproblem += '#include <cppad_ipopt_nlp.hpp>\n'
subproblem += '#include <vector>\n'
subproblem += '#include <algorithm>\n'
subproblem += '#include <tuple>\n'
subproblem += '\n'
subproblem += 'using namespace std;\n'
subproblem += '\n'
subproblem += 'namespace {\n'
subproblem += '\tusing namespace cppad_ipopt;\n'
subproblem += '\n'
subproblem += '\tclass FG_info : public cppad_ipopt_fg_info\n'
subproblem += '\t{\n'
subproblem += '\tprivate:\n'
subproblem += '\t\tbool retape_{ false };\n'
subproblem += '\t\tvector<bool> b;\n'
subproblem += '\t\tvector<double> mu;\n'
subproblem += '\tpublic:\n'
subproblem += '\t\t// derived class part of constructor\n'
subproblem += '\t\tFG_info(vector<bool> bits, vector<double> cont_penalty_coeff) : b(bits), mu(cont_penalty_coeff) {}\n'
subproblem += '\n'
subproblem += '\t\t// evaluation of the objective f(x), and constraints g(x)\n'
subproblem += '\t\t// using an algorithmic differentiation (AD) class\n'
subproblem += '\t\tADVector eval_r(size_t k, const ADVector& x_i){\n'
subproblem += '\t\t\tADVector fg(' + str(cont_constraint_count + 1) + ');\n'
subproblem += '\n'
subproblem += '\t\t\tvector<ADNumber> x;\n'


subproblem += '\t\t\tfor (size_t i = 0; i < ' + str(contvar_count + cont_constraint_count + len(cont_eqconstraint_index)) + '; ++i)\n'
subproblem += '\t\t\t\tx.push_back(x_i[i]);\n'
subproblem += '\n'

obj = cont_obj
for i in range(cont_constraint_count):
    obj += ' + mu[' + str(i) + '] * x[' + str(i + contvar_count) + ']'
for i in range(len(cont_eqconstraint_index)):
    obj += ' + mu[' + str(cont_eqconstraint_index[i]) + '] * x[' + str(i + contvar_count + cont_constraint_count) + ']'
subproblem += '\t\t\tfg[0] = ' + obj + ';\n'
subproblem += '\n'

t = 0
for i in range(cont_constraint_count):
    if t < len(cont_eqconstraint_index) and cont_eqconstraint_index[t] == i:
        subproblem += '\t\t\tfg[' + str(i + 1) + '] = ' + cont_constraints[i] + ' - x[' + str(i + contvar_count) + ']' + ' + x[' + str(t + contvar_count + cont_constraint_count) + '];\n'
        t += 1
    else:
        subproblem += '\t\t\tfg[' + str(i + 1) + '] = ' + cont_constraints[i] + ' - x[' + str(i + contvar_count) + '];\n'

subproblem += '\n'
subproblem += '\t\t\treturn fg;\n'
subproblem += '\t\t}\n'
subproblem += '\n'
subproblem += '\t\tbool retape(size_t k){\n'
subproblem += '\t\t\treturn retape_;\n'
subproblem += '\t\t}\n'
subproblem += '\t};\n'
subproblem += '}\n'
subproblem += '\n'
subproblem += 'tuple<double, bool, bool> solve_subproblem(const vector<bool>& b, const vector<double>& cont_rhs, const vector<double>& mu, vector<double>& contvar) {\n'
subproblem += '\tbool ok{ true };\n'
subproblem += '\tsize_t i;\n'
subproblem += '\tdouble temp;\n'
subproblem += '\n'
subproblem += '\n'
subproblem += '\t// number of original variables\n'
subproblem += '\tsize_t n_x = ' + str(contvar_count) + ';\n'
subproblem += '\t// number of constraints\n'
subproblem += '\tsize_t m = ' + str(cont_constraint_count) + ';\n'
subproblem += '\t// number of equality constraints\n'
subproblem += '\tsize_t me = ' + str(len(cont_eqconstraint_index)) + ';\n'
subproblem += '\t// number of original and auxiliary variables\n'
subproblem += '\tsize_t n = m + me + n_x;\n'
subproblem += '\n'
subproblem += '\n'
subproblem += '\t// initial value of variables\n'
subproblem += '\tNumberVector x(n);\n'
subproblem += '\tfor (i = 0; i < n_x; ++i)\n'
subproblem += '\t\tx[i] = contvar[i];\n'
subproblem += '\n'

t = 0
for i in range(cont_constraint_count):
    if t < len(cont_eqconstraint_index) and cont_eqconstraint_index[t] == i:
        subproblem += '\ttemp = ' + cont_constraints[i] + ' - cont_rhs[' + str(i) + '];\n'
        subproblem += '\tx[' + str(i + contvar_count) + '] = max(0.0, temp);\n'
        subproblem += '\tx[' + str(t + contvar_count + cont_constraint_count) + '] = max(0.0, -temp);\n'
        t += 1
    else:
        subproblem += '\tx[' + str(i + contvar_count) + '] = max(0.0, ' + cont_constraints[i] + ' - cont_rhs[' + str(i) + ']);\n'


subproblem += '\n'
subproblem += '\t// lower and upper limits for x\n'
subproblem += '\tNumberVector x_l(n);\n'
subproblem += '\tNumberVector x_u(n);\n'

for i in range(contvar_count):
    subproblem += '\tx_l[' + str(i) + '] = ' + contvar_bounds[i][0] + ';\t x_u[' + str(i) + '] = ' + contvar_bounds[i][1] + ';\n'


subproblem += '\tfor (i = n_x; i < n; ++i) {\n'
subproblem += '\t\tx_l[i] = 0.0;\tx_u[i] = 1.0e19;\n'
subproblem += '\t}\n'
subproblem += '\n'
subproblem += '\n'
subproblem += '\t// lower and upper limits for g\n'
subproblem += '\tNumberVector g_l(m);\n'
subproblem += '\tNumberVector g_u(m);\n'


for i in range(cont_constraint_count):
    subproblem += '\tg_l[' + str(i) + '] = ' + cont_constraint_bounds[i][0] + ';\tg_u[' + str(i) + '] = ' + cont_constraint_bounds[i][1] + ';\n'


subproblem += '\n'
subproblem += '\n'
subproblem += '\t// create the Ipopt interface\n'
subproblem += '\tFG_info fg_info(b, mu);\n'
subproblem += '\n'
subproblem += '\t// create the Ipopt interface\n'
subproblem += '\tcppad_ipopt_solution s;\n'
subproblem += '\tIpopt::SmartPtr<Ipopt::TNLP> cppad_nlp = new cppad_ipopt_nlp(\n'
subproblem += '\t\tn, m, x, x_l, x_u, g_l, g_u, &fg_info, &s\n'
subproblem += '\t);\n'
subproblem += '\n'
subproblem += '\t// create an instance of the IpoptApplication\n'
subproblem += '\tusing Ipopt::IpoptApplication;\n'
subproblem += '\tIpopt::SmartPtr<IpoptApplication> app = new IpoptApplication();\n'
subproblem += '\n'
subproblem += '\t// Ipopt options\n'
subproblem += '\tapp->Options()->SetIntegerValue(\"print_level\", 0);\n'
subproblem += '\tapp->Options()->SetStringValue(\"sb\", \"yes\");\n'
subproblem += '\tapp->Options()->SetNumericValue(\"tol\", 1e-9);\n'
subproblem += '\n'
subproblem += '\t// Initialize the IpoptApplication and process the options\n'
subproblem += '\tIpopt::ApplicationReturnStatus status = app->Initialize();\n'
subproblem += '\tok    &= status == Ipopt::Solve_Succeeded;\n'
subproblem += '\t// Run the IpoptApplication\n'
subproblem += '\tstatus = app->OptimizeTNLP(cppad_nlp);\n'
subproblem += '\tok    &= status == Ipopt::Solve_Succeeded;\n'
subproblem += '\tok    &= s.status == cppad_ipopt_solution::success;\n'
subproblem += '\n'
subproblem += '\tfor (i = 0; i < n_x; ++i)\n'
subproblem += '\t\tcontvar[i] = s.x[i];\n'
subproblem += '\n'
subproblem += '\tbool feasible{ true };\n'
subproblem += '\tfor (i = n_x; i < n; ++i) {\n'
subproblem += '\t\tif (s.x[i] >= 0.1) {\n'
subproblem += '\t\t\tfeasible = false;\n'
subproblem += '\t\t\tbreak;\n'
subproblem += '\t\t}\n'
subproblem += '\t}\n'
subproblem += '\n'
subproblem += '\tdouble energy = ' + obj.replace('x', 's.x') + ';\n'
subproblem += '\treturn make_tuple(energy, feasible, ok);\n'
subproblem += '}'


with open(os.getcwd() + '\\' + 'model.cpp', 'w') as f:
    f.write(model)

with open(os.getcwd() + '\\' + 'subproblem.cpp', 'w') as f:
    f.write(subproblem)