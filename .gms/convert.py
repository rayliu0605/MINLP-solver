from math import *

# inser a .gms file name here
filename = ''

bv = []
bv_rev = {}
bv_constraint = []

cv = []
cv_rev = {}
cv_bounds = []
cv_constraint = []

obj_bv = ''
obj_cv = ''
obj_factor = 1

# determine whether the mode is maximization or minimization
with open(filename) as f:
    for l in reversed(f.readlines()):
        if l[:5] == 'Solve':
            if 'max' in l:
                obj_factor = -1
            break

# a script that reads a .gms file and convert it to input format
with open(filename) as f:
    # read the files and convert it to the standard input format
    current_section = 0
    whole_line = ''
    for l in f:
        if l.strip() == '':
            current_section = 0
            whole_line = ''
        else:
            if l[:9] == 'Variables':
                whole_line = l.strip()
                current_section = 1

            elif l[:18] == 'Positive Variables':
                whole_line = l.strip()
                current_section = 2

            elif '..' in l:
                whole_line = l.strip()
                current_section = 3

            elif l.strip() == '* set non-default bounds':
                current_section = 4

            else:
                whole_line += l.strip()


            if ';' in whole_line:
                # replace sqr and power by pow(
                whole_line = whole_line.lower().replace('power(', 'pow(')
                while 'sqr(' in whole_line:
                    i = whole_line.index('sqr(')
                    j = i + 3
                    bracket_count = 1
                    while bracket_count > 0:
                        j += 1
                        if whole_line[j] == ')':
                            bracket_count -= 1
                        elif whole_line[j] == '(':
                            bracket_count += 1

                    whole_line = whole_line[:i] + 'pow(' + whole_line[i + 4:j] + ', 2)' + whole_line[j + 1:]

                # classify variables
                if current_section == 1:
                    var = whole_line[9:-1].strip().split(',')
                    for v in var:
                        if v != 'objvar':
                            if v[0] == 'b':
                                bv_rev[v] = len(bv)
                                bv.append(v)
                            elif v[0] == 'x':
                                cv_rev[v] = len(cv)
                                cv.append(v)
                            else:
                                print('error at the line:', whole_line)

                    for _ in range(len(cv)):
                        cv_bounds.append([-inf, inf])

                # modify lower bounds for positive continuous variables
                elif current_section == 2:
                    var = whole_line[18:-1].strip().split(',')
                    for v in var:
                        cv_bounds[cv_rev[v]][0] = 0.0

                # read objective function and constraints
                elif current_section == 3:
                    equation = whole_line[whole_line.index('..') + 2 : -1].strip()
                    label = 'bc'

                    # identify if the current equation is the objective function
                    if 'objvar' in equation:
                        label = 'o'
                        equation = equation.replace('objvar', 'o')

                    # replace variable names
                    j = 0
                    while 'b' in equation[j:]:
                        j = equation[j:].index('b') + j
                        i = j + 1
                        while equation[i].isdigit():
                            i += 1
                        
                        equation = equation[:j] + 'b[' + str(bv_rev[equation[j:i]]) + ']' + equation[i:]
                        j += 1

                    j = 0
                    while 'x' in equation[j:]:
                        if label != 'o':
                            label = 'cc'

                        j = equation[j:].index('x') + j
                        i = j + 1
                        while equation[i].isdigit():
                            i += 1

                        equation = equation[:j] + 'x[' + str(cv_rev[equation[j:i]]) + ']' + equation[i:]
                        j += 1

                    factor = 1
                    relation = ' <= '
                    if '=G=' in equation:
                        factor = -1
                    elif '=E=' in equation:
                        relation = ' = '


                    feqs_index = equation.index('=')
                    rhs = factor * float(equation[feqs_index + 3:])

                    # collect terms and signs of a equation
                    terms = []
                    if equation[0] == '+':
                        signs = [factor]
                        i = 1
                        j = 1
                    elif equation[0] == '-':
                        signs = [-factor]
                        i = 1
                        j = 1
                    else:
                        signs = [factor]
                        i = 0
                        j = 0

                    bracket_count = 0
                    while i < feqs_index:
                        if equation[i] == '(':
                            bracket_count += 1
                        elif equation[i] == ')':
                            bracket_count -= 1
                        elif equation[i] == '+' and not bracket_count:
                            terms.append(equation[j:i].strip())
                            signs.append(factor)
                            j = i + 1
                        elif equation[i] == '-' and not bracket_count:
                            terms.append(equation[j:i].strip())
                            signs.append(-factor)
                            j = i + 1

                        i += 1

                    terms.append(equation[j:feqs_index].strip())


                    # append objective function
                    if label == 'o':
                        obj_factor *= -1 * signs[terms.index('o')]
                        for k in range(len(terms)):
                            if terms[k] != 'o':
                                if 'x' not in terms[k]:
                                    if signs[k] * obj_factor == 1:
                                        obj_bv += ' + ' + terms[k]
                                    else:
                                        obj_bv += ' - ' + terms[k]
                                else:
                                    if signs[k] * obj_factor == 1:
                                        obj_cv += ' + ' + terms[k]
                                    else:
                                        obj_cv += ' - ' + terms[k]

                        if len(obj_bv) > 0:
                            if obj_bv[1] == '+':
                                obj_bv = obj_bv[3:]
                            else:
                                obj_bv = '-' + obj_bv[3:]

                        if len(obj_cv) > 0:
                            if obj_cv[1] == '+':
                                obj_cv = obj_cv[3:]
                            else:
                                obj_cv = '-' + obj_cv[3:]

                    # append a constraint involving only binary variables
                    elif label == 'bc':
                        if signs[0] == 1:
                            constraint = terms[0]
                        else:
                            constraint = '-' + terms[0]

                        for k in range(1, len(terms)):
                            if signs[k] == 1:
                                constraint += ' + ' + terms[k]
                            else:
                                constraint += ' - ' + terms[k]

                        bv_constraint.append('bc ' + constraint + relation + str(rhs))

                    # append a constraint involving continuous variables
                    else:
                        constraint_lhs = ''
                        constraint_rhs = ''
                        for k in range(len(terms)):
                            if 'x' not in terms[k]:
                                if signs[k] == 1:
                                    constraint_rhs += ' - ' + terms[k]
                                else:
                                    constraint_rhs += ' + ' + terms[k]
                            else:
                                if signs[k] == 1:
                                    constraint_lhs += ' + ' + terms[k]
                                else:
                                    constraint_lhs += ' - ' + terms[k]

                        if constraint_lhs[1] == '+':
                            constriant_lhs = constraint_lhs[3:]
                        else:
                            constraint_lhs = '-' + constraint_lhs[3:]

                        cv_constraint.append('cc ' + constriant_lhs + relation + str(rhs) + constraint_rhs)

                # modify other bounds in the bound section
                elif current_section == 4:
                    while ';' in whole_line:
                        sc_index = whole_line.index(';')
                        dot_index = whole_line.index('.')
                        cv_index = cv_rev[whole_line[:dot_index]]
                        rhs = float(whole_line[whole_line.index('=') + 1:sc_index])
                        if 'up' == whole_line[dot_index + 1 : dot_index + 3]:
                            cv_bounds[cv_index][1] = rhs
                        elif 'fx' == whole_line[dot_index + 1 : dot_index + 3]:
                            cv_bounds[cv_index][0] = rhs
                            cv_bounds[cv_index][1] = rhs
                        elif 'lo' == whole_line[dot_index + 1 : dot_index + 3]:
                            cv_bounds[cv_index][0] = rhs
                        else:
                            print('error at th eline:', whole_line)

                        whole_line = whole_line[sc_index + 1:].strip()

with open('input', 'w') as f:
    f.write('vc ' + str(len(bv)) + ' ' + str(len(cv)) + '\n')
    f.write('bo ' + obj_bv + '\n')
    f.write('co ' + obj_cv + '\n')
    for c in bv_constraint:
        f.write(c + '\n')

    for c in cv_constraint:
        f.write(c + '\n')

    f.write('cb\n')
    for b in cv_bounds:
        f.write(str(b[0]) + ' ' + str(b[1]) + '\n')