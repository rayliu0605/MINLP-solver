from math import *

# inser a .gms file name here
filename = ''

# read the corresponding solution file
b = []
x = []
with open('solution.csv') as f:
    f.readline()
    for l in f:
        terms = l.strip().split(',')
        if len(terms) == 2:
            if 'b' in terms[0]:
                b.append(float(terms[1]))
            elif 'x' in terms[0]:
                x.append(float(terms[1]))


bv = []
bv_rev = {}

cv = []
cv_rev = {}
cv_bounds = []

with open(filename) as f:
    # read the files and convert it to the standard input form
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

                #check constraints
                elif current_section == 3:
                    equation = whole_line[whole_line.index('..') + 2 : -1].strip()
                    if 'objvar' not in equation:
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
                            j = equation[j:].index('x') + j
                            i = j + 1
                            while equation[i].isdigit():
                                i += 1

                            equation = equation[:j] + 'x[' + str(cv_rev[equation[j:i]]) + ']' + equation[i:]
                            j += 1

                        feqs_index = equation.index('=')
                        feasible = False
                        if equation[feqs_index + 1] == 'L':
                            exec('feasible = ' + equation[:feqs_index] + '<=' + equation[feqs_index + 3:] + ' + 0.1')
                        elif equation[feqs_index + 1] == 'G':
                            exec('feasible = ' + equation[:feqs_index] + '>=' + equation[feqs_index + 3:] + ' - 0.1')
                        else:
                            exec('feasible = abs(' + equation[:feqs_index] + ' - ' + equation[feqs_index + 3:] + ') < 0.1')

                        if not feasible:
                            exec('print(' + equation[:feqs_index] + ', "\t" + whole_line + "\t" + equation)')

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

# check bounds
for i in range(len(cv_bounds)):
    if not (x[i] >= cv_bounds[i][0] - 0.1 and x[i] <= cv_bounds[i][1] + 0.1):
        print(i, x[i], cv_bounds[i])