from mpmath import mp
import numpy



def bisection_search(f, low:float, high:float):
    """
    A root finding method that does not rely on derivatives

    :param f: a function f: X -> R
    :param low: the lower bracket
    :param high: the upper limit bracket
    :return: the location of the root, e.g. f(mid) ~ 0
    """
    # flip high and low if out of order
    if f(high) < f(low):
        low, high = high, low

    # find mid point
    mid = .5 * (low + high)

    while True:

        # bracket up
        if f(mid) < 0:
            low = mid
        # braket down
        else:
            high = mid

        # update mid point
        mid = .5 * (high + low)

        # break if condition met
        if abs(high - low) < 10 ** (-(mp.dps / 2)):
            break

    return mid


def concave_max(f, low:float, high:float):
    """
    Forms a lambda for the approximate derivative and finds the root

    :param f: a function f: X -> R
    :param low: the lower bracket
    :param high: the upper limit bracket
    :return: the location of the root f'(mid) ~ 0
    """
    # create an approximate derivative expression
    scale = high - low

    h = mp.mpf('0.' + ''.join(['0' for i in range(int(mp.dps / 1.5))]) + '1') * scale
    df = lambda x: (f(x + h) - f(x - h)) / (2.0 * h)

    return bisection_search(df, low, high)

def chev_points(n:int, lower:float = -1, upper:float = 1):
    """
    Generates a set of chebychev points spaced in the range [lower, upper]
    :param n: number of points
    :param lower: lower limit
    :param upper: upper limit
    :return: a list of multipressison chebychev points that are in the range [lower, upper]
    """
    #generate chebeshev points on a range [-1, 1]
    index = numpy.arange(1, n+1)
    range_ = abs(upper - lower)
    return [(.5*(mp.cos((2*i-1)/(2*n)*mp.pi)+1))*range_ + lower for i in index]


def remez(func, n_degree:int, lower:float=-1, upper:float=1, max_iter:int = 10):
    """
    :param func: a function (or lambda) f: X -> R
    :param n_degree: the degree of the polynomial to approximate the function f
    :param lower: lower range of the approximation
    :param upper: upper range of the approximation
    :return: the polynomial coefficients, and an approximate maximum error associated with this approximation
    """
    # initialize the node points

    x_points = chev_points(n_degree + 2, lower, upper)

    A = mp.matrix(n_degree + 2)
    coeffs = numpy.zeros(n_degree + 2)

    # place in the E column
    mean_error = float('inf')

    for i in range(n_degree + 2):
        A[i, n_degree + 1] = (-1) ** (i + 1)

    for i in range(max_iter):

        # build the system
        vander = numpy.polynomial.chebyshev.chebvander(x_points, n_degree)

        for i in range(n_degree + 2):
            for j in range(n_degree + 1):
                A[i, j] = vander[i, j]

        b = mp.matrix([func(x) for x in x_points])
        l = mp.lu_solve(A, b)

        coeffs = l[:-1]

        # build the residual expression
        r_i = lambda x: (func(x) - numpy.polynomial.chebyshev.chebval(x, coeffs))

        interval_list = list(zip(x_points, x_points[1:]))
        #         interval_list = [[x_points[i], x_points[i+1]] for i in range(len(x_points)-1)]

        intervals = [upper]
        intervals.extend([bisection_search(r_i, *i) for i in interval_list])
        intervals.append(lower)

        extermum_interval = [[intervals[i], intervals[i + 1]] for i in range(len(intervals) - 1)]

        extremums = [concave_max(r_i, *i) for i in extermum_interval]

        extremums[0] = mp.mpf(upper)
        extremums[-1] = mp.mpf(lower)

        errors = [abs(r_i(i)) for i in extremums]
        mean_error = numpy.mean(errors)

        if numpy.max([abs(error - mean_error) for error in errors]) < 0.000001 * mean_error:
            break

        x_points = extremums

    return [float(i) for i in numpy.polynomial.chebyshev.cheb2poly(coeffs)], float(mean_error)

def c_code_gen(data_type, name, poly_coeffs, comments = None):
    method_string = f'{data_type} {name} ({data_type} x)' + '{\n'
    
    if comments is not None:
        method_string += '\t// ' + str(comments) + ' \n\n'
    
    data_type_converter = '' if data_type == 'double' else 'f'
    
    method_string += '\n'.join([f'\tconst {data_type} a_{i} = {str(val) + data_type_converter};' for i, val in enumerate(poly_coeffs)])
    
    horner = 'return a_0+'
    for i in range(len(poly_coeffs)-2):
        horner += f'x*(a_{i+1} +' 
    horner += f'x*a_{len(poly_coeffs)-1}' + ')'*(len(poly_coeffs)-2) + ';\n}'
    
    return method_string + '\n \t' + horner
