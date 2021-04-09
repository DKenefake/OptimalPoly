# OptimalPoly
Stable implementation of the Remez Algorithm using multi-precision arithmetic. Also does C/C++ code generation!

## Example of use

This is all that is needed to create a 4th order optimal polynomial approximation of e^x on [0,1].

```python
function = lambda x: mp.exp(x)
poly_coeffs, max_error = remez(f, 4, 0, 1)
```

We can inspect how well we approximated it visually with the following code snippet.  

```python
import matplotlib.pyplot as plt
x = numpy.linspace(0, 1, 50)

y_exact = numpy.polyval(poly_coeffs[::-1], x)
y_approx = numpy.array([function(x_i) for x_i in x])
plt.plot(x, y_exact)
plt.plot(x, y_approx, 'x')
plt.title(r'$f(x)$ v. $P^*_{4}(x)$')
```
![image](https://github.com/DKenefake/OptimalPoly/blob/main/assets/compare.png)

And we can also look at the distinctive equioscillation of the optimal polynomial result with the following.

```python
import matplotlib.pyplot as plt
x = numpy.linspace(0, 1, 500)

y_exact = numpy.polyval(poly_coeffs[::-1], x)
y_approx = numpy.array([function(x_i) for x_i in x])
plt.plot(x, y_exact - y_approx)
plt.title('$f(x) - P^*_{4}(x)$')
```

![image](https://github.com/DKenefake/OptimalPoly/blob/main/assets/Equioscillation.png)

To generate the C code, you need to specify the data type of either `float` or `double`. This uses horner's method to speed up the polynomial evaluation, and depending on the compiler flags used, this should compile to 4 `fmad` operations. 

```python
c_code_gen(data_type = 'float', name = 'exp_approx', poly_coeffs = poly_coeffs, comments = f'x in [0, 1]')
```

```C
float exp_approx (float x){
	// x in [0, 1] 

	float a_0 = 1.0000271624188655f;
	float a_1 = 0.9986854006378604f;
	float a_2 = 0.510139460205787f;
	float a_3 = 0.13969814854689722f;
	float a_4 = 0.0697044942307694f;
 	return a_0+x*(a_1 +x*(a_2 +x*(a_3 +x*a_4)));
}
```

## Dependancies

* numpy
* mpmath

