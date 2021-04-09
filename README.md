# OptimalPoly
Stable implementation of the Remez Algorithm using multiprecision arithmetic.


# Example of use

This is all that is needed to create a 4th order optimal polynomial approximation of e^x on [-1,1]

```python
function = lambda x: mp.exp(x)
poly_coeffs, max_error = remez(f, 4, 0, 1)
```

We can inspect how well we approximated it visually with the following code snippit. 

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


And we can also look at the distinctive equiocillation of the optimal polynomail result with the following.

```python
import matplotlib.pyplot as plt
x = numpy.linspace(0, 1, 500)

y_exact = numpy.polyval(poly_coeffs[::-1], x)
y_approx = numpy.array([function(x_i) for x_i in x])
plt.plot(x, y_exact - y_approx)
plt.title('$f(x) - P^*_{4}(x)$')
```

![image](https://github.com/DKenefake/OptimalPoly/blob/main/assets/Equioscillation.png)

# Dependancies

* numpy
* mpmath

