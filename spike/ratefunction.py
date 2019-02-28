
"""
Represent spike train  as a function of time to place an infinitely thin, infinitely tall window
to mark each spike. Use the Dirac delta function centered at t = t_hat:

1. delta(t) = 0 for t != 0
2. integral from a to b of delta(t)dt = 1 for a <= 0 <= b.
3. integral from a to b of f(t)*delta(t) = f(0) for any function f(t) and for a <= 0 <= b.


"""
