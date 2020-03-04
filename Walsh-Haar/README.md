# Walsh and Haar approximations

Input function samples and approximate the function as a linear
combination of Walsh or Haar functions. These are in the form of
square waves, and therefore go well with sampling.\
In the file input.txt, in the same directory as the project,
input the desired basic interval (the domain of the function or
a period) and the number of samples of the functions in the first
row. In the second row, input the values of the samples of the
functions. In the third row put the desired number of approximating
Walsh/Haar functions.\
For implementation reasons, the number of samples and number of
approximating functions need to be powers of two, and the latter
has to be less than or equal to the former. The code uses the
matplotlib package to plot the functions.

