import matplotlib.pyplot as plt
import math
import itertools

# returns elements of matrix of custom size with zero entries
def zero_matrix(row_size, col_size):
    return [[0] * col_size for _ in range(row_size)]

class matrix:
    def __init__(self, row_size, col_size, elements):
        self.elements = elements
        self.row_size = row_size
        self.col_size = col_size

    # multiplies two matrices
    def multiply(self, other_matrix):
        if self.col_size != other_matrix.row_size:
            print((self.row_size, self.col_size))
            print((other_matrix.row_size, other_matrix.col_size))
            raise Exception("Product of matrices is undefined (size missmatch)")
        new_matrix = matrix(self.row_size, other_matrix.col_size, zero_matrix(self.row_size, other_matrix.col_size))
        for i in range(new_matrix.row_size):
            for j in range(new_matrix.col_size):
                for k in range(self.col_size):
                    new_matrix.elements[i][j] += self.elements[i][k] * other_matrix.elements[k][j]
        return new_matrix

    def transpose(self):
        self.elements = [[self.elements[row_idx][col_idx] for row_idx in range(self.row_size)] for col_idx in range(self.col_size)]
        self.row_size, self.col_size = self.col_size, self.row_size
        return self

    def scale(self, constant):
        self.elements = [[constant*entry for entry in row] for row in self.elements]
        return self

# finds Walsh matrix for a fixed approximation accuracy
def walsh_matrix(num_approx_fcn):
    matrix_size = 1
    walsh_cur = matrix(1, 1, [[1]])

    # iterate until we achieve desired size
    while matrix_size < num_approx_fcn:
        # recursive relation W(n) = [W(n-1) W(n-1); W(n-1) -W(n-1)]
        new_elements = [[walsh_cur.elements[row_idx % matrix_size][col_idx % matrix_size] for col_idx in range(matrix_size*2)] \
                        for row_idx in range(matrix_size*2)]
        for row_offset in range(matrix_size):
            for col_offset in range(matrix_size):
                new_elements[row_offset + matrix_size][col_offset + matrix_size] *= -1
        # construct new walsh matrix
        walsh_cur = matrix(matrix_size*2, matrix_size*2, new_elements)
        matrix_size *= 2

    return walsh_cur

# helper function for obtaining rows of haar matrix
def haar_row(power_2, offset, num_approx_fcn):
    row = [0] * num_approx_fcn
    slice = num_approx_fcn // power_2
    start = slice * offset
    middle = slice * offset + slice // 2
    end = slice * (offset + 1)
    row[start:middle] = [math.sqrt(power_2) for _ in range(middle-start)]
    row[middle:end] = [-math.sqrt(power_2) for _ in range(end-middle)]
    return row


# finds Haar matrix for a fixed approximation accuracy
def haar_matrix(num_approx_fcn):
    # construct Haar matrix with this size in mind
    haar_mat = matrix(num_approx_fcn, num_approx_fcn, zero_matrix(num_approx_fcn, num_approx_fcn))
    power_2 = 1
    haar_mat.elements[0] = [1] * num_approx_fcn
    for row_idx, row in enumerate(haar_mat.elements[1:], start=1):
        if row_idx >= power_2*2:
            power_2 *= 2
        offset = row_idx - power_2
        row[:] = haar_row(power_2, offset, num_approx_fcn)
    return haar_mat

# numericly evaluates definite integral
def integral(samples, step_size):
    return sum([sample * step_size for sample in samples])

# reconstruct samples from compressed square wave representation
def reconstruct_samples(approximation_matrix, num_samples):
    approximation = approximation_matrix[0]
    num_approx_samples = len(approximation)
    sample_ratio = num_samples // num_approx_samples
    reconstructed_samples = list(itertools.chain.from_iterable([[approx_sample] * sample_ratio for approx_sample in approximation]))
    if len(reconstructed_samples) < num_samples:
        reconstructed_samples += [0] * (num_samples - len(reconstructed_samples))
    return reconstructed_samples

# turn list into column vector
def vectorize(entries):
    return [[entry] for entry in entries]

# for error vector obtains the squared error
def squared(vector):
    return [entry*entry for entry in vector]

# reads piecewise constant function from input file
# generates approximation constants for Walsh/Haar functions
if __name__ == "__main__":
    file = open("input.txt","r")
    # interval length (basic period), number of samples and the samples themselves
    # number of samples should be power of two
    period, num_samples = file.readline().split()
    period = float(period)
    num_samples = int(num_samples)
    step_size = period / num_samples
    samples = list(map(float, file.readline().split()))
    # input number of Walsh/Haar functions to be used, assume power of two
    num_approx_fcn = int(file.readline())
    # find the corresponding matrices
    W = walsh_matrix(num_approx_fcn)
    H = haar_matrix(num_approx_fcn)

    # form colon vector of areas under the curve
    area_chunk = num_samples // num_approx_fcn
    areas = [integral(samples[start*area_chunk : (start+1)*area_chunk], step_size) for start in range(num_approx_fcn)]
    area_vector = matrix(len(areas), 1, vectorize(areas))

    # calculate Walsh and Haar coefficients
    walsh_coef = (W.multiply(area_vector)).scale(1/period)
    haar_coef = (H.multiply(area_vector)).scale(1/period)

    # construct aproximations with the given coefficients
    walsh_coef.transpose()
    walsh_approx = walsh_coef.multiply(W)
    haar_coef.transpose()
    haar_approx = haar_coef.multiply(H)
    walsh_samples = reconstruct_samples(walsh_approx.elements, num_samples)
    haar_samples = reconstruct_samples(haar_approx.elements, num_samples)

    # find value of distance metric between the functions
    # the metric is the square error metric
    error_walsh = [walsh_samples[sample] - samples[sample] for sample in range(num_samples)]
    error_haar = [haar_samples[sample] - samples[sample] for sample in range(num_samples)]
    metric_walsh = integral(squared(error_walsh), step_size)
    metric_haar = integral(squared(error_haar), step_size)

    # plot the functions
    time = [sample*period / num_samples for sample in range(num_samples)]

    plt.figure()
    plt.subplot(311)
    plt.title("Original function")
    plt.step(time, samples)
    plt.subplot(312)
    plt.title("Walsh approximation")
    plt.step(time, walsh_samples)
    plt.subplot(313)
    plt.title("Haar approximation")
    plt.step(time, haar_samples)

    # plot the errors
    plt.figure()
    plt.subplot(211)
    plt.title("Walsh approximation error")
    plt.step(time, error_walsh)
    plt.subplot(212)
    plt.title("Haar approximation error")
    plt.step(time, error_haar)

    # compare which of the approximations is better
    print("The mean square error for the Walsh approximation is", metric_walsh)
    print("The mean square error for the Haar approximation is", metric_haar)

    # arbitrary small constant
    epsilon = 1e-5

    if abs(metric_walsh - metric_haar) < epsilon:
        if abs(metric_haar) > epsilon:
            print("Both approximations are equally good acording to the mean square error metric.")
        else:
            print("Both approximations coincide with the original function")
    elif metric_walsh < metric_haar:
        print("The Walsh approximation is better acording to the mean square error metric.")
    else:
        print("The Haar approximation is better acording to the mean square error metric.")

    # show the diagrams
    plt.show()





