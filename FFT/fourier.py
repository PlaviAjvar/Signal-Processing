import cmath

# pad the signal with zeroes until it's length becomes a power of two
def preprocess_signal(signal):
    num_samples = len(signal)
    # while number of samples is less than power of two
    # pad the signal with zeros at end to reach a length equal to power of two
    # the following condition is only zero for powers of two(consider bitwise representation)
    if num_samples & (num_samples - 1):
        num_samples = 1
        while num_samples < len(signal):
            num_samples *= 2
        signal += [0] * num_samples

# calculates fourier transform of given signal
# using decimation in time (DIT) FFT algorithm
def fft(signal):
    # base case only one sample, fourier transform is the same as the original signal
    if len(signal) == 1:
        return signal
    signal = preprocess_signal(signal)
    num_samples = len(signal)

    # find fft of signal formed by taking even and odd indices
    even_signal = [signal[2*time] for time in range(num_samples // 2)]
    odd_signal = [signal[2*time + 1] for time in range(num_samples // 2)]
    even_image = fft(even_signal)
    odd_image = fft(odd_signal)

    # combine the two images into one, obtaining the original signal
    image = [0] * num_samples
    base_scaler = cmath.exp(-2 * cmath.pi / num_samples)
    scaler = 1
    for frequency in range(num_samples // 2):
        image[frequency] = even_signal[frequency] + scaler * odd_signal[frequency]
        image[frequency + num_samples // 2] = even_signal[frequency] - scaler * odd_signal[frequency]
        scaler *= base_scaler
    return image

# calculates fourier transform of given signal
# using decimation in frequency (DIF) FFT algorithm
def dif_fft(signal):
    # base case only one sample, fourier transform is the same as the original signal
    if len(signal) == 1:
        return signal
    signal = preprocess_signal(signal)
    num_samples = len(signal)

    # separate first and second half
    first_half = signal[: num_samples // 2]
    second_half = signal[num_samples // 2 :]

    # find fft for even and odd frequencies
    even_image = fft([first_half[time] + second_half[time] for time in range(num_samples // 2)])
    # for efficiency sake we need to express the odd image iterativelly
    base_scaler = cmath.exp(-2 * cmath.pi / num_samples)
    odd_image = [0] * (num_samples // 2)
    scaler = 1
    for time in range(num_samples // 2):
        odd_image[time] = (first_half[time] - second_half[time]) * scaler
        scaler *= base_scaler

    # and finally combine them
    return [even_image[index] if index % 2 == 0 else odd_image[index] for index in range(num_samples)]

# computes the value of the logarithm base 2 of a number
def log_2(number):
    log_count = 0
    while number > 1:
        log_count += 1
        number /= 2
    return log_count

# reverse bits of binary number
def reverse(bin_num):
    bit_count = log_2(bin_num)
    reversed_num = 0
    for bit in range(bit_count):
        if bin_num & (1 << bit):
            reversed_num |= (1 << (bit_count-bit-1))
    return reversed_num

# reverse the bits of the indices in the signal
# e.g. a[6] -> a[3], because reversed(110) = 011 = 3
def bit_reverse(signal):
    num_samples = len(signal)
    new_signal = signal
    for index in range(num_samples):
        new_signal[reverse(index)] = signal[index]


# finds fourier transform using in-place fft algorithm
def in_place_fft(signal):
    signal = preprocess_signal(signal)
    num_samples = len(signal)

    # reverse the bits of the indices in the signal
    signal = bit_reverse(signal)
    # iterate over the table log(n) times, computing the values iterativelly
    # equivalent to starting from the base cases and building larger subproblems
    for exponent in range(log_2(num_samples)):
        slice = 2**exponent
        base_scaler = cmath.exp(-2*cmath.pi / slice)
        for start in range(0, num_samples, slice):
            scaler = 1
            for index in range(slice // 2):
                # do some in place stuff
                # completely equivalent to recursive implementation, just in place
                first = scaler * signal[start + index + slice // 2]
                second = signal[start + index]
                signal[start + index] = first + second
                signal[start + index + slice // 2] = first - second
                scaler *= base_scaler