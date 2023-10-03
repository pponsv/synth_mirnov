import numpy as np


def gaussian_profile(s, s0, amp, sigma, phase=0):
    return amp * np.exp(1j * phase) * np.exp(-((s - s0) ** 2) / sigma)


def ramp_hard(x, y, x0=0.0, x1=1.0, k=0.01):
    return np.nan_to_num((y * ((x - x0) ** k) * (x1 - x) ** k), nan=0.0)


def ramp_soft(x, y, x0, x1, k):
    return 0.5 * y * ((1 + np.tanh(k * (x - x0)) * np.tanh(k * (x1 - x))))


def gaussian_forced_zero_soft(s, s0, amp, sigma, phase=0):
    return ramp_soft(s, gaussian_profile(s, s0, amp, sigma, phase), 0.05, 0.95, 50)


def gaussian_forced_zero_hard(s, s0, amp, sigma, phase=0):
    return ramp_hard(s, gaussian_profile(s, s0, amp, sigma, phase), 0.05, 0.95, 0.01)


class potential:
    def __init__(
        self, s, m, n, freq, profile_function=gaussian_forced_zero_soft, **kwargs
    ):
        self.s = s
        self.m = m
        self.n = n
        self.freq = freq
        self.prof = profile_function(s=s, **kwargs)
