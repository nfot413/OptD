import numpy as np


class TunableMZI:
    """代表一个可调谐的马赫-曾德尔干涉仪 (MZI)。"""
    def __init__(self, theta):
        j = 1j
        self.theta = theta
        coupler_50_50 = 0.5 * np.array([[-1+j, 1+j], [1+j, -1+j]])
        phase_matrix = np.array([[np.exp(-j * self.theta), 0], [0, 1]])
        self._transfer_matrix = coupler_50_50 @ phase_matrix @ coupler_50_50

    def get_transfer_matrix(self):
        return self._transfer_matrix

class PhaseShifter:
    """代表一个应用于两个臂的并行相移器。"""
    def __init__(self, phi_t, phi_b):
        j = 1j
        self._transfer_matrix = np.array([[np.exp(-j * phi_t), 0],
                                          [0, np.exp(-j * phi_b)]])

    def get_transfer_matrix(self):
        return self._transfer_matrix

class MicroRingResonator:
    """代表单个微环谐振器 (MRR)。"""
    def __init__(self, t, k, phi_offset=np.pi):
        self.j = 1j
        self.t = t
        self.k = k
        self.phi_offset = phi_offset

    def update_coupling_coefficient(self, k):
        self.k = k
        
    def get_response(self, w):
        numerator = np.sqrt(1 - self.k) - self.t**2 * np.exp(-self.j * (2*w + self.phi_offset))
        denominator = 1 - self.t**2 * np.sqrt(1 - self.k) * np.exp(-self.j * (2*w + self.phi_offset))
        return numerator / denominator

class DelayLine:
    """代表一个简单的光延迟线。"""
    def __init__(self, t, delay, phi_c):
        self.j = 1j
        self.t = t
        self.delay = delay
        self.phi_c = phi_c
        
    def get_response(self, w):
        return self.t * np.exp(-self.j * w * self.delay - self.j * self.phi_c)
