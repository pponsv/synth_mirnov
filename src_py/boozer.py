import numpy as np
from scipy.io import netcdf_file


def make_coef_array(
    coefs, xm, xn, len_th, len_ph, deriv_order=0, deriv_dir=""
):
    """Make a 3D array of coefficients with the correct shape to be inverted."""
    coef_array = np.zeros(
        (coefs.shape[0], len_th, len_ph), dtype=np.complex128
    )
    if deriv_order == 0:
        coef_array[:, xm, -xn] = coefs
    elif deriv_order == 1:
        if deriv_dir == "th":
            coef_array[:, xm, -xn] = 1j * xm * coefs
        elif deriv_dir == "zt":
            coef_array[:, xm, -xn] = -1j * xn * coefs
        else:
            raise ValueError("Invalid deriv_dir")
    return coef_array


def invert_fourier(
    coefs_nm, xm, xn, len_th, len_ph, deriv_order=0, deriv_dir="", kind=""
):
    coef_array = make_coef_array(
        coefs_nm, xm, xn, len_th, len_ph, deriv_order, deriv_dir
    )
    if kind == "cos":
        return np.fft.ifft2(coef_array, norm="forward").real
    elif kind == "sin":
        return np.fft.ifft2(coef_array, norm="forward").imag
    else:
        return np.fft.ifft2(coef_array, norm="forward")


class Booz:
    def __init__(self, wout_path, n_th, n_ph):
        self.__wout_path = wout_path
        self._woutdata = {}
        self.n_th = n_th
        self.n_ph = n_ph
        self.th = np.linspace(0, 2 * np.pi, n_th, endpoint=False)
        self.ph = np.linspace(0, 2 * np.pi, n_ph, endpoint=False)
        self.read_wout()
        self.compute_booz()

    def read_wout(self):
        print("Reading boozer file from: ", self.__wout_path)
        with netcdf_file(self.__wout_path, "r") as wfile:
            for var in wfile.variables:
                self._woutdata[var] = wfile.variables[var].data.copy()

    def get_coefs(self):
        self.s_idx = self._woutdata["jlist"] - 1
        self.s_vmec = np.linspace(0, 1, self._woutdata["ns_b"])
        self.s = self.s_vmec[self.s_idx]
        self.iota = self._woutdata["iota_b"]
        self.phi_b_g = self._woutdata["phi_b"][-1]

        self.bmnc = self._woutdata["bmnc_b"]
        self.rmnc = self._woutdata["rmnc_b"]
        self.zmns = self._woutdata["zmns_b"]
        self.pmns = self._woutdata["pmns_b"]

        self.xm = self._woutdata["ixm_b"]
        self.xn = self._woutdata["ixn_b"]

    def compute_booz(self):
        self.rs = self.wrap_invert_fourier(self.rmnc, kind="cos")
        self.zs = self.wrap_invert_fourier(self.zmns, kind="sin")
        self.ps = self.wrap_invert_fourier(self.pmns, kind="sin")
        self.bs = self.wrap_invert_fourier(self.bmnc, kind="cos")

        self.phi_cyl = self.ph + self.ps
        # Derivatives
        self.dr_ds = self.wrap_invert_fourier(
            np.gradient(self.rmnc, np.mean(np.diff(self.s)), axis=0),
            kind="cos",
        )
        self.dr_dth = self.wrap_invert_fourier(
            self.rmnc, kind="cos", deriv_order=1, deriv_dir="th"
        )
        self.dr_dzt = self.wrap_invert_fourier(
            self.rmnc, kind="cos", deriv_order=1, deriv_dir="ph"
        )
        self.dph_ds = self.wrap_invert_fourier(
            np.gradient(self.pmns, np.mean(np.diff(self.s)), axis=0),
            kind="sin",
        )
        self.dph_dth = self.wrap_invert_fourier(
            self.pmns, kind="sin", deriv_order=1, deriv_dir="th"
        )
        self.dph_dzt = 1 - self.wrap_invert_fourier(
            self.pmns, kind="sin", deriv_order=1, deriv_dir="ph"
        )
        self.dz_ds = self.wrap_invert_fourier(
            np.gradient(self.zmns, np.mean(np.diff(self.s)), axis=0),
            kind="sin",
        )
        self.dz_dth = self.wrap_invert_fourier(
            self.zmns, kind="sin", deriv_order=1, deriv_dir="th"
        )
        self.dz_dzt = -self.wrap_invert_fourier(
            self.zmns, kind="sin", deriv_order=1, deriv_dir="ph"
        )

        # Unit vectors
        exph = np.exp(1j * self.phi_cyl)
        self.e_sub_s = np.array(
            [
                self.dr_ds * exph.real - self.rs * exph.imag * self.dph_ds,
                self.dr_ds * exph.imag + self.rs * exph.real * self.dph_ds,
                self.dz_ds,
            ]
        )
        self.e_sub_th = np.array(
            [
                self.dr_dth * exph.real - self.rs * exph.imag * self.dph_dth,
                self.dr_dth * exph.imag + self.rs * exph.real * self.dph_dth,
                self.dz_dth,
            ]
        )
        self.e_sub_zt = np.array(
            [
                self.dr_dzt * exph.real - self.rs * exph.imag * self.dph_dzt,
                self.dph_dth * exph.imag + self.rs * exph.real * self.dph_dzt,
                self.dph_dth,
            ]
        )

        self.sqrt_g = np.einsum(
            "ijkl,ijkl->jkl",
            self.e_sub_s,
            np.cross(self.e_sub_th, self.e_sub_zt, axisa=0, axisb=0, axisc=0),
        )

    def wrap_invert_fourier(self, coefs, deriv_order=0, deriv_dir="", kind=""):
        return invert_fourier(
            coefs,
            self.xm,
            self.xn,
            self.n_th,
            self.n_ph,
            deriv_order,
            deriv_dir,
            kind,
        )
