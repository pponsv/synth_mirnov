import numpy as np
from scipy.io import netcdf_file


def ispow2(n):
    return n != 0 and (n & (n - 1)) == 0


def make_coef_array(
    coefs, xm, xn, len_th, len_ph, deriv_order=0, deriv_dir=""
):
    """Make a 3D array of coefficients with the correct shape to be inverted."""
    assert deriv_dir in [
        "th",
        "ph",
        "",
    ], "Invalid derivative dir"  # "" means no derivative
    coef_array = np.zeros(
        (coefs.shape[0], len_th, len_ph), dtype=np.complex128
    )
    xm_mod = np.array(xm) % len_th
    xn_mod = -np.array(xn) % len_ph
    if deriv_order == 0:
        coefs_mod = coefs
    elif deriv_order == 1:
        if deriv_dir == "th":
            coefs_mod = 1j * xm * coefs
        elif deriv_dir == "ph":
            coefs_mod = -1j * xn * coefs
    np.add.at(coef_array, np.s_[:, xm_mod, xn_mod], coefs_mod)
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
        if (not ispow2(n_th)) or (not ispow2(n_ph)):
            raise ValueError("n_th, n_ph must be powers of 2")
        self.n_th = n_th
        self.n_ph = n_ph
        self.th = np.linspace(0, 2 * np.pi, n_th, endpoint=False)
        self.ph = np.linspace(0, 2 * np.pi, n_ph, endpoint=False)
        self.read_wout()
        self.get_coefs()
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
        self.iota_vmec = self._woutdata["iota_b"]
        self.iota = self.iota_vmec[self.s_idx]
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

        self.xs = self.rs * np.cos(self.phi_cyl)
        self.ys = self.rs * np.sin(self.phi_cyl)

        # Derivatives
        self.dr_ds = self.wrap_invert_fourier(
            np.gradient(
                self.rmnc, np.mean(np.diff(self.s)), axis=0, edge_order=2
            ),
            kind="cos",
        )
        self.dr_dth = self.wrap_invert_fourier(
            self.rmnc, kind="cos", deriv_order=1, deriv_dir="th"
        )
        self.dr_dph = self.wrap_invert_fourier(
            self.rmnc, kind="cos", deriv_order=1, deriv_dir="ph"
        )
        self.dph_ds = self.wrap_invert_fourier(
            np.gradient(
                self.pmns, np.mean(np.diff(self.s)), axis=0, edge_order=2
            ),
            kind="sin",
        )
        self.dph_dth = self.wrap_invert_fourier(
            self.pmns, kind="sin", deriv_order=1, deriv_dir="th"
        )
        self.dph_dph = 1 + self.wrap_invert_fourier(
            self.pmns, kind="sin", deriv_order=1, deriv_dir="ph"
        )
        self.dz_ds = self.wrap_invert_fourier(
            np.gradient(
                self.zmns, np.mean(np.diff(self.s)), axis=0, edge_order=2
            ),
            kind="sin",
        )
        self.dz_dth = self.wrap_invert_fourier(
            self.zmns, kind="sin", deriv_order=1, deriv_dir="th"
        )
        self.dz_dph = self.wrap_invert_fourier(
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
        self.e_sub_ph = np.array(
            [
                self.dr_dph * exph.real - self.rs * exph.imag * self.dph_dph,
                self.dr_dph * exph.imag + self.rs * exph.real * self.dph_dph,
                self.dz_dph,
            ]
        )

        self.sqrt_g = np.einsum(
            "ijkl,ijkl->jkl",
            self.e_sub_s,
            np.cross(self.e_sub_th, self.e_sub_ph, axisa=0, axisb=0, axisc=0),
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
