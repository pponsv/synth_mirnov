import numpy as np
import matplotlib.pyplot as plt
import vmec_library as vl
import libreria_calibracion as lc
from synth_mirnov import synthetic_mirnov as sm  # type: ignore
from tictoc import tic, toc

LOAD_BOOZ = True
SAVE_BOOZ = False

sm.test_meshgrid()
# exit()


def gaussian_profile(s, s0, amp, sigma):
    return amp * np.exp(-((s - s0) ** 2) / sigma)


class potential:
    def __init__(self, s, m, n, freq, profile_function=gaussian_profile, **kwargs):
        self.s = s
        self.m = m
        self.n = n
        self.freq = freq
        self.prof = profile_function(s=s, **kwargs)


torarr = lc.Mirnov_T_Array(None)
polarr = lc.Mirnov_P_Array(None)

coil_positions = [coil.xyz.flatten() for coil in torarr] + [
    coil.xyz.flatten() for coil in polarr
]
coil_positions = np.array(coil_positions).T
# coil_positions = coil_positions[:, 0]
print(coil_positions)

if LOAD_BOOZ is False:
    booz = vl.Booz(
        "./device/tjii/100_44_64/boozmn_100_44_64_0.0.nc",
        theta=np.linspace(0, 2 * np.pi, 2 + 1)[:-1],
        phi=np.linspace(0, 2 * np.pi, 3 + 1)[:-1],
    )
    booz.get_vectors()
    if SAVE_BOOZ is True:
        import pickle

        with open("./booz_pickles/booz_TJII.pkl", "wb") as tfile:
            pickle.dump(booz, tfile)
elif LOAD_BOOZ is True:
    import pickle

    with open("./booz_pickles/booz_TJII.pkl", "rb") as tfile:
        booz = pickle.load(tfile)
else:
    print("UNBOUND BOOZ")
    booz = vl.Booz(
        "./device/tjii/100_44_64/boozmn_100_44_64_0.0.nc",
        theta=np.linspace(0, 2 * np.pi, 2 + 1)[:-1],
        phi=np.linspace(0, 2 * np.pi, 4 + 1)[:-1],
    )

time = np.linspace(0, 0.1, 11)
potentials = [
    potential(
        booz.s,
        m=5,
        n=8,
        freq=150,
        profile_function=gaussian_profile,
        amp=15,
        s0=0.3,
        sigma=0.03,
    ),
    # potential(
    #     booz.s,
    #     m=4,
    #     n=7,
    #     freq=320,
    #     profile_function=gaussian_profile,
    #     amp=5 + 0.1j,
    #     s0=0.4,
    #     sigma=0.03,
    # ),
]


sm.init_booz(
    s_b=booz.s,
    th_b=booz.th,
    ph_b=booz.ph,
    b_mod_b=booz.vars["mod_b"],
    sqrt_g_b=booz.vars["sqrtg_vecs"],
    phi_b_g=booz.woutdata["phip_b"][-1],
    iota_b=booz.iota,
    x=booz.xyzs["xs"],
    y=booz.xyzs["ys"],
    z=booz.xyzs["zs"],
)
sm.init_basis(
    e_sub_s_b=booz.vecs["e_s"],
    e_sub_th_b=booz.vecs["e_th"],
    e_sub_ph_b=booz.vecs["e_ph"],
)
sm.init_pot(
    [pot.prof for pot in potentials],
    [pot.m for pot in potentials],
    [pot.n for pot in potentials],
    [pot.freq for pot in potentials],
    time=time,
)
sm.init_coils(coil_positions)

tic()
print(coil_positions.shape[1])
db = sm.run(coil_positions.shape[1], len(time))
toc()
plt.plot(time, db[0, 0, :])
plt.plot(time, db[0, 1, :])
plt.plot(time, db[0, 2, :])
plt.show()
# ths, phs = np.meshgrid(booz.th, booz.ph)
# fft_test = np.sin(3 * ths - 2 * phs).T
# print(fft_test.shape)
# print(len(booz.th), len(booz.ph))
# out = sm.test_fft_main(fft_test)
# nout = np.fft.fft2(fft_test)

# print(np.isclose(out, nout))
# print(np.all(np.isclose(out, nout)))
# print(out / nout)
# fig, ax = plt.subplots(3, 2)
# ax[0, 0].pcolor(np.abs(out))
# ax[1, 0].pcolor(np.real(out))
# ax[2, 0].pcolor(np.imag(out))
# ax[0, 1].pcolor(np.abs(nout))
# ax[1, 1].pcolor(np.real(nout))
# ax[2, 1].pcolor(np.imag(nout))

# plt.figure()
# plt.pcolor(np.abs(out - nout))
# plt.show()
