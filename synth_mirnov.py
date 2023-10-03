import numpy as np
import matplotlib.pyplot as plt
import vmec_library as vl
import libreria_calibracion as lc
from synth_mirnov import synthetic_mirnov as sm  # type: ignore
from tictoc import tic, toc

LOAD_BOOZ = True
SAVE_BOOZ = False
# LOAD_BOOZ = False
# SAVE_BOOZ = True

LEN_THETA = 64
LEN_PHI = 128

# sm.test_meshgrid()
# exit()


def gaussian_profile(s, s0, amp, sigma, phase=0):
    return amp * np.exp(1j * phase) * np.exp(-((s - s0) ** 2) / sigma)


def ramp_soft(x, y, x0, x1, k):
    return 0.5 * y * ((1 + np.tanh(k * (x - x0)) * np.tanh(k * (x1 - x))))


def ramp_hard(x, y, x0=0.0, x1=1.0, k=0.01):
    return np.nan_to_num((y * ((x - x0) ** k) * (x1 - x) ** k), nan=0.0)


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


torarr = lc.Mirnov_T_Array(None)
polarr = lc.Mirnov_P_Array(None)

# coil_positions = [coil.xyz.flatten() for coil in torarr] + [
#     coil.xyz.flatten() for coil in polarr
# ]
coil_positions = [coil.xyz.flatten() for coil in polarr]
coil_positions = np.array(coil_positions)[:]
print(coil_positions.shape)

# coil_positions = np.interp(
#     coil_positions, np.linspace(0, 1, coil_positions.shape[1]), np.linspace(0, 1)
# )
# print(coil_positions.shape)
# exit()
# coil_positions = coil_positions[:, 0]

if LOAD_BOOZ is False:
    booz = vl.Booz(
        # "./device/tjii/100_44_64/boozmn_100_44_64_0.0.nc",
        "/home/pedro/MEGA/00_doctorado/research/VMEC/TJ-II/100_44_64.0.0/boozmn_100_44_64_0.0.nc",
        theta=np.linspace(0, 2 * np.pi, LEN_THETA + 1)[:-1],
        phi=np.linspace(0, 2 * np.pi, LEN_PHI + 1)[:-1],
    )
    booz.get_vectors()
    if SAVE_BOOZ is True:
        import pickle

        with open(f"./booz_pickles/booz_TJII_{LEN_THETA}_{LEN_PHI}.pkl", "wb") as tfile:
            pickle.dump(booz, tfile)
elif LOAD_BOOZ is True:
    import pickle

    with open(f"./booz_pickles/booz_TJII_{LEN_THETA}_{LEN_PHI}.pkl", "rb") as tfile:
        booz = pickle.load(tfile)
else:
    print("UNBOUND BOOZ")
    booz = vl.Booz(
        "./device/tjii/100_44_64/boozmn_100_44_64_0.0.nc",
        theta=np.linspace(0, 2 * np.pi, 2 + 1)[:-1],
        phi=np.linspace(0, 2 * np.pi, 4 + 1)[:-1],
    )

time = np.linspace(0, 0.03, 40)
potentials = [
    potential(
        booz.s,
        m=5,
        n=8,
        freq=150,
        profile_function=gaussian_forced_zero_hard,
        # profile_function=gaussian_profile,
        amp=15,
        phase=np.angle(0),
        s0=0.4,
        sigma=0.05,
    ),
    # potential(
    #     booz.s,
    #     m=4,
    #     n=7,
    #     freq=150,
    #     profile_function=gaussian_profile,
    #     amp=15,
    #     phase=np.angle(0),
    #     s0=0.4,
    #     sigma=0.03,
    # ),
]

tic()
sm.init_booz(
    s_b=booz.s,
    th_b=booz.th,
    ph_b=booz.ph,
    b_mod_b=booz.vars["mod_b"],
    sqrt_g_b=booz.vars["sqrtg_vecs"],
    phi_b_g=booz.woutdata["phi_b"][-1],
    iota_b=booz.iota,
    x=booz.xyzs["xs"],
    y=booz.xyzs["ys"],
    z=booz.xyzs["zs"],
)
sm.init_basis(
    e_sub_s_b=np.moveaxis(booz.vecs["e_s"], 0, -1),
    e_sub_th_b=np.moveaxis(booz.vecs["e_th"], 0, -1),
    e_sub_ph_b=np.moveaxis(booz.vecs["e_ph"], 0, -1),
)
sm.init_pot(
    profiles=[pot.prof for pot in potentials],
    ms=[pot.m for pot in potentials],
    ns=[pot.n for pot in potentials],
    fs=[pot.freq for pot in potentials],
    time=time,
)
sm.init_coils(coil_positions)
t_init = toc("Initialization")

# print("Test magnetic field")
# sm.test_magnetic_field(121, 51, 41)
# exit()
sm.init_pot(
    profiles=[pot.prof for pot in potentials],
    ms=[pot.m for pot in potentials],
    ns=[pot.n for pot in potentials],
    fs=[pot.freq for pot in potentials],
    time=time,
)

tic()
db = sm.run(coil_positions.shape[0], len(time))
t_loop = toc("Main Loop")
print(f"Time elapsed\t{t_loop+t_init}\tTotal")

fig_all, axes_all = plt.subplots(
    5, 5, constrained_layout=True, sharex=True, sharey=True
)
axes_all = axes_all.flatten()
# for coil_idx, ax in enumerate(axes_all.flatten()):
for coil_idx in range(len(coil_positions)):
    axes_all[coil_idx].plot(time, db[0, coil_idx, :])
    axes_all[coil_idx].plot(time, db[1, coil_idx, :])
    axes_all[coil_idx].plot(time, db[2, coil_idx, :])
    axes_all[coil_idx].set(title=f"coil_{coil_idx}")

plt.figure()
for pot in potentials:
    plt.plot(np.abs(pot.prof), ls="-")
    plt.plot(np.real(pot.prof), ls="-.")
    plt.plot(np.imag(pot.prof), ls="--")

plt.figure()
plt.plot(time, db[0, 0, :])
plt.plot(time, db[1, 0, :])
plt.plot(time, db[2, 0, :])

plt.figure()
plt.plot(db[0, :, 0])
plt.plot(db[1, :, 0])
plt.plot(db[2, :, 0])


plt.figure()
sig_p = np.zeros((len(polarr.coils), len(time)))
for idx, coil in enumerate(polarr):
    sig_p[idx] = np.dot(coil.normal, db[:, idx])
    plt.plot(time, lc.normalize(sig_p[idx]) + 0.5 * idx, color="k", lw=1)


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
