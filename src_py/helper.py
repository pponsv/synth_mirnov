from ..synth_mirnov import synthetic_mirnov as sm
import numpy as np


def init_synth_mirnov(booz):
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


def init_potentials(potentials, time):
    sm.init_pot(
        profiles=[pot.prof for pot in potentials],
        ms=[pot.m for pot in potentials],
        ns=[pot.n for pot in potentials],
        fs=[pot.freq for pot in potentials],
        time=time,
    )


def calc_potentials():
    sm.calc_potentials()


def run(time, coil_positions):
    sm.init_coils(coil_positions)
    return sm.run(len(coil_positions), len(time))
